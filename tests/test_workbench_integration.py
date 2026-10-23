"""Real-resource tests for complete OrthoFinder research inputs and navigation."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from streamlit.testing.v1 import AppTest

from orthofinder_interrogation_app import app, launcher, workbench_page
from orthofinder_interrogation_app.models import GroupKey
from orthofinder_interrogation_app.queries import OrthoFinderQueryService
from orthofinder_interrogation_app.research.cli import main as workbench_main
from orthofinder_interrogation_app.research.common import ResearchConfig, ResearchError, file_sha256
from orthofinder_interrogation_app.research.family_inputs import import_family
from orthofinder_interrogation_app.research.jobs import artifact_path, get_job, submit_job
from orthofinder_interrogation_app.research.runner import run_next_job
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_interrogation_app.workbench_adapter import (
    select_group_inputs,
    selected_group_fasta,
)


@pytest.fixture
def workbench_service(application_resource: Path) -> OrthoFinderQueryService:
    """Open the normal test application resource through its public read-only API."""
    return OrthoFinderQueryService(resource=open_resource(path=application_resource))


@pytest.fixture
def workbench_sidecar(tmp_path: Path) -> Path:
    """Create complete protein sequences using the established motif-sidecar schema."""
    path = tmp_path / "sequences.parquet"
    pq.write_table(
        table=pa.table(
            {
                "internal_id": ["0_0", "0_1", "1_0"],
                "species_label": ["Species_A", "Species_A", "Species_B"],
                "member_id": ["alpha_1", "alpha_2", "beta_1"],
                "sequence": ["acde*", "ACDF", "ACDG"],
            }
        ),
        where=path,
    )
    return path


def test_group_selection_is_complete_scoped_and_size_ranked(workbench_service):
    """A top-N HOG selection preserves the exact hierarchy, species and full membership."""
    frame = select_group_inputs(
        service=workbench_service,
        group_type="HOG",
        hierarchy_node="N0",
        maximum_groups=1,
    )
    assert len(frame) == 3
    assert set(frame.input_id) == {"alpha_1", "alpha_2", "beta_1"}
    assert set(frame.source_group_id) == {"N0.HOG1"}
    assert set(frame.source_hierarchy_node) == {"N0"}
    assert set(frame.source_run_id) == {workbench_service.resource.run_id}
    assert set(frame.selection_rank) == {1}
    explicit = select_group_inputs(
        service=workbench_service,
        group_type="LEGACY_ORTHOGROUP",
        hierarchy_node="",
        group_ids=["OG4"],
    )
    assert set(explicit.source_hierarchy_node) == {""}
    assert len(explicit) == 1
    filtered = select_group_inputs(
        service=workbench_service,
        group_type="HOG",
        hierarchy_node="N0",
        included_species=["Species_B"],
        require_all_species=True,
        member_id_contains="beta_1",
    )
    assert len(filtered) == 3
    assert select_group_inputs(
        service=workbench_service,
        group_type="HOG",
        hierarchy_node="N0",
        maximum_groups=1,
        member_id_contains="does_not_exist",
    ).empty


@pytest.mark.parametrize(
    "changes",
    [
        {"group_ids": "N0.HOG0000001"},
        {"group_ids": [None]},
        {"group_ids": [""]},
        {"maximum_groups": 0},
        {"maximum_groups": True},
        {"maximum_inputs": 0},
        {"maximum_inputs": 2},
        {"maximum_groups": 1, "group_ids": ["A", "B"]},
    ],
)
def test_group_selection_refuses_partial_or_malformed_requests(workbench_service, changes):
    """Oversized selections fail rather than becoming unlabelled partial families."""
    arguments = {"service": workbench_service, "group_type": "HOG", "hierarchy_node": "N0"}
    with pytest.raises(ResearchError):
        select_group_inputs(**{**arguments, **changes})


def test_group_selection_checks_resource_membership_counts(workbench_service, monkeypatch):
    """Missing stored members cannot pass as a complete group selection."""
    monkeypatch.setattr(workbench_service, "get_group_members", lambda **kwargs: ())
    with pytest.raises(ResearchError, match="count disagrees"):
        select_group_inputs(service=workbench_service, group_type="HOG", hierarchy_node="N0")


def test_selected_fasta_and_ingestion_preserve_exact_context(
    workbench_service,
    workbench_sidecar,
    tmp_path,
):
    """Real sequence joins feed a successful queue job without altering the source resource."""
    source_hash = file_sha256(path=workbench_service.resource.database_path)
    inputs = select_group_inputs(
        service=workbench_service,
        group_type="HOG",
        hierarchy_node="N0",
        maximum_groups=1,
    )
    fasta = selected_group_fasta(
        service=workbench_service,
        inputs=inputs,
        sidecar=workbench_sidecar,
    )
    assert fasta == ">alpha_1\nACDE\n>alpha_2\nACDF\n>beta_1\nACDG\n"
    config = ResearchConfig(workspace=tmp_path / "work", threads=1)
    job = submit_job(
        workspace=config.workspace,
        operation="import_fasta",
        fasta=fasta,
        parameters={
            "source_context": inputs.to_dict(orient="records"),
            "source_label": "Test OrthoFinder run",
        },
    )
    assert run_next_job(config=config, job_id=job) == job
    assert get_job(workspace=config.workspace, job_id=job)["status"] == "COMPLETE"
    output = artifact_path(workspace=config.workspace, job_id=job, artifact="family.fasta")
    assert output.read_text() == fasta
    context = artifact_path(workspace=config.workspace, job_id=job, artifact="source_context.tsv")
    assert len(pd.read_csv(filepath_or_buffer=context, sep="\t")) == 3
    assert file_sha256(path=workbench_service.resource.database_path) == source_hash
    following = submit_job(
        workspace=config.workspace,
        operation="align",
        parameters={
            "source_job": job,
            "source_artifact": "family.fasta",
        },
    )
    assert get_job(workspace=config.workspace, job_id=following)["status"] == "QUEUED"


def test_sequence_transfer_refuses_missing_or_ambiguous_proteins(
    workbench_service,
    workbench_sidecar,
    tmp_path,
):
    """FASTA transfer requires one exact sequence per selected member and species."""
    inputs = select_group_inputs(
        service=workbench_service,
        group_type="HOG",
        hierarchy_node="N0",
        maximum_groups=1,
    )
    with pytest.raises(ResearchError, match="non-empty"):
        selected_group_fasta(service=workbench_service, inputs=pd.DataFrame())
    with pytest.raises(ResearchError, match="unavailable"):
        selected_group_fasta(service=workbench_service, inputs=inputs)
    duplicate = pd.concat(objs=[inputs, inputs.iloc[[0]].assign(species_label="Other")])
    with pytest.raises(ResearchError, match="repeat across species"):
        selected_group_fasta(service=workbench_service, inputs=duplicate, sidecar=workbench_sidecar)
    missing = tmp_path / "missing.parquet"
    pq.write_table(
        table=pq.read_table(source=workbench_sidecar).slice(offset=0, length=1), where=missing
    )
    with pytest.raises(ResearchError, match="Missing or duplicated"):
        selected_group_fasta(service=workbench_service, inputs=inputs, sidecar=missing)
    duplicated = tmp_path / "duplicated.parquet"
    table = pq.read_table(source=workbench_sidecar)
    pq.write_table(
        table=pa.concat_tables(tables=[table, table.slice(offset=0, length=1)]), where=duplicated
    )
    with pytest.raises(ResearchError, match="Missing or duplicated"):
        selected_group_fasta(service=workbench_service, inputs=inputs, sidecar=duplicated)


@pytest.mark.parametrize("context", ["wrong", ["wrong"], [{}], [{"input_id": "absent"}]])
def test_family_import_refuses_false_membership_context(tmp_path, context):
    """Imported FASTAs cannot acquire source-group claims for different proteins."""
    source = tmp_path / "input.fa"
    source.write_text(">a\nACD\n", encoding="utf-8")
    with pytest.raises(ResearchError):
        import_family(
            config=ResearchConfig(),
            source=source,
            output=tmp_path,
            parameters={
                "source_context": context,
            },
        )
    assert not (tmp_path / "family.fasta").exists()


def test_family_import_plain_fasta_and_overwrite_guard(tmp_path):
    """Plain families have honest local scope and cannot overwrite an earlier import."""
    source = tmp_path / "input.fa"
    source.write_text(">a\nACD\n", encoding="utf-8")
    artifacts = import_family(
        config=ResearchConfig(),
        source=source,
        output=tmp_path,
        parameters={},
    )
    assert "source_context.tsv" not in artifacts
    assert "SELECTED_LOCAL_PROTEINS" in artifacts["family_summary.tsv"].read_text()
    with pytest.raises(ResearchError, match="overwrite"):
        import_family(config=ResearchConfig(), source=source, output=tmp_path, parameters={})


def test_integrated_ui_and_current_group_handoff(
    workbench_service,
    workbench_sidecar,
    tmp_path,
    monkeypatch,
):
    """The actual app opens group inputs, queues local proteins and exposes every research page."""
    profile = tmp_path / "research.toml"
    workspace = tmp_path / "work"
    profile.write_text(f'workspace="{workspace}"\n', encoding="utf-8")
    monkeypatch.setenv(
        launcher.RESOURCE_ENVIRONMENT_VARIABLE, str(workbench_service.resource.resource_path)
    )
    monkeypatch.setenv(launcher.RESEARCH_ENVIRONMENT_VARIABLE, str(profile))
    monkeypatch.setenv(launcher.TERMINAL_MOTIF_ENVIRONMENT_VARIABLE, str(workbench_sidecar))
    test = AppTest.from_file(script_path=str(Path(app.__file__)), default_timeout=30).run()
    test.sidebar.radio[0].set_value("Cluster explorer").run()
    shortcut = next(
        item
        for item in test.button
        if item.label == "Open this group in protein research workbench"
    )
    shortcut.click().run()
    assert not test.exception
    assert test.sidebar.radio[0].value == "Protein research workbench"
    assert all(
        item.proto.help
        for controls in (
            test.radio,
            test.selectbox,
            test.text_area,
            test.text_input,
            test.number_input,
            test.checkbox,
            test.multiselect,
        )
        for item in controls
    )
    assert (
        next(item for item in test.number_input if item.label == "Largest local groups").value == 1
    )
    next(
        item for item in test.button if item.label == "Prepare complete group inputs"
    ).click().run()
    assert not test.exception and test.dataframe
    next(
        item for item in test.button if item.label == "Import selected local protein sequences"
    ).click().run()
    assert not test.exception and test.success
    config = ResearchConfig(workspace=workspace)
    assert run_next_job(config=config)
    next(item for item in test.button if item.label == "Track this job").click().run()
    assert test.radio(key="research_section").value == "Jobs"
    next(item for item in test.button if item.label == "Analyse this family").click().run()
    assert test.radio(key="research_section").value == "Family analysis"
    assert test.selectbox(key="analyse_source").value[1] == "family.fasta"
    test.radio(key="research_section").set_value("Start here").run()
    next(
        item for item in test.button if item.label == "Start from OrthoFinder groups (recommended)"
    ).click().run()
    assert test.radio(key="research_section").value == "OrthoFinder groups"
    for section in (
        "Discover",
        "Family analysis",
        "Alignments and trees",
        "Models",
        "Annotations",
        "Jobs",
        "Setup",
        "Help & glossary",
    ):
        next(item for item in test.radio if item.label == "Research section").set_value(
            section
        ).run()
        assert not test.exception, section


def test_navigation_callback_uses_exact_composite_group(monkeypatch):
    """Switching pages retains a run-scoped group key without substituting a global centroid."""
    state = {}
    monkeypatch.setattr(workbench_page.st, "session_state", state)
    key = GroupKey(run_id="run", group_type="HOG", hierarchy_node="N0", group_id="local")
    workbench_page.open_group_workbench(key=key)
    assert state["orthofinder_active_group"]["group_id"] == "local"
    assert state["app_page"] == "Protein research workbench"


def test_launcher_passes_valid_profile_and_rejects_resource_workspace(
    workbench_service,
    tmp_path,
    monkeypatch,
):
    """Named research configuration reaches Streamlit only when its workspace is external."""
    profile = tmp_path / "profile.toml"
    profile.write_text(f'workspace="{tmp_path / "research"}"\n', encoding="utf-8")
    run = Mock(return_value=Mock(returncode=0))
    monkeypatch.setattr(launcher.subprocess, "run", run)
    monkeypatch.setattr(launcher, "_port_available", lambda **kwargs: True)
    arguments = [
        "--resource-dir",
        str(workbench_service.resource.resource_path),
        "--research-config",
        str(profile),
    ]
    assert launcher.main(argv=arguments) == 0
    assert run.call_args.kwargs["env"][launcher.RESEARCH_ENVIRONMENT_VARIABLE] == str(profile)
    profile.write_text(f'workspace="{workbench_service.resource.resource_path / "research"}"\n')
    assert launcher.main(argv=arguments) == 2
    profile.write_text("unknown=true\n")
    assert launcher.main(argv=arguments) == 2


def test_workbench_can_open_with_no_orthofinder_resource(tmp_path, monkeypatch):
    """Standalone imported-protein work remains accessible before a run is selected."""
    monkeypatch.delenv(launcher.RESOURCE_ENVIRONMENT_VARIABLE, raising=False)
    profile = tmp_path / "profile.toml"
    profile.write_text(f'workspace="{tmp_path / "work"}"\n')
    monkeypatch.setenv(launcher.RESEARCH_ENVIRONMENT_VARIABLE, str(profile))
    test = AppTest.from_file(script_path=str(Path(app.__file__)), default_timeout=30).run()
    test.sidebar.radio[0].set_value("Protein research workbench").run()
    assert not test.exception
    assert any(item.value == "Protein research workbench" for item in test.header)


def test_cli_streams_large_generic_fasta_into_private_job(tmp_path, capsys):
    """Ordinary large families can bypass the web limit without asserting a global centroid."""
    source = tmp_path / "large.fasta"
    source.write_text(data=">local\n" + "A" * 20_000_001 + "\n", encoding="utf-8")
    profile = tmp_path / "profile.toml"
    workspace = tmp_path / "work"
    profile.write_text(data=f'workspace="{workspace}"\n', encoding="utf-8")
    assert (
        workbench_main(
            argv=[
                "--config",
                str(profile),
                "import-fasta",
                "--fasta",
                str(source),
                "--source-label",
                "Local export",
            ]
        )
        == 0
    )
    job = capsys.readouterr().out.strip()
    assert get_job(workspace=workspace, job_id=job)["operation"] == "import_fasta"
    snapshot = workspace / "jobs" / job / "input.fasta"
    before = file_sha256(path=snapshot)
    source.write_text(data="changed", encoding="utf-8")
    assert file_sha256(path=snapshot) == before

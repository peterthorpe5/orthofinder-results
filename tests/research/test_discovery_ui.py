"""End-to-end research discovery, mapping and family-control interactions."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from streamlit.testing.v1 import AppTest

from orthofinder_interrogation_app.research import runner, ui
from orthofinder_interrogation_app.research.catalogue import build_cluster_taxonomy
from orthofinder_interrogation_app.research.jobs import job_directory, list_jobs, submit_job
from orthofinder_interrogation_app.research.runner import run_next_job


def test_discovery_ranking_filters_mapping_and_centroid_exports(
    research_config,
    tmp_path,
    monkeypatch,
):
    """Actual DuckDB rankings, explicit mapping and per-cluster exports share the real UI."""
    config = research_config
    index_output = tmp_path / "taxonomy_index"
    index_output.mkdir()
    taxonomy = build_cluster_taxonomy(config=config, output=index_output)
    profile = tmp_path / "profile.toml"
    profile.write_text(
        data=(
            f'workspace="{config.workspace}"\nindex="{config.index}"\n'
            f'sequences="{config.sequences}"\nsequence_map="{config.sequence_map}"\n'
            f'taxonomy="{config.taxonomy}"\ncluster_taxonomy="{taxonomy}"\nthreads=1\n'
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("ORTHOFINDER_RESEARCH_CONFIG", str(profile))
    test = AppTest.from_file(
        script_path=Path(ui.__file__).with_name("standalone.py"), default_timeout=20
    ).run()
    next(item for item in test.button if item.label == "Rank global clusters").click().run()
    assert not test.exception and test.dataframe
    assert any("clusters qualify" in item.value for item in test.caption)
    next(
        item for item in test.text_input if item.label.startswith("Required species/taxon")
    ).set_value("9606")
    next(
        item for item in test.checkbox if item.label == "Require every selected species/taxon"
    ).check()
    next(
        item
        for item in test.checkbox
        if item.label == "Require at least one reviewed species assignment"
    ).check()
    next(item for item in test.button if item.label == "Rank global clusters").click().run()
    assert not test.exception and test.dataframe
    next(
        item for item in test.text_input if item.label.startswith("Required species/taxon")
    ).set_value("bad")
    next(item for item in test.button if item.label == "Rank global clusters").click().run()
    assert test.error and not test.exception
    next(item for item in test.radio if item.label == "Discover clusters").set_value(
        "Sequence identifiers"
    ).run()
    next(item for item in test.text_area if item.label.startswith("Exact sequence IDs")).set_value(
        "P1 missing multi"
    ).run()
    next(
        item for item in test.button if item.label == "Queue exact global cluster mapping"
    ).click().run()
    assert test.success and not test.exception
    run_next_job(config=config)
    test.run()
    assert not test.exception and test.dataframe
    next(
        item for item in test.text_area if item.label == "Global cluster IDs to retrieve in full"
    ).set_value("GLOBAL_A").run()
    assert any("centroid list" in item.label for item in test.get("download_button"))
    next(
        item for item in test.button if item.label == "Queue complete cluster retrieval"
    ).click().run()
    run_next_job(config=config)
    next(
        item for item in test.button if item.label == "Queue reviewed cluster/species index"
    ).click().run()
    assert test.success and len(list_jobs(workspace=config.workspace)) == 3


def test_family_controls_and_completed_subcluster_tables(
    research_config,
    protein_fasta,
    tmp_path,
    monkeypatch,
):
    """Completed subclusters render size plots and drive explicit alignment controls."""
    config = research_config

    def subcluster_fixture(*, config, source, output, parameters):
        """Publish deterministic scientific artifacts for UI contract testing."""
        shutil.copyfile(src=source, dst=output / "sequences.fasta")
        (output / "subclusters.tsv").write_text(
            "subcluster_id\tmember_count\nsub_A\t3\n",
            encoding="utf-8",
        )
        (output / "membership.tsv").write_text(
            "sequence_id\tsubcluster_id\na\tsub_A\nb\tsub_A\nc\tsub_A\n",
            encoding="utf-8",
        )
        return {
            name: output / name for name in ("sequences.fasta", "membership.tsv", "subclusters.tsv")
        }

    monkeypatch.setattr(runner, "subcluster", subcluster_fixture)
    job = submit_job(
        workspace=config.workspace, operation="subcluster", parameters={}, fasta_path=protein_fasta
    )
    run_next_job(config=config, job_id=job)
    profile = tmp_path / "profile.toml"
    profile.write_text(data=f'workspace="{config.workspace}"\n', encoding="utf-8")
    monkeypatch.setenv("ORTHOFINDER_RESEARCH_CONFIG", str(profile))
    test = AppTest.from_file(
        script_path=Path(ui.__file__).with_name("standalone.py"), default_timeout=20
    ).run()
    next(item for item in test.radio if item.label == "Research section").set_value(
        "Family analysis"
    ).run()
    assert not test.exception and test.dataframe
    assert len(test.get("download_button")) >= 2
    next(item for item in test.radio if item.label == "Family analysis").set_value(
        "Alignment and tree"
    ).run()
    next(
        item for item in test.text_input if item.label.startswith("Exact audited subcluster")
    ).set_value("sub_A")
    next(item for item in test.text_input if item.label.startswith("Reference/query")).set_value(
        "a"
    )
    next(
        item
        for item in test.checkbox
        if item.label.startswith("Use an explicit deterministic sample")
    ).check()
    next(item for item in test.number_input if item.label.startswith("Sample size")).set_value(2)
    next(item for item in test.button if item.label == "Queue family analysis").click().run()
    assert not test.exception and test.success
    pending = next(
        row for row in list_jobs(workspace=config.workspace) if row["operation"] == "align"
    )
    directory = job_directory(workspace=config.workspace, job_id=pending["job_id"])
    request = json.loads((directory / "request.json").read_text(encoding="utf-8"))
    assert request["parameters"]["sample_size"] == 2
    next(item for item in test.radio if item.label == "Family analysis").set_value(
        "Representative homology search"
    ).run()
    assert not test.exception
    next(item for item in test.button if item.label == "Queue family analysis").click().run()
    assert test.success

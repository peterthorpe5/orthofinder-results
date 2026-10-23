"""Viewer escaping, source adapters and lazy generic Streamlit integration tests."""

import hashlib
import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from orthofinder_interrogation_app.research.common import ResearchError
from orthofinder_interrogation_app.research.jobs import artifact_path, submit_job
from orthofinder_interrogation_app.research.runner import run_next_job
from orthofinder_interrogation_app.research.ui import (
    completed_manifests,
    read_table,
    split_identifiers,
)
from orthofinder_interrogation_app.research.viewers import (
    auspice_dataset_url,
    json_for_html,
    nightingale_features,
    viewer_asset,
    viewer_document,
)


def test_viewers_escape_script_content_and_ship_real_libraries():
    """Viewer documents use checksummed local bundles and inert escaped input data."""
    encoded = json_for_html(value={"x": "</script><script>alert(1)</script>&"})
    assert "</script>" not in encoded and json.loads(encoded)["x"].startswith("</script>")
    for name in ("alignment", "jbrowse", "nightingale"):
        document = viewer_document(
            viewer=name, payload={"msa": ">a\nACD", "sequence": "ACD", "features": []}
        )
        assert "https://unpkg.com" not in document and 'id="protein-data"' in document
        assert len(document) > 100_000
    assert "ReactMSAView" in viewer_asset(name="jbrowse.js")
    for name in ("../escape.js", "missing.js"):
        with pytest.raises(ResearchError):
            viewer_asset(name=name)
    with pytest.raises(ResearchError):
        viewer_document(viewer="bad", payload={})
    assert auspice_dataset_url(base_url="http://localhost:4000/proxy", job_id="a" * 32).endswith(
        "/" + "a" * 32 + "?m=div"
    )
    for url in ("javascript:alert(1)", "https://user:pass@host", "https://host/?token=secret"):
        with pytest.raises(ResearchError):
            auspice_dataset_url(base_url=url, job_id="a" * 32)


def test_nightingale_coordinates_require_exact_sequence_version():
    """Domains from another isoform cannot be silently overlaid on a model."""
    row = {
        "sequence_id": "a",
        "sequence_sha256": hashlib.sha256(b"ACDE").hexdigest(),
        "signature_accession": "PF1",
        "start": "1",
        "end": "3",
    }
    features = nightingale_features(
        sequence_id="a", sequence="ACDE", annotations=[row, {"sequence_id": "b"}]
    )
    assert features[0]["locations"][0]["fragments"][0] == {"start": 1, "end": 3}
    for changed in (
        {**row, "end": "5"},
        {**row, "start": "bad"},
        {**row, "sequence_sha256": "different"},
    ):
        with pytest.raises(ResearchError):
            nightingale_features(sequence_id="a", sequence="ACDE", annotations=[changed])
    assert (
        nightingale_features(
            sequence_id="a",
            sequence="ACDE",
            annotations=[{"sequence_id": "a", "description": "eggNOG whole-protein evidence"}],
        )
        == []
    )


def test_ui_helpers_and_corrupt_manifest(research_config):
    """UI metadata stays bounded and literal NA identifiers remain identifiers."""
    assert split_identifiers(value="a,b; a\nNA") == ["a", "b", "NA"]
    job = submit_job(
        workspace=research_config.workspace,
        operation="retrieve",
        parameters={"cluster_ids": ["GLOBAL_A"]},
    )
    run_next_job(config=research_config, job_id=job)
    manifests = completed_manifests(config=research_config, operation="retrieve")
    assert (
        job in manifests and completed_manifests(config=research_config, operation="predict") == {}
    )
    path = artifact_path(workspace=research_config.workspace, job_id=job, artifact="clusters.tsv")
    assert len(read_table(path=path, maximum_rows=1)) == 1
    manifest_path = research_config.workspace / "jobs" / job / "manifest.json"
    manifest_path.write_text("bad")
    with pytest.raises(ResearchError, match="manifest"):
        completed_manifests(config=research_config)


def test_generic_app_all_empty_sections_and_real_queue(research_config, tmp_path, monkeypatch):
    """Every generic page renders without E3 data and explicit clicks submit work."""
    config = tmp_path / "profile.toml"
    config.write_text(
        f'workspace="{research_config.workspace}"\nindex="{research_config.index}"\n'
        f'sequences="{research_config.sequences}"\nsequence_map="{research_config.sequence_map}"\n'
    )
    monkeypatch.setenv("ORTHOFINDER_RESEARCH_CONFIG", str(config))
    from orthofinder_interrogation_app.research import ui

    app = AppTest.from_file(
        str(Path(ui.__file__).with_name("standalone.py")), default_timeout=10
    ).run()
    assert not app.exception
    for section in (
        "Family analysis",
        "Alignments and trees",
        "Models",
        "Annotations",
        "Jobs",
        "Discover",
    ):
        app.radio[0].set_value(section).run()
        assert not app.exception
    app.text_area(key="research_retrieve_ids").set_value("GLOBAL_A").run()
    next(
        button for button in app.button if button.label == "Queue complete cluster retrieval"
    ).click().run()
    assert app.success and not app.exception
    app.radio[0].set_value("Jobs").run()
    assert not app.exception and len(app.dataframe) >= 1


def test_smoke_document_generator_is_reproducible(tmp_path):
    """The shipped developer script creates all three offline synthetic documents."""
    import runpy

    script = Path(__file__).parents[2] / "scripts" / "build_viewer_smoke.py"
    main = runpy.run_path(str(script))["main"]
    assert main(argv=["--output", str(tmp_path)]) == 0
    assert sorted(path.name for path in tmp_path.glob("*.html")) == [
        "alignment.html",
        "jbrowse.html",
        "nightingale.html",
    ]


def test_populated_result_pages_render_downloads_and_real_components(
    research_config, model_pdb, tmp_path, monkeypatch
):
    """Completed jobs drive alignment, model, feature and artifact result interfaces."""
    from unittest.mock import Mock

    from orthofinder_interrogation_app.research import analysis, providers, ui

    config = research_config
    profile = tmp_path / "ui.toml"
    profile.write_text(f'workspace="{config.workspace}"\n')
    monkeypatch.setenv("ORTHOFINDER_RESEARCH_CONFIG", str(profile))

    def command(*, arguments, output, config, stdout=None):
        """Return deterministic scientific tool outputs without running a provider."""
        if "--version" in arguments:
            stdout.write_text("version\n")
        elif "--auto" in arguments:
            stdout.write_text(Path(arguments[-1]).read_text())
        else:
            stdout.write_text("(p:0.1,b:0.2,c:0.3);\n")

    monkeypatch.setattr(analysis, "run_command", command)
    alignment = submit_job(
        workspace=config.workspace,
        operation="align",
        parameters={},
        fasta=">p\nAC\n>b\nAD\n>c\nAE\n",
    )
    run_next_job(config=config, job_id=alignment)
    monkeypatch.setattr(providers, "http_request", Mock(return_value=model_pdb.encode()))
    prediction = submit_job(
        workspace=config.workspace,
        operation="predict",
        parameters={"engine": "atlas", "remote_allowed": True},
        fasta=">p\nAC\n",
    )
    run_next_job(config=config, job_id=prediction)

    def remote(*, sequences, output, config, remote_allowed):
        """Return a standard exact-sequence positional InterPro fixture."""
        import hashlib

        path = output / "interpro.raw.tsv"
        md5 = hashlib.md5(b"AC", usedforsecurity=False).hexdigest()
        path.write_text(f"p\t{md5}\t2\tPfam\tPF1\tDomain\t1\t2\t1e-9\tT\tdate\n")
        return path

    monkeypatch.setattr(providers, "remote_interpro", remote)
    annotation = submit_job(
        workspace=config.workspace,
        operation="annotate",
        parameters={"engine": "interpro_remote", "remote_allowed": True},
        fasta=">p\nAC\n",
    )
    run_next_job(config=config, job_id=annotation)
    assert all(
        manifest["operation"] in {"align", "predict", "annotate"}
        for manifest in completed_manifests(config=config).values()
    )
    app = AppTest.from_file(
        str(Path(ui.__file__).with_name("standalone.py")), default_timeout=15
    ).run()
    for section in ("Alignments and trees", "Models", "Annotations", "Jobs", "Family analysis"):
        app.radio[0].set_value(section).run()
        assert not app.exception and not app.error
        if section == "Alignments and trees":
            for view in ("Alignment Viewer 2", "Statistics", "Auspice tree", "JBrowse MSA"):
                app.radio[1].set_value(view).run()
                assert not app.exception and not app.error
        if section == "Models":
            next(
                item for item in app.selectbox if item.label == "Optional annotation overlay"
            ).set_value(annotation).run()
            assert not app.exception and not app.error
        if section == "Jobs":
            next(
                item for item in app.button if item.label == "Prepare verified download"
            ).click().run()
            assert not app.exception and not app.error

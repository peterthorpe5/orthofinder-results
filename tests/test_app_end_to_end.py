"""End-to-end tests for the rendered Streamlit application."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from orthofinder_interrogation_app import app
from orthofinder_interrogation_app.evolutionary_page import (
    ACTIVE_GROUP_STATE,
    COMPARISON_STATE,
)
from orthofinder_interrogation_app.launcher import (
    CACHE_ENVIRONMENT_VARIABLE,
    EXPECTED_TAXA_ENVIRONMENT_VARIABLE,
    FOCUS_ENVIRONMENT_VARIABLE,
    RESOURCE_ENVIRONMENT_VARIABLE,
    TAXONOMY_ENVIRONMENT_VARIABLE,
)
from orthofinder_interrogation_app.report_data import load_visualisation_catalog
from orthofinder_interrogation_app.resource import open_resource


def _application_test(*, application_resource: Path, monkeypatch: pytest.MonkeyPatch) -> AppTest:
    """Run the real Streamlit script against a synthetic DuckDB resource."""

    monkeypatch.setenv(RESOURCE_ENVIRONMENT_VARIABLE, str(application_resource))
    test = AppTest.from_file(str(Path(app.__file__).resolve()), default_timeout=20)
    return test.run()


def test_overview_group_search_and_help_routes(
    application_resource: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Core navigation and DuckDB-backed group inspection render end to end."""

    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    assert not test.exception
    assert test.title[0].value == "OrthoFinder Interrogation"
    assert any(metric.label == "Group records" and metric.value == "4" for metric in test.metric)
    assert any(header.value == "Dataset summary" for header in test.header)
    assert any("Which groups contain my species?" in item.value for item in test.markdown)
    next(button for button in test.button if button.label == "Find E3 focus clusters").click()
    test.run()
    assert not test.exception
    assert any(header.value == "Focus protein clusters" for header in test.header)
    test.sidebar.radio[0].set_value("Find groups")
    test.run()
    assert not test.exception
    assert any(header.value == "Find groups" for header in test.header)
    assert any("4 matching groups" in caption.value for caption in test.caption)
    assert len(test.dataframe) >= 3
    assert any(metric.label == "Stored mean pair distance" for metric in test.metric)
    test.sidebar.radio[0].set_value("Help")
    test.run()
    assert not test.exception
    assert any(header.value == "Using the app" for header in test.header)
    assert any("Exactly the selected species set" in markdown.value for markdown in test.markdown)


def test_methods_and_glossary_routes_are_complete_and_downloadable(
    application_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Methods and glossary pages render against the validated resource end to end."""

    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Methods")
    test.run()
    assert not test.exception
    assert any(header.value == "Methods and provenance" for header in test.header)
    assert any(subheader.value == "Opened resource" for subheader in test.subheader)
    assert any("Completed OrthoFinder analysis" in expander.label for expander in test.expander)
    assert len(test.get("download_button")) >= 4

    test.sidebar.radio[0].set_value("Glossary")
    test.run()
    assert not test.exception
    assert any(header.value == "Glossary" for header in test.header)
    assert any("terms shown" in caption.value for caption in test.caption)
    assert test.dataframe
    assert len(test.get("download_button")) >= 2


def test_offline_report_and_invalid_resource_routes(
    application_resource: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Report download and visible resource errors render without blank pages."""

    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Offline report")
    test.run()
    assert not test.exception
    assert any(header.value == "Download the offline report" for header in test.header)
    assert test.button or test.get("download_button")
    test.sidebar.text_input[0].set_value(str(application_resource / "missing"))
    test.run()
    assert not test.exception
    assert test.error
    assert "could not be opened" in test.error[0].value


def test_all_distance_results_route_exports_complete_selection(
    application_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The new dataset-wide page renders every persisted distance group end to end."""

    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("All distance results")
    test.run()
    assert not test.exception
    assert any(header.value == "All distance results" for header in test.header)
    assert any(
        metric.label == "Clusters in complete export" and metric.value == "2"
        for metric in test.metric
    )
    assert test.dataframe
    assert len(test.get("download_button")) >= 2
    assert any(
        "one preferred stored result per cluster" in markdown.value for markdown in test.markdown
    )


def test_dispersion_benchmark_route_explains_older_resources(
    application_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The new page remains explicit and safe for resources without benchmarks."""

    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Dispersion benchmarks")
    test.run()
    assert not test.exception
    assert any(header.value == "Calibrated dispersion results" for header in test.header)
    assert any("predates matched-background" in item.value for item in test.info)


def test_dispersion_benchmark_route_renders_complete_inference(
    benchmark_application_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A benchmark resource renders profiles, plots, tests and paired downloads."""

    test = _application_test(
        application_resource=benchmark_application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Dispersion benchmarks")
    test.run()
    assert not test.exception
    assert any(header.value == "Calibrated dispersion results" for header in test.header)
    assert any(
        subheader.value == "Genes and proteins defining each biological profile"
        for subheader in test.subheader
    )
    assert any(subheader.value == "Compare biological backgrounds" for subheader in test.subheader)
    assert any(
        subheader.value == "Test one cluster against every background"
        for subheader in test.subheader
    )
    assert len(test.dataframe) >= 4
    assert len(test.get("plotly_chart")) >= 4
    assert len(test.get("download_button")) >= 6
    assert any("predetermined answers" in item.value for item in test.markdown)


def test_focus_protein_clusters_use_custom_authority_and_open_visuals(
    application_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A replacement focus list finds clusters and exposes the complete visual suite."""

    focus = tmp_path / "focus.tsv"
    focus.write_text(
        "protein_identifier\tprotein_name\tcategory\nalpha_1\tAlpha protein\tE3 fixture\n",
        encoding="utf-8",
    )
    monkeypatch.setenv(FOCUS_ENVIRONMENT_VARIABLE, str(focus))
    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Focus protein clusters")
    test.run()
    assert not test.exception
    assert any(header.value == "Focus protein clusters" for header in test.header)
    assert any(
        metric.label == "Matching clusters" and metric.value == "1" for metric in test.metric
    )
    assert test.dataframe
    open_button = next(
        button
        for button in test.button
        if button.label == "Show this focus cluster's complete visual suite"
    )
    open_button.click()
    test.run()
    assert not test.exception
    assert any("Highlighting configured focus protein" in info.value for info in test.info)
    assert len(test.get("plotly_chart")) >= 5


def test_selection_coverage_tree_keeps_controls_tree_and_groups_synchronised(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Headless selector changes update one last-valid tree, summary and group audit."""

    focus = tmp_path / "focus.tsv"
    focus.write_text("protein_identifier\nalpha_1\n", encoding="utf-8")
    monkeypatch.setenv(FOCUS_ENVIRONMENT_VARIABLE, str(focus))
    monkeypatch.setenv(TAXONOMY_ENVIRONMENT_VARIABLE, str(taxonomy_mapping_file))
    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Selection coverage tree")
    test.run()
    assert not test.exception
    assert any(header.value == "Selection coverage tree" for header in test.header)
    include = next(widget for widget in test.multiselect if widget.label == "Include clades")
    exclude = next(widget for widget in test.multiselect if widget.label == "Exclude exact taxa")
    include.set_value(["Target clade | unranked | ID 10"])
    exclude.set_value(["Species gamma | unranked | ID 103"])
    build = next(
        button for button in test.button if button.label == "Build selection coverage tree"
    )
    build.click()
    test.run()
    assert not test.exception
    assert any(metric.label == "Groups passing" and metric.value == "1" for metric in test.metric)
    assert test.get("plotly_chart")
    assert any(
        "Require at least one member from clade Target clade" in item.value
        for item in test.markdown
    )
    assert len(test.dataframe) >= 2
    assert len(test.get("download_button")) >= 8
    required = next(widget for widget in test.multiselect if widget.label == "Required exact taxa")
    required.set_value(["Species alpha | unranked | ID 101"])
    exclude = next(widget for widget in test.multiselect if widget.label == "Exclude exact taxa")
    exclude.set_value(
        [
            "Species gamma | unranked | ID 103",
            "Species alpha | unranked | ID 101",
        ]
    )
    next(
        button for button in test.button if button.label == "Build selection coverage tree"
    ).click()
    test.run()
    assert test.error and "same exact taxon" in test.error[-1].value
    assert any(metric.label == "Groups passing" and metric.value == "1" for metric in test.metric)
    next(
        button for button in test.button if button.label == "Open passing group in cluster explorer"
    ).click()
    test.run()
    assert not test.exception
    assert any(header.value == "Explore one cluster" for header in test.header)


def test_focus_and_coverage_pages_explain_empty_or_invalid_authorities(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Page-level failure states remain visible and distinguish absence from error."""

    monkeypatch.delenv(TAXONOMY_ENVIRONMENT_VARIABLE, raising=False)
    monkeypatch.delenv(EXPECTED_TAXA_ENVIRONMENT_VARIABLE, raising=False)
    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Selection coverage tree")
    test.run()
    assert not test.exception
    assert any("species-to-taxonomy review template" in item.value for item in test.warning)
    assert any("not an E3 protein list" in item.value for item in test.markdown)
    assert any(
        button.label == "Download this dataset's taxonomy review template"
        for button in test.get("download_button")
    )

    monkeypatch.setenv(TAXONOMY_ENVIRONMENT_VARIABLE, str(taxonomy_mapping_file))
    monkeypatch.setenv(
        EXPECTED_TAXA_ENVIRONMENT_VARIABLE,
        str(tmp_path / "missing_expected_taxa.tsv"),
    )
    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Selection coverage tree")
    test.run()
    assert not test.exception
    assert any("Expected-taxon universe could not be used" in item.value for item in test.error)

    monkeypatch.delenv(EXPECTED_TAXA_ENVIRONMENT_VARIABLE, raising=False)
    absent_focus = tmp_path / "absent_focus.tsv"
    absent_focus.write_text(
        "protein_identifier\nnot_present_in_resource\n",
        encoding="utf-8",
    )
    monkeypatch.setenv(FOCUS_ENVIRONMENT_VARIABLE, str(absent_focus))
    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Focus protein clusters")
    test.run()
    assert not test.exception
    assert any("No configured focus protein matched" in item.value for item in test.warning)


def test_protein_search_opens_focused_cluster_and_distance_table(
    application_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One canonical protein leads directly to highlighted visuals and distances."""

    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Find a protein")
    test.run()
    assert not test.exception
    assert any(header.value == "Find a gene or protein" for header in test.header)
    query_input = next(
        widget for widget in test.text_input if widget.label == "Protein or OrthoFinder internal ID"
    )
    query_input.set_value("alpha_1")
    next(button for button in test.button if button.label == "Find this protein").click()
    test.run()
    assert not test.exception
    assert any(
        metric.label == "Matching cluster records" and metric.value == "1" for metric in test.metric
    )
    assert any(
        subheader.value == "Distances from the requested protein" for subheader in test.subheader
    )
    assert any(
        metric.label == "Proteins compared" and metric.value == "2" for metric in test.metric
    )
    assert len(test.get("plotly_chart")) >= 9
    assert len(test.get("download_button")) >= 4


def test_evolutionary_and_reviewed_taxonomy_routes(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Visual and descendant-aware pages execute end to end against real fixtures."""

    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Cluster explorer")
    test.run()
    assert not test.exception
    assert any(header.value == "Explore one cluster" for header in test.header)
    assert len(test.get("plotly_chart")) >= 5

    test.sidebar.radio[0].set_value("Taxonomic search")
    test.run()
    assert not test.exception
    assert any(header.value == "Taxonomic search" for header in test.header)
    assert any("Load a reviewed mapping" in info.value for info in test.info)
    test.sidebar.text_input[2].set_value(str(taxonomy_mapping_file))
    test.run()
    assert not test.exception
    assert any(caption.value.startswith("Showing 1") for caption in test.caption)
    assert len(test.dataframe) >= 2


def test_empty_resource_path_prompts_for_input(monkeypatch: pytest.MonkeyPatch) -> None:
    """An unset launcher path produces a useful start state."""

    monkeypatch.delenv(RESOURCE_ENVIRONMENT_VARIABLE, raising=False)
    test = AppTest.from_file(str(Path(app.__file__).resolve()), default_timeout=20).run()
    assert not test.exception
    assert test.info[0].value.startswith("Provide a completed resource path")


def test_empty_and_populated_cluster_comparison_routes(
    application_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The comparison workspace gives useful empty and two-group rendered states."""

    test = _application_test(
        application_resource=application_resource,
        monkeypatch=monkeypatch,
    )
    test.sidebar.radio[0].set_value("Compare clusters")
    test.run()
    assert not test.exception
    assert any(header.value == "Compare clusters" for header in test.header)
    assert any("workspace is empty" in info.value for info in test.info)

    resource = open_resource(path=application_resource)
    assert resource.report_path is not None
    catalog = load_visualisation_catalog(
        report_path=resource.report_path,
        expected_run_id="test_run",
    )
    first = copy.deepcopy(catalog.get(label="HOG|N0|N0.HOG1"))
    second = copy.deepcopy(first)
    second["label"] = "HOG | N0 | N0.HOG3"
    payload = {
        "run": {"run_id": "test_run"},
        "networks": {
            "HOG|N0|N0.HOG1": first,
            "HOG|N0|N0.HOG3": second,
        },
        "limits": {"nearestNeighbours": 3},
    }
    resource.report_path.write_text(
        "<!doctype html><title>Comparison test</title>"
        '<script id="orthofinder-results-data" type="application/json">'
        + json.dumps(payload, separators=(",", ":"))
        + "</script>",
        encoding="utf-8",
    )
    test.session_state[COMPARISON_STATE] = [
        {
            "run_id": "test_run",
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": "N0.HOG1",
        },
        {
            "run_id": "test_run",
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": "N0.HOG3",
        },
    ]
    test.run()
    assert not test.exception
    assert len(test.get("plotly_chart")) >= 4
    assert any(subheader.value == "Comparison data and provenance" for subheader in test.subheader)


def test_schema3_pipeline_to_lazy_cluster_explorer_end_to_end(
    schema3_resource: Path,
    persistent_test_root: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A rebuilt resource drives a non-pilot HOG from portable tree through the UI."""

    cache = persistent_test_root / "e2e_analysis_cache"
    monkeypatch.setenv(CACHE_ENVIRONMENT_VARIABLE, str(cache))
    test = _application_test(
        application_resource=schema3_resource,
        monkeypatch=monkeypatch,
    )
    test.session_state[ACTIVE_GROUP_STATE] = {
        "run_id": "schema3-run",
        "group_type": "HOG",
        "hierarchy_node": "N1",
        "group_id": "N1.HOG0000001",
    }
    test.sidebar.radio[0].set_value("Cluster explorer")
    test.run()
    assert not test.exception
    assert len(test.get("plotly_chart")) >= 10
    assert any("PORTABLE_TREE_LAZY" in caption.value for caption in test.caption)
    assert tuple(cache.rglob("group_analysis_*.json.gz"))
    test.run()
    assert not test.exception
    assert any("CACHE_HIT" in caption.value for caption in test.caption)

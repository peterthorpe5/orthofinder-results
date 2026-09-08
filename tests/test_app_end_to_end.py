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
    RESOURCE_ENVIRONMENT_VARIABLE,
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
    assert any(
        metric.label == "Group records" and metric.value == "4" for metric in test.metric
    )
    assert any(header.value == "Dataset summary" for header in test.header)
    assert any("Which groups contain my species?" in item.value for item in test.markdown)
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
    assert any(header.value == "Help & glossary" for header in test.header)
    assert any("Exactly the selected species set" in markdown.value for markdown in test.markdown)


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
        "one preferred stored result per cluster" in markdown.value
        for markdown in test.markdown
    )


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
    assert any(
        subheader.value == "Comparison data and provenance"
        for subheader in test.subheader
    )


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

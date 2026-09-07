"""End-to-end tests for the rendered Streamlit application."""

from __future__ import annotations

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from orthofinder_interrogation_app import app
from orthofinder_interrogation_app.launcher import RESOURCE_ENVIRONMENT_VARIABLE


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
    assert any(metric.label == "Groups" and metric.value == "4" for metric in test.metric)
    test.sidebar.radio[0].set_value("Group search")
    test.run()
    assert not test.exception
    assert any(header.value == "Group search" for header in test.header)
    assert any("4 matching groups" in caption.value for caption in test.caption)
    assert len(test.dataframe) >= 3
    assert any(metric.label == "Mean distance" for metric in test.metric)
    test.sidebar.radio[0].set_value("Help")
    test.run()
    assert not test.exception
    assert any(header.value == "Help and interpretation" for header in test.header)
    assert any("Exactly this species set" in markdown.value for markdown in test.markdown)


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
    assert any(header.value == "Offline report" for header in test.header)
    assert test.button or test.get("download_button")
    test.sidebar.text_input[0].set_value(str(application_resource / "missing"))
    test.run()
    assert not test.exception
    assert test.error
    assert "could not be opened" in test.error[0].value


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
    test.sidebar.radio[0].set_value("Evolutionary views")
    test.run()
    assert not test.exception
    assert any(header.value == "Evolutionary views" for header in test.header)
    assert len(test.get("plotly_chart")) >= 5

    test.sidebar.radio[0].set_value("Taxonomic search")
    test.run()
    assert not test.exception
    assert any(header.value == "Taxonomic search" for header in test.header)
    assert any("Load a reviewed mapping" in info.value for info in test.info)
    test.sidebar.text_input[1].set_value(str(taxonomy_mapping_file))
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

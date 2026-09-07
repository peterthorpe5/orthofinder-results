"""Tests for visible evolutionary and taxonomy page fallback states."""

from __future__ import annotations

import copy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from orthofinder_interrogation_app import evolutionary_page, taxonomy_page
from orthofinder_interrogation_app.report_data import (
    VisualisationCatalog,
    load_visualisation_catalog,
)
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_interrogation_app.taxonomy import taxonomy_template
from orthofinder_results.errors import InputValidationError


def _entry(*, application_resource: Path) -> dict:
    """Return the synthetic report's validated visual entry."""

    resource = open_resource(path=application_resource)
    assert resource.report_path is not None
    catalog = load_visualisation_catalog(
        report_path=resource.report_path, expected_run_id=resource.run_id
    )
    return catalog.get(label=catalog.labels()[0])


def test_evolutionary_page_early_fallbacks(
    application_resource: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Absent, invalid and empty visual catalogs produce visible messages."""

    streamlit = MagicMock()
    monkeypatch.setattr(evolutionary_page, "st", streamlit)
    evolutionary_page.render_evolutionary_views(
        resource=SimpleNamespace(report_path=None, run_id="test_run")
    )
    assert streamlit.warning.called

    streamlit.reset_mock()
    evolutionary_page.render_evolutionary_views(
        resource=SimpleNamespace(
            report_path=application_resource / "missing.html", run_id="test_run"
        )
    )
    assert streamlit.error.called

    streamlit.reset_mock()
    monkeypatch.setattr(
        evolutionary_page,
        "_cached_catalog",
        lambda **kwargs: VisualisationCatalog(run_id="test_run", networks={}, nearest_neighbours=0),
    )
    report = application_resource / "report" / "orthofinder_results_summary.html"
    evolutionary_page.render_evolutionary_views(
        resource=SimpleNamespace(report_path=report, run_id="test_run")
    )
    assert streamlit.info.called


def test_evolutionary_panel_warnings_and_quality_categories(
    application_resource: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Unavailable panels and conservative POOR/MODERATE labels stay visible."""

    streamlit = MagicMock()
    streamlit.checkbox.return_value = False
    monkeypatch.setattr(evolutionary_page, "st", streamlit)
    entry = _entry(application_resource=application_resource)

    unavailable = copy.deepcopy(entry)
    unavailable["distanceProjection"] = {"status": "UNAVAILABLE", "reason": "no matrix"}
    evolutionary_page._render_pcoa(entry=unavailable, linked=frozenset())
    evolutionary_page._render_shepard(entry=unavailable)
    assert streamlit.warning.call_count >= 2

    for category in ("POOR", "MODERATE"):
        quality = copy.deepcopy(entry)
        quality["distanceProjection"]["quality_category"] = category
        quality["distanceProjection"]["quality_explanation"] = f"{category} explanation"
        evolutionary_page._render_pcoa(entry=quality, linked=frozenset())
    assert streamlit.info.called

    bad_tree = copy.deepcopy(entry)
    bad_tree["phylogram"] = {"status": "UNAVAILABLE", "reason": "no tree"}
    evolutionary_page._render_phylogram(entry=bad_tree, linked=frozenset())
    bad_matrix = copy.deepcopy(entry)
    bad_matrix["distanceMatrix"] = {"status": "UNAVAILABLE", "reason": "no distances"}
    evolutionary_page._render_matrix(entry=bad_matrix, linked=frozenset())
    bad_topology = copy.deepcopy(entry)
    bad_topology["edges"] = "invalid"
    evolutionary_page._render_topology(entry=bad_topology, linked=frozenset(), nearest_neighbours=3)
    assert streamlit.warning.call_count >= 5

    unresolved = copy.deepcopy(entry)
    unresolved["phylogram"]["unresolvedMembers"] = ["missing"]
    evolutionary_page._render_phylogram(entry=unresolved, linked=frozenset())
    assert any("Unresolved" in str(call) for call in streamlit.warning.call_args_list)


def test_evolutionary_value_and_member_validation_helpers() -> None:
    """Formatting and malformed member records retain explicit unavailable states."""

    assert evolutionary_page._format_optional(None) == "Unavailable"
    assert evolutionary_page._format_optional(0.125) == "0.125"
    assert evolutionary_page._format_percent(0.125) == "12.5%"
    assert evolutionary_page._format_percent("") == "Unavailable"
    with pytest.raises(InputValidationError, match="displayed members"):
        evolutionary_page._members(entry={})
    with pytest.raises(InputValidationError, match="malformed"):
        evolutionary_page._members(entry={"members": ["bad"]})
    with pytest.raises(InputValidationError, match="unlabelled"):
        evolutionary_page._members(entry={"members": [{"member_id": "", "species_label": "x"}]})
    with pytest.raises(InputValidationError, match="valid required"):
        evolutionary_page._required_mapping(entry={}, key="required")


def test_taxonomy_page_mapping_precedence_and_scope_messages(
    taxonomy_mapping_file: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Uploads override paths and every claim-limiting search message is exercised."""

    expected = ("Species_A", "Species_B", "Species_C", "Species_D")
    uploaded = taxonomy_page._mapping_authority(
        expected_species=expected,
        taxonomy_path_text="/definitely/missing.tsv",
        uploaded_data=taxonomy_mapping_file.read_bytes(),
    )
    assert uploaded is not None and len(uploaded.reviewed_species) == 3
    assert (
        taxonomy_page._mapping_authority(
            expected_species=expected, taxonomy_path_text="", uploaded_data=None
        )
        is None
    )

    streamlit = MagicMock()
    monkeypatch.setattr(taxonomy_page, "st", streamlit)
    option = SimpleNamespace(name="Target clade", taxon_id=10)
    for mode in ("ENRICHED", "SAMPLED_EXCLUSIVE", "NEAR_EXCLUSIVE", "CONTAINS"):
        taxonomy_page._render_search_scope(
            mode=mode,
            option=option,
            result=SimpleNamespace(
                rows=(),
                page_number=1,
                page_size=100,
                total_rows=0,
                target_species_count=2,
                outside_species_count=1,
                unresolved_species_count=1,
                tested_group_count=3,
            ),
        )
    assert streamlit.info.call_count == 3
    row = {
        "group_id": "g",
        "member_count": 2,
        "species_count": 2,
        "target_species_count": 1,
        "target_coverage": 0.5,
        "mapped_target_fraction": 0.5,
        "outside_species_count": 1,
        "outsider_species_labels": "Species_C",
        "unresolved_species_labels": "",
        "enrichment_odds_ratio": 2.0,
        "enrichment_p_value": 0.02,
        "enrichment_q_value": 0.04,
    }
    assert taxonomy_page._display_taxonomy_row(row=row)["BH q-value"] == 0.04


def test_taxonomy_page_invalid_and_unreviewed_mapping_states(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Mapping validation errors and zero reviewed records stop before searching."""

    streamlit = MagicMock()
    streamlit.file_uploader.return_value = None
    monkeypatch.setattr(taxonomy_page, "st", streamlit)
    service = MagicMock()
    service.list_species.return_value = ("Species_A",)
    service.resource.run_id = "test_run"
    taxonomy_page.render_taxonomy_search(
        service=service, taxonomy_path_text=str(tmp_path / "missing.tsv")
    )
    assert streamlit.error.called

    streamlit.reset_mock()
    streamlit.file_uploader.return_value = None
    streamlit.columns.return_value = [MagicMock() for _ in range(4)]
    mapping = tmp_path / "unreviewed.tsv"
    mapping.write_bytes(taxonomy_template(species=("Species_A",)))
    taxonomy_page.render_taxonomy_search(service=service, taxonomy_path_text=str(mapping))
    assert streamlit.warning.called
    service.search_taxonomy_groups.assert_not_called()

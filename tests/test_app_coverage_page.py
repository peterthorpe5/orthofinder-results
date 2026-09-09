"""Unit tests for selection-coverage page adapters and callbacks."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from orthofinder_interrogation_app import coverage_page
from orthofinder_interrogation_app.coverage_exports import build_selection_coverage_run
from orthofinder_interrogation_app.taxonomy_selection import (
    TaxonomyGraph,
    build_taxonomy_graph,
    default_expected_taxa,
    expected_taxa_to_tsv,
)
from orthofinder_results.errors import InputValidationError

_EXPECTED_SPECIES = ("Species_A", "Species_B", "Species_C", "Species_D")


def _graph(*, taxonomy_mapping_file: Path) -> TaxonomyGraph:
    """Return the reviewed application-fixture taxonomy graph."""

    authority = coverage_page._load_taxonomy_authority(
        expected_species=_EXPECTED_SPECIES,
        taxonomy_path_text=str(taxonomy_mapping_file),
        uploaded_data=None,
    )
    assert authority is not None
    return build_taxonomy_graph(authority=authority)


def test_coverage_authority_precedence_and_expected_universe_sources(
    taxonomy_mapping_file: Path,
    tmp_path: Path,
) -> None:
    """Uploads, configured paths and defaults retain explicit precedence."""

    uploaded = coverage_page._load_taxonomy_authority(
        expected_species=_EXPECTED_SPECIES,
        taxonomy_path_text="/missing/ignored.tsv",
        uploaded_data=taxonomy_mapping_file.read_bytes(),
    )
    assert uploaded is not None and len(uploaded.reviewed_species) == 3
    assert (
        coverage_page._load_taxonomy_authority(
            expected_species=_EXPECTED_SPECIES,
            taxonomy_path_text="",
            uploaded_data=None,
        )
        is None
    )
    graph = build_taxonomy_graph(authority=uploaded)
    default = coverage_page._load_expected_authority(
        graph=graph,
        expected_path_text="",
        uploaded_data=None,
    )
    payload = expected_taxa_to_tsv(authority=default)
    uploaded_expected = coverage_page._load_expected_authority(
        graph=graph,
        expected_path_text="/missing/ignored.tsv",
        uploaded_data=payload,
    )
    assert uploaded_expected.source_name == "uploaded_expected_taxa.tsv"
    expected_path = tmp_path / "expected.tsv"
    expected_path.write_bytes(payload)
    from_path = coverage_page._load_expected_authority(
        graph=graph,
        expected_path_text=str(expected_path),
        uploaded_data=None,
    )
    assert from_path.included_taxon_ids == default.included_taxon_ids


def test_coverage_widget_translation_addition_and_clear_callbacks(
    taxonomy_mapping_file: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keyboard actions and visible multiselect labels share one validated state."""

    graph = _graph(taxonomy_mapping_file=taxonomy_mapping_file)
    run = build_selection_coverage_run(
        graph=graph,
        expected=default_expected_taxa(graph=graph),
        selection=coverage_page.make_selection(graph=graph),
        group_species_rows=(),
        run_id="test_run",
        resource_identity="fixture",
        package_version="0.7.0",
        created_at_utc="2026-09-09T12:00:00Z",
    )
    options = coverage_page._node_options(graph=graph)
    terminal_label = next(label for label, value in options.items() if value == "101")
    clade_label = next(label for label, value in options.items() if value == "10")
    streamlit = SimpleNamespace(
        session_state={
            "coverage_required_exact": [terminal_label],
            "coverage_include_clade": [clade_label],
        }
    )
    monkeypatch.setattr(coverage_page, "st", streamlit)
    selection = coverage_page._selection_from_widgets(
        graph=graph,
        node_options=options,
    )
    assert selection.identifiers(selector_type="REQUIRED_EXACT") == ("101",)
    assert selection.identifiers(selector_type="INCLUDE_CLADE") == ("10",)

    coverage_page._add_node_selector(
        selector_type="EXCLUDE_EXACT",
        taxon_id="101",
        run=run,
    )
    coverage_page._add_node_selector(
        selector_type="EXCLUDE_EXACT",
        taxon_id="101",
        run=run,
    )
    assert streamlit.session_state["coverage_exclude_exact"] == [terminal_label]
    coverage_page._clear_node_selectors(taxon_id="101", run=run)
    assert streamlit.session_state["coverage_required_exact"] == []
    assert streamlit.session_state["coverage_exclude_exact"] == []
    assert streamlit.session_state["orthofinder_coverage_rebuild"] is True

    streamlit.session_state["coverage_include_clade"] = ["stale label"]
    with pytest.raises(InputValidationError, match="stale taxonomy label"):
        coverage_page._selection_from_widgets(graph=graph, node_options=options)


def test_coverage_resource_identity_uses_manifest_then_database(
    tmp_path: Path,
) -> None:
    """Resource identity falls back to bounded DuckDB metadata without a manifest."""

    root = tmp_path / "resource"
    database = root / "duckdb" / "orthofinder_results.duckdb"
    database.parent.mkdir(parents=True)
    database.write_bytes(b"fixture-database")
    manifest = root / "run_manifest.json"
    manifest.write_bytes(b'{"run_id":"test_run"}\n')
    resource = SimpleNamespace(
        resource_path=root,
        database_path=database,
        run_id="test_run",
        schema_version=2,
    )
    service = SimpleNamespace(resource=resource)
    assert coverage_page._resource_identity(service=service).startswith("run_manifest_sha256:")
    manifest.unlink()
    assert coverage_page._resource_identity(service=service) == ("duckdb:test_run:size=16:schema=2")


def test_coverage_display_adapters_use_plain_language_headings() -> None:
    """Coverage and predicate records expose stable, self-explanatory headings."""

    coverage = coverage_page._coverage_display(
        row={
            "taxon_id": "101",
            "taxon_name": "Species alpha",
            "taxon_rank": "species",
            "selection_state": "REQUIRED_EXACT",
            "coverage_state": "REPRESENTED_IN_INPUT",
            "represented_input_labels": "Species_A",
            "represented_terminal_count": 1,
            "expected_no_data_terminal_count": 0,
            "excluded_terminal_count": 0,
        }
    )
    assert coverage["Taxon"] == "Species alpha"
    evaluation = coverage_page._evaluation_display(
        row={
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": "N0.HOG1",
            "predicate_pass": True,
            "mapped_member_count": 2,
            "unmapped_member_count": 0,
            "represented_taxon_ids": "101;102",
            "outside_selected_scope_taxon_ids": "",
            "unmapped_species_labels": "",
            "predicate_failure_reasons": "",
            "predicate_audit": "REQUIRED_EXACT:101=PASS",
            "selection_manifest_id": "manifest",
        }
    )
    assert evaluation["Passes all predicates"] is True

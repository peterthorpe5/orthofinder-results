"""Tests for packaged RNA-seq evidence queries and visual summaries."""

from __future__ import annotations

from pathlib import Path

import pytest

from orthofinder_interrogation_app.expression import (
    expression_available,
    expression_dimensions,
    expression_group_species,
    expression_group_summaries,
    expression_heatmap_cells,
    expression_member_evidence,
)
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_interrogation_app.terminal_motif_page import (
    _expression_heatmap_figure,
    _expression_intersections,
    _expression_upset_figure,
    _rank_expression_species,
)
from orthofinder_results.errors import InputValidationError


def test_expression_queries_keep_units_contexts_and_missingness_explicit(
    expression_application_resource: Path,
) -> None:
    """Read-only expression queries return bounded, unit-safe evidence."""

    resource = open_resource(path=expression_application_resource)
    assert resource.schema_version == 5
    assert expression_available(resource=resource)
    assert expression_dimensions(resource=resource) == {
        "species": ("Species_A",),
        "units": ("TPM",),
    }
    summaries = expression_group_summaries(
        resource=resource,
        group_type="HOG",
        hierarchy_node="N0",
        group_ids=("N0.HOG1",),
    )
    assert summaries[0]["unique_mapped_member_count"] == 2
    species = expression_group_species(
        resource=resource,
        group_type="HOG",
        hierarchy_node="N0",
        group_ids=("N0.HOG1",),
    )
    assert {row["species_label"] for row in species} == {"Species_A", "Species_B"}
    beta = next(row for row in species if row["species_label"] == "Species_B")
    assert beta["expression_observed_member_count"] == 0
    cells = expression_heatmap_cells(
        resource=resource,
        group_type="HOG",
        hierarchy_node="N0",
        group_ids=("N0.HOG1",),
        context_column="organism_part",
        expression_unit="TPM",
    )
    assert {(row["context_label"], row["median_expression"]) for row in cells} == {
        ("leaf", 4.5),
        ("root", 2.0),
    }
    members = expression_member_evidence(
        resource=resource,
        group_type="HOG",
        hierarchy_node="N0",
        group_id="N0.HOG1",
    )
    not_mapped = next(row for row in members if row["member_id"] == "beta_1")
    assert not_mapped["mapping_status"] == "NOT_MAPPED"
    assert not_mapped["evidence_status"] == "NOT_MAPPED"


def test_expression_queries_reject_unbounded_or_unsupported_requests(
    expression_application_resource: Path,
) -> None:
    """Dynamic SQL identifiers and result cardinality remain allow-listed."""

    resource = open_resource(path=expression_application_resource)
    with pytest.raises(InputValidationError, match="Unsupported expression group type"):
        expression_group_summaries(
            resource=resource,
            group_type="UNSAFE",
            hierarchy_node="N0",
            group_ids=("N0.HOG1",),
        )
    with pytest.raises(InputValidationError, match="Unsupported expression context"):
        expression_heatmap_cells(
            resource=resource,
            group_type="HOG",
            hierarchy_node="N0",
            group_ids=("N0.HOG1",),
            context_column="unsafe_column",
            expression_unit="TPM",
        )
    with pytest.raises(InputValidationError, match="between 1 and 25"):
        expression_heatmap_cells(
            resource=resource,
            group_type="HOG",
            hierarchy_node="N0",
            group_ids=(),
            context_column="organism_part",
            expression_unit="TPM",
        )
    with pytest.raises(InputValidationError, match="never combined"):
        expression_heatmap_cells(
            resource=resource,
            group_type="HOG",
            hierarchy_node="N0",
            group_ids=("N0.HOG1",),
            context_column="organism_part",
            expression_unit="",
        )
    with pytest.raises(InputValidationError, match="between 1 and 50,000"):
        expression_heatmap_cells(
            resource=resource,
            group_type="HOG",
            hierarchy_node="N0",
            group_ids=("N0.HOG1",),
            context_column="organism_part",
            expression_unit="TPM",
            maximum_cells=0,
        )
    with pytest.raises(InputValidationError, match="NUL"):
        expression_group_summaries(
            resource=resource,
            group_type="HOG",
            hierarchy_node="N0",
            group_ids=("N0.HOG1\x00",),
        )


def test_expression_queries_support_legacy_groups_and_reject_old_resources(
    expression_application_resource: Path,
) -> None:
    """Legacy orthogroup queries omit hierarchy predicates safely."""

    resource = open_resource(path=expression_application_resource)
    assert expression_group_summaries(
        resource=resource,
        group_type="LEGACY_ORTHOGROUP",
        hierarchy_node="",
        group_ids=("OG4",),
    ) == ()
    legacy_species = expression_group_species(
        resource=resource,
        group_type="LEGACY_ORTHOGROUP",
        hierarchy_node="",
        group_ids=("OG4",),
    )
    assert legacy_species[0]["species_label"] == "Species_D"
    assert expression_heatmap_cells(
        resource=resource,
        group_type="LEGACY_ORTHOGROUP",
        hierarchy_node="",
        group_ids=("OG4",),
        context_column="condition",
        expression_unit="TPM",
        species=("Species_D", "", "Species_D"),
    ) == ()
    legacy_members = expression_member_evidence(
        resource=resource,
        group_type="LEGACY_ORTHOGROUP",
        hierarchy_node="",
        group_id="OG4",
    )
    assert legacy_members[0]["member_id"] == "delta_1"


def test_expression_queries_reject_old_resources(application_resource: Path) -> None:
    """Resources predating schema 5 fail with an explicit capability message."""

    old_resource = open_resource(path=application_resource)
    assert not expression_available(resource=old_resource)
    with pytest.raises(InputValidationError, match="no complete RNA-seq"):
        expression_dimensions(resource=old_resource)


def test_expression_heatmap_and_upset_use_observed_evidence_only() -> None:
    """Blank heatmap cells and exact species intersections remain distinguishable."""

    cells = (
        {
            "group_id": "HOG1",
            "species_label": "Species_A",
            "context_label": "leaf",
            "expression_unit": "TPM",
            "median_expression": 3.0,
            "mapped_member_count": 1,
            "experiment_count": 1,
            "positive_context_fraction": 1.0,
        },
    )
    heatmap = _expression_heatmap_figure(
        cells=cells,
        selected_groups=("HOG1", "HOG2"),
        log_transform=True,
    )
    assert heatmap.data[0].z[0][0] == 2.0
    assert heatmap.data[0].z[1][0] is None
    species_rows = (
        {
            "group_id": "HOG1",
            "species_label": "Species_A",
            "expression_observed_member_count": 1,
        },
        {
            "group_id": "HOG1",
            "species_label": "Species_B",
            "expression_observed_member_count": 1,
        },
        {
            "group_id": "HOG2",
            "species_label": "Species_A",
            "expression_observed_member_count": 1,
        },
    )
    assert _rank_expression_species(rows=species_rows) == ("Species_A", "Species_B")
    intersections = _expression_intersections(
        rows=species_rows,
        group_ids=("HOG1", "HOG2", "HOG3"),
        species=("Species_A", "Species_B"),
    )
    assert sum(int(row["Group count"]) for row in intersections) == 3
    assert any(row["Species with expression evidence"] == "None selected" for row in intersections)
    upset = _expression_upset_figure(
        intersections=intersections,
        species=("Species_A", "Species_B"),
    )
    assert len(upset.data) >= 4

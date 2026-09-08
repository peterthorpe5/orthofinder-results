"""Tests for protein-centred application display helpers."""

from __future__ import annotations

import pytest

from orthofinder_interrogation_app import protein_page
from orthofinder_results.errors import InputValidationError


def _raw_result() -> dict[str, object]:
    """Return one complete raw protein-to-cluster result."""

    return {
        "run_id": "run",
        "member_id": "protA",
        "internal_id": "0_1",
        "species_label": "Species_A",
        "match_source": "ORTHOFINDER_INTERNAL_ID",
        "matched_identifier": "0_1",
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "group_id": "N0.HOG1",
        "legacy_orthogroup_id": "OG1",
        "gene_tree_parent_clade": "N0",
        "member_count": 3,
        "species_count": 2,
        "max_copies_per_species": 2,
        "mean_copies_per_species": 1.5,
        "computation_status": "EXACT",
        "sampled_member_count": 3,
        "distance_pair_count": 3,
        "mean_distance": 0.2,
        "population_stddev_distance": 0.1,
    }


def test_protein_result_display_and_group_identity_are_readable() -> None:
    """Raw search fields retain identity while receiving explanatory headings."""

    row = _raw_result()
    displayed = protein_page._display_protein_result(row=row)
    assert displayed["Protein ID"] == "protA"
    assert displayed["Matched through"] == "OrthoFinder internal ID"
    assert displayed["Stored distance coverage"] == "Complete distances stored"
    key = protein_page._result_group_key(row=row)
    assert key.display_label() == "HOG | N0 | N0.HOG1"
    options = protein_page._protein_result_options(rows=(row, dict(row)))
    assert len(options) == 2
    assert tuple(options.values()) == (0, 1)


def test_protein_result_helpers_preserve_unknown_and_missing_states() -> None:
    """Unfamiliar states remain visible and empty result selection is rejected."""

    assert protein_page._distance_status_label(value=None) == "Not calculated"
    assert protein_page._distance_status_label(value="NEW_STATE") == "New State"
    with pytest.raises(InputValidationError, match="at least one"):
        protein_page._protein_result_options(rows=())

"""Unit tests for application models and TSV downloads."""

from __future__ import annotations

import pytest

from orthofinder_interrogation_app.models import (
    GroupKey,
    GroupSearchFilters,
    TaxonomySearchFilters,
)
from orthofinder_interrogation_app.tsv import records_to_tsv
from orthofinder_results.errors import InputValidationError


def test_group_key_and_search_page_controls_are_deterministic() -> None:
    """Composite labels, species selections and pagination are stable."""

    key = GroupKey(run_id="run", group_type="HOG", hierarchy_node="", group_id="HOG1")
    assert key.display_label() == "HOG | ROOT | HOG1"
    filters = GroupSearchFilters(
        included_species=(" Species_B ", "Species_A", "Species_A"),
        excluded_species=("Species_D",),
        page_size=25,
        page_number=3,
    )
    assert filters.included_species == ("Species_A", "Species_B")
    assert filters.offset == 50


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ({"include_mode": "BAD"}, "include mode"),
        ({"distance_availability": "BAD"}, "distance availability"),
        ({"sort_mode": "BAD"}, "sort mode"),
        ({"page_size": 0}, "page_size"),
        ({"page_size": 501}, "page_size"),
        ({"page_number": 0}, "page_number"),
        ({"minimum_member_count": -1}, "Minimum member count"),
        ({"maximum_member_count": -1}, "Maximum member count"),
        (
            {"minimum_member_count": 2, "maximum_member_count": 1},
            "Minimum member count",
        ),
        ({"minimum_species_count": -1}, "Minimum species count"),
        ({"maximum_species_count": -1}, "Maximum species count"),
        (
            {"minimum_species_count": 2, "maximum_species_count": 1},
            "Minimum species count",
        ),
        ({"maximum_mean_distance": -0.1}, "maximum_mean_distance"),
        ({"maximum_distance_sd": -0.1}, "maximum_distance_sd"),
        ({"included_species": ("",)}, "empty label"),
        (
            {
                "included_species": ("Species_A",),
                "excluded_species": ("Species_A",),
            },
            "both required and rejected",
        ),
    ],
)
def test_invalid_search_controls_fail_explicitly(values: dict[str, object], message: str) -> None:
    """Contradictory or unsafe controls never reach DuckDB."""

    with pytest.raises(InputValidationError, match=message):
        GroupSearchFilters(**values)


def test_records_to_tsv_preserves_types_and_missing_values() -> None:
    """TSV output retains exact fields without comma-separated output."""

    content = records_to_tsv(
        records=(
            {"identifier": "A,1", "active": True, "value": None},
            {"identifier": "B", "active": False, "value": 2.5},
        )
    ).decode("utf-8")
    assert content == ("identifier\tactive\tvalue\nA,1\ttrue\t\nB\tfalse\t2.5\n")
    assert records_to_tsv(records=(), fieldnames=("identifier",)) == b"identifier\n"


def test_taxonomy_search_controls_and_offset_are_deterministic() -> None:
    """Taxonomic searches require an exact authority and bounded page."""

    filters = TaxonomySearchFilters(
        group_type="HOG",
        hierarchy_node="N0",
        target_taxon_id=10,
        mode="NEAR_EXCLUSIVE",
        page_size=25,
        page_number=3,
    )
    assert filters.offset == 50


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ({"group_type": ""}, "exact group type"),
        ({"mode": "BAD"}, "Unsupported taxonomy"),
        ({"target_taxon_id": 0}, "positive"),
        ({"minimum_target_species_count": -1}, "must not be negative"),
        ({"maximum_outside_species_count": -1}, "must not be negative"),
        ({"maximum_unresolved_species_count": -1}, "must not be negative"),
        ({"minimum_target_coverage": -0.1}, "between zero"),
        ({"minimum_mapped_purity": 1.1}, "between zero"),
        ({"maximum_q_value": 1.1}, "between zero"),
        ({"minimum_odds_ratio": -1.0}, "must not be negative"),
        ({"page_size": 0}, "page_size"),
        ({"page_number": 0}, "page_number"),
    ],
)
def test_invalid_taxonomy_search_controls_fail_explicitly(
    values: dict[str, object], message: str
) -> None:
    """Invalid taxonomy controls fail before expensive aggregation."""

    defaults = {
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "target_taxon_id": 10,
    }
    with pytest.raises(InputValidationError, match=message):
        TaxonomySearchFilters(**(defaults | values))


@pytest.mark.parametrize(
    ("records", "fieldnames", "message"),
    [
        ((), None, "non-empty"),
        ((), ("",), "non-empty"),
        ((), ("id", "id"), "unique"),
        (({"id": "a", "extra": "b"},), ("id",), "undeclared"),
    ],
)
def test_records_to_tsv_rejects_ambiguous_schemas(
    records: tuple[dict[str, str], ...],
    fieldnames: tuple[str, ...] | None,
    message: str,
) -> None:
    """Downloads fail visibly when headings cannot represent every record."""

    with pytest.raises(InputValidationError, match=message):
        records_to_tsv(records=records, fieldnames=fieldnames)

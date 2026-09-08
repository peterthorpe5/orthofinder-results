"""Integration tests for bounded application DuckDB queries."""

from __future__ import annotations

from pathlib import Path

import pytest

from orthofinder_interrogation_app.models import GroupKey, GroupSearchFilters
from orthofinder_interrogation_app.queries import OrthoFinderQueryService
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_results.errors import InputValidationError


@pytest.fixture
def query_service(application_resource: Path) -> OrthoFinderQueryService:
    """Return a query service for the application fixture."""

    return OrthoFinderQueryService(resource=open_resource(path=application_resource))


def _group_ids(*, service: OrthoFinderQueryService, filters: GroupSearchFilters) -> list[str]:
    """Return group identifiers from one search page."""

    return [str(row["group_id"]) for row in service.search_groups(filters=filters).rows]


def test_resource_selectors_and_overview_are_complete(
    query_service: OrthoFinderQueryService,
) -> None:
    """Selector values and overview counts come from complete relations."""

    assert query_service.list_species() == (
        "Species_A",
        "Species_B",
        "Species_C",
        "Species_D",
    )
    assert query_service.list_group_types() == ("HOG", "LEGACY_ORTHOGROUP")
    assert query_service.list_hierarchy_nodes() == ("", "N0")
    assert query_service.list_hierarchy_nodes(group_type="HOG") == ("N0",)
    assert query_service.overview_counts() == {
        "group_count": 4,
        "group_species_statistic_count": 8,
        "species_count": 4,
        "distance_group_count": 2,
        "portable_tree_count": 0,
    }


@pytest.mark.parametrize(
    ("include_mode", "included", "excluded", "expected"),
    [
        ("ANY", ("Species_A", "Species_B"), (), ["N0.HOG1", "N0.HOG3", "N0.HOG2"]),
        ("ALL", ("Species_A", "Species_B"), (), ["N0.HOG1", "N0.HOG3"]),
        ("EXACT_SET", ("Species_A", "Species_B"), (), ["N0.HOG1"]),
        ("ANY", ("Species_A",), ("Species_B",), ["N0.HOG2"]),
        ("ALL", (), ("Species_C",), ["N0.HOG1", "OG4"]),
    ],
)
def test_include_and_reject_species_semantics(
    query_service: OrthoFinderQueryService,
    include_mode: str,
    included: tuple[str, ...],
    excluded: tuple[str, ...],
    expected: list[str],
) -> None:
    """ANY, ALL, exact-set and rejection filters retain their stated meanings."""

    assert _group_ids(
        service=query_service,
        filters=GroupSearchFilters(
            included_species=included,
            include_mode=include_mode,
            excluded_species=excluded,
        ),
    ) == expected


def test_identifier_ranges_pagination_and_literal_wildcards(
    query_service: OrthoFinderQueryService,
) -> None:
    """Identifier searches are literal and results remain bounded."""

    assert _group_ids(
        service=query_service,
        filters=GroupSearchFilters(member_id_contains="alpha_2"),
    ) == ["N0.HOG1"]
    assert _group_ids(
        service=query_service,
        filters=GroupSearchFilters(member_id_contains="%"),
    ) == ["N0.HOG2"]
    assert _group_ids(
        service=query_service,
        filters=GroupSearchFilters(group_id_contains="HOG_"),
    ) == []
    first = query_service.search_groups(
        filters=GroupSearchFilters(
            group_type="HOG",
            hierarchy_node="N0",
            minimum_member_count=3,
            maximum_member_count=3,
            minimum_species_count=2,
            maximum_species_count=2,
            sort_mode="GROUP_ID_ASC",
            page_size=1,
        )
    )
    assert first.total_rows == 1
    assert first.rows[0]["group_id"] == "N0.HOG1"
    assert first.page_number == 1 and first.page_size == 1
    empty_second = query_service.search_groups(
        filters=GroupSearchFilters(sort_mode="GROUP_ID_ASC", page_size=2, page_number=3)
    )
    assert empty_second.total_rows == 4 and empty_second.rows == ()
    root_only = query_service.search_groups(
        filters=GroupSearchFilters(hierarchy_node="", sort_mode="GROUP_ID_ASC")
    )
    assert [row["group_id"] for row in root_only.rows] == ["OG4"]


@pytest.mark.parametrize(
    ("filters", "expected"),
    [
        (GroupSearchFilters(distance_availability="CALCULATED"), ["N0.HOG1", "N0.HOG3"]),
        (GroupSearchFilters(distance_availability="NOT_CALCULATED"), ["N0.HOG2", "OG4"]),
        (GroupSearchFilters(maximum_mean_distance=0.2), ["N0.HOG1"]),
        (GroupSearchFilters(maximum_distance_sd=0.1), ["N0.HOG1"]),
        (
            GroupSearchFilters(sort_mode="MEAN_DISTANCE_ASC"),
            ["N0.HOG1", "N0.HOG3", "N0.HOG2", "OG4"],
        ),
        (
            GroupSearchFilters(sort_mode="MEAN_DISTANCE_DESC"),
            ["N0.HOG3", "N0.HOG1", "N0.HOG2", "OG4"],
        ),
        (
            GroupSearchFilters(sort_mode="DISTANCE_SD_ASC"),
            ["N0.HOG1", "N0.HOG3", "N0.HOG2", "OG4"],
        ),
        (
            GroupSearchFilters(sort_mode="SPECIES_COUNT_DESC"),
            ["N0.HOG3", "N0.HOG1", "N0.HOG2", "OG4"],
        ),
    ],
)
def test_distance_filters_and_sorting(
    query_service: OrthoFinderQueryService,
    filters: GroupSearchFilters,
    expected: list[str],
) -> None:
    """Compactness can be filtered and sorted without fabricating absent distances."""

    assert _group_ids(service=query_service, filters=filters) == expected


def test_group_details_memberships_and_copy_counts(
    query_service: OrthoFinderQueryService,
) -> None:
    """A selected composite key retrieves exactly its own complete records."""

    key = GroupKey(
        run_id="test_run", group_type="HOG", hierarchy_node="N0", group_id="N0.HOG1"
    )
    group = query_service.get_group(key=key)
    assert group["mean_distance"] == pytest.approx(0.1)
    assert group["population_stddev_distance"] == pytest.approx(0.02)
    species = query_service.get_group_species(key=key)
    assert [row["species_label"] for row in species] == ["Species_A", "Species_B"]
    assert species[0]["species_member_count"] == 2
    members = query_service.get_group_members(key=key)
    assert [row["member_id"] for row in members] == ["alpha_1", "alpha_2", "beta_1"]
    legacy = GroupKey(
        run_id="test_run",
        group_type="LEGACY_ORTHOGROUP",
        hierarchy_node="",
        group_id="OG4",
    )
    assert query_service.get_group_members(key=legacy)[0]["member_id"] == "delta_1"
    assert len(query_service.get_group_members(key=key, maximum=2)) == 2
    with pytest.raises(InputValidationError, match="maximum must be"):
        query_service.get_group_members(key=key, maximum=0)


def test_missing_groups_unsupported_types_and_query_errors_are_controlled(
    query_service: OrthoFinderQueryService,
) -> None:
    """Invalid detail requests and database errors remain visible."""

    missing = GroupKey(
        run_id="test_run", group_type="HOG", hierarchy_node="N0", group_id="missing"
    )
    with pytest.raises(InputValidationError, match="observed 0"):
        query_service.get_group(key=missing)
    unsupported = GroupKey(
        run_id="test_run", group_type="OTHER", hierarchy_node="", group_id="x"
    )
    with pytest.raises(InputValidationError, match="Unsupported group type"):
        query_service.get_group_members(key=unsupported)
    with pytest.raises(InputValidationError, match="Resource query failed"):
        query_service._query(sql="SELECT * FROM missing_relation", parameters=())

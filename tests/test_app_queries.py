"""Integration tests for bounded application DuckDB queries."""

from __future__ import annotations

from pathlib import Path

import duckdb
import pytest

from orthofinder_interrogation_app.models import (
    DistanceResultFilters,
    GroupKey,
    GroupSearchFilters,
    ProteinSearchFilters,
)
from orthofinder_interrogation_app.queries import OrthoFinderQueryService
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_results.errors import InputValidationError


@pytest.fixture
def query_service(application_resource: Path) -> OrthoFinderQueryService:
    """Return a query service for the application fixture."""

    return OrthoFinderQueryService(resource=open_resource(path=application_resource))


@pytest.fixture
def benchmark_query_service(
    benchmark_application_resource: Path,
) -> OrthoFinderQueryService:
    """Return a service backed by the shared benchmark application fixture."""

    return OrthoFinderQueryService(
        resource=open_resource(path=benchmark_application_resource)
    )


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


def test_benchmark_queries_are_capability_driven_and_parameterised(
    benchmark_query_service: OrthoFinderQueryService,
) -> None:
    """Profile, contrast and individual result queries retain cluster units."""

    service = benchmark_query_service
    profiles = service.benchmark_profiles()
    assert {row["profile_id"] for row in profiles} == {
        "E3_ALL",
        "HOUSEKEEPING_ALL",
        "MATCHED_NON_FOCUS",
    }
    assert all(row["group_count"] == 1 for row in profiles)
    catalogue = service.benchmark_cluster_catalogue()
    assert [row["group_id"] for row in catalogue] == ["N0.HOG1", "N0.HOG3"]
    raw = service.benchmark_distribution(
        profile_ids=("E3_ALL",),
        metric="mean_distance",
        comparison_scale="RAW",
    )
    assert float(raw[0]["metric_value"]) == pytest.approx(0.2)
    residual = service.benchmark_distribution(
        profile_ids=("E3_ALL",),
        metric="mean_distance",
        comparison_scale="MATCHED_RESIDUAL",
    )
    assert float(residual[0]["metric_value"]) == pytest.approx(0.08)
    contrasts = service.benchmark_contrasts(metric="mean_distance")
    assert float(contrasts[0]["fdr_q_value"]) == pytest.approx(0.04)
    key = GroupKey("test_run", "HOG", "N0", "N0.HOG1")
    individual = service.benchmark_individual_comparisons(
        key=key,
        metric="mean_distance",
    )
    assert len(individual) == 2
    classifications = service.benchmark_classifications()
    assert classifications[0]["status"] == "CLASSIFIED"
    assert classifications[0]["profile_ids"] == "E3_ALL"
    housekeeping = service.benchmark_profile_markers(
        profile_ids=("HOUSEKEEPING_ALL",)
    )
    assert housekeeping[0]["marker_id"] == "AT1G01010"
    assert housekeeping[0]["matched_member_id"] == "a3"
    assert housekeeping[0]["marker_name"] == "Housekeeping fixture protein"
    e3 = service.benchmark_profile_markers(profile_ids=("E3_ALL",))
    assert e3[0]["marker_id"] == "seedE3"
    assert e3[0]["matched_member_id"] == "alpha_1"
    assert service.benchmark_profile_markers(profile_ids=()) == ()
    with pytest.raises(InputValidationError, match="sequence of exact text"):
        service.benchmark_profile_markers(
            profile_ids="E3_ALL",  # type: ignore[arg-type]
        )
    with pytest.raises(InputValidationError, match="maximum must be"):
        service.benchmark_profile_markers(profile_ids=("E3_ALL",), maximum=0)
    with pytest.raises(InputValidationError, match="maximum must be an integer"):
        service.benchmark_profile_markers(profile_ids=("E3_ALL",), maximum=True)
    with pytest.raises(InputValidationError, match="contain 2 marker-to-cluster"):
        service.benchmark_profile_markers(
            profile_ids=("E3_ALL", "HOUSEKEEPING_ALL"), maximum=1
        )
    with pytest.raises(InputValidationError, match="Unsupported benchmark metric"):
        service.benchmark_distribution(
            profile_ids=("E3_ALL",),
            metric="unsafe_sql",
            comparison_scale="RAW",
        )
    with pytest.raises(InputValidationError, match="comparison scale"):
        service.benchmark_distribution(
            profile_ids=("E3_ALL",),
            metric="mean_distance",
            comparison_scale="INVALID",
        )
    assert service.benchmark_distribution(
        profile_ids=(),
        metric="mean_distance",
        comparison_scale="RAW",
    ) == ()
    with pytest.raises(InputValidationError, match="Unsupported benchmark metric"):
        service.benchmark_contrasts(metric="unsafe_sql")
    with pytest.raises(InputValidationError, match="Unsupported benchmark metric"):
        service.benchmark_individual_comparisons(key=key, metric="unsafe_sql")


def test_older_resources_return_empty_benchmark_capabilities(
    query_service: OrthoFinderQueryService,
) -> None:
    """Pre-benchmark resources remain usable and expose explicit empty results."""

    assert query_service.benchmark_profiles() == ()
    assert query_service.benchmark_cluster_catalogue() == ()
    assert query_service.benchmark_profile_markers(profile_ids=("E3_ALL",)) == ()
    assert query_service.benchmark_contrasts() == ()
    assert query_service.benchmark_classifications() == ()


def test_benchmark_marker_sources_are_independently_optional(
    benchmark_application_resource: Path,
) -> None:
    """Reference-marker browsing remains useful when an E3 relation is absent."""

    database = benchmark_application_resource / "duckdb/orthofinder_results.duckdb"
    connection = duckdb.connect(str(database))
    try:
        connection.execute("DROP TABLE e3_seed_matches")
    finally:
        connection.close()
    service = OrthoFinderQueryService(
        resource=open_resource(path=benchmark_application_resource)
    )
    markers = service.benchmark_profile_markers(profile_ids=("HOUSEKEEPING_ALL",))
    assert markers[0]["marker_id"] == "AT1G01010"


def test_benchmark_marker_query_accepts_e3_only_and_no_marker_relations(
    benchmark_application_resource: Path,
) -> None:
    """E3-only resources work and a missing marker authority returns an empty result."""

    database = benchmark_application_resource / "duckdb/orthofinder_results.duckdb"
    connection = duckdb.connect(str(database))
    try:
        connection.execute("DROP TABLE benchmark_marker_matches")
    finally:
        connection.close()
    service = OrthoFinderQueryService(
        resource=open_resource(path=benchmark_application_resource)
    )
    assert service.benchmark_profile_markers(profile_ids=("E3_ALL",))[0][
        "marker_id"
    ] == "seedE3"
    connection = duckdb.connect(str(database))
    try:
        connection.execute("DROP TABLE e3_seed_matches")
    finally:
        connection.close()
    service = OrthoFinderQueryService(
        resource=open_resource(path=benchmark_application_resource)
    )
    assert service.benchmark_profile_markers(profile_ids=("E3_ALL",)) == ()


def test_all_distance_result_facets_and_selected_statistics(
    query_service: OrthoFinderQueryService,
) -> None:
    """The complete result authority exposes filters and selected distance fields."""

    assert query_service.distance_result_facets() == {
        "group_types": ("HOG",),
        "hierarchy_nodes": ("N0",),
        "distance_methods": ("patristic_branch_length",),
        "computation_statuses": ("EXACT",),
    }
    result = query_service.distance_results(
        filters=DistanceResultFilters(),
        columns=(
            "group_id",
            "mean_distance",
            "sampling_fraction",
            "full_group_matrix",
            "interquartile_range",
            "relative_distance_spread",
        ),
    )
    assert result.total_rows == 2
    assert result.columns[0] == "group_id"
    assert [row["group_id"] for row in result.rows] == ["N0.HOG1", "N0.HOG3"]
    first = result.rows[0]
    assert first["mean_distance"] == pytest.approx(0.1)
    assert first["sampling_fraction"] == pytest.approx(1.0)
    assert first["full_group_matrix"] is True
    assert first["interquartile_range"] == pytest.approx(0.1)
    assert first["relative_distance_spread"] == pytest.approx(0.2)


def test_all_distance_result_filters_and_bounds_are_defensive(
    query_service: OrthoFinderQueryService,
) -> None:
    """Exact filters work and invalid or excessive exports fail before materialisation."""

    exact = query_service.distance_results(
        filters=DistanceResultFilters(
            group_type=" HOG ",
            hierarchy_node=" N0 ",
            distance_method=" patristic_branch_length ",
            computation_status=" EXACT ",
        ),
        columns=("group_id", "distance_method"),
    )
    assert exact.total_rows == 2
    empty = query_service.distance_results(
        filters=DistanceResultFilters(distance_method="missing"),
        columns=("group_id",),
    )
    assert empty.total_rows == 0 and empty.rows == ()
    with pytest.raises(InputValidationError, match="at least one"):
        query_service.distance_results(
            filters=DistanceResultFilters(),
            columns=(),
        )
    with pytest.raises(InputValidationError, match="must be unique"):
        query_service.distance_results(
            filters=DistanceResultFilters(),
            columns=("group_id", "group_id"),
        )
    with pytest.raises(InputValidationError, match="Unsupported all-results"):
        query_service.distance_results(
            filters=DistanceResultFilters(),
            columns=("unsafe_sql",),
        )
    with pytest.raises(InputValidationError, match="contains 2 groups"):
        query_service.distance_results(
            filters=DistanceResultFilters(),
            columns=("group_id",),
            maximum=1,
        )
    with pytest.raises(InputValidationError, match="maximum must be"):
        query_service.distance_results(
            filters=DistanceResultFilters(),
            columns=("group_id",),
            maximum=0,
        )


def test_all_distance_result_filter_models_reject_invalid_text() -> None:
    """NULs and non-text values cannot reach result-query parameters."""

    assert DistanceResultFilters(hierarchy_node="").hierarchy_node == ""
    assert DistanceResultFilters(hierarchy_node=None).hierarchy_node is None
    with pytest.raises(InputValidationError, match="group_type must be text"):
        DistanceResultFilters(group_type=1)  # type: ignore[arg-type]
    with pytest.raises(InputValidationError, match="NUL"):
        DistanceResultFilters(distance_method="bad\x00method")
    with pytest.raises(InputValidationError, match="text or None"):
        DistanceResultFilters(hierarchy_node=1)  # type: ignore[arg-type]


def test_all_distance_results_prefer_complete_then_apply_scope_filter(
    application_resource: Path,
) -> None:
    """One row per group prefers completeness unless a sampled scope is requested."""

    database = application_resource / "duckdb" / "orthofinder_results.duckdb"
    connection = duckdb.connect(str(database))
    try:
        connection.execute(
            "INSERT INTO distance_statistics SELECT run_id, group_type, hierarchy_node, "
            "group_id, distance_method, 'DETERMINISTIC_MEMBER_SAMPLE', "
            "member_identifier_resolution, total_member_count, 2, 1, "
            "unresolved_pair_count, minimum_distance, q05_distance, q25_distance, "
            "median_distance, 9.0, q75_distance, q95_distance, maximum_distance, "
            "population_stddev_distance, mean_comparable_sites, 'sample.tsv', "
            "failure_reason FROM distance_statistics WHERE group_id = 'N0.HOG1'"
        )
    finally:
        connection.close()
    service = OrthoFinderQueryService(resource=open_resource(path=application_resource))
    preferred = service.distance_results(
        filters=DistanceResultFilters(),
        columns=("group_id", "computation_status", "sampled_member_count", "mean_distance"),
    )
    first = next(row for row in preferred.rows if row["group_id"] == "N0.HOG1")
    assert first["computation_status"] == "EXACT"
    assert first["sampled_member_count"] == 3
    sampled = service.distance_results(
        filters=DistanceResultFilters(
            computation_status="DETERMINISTIC_MEMBER_SAMPLE"
        ),
        columns=("group_id", "computation_status", "mean_distance"),
    )
    assert sampled.total_rows == 1
    assert sampled.rows[0]["mean_distance"] == pytest.approx(9.0)


def test_protein_search_finds_canonical_and_internal_identifiers(
    query_service: OrthoFinderQueryService,
) -> None:
    """Protein searches retain exact clusters and resolve internal aliases."""

    canonical = query_service.search_proteins(
        filters=ProteinSearchFilters(query=" alpha_1 ")
    )
    assert canonical.total_rows == 1
    assert canonical.truncated is False
    assert canonical.rows[0]["member_id"] == "alpha_1"
    assert canonical.rows[0]["internal_id"] == "0_0"
    assert canonical.rows[0]["match_source"] == "PROTEIN_ID"
    assert canonical.rows[0]["group_id"] == "N0.HOG1"
    assert canonical.rows[0]["mean_distance"] == pytest.approx(0.1)

    alias = query_service.search_proteins(
        filters=ProteinSearchFilters(query="0_0")
    )
    assert alias.total_rows == 1
    assert alias.rows[0]["matched_identifier"] == "0_0"
    assert alias.rows[0]["member_id"] == "alpha_1"
    assert alias.rows[0]["match_source"] == "ORTHOFINDER_INTERNAL_ID"


def test_protein_contains_search_is_literal_bounded_and_filterable(
    query_service: OrthoFinderQueryService,
) -> None:
    """Literal fragments, hierarchy duplication and result bounds remain explicit."""

    contained = query_service.search_proteins(
        filters=ProteinSearchFilters(
            query="alpha",
            match_mode="CONTAINS",
            group_type="HOG",
            maximum_rows=1,
        )
    )
    assert contained.total_rows == 2
    assert len(contained.rows) == 1
    assert contained.truncated is True
    literal = query_service.search_proteins(
        filters=ProteinSearchFilters(query="%me", match_mode="CONTAINS")
    )
    assert [row["member_id"] for row in literal.rows] == ["literal%member"]
    legacy = query_service.search_proteins(
        filters=ProteinSearchFilters(
            query="delta_1",
            group_type="LEGACY_ORTHOGROUP",
        )
    )
    assert legacy.rows[0]["hierarchy_node"] == ""
    with pytest.raises(InputValidationError, match="Unsupported protein group"):
        query_service.search_proteins(
            filters=ProteinSearchFilters(query="alpha_1", group_type="OTHER")
        )


def test_protein_search_without_sequence_alias_relation(
    application_resource: Path,
) -> None:
    """Canonical identifiers remain searchable in resources without alias rows."""

    database = application_resource / "duckdb" / "orthofinder_results.duckdb"
    connection = duckdb.connect(str(database))
    try:
        connection.execute("DROP TABLE sequences")
    finally:
        connection.close()
    service = OrthoFinderQueryService(resource=open_resource(path=application_resource))
    canonical = service.search_proteins(
        filters=ProteinSearchFilters(query="alpha_1")
    )
    assert canonical.total_rows == 1
    assert canonical.rows[0]["internal_id"] == ""
    alias = service.search_proteins(filters=ProteinSearchFilters(query="0_0"))
    assert alias.total_rows == 0


@pytest.mark.parametrize(
    ("query", "match_mode", "expected_source", "expected_identifier"),
    [
        ("Q9SA03", "EXACT", "UNIPROT_ACCESSION", "Q9SA03"),
        ("FB27_ARATH", "EXACT", "UNIPROT_ENTRY", "FB27_ARATH"),
        ("q9s", "CONTAINS", "UNIPROT_ACCESSION", "Q9SA03"),
        ("arath", "CONTAINS", "UNIPROT_ENTRY", "FB27_ARATH"),
    ],
)
def test_protein_search_resolves_controlled_uniprot_aliases(
    application_resource: Path,
    query: str,
    match_mode: str,
    expected_source: str,
    expected_identifier: str,
) -> None:
    """UniProt pipe identifiers expose bounded accession and entry aliases."""

    database = application_resource / "duckdb" / "orthofinder_results.duckdb"
    connection = duckdb.connect(str(database))
    try:
        connection.execute(
            "UPDATE hog_memberships SET member_id = 'sp|Q9SA03|FB27_ARATH' "
            "WHERE member_id = 'alpha_1'"
        )
        connection.execute(
            "UPDATE sequences SET member_id = 'sp|Q9SA03|FB27_ARATH' "
            "WHERE member_id = 'alpha_1'"
        )
    finally:
        connection.close()
    service = OrthoFinderQueryService(resource=open_resource(path=application_resource))
    result = service.search_proteins(
        filters=ProteinSearchFilters(query=query, match_mode=match_mode)
    )
    assert result.total_rows == 1
    assert result.rows[0]["member_id"] == "sp|Q9SA03|FB27_ARATH"
    assert result.rows[0]["match_source"] == expected_source
    assert result.rows[0]["matched_identifier"] == expected_identifier


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"query": ""}, "Enter a protein"),
        ({"query": "ab", "match_mode": "CONTAINS"}, "at least three"),
        ({"query": "a", "match_mode": "UNKNOWN"}, "Unsupported protein"),
        ({"query": "bad\nvalue"}, "control character"),
        ({"query": "x" * 513}, "512"),
        ({"query": "a", "group_type": 1}, "group_type must be text"),
        ({"query": "a", "maximum_rows": 0}, "between 1 and 1,000"),
        ({"query": "a", "maximum_rows": True}, "must be an integer"),
    ],
)
def test_protein_search_filter_validation(
    kwargs: dict[str, object], message: str
) -> None:
    """Malformed and unbounded protein searches fail before reaching DuckDB."""

    with pytest.raises(InputValidationError, match=message):
        ProteinSearchFilters(**kwargs)  # type: ignore[arg-type]


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

"""Bounded, parameterised DuckDB queries for the interactive application."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import Any

import duckdb

from orthofinder_results.errors import InputValidationError

from .models import (
    GroupKey,
    GroupSearchFilters,
    ResourceIdentity,
    SearchPage,
    TaxonomySearchFilters,
    TaxonomySearchPage,
)
from .resource import connect_read_only
from .taxonomy import TaxonomyAuthority, add_fisher_enrichment

_LOGGER = logging.getLogger("orthofinder_interrogation_app.queries")
MAX_GROUP_MEMBER_ROWS = 50_000
MAX_TAXONOMY_TEST_GROUPS = 250_000
_MEMBERSHIP_RELATIONS = {
    "HOG": "hog_memberships",
    "LEGACY_ORTHOGROUP": "legacy_orthogroup_memberships",
}
_SORT_SQL = {
    "MEMBER_COUNT_DESC": "member_count DESC, group_type, hierarchy_node, group_id",
    "SPECIES_COUNT_DESC": "species_count DESC, member_count DESC, group_id",
    "GROUP_ID_ASC": "group_type, hierarchy_node, group_id",
    "MEAN_DISTANCE_ASC": (
        "mean_distance NULLS LAST, member_count DESC, group_type, hierarchy_node, group_id"
    ),
    "MEAN_DISTANCE_DESC": (
        "mean_distance DESC NULLS LAST, member_count DESC, group_type, hierarchy_node, group_id"
    ),
    "DISTANCE_SD_ASC": (
        "population_stddev_distance NULLS LAST, member_count DESC, group_type, "
        "hierarchy_node, group_id"
    ),
}
_GROUP_COLUMNS = (
    "run_id",
    "group_type",
    "hierarchy_node",
    "group_id",
    "legacy_orthogroup_id",
    "gene_tree_parent_clade",
    "member_count",
    "species_count",
    "single_copy_species_count",
    "max_copies_per_species",
    "mean_copies_per_species",
    "is_singleton",
    "distance_method",
    "computation_status",
    "total_member_count",
    "sampled_member_count",
    "distance_pair_count",
    "minimum_distance",
    "q25_distance",
    "median_distance",
    "mean_distance",
    "q75_distance",
    "maximum_distance",
    "population_stddev_distance",
    "failure_reason",
)


class OrthoFinderQueryService:
    """Execute safe read-only queries against one validated resource."""

    def __init__(self, *, resource: ResourceIdentity) -> None:
        """Initialise the service without retaining an open connection.

        Args:
            resource: Validated immutable resource identity.
        """

        self.resource = resource

    def list_species(self) -> tuple[str, ...]:
        """Return every exact stored species label in lexical order."""

        rows = self._query(
            sql="SELECT DISTINCT species_label FROM species ORDER BY species_label",
            parameters=(),
        )
        return tuple(str(row["species_label"]) for row in rows)

    def list_group_types(self) -> tuple[str, ...]:
        """Return available group authority types."""

        rows = self._query(
            sql="SELECT DISTINCT group_type FROM group_statistics ORDER BY group_type",
            parameters=(),
        )
        return tuple(str(row["group_type"]) for row in rows)

    def list_hierarchy_nodes(self, *, group_type: str = "") -> tuple[str, ...]:
        """Return available hierarchy nodes, optionally for one group type.

        Args:
            group_type: Optional exact group type.

        Returns:
            Deterministically ordered hierarchy nodes.
        """

        if group_type:
            rows = self._query(
                sql=(
                    "SELECT DISTINCT hierarchy_node FROM group_statistics "
                    "WHERE group_type = ? ORDER BY hierarchy_node"
                ),
                parameters=(group_type,),
            )
        else:
            rows = self._query(
                sql=(
                    "SELECT DISTINCT hierarchy_node FROM group_statistics ORDER BY hierarchy_node"
                ),
                parameters=(),
            )
        return tuple(str(row["hierarchy_node"]) for row in rows)

    def overview_counts(self) -> dict[str, int]:
        """Return complete-authority counts for the overview page."""

        rows = self._query(
            sql=(
                "SELECT (SELECT count(*) FROM group_statistics) AS group_count, "
                "(SELECT count(*) FROM group_species_statistics) "
                "AS group_species_statistic_count, "
                "(SELECT count(*) FROM species) AS species_count, "
                "(SELECT count(*) FROM distance_statistics "
                "WHERE distance_pair_count > 0) AS distance_group_count"
            ),
            parameters=(),
        )
        return {key: int(value) for key, value in rows[0].items()}

    def search_groups(self, *, filters: GroupSearchFilters) -> SearchPage:
        """Return one bounded page matching exact search semantics.

        Args:
            filters: Validated group-search controls.

        Returns:
            Search rows with the complete matching row count.
        """

        query_sql, parameters = _group_search_query(filters=filters)
        count_rows = self._query(
            sql=f"SELECT count(*) AS total_rows FROM ({query_sql}) AS matched_groups",
            parameters=parameters,
        )
        order_sql = _SORT_SQL[filters.sort_mode]
        page_rows = self._query(
            sql=(
                f"SELECT * FROM ({query_sql}) AS matched_groups ORDER BY {order_sql} "
                "LIMIT ? OFFSET ?"
            ),
            parameters=(*parameters, filters.page_size, filters.offset),
        )
        total_rows = int(count_rows[0]["total_rows"])
        _LOGGER.info(
            "Group search completed: run=%s, matches=%s, returned=%s, page=%s",
            self.resource.run_id,
            total_rows,
            len(page_rows),
            filters.page_number,
        )
        return SearchPage(
            rows=tuple(page_rows),
            total_rows=total_rows,
            page_number=filters.page_number,
            page_size=filters.page_size,
        )

    def get_group(self, *, key: GroupKey) -> dict[str, Any]:
        """Return complete group and distance-summary fields for one key.

        Args:
            key: Composite run-scoped group identity.

        Returns:
            Exactly one group record.

        Raises:
            InputValidationError: If the group does not exist or is duplicated.
        """

        filters = GroupSearchFilters(
            group_type=key.group_type,
            hierarchy_node=key.hierarchy_node,
            group_id_contains=key.group_id,
            sort_mode="GROUP_ID_ASC",
            page_size=2,
        )
        query_sql, parameters = _group_search_query(
            filters=filters,
            exact_group_id=key.group_id,
            run_id=key.run_id,
        )
        rows = self._query(sql=query_sql, parameters=parameters)
        if len(rows) != 1:
            raise InputValidationError(
                f"Expected exactly one group for {key.display_label()}; observed {len(rows)}."
            )
        return rows[0]

    def get_group_species(self, *, key: GroupKey) -> tuple[dict[str, Any], ...]:
        """Return complete per-species copy counts for one group."""

        rows = self._query(
            sql=(
                "SELECT species_label, species_member_count, member_fraction "
                "FROM group_species_statistics WHERE run_id = ? AND group_type = ? "
                "AND hierarchy_node = ? AND group_id = ? "
                "ORDER BY species_member_count DESC, species_label"
            ),
            parameters=_key_parameters(key=key),
        )
        return tuple(rows)

    def get_group_members(
        self, *, key: GroupKey, maximum: int = MAX_GROUP_MEMBER_ROWS
    ) -> tuple[dict[str, Any], ...]:
        """Return a bounded, deterministic set of exact group memberships.

        Args:
            key: Composite run-scoped group identity.
            maximum: Maximum rows to materialise in the browser process.

        Returns:
            At most ``maximum`` membership records.

        Raises:
            InputValidationError: If the group type or row bound is unsupported.
        """

        relation = _MEMBERSHIP_RELATIONS.get(key.group_type)
        if relation is None:
            raise InputValidationError(f"Unsupported group type: {key.group_type}")
        if not 1 <= maximum <= MAX_GROUP_MEMBER_ROWS:
            raise InputValidationError(f"maximum must be between 1 and {MAX_GROUP_MEMBER_ROWS:,}.")
        rows = self._query(
            sql=(
                "SELECT species_label, member_id, legacy_orthogroup_id, "
                "gene_tree_parent_clade FROM "
                f"{relation} WHERE run_id = ? AND group_type = ? "
                "AND hierarchy_node = ? AND group_id = ? "
                "ORDER BY species_label, member_id LIMIT ?"
            ),
            parameters=(*_key_parameters(key=key), maximum),
        )
        return tuple(rows)

    def search_taxonomy_groups(
        self,
        *,
        authority: TaxonomyAuthority,
        filters: TaxonomySearchFilters,
    ) -> TaxonomySearchPage:
        """Search groups using reviewed descendant membership and explicit scope.

        Args:
            authority: Mapping decisions validated against this resource's species.
            filters: Exact group authority, target taxon and search semantics.

        Returns:
            Bounded result page with mapping and statistical-universe counts.
        """

        resource_species = self.list_species()
        if authority.expected_species != resource_species:
            raise InputValidationError(
                "Taxonomy mapping was not validated against this resource's exact species set."
            )
        target_species = authority.target_species(taxon_id=filters.target_taxon_id)
        if not target_species:
            raise InputValidationError(
                f"No reviewed sampled species descends from NCBI taxon {filters.target_taxon_id}."
            )
        reviewed = set(authority.reviewed_species)
        target = set(target_species)
        outside_species = tuple(sorted(reviewed.difference(target)))
        unresolved_species = authority.unresolved_species
        categories = tuple(
            [(species, "TARGET") for species in sorted(target)]
            + [(species, "OUTSIDE") for species in outside_species]
            + [(species, "UNRESOLVED") for species in unresolved_species]
        )
        summary_sql, summary_parameters = _taxonomy_summary_query(
            run_id=self.resource.run_id,
            filters=filters,
            categories=categories,
            target_species_count=len(target_species),
        )
        if filters.mode == "ENRICHED":
            page_rows, total_rows, tested_count = self._search_enriched_groups(
                summary_sql=summary_sql,
                summary_parameters=summary_parameters,
                filters=filters,
                total_target_species=len(target_species),
                total_outside_species=len(outside_species),
            )
        else:
            condition_sql, condition_parameters = _taxonomy_filter_condition(filters=filters)
            filtered_sql = (
                f"SELECT * FROM ({summary_sql}) AS taxonomy_summary WHERE {condition_sql}"
            )
            count_rows = self._query(
                sql=f"SELECT count(*) AS total_rows FROM ({filtered_sql}) AS matched",
                parameters=(*summary_parameters, *condition_parameters),
            )
            total_rows = int(count_rows[0]["total_rows"])
            page_rows = tuple(
                self._query(
                    sql=(
                        f"{filtered_sql} ORDER BY target_coverage DESC, "
                        "mapped_target_fraction DESC, target_species_count DESC, "
                        "member_count DESC, group_id LIMIT ? OFFSET ?"
                    ),
                    parameters=(
                        *summary_parameters,
                        *condition_parameters,
                        filters.page_size,
                        filters.offset,
                    ),
                )
            )
            tested_count = 0
        _LOGGER.info(
            "Taxonomy search completed: run=%s, mode=%s, taxon=%s, matches=%s, "
            "returned=%s, tested=%s",
            self.resource.run_id,
            filters.mode,
            filters.target_taxon_id,
            total_rows,
            len(page_rows),
            tested_count,
        )
        return TaxonomySearchPage(
            rows=page_rows,
            total_rows=total_rows,
            tested_group_count=tested_count,
            target_species_count=len(target_species),
            outside_species_count=len(outside_species),
            unresolved_species_count=len(unresolved_species),
            page_number=filters.page_number,
            page_size=filters.page_size,
        )

    def _search_enriched_groups(
        self,
        *,
        summary_sql: str,
        summary_parameters: tuple[object, ...],
        filters: TaxonomySearchFilters,
        total_target_species: int,
        total_outside_species: int,
    ) -> tuple[tuple[dict[str, Any], ...], int, int]:
        """Calculate bounded Fisher tests and adjusted enrichment results."""

        count_rows = self._query(
            sql=f"SELECT count(*) AS total_rows FROM ({summary_sql}) AS tested",
            parameters=summary_parameters,
        )
        tested_count = int(count_rows[0]["total_rows"])
        if tested_count > MAX_TAXONOMY_TEST_GROUPS:
            raise InputValidationError(
                f"Enrichment authority contains {tested_count:,} groups; the interactive "
                f"limit is {MAX_TAXONOMY_TEST_GROUPS:,}. Select a narrower hierarchy node."
            )
        rows = tuple(
            self._query(
                sql=f"{summary_sql} ORDER BY group_id",
                parameters=summary_parameters,
            )
        )
        tested = add_fisher_enrichment(
            rows=rows,
            total_target_species=total_target_species,
            total_outside_species=total_outside_species,
        )
        matched = [
            row
            for row in tested
            if int(row["target_species_count"]) >= filters.minimum_target_species_count
            and float(row["target_coverage"]) >= filters.minimum_target_coverage
            and float(row["mapped_target_fraction"]) >= filters.minimum_mapped_purity
            and float(row["enrichment_q_value"]) <= filters.maximum_q_value
            and float(row["enrichment_odds_ratio"]) >= filters.minimum_odds_ratio
        ]
        matched.sort(
            key=lambda row: (
                float(row["enrichment_q_value"]),
                float(row["enrichment_p_value"]),
                -float(row["enrichment_odds_ratio"]),
                str(row["group_id"]),
            )
        )
        start, finish = filters.offset, filters.offset + filters.page_size
        return tuple(matched[start:finish]), len(matched), tested_count

    def _query(self, *, sql: str, parameters: Sequence[object]) -> list[dict[str, Any]]:
        """Execute one parameterised query and close its read-only connection."""

        connection = connect_read_only(database_path=self.resource.database_path)
        try:
            cursor = connection.execute(sql, list(parameters))
            headings = tuple(column[0] for column in cursor.description)
            return [dict(zip(headings, row, strict=True)) for row in cursor.fetchall()]
        except duckdb.Error as error:
            _LOGGER.exception("Read-only DuckDB query failed")
            raise InputValidationError(f"Resource query failed: {error}") from error
        finally:
            connection.close()


def _group_search_query(
    *,
    filters: GroupSearchFilters,
    exact_group_id: str = "",
    run_id: str = "",
) -> tuple[str, tuple[object, ...]]:
    """Build one parameterised, one-row-per-group search relation."""

    conditions: list[str] = []
    parameters: list[object] = []
    if run_id:
        conditions.append("g.run_id = ?")
        parameters.append(run_id)
    if filters.group_type:
        conditions.append("g.group_type = ?")
        parameters.append(filters.group_type)
    if filters.hierarchy_node is not None:
        conditions.append("g.hierarchy_node = ?")
        parameters.append(filters.hierarchy_node)
    if exact_group_id:
        conditions.append("g.group_id = ?")
        parameters.append(exact_group_id)
    elif filters.group_id_contains.strip():
        conditions.append("g.group_id ILIKE ? ESCAPE '\\'")
        parameters.append(_literal_contains(value=filters.group_id_contains))
    if filters.member_id_contains.strip():
        conditions.append(_member_exists_condition())
        parameters.extend(
            (
                _literal_contains(value=filters.member_id_contains),
                _literal_contains(value=filters.member_id_contains),
            )
        )
    _append_species_conditions(
        conditions=conditions,
        parameters=parameters,
        filters=filters,
    )
    _append_range_condition(
        conditions=conditions,
        parameters=parameters,
        column="g.member_count",
        minimum=filters.minimum_member_count,
        maximum=filters.maximum_member_count,
    )
    _append_range_condition(
        conditions=conditions,
        parameters=parameters,
        column="g.species_count",
        minimum=filters.minimum_species_count,
        maximum=filters.maximum_species_count,
    )
    if filters.distance_availability == "CALCULATED":
        conditions.append("d.distance_pair_count > 0")
    elif filters.distance_availability == "NOT_CALCULATED":
        conditions.append("coalesce(d.distance_pair_count, 0) = 0")
    if filters.maximum_mean_distance is not None:
        conditions.append("d.mean_distance <= ?")
        parameters.append(filters.maximum_mean_distance)
    if filters.maximum_distance_sd is not None:
        conditions.append("d.population_stddev_distance <= ?")
        parameters.append(filters.maximum_distance_sd)
    where_sql = " AND ".join(conditions) if conditions else "TRUE"
    selected = ", ".join(
        f"g.{column}" if column in _GROUP_COLUMNS[:12] else f"d.{column}"
        for column in _GROUP_COLUMNS
    )
    sql = (
        "WITH preferred_distance AS ("
        "SELECT *, row_number() OVER (PARTITION BY run_id, group_type, hierarchy_node, "
        "group_id ORDER BY CASE WHEN distance_pair_count > 0 THEN 0 ELSE 1 END, "
        "distance_method) AS app_distance_rank FROM distance_statistics"
        ") SELECT "
        f"{selected} FROM group_statistics AS g LEFT JOIN preferred_distance AS d "
        "ON d.run_id = g.run_id AND d.group_type = g.group_type "
        "AND d.hierarchy_node = g.hierarchy_node AND d.group_id = g.group_id "
        f"AND d.app_distance_rank = 1 WHERE {where_sql}"
    )
    return sql, tuple(parameters)


def _append_species_conditions(
    *,
    conditions: list[str],
    parameters: list[object],
    filters: GroupSearchFilters,
) -> None:
    """Append exact species-set containment and rejection conditions."""

    if filters.included_species:
        placeholders = ", ".join("?" for _ in filters.included_species)
        count_sql = (
            "(SELECT count(DISTINCT inc.species_label) FROM group_species_statistics AS inc "
            "WHERE inc.run_id = g.run_id AND inc.group_type = g.group_type "
            "AND inc.hierarchy_node = g.hierarchy_node AND inc.group_id = g.group_id "
            f"AND inc.species_label IN ({placeholders}))"
        )
        if filters.include_mode == "ANY":
            conditions.append(f"{count_sql} >= 1")
        else:
            conditions.append(f"{count_sql} = ?")
        parameters.extend(filters.included_species)
        if filters.include_mode != "ANY":
            parameters.append(len(filters.included_species))
        if filters.include_mode == "EXACT_SET":
            conditions.append("g.species_count = ?")
            parameters.append(len(filters.included_species))
    if filters.excluded_species:
        placeholders = ", ".join("?" for _ in filters.excluded_species)
        conditions.append(
            "NOT EXISTS (SELECT 1 FROM group_species_statistics AS exc "
            "WHERE exc.run_id = g.run_id AND exc.group_type = g.group_type "
            "AND exc.hierarchy_node = g.hierarchy_node AND exc.group_id = g.group_id "
            f"AND exc.species_label IN ({placeholders}))"
        )
        parameters.extend(filters.excluded_species)


def _member_exists_condition() -> str:
    """Return the parameterised membership identifier search condition."""

    return (
        "((g.group_type = 'HOG' AND EXISTS (SELECT 1 FROM hog_memberships AS hm "
        "WHERE hm.run_id = g.run_id AND hm.group_type = g.group_type "
        "AND hm.hierarchy_node = g.hierarchy_node AND hm.group_id = g.group_id "
        "AND hm.member_id ILIKE ? ESCAPE '\\')) OR "
        "(g.group_type = 'LEGACY_ORTHOGROUP' AND EXISTS (SELECT 1 "
        "FROM legacy_orthogroup_memberships AS lm WHERE lm.run_id = g.run_id "
        "AND lm.group_type = g.group_type AND lm.hierarchy_node = g.hierarchy_node "
        "AND lm.group_id = g.group_id AND lm.member_id ILIKE ? ESCAPE '\\')))"
    )


def _append_range_condition(
    *,
    conditions: list[str],
    parameters: list[object],
    column: str,
    minimum: int | None,
    maximum: int | None,
) -> None:
    """Append optional inclusive SQL bounds."""

    if minimum is not None:
        conditions.append(f"{column} >= ?")
        parameters.append(minimum)
    if maximum is not None:
        conditions.append(f"{column} <= ?")
        parameters.append(maximum)


def _literal_contains(*, value: str) -> str:
    """Escape SQL wildcard characters and return a literal contains pattern."""

    escaped = value.strip().replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
    return f"%{escaped}%"


def _key_parameters(*, key: GroupKey) -> tuple[str, str, str, str]:
    """Return one composite key in stable SQL parameter order."""

    return key.run_id, key.group_type, key.hierarchy_node, key.group_id


def _taxonomy_summary_query(
    *,
    run_id: str,
    filters: TaxonomySearchFilters,
    categories: tuple[tuple[str, str], ...],
    target_species_count: int,
) -> tuple[str, tuple[object, ...]]:
    """Build one parameterised per-group taxonomic composition relation."""

    if not categories or target_species_count < 1:
        raise InputValidationError("Taxonomy search requires mapped sampled species.")
    values_sql = ", ".join("(?, ?)" for _ in categories)
    parameters: list[object] = [value for row in categories for value in row]
    parameters.extend((run_id, filters.group_type, filters.hierarchy_node, target_species_count))
    category = "coalesce(tm.app_category, 'UNRESOLVED')"
    sql = (
        f"WITH taxonomy_map(species_label, app_category) AS (VALUES {values_sql}), "
        "summaries AS (SELECT g.run_id, g.group_type, g.hierarchy_node, g.group_id, "
        "g.legacy_orthogroup_id, g.member_count, g.species_count, "
        f"count(*) FILTER (WHERE {category} = 'TARGET') AS target_species_count, "
        f"count(*) FILTER (WHERE {category} = 'OUTSIDE') AS outside_species_count, "
        f"count(*) FILTER (WHERE {category} = 'UNRESOLVED') "
        "AS unresolved_species_count, "
        f"sum(gs.species_member_count) FILTER (WHERE {category} = 'TARGET') "
        "AS target_member_count, "
        f"sum(gs.species_member_count) FILTER (WHERE {category} = 'OUTSIDE') "
        "AS outside_member_count, "
        f"sum(gs.species_member_count) FILTER (WHERE {category} = 'UNRESOLVED') "
        "AS unresolved_member_count, "
        f"string_agg(gs.species_label, '; ' ORDER BY gs.species_label) "
        f"FILTER (WHERE {category} = 'OUTSIDE') AS outsider_species, "
        f"string_agg(gs.species_label, '; ' ORDER BY gs.species_label) "
        f"FILTER (WHERE {category} = 'UNRESOLVED') AS unresolved_species "
        "FROM group_statistics AS g JOIN group_species_statistics AS gs "
        "ON gs.run_id = g.run_id AND gs.group_type = g.group_type "
        "AND gs.hierarchy_node = g.hierarchy_node AND gs.group_id = g.group_id "
        "LEFT JOIN taxonomy_map AS tm ON tm.species_label = gs.species_label "
        "WHERE g.run_id = ? AND g.group_type = ? AND g.hierarchy_node = ? "
        "GROUP BY g.run_id, g.group_type, g.hierarchy_node, g.group_id, "
        "g.legacy_orthogroup_id, g.member_count, g.species_count) "
        "SELECT *, target_species_count::DOUBLE / ? AS target_coverage, "
        "CASE WHEN target_species_count + outside_species_count > 0 THEN "
        "target_species_count::DOUBLE / (target_species_count + outside_species_count) "
        "ELSE 0.0 END AS mapped_target_fraction, "
        "coalesce(outsider_species, '') AS outsider_species_labels, "
        "coalesce(unresolved_species, '') AS unresolved_species_labels "
        "FROM summaries"
    )
    return sql, tuple(parameters)


def _taxonomy_filter_condition(*, filters: TaxonomySearchFilters) -> tuple[str, tuple[object, ...]]:
    """Return explicit contains/exclusive/near-exclusive search semantics."""

    conditions = [
        "target_species_count >= ?",
        "target_coverage >= ?",
        "mapped_target_fraction >= ?",
    ]
    parameters: list[object] = [
        filters.minimum_target_species_count,
        filters.minimum_target_coverage,
        filters.minimum_mapped_purity,
    ]
    if filters.mode == "SAMPLED_EXCLUSIVE":
        conditions.extend(("outside_species_count = 0", "unresolved_species_count = 0"))
    elif filters.mode == "NEAR_EXCLUSIVE":
        conditions.extend(("outside_species_count <= ?", "unresolved_species_count <= ?"))
        parameters.extend(
            (
                filters.maximum_outside_species_count,
                filters.maximum_unresolved_species_count,
            )
        )
    elif filters.mode != "CONTAINS":
        raise InputValidationError(f"SQL taxonomy filtering does not support mode {filters.mode}.")
    return " AND ".join(conditions), tuple(parameters)

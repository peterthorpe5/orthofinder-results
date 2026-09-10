"""Bounded, parameterised DuckDB queries for the interactive application."""

from __future__ import annotations

import logging
import re
from collections.abc import Sequence
from typing import Any

import duckdb

from orthofinder_results.errors import InputValidationError

from .models import (
    DistanceResultFilters,
    DistanceResultSet,
    FocusClusterFilters,
    FocusClusterPage,
    GroupKey,
    GroupSearchFilters,
    ProteinSearchFilters,
    ProteinSearchPage,
    ResourceIdentity,
    SearchPage,
    TaxonomySearchFilters,
    TaxonomySearchPage,
)
from .resource import connect_read_only
from .taxonomy import TaxonomyAuthority, add_fisher_enrichment

_LOGGER = logging.getLogger("orthofinder_interrogation_app.queries")
MAX_GROUP_MEMBER_ROWS = 50_000
MAX_GROUP_DISTANCE_ROWS = 124_750
MAX_TAXONOMY_TEST_GROUPS = 250_000
MAX_ALL_DISTANCE_RESULTS = 250_000
MAX_PROTEIN_SEARCH_ROWS = 1_000
MAX_FOCUS_CLUSTER_ROWS = 250_000
MAX_GROUP_SPECIES_COLLECTION_ROWS = 5_000_000
BENCHMARK_METRICS = (
    "mean_distance",
    "median_distance",
    "population_stddev_distance",
    "distance_interquartile_range",
    "distance_coefficient_of_variation",
)
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
DISTANCE_RESULT_COLUMNS = (
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
    "distance_method",
    "computation_status",
    "member_identifier_resolution",
    "total_member_count",
    "sampled_member_count",
    "sampling_fraction",
    "full_group_matrix",
    "distance_pair_count",
    "unresolved_pair_count",
    "minimum_distance",
    "q05_distance",
    "q25_distance",
    "median_distance",
    "mean_distance",
    "q75_distance",
    "q95_distance",
    "maximum_distance",
    "population_stddev_distance",
    "interquartile_range",
    "relative_distance_spread",
    "mean_comparable_sites",
    "source_file",
    "failure_reason",
)
_DISTANCE_RESULT_EXPRESSIONS = {
    "run_id": "g.run_id",
    "group_type": "g.group_type",
    "hierarchy_node": "g.hierarchy_node",
    "group_id": "g.group_id",
    "legacy_orthogroup_id": "g.legacy_orthogroup_id",
    "gene_tree_parent_clade": "g.gene_tree_parent_clade",
    "member_count": "g.member_count",
    "species_count": "g.species_count",
    "single_copy_species_count": "g.single_copy_species_count",
    "max_copies_per_species": "g.max_copies_per_species",
    "mean_copies_per_species": "g.mean_copies_per_species",
    "distance_method": "d.distance_method",
    "computation_status": "d.computation_status",
    "member_identifier_resolution": "d.member_identifier_resolution",
    "total_member_count": "d.total_member_count",
    "sampled_member_count": "d.sampled_member_count",
    "sampling_fraction": (
        "CASE WHEN d.total_member_count > 0 THEN "
        "d.sampled_member_count::DOUBLE / d.total_member_count ELSE NULL END"
    ),
    "full_group_matrix": "d.sampled_member_count = d.total_member_count",
    "distance_pair_count": "d.distance_pair_count",
    "unresolved_pair_count": "d.unresolved_pair_count",
    "minimum_distance": "d.minimum_distance",
    "q05_distance": "d.q05_distance",
    "q25_distance": "d.q25_distance",
    "median_distance": "d.median_distance",
    "mean_distance": "d.mean_distance",
    "q75_distance": "d.q75_distance",
    "q95_distance": "d.q95_distance",
    "maximum_distance": "d.maximum_distance",
    "population_stddev_distance": "d.population_stddev_distance",
    "interquartile_range": "d.q75_distance - d.q25_distance",
    "relative_distance_spread": (
        "CASE WHEN d.mean_distance > 0 THEN "
        "d.population_stddev_distance / d.mean_distance ELSE NULL END"
    ),
    "mean_comparable_sites": "d.mean_comparable_sites",
    "source_file": "d.source_file",
    "failure_reason": "d.failure_reason",
}


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

    def has_relation(self, *, relation: str) -> bool:
        """Return whether a validated resource exposes one exact relation.

        Args:
            relation: Unquoted DuckDB relation name.

        Returns:
            ``True`` only for relations observed while opening the resource.
        """

        return relation in self.resource.relations

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
        counts = {key: int(value) for key, value in rows[0].items()}
        counts["portable_tree_count"] = 0
        if self.has_relation(relation="tree_payloads"):
            counts["portable_tree_count"] = int(
                self._query(
                    sql="SELECT count(*) AS tree_count FROM tree_payloads",
                    parameters=(),
                )[0]["tree_count"]
            )
        return counts

    def overview_authorities(self) -> tuple[dict[str, Any], ...]:
        """Return compact run-wide group authority summaries.

        Returns:
            One deterministic record per group type and hierarchy node.
        """

        rows = self._query(
            sql=(
                "SELECT group_type, hierarchy_node, count(*) AS group_count, "
                "min(member_count) AS minimum_members, "
                "median(member_count) AS median_members, "
                "max(member_count) AS maximum_members, "
                "avg(species_count) AS mean_species, "
                "max(species_count) AS maximum_species FROM group_statistics "
                "GROUP BY group_type, hierarchy_node ORDER BY group_type, hierarchy_node"
            ),
            parameters=(),
        )
        return tuple(rows)

    def benchmark_profiles(self) -> tuple[dict[str, Any], ...]:
        """Return available biological dispersion profiles and group counts.

        Returns:
            One row per stored profile, or an empty tuple for older resources.
        """

        if not self.has_relation(relation="benchmark_group_profiles"):
            return ()
        rows = self._query(
            sql=(
                "SELECT profile_id, profile_class, profile_subclass, "
                "membership_role, count(DISTINCT (group_type, hierarchy_node, "
                "group_id)) AS group_count FROM benchmark_group_profiles "
                "GROUP BY profile_id, profile_class, profile_subclass, "
                "membership_role ORDER BY profile_class, profile_subclass, profile_id"
            ),
            parameters=(),
        )
        return tuple(rows)

    def benchmark_cluster_catalogue(self) -> tuple[dict[str, Any], ...]:
        """Return biological target clusters available for individual tests."""

        if not self.has_relation(relation="benchmark_cluster_results"):
            return ()
        rows = self._query(
            sql=(
                "SELECT run_id, group_type, hierarchy_node, group_id, "
                "profile_ids, profile_classes, profile_subclasses, member_count, "
                "species_count, computation_status FROM benchmark_cluster_results "
                "WHERE membership_roles LIKE '%TARGET%' "
                "ORDER BY profile_classes, group_type, hierarchy_node, group_id"
            ),
            parameters=(),
        )
        return tuple(rows)

    def benchmark_distribution(
        self,
        *,
        profile_ids: tuple[str, ...],
        metric: str,
        comparison_scale: str,
    ) -> tuple[dict[str, Any], ...]:
        """Return cluster-level values for selected benchmark profiles.

        Args:
            profile_ids: Exact stored biological profile identifiers.
            metric: Whitelisted cluster dispersion statistic.
            comparison_scale: ``RAW`` or ``MATCHED_RESIDUAL``.

        Returns:
            One cluster-level observation per selected profile membership.

        Raises:
            InputValidationError: If controls are unsupported or unsafe.
        """

        if metric not in BENCHMARK_METRICS:
            raise InputValidationError(f"Unsupported benchmark metric: {metric}")
        if comparison_scale not in {"RAW", "MATCHED_RESIDUAL"}:
            raise InputValidationError(
                f"Unsupported benchmark comparison scale: {comparison_scale}"
            )
        selected = tuple(dict.fromkeys(profile_ids))
        if not selected:
            return ()
        if not self.has_relation(relation="benchmark_group_profiles"):
            return ()
        placeholders = ", ".join("?" for _ in selected)
        if comparison_scale == "RAW":
            sql = (
                "SELECT p.profile_id, p.profile_class, p.profile_subclass, "
                "c.group_type, c.hierarchy_node, c.group_id, "
                f"c.{metric} AS metric_value FROM benchmark_group_profiles AS p "
                "JOIN benchmark_cluster_results AS c USING "
                "(run_id, group_type, hierarchy_node, group_id) "
                f"WHERE p.profile_id IN ({placeholders}) AND c.{metric} IS NOT NULL "
                "ORDER BY p.profile_id, c.group_type, c.hierarchy_node, c.group_id"
            )
        else:
            sql = (
                "WITH residual AS (SELECT run_id, group_type, hierarchy_node, "
                "group_id, observed_value, row_number() OVER (PARTITION BY run_id, "
                "group_type, hierarchy_node, group_id ORDER BY background_profile_id) "
                "AS residual_rank FROM benchmark_individual_comparisons WHERE "
                "comparison_scale = 'MATCHED_RESIDUAL' AND metric = ?) "
                "SELECT p.profile_id, p.profile_class, p.profile_subclass, "
                "p.group_type, p.hierarchy_node, p.group_id, "
                "r.observed_value AS metric_value FROM benchmark_group_profiles AS p "
                "JOIN residual AS r USING (run_id, group_type, hierarchy_node, group_id) "
                f"WHERE r.residual_rank = 1 AND p.profile_id IN ({placeholders}) "
                "ORDER BY p.profile_id, p.group_type, p.hierarchy_node, p.group_id"
            )
            selected = (metric, *selected)
        return tuple(self._query(sql=sql, parameters=selected))

    def benchmark_contrasts(
        self, *, metric: str = ""
    ) -> tuple[dict[str, Any], ...]:
        """Return precomputed profile contrasts, optionally for one metric."""

        if metric and metric not in BENCHMARK_METRICS:
            raise InputValidationError(f"Unsupported benchmark metric: {metric}")
        if not self.has_relation(relation="benchmark_contrasts"):
            return ()
        condition = "WHERE metric = ?" if metric else ""
        parameters: tuple[object, ...] = (metric,) if metric else ()
        rows = self._query(
            sql=(
                "SELECT * FROM benchmark_contrasts "
                f"{condition} ORDER BY metric, target_profile_id, reference_profile_id"
            ),
            parameters=parameters,
        )
        return tuple(rows)

    def benchmark_individual_comparisons(
        self, *, key: GroupKey, metric: str = ""
    ) -> tuple[dict[str, Any], ...]:
        """Return all stored cluster-to-background comparisons for one group."""

        if metric and metric not in BENCHMARK_METRICS:
            raise InputValidationError(f"Unsupported benchmark metric: {metric}")
        if not self.has_relation(relation="benchmark_individual_comparisons"):
            return ()
        condition = (
            "run_id = ? AND group_type = ? AND hierarchy_node = ? AND group_id = ?"
        )
        parameters: tuple[object, ...] = _key_parameters(key=key)
        if metric:
            condition += " AND metric = ?"
            parameters = (*parameters, metric)
        rows = self._query(
            sql=(
                "SELECT * FROM benchmark_individual_comparisons WHERE "
                f"{condition} ORDER BY comparison_scale, metric, background_profile_id"
            ),
            parameters=parameters,
        )
        return tuple(rows)

    def benchmark_classifications(self) -> tuple[dict[str, Any], ...]:
        """Return matched-control central-divergence and spread classes."""

        if not self.has_relation(relation="benchmark_cluster_classifications"):
            return ()
        return tuple(
            self._query(
                sql=(
                    "SELECT * FROM benchmark_cluster_classifications "
                    "ORDER BY group_type, hierarchy_node, group_id"
                ),
                parameters=(),
            )
        )

    def distance_result_facets(self) -> dict[str, tuple[str, ...]]:
        """Return exact stored-distance values available for export filters.

        Returns:
            Group systems, hierarchy nodes, methods and calculation statuses
            represented by successful persisted distance summaries.
        """

        dimensions = {
            "group_types": "group_type",
            "hierarchy_nodes": "hierarchy_node",
            "distance_methods": "distance_method",
            "computation_statuses": "computation_status",
        }
        facets: dict[str, tuple[str, ...]] = {}
        for label, column in dimensions.items():
            rows = self._query(
                sql=(
                    f"SELECT DISTINCT {column} AS value FROM distance_statistics "
                    "WHERE distance_pair_count > 0 ORDER BY value"
                ),
                parameters=(),
            )
            facets[label] = tuple(str(row["value"] or "") for row in rows)
        return facets

    def distance_results(
        self,
        *,
        filters: DistanceResultFilters,
        columns: tuple[str, ...],
        maximum: int = MAX_ALL_DISTANCE_RESULTS,
    ) -> DistanceResultSet:
        """Return every bounded persisted group-distance summary after filtering.

        Exactly one preferred successful distance summary is returned per group.
        If a method or calculation status is selected, preference is evaluated
        only within that exact subset.

        Args:
            filters: Exact result-scope filters.
            columns: Whitelisted fields to materialise, in requested order.
            maximum: Maximum complete rows allowed in the application process.

        Returns:
            Complete result rows, count and selected column order.

        Raises:
            InputValidationError: If columns or the complete result size are unsafe.
        """

        if not columns:
            raise InputValidationError("Select at least one all-results column.")
        if len(set(columns)) != len(columns):
            raise InputValidationError("All-results columns must be unique.")
        unsupported = tuple(column for column in columns if column not in DISTANCE_RESULT_COLUMNS)
        if unsupported:
            raise InputValidationError(
                "Unsupported all-results columns: " + "; ".join(unsupported)
            )
        if not 1 <= maximum <= MAX_ALL_DISTANCE_RESULTS:
            raise InputValidationError(
                f"maximum must be between 1 and {MAX_ALL_DISTANCE_RESULTS:,}."
            )
        count_sql, parameters = _distance_result_query(
            run_id=self.resource.run_id,
            filters=filters,
            selected_sql="count(*) AS result_count",
            ordered=False,
        )
        count = int(
            self._query(
                sql=count_sql,
                parameters=parameters,
            )[0]["result_count"]
        )
        if count > maximum:
            raise InputValidationError(
                f"The selected distance authority contains {count:,} groups; the complete "
                f"interactive export limit is {maximum:,}. Narrow the group system, "
                "species-tree level, method or calculation scope."
            )
        selected = ", ".join(
            f"{_DISTANCE_RESULT_EXPRESSIONS[column]} AS {column}" for column in columns
        )
        result_sql, result_parameters = _distance_result_query(
            run_id=self.resource.run_id,
            filters=filters,
            selected_sql=selected,
        )
        rows = self._query(
            sql=result_sql,
            parameters=result_parameters,
        )
        _LOGGER.info(
            "Complete distance result export loaded: run=%s, rows=%s, columns=%s",
            self.resource.run_id,
            len(rows),
            len(columns),
        )
        return DistanceResultSet(rows=tuple(rows), total_rows=count, columns=columns)

    def search_proteins(self, *, filters: ProteinSearchFilters) -> ProteinSearchPage:
        """Find exact protein memberships and every associated cluster.

        Canonical protein identifiers are always searched. OrthoFinder internal
        identifiers are searched as aliases when the resource contains the
        ``sequences`` relation. Results retain each HOG hierarchy level and legacy
        orthogroup separately because these are distinct biological group records.

        Args:
            filters: Validated identifier, matching and result-bound controls.

        Returns:
            Bounded matching memberships with complete pre-limit result count.

        Raises:
            InputValidationError: If the requested group system is unsupported.
        """

        if filters.maximum_rows > MAX_PROTEIN_SEARCH_ROWS:
            raise InputValidationError(
                f"Protein searches return at most {MAX_PROTEIN_SEARCH_ROWS:,} rows."
            )
        sql, parameters = _protein_search_query(
            run_id=self.resource.run_id,
            filters=filters,
            include_sequence_aliases=self.has_relation(relation="sequences"),
        )
        raw_rows = self._query(sql=sql, parameters=parameters)
        total_rows = int(raw_rows[0]["_complete_match_count"]) if raw_rows else 0
        rows = tuple(
            {
                key: value
                for key, value in row.items()
                if key != "_complete_match_count"
            }
            for row in raw_rows
        )
        _LOGGER.info(
            "Protein search completed: run=%s, mode=%s, group_type=%s, "
            "matches=%s, returned=%s",
            self.resource.run_id,
            filters.match_mode,
            filters.group_type or "ALL",
            total_rows,
            len(rows),
        )
        return ProteinSearchPage(
            rows=rows,
            total_rows=total_rows,
            maximum_rows=filters.maximum_rows,
        )

    def search_focus_clusters(self, *, filters: FocusClusterFilters) -> FocusClusterPage:
        """Map an exact focus-protein authority to aggregated OrthoFinder groups.

        Canonical membership IDs, OrthoFinder internal IDs and controlled UniProt
        accession/entry aliases are matched exactly. The supplied focus list is never
        interpreted as functional proof for unlisted cluster members.

        Args:
            filters: Validated focus identifiers and exact group authority.

        Returns:
            Bounded cluster summaries and complete matching counts.

        Raises:
            InputValidationError: If the group authority is unsupported.
        """

        if filters.maximum_rows > MAX_FOCUS_CLUSTER_ROWS:
            raise InputValidationError(
                f"Focus searches return at most {MAX_FOCUS_CLUSTER_ROWS:,} clusters."
            )
        sql, parameters = _focus_cluster_query(
            run_id=self.resource.run_id,
            filters=filters,
            include_sequence_aliases=self.has_relation(relation="sequences"),
        )
        raw_rows = self._query(sql=sql, parameters=parameters)
        total_rows = int(raw_rows[0]["_complete_cluster_count"]) if raw_rows else 0
        matched_focus = int(raw_rows[0]["_matched_focus_total"]) if raw_rows else 0
        rows = tuple(
            {
                name: value
                for name, value in row.items()
                if name not in {"_complete_cluster_count", "_matched_focus_total"}
            }
            for row in raw_rows
        )
        _LOGGER.info(
            "Focus cluster search completed: run=%s, group_type=%s, node=%s, "
            "submitted=%s, matched_focus=%s, clusters=%s, returned=%s",
            self.resource.run_id,
            filters.group_type,
            filters.hierarchy_node or "ROOT",
            len(filters.protein_identifiers),
            matched_focus,
            total_rows,
            len(rows),
        )
        return FocusClusterPage(
            rows=rows,
            total_rows=total_rows,
            matched_focus_identifiers=matched_focus,
            submitted_focus_identifiers=len(filters.protein_identifiers),
        )

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

    def get_group_species_collection(
        self,
        *,
        group_type: str,
        hierarchy_node: str,
        group_ids: Sequence[str] = (),
        maximum_rows: int = MAX_GROUP_SPECIES_COLLECTION_ROWS,
    ) -> tuple[dict[str, Any], ...]:
        """Return a bounded exact group-by-species collection for selection audits.

        Args:
            group_type: Exact HOG or legacy group authority.
            hierarchy_node: Exact HOG level or empty legacy root.
            group_ids: Optional exact focus-cluster identifiers.
            maximum_rows: Maximum group/species rows materialised from DuckDB.

        Returns:
            Complete rows within the declared bound, in deterministic key order.

        Raises:
            InputValidationError: If authority, identifiers or size are unsafe.
        """

        if group_type not in _MEMBERSHIP_RELATIONS:
            raise InputValidationError(f"Unsupported group collection type: {group_type}")
        if group_type == "LEGACY_ORTHOGROUP" and hierarchy_node:
            raise InputValidationError(
                "Legacy orthogroup collections require the ROOT hierarchy."
            )
        if not isinstance(maximum_rows, int) or isinstance(maximum_rows, bool):
            raise InputValidationError("maximum_rows must be an integer.")
        if not 1 <= maximum_rows <= MAX_GROUP_SPECIES_COLLECTION_ROWS:
            raise InputValidationError(
                "maximum_rows must be between 1 and "
                f"{MAX_GROUP_SPECIES_COLLECTION_ROWS:,}."
            )
        if isinstance(group_ids, (str, bytes)):
            raise InputValidationError("group_ids must be a sequence of exact identifiers.")
        identifiers = tuple(value.strip() for value in group_ids)
        if any(not value or len(value) > 512 or "\x00" in value for value in identifiers):
            raise InputValidationError("A group collection identifier is empty or unsafe.")
        if len(identifiers) != len(set(identifiers)):
            raise InputValidationError("Group collection identifiers must be unique.")
        if len(identifiers) > 100_000:
            raise InputValidationError("A group collection accepts at most 100,000 IDs.")
        conditions = ["run_id = ?", "group_type = ?", "hierarchy_node = ?"]
        parameters: list[object] = [self.resource.run_id, group_type, hierarchy_node]
        if identifiers:
            conditions.append("group_id IN (SELECT unnest(?::VARCHAR[]))")
            parameters.append(list(identifiers))
        where_sql = " AND ".join(conditions)
        count = int(
            self._query(
                sql=(
                    "SELECT count(*) AS row_count FROM group_species_statistics WHERE "
                    + where_sql
                ),
                parameters=parameters,
            )[0]["row_count"]
        )
        if count > maximum_rows:
            raise InputValidationError(
                f"Selected group collection contains {count:,} group/species rows; "
                f"limit is {maximum_rows:,}. Restrict the hierarchy or focus clusters."
            )
        rows = self._query(
            sql=(
                "SELECT run_id, group_type, hierarchy_node, group_id, species_label, "
                "species_member_count FROM group_species_statistics WHERE "
                f"{where_sql} ORDER BY run_id, group_type, hierarchy_node, group_id, "
                "species_label"
            ),
            parameters=parameters,
        )
        _LOGGER.info(
            "Group/species collection loaded: run=%s, type=%s, node=%s, groups=%s, "
            "rows=%s",
            self.resource.run_id,
            group_type,
            hierarchy_node or "ROOT",
            len(identifiers) if identifiers else "ALL",
            len(rows),
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
                "SELECT run_id, group_type, hierarchy_node, group_id, species_label, "
                "member_id, legacy_orthogroup_id, gene_tree_parent_clade FROM "
                f"{relation} WHERE run_id = ? AND group_type = ? "
                "AND hierarchy_node = ? AND group_id = ? "
                "ORDER BY species_label, member_id LIMIT ?"
            ),
            parameters=(*_key_parameters(key=key), maximum),
        )
        return tuple(rows)

    def get_group_distances(
        self,
        *,
        key: GroupKey,
        distance_method: str = "",
        maximum: int = MAX_GROUP_DISTANCE_ROWS,
    ) -> tuple[dict[str, Any], ...]:
        """Return bounded persisted member-to-member distance records.

        Args:
            key: Composite run-scoped group identity.
            distance_method: Optional exact calculation method.
            maximum: Maximum rows permitted in the application process.

        Returns:
            Complete persisted rows when their count is within ``maximum``.

        Raises:
            InputValidationError: If the relation is absent or exceeds the safe bound.
        """

        if not self.has_relation(relation="pairwise_distances"):
            return ()
        if not 1 <= maximum <= MAX_GROUP_DISTANCE_ROWS:
            raise InputValidationError(
                f"maximum must be between 1 and {MAX_GROUP_DISTANCE_ROWS:,}."
            )
        condition = (
            "run_id = ? AND group_type = ? AND hierarchy_node = ? AND group_id = ?"
        )
        parameters: tuple[object, ...] = _key_parameters(key=key)
        if distance_method:
            condition += " AND distance_method = ?"
            parameters = (*parameters, distance_method)
        count = int(
            self._query(
                sql=f"SELECT count(*) AS row_count FROM pairwise_distances WHERE {condition}",
                parameters=parameters,
            )[0]["row_count"]
        )
        if count > maximum:
            raise InputValidationError(
                f"Persisted distances for {key.display_label()} contain {count:,} rows; "
                f"the interactive limit is {maximum:,}."
            )
        rows = self._query(
            sql=(
                "SELECT run_id, group_type, hierarchy_node, group_id, member_a, "
                "member_b, distance_method, distance, comparable_sites, mismatch_sites, "
                "computation_status, source_file FROM pairwise_distances WHERE "
                f"{condition} ORDER BY member_a, member_b"
            ),
            parameters=parameters,
        )
        return tuple(rows)

    def get_group_sequence_aliases(self, *, key: GroupKey) -> tuple[dict[str, Any], ...]:
        """Return canonical-to-internal aliases for members of one group.

        Args:
            key: Composite run-scoped group identity.

        Returns:
            Exact member, species and internal identifier triples. An empty tuple is
            valid when ``SequenceIDs.txt`` was unavailable during resource creation.
        """

        relation = _MEMBERSHIP_RELATIONS.get(key.group_type)
        if relation is None:
            raise InputValidationError(f"Unsupported group type: {key.group_type}")
        if not self.has_relation(relation="sequences"):
            return ()
        rows = self._query(
            sql=(
                "SELECT DISTINCT s.member_id, s.species_label, s.internal_id "
                f"FROM sequences AS s JOIN {relation} AS m "
                "ON m.run_id = s.run_id AND m.member_id = s.member_id "
                "AND m.species_label = s.species_label WHERE m.run_id = ? "
                "AND m.group_type = ? AND m.hierarchy_node = ? AND m.group_id = ? "
                "ORDER BY s.member_id, s.species_label, s.internal_id"
            ),
            parameters=_key_parameters(key=key),
        )
        return tuple(rows)

    def get_portable_tree(
        self, *, key: GroupKey, legacy_orthogroup_id: str = ""
    ) -> dict[str, Any] | None:
        """Return the preferred portable gene tree for one group.

        Args:
            key: Composite run-scoped group identity.
            legacy_orthogroup_id: Optional parent orthogroup tree identifier for a HOG.

        Returns:
            Checksum-bound compressed tree record, or ``None`` for schema-2 resources.

        Raises:
            InputValidationError: If duplicate preferred payloads make authority ambiguous.
        """

        if not self.has_relation(relation="tree_payloads"):
            return None
        derived_tree_id = _tree_id_from_hog(group_id=key.group_id)
        candidates = tuple(
            dict.fromkeys(
                value
                for value in (
                    legacy_orthogroup_id.strip(),
                    derived_tree_id,
                    key.group_id,
                )
                if value
            )
        )
        if not candidates:
            return None
        placeholders = ", ".join("?" for _ in candidates)
        rows = self._query(
            sql=(
                "SELECT run_id, tree_type, tree_id, group_id, source_path, "
                "source_size_bytes, source_sha256, payload_encoding, newick_payload "
                "FROM tree_payloads WHERE run_id = ? AND tree_id IN ("
                f"{placeholders}) ORDER BY CASE tree_type "
                "WHEN 'RESOLVED_GENE_TREE' THEN 0 WHEN 'GENE_TREE' THEN 1 ELSE 2 END, "
                "tree_id"
            ),
            parameters=(self.resource.run_id, *candidates),
        )
        if not rows:
            return None
        rows.sort(
            key=lambda row: _tree_priority(
                tree_type=str(row["tree_type"]),
                tree_id=str(row["tree_id"]),
                candidate_ids=candidates,
            )
        )
        first = rows[0]
        first_priority = _tree_priority(
            tree_type=str(first["tree_type"]),
            tree_id=str(first["tree_id"]),
            candidate_ids=candidates,
        )
        if len(rows) > 1:
            second = rows[1]
            second_priority = _tree_priority(
                tree_type=str(second["tree_type"]),
                tree_id=str(second["tree_id"]),
                candidate_ids=candidates,
            )
            if first_priority == second_priority:
                raise InputValidationError(
                    f"Multiple equally preferred portable trees match {key.display_label()}."
                )
        return first

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


def _distance_result_query(
    *,
    run_id: str,
    filters: DistanceResultFilters,
    selected_sql: str,
    ordered: bool = True,
) -> tuple[str, tuple[object, ...]]:
    """Build a one-preferred-summary-per-group distance query."""

    distance_conditions = ["d.distance_pair_count > 0"]
    group_conditions = ["g.run_id = ?"]
    distance_parameters: list[object] = []
    if filters.distance_method:
        distance_conditions.append("d.distance_method = ?")
        distance_parameters.append(filters.distance_method)
    if filters.computation_status:
        distance_conditions.append("d.computation_status = ?")
        distance_parameters.append(filters.computation_status)
    distance_where = " AND ".join(distance_conditions)
    group_parameters: list[object] = [run_id]
    if filters.group_type:
        group_conditions.append("g.group_type = ?")
        group_parameters.append(filters.group_type)
    if filters.hierarchy_node is not None:
        group_conditions.append("g.hierarchy_node = ?")
        group_parameters.append(filters.hierarchy_node)
    group_where = " AND ".join(group_conditions)
    order_sql = " ORDER BY g.group_type, g.hierarchy_node, g.group_id" if ordered else ""
    sql = (
        "WITH ranked_distance AS (SELECT d.*, row_number() OVER (PARTITION BY "
        "d.run_id, d.group_type, d.hierarchy_node, d.group_id ORDER BY "
        "CASE WHEN d.sampled_member_count = d.total_member_count THEN 0 ELSE 1 END, "
        "d.sampled_member_count DESC, d.distance_method, d.source_file) "
        f"AS app_distance_rank FROM distance_statistics AS d WHERE {distance_where}) "
        f"SELECT {selected_sql} FROM group_statistics AS g JOIN ranked_distance AS d "
        "ON d.run_id = g.run_id AND d.group_type = g.group_type "
        "AND d.hierarchy_node = g.hierarchy_node AND d.group_id = g.group_id "
        f"AND d.app_distance_rank = 1 WHERE {group_where}{order_sql}"
    )
    return sql, tuple((*distance_parameters, *group_parameters))


def _protein_search_query(
    *,
    run_id: str,
    filters: ProteinSearchFilters,
    include_sequence_aliases: bool,
) -> tuple[str, tuple[object, ...]]:
    """Build one parameterised protein-to-cluster membership query.

    Args:
        run_id: Immutable resource run identifier.
        filters: Validated protein search controls.
        include_sequence_aliases: Whether the ``sequences`` relation is available.

    Returns:
        DuckDB SQL and ordered bound parameters.

    Raises:
        InputValidationError: If an exact group-system filter is unsupported.
    """

    if filters.group_type:
        relation = _MEMBERSHIP_RELATIONS.get(filters.group_type)
        if relation is None:
            raise InputValidationError(
                f"Unsupported protein group type: {filters.group_type}"
            )
        relations = (relation,)
    else:
        relations = tuple(_MEMBERSHIP_RELATIONS.values())
    exact_mode = filters.match_mode == "EXACT"
    if exact_mode:
        search_value = filters.query
        member_match = "m.member_id = q.value"
        alias_match = "s.internal_id = q.value"
        accession_match = "split_part(m.member_id, '|', 2) = q.value"
        entry_match = "split_part(m.member_id, '|', 3) = q.value"
    else:
        search_value = _literal_contains(value=filters.query)
        member_match = "m.member_id ILIKE q.value ESCAPE '\\'"
        alias_match = "s.internal_id ILIKE q.value ESCAPE '\\'"
        accession_match = "split_part(m.member_id, '|', 2) ILIKE q.value ESCAPE '\\'"
        entry_match = "split_part(m.member_id, '|', 3) ILIKE q.value ESCAPE '\\'"
    alias_pattern = "^(sp|tr)\\|[^|[:space:]]+\\|[^|[:space:]]+$"
    pipe_identifier = f"regexp_full_match(m.member_id, '{alias_pattern}')"
    accession_match = f"({pipe_identifier} AND {accession_match})"
    entry_match = f"({pipe_identifier} AND {entry_match})"

    branches = []
    parameters: list[object] = [search_value]
    for relation in relations:
        if include_sequence_aliases:
            alias_join = (
                "LEFT JOIN sequences AS s ON s.run_id = m.run_id "
                "AND s.member_id = m.member_id AND s.species_label = m.species_label "
            )
            internal_identifier = "coalesce(s.internal_id, '')"
            combined_match = (
                f"({member_match} OR {alias_match} OR {accession_match} OR {entry_match})"
            )
            if exact_mode:
                match_source = (
                    f"CASE WHEN {member_match} AND {alias_match} THEN "
                    "'PROTEIN_AND_INTERNAL_ID' "
                    f"WHEN {member_match} THEN 'PROTEIN_ID' "
                    f"WHEN {alias_match} THEN 'ORTHOFINDER_INTERNAL_ID' "
                    f"WHEN {accession_match} THEN 'UNIPROT_ACCESSION' "
                    "ELSE 'UNIPROT_ENTRY' END"
                )
                matched_identifier = (
                    f"CASE WHEN {member_match} THEN m.member_id "
                    f"WHEN {alias_match} THEN s.internal_id "
                    f"WHEN {accession_match} THEN split_part(m.member_id, '|', 2) "
                    "ELSE split_part(m.member_id, '|', 3) END"
                )
            else:
                match_source = (
                    f"CASE WHEN {alias_match} THEN 'ORTHOFINDER_INTERNAL_ID' "
                    f"WHEN {accession_match} THEN 'UNIPROT_ACCESSION' "
                    f"WHEN {entry_match} THEN 'UNIPROT_ENTRY' "
                    "ELSE 'PROTEIN_ID' END"
                )
                matched_identifier = (
                    f"CASE WHEN {alias_match} THEN s.internal_id "
                    f"WHEN {accession_match} THEN split_part(m.member_id, '|', 2) "
                    f"WHEN {entry_match} THEN split_part(m.member_id, '|', 3) "
                    "ELSE m.member_id END"
                )
        else:
            alias_join = ""
            internal_identifier = "''"
            combined_match = f"({member_match} OR {accession_match} OR {entry_match})"
            if exact_mode:
                match_source = (
                    f"CASE WHEN {member_match} THEN 'PROTEIN_ID' "
                    f"WHEN {accession_match} THEN 'UNIPROT_ACCESSION' "
                    "ELSE 'UNIPROT_ENTRY' END"
                )
                matched_identifier = (
                    f"CASE WHEN {member_match} THEN m.member_id "
                    f"WHEN {accession_match} THEN split_part(m.member_id, '|', 2) "
                    "ELSE split_part(m.member_id, '|', 3) END"
                )
            else:
                match_source = (
                    f"CASE WHEN {accession_match} THEN 'UNIPROT_ACCESSION' "
                    f"WHEN {entry_match} THEN 'UNIPROT_ENTRY' "
                    "ELSE 'PROTEIN_ID' END"
                )
                matched_identifier = (
                    f"CASE WHEN {accession_match} THEN split_part(m.member_id, '|', 2) "
                    f"WHEN {entry_match} THEN split_part(m.member_id, '|', 3) "
                    "ELSE m.member_id END"
                )
        branches.append(
            "SELECT DISTINCT m.run_id, m.group_type, m.hierarchy_node, m.group_id, "
            "m.legacy_orthogroup_id, m.gene_tree_parent_clade, m.species_label, "
            f"m.member_id, {internal_identifier} AS internal_id, "
            f"{match_source} AS match_source, "
            f"{matched_identifier} AS matched_identifier FROM {relation} AS m "
            f"CROSS JOIN search_term AS q {alias_join}"
            f"WHERE m.run_id = ? AND {combined_match}"
        )
        parameters.append(run_id)
    membership_sql = " UNION ALL ".join(branches)
    sql = (
        "WITH search_term AS (SELECT ?::VARCHAR AS value), matched_members AS ("
        f"{membership_sql}), preferred_distance AS (SELECT *, row_number() OVER ("
        "PARTITION BY run_id, group_type, hierarchy_node, group_id ORDER BY "
        "CASE WHEN sampled_member_count = total_member_count THEN 0 ELSE 1 END, "
        "sampled_member_count DESC, distance_method, source_file) "
        "AS app_distance_rank FROM distance_statistics) SELECT mm.run_id, mm.group_type, "
        "mm.hierarchy_node, mm.group_id, mm.legacy_orthogroup_id, "
        "mm.gene_tree_parent_clade, mm.species_label, mm.member_id, mm.internal_id, "
        "mm.match_source, mm.matched_identifier, g.member_count, g.species_count, "
        "g.single_copy_species_count, g.max_copies_per_species, "
        "g.mean_copies_per_species, d.distance_method, d.computation_status, "
        "d.sampled_member_count, d.distance_pair_count, d.mean_distance, "
        "d.population_stddev_distance, count(*) OVER () AS _complete_match_count "
        "FROM matched_members AS mm JOIN group_statistics AS g ON g.run_id = mm.run_id "
        "AND g.group_type = mm.group_type AND g.hierarchy_node = mm.hierarchy_node "
        "AND g.group_id = mm.group_id LEFT JOIN preferred_distance AS d "
        "ON d.run_id = mm.run_id AND d.group_type = mm.group_type "
        "AND d.hierarchy_node = mm.hierarchy_node AND d.group_id = mm.group_id "
        "AND d.app_distance_rank = 1 ORDER BY CASE mm.match_source "
        "WHEN 'PROTEIN_ID' THEN 0 WHEN 'PROTEIN_AND_INTERNAL_ID' THEN 1 "
        "WHEN 'ORTHOFINDER_INTERNAL_ID' THEN 2 WHEN 'UNIPROT_ACCESSION' THEN 3 "
        "ELSE 4 END, "
        "mm.group_type, mm.hierarchy_node, mm.group_id, mm.species_label, mm.member_id "
        "LIMIT ?"
    )
    parameters.append(filters.maximum_rows)
    return sql, tuple(parameters)


def _focus_cluster_query(
    *,
    run_id: str,
    filters: FocusClusterFilters,
    include_sequence_aliases: bool,
) -> tuple[str, tuple[object, ...]]:
    """Build a parameterised focus-authority-to-cluster aggregation query.

    Args:
        run_id: Immutable resource run identifier.
        filters: Exact protein and group-authority controls.
        include_sequence_aliases: Whether SequenceIDs-derived aliases are available.

    Returns:
        DuckDB SQL and ordered bound parameters.

    Raises:
        InputValidationError: If the group system is unsupported.
    """

    relation = _MEMBERSHIP_RELATIONS.get(filters.group_type)
    if relation is None:
        raise InputValidationError(
            f"Unsupported focus group type: {filters.group_type}"
        )
    if filters.group_type == "LEGACY_ORTHOGROUP" and filters.hierarchy_node:
        raise InputValidationError(
            "Legacy orthogroup focus searches require the ROOT hierarchy."
        )
    alias_pattern = "^(sp|tr)\\|[^|[:space:]]+\\|[^|[:space:]]+$"
    parameters: list[object] = [list(filters.protein_identifiers)]
    if include_sequence_aliases:
        alias_sql = (
            "identifier_aliases AS ("
            "SELECT run_id, member_id, species_label, internal_id, member_id AS alias, "
            "'PROTEIN_ID' AS match_authority FROM sequences WHERE run_id = ? UNION ALL "
            "SELECT run_id, member_id, species_label, internal_id, internal_id AS alias, "
            "'ORTHOFINDER_INTERNAL_ID' AS match_authority FROM sequences "
            "WHERE run_id = ? AND internal_id <> '' UNION ALL "
            "SELECT run_id, member_id, species_label, internal_id, "
            "split_part(member_id, '|', 2) AS alias, "
            "'UNIPROT_ACCESSION' AS match_authority FROM sequences "
            f"WHERE run_id = ? AND regexp_full_match(member_id, '{alias_pattern}') UNION ALL "
            "SELECT run_id, member_id, species_label, internal_id, "
            "split_part(member_id, '|', 3) AS alias, "
            "'UNIPROT_ENTRY' AS match_authority FROM sequences "
            f"WHERE run_id = ? AND regexp_full_match(member_id, '{alias_pattern}')), "
            "matched_sequences AS (SELECT DISTINCT aliases.run_id, aliases.member_id, "
            "aliases.species_label, aliases.internal_id, focus.focus_identifier, "
            "aliases.match_authority FROM identifier_aliases AS aliases JOIN focus_terms "
            "AS focus ON focus.focus_identifier = aliases.alias), "
            "matched_memberships AS (SELECT DISTINCT m.run_id, m.group_type, "
            "m.hierarchy_node, m.group_id, m.legacy_orthogroup_id, "
            "m.gene_tree_parent_clade, m.member_id, m.species_label, "
            "matches.focus_identifier, matches.match_authority FROM "
            f"{relation} AS m JOIN matched_sequences AS matches ON "
            "matches.run_id = m.run_id AND matches.member_id = m.member_id "
            "AND matches.species_label = m.species_label WHERE m.run_id = ? "
            "AND m.hierarchy_node = ?), "
        )
        parameters.extend((run_id, run_id, run_id, run_id, run_id, filters.hierarchy_node))
    else:
        alias_sql = (
            "scoped_memberships AS (SELECT * FROM "
            f"{relation} WHERE run_id = ? AND hierarchy_node = ?), "
            "identifier_aliases AS (SELECT *, member_id AS alias, "
            "'PROTEIN_ID' AS match_authority FROM scoped_memberships UNION ALL "
            "SELECT *, split_part(member_id, '|', 2) AS alias, "
            "'UNIPROT_ACCESSION' AS match_authority FROM scoped_memberships "
            f"WHERE regexp_full_match(member_id, '{alias_pattern}') UNION ALL "
            "SELECT *, split_part(member_id, '|', 3) AS alias, "
            "'UNIPROT_ENTRY' AS match_authority FROM scoped_memberships "
            f"WHERE regexp_full_match(member_id, '{alias_pattern}')), "
            "matched_memberships AS (SELECT DISTINCT aliases.run_id, "
            "aliases.group_type, aliases.hierarchy_node, aliases.group_id, "
            "aliases.legacy_orthogroup_id, aliases.gene_tree_parent_clade, "
            "aliases.member_id, aliases.species_label, focus.focus_identifier, "
            "aliases.match_authority FROM identifier_aliases AS aliases JOIN "
            "focus_terms AS focus ON focus.focus_identifier = aliases.alias), "
        )
        parameters.extend((run_id, filters.hierarchy_node))
    sql = (
        "WITH focus_terms AS (SELECT DISTINCT unnest(?::VARCHAR[]) AS focus_identifier), "
        f"{alias_sql} aggregated AS (SELECT run_id, group_type, hierarchy_node, group_id, "
        "min(legacy_orthogroup_id) AS legacy_orthogroup_id, "
        "min(gene_tree_parent_clade) AS gene_tree_parent_clade, "
        "count(DISTINCT focus_identifier) AS matched_focus_count, "
        "count(DISTINCT member_id) AS matched_protein_count, "
        "string_agg(DISTINCT focus_identifier, ';' ORDER BY focus_identifier) "
        "AS matched_focus_identifiers, "
        "string_agg(DISTINCT member_id, ';' ORDER BY member_id) AS matched_member_ids, "
        "string_agg(DISTINCT species_label, ';' ORDER BY species_label) "
        "AS matched_species_labels, "
        "string_agg(DISTINCT match_authority, ';' ORDER BY match_authority) "
        "AS match_authorities FROM matched_memberships GROUP BY run_id, group_type, "
        "hierarchy_node, group_id), preferred_distance AS (SELECT *, row_number() OVER ("
        "PARTITION BY run_id, group_type, hierarchy_node, group_id ORDER BY "
        "CASE WHEN sampled_member_count = total_member_count THEN 0 ELSE 1 END, "
        "sampled_member_count DESC, distance_method, source_file) AS app_distance_rank "
        "FROM distance_statistics), results AS (SELECT matches.*, stats.member_count, "
        "stats.species_count, stats.max_copies_per_species, stats.mean_copies_per_species, "
        "distance.distance_method, distance.computation_status, "
        "distance.sampled_member_count, distance.distance_pair_count, "
        "distance.mean_distance, distance.population_stddev_distance FROM aggregated "
        "AS matches JOIN group_statistics AS stats ON stats.run_id = matches.run_id "
        "AND stats.group_type = matches.group_type "
        "AND stats.hierarchy_node = matches.hierarchy_node "
        "AND stats.group_id = matches.group_id LEFT JOIN preferred_distance AS distance "
        "ON distance.run_id = matches.run_id AND distance.group_type = matches.group_type "
        "AND distance.hierarchy_node = matches.hierarchy_node "
        "AND distance.group_id = matches.group_id AND distance.app_distance_rank = 1), "
        "matched_total AS (SELECT count(DISTINCT focus_identifier) AS value "
        "FROM matched_memberships) SELECT results.*, count(*) OVER () "
        "AS _complete_cluster_count, matched_total.value AS _matched_focus_total "
        "FROM results CROSS JOIN matched_total ORDER BY matched_focus_count DESC, "
        "matched_protein_count DESC, member_count DESC, group_id LIMIT ?"
    )
    parameters.append(filters.maximum_rows)
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


def _tree_priority(
    *, tree_type: str, tree_id: str, candidate_ids: tuple[str, ...]
) -> tuple[int, int]:
    """Return deterministic portable-tree selection priority.

    Args:
        tree_type: Stored tree authority type.
        tree_id: Stored tree identifier.
        candidate_ids: Ordered group-derived candidate identifiers.

    Returns:
        Resolved-tree and preferred-identifier sort ranks.
    """

    type_rank = {"RESOLVED_GENE_TREE": 0, "GENE_TREE": 1}.get(tree_type, 2)
    try:
        identifier_rank = candidate_ids.index(tree_id)
    except ValueError:
        identifier_rank = len(candidate_ids)
    return type_rank, identifier_rank


def _tree_id_from_hog(*, group_id: str) -> str:
    """Return OrthoFinder's corresponding OG tree ID for a canonical HOG ID.

    Args:
        group_id: Exact HOG identifier such as ``N0.HOG0000001``.

    Returns:
        Corresponding ``OG`` identifier, or an empty string for other formats.
    """

    match = re.fullmatch(r"N\d+\.HOG(\d+)", group_id)
    return f"OG{match.group(1)}" if match is not None else ""


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
        "g.legacy_orthogroup_id, g.member_count, g.species_count), "
        "preferred_distance AS (SELECT *, row_number() OVER (PARTITION BY run_id, "
        "group_type, hierarchy_node, group_id ORDER BY CASE WHEN distance_pair_count > 0 "
        "THEN 0 ELSE 1 END, distance_method) AS app_distance_rank "
        "FROM distance_statistics) "
        "SELECT s.*, d.distance_method, d.computation_status, d.sampled_member_count, "
        "d.distance_pair_count, d.mean_distance, d.median_distance, "
        "d.population_stddev_distance, s.target_species_count::DOUBLE / ? "
        "AS target_coverage, "
        "CASE WHEN s.target_species_count + s.outside_species_count > 0 THEN "
        "s.target_species_count::DOUBLE / (s.target_species_count + s.outside_species_count) "
        "ELSE 0.0 END AS mapped_target_fraction, "
        "coalesce(s.outsider_species, '') AS outsider_species_labels, "
        "coalesce(s.unresolved_species, '') AS unresolved_species_labels "
        "FROM summaries AS s LEFT JOIN preferred_distance AS d ON d.run_id = s.run_id "
        "AND d.group_type = s.group_type AND d.hierarchy_node = s.hierarchy_node "
        "AND d.group_id = s.group_id AND d.app_distance_rank = 1"
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

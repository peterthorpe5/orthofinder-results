"""Bounded read-only queries for packaged RNA-seq expression evidence."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Sequence

from orthofinder_results.errors import InputValidationError

from .models import ResourceIdentity
from .resource import connect_read_only

EXPRESSION_RELATIONS = frozenset(
    {
        "expression_member_mapping",
        "expression_member_summary",
        "expression_context",
        "expression_group_summary",
        "expression_import_audit",
    }
)
CONTEXT_COLUMNS = {
    "Expression context": "expression_context",
    "Organism part / tissue": "organism_part",
    "Developmental stage": "developmental_stage",
    "Condition": "condition",
}
_GROUP_TABLES = {
    "HOG": "hog_memberships",
    "LEGACY_ORTHOGROUP": "legacy_orthogroup_memberships",
}


@dataclass(frozen=True)
class BoundedExpressionRows:
    """Bounded expression-query rows and an explicit truncation state.

    Attributes:
        rows: Exact rows retained for display and export.
        truncated: Whether at least one additional matching row existed.
        maximum_rows: User-visible row limit applied to the query.
    """

    rows: tuple[dict[str, Any], ...]
    truncated: bool
    maximum_rows: int


def expression_available(*, resource: ResourceIdentity) -> bool:
    """Return whether the opened resource has the complete expression contract."""

    return EXPRESSION_RELATIONS.issubset(resource.relations)


def expression_dimensions(*, resource: ResourceIdentity) -> dict[str, tuple[str, ...]]:
    """Return exact packaged species, units and biological contexts."""

    _require_expression(resource=resource)
    stat = resource.database_path.stat()
    return _expression_dimensions_cached(
        database_path=str(resource.database_path),
        size_bytes=stat.st_size,
        modified_ns=stat.st_mtime_ns,
    )


@lru_cache(maxsize=8)
def _expression_dimensions_cached(
    *, database_path: str, size_bytes: int, modified_ns: int
) -> dict[str, tuple[str, ...]]:
    """Load low-cardinality expression selectors once per unchanged resource."""

    del size_bytes, modified_ns
    connection = connect_read_only(database_path=Path(database_path))
    try:
        species = tuple(
            str(row[0])
            for row in connection.execute(
                "SELECT DISTINCT species_label FROM expression_context "
                "WHERE species_label <> '' ORDER BY species_label"
            ).fetchall()
        )
        units = tuple(
            str(row[0])
            for row in connection.execute(
                "SELECT DISTINCT expression_unit FROM expression_context "
                "WHERE expression_unit <> '' ORDER BY expression_unit"
            ).fetchall()
        )
    finally:
        connection.close()
    return {"species": species, "units": units}


def expression_group_summaries(
    *,
    resource: ResourceIdentity,
    group_type: str,
    hierarchy_node: str,
    group_ids: Sequence[str],
) -> tuple[dict[str, Any], ...]:
    """Return packaged mapping and expression coverage for selected groups."""

    selected = _validate_group_ids(group_ids=group_ids, maximum=2_000)
    _require_group_type(group_type=group_type)
    _require_expression(resource=resource)
    connection = connect_read_only(database_path=resource.database_path)
    try:
        rows = connection.execute(
            "SELECT * FROM expression_group_summary WHERE run_id = ? AND group_type = ? "
            "AND hierarchy_node = ? AND group_id IN (SELECT unnest(?)) "
            "ORDER BY expression_observed_species_count DESC, mapping_fraction DESC, group_id",
            [resource.run_id, group_type, hierarchy_node, list(selected)],
        ).fetchall()
        columns = tuple(str(column[0]) for column in connection.description)
    finally:
        connection.close()
    return tuple(dict(zip(columns, row, strict=True)) for row in rows)


def expression_group_species(
    *,
    resource: ResourceIdentity,
    group_type: str,
    hierarchy_node: str,
    group_ids: Sequence[str],
    species: Sequence[str] = (),
) -> tuple[dict[str, Any], ...]:
    """Summarise mapping and measured contexts per selected group and species.

    Args:
        resource: Validated immutable OrthoFinder resource.
        group_type: Allow-listed HOG or original-orthogroup authority.
        hierarchy_node: Exact HOG hierarchy node, or empty for original groups.
        group_ids: Bounded group identifiers to summarise.
        species: Optional exact species labels; empty selects every species.

    Returns:
        Ordered group-by-species mapping and expression summaries.
    """

    selected = _validate_group_ids(group_ids=group_ids, maximum=500)
    selected_species = _normalise_species(species=species)
    table = _require_group_type(group_type=group_type)
    _require_expression(resource=resource)
    hierarchy_sql = "AND m.hierarchy_node = ?" if group_type == "HOG" else ""
    species_sql = (
        "AND m.species_label IN (SELECT unnest(?))" if selected_species else ""
    )
    parameters: list[Any] = [resource.run_id]
    if group_type == "HOG":
        parameters.append(hierarchy_node)
    parameters.append(list(selected))
    if selected_species:
        parameters.append(list(selected_species))
    connection = connect_read_only(database_path=resource.database_path)
    try:
        rows = connection.execute(
            f"SELECT m.group_id, m.species_label, count(*) member_count, "
            "count(*) FILTER (WHERE s.mapping_status = 'MAPPED_UNIQUE') mapped_member_count, "
            "count(*) FILTER (WHERE s.context_count > 0) expression_observed_member_count, "
            "count(*) FILTER (WHERE s.broad_expression_supported) broad_expression_member_count, "
            "max(COALESCE(s.context_count, 0)) maximum_member_context_count "
            f"FROM {table} m LEFT JOIN expression_member_summary s USING "
            "(run_id, species_label, member_id) WHERE m.run_id = ? "
            f"{hierarchy_sql} AND m.group_id IN (SELECT unnest(?)) {species_sql} "
            "GROUP BY m.group_id, m.species_label ORDER BY m.group_id, m.species_label",
            parameters,
        ).fetchall()
        columns = tuple(str(column[0]) for column in connection.description)
    finally:
        connection.close()
    return tuple(dict(zip(columns, row, strict=True)) for row in rows)


def expression_selected_member_evidence(
    *,
    resource: ResourceIdentity,
    group_type: str,
    hierarchy_node: str,
    group_ids: Sequence[str],
    species: Sequence[str] = (),
    maximum_rows: int = 100_000,
) -> BoundedExpressionRows:
    """Return mapping and missingness states for selected groups and species.

    Args:
        resource: Validated immutable OrthoFinder resource.
        group_type: Allow-listed HOG or original-orthogroup authority.
        hierarchy_node: Exact HOG hierarchy node, or empty for original groups.
        group_ids: One to 25 selected group identifiers.
        species: Optional exact species labels; empty selects every species.
        maximum_rows: Maximum member rows retained for display and export.

    Returns:
        Bounded member rows with an explicit truncation flag.
    """

    selected = _validate_group_ids(group_ids=group_ids, maximum=25)
    selected_species = _normalise_species(species=species)
    row_limit = _validate_maximum_rows(maximum_rows=maximum_rows)
    table = _require_group_type(group_type=group_type)
    _require_expression(resource=resource)
    hierarchy_sql = "AND m.hierarchy_node = ?" if group_type == "HOG" else ""
    species_sql = (
        "AND m.species_label IN (SELECT unnest(?))" if selected_species else ""
    )
    parameters: list[Any] = [resource.run_id]
    if group_type == "HOG":
        parameters.append(hierarchy_node)
    parameters.append(list(selected))
    if selected_species:
        parameters.append(list(selected_species))
    connection = connect_read_only(database_path=resource.database_path)
    try:
        rows = connection.execute(
            f"SELECT m.group_id, m.species_label, m.member_id, map.mapping_status, "
            "map.mapping_tier, map.matched_gene_count, map.matched_gene_ids, "
            "map.matched_gene_names, map.matched_aliases, map.reason mapping_reason, "
            "s.gene_id, s.gene_name, s.experiment_count, s.selected_expression_units, "
            "s.expression_unit_count, s.context_count, s.positive_context_count, "
            "s.positive_context_fraction, s.minimum_context_expression_value, "
            "s.maximum_context_expression_value, s.median_context_expression_value, "
            "s.broad_expression_supported, s.evidence_status "
            f"FROM {table} m LEFT JOIN expression_member_mapping map USING "
            "(run_id, species_label, member_id) LEFT JOIN expression_member_summary s USING "
            "(run_id, species_label, member_id) WHERE m.run_id = ? "
            f"{hierarchy_sql} AND m.group_id IN (SELECT unnest(?)) {species_sql} "
            "ORDER BY m.group_id, m.species_label, m.member_id "
            f"LIMIT {row_limit + 1}",
            parameters,
        ).fetchall()
        columns = tuple(str(column[0]) for column in connection.description)
    finally:
        connection.close()
    records = tuple(dict(zip(columns, row, strict=True)) for row in rows)
    return BoundedExpressionRows(
        rows=records[:row_limit],
        truncated=len(records) > row_limit,
        maximum_rows=row_limit,
    )


def expression_heatmap_cells(
    *,
    resource: ResourceIdentity,
    group_type: str,
    hierarchy_node: str,
    group_ids: Sequence[str],
    context_column: str,
    expression_unit: str,
    species: Sequence[str] = (),
    maximum_cells: int = 50_000,
) -> tuple[dict[str, Any], ...]:
    """Aggregate selected groups into unit-safe group/species/context cells."""

    selected = _validate_group_ids(group_ids=group_ids, maximum=25)
    table = _require_group_type(group_type=group_type)
    _require_expression(resource=resource)
    if context_column not in CONTEXT_COLUMNS.values():
        raise InputValidationError(f"Unsupported expression context: {context_column}")
    if not expression_unit.strip():
        raise InputValidationError("Select one expression unit; units are never combined.")
    if not 1 <= maximum_cells <= 50_000:
        raise InputValidationError("Expression heatmap cells must be between 1 and 50,000.")
    selected_species = tuple(
        sorted({str(value).strip() for value in species if str(value).strip()})
    )
    hierarchy_sql = "AND m.hierarchy_node = ?" if group_type == "HOG" else ""
    species_sql = (
        "AND e.species_label IN (SELECT unnest(?))" if selected_species else ""
    )
    parameters: list[Any] = [resource.run_id]
    if group_type == "HOG":
        parameters.append(hierarchy_node)
    parameters.extend((list(selected), expression_unit))
    if selected_species:
        parameters.append(list(selected_species))
    connection = connect_read_only(database_path=resource.database_path)
    try:
        rows = connection.execute(
            f"SELECT m.group_id, e.species_label, "
            f"COALESCE(NULLIF(trim(CAST(e.{context_column} AS VARCHAR)), ''), 'Unknown') "
            "context_label, e.expression_unit, median(e.expression_value) median_expression, "
            "min(e.expression_value) minimum_expression, max(e.expression_value) "
            "maximum_expression, count(*) context_row_count, "
            "count(DISTINCT e.member_id) mapped_member_count, "
            "count(DISTINCT e.experiment_accession) experiment_count, "
            "avg(CASE WHEN e.expression_positive THEN 1.0 ELSE 0.0 END) "
            "positive_context_fraction "
            f"FROM {table} m JOIN expression_context e USING "
            "(run_id, species_label, member_id) WHERE m.run_id = ? "
            f"{hierarchy_sql} AND m.group_id IN (SELECT unnest(?)) "
            f"AND e.expression_unit = ? {species_sql} "
            "GROUP BY m.group_id, e.species_label, context_label, e.expression_unit "
            "ORDER BY m.group_id, e.species_label, context_label "
            f"LIMIT {int(maximum_cells)}",
            parameters,
        ).fetchall()
        columns = tuple(str(column[0]) for column in connection.description)
    finally:
        connection.close()
    return tuple(dict(zip(columns, row, strict=True)) for row in rows)


def expression_context_records(
    *,
    resource: ResourceIdentity,
    group_type: str,
    hierarchy_node: str,
    group_ids: Sequence[str],
    expression_unit: str,
    species: Sequence[str] = (),
    maximum_rows: int = 50_000,
) -> BoundedExpressionRows:
    """Return bounded underlying RNA-seq contexts for selected groups and species.

    Args:
        resource: Validated immutable OrthoFinder resource.
        group_type: Allow-listed HOG or original-orthogroup authority.
        hierarchy_node: Exact HOG hierarchy node, or empty for original groups.
        group_ids: One to 25 selected group identifiers.
        expression_unit: Exact TPM or FPKM-style unit retained in the resource.
        species: Optional exact species labels; empty selects every species.
        maximum_rows: Maximum context rows retained for display and export.

    Returns:
        Bounded context-level evidence with an explicit truncation flag.

    Raises:
        InputValidationError: If a unit, species, group or row bound is unsafe.
    """

    selected = _validate_group_ids(group_ids=group_ids, maximum=25)
    selected_species = _normalise_species(species=species)
    row_limit = _validate_maximum_rows(maximum_rows=maximum_rows)
    table = _require_group_type(group_type=group_type)
    _require_expression(resource=resource)
    if not expression_unit.strip():
        raise InputValidationError("Select one expression unit; units are never combined.")
    hierarchy_sql = "AND m.hierarchy_node = ?" if group_type == "HOG" else ""
    species_sql = (
        "AND e.species_label IN (SELECT unnest(?))" if selected_species else ""
    )
    parameters: list[Any] = [resource.run_id]
    if group_type == "HOG":
        parameters.append(hierarchy_node)
    parameters.extend((list(selected), expression_unit))
    if selected_species:
        parameters.append(list(selected_species))
    connection = connect_read_only(database_path=resource.database_path)
    try:
        rows = connection.execute(
            f"SELECT m.group_id, e.species_label, e.member_id, e.gene_id, e.gene_name, "
            "e.source_database, e.experiment_accession, e.expression_unit, "
            "e.sample_or_condition, e.atlas_group_label, e.assay_ids, e.assay_count, "
            "e.organism_part, e.developmental_stage, e.genotype, e.cultivar, "
            "e.treatment, e.condition, e.expression_context, e.metadata_status, "
            "e.expression_value_statistic, e.expression_summary_type, "
            "e.expression_value, e.expression_minimum, e.expression_lower_quartile, "
            "e.expression_median, e.expression_upper_quartile, e.expression_maximum, "
            "e.expression_positive, e.expression_source_file, "
            "e.expression_source_file_sha256, e.metadata_source_file, "
            "e.metadata_source_file_sha256 "
            f"FROM {table} m JOIN expression_context e USING "
            "(run_id, species_label, member_id) WHERE m.run_id = ? "
            f"{hierarchy_sql} AND m.group_id IN (SELECT unnest(?)) "
            f"AND e.expression_unit = ? {species_sql} "
            "ORDER BY m.group_id, e.species_label, e.member_id, "
            "e.experiment_accession, e.expression_context, e.sample_or_condition "
            f"LIMIT {row_limit + 1}",
            parameters,
        ).fetchall()
        columns = tuple(str(column[0]) for column in connection.description)
    finally:
        connection.close()
    records = tuple(dict(zip(columns, row, strict=True)) for row in rows)
    return BoundedExpressionRows(
        rows=records[:row_limit],
        truncated=len(records) > row_limit,
        maximum_rows=row_limit,
    )


def expression_member_evidence(
    *,
    resource: ResourceIdentity,
    group_type: str,
    hierarchy_node: str,
    group_id: str,
) -> tuple[dict[str, Any], ...]:
    """Return exact mapping and aggregate expression states for one group."""

    selected = _validate_group_ids(group_ids=(group_id,), maximum=1)[0]
    table = _require_group_type(group_type=group_type)
    _require_expression(resource=resource)
    hierarchy_sql = "AND m.hierarchy_node = ?" if group_type == "HOG" else ""
    parameters: list[Any] = [resource.run_id]
    if group_type == "HOG":
        parameters.append(hierarchy_node)
    parameters.append(selected)
    connection = connect_read_only(database_path=resource.database_path)
    try:
        rows = connection.execute(
            f"SELECT m.species_label, m.member_id, map.mapping_status, map.mapping_tier, "
            "map.matched_gene_ids, map.matched_gene_names, map.matched_aliases, map.reason, "
            "s.gene_id, s.gene_name, s.experiment_count, s.selected_expression_units, "
            "s.context_count, s.positive_context_count, s.positive_context_fraction, "
            "s.minimum_context_expression_value, s.maximum_context_expression_value, "
            "s.median_context_expression_value, s.broad_expression_supported, "
            f"s.evidence_status FROM {table} m LEFT JOIN expression_member_mapping map USING "
            "(run_id, species_label, member_id) LEFT JOIN expression_member_summary s USING "
            "(run_id, species_label, member_id) WHERE m.run_id = ? "
            f"{hierarchy_sql} AND m.group_id = ? ORDER BY m.species_label, m.member_id",
            parameters,
        ).fetchall()
        columns = tuple(str(column[0]) for column in connection.description)
    finally:
        connection.close()
    return tuple(dict(zip(columns, row, strict=True)) for row in rows)


def _require_expression(*, resource: ResourceIdentity) -> None:
    """Reject resources lacking the complete schema-5 expression contract."""

    missing = sorted(EXPRESSION_RELATIONS.difference(resource.relations))
    if missing:
        raise InputValidationError(
            "This resource has no complete RNA-seq expression evidence: "
            + "; ".join(missing)
        )


def _require_group_type(*, group_type: str) -> str:
    """Return the allow-listed membership table for one group authority."""

    try:
        return _GROUP_TABLES[group_type]
    except KeyError as error:
        raise InputValidationError(
            f"Unsupported expression group type: {group_type}"
        ) from error


def _validate_group_ids(*, group_ids: Sequence[str], maximum: int) -> tuple[str, ...]:
    """Return a bounded, unique, non-empty group selection."""

    selected = tuple(
        dict.fromkeys(str(value).strip() for value in group_ids if str(value).strip())
    )
    if not 1 <= len(selected) <= maximum:
        raise InputValidationError(f"Select between 1 and {maximum:,} groups.")
    if any("\x00" in value for value in selected):
        raise InputValidationError("Group identifiers must not contain NUL characters.")
    return selected


def _normalise_species(*, species: Sequence[str]) -> tuple[str, ...]:
    """Return bounded, unique species labels safe for parameterised queries."""

    selected = tuple(
        sorted({str(value).strip() for value in species if str(value).strip()})
    )
    if len(selected) > 500:
        raise InputValidationError("Select no more than 500 expression species.")
    if any("\x00" in value for value in selected):
        raise InputValidationError("Species labels must not contain NUL characters.")
    return selected


def _validate_maximum_rows(*, maximum_rows: int) -> int:
    """Return a defensive bounded row limit for detailed expression exports."""

    if not isinstance(maximum_rows, int) or isinstance(maximum_rows, bool):
        raise InputValidationError("Expression row limit must be an integer.")
    if not 1 <= maximum_rows <= 100_000:
        raise InputValidationError("Expression row limit must be between 1 and 100,000.")
    return maximum_rows

"""Validated protein-sequence motif queries over complete OrthoFinder groups."""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Sequence

import duckdb

from orthofinder_results.errors import InputValidationError

from .models import ResourceIdentity
from .resource import connect_read_only

DEFAULT_MOTIF = "N"
DEFAULT_THRESHOLD = 0.80
EXACT_C_TERMINAL = "EXACT_C_TERMINAL"
REGEX_ANYWHERE = "REGEX_ANYWHERE"
REGEX_C_TERMINAL = "REGEX_C_TERMINAL"
SEARCH_MODES = frozenset({EXACT_C_TERMINAL, REGEX_ANYWHERE, REGEX_C_TERMINAL})
GROUP_TYPES = frozenset({"HOG", "LEGACY_ORTHOGROUP"})
MAX_MOTIF_LENGTH = 100
MAX_REGEX_LENGTH = 500
MAX_RESULT_ROWS = 20_000
_CANONICAL_MOTIF = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
_CONTROL_CHARACTERS = re.compile(r"[\x00-\x1f\x7f]")
REQUIRED_COLUMNS = frozenset({"internal_id", "species_label", "member_id", "sequence"})


@dataclass(frozen=True)
class SequenceSearch:
    """One validated exact or regular-expression protein search."""

    mode: str
    expression: str
    query_expression: str

    @property
    def display_label(self) -> str:
        """Return a concise description for plots and exports."""

        labels = {
            EXACT_C_TERMINAL: "exact C-terminal motif",
            REGEX_ANYWHERE: "regular expression anywhere",
            REGEX_C_TERMINAL: "C-terminal regular expression",
        }
        return f"{labels[self.mode]} {self.expression}"


def validate_motif(*, motif: str) -> str:
    """Return a normalised exact amino-acid suffix.

    Args:
        motif: User-supplied one-letter amino-acid sequence.

    Returns:
        Upper-case canonical motif.

    Raises:
        InputValidationError: If the motif is empty, too long or non-canonical.
    """

    if not isinstance(motif, str):
        raise InputValidationError("The C-terminal motif must be text.")
    normalised = "".join(motif.split()).upper()
    if not normalised:
        raise InputValidationError("Enter at least one terminal amino-acid residue.")
    if len(normalised) > MAX_MOTIF_LENGTH:
        raise InputValidationError(f"Terminal motifs are limited to {MAX_MOTIF_LENGTH} residues.")
    if _CANONICAL_MOTIF.fullmatch(normalised) is None:
        raise InputValidationError(
            "Use canonical one-letter amino-acid codes only: ACDEFGHIKLMNPQRSTVWY."
        )
    return normalised


def validate_sequence_search(*, expression: str, mode: str) -> SequenceSearch:
    """Validate an exact suffix or an explicitly enabled regular expression.

    Args:
        expression: Exact amino-acid suffix or regular expression.
        mode: One of the supported sequence-search modes.

    Returns:
        Immutable validated search definition.

    Raises:
        InputValidationError: If the mode or expression is invalid.
    """

    if mode not in SEARCH_MODES:
        raise InputValidationError(f"Unsupported sequence-search mode: {mode}")
    if mode == EXACT_C_TERMINAL:
        motif = validate_motif(motif=expression)
        return SequenceSearch(mode=mode, expression=motif, query_expression=motif)
    if not isinstance(expression, str):
        raise InputValidationError("The regular expression must be text.")
    pattern = expression.strip()
    if not pattern:
        raise InputValidationError("Enter a regular expression after enabling regex search.")
    if len(pattern) > MAX_REGEX_LENGTH:
        raise InputValidationError(
            f"Regular expressions are limited to {MAX_REGEX_LENGTH} characters."
        )
    if _CONTROL_CHARACTERS.search(pattern) is not None:
        raise InputValidationError("Regular expressions must not contain control characters.")
    terminal_anchor_present = pattern.endswith("$") or pattern.endswith(r"\Z")
    query_expression = (
        pattern
        if mode == REGEX_ANYWHERE or terminal_anchor_present
        else f"(?:{pattern})$"
    )
    try:
        re.compile(query_expression)
    except re.error as error:
        raise InputValidationError(f"Invalid regular expression: {error}") from error
    return SequenceSearch(
        mode=mode,
        expression=pattern,
        query_expression=query_expression,
    )


def validate_sequence_sidecar(*, path: Path) -> Path:
    """Validate the sequence Parquet schema without loading its sequences."""

    source = Path(path).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise InputValidationError(f"Sequence sidecar is missing or empty: {source}")
    details = source.stat()
    _validate_sequence_sidecar_cached(
        source_text=str(source),
        size_bytes=details.st_size,
        modified_ns=details.st_mtime_ns,
    )
    return source


@lru_cache(maxsize=8)
def _validate_sequence_sidecar_cached(
    *, source_text: str, size_bytes: int, modified_ns: int
) -> tuple[str, ...]:
    """Inspect one unchanged Parquet schema once per application process."""

    del size_bytes, modified_ns
    connection = duckdb.connect()
    try:
        description = connection.execute(
            "DESCRIBE SELECT * FROM read_parquet(?)", [source_text]
        ).fetchall()
    except duckdb.Error as error:
        raise InputValidationError(
            f"Sequence sidecar is not readable Parquet: {source_text}: {error}"
        ) from error
    finally:
        connection.close()
    columns = {str(row[0]) for row in description}
    missing = sorted(REQUIRED_COLUMNS - columns)
    if missing:
        raise InputValidationError("Sequence sidecar lacks required columns: " + "; ".join(missing))
    return tuple(sorted(columns))


def motif_species(*, sidecar_path: Path) -> tuple[str, ...]:
    """Return exact species labels represented by the sequence sidecar."""

    source = validate_sequence_sidecar(path=sidecar_path)
    details = source.stat()
    return _motif_species_cached(
        source_text=str(source),
        size_bytes=details.st_size,
        modified_ns=details.st_mtime_ns,
    )


@lru_cache(maxsize=8)
def _motif_species_cached(
    *, source_text: str, size_bytes: int, modified_ns: int
) -> tuple[str, ...]:
    """Return species once per unchanged sidecar."""

    del size_bytes, modified_ns
    connection = duckdb.connect()
    try:
        rows = connection.execute(
            "SELECT DISTINCT species_label FROM read_parquet(?) "
            "WHERE species_label <> '' ORDER BY species_label",
            [source_text],
        ).fetchall()
    finally:
        connection.close()
    return tuple(str(row[0]) for row in rows)


def motif_hierarchy_nodes(*, resource: ResourceIdentity) -> tuple[str, ...]:
    """Return HOG hierarchy nodes present in the completed resource."""

    details = resource.database_path.stat()
    return _motif_hierarchy_nodes_cached(
        database_path=str(resource.database_path),
        size_bytes=details.st_size,
        modified_ns=details.st_mtime_ns,
        run_id=resource.run_id,
    )


@lru_cache(maxsize=8)
def _motif_hierarchy_nodes_cached(
    *, database_path: str, size_bytes: int, modified_ns: int, run_id: str
) -> tuple[str, ...]:
    """Return hierarchy nodes once per unchanged resource."""

    del size_bytes, modified_ns
    connection = connect_read_only(database_path=Path(database_path))
    try:
        rows = connection.execute(
            "SELECT DISTINCT hierarchy_node FROM group_statistics "
            "WHERE run_id = ? AND group_type = 'HOG' AND hierarchy_node <> '' "
            "ORDER BY hierarchy_node",
            [run_id],
        ).fetchall()
    except duckdb.Error as error:
        raise InputValidationError(f"Could not retrieve HOG hierarchy levels: {error}") from error
    finally:
        connection.close()
    return tuple(str(row[0]) for row in rows)


def motif_group_summary(
    *,
    resource: ResourceIdentity,
    sidecar_path: Path,
    motif: str,
    threshold: float = DEFAULT_THRESHOLD,
    minimum_species: int = 2,
    minimum_proteins: int = 2,
    analysis_species: Sequence[str] = (),
    required_species: Sequence[str] = (),
    group_type: str = "HOG",
    hierarchy_node: str = "N0",
    maximum_rows: int = 5_000,
    search_mode: str = EXACT_C_TERMINAL,
) -> tuple[dict[str, Any], ...]:
    """Summarise sequence-pattern conservation for OrthoFinder groups.

    Qualification uses sequence-bearing proteins from ``analysis_species``.
    All sampled members remain in the returned evidence so independent
    comparison lineages can be inspected without entering the denominator.
    """

    search = validate_sequence_search(expression=motif, mode=search_mode)
    source = validate_sequence_sidecar(path=sidecar_path)
    _validate_summary_filters(
        threshold=threshold,
        minimum_species=minimum_species,
        minimum_proteins=minimum_proteins,
        maximum_rows=maximum_rows,
        group_type=group_type,
    )
    selected_species = _normalise_species(
        values=analysis_species,
        fallback=motif_species(sidecar_path=source),
    )
    required = _normalise_species(values=required_species, fallback=(), allow_empty=True)
    table_name = _membership_table(group_type=group_type)
    hierarchy_clause = " AND hierarchy_node = ?" if group_type == "HOG" else ""
    predicate_sql = _match_predicate_sql(search=search)
    connection = connect_read_only(database_path=resource.database_path)
    parameters: list[Any] = [resource.run_id]
    if group_type == "HOG":
        parameters.append(hierarchy_node.strip())
    parameters.extend(
        [
            str(source),
            search.query_expression,
            list(selected_species),
            minimum_species,
            minimum_proteins,
            threshold,
            len(required),
            MAX_RESULT_ROWS if required else maximum_rows,
        ]
    )
    try:
        rows = connection.execute(
            f"WITH membership AS ("
            f" SELECT DISTINCT group_id, species_label, member_id FROM {table_name}"
            f" WHERE run_id = ?{hierarchy_clause}"
            "), joined AS ("
            " SELECT m.group_id, m.species_label, m.member_id,"
            "   upper(rtrim(s.sequence, '*')) AS clean_sequence,"
            "   s.sequence IS NOT NULL AND length(rtrim(s.sequence, '*')) > 0"
            "     AS sequence_available"
            " FROM membership AS m LEFT JOIN read_parquet(?) AS s"
            " ON s.species_label = m.species_label AND s.member_id = m.member_id"
            "), classified AS ("
            " SELECT *, CASE WHEN sequence_available THEN "
            f"   {predicate_sql} ELSE FALSE END AS motif_match,"
            "   species_label IN (SELECT unnest(?)) AS analysis_species"
            " FROM joined"
            "), summarised AS ("
            " SELECT group_id, count(*) AS total_member_count,"
            "   count(*) FILTER (WHERE sequence_available) AS sequence_count,"
            "   count(*) FILTER (WHERE NOT sequence_available) AS unavailable_sequence_count,"
            "   count(DISTINCT species_label) AS species_count,"
            "   count(DISTINCT species_label) FILTER (WHERE sequence_available)"
            "     AS assessed_species_count,"
            "   count(*) FILTER (WHERE motif_match) AS matching_sequence_count,"
            "   count(DISTINCT species_label) FILTER (WHERE motif_match)"
            "     AS matching_species_count,"
            "   count(*) FILTER (WHERE analysis_species) AS analysis_member_count,"
            "   count(*) FILTER (WHERE analysis_species AND sequence_available)"
            "     AS analysis_sequence_count,"
            "   count(*) FILTER (WHERE analysis_species AND NOT sequence_available)"
            "     AS analysis_unavailable_sequence_count,"
            "   count(DISTINCT species_label) FILTER (WHERE analysis_species)"
            "     AS analysis_species_count,"
            "   count(DISTINCT species_label)"
            "     FILTER (WHERE analysis_species AND sequence_available)"
            "     AS analysis_assessed_species_count,"
            "   count(*) FILTER (WHERE analysis_species AND motif_match)"
            "     AS analysis_matching_sequence_count,"
            "   count(DISTINCT species_label)"
            "     FILTER (WHERE analysis_species AND motif_match)"
            "     AS analysis_matching_species_count,"
            "   string_agg(DISTINCT species_label, '; ' ORDER BY species_label)"
            "     AS represented_species,"
            "   string_agg(DISTINCT species_label, '; ' ORDER BY species_label)"
            "     FILTER (WHERE sequence_available) AS assessed_species,"
            "   string_agg(DISTINCT species_label, '; ' ORDER BY species_label)"
            "     FILTER (WHERE motif_match) AS matching_species,"
            "   string_agg(DISTINCT species_label, '; ' ORDER BY species_label)"
            "     FILTER (WHERE analysis_species AND sequence_available)"
            "     AS analysis_assessed_species,"
            "   string_agg(DISTINCT species_label, '; ' ORDER BY species_label)"
            "     FILTER (WHERE analysis_species AND motif_match)"
            "     AS analysis_matching_species"
            " FROM classified GROUP BY group_id"
            ") SELECT *,"
            " sequence_count::DOUBLE / total_member_count AS sequence_coverage,"
            " matching_sequence_count::DOUBLE / nullif(sequence_count, 0)"
            "   AS matching_fraction,"
            " analysis_sequence_count::DOUBLE / nullif(analysis_member_count, 0)"
            "   AS analysis_sequence_coverage,"
            " analysis_matching_sequence_count::DOUBLE / nullif(analysis_sequence_count, 0)"
            "   AS analysis_matching_fraction"
            " FROM summarised WHERE analysis_assessed_species_count >= ?"
            " AND analysis_sequence_count >= ?"
            " AND analysis_matching_sequence_count::DOUBLE"
            "       / nullif(analysis_sequence_count, 0) >= ?"
            " AND matching_species_count >= ?"
            " ORDER BY analysis_matching_species_count DESC,"
            " analysis_matching_fraction DESC, analysis_sequence_count DESC, group_id"
            " LIMIT ?",
            parameters,
        ).fetchall()
        columns = tuple(item[0] for item in connection.description)
    except duckdb.Error as error:
        raise InputValidationError(f"Sequence motif query failed: {error}") from error
    finally:
        connection.close()
    records = tuple(dict(zip(columns, row, strict=True)) for row in rows)
    if required:
        required_set = set(required)
        records = tuple(
            row
            for row in records
            if required_set.issubset(_split_species(row.get("matching_species")))
        )[:maximum_rows]
    return records


def motif_group_members(
    *,
    resource: ResourceIdentity,
    sidecar_path: Path,
    motif: str,
    group_id: str,
    group_type: str = "HOG",
    hierarchy_node: str = "N0",
    search_mode: str = EXACT_C_TERMINAL,
) -> tuple[dict[str, Any], ...]:
    """Return sequence-level pattern calls for one selected group."""

    search = validate_sequence_search(expression=motif, mode=search_mode)
    source = validate_sequence_sidecar(path=sidecar_path)
    details = source.stat()
    sidecar_columns = set(
        _validate_sequence_sidecar_cached(
            source_text=str(source),
            size_bytes=details.st_size,
            modified_ns=details.st_mtime_ns,
        )
    )
    raw_header_sql = (
        "COALESCE(s.raw_header, '')" if "raw_header" in sidecar_columns else "''"
    )
    source_fasta_sql = (
        "COALESCE(s.source_fasta, '')" if "source_fasta" in sidecar_columns else "''"
    )
    table_name = _membership_table(group_type=group_type)
    hierarchy_clause = " AND m.hierarchy_node = ?" if group_type == "HOG" else ""
    predicate_sql = _match_predicate_sql(search=search)
    match_text_sql = _matched_text_sql(search=search)
    parameters: list[Any] = [str(source), resource.run_id]
    if group_type == "HOG":
        parameters.append(hierarchy_node.strip())
    parameters.extend([group_id, search.query_expression])
    parameters.extend(_matched_text_parameters(search=search))
    connection = connect_read_only(database_path=resource.database_path)
    try:
        rows = connection.execute(
            "WITH joined AS ("
            " SELECT DISTINCT m.species_label, m.member_id, s.internal_id,"
            f"   {raw_header_sql} AS raw_header,"
            f"   {source_fasta_sql} AS source_fasta,"
            "   upper(rtrim(s.sequence, '*')) AS clean_sequence,"
            "   s.sequence IS NOT NULL AND length(rtrim(s.sequence, '*')) > 0"
            "     AS sequence_available"
            f" FROM {table_name} AS m LEFT JOIN read_parquet(?) AS s"
            " ON s.species_label = m.species_label AND s.member_id = m.member_id"
            f" WHERE m.run_id = ?{hierarchy_clause} AND m.group_id = ?"
            ") SELECT species_label, member_id, internal_id, raw_header, source_fasta,"
            " CASE WHEN sequence_available THEN length(clean_sequence) END AS sequence_length,"
            " sequence_available, CASE WHEN sequence_available THEN "
            f" {predicate_sql} ELSE FALSE END AS motif_match,"
            " CASE WHEN sequence_available THEN "
            f" {match_text_sql} ELSE '' END AS matched_sequence,"
            " clean_sequence AS sequence"
            " FROM joined ORDER BY motif_match DESC, species_label, member_id",
            parameters,
        ).fetchall()
        columns = tuple(item[0] for item in connection.description)
    except duckdb.Error as error:
        raise InputValidationError(f"Could not retrieve group motif members: {error}") from error
    finally:
        connection.close()
    return tuple(dict(zip(columns, row, strict=True)) for row in rows)


def _validate_summary_filters(
    *,
    threshold: float,
    minimum_species: int,
    minimum_proteins: int,
    maximum_rows: int,
    group_type: str,
) -> None:
    """Reject unsafe or contradictory sequence-summary controls."""

    if not 0.0 <= threshold <= 1.0:
        raise InputValidationError("Matching threshold must be between zero and one.")
    if not 1 <= minimum_species <= 100_000:
        raise InputValidationError("Minimum species must be between 1 and 100,000.")
    if not 1 <= minimum_proteins <= 10_000_000:
        raise InputValidationError("Minimum proteins must be between 1 and 10,000,000.")
    if not 1 <= maximum_rows <= MAX_RESULT_ROWS:
        raise InputValidationError(
            f"Maximum result rows must be between 1 and {MAX_RESULT_ROWS:,}."
        )
    if group_type not in GROUP_TYPES:
        raise InputValidationError(f"Unsupported OrthoFinder group type: {group_type}")


def _normalise_species(
    *, values: Sequence[str], fallback: Sequence[str], allow_empty: bool = False
) -> tuple[str, ...]:
    """Return unique non-empty species labels or a validated fallback."""

    normalised = tuple(sorted({str(value).strip() for value in values if str(value).strip()}))
    if normalised:
        return normalised
    fallback_values = tuple(
        sorted({str(value).strip() for value in fallback if str(value).strip()})
    )
    if fallback_values or allow_empty:
        return fallback_values
    raise InputValidationError("At least one analysis species is required.")


def _membership_table(*, group_type: str) -> str:
    """Return the allow-listed membership relation for a group type."""

    tables = {
        "HOG": "hog_memberships",
        "LEGACY_ORTHOGROUP": "legacy_orthogroup_memberships",
    }
    try:
        return tables[group_type]
    except KeyError as error:
        raise InputValidationError(f"Unsupported OrthoFinder group type: {group_type}") from error


def _match_predicate_sql(*, search: SequenceSearch) -> str:
    """Return a fixed SQL fragment for the validated search mode."""

    if search.mode == EXACT_C_TERMINAL:
        return "ends_with(clean_sequence, ?)"
    return "regexp_matches(clean_sequence, ?)"


def _matched_text_sql(*, search: SequenceSearch) -> str:
    """Return a fixed SQL fragment extracting the first sequence match."""

    if search.mode == EXACT_C_TERMINAL:
        return "CASE WHEN ends_with(clean_sequence, ?) THEN right(clean_sequence, ?) ELSE '' END"
    return "regexp_extract(clean_sequence, ?, 0)"


def _matched_text_parameters(*, search: SequenceSearch) -> list[Any]:
    """Return parameters corresponding to the fixed match-extraction SQL."""

    if search.mode == EXACT_C_TERMINAL:
        return [search.query_expression, len(search.query_expression)]
    return [search.query_expression]


def _split_species(value: Any) -> set[str]:
    """Parse a semicolon-delimited species aggregation."""

    return {item for item in str(value or "").split("; ") if item}

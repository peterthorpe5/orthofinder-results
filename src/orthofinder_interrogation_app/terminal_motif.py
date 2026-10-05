"""Validated C-terminal motif queries over complete OrthoFinder sequences."""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Sequence

import duckdb

from orthofinder_results.errors import InputValidationError

from .models import ResourceIdentity
from .resource import connect_read_only

DEFAULT_MOTIF = "N"
DEFAULT_THRESHOLD = 0.80
MAX_MOTIF_LENGTH = 100
MAX_RESULT_ROWS = 20_000
_CANONICAL_MOTIF = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
REQUIRED_COLUMNS = frozenset({"internal_id", "species_label", "member_id", "sequence"})


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


def validate_sequence_sidecar(*, path: Path) -> Path:
    """Validate the terminal-motif Parquet schema without loading its sequences."""

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
) -> None:
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


def motif_species(*, sidecar_path: Path) -> tuple[str, ...]:
    """Return exact species labels represented by the sequence sidecar."""

    source = validate_sequence_sidecar(path=sidecar_path)
    connection = duckdb.connect()
    try:
        rows = connection.execute(
            "SELECT DISTINCT species_label FROM read_parquet(?) "
            "WHERE species_label <> '' ORDER BY species_label",
            [str(source)],
        ).fetchall()
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
    required_species: Sequence[str] = (),
    hierarchy_node: str = "N0",
    maximum_rows: int = 5_000,
) -> tuple[dict[str, Any], ...]:
    """Summarise exact C-terminal motif conservation for HOGs.

    The denominator is every sequence-bearing HOG member. A required species
    must contain at least one matching protein, not merely occur in the HOG.
    """

    suffix = validate_motif(motif=motif)
    source = validate_sequence_sidecar(path=sidecar_path)
    if not 0.0 <= threshold <= 1.0:
        raise InputValidationError("Matching threshold must be between zero and one.")
    if not 1 <= minimum_species <= 100_000:
        raise InputValidationError("Minimum species must be between 1 and 100,000.")
    if not 1 <= maximum_rows <= MAX_RESULT_ROWS:
        raise InputValidationError(
            f"Maximum result rows must be between 1 and {MAX_RESULT_ROWS:,}."
        )
    required = tuple(
        sorted({str(value).strip() for value in required_species if str(value).strip()})
    )
    connection = connect_read_only(database_path=resource.database_path)
    try:
        rows = connection.execute(
            "WITH joined AS ("
            " SELECT DISTINCT m.group_id, m.species_label, m.member_id, s.sequence,"
            "   ends_with(upper(rtrim(s.sequence, '*')), ?) AS motif_match"
            " FROM hog_memberships AS m"
            " JOIN read_parquet(?) AS s"
            "   ON s.species_label = m.species_label AND s.member_id = m.member_id"
            " WHERE m.run_id = ? AND m.hierarchy_node = ? AND length(s.sequence) > 0"
            "), summarised AS ("
            " SELECT group_id, count(*) AS sequence_count,"
            "   count(DISTINCT species_label) AS species_count,"
            "   count(*) FILTER (WHERE motif_match) AS matching_sequence_count,"
            "   count(DISTINCT species_label) FILTER (WHERE motif_match) AS matching_species_count,"
            "   string_agg(DISTINCT species_label, '; ' ORDER BY species_label)"
            "     AS represented_species,"
            "   string_agg(DISTINCT species_label, '; ' ORDER BY species_label)"
            "     FILTER (WHERE motif_match) AS matching_species"
            " FROM joined GROUP BY group_id"
            ") SELECT *, matching_sequence_count::DOUBLE / sequence_count AS matching_fraction"
            " FROM summarised WHERE species_count >= ?"
            " AND matching_sequence_count::DOUBLE / sequence_count >= ?"
            " AND matching_species_count >= ?"
            " ORDER BY matching_species_count DESC, matching_fraction DESC,"
            " sequence_count DESC, group_id LIMIT ?",
            [
                suffix,
                str(source),
                resource.run_id,
                hierarchy_node.strip(),
                minimum_species,
                threshold,
                len(required),
                MAX_RESULT_ROWS if required else maximum_rows,
            ],
        ).fetchall()
        columns = tuple(item[0] for item in connection.description)
    except duckdb.Error as error:
        raise InputValidationError(f"C-terminal motif query failed: {error}") from error
    finally:
        connection.close()
    records = tuple(dict(zip(columns, row, strict=True)) for row in rows)
    if required:
        required_set = set(required)
        records = tuple(
            row
            for row in records
            if required_set.issubset(set(str(row.get("matching_species") or "").split("; ")))
        )[:maximum_rows]
    return records


def motif_group_members(
    *,
    resource: ResourceIdentity,
    sidecar_path: Path,
    motif: str,
    group_id: str,
    hierarchy_node: str = "N0",
) -> tuple[dict[str, Any], ...]:
    """Return exact sequence-level motif calls for one selected HOG."""

    suffix = validate_motif(motif=motif)
    source = validate_sequence_sidecar(path=sidecar_path)
    connection = connect_read_only(database_path=resource.database_path)
    try:
        rows = connection.execute(
            "SELECT DISTINCT m.species_label, m.member_id, s.internal_id,"
            " length(rtrim(s.sequence, '*')) AS sequence_length,"
            " right(upper(rtrim(s.sequence, '*')), ?) AS observed_terminus,"
            " ends_with(upper(rtrim(s.sequence, '*')), ?) AS motif_match, s.sequence"
            " FROM hog_memberships AS m JOIN read_parquet(?) AS s"
            " ON s.species_label = m.species_label AND s.member_id = m.member_id"
            " WHERE m.run_id = ? AND m.hierarchy_node = ? AND m.group_id = ?"
            " ORDER BY motif_match DESC, m.species_label, m.member_id",
            [len(suffix), suffix, str(source), resource.run_id, hierarchy_node, group_id],
        ).fetchall()
        columns = tuple(item[0] for item in connection.description)
    except duckdb.Error as error:
        raise InputValidationError(f"Could not retrieve HOG motif members: {error}") from error
    finally:
        connection.close()
    return tuple(dict(zip(columns, row, strict=True)) for row in rows)

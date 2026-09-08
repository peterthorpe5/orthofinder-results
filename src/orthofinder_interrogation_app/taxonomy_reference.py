"""Review-required taxonomy candidates from a local NCBI taxdump."""

from __future__ import annotations

import csv
import io
import logging
from collections import defaultdict
from collections.abc import Iterator
from datetime import date
from pathlib import Path

from orthofinder_results.errors import InputValidationError

from .taxonomy import TAXONOMY_COLUMNS

_LOGGER = logging.getLogger("orthofinder_interrogation_app.taxonomy_reference")
NAMES_FILENAME = "names.dmp"
NODES_FILENAME = "nodes.dmp"


def build_taxonomy_candidates(
    *,
    species: tuple[str, ...],
    taxdump_dir: Path,
    source_date: str,
    source_version: str,
) -> bytes:
    """Build exact-name candidate rows for arbitrary resource species labels.

    No candidate is marked reviewed. Unique exact NCBI name matches are emitted as
    ``PENDING_REVIEW``; multiple matches remain ``AMBIGUOUS`` and absent matches remain
    ``UNMAPPED``. This preserves a human decision boundary before descendant searches.

    Args:
        species: Exact species labels stored in one OrthoFinder resource.
        taxdump_dir: Extracted NCBI taxdump directory containing names.dmp and nodes.dmp.
        source_date: ISO date associated with the downloaded taxonomy source.
        source_version: User-declared taxdump release identifier.

    Returns:
        UTF-8 review TSV bytes covering every resource species label.

    Raises:
        InputValidationError: If labels, metadata or taxdump authorities are invalid.
    """

    labels = tuple(sorted(species))
    if not labels or len(labels) != len(set(labels)) or any(not label.strip() for label in labels):
        raise InputValidationError("Resource species labels must be unique and non-empty.")
    try:
        date.fromisoformat(source_date)
    except ValueError as error:
        raise InputValidationError("Taxdump source date must use ISO YYYY-MM-DD format.") from error
    if not source_version.strip():
        raise InputValidationError("Taxdump source version must not be empty.")
    root = Path(taxdump_dir).expanduser().resolve()
    names_path = _required_taxdump_file(root=root, name=NAMES_FILENAME)
    nodes_path = _required_taxdump_file(root=root, name=NODES_FILENAME)
    source_by_label = {label: _source_species_name(label=label) for label in labels}
    targets_by_normalised: dict[str, list[str]] = defaultdict(list)
    for label, source_name in source_by_label.items():
        targets_by_normalised[_normalise_name(value=source_name)].append(label)
    matches: dict[str, dict[int, set[str]]] = {
        label: defaultdict(set) for label in labels
    }
    for taxon_id, name, name_class in _iter_names(path=names_path):
        for label in targets_by_normalised.get(_normalise_name(value=name), ()):
            matches[label][taxon_id].add(name_class)
    candidate_ids = {
        taxon_id for by_taxon in matches.values() for taxon_id in by_taxon
    }
    parent_by_id = _read_parent_map(path=nodes_path)
    lineages = {
        taxon_id: _lineage_ids(taxon_id=taxon_id, parent_by_id=parent_by_id)
        for taxon_id in candidate_ids
    }
    required_names = set(candidate_ids)
    for lineage in lineages.values():
        required_names.update(lineage)
    scientific_names = _read_scientific_names(
        path=names_path,
        required_taxon_ids=required_names,
    )
    missing_names = sorted(required_names.difference(scientific_names))
    if missing_names:
        raise InputValidationError(
            "NCBI names.dmp lacks scientific names for required taxon IDs: "
            + "; ".join(str(value) for value in missing_names[:20])
        )
    rows = [
        _candidate_row(
            label=label,
            source_name=source_by_label[label],
            matches=matches[label],
            lineages=lineages,
            scientific_names=scientific_names,
            source_date=source_date,
            source_version=source_version,
        )
        for label in labels
    ]
    _LOGGER.info(
        "Taxonomy candidates prepared: species=%s, pending=%s, ambiguous=%s, unmapped=%s",
        len(rows),
        sum(row["mapping_status"] == "PENDING_REVIEW" for row in rows),
        sum(row["mapping_status"] == "AMBIGUOUS" for row in rows),
        sum(row["mapping_status"] == "UNMAPPED" for row in rows),
    )
    return _serialise_rows(rows=rows)


def _candidate_row(
    *,
    label: str,
    source_name: str,
    matches: dict[int, set[str]],
    lineages: dict[int, tuple[int, ...]],
    scientific_names: dict[int, str],
    source_date: str,
    source_version: str,
) -> dict[str, object]:
    """Return one pending, ambiguous or unmapped candidate record."""

    row: dict[str, object] = {column: "" for column in TAXONOMY_COLUMNS}
    row.update(
        {
            "workflow_species_label": label,
            "source_species_name": source_name,
            "mapping_source": "NCBI Taxonomy taxdump",
            "source_date": source_date,
            "source_version": source_version,
        }
    )
    if not matches:
        row.update(
            {
                "mapping_status": "UNMAPPED",
                "mapping_method": "NO_EXACT_NCBI_NAME_MATCH",
                "review_note": "No exact taxdump name match; manual review required.",
            }
        )
        return row
    if len(matches) > 1:
        candidate_text = "; ".join(
            f"{taxon_id}:{scientific_names[taxon_id]}"
            for taxon_id in sorted(matches)
        )
        row.update(
            {
                "mapping_status": "AMBIGUOUS",
                "mapping_method": "MULTIPLE_EXACT_NCBI_NAME_MATCHES",
                "review_note": f"Candidate taxon IDs: {candidate_text}",
            }
        )
        return row
    taxon_id = next(iter(matches))
    classes = matches[taxon_id]
    lineage = lineages[taxon_id]
    parent_id = lineage[-1] if lineage else None
    method = (
        "EXACT_SCIENTIFIC_NAME_CANDIDATE"
        if "scientific name" in classes
        else "EXACT_NCBI_NAME_CANDIDATE"
    )
    row.update(
        {
            "accepted_species_name": scientific_names[taxon_id],
            "ncbi_taxon_id": taxon_id,
            "parent_taxon_id": parent_id or "",
            "parent_taxon_name": scientific_names[parent_id] if parent_id else "",
            "lineage_taxon_ids": ";".join(str(value) for value in lineage),
            "lineage_names": ";".join(scientific_names[value] for value in lineage),
            "mapping_status": "PENDING_REVIEW",
            "mapping_method": method,
            "review_note": (
                "Exact-name candidate only. Verify the taxon, then change mapping_status "
                "to REVIEWED and complete reviewer provenance."
            ),
        }
    )
    return row


def _required_taxdump_file(*, root: Path, name: str) -> Path:
    """Return one non-empty extracted NCBI taxdump authority file."""

    path = root / name
    if not path.is_file() or path.stat().st_size == 0:
        raise InputValidationError(f"Taxdump directory lacks non-empty {name}: {root}")
    return path


def _iter_names(*, path: Path) -> Iterator[tuple[int, str, str]]:
    """Yield taxon ID, name text and name class from NCBI names.dmp."""

    try:
        with path.open(mode="r", encoding="utf-8", errors="strict") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = _dmp_fields(line=line)
                if len(fields) < 4:
                    raise InputValidationError(
                        f"Malformed names.dmp line {line_number}: expected four fields."
                    )
                try:
                    taxon_id = int(fields[0])
                except ValueError as error:
                    raise InputValidationError(
                        f"Malformed names.dmp taxon ID at line {line_number}."
                    ) from error
                if taxon_id <= 0 or not fields[1] or not fields[3]:
                    raise InputValidationError(
                        f"Malformed names.dmp record at line {line_number}."
                    )
                yield taxon_id, fields[1], fields[3]
    except (OSError, UnicodeError) as error:
        raise InputValidationError(f"Could not read NCBI names.dmp: {error}") from error


def _read_parent_map(*, path: Path) -> dict[int, int]:
    """Read the complete NCBI node-parent authority with duplicate checks."""

    parents: dict[int, int] = {}
    try:
        with path.open(mode="r", encoding="utf-8", errors="strict") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = _dmp_fields(line=line)
                if len(fields) < 2:
                    raise InputValidationError(
                        f"Malformed nodes.dmp line {line_number}: expected two fields."
                    )
                try:
                    taxon_id, parent_id = int(fields[0]), int(fields[1])
                except ValueError as error:
                    raise InputValidationError(
                        f"Malformed nodes.dmp taxon ID at line {line_number}."
                    ) from error
                if taxon_id <= 0 or parent_id <= 0:
                    raise InputValidationError(
                        f"Non-positive nodes.dmp taxon ID at line {line_number}."
                    )
                if taxon_id in parents:
                    raise InputValidationError(
                        f"Duplicate nodes.dmp taxon ID {taxon_id} at line {line_number}."
                    )
                parents[taxon_id] = parent_id
    except (OSError, UnicodeError) as error:
        raise InputValidationError(f"Could not read NCBI nodes.dmp: {error}") from error
    if not parents:
        raise InputValidationError("NCBI nodes.dmp contains no nodes.")
    return parents


def _lineage_ids(*, taxon_id: int, parent_by_id: dict[int, int]) -> tuple[int, ...]:
    """Return root-to-parent lineage IDs for one candidate taxon."""

    if taxon_id not in parent_by_id:
        raise InputValidationError(f"NCBI nodes.dmp lacks candidate taxon {taxon_id}.")
    reverse: list[int] = []
    seen = {taxon_id}
    current = taxon_id
    while True:
        parent = parent_by_id.get(current)
        if parent is None:
            raise InputValidationError(
                f"NCBI nodes.dmp lineage for taxon {taxon_id} stops at missing node {current}."
            )
        if parent == current:
            break
        if parent in seen:
            raise InputValidationError(f"NCBI nodes.dmp lineage cycle includes taxon {parent}.")
        reverse.append(parent)
        seen.add(parent)
        current = parent
    return tuple(reversed(reverse))


def _read_scientific_names(
    *, path: Path, required_taxon_ids: set[int]
) -> dict[int, str]:
    """Return exact scientific names only for candidate and lineage taxa."""

    names: dict[int, str] = {}
    if not required_taxon_ids:
        return names
    for taxon_id, name, name_class in _iter_names(path=path):
        if taxon_id in required_taxon_ids and name_class == "scientific name":
            if taxon_id in names:
                raise InputValidationError(
                    f"NCBI names.dmp contains multiple scientific names for taxon {taxon_id}."
                )
            names[taxon_id] = name
    return names


def _source_species_name(*, label: str) -> str:
    """Convert common underscore separators without changing other label content."""

    return " ".join(label.replace("_", " ").split())


def _normalise_name(*, value: str) -> str:
    """Return case-insensitive whitespace-normalised exact-match text."""

    return " ".join(value.split()).casefold()


def _dmp_fields(*, line: str) -> tuple[str, ...]:
    """Split one pipe-delimited NCBI taxdump line without interpreting escapes."""

    fields = tuple(field.strip() for field in line.rstrip("\r\n").split("|"))
    return fields[:-1] if fields and not fields[-1] else fields


def _serialise_rows(*, rows: list[dict[str, object]]) -> bytes:
    """Serialise complete candidate rows as deterministic UTF-8 TSV bytes."""

    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=TAXONOMY_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")

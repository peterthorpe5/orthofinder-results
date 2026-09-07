"""Reviewed taxonomy mapping, descendant resolution and enrichment statistics."""

from __future__ import annotations

import csv
import io
import logging
import math
from collections import Counter
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any

from scipy.stats import fisher_exact

from orthofinder_results.errors import InputValidationError

_LOGGER = logging.getLogger("orthofinder_interrogation_app.taxonomy")
MAX_TAXONOMY_BYTES = 5 * 1024 * 1024
MAPPING_STATUSES = frozenset({"REVIEWED", "UNMAPPED", "AMBIGUOUS"})
TAXONOMY_COLUMNS = (
    "workflow_species_label",
    "source_species_name",
    "accepted_species_name",
    "ncbi_taxon_id",
    "parent_taxon_id",
    "parent_taxon_name",
    "lineage_taxon_ids",
    "lineage_names",
    "mapping_status",
    "mapping_method",
    "mapping_source",
    "source_date",
    "source_version",
    "reviewed_by",
    "reviewed_at_utc",
    "review_note",
)


@dataclass(frozen=True)
class TaxonomyRecord:
    """One explicit workflow-label-to-taxonomy mapping decision."""

    workflow_species_label: str
    source_species_name: str
    accepted_species_name: str
    ncbi_taxon_id: int | None
    parent_taxon_id: int | None
    parent_taxon_name: str
    lineage_taxon_ids: tuple[int, ...]
    lineage_names: tuple[str, ...]
    mapping_status: str
    mapping_method: str
    mapping_source: str
    source_date: str
    source_version: str
    reviewed_by: str
    reviewed_at_utc: str
    review_note: str

    def is_descendant_of(self, *, taxon_id: int) -> bool:
        """Return whether this reviewed record is the target or its descendant."""

        return self.mapping_status == "REVIEWED" and (
            self.ncbi_taxon_id == taxon_id or taxon_id in self.lineage_taxon_ids
        )


@dataclass(frozen=True)
class TaxonOption:
    """One taxon identifier available as a descendant-search target."""

    taxon_id: int
    name: str

    def display_label(self) -> str:
        """Return an unambiguous selector label."""

        return f"{self.name} | NCBI taxon {self.taxon_id}"


@dataclass(frozen=True)
class TaxonomyAuthority:
    """Validated mapping decisions scoped to one resource's sampled species."""

    records: tuple[TaxonomyRecord, ...]
    expected_species: tuple[str, ...]

    @property
    def reviewed_records(self) -> tuple[TaxonomyRecord, ...]:
        """Return records explicitly marked reviewed."""

        return tuple(row for row in self.records if row.mapping_status == "REVIEWED")

    @property
    def reviewed_species(self) -> tuple[str, ...]:
        """Return exact workflow labels eligible for taxonomic inference."""

        return tuple(sorted(row.workflow_species_label for row in self.reviewed_records))

    @property
    def unresolved_species(self) -> tuple[str, ...]:
        """Return missing, unmapped and ambiguous workflow labels."""

        reviewed = set(self.reviewed_species)
        return tuple(sorted(set(self.expected_species).difference(reviewed)))

    def summary(self) -> dict[str, int]:
        """Return complete mapping-status counts including missing labels."""

        counts = {"REVIEWED": 0, "UNMAPPED": 0, "AMBIGUOUS": 0, "MISSING": 0}
        by_label = {row.workflow_species_label: row for row in self.records}
        for species in self.expected_species:
            row = by_label.get(species)
            counts[row.mapping_status if row is not None else "MISSING"] += 1
        return counts

    def target_species(self, *, taxon_id: int) -> tuple[str, ...]:
        """Return reviewed sampled species at or below one NCBI taxon ID."""

        if taxon_id <= 0:
            raise InputValidationError("Target NCBI taxon ID must be a positive integer.")
        return tuple(
            sorted(
                row.workflow_species_label
                for row in self.reviewed_records
                if row.is_descendant_of(taxon_id=taxon_id)
            )
        )

    def taxon_options(self) -> tuple[TaxonOption, ...]:
        """Return consistent self and lineage taxa represented by reviewed rows."""

        names_by_id: dict[int, set[str]] = {}
        for row in self.reviewed_records:
            if row.ncbi_taxon_id is not None:
                names_by_id.setdefault(row.ncbi_taxon_id, set()).add(row.accepted_species_name)
            for taxon_id, name in zip(row.lineage_taxon_ids, row.lineage_names, strict=True):
                names_by_id.setdefault(taxon_id, set()).add(name)
        conflicts = {
            taxon_id: sorted(names) for taxon_id, names in names_by_id.items() if len(names) != 1
        }
        if conflicts:
            first_id = sorted(conflicts)[0]
            raise InputValidationError(
                f"Taxonomy mapping uses conflicting names for taxon {first_id}: "
                + "; ".join(conflicts[first_id])
            )
        options = [
            TaxonOption(taxon_id=taxon_id, name=next(iter(names)))
            for taxon_id, names in names_by_id.items()
        ]
        return tuple(sorted(options, key=lambda value: (value.name.casefold(), value.taxon_id)))


def read_taxonomy_mapping(*, path: Path, expected_species: tuple[str, ...]) -> TaxonomyAuthority:
    """Read a bounded taxonomy TSV from a user-selected path."""

    source = Path(path).expanduser().resolve()
    try:
        size = source.stat().st_size
    except OSError as error:
        raise InputValidationError(f"Taxonomy mapping is unavailable: {source}: {error}") from error
    if not 1 <= size <= MAX_TAXONOMY_BYTES:
        raise InputValidationError(
            f"Taxonomy mapping size must be between 1 and {MAX_TAXONOMY_BYTES:,} bytes; "
            f"observed {size:,}."
        )
    try:
        content = source.read_bytes()
    except OSError as error:
        raise InputValidationError(
            f"Taxonomy mapping could not be read: {source}: {error}"
        ) from error
    authority = parse_taxonomy_mapping(data=content, expected_species=expected_species)
    _LOGGER.info(
        "Validated taxonomy mapping: path=%s, reviewed=%s, unresolved=%s",
        source,
        len(authority.reviewed_species),
        len(authority.unresolved_species),
    )
    return authority


def parse_taxonomy_mapping(*, data: bytes, expected_species: tuple[str, ...]) -> TaxonomyAuthority:
    """Parse and validate one versioned reviewed-taxonomy TSV payload."""

    if not data or len(data) > MAX_TAXONOMY_BYTES:
        raise InputValidationError(f"Taxonomy TSV must contain 1–{MAX_TAXONOMY_BYTES:,} bytes.")
    if b"\x00" in data:
        raise InputValidationError("Taxonomy TSV must not contain NUL bytes.")
    try:
        text = data.decode("utf-8-sig")
    except UnicodeError as error:
        raise InputValidationError("Taxonomy TSV must be valid UTF-8.") from error
    reader = csv.DictReader(io.StringIO(text), delimiter="\t", strict=True)
    headings = tuple(reader.fieldnames or ())
    missing = [column for column in TAXONOMY_COLUMNS if column not in headings]
    if missing:
        raise InputValidationError("Taxonomy TSV lacks required columns: " + "; ".join(missing))
    if len(headings) != len(set(headings)):
        raise InputValidationError("Taxonomy TSV contains duplicate column headings.")
    records: list[TaxonomyRecord] = []
    try:
        for line_number, raw_row in enumerate(reader, start=2):
            if None in raw_row:
                raise InputValidationError(
                    f"Taxonomy TSV line {line_number} contains extra tab-separated fields."
                )
            records.append(_parse_record(raw_row=raw_row, line_number=line_number))
    except csv.Error as error:
        raise InputValidationError(f"Taxonomy TSV could not be parsed: {error}") from error
    labels = [row.workflow_species_label for row in records]
    label_counts = Counter(labels)
    duplicates = sorted(label for label, count in label_counts.items() if count > 1)
    if duplicates:
        raise InputValidationError(
            "Taxonomy TSV contains duplicate workflow species labels: " + "; ".join(duplicates)
        )
    expected = tuple(sorted(set(expected_species)))
    if len(expected) != len(expected_species) or any(not value.strip() for value in expected):
        raise InputValidationError("Expected resource species labels must be unique and non-empty.")
    unexpected = sorted(set(labels).difference(expected))
    if unexpected:
        raise InputValidationError(
            "Taxonomy TSV contains species absent from this resource: " + "; ".join(unexpected)
        )
    return TaxonomyAuthority(records=tuple(records), expected_species=expected)


def taxonomy_template(*, species: tuple[str, ...]) -> bytes:
    """Return a review-ready TSV template for exact resource species labels."""

    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=TAXONOMY_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for label in sorted(set(species)):
        row = {column: "" for column in TAXONOMY_COLUMNS}
        row.update(
            {
                "workflow_species_label": label,
                "source_species_name": label.replace("_", " "),
                "mapping_status": "UNMAPPED",
                "review_note": "Review required before descendant filtering.",
            }
        )
        writer.writerow(row)
    return output.getvalue().encode("utf-8")


def taxonomy_audit_rows(*, authority: TaxonomyAuthority) -> tuple[dict[str, Any], ...]:
    """Return complete user-facing mapping rows, including missing decisions."""

    by_label = {row.workflow_species_label: row for row in authority.records}
    rows: list[dict[str, Any]] = []
    for label in authority.expected_species:
        record = by_label.get(label)
        if record is None:
            rows.append(
                {
                    "workflow_species_label": label,
                    "mapping_status": "MISSING",
                    "accepted_species_name": "",
                    "ncbi_taxon_id": None,
                    "parent_taxon_name": "",
                    "parent_taxon_id": None,
                    "mapping_source": "",
                    "source_version": "",
                    "reviewed_by": "",
                    "reviewed_at_utc": "",
                    "review_note": "No mapping row was supplied.",
                }
            )
            continue
        rows.append(
            {
                "workflow_species_label": label,
                "mapping_status": record.mapping_status,
                "accepted_species_name": record.accepted_species_name,
                "ncbi_taxon_id": record.ncbi_taxon_id,
                "parent_taxon_name": record.parent_taxon_name,
                "parent_taxon_id": record.parent_taxon_id,
                "mapping_source": record.mapping_source,
                "source_version": record.source_version,
                "reviewed_by": record.reviewed_by,
                "reviewed_at_utc": record.reviewed_at_utc,
                "review_note": record.review_note,
            }
        )
    return tuple(rows)


def add_fisher_enrichment(
    *,
    rows: tuple[dict[str, Any], ...],
    total_target_species: int,
    total_outside_species: int,
) -> tuple[dict[str, Any], ...]:
    """Add one-sided Fisher exact tests and Benjamini-Hochberg q-values.

    Species presence, rather than copy number, is tested against the reviewed
    sampled species universe. This makes the statistical scope explicit and
    avoids treating unmapped labels as known out-clade absences.
    """

    if total_target_species < 1 or total_outside_species < 1:
        raise InputValidationError(
            "Enrichment requires at least one reviewed target and one reviewed "
            "outside species in the sampled analysis."
        )
    tested: list[dict[str, Any]] = []
    p_values: list[float] = []
    for raw_row in rows:
        row = dict(raw_row)
        target = int(row.get("target_species_count", 0))
        outside = int(row.get("outside_species_count", 0))
        if not 0 <= target <= total_target_species:
            raise InputValidationError("A group target-species count exceeds the universe.")
        if not 0 <= outside <= total_outside_species:
            raise InputValidationError("A group outside-species count exceeds the universe.")
        result = fisher_exact(
            [
                [target, outside],
                [total_target_species - target, total_outside_species - outside],
            ],
            alternative="greater",
        )
        odds_ratio, p_value = float(result.statistic), float(result.pvalue)
        if math.isnan(odds_ratio):
            odds_ratio = 0.0
        row["enrichment_odds_ratio"] = odds_ratio
        row["enrichment_p_value"] = p_value
        tested.append(row)
        p_values.append(p_value)
    q_values = benjamini_hochberg(p_values=tuple(p_values))
    for row, q_value in zip(tested, q_values, strict=True):
        row["enrichment_q_value"] = q_value
    return tuple(tested)


def benjamini_hochberg(*, p_values: tuple[float, ...]) -> tuple[float, ...]:
    """Return stable Benjamini-Hochberg adjusted p-values."""

    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in p_values):
        raise InputValidationError("Enrichment p-values must be finite values from zero to one.")
    count = len(p_values)
    if count == 0:
        return ()
    ordered = sorted(range(count), key=lambda index: (p_values[index], index))
    adjusted = [1.0] * count
    running = 1.0
    for reverse_rank, index in enumerate(reversed(ordered), start=1):
        rank = count - reverse_rank + 1
        running = min(running, p_values[index] * count / rank)
        adjusted[index] = min(1.0, running)
    return tuple(adjusted)


def _parse_record(*, raw_row: dict[str | None, str | None], line_number: int) -> TaxonomyRecord:
    """Parse and validate one mapping decision."""

    row = {str(key): str(value or "").strip() for key, value in raw_row.items()}
    label = row["workflow_species_label"]
    status = row["mapping_status"].upper()
    if not label:
        raise InputValidationError(
            f"Taxonomy TSV line {line_number} has an empty workflow species label."
        )
    if status not in MAPPING_STATUSES:
        raise InputValidationError(
            f"Taxonomy TSV line {line_number} has unsupported mapping_status {status!r}."
        )
    taxon_id = _optional_taxon_id(value=row["ncbi_taxon_id"], line_number=line_number)
    parent_id = _optional_taxon_id(value=row["parent_taxon_id"], line_number=line_number)
    lineage_ids = _taxon_id_list(value=row["lineage_taxon_ids"], line_number=line_number)
    lineage_names = _text_list(value=row["lineage_names"])
    if len(lineage_ids) != len(lineage_names):
        raise InputValidationError(
            f"Taxonomy TSV line {line_number} has unequal lineage ID/name counts."
        )
    if status == "REVIEWED":
        required = (
            "source_species_name",
            "accepted_species_name",
            "mapping_method",
            "mapping_source",
            "source_date",
            "source_version",
            "reviewed_by",
            "reviewed_at_utc",
        )
        absent = [column for column in required if not row[column]]
        if taxon_id is None:
            absent.append("ncbi_taxon_id")
        if absent:
            raise InputValidationError(
                f"Reviewed taxonomy row {line_number} lacks: " + "; ".join(absent)
            )
        _validate_dates(
            source_date=row["source_date"],
            reviewed_at_utc=row["reviewed_at_utc"],
            line_number=line_number,
        )
        if taxon_id in lineage_ids:
            raise InputValidationError(
                f"Taxonomy TSV line {line_number} repeats its own taxon ID in its lineage."
            )
        if parent_id is not None and (not lineage_ids or lineage_ids[-1] != parent_id):
            raise InputValidationError(
                f"Taxonomy TSV line {line_number} parent taxon must be the final lineage taxon."
            )
        if (parent_id is None) != (not row["parent_taxon_name"]):
            raise InputValidationError(
                f"Taxonomy TSV line {line_number} must supply parent ID and name together."
            )
    return TaxonomyRecord(
        workflow_species_label=label,
        source_species_name=row["source_species_name"],
        accepted_species_name=row["accepted_species_name"],
        ncbi_taxon_id=taxon_id,
        parent_taxon_id=parent_id,
        parent_taxon_name=row["parent_taxon_name"],
        lineage_taxon_ids=lineage_ids,
        lineage_names=lineage_names,
        mapping_status=status,
        mapping_method=row["mapping_method"],
        mapping_source=row["mapping_source"],
        source_date=row["source_date"],
        source_version=row["source_version"],
        reviewed_by=row["reviewed_by"],
        reviewed_at_utc=row["reviewed_at_utc"],
        review_note=row["review_note"],
    )


def _optional_taxon_id(*, value: str, line_number: int) -> int | None:
    """Parse an optional positive NCBI taxon identifier."""

    if not value:
        return None
    try:
        taxon_id = int(value)
    except ValueError as error:
        raise InputValidationError(
            f"Taxonomy TSV line {line_number} contains a non-integer taxon ID: {value!r}."
        ) from error
    if taxon_id <= 0:
        raise InputValidationError(
            f"Taxonomy TSV line {line_number} contains a non-positive taxon ID."
        )
    return taxon_id


def _taxon_id_list(*, value: str, line_number: int) -> tuple[int, ...]:
    """Parse a semicolon-delimited lineage ID list."""

    if not value:
        return ()
    parsed = tuple(
        _optional_taxon_id(value=part.strip(), line_number=line_number) for part in value.split(";")
    )
    if any(item is None for item in parsed):
        raise InputValidationError(
            f"Taxonomy TSV line {line_number} contains an empty lineage taxon ID."
        )
    result = tuple(int(item) for item in parsed if item is not None)
    if len(result) != len(set(result)):
        raise InputValidationError(
            f"Taxonomy TSV line {line_number} contains duplicate lineage taxon IDs."
        )
    return result


def _text_list(*, value: str) -> tuple[str, ...]:
    """Parse a semicolon-delimited lineage name list."""

    if not value:
        return ()
    return tuple(part.strip() for part in value.split(";") if part.strip())


def _validate_dates(*, source_date: str, reviewed_at_utc: str, line_number: int) -> None:
    """Validate review provenance date formats."""

    try:
        date.fromisoformat(source_date)
        parsed_review = datetime.fromisoformat(reviewed_at_utc.replace("Z", "+00:00"))
    except ValueError as error:
        raise InputValidationError(
            f"Taxonomy TSV line {line_number} contains an invalid provenance date."
        ) from error
    if parsed_review.tzinfo is None:
        raise InputValidationError(
            f"Taxonomy TSV line {line_number} reviewed_at_utc must include a timezone."
        )

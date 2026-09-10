"""Versioned protein-marker authorities for dispersion benchmarking."""

from __future__ import annotations

import csv
import hashlib
import logging
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from orthofinder_interrogation_app.focus import FocusProtein, FocusProteinAuthority

from .errors import InputValidationError

_LOGGER = logging.getLogger("orthofinder_results.benchmark_authority")
BENCHMARK_FILENAME = "arabidopsis_dispersion_benchmarks.tsv"
BENCHMARK_FIELDS = (
    "marker_id",
    "protein_identifier",
    "protein_entry",
    "marker_name",
    "benchmark_class",
    "benchmark_subclass",
    "domain_architecture",
    "evidence_type",
    "organism",
    "taxon_id",
    "source_title",
    "source_doi",
    "source_table",
    "source_version",
    "enabled",
    "note",
)
MAX_BENCHMARK_BYTES = 8 * 1024 * 1024
MAX_BENCHMARK_RECORDS = 10_000
MAX_FIELD_CHARACTERS = 4_096


@dataclass(frozen=True)
class BenchmarkMarker:
    """One exact protein marker with biological class provenance."""

    marker_id: str
    protein_identifier: str
    protein_entry: str
    marker_name: str
    benchmark_class: str
    benchmark_subclass: str
    domain_architecture: str
    evidence_type: str
    organism: str
    taxon_id: str
    source_title: str
    source_doi: str
    source_table: str
    source_version: str
    enabled: bool
    note: str


@dataclass(frozen=True)
class BenchmarkAuthority:
    """Validated marker records and checksum-bound source identity."""

    records: tuple[BenchmarkMarker, ...]
    source_name: str
    sha256: str

    @property
    def enabled_records(self) -> tuple[BenchmarkMarker, ...]:
        """Return enabled records in deterministic order."""

        return tuple(record for record in self.records if record.enabled)

    def to_focus_authority(self) -> FocusProteinAuthority:
        """Return an exact-identifier authority accepted by group selection."""

        return FocusProteinAuthority(
            records=tuple(
                FocusProtein(
                    identifier=record.protein_identifier,
                    protein_name=record.marker_name,
                    category=(
                        f"{record.benchmark_class}::{record.benchmark_subclass}"
                    ),
                    evidence_type=record.evidence_type,
                    organism=record.organism,
                    source=_source_label(record=record),
                    enabled=record.enabled,
                    note=record.note,
                    taxon_id=record.taxon_id,
                    annotation_scope="benchmark protein marker",
                    catalogue_source=record.source_version,
                )
                for record in self.records
            ),
            source_name=self.source_name,
            sha256=self.sha256,
        )

    def by_identifier(self) -> dict[str, BenchmarkMarker]:
        """Return exact protein identifier to marker metadata."""

        return {record.protein_identifier: record for record in self.records}


def bundled_benchmark_path() -> Path:
    """Return the installed Arabidopsis benchmark authority path.

    Returns:
        Resolved package-data path.

    Raises:
        InputValidationError: If package data are unavailable.
    """

    candidate = files("orthofinder_interrogation_app").joinpath(
        "data", BENCHMARK_FILENAME
    )
    path = Path(str(candidate)).resolve()
    if not path.is_file():
        raise InputValidationError(
            f"Bundled dispersion benchmark authority is unavailable: {path}"
        )
    return path


def read_benchmark_authority(*, path: Path) -> BenchmarkAuthority:
    """Read and validate a plain UTF-8 benchmark TSV.

    Args:
        path: Marker-authority TSV.

    Returns:
        Checksum-bound benchmark authority.

    Raises:
        InputValidationError: If the authority is missing or malformed.
    """

    source = Path(path).expanduser().resolve()
    try:
        size = source.stat().st_size
    except OSError as error:
        raise InputValidationError(
            f"Benchmark authority is unavailable: {source}"
        ) from error
    if not 1 <= size <= MAX_BENCHMARK_BYTES:
        raise InputValidationError(
            f"Benchmark authority must contain 1–{MAX_BENCHMARK_BYTES:,} bytes."
        )
    try:
        raw = source.read_bytes()
        text = raw.decode("utf-8-sig")
    except (OSError, UnicodeError) as error:
        raise InputValidationError(
            f"Benchmark authority must be readable UTF-8: {source}"
        ) from error
    if b"\x00" in raw:
        raise InputValidationError("Benchmark authority contains a NUL byte.")
    reader = csv.DictReader(text.splitlines(), delimiter="\t", strict=True)
    if tuple(reader.fieldnames or ()) != BENCHMARK_FIELDS:
        raise InputValidationError(
            "Benchmark authority headings do not match the required schema."
        )
    records: list[BenchmarkMarker] = []
    try:
        for line_number, raw_row in enumerate(reader, start=2):
            if None in raw_row:
                raise InputValidationError(
                    f"Benchmark authority line {line_number} has extra fields."
                )
            row = {key: str(value or "").strip() for key, value in raw_row.items()}
            _validate_row(row=row, line_number=line_number)
            records.append(
                BenchmarkMarker(
                    marker_id=row["marker_id"],
                    protein_identifier=row["protein_identifier"],
                    protein_entry=row["protein_entry"],
                    marker_name=row["marker_name"],
                    benchmark_class=row["benchmark_class"],
                    benchmark_subclass=row["benchmark_subclass"],
                    domain_architecture=row["domain_architecture"],
                    evidence_type=row["evidence_type"],
                    organism=row["organism"],
                    taxon_id=row["taxon_id"],
                    source_title=row["source_title"],
                    source_doi=row["source_doi"],
                    source_table=row["source_table"],
                    source_version=row["source_version"],
                    enabled=_enabled(value=row["enabled"], line_number=line_number),
                    note=row["note"],
                )
            )
    except csv.Error as error:
        raise InputValidationError(
            f"Benchmark authority could not be parsed: {error}"
        ) from error
    if not records or len(records) > MAX_BENCHMARK_RECORDS:
        raise InputValidationError(
            f"Benchmark authority must contain 1–{MAX_BENCHMARK_RECORDS:,} records."
        )
    identifiers = [record.protein_identifier for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise InputValidationError(
            "Benchmark authority contains duplicate protein identifiers."
        )
    marker_classes: dict[str, tuple[str, str]] = {}
    for record in records:
        classification = (record.benchmark_class, record.benchmark_subclass)
        previous = marker_classes.setdefault(record.marker_id, classification)
        if previous != classification:
            raise InputValidationError(
                f"Benchmark marker {record.marker_id!r} has conflicting classes."
            )
    authority = BenchmarkAuthority(
        records=tuple(records),
        source_name=source.name,
        sha256=hashlib.sha256(raw).hexdigest(),
    )
    _LOGGER.info(
        "Validated dispersion benchmark authority: source=%s, records=%s, "
        "enabled=%s, classes=%s, sha256=%s",
        source,
        len(records),
        len(authority.enabled_records),
        len({record.benchmark_class for record in authority.enabled_records}),
        authority.sha256,
    )
    return authority


def _validate_row(*, row: dict[str, str], line_number: int) -> None:
    """Validate required text and controlled class fields for one row."""

    required = (
        "marker_id",
        "protein_identifier",
        "protein_entry",
        "marker_name",
        "benchmark_class",
        "benchmark_subclass",
        "evidence_type",
        "organism",
        "taxon_id",
        "source_title",
        "source_version",
    )
    missing = [field for field in required if not row[field]]
    if missing:
        raise InputValidationError(
            f"Benchmark authority line {line_number} lacks: {', '.join(missing)}."
        )
    if row["benchmark_class"] not in {"HOUSEKEEPING", "R_NLR"}:
        raise InputValidationError(
            f"Benchmark authority line {line_number} has unsupported class "
            f"{row['benchmark_class']!r}."
        )
    if row["taxon_id"] != "3702":
        raise InputValidationError(
            f"Benchmark authority line {line_number} must describe taxon 3702."
        )
    if any(
        len(value) > MAX_FIELD_CHARACTERS
        or "\x00" in value
        or any(ord(character) < 32 for character in value)
        for value in row.values()
    ):
        raise InputValidationError(
            f"Benchmark authority line {line_number} contains unsafe text."
        )


def _enabled(*, value: str, line_number: int) -> bool:
    """Return a controlled boolean value from one authority row."""

    folded = value.casefold()
    if folded in {"true", "1", "yes"}:
        return True
    if folded in {"false", "0", "no"}:
        return False
    raise InputValidationError(
        f"Benchmark authority line {line_number} has invalid enabled value {value!r}."
    )


def _source_label(*, record: BenchmarkMarker) -> str:
    """Return a compact source label without losing DOI provenance."""

    parts = [record.source_title]
    if record.source_doi:
        parts.append(f"doi:{record.source_doi}")
    if record.source_table:
        parts.append(record.source_table)
    parts.append(record.source_version)
    return "; ".join(parts)

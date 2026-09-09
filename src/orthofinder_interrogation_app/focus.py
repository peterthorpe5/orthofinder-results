"""Versioned, replaceable protein-focus authorities for cluster discovery."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import logging
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from orthofinder_results.errors import InputValidationError

_LOGGER = logging.getLogger("orthofinder_interrogation_app.focus")
MAX_FOCUS_BYTES = 32 * 1024 * 1024
MAX_FOCUS_RECORDS = 100_000
MAX_FOCUS_FIELD_CHARACTERS = 8_192
DEFAULT_FOCUS_FILENAME = "e3_seed_catalogue.tsv"
FOCUS_TEMPLATE_COLUMNS = (
    "protein_identifier",
    "protein_name",
    "category",
    "evidence_type",
    "organism",
    "source",
    "enabled",
    "note",
)
_IDENTIFIER_ALIASES = (
    "protein_identifier",
    "accession",
    "protein_id",
    "gene_id",
)
_NAME_ALIASES = ("protein_name", "gene_name", "name")
_CATEGORY_ALIASES = ("category", "e3_category")
_EVIDENCE_TYPE_ALIASES = ("evidence_type",)
_ORGANISM_ALIASES = ("organism", "species")
_SOURCE_ALIASES = ("source", "mapping_source")
_NOTE_ALIASES = ("note", "review_note")


@dataclass(frozen=True)
class FocusProtein:
    """One exact protein identifier in a user-editable focus authority."""

    identifier: str
    protein_name: str
    category: str
    evidence_type: str
    organism: str
    source: str
    enabled: bool
    note: str
    review_status: str = ""
    ubiquitin_go_status: str = ""
    exclusion_go_term: str = ""
    taxon_id: str = ""
    sequence_md5: str = ""
    sequence_available: str = ""
    sequence_match_count: str = ""
    distinct_sequence_count: str = ""
    sequence_species: str = ""
    sequence_identifiers: str = ""
    protein_sequence_length: str = ""
    annotation_scope: str = ""
    catalogue_source: str = ""


@dataclass(frozen=True)
class FocusProteinAuthority:
    """Validated focus records and their immutable input provenance."""

    records: tuple[FocusProtein, ...]
    source_name: str
    sha256: str

    @property
    def enabled_records(self) -> tuple[FocusProtein, ...]:
        """Return enabled focus proteins in deterministic identifier order."""

        return tuple(row for row in self.records if row.enabled)

    @property
    def identifiers(self) -> tuple[str, ...]:
        """Return exact enabled identifiers for a bounded cluster query."""

        return tuple(row.identifier for row in self.enabled_records)

    def summary(self) -> dict[str, int]:
        """Return complete, enabled and annotated focus-record counts."""

        enabled = self.enabled_records
        return {
            "records": len(self.records),
            "enabled": len(enabled),
            "named": sum(bool(row.protein_name) for row in enabled),
            "categories": len({row.category for row in enabled if row.category}),
            "evidence_types": len(
                {row.evidence_type for row in enabled if row.evidence_type}
            ),
            "organisms": len({row.organism for row in enabled if row.organism}),
        }


def bundled_focus_path() -> Path:
    """Return the installed default E3 focus authority path.

    Returns:
        Path to the packaged, versioned TSV.

    Raises:
        InputValidationError: If package data are unavailable.
    """

    candidate = files("orthofinder_interrogation_app").joinpath(
        "data", DEFAULT_FOCUS_FILENAME
    )
    path = Path(str(candidate)).resolve()
    if not path.is_file():
        raise InputValidationError(f"Bundled focus authority is unavailable: {path}")
    return path


def read_focus_proteins(*, path: Path) -> FocusProteinAuthority:
    """Read a bounded plain or gzip-compressed focus TSV.

    Args:
        path: User-selected focus authority.

    Returns:
        Validated focus authority with checksum provenance.

    Raises:
        InputValidationError: If the file is unavailable, unsafe or malformed.
    """

    source = Path(path).expanduser().resolve()
    try:
        size = source.stat().st_size
    except OSError as error:
        raise InputValidationError(f"Focus protein TSV is unavailable: {source}") from error
    if not 1 <= size <= MAX_FOCUS_BYTES:
        raise InputValidationError(
            f"Focus protein file must contain 1–{MAX_FOCUS_BYTES:,} bytes."
        )
    try:
        raw = source.read_bytes()
        if source.suffix.casefold() == ".gz":
            with gzip.GzipFile(fileobj=io.BytesIO(raw), mode="rb") as handle:
                data = handle.read(MAX_FOCUS_BYTES + 1)
        else:
            data = raw
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise InputValidationError(f"Focus protein file could not be read: {source}") from error
    if len(data) > MAX_FOCUS_BYTES:
        raise InputValidationError("Expanded focus protein TSV exceeds the safety limit.")
    authority = parse_focus_proteins(
        data=data,
        source_name=source.name,
        source_sha256=hashlib.sha256(raw).hexdigest(),
    )
    _LOGGER.info(
        "Validated focus authority: source=%s, records=%s, enabled=%s, sha256=%s",
        source,
        len(authority.records),
        len(authority.enabled_records),
        authority.sha256,
    )
    return authority


def parse_focus_proteins(
    *, data: bytes, source_name: str, source_sha256: str = ""
) -> FocusProteinAuthority:
    """Parse an in-memory focus-protein TSV with controlled header aliases.

    Args:
        data: Decompressed UTF-8 TSV bytes.
        source_name: Human-readable input authority name.
        source_sha256: Optional checksum of the original physical input.

    Returns:
        Validated deterministic focus authority.

    Raises:
        InputValidationError: If text, headings or records are ambiguous.
    """

    if not data or len(data) > MAX_FOCUS_BYTES or b"\x00" in data:
        raise InputValidationError("Focus protein TSV is empty, overlong or contains NUL.")
    try:
        text = data.decode("utf-8-sig")
    except UnicodeError as error:
        raise InputValidationError("Focus protein TSV must be valid UTF-8.") from error
    reader = csv.DictReader(io.StringIO(text), delimiter="\t", strict=True)
    headings = tuple(reader.fieldnames or ())
    if not headings or len(headings) != len(set(headings)):
        raise InputValidationError("Focus protein TSV headings are empty or duplicated.")
    catalogue_schema = "seed_id" in headings
    identifier_column = (
        "seed_id"
        if catalogue_schema
        else _one_alias(
            headings=headings,
            aliases=_IDENTIFIER_ALIASES,
            required=True,
            meaning="protein identifier",
        )
    )
    name_column = _optional_generic_column(
        headings=headings,
        aliases=_NAME_ALIASES,
        meaning="protein name",
        catalogue_schema=catalogue_schema,
    )
    category_column = _optional_generic_column(
        headings=headings,
        aliases=_CATEGORY_ALIASES,
        meaning="category",
        catalogue_schema=catalogue_schema,
    )
    evidence_type_column = _optional_generic_column(
        headings=headings,
        aliases=_EVIDENCE_TYPE_ALIASES,
        meaning="evidence type",
        catalogue_schema=catalogue_schema,
    )
    organism_column = _optional_generic_column(
        headings=headings,
        aliases=_ORGANISM_ALIASES,
        meaning="organism",
        catalogue_schema=catalogue_schema,
    )
    source_column = _optional_generic_column(
        headings=headings,
        aliases=_SOURCE_ALIASES,
        meaning="source",
        catalogue_schema=catalogue_schema,
    )
    note_column = _optional_generic_column(
        headings=headings,
        aliases=_NOTE_ALIASES,
        meaning="note",
        catalogue_schema=catalogue_schema,
    )
    enabled_column = "enabled" if "enabled" in headings else ""
    records: list[FocusProtein] = []
    try:
        for line_number, raw_row in enumerate(reader, start=2):
            if None in raw_row:
                raise InputValidationError(
                    f"Focus protein TSV line {line_number} has extra fields."
                )
            row = {key: str(value or "").strip() for key, value in raw_row.items()}
            identifier = row[identifier_column]
            if catalogue_schema:
                protein_name = _first_populated(
                    row=row,
                    columns=("seed_protein_names", "associated_seed_protein_names"),
                )
                category = _first_populated(
                    row=row,
                    columns=("seed_category", "associated_seed_categories"),
                )
                evidence_type = row.get("seed_evidence_type", "")
                organism = _first_populated(
                    row=row,
                    columns=("seed_organism", "associated_seed_organisms"),
                )
                source = _first_populated(
                    row=row,
                    columns=("seed_source", "catalogue_source"),
                )
                note = row.get("source_value", "")
            else:
                protein_name = row.get(name_column, "")
                category = row.get(category_column, "")
                evidence_type = row.get(evidence_type_column, "")
                organism = row.get(organism_column, "")
                source = row.get(source_column, "")
                note = row.get(note_column, "")
            values = (
                identifier,
                protein_name,
                category,
                evidence_type,
                organism,
                source,
                note,
            )
            if not identifier:
                raise InputValidationError(
                    f"Focus protein TSV line {line_number} has an empty identifier."
                )
            if any(
                len(value) > MAX_FOCUS_FIELD_CHARACTERS
                or "\x00" in value
                or any(ord(character) < 32 for character in value)
                for value in values
            ):
                raise InputValidationError(
                    f"Focus protein TSV line {line_number} contains unsafe text."
                )
            records.append(
                FocusProtein(
                    identifier=identifier,
                    protein_name=protein_name,
                    category=category,
                    evidence_type=evidence_type,
                    organism=organism,
                    source=source,
                    enabled=_parse_enabled(
                        value=row.get(enabled_column, "true"),
                        line_number=line_number,
                    ),
                    note=note,
                    review_status=_first_populated(
                        row=row,
                        columns=(
                            "seed_review_status",
                            "associated_seed_review_statuses",
                        ),
                    ),
                    ubiquitin_go_status=_first_populated(
                        row=row,
                        columns=(
                            "seed_ubiquitin_go_status",
                            "associated_seed_ubiquitin_go_statuses",
                        ),
                    ),
                    exclusion_go_term=row.get("seed_exclusion_go_term", ""),
                    taxon_id=row.get("seed_taxon_id", ""),
                    sequence_md5=row.get("seed_sequence_md5", ""),
                    sequence_available=row.get("sequence_available", ""),
                    sequence_match_count=row.get("sequence_match_count", ""),
                    distinct_sequence_count=row.get("distinct_sequence_count", ""),
                    sequence_species=row.get("sequence_species", ""),
                    sequence_identifiers=row.get("sequence_identifiers", ""),
                    protein_sequence_length=row.get("protein_sequence_length", ""),
                    annotation_scope=row.get("annotation_scope", ""),
                    catalogue_source=row.get("catalogue_source", ""),
                )
            )
    except csv.Error as error:
        raise InputValidationError(f"Focus protein TSV could not be parsed: {error}") from error
    if not records:
        raise InputValidationError("Focus protein TSV contains no data records.")
    if len(records) > MAX_FOCUS_RECORDS:
        raise InputValidationError(
            f"Focus protein TSV exceeds {MAX_FOCUS_RECORDS:,} records."
        )
    identifiers = [row.identifier for row in records]
    if len(identifiers) != len(set(identifiers)):
        raise InputValidationError("Focus protein TSV contains duplicate identifiers.")
    folded: dict[str, list[str]] = {}
    for identifier in identifiers:
        folded.setdefault(identifier.casefold(), []).append(identifier)
    collisions = [values for values in folded.values() if len(set(values)) > 1]
    if collisions:
        raise InputValidationError(
            "Focus protein TSV contains case-colliding identifiers: "
            + "; ".join(" / ".join(sorted(values)) for values in collisions[:20])
        )
    ordered = tuple(sorted(records, key=lambda row: row.identifier))
    if not any(row.enabled for row in ordered):
        raise InputValidationError("Focus protein TSV has no enabled records.")
    digest = source_sha256 or hashlib.sha256(data).hexdigest()
    return FocusProteinAuthority(
        records=ordered,
        source_name=source_name.strip() or "uploaded_focus_proteins.tsv",
        sha256=digest,
    )


def focus_template() -> bytes:
    """Return a small editable TSV demonstrating the custom focus contract."""

    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=FOCUS_TEMPLATE_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerow(
        {
            "protein_identifier": "Q9SA03",
            "protein_name": "FB27_ARATH",
            "category": "F-box E3-associated protein",
            "evidence_type": "reviewed project seed",
            "organism": "Arabidopsis thaliana",
            "source": "Replace with reviewed project authority",
            "enabled": "true",
            "note": "Example only; add, remove or disable rows for your analysis.",
        }
    )
    return output.getvalue().encode("utf-8")


def _one_alias(
    *, headings: tuple[str, ...], aliases: tuple[str, ...], required: bool, meaning: str
) -> str:
    """Resolve zero or one accepted header without silent precedence.

    Args:
        headings: Exact input headings.
        aliases: Accepted alternatives for one semantic field.
        required: Whether one alternative must be present.
        meaning: Human-readable field meaning for errors.

    Returns:
        Selected heading or an empty string for an absent optional field.

    Raises:
        InputValidationError: If the field is missing or ambiguously duplicated.
    """

    present = tuple(alias for alias in aliases if alias in headings)
    if len(present) > 1:
        raise InputValidationError(
            f"Focus protein TSV supplies multiple aliases for {meaning}: "
            + "; ".join(present)
        )
    if required and not present:
        raise InputValidationError(
            f"Focus protein TSV requires one {meaning} heading: " + "; ".join(aliases)
        )
    return present[0] if present else ""


def _optional_generic_column(
    *,
    headings: tuple[str, ...],
    aliases: tuple[str, ...],
    meaning: str,
    catalogue_schema: bool,
) -> str:
    """Resolve an optional generic field outside the rich catalogue schema.

    Args:
        headings: Exact input headings.
        aliases: Accepted alternatives for one semantic field.
        meaning: Human-readable field meaning for errors.
        catalogue_schema: Whether the explicit E3 seed schema was detected.

    Returns:
        Selected generic heading or an empty string.
    """

    if catalogue_schema:
        return ""
    return _one_alias(
        headings=headings,
        aliases=aliases,
        required=False,
        meaning=meaning,
    )


def _first_populated(*, row: dict[str, str], columns: tuple[str, ...]) -> str:
    """Return the first non-empty value from explicitly ordered columns.

    Args:
        row: Normalised TSV row.
        columns: Ordered primary then fallback headings.

    Returns:
        First populated value or an empty string.
    """

    return next((row.get(column, "") for column in columns if row.get(column, "")), "")


def _parse_enabled(*, value: str, line_number: int) -> bool:
    """Parse an optional explicit Boolean enable flag."""

    normalised = value.strip().casefold()
    if normalised in {"true", "yes", "1", "enabled", ""}:
        return True
    if normalised in {"false", "no", "0", "disabled"}:
        return False
    raise InputValidationError(
        f"Focus protein TSV line {line_number} has invalid enabled value {value!r}."
    )

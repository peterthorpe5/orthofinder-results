"""Audited Expression Atlas integration for generic OrthoFinder members.

The integration is deliberately protein-family agnostic.  It derives exact
identifier aliases from the OrthoFinder ``SequenceIDs.txt`` authority, scopes
every match by species, chooses one expression unit per experiment, and keeps
unmapped or ambiguous proteins distinct from measured zero expression.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
import resource
import shutil
import subprocess
import sys
import time
import uuid
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import duckdb

from .errors import InputValidationError, OrthoFinderResultsError, PublicationError
from .io_utils import atomic_write_json, read_tsv, sha256_file, write_tsv

_LOGGER = logging.getLogger("orthofinder_results.expression_analysis")
_WORKER_METADATA_VERSION = 1

EXPRESSION_RESOURCE_TYPES = frozenset(
    {
        "atlas_expression_long",
        "atlas_sample_metadata_long",
        "atlas_sample_metadata_wide",
    }
)
MANIFEST_FIELDS = (
    "resource_id",
    "resource_type",
    "species_column",
    "dataset",
    "path",
    "sha256",
    "include",
)
ALIAS_FIELDS = (
    "run_id",
    "species_label",
    "internal_id",
    "member_id",
    "alias_type",
    "alias_value",
    "mapping_tier",
    "source",
)
MAPPING_FIELDS = (
    "run_id",
    "species_label",
    "member_id",
    "mapping_status",
    "mapping_tier",
    "matched_gene_count",
    "matched_gene_ids",
    "matched_gene_names",
    "matched_aliases",
    "reason",
)
SUMMARY_FIELDS = (
    "run_id",
    "species_label",
    "member_id",
    "mapping_status",
    "gene_id",
    "gene_name",
    "experiment_count",
    "selected_expression_units",
    "expression_unit_count",
    "context_count",
    "positive_context_count",
    "positive_context_fraction",
    "minimum_context_expression_value",
    "maximum_context_expression_value",
    "median_context_expression_value",
    "broad_expression_supported",
    "evidence_status",
)
CONTEXT_FIELDS = (
    "run_id",
    "species_label",
    "member_id",
    "gene_id",
    "gene_name",
    "source_database",
    "experiment_accession",
    "expression_unit",
    "sample_or_condition",
    "atlas_group_label",
    "assay_ids",
    "assay_count",
    "organism_part",
    "developmental_stage",
    "genotype",
    "cultivar",
    "treatment",
    "condition",
    "expression_context",
    "metadata_status",
    "expression_value_statistic",
    "expression_summary_type",
    "expression_value",
    "expression_minimum",
    "expression_lower_quartile",
    "expression_median",
    "expression_upper_quartile",
    "expression_maximum",
    "expression_positive",
    "expression_source_file",
    "expression_source_file_sha256",
    "metadata_source_file",
    "metadata_source_file_sha256",
)
GROUP_SUMMARY_FIELDS = (
    "run_id",
    "group_type",
    "hierarchy_node",
    "group_id",
    "member_count",
    "unique_mapped_member_count",
    "ambiguous_member_count",
    "not_mapped_member_count",
    "expression_observed_member_count",
    "broad_expression_member_count",
    "mapped_species_count",
    "expression_observed_species_count",
    "mapping_fraction",
    "expression_observed_fraction",
    "selected_expression_units",
)
AUDIT_FIELDS = (
    "run_id",
    "expression_manifest",
    "expression_manifest_sha256",
    "expression_partition_count",
    "metadata_partition_count",
    "resource_species_count",
    "orthofinder_species_count",
    "scanned_expression_partition_count",
    "scanned_metadata_partition_count",
    "alias_count",
    "member_count",
    "unique_mapped_member_count",
    "ambiguous_member_count",
    "not_mapped_member_count",
    "expression_context_count",
    "expression_group_count",
    "expression_species",
    "expression_worker_peak_rss_mib",
    "expression_worker_database_bytes",
    "expression_worker_spill_bytes",
    "minimum_expression_value",
    "broad_positive_fraction",
    "unit_selection_policy",
    "mapping_policy",
    "missing_data_policy",
)
ADDITIONAL_ALIAS_FIELDS = (
    "species_label",
    "member_id",
    "alias_type",
    "alias_value",
    "mapping_tier",
    "source",
)
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_PIPE_IDENTIFIER = re.compile(r"^(?:sp|tr)\|([^|]+)\|([^|\s]+)", re.IGNORECASE)
_HEADER_ALIAS = re.compile(r"(?:^|\s)(GN|gene|gene_id|locus|locus_tag)=([^\s;]+)", re.IGNORECASE)


@dataclass(frozen=True)
class ExpressionResource:
    """One checksum-verified Expression Atlas resource partition."""

    resource_id: str
    resource_type: str
    species_column: str
    dataset: str
    path: Path
    sha256: str


@dataclass(frozen=True)
class ExpressionAuthority:
    """Versioned resource manifest and its verified included partitions."""

    manifest_path: Path
    manifest_sha256: str
    records: tuple[ExpressionResource, ...]

    def paths(self, *, resource_type: str, species: Iterable[str]) -> tuple[Path, ...]:
        """Return selected partitions matching exact case-insensitive species labels."""

        selected = {str(value).strip().upper() for value in species if str(value).strip()}
        return tuple(
            sorted(
                {
                    record.path
                    for record in self.records
                    if record.resource_type == resource_type
                    and record.species_column.upper() in selected
                }
            )
        )


@dataclass(frozen=True)
class ExpressionPublication:
    """Counts and metadata produced by one expression integration stage."""

    counts: Mapping[str, int]
    scanned_expression_paths: tuple[Path, ...]
    scanned_metadata_paths: tuple[Path, ...]


def read_expression_manifest(*, path: Path, verify_checksums: bool = True) -> ExpressionAuthority:
    """Read and verify one corrected Expression Atlas resource manifest.

    Args:
        path: Manifest containing checksum-bound Parquet partitions.
        verify_checksums: Recalculate every included partition digest.

    Returns:
        Immutable verified authority.

    Raises:
        InputValidationError: If the manifest or any included resource is invalid.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise InputValidationError(f"Expression resource manifest is missing or empty: {source}")
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise InputValidationError(f"Expression resource manifest has no header: {source}")
        missing = sorted(set(MANIFEST_FIELDS).difference(reader.fieldnames))
        if missing:
            raise InputValidationError(
                "Expression resource manifest lacks fields: " + "; ".join(missing)
            )
        rows = tuple(reader)
    records: list[ExpressionResource] = []
    identifiers: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        include = str(row.get("include", "")).strip().casefold()
        if include not in {"true", "false", "1", "0", "yes", "no"}:
            raise InputValidationError(
                f"Invalid include value at expression manifest line {line_number}: {include!r}"
            )
        if include in {"false", "0", "no"}:
            continue
        resource_id = str(row.get("resource_id", "")).strip()
        resource_type = str(row.get("resource_type", "")).strip()
        species_column = str(row.get("species_column", "")).strip()
        digest = str(row.get("sha256", "")).strip().lower()
        if not resource_id or resource_id in identifiers:
            raise InputValidationError(
                f"Empty or duplicate expression resource_id at line {line_number}: {resource_id!r}"
            )
        if resource_type not in EXPRESSION_RESOURCE_TYPES:
            raise InputValidationError(
                f"Unsupported expression resource_type at line {line_number}: {resource_type!r}"
            )
        if not species_column:
            raise InputValidationError(
                f"Expression resource lacks species_column at line {line_number}."
            )
        resource_path = Path(str(row.get("path", ""))).expanduser().resolve()
        if not resource_path.is_file() or resource_path.stat().st_size == 0:
            raise InputValidationError(
                f"Expression resource is missing or empty at line {line_number}: {resource_path}"
            )
        if _SHA256_PATTERN.fullmatch(digest) is None:
            raise InputValidationError(
                f"Expression resource has an invalid SHA-256 at line {line_number}."
            )
        if verify_checksums and sha256_file(path=resource_path) != digest:
            raise InputValidationError(
                f"Expression resource checksum differs at line {line_number}: {resource_path}"
            )
        identifiers.add(resource_id)
        records.append(
            ExpressionResource(
                resource_id=resource_id,
                resource_type=resource_type,
                species_column=species_column,
                dataset=str(row.get("dataset", "")).strip(),
                path=resource_path,
                sha256=digest,
            )
        )
    if not any(record.resource_type == "atlas_expression_long" for record in records):
        raise InputValidationError(
            "Expression resource manifest contains no included atlas_expression_long partition."
        )
    return ExpressionAuthority(
        manifest_path=source,
        manifest_sha256=sha256_file(path=source),
        records=tuple(records),
    )


def publish_expression_evidence(
    *,
    tables_dir: Path,
    work_dir: Path,
    run_id: str,
    authority: ExpressionAuthority,
    additional_aliases_path: Path | None = None,
    minimum_expression_value: float = 0.5,
    broad_positive_fraction: float = 0.5,
    threads: int = 4,
    memory_limit_mb: int = 8_192,
) -> ExpressionPublication:
    """Map every eligible OrthoFinder member to versioned RNA-seq evidence.

    Args:
        tables_dir: Staging directory containing sequence and membership TSVs.
        work_dir: Scratch directory available for DuckDB spill files.
        run_id: Immutable OrthoFinder resource identifier.
        authority: Verified Expression Atlas manifest.
        additional_aliases_path: Optional reviewed exact-alias TSV.
        minimum_expression_value: Context value treated as positive evidence.
        broad_positive_fraction: Positive-context fraction called broad.
        threads: Maximum analytical DuckDB threads.
        memory_limit_mb: DuckDB memory bound before scratch spilling.

    Returns:
        Published relation counts and selected input paths.

    Raises:
        InputValidationError: If controls or source schemas are invalid.
        PublicationError: If expression integration fails.
    """

    _validate_controls(
        minimum_expression_value=minimum_expression_value,
        broad_positive_fraction=broad_positive_fraction,
        threads=threads,
        memory_limit_mb=memory_limit_mb,
    )
    sequence_path = tables_dir / "sequences.tsv.gz"
    species = tuple(
        sorted(
            {
                str(row.get("species_label", "")).strip()
                for row in read_tsv(path=sequence_path)
                if str(row.get("species_label", "")).strip()
            }
        )
    )
    expression_paths = authority.paths(resource_type="atlas_expression_long", species=species)
    metadata_paths = authority.paths(resource_type="atlas_sample_metadata_wide", species=species)
    alias_path = tables_dir / "expression_identifier_aliases.tsv.gz"
    alias_count = write_tsv(
        path=alias_path,
        fieldnames=ALIAS_FIELDS,
        records=_iter_aliases(
            sequence_path=sequence_path,
            run_id=run_id,
            additional_aliases_path=additional_aliases_path,
            allowed_species={
                record.species_column
                for record in authority.records
                if record.path in set(expression_paths)
            },
        ),
    )
    if alias_count == 0:
        raise InputValidationError(
            "No exact member identifiers were available for RNA-seq expression mapping."
        )
    worker_root = Path(work_dir).expanduser().resolve() / f"expression_worker.{uuid.uuid4().hex}"
    worker_root.mkdir(parents=True, exist_ok=False)
    request_path = worker_root / "request.json"
    metadata_path = worker_root / "metadata.json"
    atomic_write_json(
        path=request_path,
        record={
            "metadata_version": _WORKER_METADATA_VERSION,
            "tables_dir": str(Path(tables_dir).expanduser().resolve()),
            "work_dir": str(worker_root),
            "run_id": run_id,
            "expression_paths": [str(path) for path in expression_paths],
            "all_expression_paths": [
                str(record.path)
                for record in authority.records
                if record.resource_type == "atlas_expression_long"
            ],
            "metadata_paths": [str(path) for path in metadata_paths],
            "minimum_expression_value": minimum_expression_value,
            "broad_positive_fraction": broad_positive_fraction,
            "threads": threads,
            "memory_limit_mb": memory_limit_mb,
        },
    )
    worker_metadata = _run_expression_worker(
        request_path=request_path,
        metadata_path=metadata_path,
    )
    counts = dict(worker_metadata["counts"])
    counts["expression_alias_count"] = alias_count
    shutil.rmtree(worker_root)
    audit = {
        "run_id": run_id,
        "expression_manifest": str(authority.manifest_path),
        "expression_manifest_sha256": authority.manifest_sha256,
        "expression_partition_count": sum(
            record.resource_type == "atlas_expression_long" for record in authority.records
        ),
        "metadata_partition_count": sum(
            record.resource_type == "atlas_sample_metadata_wide" for record in authority.records
        ),
        "resource_species_count": len({record.species_column for record in authority.records}),
        "orthofinder_species_count": len(species),
        "scanned_expression_partition_count": len(expression_paths),
        "scanned_metadata_partition_count": len(metadata_paths),
        "alias_count": alias_count,
        "member_count": counts["expression_member_count"],
        "unique_mapped_member_count": counts["expression_unique_mapping_count"],
        "ambiguous_member_count": counts["expression_ambiguous_mapping_count"],
        "not_mapped_member_count": counts["expression_not_mapped_count"],
        "expression_context_count": counts["expression_context_count"],
        "expression_group_count": counts["expression_group_count"],
        "expression_species": ";".join(
            sorted(
                {
                    record.species_column
                    for record in authority.records
                    if record.path in set(expression_paths)
                }
            )
        ),
        "expression_worker_peak_rss_mib": worker_metadata["peak_rss_mib"],
        "expression_worker_database_bytes": worker_metadata["database_bytes"],
        "expression_worker_spill_bytes": worker_metadata["spill_bytes"],
        "minimum_expression_value": minimum_expression_value,
        "broad_positive_fraction": broad_positive_fraction,
        "unit_selection_policy": "TPM when present per species/experiment; otherwise FPKM",
        "mapping_policy": (
            "case-insensitive exact gene_id or gene_name match within exact species; "
            "lowest mapping tier retained"
        ),
        "missing_data_policy": (
            "unmapped, ambiguous and no-expression-record states remain unavailable; "
            "they are never converted to measured zero"
        ),
    }
    write_tsv(
        path=tables_dir / "expression_import_audit.tsv.gz",
        fieldnames=AUDIT_FIELDS,
        records=(audit,),
    )
    return ExpressionPublication(
        counts=counts,
        scanned_expression_paths=expression_paths,
        scanned_metadata_paths=metadata_paths,
    )


def _run_expression_worker(*, request_path: Path, metadata_path: Path) -> Mapping[str, Any]:
    """Run the native DuckDB workload in a fresh, diagnosable process."""

    environment = dict(os.environ)
    environment["MALLOC_ARENA_MAX"] = "2"
    command = (
        str(Path(sys.executable).absolute()),
        "-m",
        "orthofinder_results.expression_analysis",
        "--worker-request",
        str(request_path),
        "--worker-metadata",
        str(metadata_path),
    )
    completed = subprocess.run(command, check=False, env=environment)
    record = _read_expression_worker_metadata(path=metadata_path)
    if record is not None and record.get("status") == "failed":
        message = str(record.get("message", "RNA-seq worker failed without a message."))
        if record.get("error_type") == "InputValidationError":
            raise InputValidationError(message)
        raise PublicationError(message)
    if completed.returncode != 0:
        if completed.returncode < 0:
            status = f"terminated by signal {-completed.returncode}"
        else:
            status = f"exited with status {completed.returncode}"
        raise PublicationError(
            "RNA-seq expression worker "
            f"{status}. Inspect its last logged substage and scheduler memory accounting."
        )
    if record is None or record.get("status") != "complete":
        raise PublicationError(
            f"RNA-seq expression worker did not publish valid completion metadata: {metadata_path}"
        )
    expected_counts = {
        "expression_member_count",
        "expression_unique_mapping_count",
        "expression_ambiguous_mapping_count",
        "expression_not_mapped_count",
        "expression_context_count",
        "expression_group_count",
    }
    counts = record.get("counts")
    if not isinstance(counts, dict) or set(counts) != expected_counts:
        raise PublicationError("RNA-seq worker completion metadata has invalid counts.")
    if any(
        isinstance(value, bool) or not isinstance(value, int) or value < 0
        for value in counts.values()
    ):
        raise PublicationError("RNA-seq worker completion counts must be non-negative integers.")
    for field in ("worker_pid", "database_bytes", "spill_bytes"):
        value = record.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise PublicationError(
                f"RNA-seq worker completion metadata has invalid {field}: {value!r}"
            )
    peak = record.get("peak_rss_mib")
    if isinstance(peak, bool) or not isinstance(peak, (int, float)) or peak < 0:
        raise PublicationError(
            f"RNA-seq worker completion metadata has invalid peak_rss_mib: {peak!r}"
        )
    return record


def _read_expression_worker_metadata(*, path: Path) -> Mapping[str, Any] | None:
    """Read an optional success or controlled-failure worker record."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        return None
    try:
        record = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PublicationError(f"Invalid RNA-seq worker metadata: {source}") from error
    if not isinstance(record, dict):
        raise PublicationError(f"RNA-seq worker metadata is not an object: {source}")
    if record.get("metadata_version") != _WORKER_METADATA_VERSION:
        raise PublicationError("RNA-seq worker metadata version is unsupported.")
    return record


def _publish_expression_worker(*, request: Mapping[str, Any]) -> Mapping[str, Any]:
    """Execute one file-backed RNA-seq integration worker."""

    tables_dir = Path(str(request["tables_dir"])).expanduser().resolve()
    work_dir = Path(str(request["work_dir"])).expanduser().resolve()
    run_id = str(request["run_id"])
    expression_paths = tuple(Path(str(path)).resolve() for path in request["expression_paths"])
    all_expression_paths = tuple(
        Path(str(path)).resolve() for path in request["all_expression_paths"]
    )
    metadata_paths = tuple(Path(str(path)).resolve() for path in request["metadata_paths"])
    minimum_expression_value = float(request["minimum_expression_value"])
    broad_positive_fraction = float(request["broad_positive_fraction"])
    threads = int(request["threads"])
    memory_limit_mb = int(request["memory_limit_mb"])
    _validate_controls(
        minimum_expression_value=minimum_expression_value,
        broad_positive_fraction=broad_positive_fraction,
        threads=threads,
        memory_limit_mb=memory_limit_mb,
    )
    if not tables_dir.is_dir() or not work_dir.is_dir():
        raise InputValidationError("RNA-seq worker input directories are unavailable.")

    sequence_path = tables_dir / "sequences.tsv.gz"
    alias_path = tables_dir / "expression_identifier_aliases.tsv.gz"
    hog_path = tables_dir / "hog_memberships.tsv.gz"
    legacy_path = tables_dir / "legacy_orthogroup_memberships.tsv.gz"
    database_path = work_dir / "expression_work.duckdb"
    spill = work_dir / "spill"
    spill.mkdir(parents=False, exist_ok=False)
    context_path = tables_dir / "expression_context.parquet"
    connection = duckdb.connect(str(database_path))
    try:
        connection.execute(f"SET threads = {threads}")
        connection.execute(f"SET memory_limit = '{memory_limit_mb}MB'")
        connection.execute("SET temp_directory = ?", [str(spill)])
        connection.execute("SET preserve_insertion_order = false")
        with _expression_substage(
            name="source_views",
            database_path=database_path,
            spill_path=spill,
        ):
            _create_source_views(
                connection=connection,
                expression_paths=expression_paths,
                all_expression_paths=all_expression_paths,
                metadata_paths=metadata_paths,
            )
        with _expression_substage(
            name="source_validation",
            database_path=database_path,
            spill_path=spill,
        ):
            _validate_atlas_rows(connection=connection)
        _create_expression_tables(
            connection=connection,
            sequence_path=sequence_path,
            alias_path=alias_path,
            hog_path=hog_path,
            legacy_path=legacy_path,
            context_path=context_path,
            minimum_expression_value=minimum_expression_value,
            broad_positive_fraction=broad_positive_fraction,
            database_path=database_path,
            spill_path=spill,
        )
        counts = {
            "expression_member_count": _row_count(connection, "member_mapping"),
            "expression_unique_mapping_count": _row_count(
                connection,
                "member_mapping",
                "mapping_status = 'MAPPED_UNIQUE'",
            ),
            "expression_ambiguous_mapping_count": _row_count(
                connection,
                "member_mapping",
                "mapping_status = 'AMBIGUOUS'",
            ),
            "expression_not_mapped_count": _row_count(
                connection,
                "member_mapping",
                "mapping_status = 'NOT_MAPPED'",
            ),
            "expression_context_count": _scalar_int(
                connection,
                "SELECT count(*) FROM read_parquet(" + _sql_literal(str(context_path)) + ")",
            ),
            "expression_group_count": _row_count(connection, "group_expression_summary"),
        }
        with _expression_substage(
            name="mapping_tsv",
            database_path=database_path,
            spill_path=spill,
        ):
            _write_query(
                connection=connection,
                query="SELECT * FROM member_mapping ORDER BY species_label, member_id",
                path=tables_dir / "expression_member_mapping.tsv.gz",
                fieldnames=MAPPING_FIELDS,
            )
        with _expression_substage(
            name="member_summary_tsv",
            database_path=database_path,
            spill_path=spill,
        ):
            _write_query(
                connection=connection,
                query="SELECT * FROM member_summary ORDER BY species_label, member_id",
                path=tables_dir / "expression_member_summary.tsv.gz",
                fieldnames=SUMMARY_FIELDS,
            )
        with _expression_substage(
            name="group_summary_tsv",
            database_path=database_path,
            spill_path=spill,
        ):
            _write_query(
                connection=connection,
                query=(
                    "SELECT * FROM group_expression_summary ORDER BY group_type, "
                    "hierarchy_node, group_id"
                ),
                path=tables_dir / "expression_group_summary.tsv.gz",
                fieldnames=GROUP_SUMMARY_FIELDS,
            )
        connection.execute("CHECKPOINT")
    except duckdb.Error as error:
        raise PublicationError(f"RNA-seq expression integration failed: {error}") from error
    finally:
        connection.close()
    database_bytes = database_path.stat().st_size if database_path.is_file() else 0
    spill_bytes = _directory_size_bytes(path=spill)
    _LOGGER.info(
        "RNA-seq worker complete: contexts=%s, groups=%s, peak_rss_mib=%.1f, "
        "database_bytes=%s, spill_bytes=%s",
        f"{counts['expression_context_count']:,}",
        f"{counts['expression_group_count']:,}",
        peak_rss_mib(),
        f"{database_bytes:,}",
        f"{spill_bytes:,}",
    )
    return {
        "metadata_version": _WORKER_METADATA_VERSION,
        "status": "complete",
        "run_id": run_id,
        "counts": counts,
        "worker_pid": os.getpid(),
        "peak_rss_mib": peak_rss_mib(),
        "database_bytes": database_bytes,
        "spill_bytes": spill_bytes,
    }


@contextmanager
def _expression_substage(*, name: str, database_path: Path, spill_path: Path) -> Iterable[None]:
    """Log elapsed time, resident memory and disk use around one large operation."""

    started = time.perf_counter()
    _LOGGER.info(
        "RNA-seq substage started: %s | rss_mib=%.1f | database_bytes=%s | spill_bytes=%s",
        name,
        peak_rss_mib(),
        f"{database_path.stat().st_size if database_path.is_file() else 0:,}",
        f"{_directory_size_bytes(path=spill_path):,}",
    )
    try:
        yield
    except Exception:
        _LOGGER.exception(
            "RNA-seq substage failed: %s | rss_mib=%.1f | elapsed_seconds=%.3f",
            name,
            peak_rss_mib(),
            time.perf_counter() - started,
        )
        raise
    _LOGGER.info(
        "RNA-seq substage finished: %s | rss_mib=%.1f | database_bytes=%s | "
        "spill_bytes=%s | elapsed_seconds=%.3f",
        name,
        peak_rss_mib(),
        f"{database_path.stat().st_size if database_path.is_file() else 0:,}",
        f"{_directory_size_bytes(path=spill_path):,}",
        time.perf_counter() - started,
    )


def _directory_size_bytes(*, path: Path) -> int:
    """Return the total size of regular files below a small work directory."""

    source = Path(path)
    if not source.is_dir():
        return 0
    return sum(item.stat().st_size for item in source.rglob("*") if item.is_file())


def peak_rss_mib() -> float:
    """Return this process's peak resident set size in MiB."""

    maximum = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    divisor = 1024.0 * 1024.0 if sys.platform == "darwin" else 1024.0
    return maximum / divisor


def _validate_controls(
    *,
    minimum_expression_value: float,
    broad_positive_fraction: float,
    threads: int,
    memory_limit_mb: int,
) -> None:
    """Validate bounded expression controls before opening large resources."""

    if not math.isfinite(minimum_expression_value) or minimum_expression_value < 0:
        raise InputValidationError("Minimum expression value must be finite and non-negative.")
    if not math.isfinite(broad_positive_fraction) or not 0 <= broad_positive_fraction <= 1:
        raise InputValidationError("Broad-expression fraction must be between zero and one.")
    if not 1 <= threads <= 256:
        raise InputValidationError("Expression threads must be between 1 and 256.")
    if not 256 <= memory_limit_mb <= 1_048_576:
        raise InputValidationError("Expression memory limit must be between 256 and 1,048,576 MiB.")


def _iter_aliases(
    *,
    sequence_path: Path,
    run_id: str,
    additional_aliases_path: Path | None,
    allowed_species: set[str] | None = None,
) -> Iterable[Mapping[str, Any]]:
    """Yield deterministic exact aliases for each OrthoFinder member."""

    selected_species = (
        {value.upper() for value in allowed_species} if allowed_species is not None else None
    )
    additional: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    if additional_aliases_path is not None:
        source = Path(additional_aliases_path).expanduser().resolve()
        rows = tuple(read_tsv(path=source))
        if rows:
            missing = sorted(set(ADDITIONAL_ALIAS_FIELDS).difference(rows[0]))
            if missing:
                raise InputValidationError(
                    "Expression alias TSV lacks fields: " + "; ".join(missing)
                )
        for row in rows:
            species = str(row.get("species_label", "")).strip()
            if selected_species is not None and species.upper() not in selected_species:
                continue
            member = str(row.get("member_id", "")).strip()
            alias = str(row.get("alias_value", "")).strip()
            try:
                tier = int(str(row.get("mapping_tier", "")))
            except ValueError as error:
                raise InputValidationError(
                    f"Expression alias mapping_tier is not an integer: {row.get('mapping_tier')!r}"
                ) from error
            if not species or not member or not alias or not 1 <= tier <= 100:
                raise InputValidationError(
                    "Expression aliases require species_label, member_id, alias_value and "
                    "mapping_tier between 1 and 100."
                )
            additional[(species, member)].append(dict(row))
    for row in read_tsv(path=sequence_path):
        species = str(row.get("species_label", "")).strip()
        if selected_species is not None and species.upper() not in selected_species:
            continue
        member = str(row.get("member_id", "")).strip()
        internal = str(row.get("internal_id", "")).strip()
        raw_header = str(row.get("raw_header", "")).strip()
        candidates = _header_aliases(member_id=member, raw_header=raw_header)
        for extra in additional.get((species, member), ()):
            candidates.append(
                (
                    str(extra["alias_type"]).strip() or "reviewed_alias",
                    str(extra["alias_value"]).strip(),
                    int(extra["mapping_tier"]),
                    str(extra["source"]).strip() or "reviewed_alias_tsv",
                )
            )
        unique: dict[str, tuple[str, str, int, str]] = {}
        for alias_type, alias_value, tier, source in candidates:
            value = alias_value.strip()
            if not value or "\x00" in value:
                continue
            key = value.upper()
            existing = unique.get(key)
            candidate = (alias_type, value, tier, source)
            if existing is None or (tier, alias_type, source) < (
                existing[2],
                existing[0],
                existing[3],
            ):
                unique[key] = candidate
        for key in sorted(unique):
            alias_type, alias_value, tier, source = unique[key]
            yield {
                "run_id": run_id,
                "species_label": species,
                "internal_id": internal,
                "member_id": member,
                "alias_type": alias_type,
                "alias_value": alias_value,
                "mapping_tier": tier,
                "source": source,
            }


def _header_aliases(*, member_id: str, raw_header: str) -> list[tuple[str, str, int, str]]:
    """Return exact identifiers encoded by one source protein header."""

    aliases: list[tuple[str, str, int, str]] = []
    if member_id:
        aliases.append(("member_id", member_id, 1, "SequenceIDs.txt"))
    primary = raw_header.split(maxsplit=1)[0] if raw_header else ""
    if primary and primary != member_id:
        aliases.append(("raw_primary_token", primary, 1, "SequenceIDs.txt"))
    pipe = _PIPE_IDENTIFIER.match(primary or member_id)
    if pipe is not None:
        aliases.extend(
            (
                ("uniprot_accession", pipe.group(1), 1, "SequenceIDs.txt"),
                ("uniprot_entry", pipe.group(2), 2, "SequenceIDs.txt"),
            )
        )
    for field, value in _HEADER_ALIAS.findall(raw_header):
        aliases.append((field.lower(), value, 1, "SequenceIDs.txt"))
    for alias_type, value, tier, source in tuple(aliases):
        if re.search(r"\.\d+$", value):
            aliases.append(
                (
                    f"{alias_type}_versionless",
                    re.sub(r"\.\d+$", "", value),
                    tier + 1,
                    source,
                )
            )
    return aliases


def _create_source_views(
    *,
    connection: duckdb.DuckDBPyConnection,
    expression_paths: Sequence[Path],
    all_expression_paths: Sequence[Path],
    metadata_paths: Sequence[Path],
) -> None:
    """Validate source schemas and create pruned temporary views."""

    if not all_expression_paths:
        raise InputValidationError("Expression authority contains no expression partitions.")
    schema_query = _parquet_query(paths=all_expression_paths)
    expression_columns = {
        str(row[0]) for row in connection.execute(f"DESCRIBE {schema_query}").fetchall()
    }
    required = {
        "source_database",
        "experiment_accession",
        "species_column",
        "gene_id",
        "gene_name",
        "sample_or_condition",
        "expression_value",
        "expression_minimum",
        "expression_lower_quartile",
        "expression_median",
        "expression_upper_quartile",
        "expression_maximum",
        "expression_value_statistic",
        "expression_summary_type",
        "expression_unit",
        "source_file",
        "source_file_sha256",
    }
    missing = sorted(required.difference(expression_columns))
    if missing:
        raise InputValidationError("Expression Parquet lacks fields: " + "; ".join(missing))
    if expression_paths:
        connection.execute(
            f"CREATE TEMP VIEW atlas_expression AS {_parquet_query(paths=expression_paths)}"
        )
    else:
        connection.execute(
            "CREATE TEMP VIEW atlas_expression AS SELECT "
            "NULL::VARCHAR source_database, NULL::VARCHAR experiment_accession, "
            "NULL::VARCHAR species_column, NULL::VARCHAR gene_id, NULL::VARCHAR gene_name, "
            "NULL::VARCHAR sample_or_condition, NULL::DOUBLE expression_value, "
            "NULL::DOUBLE expression_minimum, NULL::DOUBLE expression_lower_quartile, "
            "NULL::DOUBLE expression_median, NULL::DOUBLE expression_upper_quartile, "
            "NULL::DOUBLE expression_maximum, NULL::VARCHAR expression_value_statistic, "
            "NULL::VARCHAR expression_summary_type, NULL::VARCHAR expression_unit, "
            "NULL::VARCHAR source_file, NULL::VARCHAR source_file_sha256 WHERE FALSE"
        )
    if metadata_paths:
        metadata_query = _parquet_query(paths=metadata_paths)
        metadata_columns = {
            str(row[0]) for row in connection.execute(f"DESCRIBE {metadata_query}").fetchall()
        }
        required_metadata = {
            "experiment_accession",
            "species_column",
            "sample_or_condition",
            "atlas_group_label",
            "assay_ids",
            "assay_count",
            "organism_part",
            "developmental_stage",
            "genotype",
            "cultivar",
            "treatment",
            "condition",
            "source_file",
            "source_file_sha256",
            "expression_file_sha256",
        }
        missing_metadata = sorted(required_metadata.difference(metadata_columns))
        if missing_metadata:
            raise InputValidationError(
                "Expression metadata Parquet lacks fields: " + "; ".join(missing_metadata)
            )
        connection.execute(f"CREATE TEMP VIEW atlas_metadata AS {metadata_query}")
    else:
        connection.execute(
            "CREATE TEMP VIEW atlas_metadata AS SELECT "
            "NULL::VARCHAR experiment_accession, NULL::VARCHAR species_column, "
            "NULL::VARCHAR sample_or_condition, NULL::VARCHAR atlas_group_label, "
            "NULL::VARCHAR assay_ids, NULL::BIGINT assay_count, "
            "NULL::VARCHAR organism_part, NULL::VARCHAR developmental_stage, "
            "NULL::VARCHAR genotype, NULL::VARCHAR cultivar, NULL::VARCHAR treatment, "
            "NULL::VARCHAR condition, NULL::VARCHAR source_file, "
            "NULL::VARCHAR source_file_sha256, NULL::VARCHAR expression_file_sha256 "
            "WHERE FALSE"
        )


def _validate_atlas_rows(*, connection: duckdb.DuckDBPyConnection) -> None:
    """Reject mixed, duplicated or scientifically invalid expression records."""

    invalid_units = _scalar_int(
        connection,
        "SELECT count(*) FROM atlas_expression WHERE "
        "upper(CAST(expression_unit AS VARCHAR)) NOT IN ('TPM', 'FPKM')",
    )
    if invalid_units:
        raise InputValidationError(
            f"Expression input contains {invalid_units:,} rows outside TPM/FPKM."
        )
    invalid_values = _scalar_int(
        connection,
        "SELECT count(*) FROM atlas_expression WHERE expression_value < 0 "
        "OR NOT isfinite(expression_value)",
    )
    if invalid_values:
        raise InputValidationError(
            f"Expression input contains {invalid_values:,} negative or non-finite values."
        )
    duplicate_metadata = _scalar_int(
        connection,
        "SELECT count(*) FROM (SELECT species_column, experiment_accession, "
        "sample_or_condition, count(*) n FROM atlas_metadata GROUP BY ALL HAVING n <> 1)",
    )
    if duplicate_metadata:
        raise InputValidationError(
            f"Expression metadata contains {duplicate_metadata:,} duplicated context keys."
        )


def _create_expression_tables(
    *,
    connection: duckdb.DuckDBPyConnection,
    sequence_path: Path,
    alias_path: Path,
    hog_path: Path,
    legacy_path: Path,
    context_path: Path,
    minimum_expression_value: float,
    broad_positive_fraction: float,
    database_path: Path,
    spill_path: Path,
) -> None:
    """Create normalised mappings and summaries with disk-backed intermediates."""

    sequences = _tsv_query(path=sequence_path)
    aliases = _tsv_query(path=alias_path)
    hogs = _tsv_query(path=hog_path)
    legacy = _tsv_query(path=legacy_path)
    _execute_expression_sql(
        connection=connection,
        name="unit_selection",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            "CREATE TABLE unit_selection AS SELECT "
            "upper(CAST(species_column AS VARCHAR)) species_key, "
            "CAST(experiment_accession AS VARCHAR) experiment_accession, "
            "CASE WHEN bool_or(upper(CAST(expression_unit AS VARCHAR)) = 'TPM') THEN 'TPM' "
            "WHEN bool_or(upper(CAST(expression_unit AS VARCHAR)) = 'FPKM') THEN 'FPKM' END "
            "selected_expression_unit FROM atlas_expression GROUP BY ALL"
        ),
    )
    _execute_expression_sql(
        connection=connection,
        name="preferred_expression_view",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            "CREATE TEMP VIEW expression_preferred AS SELECT e.*, u.selected_expression_unit "
            "FROM atlas_expression e JOIN unit_selection u ON "
            "upper(CAST(e.species_column AS VARCHAR)) = u.species_key AND "
            "CAST(e.experiment_accession AS VARCHAR) = u.experiment_accession "
            "WHERE upper(CAST(e.expression_unit AS VARCHAR)) = u.selected_expression_unit"
        ),
    )
    _execute_expression_sql(
        connection=connection,
        name="atlas_gene_catalogue",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            "CREATE TABLE atlas_gene_catalogue AS SELECT DISTINCT "
            "upper(CAST(species_column AS VARCHAR)) species_key, "
            "CAST(gene_id AS VARCHAR) gene_id, "
            "upper(CAST(gene_id AS VARCHAR)) gene_id_key, "
            "CAST(gene_name AS VARCHAR) gene_name, "
            "upper(COALESCE(CAST(gene_name AS VARCHAR), '')) gene_name_key "
            "FROM atlas_expression WHERE gene_id IS NOT NULL "
            "AND CAST(gene_id AS VARCHAR) <> ''"
        ),
    )
    _execute_expression_sql(
        connection=connection,
        name="candidate_members",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            f"CREATE TABLE candidates AS SELECT DISTINCT s.run_id, s.species_label, "
            f"s.member_id FROM {sequences} s JOIN "
            "(SELECT DISTINCT species_key FROM atlas_gene_catalogue) scoped "
            "ON upper(s.species_label) = scoped.species_key"
        ),
    )
    _execute_expression_sql(
        connection=connection,
        name="candidate_aliases",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            f"CREATE TABLE aliases AS SELECT run_id, species_label, member_id, "
            "alias_type, alias_value, TRY_CAST(mapping_tier AS INTEGER) "
            f"mapping_tier FROM {aliases}"
        ),
    )
    _execute_expression_sql(
        connection=connection,
        name="matched_genes",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            "CREATE TABLE matched_genes AS SELECT DISTINCT a.run_id, a.species_label, "
            "a.member_id, a.alias_type, a.alias_value, a.mapping_tier, "
            "g.gene_id, g.gene_name FROM aliases a JOIN atlas_gene_catalogue g ON "
            "upper(a.species_label) = g.species_key AND "
            "upper(a.alias_value) = g.gene_id_key "
            "UNION SELECT DISTINCT a.run_id, a.species_label, a.member_id, a.alias_type, "
            "a.alias_value, a.mapping_tier, g.gene_id, g.gene_name "
            "FROM aliases a JOIN atlas_gene_catalogue g ON "
            "upper(a.species_label) = g.species_key AND "
            "upper(a.alias_value) = g.gene_name_key WHERE g.gene_name_key <> ''"
        ),
    )
    _execute_expression_sql(
        connection=connection,
        name="member_mapping",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            "CREATE TABLE member_mapping AS WITH best_tier AS (SELECT run_id, "
            "species_label, member_id, min(mapping_tier) mapping_tier FROM matched_genes "
            "GROUP BY ALL), best AS (SELECT m.* FROM matched_genes m JOIN best_tier b USING "
            "(run_id, species_label, member_id, mapping_tier)), summaries AS (SELECT run_id, "
            "species_label, member_id, min(mapping_tier) mapping_tier, "
            "count(DISTINCT gene_id) matched_gene_count, "
            "string_agg(DISTINCT gene_id, ';' ORDER BY gene_id) matched_gene_ids, "
            "string_agg(DISTINCT COALESCE(gene_name, ''), ';' ORDER BY COALESCE(gene_name, '')) "
            "matched_gene_names, string_agg(DISTINCT alias_value, ';' ORDER BY alias_value) "
            "matched_aliases FROM best GROUP BY run_id, species_label, member_id) "
            "SELECT c.run_id, c.species_label, c.member_id, CASE WHEN s.matched_gene_count IS NULL "
            "THEN 'NOT_MAPPED' WHEN s.matched_gene_count = 1 THEN 'MAPPED_UNIQUE' ELSE 'AMBIGUOUS' "
            "END mapping_status, COALESCE(CAST(s.mapping_tier AS VARCHAR), '') mapping_tier, "
            "COALESCE(s.matched_gene_count, 0) matched_gene_count, "
            "COALESCE(s.matched_gene_ids, '') matched_gene_ids, "
            "COALESCE(s.matched_gene_names, '') matched_gene_names, "
            "COALESCE(s.matched_aliases, '') matched_aliases, CASE "
            "WHEN s.matched_gene_count IS NULL THEN 'no_exact_species_scoped_identifier_match' "
            "WHEN s.matched_gene_count = 1 THEN 'unique_best_tier_exact_match' "
            "ELSE 'multiple_genes_at_best_mapping_tier' END reason FROM candidates c "
            "LEFT JOIN summaries s USING (run_id, species_label, member_id)"
        ),
    )
    threshold = float(minimum_expression_value)
    broad = float(broad_positive_fraction)
    context_query = (
        "WITH unique_mapping AS (SELECT run_id, "
        "species_label, member_id, matched_gene_ids gene_id FROM member_mapping "
        "WHERE mapping_status = 'MAPPED_UNIQUE') SELECT m.run_id, m.species_label, "
        "m.member_id, m.gene_id, CAST(e.gene_name AS VARCHAR) gene_name, "
        "CAST(e.source_database AS VARCHAR) source_database, "
        "CAST(e.experiment_accession AS VARCHAR) experiment_accession, "
        "CAST(e.selected_expression_unit AS VARCHAR) expression_unit, "
        "CAST(e.sample_or_condition AS VARCHAR) sample_or_condition, "
        "COALESCE(CAST(meta.atlas_group_label AS VARCHAR), '') atlas_group_label, "
        "COALESCE(CAST(meta.assay_ids AS VARCHAR), '') assay_ids, "
        "COALESCE(TRY_CAST(meta.assay_count AS BIGINT), 0) assay_count, "
        "COALESCE(CAST(meta.organism_part AS VARCHAR), '') organism_part, "
        "COALESCE(CAST(meta.developmental_stage AS VARCHAR), '') developmental_stage, "
        "COALESCE(CAST(meta.genotype AS VARCHAR), '') genotype, "
        "COALESCE(CAST(meta.cultivar AS VARCHAR), '') cultivar, "
        "COALESCE(CAST(meta.treatment AS VARCHAR), '') treatment, "
        "COALESCE(CAST(meta.condition AS VARCHAR), '') condition, "
        "COALESCE(NULLIF(CAST(meta.organism_part AS VARCHAR), ''), "
        "NULLIF(CAST(meta.condition AS VARCHAR), ''), "
        "NULLIF(CAST(meta.atlas_group_label AS VARCHAR), ''), "
        "CAST(e.sample_or_condition AS VARCHAR), 'UNSPECIFIED') expression_context, "
        "CASE WHEN meta.sample_or_condition IS NULL THEN 'METADATA_NOT_MAPPED' "
        "WHEN COALESCE(CAST(meta.organism_part AS VARCHAR), '') = '' THEN "
        "'MAPPED_WITHOUT_TISSUE' ELSE 'MAPPED_WITH_TISSUE' END metadata_status, "
        "CAST(e.expression_value_statistic AS VARCHAR) expression_value_statistic, "
        "CAST(e.expression_summary_type AS VARCHAR) expression_summary_type, "
        "CAST(e.expression_value AS DOUBLE) expression_value, "
        "CAST(e.expression_minimum AS DOUBLE) expression_minimum, "
        "CAST(e.expression_lower_quartile AS DOUBLE) expression_lower_quartile, "
        "CAST(e.expression_median AS DOUBLE) expression_median, "
        "CAST(e.expression_upper_quartile AS DOUBLE) expression_upper_quartile, "
        "CAST(e.expression_maximum AS DOUBLE) expression_maximum, "
        f"CAST(e.expression_value AS DOUBLE) >= {threshold!r} expression_positive, "
        "CAST(e.source_file AS VARCHAR) expression_source_file, "
        "CAST(e.source_file_sha256 AS VARCHAR) expression_source_file_sha256, "
        "COALESCE(CAST(meta.source_file AS VARCHAR), '') metadata_source_file, "
        "COALESCE(CAST(meta.source_file_sha256 AS VARCHAR), '') metadata_source_file_sha256 "
        "FROM unique_mapping m JOIN expression_preferred e ON "
        "upper(CAST(e.species_column AS VARCHAR)) = upper(m.species_label) AND "
        "upper(CAST(e.gene_id AS VARCHAR)) = upper(m.gene_id) LEFT JOIN atlas_metadata meta "
        "ON upper(CAST(meta.species_column AS VARCHAR)) = upper(m.species_label) AND "
        "CAST(meta.experiment_accession AS VARCHAR) = CAST(e.experiment_accession AS VARCHAR) "
        "AND CAST(meta.sample_or_condition AS VARCHAR) = CAST(e.sample_or_condition AS VARCHAR) "
        "AND CAST(meta.expression_file_sha256 AS VARCHAR) = CAST(e.source_file_sha256 AS VARCHAR)"
    )
    with _expression_substage(
        name="context_parquet",
        database_path=database_path,
        spill_path=spill_path,
    ):
        _write_query_parquet(
            connection=connection,
            query=context_query,
            path=context_path,
            fieldnames=CONTEXT_FIELDS,
        )
    _execute_expression_sql(
        connection=connection,
        name="context_view",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            "CREATE TEMP VIEW member_context AS SELECT * FROM read_parquet("
            + _sql_literal(str(Path(context_path).resolve()))
            + ")"
        ),
    )
    _execute_expression_sql(
        connection=connection,
        name="member_summary",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            "CREATE TABLE member_summary AS WITH evidence AS (SELECT run_id, species_label, "
            "member_id, max(gene_id) gene_id, max(gene_name) gene_name, "
            "count(DISTINCT experiment_accession) experiment_count, "
            "string_agg(DISTINCT expression_unit, ';' ORDER BY expression_unit) "
            "selected_expression_units, count(DISTINCT expression_unit) expression_unit_count, "
            "count(*) context_count, count(*) FILTER (WHERE expression_positive) "
            "positive_context_count, avg(CASE WHEN expression_positive THEN 1.0 ELSE 0.0 END) "
            "positive_context_fraction, min(expression_value) minimum_context_expression_value, "
            "max(expression_value) maximum_context_expression_value, median(expression_value) "
            "median_context_expression_value FROM member_context GROUP BY run_id, species_label, "
            "member_id) SELECT m.run_id, m.species_label, m.member_id, m.mapping_status, "
            "COALESCE(e.gene_id, '') gene_id, COALESCE(e.gene_name, '') gene_name, "
            "CASE WHEN m.mapping_status = 'MAPPED_UNIQUE' THEN COALESCE(e.experiment_count, 0) END "
            "experiment_count, CASE WHEN m.mapping_status = 'MAPPED_UNIQUE' THEN "
            "COALESCE(e.selected_expression_units, '') END selected_expression_units, "
            "CASE WHEN m.mapping_status = 'MAPPED_UNIQUE' THEN "
            "COALESCE(e.expression_unit_count, 0) "
            "END expression_unit_count, CASE WHEN m.mapping_status = 'MAPPED_UNIQUE' THEN "
            "COALESCE(e.context_count, 0) END context_count, CASE WHEN m.mapping_status = "
            "'MAPPED_UNIQUE' THEN COALESCE(e.positive_context_count, 0) END "
            "positive_context_count, "
            "e.positive_context_fraction, e.minimum_context_expression_value, "
            "e.maximum_context_expression_value, e.median_context_expression_value, CASE "
            f"WHEN e.positive_context_fraction >= {broad!r} THEN true "
            "WHEN e.context_count IS NOT NULL THEN false END broad_expression_supported, CASE "
            "WHEN m.mapping_status <> 'MAPPED_UNIQUE' THEN m.mapping_status "
            "WHEN e.context_count IS NULL THEN 'NO_EXPRESSION_RECORDS' "
            f"WHEN e.positive_context_fraction >= {broad!r} THEN 'BROAD_EXPRESSION_SUPPORTED' "
            "ELSE 'LIMITED_OR_ZERO_EXPRESSION' END evidence_status FROM member_mapping m "
            "LEFT JOIN evidence e USING (run_id, species_label, member_id)"
        ),
    )
    _execute_expression_sql(
        connection=connection,
        name="scoped_membership_view",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            f"CREATE TEMP VIEW all_memberships AS SELECT m.run_id, m.group_type, "
            f"m.hierarchy_node, m.group_id, m.species_label, m.member_id FROM {hogs} m "
            "JOIN (SELECT DISTINCT species_key FROM atlas_gene_catalogue) scoped "
            "ON upper(m.species_label) = scoped.species_key UNION ALL SELECT m.run_id, "
            "m.group_type, m.hierarchy_node, m.group_id, m.species_label, "
            f"m.member_id FROM {legacy} m "
            "JOIN (SELECT DISTINCT species_key FROM atlas_gene_catalogue) scoped "
            "ON upper(m.species_label) = scoped.species_key"
        ),
    )
    _execute_expression_sql(
        connection=connection,
        name="group_expression_summary",
        database_path=database_path,
        spill_path=spill_path,
        query=(
            "CREATE TABLE group_expression_summary AS SELECT g.run_id, g.group_type, "
            "g.hierarchy_node, g.group_id, count(*) member_count, "
            "count(*) FILTER (WHERE s.mapping_status = 'MAPPED_UNIQUE') "
            "unique_mapped_member_count, "
            "count(*) FILTER (WHERE s.mapping_status = 'AMBIGUOUS') ambiguous_member_count, "
            "count(*) FILTER (WHERE s.mapping_status IS NULL OR s.mapping_status = 'NOT_MAPPED') "
            "not_mapped_member_count, count(*) FILTER (WHERE s.context_count > 0) "
            "expression_observed_member_count, count(*) FILTER (WHERE "
            "s.broad_expression_supported) broad_expression_member_count, "
            "count(DISTINCT g.species_label) FILTER (WHERE s.mapping_status = 'MAPPED_UNIQUE') "
            "mapped_species_count, count(DISTINCT g.species_label) FILTER "
            "(WHERE s.context_count > 0) "
            "expression_observed_species_count, unique_mapped_member_count::DOUBLE / "
            "nullif(member_count, 0) mapping_fraction, expression_observed_member_count::DOUBLE / "
            "nullif(unique_mapped_member_count, 0) expression_observed_fraction, "
            "COALESCE(string_agg(DISTINCT s.selected_expression_units, ';' ORDER BY "
            "s.selected_expression_units) FILTER (WHERE s.selected_expression_units <> ''), '') "
            "selected_expression_units FROM all_memberships g LEFT JOIN member_summary s USING "
            "(run_id, species_label, member_id) GROUP BY g.run_id, g.group_type, "
            "g.hierarchy_node, g.group_id"
        ),
    )


def _execute_expression_sql(
    *,
    connection: duckdb.DuckDBPyConnection,
    name: str,
    query: str,
    database_path: Path,
    spill_path: Path,
) -> None:
    """Execute one trusted expression SQL statement with resource diagnostics."""

    with _expression_substage(
        name=name,
        database_path=database_path,
        spill_path=spill_path,
    ):
        connection.execute(query)


def _parquet_query(*, paths: Sequence[Path]) -> str:
    """Return a controlled DuckDB Parquet scan over verified paths."""

    if not paths:
        raise InputValidationError("At least one Expression Atlas Parquet path is required.")
    literals = ", ".join("'" + str(path).replace("'", "''") + "'" for path in paths)
    return f"SELECT * FROM read_parquet([{literals}], union_by_name=true, hive_partitioning=true)"


def _tsv_query(*, path: Path) -> str:
    """Return a controlled all-VARCHAR DuckDB TSV scan."""

    literal = str(Path(path).resolve()).replace("'", "''")
    return (
        f"read_csv('{literal}', delim='\\t', header=true, all_varchar=true, "
        "compression='gzip', quote='\"', escape='\"')"
    )


def _write_query(
    *,
    connection: duckdb.DuckDBPyConnection,
    query: str,
    path: Path,
    fieldnames: tuple[str, ...],
) -> int:
    """Stream one DuckDB result into the package's unquoted TSV authority."""

    cursor = connection.execute(query)
    observed = tuple(str(column[0]) for column in cursor.description)
    if observed != fieldnames:
        raise PublicationError(
            f"Expression query columns differ: observed={observed}; expected={fieldnames}"
        )
    batches = cursor.to_arrow_reader(65_536)

    def records() -> Iterable[Mapping[str, Any]]:
        for batch in batches:
            yield from batch.to_pylist()

    return write_tsv(path=path, fieldnames=fieldnames, records=records())


def _write_query_parquet(
    *,
    connection: duckdb.DuckDBPyConnection,
    query: str,
    path: Path,
    fieldnames: tuple[str, ...],
) -> int:
    """Materialise one large typed query directly as compressed Parquet.

    This path deliberately bypasses TSV for high-volume relations. Text fields may
    legitimately contain tabs, quotes or newlines, all of which remain scalar Parquet
    values rather than becoming structural delimiters during a second parse.

    Args:
        connection: Connection owning the trusted internal query relations.
        query: Trusted internal SELECT statement with deterministic column order.
        path: Destination Parquet path inside the staging resource.
        fieldnames: Required ordered output headings.

    Returns:
        Number of materialised rows.

    Raises:
        PublicationError: If the query schema differs or publication is incomplete.
    """

    cursor = connection.execute(f"SELECT * FROM ({query}) AS source LIMIT 0")
    observed = tuple(str(column[0]) for column in cursor.description)
    if observed != fieldnames:
        raise PublicationError(
            f"Expression query columns differ: observed={observed}; expected={fieldnames}"
        )
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        destination.unlink()
    connection.execute(
        f"COPY ({query}) TO ? (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880)",
        [str(destination)],
    )
    if not destination.is_file() or destination.stat().st_size == 0:
        raise PublicationError(f"Expression Parquet publication is empty: {destination}")
    return _scalar_int(
        connection,
        "SELECT count(*) FROM read_parquet(" + _sql_literal(str(destination)) + ")",
    )


def _sql_literal(value: str) -> str:
    """Return one safely quoted SQL string literal for a trusted file path."""

    return "'" + str(value).replace("'", "''") + "'"


def _row_count(connection: duckdb.DuckDBPyConnection, relation: str, condition: str = "") -> int:
    """Return one exact temporary-relation row count."""

    where = f" WHERE {condition}" if condition else ""
    return _scalar_int(connection, f"SELECT count(*) FROM {relation}{where}")


def _scalar_int(connection: duckdb.DuckDBPyConnection, query: str) -> int:
    """Execute a trusted internal scalar query and return an integer."""

    return int(connection.execute(query).fetchone()[0])


def _remove_empty_directory(*, path: Path) -> None:
    """Remove a DuckDB spill directory only when it is empty."""

    try:
        path.rmdir()
    except OSError:
        pass


def _load_worker_request(*, path: Path) -> Mapping[str, Any]:
    """Read and minimally validate one private RNA-seq worker request."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise InputValidationError(f"RNA-seq worker request is unavailable: {source}")
    try:
        request = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InputValidationError(f"RNA-seq worker request is invalid: {source}") from error
    if not isinstance(request, dict):
        raise InputValidationError("RNA-seq worker request must be a JSON object.")
    if request.get("metadata_version") != _WORKER_METADATA_VERSION:
        raise InputValidationError("RNA-seq worker request version is unsupported.")
    required = {
        "tables_dir",
        "work_dir",
        "run_id",
        "expression_paths",
        "all_expression_paths",
        "metadata_paths",
        "minimum_expression_value",
        "broad_positive_fraction",
        "threads",
        "memory_limit_mb",
    }
    missing = sorted(required.difference(request))
    if missing:
        raise InputValidationError("RNA-seq worker request lacks fields: " + "; ".join(missing))
    return request


def build_worker_parser() -> argparse.ArgumentParser:
    """Return the named-option parser for the private RNA-seq worker."""

    parser = argparse.ArgumentParser(
        prog="python -m orthofinder_results.expression_analysis",
        description="Run one process-isolated, disk-backed RNA-seq integration.",
    )
    parser.add_argument("--worker-request", required=True, type=Path)
    parser.add_argument("--worker-metadata", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the private expression worker and return its process status."""

    arguments = build_worker_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )
    metadata_path = Path(arguments.worker_metadata).expanduser().resolve()
    try:
        request = _load_worker_request(path=arguments.worker_request)
        record = _publish_expression_worker(request=request)
    except (
        OrthoFinderResultsError,
        OSError,
        UnicodeError,
        ValueError,
        TypeError,
        KeyError,
    ) as error:
        _LOGGER.exception("RNA-seq expression worker failed: %s", error)
        atomic_write_json(
            path=metadata_path,
            record={
                "metadata_version": _WORKER_METADATA_VERSION,
                "status": "failed",
                "error_type": type(error).__name__,
                "message": str(error),
                "worker_pid": os.getpid(),
                "peak_rss_mib": peak_rss_mib(),
            },
        )
        return 2
    atomic_write_json(path=metadata_path, record=record)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

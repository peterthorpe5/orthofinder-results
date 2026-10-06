"""Build a complete-proteome sequence sidecar for terminal-motif queries."""

from __future__ import annotations

import argparse
import logging
import os
import sys
import uuid
from collections.abc import Iterator, Sequence
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from orthofinder_results.errors import InputValidationError, OrthoFinderResultsError
from orthofinder_results.io_utils import configure_logging
from orthofinder_results.parsers import iter_sequence_ids, read_species_ids

_LOGGER = logging.getLogger("orthofinder_interrogation_app.terminal_motif_cli")
_FIELDS = (
    "internal_id",
    "species_label",
    "source_fasta",
    "member_id",
    "raw_header",
    "sequence",
)
_DEFAULT_BATCH_SIZE = 25_000
_PROGRESS_INTERVAL = 250_000


def build_parser() -> argparse.ArgumentParser:
    """Return the all-named-option command-line parser."""

    parser = argparse.ArgumentParser(
        prog="orthofinder-terminal-motif-build",
        description="Build an atomic Parquet sequence sidecar from an OrthoFinder run.",
    )
    parser.add_argument("--orthofinder-results-dir", required=True, type=Path)
    parser.add_argument("--output-parquet", required=True, type=Path)
    parser.add_argument("--log-file", type=Path)
    parser.add_argument("--verbose", action="store_true")
    return parser


def build_sequence_sidecar(
    *,
    results_dir: Path,
    output_path: Path,
    batch_size: int = _DEFAULT_BATCH_SIZE,
) -> int:
    """Reconcile OrthoFinder identifiers and FASTA sequences into Parquet.

    Args:
        results_dir: Completed OrthoFinder ``Results_*`` directory or its
            ``WorkingDirectory``.
        output_path: New Parquet destination.
        batch_size: Maximum protein records converted to Arrow at once.

    Returns:
        Number of published protein sequences.

    Raises:
        InputValidationError: If identifiers or files cannot be reconciled.
    """

    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size < 1:
        raise InputValidationError("Sequence publication batch size must be a positive integer.")
    root = Path(results_dir).expanduser().resolve()
    working = root if root.name == "WorkingDirectory" else root / "WorkingDirectory"
    species_path = working / "SpeciesIDs.txt"
    sequence_path = working / "SequenceIDs.txt"
    species_files, _ = read_species_ids(path=species_path, run_id="sequence_sidecar")
    by_internal: dict[str, tuple[str, str, str, str]] = {}
    for row in iter_sequence_ids(
        path=sequence_path,
        run_id="sequence_sidecar",
        species_by_index=species_files,
    ):
        internal_id = str(row["internal_id"])
        if internal_id in by_internal:
            raise InputValidationError(
                "SequenceIDs.txt contains duplicate internal identifiers."
            )
        by_internal[internal_id] = (
            str(row["species_label"]),
            str(row["source_fasta"]),
            str(row["member_id"]),
            str(row["raw_header"]),
        )
    expected_count = len(by_internal)
    destination = Path(output_path).expanduser().resolve()
    if destination.suffix.lower() != ".parquet":
        raise InputValidationError("--output-parquet must end in .parquet.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    incoming = destination.with_name(f".{destination.name}.incoming.{uuid.uuid4().hex}")
    writer: pq.ParquetWriter | None = None
    published_ids: set[str] = set()
    unexpected_count = 0
    published_count = 0
    batch: list[dict[str, str]] = []
    try:
        schema = pa.schema([(name, pa.string()) for name in _FIELDS])
        writer = pq.ParquetWriter(incoming, schema=schema, compression="zstd")
        for species_index, fasta_name in sorted(species_files.items()):
            fasta_path = working / f"Species{species_index}.fa"
            if not fasta_path.is_file():
                alternative = working / fasta_name
                fasta_path = alternative if alternative.is_file() else fasta_path
            for internal_id, sequence in _iter_fasta_records(path=fasta_path):
                if internal_id in published_ids:
                    raise InputValidationError(f"Duplicate FASTA identifier: {internal_id}")
                published_ids.add(internal_id)
                metadata = by_internal.pop(internal_id, None)
                if metadata is None:
                    unexpected_count += 1
                    continue
                species_label, source_fasta, member_id, raw_header = metadata
                batch.append(
                    {
                        "internal_id": internal_id,
                        "species_label": species_label,
                        "source_fasta": source_fasta,
                        "member_id": member_id,
                        "raw_header": raw_header,
                        "sequence": sequence.rstrip("*"),
                    }
                )
                if len(batch) >= batch_size:
                    published_count += _write_sequence_batch(
                        writer=writer,
                        schema=schema,
                        records=batch,
                    )
                    batch.clear()
                    if published_count % _PROGRESS_INTERVAL == 0:
                        _LOGGER.info(
                            "Sequence sidecar progress: proteins=%s",
                            f"{published_count:,}",
                        )
        if by_internal or unexpected_count:
            raise InputValidationError(
                "Sequence/identifier reconciliation failed: "
                f"missing={len(by_internal):,}; unexpected={unexpected_count:,}."
            )
        if batch:
            published_count += _write_sequence_batch(
                writer=writer,
                schema=schema,
                records=batch,
            )
            batch.clear()
        if published_count != expected_count or published_count == 0:
            raise InputValidationError(
                "Sequence/identifier reconciliation failed: "
                f"expected={expected_count:,}; published={published_count:,}."
            )
        writer.close()
        writer = None
        os.replace(incoming, destination)
    finally:
        if writer is not None:
            writer.close()
        if incoming.exists():
            incoming.unlink()
    _LOGGER.info(
        "Published sequence sidecar: proteins=%s path=%s",
        published_count,
        destination,
    )
    return published_count


def _iter_fasta_records(*, path: Path) -> Iterator[tuple[str, str]]:
    """Yield validated FASTA records without retaining a complete proteome in memory."""

    source = Path(path).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise InputValidationError(f"Missing or empty FASTA file: {source}")
    current: str | None = None
    sequence_parts: list[str] = []
    observed: set[str] = set()
    record_count = 0
    with source.open(mode="r", encoding="utf-8", errors="strict") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    sequence = "".join(sequence_parts)
                    if not sequence:
                        raise InputValidationError(
                            f"FASTA contains an empty sequence for {current!r}: {source}"
                        )
                    yield current, sequence
                    record_count += 1
                identifier_text = line[1:].strip()
                identifier = identifier_text.split(maxsplit=1)[0] if identifier_text else ""
                if not identifier:
                    raise InputValidationError(
                        f"Empty FASTA identifier at line {line_number} in {source}."
                    )
                if identifier in observed:
                    raise InputValidationError(
                        f"Duplicate FASTA identifier {identifier!r} in {source}."
                    )
                observed.add(identifier)
                current = identifier
                sequence_parts = []
                continue
            if current is None:
                raise InputValidationError(
                    "Sequence precedes the first FASTA header at line "
                    f"{line_number} in {source}."
                )
            sequence_parts.append("".join(line.split()).upper())
    if current is not None:
        sequence = "".join(sequence_parts)
        if not sequence:
            raise InputValidationError(
                f"FASTA contains an empty sequence for {current!r}: {source}"
            )
        yield current, sequence
        record_count += 1
    if record_count == 0:
        raise InputValidationError(f"FASTA contains no complete sequences: {source}")


def _write_sequence_batch(
    *,
    writer: pq.ParquetWriter,
    schema: pa.Schema,
    records: list[dict[str, str]],
) -> int:
    """Write one bounded protein batch and return its record count."""

    writer.write_table(pa.Table.from_pylist(records, schema=schema))
    return len(records)


def main(argv: Sequence[str] | None = None) -> int:
    """Build the sidecar and return a process status."""

    arguments = build_parser().parse_args(argv)
    configure_logging(log_path=arguments.log_file, verbose=arguments.verbose)
    try:
        build_sequence_sidecar(
            results_dir=arguments.orthofinder_results_dir,
            output_path=arguments.output_parquet,
        )
    except (OrthoFinderResultsError, OSError, ValueError) as error:
        _LOGGER.error("Terminal-motif sidecar build failed: %s", error)
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

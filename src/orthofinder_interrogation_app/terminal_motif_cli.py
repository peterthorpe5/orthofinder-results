"""Build a complete-proteome sequence sidecar for terminal-motif queries."""

from __future__ import annotations

import argparse
import logging
import os
import sys
import uuid
from collections.abc import Sequence
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from orthofinder_results.distances import read_fasta
from orthofinder_results.errors import InputValidationError, OrthoFinderResultsError
from orthofinder_results.io_utils import configure_logging
from orthofinder_results.parsers import iter_sequence_ids, read_species_ids

_LOGGER = logging.getLogger("orthofinder_interrogation_app.terminal_motif_cli")
_FIELDS = ("internal_id", "species_label", "member_id", "sequence")


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


def build_sequence_sidecar(*, results_dir: Path, output_path: Path) -> int:
    """Reconcile OrthoFinder identifiers and FASTA sequences into Parquet.

    Args:
        results_dir: Completed OrthoFinder ``Results_*`` directory or its
            ``WorkingDirectory``.
        output_path: New Parquet destination.

    Returns:
        Number of published protein sequences.

    Raises:
        InputValidationError: If identifiers or files cannot be reconciled.
    """

    root = Path(results_dir).expanduser().resolve()
    working = root if root.name == "WorkingDirectory" else root / "WorkingDirectory"
    species_path = working / "SpeciesIDs.txt"
    sequence_path = working / "SequenceIDs.txt"
    species_files, _ = read_species_ids(path=species_path, run_id="sequence_sidecar")
    aliases = tuple(
        iter_sequence_ids(
            path=sequence_path,
            run_id="sequence_sidecar",
            species_by_index=species_files,
        )
    )
    by_internal = {str(row["internal_id"]): row for row in aliases}
    if len(by_internal) != len(aliases):
        raise InputValidationError("SequenceIDs.txt contains duplicate internal identifiers.")
    sequences: dict[str, str] = {}
    for species_index, fasta_name in sorted(species_files.items()):
        fasta_path = working / f"Species{species_index}.fa"
        if not fasta_path.is_file():
            alternative = working / fasta_name
            fasta_path = alternative if alternative.is_file() else fasta_path
        for internal_id, sequence in read_fasta(path=fasta_path).items():
            if internal_id in sequences:
                raise InputValidationError(f"Duplicate FASTA identifier: {internal_id}")
            sequences[internal_id] = sequence
    missing = sorted(set(by_internal).difference(sequences))
    unexpected = sorted(set(sequences).difference(by_internal))
    if missing or unexpected:
        raise InputValidationError(
            "Sequence/identifier reconciliation failed: "
            f"missing={len(missing):,}; unexpected={len(unexpected):,}."
        )
    records = [
        {
            "internal_id": internal_id,
            "species_label": str(by_internal[internal_id]["species_label"]),
            "member_id": str(by_internal[internal_id]["member_id"]),
            "sequence": sequences[internal_id].upper().rstrip("*"),
        }
        for internal_id in sorted(by_internal)
    ]
    destination = Path(output_path).expanduser().resolve()
    if destination.suffix.lower() != ".parquet":
        raise InputValidationError("--output-parquet must end in .parquet.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    incoming = destination.with_name(f".{destination.name}.incoming.{uuid.uuid4().hex}")
    try:
        schema = pa.schema([(name, pa.string()) for name in _FIELDS])
        table = pa.Table.from_pylist(records, schema=schema)
        pq.write_table(table, incoming, compression="zstd")
        os.replace(incoming, destination)
    finally:
        if incoming.exists():
            incoming.unlink()
    _LOGGER.info("Published sequence sidecar: proteins=%s path=%s", len(records), destination)
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

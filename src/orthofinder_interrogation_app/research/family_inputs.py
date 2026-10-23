"""Audited ingestion of selected local proteins into a reusable research job."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Mapping

from .analysis import index_fasta
from .common import ResearchConfig, ResearchError, write_tsv

LOGGER = logging.getLogger(__name__)


def import_family(
    *, config: ResearchConfig, source: Path, output: Path, parameters: Mapping[str, object]
) -> dict[str, Path]:
    """Validate a complete selected FASTA and retain its OrthoFinder membership context.

    Args:
        config: Worker profile and sequence limit.
        source: Private checksummed input snapshot.
        output: Current owned job directory.
        parameters: Optional source_context records and source_label.

    Returns:
        Normalised FASTA, summary, and optional exact source membership table.

    Raises:
        ResearchError: If source identifiers repeat or context disagrees with the FASTA.
    """
    destination = output / "family.fasta"
    if destination.exists() or (output / "family.sqlite").exists():
        raise ResearchError("Refusing to overwrite an existing family import")
    context = parameters.get("source_context", [])
    if not isinstance(context, list) or any(not isinstance(row, dict) for row in context):
        raise ResearchError("source_context must be a list of membership records")
    connection, count = index_fasta(
        source=source, destination=output / "family.sqlite", maximum=config.max_sequences
    )
    try:
        if context:
            if any(not isinstance(row.get("input_id"), str) for row in context):
                raise ResearchError("Source context requires exact input_id values")
            identifiers = {row[0] for row in connection.execute("SELECT id FROM sequences")}
            if identifiers != {row["input_id"] for row in context}:
                raise ResearchError("Source context and imported protein IDs disagree")
        with destination.open(mode="x", encoding="utf-8") as handle:
            for identifier, sequence in connection.execute(
                "SELECT id,sequence FROM sequences ORDER BY id"
            ):
                handle.write(f">{identifier}\n{sequence}\n")
        minimum, maximum, mean = connection.execute(
            "SELECT min(length(sequence)),max(length(sequence)),avg(length(sequence)) "
            "FROM sequences"
        ).fetchone()
    finally:
        connection.close()
    summary = {
        "sequence_count": count,
        "minimum_length": minimum,
        "maximum_length": maximum,
        "mean_length": mean,
        "scope": "SELECTED_LOCAL_PROTEINS",
        "source_label": str(parameters.get("source_label", "Uploaded protein family")),
    }
    path = output / "family_summary.tsv"
    write_tsv(path=path, rows=[summary], columns=list(summary))
    artifacts = {destination.name: destination, path.name: path}
    if context:
        path = output / "source_context.tsv"
        columns = sorted({name for row in context for name in row})
        write_tsv(path=path, rows=context, columns=columns)
        artifacts[path.name] = path
    LOGGER.info("Imported %d complete local proteins", count)
    return artifacts

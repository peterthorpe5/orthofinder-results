"""Exact tab-separated downloads for bounded application results."""

from __future__ import annotations

import csv
import io
from collections.abc import Mapping, Sequence
from typing import Any

from orthofinder_results.errors import InputValidationError


def records_to_tsv(
    *, records: Sequence[Mapping[str, Any]], fieldnames: Sequence[str] | None = None
) -> bytes:
    """Serialise records as UTF-8 tab-separated text.

    Args:
        records: Ordered tabular records.
        fieldnames: Optional explicit ordered headings. The first record defines
            headings when omitted.

    Returns:
        UTF-8 encoded TSV content, including a header.

    Raises:
        InputValidationError: If headings are unavailable or inconsistent.
    """

    headings = tuple(fieldnames or (tuple(records[0]) if records else ()))
    if not headings or any(not heading for heading in headings):
        raise InputValidationError("TSV output requires non-empty field names.")
    if len(set(headings)) != len(headings):
        raise InputValidationError("TSV output field names must be unique.")
    expected = set(headings)
    for index, record in enumerate(records, start=1):
        unknown = sorted(set(record) - expected)
        if unknown:
            raise InputValidationError(
                f"TSV row {index} contains undeclared fields: {'; '.join(unknown)}"
            )
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=headings,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for record in records:
        writer.writerow({heading: _scalar(record.get(heading)) for heading in headings})
    return buffer.getvalue().encode("utf-8")


def _scalar(value: Any) -> Any:
    """Convert optional and Boolean values to unambiguous TSV scalars."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return value

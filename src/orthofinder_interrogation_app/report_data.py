"""Validated access to bounded visual data embedded in an offline report."""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from orthofinder_results.errors import InputValidationError

_LOGGER = logging.getLogger("orthofinder_interrogation_app.report_data")
_START_MARKER = b'<script id="orthofinder-results-data" type="application/json">'
_END_MARKER = b"</script>"
MAX_REPORT_BYTES = 250 * 1024 * 1024
MAX_VISUAL_GROUPS = 100
MAX_VISUAL_MEMBERS = 1_000
MAX_VISUAL_EDGES = 100_000


@dataclass(frozen=True)
class VisualisationCatalog:
    """Bounded visual records from one checksum-bound offline report.

    Attributes:
        run_id: Run identifier declared by the report payload.
        networks: Visual records keyed by composite group identity.
        nearest_neighbours: Requested nearest neighbours per displayed member.
    """

    run_id: str
    networks: dict[str, dict[str, Any]]
    nearest_neighbours: int

    def labels(self) -> tuple[str, ...]:
        """Return deterministic visual-group labels."""

        return tuple(sorted(self.networks))

    def get(self, *, label: str) -> dict[str, Any]:
        """Return one validated visual entry.

        Args:
            label: Composite report group label.

        Returns:
            Visual record for the requested group.

        Raises:
            InputValidationError: If the group is not in the bounded catalog.
        """

        try:
            return self.networks[label]
        except KeyError as error:
            raise InputValidationError(
                f"Visual group is not available in this report: {label}"
            ) from error


def load_visualisation_catalog(*, report_path: Path, expected_run_id: str) -> VisualisationCatalog:
    """Load and validate the bounded JSON payload from an offline HTML report.

    The embedded payload is a compatibility authority for visual records that
    schema-2 resources did not persist separately, most notably pruned resolved
    gene-tree phylograms. It never replaces complete DuckDB group searching.

    Args:
        report_path: Existing self-contained report.
        expected_run_id: Run identifier already validated from DuckDB.

    Returns:
        Validated bounded visualisation catalog.

    Raises:
        InputValidationError: If the report is absent, unsafe or inconsistent.
    """

    source = Path(report_path).expanduser().resolve()
    try:
        size = source.stat().st_size
    except OSError as error:
        raise InputValidationError(f"Offline report is unavailable: {source}: {error}") from error
    if not 1 <= size <= MAX_REPORT_BYTES:
        raise InputValidationError(
            f"Offline report size must be between 1 and {MAX_REPORT_BYTES:,} bytes; "
            f"observed {size:,}: {source}"
        )
    try:
        document = source.read_bytes()
    except OSError as error:
        raise InputValidationError(
            f"Offline report could not be read: {source}: {error}"
        ) from error
    start = document.find(_START_MARKER)
    if start < 0:
        raise InputValidationError(
            "Offline report lacks the orthofinder-results visual data payload."
        )
    payload_start = start + len(_START_MARKER)
    payload_end = document.find(_END_MARKER, payload_start)
    if payload_end < 0:
        raise InputValidationError("Offline report visual data payload is not terminated.")
    try:
        payload = json.loads(document[payload_start:payload_end].decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise InputValidationError(
            f"Offline report visual data payload is invalid JSON: {error}"
        ) from error
    catalog = _validate_payload(payload=payload, expected_run_id=expected_run_id)
    _LOGGER.info(
        "Validated bounded report visuals: run=%s, groups=%s, report=%s",
        catalog.run_id,
        len(catalog.networks),
        source,
    )
    return catalog


def _validate_payload(*, payload: object, expected_run_id: str) -> VisualisationCatalog:
    """Validate top-level report identity and bounded network collections."""

    if not isinstance(payload, dict):
        raise InputValidationError("Offline report visual payload must be a JSON object.")
    run = payload.get("run")
    networks = payload.get("networks")
    limits = payload.get("limits", {})
    if not isinstance(run, dict) or not isinstance(networks, dict):
        raise InputValidationError("Offline report visual payload lacks run or network data.")
    run_id = str(run.get("run_id", ""))
    if run_id != expected_run_id:
        raise InputValidationError(
            f"Offline report and DuckDB run identifiers disagree: {run_id!r} versus "
            f"{expected_run_id!r}."
        )
    if len(networks) > MAX_VISUAL_GROUPS:
        raise InputValidationError(
            f"Offline report contains {len(networks):,} visual groups; the safe limit is "
            f"{MAX_VISUAL_GROUPS:,}."
        )
    validated: dict[str, dict[str, Any]] = {}
    for raw_label, raw_entry in networks.items():
        label = str(raw_label)
        if not label or not isinstance(raw_entry, dict):
            raise InputValidationError("Every report visual group must have a label and object.")
        _validate_entry(label=label, entry=raw_entry)
        validated[label] = raw_entry
    nearest_neighbours = 0
    if isinstance(limits, dict):
        try:
            nearest_neighbours = int(limits.get("nearestNeighbours", 0))
        except (TypeError, ValueError):
            nearest_neighbours = 0
    if nearest_neighbours < 0 or nearest_neighbours > 100:
        raise InputValidationError(
            f"Report nearest-neighbour limit is unsafe: {nearest_neighbours}."
        )
    return VisualisationCatalog(
        run_id=run_id,
        networks=validated,
        nearest_neighbours=nearest_neighbours,
    )


def _validate_entry(*, label: str, entry: dict[str, Any]) -> None:
    """Reject malformed or unexpectedly large per-group visual records."""

    members = entry.get("members")
    nodes = entry.get("nodes")
    edges = entry.get("edges")
    if not isinstance(members, list) or not isinstance(nodes, list) or not isinstance(edges, list):
        raise InputValidationError(f"Visual group {label} lacks member, node or edge arrays.")
    if len(members) > MAX_VISUAL_MEMBERS or len(nodes) > MAX_VISUAL_MEMBERS:
        raise InputValidationError(
            f"Visual group {label} exceeds the {MAX_VISUAL_MEMBERS:,}-member safety limit."
        )
    if len(edges) > MAX_VISUAL_EDGES:
        raise InputValidationError(
            f"Visual group {label} exceeds the {MAX_VISUAL_EDGES:,}-edge safety limit."
        )
    member_ids = {
        str(row.get("member_id", ""))
        for row in members
        if isinstance(row, dict) and row.get("member_id")
    }
    node_ids = {str(row.get("id", "")) for row in nodes if isinstance(row, dict) and row.get("id")}
    if len(member_ids) != len(members) or len(node_ids) != len(nodes):
        raise InputValidationError(
            f"Visual group {label} contains duplicate or malformed members/nodes."
        )

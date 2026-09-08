"""Lazy, cacheable distance analysis for application-selected groups."""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import sys
import uuid
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Mapping

from orthofinder_results.distances import (
    calculate_patristic_distances_from_newick,
    summarise_distances,
)
from orthofinder_results.errors import InputValidationError
from orthofinder_results.report import build_group_visualisation
from orthofinder_results.trees import decode_newick_payload, normalise_newick_text

from .models import GroupKey
from .queries import (
    MAX_GROUP_DISTANCE_ROWS,
    MAX_GROUP_MEMBER_ROWS,
    OrthoFinderQueryService,
)
from .report_data import VisualisationCatalog

_LOGGER = logging.getLogger("orthofinder_interrogation_app.distance_data")
DEFAULT_ANALYSIS_MEMBERS = 250
MAX_ANALYSIS_MEMBERS = 500
DEFAULT_NEAREST_NEIGHBOURS = 3
MAX_NEAREST_NEIGHBOURS = 20
MAX_CACHE_BYTES = 250 * 1024 * 1024
CACHE_FORMAT_VERSION = 1


@dataclass(frozen=True)
class GroupAnalysis:
    """Complete bounded analysis used by all linked cluster views.

    Attributes:
        key: Composite group identity.
        group: Complete group-statistics record.
        members: Displayed canonical members and species.
        distances: Exact displayed member-to-member rows, including species labels.
        summary: Distribution and sampling summary for those distances.
        visual_entry: Shared PCoA, matrix, phylogram and neighbour payload.
        source: Explicit persisted-report, persisted-DuckDB or lazy-tree authority.
        tree_authority: Resolved or original gene-tree source, when available.
        cache_status: Whether a lazy result was calculated or loaded from sidecar cache.
    """

    key: GroupKey
    group: dict[str, Any]
    members: tuple[dict[str, Any], ...]
    distances: tuple[dict[str, Any], ...]
    summary: dict[str, Any]
    visual_entry: dict[str, Any]
    source: str
    tree_authority: str
    cache_status: str


class DistanceAnalysisProvider:
    """Resolve persisted or portable-tree distances without modifying a resource."""

    def __init__(
        self,
        *,
        service: OrthoFinderQueryService,
        cache_dir: Path,
        report_catalog: VisualisationCatalog | None = None,
    ) -> None:
        """Initialise a provider for one validated immutable resource.

        Args:
            service: Read-only query service.
            cache_dir: Writable sidecar location outside the completed resource.
            report_catalog: Optional schema-2-compatible embedded visual catalogue.

        Raises:
            InputValidationError: If the cache would be created inside the resource.
        """

        self.service = service
        self.cache_dir = _validate_cache_directory(
            cache_dir=cache_dir,
            resource_path=service.resource.resource_path,
        )
        self.report_catalog = report_catalog

    def analyse(
        self,
        *,
        key: GroupKey,
        max_members: int = DEFAULT_ANALYSIS_MEMBERS,
        nearest_neighbours: int = DEFAULT_NEAREST_NEIGHBOURS,
        force_recompute: bool = False,
    ) -> GroupAnalysis:
        """Return one bounded analysis from the best available authority.

        Persisted report matrices preserve the original schema-2 pilot. Persisted
        DuckDB pairs are used next. Schema-3 portable trees then permit deterministic
        on-demand calculation for groups whose distances were not precomputed.

        Args:
            key: Composite run-scoped group identity.
            max_members: Maximum exact or deterministic sampled members.
            nearest_neighbours: Neighbour edges retained per displayed member.
            force_recompute: Ignore a valid sidecar cache entry.

        Returns:
            Complete linked analysis record.

        Raises:
            InputValidationError: If controls or resource capabilities are insufficient.
        """

        _validate_analysis_controls(
            max_members=max_members,
            nearest_neighbours=nearest_neighbours,
        )
        if key.run_id != self.service.resource.run_id:
            raise InputValidationError(
                f"Group run {key.run_id!r} does not match resource run "
                f"{self.service.resource.run_id!r}."
            )
        group = self.service.get_group(key=key)
        report_entry = self._report_entry(key=key)
        report_matrix = report_entry.get("distanceMatrix") if report_entry is not None else None
        report_is_exact = (
            report_entry is not None
            and isinstance(report_matrix, dict)
            and report_matrix.get("status") == "EXACT_COMPLETE_DISPLAYED_MATRIX"
        )
        if self.service.resource.schema_version < 3 and report_is_exact:
            return _analysis_from_report(
                key=key,
                group=group,
                entry=report_entry,
            )

        method = str(group.get("distance_method", ""))
        persisted = self.service.get_group_distances(
            key=key,
            distance_method=method,
            maximum=MAX_GROUP_DISTANCE_ROWS,
        )
        if persisted:
            return self._analysis_from_rows(
                key=key,
                group=group,
                rows=persisted,
                summary=group,
                nearest_neighbours=nearest_neighbours,
                source="PERSISTED_DUCKDB_DISTANCES",
                tree_record=None,
                newick_text=None,
            )
        if report_is_exact:
            return _analysis_from_report(
                key=key,
                group=group,
                entry=report_entry,
            )

        tree_record = self.service.get_portable_tree(
            key=key,
            legacy_orthogroup_id=str(group.get("legacy_orthogroup_id", "")),
        )
        if tree_record is None:
            if self.service.resource.schema_version < 3:
                reason = (
                    "This schema-2 resource has no portable tree for on-demand distances. "
                    "Rebuild it with orthofinder-results 0.4 or later; the existing 25 "
                    "pilot groups remain available without rebuilding."
                )
            else:
                reason = "No portable gene tree matches this group or its parent orthogroup."
            raise InputValidationError(reason)
        cache_path = self._cache_path(
            key=key,
            tree_record=tree_record,
            max_members=max_members,
            nearest_neighbours=nearest_neighbours,
        )
        if not force_recompute:
            cached = _read_cached_analysis(path=cache_path, expected_key=key)
            if cached is not None:
                _LOGGER.info("Loaded lazy group analysis from cache: %s", cache_path)
                return replace(cached, cache_status="CACHE_HIT")

        newick_text = decode_newick_payload(
            payload=str(tree_record["newick_payload"]),
            payload_encoding=str(tree_record["payload_encoding"]),
            expected_size=int(tree_record["source_size_bytes"]),
            expected_sha256=str(tree_record["source_sha256"]),
        )
        declared_member_count = int(group.get("member_count", 0))
        if declared_member_count > MAX_GROUP_MEMBER_ROWS:
            raise InputValidationError(
                f"Group {key.display_label()} contains {declared_member_count:,} members, "
                f"above the safe browser-analysis limit of {MAX_GROUP_MEMBER_ROWS:,}. "
                "Precompute this group in the cluster pipeline before opening it in the app."
            )
        all_members = self.service.get_group_members(key=key)
        if declared_member_count and len(all_members) != declared_member_count:
            raise InputValidationError(
                f"Membership authority is incomplete for {key.display_label()}: expected "
                f"{declared_member_count:,} rows but loaded {len(all_members):,}."
            )
        member_ids, aliases = _member_aliases(
            members=all_members,
            sequence_aliases=self.service.get_group_sequence_aliases(key=key),
        )
        _LOGGER.info(
            "Lazy patristic calculation started: group=%s, members=%s, maximum=%s, tree=%s",
            key.display_label(),
            len(member_ids),
            max_members,
            tree_record["tree_id"],
        )
        rows, summary = calculate_patristic_distances_from_newick(
            newick_text=newick_text,
            run_id=key.run_id,
            group_type=key.group_type,
            hierarchy_node=key.hierarchy_node,
            group_id=key.group_id,
            max_members=max_members,
            member_ids=member_ids,
            member_aliases=aliases,
            source_file=(
                "portable-tree://"
                f"{tree_record['tree_type']}/{tree_record['tree_id']}"
                f"@sha256:{tree_record['source_sha256']}"
            ),
        )
        analysis = self._analysis_from_rows(
            key=key,
            group=group,
            rows=tuple(rows),
            summary=summary,
            nearest_neighbours=nearest_neighbours,
            source="PORTABLE_TREE_LAZY",
            tree_record=tree_record,
            newick_text=newick_text,
        )
        _write_cached_analysis(path=cache_path, analysis=analysis)
        _LOGGER.info(
            "Lazy patristic calculation finished: group=%s, sampled_members=%s, pairs=%s, "
            "cache=%s",
            key.display_label(),
            summary["sampled_member_count"],
            summary["distance_pair_count"],
            cache_path,
        )
        return analysis

    def _report_entry(self, *, key: GroupKey) -> dict[str, Any] | None:
        """Return a report visual for one exact group key when present."""

        if self.report_catalog is None:
            return None
        return self.report_catalog.networks.get(_catalog_key(key=key))

    def _analysis_from_rows(
        self,
        *,
        key: GroupKey,
        group: Mapping[str, Any],
        rows: tuple[dict[str, Any], ...],
        summary: Mapping[str, Any],
        nearest_neighbours: int,
        source: str,
        tree_record: Mapping[str, Any] | None,
        newick_text: str | None,
    ) -> GroupAnalysis:
        """Build one shared visual and table record from exact pair rows."""

        displayed_ids = {
            str(row[field]) for row in rows for field in ("member_a", "member_b")
        }
        all_members = self.service.get_group_members(key=key)
        members = tuple(
            row for row in all_members if str(row.get("member_id", "")) in displayed_ids
        )
        if len({str(row["member_id"]) for row in members}) != len(displayed_ids):
            raise InputValidationError(
                "Distance members could not be resolved uniquely to group memberships."
            )
        sequence_aliases = self.service.get_group_sequence_aliases(key=key)
        tree_nodes: tuple[dict[str, Any], ...] = ()
        tree_edges: tuple[dict[str, Any], ...] = ()
        tree_authority = ""
        if tree_record is not None and newick_text is not None:
            tree_authority = str(tree_record["tree_type"])
            raw_nodes, raw_edges = normalise_newick_text(
                newick_text=newick_text,
                run_id=key.run_id,
                tree_type="RESOLVED_GENE_TREE",
                tree_id=str(tree_record["tree_id"]),
                source_label=str(tree_record["source_path"]),
            )
            tree_nodes, tree_edges = tuple(raw_nodes), tuple(raw_edges)
        entry = build_group_visualisation(
            group_statistic=group,
            memberships=members,
            distances=rows,
            distance_statistic=summary,
            max_members=max(2, len(displayed_ids)),
            nearest_neighbours=nearest_neighbours,
            tree_nodes=tree_nodes,
            tree_edges=tree_edges,
            sequence_identifiers=sequence_aliases,
        )
        if tree_authority:
            phylogram = entry.get("phylogram")
            if isinstance(phylogram, dict):
                phylogram["treeAuthority"] = tree_authority
                if tree_authority == "GENE_TREE":
                    phylogram["method"] = "GENE_TREE_BRANCH_LENGTH"
        decorated = _decorate_distance_rows(rows=rows, members=members)
        return GroupAnalysis(
            key=key,
            group=dict(group),
            members=tuple(dict(row) for row in members),
            distances=decorated,
            summary=dict(summary),
            visual_entry=entry,
            source=source,
            tree_authority=tree_authority,
            cache_status=("CACHE_WRITE" if source == "PORTABLE_TREE_LAZY" else "NOT_APPLICABLE"),
        )

    def _cache_path(
        self,
        *,
        key: GroupKey,
        tree_record: Mapping[str, Any],
        max_members: int,
        nearest_neighbours: int,
    ) -> Path:
        """Return a content-addressed sidecar cache path."""

        identity = {
            "cache_format": CACHE_FORMAT_VERSION,
            "run_id": key.run_id,
            "group_type": key.group_type,
            "hierarchy_node": key.hierarchy_node,
            "group_id": key.group_id,
            "tree_sha256": str(tree_record["source_sha256"]),
            "max_members": max_members,
            "nearest_neighbours": nearest_neighbours,
        }
        digest = hashlib.sha256(
            json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        return self.cache_dir / key.run_id / f"group_analysis_{digest}.json.gz"


def default_cache_directory() -> Path:
    """Return a persistent per-user cache path without relying on ``/tmp``.

    Returns:
        macOS Library cache or the XDG/Linux user cache equivalent.
    """

    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "orthofinder-results"
    xdg = os.environ.get("XDG_CACHE_HOME", "").strip()
    root = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return root / "orthofinder-results"


def _validate_analysis_controls(*, max_members: int, nearest_neighbours: int) -> None:
    """Validate bounded analysis controls."""

    if not 2 <= max_members <= MAX_ANALYSIS_MEMBERS:
        raise InputValidationError(
            f"Analysis members must be between 2 and {MAX_ANALYSIS_MEMBERS:,}."
        )
    if not 1 <= nearest_neighbours <= MAX_NEAREST_NEIGHBOURS:
        raise InputValidationError(
            f"Nearest neighbours must be between 1 and {MAX_NEAREST_NEIGHBOURS:,}."
        )


def _validate_cache_directory(*, cache_dir: Path, resource_path: Path) -> Path:
    """Resolve a sidecar cache path and reject resource-internal writes."""

    cache = Path(cache_dir).expanduser().resolve()
    resource = Path(resource_path).expanduser().resolve()
    resource_root = resource if resource.is_dir() else resource.parent
    if cache == resource_root or resource_root in cache.parents:
        raise InputValidationError(
            "The analysis cache must be outside the immutable completed resource."
        )
    return cache


def _member_aliases(
    *,
    members: tuple[dict[str, Any], ...],
    sequence_aliases: tuple[dict[str, Any], ...],
) -> tuple[tuple[str, ...], dict[str, dict[str, str]]]:
    """Build unambiguous canonical member aliases for tree-leaf resolution."""

    species_by_member: dict[str, set[str]] = defaultdict(set)
    for row in members:
        species_by_member[str(row["member_id"])].add(str(row["species_label"]))
    ambiguous = {
        member: species for member, species in species_by_member.items() if len(species) != 1
    }
    if ambiguous:
        member, species = sorted(ambiguous.items())[0]
        raise InputValidationError(
            f"Canonical member {member!r} occurs under multiple species: "
            + ";".join(sorted(species))
        )
    internal_by_pair: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in sequence_aliases:
        internal_by_pair[
            (str(row["member_id"]), str(row["species_label"]))
        ].add(str(row["internal_id"]))
    aliases: dict[str, dict[str, str]] = {}
    for member, species_values in species_by_member.items():
        species = next(iter(species_values))
        member_aliases = {
            f"{species}_{member}": "SPECIES_PREFIXED_MEMBER_ID",
        }
        member_aliases.update(
            {
                internal: "ORTHOFINDER_INTERNAL_ID"
                for internal in sorted(internal_by_pair.get((member, species), set()))
                if internal
            }
        )
        aliases[member] = member_aliases
    return tuple(sorted(species_by_member)), aliases


def _decorate_distance_rows(
    *, rows: tuple[dict[str, Any], ...], members: tuple[dict[str, Any], ...]
) -> tuple[dict[str, Any], ...]:
    """Add exact endpoint species to long-form distance rows."""

    species_by_member: dict[str, str] = {}
    for row in members:
        member = str(row["member_id"])
        species = str(row["species_label"])
        existing = species_by_member.get(member)
        if existing is not None and existing != species:
            raise InputValidationError(
                f"Member {member!r} maps to more than one species in this group."
            )
        species_by_member[member] = species
    decorated = []
    for row in rows:
        left, right = str(row["member_a"]), str(row["member_b"])
        if left not in species_by_member or right not in species_by_member:
            raise InputValidationError(
                "A persisted distance endpoint is absent from the group membership authority."
            )
        decorated.append(
            {
                **dict(row),
                "species_a": species_by_member[left],
                "species_b": species_by_member[right],
            }
        )
    return tuple(decorated)


def _analysis_from_report(
    *, key: GroupKey, group: Mapping[str, Any], entry: Mapping[str, Any]
) -> GroupAnalysis:
    """Recover exact displayed rows from a validated schema-2 report matrix."""

    raw_members = entry.get("members")
    matrix = entry.get("distanceMatrix")
    raw_summary = entry.get("distanceSummary")
    if not isinstance(raw_members, list) or not isinstance(matrix, dict):
        raise InputValidationError("The report group lacks members or an exact matrix.")
    if matrix.get("status") != "EXACT_COMPLETE_DISPLAYED_MATRIX":
        raise InputValidationError(
            "The report group does not contain a complete displayed distance matrix."
        )
    order, values = matrix.get("memberOrder"), matrix.get("upperTriangle")
    if not isinstance(order, list) or not isinstance(values, list):
        raise InputValidationError("The report distance matrix is malformed.")
    expected = len(order) * (len(order) - 1) // 2
    if len(values) != expected:
        raise InputValidationError(
            f"The report matrix requires {expected:,} pairs; observed {len(values):,}."
        )
    members = tuple(dict(row) for row in raw_members if isinstance(row, dict))
    method = (
        str(raw_summary.get("distance_method", ""))
        if isinstance(raw_summary, dict)
        else ""
    ) or "report_embedded_pairwise_distance"
    status = (
        str(raw_summary.get("computation_status", ""))
        if isinstance(raw_summary, dict)
        else ""
    ) or "EXACT"
    rows: list[dict[str, Any]] = []
    value_index = 0
    for left_index, left in enumerate(order):
        for right in order[left_index + 1 :]:
            rows.append(
                {
                    "run_id": key.run_id,
                    "group_type": key.group_type,
                    "hierarchy_node": key.hierarchy_node,
                    "group_id": key.group_id,
                    "member_a": str(left),
                    "member_b": str(right),
                    "distance_method": method,
                    "distance": float(values[value_index]),
                    "comparable_sites": "",
                    "mismatch_sites": "",
                    "computation_status": status,
                    "source_file": "offline report exact displayed matrix",
                }
            )
            value_index += 1
    total_members = int(group.get("member_count", len(order)))
    if isinstance(raw_summary, dict):
        total_members = int(raw_summary.get("total_member_count", total_members))
    summary = summarise_distances(
        rows=rows,
        run_id=key.run_id,
        group_type=key.group_type,
        hierarchy_node=key.hierarchy_node,
        group_id=key.group_id,
        method=method,
        status=status,
        total_member_count=total_members,
        sampled_member_count=len(order),
        member_identifier_resolution="REPORT_DISPLAYED_MEMBER_ID",
        source_file="offline report exact displayed matrix",
    )
    return GroupAnalysis(
        key=key,
        group=dict(group),
        members=members,
        distances=_decorate_distance_rows(rows=tuple(rows), members=members),
        summary=summary,
        visual_entry=dict(entry),
        source="PERSISTED_REPORT_MATRIX",
        tree_authority="RESOLVED_GENE_TREE",
        cache_status="NOT_APPLICABLE",
    )


def _catalog_key(*, key: GroupKey) -> str:
    """Return the report's collision-safe group key."""

    return "|".join((key.group_type, key.hierarchy_node, key.group_id))


def _analysis_record(*, analysis: GroupAnalysis) -> dict[str, Any]:
    """Serialise one analysis for a content-addressed sidecar cache."""

    return {
        "cache_format": CACHE_FORMAT_VERSION,
        "key": {
            "run_id": analysis.key.run_id,
            "group_type": analysis.key.group_type,
            "hierarchy_node": analysis.key.hierarchy_node,
            "group_id": analysis.key.group_id,
        },
        "group": analysis.group,
        "members": list(analysis.members),
        "distances": list(analysis.distances),
        "summary": analysis.summary,
        "visual_entry": analysis.visual_entry,
        "source": analysis.source,
        "tree_authority": analysis.tree_authority,
    }


def _read_cached_analysis(*, path: Path, expected_key: GroupKey) -> GroupAnalysis | None:
    """Return a validated cached analysis or ignore a corrupt sidecar safely."""

    if not path.is_file():
        return None
    try:
        if not 1 <= path.stat().st_size <= MAX_CACHE_BYTES:
            raise InputValidationError("Cached analysis has an unsafe compressed size.")
        with gzip.open(path, mode="rb") as handle:
            raw = handle.read(MAX_CACHE_BYTES + 1)
        if len(raw) > MAX_CACHE_BYTES:
            raise InputValidationError("Cached analysis exceeds the safe decoded size.")
        record = json.loads(raw.decode("utf-8"))
        return _cached_analysis_from_record(record=record, expected_key=expected_key)
    except (OSError, UnicodeError, json.JSONDecodeError, InputValidationError) as error:
        _LOGGER.warning("Ignoring invalid analysis cache %s: %s", path, error)
        return None


def _cached_analysis_from_record(
    *, record: object, expected_key: GroupKey
) -> GroupAnalysis:
    """Validate a decoded cache object and reconstruct its immutable model."""

    if not isinstance(record, dict) or record.get("cache_format") != CACHE_FORMAT_VERSION:
        raise InputValidationError("Cached analysis format is unsupported.")
    key_record = record.get("key")
    if not isinstance(key_record, dict):
        raise InputValidationError("Cached analysis lacks a group key.")
    observed_key = GroupKey(
        run_id=str(key_record.get("run_id", "")),
        group_type=str(key_record.get("group_type", "")),
        hierarchy_node=str(key_record.get("hierarchy_node", "")),
        group_id=str(key_record.get("group_id", "")),
    )
    if observed_key != expected_key:
        raise InputValidationError("Cached analysis group identity does not match its request.")
    required_mappings = ("group", "summary", "visual_entry")
    if any(not isinstance(record.get(field), dict) for field in required_mappings):
        raise InputValidationError("Cached analysis lacks a required object.")
    raw_members, raw_distances = record.get("members"), record.get("distances")
    if not isinstance(raw_members, list) or not isinstance(raw_distances, list):
        raise InputValidationError("Cached analysis lacks member or distance rows.")
    if any(not isinstance(row, dict) for row in (*raw_members, *raw_distances)):
        raise InputValidationError("Cached analysis contains a malformed row.")
    if len(raw_distances) > MAX_GROUP_DISTANCE_ROWS:
        raise InputValidationError("Cached analysis contains too many distance rows.")
    return GroupAnalysis(
        key=observed_key,
        group=dict(record["group"]),
        members=tuple(dict(row) for row in raw_members),
        distances=tuple(dict(row) for row in raw_distances),
        summary=dict(record["summary"]),
        visual_entry=dict(record["visual_entry"]),
        source=str(record.get("source", "PORTABLE_TREE_LAZY")),
        tree_authority=str(record.get("tree_authority", "")),
        cache_status="CACHE_HIT",
    )


def _write_cached_analysis(*, path: Path, analysis: GroupAnalysis) -> None:
    """Atomically publish one compressed sidecar cache record."""

    path.parent.mkdir(parents=True, exist_ok=True)
    incoming = path.with_name(f".{path.name}.incoming.{uuid.uuid4().hex}")
    encoded = json.dumps(
        _analysis_record(analysis=analysis),
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    if len(encoded) > MAX_CACHE_BYTES:
        raise InputValidationError(
            f"Calculated analysis exceeds the {MAX_CACHE_BYTES:,}-byte cache limit."
        )
    try:
        with incoming.open(mode="wb") as raw_handle:
            with gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0) as handle:
                handle.write(encoded)
            raw_handle.flush()
            os.fsync(raw_handle.fileno())
        os.chmod(incoming, 0o600)
        os.replace(incoming, path)
    except OSError as error:
        _LOGGER.exception("Could not publish analysis cache: %s", path)
        try:
            incoming.unlink(missing_ok=True)
        except OSError:
            pass
        raise InputValidationError(f"Could not publish analysis cache {path}: {error}") from error

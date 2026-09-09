"""Reproducible selection-coverage orchestration and atomic publication."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import uuid
import zipfile
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from io import BytesIO
from pathlib import Path
from typing import Any, Mapping, Sequence

from orthofinder_results.errors import InputValidationError

from .coverage_tree import (
    CoverageTree,
    build_coverage_tree,
    coverage_records,
    tree_edge_records,
    tree_to_newick,
)
from .coverage_tree_render import style_records, tree_to_pdf, tree_to_svg
from .taxonomy_selection import (
    ExpectedTaxaAuthority,
    GroupTaxonEvaluation,
    SelectionSpec,
    TaxonomyGraph,
    evaluate_group_rows,
    expected_taxa_to_tsv,
    input_taxon_mapping_records,
    selection_predicate_records,
    unmapped_input_records,
)
from .tsv import records_to_tsv

_LOGGER = logging.getLogger("orthofinder_interrogation_app.coverage_exports")
EXPORT_FILENAMES = (
    "taxonomy_nodes.tsv",
    "input_taxon_mapping.tsv",
    "expected_taxa.tsv",
    "selection_predicates.tsv",
    "taxon_coverage.tsv",
    "group_taxon_evaluation.tsv",
    "tree_edges.tsv",
    "unmapped_input_taxa.tsv",
    "selection_coverage_tree.newick",
    "selection_coverage_tree_styles.tsv",
    "selection_coverage_tree.svg",
    "selection_coverage_tree.pdf",
    "selection_manifest.json",
    "checksums.sha256",
)
GROUP_TAXON_EVALUATION_COLUMNS = (
    "run_id",
    "group_type",
    "hierarchy_node",
    "group_id",
    "mapped_member_count",
    "unmapped_member_count",
    "represented_taxon_ids",
    "outside_selected_scope_taxon_ids",
    "unmapped_species_labels",
    "predicate_pass",
    "predicate_failure_reasons",
    "predicate_audit",
    "selection_manifest_id",
)


@dataclass(frozen=True)
class SelectionCoverageRun:
    """One reconciled tree, group audit and immutable provenance manifest."""

    graph: TaxonomyGraph
    expected: ExpectedTaxaAuthority
    selection: SelectionSpec
    evaluations: tuple[GroupTaxonEvaluation, ...]
    tree: CoverageTree
    unmapped_rows: tuple[dict[str, Any], ...]
    run_id: str
    resource_identity: str
    package_version: str
    created_at_utc: str
    selection_manifest_id: str
    focus_authority_name: str
    focus_authority_sha256: str
    focus_limited: bool
    maximum_groups: int
    maximum_nodes: int

    @property
    def passing_evaluations(self) -> tuple[GroupTaxonEvaluation, ...]:
        """Return groups satisfying every active AND-composed predicate."""

        return tuple(row for row in self.evaluations if row.predicate_pass)

    def summary(self) -> dict[str, int]:
        """Return high-level tree, mapping, expected and group counts."""

        represented = sum(node.is_represented for node in self.tree.nodes)
        expected_no_data = sum(
            node.coverage_state == "EXPECTED_NO_DATA" for node in self.tree.nodes
        )
        outside = sum(node.is_outside_selected_scope for node in self.tree.nodes)
        return {
            "displayed_nodes": len(self.tree.nodes),
            "displayed_edges": len(self.tree.edges),
            "represented_terminal_taxa": represented,
            "expected_no_data_taxa": expected_no_data,
            "outside_scope_taxa_with_hits": outside,
            "unmapped_input_labels": len(self.unmapped_rows),
            "evaluated_groups": len(self.evaluations),
            "passing_groups": len(self.passing_evaluations),
        }


def build_selection_coverage_run(
    *,
    graph: TaxonomyGraph,
    expected: ExpectedTaxaAuthority,
    selection: SelectionSpec,
    group_species_rows: Sequence[Mapping[str, Any]],
    run_id: str,
    resource_identity: str,
    package_version: str,
    compact: bool = False,
    maximum_groups: int = 250_000,
    maximum_nodes: int = 5_000,
    focus_authority_name: str = "",
    focus_authority_sha256: str = "",
    focus_limited: bool = False,
    created_at_utc: str = "",
) -> SelectionCoverageRun:
    """Evaluate groups and build one internally reconciled coverage run.

    Args:
        graph: Reviewed bounded taxonomy graph.
        expected: Explicit expected-taxon universe.
        selection: Validated AND-composed predicates.
        group_species_rows: Immutable group-by-species count records.
        run_id: Source OrthoFinder resource identifier.
        resource_identity: Checksum or immutable source identity string.
        package_version: Application package version.
        compact: Whether unary neutral internal nodes are collapsed.
        maximum_groups: Defensive group evaluation bound.
        maximum_nodes: Defensive displayed-node bound.
        focus_authority_name: Optional protein-focus authority name.
        focus_authority_sha256: Optional focus authority checksum.
        focus_limited: Whether group evaluation was restricted to focus clusters.
        created_at_utc: Optional reproducible UTC timestamp.

    Returns:
        Complete selection state used by app and command-line exports.

    Raises:
        InputValidationError: If provenance or timestamps are incomplete.
    """

    for label, value in (
        ("run_id", run_id),
        ("resource_identity", resource_identity),
        ("package_version", package_version),
    ):
        if not isinstance(value, str) or not value.strip():
            raise InputValidationError(f"Selection coverage {label} must be non-empty text.")
    timestamp = created_at_utc or datetime.now(tz=UTC).isoformat().replace("+00:00", "Z")
    _validate_timestamp(value=timestamp)
    evaluations = evaluate_group_rows(
        graph=graph,
        selection=selection,
        group_species_rows=group_species_rows,
        maximum_groups=maximum_groups,
    )
    outside_hits = {
        taxon_id
        for evaluation in evaluations
        if evaluation.predicate_pass
        for taxon_id in evaluation.outside_selected_scope_taxon_ids
    }
    tree = build_coverage_tree(
        graph=graph,
        selection=selection,
        expected=expected,
        outside_hit_taxon_ids=outside_hits,
        compact=compact,
        maximum_nodes=maximum_nodes,
    )
    identity_payload = {
        "run_id": run_id,
        "resource_identity": resource_identity,
        "mapping_sha256": graph.mapping_sha256,
        "expected_sha256": expected.sha256,
        "selection_semantics_version": selection.semantics_version,
        "predicates": selection_predicate_records(selection=selection),
        "compact": compact,
        "focus_authority_sha256": focus_authority_sha256,
        "focus_limited": focus_limited,
    }
    manifest_id = hashlib.sha256(
        json.dumps(
            identity_payload,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()
    evaluations = tuple(
        replace(row, selection_manifest_id=manifest_id) for row in evaluations
    )
    unmapped = unmapped_input_records(
        graph=graph,
        group_species_rows=group_species_rows,
    )
    result = SelectionCoverageRun(
        graph=graph,
        expected=expected,
        selection=selection,
        evaluations=evaluations,
        tree=tree,
        unmapped_rows=unmapped,
        run_id=run_id,
        resource_identity=resource_identity,
        package_version=package_version,
        created_at_utc=timestamp,
        selection_manifest_id=manifest_id,
        focus_authority_name=focus_authority_name,
        focus_authority_sha256=focus_authority_sha256,
        focus_limited=focus_limited,
        maximum_groups=maximum_groups,
        maximum_nodes=maximum_nodes,
    )
    _LOGGER.info(
        "Selection coverage run built: run=%s, manifest=%s, groups=%s, passing=%s, "
        "nodes=%s, focus_limited=%s",
        run_id,
        manifest_id,
        len(evaluations),
        len(result.passing_evaluations),
        len(tree.nodes),
        focus_limited,
    )
    return result


def coverage_export_files(*, run: SelectionCoverageRun) -> dict[str, bytes]:
    """Return every reconciled download file and checksum manifest in memory."""

    files: dict[str, bytes] = {
        "taxonomy_nodes.tsv": records_to_tsv(
            records=_taxonomy_node_records(graph=run.graph)
        ),
        "input_taxon_mapping.tsv": records_to_tsv(
            records=input_taxon_mapping_records(graph=run.graph)
        ),
        "expected_taxa.tsv": expected_taxa_to_tsv(authority=run.expected),
        "selection_predicates.tsv": records_to_tsv(
            records=selection_predicate_records(selection=run.selection),
            fieldnames=(
                "selector_type",
                "taxon_id",
                "taxon_name",
                "taxon_rank",
                "predicate_order",
                "normalised_logical_mode",
                "semantics_version",
            ),
        ),
        "taxon_coverage.tsv": records_to_tsv(records=coverage_records(tree=run.tree)),
        "group_taxon_evaluation.tsv": records_to_tsv(
            records=tuple(row.as_record() for row in run.evaluations),
            fieldnames=GROUP_TAXON_EVALUATION_COLUMNS,
        ),
        "tree_edges.tsv": records_to_tsv(records=tree_edge_records(tree=run.tree)),
        "unmapped_input_taxa.tsv": records_to_tsv(
            records=run.unmapped_rows,
            fieldnames=(
                "workflow_species_label",
                "input_label_count",
                "evaluated_group_member_occurrences",
                "mapping_status",
                "mapping_reason",
                "source_name_original",
            ),
        ),
        "selection_coverage_tree.newick": tree_to_newick(tree=run.tree).encode("utf-8"),
        "selection_coverage_tree_styles.tsv": records_to_tsv(
            records=style_records(tree=run.tree)
        ),
        "selection_coverage_tree.svg": tree_to_svg(tree=run.tree),
        "selection_coverage_tree.pdf": tree_to_pdf(tree=run.tree),
    }
    file_checksums = {
        name: hashlib.sha256(payload).hexdigest()
        for name, payload in sorted(files.items())
    }
    manifest = _selection_manifest(run=run, file_checksums=file_checksums)
    files["selection_manifest.json"] = (
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    checksum_rows = [
        f"{hashlib.sha256(payload).hexdigest()}  {name}"
        for name, payload in sorted(files.items())
    ]
    files["checksums.sha256"] = ("\n".join(checksum_rows) + "\n").encode("utf-8")
    if tuple(sorted(files)) != tuple(sorted(EXPORT_FILENAMES)):
        raise InputValidationError("Selection export file contract is incomplete.")
    return files


def coverage_export_zip(*, run: SelectionCoverageRun) -> bytes:
    """Return a deterministic ZIP containing every selection export."""

    files = coverage_export_files(run=run)
    output = BytesIO()
    with zipfile.ZipFile(output, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, payload in sorted(files.items()):
            info = zipfile.ZipInfo(filename=name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, payload)
    return output.getvalue()


def publish_selection_coverage_run(
    *, run: SelectionCoverageRun, output_dir: Path, dry_run: bool = False
) -> dict[str, Any]:
    """Publish exports beside a persistent target using verified atomic rename.

    Args:
        run: Complete reconciled selection state.
        output_dir: New persistent output directory.
        dry_run: Validate and return the manifest without writing files.

    Returns:
        Complete selection manifest.

    Raises:
        InputValidationError: If the target exists or publication verification fails.
    """

    target = Path(output_dir).expanduser().resolve()
    files = coverage_export_files(run=run)
    manifest = json.loads(files["selection_manifest.json"].decode("utf-8"))
    if dry_run:
        _LOGGER.info("Selection coverage dry run passed: target=%s", target)
        return manifest
    if target.exists():
        raise InputValidationError(
            f"Selection coverage output already exists; choose a new directory: {target}"
        )
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise InputValidationError(
            f"Selection coverage output parent cannot be created: {target.parent}"
        ) from error
    staging = target.parent / f".{target.name}.incoming.{uuid.uuid4().hex}"
    try:
        staging.mkdir(mode=0o700)
        for name, payload in sorted(files.items()):
            path = staging / name
            path.write_bytes(payload)
            if hashlib.sha256(path.read_bytes()).digest() != hashlib.sha256(payload).digest():
                raise InputValidationError(f"Checksum verification failed while writing {name}.")
        os.replace(staging, target)
    except (OSError, InputValidationError) as error:
        shutil.rmtree(staging, ignore_errors=True)
        if isinstance(error, InputValidationError):
            raise
        raise InputValidationError(
            f"Selection coverage publication failed for {target}: {error}"
        ) from error
    _LOGGER.info(
        "Selection coverage exports published atomically: target=%s, files=%s",
        target,
        len(files),
    )
    return manifest


def _taxonomy_node_records(*, graph: TaxonomyGraph) -> tuple[dict[str, Any], ...]:
    """Return stable derived ``taxonomy_nodes`` rows."""

    return tuple(
        {
            "taxon_id": node.taxon_id,
            "name": node.name,
            "rank": node.rank,
            "parent_taxon_id": node.parent_taxon_id or "",
            "lineage_depth": node.lineage_depth,
            "authority": node.authority,
            "release": node.release,
            "is_reviewed_terminal": node.is_terminal,
        }
        for node in graph.nodes
    )


def _selection_manifest(
    *, run: SelectionCoverageRun, file_checksums: Mapping[str, str]
) -> dict[str, Any]:
    """Build the canonical JSON manifest before its own checksum is available."""

    return {
        "manifest_schema_version": 1,
        "selection_manifest_id": run.selection_manifest_id,
        "created_at_utc": run.created_at_utc,
        "run_id": run.run_id,
        "resource_identity": run.resource_identity,
        "package_version": run.package_version,
        "taxonomy": {
            "authority": run.graph.authority,
            "release": run.graph.release,
            "mapping_sha256": run.graph.mapping_sha256,
        },
        "expected_universe": {
            "source_name": run.expected.source_name,
            "sha256": run.expected.sha256,
            "included_taxa": len(run.expected.included_taxon_ids),
        },
        "focus_authority": {
            "limited_to_focus_clusters": run.focus_limited,
            "source_name": run.focus_authority_name,
            "sha256": run.focus_authority_sha256,
        },
        "selection": {
            "semantics_version": run.selection.semantics_version,
            "logical_mode": "AND",
            "only_in_multiple_clade_mode": "INTERSECTION",
            "predicates": selection_predicate_records(selection=run.selection),
        },
        "tree": {
            "displayed_root_taxon_id": run.tree.root_taxon_id,
            "compact_unary_neutral_nodes": run.tree.compact,
            "displayed_nodes": len(run.tree.nodes),
            "displayed_edges": len(run.tree.edges),
        },
        "limits": {
            "maximum_groups": run.maximum_groups,
            "maximum_nodes": run.maximum_nodes,
        },
        "summary": run.summary(),
        "files": dict(sorted(file_checksums.items())),
        "interpretation_caution": (
            "EXPECTED_NO_DATA means not represented in this dataset; it is not evidence "
            "of biological absence. Unmapped input labels are never placed speculatively."
        ),
    }


def _validate_timestamp(*, value: str) -> None:
    """Validate one timezone-aware UTC creation timestamp."""

    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise InputValidationError("Selection creation timestamp is not valid ISO-8601.") from error
    if parsed.tzinfo is None or parsed.utcoffset() != UTC.utcoffset(parsed):
        raise InputValidationError("Selection creation timestamp must be UTC and timezone-aware.")

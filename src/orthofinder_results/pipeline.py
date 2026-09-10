"""Atomic production pipeline for generic OrthoFinder result resources."""

from __future__ import annotations

import csv
import hashlib
import heapq
import json
import logging
import math
import os
import re
import shutil
import time
import uuid
from collections import defaultdict
from collections.abc import Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

from orthofinder_interrogation_app.focus import (
    FocusProteinAuthority,
    read_focus_proteins,
)

from . import __schema_version__, __version__
from .benchmark_analysis import (
    BenchmarkPlan,
    build_benchmark_plan,
    publish_benchmark_inputs,
    publish_benchmark_results,
)
from .benchmark_authority import BenchmarkAuthority, read_benchmark_authority
from .distances import (
    DISTANCE_FIELDS,
    DISTANCE_STATISTIC_FIELDS,
    calculate_alignment_distances,
    calculate_patristic_distances,
    read_fasta,
    summarise_distances,
)
from .errors import DistanceCalculationError, InputValidationError, PublicationError
from .focus_analysis import (
    FocusSelection,
    publish_focus_cluster_results,
    publish_focus_selection,
    select_focus_groups,
)
from .io_utils import (
    atomic_write_json,
    configure_logging,
    create_duckdb,
    file_record,
    open_text,
    read_tsv,
    sha256_file,
    tsv_to_parquet,
    utc_now_iso,
    validate_persistent_path,
    write_tsv,
)
from .layout import discover_layout
from .models import ResultLayout
from .parsers import (
    MEMBERSHIP_FIELDS,
    SEQUENCE_FIELDS,
    SPECIES_FIELDS,
    iter_memberships,
    iter_sequence_ids,
    read_species_ids,
)
from .report import build_interactive_report
from .statistics import (
    GROUP_SPECIES_STATISTIC_FIELDS,
    GROUP_STATISTIC_FIELDS,
    GroupAccumulator,
)
from .trees import (
    TREE_EDGE_FIELDS,
    TREE_INVENTORY_FIELDS,
    TREE_NODE_FIELDS,
    TREE_PAYLOAD_FIELDS,
    iter_portable_tree_payloads,
    iter_tree_inventory,
    normalise_newick_tree,
    tree_id_from_path,
)

_LOGGER = logging.getLogger("orthofinder_results.pipeline")
_RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_ALIGNMENT_SUFFIXES = (".fa", ".faa", ".fasta", ".fas", ".aln")
_MAX_REPORT_STATISTIC_ROWS = 50_000
_MAX_REPORT_DISTANCE_PAIRS = 1_000_000
_STAGE_FIELDS = (
    "stage",
    "status",
    "started_at_utc",
    "finished_at_utc",
    "elapsed_seconds",
    "details",
)

GROUP_TYPES = {
    "group_statistics": {
        "member_count": "int64",
        "species_count": "int64",
        "single_copy_species_count": "int64",
        "max_copies_per_species": "int64",
        "mean_copies_per_species": "float64",
        "is_singleton": "bool",
    },
    "group_species_statistics": {
        "species_member_count": "int64",
        "member_fraction": "float64",
    },
    "species": {"source_line": "int64"},
    "sequences": {"source_line": "int64"},
    "legacy_orthogroup_memberships": {"source_row": "int64"},
    "hog_memberships": {"source_row": "int64"},
    "tree_inventory": {"size_bytes": "int64"},
    "tree_payloads": {"source_size_bytes": "int64"},
    "tree_nodes": {
        "is_leaf": "bool",
        "branch_length": "float64",
        "confidence": "float64",
        "descendant_leaf_count": "int64",
    },
    "tree_edges": {"branch_length": "float64"},
    "pairwise_distances": {
        "distance": "float64",
        "comparable_sites": "int64",
        "mismatch_sites": "int64",
    },
    "distance_statistics": {
        "total_member_count": "int64",
        "sampled_member_count": "int64",
        "distance_pair_count": "int64",
        "unresolved_pair_count": "int64",
        "minimum_distance": "float64",
        "q05_distance": "float64",
        "q25_distance": "float64",
        "median_distance": "float64",
        "mean_distance": "float64",
        "q75_distance": "float64",
        "q95_distance": "float64",
        "maximum_distance": "float64",
        "population_stddev_distance": "float64",
        "mean_comparable_sites": "float64",
    },
    "e3_seed_catalogue_audit": {
        "enabled": "bool",
        "matched_group_count": "int64",
        "matched_member_count": "int64",
        "matched_species_count": "int64",
        "seed_protein_sequence_length": "int64",
    },
    "e3_seed_matches": {
        "seed_protein_sequence_length": "int64",
    },
    "e3_cluster_results": {
        "member_count": "int64",
        "species_count": "int64",
        "single_copy_species_count": "int64",
        "max_copies_per_species": "int64",
        "mean_copies_per_species": "float64",
        "is_singleton": "bool",
        "matched_seed_count": "int64",
        "matched_e3_member_count": "int64",
        "matched_e3_species_count": "int64",
        "total_member_count": "int64",
        "sampled_member_count": "int64",
        "distance_pair_count": "int64",
        "unresolved_pair_count": "int64",
        "minimum_distance": "float64",
        "q05_distance": "float64",
        "q25_distance": "float64",
        "median_distance": "float64",
        "mean_distance": "float64",
        "q75_distance": "float64",
        "q95_distance": "float64",
        "maximum_distance": "float64",
        "population_stddev_distance": "float64",
        "mean_comparable_sites": "float64",
        "distance_sampling_fraction": "float64",
        "expected_sample_pair_count": "int64",
        "resolved_pair_fraction": "float64",
        "full_group_distance_matrix": "bool",
        "distance_interquartile_range": "float64",
        "distance_range": "float64",
        "distance_coefficient_of_variation": "float64",
    },
    "benchmark_marker_audit": {
        "enabled": "bool",
        "matched_group_count": "int64",
        "matched_member_count": "int64",
    },
    "benchmark_group_profiles": {
        "marker_count": "int64",
        "matched_member_count": "int64",
    },
    "benchmark_matched_controls": {
        "control_rank": "int64",
        "matching_score": "float64",
        "anchor_member_count": "int64",
        "control_member_count": "int64",
        "anchor_species_count": "int64",
        "control_species_count": "int64",
        "anchor_mean_copies_per_species": "float64",
        "control_mean_copies_per_species": "float64",
        "anchor_single_copy_fraction": "float64",
        "control_single_copy_fraction": "float64",
        "control_reused": "bool",
    },
    "benchmark_cluster_results": {
        "member_count": "int64",
        "species_count": "int64",
        "single_copy_species_count": "int64",
        "max_copies_per_species": "int64",
        "mean_copies_per_species": "float64",
        "total_member_count": "int64",
        "sampled_member_count": "int64",
        "distance_pair_count": "int64",
        "unresolved_pair_count": "int64",
        "minimum_distance": "float64",
        "q05_distance": "float64",
        "q25_distance": "float64",
        "median_distance": "float64",
        "mean_distance": "float64",
        "q75_distance": "float64",
        "q95_distance": "float64",
        "maximum_distance": "float64",
        "population_stddev_distance": "float64",
        "distance_interquartile_range": "float64",
        "distance_coefficient_of_variation": "float64",
        "sampling_fraction": "float64",
        "pairwise_rows_persisted": "bool",
    },
    "benchmark_background_statistics": {
        "eligible_group_count": "int64",
        "minimum_value": "float64",
        "q25_value": "float64",
        "median_value": "float64",
        "mean_value": "float64",
        "q75_value": "float64",
        "maximum_value": "float64",
        "population_stddev_value": "float64",
    },
    "benchmark_contrasts": {
        "target_group_count": "int64",
        "reference_group_count": "int64",
        "target_median": "float64",
        "reference_median": "float64",
        "median_difference": "float64",
        "median_difference_ci_low": "float64",
        "median_difference_ci_high": "float64",
        "cliffs_delta": "float64",
        "mann_whitney_u": "float64",
        "p_value_two_sided": "float64",
        "fdr_q_value": "float64",
    },
    "benchmark_individual_comparisons": {
        "observed_value": "float64",
        "background_group_count": "int64",
        "background_median": "float64",
        "difference_from_background_median": "float64",
        "empirical_percentile": "float64",
        "lower_tail_p_value": "float64",
        "upper_tail_p_value": "float64",
        "two_sided_p_value": "float64",
        "fdr_q_value": "float64",
        "leave_one_out": "bool",
    },
    "benchmark_cluster_classifications": {
        "matched_control_count": "int64",
        "mean_distance": "float64",
        "mean_distance_control_median": "float64",
        "mean_distance_control_percentile": "float64",
        "distance_sd": "float64",
        "distance_sd_control_median": "float64",
        "distance_sd_control_percentile": "float64",
    },
}


class _StageRecorder:
    """Persist machine-readable stage timing while also logging progress."""

    def __init__(self, *, path: Path) -> None:
        """Initialise an empty stage recorder.

        Args:
            path: Persistent TSV destination within the staging resource.
        """

        self.path = path
        self.rows: list[dict[str, Any]] = []

    @contextmanager
    def record(self, *, stage: str) -> Iterator[dict[str, Any]]:
        """Time one stage and persist its completion or failure.

        Args:
            stage: Stable machine-readable stage name.

        Yields:
            Mutable row whose ``details`` value may be populated by the caller.
        """

        started_at = utc_now_iso()
        started = time.perf_counter()
        row: dict[str, Any] = {
            "stage": stage,
            "status": "RUNNING",
            "started_at_utc": started_at,
            "finished_at_utc": "",
            "elapsed_seconds": "",
            "details": "",
        }
        _LOGGER.info("Stage started: %s", stage)
        try:
            yield row
        except Exception:
            row["status"] = "FAILED"
            raise
        else:
            row["status"] = "PASS"
        finally:
            row["finished_at_utc"] = utc_now_iso()
            row["elapsed_seconds"] = f"{time.perf_counter() - started:.3f}"
            self.rows.append(row)
            write_tsv(path=self.path, fieldnames=_STAGE_FIELDS, records=self.rows)
            _LOGGER.info(
                "Stage finished: %s | status=%s | elapsed_seconds=%s | %s",
                stage,
                row["status"],
                row["elapsed_seconds"],
                row["details"],
            )


def run_pipeline(
    *,
    results_dir: Path,
    output_dir: Path,
    run_id: str,
    work_dir: Path | None,
    alignment_dir: Path | None,
    distance_source: str,
    distance_group_type: str,
    distance_hierarchy_node: str,
    distance_max_groups: int,
    distance_max_members: int,
    parse_gene_trees: bool,
    report_max_statistic_rows: int,
    report_max_groups: int,
    report_max_members: int,
    report_nearest_neighbours: int,
    resume: bool,
    force: bool,
    verbose: bool,
    keep_failed_work: bool = False,
    focus_proteins_path: Path | None = None,
    focus_group_type: str = "HOG",
    focus_hierarchy_node: str = "N0",
    benchmark_proteins_path: Path | None = None,
    benchmark_controls_per_group: int = 3,
    benchmark_bootstrap_resamples: int = 1_000,
) -> dict[str, Any]:
    """Build a complete versioned result resource without mutating its authority.

    Args:
        results_dir: Read-only completed OrthoFinder result directory.
        output_dir: New formal run directory.
        run_id: Stable identifier unique to this OrthoFinder run.
        work_dir: Optional staging root. Slurm jobs should use node-local scratch;
            other runs default beside ``output_dir``.
        alignment_dir: Optional aligned FASTA directory for distance calculation.
        distance_source: ``AUTO``, aligned sequence, resolved tree or disabled.
        distance_group_type: ``AUTO``, ``HOG`` or ``LEGACY_ORTHOGROUP``.
        distance_hierarchy_node: HOG node assigned to aligned FASTA filenames.
        distance_max_groups: Maximum alignment groups; zero means unlimited.
        distance_max_members: Exact/sample distance member limit per group.
        parse_gene_trees: Normalise all available gene-tree nodes and edges.
        report_max_statistic_rows: Browser-safe maximum embedded summary rows.
        report_max_groups: Maximum interactive report networks.
        report_max_members: Maximum rendered protein nodes per network.
        report_nearest_neighbours: Nearest-neighbour edges retained per node.
        resume: Reuse an exactly matching completed formal output.
        force: Supersede an existing non-matching output.
        verbose: Enable debug logging.
        keep_failed_work: Retain partial staging/copy directories after a failure.
        focus_proteins_path: Optional E3/focus authority restricting groups and distances.
        focus_group_type: Exact group collection used for focus selection.
        focus_hierarchy_node: Exact HOG hierarchy node used for focus selection.
        benchmark_proteins_path: Optional housekeeping/R marker authority.
        benchmark_controls_per_group: Unique distance-blind controls per target.
        benchmark_bootstrap_resamples: Contrast confidence-interval iterations.

    Returns:
        Completed run manifest.

    Raises:
        InputValidationError: If inputs or named controls are invalid.
        PublicationError: If verified publication cannot be completed.
    """

    _validate_controls(
        run_id=run_id,
        distance_source=distance_source,
        distance_group_type=distance_group_type,
        distance_hierarchy_node=distance_hierarchy_node,
        distance_max_groups=distance_max_groups,
        distance_max_members=distance_max_members,
        report_max_statistic_rows=report_max_statistic_rows,
        report_max_groups=report_max_groups,
        report_max_members=report_max_members,
        report_nearest_neighbours=report_nearest_neighbours,
        resume=resume,
        force=force,
        focus_enabled=focus_proteins_path is not None,
        focus_group_type=focus_group_type,
        focus_hierarchy_node=focus_hierarchy_node,
        benchmark_enabled=benchmark_proteins_path is not None,
        benchmark_controls_per_group=benchmark_controls_per_group,
        benchmark_bootstrap_resamples=benchmark_bootstrap_resamples,
    )
    focus_path = (
        Path(focus_proteins_path).expanduser().resolve()
        if focus_proteins_path is not None
        else None
    )
    focus_authority = (
        read_focus_proteins(path=focus_path) if focus_path is not None else None
    )
    benchmark_path = (
        Path(benchmark_proteins_path).expanduser().resolve()
        if benchmark_proteins_path is not None
        else None
    )
    benchmark_authority = (
        read_benchmark_authority(path=benchmark_path)
        if benchmark_path is not None
        else None
    )
    layout = discover_layout(results_dir=results_dir)
    output = validate_persistent_path(path=output_dir, role="output_dir")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging_root = Path(
        work_dir if work_dir is not None else output.parent / ".orthofinder_results_work"
    ).expanduser().resolve()
    if staging_root == output:
        raise InputValidationError("work_dir must not be the formal output directory.")
    staging_root.mkdir(parents=True, exist_ok=True)
    publication_method = (
        "ATOMIC_RENAME"
        if _same_filesystem(first=staging_root, second=output.parent)
        else "VERIFIED_COPY_THEN_ATOMIC_RENAME"
    )
    resolved_alignment_dir = _resolve_pipeline_alignment_dir(
        requested=alignment_dir,
        discovered=layout.alignments_dir,
        distance_source=distance_source,
    )
    _LOGGER.info("Source inventory started: %s", layout.results_dir)
    inventory_started = time.perf_counter()
    source_inventory = _build_source_inventory(
        layout=layout,
        alignment_dir=resolved_alignment_dir,
        focus_proteins_path=focus_path,
        benchmark_proteins_path=benchmark_path,
    )
    _LOGGER.info(
        "Source inventory finished: files=%s, elapsed_seconds=%.3f",
        f"{len(source_inventory):,}",
        time.perf_counter() - inventory_started,
    )
    input_digest = _inventory_digest(records=source_inventory)
    reusable = _resolve_existing_output(
        output_dir=output,
        run_id=run_id,
        input_digest=input_digest,
        resume=resume,
        force=force,
    )
    if reusable is not None:
        return reusable

    staging = staging_root / f"{output.name}.staging.{uuid.uuid4().hex}"
    staging.mkdir(parents=False, exist_ok=False)
    configure_logging(log_path=staging / "logs" / "run.log", verbose=verbose)
    _LOGGER.info("Starting orthofinder-results %s for run %s", __version__, run_id)
    _LOGGER.info("Read-only authority: %s", layout.results_dir)
    _LOGGER.info("Formal output: %s", output)
    _LOGGER.info("Staging root: %s", staging_root)
    _LOGGER.info("Publication method: %s", publication_method)
    started_at = utc_now_iso()
    try:
        manifest = _build_resource(
            staging=staging,
            layout=layout,
            run_id=run_id,
            source_inventory=source_inventory,
            input_digest=input_digest,
            alignment_dir=resolved_alignment_dir,
            distance_source=distance_source,
            distance_group_type=distance_group_type,
            distance_hierarchy_node=distance_hierarchy_node,
            distance_max_groups=distance_max_groups,
            distance_max_members=distance_max_members,
            parse_gene_trees=parse_gene_trees,
            report_max_statistic_rows=report_max_statistic_rows,
            report_max_groups=report_max_groups,
            report_max_members=report_max_members,
            report_nearest_neighbours=report_nearest_neighbours,
            started_at=started_at,
            staging_root=staging_root,
            publication_method=publication_method,
            focus_authority=focus_authority,
            focus_proteins_path=focus_path,
            focus_group_type=focus_group_type,
            focus_hierarchy_node=focus_hierarchy_node,
            benchmark_authority=benchmark_authority,
            benchmark_proteins_path=benchmark_path,
            benchmark_controls_per_group=benchmark_controls_per_group,
            benchmark_bootstrap_resamples=benchmark_bootstrap_resamples,
        )
        # Close the staging file handler before checksums are verified or files
        # cross filesystems. Subsequent CLI messages remain console-only.
        configure_logging(verbose=verbose)
        _publish_completed_resource(
            staging=staging,
            output=output,
            keep_failed_work=keep_failed_work,
        )
        _LOGGER.info("Published completed resource: %s", output)
        return manifest
    except Exception:
        _LOGGER.exception("Run failed before formal publication.")
        configure_logging(verbose=verbose)
        if staging.exists():
            if keep_failed_work:
                failed = staging.with_name(staging.name.replace(".staging.", ".failed."))
                os.replace(staging, failed)
                _LOGGER.error("Retained diagnostic staging directory: %s", failed)
            else:
                shutil.rmtree(staging)
                _LOGGER.info("Removed partial staging directory: %s", staging)
        raise


def _same_filesystem(*, first: Path, second: Path) -> bool:
    """Return whether two existing paths use the same filesystem.

    Args:
        first: First existing path.
        second: Second existing path.

    Returns:
        ``True`` when both device identifiers match.
    """

    return first.stat().st_dev == second.stat().st_dev


def _publish_completed_resource(
    *, staging: Path, output: Path, keep_failed_work: bool
) -> None:
    """Publish a validated staging tree to persistent storage.

    Same-filesystem runs use a direct atomic rename. Cross-filesystem runs copy
    into a hidden directory beside the formal output, verify every manifested
    file, and then atomically rename that verified copy.

    Args:
        staging: Completed resource staging directory.
        output: Formal persistent output directory.
        keep_failed_work: Retain an incomplete cross-filesystem copy for diagnosis.

    Raises:
        PublicationError: If copying, validation or final publication fails.
    """

    publication_started = time.perf_counter()
    if _same_filesystem(first=staging, second=output.parent):
        _LOGGER.info("Atomic same-filesystem publication started: %s", output)
        os.replace(staging, output)
        _LOGGER.info(
            "Atomic publication finished in %.3f seconds.",
            time.perf_counter() - publication_started,
        )
        return

    token = uuid.uuid4().hex
    incoming = output.parent / f".{output.name}.incoming.{token}"
    try:
        _LOGGER.info("Cross-filesystem copy started: %s", incoming)
        shutil.copytree(staging, incoming, copy_function=shutil.copy2)
        _LOGGER.info(
            "Cross-filesystem copy finished in %.3f seconds; checksum validation started.",
            time.perf_counter() - publication_started,
        )
        _validate_published_copy(source=staging, destination=incoming)
        os.replace(incoming, output)
        _LOGGER.info(
            "Verified publication finished in %.3f seconds.",
            time.perf_counter() - publication_started,
        )
    except Exception as error:
        if incoming.exists():
            if keep_failed_work:
                failed_copy = output.parent / f".{output.name}.copy_failed.{token}"
                os.replace(incoming, failed_copy)
                _LOGGER.error("Retained incomplete persistent copy for diagnosis: %s", failed_copy)
            else:
                shutil.rmtree(incoming)
                _LOGGER.info("Removed incomplete persistent copy: %s", incoming)
        if isinstance(error, PublicationError):
            raise
        raise PublicationError(
            f"Could not publish verified resource to {output}: {error}"
        ) from error

    try:
        shutil.rmtree(staging)
    except OSError as error:
        _LOGGER.warning(
            "Published output but could not remove scratch staging %s: %s", staging, error
        )


def _validate_published_copy(*, source: Path, destination: Path) -> None:
    """Verify a cross-filesystem copy against its completed manifest.

    Args:
        source: Completed scratch resource.
        destination: Persistent incoming copy.

    Raises:
        PublicationError: If manifests, file sets, sizes or checksums differ.
    """

    source_manifest = source / "run_manifest.json"
    destination_manifest = destination / "run_manifest.json"
    if not source_manifest.is_file() or not destination_manifest.is_file():
        raise PublicationError("Completed resource copy is missing run_manifest.json.")
    if sha256_file(path=source_manifest) != sha256_file(path=destination_manifest):
        raise PublicationError(
            "Copied run_manifest.json checksum does not match scratch authority."
        )

    try:
        manifest = json.loads(destination_manifest.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as error:
        raise PublicationError(f"Copied run manifest is unreadable: {error}") from error
    if manifest.get("status") != "complete":
        raise PublicationError("Copied run manifest is not marked complete.")

    expected_paths: set[Path] = set()
    for record in manifest.get("outputs", []):
        relative_path = Path(str(record.get("path", "")))
        if (
            not relative_path.parts
            or relative_path.is_absolute()
            or ".." in relative_path.parts
        ):
            raise PublicationError(f"Unsafe output path in run manifest: {relative_path}")
        expected_paths.add(relative_path)
        copied_path = destination / relative_path
        if not copied_path.is_file():
            raise PublicationError(f"Copied resource is missing manifested file: {relative_path}")
        if copied_path.stat().st_size != int(record["size_bytes"]):
            raise PublicationError(f"Copied file size differs for: {relative_path}")
        if sha256_file(path=copied_path) != record["sha256"]:
            raise PublicationError(f"Copied file checksum differs for: {relative_path}")

    actual_paths = {
        path.relative_to(destination)
        for path in destination.rglob("*")
        if path.is_file() and path.name != "run_manifest.json"
    }
    if actual_paths != expected_paths:
        missing = sorted(str(path) for path in expected_paths - actual_paths)
        unexpected = sorted(str(path) for path in actual_paths - expected_paths)
        raise PublicationError(
            "Copied resource file set differs from its manifest; "
            f"missing={missing}, unexpected={unexpected}."
        )


def inspect_results(*, results_dir: Path) -> dict[str, Any]:
    """Return a read-only version and capability inspection.

    Args:
        results_dir: Completed OrthoFinder results directory.

    Returns:
        JSON-safe layout record.
    """

    return discover_layout(results_dir=results_dir).to_record()


def regenerate_report(
    *,
    resource_dir: Path,
    output_path: Path,
    work_dir: Path | None,
    report_max_statistic_rows: int,
    report_max_groups: int,
    report_max_members: int,
    report_nearest_neighbours: int,
    force: bool,
) -> dict[str, Any]:
    """Regenerate only the HTML from a completed analytical resource.

    The completed resource is read-only. The new report is written to a separate
    persistent path, avoiding membership parsing, distance calculation, Parquet
    conversion and DuckDB rebuilding.

    Args:
        resource_dir: Completed orthofinder-results resource directory.
        output_path: New standalone HTML output.
        work_dir: Optional temporary root for node-local copies of report inputs.
        report_max_statistic_rows: Maximum embedded group summaries.
        report_max_groups: Maximum interactive networks.
        report_max_members: Maximum members per network.
        report_nearest_neighbours: Nearest-neighbour edges retained per node.
        force: Permit atomic replacement of an existing standalone report.

    Returns:
        Output path, size and SHA-256 record.

    Raises:
        InputValidationError: If the resource or controls are invalid.
        PublicationError: If the resource is incomplete or the output exists.
    """

    _validate_report_controls(
        report_max_statistic_rows=report_max_statistic_rows,
        report_max_groups=report_max_groups,
        report_max_members=report_max_members,
        report_nearest_neighbours=report_nearest_neighbours,
    )
    resource = validate_persistent_path(path=resource_dir, role="resource_dir")
    destination = validate_persistent_path(path=output_path, role="report_output")
    if not resource.is_dir():
        raise InputValidationError(f"Resource directory does not exist: {resource}")
    if resource == destination or resource in destination.parents:
        raise InputValidationError(
            "Standalone report output must be outside the immutable completed resource."
        )
    if work_dir is not None:
        resolved_work = Path(work_dir).expanduser().resolve()
        if resource == resolved_work or resource in resolved_work.parents:
            raise InputValidationError(
                "Report work directory must be outside the immutable completed resource."
            )
    manifest_path = resource / "run_manifest.json"
    if not manifest_path.is_file():
        raise InputValidationError(
            f"Completed resource lacks run_manifest.json: {resource}"
        )
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as error:
        raise InputValidationError(f"Resource manifest is unreadable: {error}") from error
    if manifest.get("status") != "complete":
        raise PublicationError("Report regeneration requires a complete resource manifest.")
    if destination.exists() and not force:
        raise PublicationError(
            f"Report output already exists: {destination}. Use --force to replace it atomically."
        )
    tables = resource / "tables"
    required_relations = (
        "distance_statistics",
        "group_species_statistics",
        "group_statistics",
        "hog_memberships",
        "legacy_orthogroup_memberships",
        "pairwise_distances",
        "sequences",
        "tree_edges",
        "tree_inventory",
        "tree_nodes",
    )
    missing = [
        relation
        for relation in required_relations
        if not _table_path(tables_dir=tables, relation=relation).is_file()
    ]
    if missing:
        raise InputValidationError(
            f"Resource lacks required compressed TSV relations: {', '.join(missing)}"
        )

    _LOGGER.info("Starting report-only regeneration with orthofinder-results %s.", __version__)
    _LOGGER.info("Read-only completed resource: %s", resource)
    _LOGGER.info("Standalone report output: %s", destination)
    with _report_table_source(
        source_tables=tables,
        work_dir=work_dir,
        relations=required_relations,
    ) as report_tables:
        return _build_standalone_report(
            tables=report_tables,
            destination=destination,
            manifest=manifest,
            report_max_statistic_rows=report_max_statistic_rows,
            report_max_groups=report_max_groups,
            report_max_members=report_max_members,
            report_nearest_neighbours=report_nearest_neighbours,
        )


@contextmanager
def _report_table_source(
    *, source_tables: Path, work_dir: Path | None, relations: Sequence[str]
) -> Iterator[Path]:
    """Yield source tables or node-local copies and always clean temporary data.

    Args:
        source_tables: Persistent completed-resource table directory.
        work_dir: Optional scheduler-provided temporary root.
        relations: Compressed relations required for report construction.

    Yields:
        Directory containing report input tables.
    """

    if work_dir is None:
        _LOGGER.info("Report input mode: direct persistent reads.")
        yield source_tables
        return
    root = Path(work_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    staging = root / f"report_inputs_{uuid.uuid4().hex}"
    staging.mkdir(parents=False, exist_ok=False)
    started = time.perf_counter()
    copied_bytes = 0
    try:
        _LOGGER.info("Node-local report input staging started: %s", staging)
        for relation in relations:
            source = _table_path(tables_dir=source_tables, relation=relation)
            destination = _table_path(tables_dir=staging, relation=relation)
            shutil.copy2(source, destination)
            copied_bytes += destination.stat().st_size
            _LOGGER.info(
                "Node-local report input copied: relation=%s, size_bytes=%s",
                relation,
                f"{destination.stat().st_size:,}",
            )
        _LOGGER.info(
            "Node-local report input staging finished: files=%s, bytes=%s, "
            "elapsed_seconds=%.3f",
            len(relations),
            f"{copied_bytes:,}",
            time.perf_counter() - started,
        )
        yield staging
    finally:
        if staging.exists():
            shutil.rmtree(staging)
            _LOGGER.info("Removed node-local report inputs: %s", staging)


def _build_standalone_report(
    *,
    tables: Path,
    destination: Path,
    manifest: Mapping[str, Any],
    report_max_statistic_rows: int,
    report_max_groups: int,
    report_max_members: int,
    report_nearest_neighbours: int,
) -> dict[str, Any]:
    """Build a standalone report from validated compressed analytical tables."""

    started = time.perf_counter()
    distance_summaries = list(
        read_tsv(path=_table_path(tables_dir=tables, relation="distance_statistics"))
    )
    group_rows, overview_statistics = _load_report_group_statistics_and_aggregates(
        path=_table_path(tables_dir=tables, relation="group_statistics"),
        maximum=report_max_statistic_rows,
    )
    _LOGGER.info("Loaded %s bounded group summaries.", f"{len(group_rows):,}")
    network_rows, memberships, distances = _load_report_network_data(
        tables_dir=tables,
        group_statistics_path=_table_path(tables_dir=tables, relation="group_statistics"),
        distance_summaries=distance_summaries,
        maximum_groups=report_max_groups,
        maximum_members=report_max_members,
    )
    group_species = _load_report_group_species_data(
        tables_dir=tables,
        memberships=memberships,
    )
    tree_nodes, tree_edges, sequence_identifiers = _load_report_tree_data(
        tables_dir=tables,
        group_statistics=network_rows,
        memberships=memberships,
        distance_summaries=distance_summaries,
    )
    counts = dict(manifest.get("counts", {}))
    run_metadata = {
        key: manifest.get(key)
        for key in (
            "run_id",
            "orthofinder_version",
            "adapter_name",
            "primary_group_authority",
            "schema_version",
            "distance_source_requested",
            "capabilities",
            "publication",
        )
    }
    run_metadata["resource_package_version"] = manifest.get("package_version", "")
    run_metadata["package_version"] = __version__
    run_metadata["counts"] = counts
    build_interactive_report(
        output_path=destination,
        run_metadata=run_metadata,
        group_statistics=group_rows,
        network_group_statistics=network_rows,
        memberships=memberships,
        group_species_statistics=group_species,
        distances=distances,
        distance_statistics=distance_summaries,
        total_group_statistic_count=int(counts.get("group_count", 0)),
        total_membership_count=(
            int(counts.get("hog_membership_count", 0))
            + int(counts.get("legacy_orthogroup_membership_count", 0))
        ),
        max_network_groups=report_max_groups,
        max_network_members=report_max_members,
        nearest_neighbours=report_nearest_neighbours,
        overview_statistics=overview_statistics,
        tree_nodes=tree_nodes,
        tree_edges=tree_edges,
        sequence_identifiers=sequence_identifiers,
    )
    elapsed = time.perf_counter() - started
    record = file_record(path=destination)
    _LOGGER.info(
        "Finished report-only regeneration in %.2f seconds: networks=%s, "
        "memberships=%s, distance_rows=%s, size_bytes=%s.",
        elapsed,
        f"{len(network_rows):,}",
        f"{len(memberships):,}",
        f"{len(distances):,}",
        f"{record['size_bytes']:,}",
    )
    return record


def _build_resource(
    *,
    staging: Path,
    layout: ResultLayout,
    run_id: str,
    source_inventory: Sequence[Mapping[str, Any]],
    input_digest: str,
    alignment_dir: Path | None,
    distance_source: str,
    distance_group_type: str,
    distance_hierarchy_node: str,
    distance_max_groups: int,
    distance_max_members: int,
    parse_gene_trees: bool,
    report_max_statistic_rows: int,
    report_max_groups: int,
    report_max_members: int,
    report_nearest_neighbours: int,
    started_at: str,
    staging_root: Path,
    publication_method: str,
    focus_authority: FocusProteinAuthority | None,
    focus_proteins_path: Path | None,
    focus_group_type: str,
    focus_hierarchy_node: str,
    benchmark_authority: BenchmarkAuthority | None,
    benchmark_proteins_path: Path | None,
    benchmark_controls_per_group: int,
    benchmark_bootstrap_resamples: int,
) -> dict[str, Any]:
    """Populate one staging directory and return its complete manifest.

    Args:
        staging: Unique incomplete resource directory.
        layout: Validated completed OrthoFinder layout.
        run_id: Immutable output run identifier.
        source_inventory: Checksum-bound input inventory.
        input_digest: Digest of the complete input inventory.
        alignment_dir: Optional alignment distance authority.
        distance_source: Requested distance authority.
        distance_group_type: Requested distance group collection.
        distance_hierarchy_node: Exact distance hierarchy node.
        distance_max_groups: Maximum ordinary-run distance groups; zero is unlimited.
        distance_max_members: Per-group exact or deterministic member bound.
        parse_gene_trees: Whether selected gene trees should be normalised.
        report_max_statistic_rows: Browser-safe group-summary bound.
        report_max_groups: Browser-safe network-group bound.
        report_max_members: Browser-safe network-member bound.
        report_nearest_neighbours: Network neighbour count.
        started_at: Run start time in UTC.
        staging_root: Parent temporary directory used for publication metadata.
        publication_method: Atomic publication strategy.
        focus_authority: Optional validated E3/focus authority.
        focus_proteins_path: Physical focus input copied into provenance.
        focus_group_type: Exact focus group collection.
        focus_hierarchy_node: Exact focus HOG hierarchy node.
        benchmark_authority: Optional reviewed housekeeping/R marker authority.
        benchmark_proteins_path: Physical benchmark authority copied to provenance.
        benchmark_controls_per_group: Unique non-focus controls per target.
        benchmark_bootstrap_resamples: Deterministic contrast bootstrap iterations.

    Returns:
        Complete resource manifest.
    """

    tables = staging / "tables"
    provenance = staging / "provenance"
    qc_dir = staging / "qc"
    report_dir = staging / "report"
    database_dir = staging / "duckdb"
    for directory in (tables, provenance, qc_dir, report_dir, database_dir):
        directory.mkdir(parents=True, exist_ok=True)
    stages = _StageRecorder(path=staging / "logs" / "stage_metrics.tsv")
    with stages.record(stage="provenance") as stage:
        atomic_write_json(path=provenance / "resolved_layout.json", record=layout.to_record())
        write_tsv(
            path=provenance / "input_inventory.tsv",
            fieldnames=("role", "path", "size_bytes", "sha256"),
            records=source_inventory,
        )
        if focus_proteins_path is not None:
            suffix = ".tsv.gz" if focus_proteins_path.name.endswith(".tsv.gz") else ".tsv"
            shutil.copy2(
                focus_proteins_path,
                provenance / f"focus_protein_authority{suffix}",
            )
        if benchmark_proteins_path is not None:
            shutil.copy2(
                benchmark_proteins_path,
                provenance / "dispersion_benchmark_authority.tsv",
            )
        stage["details"] = f"input_files={len(source_inventory)}"

    with stages.record(stage="memberships_and_group_statistics") as stage:
        (
            membership_counts,
            group_count,
            group_species_statistic_count,
            species_from_groups,
        ) = _publish_memberships(
            tables_dir=tables,
            layout=layout,
            run_id=run_id,
        )
        stage["details"] = (
            f"memberships={sum(membership_counts.values())};groups={group_count};"
            f"group_species_rows={group_species_statistic_count}"
        )
    with stages.record(stage="species_and_sequence_identifiers") as stage:
        species_count, sequence_count = _publish_identifiers(
            tables_dir=tables,
            layout=layout,
            run_id=run_id,
            species_from_groups=species_from_groups,
        )
        stage["details"] = f"species={species_count};sequences={sequence_count}"
    focus_selection: FocusSelection | None = None
    if focus_authority is not None:
        with stages.record(stage="e3_focus_selection") as stage:
            focus_selection = select_focus_groups(
                tables_dir=tables,
                run_id=run_id,
                authority=focus_authority,
                group_type=focus_group_type,
                hierarchy_node=focus_hierarchy_node,
            )
            if not focus_selection.group_statistics:
                raise InputValidationError(
                    "No enabled focus protein matched the selected group collection."
                )
            publish_focus_selection(
                tables_dir=tables,
                selection=focus_selection,
            )
            stage["details"] = (
                f"authority_records={len(focus_authority.records)};"
                f"matched_seeds={focus_selection.matched_seed_count};"
                f"groups={len(focus_selection.group_statistics)};"
                f"matches={len(focus_selection.match_rows)}"
            )
    benchmark_plan: BenchmarkPlan | None = None
    if benchmark_authority is not None:
        if focus_selection is None:
            raise InputValidationError(
                "Dispersion benchmarking requires an E3/focus authority."
            )
        with stages.record(stage="dispersion_benchmark_selection") as stage:
            marker_selection = select_focus_groups(
                tables_dir=tables,
                run_id=run_id,
                authority=benchmark_authority.to_focus_authority(),
                group_type=focus_group_type,
                hierarchy_node=focus_hierarchy_node,
            )
            if not marker_selection.group_statistics:
                raise InputValidationError(
                    "No enabled dispersion benchmark marker matched the group collection."
                )
            benchmark_plan = build_benchmark_plan(
                tables_dir=tables,
                run_id=run_id,
                e3_selection=focus_selection,
                marker_selection=marker_selection,
                marker_authority=benchmark_authority,
                controls_per_group=benchmark_controls_per_group,
            )
            publish_benchmark_inputs(tables_dir=tables, plan=benchmark_plan)
            stage["details"] = (
                f"marker_records={len(benchmark_authority.records)};"
                f"matched_markers={marker_selection.matched_seed_count};"
                f"target_groups={len(benchmark_plan.target_group_keys)};"
                f"matched_controls={len(benchmark_plan.control_group_keys)}"
            )
    with stages.record(stage="tree_inventory_and_normalisation") as stage:
        (
            tree_inventory,
            tree_payload_count,
            tree_node_count,
            tree_edge_count,
        ) = _publish_trees(
            tables_dir=tables,
            layout=layout,
            run_id=run_id,
            parse_gene_trees=parse_gene_trees,
            selected_gene_tree_ids=(
                benchmark_plan.selected_tree_ids
                if benchmark_plan is not None
                else focus_selection.selected_tree_ids
                if focus_selection is not None
                else None
            ),
        )
        stage["details"] = (
            f"tree_files={len(tree_inventory)};portable_trees={tree_payload_count};"
            f"nodes={tree_node_count};edges={tree_edge_count}"
        )
    with stages.record(stage="pairwise_distances") as stage:
        distance_count, distance_summaries = _publish_distances(
            tables_dir=tables,
            layout=layout,
            run_id=run_id,
            alignment_dir=alignment_dir,
            distance_source=distance_source,
            distance_group_type=distance_group_type,
            distance_hierarchy_node=distance_hierarchy_node,
            distance_max_groups=distance_max_groups,
            distance_max_members=distance_max_members,
            selected_group_keys=(
                benchmark_plan.selected_group_keys
                if benchmark_plan is not None
                else focus_selection.group_keys
                if focus_selection is not None
                else None
            ),
            required_members_by_group=(
                benchmark_plan.required_members_by_group
                if benchmark_plan is not None
                else focus_selection.required_members_by_group
                if focus_selection is not None
                else None
            ),
            persist_pairwise_group_keys=(
                benchmark_plan.target_group_keys
                if benchmark_plan is not None
                else None
            ),
        )
        stage["details"] = (
            f"groups={len(distance_summaries)};pairs={distance_count};requested={distance_source}"
        )
    focus_result_rows: tuple[dict[str, Any], ...] = ()
    if focus_selection is not None:
        with stages.record(stage="e3_cluster_results") as stage:
            focus_result_rows = publish_focus_cluster_results(
                tables_dir=tables,
                selection=focus_selection,
                distance_summaries=distance_summaries,
            )
            unavailable = sum(
                row["computation_status"] == "UNAVAILABLE"
                for row in focus_result_rows
            )
            stage["details"] = (
                f"groups={len(focus_result_rows)};"
                f"unavailable_distance_groups={unavailable}"
            )
    benchmark_counts: dict[str, int] = {}
    if benchmark_plan is not None:
        with stages.record(stage="dispersion_benchmark_statistics") as stage:
            benchmark_counts = publish_benchmark_results(
                tables_dir=tables,
                plan=benchmark_plan,
                distance_summaries=distance_summaries,
                bootstrap_resamples=benchmark_bootstrap_resamples,
            )
            stage["details"] = ";".join(
                f"{key}={value}" for key, value in sorted(benchmark_counts.items())
            )

    with stages.record(stage="parquet_publication") as stage:
        parquet_tables = _publish_parquet(tables_dir=tables)
        stage["details"] = f"relations={len(parquet_tables)}"
    with stages.record(stage="duckdb_publication") as stage:
        database_path = database_dir / "orthofinder_results.duckdb"
        create_duckdb(
            database_path=database_path,
            parquet_tables=parquet_tables,
        )
        stage["details"] = f"size_bytes={database_path.stat().st_size}"
    group_rows, overview_statistics = _load_report_group_statistics_and_aggregates(
        path=_table_path(tables_dir=tables, relation="group_statistics"),
        maximum=report_max_statistic_rows,
    )
    network_group_rows, report_memberships, report_distances = _load_report_network_data(
        tables_dir=tables,
        group_statistics_path=_table_path(
            tables_dir=tables,
            relation="group_statistics",
        ),
        distance_summaries=distance_summaries,
        maximum_groups=report_max_groups,
        maximum_members=report_max_members,
    )
    report_group_species = _load_report_group_species_data(
        tables_dir=tables,
        memberships=report_memberships,
    )
    report_tree_nodes, report_tree_edges, report_sequence_identifiers = (
        _load_report_tree_data(
            tables_dir=tables,
            group_statistics=network_group_rows,
            memberships=report_memberships,
            distance_summaries=distance_summaries,
        )
    )
    run_metadata = {
        "run_id": run_id,
        "orthofinder_version": layout.orthofinder_version,
        "adapter_name": layout.adapter_name,
        "primary_group_authority": layout.primary_group_authority,
        "package_version": __version__,
        "schema_version": __schema_version__,
        "distance_source_requested": distance_source,
        "capabilities": layout.capabilities.to_record(),
        "publication": {
            "method": publication_method,
            "staging_root": str(staging_root),
            "copy_verified": publication_method == "VERIFIED_COPY_THEN_ATOMIC_RENAME",
        },
        "counts": {
            "group_count": group_count,
            "species_count": species_count,
            "tree_inventory_count": len(tree_inventory),
            "tree_payload_count": tree_payload_count,
            "tree_node_count": tree_node_count,
            "distance_group_count": len(distance_summaries),
            "distance_pair_count": distance_count,
            "focus_seed_count": (
                len(focus_authority.records) if focus_authority is not None else 0
            ),
            "matched_focus_seed_count": (
                focus_selection.matched_seed_count
                if focus_selection is not None
                else 0
            ),
            "focus_group_count": len(focus_result_rows),
            **benchmark_counts,
        },
    }
    if focus_selection is not None:
        run_metadata["focus_analysis"] = {
            "authority_name": focus_selection.authority.source_name,
            "authority_sha256": focus_selection.authority.sha256,
            "group_type": focus_selection.group_type,
            "hierarchy_node": focus_selection.hierarchy_node,
            "cluster_results": "tables/e3_cluster_results.tsv.gz",
            "seed_matches": "tables/e3_seed_matches.tsv.gz",
            "seed_audit": "tables/e3_seed_catalogue_audit.tsv.gz",
            "pairwise_distances": "tables/pairwise_distances.tsv.gz",
        }
    if benchmark_plan is not None:
        run_metadata["dispersion_benchmark"] = {
            "marker_authority_name": benchmark_plan.marker_authority.source_name,
            "marker_authority_sha256": benchmark_plan.marker_authority.sha256,
            "controls_per_target": benchmark_controls_per_group,
            "bootstrap_resamples": benchmark_bootstrap_resamples,
            "statistical_unit": "ORTHOFINDER_CLUSTER",
            "matching_variables": (
                "member_count;species_count;mean_copies_per_species;"
                "single_copy_species_fraction"
            ),
            "cluster_results": "tables/benchmark_cluster_results.tsv.gz",
            "background_statistics": (
                "tables/benchmark_background_statistics.tsv.gz"
            ),
            "profile_contrasts": "tables/benchmark_contrasts.tsv.gz",
            "individual_comparisons": (
                "tables/benchmark_individual_comparisons.tsv.gz"
            ),
            "classifications": (
                "tables/benchmark_cluster_classifications.tsv.gz"
            ),
            "matched_controls": "tables/benchmark_matched_controls.tsv.gz",
        }
    report_path = report_dir / "orthofinder_results_summary.html"
    with stages.record(stage="offline_html_report") as stage:
        build_interactive_report(
            output_path=report_path,
            run_metadata=run_metadata,
            group_statistics=group_rows,
            network_group_statistics=network_group_rows,
            memberships=report_memberships,
            group_species_statistics=report_group_species,
            distances=report_distances,
            distance_statistics=distance_summaries,
            total_group_statistic_count=group_count,
            total_membership_count=sum(membership_counts.values()),
            max_network_groups=report_max_groups,
            max_network_members=report_max_members,
            nearest_neighbours=report_nearest_neighbours,
            overview_statistics=overview_statistics,
            tree_nodes=report_tree_nodes,
            tree_edges=report_tree_edges,
            sequence_identifiers=report_sequence_identifiers,
        )
        stage["details"] = (
            f"group_rows={len(group_rows)};network_groups={len(network_group_rows)};"
            f"members={len(report_memberships)};distance_rows={len(report_distances)};"
            f"size_bytes={report_path.stat().st_size}"
        )
    report_text = report_path.read_text(encoding="utf-8")
    offline_report = not re.search(r"(?:src|href)=[\"']https?://", report_text, re.IGNORECASE)
    qc_rows = _qc_rows(
        layout=layout,
        membership_counts=membership_counts,
        group_count=group_count,
        species_count=species_count,
        sequence_count=sequence_count,
        tree_inventory_count=len(tree_inventory),
        tree_payload_count=tree_payload_count,
        tree_node_count=tree_node_count,
        distance_count=distance_count,
        offline_report=offline_report,
        focus_selection=focus_selection,
        focus_result_rows=focus_result_rows,
        distance_summaries=distance_summaries,
        benchmark_plan=benchmark_plan,
        benchmark_counts=benchmark_counts,
    )
    with stages.record(stage="quality_control") as stage:
        write_tsv(
            path=qc_dir / "validation_checks.tsv",
            fieldnames=("check_name", "status", "observed_value", "expected_value", "details"),
            records=qc_rows,
        )
        failed_checks = [row["check_name"] for row in qc_rows if row["status"] == "FAIL"]
        stage["details"] = (
            f"checks={len(qc_rows)};failed={';'.join(failed_checks) if failed_checks else 'none'}"
        )
        if failed_checks:
            raise PublicationError("One or more required validation checks failed.")

    output_inventory = [
        file_record(path=path, relative_to=staging)
        for path in sorted(staging.rglob("*"))
        if path.is_file() and path.name != "run_manifest.json"
    ]
    manifest: dict[str, Any] = {
        **run_metadata,
        "status": "complete",
        "started_at_utc": started_at,
        "finished_at_utc": utc_now_iso(),
        "input_digest": input_digest,
        "source_results_dir": str(layout.results_dir),
        "counts": {
            **membership_counts,
            "group_count": group_count,
            "group_species_statistic_count": group_species_statistic_count,
            "species_count": species_count,
            "sequence_count": sequence_count,
            "tree_file_count": len(tree_inventory),
            "tree_payload_count": tree_payload_count,
            "tree_node_count": tree_node_count,
            "tree_edge_count": tree_edge_count,
            "distance_pair_count": distance_count,
            "distance_group_count": len(distance_summaries),
            "focus_seed_count": (
                len(focus_authority.records) if focus_authority is not None else 0
            ),
            "matched_focus_seed_count": (
                focus_selection.matched_seed_count
                if focus_selection is not None
                else 0
            ),
            "focus_group_count": len(focus_result_rows),
            **benchmark_counts,
        },
        "report_limits": {
            "maximum_statistic_rows": report_max_statistic_rows,
            "maximum_network_groups": report_max_groups,
            "maximum_network_members": report_max_members,
            "nearest_neighbours": report_nearest_neighbours,
        },
        "outputs": output_inventory,
        "scientific_limitations": [
            "Group membership does not by itself prove orthology or biological function.",
            "Distance method and exact-versus-sampled status must accompany every interpretation.",
            (
                "Interactive networks are bounded exploratory views; analytical tables "
                "are authoritative."
            ),
            (
                "OrthoFinder group identifiers are scoped to this run and must not be "
                "joined across runs by label alone."
            ),
            (
                "Housekeeping and R/NLR panels are empirical benchmark candidates, "
                "not assumed compact or dispersed truths."
            ),
            (
                "Dispersion tests use clusters as independent units and matched-control "
                "residuals; individual protein-pair rows are not statistical replicates."
            ),
            (
                "Matched non-focus controls reduce measured structural confounding but "
                "cannot remove unmeasured biological or annotation confounding."
            ),
        ],
    }
    atomic_write_json(path=staging / "run_manifest.json", record=manifest)
    return manifest


def _publish_memberships(
    *, tables_dir: Path, layout: ResultLayout, run_id: str
) -> tuple[dict[str, int], int, int, set[str]]:
    """Publish long-form membership tables and streaming group statistics."""

    statistics_path = _table_path(tables_dir=tables_dir, relation="group_statistics")
    counts = {"legacy_orthogroup_membership_count": 0, "hog_membership_count": 0}
    group_count = 0
    group_species_count = 0
    species: set[str] = set()
    species_statistics_path = _table_path(
        tables_dir=tables_dir,
        relation="group_species_statistics",
    )
    with (
        open_text(path=statistics_path, mode="w") as stats_handle,
        open_text(path=species_statistics_path, mode="w") as species_stats_handle,
    ):
        stats_writer = csv.DictWriter(
            stats_handle,
            fieldnames=GROUP_STATISTIC_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        stats_writer.writeheader()
        species_stats_writer = csv.DictWriter(
            species_stats_handle,
            fieldnames=GROUP_SPECIES_STATISTIC_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        species_stats_writer.writeheader()
        legacy_sources = (
            []
            if layout.orthogroups_path is None
            else [(layout.orthogroups_path, "LEGACY_ORTHOGROUP", "")]
        )
        count, groups, group_species_rows = _write_membership_authority(
            path=_table_path(
                tables_dir=tables_dir,
                relation="legacy_orthogroup_memberships",
            ),
            sources=legacy_sources,
            run_id=run_id,
            statistics_writer=stats_writer,
            species_statistics_writer=species_stats_writer,
            species=species,
        )
        counts["legacy_orthogroup_membership_count"] = count
        group_count += groups
        group_species_count += group_species_rows
        hog_sources = [(path, "HOG", path.stem) for path in layout.hog_paths]
        count, groups, group_species_rows = _write_membership_authority(
            path=_table_path(tables_dir=tables_dir, relation="hog_memberships"),
            sources=hog_sources,
            run_id=run_id,
            statistics_writer=stats_writer,
            species_statistics_writer=species_stats_writer,
            species=species,
        )
        counts["hog_membership_count"] = count
        group_count += groups
        group_species_count += group_species_rows
    _LOGGER.info(
        "Published %s memberships across %s run-scoped groups.",
        f"{sum(counts.values()):,}",
        f"{group_count:,}",
    )
    return counts, group_count, group_species_count, species


def _write_membership_authority(
    *,
    path: Path,
    sources: Sequence[tuple[Path, str, str]],
    run_id: str,
    statistics_writer: csv.DictWriter,
    species_statistics_writer: csv.DictWriter,
    species: set[str],
) -> tuple[int, int, int]:
    """Stream related source tables to one membership authority and statistics sink."""

    member_count = 0
    group_count = 0
    group_species_count = 0
    with open_text(path=path, mode="w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MEMBERSHIP_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for source, group_type, hierarchy_node in sources:
            source_started = time.perf_counter()
            source_member_start = member_count
            _LOGGER.info(
                "Membership source started: type=%s, node=%s, file=%s",
                group_type,
                hierarchy_node or "ROOT",
                source,
            )
            current_key: tuple[str, str, str] | None = None
            accumulator: GroupAccumulator | None = None
            for row in iter_memberships(
                path=source,
                run_id=run_id,
                group_type=group_type,
                hierarchy_node=hierarchy_node,
            ):
                key = (group_type, hierarchy_node, str(row["group_id"]))
                if key != current_key:
                    if accumulator is not None:
                        _write_accumulator_statistics(
                            accumulator=accumulator,
                            statistics_writer=statistics_writer,
                            species_statistics_writer=species_statistics_writer,
                        )
                        group_count += 1
                        group_species_count += len(accumulator.species_counts)
                    accumulator = GroupAccumulator(
                        run_id=run_id,
                        group_type=group_type,
                        hierarchy_node=hierarchy_node,
                        group_id=str(row["group_id"]),
                        legacy_orthogroup_id=str(row["legacy_orthogroup_id"]),
                        gene_tree_parent_clade=str(row["gene_tree_parent_clade"]),
                        source_file=str(row["source_file"]),
                    )
                    current_key = key
                if accumulator is None:  # pragma: no cover - guarded by key transition
                    raise AssertionError("Membership accumulator was not initialised.")
                accumulator.add_member(species_label=str(row["species_label"]))
                species.add(str(row["species_label"]))
                writer.writerow(row)
                member_count += 1
                source_members = member_count - source_member_start
                if source_members % 1_000_000 == 0:
                    elapsed = max(time.perf_counter() - source_started, 0.001)
                    _LOGGER.info(
                        "Membership source progress: type=%s, node=%s, rows=%s, "
                        "rows_per_second=%.1f",
                        group_type,
                        hierarchy_node or "ROOT",
                        f"{source_members:,}",
                        source_members / elapsed,
                    )
            if accumulator is not None:
                _write_accumulator_statistics(
                    accumulator=accumulator,
                    statistics_writer=statistics_writer,
                    species_statistics_writer=species_statistics_writer,
                )
                group_count += 1
                group_species_count += len(accumulator.species_counts)
            source_members = member_count - source_member_start
            _LOGGER.info(
                "Membership source finished: type=%s, node=%s, rows=%s, "
                "elapsed_seconds=%.3f",
                group_type,
                hierarchy_node or "ROOT",
                f"{source_members:,}",
                time.perf_counter() - source_started,
            )
    return member_count, group_count, group_species_count


def _write_accumulator_statistics(
    *,
    accumulator: GroupAccumulator,
    statistics_writer: csv.DictWriter,
    species_statistics_writer: csv.DictWriter,
) -> None:
    """Write group-wide and per-species statistics for one completed group."""

    statistics_writer.writerow(accumulator.to_record())
    species_statistics_writer.writerows(accumulator.to_species_records())


def _publish_identifiers(
    *,
    tables_dir: Path,
    layout: ResultLayout,
    run_id: str,
    species_from_groups: set[str],
) -> tuple[int, int]:
    """Publish source species and optional sequence-identifier authorities."""

    species_lookup: dict[str, str] = {}
    if layout.species_ids_path is not None:
        species_lookup, species_rows = read_species_ids(path=layout.species_ids_path, run_id=run_id)
    else:
        species_rows = [
            {
                "run_id": run_id,
                "species_index": "",
                "species_label": label,
                "source_fasta": "",
                "source_file": "derived_from_group_table_headings",
                "source_line": index,
            }
            for index, label in enumerate(sorted(species_from_groups), start=1)
        ]
    write_tsv(
        path=_table_path(tables_dir=tables_dir, relation="species"),
        fieldnames=SPECIES_FIELDS,
        records=species_rows,
    )
    if layout.sequence_ids_path is not None and species_lookup:
        sequence_count = write_tsv(
            path=_table_path(tables_dir=tables_dir, relation="sequences"),
            fieldnames=SEQUENCE_FIELDS,
            records=iter_sequence_ids(
                path=layout.sequence_ids_path,
                run_id=run_id,
                species_by_index=species_lookup,
            ),
        )
    else:
        sequence_count = write_tsv(
            path=_table_path(tables_dir=tables_dir, relation="sequences"),
            fieldnames=SEQUENCE_FIELDS,
            records=(),
        )
    return len(species_rows), sequence_count


def _publish_trees(
    *,
    tables_dir: Path,
    layout: ResultLayout,
    run_id: str,
    parse_gene_trees: bool,
    selected_gene_tree_ids: frozenset[str] | None = None,
) -> tuple[list[dict[str, Any]], int, int, int]:
    """Publish tree provenance and optional selected portable representations.

    Args:
        tables_dir: Analytical output directory.
        layout: Validated completed OrthoFinder layout.
        run_id: Immutable resource run identifier.
        parse_gene_trees: Whether gene-tree nodes and edges should be normalised.
        selected_gene_tree_ids: Optional exact focus-related tree identifiers.

    Returns:
        Inventory plus portable payload, node and edge counts.
    """

    _LOGGER.info("Tree checksum inventory started.")
    inventory_started = time.perf_counter()
    inventory = list(iter_tree_inventory(layout=layout, run_id=run_id))
    _LOGGER.info(
        "Tree checksum inventory finished: files=%s, elapsed_seconds=%.3f",
        f"{len(inventory):,}",
        time.perf_counter() - inventory_started,
    )
    write_tsv(
        path=_table_path(tables_dir=tables_dir, relation="tree_inventory"),
        fieldnames=TREE_INVENTORY_FIELDS,
        records=inventory,
    )
    payload_count = write_tsv(
        path=_table_path(tables_dir=tables_dir, relation="tree_payloads"),
        fieldnames=TREE_PAYLOAD_FIELDS,
        records=iter_portable_tree_payloads(
            layout=layout,
            run_id=run_id,
            inventory=inventory,
            selected_tree_ids=selected_gene_tree_ids,
        ),
    )
    _LOGGER.info("Portable tree publication finished: trees=%s", f"{payload_count:,}")
    node_count = 0
    edge_count = 0
    with (
        open_text(
            path=_table_path(tables_dir=tables_dir, relation="tree_nodes"),
            mode="w",
        ) as node_handle,
        open_text(
            path=_table_path(tables_dir=tables_dir, relation="tree_edges"),
            mode="w",
        ) as edge_handle,
    ):
        node_writer = csv.DictWriter(
            node_handle, fieldnames=TREE_NODE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        edge_writer = csv.DictWriter(
            edge_handle, fieldnames=TREE_EDGE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        node_writer.writeheader()
        edge_writer.writeheader()
        for record in inventory:
            if record["tree_type"] != "SPECIES_TREE" and not parse_gene_trees:
                continue
            if (
                record["tree_type"] != "SPECIES_TREE"
                and selected_gene_tree_ids is not None
                and str(record["tree_id"]) not in selected_gene_tree_ids
            ):
                continue
            source = layout.results_dir / str(record["path"])
            nodes, edges = normalise_newick_tree(
                path=source,
                run_id=run_id,
                tree_type=str(record["tree_type"]),
                tree_id=str(record["tree_id"]),
            )
            node_writer.writerows(nodes)
            edge_writer.writerows(edges)
            node_count += len(nodes)
            edge_count += len(edges)
    return inventory, payload_count, node_count, edge_count


def _publish_distances(
    *,
    tables_dir: Path,
    layout: ResultLayout,
    run_id: str,
    alignment_dir: Path | None,
    distance_source: str,
    distance_group_type: str,
    distance_hierarchy_node: str,
    distance_max_groups: int,
    distance_max_members: int,
    selected_group_keys: frozenset[str] | None = None,
    required_members_by_group: Mapping[str, tuple[str, ...]] | None = None,
    persist_pairwise_group_keys: frozenset[str] | None = None,
) -> tuple[int, list[dict[str, Any]]]:
    """Publish optional aligned-sequence or resolved-tree distances.

    Args:
        tables_dir: Analytical output directory.
        layout: Validated completed OrthoFinder layout.
        run_id: Immutable resource run identifier.
        alignment_dir: Optional aligned-sequence authority.
        distance_source: Requested distance authority.
        distance_group_type: Requested distance group collection.
        distance_hierarchy_node: Exact HOG hierarchy node.
        distance_max_groups: Maximum groups for ordinary runs; zero is unlimited.
        distance_max_members: Per-group member bound.
        selected_group_keys: Optional exact focus-selected group keys.
        required_members_by_group: Focus members that bounded samples must retain.
        persist_pairwise_group_keys: Optional subset whose individual pairs are
            published. Every selected group still receives a distance summary.

    Returns:
        Pairwise row count and one summary per attempted group.
    """

    resolved_group_type = (
        ("HOG" if layout.primary_group_authority == "HOG" else "LEGACY_ORTHOGROUP")
        if distance_group_type == "AUTO"
        else distance_group_type
    )
    resolved_source = _resolve_distance_source(
        requested=distance_source,
        alignment_dir=alignment_dir,
        resolved_tree_dir=layout.resolved_gene_trees_dir,
    )
    oversized_required = sorted(
        (key, len(members))
        for key, members in (required_members_by_group or {}).items()
        if len(members) > distance_max_members
    )
    if oversized_required:
        key, count = oversized_required[0]
        raise InputValidationError(
            f"Focus group {key} contains {count:,} matched focus members, exceeding "
            f"distance_max_members={distance_max_members:,}. Increase the member bound "
            "so every focus member remains in the deterministic matrix."
        )
    summaries: list[dict[str, Any]] = []
    pair_count = 0
    with open_text(
        path=_table_path(tables_dir=tables_dir, relation="pairwise_distances"),
        mode="w",
    ) as distance_handle:
        writer = csv.DictWriter(
            distance_handle,
            fieldnames=DISTANCE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        if resolved_source == "ALIGNED_SEQUENCE":
            pair_count = _write_alignment_distances(
                writer=writer,
                summaries=summaries,
                alignment_dir=alignment_dir,
                run_id=run_id,
                group_type=resolved_group_type,
                hierarchy_node=distance_hierarchy_node,
                maximum_groups=distance_max_groups,
                maximum_members=distance_max_members,
                selected_group_keys=selected_group_keys,
                required_members_by_group=required_members_by_group,
                persist_pairwise_group_keys=persist_pairwise_group_keys,
            )
        elif resolved_source == "RESOLVED_GENE_TREE":
            pair_count = _write_tree_distances(
                writer=writer,
                summaries=summaries,
                tables_dir=tables_dir,
                layout=layout,
                run_id=run_id,
                group_type=resolved_group_type,
                hierarchy_node=distance_hierarchy_node,
                maximum_groups=distance_max_groups,
                maximum_members=distance_max_members,
                selected_group_keys=selected_group_keys,
                required_members_by_group=required_members_by_group,
                persist_pairwise_group_keys=persist_pairwise_group_keys,
            )
    write_tsv(
        path=_table_path(tables_dir=tables_dir, relation="distance_statistics"),
        fieldnames=DISTANCE_STATISTIC_FIELDS,
        records=summaries,
    )
    return pair_count, summaries


def _write_alignment_distances(
    *,
    writer: csv.DictWriter,
    summaries: list[dict[str, Any]],
    alignment_dir: Path | None,
    run_id: str,
    group_type: str,
    hierarchy_node: str,
    maximum_groups: int,
    maximum_members: int,
    selected_group_keys: frozenset[str] | None = None,
    required_members_by_group: Mapping[str, tuple[str, ...]] | None = None,
    persist_pairwise_group_keys: frozenset[str] | None = None,
) -> int:
    """Write selected aligned-sequence distances and return their pair count.

    Args:
        writer: Open pairwise-distance TSV writer.
        summaries: Mutable distance-summary sink.
        alignment_dir: Required aligned FASTA directory.
        run_id: Immutable resource run identifier.
        group_type: Exact distance group collection.
        hierarchy_node: Exact HOG hierarchy node.
        maximum_groups: Maximum ordinary-run groups; zero is unlimited.
        maximum_members: Per-group member calculation bound.
        selected_group_keys: Optional exact focus-selected groups.
        required_members_by_group: Focus members retained in bounded samples.
        persist_pairwise_group_keys: Optional selected groups whose pair rows
            are retained; summaries are retained for every calculated group.

    Returns:
        Published pairwise-distance row count.
    """

    if alignment_dir is None:  # pragma: no cover - guarded by source resolution
        raise AssertionError("Aligned-sequence distance source lacks an alignment directory.")
    paths = _alignment_paths(directory=alignment_dir)
    if selected_group_keys is not None:
        paths = [
            path
            for path in paths
            if _group_key(
                {
                    "group_type": group_type,
                    "hierarchy_node": hierarchy_node if group_type == "HOG" else "",
                    "group_id": _group_id_from_alignment(path=path),
                }
            )
            in selected_group_keys
        ]
    if maximum_groups and selected_group_keys is not None and len(paths) > maximum_groups:
        raise InputValidationError(
            "Focus selection exceeds distance_max_groups; focus runs are never truncated."
        )
    if maximum_groups and selected_group_keys is None:
        paths = paths[:maximum_groups]
    pair_count = 0
    for index, path in enumerate(paths, start=1):
        group_id = _group_id_from_alignment(path=path)
        key = _group_key(
            {
                "group_type": group_type,
                "hierarchy_node": hierarchy_node if group_type == "HOG" else "",
                "group_id": group_id,
            }
        )
        group_started = time.perf_counter()
        _LOGGER.info(
            "Distance group started: %s/%s, group=%s, source=alignment",
            index,
            len(paths),
            group_id,
        )
        rows, summary = calculate_alignment_distances(
            sequences=read_fasta(path=path),
            run_id=run_id,
            group_type=group_type,
            hierarchy_node=hierarchy_node if group_type == "HOG" else "",
            group_id=group_id,
            max_members=maximum_members,
            required_member_ids=(required_members_by_group or {}).get(key, ()),
            source_file=str(path),
        )
        persisted = persist_pairwise_group_keys is None or key in (
            persist_pairwise_group_keys
        )
        if persisted:
            writer.writerows(rows)
        summaries.append(summary)
        pair_count += len(rows) if persisted else 0
        _LOGGER.info(
            "Distance group finished: %s/%s, group=%s, status=%s, total_members=%s, "
            "sampled_members=%s, calculated_pairs=%s, persisted=%s, "
            "elapsed_seconds=%.3f",
            index,
            len(paths),
            group_id,
            summary["computation_status"],
            summary["total_member_count"],
            summary["sampled_member_count"],
            f"{len(rows):,}",
            persisted,
            time.perf_counter() - group_started,
        )
    return pair_count


def _write_tree_distances(
    *,
    writer: csv.DictWriter,
    summaries: list[dict[str, Any]],
    tables_dir: Path,
    layout: ResultLayout,
    run_id: str,
    group_type: str,
    hierarchy_node: str,
    maximum_groups: int,
    maximum_members: int,
    selected_group_keys: frozenset[str] | None = None,
    required_members_by_group: Mapping[str, tuple[str, ...]] | None = None,
    persist_pairwise_group_keys: frozenset[str] | None = None,
) -> int:
    """Write selected HOG/orthogroup patristic distances and summaries.

    Args:
        writer: Open pairwise-distance TSV writer.
        summaries: Mutable distance-summary sink.
        tables_dir: Published analytical TSV directory.
        layout: Validated completed OrthoFinder layout.
        run_id: Immutable resource run identifier.
        group_type: Exact distance group collection.
        hierarchy_node: Exact HOG hierarchy node.
        maximum_groups: Maximum ordinary-run groups; zero is unlimited.
        maximum_members: Per-group member calculation bound.
        selected_group_keys: Optional exact focus-selected groups.
        required_members_by_group: Focus members retained in bounded samples.
        persist_pairwise_group_keys: Optional selected groups whose pair rows
            are retained; summaries are retained for every calculated group.

    Returns:
        Published pairwise-distance row count.
    """

    tree_dir = layout.resolved_gene_trees_dir
    if tree_dir is None:  # pragma: no cover - guarded by source resolution
        raise AssertionError("Resolved-tree distance source lacks a tree directory.")
    statistics = [
        row
        for row in read_tsv(
            path=_table_path(tables_dir=tables_dir, relation="group_statistics")
        )
        if row["group_type"] == group_type
        and (group_type != "HOG" or row["hierarchy_node"] == hierarchy_node)
        and (
            selected_group_keys is not None
            or int(row["member_count"]) >= 2
        )
        and (
            selected_group_keys is None
            or _group_key(row) in selected_group_keys
        )
    ]
    statistics.sort(key=lambda row: (-int(row["member_count"]), row["group_id"]))
    if (
        maximum_groups
        and selected_group_keys is not None
        and len(statistics) > maximum_groups
    ):
        raise InputValidationError(
            "Focus selection exceeds distance_max_groups; focus runs are never truncated."
        )
    if maximum_groups and selected_group_keys is None:
        statistics = statistics[:maximum_groups]
    selected_keys = {_group_key(row) for row in statistics}
    membership_path = _table_path(
        tables_dir=tables_dir,
        relation=(
            "hog_memberships" if group_type == "HOG" else "legacy_orthogroup_memberships"
        ),
    )
    species_by_member_by_key: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for row in read_tsv(path=membership_path):
        key = _group_key(row)
        if key in selected_keys:
            species_by_member_by_key[key][row["member_id"]].add(row["species_label"])
    selected_member_species = {
        (member_id, species)
        for by_member in species_by_member_by_key.values()
        for member_id, species_labels in by_member.items()
        for species in species_labels
    }
    internal_ids_by_member_species: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in read_tsv(path=_table_path(tables_dir=tables_dir, relation="sequences")):
        member_species = (row["member_id"], row["species_label"])
        if member_species in selected_member_species:
            internal_ids_by_member_species[member_species].add(row["internal_id"])
    tree_paths = {
        tree_id_from_path(path=path, tree_type="RESOLVED_GENE_TREE"): path.resolve()
        for path in sorted(tree_dir.iterdir())
        if path.is_file() and not path.name.startswith("._")
    }
    pair_count = 0
    for index, statistic in enumerate(statistics, start=1):
        key = _group_key(statistic)
        group_id = statistic["group_id"]
        group_started = time.perf_counter()
        _LOGGER.info(
            "Distance group started: %s/%s, group=%s, total_members=%s, "
            "source=resolved_gene_tree",
            index,
            len(statistics),
            group_id,
            statistic["member_count"],
        )
        species_by_member = species_by_member_by_key.get(key, {})
        members = tuple(sorted(species_by_member))
        if len(members) < 2:
            summaries.append(
                _unavailable_distance_summary(
                    run_id=run_id,
                    group_type=group_type,
                    hierarchy_node=hierarchy_node if group_type == "HOG" else "",
                    group_id=group_id,
                    member_count=len(members),
                    reason="The selected group contains fewer than two members.",
                )
            )
            _LOGGER.info(
                "Distance group finished: %s/%s, group=%s, status=UNAVAILABLE, "
                "reason=fewer_than_two_members, elapsed_seconds=%.3f",
                index,
                len(statistics),
                group_id,
                time.perf_counter() - group_started,
            )
            continue
        tree_candidates = (
            statistic.get("legacy_orthogroup_id", ""),
            group_id,
        )
        tree_path = next(
            (tree_paths[candidate] for candidate in tree_candidates if candidate in tree_paths),
            None,
        )
        if tree_path is None:
            summaries.append(
                _unavailable_distance_summary(
                    run_id=run_id,
                    group_type=group_type,
                    hierarchy_node=hierarchy_node if group_type == "HOG" else "",
                    group_id=group_id,
                    member_count=len(members),
                    reason="No resolved gene tree matched the group or legacy orthogroup ID.",
                )
            )
            _LOGGER.info(
                "Distance group finished: %s/%s, group=%s, status=UNAVAILABLE, "
                "reason=no_matching_tree, elapsed_seconds=%.3f",
                index,
                len(statistics),
                group_id,
                time.perf_counter() - group_started,
            )
            continue
        ambiguous_members = {
            member_id: species_labels
            for member_id, species_labels in species_by_member.items()
            if len(species_labels) != 1
        }
        if ambiguous_members:
            member_id, species_labels = sorted(ambiguous_members.items())[0]
            summaries.append(
                _unavailable_distance_summary(
                    run_id=run_id,
                    group_type=group_type,
                    hierarchy_node=hierarchy_node if group_type == "HOG" else "",
                    group_id=group_id,
                    member_count=len(members),
                    reason=(
                        f"Canonical member {member_id!r} occurs under multiple species: "
                        f"{';'.join(sorted(species_labels))}"
                    ),
                    source_file=str(tree_path),
                )
            )
            _LOGGER.info(
                "Distance group finished: %s/%s, group=%s, status=UNAVAILABLE, "
                "reason=ambiguous_membership, elapsed_seconds=%.3f",
                index,
                len(statistics),
                group_id,
                time.perf_counter() - group_started,
            )
            continue
        member_aliases: dict[str, dict[str, str]] = {}
        for member_id, species_labels in species_by_member.items():
            species_label = next(iter(species_labels))
            aliases = {
                f"{species_label}_{member_id}": "SPECIES_PREFIXED_MEMBER_ID",
            }
            aliases.update(
                {
                    internal_id: "ORTHOFINDER_INTERNAL_ID"
                    for internal_id in internal_ids_by_member_species.get(
                        (member_id, species_label), set()
                    )
                }
            )
            member_aliases[member_id] = aliases
        try:
            rows, summary = calculate_patristic_distances(
                tree_path=tree_path,
                run_id=run_id,
                group_type=group_type,
                hierarchy_node=hierarchy_node if group_type == "HOG" else "",
                group_id=group_id,
                max_members=maximum_members,
                member_ids=members,
                member_aliases=member_aliases,
                required_member_ids=(required_members_by_group or {}).get(key, ()),
                source_file=str(tree_path),
            )
        except DistanceCalculationError as error:
            _LOGGER.warning("Distance unavailable for %s: %s", group_id, error)
            summaries.append(
                _unavailable_distance_summary(
                    run_id=run_id,
                    group_type=group_type,
                    hierarchy_node=hierarchy_node if group_type == "HOG" else "",
                    group_id=group_id,
                    member_count=len(members),
                    reason=str(error),
                    source_file=str(tree_path),
                )
            )
            _LOGGER.info(
                "Distance group finished: %s/%s, group=%s, status=UNAVAILABLE, "
                "reason=calculation_error, elapsed_seconds=%.3f",
                index,
                len(statistics),
                group_id,
                time.perf_counter() - group_started,
            )
            continue
        persisted = persist_pairwise_group_keys is None or key in (
            persist_pairwise_group_keys
        )
        if persisted:
            writer.writerows(rows)
        summaries.append(summary)
        pair_count += len(rows) if persisted else 0
        _LOGGER.info(
            "Distance group finished: %s/%s, group=%s, status=%s, "
            "identifier_resolution=%s, sampled_members=%s, calculated_pairs=%s, "
            "persisted=%s, elapsed_seconds=%.3f",
            index,
            len(statistics),
            group_id,
            summary["computation_status"],
            summary["member_identifier_resolution"],
            summary["sampled_member_count"],
            f"{len(rows):,}",
            persisted,
            time.perf_counter() - group_started,
        )
    return pair_count


def _unavailable_distance_summary(
    *,
    run_id: str,
    group_type: str,
    hierarchy_node: str,
    group_id: str,
    member_count: int,
    reason: str,
    source_file: str = "",
) -> dict[str, Any]:
    """Return a schema-complete unavailable distance record."""

    return summarise_distances(
        rows=(),
        run_id=run_id,
        group_type=group_type,
        hierarchy_node=hierarchy_node,
        group_id=group_id,
        method="patristic_branch_length",
        status="UNAVAILABLE",
        total_member_count=member_count,
        sampled_member_count=0,
        source_file=source_file,
        failure_reason=reason,
    )


def _publish_parquet(*, tables_dir: Path) -> dict[str, Path]:
    """Convert every TSV analytical authority into typed Parquet."""

    parquet_tables: dict[str, Path] = {}
    for tsv_path in sorted(tables_dir.glob("*.tsv.gz")):
        relation = tsv_path.name.removesuffix(".tsv.gz")
        parquet_path = tables_dir / f"{relation}.parquet"
        started = time.perf_counter()
        _LOGGER.info(
            "Parquet conversion started: relation=%s, compressed_tsv_bytes=%s",
            relation,
            f"{tsv_path.stat().st_size:,}",
        )
        row_count = tsv_to_parquet(
            tsv_path=tsv_path,
            parquet_path=parquet_path,
            column_types=GROUP_TYPES.get(relation),
        )
        _LOGGER.info(
            "Parquet conversion finished: relation=%s, rows=%s, parquet_bytes=%s, "
            "elapsed_seconds=%.3f",
            relation,
            f"{row_count:,}",
            f"{parquet_path.stat().st_size:,}",
            time.perf_counter() - started,
        )
        parquet_tables[relation] = parquet_path
    return parquet_tables


def _table_path(*, tables_dir: Path, relation: str) -> Path:
    """Return the compressed TSV authority path for an analytical relation.

    Args:
        tables_dir: Analytical table directory.
        relation: Safe relation/file stem.

    Returns:
        Path ending in ``.tsv.gz``.

    Raises:
        ValueError: If the relation cannot be used as a safe table name.
    """

    if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", relation) is None:
        raise ValueError(f"Unsafe analytical relation name: {relation!r}")
    return tables_dir / f"{relation}.tsv.gz"


def _load_report_group_statistics_and_aggregates(
    *, path: Path, maximum: int
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Load bounded rows and exact compact aggregates for every group level.

    Args:
        path: Complete compressed group-statistics authority.
        maximum: Maximum individual rows embedded in the report.

    Returns:
        Deterministically sampled rows and full-authority histogram aggregates.
    """

    counts: dict[str, int] = defaultdict(int)
    aggregates: dict[str, dict[str, Any]] = {}
    for row in read_tsv(path=path):
        key = _group_level_key(row)
        counts[key] += 1
        accumulator = aggregates.setdefault(
            key,
            {
                "groupType": str(row.get("group_type", "")),
                "hierarchyNode": str(row.get("hierarchy_node", "")),
                "groupCount": 0,
                "singleCopyGroupCount": 0,
                "multicopyGroupCount": 0,
                "memberCount": _new_report_metric(logarithmic=True),
                "speciesCount": _new_report_metric(logarithmic=False),
                "maximumCopiesPerSpecies": _new_report_metric(logarithmic=True),
            },
        )
        accumulator["groupCount"] += 1
        _update_report_metric(
            metric=accumulator["memberCount"], value=row.get("member_count")
        )
        _update_report_metric(
            metric=accumulator["speciesCount"], value=row.get("species_count")
        )
        maximum_copies = _optional_report_int(row.get("max_copies_per_species"))
        _update_report_metric(
            metric=accumulator["maximumCopiesPerSpecies"], value=maximum_copies
        )
        if maximum_copies is not None:
            classification = (
                "singleCopyGroupCount" if maximum_copies <= 1 else "multicopyGroupCount"
            )
            accumulator[classification] += 1
    total = sum(counts.values())
    if maximum == 0 or total <= maximum:
        rows = list(read_tsv(path=path))
    else:
        quotas = _stratified_quotas(counts=counts, maximum=maximum)
        samplers = {
            key: _HashRowSampler(
                maximum=quota,
                salt=f"report_group_statistics\0{key}",
                identifier_field="group_id",
            )
            for key, quota in quotas.items()
            if quota
        }
        for row in read_tsv(path=path):
            sampler = samplers.get(_group_level_key(row))
            if sampler is not None:
                sampler.add(row=row)
        rows = [row for key in sorted(samplers) for row in samplers[key].rows()]
    final_aggregates = {
        key: {
            **{
                field: value
                for field, value in record.items()
                if field
                not in {"memberCount", "speciesCount", "maximumCopiesPerSpecies"}
            },
            "memberCount": _finalise_report_metric(metric=record["memberCount"]),
            "speciesCount": _finalise_report_metric(metric=record["speciesCount"]),
            "maximumCopiesPerSpecies": _finalise_report_metric(
                metric=record["maximumCopiesPerSpecies"]
            ),
        }
        for key, record in sorted(aggregates.items())
    }
    _LOGGER.info(
        "Report group summaries prepared: embedded=%s, complete=%s, levels=%s.",
        f"{len(rows):,}",
        f"{total:,}",
        f"{len(final_aggregates):,}",
    )
    return rows, final_aggregates


def _load_report_group_statistics(*, path: Path, maximum: int) -> list[dict[str, str]]:
    """Load a deterministic hierarchy-stratified subset for HTML embedding."""

    rows, _ = _load_report_group_statistics_and_aggregates(
        path=path,
        maximum=maximum,
    )
    return rows


def _new_report_metric(*, logarithmic: bool) -> dict[str, Any]:
    """Return an empty streaming accumulator for a report overview metric."""

    return {
        "binning": "LOG2_INTEGER" if logarithmic else "INTEGER",
        "count": 0,
        "sum": 0,
        "minimum": None,
        "maximum": None,
        "bins": defaultdict(int),
    }


def _optional_report_int(value: Any) -> int | None:
    """Return a non-negative integer or ``None`` for an unavailable value."""

    if value in {None, ""}:
        return None
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed >= 0 else None


def _update_report_metric(*, metric: dict[str, Any], value: Any) -> None:
    """Update an exact streaming report aggregate when a value is available."""

    parsed = _optional_report_int(value)
    if parsed is None:
        return
    metric["count"] += 1
    metric["sum"] += parsed
    metric["minimum"] = parsed if metric["minimum"] is None else min(metric["minimum"], parsed)
    metric["maximum"] = parsed if metric["maximum"] is None else max(metric["maximum"], parsed)
    bin_index = (
        int(math.log2(max(1, parsed)))
        if metric["binning"] == "LOG2_INTEGER"
        else parsed
    )
    metric["bins"][bin_index] += 1


def _finalise_report_metric(*, metric: Mapping[str, Any]) -> dict[str, Any]:
    """Convert one streaming metric accumulator into a compact JSON record."""

    binning = str(metric["binning"])
    bins = dict(metric["bins"])
    ordered_indices = (
        list(range(min(bins), max(bins) + 1)) if bins else []
    )
    labels = [
        (
            "1"
            if index == 0 and binning == "LOG2_INTEGER"
            else (
                f"{2**index}–{2 ** (index + 1) - 1}"
                if binning == "LOG2_INTEGER"
                else str(index)
            )
        )
        for index in ordered_indices
    ]
    count = int(metric["count"])
    return {
        "binning": binning,
        "count": count,
        "minimum": metric["minimum"],
        "maximum": metric["maximum"],
        "mean": (float(metric["sum"]) / count) if count else None,
        "labels": labels,
        "counts": [int(bins.get(index, 0)) for index in ordered_indices],
    }


def _stratified_quotas(*, counts: Mapping[str, int], maximum: int) -> dict[str, int]:
    """Allocate a bounded row budget across hierarchy levels proportionally.

    Args:
        counts: Available rows keyed by group type and hierarchy node.
        maximum: Total row budget.

    Returns:
        Per-level quotas summing to at most the requested maximum.
    """

    ordered = sorted(counts)
    if maximum < len(ordered):
        retained = sorted(ordered, key=lambda key: (-counts[key], key))[:maximum]
        return {key: int(key in retained) for key in ordered}
    quotas = {key: min(1, counts[key]) for key in ordered}
    remaining = maximum - sum(quotas.values())
    capacity = {key: max(0, counts[key] - quotas[key]) for key in ordered}
    total_capacity = sum(capacity.values())
    if not remaining or not total_capacity:
        return quotas
    targets = {key: remaining * capacity[key] / total_capacity for key in ordered}
    for key in ordered:
        addition = min(capacity[key], int(targets[key]))
        quotas[key] += addition
        remaining -= addition
    for key in sorted(
        ordered,
        key=lambda item: (-(targets[item] - int(targets[item])), item),
    ):
        if not remaining:
            break
        if quotas[key] < counts[key]:
            quotas[key] += 1
            remaining -= 1
    return quotas


def _group_level_key(row: Mapping[str, Any]) -> str:
    """Return the group type and hierarchy node used for report stratification."""

    return "|".join(
        (str(row.get("group_type", "")), str(row.get("hierarchy_node", "")))
    )


def _load_report_network_data(
    *,
    tables_dir: Path,
    group_statistics_path: Path,
    distance_summaries: Sequence[Mapping[str, Any]],
    maximum_groups: int,
    maximum_members: int,
) -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    """Load bounded network members and pair distances for the interactive report."""

    distance_keys = {_group_key(row) for row in distance_summaries}
    selected_statistics = heapq.nlargest(
        maximum_groups,
        read_tsv(path=group_statistics_path),
        key=lambda row: (
            _group_key(row) in distance_keys,
            int(row.get("member_count", 0)),
            _group_key(row),
        ),
    )
    selected_keys = {_group_key(row) for row in selected_statistics}
    report_distances = [
        row
        for row in read_tsv(
            path=_table_path(tables_dir=tables_dir, relation="pairwise_distances")
        )
        if _group_key(row) in selected_keys
    ]
    distance_members: dict[str, set[str]] = defaultdict(set)
    for row in report_distances:
        key = _group_key(row)
        distance_members[key].update((row["member_a"], row["member_b"]))

    required_ids = {
        key: set(_sample_member_ids(member_ids=member_ids, maximum=maximum_members, salt=key))
        for key, member_ids in distance_members.items()
    }
    samplers = {
        key: _HashRowSampler(maximum=maximum_members, salt=key)
        for key in selected_keys
        if key not in required_ids
    }
    required_rows: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for membership_path in (
        _table_path(tables_dir=tables_dir, relation="legacy_orthogroup_memberships"),
        _table_path(tables_dir=tables_dir, relation="hog_memberships"),
    ):
        for row in read_tsv(path=membership_path):
            key = _group_key(row)
            if key not in selected_keys:
                continue
            if key in required_ids:
                if row["member_id"] in required_ids[key]:
                    required_rows[key].setdefault(row["member_id"], row)
            else:
                samplers[key].add(row=row)
    report_members = []
    for key in sorted(selected_keys):
        if key not in required_ids:
            report_members.extend(samplers[key].rows())
            continue
        selected = required_rows[key]
        # Distance identifiers should be canonical membership identifiers. If a
        # future adapter supplies a different scheme, retain every endpoint as an
        # explicit unresolved node rather than silently constructing another sample.
        for member_id in sorted(required_ids[key]):
            report_members.append(
                selected.get(
                    member_id,
                    {
                        "run_id": "",
                        "group_type": key.split("|", maxsplit=2)[0],
                        "hierarchy_node": key.split("|", maxsplit=2)[1],
                        "group_id": key.split("|", maxsplit=2)[2],
                        "legacy_orthogroup_id": "",
                        "gene_tree_parent_clade": "",
                        "species_label": "unresolved_distance_identifier",
                        "member_id": member_id,
                        "source_file": "distance_table",
                        "source_row": "",
                    },
                )
            )
    return selected_statistics, report_members, report_distances


def _sample_member_ids(
    *, member_ids: set[str], maximum: int, salt: str
) -> list[str]:
    """Select a deterministic bounded identifier subset.

    Args:
        member_ids: Unique candidate member identifiers.
        maximum: Maximum identifiers retained.
        salt: Stable group-specific hash salt.

    Returns:
        Selected member identifiers in lexical order.
    """

    ranked = sorted(
        member_ids,
        key=lambda member_id: hashlib.sha256(
            f"{salt}\0{member_id}".encode("utf-8")
        ).hexdigest(),
    )
    return sorted(ranked[:maximum])


def _load_report_group_species_data(
    *,
    tables_dir: Path,
    memberships: Sequence[Mapping[str, Any]],
) -> list[dict[str, str]]:
    """Load full species copy counts for the bounded network groups.

    Args:
        tables_dir: Analytical table directory.
        memberships: Bounded report memberships identifying selected groups.

    Returns:
        Species-level rows for the selected group keys.
    """

    selected_keys = {_group_key(row) for row in memberships}
    if not selected_keys:
        return []
    return [
        row
        for row in read_tsv(
            path=_table_path(tables_dir=tables_dir, relation="group_species_statistics")
        )
        if _group_key(row) in selected_keys
    ]


def _load_report_tree_data(
    *,
    tables_dir: Path,
    group_statistics: Sequence[Mapping[str, Any]],
    memberships: Sequence[Mapping[str, Any]],
    distance_summaries: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, str]]]:
    """Load only tree and identifier records needed by report-selected groups.

    Normalised resource tables are preferred. When an older completed resource
    checksum-inventoried trees but did not expand every gene tree, the exact
    distance-summary source may be parsed only after its SHA-256 matches the
    immutable tree inventory. This fallback changes no analytical authority.

    Args:
        tables_dir: Completed resource table directory or node-local copy.
        group_statistics: Exact groups selected for interactive views.
        memberships: Exact bounded members selected for those groups.
        distance_summaries: Complete distance summaries with source provenance.

    Returns:
        Selected normalised tree nodes, edges and sequence-identifier aliases.
    """

    if not group_statistics or not memberships:
        return [], [], []
    candidate_ids: set[str] = set()
    candidates_by_key: dict[str, list[str]] = {}
    for row in group_statistics:
        candidates = list(
            dict.fromkeys(
                str(value)
                for value in (
                    row.get("legacy_orthogroup_id", ""),
                    row.get("group_id", ""),
                )
                if value
            )
        )
        candidates_by_key[_group_key(row)] = candidates
        candidate_ids.update(candidates)
    selected_nodes = [
        dict(row)
        for row in read_tsv(
            path=_table_path(tables_dir=tables_dir, relation="tree_nodes")
        )
        if row.get("tree_type") == "RESOLVED_GENE_TREE"
        and row.get("tree_id") in candidate_ids
    ]
    selected_edges = [
        dict(row)
        for row in read_tsv(
            path=_table_path(tables_dir=tables_dir, relation="tree_edges")
        )
        if row.get("tree_type") == "RESOLVED_GENE_TREE"
        and row.get("tree_id") in candidate_ids
    ]
    available_tree_ids = {str(row.get("tree_id", "")) for row in selected_nodes}
    summaries_by_key = {_group_key(row): row for row in distance_summaries}
    inventory = [
        row
        for row in read_tsv(
            path=_table_path(tables_dir=tables_dir, relation="tree_inventory")
        )
        if row.get("tree_type") == "RESOLVED_GENE_TREE"
        and row.get("tree_id") in candidate_ids
    ]
    inventory_by_id = {str(row.get("tree_id", "")): row for row in inventory}
    for group_key, candidates in candidates_by_key.items():
        if any(candidate in available_tree_ids for candidate in candidates):
            continue
        summary = summaries_by_key.get(group_key, {})
        source_text = str(summary.get("source_file", ""))
        tree_id = next(
            (candidate for candidate in candidates if candidate in inventory_by_id),
            "",
        )
        if not source_text or not tree_id:
            continue
        source = Path(source_text).expanduser().resolve()
        inventory_row = inventory_by_id[tree_id]
        if not source.is_file():
            _LOGGER.warning(
                "Report phylogram unavailable: checksum-inventoried source is missing: %s",
                source,
            )
            continue
        observed_digest = sha256_file(path=source)
        expected_digest = str(inventory_row.get("sha256", ""))
        if not expected_digest or observed_digest != expected_digest:
            _LOGGER.warning(
                "Report phylogram source rejected after SHA-256 verification: %s",
                source,
            )
            continue
        nodes, edges = normalise_newick_tree(
            path=source,
            run_id=str(summary.get("run_id", "")),
            tree_type="RESOLVED_GENE_TREE",
            tree_id=tree_id,
        )
        selected_nodes.extend(nodes)
        selected_edges.extend(edges)
        available_tree_ids.add(tree_id)
        _LOGGER.info(
            "Report phylogram source expanded after checksum verification: "
            "tree=%s, nodes=%s, edges=%s.",
            tree_id,
            f"{len(nodes):,}",
            f"{len(edges):,}",
        )

    selected_member_species = {
        (str(row.get("member_id", "")), str(row.get("species_label", "")))
        for row in memberships
    }
    sequence_identifiers = [
        row
        for row in read_tsv(
            path=_table_path(tables_dir=tables_dir, relation="sequences")
        )
        if (str(row.get("member_id", "")), str(row.get("species_label", "")))
        in selected_member_species
    ]
    _LOGGER.info(
        "Report tree data loaded: candidate_trees=%s, available_trees=%s, "
        "nodes=%s, edges=%s, identifier_aliases=%s.",
        f"{len(candidate_ids):,}",
        f"{len(available_tree_ids):,}",
        f"{len(selected_nodes):,}",
        f"{len(selected_edges):,}",
        f"{len(sequence_identifiers):,}",
    )
    return selected_nodes, selected_edges, sequence_identifiers


class _HashRowSampler:
    """Deterministic bounded row sampler independent of input order."""

    def __init__(
        self, *, maximum: int, salt: str, identifier_field: str = "member_id"
    ) -> None:
        """Initialise a sampler.

        Args:
            maximum: Maximum retained unique rows.
            salt: Stable group-specific hash salt.
            identifier_field: Row field defining uniqueness and stable order.
        """

        self.maximum = maximum
        self.salt = salt
        self.identifier_field = identifier_field
        self._rows: dict[str, tuple[str, dict[str, str]]] = {}

    def add(self, *, row: Mapping[str, str]) -> None:
        """Consider one membership row for retention.

        Args:
            row: Record containing the configured identifier field.
        """

        identifier = row[self.identifier_field]
        if identifier in self._rows:
            return
        rank = hashlib.sha256(f"{self.salt}\0{identifier}".encode("utf-8")).hexdigest()
        self._rows[identifier] = (rank, dict(row))
        if len(self._rows) > self.maximum:
            worst = max(self._rows, key=lambda key: self._rows[key][0])
            del self._rows[worst]

    def rows(self) -> list[dict[str, str]]:
        """Return retained rows ordered by member identifier.

        Returns:
            Deterministically selected row mappings.
        """

        return [self._rows[key][1] for key in sorted(self._rows)]


def _build_source_inventory(
    *,
    layout: ResultLayout,
    alignment_dir: Path | None,
    focus_proteins_path: Path | None = None,
    benchmark_proteins_path: Path | None = None,
) -> list[dict[str, Any]]:
    """Checksum every source that can alter the requested analytical result.

    Args:
        layout: Validated completed OrthoFinder layout.
        alignment_dir: Optional alignment authority.
        focus_proteins_path: Optional E3/focus selection authority.
        benchmark_proteins_path: Optional housekeeping/R marker authority.

    Returns:
        Complete input file inventory in deterministic role/path order.
    """

    roles: list[tuple[str, Path]] = [("orthofinder_log", layout.log_path)]
    for role, path in (
        ("species_ids", layout.species_ids_path),
        ("sequence_ids", layout.sequence_ids_path),
        ("legacy_orthogroups", layout.orthogroups_path),
        ("species_tree", layout.species_tree_path),
    ):
        if path is not None:
            roles.append((role, path))
    roles.extend(("hog_table", path) for path in layout.hog_paths)
    for role, directory in (
        ("gene_tree", layout.gene_trees_dir),
        ("resolved_gene_tree", layout.resolved_gene_trees_dir),
    ):
        if directory is not None:
            roles.extend(
                (role, path)
                for path in sorted(directory.iterdir())
                if path.is_file() and not path.name.startswith("._")
            )
    if alignment_dir is not None:
        roles.extend(("alignment", path) for path in _alignment_paths(directory=alignment_dir))
    if focus_proteins_path is not None:
        roles.append(("focus_protein_authority", focus_proteins_path))
    if benchmark_proteins_path is not None:
        roles.append(("dispersion_benchmark_authority", benchmark_proteins_path))
    records = []
    for role, path in roles:
        record = file_record(path=path)
        records.append({"role": role, **record})
    return records


def _inventory_digest(*, records: Sequence[Mapping[str, Any]]) -> str:
    """Calculate a deterministic digest over an input inventory."""

    digest = hashlib.sha256()
    for record in sorted(records, key=lambda row: (str(row["role"]), str(row["path"]))):
        digest.update(
            json.dumps(dict(record), sort_keys=True, separators=(",", ":")).encode("utf-8")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _resolve_existing_output(
    *,
    output_dir: Path,
    run_id: str,
    input_digest: str,
    resume: bool,
    force: bool,
) -> dict[str, Any] | None:
    """Reuse, reject or supersede an existing formal output."""

    if not output_dir.exists():
        return None
    manifest_path = output_dir / "run_manifest.json"
    if resume and manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            manifest.get("status") == "complete"
            and manifest.get("run_id") == run_id
            and manifest.get("input_digest") == input_digest
            and manifest.get("package_version") == __version__
        ):
            _LOGGER.info("Reusing checksum-matched completed output: %s", output_dir)
            return manifest
    if not force:
        raise PublicationError(
            f"Output directory already exists and was not reusable: {output_dir}. "
            "Use --force only after reviewing it."
        )
    timestamp = utc_now_iso().replace(":", "").replace("-", "")
    superseded = output_dir.with_name(f"{output_dir.name}.superseded.{timestamp}")
    os.replace(output_dir, superseded)
    _LOGGER.warning("Moved existing output to recoverable location: %s", superseded)
    return None


def _resolve_pipeline_alignment_dir(
    *,
    requested: Path | None,
    discovered: Path | None,
    distance_source: str,
) -> Path | None:
    """Resolve alignments only when the selected distance policy can use them."""

    if distance_source in {"NONE", "RESOLVED_GENE_TREE"}:
        return None
    return _resolve_alignment_dir(requested=requested, discovered=discovered)


def _resolve_distance_source(
    *,
    requested: str,
    alignment_dir: Path | None,
    resolved_tree_dir: Path | None,
) -> str:
    """Resolve an explicit or capability-driven distance authority."""

    has_resolved_trees = resolved_tree_dir is not None and any(
        path.is_file() and not path.name.startswith("._") for path in resolved_tree_dir.iterdir()
    )
    if requested == "AUTO":
        if alignment_dir is not None:
            return "ALIGNED_SEQUENCE"
        return "RESOLVED_GENE_TREE" if has_resolved_trees else "NONE"
    if requested == "ALIGNED_SEQUENCE" and alignment_dir is None:
        raise InputValidationError(
            "ALIGNED_SEQUENCE distance source requires a recognised alignment directory."
        )
    if requested == "RESOLVED_GENE_TREE" and not has_resolved_trees:
        raise InputValidationError(
            "RESOLVED_GENE_TREE distance source requires resolved gene-tree files."
        )
    return requested


def _resolve_alignment_dir(*, requested: Path | None, discovered: Path | None) -> Path | None:
    """Resolve an explicit or OrthoFinder-supplied alignment directory."""

    candidate = requested if requested is not None else discovered
    if candidate is None:
        return None
    resolved = Path(candidate).expanduser().resolve()
    if not resolved.is_dir():
        raise InputValidationError(f"Alignment directory does not exist: {resolved}")
    if not _alignment_paths(directory=resolved):
        raise InputValidationError(
            f"Alignment directory contains no recognised FASTA files: {resolved}"
        )
    return resolved


def _alignment_paths(*, directory: Path) -> list[Path]:
    """Return recognised non-sidecar alignment files in deterministic order."""

    return sorted(
        path.resolve()
        for path in directory.iterdir()
        if path.is_file()
        and not path.name.startswith("._")
        and path.suffix.lower() in _ALIGNMENT_SUFFIXES
    )


def _group_id_from_alignment(*, path: Path) -> str:
    """Return an exact filename-derived group identifier."""

    group_id = path.stem
    if not group_id:
        raise InputValidationError(f"Could not derive group identifier from alignment: {path}")
    return group_id


def _group_key(row: Mapping[str, Any]) -> str:
    """Return the report/database collision-safe group key."""

    return "|".join(
        (
            str(row.get("group_type", "")),
            str(row.get("hierarchy_node", "")),
            str(row.get("group_id", "")),
        )
    )


def _qc_rows(
    *,
    layout: ResultLayout,
    membership_counts: Mapping[str, int],
    group_count: int,
    species_count: int,
    sequence_count: int,
    tree_inventory_count: int,
    tree_payload_count: int,
    tree_node_count: int,
    distance_count: int,
    offline_report: bool,
    focus_selection: FocusSelection | None = None,
    focus_result_rows: Sequence[Mapping[str, Any]] = (),
    distance_summaries: Sequence[Mapping[str, Any]] = (),
    benchmark_plan: BenchmarkPlan | None = None,
    benchmark_counts: Mapping[str, int] | None = None,
) -> list[dict[str, Any]]:
    """Build explicit run validation checks.

    Args:
        layout: Validated source layout.
        membership_counts: Published membership counts by authority.
        group_count: Published group-statistics row count.
        species_count: Published species count.
        sequence_count: Published sequence-identifier count.
        tree_inventory_count: Checksum-inventoried tree count.
        tree_payload_count: Portable selected gene-tree count.
        tree_node_count: Normalised selected tree-node count.
        distance_count: Published pairwise-distance count.
        offline_report: Whether the report contains no HTTP dependencies.
        focus_selection: Optional complete E3/focus selection.
        focus_result_rows: Flat selected-cluster result rows.
        distance_summaries: Explicit selected distance summaries.
        benchmark_plan: Optional complete matched-background selection.
        benchmark_counts: Published benchmark relation row counts.

    Returns:
        Required validation rows, including focus reconciliation when enabled.
    """

    rows = [
        _qc(
            "supported_adapter",
            layout.adapter_name in {"orthofinder_2", "orthofinder_3"},
            layout.adapter_name,
            "orthofinder_2|orthofinder_3",
            "Detected version-specific adapter.",
        ),
        _qc(
            "group_memberships_present",
            sum(membership_counts.values()) > 0,
            sum(membership_counts.values()),
            ">0",
            "At least one group membership is required.",
        ),
        _qc(
            "group_statistics_present",
            group_count > 0,
            group_count,
            ">0",
            "At least one run-scoped group is required.",
        ),
        _qc(
            "species_present",
            species_count > 0,
            species_count,
            ">0",
            "Species are sourced from SpeciesIDs or group headings.",
        ),
        _qc(
            "sequence_identifier_status",
            sequence_count > 0 or not layout.capabilities.has_sequence_ids,
            sequence_count,
            ">0 when SequenceIDs and SpeciesIDs are available",
            "SequenceIDs are optional in portable completed results.",
        ),
        _qc(
            "tree_inventory_status",
            tree_inventory_count > 0
            or not (layout.capabilities.has_species_tree or layout.capabilities.has_gene_trees),
            tree_inventory_count,
            ">0 when trees are available",
            "Every discovered tree is checksum inventoried.",
        ),
        _qc(
            "tree_normalisation_status",
            tree_node_count > 0 or not layout.capabilities.has_species_tree,
            tree_node_count,
            ">0 when a species tree is available",
            "Species tree nodes are always normalised; gene trees are optional.",
        ),
        _qc(
            "portable_gene_tree_status",
            tree_payload_count > 0
            or not (
                layout.capabilities.has_resolved_gene_trees
                or layout.capabilities.has_gene_trees
            ),
            tree_payload_count,
            ">0 when gene trees are available",
            "Preferred gene-tree Newick is embedded for portable lazy distances.",
        ),
        _qc(
            "distance_status",
            True,
            distance_count,
            ">=0",
            "Distance absence is explicit when no alignment directory is supplied.",
        ),
        _qc(
            "offline_html_report",
            offline_report,
            str(offline_report).lower(),
            "true",
            "The HTML report must not depend on HTTP resources.",
        ),
    ]
    if focus_selection is not None:
        selected_count = len(focus_selection.group_statistics)
        summary_keys = {_group_key(row) for row in distance_summaries}
        rows.extend(
            (
                _qc(
                    "focus_seed_matches_present",
                    focus_selection.matched_seed_count > 0,
                    focus_selection.matched_seed_count,
                    ">0",
                    "At least one enabled E3/focus seed must match exactly.",
                ),
                _qc(
                    "focus_seed_audit_complete",
                    len(focus_selection.audit_rows)
                    == len(focus_selection.authority.records),
                    len(focus_selection.audit_rows),
                    len(focus_selection.authority.records),
                    "Every configured seed has a matched, unmatched or disabled audit row.",
                ),
                _qc(
                    "focus_cluster_export_complete",
                    len(focus_result_rows) == selected_count,
                    len(focus_result_rows),
                    selected_count,
                    "The compressed cluster result contains every matched focus group.",
                ),
                _qc(
                    "focus_distance_summary_complete",
                    focus_selection.group_keys.issubset(summary_keys),
                    len(focus_selection.group_keys.intersection(summary_keys)),
                    selected_count,
                    "Every focus group has an explicit distance result or reason.",
                ),
            )
        )
    if benchmark_plan is not None:
        summary_keys = {_group_key(row) for row in distance_summaries}
        expected_groups = len(benchmark_plan.selected_group_keys)
        target_groups = len(benchmark_plan.target_group_keys)
        controls = len(benchmark_plan.control_group_keys)
        counts = dict(benchmark_counts or {})
        linked_anchors = {
            _group_key(
                {
                    "group_type": row["anchor_group_type"],
                    "hierarchy_node": row["anchor_hierarchy_node"],
                    "group_id": row["anchor_group_id"],
                }
            )
            for row in benchmark_plan.control_links
        }
        rows.extend(
            (
                _qc(
                    "benchmark_marker_audit_complete",
                    len(benchmark_plan.marker_selection.audit_rows)
                    == len(benchmark_plan.marker_authority.records),
                    len(benchmark_plan.marker_selection.audit_rows),
                    len(benchmark_plan.marker_authority.records),
                    "Every marker has a matched, unmatched or disabled audit state.",
                ),
                _qc(
                    "benchmark_matched_controls_complete",
                    benchmark_plan.target_group_keys.issubset(linked_anchors),
                    len(benchmark_plan.target_group_keys.intersection(linked_anchors)),
                    target_groups,
                    "Every target has one or more distance-blind matched controls.",
                ),
                _qc(
                    "benchmark_distance_summaries_complete",
                    benchmark_plan.selected_group_keys.issubset(summary_keys),
                    len(benchmark_plan.selected_group_keys.intersection(summary_keys)),
                    expected_groups,
                    "Every target and matched control has a result or explicit reason.",
                ),
                _qc(
                    "benchmark_cluster_export_complete",
                    counts.get("benchmark_cluster_count", -1) == expected_groups,
                    counts.get("benchmark_cluster_count", -1),
                    expected_groups,
                    "The cluster export contains every target and matched control.",
                ),
                _qc(
                    "benchmark_target_and_control_groups_present",
                    target_groups > 0 and controls > 0,
                    f"targets={target_groups};controls={controls}",
                    "targets>0;controls>0",
                    "Benchmarking requires biological targets and non-focus controls.",
                ),
            )
        )
    return rows


def _qc(name: str, passed: bool, observed: Any, expected: Any, details: str) -> dict[str, Any]:
    """Return one standard validation row."""

    return {
        "check_name": name,
        "status": "PASS" if passed else "FAIL",
        "observed_value": observed,
        "expected_value": expected,
        "details": details,
    }


def _validate_controls(
    *,
    run_id: str,
    distance_source: str,
    distance_group_type: str,
    distance_max_groups: int,
    distance_max_members: int,
    report_max_statistic_rows: int,
    report_max_groups: int,
    report_max_members: int,
    report_nearest_neighbours: int,
    resume: bool,
    force: bool,
    focus_enabled: bool = False,
    focus_group_type: str = "HOG",
    focus_hierarchy_node: str = "N0",
    distance_hierarchy_node: str = "N0",
    benchmark_enabled: bool = False,
    benchmark_controls_per_group: int = 3,
    benchmark_bootstrap_resamples: int = 1_000,
) -> None:
    """Validate named execution controls before filesystem mutation.

    Args:
        run_id: Immutable output run identifier.
        distance_source: Requested distance authority.
        distance_group_type: Requested distance group collection.
        distance_max_groups: Maximum ordinary-run groups or zero for unlimited.
        distance_max_members: Per-group member calculation bound.
        report_max_statistic_rows: Browser-safe summary-row bound.
        report_max_groups: Browser-safe network-group bound.
        report_max_members: Browser-safe network-member bound.
        report_nearest_neighbours: Retained neighbour count.
        resume: Whether an exact completed output may be reused.
        force: Whether an existing output may be superseded.
        focus_enabled: Whether the run is restricted to a focus authority.
        focus_group_type: Exact focus group collection.
        focus_hierarchy_node: Exact focus hierarchy node.
        distance_hierarchy_node: Exact distance hierarchy node.
        benchmark_enabled: Whether matched-background benchmarking is enabled.
        benchmark_controls_per_group: Unique controls selected per target.
        benchmark_bootstrap_resamples: Deterministic confidence-interval iterations.

    Raises:
        InputValidationError: If any controls are unsafe or inconsistent.
    """

    if _RUN_ID_PATTERN.fullmatch(run_id) is None:
        raise InputValidationError(
            "run_id must begin with an alphanumeric character and contain only "
            "letters, numbers, underscores, dots and hyphens."
        )
    if distance_source not in {
        "AUTO",
        "ALIGNED_SEQUENCE",
        "RESOLVED_GENE_TREE",
        "NONE",
    }:
        raise InputValidationError(f"Unsupported distance_source: {distance_source}")
    if distance_group_type not in {"AUTO", "HOG", "LEGACY_ORTHOGROUP"}:
        raise InputValidationError(f"Unsupported distance_group_type: {distance_group_type}")
    if focus_enabled:
        if distance_source != "RESOLVED_GENE_TREE":
            raise InputValidationError(
                "Focus precursor runs require RESOLVED_GENE_TREE distances."
            )
        if focus_group_type not in {"HOG", "LEGACY_ORTHOGROUP"}:
            raise InputValidationError(
                f"Unsupported focus_group_type: {focus_group_type}"
            )
        if focus_group_type == "LEGACY_ORTHOGROUP" and focus_hierarchy_node:
            raise InputValidationError(
                "Legacy focus runs require an empty ROOT hierarchy node."
            )
        resolved_distance_type = (
            focus_group_type if distance_group_type == "AUTO" else distance_group_type
        )
        if resolved_distance_type != focus_group_type:
            raise InputValidationError(
                "Focus and distance group types must identify the same collection."
            )
        if distance_hierarchy_node != focus_hierarchy_node:
            raise InputValidationError(
                "Focus and distance hierarchy nodes must be identical."
            )
    if benchmark_enabled:
        if not focus_enabled:
            raise InputValidationError(
                "Dispersion benchmarking requires an E3/focus authority."
            )
        if focus_group_type != "HOG" or focus_hierarchy_node != "N0":
            raise InputValidationError(
                "Dispersion benchmarking currently requires HOGs at hierarchy N0."
            )
        if not 1 <= benchmark_controls_per_group <= 50:
            raise InputValidationError(
                "benchmark_controls_per_group must be between 1 and 50."
            )
        if not 100 <= benchmark_bootstrap_resamples <= 100_000:
            raise InputValidationError(
                "benchmark_bootstrap_resamples must be between 100 and 100,000."
            )
    _validate_report_controls(
        report_max_statistic_rows=report_max_statistic_rows,
        report_max_groups=report_max_groups,
        report_max_members=report_max_members,
        report_nearest_neighbours=report_nearest_neighbours,
    )
    non_negative = {
        "distance_max_groups": distance_max_groups,
    }
    for name, value in non_negative.items():
        if value < 0:
            raise InputValidationError(f"{name} must not be negative.")
    positive = {"distance_max_members": distance_max_members}
    for name, value in positive.items():
        if value <= 0:
            raise InputValidationError(f"{name} must be positive.")
    if distance_max_members < 2:
        raise InputValidationError("Distance member limit must be at least two.")
    if resume and force:
        raise InputValidationError("--resume and --force are mutually exclusive.")


def _validate_report_controls(
    *,
    report_max_statistic_rows: int,
    report_max_groups: int,
    report_max_members: int,
    report_nearest_neighbours: int,
) -> None:
    """Validate controls shared by full runs and report-only regeneration."""

    if not 1 <= report_max_statistic_rows <= _MAX_REPORT_STATISTIC_ROWS:
        raise InputValidationError(
            "report_max_statistic_rows must be between 1 and "
            f"{_MAX_REPORT_STATISTIC_ROWS:,} for browser safety."
        )
    positive = {
        "report_max_groups": report_max_groups,
        "report_max_members": report_max_members,
        "report_nearest_neighbours": report_nearest_neighbours,
    }
    for name, value in positive.items():
        if value <= 0:
            raise InputValidationError(f"{name} must be positive.")
    if report_max_members < 2:
        raise InputValidationError("Report member limit must be at least two.")
    pair_budget = (
        report_max_groups * report_max_members * (report_max_members - 1) // 2
    )
    if pair_budget > _MAX_REPORT_DISTANCE_PAIRS:
        raise InputValidationError(
            "The requested report network bounds could embed up to "
            f"{pair_budget:,} pair distances; the browser-safety limit is "
            f"{_MAX_REPORT_DISTANCE_PAIRS:,}. Reduce --report-max-groups or "
            "--report-max-members."
        )

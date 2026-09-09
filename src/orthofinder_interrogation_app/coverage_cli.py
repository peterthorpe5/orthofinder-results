"""Command-line orchestration for reproducible selection coverage trees."""

from __future__ import annotations

import argparse
import hashlib
import logging
from pathlib import Path
from typing import Any

from orthofinder_results import __version__
from orthofinder_results.errors import InputValidationError

from .coverage_exports import build_selection_coverage_run, publish_selection_coverage_run
from .focus import bundled_focus_path, read_focus_proteins
from .models import FocusClusterFilters
from .queries import OrthoFinderQueryService
from .resource import open_resource
from .taxonomy import read_taxonomy_mapping
from .taxonomy_selection import (
    build_taxonomy_graph,
    default_expected_taxa,
    make_selection,
    read_expected_taxa,
    selection_summary,
)

_LOGGER = logging.getLogger("orthofinder_interrogation_app.coverage_cli")


def run_coverage_tree_action(*, arguments: argparse.Namespace) -> dict[str, Any]:
    """Validate, dry-run or atomically publish one command-line selection.

    Args:
        arguments: Parsed named command-line options from ``orthofinder-results``.

    Returns:
        Complete selection manifest or validation summary.

    Raises:
        InputValidationError: If inputs, selectors, bounds or targets are invalid.
    """

    resource = open_resource(path=arguments.resource_dir)
    service = OrthoFinderQueryService(resource=resource)
    authority = read_taxonomy_mapping(
        path=arguments.taxonomy_map,
        expected_species=service.list_species(),
    )
    graph = build_taxonomy_graph(authority=authority)
    expected = (
        read_expected_taxa(path=arguments.expected_taxa, graph=graph)
        if arguments.expected_taxa is not None
        else default_expected_taxa(graph=graph)
    )
    selection = make_selection(
        graph=graph,
        required_exact=tuple(arguments.require_exact_tax_id),
        include_clade=tuple(arguments.include_clade_tax_id),
        only_in_clade=tuple(arguments.only_in_clade_tax_id),
        exclude_exact=tuple(arguments.exclude_exact_tax_id),
        exclude_clade=tuple(arguments.exclude_clade_tax_id),
    )
    focus_authority = None
    group_rows: tuple[dict[str, Any], ...] = ()
    if not arguments.coverage_all_groups:
        focus_authority = read_focus_proteins(
            path=(
                arguments.focus_proteins
                if arguments.focus_proteins is not None
                else bundled_focus_path()
            )
        )
    if not arguments.validate_only:
        group_ids: tuple[str, ...] = ()
        if focus_authority is not None:
            focus_result = service.search_focus_clusters(
                filters=FocusClusterFilters(
                    protein_identifiers=focus_authority.identifiers,
                    group_type=arguments.coverage_group_type,
                    hierarchy_node=arguments.coverage_hierarchy_node,
                    maximum_rows=arguments.coverage_max_groups,
                )
            )
            if focus_result.truncated:
                raise InputValidationError(
                    f"Focus authority matches {focus_result.total_rows:,} groups, exceeding "
                    f"--coverage-max-groups={arguments.coverage_max_groups:,}."
                )
            group_ids = tuple(str(row["group_id"]) for row in focus_result.rows)
        if group_ids or focus_authority is None:
            group_rows = service.get_group_species_collection(
                group_type=arguments.coverage_group_type,
                hierarchy_node=arguments.coverage_hierarchy_node,
                group_ids=group_ids,
            )
    run = build_selection_coverage_run(
        graph=graph,
        expected=expected,
        selection=selection,
        group_species_rows=group_rows,
        run_id=resource.run_id,
        resource_identity=_resource_identity(resource_path=resource.resource_path, service=service),
        package_version=__version__,
        compact=arguments.coverage_compact,
        maximum_groups=arguments.coverage_max_groups,
        maximum_nodes=arguments.coverage_max_nodes,
        focus_authority_name=(
            focus_authority.source_name if focus_authority is not None else ""
        ),
        focus_authority_sha256=(
            focus_authority.sha256 if focus_authority is not None else ""
        ),
        focus_limited=focus_authority is not None,
        created_at_utc=arguments.coverage_created_at_utc,
    )
    if arguments.validate_only:
        _LOGGER.info(
            "Selection coverage validation passed: nodes=%s, predicates=%s, summary=%s",
            len(run.tree.nodes),
            len(selection.predicates),
            " ".join(selection_summary(selection=selection)),
        )
        return {
            "status": "validated",
            "selection_manifest_id": run.selection_manifest_id,
            "summary": run.summary(),
            "predicates": selection_summary(selection=selection),
        }
    if arguments.output_dir is None and not arguments.dry_run:
        raise InputValidationError(
            "--output-dir is required for coverage-tree publication."
        )
    target = (
        arguments.output_dir
        if arguments.output_dir is not None
        else resource.resource_path.parent / "coverage_tree_dry_run"
    )
    _validate_output_outside_resource(
        output_dir=target,
        resource_path=resource.resource_path,
    )
    return publish_selection_coverage_run(
        run=run,
        output_dir=target,
        dry_run=arguments.dry_run,
    )


def _resource_identity(
    *, resource_path: Path, service: OrthoFinderQueryService
) -> str:
    """Return a run-manifest checksum or explicit DuckDB identity."""

    root = resource_path if resource_path.is_dir() else resource_path.parent.parent
    manifest = root / "run_manifest.json"
    if manifest.is_file():
        return "run_manifest_sha256:" + hashlib.sha256(manifest.read_bytes()).hexdigest()
    stat = service.resource.database_path.stat()
    return (
        f"duckdb:{service.resource.run_id}:size={stat.st_size}:"
        f"schema={service.resource.schema_version}"
    )


def _validate_output_outside_resource(*, output_dir: Path, resource_path: Path) -> None:
    """Reject output paths that could mutate the completed read-only resource."""

    output = Path(output_dir).expanduser().resolve()
    resource = Path(resource_path).expanduser().resolve()
    resource_root = resource if resource.is_dir() else resource.parent.parent
    if output == resource_root or resource_root in output.parents:
        raise InputValidationError(
            "Selection coverage output must remain outside the immutable resource."
        )

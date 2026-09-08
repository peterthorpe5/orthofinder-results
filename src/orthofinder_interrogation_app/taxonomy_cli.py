"""Command-line creation of review-required taxonomy candidate sidecars."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

from orthofinder_results.errors import InputValidationError, OrthoFinderResultsError
from orthofinder_results.io_utils import (
    atomic_write_text,
    configure_logging,
    validate_persistent_path,
)

from .queries import OrthoFinderQueryService
from .resource import open_resource
from .taxonomy_reference import build_taxonomy_candidates

_LOGGER = logging.getLogger("orthofinder_interrogation_app.taxonomy_cli")


def build_parser() -> argparse.ArgumentParser:
    """Build the all-named-option taxonomy candidate parser."""

    parser = argparse.ArgumentParser(
        prog="orthofinder-taxonomy-map",
        description=(
            "Create exact-name NCBI taxonomy candidates for every species in a completed "
            "resource. Candidates remain pending until human review."
        ),
    )
    parser.add_argument(
        "--resource-dir",
        required=True,
        type=Path,
        help="Completed resource directory or its orthofinder_results.duckdb file.",
    )
    parser.add_argument(
        "--taxdump-dir",
        required=True,
        type=Path,
        help="Extracted NCBI taxdump directory containing names.dmp and nodes.dmp.",
    )
    parser.add_argument("--source-date", required=True, help="Taxdump date as YYYY-MM-DD.")
    parser.add_argument(
        "--source-version",
        required=True,
        help="Exact downloaded taxdump release or checksum label.",
    )
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--log-file", type=Path)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing sidecar TSV after all validation succeeds.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Create one generic review-required taxonomy candidate mapping.

    Args:
        argv: Optional named command-line arguments.

    Returns:
        Zero on success or two for a controlled validation failure.
    """

    parser = build_parser()
    arguments = parser.parse_args(argv)
    configure_logging(log_path=arguments.log_file, verbose=arguments.verbose)
    try:
        resource = open_resource(path=arguments.resource_dir)
        output = validate_persistent_path(
            path=arguments.output_tsv,
            role="Taxonomy mapping output",
        )
        _validate_sidecar_location(output=output, resource_path=resource.resource_path)
        if output.exists() and not arguments.force:
            parser.error(f"--output-tsv already exists; use --force to replace it: {output}")
        service = OrthoFinderQueryService(resource=resource)
        content = build_taxonomy_candidates(
            species=service.list_species(),
            taxdump_dir=arguments.taxdump_dir,
            source_date=arguments.source_date,
            source_version=arguments.source_version,
        )
        atomic_write_text(path=output, text=content.decode("utf-8"))
    except OrthoFinderResultsError as error:
        _LOGGER.error("Taxonomy candidate mapping failed: %s", error)
        return 2
    _LOGGER.info(
        "Taxonomy candidate mapping written: run=%s, output=%s, bytes=%s",
        resource.run_id,
        output,
        len(content),
    )
    return 0


def _validate_sidecar_location(*, output: Path, resource_path: Path) -> None:
    """Reject mapping writes inside the immutable completed resource."""

    resource = Path(resource_path).expanduser().resolve()
    resource_root = resource if resource.is_dir() else resource.parent
    if output == resource_root or resource_root in output.parents:
        raise InputValidationError(
            "Taxonomy mapping output must remain outside the immutable resource."
        )


if __name__ == "__main__":
    raise SystemExit(main())

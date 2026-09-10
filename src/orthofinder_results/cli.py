"""Named command-line interface for OrthoFinder result interrogation."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Sequence

from . import __version__
from .errors import OrthoFinderResultsError
from .io_utils import atomic_write_text, configure_logging, validate_persistent_path
from .pipeline import inspect_results, regenerate_report, run_pipeline

_LOGGER = logging.getLogger("orthofinder_results.cli")


def build_parser() -> argparse.ArgumentParser:
    """Build the all-named-option command-line parser.

    Returns:
        Configured argument parser.
    """

    parser = argparse.ArgumentParser(
        prog="orthofinder-results",
        description=(
            "Version-aware interrogation of OrthoFinder 2 and 3 results with "
            "compressed TSV, Parquet, DuckDB and offline HTML publication."
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument(
        "--action",
        required=True,
        choices=(
            "inspect",
            "run",
            "e3-precursor",
            "dispersion-benchmark",
            "report",
            "coverage-tree",
        ),
        help=(
            "Read-only inspection, complete resource publication, report-only "
            "regeneration, E3-focus precursor publication, matched-background "
            "dispersion benchmarking, or selection-coverage export."
        ),
    )
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument(
        "--inspection-output",
        type=Path,
        help="Required persistent JSON output for --action inspect.",
    )
    parser.add_argument("--output-dir", type=Path, help="Formal output for --action run.")
    parser.add_argument("--run-id", help="Immutable run identifier for --action run.")
    parser.add_argument(
        "--resource-dir",
        type=Path,
        help="Completed resource for --action report or --action coverage-tree.",
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        help="New standalone HTML file for --action report.",
    )
    parser.add_argument(
        "--log-output",
        type=Path,
        help="Persistent log file for --action report; defaults beside the HTML.",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        help=(
            "Optional staging root; Slurm uses node-local temporary storage by "
            "default, while direct runs default beside --output-dir."
        ),
    )
    parser.add_argument(
        "--alignment-dir",
        type=Path,
        help="Optional aligned FASTA directory; discovered OrthoFinder MSAs are used otherwise.",
    )
    parser.add_argument(
        "--distance-source",
        choices=("AUTO", "ALIGNED_SEQUENCE", "RESOLVED_GENE_TREE", "NONE"),
        default="AUTO",
        help="Distance authority; AUTO prefers alignments then resolved gene trees.",
    )
    parser.add_argument(
        "--distance-group-type",
        choices=("AUTO", "HOG", "LEGACY_ORTHOGROUP"),
        default="AUTO",
    )
    parser.add_argument("--distance-hierarchy-node", default="N0")
    parser.add_argument(
        "--distance-max-groups",
        type=int,
        default=0,
        help="Maximum groups for distance computation; zero means all eligible groups.",
    )
    parser.add_argument("--distance-max-members", type=int, default=250)
    parser.add_argument("--parse-gene-trees", action="store_true")
    parser.add_argument(
        "--report-max-statistic-rows",
        type=int,
        default=20000,
        help="Maximum embedded group rows (1-50,000; default 20,000).",
    )
    parser.add_argument(
        "--report-max-groups",
        type=int,
        default=25,
        help="Maximum interactive groups; combined distance-matrix budget is one million pairs.",
    )
    parser.add_argument(
        "--report-max-members",
        type=int,
        default=250,
        help="Maximum members per interactive group; default 250.",
    )
    parser.add_argument("--report-nearest-neighbours", type=int, default=3)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--keep-failed-work",
        action="store_true",
        help=(
            "Retain failed staging/copy directories for diagnosis; by default "
            "tracebacks are logged and partial data are removed."
        ),
    )
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--taxonomy-map",
        type=Path,
        help="Reviewed offline taxonomy mapping for --action coverage-tree.",
    )
    parser.add_argument(
        "--expected-taxa",
        type=Path,
        help="Optional explicit expected-taxon TSV for --action coverage-tree.",
    )
    parser.add_argument(
        "--focus-proteins",
        type=Path,
        help="Optional protein focus TSV; packaged E3 seeds are the default.",
    )
    parser.add_argument(
        "--benchmark-proteins",
        type=Path,
        help=(
            "Optional reviewed housekeeping/R marker TSV; the packaged "
            "Arabidopsis authority is the dispersion-benchmark default."
        ),
    )
    parser.add_argument(
        "--benchmark-controls-per-group",
        type=int,
        default=3,
        help="Unique distance-blind non-focus controls per benchmark target.",
    )
    parser.add_argument(
        "--benchmark-bootstrap-resamples",
        type=int,
        default=1_000,
        help="Deterministic bootstrap iterations for median-difference intervals.",
    )
    parser.add_argument("--require-exact-tax-id", action="append", default=[])
    parser.add_argument("--include-clade-tax-id", action="append", default=[])
    parser.add_argument("--only-in-clade-tax-id", action="append", default=[])
    parser.add_argument("--exclude-exact-tax-id", action="append", default=[])
    parser.add_argument("--exclude-clade-tax-id", action="append", default=[])
    parser.add_argument(
        "--coverage-group-type",
        choices=("HOG", "LEGACY_ORTHOGROUP"),
        default="HOG",
    )
    parser.add_argument("--coverage-hierarchy-node", default="N0")
    parser.add_argument(
        "--coverage-all-groups",
        action="store_true",
        help="Evaluate all selected-authority groups instead of E3 focus clusters.",
    )
    parser.add_argument("--coverage-compact", action="store_true")
    parser.add_argument("--coverage-max-groups", type=int, default=10_000)
    parser.add_argument("--coverage-max-nodes", type=int, default=5_000)
    parser.add_argument(
        "--coverage-created-at-utc",
        default="",
        help=(
            "Optional fixed ISO-8601 UTC manifest timestamp for byte-reproducible exports."
        ),
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Execute a named package action.

    Args:
        argv: Optional argument sequence excluding the program name.

    Returns:
        Process exit status.
    """

    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(verbose=args.verbose)
    try:
        if args.action == "inspect":
            _validate_inspect_arguments(parser=parser, args=args)
            record = inspect_results(results_dir=args.results_dir)
            atomic_write_text(
                path=args.inspection_output,
                text=json.dumps(record, indent=2, sort_keys=True) + "\n",
            )
            _LOGGER.info("Inspection written to %s", args.inspection_output.resolve())
            return 0
        if args.action == "report":
            _validate_report_arguments(parser=parser, args=args)
            log_output = validate_persistent_path(
                path=args.log_output or args.report_output.with_suffix(".log"),
                role="log_output",
            )
            configure_logging(log_path=log_output, verbose=args.verbose)
            record = regenerate_report(
                resource_dir=args.resource_dir,
                output_path=args.report_output,
                work_dir=args.work_dir,
                report_max_statistic_rows=args.report_max_statistic_rows,
                report_max_groups=args.report_max_groups,
                report_max_members=args.report_max_members,
                report_nearest_neighbours=args.report_nearest_neighbours,
                force=args.force,
            )
            _LOGGER.info(
                "Report regenerated: %s bytes at %s",
                f"{record['size_bytes']:,}",
                record["path"],
            )
            return 0
        if args.action == "coverage-tree":
            _validate_coverage_arguments(parser=parser, args=args)
            from orthofinder_interrogation_app.coverage_cli import (
                run_coverage_tree_action,
            )

            manifest = run_coverage_tree_action(arguments=args)
            print(json.dumps(manifest, indent=2, sort_keys=True))
            return 0
        if args.action == "e3-precursor":
            _validate_e3_precursor_arguments(parser=parser, args=args)
            from orthofinder_interrogation_app.focus import bundled_focus_path

            focus_path = (
                args.focus_proteins
                if args.focus_proteins is not None
                else bundled_focus_path()
            )
            manifest = run_pipeline(
                results_dir=args.results_dir,
                output_dir=args.output_dir,
                run_id=args.run_id,
                work_dir=args.work_dir,
                alignment_dir=None,
                distance_source="RESOLVED_GENE_TREE",
                distance_group_type="HOG",
                distance_hierarchy_node="N0",
                distance_max_groups=0,
                distance_max_members=args.distance_max_members,
                parse_gene_trees=True,
                report_max_statistic_rows=args.report_max_statistic_rows,
                report_max_groups=args.report_max_groups,
                report_max_members=args.report_max_members,
                report_nearest_neighbours=args.report_nearest_neighbours,
                resume=args.resume,
                force=args.force,
                keep_failed_work=args.keep_failed_work,
                verbose=args.verbose,
                focus_proteins_path=focus_path,
                focus_group_type="HOG",
                focus_hierarchy_node="N0",
            )
            _LOGGER.info(
                "Completed E3 precursor run %s with status %s",
                manifest["run_id"],
                manifest["status"],
            )
            return 0
        if args.action == "dispersion-benchmark":
            _validate_dispersion_benchmark_arguments(parser=parser, args=args)
            from orthofinder_interrogation_app.focus import bundled_focus_path

            from .benchmark_authority import bundled_benchmark_path

            focus_path = args.focus_proteins or bundled_focus_path()
            benchmark_path = args.benchmark_proteins or bundled_benchmark_path()
            manifest = run_pipeline(
                results_dir=args.results_dir,
                output_dir=args.output_dir,
                run_id=args.run_id,
                work_dir=args.work_dir,
                alignment_dir=None,
                distance_source="RESOLVED_GENE_TREE",
                distance_group_type="HOG",
                distance_hierarchy_node="N0",
                distance_max_groups=0,
                distance_max_members=args.distance_max_members,
                parse_gene_trees=True,
                report_max_statistic_rows=args.report_max_statistic_rows,
                report_max_groups=args.report_max_groups,
                report_max_members=args.report_max_members,
                report_nearest_neighbours=args.report_nearest_neighbours,
                resume=args.resume,
                force=args.force,
                keep_failed_work=args.keep_failed_work,
                verbose=args.verbose,
                focus_proteins_path=focus_path,
                focus_group_type="HOG",
                focus_hierarchy_node="N0",
                benchmark_proteins_path=benchmark_path,
                benchmark_controls_per_group=args.benchmark_controls_per_group,
                benchmark_bootstrap_resamples=(
                    args.benchmark_bootstrap_resamples
                ),
            )
            _LOGGER.info(
                "Completed dispersion benchmark run %s with status %s",
                manifest["run_id"],
                manifest["status"],
            )
            return 0
        _validate_run_arguments(parser=parser, args=args)
        manifest = run_pipeline(
            results_dir=args.results_dir,
            output_dir=args.output_dir,
            run_id=args.run_id,
            work_dir=args.work_dir,
            alignment_dir=args.alignment_dir,
            distance_source=args.distance_source,
            distance_group_type=args.distance_group_type,
            distance_hierarchy_node=args.distance_hierarchy_node,
            distance_max_groups=args.distance_max_groups,
            distance_max_members=args.distance_max_members,
            parse_gene_trees=args.parse_gene_trees,
            report_max_statistic_rows=args.report_max_statistic_rows,
            report_max_groups=args.report_max_groups,
            report_max_members=args.report_max_members,
            report_nearest_neighbours=args.report_nearest_neighbours,
            resume=args.resume,
            force=args.force,
            keep_failed_work=args.keep_failed_work,
            verbose=args.verbose,
        )
        _LOGGER.info("Completed run %s with status %s", manifest["run_id"], manifest["status"])
        return 0
    except OrthoFinderResultsError as error:
        _LOGGER.error("%s", error)
        return 2


def _validate_inspect_arguments(
    *, parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    """Validate arguments specific to layout inspection.

    Args:
        parser: Parser used for a controlled error.
        args: Parsed arguments.
    """

    if args.inspection_output is None:
        parser.error("--inspection-output is required for --action inspect.")
    if args.results_dir is None:
        parser.error("--results-dir is required for --action inspect.")
    if any(
        value is not None
        for value in (args.output_dir, args.run_id, args.resource_dir, args.report_output)
    ):
        parser.error("Run and report output arguments are not valid for --action inspect.")


def _validate_run_arguments(*, parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate arguments specific to resource publication.

    Args:
        parser: Parser used for a controlled error.
        args: Parsed arguments.
    """

    missing = [
        name
        for name, value in (
            ("--results-dir", args.results_dir),
            ("--output-dir", args.output_dir),
            ("--run-id", args.run_id),
        )
        if value is None
    ]
    if missing:
        parser.error(f"{' and '.join(missing)} required for --action run.")
    if args.inspection_output is not None:
        parser.error("--inspection-output is not valid for --action run.")
    if any(value is not None for value in (args.resource_dir, args.report_output, args.log_output)):
        parser.error("Report-only arguments are not valid for --action run.")
    if args.focus_proteins is not None:
        parser.error(
            "--focus-proteins requires --action e3-precursor, "
            "--action dispersion-benchmark or --action coverage-tree."
        )
    if args.benchmark_proteins is not None:
        parser.error(
            "--benchmark-proteins requires --action dispersion-benchmark."
        )


def _validate_e3_precursor_arguments(
    *, parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    """Validate the fixed HOG N0 E3 precursor publication contract.

    Args:
        parser: Parser used for controlled named-option errors.
        args: Parsed command-line arguments.
    """

    missing = [
        name
        for name, value in (
            ("--results-dir", args.results_dir),
            ("--output-dir", args.output_dir),
            ("--run-id", args.run_id),
        )
        if value is None
    ]
    if missing:
        parser.error(
            f"{' and '.join(missing)} required for --action e3-precursor."
        )
    if any(
        value is not None
        for value in (
            args.inspection_output,
            args.resource_dir,
            args.report_output,
            args.log_output,
            args.alignment_dir,
        )
    ):
        parser.error(
            "Inspection, report, resource and alignment arguments are not valid for "
            "--action e3-precursor."
        )
    if args.distance_source not in {"AUTO", "RESOLVED_GENE_TREE"}:
        parser.error(
            "--action e3-precursor uses RESOLVED_GENE_TREE distances."
        )
    if args.distance_group_type not in {"AUTO", "HOG"}:
        parser.error("--action e3-precursor uses the HOG group collection.")
    if args.distance_hierarchy_node != "N0":
        parser.error("--action e3-precursor uses the N0 HOG hierarchy.")
    if args.distance_max_groups != 0:
        parser.error(
            "--action e3-precursor requires --distance-max-groups 0 so matching "
            "clusters are never truncated."
        )
    if args.benchmark_proteins is not None:
        parser.error(
            "--benchmark-proteins requires --action dispersion-benchmark."
        )


def _validate_dispersion_benchmark_arguments(
    *, parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    """Validate the fixed HOG-N0 matched-background benchmark contract.

    Args:
        parser: Parser used for controlled named-option errors.
        args: Parsed command-line arguments.
    """

    missing = [
        name
        for name, value in (
            ("--results-dir", args.results_dir),
            ("--output-dir", args.output_dir),
            ("--run-id", args.run_id),
        )
        if value is None
    ]
    if missing:
        parser.error(
            f"{' and '.join(missing)} required for --action dispersion-benchmark."
        )
    if any(
        value is not None
        for value in (
            args.inspection_output,
            args.resource_dir,
            args.report_output,
            args.log_output,
            args.alignment_dir,
        )
    ):
        parser.error(
            "Inspection, report, resource and alignment arguments are not valid "
            "for --action dispersion-benchmark."
        )
    if args.distance_source not in {"AUTO", "RESOLVED_GENE_TREE"}:
        parser.error(
            "--action dispersion-benchmark uses RESOLVED_GENE_TREE distances."
        )
    if args.distance_group_type not in {"AUTO", "HOG"}:
        parser.error("--action dispersion-benchmark uses the HOG collection.")
    if args.distance_hierarchy_node != "N0":
        parser.error("--action dispersion-benchmark uses the N0 HOG hierarchy.")
    if args.distance_max_groups != 0:
        parser.error(
            "--action dispersion-benchmark requires --distance-max-groups 0."
        )
    if not 1 <= args.benchmark_controls_per_group <= 50:
        parser.error("--benchmark-controls-per-group must be between 1 and 50.")
    if not 100 <= args.benchmark_bootstrap_resamples <= 100_000:
        parser.error(
            "--benchmark-bootstrap-resamples must be between 100 and 100,000."
        )


def _validate_report_arguments(
    *, parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    """Validate arguments specific to standalone report regeneration.

    Args:
        parser: Parser used for a controlled error.
        args: Parsed arguments.
    """

    missing = [
        name
        for name, value in (
            ("--resource-dir", args.resource_dir),
            ("--report-output", args.report_output),
        )
        if value is None
    ]
    if missing:
        parser.error(f"{' and '.join(missing)} required for --action report.")
    if any(
        value is not None
        for value in (args.results_dir, args.inspection_output, args.output_dir, args.run_id)
    ):
        parser.error("Inspection and run arguments are not valid for --action report.")
    resource = args.resource_dir.expanduser().resolve()
    report_paths = [args.report_output]
    if args.log_output is not None:
        report_paths.append(args.log_output)
    if any(
        resource == path.expanduser().resolve()
        or resource in path.expanduser().resolve().parents
        for path in report_paths
    ):
        parser.error(
            "--report-output and --log-output must be outside the immutable --resource-dir."
        )


def _validate_coverage_arguments(
    *, parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    """Validate arguments specific to selection coverage publication.

    Args:
        parser: Parser used for controlled named-option errors.
        args: Parsed command-line arguments.
    """

    if args.resource_dir is None or args.taxonomy_map is None:
        parser.error(
            "--resource-dir and --taxonomy-map are required for --action coverage-tree."
        )
    if args.validate_only and args.dry_run:
        parser.error("--validate-only and --dry-run are mutually exclusive.")
    if not 1 <= args.coverage_max_groups <= 250_000:
        parser.error("--coverage-max-groups must be between 1 and 250,000.")
    if not 1 <= args.coverage_max_nodes <= 5_000:
        parser.error("--coverage-max-nodes must be between 1 and 5,000.")
    if args.coverage_group_type == "LEGACY_ORTHOGROUP":
        if args.coverage_hierarchy_node in {"ROOT", "root"}:
            args.coverage_hierarchy_node = ""
        elif args.coverage_hierarchy_node:
            parser.error(
                "--coverage-hierarchy-node must be ROOT for LEGACY_ORTHOGROUP."
            )
    if any(
        value is not None
        for value in (
            args.results_dir,
            args.inspection_output,
            args.report_output,
            args.log_output,
        )
    ):
        parser.error(
            "Inspection, raw-run and report-only paths are not valid for coverage-tree."
        )

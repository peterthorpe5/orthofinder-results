"""Process-isolated publication of complete OrthoFinder memberships."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import resource
import subprocess
import sys
import tempfile
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Any, TextIO

from .errors import InputValidationError, OrthoFinderResultsError, PublicationError
from .io_utils import atomic_write_json, open_text
from .parsers import MEMBERSHIP_FIELDS, iter_memberships
from .statistics import (
    GROUP_SPECIES_STATISTIC_FIELDS,
    GROUP_STATISTIC_FIELDS,
    GroupAccumulator,
)

_LOGGER = logging.getLogger("orthofinder_results.membership_publication")
_FRAGMENT_CHUNK_CHARACTERS = 1024 * 1024
_METADATA_VERSION = 1


def publish_membership_tables(
    *,
    tables_dir: Path,
    legacy_sources: Sequence[tuple[Path, str, str]],
    hog_sources: Sequence[tuple[Path, str, str]],
    run_id: str,
    worker_python: Path | None = None,
) -> tuple[dict[str, int], int, int, set[str]]:
    """Publish membership authorities with one fresh process per source file.

    Args:
        tables_dir: Staging directory for compressed analytical tables.
        legacy_sources: Legacy orthogroup source specifications.
        hog_sources: Hierarchical orthogroup source specifications.
        run_id: Immutable resource identifier.
        worker_python: Optional interpreter override used by tests.

    Returns:
        Membership counts, group count, group-species count and species labels.

    Raises:
        InputValidationError: If the staging directory or source specification is invalid.
        PublicationError: If an isolated worker or fragment validation fails.
    """

    destination = Path(tables_dir).expanduser().resolve()
    if not destination.is_dir():
        raise InputValidationError(
            f"Membership table directory does not exist: {destination}"
        )
    python = (
        Path(worker_python).expanduser().absolute()
        if worker_python is not None
        else Path(sys.executable).absolute()
    )
    if not python.is_file():
        raise InputValidationError(f"Membership worker interpreter is unavailable: {python}")

    statistics_path = destination / "group_statistics.tsv.gz"
    species_statistics_path = destination / "group_species_statistics.tsv.gz"
    counts = {
        "legacy_orthogroup_membership_count": 0,
        "hog_membership_count": 0,
    }
    group_count = 0
    group_species_count = 0
    species: set[str] = set()
    with (
        open_text(path=statistics_path, mode="w") as statistics_handle,
        open_text(path=species_statistics_path, mode="w") as species_statistics_handle,
    ):
        _write_header(handle=statistics_handle, fields=GROUP_STATISTIC_FIELDS)
        _write_header(
            handle=species_statistics_handle,
            fields=GROUP_SPECIES_STATISTIC_FIELDS,
        )
        legacy_result = _write_membership_authority(
            path=destination / "legacy_orthogroup_memberships.tsv.gz",
            sources=legacy_sources,
            run_id=run_id,
            statistics_handle=statistics_handle,
            species_statistics_handle=species_statistics_handle,
            worker_python=python,
        )
        counts["legacy_orthogroup_membership_count"] = legacy_result[0]
        group_count += legacy_result[1]
        group_species_count += legacy_result[2]
        species.update(legacy_result[3])

        hog_result = _write_membership_authority(
            path=destination / "hog_memberships.tsv.gz",
            sources=hog_sources,
            run_id=run_id,
            statistics_handle=statistics_handle,
            species_statistics_handle=species_statistics_handle,
            worker_python=python,
        )
        counts["hog_membership_count"] = hog_result[0]
        group_count += hog_result[1]
        group_species_count += hog_result[2]
        species.update(hog_result[3])

    _LOGGER.info(
        "Published %s memberships across %s run-scoped groups with isolated workers.",
        f"{sum(counts.values()):,}",
        f"{group_count:,}",
    )
    return counts, group_count, group_species_count, species


def _write_membership_authority(
    *,
    path: Path,
    sources: Sequence[tuple[Path, str, str]],
    run_id: str,
    statistics_handle: TextIO,
    species_statistics_handle: TextIO,
    worker_python: Path,
) -> tuple[int, int, int, set[str]]:
    """Merge validated source fragments into one compressed authority."""

    member_count = 0
    group_count = 0
    group_species_count = 0
    species: set[str] = set()
    with open_text(path=path, mode="w") as membership_handle:
        _write_header(handle=membership_handle, fields=MEMBERSHIP_FIELDS)
        for source, group_type, hierarchy_node in sources:
            started = time.perf_counter()
            source_path = Path(source).expanduser().resolve()
            _LOGGER.info(
                "Membership source started: type=%s, node=%s, file=%s, "
                "parent_peak_rss_mib=%.1f",
                group_type,
                hierarchy_node or "ROOT",
                source_path,
                peak_rss_mib(),
            )
            with tempfile.TemporaryDirectory(
                prefix=".membership_source.",
                dir=path.parent.parent,
            ) as temporary_name:
                temporary = Path(temporary_name)
                membership_fragment = temporary / "memberships.tsv"
                statistics_fragment = temporary / "group_statistics.tsv"
                species_statistics_fragment = temporary / "group_species_statistics.tsv"
                metadata_path = temporary / "metadata.json"
                _run_worker(
                    worker_python=worker_python,
                    source=source_path,
                    run_id=run_id,
                    group_type=group_type,
                    hierarchy_node=hierarchy_node,
                    membership_output=membership_fragment,
                    statistics_output=statistics_fragment,
                    species_statistics_output=species_statistics_fragment,
                    metadata_output=metadata_path,
                )
                metadata = _read_worker_metadata(
                    path=metadata_path,
                    source=source_path,
                    group_type=group_type,
                    hierarchy_node=hierarchy_node,
                )
                _append_fragment(
                    source=membership_fragment,
                    target=membership_handle,
                    expected_rows=metadata["member_count"],
                    role="membership",
                )
                _append_fragment(
                    source=statistics_fragment,
                    target=statistics_handle,
                    expected_rows=metadata["group_count"],
                    role="group-statistics",
                )
                _append_fragment(
                    source=species_statistics_fragment,
                    target=species_statistics_handle,
                    expected_rows=metadata["group_species_count"],
                    role="group-species-statistics",
                )
            member_count += metadata["member_count"]
            group_count += metadata["group_count"]
            group_species_count += metadata["group_species_count"]
            species.update(metadata["species"])
            _LOGGER.info(
                "Membership source finished: type=%s, node=%s, rows=%s, "
                "groups=%s, worker_pid=%s, worker_peak_rss_mib=%.1f, "
                "parent_peak_rss_mib=%.1f, elapsed_seconds=%.3f",
                group_type,
                hierarchy_node or "ROOT",
                f"{metadata['member_count']:,}",
                f"{metadata['group_count']:,}",
                metadata["worker_pid"],
                metadata["peak_rss_mib"],
                peak_rss_mib(),
                time.perf_counter() - started,
            )
    return member_count, group_count, group_species_count, species


def _run_worker(
    *,
    worker_python: Path,
    source: Path,
    run_id: str,
    group_type: str,
    hierarchy_node: str,
    membership_output: Path,
    statistics_output: Path,
    species_statistics_output: Path,
    metadata_output: Path,
) -> None:
    """Execute one lean membership worker and diagnose abnormal termination."""

    environment = dict(os.environ)
    environment["MALLOC_ARENA_MAX"] = "2"
    command = (
        str(worker_python),
        "-m",
        "orthofinder_results.membership_publication",
        "--source",
        str(source),
        "--run-id",
        run_id,
        "--group-type",
        group_type,
        "--hierarchy-node",
        hierarchy_node,
        "--membership-output",
        str(membership_output),
        "--statistics-output",
        str(statistics_output),
        "--species-statistics-output",
        str(species_statistics_output),
        "--metadata-output",
        str(metadata_output),
    )
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    if completed.returncode == 0:
        return
    detail = completed.stderr.strip() or completed.stdout.strip() or "no worker output"
    if completed.returncode < 0:
        status = f"terminated by signal {-completed.returncode}"
    else:
        status = f"exited with status {completed.returncode}"
    raise PublicationError(
        f"Membership worker {status} for {group_type} {hierarchy_node or 'ROOT'} "
        f"from {source}: {detail[-4000:]}"
    )


def _read_worker_metadata(
    *,
    path: Path,
    source: Path,
    group_type: str,
    hierarchy_node: str,
) -> dict[str, Any]:
    """Read and validate the small completion record from one worker."""

    if not path.is_file():
        raise PublicationError(f"Membership worker did not publish metadata: {path}")
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PublicationError(f"Invalid membership worker metadata: {path}") from error
    expected_scalars = {
        "metadata_version": _METADATA_VERSION,
        "source": str(source),
        "group_type": group_type,
        "hierarchy_node": hierarchy_node,
    }
    for field, expected in expected_scalars.items():
        if record.get(field) != expected:
            raise PublicationError(
                f"Membership worker metadata field {field!r} was "
                f"{record.get(field)!r}; expected {expected!r}."
            )
    for field in ("member_count", "group_count", "group_species_count"):
        value = record.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise PublicationError(
                f"Membership worker metadata has invalid {field}: {value!r}"
            )
    worker_pid = record.get("worker_pid")
    if (
        isinstance(worker_pid, bool)
        or not isinstance(worker_pid, int)
        or worker_pid < 1
    ):
        raise PublicationError(
            f"Membership worker metadata has invalid worker_pid: {worker_pid!r}"
        )
    peak = record.get("peak_rss_mib")
    if isinstance(peak, bool) or not isinstance(peak, (int, float)) or peak < 0:
        raise PublicationError(
            f"Membership worker metadata has invalid peak_rss_mib: {peak!r}"
        )
    labels = record.get("species")
    if not isinstance(labels, list) or any(
        not isinstance(label, str) or not label for label in labels
    ):
        raise PublicationError("Membership worker metadata has invalid species labels.")
    if len(labels) != len(set(labels)) or labels != sorted(labels):
        raise PublicationError(
            "Membership worker species labels must be unique and sorted."
        )
    return record


def _append_fragment(
    *,
    source: Path,
    target: TextIO,
    expected_rows: int,
    role: str,
) -> None:
    """Append one bounded text fragment while validating its row count."""

    if not source.is_file():
        raise PublicationError(f"Membership worker lacks {role} fragment: {source}")
    row_count = 0
    with source.open(mode="r", encoding="utf-8", newline="") as handle:
        while chunk := handle.read(_FRAGMENT_CHUNK_CHARACTERS):
            row_count += chunk.count("\n")
            target.write(chunk)
    if row_count != expected_rows:
        raise PublicationError(
            f"Membership worker {role} fragment has {row_count:,} rows; "
            f"expected {expected_rows:,}: {source}"
        )


def _write_header(*, handle: TextIO, fields: Sequence[str]) -> None:
    """Write a deterministic TSV heading without allocating row objects."""

    if not fields or any(not field or "\t" in field or "\n" in field for field in fields):
        raise ValueError("TSV heading fields must be non-empty single-line values.")
    handle.write("\t".join(fields) + "\n")


def publish_membership_source(
    *,
    source: Path,
    run_id: str,
    group_type: str,
    hierarchy_node: str,
    membership_output: Path,
    statistics_output: Path,
    species_statistics_output: Path,
    metadata_output: Path,
) -> dict[str, Any]:
    """Publish headerless fragments for exactly one OrthoFinder source table."""

    outputs = tuple(
        Path(path).expanduser().resolve()
        for path in (
            membership_output,
            statistics_output,
            species_statistics_output,
            metadata_output,
        )
    )
    if len(set(outputs)) != len(outputs):
        raise InputValidationError("Membership worker output paths must be distinct.")
    if any(path.exists() for path in outputs):
        raise InputValidationError("Membership worker refuses to replace an output file.")
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)

    resolved_source = Path(source).expanduser().resolve()
    member_count = 0
    group_count = 0
    group_species_count = 0
    species: set[str] = set()
    with (
        outputs[0].open(mode="w", encoding="utf-8", newline="") as membership_handle,
        outputs[1].open(mode="w", encoding="utf-8", newline="") as statistics_handle,
        outputs[2].open(
            mode="w", encoding="utf-8", newline=""
        ) as species_statistics_handle,
    ):
        membership_writer = csv.DictWriter(
            membership_handle,
            fieldnames=MEMBERSHIP_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        statistics_writer = csv.DictWriter(
            statistics_handle,
            fieldnames=GROUP_STATISTIC_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        species_statistics_writer = csv.DictWriter(
            species_statistics_handle,
            fieldnames=GROUP_SPECIES_STATISTIC_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        current_key: tuple[str, str, str] | None = None
        accumulator: GroupAccumulator | None = None
        for row in iter_memberships(
            path=resolved_source,
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
            species_label = str(row["species_label"])
            accumulator.add_member(species_label=species_label)
            species.add(species_label)
            membership_writer.writerow(row)
            member_count += 1
        if accumulator is not None:
            _write_accumulator_statistics(
                accumulator=accumulator,
                statistics_writer=statistics_writer,
                species_statistics_writer=species_statistics_writer,
            )
            group_count += 1
            group_species_count += len(accumulator.species_counts)

    metadata = {
        "metadata_version": _METADATA_VERSION,
        "source": str(resolved_source),
        "group_type": group_type,
        "hierarchy_node": hierarchy_node,
        "member_count": member_count,
        "group_count": group_count,
        "group_species_count": group_species_count,
        "species": sorted(species),
        "peak_rss_mib": peak_rss_mib(),
        "worker_pid": os.getpid(),
    }
    atomic_write_json(path=outputs[3], record=metadata)
    return metadata


def _write_accumulator_statistics(
    *,
    accumulator: GroupAccumulator,
    statistics_writer: csv.DictWriter,
    species_statistics_writer: csv.DictWriter,
) -> None:
    """Write one group's summary and represented-species rows."""

    statistics_writer.writerow(accumulator.to_record())
    species_statistics_writer.writerows(accumulator.to_species_records())


def peak_rss_mib() -> float:
    """Return the current process's peak resident set size in MiB."""

    maximum = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    divisor = 1024.0 * 1024.0 if sys.platform == "darwin" else 1024.0
    return maximum / divisor


def build_parser() -> argparse.ArgumentParser:
    """Return the named-option parser for one private membership worker."""

    parser = argparse.ArgumentParser(
        prog="python -m orthofinder_results.membership_publication",
        description="Publish one process-isolated OrthoFinder membership fragment.",
    )
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument(
        "--group-type",
        required=True,
        choices=("LEGACY_ORTHOGROUP", "HOG"),
    )
    parser.add_argument("--hierarchy-node", required=True)
    parser.add_argument("--membership-output", required=True, type=Path)
    parser.add_argument("--statistics-output", required=True, type=Path)
    parser.add_argument("--species-statistics-output", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run one membership worker and return its process status."""

    arguments = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )
    try:
        publish_membership_source(
            source=arguments.source,
            run_id=arguments.run_id,
            group_type=arguments.group_type,
            hierarchy_node=arguments.hierarchy_node,
            membership_output=arguments.membership_output,
            statistics_output=arguments.statistics_output,
            species_statistics_output=arguments.species_statistics_output,
            metadata_output=arguments.metadata_output,
        )
    except (OrthoFinderResultsError, OSError, UnicodeError, ValueError) as error:
        _LOGGER.error("Membership worker failed: %s", error)
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

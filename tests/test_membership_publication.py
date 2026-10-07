"""Tests for process-isolated, bounded-memory membership publication."""

from __future__ import annotations

import csv
import io
import json
import logging
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

import orthofinder_results.membership_publication as membership_publication
from orthofinder_results.errors import InputValidationError, PublicationError
from orthofinder_results.io_utils import read_tsv
from orthofinder_results.membership_publication import (
    _append_fragment,
    _read_worker_metadata,
    _run_worker,
    _write_header,
    main,
    publish_membership_source,
    publish_membership_tables,
)


def _write_hog(*, path: Path, node: str, member_suffix: str = "") -> None:
    """Write one small two-group hierarchical orthogroup authority."""

    path.write_text(
        "HOG\tOG\tGene Tree Parent Clade\tSpecies_A\tSpecies_B\n"
        f"{node}.HOG0000001\tOG0000001\t{node}\tprotA{member_suffix}\tprotB{member_suffix}\n"
        f"{node}.HOG0000002\tOG0000002\t{node}\tprotA2{member_suffix}\t\n",
        encoding="utf-8",
    )


def test_membership_source_publishes_headerless_validated_fragments(
    tmp_path: Path,
) -> None:
    """One worker retains exact rows and compact per-group statistics."""

    source = tmp_path / "N0.tsv"
    _write_hog(path=source, node="N0")
    membership = tmp_path / "memberships.tsv"
    statistics = tmp_path / "statistics.tsv"
    species_statistics = tmp_path / "species_statistics.tsv"
    metadata = tmp_path / "metadata.json"
    record = publish_membership_source(
        source=source,
        run_id="test-run",
        group_type="HOG",
        hierarchy_node="N0",
        membership_output=membership,
        statistics_output=statistics,
        species_statistics_output=species_statistics,
        metadata_output=metadata,
    )
    assert record["member_count"] == 3
    assert record["group_count"] == 2
    assert record["group_species_count"] == 3
    assert record["species"] == ["Species_A", "Species_B"]
    assert record["peak_rss_mib"] > 0
    assert record["worker_pid"] > 0
    assert metadata.is_file()
    assert len(membership.read_text(encoding="utf-8").splitlines()) == 3
    assert len(statistics.read_text(encoding="utf-8").splitlines()) == 2
    assert len(species_statistics.read_text(encoding="utf-8").splitlines()) == 3
    assert not membership.read_text(encoding="utf-8").startswith("run_id\t")


def test_membership_tables_use_isolated_workers_and_merge_exact_rows(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Independent source workers merge into the unchanged compressed authority."""

    tables = tmp_path / "staging" / "tables"
    tables.mkdir(parents=True)
    legacy = tmp_path / "Orthogroups.tsv"
    legacy.write_text(
        "Orthogroup\tSpecies_A\tSpecies_B\n"
        "OG0000001\tlegacyA\tlegacyB\n",
        encoding="utf-8",
    )
    first = tmp_path / "N0.tsv"
    second = tmp_path / "N1.tsv"
    _write_hog(path=first, node="N0")
    _write_hog(path=second, node="N1", member_suffix="x")
    with caplog.at_level(logging.INFO):
        counts, group_count, group_species_count, species = publish_membership_tables(
            tables_dir=tables,
            legacy_sources=((legacy, "LEGACY_ORTHOGROUP", ""),),
            hog_sources=(
                (first, "HOG", "N0"),
                (second, "HOG", "N1"),
            ),
            run_id="isolated-run",
        )
    assert counts == {
        "legacy_orthogroup_membership_count": 2,
        "hog_membership_count": 6,
    }
    assert group_count == 5
    assert group_species_count == 8
    assert species == {"Species_A", "Species_B"}
    legacy_rows = list(read_tsv(path=tables / "legacy_orthogroup_memberships.tsv.gz"))
    hog_rows = list(read_tsv(path=tables / "hog_memberships.tsv.gz"))
    group_rows = list(read_tsv(path=tables / "group_statistics.tsv.gz"))
    species_rows = list(read_tsv(path=tables / "group_species_statistics.tsv.gz"))
    assert len(legacy_rows) == 2
    assert len(hog_rows) == 6
    assert len(group_rows) == 5
    assert len(species_rows) == 8
    assert [row["hierarchy_node"] for row in hog_rows] == [
        "N0",
        "N0",
        "N0",
        "N1",
        "N1",
        "N1",
    ]
    worker_pids = re.findall(r"worker_pid=(\d+)", caplog.text)
    assert len(worker_pids) == 3
    assert len(set(worker_pids)) == 3
    assert not tuple((tmp_path / "staging").glob(".membership_source.*"))


def test_worker_main_and_output_safety_fail_closed(tmp_path: Path) -> None:
    """The private worker returns controlled errors and never replaces fragments."""

    outputs = {
        "membership_output": tmp_path / "members.tsv",
        "statistics_output": tmp_path / "statistics.tsv",
        "species_statistics_output": tmp_path / "species.tsv",
        "metadata_output": tmp_path / "metadata.json",
    }
    arguments = [
        "--source",
        str(tmp_path / "missing.tsv"),
        "--run-id",
        "run",
        "--group-type",
        "HOG",
        "--hierarchy-node",
        "N0",
        "--membership-output",
        str(outputs["membership_output"]),
        "--statistics-output",
        str(outputs["statistics_output"]),
        "--species-statistics-output",
        str(outputs["species_statistics_output"]),
        "--metadata-output",
        str(outputs["metadata_output"]),
    ]
    assert main(arguments) == 2

    source = tmp_path / "N0.tsv"
    _write_hog(path=source, node="N0")
    with pytest.raises(InputValidationError, match="distinct"):
        publish_membership_source(
            source=source,
            run_id="run",
            group_type="HOG",
            hierarchy_node="N0",
            membership_output=outputs["membership_output"],
            statistics_output=outputs["membership_output"],
            species_statistics_output=outputs["species_statistics_output"],
            metadata_output=outputs["metadata_output"],
        )
    outputs["membership_output"].write_text("retain\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="refuses to replace"):
        publish_membership_source(
            source=source,
            run_id="run",
            group_type="HOG",
            hierarchy_node="N0",
            **outputs,
        )


def test_worker_failures_and_fragment_damage_are_diagnostic(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Signals, invalid metadata and incomplete fragments remain explicit."""

    observed_environment: dict[str, str] = {}

    def fail_worker(*args: object, **kwargs: object) -> SimpleNamespace:
        environment = kwargs["env"]
        assert isinstance(environment, dict)
        observed_environment.update(environment)
        return SimpleNamespace(
            returncode=-11,
            stderr="native failure",
            stdout="",
        )

    monkeypatch.setenv("MALLOC_ARENA_MAX", "64")
    monkeypatch.setattr(membership_publication.subprocess, "run", fail_worker)
    with pytest.raises(PublicationError, match="signal 11"):
        _run_worker(
            worker_python=Path(sys.executable),
            source=tmp_path / "N0.tsv",
            run_id="run",
            group_type="HOG",
            hierarchy_node="N0",
            membership_output=tmp_path / "members.tsv",
            statistics_output=tmp_path / "statistics.tsv",
            species_statistics_output=tmp_path / "species.tsv",
            metadata_output=tmp_path / "metadata.json",
        )
    assert observed_environment["MALLOC_ARENA_MAX"] == "2"

    fragment = tmp_path / "fragment.tsv"
    fragment.write_text("one\ntwo\n", encoding="utf-8")
    with pytest.raises(PublicationError, match="has 2 rows; expected 3"):
        _append_fragment(
            source=fragment,
            target=io.StringIO(),
            expected_rows=3,
            role="test",
        )
    with pytest.raises(PublicationError, match="lacks test fragment"):
        _append_fragment(
            source=tmp_path / "absent.tsv",
            target=io.StringIO(),
            expected_rows=0,
            role="test",
        )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("metadata_version", 2, "metadata_version"),
        ("member_count", -1, "member_count"),
        ("worker_pid", 0, "worker_pid"),
        ("peak_rss_mib", -1.0, "peak_rss_mib"),
        ("species", ["Species_B", "Species_A"], "unique and sorted"),
        ("species", [""], "invalid species"),
    ],
)
def test_worker_metadata_validation(
    tmp_path: Path,
    field: str,
    value: object,
    message: str,
) -> None:
    """Every worker completion field is validated before fragments are merged."""

    source = (tmp_path / "N0.tsv").resolve()
    record = {
        "metadata_version": 1,
        "source": str(source),
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "member_count": 1,
        "group_count": 1,
        "group_species_count": 1,
        "species": ["Species_A"],
        "peak_rss_mib": 1.0,
        "worker_pid": 1,
    }
    record[field] = value
    metadata = tmp_path / f"metadata_{field}.json"
    metadata.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(PublicationError, match=message):
        _read_worker_metadata(
            path=metadata,
            source=source,
            group_type="HOG",
            hierarchy_node="N0",
        )


def test_metadata_headers_interpreter_and_platform_branches(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Missing metadata, unsafe headings and unavailable interpreters fail early."""

    with pytest.raises(PublicationError, match="did not publish metadata"):
        _read_worker_metadata(
            path=tmp_path / "missing.json",
            source=tmp_path / "N0.tsv",
            group_type="HOG",
            hierarchy_node="N0",
        )
    malformed = tmp_path / "malformed.json"
    malformed.write_text("{", encoding="utf-8")
    with pytest.raises(PublicationError, match="Invalid membership worker metadata"):
        _read_worker_metadata(
            path=malformed,
            source=tmp_path / "N0.tsv",
            group_type="HOG",
            hierarchy_node="N0",
        )
    with pytest.raises(ValueError, match="heading fields"):
        _write_header(handle=io.StringIO(), fields=("valid", "bad\nfield"))
    with pytest.raises(InputValidationError, match="does not exist"):
        publish_membership_tables(
            tables_dir=tmp_path / "absent",
            legacy_sources=(),
            hog_sources=(),
            run_id="run",
        )
    tables = tmp_path / "tables"
    tables.mkdir()
    with pytest.raises(InputValidationError, match="interpreter is unavailable"):
        publish_membership_tables(
            tables_dir=tables,
            legacy_sources=(),
            hog_sources=(),
            run_id="run",
            worker_python=tmp_path / "missing-python",
        )

    monkeypatch.setattr(membership_publication.sys, "platform", "darwin")
    monkeypatch.setattr(
        membership_publication.resource,
        "getrusage",
        lambda _: SimpleNamespace(ru_maxrss=2 * 1024 * 1024),
    )
    assert membership_publication.peak_rss_mib() == 2.0


def test_header_writer_matches_csv_heading() -> None:
    """Raw heading publication retains the existing DictWriter contract."""

    raw = io.StringIO()
    _write_header(handle=raw, fields=("first", "second"))
    expected = io.StringIO()
    csv.DictWriter(
        expected,
        fieldnames=("first", "second"),
        delimiter="\t",
        lineterminator="\n",
    ).writeheader()
    assert raw.getvalue() == expected.getvalue()

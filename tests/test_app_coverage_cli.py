"""End-to-end command-line tests for selection coverage publication."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from orthofinder_interrogation_app import coverage_cli
from orthofinder_interrogation_app.coverage_exports import EXPORT_FILENAMES
from orthofinder_interrogation_app.models import FocusClusterPage
from orthofinder_results.cli import main


def _base_arguments(*, application_resource: Path, taxonomy_mapping_file: Path) -> list[str]:
    """Return named coverage-tree arguments shared by command-line tests."""

    return [
        "--action",
        "coverage-tree",
        "--resource-dir",
        str(application_resource),
        "--taxonomy-map",
        str(taxonomy_mapping_file),
        "--coverage-group-type",
        "HOG",
        "--coverage-hierarchy-node",
        "N0",
        "--include-clade-tax-id",
        "10",
        "--exclude-exact-tax-id",
        "103",
        "--coverage-created-at-utc",
        "2026-09-09T12:00:00Z",
    ]


def test_coverage_cli_validate_dry_run_and_atomic_publication(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Named CLI modes validate, dry-run and publish without touching the resource."""

    base = _base_arguments(
        application_resource=application_resource,
        taxonomy_mapping_file=taxonomy_mapping_file,
    )
    assert main([*base, "--coverage-all-groups", "--validate-only"]) == 0
    validated = json.loads(capsys.readouterr().out)
    assert validated["status"] == "validated"
    assert main([*base, "--coverage-all-groups", "--dry-run"]) == 0
    dry = json.loads(capsys.readouterr().out)
    assert dry["summary"]["evaluated_groups"] == 3
    first = tmp_path / "coverage_first"
    assert main([*base, "--coverage-all-groups", "--output-dir", str(first)]) == 0
    first_manifest = json.loads(capsys.readouterr().out)
    assert first_manifest["summary"]["passing_groups"] == 1
    assert set(path.name for path in first.iterdir()) == set(EXPORT_FILENAMES)
    second = tmp_path / "coverage_second"
    assert main([*base, "--coverage-all-groups", "--output-dir", str(second)]) == 0
    capsys.readouterr()
    assert {path.name: path.read_bytes() for path in first.iterdir()} == {
        path.name: path.read_bytes() for path in second.iterdir()
    }


def test_coverage_cli_defaults_to_replaceable_e3_focus_scope(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """A custom focus file replaces E3 seeds while retaining focus-limited semantics."""

    focus = tmp_path / "focus.tsv"
    focus.write_text(
        "protein_identifier\tprotein_name\tcategory\nalpha_1\tAlpha\tE3 test seed\n",
        encoding="utf-8",
    )
    output = tmp_path / "focus_coverage"
    arguments = _base_arguments(
        application_resource=application_resource,
        taxonomy_mapping_file=taxonomy_mapping_file,
    )
    assert (
        main(
            [
                *arguments,
                "--focus-proteins",
                str(focus),
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    manifest = json.loads(capsys.readouterr().out)
    assert manifest["focus_authority"]["limited_to_focus_clusters"]
    assert manifest["summary"]["evaluated_groups"] == 1
    audit = (output / "group_taxon_evaluation.tsv").read_text(encoding="utf-8")
    assert "N0.HOG1" in audit and "N0.HOG2" not in audit


def test_coverage_cli_rejects_conflicts_bad_modes_and_resource_mutation(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    tmp_path: Path,
) -> None:
    """Parser and domain validation fail before contradictory or unsafe publication."""

    base = _base_arguments(
        application_resource=application_resource,
        taxonomy_mapping_file=taxonomy_mapping_file,
    )
    with pytest.raises(SystemExit) as incompatible:
        main([*base, "--validate-only", "--dry-run"])
    assert incompatible.value.code == 2
    assert (
        main(
            [
                *base,
                "--require-exact-tax-id",
                "103",
                "--coverage-all-groups",
                "--validate-only",
            ]
        )
        == 2
    )
    assert (
        main(
            [
                *base,
                "--coverage-all-groups",
                "--output-dir",
                str(application_resource / "forbidden"),
            ]
        )
        == 2
    )
    with pytest.raises(SystemExit) as legacy:
        main(
            [
                "--action",
                "coverage-tree",
                "--resource-dir",
                str(application_resource),
                "--taxonomy-map",
                str(taxonomy_mapping_file),
                "--coverage-group-type",
                "LEGACY_ORTHOGROUP",
                "--coverage-hierarchy-node",
                "N0",
                "--validate-only",
            ]
        )
    assert legacy.value.code == 2
    assert not (tmp_path / "unused").exists()


def test_coverage_cli_requires_output_and_handles_zero_focus_matches(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Publication needs a target while a complete zero-match focus audit remains valid."""

    base = _base_arguments(
        application_resource=application_resource,
        taxonomy_mapping_file=taxonomy_mapping_file,
    )
    assert main([*base, "--coverage-all-groups"]) == 2
    focus = tmp_path / "absent_focus.tsv"
    focus.write_text("protein_identifier\nnot_present_in_resource\n", encoding="utf-8")
    output = tmp_path / "zero_focus_coverage"
    assert (
        main(
            [
                *base,
                "--focus-proteins",
                str(focus),
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    manifest = json.loads(capsys.readouterr().out)
    assert manifest["summary"]["evaluated_groups"] == 0
    assert manifest["focus_authority"]["limited_to_focus_clusters"] is True


def test_coverage_cli_rejects_truncated_focus_scope(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A partial focus-to-cluster query cannot be published as a complete audit."""

    focus = tmp_path / "focus.tsv"
    focus.write_text("protein_identifier\nalpha_1\n", encoding="utf-8")
    monkeypatch.setattr(
        coverage_cli.OrthoFinderQueryService,
        "search_focus_clusters",
        lambda self, *, filters: FocusClusterPage(
            rows=({"group_id": "N0.HOG1"},),
            total_rows=2,
            matched_focus_identifiers=1,
            submitted_focus_identifiers=1,
        ),
    )
    assert (
        main(
            [
                *_base_arguments(
                    application_resource=application_resource,
                    taxonomy_mapping_file=taxonomy_mapping_file,
                ),
                "--focus-proteins",
                str(focus),
                "--output-dir",
                str(tmp_path / "must_not_publish"),
            ]
        )
        == 2
    )
    assert not (tmp_path / "must_not_publish").exists()


def test_coverage_cli_resource_identity_falls_back_to_database_metadata(
    tmp_path: Path,
) -> None:
    """A standalone DuckDB receives an explicit size, run and schema identity."""

    database = tmp_path / "duckdb" / "orthofinder_results.duckdb"
    database.parent.mkdir()
    database.write_bytes(b"database")
    service = SimpleNamespace(
        resource=SimpleNamespace(
            database_path=database,
            run_id="standalone",
            schema_version=3,
        )
    )
    identity = coverage_cli._resource_identity(
        resource_path=database,
        service=service,
    )
    assert identity == "duckdb:standalone:size=8:schema=3"

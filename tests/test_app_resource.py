"""Tests for immutable resource opening and validation."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import duckdb
import pytest

from orthofinder_interrogation_app.resource import connect_read_only, open_resource
from orthofinder_results.errors import InputValidationError


def test_completed_resource_and_direct_database_open_read_only(
    application_resource: Path, tmp_path: Path
) -> None:
    """Directory and DuckDB modes return the same scientific run identity."""

    identity = open_resource(path=application_resource)
    assert identity.run_id == "test_run"
    assert identity.schema_version == 2
    assert identity.resource_package_version == "0.1.5"
    assert identity.counts["group_count"] == 4
    assert identity.report_path is not None
    direct_database = tmp_path / "standalone.duckdb"
    shutil.copyfile(identity.database_path, direct_database)
    direct = open_resource(path=direct_database)
    assert direct.run_id == identity.run_id
    assert direct.resource_path == direct_database
    assert direct.report_path is None
    assert direct.resource_package_version == "unknown"
    connection = connect_read_only(database_path=direct_database)
    try:
        with pytest.raises(duckdb.Error):
            connection.execute("CREATE TABLE forbidden(value INTEGER)")
    finally:
        connection.close()


@pytest.mark.parametrize("value", ("missing", "file.txt"))
def test_unsupported_resource_paths_fail(tmp_path: Path, value: str) -> None:
    """Only a directory or DuckDB file can be opened."""

    path = tmp_path / value
    if path.suffix:
        path.write_text("not a database", encoding="utf-8")
    with pytest.raises(InputValidationError, match="directory or a .duckdb"):
        open_resource(path=path)


def test_missing_or_empty_database_fails(tmp_path: Path) -> None:
    """A nominal resource directory must contain a non-empty database."""

    resource = tmp_path / "resource"
    (resource / "duckdb").mkdir(parents=True)
    with pytest.raises(InputValidationError, match="missing or empty"):
        open_resource(path=resource)
    (resource / "duckdb" / "orthofinder_results.duckdb").touch()
    with pytest.raises(InputValidationError, match="missing or empty"):
        open_resource(path=resource)


@pytest.mark.parametrize(
    ("manifest_text", "message"),
    [
        ("{\n", "unreadable"),
        ("[]\n", "JSON object"),
        ('{"status":"running"}\n', "not marked complete"),
    ],
)
def test_invalid_manifests_fail(
    application_resource: Path, manifest_text: str, message: str
) -> None:
    """Malformed and incomplete manifests cannot be presented as resources."""

    (application_resource / "run_manifest.json").write_text(manifest_text, encoding="utf-8")
    with pytest.raises(InputValidationError, match=message):
        open_resource(path=application_resource)


def test_missing_manifest_fails(application_resource: Path) -> None:
    """Resource-directory mode requires the completed manifest."""

    (application_resource / "run_manifest.json").unlink()
    with pytest.raises(InputValidationError, match="lacks run_manifest"):
        open_resource(path=application_resource)


@pytest.mark.parametrize(
    ("manifest_change", "database_sql", "message"),
    [
        ({"run_id": "other"}, "", "identifiers disagree"),
        ({"schema_version": 3}, "", "schema versions disagree"),
        ({"schema_version": "invalid"}, "", "not an integer"),
        ({}, "UPDATE resource_metadata SET schema_version = 3", "Unsupported resource schema"),
        ({}, "DELETE FROM group_statistics", "exactly one run_id"),
        (
            {},
            "INSERT INTO group_statistics SELECT 'other', * EXCLUDE (run_id) "
            "FROM group_statistics LIMIT 1",
            "exactly one run_id",
        ),
        ({}, "DROP TABLE species", "lacks required relations"),
    ],
)
def test_manifest_and_database_inconsistencies_fail(
    application_resource: Path,
    manifest_change: dict[str, object],
    database_sql: str,
    message: str,
) -> None:
    """Schema, relation and run mismatches are rejected before UI queries."""

    manifest_path = application_resource / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update(manifest_change)
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    if database_sql:
        connection = duckdb.connect(
            str(application_resource / "duckdb" / "orthofinder_results.duckdb")
        )
        try:
            connection.execute(database_sql)
            connection.execute("CHECKPOINT")
        finally:
            connection.close()
    with pytest.raises(InputValidationError, match=message):
        open_resource(path=application_resource)


def test_database_open_and_inspection_failures_are_controlled(tmp_path: Path) -> None:
    """Unreadable and structurally invalid databases produce domain errors."""

    invalid = tmp_path / "invalid.duckdb"
    invalid.write_text("not duckdb", encoding="utf-8")
    with pytest.raises(InputValidationError, match="Could not open"):
        open_resource(path=invalid)
    empty = tmp_path / "empty.duckdb"
    connection = duckdb.connect(str(empty))
    connection.close()
    with pytest.raises(InputValidationError, match="required relations"):
        open_resource(path=empty)

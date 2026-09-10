"""Defensive read-only opening of completed OrthoFinder resources."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import duckdb

from orthofinder_results.errors import InputValidationError

from .models import ResourceIdentity

_LOGGER = logging.getLogger("orthofinder_interrogation_app.resource")
SUPPORTED_SCHEMA_VERSIONS = frozenset({2, 3, 4})
REQUIRED_RELATIONS = frozenset(
    {
        "distance_statistics",
        "group_species_statistics",
        "group_statistics",
        "hog_memberships",
        "legacy_orthogroup_memberships",
        "resource_metadata",
        "species",
    }
)


def open_resource(*, path: Path) -> ResourceIdentity:
    """Validate and describe a completed resource without modifying it.

    Args:
        path: Completed resource directory or its physical DuckDB file.

    Returns:
        Validated resource identity.

    Raises:
        InputValidationError: If the resource is missing, incomplete or incompatible.
    """

    source = Path(path).expanduser().resolve()
    if source.is_dir():
        resource_dir = source
        database_path = resource_dir / "duckdb" / "orthofinder_results.duckdb"
        manifest_path = resource_dir / "run_manifest.json"
    elif source.is_file() and source.suffix.lower() == ".duckdb":
        resource_dir = _resource_root_for_database(database_path=source)
        database_path = source
        manifest_path = resource_dir / "run_manifest.json" if resource_dir else None
    else:
        raise InputValidationError(
            "Resource path must be a completed resource directory or a .duckdb file: "
            f"{source}"
        )
    if not database_path.is_file() or database_path.stat().st_size == 0:
        raise InputValidationError(f"Resource DuckDB is missing or empty: {database_path}")

    manifest = _read_manifest(path=manifest_path) if manifest_path is not None else {}
    schema_version, relations, run_ids = _inspect_database(path=database_path)
    missing = sorted(REQUIRED_RELATIONS - relations)
    if missing:
        raise InputValidationError(
            "Resource DuckDB lacks required relations: " + "; ".join(missing)
        )
    if schema_version not in SUPPORTED_SCHEMA_VERSIONS:
        supported = ", ".join(str(value) for value in sorted(SUPPORTED_SCHEMA_VERSIONS))
        raise InputValidationError(
            f"Unsupported resource schema {schema_version}; supported schemas: {supported}."
        )
    if schema_version >= 3 and "tree_payloads" not in relations:
        raise InputValidationError(
            "Schema-3-or-newer resource DuckDB lacks required relation: tree_payloads"
        )
    if len(run_ids) != 1:
        raise InputValidationError(
            "Resource must contain exactly one run_id; observed: "
            + ("; ".join(run_ids) if run_ids else "none")
        )
    run_id = run_ids[0]
    if manifest and str(manifest.get("run_id", "")) != run_id:
        raise InputValidationError(
            "Manifest and DuckDB run identifiers disagree: "
            f"{manifest.get('run_id')!r} versus {run_id!r}."
        )
    manifest_schema = manifest.get("schema_version")
    if manifest_schema is not None:
        try:
            parsed_manifest_schema = int(manifest_schema)
        except (TypeError, ValueError) as error:
            raise InputValidationError(
                f"Manifest schema_version is not an integer: {manifest_schema!r}."
            ) from error
        if parsed_manifest_schema != schema_version:
            raise InputValidationError(
                "Manifest and DuckDB schema versions disagree: "
                f"{manifest_schema!r} versus {schema_version!r}."
            )
    report_path = (
        resource_dir / "report" / "orthofinder_results_summary.html"
        if resource_dir is not None
        else None
    )
    identity = ResourceIdentity(
        resource_path=resource_dir or database_path,
        database_path=database_path,
        report_path=report_path if report_path is not None and report_path.is_file() else None,
        run_id=run_id,
        schema_version=schema_version,
        resource_package_version=str(manifest.get("package_version", "unknown")),
        orthofinder_version=str(manifest.get("orthofinder_version", "unknown")),
        adapter_name=str(manifest.get("adapter_name", "unknown")),
        primary_group_authority=str(manifest.get("primary_group_authority", "unknown")),
        counts=_integer_counts(record=manifest.get("counts", {})),
        relations=frozenset(relations),
    )
    _LOGGER.info(
        "Validated read-only resource: run=%s, schema=%s, database=%s",
        identity.run_id,
        identity.schema_version,
        identity.database_path,
    )
    return identity


def connect_read_only(*, database_path: Path) -> duckdb.DuckDBPyConnection:
    """Open one DuckDB connection with physical writes disabled.

    Args:
        database_path: Existing DuckDB resource.

    Returns:
        Caller-owned read-only DuckDB connection.

    Raises:
        InputValidationError: If DuckDB cannot be opened read-only.
    """

    source = Path(database_path).expanduser().resolve()
    try:
        return duckdb.connect(str(source), read_only=True)
    except (duckdb.Error, OSError) as error:
        raise InputValidationError(
            f"Could not open resource DuckDB read-only: {source}: {error}"
        ) from error


def _resource_root_for_database(*, database_path: Path) -> Path | None:
    """Find the conventional resource root for a direct DuckDB path."""

    parent = database_path.parent
    candidate = parent.parent if parent.name == "duckdb" else None
    if candidate is not None and (candidate / "run_manifest.json").is_file():
        return candidate
    return None


def _read_manifest(*, path: Path) -> dict[str, Any]:
    """Read and validate a completed resource manifest."""

    if not path.is_file():
        raise InputValidationError(f"Completed resource lacks run_manifest.json: {path}")
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InputValidationError(f"Resource manifest is unreadable: {path}: {error}") from error
    if not isinstance(record, dict):
        raise InputValidationError(f"Resource manifest must contain a JSON object: {path}")
    if record.get("status") != "complete":
        raise InputValidationError(f"Resource manifest is not marked complete: {path}")
    return record


def _inspect_database(*, path: Path) -> tuple[int, set[str], tuple[str, ...]]:
    """Read schema, relation and run identities from a DuckDB file."""

    connection = connect_read_only(database_path=path)
    try:
        relations = {str(row[0]) for row in connection.execute("SHOW TABLES").fetchall()}
        if "resource_metadata" not in relations or "group_statistics" not in relations:
            return -1, relations, ()
        schema_rows = connection.execute(
            "SELECT DISTINCT schema_version FROM resource_metadata"
        ).fetchall()
        run_rows = connection.execute(
            "SELECT DISTINCT run_id FROM group_statistics ORDER BY run_id"
        ).fetchall()
    except duckdb.Error as error:
        raise InputValidationError(f"Could not inspect resource DuckDB {path}: {error}") from error
    finally:
        connection.close()
    if len(schema_rows) != 1:
        raise InputValidationError(
            "Resource metadata must contain exactly one schema version."
        )
    try:
        schema_version = int(schema_rows[0][0])
    except (TypeError, ValueError) as error:
        raise InputValidationError(
            f"Resource schema version is not an integer: {schema_rows[0][0]!r}."
        ) from error
    return schema_version, relations, tuple(str(row[0]) for row in run_rows)


def _integer_counts(*, record: object) -> dict[str, int]:
    """Retain valid manifest integer counts without fabricating missing values."""

    if not isinstance(record, dict):
        return {}
    return {
        str(key): value
        for key, value in record.items()
        if isinstance(value, int) and not isinstance(value, bool)
    }

"""Tests for schema-3 portable trees and lazy distance sidecars."""

from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path

import duckdb
import pytest

import orthofinder_interrogation_app.distance_data as distance_data
from orthofinder_interrogation_app.distance_data import DistanceAnalysisProvider
from orthofinder_interrogation_app.models import GroupKey
from orthofinder_interrogation_app.queries import OrthoFinderQueryService
from orthofinder_interrogation_app.report_data import load_visualisation_catalog
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_results.errors import InputValidationError


def _service(resource_path: Path) -> OrthoFinderQueryService:
    """Open one test resource and return its query service."""

    return OrthoFinderQueryService(resource=open_resource(path=resource_path))


def _key(*, run_id: str, node: str = "N0") -> GroupKey:
    """Return the synthetic HOG key at one hierarchy node."""

    return GroupKey(
        run_id=run_id,
        group_type="HOG",
        hierarchy_node=node,
        group_id=f"{node}.HOG0000001",
    )


def test_schema3_lazy_analysis_is_exact_cached_and_forceable(
    schema3_resource: Path,
    persistent_test_root: Path,
) -> None:
    """A non-precomputed HOG resolves its OG tree and writes a reusable sidecar."""

    service = _service(schema3_resource)
    provider = DistanceAnalysisProvider(
        service=service,
        cache_dir=persistent_test_root / "analysis_cache",
    )
    key = _key(run_id="schema3-run", node="N1")
    first = provider.analyse(key=key, max_members=3, nearest_neighbours=2)
    assert first.source == "PORTABLE_TREE_LAZY"
    assert first.tree_authority == "RESOLVED_GENE_TREE"
    assert first.cache_status == "CACHE_WRITE"
    assert len(first.members) == 3
    assert len(first.distances) == 3
    assert {row["species_a"] for row in first.distances} <= {"Species_A", "Species_B"}
    assert first.summary["distance_pair_count"] == 3
    assert first.visual_entry["distanceMatrix"]["status"] == (
        "EXACT_COMPLETE_DISPLAYED_MATRIX"
    )
    second = provider.analyse(key=key, max_members=3, nearest_neighbours=2)
    assert second.cache_status == "CACHE_HIT"
    assert second.distances == first.distances
    forced = provider.analyse(
        key=key,
        max_members=3,
        nearest_neighbours=2,
        force_recompute=True,
    )
    assert forced.cache_status == "CACHE_WRITE"


def test_schema3_lazy_analysis_retains_a_required_focus_protein(
    schema3_resource: Path,
    persistent_test_root: Path,
) -> None:
    """A protein lookup cannot silently sample its requested protein away."""

    service = _service(schema3_resource)
    provider = DistanceAnalysisProvider(
        service=service,
        cache_dir=persistent_test_root / "focused_analysis_cache",
    )
    key = _key(run_id="schema3-run", node="N1")
    focused = provider.analyse(
        key=key,
        max_members=2,
        nearest_neighbours=2,
        required_members=("protB",),
    )
    assert focused.summary["sampled_member_count"] == 2
    assert "protB" in {str(row["member_id"]) for row in focused.members}
    assert any(
        "protB" in (row["member_a"], row["member_b"])
        for row in focused.distances
    )
    repeated = provider.analyse(
        key=key,
        max_members=2,
        nearest_neighbours=2,
        required_members=("protB",),
    )
    assert repeated.cache_status == "CACHE_HIT"


def test_provider_prefers_persisted_pairs_and_schema2_report(
    schema3_resource: Path,
    application_resource: Path,
    persistent_test_root: Path,
) -> None:
    """Existing exact pair rows and the old pilot matrix remain first-class inputs."""

    schema3_service = _service(schema3_resource)
    schema3_report_path = schema3_service.resource.report_path
    assert schema3_report_path is not None
    schema3_catalog = load_visualisation_catalog(
        report_path=schema3_report_path,
        expected_run_id="schema3-run",
    )
    persisted = DistanceAnalysisProvider(
        service=schema3_service,
        cache_dir=persistent_test_root / "schema3_cache",
        report_catalog=schema3_catalog,
    ).analyse(key=_key(run_id="schema3-run"), max_members=3)
    assert persisted.source == "PERSISTED_DUCKDB_DISTANCES"
    assert persisted.cache_status == "NOT_APPLICABLE"
    assert persisted.tree_authority == ""

    schema2_service = _service(application_resource)
    report_path = schema2_service.resource.report_path
    assert report_path is not None
    catalog = load_visualisation_catalog(
        report_path=report_path,
        expected_run_id="test_run",
    )
    report = DistanceAnalysisProvider(
        service=schema2_service,
        cache_dir=persistent_test_root / "schema2_cache",
        report_catalog=catalog,
    ).analyse(
        key=GroupKey("test_run", "HOG", "N0", "N0.HOG1"),
        max_members=3,
    )
    assert report.source == "PERSISTED_REPORT_MATRIX"
    assert report.summary["mean_distance"] == pytest.approx(0.2)
    assert [row["distance"] for row in report.distances] == [0.1, 0.2, 0.3]


def test_provider_reports_controls_identity_and_missing_capabilities(
    application_resource: Path,
    persistent_test_root: Path,
) -> None:
    """Unsafe controls, cross-run keys and schema-2 gaps fail explicitly."""

    service = _service(application_resource)
    provider = DistanceAnalysisProvider(
        service=service,
        cache_dir=persistent_test_root / "cache",
    )
    key = GroupKey("test_run", "HOG", "N0", "N0.HOG3")
    with pytest.raises(InputValidationError, match="between 2"):
        provider.analyse(key=key, max_members=1)
    with pytest.raises(InputValidationError, match="Nearest neighbours"):
        provider.analyse(key=key, nearest_neighbours=0)
    with pytest.raises(InputValidationError, match="does not match"):
        provider.analyse(key=GroupKey("other", "HOG", "N0", "N0.HOG3"))
    with pytest.raises(InputValidationError, match="schema-2 resource"):
        provider.analyse(key=key)
    with pytest.raises(InputValidationError, match="exactly one group"):
        provider.analyse(key=GroupKey("test_run", "HOG", "N0", "missing"))
    with pytest.raises(InputValidationError, match="outside"):
        DistanceAnalysisProvider(
            service=service,
            cache_dir=application_resource / "unsafe-cache",
        )


def test_lazy_analysis_rejects_truncated_membership_authority(
    schema3_resource: Path,
    persistent_test_root: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Lazy analysis must never treat a browser row limit as the complete group."""

    service = _service(schema3_resource)
    key = _key(run_id="schema3-run", node="N1")
    group = service.get_group(key=key)
    provider = DistanceAnalysisProvider(
        service=service,
        cache_dir=persistent_test_root / "guard_cache",
    )
    with monkeypatch.context() as patch:
        patch.setattr(
            service,
            "get_group",
            lambda *, key: {
                **group,
                "member_count": distance_data.MAX_GROUP_MEMBER_ROWS + 1,
            },
        )
        with pytest.raises(InputValidationError, match="safe browser-analysis limit"):
            provider.analyse(key=key)
    with monkeypatch.context() as patch:
        patch.setattr(service, "get_group_members", lambda *, key: ())
        with pytest.raises(InputValidationError, match="Membership authority is incomplete"):
            provider.analyse(key=key)


def test_query_service_exposes_portable_tree_pairs_and_aliases(
    schema3_resource: Path,
) -> None:
    """Schema-3 query methods preserve bounded relation and tree semantics."""

    service = _service(schema3_resource)
    n0 = _key(run_id="schema3-run")
    n1 = _key(run_id="schema3-run", node="N1")
    assert service.overview_counts()["portable_tree_count"] == 1
    assert service.has_relation(relation="tree_payloads")
    assert service.overview_authorities()
    aliases = service.get_group_sequence_aliases(key=n1)
    assert {(row["member_id"], row["internal_id"]) for row in aliases} == {
        ("protA", "0_0"),
        ("protA2", "0_1"),
        ("protB", "1_0"),
    }
    tree = service.get_portable_tree(key=n1, legacy_orthogroup_id="OG0000001")
    assert tree is not None
    assert tree["tree_id"] == "OG0000001"
    assert tree["tree_type"] == "RESOLVED_GENE_TREE"
    assert len(service.get_group_distances(key=n0)) == 3
    with pytest.raises(InputValidationError, match="contain 3 rows"):
        service.get_group_distances(key=n0, maximum=2)
    with pytest.raises(InputValidationError, match="maximum"):
        service.get_group_distances(key=n0, maximum=0)
    assert service.get_portable_tree(
        key=GroupKey("schema3-run", "HOG", "N0", "unrelated")
    ) is None
    with pytest.raises(InputValidationError, match="Unsupported group type"):
        service.get_group_sequence_aliases(
            key=GroupKey("schema3-run", "UNKNOWN", "", "group")
        )


def test_duplicate_portable_tree_authority_is_rejected(
    schema3_resource: Path,
) -> None:
    """Two equally ranked payload records cannot be silently selected."""

    database = schema3_resource / "duckdb/orthofinder_results.duckdb"
    connection = duckdb.connect(str(database))
    try:
        connection.execute(
            "INSERT INTO tree_payloads SELECT * FROM tree_payloads "
            "WHERE tree_id = 'OG0000001'"
        )
        connection.execute("CHECKPOINT")
    finally:
        connection.close()
    service = _service(schema3_resource)
    with pytest.raises(InputValidationError, match="Multiple equally preferred"):
        service.get_portable_tree(
            key=_key(run_id="schema3-run", node="N1"),
            legacy_orthogroup_id="OG0000001",
        )


def test_member_alias_and_distance_decoration_defences() -> None:
    """Canonical aliases and exact endpoint species reject ambiguous input."""

    members = (
        {"member_id": "a", "species_label": "Species_A"},
        {"member_id": "b", "species_label": "Species_B"},
    )
    aliases = ({"member_id": "a", "species_label": "Species_A", "internal_id": "0_0"},)
    member_ids, mapping = distance_data._member_aliases(
        members=members,
        sequence_aliases=aliases,
    )
    assert member_ids == ("a", "b")
    assert mapping["a"]["0_0"] == "ORTHOFINDER_INTERNAL_ID"
    rows = ({"member_a": "a", "member_b": "b", "distance": 1.0},)
    decorated = distance_data._decorate_distance_rows(rows=rows, members=members)
    assert decorated[0]["species_b"] == "Species_B"
    with pytest.raises(InputValidationError, match="multiple species"):
        distance_data._member_aliases(
            members=members + ({"member_id": "a", "species_label": "Species_C"},),
            sequence_aliases=(),
        )
    with pytest.raises(InputValidationError, match="more than one species"):
        distance_data._decorate_distance_rows(
            rows=rows,
            members=members + ({"member_id": "a", "species_label": "Species_C"},),
        )
    with pytest.raises(InputValidationError, match="absent"):
        distance_data._decorate_distance_rows(
            rows=({"member_a": "a", "member_b": "missing", "distance": 1.0},),
            members=members,
        )


def test_report_matrix_and_cache_records_are_validated(
    application_resource: Path,
    persistent_test_root: Path,
) -> None:
    """Malformed report matrices and corrupt sidecars are controlled failures."""

    service = _service(application_resource)
    report_path = service.resource.report_path
    assert report_path is not None
    catalog = load_visualisation_catalog(
        report_path=report_path,
        expected_run_id="test_run",
    )
    key = GroupKey("test_run", "HOG", "N0", "N0.HOG1")
    entry = catalog.get(label="HOG|N0|N0.HOG1")
    group = service.get_group(key=key)
    malformed = copy.deepcopy(entry)
    malformed["distanceMatrix"]["upperTriangle"] = [0.1]
    with pytest.raises(InputValidationError, match="requires 3"):
        distance_data._analysis_from_report(key=key, group=group, entry=malformed)
    malformed = copy.deepcopy(entry)
    malformed["members"] = "wrong"
    with pytest.raises(InputValidationError, match="lacks members"):
        distance_data._analysis_from_report(key=key, group=group, entry=malformed)

    analysis = distance_data._analysis_from_report(key=key, group=group, entry=entry)
    record = distance_data._analysis_record(analysis=analysis)
    restored = distance_data._cached_analysis_from_record(
        record=record,
        expected_key=key,
    )
    assert restored.distances == analysis.distances
    bad_records = (
        [],
        {**record, "cache_format": 99},
        {**record, "key": "wrong"},
        {**record, "group": "wrong"},
        {**record, "members": "wrong"},
        {**record, "members": ["wrong"]},
    )
    for bad in bad_records:
        with pytest.raises(InputValidationError):
            distance_data._cached_analysis_from_record(record=bad, expected_key=key)
    wrong_key = copy.deepcopy(record)
    wrong_key["key"]["group_id"] = "other"
    with pytest.raises(InputValidationError, match="identity"):
        distance_data._cached_analysis_from_record(record=wrong_key, expected_key=key)

    cache_path = persistent_test_root / "records" / "analysis.json.gz"
    distance_data._write_cached_analysis(path=cache_path, analysis=analysis)
    assert cache_path.stat().st_mode & 0o777 == 0o600
    assert distance_data._read_cached_analysis(path=cache_path, expected_key=key) is not None
    cache_path.write_bytes(b"not gzip")
    assert distance_data._read_cached_analysis(path=cache_path, expected_key=key) is None
    with gzip.open(cache_path, "wt", encoding="utf-8") as handle:
        json.dump(wrong_key, handle)
    assert distance_data._read_cached_analysis(path=cache_path, expected_key=key) is None
    assert distance_data._read_cached_analysis(
        path=persistent_test_root / "missing.json.gz",
        expected_key=key,
    ) is None


def test_cache_defaults_and_size_limits(
    monkeypatch: pytest.MonkeyPatch,
    application_resource: Path,
    persistent_test_root: Path,
) -> None:
    """Default paths avoid system temp and cache byte limits fail before publication."""

    monkeypatch.setattr(distance_data.sys, "platform", "darwin")
    monkeypatch.setenv("HOME", "/Users/Tester")
    assert str(distance_data.default_cache_directory()).endswith(
        "Library/Caches/orthofinder-results"
    )
    monkeypatch.setattr(distance_data.sys, "platform", "linux")
    monkeypatch.setenv("XDG_CACHE_HOME", str(persistent_test_root / "xdg"))
    assert distance_data.default_cache_directory() == (
        persistent_test_root / "xdg/orthofinder-results"
    )

    service = _service(application_resource)
    report_path = service.resource.report_path
    assert report_path is not None
    entry = load_visualisation_catalog(
        report_path=report_path,
        expected_run_id="test_run",
    ).get(label="HOG|N0|N0.HOG1")
    key = GroupKey("test_run", "HOG", "N0", "N0.HOG1")
    analysis = distance_data._analysis_from_report(
        key=key,
        group=service.get_group(key=key),
        entry=entry,
    )
    monkeypatch.setattr(distance_data, "MAX_CACHE_BYTES", 10)
    with pytest.raises(InputValidationError, match="cache limit"):
        distance_data._write_cached_analysis(
            path=persistent_test_root / "too_large.json.gz",
            analysis=analysis,
        )


def test_tree_identifier_helpers_are_deterministic() -> None:
    """HOG-to-OG conversion and tree ranks remain explicit and stable."""

    from orthofinder_interrogation_app import queries

    assert queries._tree_id_from_hog(group_id="N42.HOG0000123") == "OG0000123"
    assert queries._tree_id_from_hog(group_id="OG0000123") == ""
    candidates = ("OG1", "N0.HOG1")
    assert queries._tree_priority(
        tree_type="RESOLVED_GENE_TREE",
        tree_id="OG1",
        candidate_ids=candidates,
    ) == (0, 0)
    assert queries._tree_priority(
        tree_type="OTHER",
        tree_id="missing",
        candidate_ids=candidates,
    ) == (2, 2)

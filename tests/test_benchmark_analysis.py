"""Tests for cluster-level matched-background dispersion statistics."""

from __future__ import annotations

import hashlib
from pathlib import Path

import duckdb
import pytest
from conftest import make_results

from orthofinder_interrogation_app.focus import FocusProtein, FocusProteinAuthority
from orthofinder_results.benchmark_analysis import (
    BenchmarkPlan,
    _classification,
    _match_controls,
    _optional_float,
    _quantile,
    benjamini_hochberg,
    build_profile_rows,
    mann_whitney_test,
    publish_benchmark_results,
)
from orthofinder_results.benchmark_authority import (
    BENCHMARK_FIELDS,
    BenchmarkAuthority,
    BenchmarkMarker,
)
from orthofinder_results.cli import main
from orthofinder_results.focus_analysis import FocusSelection
from orthofinder_results.io_utils import read_tsv


def _group(*, identifier: str, member_count: int = 6) -> dict[str, str]:
    """Return one structurally complete synthetic N0 HOG record."""

    return {
        "run_id": "benchmark-run",
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "group_id": identifier,
        "legacy_orthogroup_id": identifier.replace("N0.HOG", "OG"),
        "gene_tree_parent_clade": "N0",
        "member_count": str(member_count),
        "species_count": "3",
        "single_copy_species_count": "1",
        "max_copies_per_species": str(member_count - 2),
        "mean_copies_per_species": str(member_count / 3),
        "is_singleton": "False",
        "species_labels": "Species_A;Species_B;Species_C",
        "source_file": "N0.tsv",
    }


def test_statistical_boundary_helpers_keep_empty_and_tied_results_explicit() -> None:
    """Empty corrections, all-tie tests and central classifications are defined."""

    assert benjamini_hochberg(p_values=()) == ()
    _, p_value, effect = mann_whitney_test(
        target=(1.0, 1.0, 1.0),
        reference=(1.0, 1.0, 1.0),
    )
    assert p_value == pytest.approx(1.0)
    assert effect == pytest.approx(0.0)
    assert _classification(observed=2.0, background=(1.0, 2.0, 3.0)) == (
        2.0,
        0.5,
        "TYPICAL",
    )
    assert _optional_float(value=None) is None
    with pytest.raises(ValueError, match="must not be empty"):
        _quantile(values=(), fraction=0.5)


def _focus_selection(*, prefix: str, means: tuple[float, ...]) -> FocusSelection:
    """Return one three-cluster synthetic focus selection."""

    records = tuple(
        FocusProtein(
            identifier=f"{prefix}{index}",
            protein_name=f"{prefix} protein {index}",
            category="RING" if prefix == "e3" else prefix.upper(),
            evidence_type="reviewed",
            organism="Arabidopsis thaliana",
            source="test authority",
            enabled=True,
            note="",
        )
        for index, _ in enumerate(means, start=1)
    )
    authority = FocusProteinAuthority(
        records=records,
        source_name=f"{prefix}.tsv",
        sha256=hashlib.sha256(prefix.encode()).hexdigest(),
    )
    groups = tuple(
        _group(identifier=f"N0.HOG{prefix.upper()}{index}")
        for index, _ in enumerate(means, start=1)
    )
    matches = tuple(
        {
            "run_id": "benchmark-run",
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": group["group_id"],
            "legacy_orthogroup_id": group["legacy_orthogroup_id"],
            "seed_id": record.identifier,
            "seed_categories": record.category,
            "matched_member_id": record.identifier,
            "matched_species_label": "Arabidopsis_thaliana",
        }
        for group, record in zip(groups, records, strict=True)
    )
    return FocusSelection(
        authority=authority,
        group_type="HOG",
        hierarchy_node="N0",
        match_rows=matches,
        audit_rows=tuple({"match_status": "MATCHED"} for _ in records),
        group_statistics=groups,
        required_members_by_group={
            f"HOG|N0|{group['group_id']}": (record.identifier,)
            for group, record in zip(groups, records, strict=True)
        },
    )


def _marker_authority_and_selection() -> tuple[BenchmarkAuthority, FocusSelection]:
    """Return three housekeeping and three R/NLR target selections."""

    markers = []
    groups = []
    matches = []
    for marker_class, prefix, subclass in (
        ("HOUSEKEEPING", "hk", "STABLE_EXPRESSION_REFERENCE"),
        ("R_NLR", "nlr", "TNL"),
    ):
        for index in range(1, 4):
            identifier = f"{prefix}{index}"
            marker = BenchmarkMarker(
                marker_id=f"AT{index}G{index:05d}",
                protein_identifier=identifier,
                protein_entry=f"{identifier.upper()}_ARATH",
                marker_name=f"{marker_class} marker {index}",
                benchmark_class=marker_class,
                benchmark_subclass=subclass,
                domain_architecture="TIR,NB,LRR" if marker_class == "R_NLR" else "",
                evidence_type="reviewed",
                organism="Arabidopsis thaliana",
                taxon_id="3702",
                source_title="Test source",
                source_doi="10.0000/test",
                source_table="Table 1",
                source_version="1",
                enabled=True,
                note="Test marker.",
            )
            group = _group(identifier=f"N0.HOG{prefix.upper()}{index}")
            markers.append(marker)
            groups.append(group)
            matches.append(
                {
                    "run_id": "benchmark-run",
                    "group_type": "HOG",
                    "hierarchy_node": "N0",
                    "group_id": group["group_id"],
                    "legacy_orthogroup_id": group["legacy_orthogroup_id"],
                    "seed_id": identifier,
                    "seed_categories": f"{marker_class}::{subclass}",
                    "matched_member_id": identifier,
                    "matched_species_label": "Arabidopsis_thaliana",
                    "match_authority": "UNIPROT_ACCESSION",
                }
            )
    authority = BenchmarkAuthority(
        records=tuple(markers),
        source_name="benchmarks.tsv",
        sha256=hashlib.sha256(b"benchmarks").hexdigest(),
    )
    selection = FocusSelection(
        authority=authority.to_focus_authority(),
        group_type="HOG",
        hierarchy_node="N0",
        match_rows=tuple(matches),
        audit_rows=tuple({"match_status": "MATCHED"} for _ in markers),
        group_statistics=tuple(groups),
        required_members_by_group={
            f"HOG|N0|{group['group_id']}": (marker.protein_identifier,)
            for group, marker in zip(groups, markers, strict=True)
        },
    )
    return authority, selection


def _plan() -> tuple[BenchmarkPlan, dict[str, float]]:
    """Return a complete nine-target, 27-control benchmark plan."""

    e3 = _focus_selection(prefix="e3", means=(4.0, 5.0, 6.0))
    authority, marker_selection = _marker_authority_and_selection()
    targets = (*e3.group_statistics, *marker_selection.group_statistics)
    candidates = tuple(
        _group(identifier=f"N0.HOGCONTROL{index}", member_count=5 + index % 3)
        for index in range(1, 28)
    )
    controls, links = _match_controls(
        anchors=targets,
        candidates=candidates,
        controls_per_group=3,
    )
    plan = BenchmarkPlan(
        run_id="benchmark-run",
        e3_selection=e3,
        marker_selection=marker_selection,
        marker_authority=authority,
        target_statistics=tuple(targets),
        control_statistics=controls,
        control_links=links,
    )
    means = {
        **{
            f"HOG|N0|N0.HOGE3{index}": value
            for index, value in enumerate((4.0, 5.0, 6.0), start=1)
        },
        **{
            f"HOG|N0|N0.HOGHK{index}": value
            for index, value in enumerate((1.0, 1.2, 1.4), start=1)
        },
        **{
            f"HOG|N0|N0.HOGNLR{index}": value
            for index, value in enumerate((7.0, 8.0, 9.0), start=1)
        },
        **{
            f"HOG|N0|N0.HOGCONTROL{index}": 2.0 + (index % 3) * 0.1
            for index in range(1, 28)
        },
    }
    return plan, means


def _summary(*, key: str, mean: float) -> dict[str, object]:
    """Return one complete successful synthetic distance summary."""

    group_type, hierarchy, group_id = key.split("|")
    return {
        "run_id": "benchmark-run",
        "group_type": group_type,
        "hierarchy_node": hierarchy,
        "group_id": group_id,
        "distance_method": "patristic_branch_length",
        "computation_status": "EXACT",
        "member_identifier_resolution": "EXACT_MEMBER_ID",
        "total_member_count": 3,
        "sampled_member_count": 3,
        "distance_pair_count": 3,
        "unresolved_pair_count": 0,
        "minimum_distance": mean - 0.2,
        "q05_distance": mean - 0.18,
        "q25_distance": mean - 0.1,
        "median_distance": mean,
        "mean_distance": mean,
        "q75_distance": mean + 0.1,
        "q95_distance": mean + 0.18,
        "maximum_distance": mean + 0.2,
        "population_stddev_distance": mean / 10,
        "failure_reason": "",
    }


def test_fdr_and_mann_whitney_are_cluster_level_and_tie_corrected() -> None:
    """Core non-parametric helpers return known monotonic finite results."""

    assert benjamini_hochberg(p_values=(0.01, 0.04, 0.03)) == pytest.approx(
        (0.03, 0.04, 0.04)
    )
    u_value, p_value, delta = mann_whitney_test(
        target=(4.0, 5.0, 6.0),
        reference=(1.0, 2.0, 3.0),
    )
    assert u_value == pytest.approx(9.0)
    assert 0.0 <= p_value <= 1.0
    assert delta == pytest.approx(1.0)
    tied_u, tied_p, tied_delta = mann_whitney_test(
        target=(1.0, 1.0, 2.0),
        reference=(1.0, 2.0, 2.0),
    )
    assert tied_u == pytest.approx(3.0)
    assert 0.0 <= tied_p <= 1.0
    assert -1.0 <= tied_delta <= 1.0
    with pytest.raises(ValueError, match="finite"):
        mann_whitney_test(target=(float("nan"),), reference=(1.0,))
    with pytest.raises(ValueError, match="between zero and one"):
        benjamini_hochberg(p_values=(1.2,))


def test_benchmark_results_cover_profiles_contrasts_and_individual_clusters(
    tmp_path: Path,
) -> None:
    """Every selected cluster reaches calibrated, FDR-controlled result tables."""

    plan, means = _plan()
    profiles = build_profile_rows(plan=plan)
    assert {row["profile_id"] for row in profiles}.issuperset(
        {"E3_ALL", "E3_CATEGORY::RING", "HOUSEKEEPING_ALL", "R_NLR_ALL"}
    )
    summaries = tuple(_summary(key=key, mean=value) for key, value in means.items())
    counts = publish_benchmark_results(
        tables_dir=tmp_path,
        plan=plan,
        distance_summaries=summaries,
        bootstrap_resamples=100,
    )
    assert counts["benchmark_cluster_count"] == 36
    clusters = tuple(
        read_tsv(path=tmp_path / "benchmark_cluster_results.tsv.gz")
    )
    assert len(clusters) == 36
    assert sum(row["pairwise_rows_persisted"] == "true" for row in clusters) == 9

    contrasts = tuple(read_tsv(path=tmp_path / "benchmark_contrasts.tsv.gz"))
    tested = [row for row in contrasts if row["status"] == "TESTED"]
    assert tested
    assert all(0.0 <= float(row["fdr_q_value"]) <= 1.0 for row in tested)
    e3_housekeeping = next(
        row
        for row in tested
        if row["contrast_id"] == "E3_ALL__VS__HOUSEKEEPING_ALL"
        and row["metric"] == "mean_distance"
    )
    assert float(e3_housekeeping["median_difference"]) > 0

    individual = tuple(
        read_tsv(path=tmp_path / "benchmark_individual_comparisons.tsv.gz")
    )
    own_control = [
        row
        for row in individual
        if row["background_profile_id"] == "MATCHED_NON_FOCUS"
    ]
    assert own_control
    assert all(row["background_group_count"] == "3" for row in own_control)
    assert any(row["leave_one_out"] == "true" for row in individual)
    classifications = tuple(
        read_tsv(path=tmp_path / "benchmark_cluster_classifications.tsv.gz")
    )
    assert len(classifications) == 9
    assert all(row["status"] == "CLASSIFIED" for row in classifications)


def _extend_results_for_benchmark(*, results: Path) -> None:
    """Add marker targets and three eligible controls to a tiny result set."""

    sequence_path = results / "WorkingDirectory/SequenceIDs.txt"
    sequence_path.write_text(
        sequence_path.read_text(encoding="utf-8")
        + "0_2: hk1\n1_1: hk1b\n0_3: nlr1\n1_2: nlr1b\n"
        + "0_4: control1a\n1_3: control1b\n0_5: control2a\n"
        + "1_4: control2b\n0_6: control3a\n1_5: control3b\n",
        encoding="utf-8",
    )
    hog_path = results / "Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
    hog_path.write_text(
        hog_path.read_text(encoding="utf-8")
        + "N0.HOG0000002\tOG0000002\tN0\thk1\thk1b\n"
        + "N0.HOG0000003\tOG0000003\tN0\tnlr1\tnlr1b\n"
        + "N0.HOG0000004\tOG0000004\tN0\tcontrol1a\tcontrol1b\n"
        + "N0.HOG0000005\tOG0000005\tN0\tcontrol2a\tcontrol2b\n"
        + "N0.HOG0000006\tOG0000006\tN0\tcontrol3a\tcontrol3b\n",
        encoding="utf-8",
    )
    orthogroup_path = results / "Orthogroups/Orthogroups.tsv"
    orthogroup_path.write_text(
        orthogroup_path.read_text(encoding="utf-8")
        + "OG0000002\thk1\thk1b\n"
        + "OG0000003\tnlr1\tnlr1b\n"
        + "OG0000004\tcontrol1a\tcontrol1b\n"
        + "OG0000005\tcontrol2a\tcontrol2b\n"
        + "OG0000006\tcontrol3a\tcontrol3b\n",
        encoding="utf-8",
    )
    tree_members = {
        "OG0000002": ("hk1", "hk1b", 0.2),
        "OG0000003": ("nlr1", "nlr1b", 0.8),
        "OG0000004": ("control1a", "control1b", 0.3),
        "OG0000005": ("control2a", "control2b", 0.4),
        "OG0000006": ("control3a", "control3b", 0.5),
    }
    for tree_id, (left, right, branch) in tree_members.items():
        tree = f"({left}:{branch},{right}:{branch})N0:0.0;\n"
        for directory in ("Gene_Trees", "Resolved_Gene_Trees"):
            (results / directory / f"{tree_id}_tree.txt").write_text(
                tree,
                encoding="utf-8",
            )


def _write_benchmark_authority(*, path: Path) -> None:
    """Write one housekeeping and one R/NLR exact test marker."""

    common = {
        "protein_entry": "TEST_ARATH",
        "domain_architecture": "",
        "evidence_type": "reviewed",
        "organism": "Arabidopsis thaliana",
        "taxon_id": "3702",
        "source_title": "Test marker source",
        "source_doi": "10.0000/test",
        "source_table": "Table 1",
        "source_version": "1",
        "enabled": "true",
        "note": "Evaluate rather than assume.",
    }
    rows = (
        {
            **common,
            "marker_id": "AT1G00001",
            "protein_identifier": "hk1",
            "marker_name": "Housekeeping candidate",
            "benchmark_class": "HOUSEKEEPING",
            "benchmark_subclass": "STABLE_EXPRESSION_REFERENCE",
        },
        {
            **common,
            "marker_id": "AT1G00002",
            "protein_identifier": "nlr1",
            "marker_name": "R/NLR candidate",
            "benchmark_class": "R_NLR",
            "benchmark_subclass": "TNL",
            "domain_architecture": "TIR,NB,LRR",
        },
    )
    text = "\t".join(BENCHMARK_FIELDS) + "\n"
    text += "".join(
        "\t".join(str(row[field]) for field in BENCHMARK_FIELDS) + "\n"
        for row in rows
    )
    path.write_text(text, encoding="utf-8")


def test_dispersion_benchmark_cli_publishes_compact_controls_and_target_pairs(
    persistent_test_root: Path,
    tmp_path: Path,
) -> None:
    """The end-to-end action persists target matrices and control summaries."""

    results = make_results(root=tmp_path)
    _extend_results_for_benchmark(results=results)
    focus = tmp_path / "focus.tsv"
    focus.write_text("protein_identifier\nprotA\n", encoding="utf-8")
    benchmarks = tmp_path / "benchmarks.tsv"
    _write_benchmark_authority(path=benchmarks)
    output = persistent_test_root / "dispersion_benchmark"
    status = main(
        [
            "--action",
            "dispersion-benchmark",
            "--results-dir",
            str(results),
            "--output-dir",
            str(output),
            "--run-id",
            "benchmark-cli-test",
            "--work-dir",
            str(persistent_test_root / "benchmark_work"),
            "--focus-proteins",
            str(focus),
            "--benchmark-proteins",
            str(benchmarks),
            "--benchmark-controls-per-group",
            "1",
            "--benchmark-bootstrap-resamples",
            "100",
            "--distance-max-members",
            "3",
        ]
    )
    assert status == 0
    cluster_rows = tuple(
        read_tsv(path=output / "tables/benchmark_cluster_results.tsv.gz")
    )
    assert len(cluster_rows) == 6
    assert sum(row["pairwise_rows_persisted"] == "true" for row in cluster_rows) == 3
    pair_rows = tuple(read_tsv(path=output / "tables/pairwise_distances.tsv.gz"))
    assert len(pair_rows) == 5
    manifest_checks = {
        row["check_name"]: row["status"]
        for row in read_tsv(path=output / "qc/validation_checks.tsv")
    }
    assert manifest_checks["benchmark_cluster_export_complete"] == "PASS"
    assert manifest_checks["benchmark_distance_summaries_complete"] == "PASS"
    connection = duckdb.connect(
        str(output / "duckdb/orthofinder_results.duckdb"),
        read_only=True,
    )
    try:
        assert connection.execute(
            "SELECT count(*) FROM benchmark_matched_controls"
        ).fetchone()[0] == 3
        assert connection.execute(
            "SELECT count(*) FROM benchmark_contrasts"
        ).fetchone()[0] > 0
    finally:
        connection.close()

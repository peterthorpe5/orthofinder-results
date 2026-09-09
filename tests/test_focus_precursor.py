"""Tests for complete E3-focus precursor selection and publication."""

from __future__ import annotations

import hashlib
from pathlib import Path

import duckdb
import pytest
from conftest import make_results

from orthofinder_interrogation_app.focus import FocusProtein, FocusProteinAuthority
from orthofinder_results.cli import main
from orthofinder_results.errors import InputValidationError
from orthofinder_results.focus_analysis import (
    FocusSelection,
    _matched_sequence_aliases,
    _matching_aliases,
    _selected_group_statistics,
    publish_focus_cluster_results,
    select_focus_groups,
)
from orthofinder_results.io_utils import read_tsv, sha256_file, write_tsv
from orthofinder_results.parsers import SEQUENCE_FIELDS
from orthofinder_results.pipeline import _validate_controls, run_pipeline
from orthofinder_results.statistics import GROUP_STATISTIC_FIELDS


def _precursor_arguments(
    *, results: Path, output: Path, focus_path: Path
) -> dict[str, object]:
    """Return a complete small E3 precursor pipeline invocation.

    Args:
        results: Synthetic OrthoFinder result directory.
        output: Formal resource output directory.
        focus_path: Custom seed authority.

    Returns:
        Named pipeline arguments.
    """

    return {
        "results_dir": results,
        "output_dir": output,
        "run_id": "e3-precursor-test",
        "work_dir": output.parent / "work",
        "alignment_dir": None,
        "distance_source": "RESOLVED_GENE_TREE",
        "distance_group_type": "HOG",
        "distance_hierarchy_node": "N0",
        "distance_max_groups": 0,
        "distance_max_members": 2,
        "parse_gene_trees": True,
        "report_max_statistic_rows": 100,
        "report_max_groups": 10,
        "report_max_members": 10,
        "report_nearest_neighbours": 2,
        "resume": False,
        "force": False,
        "verbose": False,
        "focus_proteins_path": focus_path,
        "focus_group_type": "HOG",
        "focus_hierarchy_node": "N0",
    }


def _add_singleton_focus_group(*, results: Path) -> None:
    """Add a focus-matched singleton without a resolved gene tree.

    Args:
        results: Synthetic OrthoFinder result directory to extend.
    """

    sequence_ids = results / "WorkingDirectory/SequenceIDs.txt"
    sequence_ids.write_text(
        sequence_ids.read_text(encoding="utf-8") + "0_2: singleE3 singleton\n",
        encoding="utf-8",
    )
    hogs = results / "Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
    hogs.write_text(
        hogs.read_text(encoding="utf-8")
        + "N0.HOG0000002\tOG0000002\tN0\tsingleE3\t\n",
        encoding="utf-8",
    )


def _write_focus_catalogue(*, path: Path) -> None:
    """Write a rich catalogue with matched and unmatched seeds.

    Args:
        path: New focus authority destination.
    """

    path.write_text(
        "seed_id\tassociated_seed_protein_names\tassociated_seed_categories\t"
        "associated_seed_review_statuses\tassociated_seed_ubiquitin_go_statuses\t"
        "associated_seed_organisms\tsequence_identifiers\tsequence_species\t"
        "protein_sequence_length\tannotation_scope\tcatalogue_source\n"
        "protA\tAlpha E3\tRING\treviewed\tUbiquitin GO term\tSpecies A\t"
        "protA\tSpecies_A\t4\texact seed\treviewed catalogue\n"
        "singleE3\tSingleton E3\tHECT\treviewed\tUbiquitin GO term\tSpecies A\t"
        "singleE3\tSpecies_A\t4\texact seed\treviewed catalogue\n"
        "ABSENT\tAbsent E3\tRBR\treviewed\tUbiquitin GO term\tSpecies Z\t"
        "ABSENT\tSpecies_Z\t4\texact seed\treviewed catalogue\n",
        encoding="utf-8",
    )


def test_e3_precursor_publishes_every_matched_cluster_and_distance_state(
    persistent_test_root: Path,
    tmp_path: Path,
) -> None:
    """Matched multi-member and singleton groups both reach the flat gzip export."""

    results = make_results(root=tmp_path)
    _add_singleton_focus_group(results=results)
    focus_path = tmp_path / "e3_catalogue.tsv"
    _write_focus_catalogue(path=focus_path)
    output = persistent_test_root / "e3_precursor"

    manifest = run_pipeline(
        **_precursor_arguments(
            results=results,
            output=output,
            focus_path=focus_path,
        )
    )

    assert manifest["status"] == "complete"
    assert manifest["counts"]["focus_seed_count"] == 3
    assert manifest["counts"]["matched_focus_seed_count"] == 2
    assert manifest["counts"]["focus_group_count"] == 2
    assert manifest["counts"]["distance_group_count"] == 2
    assert manifest["focus_analysis"]["cluster_results"] == (
        "tables/e3_cluster_results.tsv.gz"
    )
    results_rows = list(read_tsv(path=output / "tables/e3_cluster_results.tsv.gz"))
    assert [row["group_id"] for row in results_rows] == [
        "N0.HOG0000001",
        "N0.HOG0000002",
    ]
    multi, singleton = results_rows
    assert multi["matched_seed_ids"] == "protA"
    assert multi["computation_status"] == "DETERMINISTIC_MEMBER_SAMPLE"
    assert multi["sampled_member_count"] == "2"
    assert multi["distance_pair_count"] == "1"
    assert multi["maximum_distance"]
    assert multi["distance_sampling_fraction"] == str(2 / 3)
    assert singleton["matched_seed_ids"] == "singleE3"
    assert singleton["computation_status"] == "UNAVAILABLE"
    assert singleton["distance_pair_count"] == "0"
    assert singleton["mean_distance"] == ""
    assert "fewer than two" in singleton["failure_reason"]

    pair_rows = list(read_tsv(path=output / "tables/pairwise_distances.tsv.gz"))
    assert len(pair_rows) == 1
    assert "protA" in {pair_rows[0]["member_a"], pair_rows[0]["member_b"]}
    audit = {
        row["seed_id"]: row
        for row in read_tsv(path=output / "tables/e3_seed_catalogue_audit.tsv.gz")
    }
    assert audit["protA"]["match_status"] == "MATCHED"
    assert audit["singleE3"]["match_status"] == "MATCHED"
    assert audit["ABSENT"]["match_status"] == "UNMATCHED"
    assert len(list(read_tsv(path=output / "tables/e3_seed_matches.tsv.gz"))) == 2
    provenance = output / "provenance/focus_protein_authority.tsv"
    assert sha256_file(path=provenance) == sha256_file(path=focus_path)
    checks = {
        row["check_name"]: row["status"]
        for row in read_tsv(path=output / "qc/validation_checks.tsv")
    }
    assert checks["focus_cluster_export_complete"] == "PASS"
    assert checks["focus_distance_summary_complete"] == "PASS"

    connection = duckdb.connect(
        str(output / "duckdb/orthofinder_results.duckdb"), read_only=True
    )
    try:
        count, maximum = connection.execute(
            "SELECT count(*), max(maximum_distance) FROM e3_cluster_results"
        ).fetchone()
        assert count == 2
        assert maximum is not None
        assert connection.execute("SELECT count(*) FROM tree_payloads").fetchone()[0] == 1
    finally:
        connection.close()


def test_controlled_focus_aliases_are_exact_and_auditable() -> None:
    """Canonical, internal and structured UniProt aliases never use substring matching."""

    identifiers = frozenset({"sp|Q1|ENTRY", "Q1", "ENTRY", "0_1", "Q"})
    assert _matching_aliases(
        member_id="sp|Q1|ENTRY",
        internal_id="0_1",
        focus_identifiers=identifiers,
    ) == (
        ("0_1", "ORTHOFINDER_INTERNAL_ID"),
        ("ENTRY", "UNIPROT_ENTRY"),
        ("Q1", "UNIPROT_ACCESSION"),
        ("sp|Q1|ENTRY", "PROTEIN_ID"),
    )
    assert not any(
        identifier == "Q"
        for identifier, _ in _matching_aliases(
            member_id="sp|Q1|ENTRY",
            internal_id="0_1",
            focus_identifiers=identifiers,
        )
    )


def test_e3_precursor_cli_runs_all_groups_and_rejects_truncation(
    persistent_test_root: Path,
    tmp_path: Path,
) -> None:
    """The dedicated CLI fixes scientific scope while retaining a custom seed option."""

    results = make_results(root=tmp_path)
    focus_path = tmp_path / "focus.tsv"
    focus_path.write_text("protein_identifier\nprotA\n", encoding="utf-8")
    output = persistent_test_root / "cli_e3"
    base = [
        "--action",
        "e3-precursor",
        "--results-dir",
        str(results),
        "--output-dir",
        str(output),
        "--run-id",
        "cli-e3",
        "--work-dir",
        str(persistent_test_root / "cli_e3_work"),
        "--focus-proteins",
        str(focus_path),
        "--distance-max-members",
        "2",
    ]
    assert main(base) == 0
    assert (output / "tables/e3_cluster_results.tsv.gz").is_file()
    with pytest.raises(SystemExit):
        main([*base, "--distance-max-groups", "1"])
    with pytest.raises(SystemExit):
        main([*base, "--distance-source", "NONE"])


def test_focus_cluster_export_rejects_missing_or_duplicate_summaries(
    tmp_path: Path,
) -> None:
    """A selected cluster cannot disappear or gain an ambiguous distance result."""

    seed = FocusProtein(
        identifier="A",
        protein_name="Alpha",
        category="RING",
        evidence_type="reviewed",
        organism="Species",
        source="source",
        enabled=True,
        note="",
    )
    authority = FocusProteinAuthority(
        records=(seed,),
        source_name="focus.tsv",
        sha256=hashlib.sha256(b"focus").hexdigest(),
    )
    statistic = {
        "run_id": "run",
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "group_id": "N0.HOG1",
        "legacy_orthogroup_id": "OG1",
        "gene_tree_parent_clade": "N0",
        "member_count": "2",
        "species_count": "2",
        "single_copy_species_count": "2",
        "max_copies_per_species": "1",
        "mean_copies_per_species": "1",
        "is_singleton": "False",
        "species_labels": "A;B",
        "source_file": "N0.tsv",
    }
    selection = FocusSelection(
        authority=authority,
        group_type="HOG",
        hierarchy_node="N0",
        match_rows=(),
        audit_rows=(),
        group_statistics=(statistic,),
        required_members_by_group={},
    )
    with pytest.raises(InputValidationError, match="lack explicit"):
        publish_focus_cluster_results(
            tables_dir=tmp_path,
            selection=selection,
            distance_summaries=(),
        )
    summary = {
        "run_id": "run",
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "group_id": "N0.HOG1",
    }
    with pytest.raises(InputValidationError, match="multiple"):
        publish_focus_cluster_results(
            tables_dir=tmp_path,
            selection=selection,
            distance_summaries=(summary, summary),
        )


@pytest.mark.parametrize(
    ("group_type", "hierarchy_node", "message"),
    (
        ("UNSUPPORTED", "N0", "Unsupported focus group type"),
        ("LEGACY_ORTHOGROUP", "ROOT", "requires the ROOT hierarchy"),
    ),
)
def test_focus_selection_rejects_ambiguous_group_controls(
    tmp_path: Path,
    group_type: str,
    hierarchy_node: str,
    message: str,
) -> None:
    """Focus matching rejects an unsupported or ambiguous group collection."""

    authority = FocusProteinAuthority(
        records=(),
        source_name="focus.tsv",
        sha256=hashlib.sha256(b"focus").hexdigest(),
    )
    with pytest.raises(InputValidationError, match=message):
        select_focus_groups(
            tables_dir=tmp_path,
            run_id="run",
            authority=authority,
            group_type=group_type,
            hierarchy_node=hierarchy_node,
        )


def test_sequence_alias_index_is_run_scoped_and_omits_blank_internal_ids(
    tmp_path: Path,
) -> None:
    """Sequence aliases never leak between runs or invent blank identifiers."""

    path = tmp_path / "sequences.tsv.gz"
    common = {
        "species_index": "0",
        "species_label": "Species_A",
        "source_fasta": "Species_A.fa",
        "raw_header": "protein A",
        "source_file": "SequenceIDs.txt",
        "source_line": "1",
    }
    write_tsv(
        path=path,
        fieldnames=SEQUENCE_FIELDS,
        records=(
            {
                **common,
                "run_id": "run",
                "internal_id": "0_1",
                "member_id": "sp|Q1|ENTRY",
            },
            {
                **common,
                "run_id": "other-run",
                "internal_id": "0_2",
                "member_id": "ABSENT",
            },
            {
                **common,
                "run_id": "run",
                "internal_id": "",
                "member_id": "plain-id",
            },
        ),
    )

    matches, internal_ids = _matched_sequence_aliases(
        path=path,
        run_id="run",
        focus_identifiers=frozenset({"Q1", "ABSENT", "plain-id"}),
    )

    assert matches[("sp|Q1|ENTRY", "Species_A")] == (
        ("Q1", "UNIPROT_ACCESSION"),
    )
    assert matches[("plain-id", "Species_A")] == (("plain-id", "PROTEIN_ID"),)
    assert ("ABSENT", "Species_A") not in matches
    assert internal_ids == {("sp|Q1|ENTRY", "Species_A"): ("0_1",)}


def test_selected_group_statistics_rejects_duplicates_and_missing_rows(
    tmp_path: Path,
) -> None:
    """Selected focus memberships require one and only one statistics row."""

    path = tmp_path / "group_statistics.tsv.gz"
    row = {
        "run_id": "run",
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "group_id": "N0.HOG1",
        "legacy_orthogroup_id": "OG1",
        "gene_tree_parent_clade": "N0",
        "member_count": "2",
        "species_count": "2",
        "single_copy_species_count": "2",
        "max_copies_per_species": "1",
        "mean_copies_per_species": "1",
        "is_singleton": "False",
        "species_labels": "A;B",
        "source_file": "N0.tsv",
    }
    key = "HOG|N0|N0.HOG1"
    write_tsv(
        path=path,
        fieldnames=GROUP_STATISTIC_FIELDS,
        records=(row, row),
    )
    with pytest.raises(InputValidationError, match="Duplicate group-statistics"):
        _selected_group_statistics(path=path, selected_keys=frozenset({key}))

    write_tsv(path=path, fieldnames=GROUP_STATISTIC_FIELDS, records=(row,))
    with pytest.raises(InputValidationError, match="lack group statistics"):
        _selected_group_statistics(
            path=path,
            selected_keys=frozenset({key, "HOG|N0|N0.HOG2"}),
        )


@pytest.mark.parametrize(
    ("overrides", "message"),
    (
        ({"distance_source": "NONE"}, "require RESOLVED_GENE_TREE"),
        ({"focus_group_type": "BAD"}, "Unsupported focus_group_type"),
        (
            {
                "focus_group_type": "LEGACY_ORTHOGROUP",
                "focus_hierarchy_node": "ROOT",
            },
            "require an empty ROOT hierarchy",
        ),
        (
            {"distance_group_type": "LEGACY_ORTHOGROUP"},
            "same collection",
        ),
        ({"distance_hierarchy_node": "N1"}, "nodes must be identical"),
    ),
)
def test_focus_pipeline_controls_fail_closed(
    overrides: dict[str, object],
    message: str,
) -> None:
    """Focus runs reject scientific controls that would change their scope."""

    controls: dict[str, object] = {
        "run_id": "focus-run",
        "distance_source": "RESOLVED_GENE_TREE",
        "distance_group_type": "HOG",
        "distance_max_groups": 0,
        "distance_max_members": 250,
        "report_max_statistic_rows": 100,
        "report_max_groups": 25,
        "report_max_members": 250,
        "report_nearest_neighbours": 3,
        "resume": False,
        "force": False,
        "focus_enabled": True,
        "focus_group_type": "HOG",
        "focus_hierarchy_node": "N0",
        "distance_hierarchy_node": "N0",
    }
    controls.update(overrides)
    with pytest.raises(InputValidationError, match=message):
        _validate_controls(**controls)

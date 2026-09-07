"""Shared synthetic OrthoFinder result fixtures."""

from __future__ import annotations

import csv
import io
import json
import shutil
import uuid
from pathlib import Path

import duckdb
import pytest

from orthofinder_interrogation_app.taxonomy import TAXONOMY_COLUMNS


@pytest.fixture
def persistent_test_root() -> Path:
    """Provide a non-system-temporary root for publication policy tests."""

    root = Path.cwd() / ".pytest_work" / uuid.uuid4().hex
    root.mkdir(parents=True)
    try:
        yield root
    finally:
        shutil.rmtree(root)


def make_results(
    root: Path,
    *,
    version: str = "2.5.5",
    include_legacy: bool = True,
    include_alignments: bool = True,
) -> Path:
    """Create a minimal completed OrthoFinder result directory."""

    result = root / f"Results_{version.replace('.', '_')}"
    working = result / "WorkingDirectory"
    hogs = result / "Phylogenetic_Hierarchical_Orthogroups"
    species_trees = result / "Species_Tree"
    gene_trees = result / "Gene_Trees"
    resolved_trees = result / "Resolved_Gene_Trees"
    sequences = result / "Orthogroup_Sequences"
    for directory in (
        working,
        hogs,
        species_trees,
        gene_trees,
        resolved_trees,
        sequences,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    (result / "Log.txt").write_text(
        f"OrthoFinder version {version}\nCommand Line: orthofinder -f Example\n",
        encoding="utf-8",
    )
    (working / "SpeciesIDs.txt").write_text("0: Species_A.fa\n1: Species_B.faa\n", encoding="utf-8")
    (working / "SequenceIDs.txt").write_text(
        "0_0: protA first protein\n0_1: protA2 second protein\n1_0: protB\n",
        encoding="utf-8",
    )
    if include_legacy:
        orthogroups = result / "Orthogroups"
        orthogroups.mkdir()
        (orthogroups / "Orthogroups.tsv").write_text(
            "Orthogroup\tSpecies_A\tSpecies_B\nOG0000001\tprotA, protA2\tprotB\n",
            encoding="utf-8",
        )
    if version.startswith("3"):
        hog_header = "HOG\tGene Tree Parent Clade\tSpecies_A\tSpecies_B\n"
        hog_row = "N0.HOG0000001\tN0\tprotA, protA2\tprotB\n"
    else:
        hog_header = "HOG\tOG\tGene Tree Parent Clade\tSpecies_A\tSpecies_B\n"
        hog_row = "N0.HOG0000001\tOG0000001\tN0\tprotA, protA2\tprotB\n"
    (hogs / "N0.tsv").write_text(hog_header + hog_row, encoding="utf-8")
    (hogs / "N1.tsv").write_text(
        hog_header.replace("N0", "N1")
        + hog_row.replace("N0.HOG", "N1.HOG").replace("\tN0\t", "\tN1\t"),
        encoding="utf-8",
    )
    (species_trees / "SpeciesTree_rooted_node_labels.txt").write_text(
        "(Species_A:0.1,Species_B:0.2)N0:0.0;\n", encoding="utf-8"
    )
    tree = "((protA:0.1,protA2:0.2):0.1,protB:0.3)N0:0.0;\n"
    (gene_trees / "OG0000001_tree.txt").write_text(tree, encoding="utf-8")
    (resolved_trees / "OG0000001_tree.txt").write_text(tree, encoding="utf-8")
    (sequences / "OG0000001.fa").write_text(
        ">protA\nAAAA\n>protA2\nAATA\n>protB\nTTTT\n", encoding="utf-8"
    )
    if include_alignments:
        alignments = result / "MultipleSequenceAlignments"
        alignments.mkdir()
        (alignments / "N0.HOG0000001.fa").write_text(
            ">protA\nAAAA\n>protA2\nAATA\n>protB\nTTTT\n", encoding="utf-8"
        )
    return result


@pytest.fixture
def orthofinder2_results(tmp_path: Path) -> Path:
    """Return a compact OrthoFinder 2 fixture."""

    return make_results(tmp_path)


@pytest.fixture
def orthofinder3_results(tmp_path: Path) -> Path:
    """Return a compact HOG-only OrthoFinder 3 fixture."""

    return make_results(tmp_path, version="3.1.0", include_legacy=False)


@pytest.fixture
def application_resource(tmp_path: Path) -> Path:
    """Create a compact schema-2 resource for application tests."""

    resource = tmp_path / "application_resource"
    database_dir = resource / "duckdb"
    report_dir = resource / "report"
    database_dir.mkdir(parents=True)
    report_dir.mkdir()
    database = database_dir / "orthofinder_results.duckdb"
    connection = duckdb.connect(str(database))
    try:
        connection.execute(
            "CREATE TABLE resource_metadata AS SELECT '2026-09-07T00:00:00Z' "
            "AS created_at_utc, 2::INTEGER AS schema_version"
        )
        connection.execute(
            "CREATE TABLE species(run_id VARCHAR, species_index VARCHAR, "
            "species_label VARCHAR, source_fasta VARCHAR, source_file VARCHAR, "
            "source_line BIGINT)"
        )
        connection.executemany(
            "INSERT INTO species VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("test_run", str(index), species, f"{species}.fa", "SpeciesIDs.txt", index)
                for index, species in enumerate(
                    ("Species_A", "Species_B", "Species_C", "Species_D")
                )
            ],
        )
        connection.execute(
            "CREATE TABLE group_statistics(run_id VARCHAR, group_type VARCHAR, "
            "hierarchy_node VARCHAR, group_id VARCHAR, legacy_orthogroup_id VARCHAR, "
            "gene_tree_parent_clade VARCHAR, member_count BIGINT, species_count BIGINT, "
            "single_copy_species_count BIGINT, max_copies_per_species BIGINT, "
            "mean_copies_per_species DOUBLE, is_singleton BOOLEAN, species_labels VARCHAR, "
            "source_file VARCHAR)"
        )
        group_rows = [
            ("HOG", "N0", "N0.HOG1", "OG1", 3, 2, 1, 2, 1.5, False, "Species_A;Species_B"),
            ("HOG", "N0", "N0.HOG2", "OG2", 2, 2, 2, 1, 1.0, False, "Species_A;Species_C"),
            (
                "HOG",
                "N0",
                "N0.HOG3",
                "OG3",
                3,
                3,
                3,
                1,
                1.0,
                False,
                "Species_A;Species_B;Species_C",
            ),
            (
                "LEGACY_ORTHOGROUP",
                "",
                "OG4",
                "OG4",
                1,
                1,
                1,
                1,
                1.0,
                True,
                "Species_D",
            ),
        ]
        connection.executemany(
            "INSERT INTO group_statistics VALUES "
            "('test_run', ?, ?, ?, ?, 'N0', ?, ?, ?, ?, ?, ?, ?, 'groups.tsv')",
            group_rows,
        )
        connection.execute(
            "CREATE TABLE group_species_statistics(run_id VARCHAR, group_type VARCHAR, "
            "hierarchy_node VARCHAR, group_id VARCHAR, legacy_orthogroup_id VARCHAR, "
            "gene_tree_parent_clade VARCHAR, species_label VARCHAR, "
            "species_member_count BIGINT, member_fraction DOUBLE, source_file VARCHAR)"
        )
        species_rows = [
            ("HOG", "N0", "N0.HOG1", "OG1", "Species_A", 2, 2 / 3),
            ("HOG", "N0", "N0.HOG1", "OG1", "Species_B", 1, 1 / 3),
            ("HOG", "N0", "N0.HOG2", "OG2", "Species_A", 1, 0.5),
            ("HOG", "N0", "N0.HOG2", "OG2", "Species_C", 1, 0.5),
            ("HOG", "N0", "N0.HOG3", "OG3", "Species_A", 1, 1 / 3),
            ("HOG", "N0", "N0.HOG3", "OG3", "Species_B", 1, 1 / 3),
            ("HOG", "N0", "N0.HOG3", "OG3", "Species_C", 1, 1 / 3),
            ("LEGACY_ORTHOGROUP", "", "OG4", "OG4", "Species_D", 1, 1.0),
        ]
        connection.executemany(
            "INSERT INTO group_species_statistics VALUES "
            "('test_run', ?, ?, ?, ?, 'N0', ?, ?, ?, 'group_species.tsv')",
            species_rows,
        )
        membership_schema = (
            "(run_id VARCHAR, group_type VARCHAR, hierarchy_node VARCHAR, "
            "group_id VARCHAR, legacy_orthogroup_id VARCHAR, "
            "gene_tree_parent_clade VARCHAR, species_label VARCHAR, member_id VARCHAR, "
            "source_file VARCHAR, source_row BIGINT)"
        )
        connection.execute(f"CREATE TABLE hog_memberships{membership_schema}")
        connection.execute(f"CREATE TABLE legacy_orthogroup_memberships{membership_schema}")
        connection.executemany(
            "INSERT INTO hog_memberships VALUES "
            "('test_run', 'HOG', 'N0', ?, ?, 'N0', ?, ?, 'hogs.tsv', ?)",
            [
                ("N0.HOG1", "OG1", "Species_A", "alpha_1", 1),
                ("N0.HOG1", "OG1", "Species_A", "alpha_2", 1),
                ("N0.HOG1", "OG1", "Species_B", "beta_1", 1),
                ("N0.HOG2", "OG2", "Species_A", "literal%member", 2),
                ("N0.HOG2", "OG2", "Species_C", "gamma_1", 2),
                ("N0.HOG3", "OG3", "Species_A", "a3", 3),
                ("N0.HOG3", "OG3", "Species_B", "b3", 3),
                ("N0.HOG3", "OG3", "Species_C", "c3", 3),
            ],
        )
        connection.execute(
            "INSERT INTO legacy_orthogroup_memberships VALUES "
            "('test_run', 'LEGACY_ORTHOGROUP', '', 'OG4', 'OG4', 'N0', "
            "'Species_D', 'delta_1', 'orthogroups.tsv', 4)"
        )
        connection.execute(
            "CREATE TABLE distance_statistics(run_id VARCHAR, group_type VARCHAR, "
            "hierarchy_node VARCHAR, group_id VARCHAR, distance_method VARCHAR, "
            "computation_status VARCHAR, member_identifier_resolution VARCHAR, "
            "total_member_count BIGINT, sampled_member_count BIGINT, "
            "distance_pair_count BIGINT, unresolved_pair_count BIGINT, "
            "minimum_distance DOUBLE, q05_distance DOUBLE, q25_distance DOUBLE, "
            "median_distance DOUBLE, mean_distance DOUBLE, q75_distance DOUBLE, "
            "q95_distance DOUBLE, maximum_distance DOUBLE, "
            "population_stddev_distance DOUBLE, mean_comparable_sites DOUBLE, "
            "source_file VARCHAR, failure_reason VARCHAR)"
        )
        distance_rows = [
            ("N0.HOG1", 3, 3, 0.05, 0.1, 0.15, 0.02),
            ("N0.HOG3", 3, 3, 0.2, 0.8, 1.4, 0.4),
        ]
        connection.executemany(
            "INSERT INTO distance_statistics VALUES "
            "('test_run', 'HOG', 'N0', ?, 'patristic_branch_length', 'EXACT', "
            "'EXACT_MEMBER_ID', ?, ?, 3, 0, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0.0, "
            "'tree.txt', '')",
            [
                (
                    group_id,
                    total,
                    sampled,
                    minimum,
                    minimum,
                    minimum,
                    median,
                    mean,
                    maximum,
                    maximum,
                    maximum,
                    stddev,
                )
                for group_id, total, sampled, minimum, mean, maximum, stddev in distance_rows
                for median in ((minimum + maximum) / 2,)
            ],
        )
        connection.execute("CHECKPOINT")
    finally:
        connection.close()
    manifest = {
        "status": "complete",
        "run_id": "test_run",
        "package_version": "0.1.5",
        "schema_version": 2,
        "orthofinder_version": "2.5.5",
        "adapter_name": "orthofinder_2",
        "primary_group_authority": "HOG",
        "counts": {
            "group_count": 4,
            "group_species_statistic_count": 8,
            "species_count": 4,
            "distance_group_count": 2,
        },
    }
    (resource / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    visual_payload = _application_visual_payload()
    (report_dir / "orthofinder_results_summary.html").write_text(
        "<!doctype html><title>Test report</title>"
        '<script id="orthofinder-results-data" type="application/json">'
        + json.dumps(visual_payload, separators=(",", ":"))
        + "</script>",
        encoding="utf-8",
    )
    return resource


def _application_visual_payload() -> dict[str, object]:
    """Return one complete three-member report visual for application tests."""

    member_species = {
        "alpha_1": "Species_A",
        "alpha_2": "Species_A",
        "beta_1": "Species_B",
    }
    colours = {"Species_A": "#2c7fb8", "Species_B": "#d95f0e"}
    coordinates = {
        "alpha_1": (-0.1, 0.0),
        "alpha_2": (0.0, 0.1),
        "beta_1": (0.2, -0.1),
    }
    nodes = [
        {
            "id": member,
            "species": species,
            "speciesColour": colours[species],
            "projectionX": coordinates[member][0],
            "projectionY": coordinates[member][1],
            "isMedoid": member == "alpha_1",
        }
        for member, species in member_species.items()
    ]
    tree_nodes = [
        {
            "id": "root",
            "x": 0.0,
            "y": 1.0,
            "isLeaf": False,
            "memberId": "",
            "species": "",
            "colour": "#555555",
            "isMedoid": False,
        },
        *[
            {
                "id": f"leaf_{index}",
                "x": distance,
                "y": float(index),
                "isLeaf": True,
                "memberId": member,
                "species": member_species[member],
                "colour": colours[member_species[member]],
                "isMedoid": member == "alpha_1",
            }
            for index, (member, distance) in enumerate(
                (("alpha_1", 0.1), ("alpha_2", 0.2), ("beta_1", 0.3))
            )
        ],
    ]
    entry = {
        "label": "HOG | N0 | N0.HOG1",
        "nodes": nodes,
        "edges": [
            {
                "from": "alpha_1",
                "to": "alpha_2",
                "distance": 0.1,
                "edgeType": "NEAREST_NEIGHBOUR",
            },
            {
                "from": "alpha_2",
                "to": "beta_1",
                "distance": 0.3,
                "edgeType": "COMPONENT_CONNECTOR",
            },
        ],
        "members": [
            {"member_id": member, "species_label": species}
            for member, species in member_species.items()
        ],
        "medoid": {
            "status": "EXACT_WITHIN_RENDERED_SAMPLE",
            "member_id": "alpha_1",
            "mean_distance": 0.15,
            "sample_size": 3,
        },
        "distanceProjection": {
            "status": "COMPLETE_DISTANCE_PCOA_2D",
            "axis_1_positive_inertia_fraction": 0.55,
            "axis_2_positive_inertia_fraction": 0.30,
            "two_axis_positive_inertia_fraction": 0.85,
            "distance_correlation": 0.95,
            "normalised_stress": 0.1,
            "negative_inertia_fraction": 0.0,
            "quality_category": "BETTER",
            "quality_explanation": "Better diagnostic fit for this fixture.",
            "shepard_points": [[0.1, 0.09], [0.2, 0.21], [0.3, 0.29]],
            "shepard_point_count": 3,
            "shepard_total_pair_count": 3,
        },
        "phylogram": {
            "status": "COMPLETE_PRUNED_PHYLOGRAM",
            "treeId": "OG1",
            "displayedLeafCount": 3,
            "requestedMemberCount": 3,
            "maximumRootDistance": 0.3,
            "unresolvedMembers": [],
            "nodes": tree_nodes,
            "edges": [
                {"parentId": "root", "childId": f"leaf_{index}", "branchLength": value}
                for index, value in enumerate((0.1, 0.2, 0.3))
            ],
            "memberOrder": ["alpha_1", "alpha_2", "beta_1"],
        },
        "distanceMatrix": {
            "status": "EXACT_COMPLETE_DISPLAYED_MATRIX",
            "memberOrder": ["alpha_1", "alpha_2", "beta_1"],
            "orderMethod": "PRUNED_GENE_TREE_LEAF_ORDER",
            "pairCount": 3,
            "minimum": 0.1,
            "maximum": 0.3,
            "upperTriangle": [0.1, 0.2, 0.3],
        },
        "networkMetrics": {
            "rawComponentCount": 2,
            "rawIsolateCount": 0,
            "nearestNeighbourEdgeCount": 1,
            "connectorCount": 1,
        },
        "distanceSummary": {
            "distance_method": "patristic_branch_length",
            "computation_status": "EXACT",
            "total_member_count": 3,
            "sampled_member_count": 3,
            "distance_pair_count": 3,
            "mean_distance": 0.2,
            "population_stddev_distance": 0.08,
        },
    }
    return {
        "run": {"run_id": "test_run"},
        "networks": {"HOG|N0|N0.HOG1": entry},
        "limits": {"nearestNeighbours": 3},
    }


@pytest.fixture
def taxonomy_mapping_file(tmp_path: Path) -> Path:
    """Write a reviewed three-species mapping plus one ambiguous decision."""

    path = tmp_path / "reviewed_taxonomy.tsv"
    path.write_bytes(_application_taxonomy_mapping())
    return path


def _application_taxonomy_mapping() -> bytes:
    """Return a complete taxonomy mapping for the application fixture."""

    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=TAXONOMY_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    reviewed = (
        ("Species_A", "Species alpha", 101, 10, "Target clade"),
        ("Species_B", "Species beta", 102, 10, "Target clade"),
        ("Species_C", "Species gamma", 103, 20, "Outside clade"),
    )
    for label, accepted, taxon_id, parent_id, parent_name in reviewed:
        writer.writerow(
            {
                "workflow_species_label": label,
                "source_species_name": accepted,
                "accepted_species_name": accepted,
                "ncbi_taxon_id": taxon_id,
                "parent_taxon_id": parent_id,
                "parent_taxon_name": parent_name,
                "lineage_taxon_ids": f"1;{parent_id}",
                "lineage_names": f"cellular organisms;{parent_name}",
                "mapping_status": "REVIEWED",
                "mapping_method": "EXACT_SCIENTIFIC_NAME",
                "mapping_source": "NCBI Taxonomy",
                "source_date": "2026-09-07",
                "source_version": "fixture-2026-09-07",
                "reviewed_by": "test reviewer",
                "reviewed_at_utc": "2026-09-07T12:00:00Z",
                "review_note": "Synthetic fixture decision.",
            }
        )
    writer.writerow(
        {
            "workflow_species_label": "Species_D",
            "source_species_name": "Species delta",
            "mapping_status": "AMBIGUOUS",
            "review_note": "Two source taxa remain possible.",
        }
    )
    return output.getvalue().encode("utf-8")

"""Tests for checksum-bound portable Newick publication and decoding."""

from __future__ import annotations

import base64
import hashlib
import zlib
from pathlib import Path

import pytest

import orthofinder_results.trees as trees
from orthofinder_results.distances import (
    _binary_ancestor_table,
    _indexed_patristic_distance,
    calculate_patristic_distances_from_newick,
)
from orthofinder_results.errors import DistanceCalculationError, InputValidationError
from orthofinder_results.layout import discover_layout
from orthofinder_results.trees import (
    decode_newick_payload,
    encode_newick_payload,
    iter_portable_tree_payloads,
    iter_tree_inventory,
    normalise_newick_text,
    normalise_newick_tree,
    tree_id_from_path,
)


def _encoded(*, path: Path) -> dict[str, object]:
    """Encode one fixture tree with its current checksum."""

    raw = path.read_bytes()
    return encode_newick_payload(
        path=path,
        run_id="run",
        tree_type="RESOLVED_GENE_TREE",
        tree_id="OG1",
        group_id="OG1",
        source_path="Resolved_Gene_Trees/OG1_tree.txt",
        expected_sha256=hashlib.sha256(raw).hexdigest(),
    )


def _compressed(raw: bytes) -> str:
    """Return the portable base64/zlib representation of arbitrary bytes."""

    return base64.b64encode(zlib.compress(raw)).decode("ascii")


def test_inventory_prefers_resolved_payload_and_falls_back_to_gene_tree(
    orthofinder2_results: Path,
) -> None:
    """A resource stores one tree family, preferring resolved OrthoFinder trees."""

    layout = discover_layout(results_dir=orthofinder2_results)
    inventory = tuple(iter_tree_inventory(layout=layout, run_id="run"))
    assert {row["tree_type"] for row in inventory} == {
        "SPECIES_TREE",
        "GENE_TREE",
        "RESOLVED_GENE_TREE",
    }
    resolved = tuple(
        iter_portable_tree_payloads(layout=layout, run_id="run", inventory=inventory)
    )
    assert len(resolved) == 1
    assert resolved[0]["tree_type"] == "RESOLVED_GENE_TREE"
    assert decode_newick_payload(
        payload=str(resolved[0]["newick_payload"]),
        payload_encoding=str(resolved[0]["payload_encoding"]),
        expected_size=int(resolved[0]["source_size_bytes"]),
        expected_sha256=str(resolved[0]["source_sha256"]),
    ).endswith(";\n")

    fallback_inventory = tuple(
        row for row in inventory if row["tree_type"] != "RESOLVED_GENE_TREE"
    )
    original = tuple(
        iter_portable_tree_payloads(
            layout=layout,
            run_id="run",
            inventory=fallback_inventory,
        )
    )
    assert len(original) == 1
    assert original[0]["tree_type"] == "GENE_TREE"

    gene_tree = orthofinder2_results / "Gene_Trees" / "OG0000002_tree.txt"
    gene_tree.write_text("(protA:0.1,protB:0.2);\n", encoding="utf-8")
    mixed_inventory = tuple(iter_tree_inventory(layout=layout, run_id="run"))
    mixed = tuple(
        iter_portable_tree_payloads(
            layout=layout,
            run_id="run",
            inventory=mixed_inventory,
        )
    )
    assert [(row["tree_id"], row["tree_type"]) for row in mixed] == [
        ("OG0000001", "RESOLVED_GENE_TREE"),
        ("OG0000002", "GENE_TREE"),
    ]


def test_duplicate_portable_tree_authority_is_rejected(
    orthofinder2_results: Path,
) -> None:
    """Two same-type files resolving to one tree ID cannot be selected silently."""

    layout = discover_layout(results_dir=orthofinder2_results)
    inventory = list(iter_tree_inventory(layout=layout, run_id="run"))
    duplicate = next(row for row in inventory if row["tree_type"] == "GENE_TREE")
    inventory.append(dict(duplicate))
    with pytest.raises(InputValidationError, match="Duplicate GENE_TREE authority"):
        tuple(iter_portable_tree_payloads(layout=layout, run_id="run", inventory=inventory))


def test_portable_encoder_validates_source_size_utf8_and_checksum(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Only bounded, readable and unchanged UTF-8 tree bytes are embedded."""

    tree_path = tmp_path / "OG1_tree.txt"
    tree_path.write_text("(a:0.1,b:0.2);\n", encoding="utf-8")
    record = _encoded(path=tree_path)
    assert record["source_size_bytes"] == tree_path.stat().st_size
    assert record["payload_encoding"] == trees.TREE_PAYLOAD_ENCODING
    with pytest.raises(InputValidationError, match="unavailable"):
        encode_newick_payload(
            path=tmp_path / "missing.txt",
            run_id="run",
            tree_type="GENE_TREE",
            tree_id="OG1",
            group_id="OG1",
            source_path="missing.txt",
            expected_sha256="0" * 64,
        )
    empty = tmp_path / "empty.txt"
    empty.write_bytes(b"")
    with pytest.raises(InputValidationError, match="must contain"):
        _encoded(path=empty)
    invalid_utf8 = tmp_path / "invalid.txt"
    invalid_utf8.write_bytes(b"\xff")
    with pytest.raises(InputValidationError, match="not readable UTF-8"):
        _encoded(path=invalid_utf8)
    with pytest.raises(InputValidationError, match="checksum changed"):
        encode_newick_payload(
            path=tree_path,
            run_id="run",
            tree_type="GENE_TREE",
            tree_id="OG1",
            group_id="OG1",
            source_path="tree",
            expected_sha256="0" * 64,
        )
    monkeypatch.setattr(trees, "MAX_TREE_SOURCE_BYTES", 2)
    with pytest.raises(InputValidationError, match="must contain"):
        _encoded(path=tree_path)
    monkeypatch.setattr(trees, "MAX_TREE_PAYLOAD_CHARACTERS", 2)
    with pytest.raises(InputValidationError, match="encoded size"):
        decode_newick_payload(
            payload=_compressed(b"x"),
            payload_encoding=trees.TREE_PAYLOAD_ENCODING,
            expected_size=1,
            expected_sha256=hashlib.sha256(b"x").hexdigest(),
        )


@pytest.mark.parametrize(
    ("payload", "encoding", "size", "sha256", "message"),
    [
        (_compressed(b"tree"), "OTHER", 4, hashlib.sha256(b"tree").hexdigest(), "encoding"),
        (_compressed(b"tree"), trees.TREE_PAYLOAD_ENCODING, 0, "0" * 64, "unsafe"),
        ("", trees.TREE_PAYLOAD_ENCODING, 4, "0" * 64, "encoded size"),
        ("not-base64", trees.TREE_PAYLOAD_ENCODING, 4, "0" * 64, "could not be decoded"),
        (
            base64.b64encode(zlib.compress(b"tree") + zlib.compress(b"extra")).decode(),
            trees.TREE_PAYLOAD_ENCODING,
            4,
            hashlib.sha256(b"tree").hexdigest(),
            "trailing compressed data",
        ),
        (_compressed(b"tree"), trees.TREE_PAYLOAD_ENCODING, 5, "0" * 64, "size mismatch"),
        (_compressed(b"tree"), trees.TREE_PAYLOAD_ENCODING, 4, "0" * 64, "checksum"),
        (
            _compressed(b"\xff"),
            trees.TREE_PAYLOAD_ENCODING,
            1,
            hashlib.sha256(b"\xff").hexdigest(),
            "valid UTF-8",
        ),
        (
            _compressed(b" "),
            trees.TREE_PAYLOAD_ENCODING,
            1,
            hashlib.sha256(b" ").hexdigest(),
            "empty",
        ),
    ],
)
def test_portable_decoder_rejects_tampering_and_unsafe_data(
    payload: str,
    encoding: str,
    size: int,
    sha256: str,
    message: str,
) -> None:
    """Every portable payload boundary has a controlled provenance error."""

    with pytest.raises(InputValidationError, match=message):
        decode_newick_payload(
            payload=payload,
            payload_encoding=encoding,
            expected_size=size,
            expected_sha256=sha256,
        )


def test_newick_text_normalisation_and_filename_defences(tmp_path: Path) -> None:
    """In-memory and file trees share stable node identities and parse failures."""

    text = "((a:0.1,b:0.2)90:0.3,c:0.4)root;\n"
    nodes, edges = normalise_newick_text(
        newick_text=text,
        run_id="run",
        tree_type="GENE_TREE",
        tree_id="OG1",
    )
    assert len(nodes) == 5
    assert len(edges) == 4
    assert sum(bool(row["is_leaf"]) for row in nodes) == 3
    tree_path = tmp_path / "OG1_tree.newick"
    tree_path.write_text(text, encoding="utf-8")
    file_nodes, file_edges = normalise_newick_tree(
        path=tree_path,
        run_id="run",
        tree_type="GENE_TREE",
        tree_id="OG1",
    )
    assert file_nodes == nodes
    assert file_edges == edges
    assert tree_id_from_path(path=tree_path, tree_type="GENE_TREE") == "OG1"
    assert tree_id_from_path(path=tree_path, tree_type="SPECIES_TREE") == "SPECIES_TREE"
    with pytest.raises(ValueError, match="Unsupported"):
        tree_id_from_path(path=tree_path, tree_type="OTHER")
    with pytest.raises(InputValidationError, match="derive"):
        tree_id_from_path(path=tmp_path / "_tree.txt", tree_type="GENE_TREE")
    with pytest.raises(InputValidationError, match="Missing or empty"):
        normalise_newick_text(
            newick_text=" ",
            run_id="run",
            tree_type="GENE_TREE",
            tree_id="OG1",
        )
    with pytest.raises(InputValidationError, match="Could not parse"):
        normalise_newick_text(
            newick_text="(a,b",
            run_id="run",
            tree_type="GENE_TREE",
            tree_id="OG1",
        )
    invalid = tmp_path / "invalid_utf8.tree"
    invalid.write_bytes(b"\xff")
    with pytest.raises(InputValidationError, match="Could not read"):
        normalise_newick_tree(
            path=invalid,
            run_id="run",
            tree_type="GENE_TREE",
            tree_id="OG1",
        )


def test_portable_patristic_errors_and_binary_lifting_branches() -> None:
    """Portable calculation rejects parse, alias and branch-length ambiguities."""

    common = {
        "run_id": "run",
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "group_id": "N0.HOG1",
        "max_members": 10,
    }
    with pytest.raises(DistanceCalculationError, match="Missing or empty"):
        calculate_patristic_distances_from_newick(newick_text="", **common)
    with pytest.raises(DistanceCalculationError, match="Could not parse"):
        calculate_patristic_distances_from_newick(newick_text="(a,b", **common)
    with pytest.raises(DistanceCalculationError, match="ambiguous labels"):
        calculate_patristic_distances_from_newick(
            newick_text="(a:1,alias:1,b:1);",
            member_ids=("a", "b"),
            member_aliases={"a": {"alias": "ALIAS"}},
            **common,
        )
    with pytest.raises(DistanceCalculationError, match="invalid branch length"):
        calculate_patristic_distances_from_newick(
            newick_text="(a:-1,b:1);",
            **common,
        )
    assert _binary_ancestor_table(parent_indices=()) == ()
    ancestors = _binary_ancestor_table(parent_indices=(-1, 0, 1, 0))
    distance = _indexed_patristic_distance(
        left_index=3,
        right_index=2,
        levels=(0, 1, 2, 1),
        root_distances=(0.0, 1.0, 2.0, 1.0),
        ancestors=ancestors,
    )
    assert distance == pytest.approx(3.0)

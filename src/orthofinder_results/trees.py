"""Tree inventory and normalisation for OrthoFinder phylogenies."""

from __future__ import annotations

import base64
import hashlib
import io
import re
import zlib
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from .errors import InputValidationError
from .io_utils import file_record
from .models import ResultLayout

TREE_INVENTORY_FIELDS = (
    "run_id",
    "tree_type",
    "tree_id",
    "group_id",
    "path",
    "size_bytes",
    "sha256",
)
TREE_NODE_FIELDS = (
    "run_id",
    "tree_type",
    "tree_id",
    "node_id",
    "parent_node_id",
    "node_name",
    "is_leaf",
    "branch_length",
    "confidence",
    "descendant_leaf_count",
)
TREE_EDGE_FIELDS = (
    "run_id",
    "tree_type",
    "tree_id",
    "parent_node_id",
    "child_node_id",
    "branch_length",
)
TREE_PAYLOAD_FIELDS = (
    "run_id",
    "tree_type",
    "tree_id",
    "group_id",
    "source_path",
    "source_size_bytes",
    "source_sha256",
    "payload_encoding",
    "newick_payload",
)
TREE_PAYLOAD_ENCODING = "ZLIB_BASE64_UTF8"
MAX_TREE_SOURCE_BYTES = 64 * 1024 * 1024
MAX_TREE_PAYLOAD_CHARACTERS = 4 * ((MAX_TREE_SOURCE_BYTES + 64 * 1024 + 2) // 3)


def iter_tree_inventory(*, layout: ResultLayout, run_id: str) -> Iterator[dict[str, Any]]:
    """Yield checksum-bound records for all available tree files.

    Args:
        layout: Discovered OrthoFinder result layout.
        run_id: Immutable run identifier.

    Yields:
        Species, gene and resolved-gene tree file records.
    """

    sources: list[tuple[str, Path]] = []
    if layout.species_tree_path is not None:
        sources.append(("SPECIES_TREE", layout.species_tree_path))
    for tree_type, directory in (
        ("GENE_TREE", layout.gene_trees_dir),
        ("RESOLVED_GENE_TREE", layout.resolved_gene_trees_dir),
    ):
        if directory is None:
            continue
        sources.extend(
            (tree_type, path)
            for path in sorted(directory.iterdir())
            if path.is_file() and not path.name.startswith("._")
        )
    for tree_type, path in sources:
        inventory = file_record(path=path, relative_to=layout.results_dir)
        tree_id = tree_id_from_path(path=path, tree_type=tree_type)
        yield {
            "run_id": run_id,
            "tree_type": tree_type,
            "tree_id": tree_id,
            "group_id": "" if tree_type == "SPECIES_TREE" else tree_id,
            **inventory,
        }


def tree_id_from_path(*, path: Path, tree_type: str) -> str:
    """Derive a run-scoped tree identifier from an OrthoFinder filename.

    Args:
        path: Tree file path.
        tree_type: Declared tree type.

    Returns:
        ``SPECIES_TREE`` or the filename-derived group identifier.

    Raises:
        ValueError: If the tree type is unsupported.
    """

    if tree_type == "SPECIES_TREE":
        return "SPECIES_TREE"
    if tree_type not in {"GENE_TREE", "RESOLVED_GENE_TREE"}:
        raise ValueError(f"Unsupported tree_type: {tree_type!r}")
    name = path.name
    for suffix in (".txt", ".tree", ".nwk", ".newick"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = re.sub(r"_tree(?:_id)?$", "", name, flags=re.IGNORECASE)
    if not name:
        raise InputValidationError(f"Could not derive a tree identifier from {path}.")
    return name


def normalise_newick_tree(
    *, path: Path, run_id: str, tree_type: str, tree_id: str
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Convert one Newick tree into stable node and edge records.

    Args:
        path: Newick tree file.
        run_id: Immutable run identifier.
        tree_type: Species, gene or resolved-gene tree type.
        tree_id: Run-scoped tree identifier.

    Returns:
        Node records and directed parent-to-child edge records.

    Raises:
        InputValidationError: If Biopython cannot parse exactly one tree.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_file() or source.stat().st_size == 0:
        raise InputValidationError(f"Missing or empty Newick tree: {source}")
    try:
        newick_text = source.read_text(encoding="utf-8", errors="strict")
    except (OSError, UnicodeError) as error:
        raise InputValidationError(f"Could not read Newick tree {source}: {error}") from error
    return normalise_newick_text(
        newick_text=newick_text,
        run_id=run_id,
        tree_type=tree_type,
        tree_id=tree_id,
        source_label=str(source),
    )


def normalise_newick_text(
    *,
    newick_text: str,
    run_id: str,
    tree_type: str,
    tree_id: str,
    source_label: str = "portable tree payload",
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Convert Newick text into stable node and edge records.

    Args:
        newick_text: Complete UTF-8 Newick document.
        run_id: Immutable run identifier.
        tree_type: Species, gene or resolved-gene tree type.
        tree_id: Run-scoped tree identifier.
        source_label: Human-readable provenance used in controlled errors.

    Returns:
        Node records and directed parent-to-child edge records.

    Raises:
        InputValidationError: If the text is empty or cannot be parsed.
    """

    from Bio import Phylo
    from Bio.Phylo.NewickIO import NewickError

    if not newick_text.strip():
        raise InputValidationError(f"Missing or empty Newick tree: {source_label}")
    try:
        tree = Phylo.read(io.StringIO(newick_text), "newick")
    except (NewickError, ValueError, OSError) as error:
        raise InputValidationError(
            f"Could not parse Newick tree {source_label}: {error}"
        ) from error

    clades = list(tree.find_clades(order="preorder"))
    node_ids = {id(clade): f"{tree_id}:n{index}" for index, clade in enumerate(clades)}
    parents: dict[int, str] = {}
    for parent in clades:
        for child in parent.clades:
            parents[id(child)] = node_ids[id(parent)]
    descendant_counts: dict[int, int] = {}
    for clade in tree.find_clades(order="postorder"):
        descendant_counts[id(clade)] = (
            1
            if clade.is_terminal()
            else sum(descendant_counts[id(child)] for child in clade.clades)
        )

    nodes: list[dict[str, Any]] = []
    edges: list[dict[str, Any]] = []
    for clade in clades:
        node_id = node_ids[id(clade)]
        parent_id = parents.get(id(clade), "")
        branch_length = clade.branch_length
        confidence = clade.confidence
        nodes.append(
            {
                "run_id": run_id,
                "tree_type": tree_type,
                "tree_id": tree_id,
                "node_id": node_id,
                "parent_node_id": parent_id,
                "node_name": clade.name or "",
                "is_leaf": clade.is_terminal(),
                "branch_length": "" if branch_length is None else float(branch_length),
                "confidence": "" if confidence is None else float(confidence),
                "descendant_leaf_count": descendant_counts[id(clade)],
            }
        )
        if parent_id:
            edges.append(
                {
                    "run_id": run_id,
                    "tree_type": tree_type,
                    "tree_id": tree_id,
                    "parent_node_id": parent_id,
                    "child_node_id": node_id,
                    "branch_length": "" if branch_length is None else float(branch_length),
                }
            )
    return nodes, edges


def iter_portable_tree_payloads(
    *,
    layout: ResultLayout,
    run_id: str,
    inventory: Iterable[dict[str, Any]],
    selected_tree_ids: frozenset[str] | None = None,
) -> Iterator[dict[str, Any]]:
    """Yield compressed portable gene-tree payloads for lazy app calculations.

    Resolved gene-tree records are preferred independently for each tree identifier.
    An original gene tree fills any identifier absent from the resolved directory,
    avoiding redundant copies without losing partially resolved datasets.

    Args:
        layout: Discovered OrthoFinder result layout.
        run_id: Immutable run identifier.
        inventory: Checksum-bound tree inventory records.
        selected_tree_ids: Optional exact identifiers restricting portable payloads.

    Yields:
        One compressed, checksum-bound Newick payload per selected gene tree.
    """

    selected: dict[str, dict[str, Any]] = {}
    seen_authorities: set[tuple[str, str]] = set()
    for record in inventory:
        tree_type = str(record.get("tree_type", ""))
        if tree_type not in {"GENE_TREE", "RESOLVED_GENE_TREE"}:
            continue
        tree_id = str(record.get("tree_id", ""))
        if selected_tree_ids is not None and tree_id not in selected_tree_ids:
            continue
        authority = (tree_type, tree_id)
        if authority in seen_authorities:
            raise InputValidationError(
                f"Duplicate {tree_type} authority for portable tree {tree_id!r}."
            )
        seen_authorities.add(authority)
        current = selected.get(tree_id)
        if current is None or tree_type == "RESOLVED_GENE_TREE":
            selected[tree_id] = record
    for tree_id in sorted(selected):
        record = selected[tree_id]
        tree_type = str(record["tree_type"])
        source = layout.results_dir / str(record["path"])
        yield encode_newick_payload(
            path=source,
            run_id=run_id,
            tree_type=tree_type,
            tree_id=tree_id,
            group_id=str(record["group_id"]),
            source_path=str(record["path"]),
            expected_sha256=str(record["sha256"]),
        )


def encode_newick_payload(
    *,
    path: Path,
    run_id: str,
    tree_type: str,
    tree_id: str,
    group_id: str,
    source_path: str,
    expected_sha256: str,
) -> dict[str, Any]:
    """Encode one checksum-verified Newick file for portable publication.

    Args:
        path: Source Newick file.
        run_id: Immutable run identifier.
        tree_type: Gene-tree authority type.
        tree_id: Run-scoped tree identifier.
        group_id: OrthoFinder group associated with the tree.
        source_path: Portable source-relative provenance path.
        expected_sha256: Digest recorded by the tree inventory.

    Returns:
        Schema-complete compressed payload record.

    Raises:
        InputValidationError: If size, encoding or checksum validation fails.
    """

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise InputValidationError(f"Portable tree source is unavailable: {source}")
    size = source.stat().st_size
    if not 1 <= size <= MAX_TREE_SOURCE_BYTES:
        raise InputValidationError(
            f"Portable tree source must contain 1–{MAX_TREE_SOURCE_BYTES:,} bytes; "
            f"observed {size:,}: {source}"
        )
    try:
        raw = source.read_bytes()
        raw.decode("utf-8", errors="strict")
    except (OSError, UnicodeError) as error:
        raise InputValidationError(f"Portable tree is not readable UTF-8: {source}") from error
    observed_sha256 = hashlib.sha256(raw).hexdigest()
    if observed_sha256 != expected_sha256:
        raise InputValidationError(
            f"Portable tree checksum changed during publication: {source}"
        )
    payload = base64.b64encode(zlib.compress(raw, level=9)).decode("ascii")
    return {
        "run_id": run_id,
        "tree_type": tree_type,
        "tree_id": tree_id,
        "group_id": group_id,
        "source_path": source_path,
        "source_size_bytes": size,
        "source_sha256": observed_sha256,
        "payload_encoding": TREE_PAYLOAD_ENCODING,
        "newick_payload": payload,
    }


def decode_newick_payload(
    *,
    payload: str,
    payload_encoding: str,
    expected_size: int,
    expected_sha256: str,
) -> str:
    """Decode and validate one bounded portable Newick payload.

    Args:
        payload: Base64-encoded compressed bytes.
        payload_encoding: Declared encoding contract.
        expected_size: Uncompressed byte count from publication.
        expected_sha256: SHA-256 digest of the uncompressed source.

    Returns:
        Verified UTF-8 Newick text.

    Raises:
        InputValidationError: If decoding, bounds or provenance validation fails.
    """

    if payload_encoding != TREE_PAYLOAD_ENCODING:
        raise InputValidationError(f"Unsupported portable tree encoding: {payload_encoding!r}.")
    if not 1 <= expected_size <= MAX_TREE_SOURCE_BYTES:
        raise InputValidationError(
            f"Portable tree declares an unsafe source size: {expected_size:,}."
        )
    if not 1 <= len(payload) <= MAX_TREE_PAYLOAD_CHARACTERS:
        raise InputValidationError("Portable tree payload has an unsafe encoded size.")
    try:
        compressed = base64.b64decode(payload.encode("ascii"), validate=True)
        decompressor = zlib.decompressobj()
        raw = decompressor.decompress(compressed, expected_size + 1)
    except (UnicodeError, ValueError, zlib.error) as error:
        raise InputValidationError("Portable tree payload could not be decoded.") from error
    if decompressor.unused_data or decompressor.unconsumed_tail or not decompressor.eof:
        raise InputValidationError("Portable tree payload contains trailing compressed data.")
    if len(raw) != expected_size:
        raise InputValidationError(
            f"Portable tree size mismatch: expected {expected_size:,}, observed {len(raw):,}."
        )
    if hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise InputValidationError("Portable tree payload checksum validation failed.")
    try:
        text = raw.decode("utf-8", errors="strict")
    except UnicodeError as error:
        raise InputValidationError("Portable tree payload is not valid UTF-8.") from error
    if not text.strip():
        raise InputValidationError("Portable tree payload is empty.")
    return text

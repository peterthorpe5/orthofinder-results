"""Deterministic minimal taxonomy coverage-tree construction."""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Iterable

from orthofinder_results.errors import InputValidationError

from .taxonomy_selection import (
    ExpectedTaxaAuthority,
    SelectionSpec,
    TaxonomyGraph,
)

_LOGGER = logging.getLogger("orthofinder_interrogation_app.coverage_tree")
MAX_COVERAGE_TREE_NODES = 5_000


@dataclass(frozen=True)
class TaxonCoverage:
    """Independent selection and dataset-coverage states for one displayed node."""

    taxon_id: str
    taxon_name: str
    taxon_rank: str
    parent_taxon_id: str | None
    lineage_depth: int
    taxonomy_authority: str
    taxonomy_release: str
    selection_state: str
    coverage_state: str
    is_selected: bool
    is_excluded: bool
    is_represented: bool
    is_expected: bool
    has_represented_descendant: bool
    has_expected_no_data_descendant: bool
    is_outside_selected_scope: bool
    represented_input_labels: tuple[str, ...]
    expected_reasons: tuple[str, ...]
    expected_sources: tuple[str, ...]
    represented_terminal_count: int
    expected_no_data_terminal_count: int
    excluded_terminal_count: int

    def as_record(self) -> dict[str, Any]:
        """Return a stable ``taxon_coverage`` export record."""

        return {
            "taxon_id": self.taxon_id,
            "taxon_name": self.taxon_name,
            "taxon_rank": self.taxon_rank,
            "parent_taxon_id": self.parent_taxon_id or "",
            "lineage_depth": self.lineage_depth,
            "taxonomy_authority": self.taxonomy_authority,
            "taxonomy_release": self.taxonomy_release,
            "selection_state": self.selection_state,
            "coverage_state": self.coverage_state,
            "is_selected": self.is_selected,
            "is_excluded": self.is_excluded,
            "is_represented": self.is_represented,
            "is_expected": self.is_expected,
            "has_represented_descendant": self.has_represented_descendant,
            "has_expected_no_data_descendant": self.has_expected_no_data_descendant,
            "is_outside_selected_scope": self.is_outside_selected_scope,
            "represented_input_labels": ";".join(self.represented_input_labels),
            "expected_reasons": ";".join(self.expected_reasons),
            "expected_sources": ";".join(self.expected_sources),
            "represented_terminal_count": self.represented_terminal_count,
            "expected_no_data_terminal_count": self.expected_no_data_terminal_count,
            "excluded_terminal_count": self.excluded_terminal_count,
        }


@dataclass(frozen=True)
class TreeEdge:
    """One parent-to-child edge retained by the pruned displayed tree."""

    parent_taxon_id: str
    child_taxon_id: str

    def as_record(self) -> dict[str, str]:
        """Return a stable ``tree_edges`` export record."""

        return {
            "parent_taxon_id": self.parent_taxon_id,
            "child_taxon_id": self.child_taxon_id,
        }


@dataclass(frozen=True)
class CoverageTree:
    """Reconciled nodes and edges for one displayed selection/coverage state."""

    root_taxon_id: str
    nodes: tuple[TaxonCoverage, ...]
    edges: tuple[TreeEdge, ...]
    compact: bool

    @property
    def node_by_id(self) -> dict[str, TaxonCoverage]:
        """Return displayed nodes indexed by stable identifier."""

        return {node.taxon_id: node for node in self.nodes}

    @property
    def children_by_id(self) -> dict[str, tuple[str, ...]]:
        """Return deterministically ordered displayed children."""

        children: dict[str, list[str]] = defaultdict(list)
        node_by_id = self.node_by_id
        for edge in self.edges:
            children[edge.parent_taxon_id].append(edge.child_taxon_id)
        return {
            parent: tuple(
                sorted(
                    values,
                    key=lambda value: (
                        node_by_id[value].taxon_name.casefold(),
                        node_by_id[value].taxon_name,
                        value,
                    ),
                )
            )
            for parent, values in children.items()
        }


def build_coverage_tree(
    *,
    graph: TaxonomyGraph,
    selection: SelectionSpec,
    expected: ExpectedTaxaAuthority,
    outside_hit_taxon_ids: Iterable[str] = (),
    compact: bool = False,
    maximum_nodes: int = MAX_COVERAGE_TREE_NODES,
) -> CoverageTree:
    """Build the minimal ancestor closure for selected and covered terminals.

    Args:
        graph: Complete reviewed bounded taxonomy graph.
        selection: Validated taxonomy predicates.
        expected: Explicit expected-taxon universe.
        outside_hit_taxon_ids: Represented taxa observed outside include scope.
        compact: Collapse unary neutral internal nodes when true.
        maximum_nodes: Defensive displayed-node bound.

    Returns:
        Reconciled pruned tree with independent selection and coverage states.

    Raises:
        InputValidationError: If identifiers or bounds cannot be reconciled.
    """

    if not isinstance(maximum_nodes, int) or isinstance(maximum_nodes, bool):
        raise InputValidationError("maximum_nodes must be an integer.")
    if not 1 <= maximum_nodes <= MAX_COVERAGE_TREE_NODES:
        raise InputValidationError(
            f"maximum_nodes must be between 1 and {MAX_COVERAGE_TREE_NODES:,}."
        )
    expected_ids = set(expected.included_taxon_ids)
    unknown_expected = expected_ids.difference(graph.terminal_ids)
    if unknown_expected:
        raise InputValidationError(
            "Expected universe contains unknown or non-terminal taxa: "
            + "; ".join(sorted(unknown_expected))
        )
    outside_hits = set(outside_hit_taxon_ids)
    if not outside_hits.issubset(graph.represented_taxon_ids):
        raise InputValidationError(
            "Outside-scope hit identifiers must be represented reviewed input taxa."
        )
    terminals = set(graph.represented_taxon_ids).union(expected_ids)
    clade_types = {"INCLUDE_CLADE", "ONLY_IN_CLADE", "EXCLUDE_CLADE"}
    for predicate in selection.predicates:
        if predicate.selector_type in clade_types:
            terminals.update(
                graph.descendants(taxon_id=predicate.taxon_id, terminals_only=True)
            )
        else:
            terminals.add(predicate.taxon_id)
    if not terminals:
        raise InputValidationError("The displayed terminal universe is empty.")
    paths = {taxon_id: graph.ancestors(taxon_id=taxon_id) for taxon_id in terminals}
    root = _least_common_root(paths=tuple(paths.values()))
    displayed: set[str] = set()
    for path in paths.values():
        root_index = path.index(root)
        displayed.update(path[root_index:])
    if len(displayed) > maximum_nodes:
        raise InputValidationError(
            f"Selection coverage tree requires {len(displayed):,} nodes; "
            f"limit is {maximum_nodes:,}. Narrow the reviewed expected universe."
        )
    parents = {
        taxon_id: (
            None if taxon_id == root else graph.node_by_id[taxon_id].parent_taxon_id
        )
        for taxon_id in displayed
    }
    if compact:
        parents, displayed = _collapse_unary_neutral_nodes(
            graph=graph,
            selection=selection,
            root_taxon_id=root,
            displayed=displayed,
            parents=parents,
        )
    represented = set(graph.represented_taxon_ids)
    expected_no_data = expected_ids.difference(represented)
    excluded = _excluded_terminals(graph=graph, selection=selection)
    direct_states = _direct_selection_states(selection=selection)
    label_mappings: dict[str, list[str]] = defaultdict(list)
    for mapping in graph.input_mappings:
        if (
            mapping.workflow_species_label in graph.represented_labels
            and mapping.mapping_status == "REVIEWED"
            and mapping.taxon_id
        ):
            label_mappings[mapping.taxon_id].append(mapping.workflow_species_label)
    expected_rows = {row.taxon_id: row for row in expected.records if row.included}
    displayed_terminals = terminals.intersection(displayed)
    rows: list[TaxonCoverage] = []
    for taxon_id in sorted(
        displayed,
        key=lambda value: (
            graph.node_by_id[value].lineage_depth,
            graph.node_by_id[value].name.casefold(),
            graph.node_by_id[value].name,
            value,
        ),
    ):
        source = graph.node_by_id[taxon_id]
        descendants = set(graph.descendants(taxon_id=taxon_id, terminals_only=True))
        descendants.intersection_update(displayed_terminals)
        represented_descendants = descendants.intersection(represented)
        expected_no_data_descendants = descendants.intersection(expected_no_data)
        excluded_descendants = descendants.intersection(excluded)
        is_represented = taxon_id in represented
        is_expected = taxon_id in expected_ids
        coverage_state = _coverage_state(
            taxon_id=taxon_id,
            represented=represented,
            expected_no_data=expected_no_data,
            outside_hits=outside_hits,
            represented_descendants=represented_descendants,
        )
        selection_state, is_selected, is_excluded = _selection_state(
            graph=graph,
            selection=selection,
            taxon_id=taxon_id,
            direct_states=direct_states,
        )
        expected_row = expected_rows.get(taxon_id)
        rows.append(
            TaxonCoverage(
                taxon_id=taxon_id,
                taxon_name=source.name,
                taxon_rank=source.rank,
                parent_taxon_id=parents[taxon_id],
                lineage_depth=source.lineage_depth,
                taxonomy_authority=source.authority,
                taxonomy_release=source.release,
                selection_state=selection_state,
                coverage_state=coverage_state,
                is_selected=is_selected,
                is_excluded=is_excluded,
                is_represented=is_represented,
                is_expected=is_expected,
                has_represented_descendant=bool(represented_descendants),
                has_expected_no_data_descendant=bool(expected_no_data_descendants),
                is_outside_selected_scope=taxon_id in outside_hits,
                represented_input_labels=tuple(sorted(label_mappings[taxon_id])),
                expected_reasons=(expected_row.reason,) if expected_row is not None else (),
                expected_sources=(expected_row.source,) if expected_row is not None else (),
                represented_terminal_count=len(represented_descendants),
                expected_no_data_terminal_count=len(expected_no_data_descendants),
                excluded_terminal_count=len(excluded_descendants),
            )
        )
    edges = tuple(
        TreeEdge(parent_taxon_id=parent, child_taxon_id=child)
        for child, parent in sorted(
            parents.items(),
            key=lambda item: (
                "" if item[1] is None else item[1],
                graph.node_by_id[item[0]].name.casefold(),
                item[0],
            ),
        )
        if parent is not None
    )
    tree = CoverageTree(root_taxon_id=root, nodes=tuple(rows), edges=edges, compact=compact)
    _reconcile_tree(tree=tree, displayed_terminals=displayed_terminals)
    _LOGGER.info(
        "Selection coverage tree built: root=%s, nodes=%s, edges=%s, compact=%s",
        root,
        len(tree.nodes),
        len(tree.edges),
        compact,
    )
    return tree


def tree_edge_records(*, tree: CoverageTree) -> tuple[dict[str, str], ...]:
    """Return the complete derived ``tree_edges`` table."""

    return tuple(edge.as_record() for edge in tree.edges)


def coverage_records(*, tree: CoverageTree) -> tuple[dict[str, Any], ...]:
    """Return the complete derived ``taxon_coverage`` table."""

    return tuple(node.as_record() for node in tree.nodes)


def tree_to_newick(*, tree: CoverageTree) -> str:
    """Serialise the displayed taxonomy tree as deterministic escaped Newick."""

    nodes = tree.node_by_id
    children = tree.children_by_id

    def serialise(taxon_id: str) -> str:
        """Serialise one node and its displayed descendants."""

        child_text = ""
        if children.get(taxon_id):
            child_text = "(" + ",".join(
                serialise(child_id) for child_id in children[taxon_id]
            ) + ")"
        node = nodes[taxon_id]
        label = f"{node.taxon_name} [{node.taxon_rank}; ID={node.taxon_id}]"
        return child_text + "'" + label.replace("'", "''") + "'"

    return serialise(tree.root_taxon_id) + ";\n"


def _least_common_root(*, paths: tuple[tuple[str, ...], ...]) -> str:
    """Return the deepest common identifier across root-to-terminal paths."""

    if not paths or any(not path for path in paths):
        raise InputValidationError("Cannot identify a displayed root from empty paths.")
    common: list[str] = []
    for values in zip(*paths):
        if len(set(values)) != 1:
            break
        common.append(values[0])
    if not common:
        raise InputValidationError("Displayed taxonomy terminals have no common root.")
    return common[-1]


def _collapse_unary_neutral_nodes(
    *,
    graph: TaxonomyGraph,
    selection: SelectionSpec,
    root_taxon_id: str,
    displayed: set[str],
    parents: dict[str, str | None],
) -> tuple[dict[str, str | None], set[str]]:
    """Collapse only unary, unselected internal nodes in labelled compact mode."""

    selected_ids = {row.taxon_id for row in selection.predicates}
    retained = set(displayed)
    changed = True
    while changed:
        changed = False
        children: dict[str, list[str]] = defaultdict(list)
        for child, parent in parents.items():
            if child in retained and parent is not None and parent in retained:
                children[parent].append(child)
        candidates = sorted(
            (
                taxon_id
                for taxon_id in retained
                if taxon_id != root_taxon_id
                and taxon_id not in selected_ids
                and taxon_id not in graph.terminal_ids
                and len(children.get(taxon_id, ())) == 1
            ),
            key=lambda value: (-graph.node_by_id[value].lineage_depth, value),
        )
        if not candidates:
            continue
        for taxon_id in candidates:
            child_ids = [
                child
                for child, parent in parents.items()
                if child in retained and parent == taxon_id
            ]
            if len(child_ids) != 1:
                continue
            child = child_ids[0]
            parents[child] = parents[taxon_id]
            retained.remove(taxon_id)
            changed = True
    return {key: value for key, value in parents.items() if key in retained}, retained


def _direct_selection_states(*, selection: SelectionSpec) -> dict[str, tuple[str, ...]]:
    """Return all direct selection states indexed by taxon identifier."""

    names = {
        "REQUIRED_EXACT": "REQUIRED_EXACT",
        "INCLUDE_CLADE": "INCLUDE_CLADE",
        "ONLY_IN_CLADE": "ONLY_IN_SCOPE",
        "EXCLUDE_EXACT": "EXCLUDED_EXACT",
        "EXCLUDE_CLADE": "EXCLUDED_CLADE",
    }
    states: dict[str, list[str]] = defaultdict(list)
    for row in selection.predicates:
        states[row.taxon_id].append(names[row.selector_type])
    return {key: tuple(values) for key, values in states.items()}


def _selection_state(
    *,
    graph: TaxonomyGraph,
    selection: SelectionSpec,
    taxon_id: str,
    direct_states: dict[str, tuple[str, ...]],
) -> tuple[str, bool, bool]:
    """Resolve direct and inherited display selection state for one node."""

    if taxon_id in direct_states:
        states = direct_states[taxon_id]
        return (
            "|".join(states),
            any(not state.startswith("EXCLUDED") for state in states),
            any(state.startswith("EXCLUDED") for state in states),
        )
    excluded_clades = selection.identifiers(selector_type="EXCLUDE_CLADE")
    if any(
        graph.is_descendant(taxon_id=taxon_id, ancestor_taxon_id=ancestor)
        for ancestor in excluded_clades
    ):
        return "EXCLUDED_BY_CLADE", False, True
    only_clades = selection.identifiers(selector_type="ONLY_IN_CLADE")
    if only_clades and all(
        graph.is_descendant(taxon_id=taxon_id, ancestor_taxon_id=ancestor)
        for ancestor in only_clades
    ):
        return "WITHIN_ONLY_IN_SCOPE", False, False
    return "NEUTRAL", False, False


def _coverage_state(
    *,
    taxon_id: str,
    represented: set[str],
    expected_no_data: set[str],
    outside_hits: set[str],
    represented_descendants: set[str],
) -> str:
    """Return one primary coverage display state without changing audit Booleans."""

    if taxon_id in outside_hits:
        return "OUTSIDE_SELECTED_SCOPE_WITH_HITS"
    if taxon_id in represented:
        return "REPRESENTED_IN_INPUT"
    if taxon_id in expected_no_data:
        return "EXPECTED_NO_DATA"
    if represented_descendants:
        return "ANCESTOR_OF_REPRESENTED"
    return "NOT_IN_EXPECTED_UNIVERSE"


def _excluded_terminals(*, graph: TaxonomyGraph, selection: SelectionSpec) -> set[str]:
    """Return bounded terminals covered by exact or clade exclusions."""

    excluded = set(selection.identifiers(selector_type="EXCLUDE_EXACT"))
    for taxon_id in selection.identifiers(selector_type="EXCLUDE_CLADE"):
        excluded.update(graph.descendants(taxon_id=taxon_id, terminals_only=True))
    return excluded


def _reconcile_tree(*, tree: CoverageTree, displayed_terminals: set[str]) -> None:
    """Verify topology, terminal inventory and all descendant audit counts."""

    node_ids = {node.taxon_id for node in tree.nodes}
    if tree.root_taxon_id not in node_ids:
        raise InputValidationError("Coverage tree root is missing from displayed nodes.")
    if len(tree.edges) != len(node_ids) - 1:
        raise InputValidationError("Coverage tree is disconnected or contains duplicate edges.")
    edge_pairs = {(edge.parent_taxon_id, edge.child_taxon_id) for edge in tree.edges}
    if len(edge_pairs) != len(tree.edges):
        raise InputValidationError("Coverage tree contains duplicate edges.")
    if any(parent not in node_ids or child not in node_ids for parent, child in edge_pairs):
        raise InputValidationError("Coverage tree edge references an undisplayed node.")
    node_by_id = tree.node_by_id
    if any(
        node_by_id[child].parent_taxon_id != parent
        for parent, child in edge_pairs
    ):
        raise InputValidationError("Coverage tree edges disagree with node parent identifiers.")
    roots = [node.taxon_id for node in tree.nodes if node.parent_taxon_id is None]
    if roots != [tree.root_taxon_id]:
        raise InputValidationError("Coverage tree must contain exactly one displayed root.")
    if not displayed_terminals.issubset(node_ids):
        raise InputValidationError("Coverage tree omitted a required displayed terminal.")
    children = tree.children_by_id
    visited: set[str] = set()

    terminal_sets: dict[str, set[str]] = {}

    def terminal_descendants(taxon_id: str, active: frozenset[str]) -> set[str]:
        """Return reviewed terminals while detecting cycles and repeated visits."""

        if taxon_id in active:
            raise InputValidationError("Coverage tree contains a directed cycle.")
        if taxon_id in visited:
            raise InputValidationError("Coverage tree contains a multiply parented node.")
        visited.add(taxon_id)
        child_ids = children.get(taxon_id, ())
        next_active = active.union({taxon_id})
        result = {taxon_id} if taxon_id in displayed_terminals else set()
        for child_id in child_ids:
            result.update(terminal_descendants(child_id, next_active))
        terminal_sets[taxon_id] = result
        return result

    observed_terminals = terminal_descendants(tree.root_taxon_id, frozenset())
    if visited != node_ids:
        raise InputValidationError("Coverage tree contains nodes disconnected from the root.")
    if observed_terminals != displayed_terminals:
        raise InputValidationError("Coverage tree terminal inventory does not reconcile.")
    for node in tree.nodes:
        terminals = terminal_sets[node.taxon_id]
        represented_count = sum(node_by_id[value].is_represented for value in terminals)
        expected_no_data_count = sum(
            node_by_id[value].coverage_state == "EXPECTED_NO_DATA"
            for value in terminals
        )
        excluded_count = sum(node_by_id[value].is_excluded for value in terminals)
        if (
            node.represented_terminal_count != represented_count
            or node.expected_no_data_terminal_count != expected_no_data_count
            or node.excluded_terminal_count != excluded_count
        ):
            raise InputValidationError(
                f"Coverage tree descendant counts do not reconcile for taxon {node.taxon_id}."
            )

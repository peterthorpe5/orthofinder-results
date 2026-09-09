"""Offline taxonomy selection semantics and auditable group evaluation."""

from __future__ import annotations

import csv
import hashlib
import io
import logging
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from orthofinder_results.errors import InputValidationError

from .taxonomy import TaxonomyAuthority, TaxonomyRecord

_LOGGER = logging.getLogger("orthofinder_interrogation_app.taxonomy_selection")
SELECTION_SEMANTICS_VERSION = "1.0"
PREDICATE_TYPES = (
    "REQUIRED_EXACT",
    "INCLUDE_CLADE",
    "ONLY_IN_CLADE",
    "EXCLUDE_EXACT",
    "EXCLUDE_CLADE",
)
MAX_EXPECTED_TAXA_BYTES = 5 * 1024 * 1024
MAX_EXPECTED_TAXA = 100_000
MAX_SELECTION_GROUPS = 250_000
_SAFE_TEXT = re.compile(r"^[^\x00-\x1f\x7f]{1,2048}$")


@dataclass(frozen=True)
class TaxonomyNode:
    """One normalised node in a pinned, reviewed taxonomy graph."""

    taxon_id: str
    name: str
    rank: str
    parent_taxon_id: str | None
    lineage_depth: int
    authority: str
    release: str
    is_terminal: bool


@dataclass(frozen=True)
class InputTaxonMapping:
    """One exact workflow label and its reviewed or unresolved mapping state."""

    workflow_species_label: str
    taxon_id: str
    mapping_status: str
    mapping_method: str
    mapping_reason: str
    source_name_original: str
    authority: str
    release: str
    role: str
    source_species_name: str = ""
    accepted_species_name: str = ""
    ncbi_taxon_id: int | None = None
    parent_taxon_id: int | None = None
    parent_taxon_name: str = ""
    lineage_taxon_ids: tuple[int, ...] = ()
    lineage_names: tuple[str, ...] = ()
    lineage_ranks: tuple[str, ...] = ()
    taxon_rank: str = ""
    authority_taxon_id: str = ""
    mapping_source: str = ""
    source_date: str = ""
    source_version: str = ""
    reviewed_by: str = ""
    reviewed_at_utc: str = ""
    review_note: str = ""


@dataclass(frozen=True)
class TaxonomyGraph:
    """Validated taxonomy nodes plus exact workflow-label mappings."""

    nodes: tuple[TaxonomyNode, ...]
    input_mappings: tuple[InputTaxonMapping, ...]
    represented_labels: tuple[str, ...]
    authority: str
    release: str
    mapping_sha256: str

    @cached_property
    def node_by_id(self) -> dict[str, TaxonomyNode]:
        """Return nodes indexed by stable authority taxon identifier."""

        return {node.taxon_id: node for node in self.nodes}

    @cached_property
    def children_by_id(self) -> dict[str, tuple[str, ...]]:
        """Return deterministically ordered direct children for every node."""

        children: dict[str, list[str]] = defaultdict(list)
        for node in self.nodes:
            if node.parent_taxon_id is not None:
                children[node.parent_taxon_id].append(node.taxon_id)
        result: dict[str, tuple[str, ...]] = {}
        for parent, identifiers in children.items():
            result[parent] = tuple(
                sorted(
                    identifiers,
                    key=lambda value: (
                        self.node_by_id[value].name.casefold(),
                        self.node_by_id[value].name,
                        value,
                    ),
                )
            )
        return result

    @cached_property
    def mapping_by_label(self) -> dict[str, InputTaxonMapping]:
        """Return exact workflow-label mappings, including unresolved rows."""

        return {row.workflow_species_label: row for row in self.input_mappings}

    @cached_property
    def terminal_ids(self) -> frozenset[str]:
        """Return all reviewed taxa that are terminals in the bounded universe."""

        return frozenset(node.taxon_id for node in self.nodes if node.is_terminal)

    @cached_property
    def represented_taxon_ids(self) -> frozenset[str]:
        """Return reviewed terminal taxa represented in the OrthoFinder input."""

        represented = set(self.represented_labels)
        return frozenset(
            row.taxon_id
            for row in self.input_mappings
            if row.workflow_species_label in represented
            and row.mapping_status == "REVIEWED"
            and row.taxon_id
        )

    def descendants(self, *, taxon_id: str, terminals_only: bool = False) -> frozenset[str]:
        """Return a node and all descendants in the reviewed bounded graph.

        Args:
            taxon_id: Stable selected-authority identifier.
            terminals_only: Retain only reviewed bounded terminal taxa.

        Returns:
            Descendant identifiers, including ``taxon_id`` when eligible.

        Raises:
            InputValidationError: If the identifier is unknown.
        """

        if taxon_id not in self.node_by_id:
            raise InputValidationError(f"Unknown taxonomy identifier: {taxon_id}")
        pending = [taxon_id]
        observed: set[str] = set()
        while pending:
            current = pending.pop()
            if current in observed:
                continue
            observed.add(current)
            pending.extend(self.children_by_id.get(current, ()))
        if terminals_only:
            observed.intersection_update(self.terminal_ids)
        return frozenset(observed)

    def ancestors(self, *, taxon_id: str) -> tuple[str, ...]:
        """Return root-to-self ancestors for one known node."""

        if taxon_id not in self.node_by_id:
            raise InputValidationError(f"Unknown taxonomy identifier: {taxon_id}")
        reversed_path: list[str] = []
        current: str | None = taxon_id
        while current is not None:
            if current in reversed_path:
                raise InputValidationError("Taxonomy graph contains a parent cycle.")
            reversed_path.append(current)
            current = self.node_by_id[current].parent_taxon_id
        return tuple(reversed(reversed_path))

    def is_descendant(self, *, taxon_id: str, ancestor_taxon_id: str) -> bool:
        """Return whether one known node is at or below another."""

        return ancestor_taxon_id in self.ancestors(taxon_id=taxon_id)


@dataclass(frozen=True)
class ExpectedTaxon:
    """One explicitly expected terminal taxon in the dataset coverage universe."""

    taxon_id: str
    reason: str
    source: str
    included: bool = True


@dataclass(frozen=True)
class ExpectedTaxaAuthority:
    """Versioned, bounded expected-taxon universe with checksum provenance."""

    records: tuple[ExpectedTaxon, ...]
    source_name: str
    sha256: str

    @property
    def included_taxon_ids(self) -> tuple[str, ...]:
        """Return deterministic included terminal identifiers."""

        return tuple(row.taxon_id for row in self.records if row.included)


@dataclass(frozen=True)
class SelectionPredicate:
    """One validated, ordered taxonomy predicate."""

    selector_type: str
    taxon_id: str
    taxon_name: str
    taxon_rank: str
    order: int
    logical_mode: str = "AND"


@dataclass(frozen=True)
class SelectionSpec:
    """Normalised taxonomy predicates with stable AND semantics."""

    predicates: tuple[SelectionPredicate, ...]
    semantics_version: str = SELECTION_SEMANTICS_VERSION

    def identifiers(self, *, selector_type: str) -> tuple[str, ...]:
        """Return identifiers for one exact selector type."""

        if selector_type not in PREDICATE_TYPES:
            raise InputValidationError(f"Unknown selector type: {selector_type}")
        return tuple(
            row.taxon_id
            for row in self.predicates
            if row.selector_type == selector_type
        )


@dataclass(frozen=True)
class GroupTaxonEvaluation:
    """Compact audit of all taxonomy predicates for one OrthoFinder group."""

    run_id: str
    group_type: str
    hierarchy_node: str
    group_id: str
    mapped_member_count: int
    unmapped_member_count: int
    represented_taxon_ids: tuple[str, ...]
    outside_selected_scope_taxon_ids: tuple[str, ...]
    unmapped_species_labels: tuple[str, ...]
    predicate_pass: bool
    predicate_failure_reasons: tuple[str, ...]
    predicate_audit: tuple[str, ...]
    selection_manifest_id: str = ""

    def as_record(self) -> dict[str, Any]:
        """Return the stable bridge/export record for this group."""

        return {
            "run_id": self.run_id,
            "group_type": self.group_type,
            "hierarchy_node": self.hierarchy_node,
            "group_id": self.group_id,
            "mapped_member_count": self.mapped_member_count,
            "unmapped_member_count": self.unmapped_member_count,
            "represented_taxon_ids": ";".join(self.represented_taxon_ids),
            "outside_selected_scope_taxon_ids": ";".join(
                self.outside_selected_scope_taxon_ids
            ),
            "unmapped_species_labels": ";".join(self.unmapped_species_labels),
            "predicate_pass": self.predicate_pass,
            "predicate_failure_reasons": ";".join(self.predicate_failure_reasons),
            "predicate_audit": ";".join(self.predicate_audit),
            "selection_manifest_id": self.selection_manifest_id,
        }


def build_taxonomy_graph(*, authority: TaxonomyAuthority) -> TaxonomyGraph:
    """Build and reconcile an offline graph from reviewed mapping lineages.

    Args:
        authority: Mapping parsed and validated against exact resource labels.

    Returns:
        Consistent nodes and complete input-mapping inventory.

    Raises:
        InputValidationError: If names, ranks, parents, depths or releases conflict.
    """

    reviewed = authority.all_reviewed_records
    if not reviewed:
        raise InputValidationError("A selection coverage tree requires reviewed mappings.")
    mutable: dict[str, dict[str, Any]] = {}
    terminal_ids: set[str] = set()
    authority_releases: dict[str, set[str]] = defaultdict(set)
    for record in sorted(reviewed, key=lambda row: row.workflow_species_label):
        authority_releases[record.taxonomy_authority].add(record.taxonomy_release)
        parent: str | None = None
        for depth, (taxon_id, name, rank) in enumerate(
            zip(
                record.lineage_taxon_ids,
                record.lineage_names,
                record.lineage_ranks,
                strict=True,
            )
        ):
            key = str(taxon_id)
            _merge_node(
                nodes=mutable,
                taxon_id=key,
                name=name,
                rank=rank,
                parent_taxon_id=parent,
                depth=depth,
                authority="NCBI Taxonomy",
                release=record.taxonomy_release,
            )
            parent = key
        terminal_key = record.taxon_key
        if not terminal_key:
            raise InputValidationError(
                f"Reviewed label {record.workflow_species_label!r} lacks a terminal ID."
            )
        _merge_node(
            nodes=mutable,
            taxon_id=terminal_key,
            name=record.accepted_species_name,
            rank=record.taxon_rank,
            parent_taxon_id=parent,
            depth=len(record.lineage_taxon_ids),
            authority=(
                "NCBI Taxonomy"
                if record.ncbi_taxon_id is not None
                else record.taxonomy_authority
            ),
            release=record.taxonomy_release,
        )
        terminal_ids.add(terminal_key)
    inconsistent_releases = {
        name: releases
        for name, releases in authority_releases.items()
        if len(releases) > 1
    }
    if inconsistent_releases:
        name = sorted(inconsistent_releases)[0]
        raise InputValidationError(
            f"Taxonomy authority {name!r} mixes releases: "
            + "; ".join(sorted(inconsistent_releases[name]))
        )
    roots = sorted(key for key, row in mutable.items() if row["parent_taxon_id"] is None)
    if len(roots) != 1:
        raise InputValidationError(
            "Reviewed taxonomy lineages must share one displayed root; observed: "
            + "; ".join(roots)
        )
    nodes = tuple(
        TaxonomyNode(
            taxon_id=taxon_id,
            name=str(row["name"]),
            rank=str(row["rank"]),
            parent_taxon_id=row["parent_taxon_id"],
            lineage_depth=int(row["depth"]),
            authority=str(row["authority"]),
            release=str(row["release"]),
            is_terminal=taxon_id in terminal_ids,
        )
        for taxon_id, row in sorted(
            mutable.items(),
            key=lambda item: (
                int(item[1]["depth"]),
                str(item[1]["name"]).casefold(),
                item[0],
            ),
        )
    )
    records_by_label = {row.workflow_species_label: row for row in authority.records}
    mappings: list[InputTaxonMapping] = []
    all_labels = sorted(set(records_by_label).union(authority.expected_species))
    for label in all_labels:
        record = records_by_label.get(label)
        mappings.append(_mapping_record(label=label, record=record))
    graph_authorities = tuple(sorted(authority_releases))
    graph_releases = tuple(
        sorted({release for values in authority_releases.values() for release in values})
    )
    graph = TaxonomyGraph(
        nodes=nodes,
        input_mappings=tuple(mappings),
        represented_labels=authority.expected_species,
        authority="; ".join(graph_authorities),
        release="; ".join(graph_releases),
        mapping_sha256=authority.mapping_sha256,
    )
    for node in graph.nodes:
        graph.ancestors(taxon_id=node.taxon_id)
    _LOGGER.info(
        "Taxonomy graph built: nodes=%s, terminals=%s, represented=%s, unresolved=%s",
        len(graph.nodes),
        len(graph.terminal_ids),
        len(graph.represented_taxon_ids),
        sum(row.mapping_status != "REVIEWED" for row in graph.input_mappings),
    )
    return graph


def default_expected_taxa(*, graph: TaxonomyGraph) -> ExpectedTaxaAuthority:
    """Use reviewed OrthoFinder input taxa as the explicit default expectation."""

    records = tuple(
        ExpectedTaxon(
            taxon_id=taxon_id,
            reason="Exact reviewed species list supplied to this OrthoFinder run.",
            source="OrthoFinder input species authority",
            included=True,
        )
        for taxon_id in sorted(graph.represented_taxon_ids)
    )
    payload = expected_taxa_to_tsv(authority=ExpectedTaxaAuthority(records, "", ""))
    return ExpectedTaxaAuthority(
        records=records,
        source_name="resource_input_species",
        sha256=hashlib.sha256(payload).hexdigest(),
    )


def read_expected_taxa(*, path: Path, graph: TaxonomyGraph) -> ExpectedTaxaAuthority:
    """Read a bounded expected-taxon TSV from a user-selected path."""

    source = Path(path).expanduser().resolve()
    try:
        data = source.read_bytes()
    except OSError as error:
        raise InputValidationError(f"Expected-taxon TSV is unavailable: {source}") from error
    authority = parse_expected_taxa(data=data, source_name=source.name, graph=graph)
    _LOGGER.info(
        "Expected-taxon universe loaded: source=%s, included=%s, sha256=%s",
        source,
        len(authority.included_taxon_ids),
        authority.sha256,
    )
    return authority


def parse_expected_taxa(
    *, data: bytes, source_name: str, graph: TaxonomyGraph
) -> ExpectedTaxaAuthority:
    """Parse and validate an explicit expected-taxon universe.

    Args:
        data: UTF-8 TSV containing taxon_id, reason, source and included.
        source_name: Human-readable authority filename.
        graph: Reviewed bounded taxonomy graph.

    Returns:
        Deterministically ordered expected-taxon authority.

    Raises:
        InputValidationError: If rows are malformed, unknown or non-terminal.
    """

    if not data or len(data) > MAX_EXPECTED_TAXA_BYTES or b"\x00" in data:
        raise InputValidationError("Expected-taxon TSV is empty, unsafe or overlong.")
    try:
        text = data.decode("utf-8-sig")
    except UnicodeError as error:
        raise InputValidationError("Expected-taxon TSV must be valid UTF-8.") from error
    reader = csv.DictReader(io.StringIO(text), delimiter="\t", strict=True)
    headings = tuple(reader.fieldnames or ())
    required = ("taxon_id", "reason", "source", "included")
    missing = tuple(column for column in required if column not in headings)
    if missing or len(headings) != len(set(headings)):
        raise InputValidationError(
            "Expected-taxon TSV requires unique columns: " + "; ".join(required)
        )
    records: list[ExpectedTaxon] = []
    try:
        for line_number, raw_row in enumerate(reader, start=2):
            if None in raw_row:
                raise InputValidationError(
                    f"Expected-taxon TSV line {line_number} has extra fields."
                )
            row = {key: str(value or "").strip() for key, value in raw_row.items()}
            taxon_id = _normalise_taxon_id(value=row["taxon_id"])
            reason, source = row["reason"], row["source"]
            if not reason or not source or any(
                _SAFE_TEXT.fullmatch(value) is None for value in (reason, source)
            ):
                raise InputValidationError(
                    f"Expected-taxon TSV line {line_number} has missing or unsafe provenance."
                )
            included = _parse_boolean(value=row["included"], line_number=line_number)
            if taxon_id not in graph.node_by_id:
                raise InputValidationError(f"Unknown expected taxon identifier: {taxon_id}")
            if taxon_id not in graph.terminal_ids:
                raise InputValidationError(
                    f"Expected taxon {taxon_id} is not a reviewed terminal taxon."
                )
            records.append(ExpectedTaxon(taxon_id, reason, source, included))
    except csv.Error as error:
        raise InputValidationError(f"Expected-taxon TSV could not be parsed: {error}") from error
    if not records or len(records) > MAX_EXPECTED_TAXA:
        raise InputValidationError(
            f"Expected-taxon TSV must contain 1–{MAX_EXPECTED_TAXA:,} rows."
        )
    identifiers = [row.taxon_id for row in records]
    if len(identifiers) != len(set(identifiers)):
        raise InputValidationError("Expected-taxon TSV contains duplicate taxon IDs.")
    if not any(row.included for row in records):
        raise InputValidationError("Expected-taxon universe has no included taxa.")
    return ExpectedTaxaAuthority(
        records=tuple(sorted(records, key=lambda row: row.taxon_id)),
        source_name=source_name.strip() or "expected_taxa.tsv",
        sha256=hashlib.sha256(data).hexdigest(),
    )


def expected_taxa_to_tsv(*, authority: ExpectedTaxaAuthority) -> bytes:
    """Serialise an expected-taxon authority as deterministic TSV."""

    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=("taxon_id", "reason", "source", "included"),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for row in authority.records:
        writer.writerow(
            {
                "taxon_id": row.taxon_id,
                "reason": row.reason,
                "source": row.source,
                "included": str(row.included).lower(),
            }
        )
    return output.getvalue().encode("utf-8")


def make_selection(
    *,
    graph: TaxonomyGraph,
    required_exact: Sequence[str] = (),
    include_clade: Sequence[str] = (),
    only_in_clade: Sequence[str] = (),
    exclude_exact: Sequence[str] = (),
    exclude_clade: Sequence[str] = (),
) -> SelectionSpec:
    """Normalise selectors and reject contradictions before group evaluation."""

    supplied = {
        "REQUIRED_EXACT": required_exact,
        "INCLUDE_CLADE": include_clade,
        "ONLY_IN_CLADE": only_in_clade,
        "EXCLUDE_EXACT": exclude_exact,
        "EXCLUDE_CLADE": exclude_clade,
    }
    predicates: list[SelectionPredicate] = []
    order = 0
    for selector_type in PREDICATE_TYPES:
        values = _normalise_identifiers(values=supplied[selector_type])
        for taxon_id in values:
            node = graph.node_by_id.get(taxon_id)
            if node is None:
                raise InputValidationError(f"Unknown taxonomy identifier: {taxon_id}")
            if selector_type in {"REQUIRED_EXACT", "EXCLUDE_EXACT"} and not node.is_terminal:
                raise InputValidationError(
                    f"{selector_type} requires a reviewed terminal taxon: {taxon_id}"
                )
            predicates.append(
                SelectionPredicate(
                    selector_type=selector_type,
                    taxon_id=taxon_id,
                    taxon_name=node.name,
                    taxon_rank=node.rank,
                    order=order,
                )
            )
            order += 1
    selection = SelectionSpec(predicates=tuple(predicates))
    _validate_selection_contradictions(graph=graph, selection=selection)
    _LOGGER.info("Taxonomy selection validated: predicates=%s", len(predicates))
    return selection


def selection_summary(*, selection: SelectionSpec) -> tuple[str, ...]:
    """Return synchronised plain-English statements for selected predicates."""

    descriptions = {
        "REQUIRED_EXACT": "Require exact taxon",
        "INCLUDE_CLADE": "Require at least one member from clade",
        "ONLY_IN_CLADE": "Restrict every mapped member to clade",
        "EXCLUDE_EXACT": "Reject exact taxon",
        "EXCLUDE_CLADE": "Reject every descendant of clade",
    }
    if not selection.predicates:
        return ("No taxonomic predicates are active; all evaluated groups pass.",)
    return tuple(
        f"{descriptions[row.selector_type]} {row.taxon_name} "
        f"({row.taxon_rank}, ID {row.taxon_id})."
        for row in selection.predicates
    )


def selection_predicate_records(*, selection: SelectionSpec) -> tuple[dict[str, Any], ...]:
    """Return stable derived ``selection_predicates`` table rows."""

    return tuple(
        {
            "selector_type": row.selector_type,
            "taxon_id": row.taxon_id,
            "taxon_name": row.taxon_name,
            "taxon_rank": row.taxon_rank,
            "predicate_order": row.order,
            "normalised_logical_mode": row.logical_mode,
            "semantics_version": selection.semantics_version,
        }
        for row in selection.predicates
    )


def input_taxon_mapping_records(*, graph: TaxonomyGraph) -> tuple[dict[str, Any], ...]:
    """Return the complete derived input-mapping table."""

    represented = set(graph.represented_labels)
    return tuple(
        {
            "workflow_species_label": row.workflow_species_label,
            "source_species_name": row.source_species_name,
            "accepted_species_name": row.accepted_species_name,
            "ncbi_taxon_id": row.ncbi_taxon_id if row.ncbi_taxon_id is not None else "",
            "parent_taxon_id": (
                row.parent_taxon_id if row.parent_taxon_id is not None else ""
            ),
            "parent_taxon_name": row.parent_taxon_name,
            "lineage_taxon_ids": ";".join(str(value) for value in row.lineage_taxon_ids),
            "lineage_names": ";".join(row.lineage_names),
            "lineage_ranks": ";".join(row.lineage_ranks),
            "taxon_rank": row.taxon_rank,
            "mapping_status": (
                "UNMAPPED" if row.mapping_status == "MISSING" else row.mapping_status
            ),
            "derived_mapping_status": row.mapping_status,
            "role": row.role,
            "reviewed_terminal_taxon_id": row.taxon_id,
            "authority_taxon_id": row.authority_taxon_id,
            "mapping_method": row.mapping_method,
            "mapping_reason": row.mapping_reason,
            "mapping_source": row.mapping_source,
            "source_date": row.source_date,
            "source_version": row.source_version,
            "reviewed_by": row.reviewed_by,
            "reviewed_at_utc": row.reviewed_at_utc,
            "review_note": row.review_note,
            "source_name_original": row.source_name_original,
            "taxonomy_authority": row.authority,
            "taxonomy_release": row.release,
            "represented_in_input": row.workflow_species_label in represented,
        }
        for row in graph.input_mappings
    )


def unmapped_input_records(
    *, graph: TaxonomyGraph, group_species_rows: Sequence[Mapping[str, Any]] = ()
) -> tuple[dict[str, Any], ...]:
    """Return unresolved exact input labels without attaching them to the tree."""

    occurrences: Counter[str] = Counter()
    for row in group_species_rows:
        label = str(row.get("species_label", ""))
        count = _positive_count(value=row.get("species_member_count", 0), allow_zero=True)
        occurrences[label] += count
    represented = set(graph.represented_labels)
    rows = []
    for mapping in graph.input_mappings:
        if (
            mapping.workflow_species_label in represented
            and mapping.mapping_status != "REVIEWED"
        ):
            rows.append(
                {
                    "workflow_species_label": mapping.workflow_species_label,
                    "input_label_count": 1,
                    "evaluated_group_member_occurrences": occurrences[
                        mapping.workflow_species_label
                    ],
                    "mapping_status": mapping.mapping_status,
                    "mapping_reason": mapping.mapping_reason,
                    "source_name_original": mapping.source_name_original,
                }
            )
    return tuple(sorted(rows, key=lambda row: str(row["workflow_species_label"])))


def evaluate_group_taxonomy(
    *,
    graph: TaxonomyGraph,
    selection: SelectionSpec,
    run_id: str,
    group_type: str,
    hierarchy_node: str,
    group_id: str,
    species_counts: Mapping[str, int],
) -> GroupTaxonEvaluation:
    """Evaluate all AND-composed selectors for one immutable group composition."""

    if not all(isinstance(value, str) for value in (run_id, group_type, hierarchy_node, group_id)):
        raise InputValidationError("Group identity fields must be text.")
    if not group_id or not group_type:
        raise InputValidationError("Group evaluation requires group type and identifier.")
    mapped_counts: Counter[str] = Counter()
    unmapped_counts: Counter[str] = Counter()
    for species_label, raw_count in species_counts.items():
        if not isinstance(species_label, str) or not species_label:
            raise InputValidationError("Group species labels must be non-empty text.")
        count = _positive_count(value=raw_count)
        mapping = graph.mapping_by_label.get(species_label)
        if mapping is not None and mapping.mapping_status == "REVIEWED" and mapping.taxon_id:
            mapped_counts[mapping.taxon_id] += count
        else:
            unmapped_counts[species_label] += count
    present = frozenset(mapped_counts)
    failures: list[str] = []
    audit: list[str] = []
    for predicate in selection.predicates:
        matched = _predicate_matches_group(
            graph=graph,
            predicate=predicate,
            present_taxa=present,
            has_unmapped=bool(unmapped_counts),
            selection=selection,
        )
        audit.append(
            f"{predicate.selector_type}:{predicate.taxon_id}="
            f"{'PASS' if matched else 'FAIL'}"
        )
        if not matched:
            failures.append(
                _predicate_failure(
                    predicate=predicate,
                    has_unmapped=bool(unmapped_counts),
                )
            )
    outside = _outside_scope_taxa(graph=graph, selection=selection, present_taxa=present)
    return GroupTaxonEvaluation(
        run_id=run_id,
        group_type=group_type,
        hierarchy_node=hierarchy_node,
        group_id=group_id,
        mapped_member_count=sum(mapped_counts.values()),
        unmapped_member_count=sum(unmapped_counts.values()),
        represented_taxon_ids=tuple(sorted(present)),
        outside_selected_scope_taxon_ids=tuple(sorted(outside)),
        unmapped_species_labels=tuple(sorted(unmapped_counts)),
        predicate_pass=not failures,
        predicate_failure_reasons=tuple(failures),
        predicate_audit=tuple(audit),
    )


def evaluate_group_rows(
    *,
    graph: TaxonomyGraph,
    selection: SelectionSpec,
    group_species_rows: Sequence[Mapping[str, Any]],
    maximum_groups: int = MAX_SELECTION_GROUPS,
) -> tuple[GroupTaxonEvaluation, ...]:
    """Aggregate group/species rows and evaluate a bounded group collection."""

    if not isinstance(maximum_groups, int) or isinstance(maximum_groups, bool):
        raise InputValidationError("maximum_groups must be an integer.")
    if not 1 <= maximum_groups <= MAX_SELECTION_GROUPS:
        raise InputValidationError(
            f"maximum_groups must be between 1 and {MAX_SELECTION_GROUPS:,}."
        )
    grouped: dict[tuple[str, str, str, str], dict[str, int]] = defaultdict(dict)
    for row in group_species_rows:
        key = (
            str(row.get("run_id", "")),
            str(row.get("group_type", "")),
            str(row.get("hierarchy_node", "")),
            str(row.get("group_id", "")),
        )
        species = str(row.get("species_label", ""))
        if species in grouped[key]:
            raise InputValidationError(
                f"Group/species input contains duplicate row for {key[3]} and {species}."
            )
        grouped[key][species] = _positive_count(value=row.get("species_member_count", 0))
    if len(grouped) > maximum_groups:
        raise InputValidationError(
            f"Selection contains {len(grouped):,} groups; limit is {maximum_groups:,}."
        )
    return tuple(
        evaluate_group_taxonomy(
            graph=graph,
            selection=selection,
            run_id=key[0],
            group_type=key[1],
            hierarchy_node=key[2],
            group_id=key[3],
            species_counts=grouped[key],
        )
        for key in sorted(grouped)
    )


def _merge_node(
    *,
    nodes: dict[str, dict[str, Any]],
    taxon_id: str,
    name: str,
    rank: str,
    parent_taxon_id: str | None,
    depth: int,
    authority: str,
    release: str,
) -> None:
    """Add one node or reject inconsistent metadata for an observed identifier."""

    if not taxon_id or not name or not rank or depth < 0:
        raise InputValidationError("Taxonomy node identifiers, names and ranks are required.")
    candidate = {
        "name": name,
        "rank": rank,
        "parent_taxon_id": parent_taxon_id,
        "depth": depth,
        "authority": authority,
        "release": release,
    }
    existing = nodes.get(taxon_id)
    if existing is None:
        nodes[taxon_id] = candidate
        return
    compared = ("name", "rank", "parent_taxon_id", "depth")
    differences = tuple(key for key in compared if existing[key] != candidate[key])
    if differences:
        raise InputValidationError(
            f"Taxonomy identifier {taxon_id} has inconsistent " + "; ".join(differences)
        )
    existing["authority"] = "; ".join(
        sorted(set(str(existing["authority"]).split("; ")).union({authority}))
    )
    existing["release"] = "; ".join(
        sorted(set(str(existing["release"]).split("; ")).union({release}))
    )


def _mapping_record(*, label: str, record: TaxonomyRecord | None) -> InputTaxonMapping:
    """Convert one supplied or missing taxonomy decision to its derived record."""

    if record is None:
        return InputTaxonMapping(
            workflow_species_label=label,
            taxon_id="",
            mapping_status="MISSING",
            mapping_method="",
            mapping_reason="No reviewed mapping row was supplied.",
            source_name_original=label,
            authority="",
            release="",
            role="input",
        )
    return InputTaxonMapping(
        workflow_species_label=label,
        taxon_id=record.taxon_key if record.mapping_status == "REVIEWED" else "",
        mapping_status=record.mapping_status,
        mapping_method=record.mapping_method,
        mapping_reason=record.review_note,
        source_name_original=record.source_name_original,
        authority=record.taxonomy_authority,
        release=record.taxonomy_release,
        role=record.role,
        source_species_name=record.source_species_name,
        accepted_species_name=record.accepted_species_name,
        ncbi_taxon_id=record.ncbi_taxon_id,
        parent_taxon_id=record.parent_taxon_id,
        parent_taxon_name=record.parent_taxon_name,
        lineage_taxon_ids=record.lineage_taxon_ids,
        lineage_names=record.lineage_names,
        lineage_ranks=record.lineage_ranks,
        taxon_rank=record.taxon_rank,
        authority_taxon_id=record.authority_taxon_id,
        mapping_source=record.mapping_source,
        source_date=record.source_date,
        source_version=record.source_version,
        reviewed_by=record.reviewed_by,
        reviewed_at_utc=record.reviewed_at_utc,
        review_note=record.review_note,
    )


def _normalise_taxon_id(*, value: object) -> str:
    """Return one safe numeric or authority-qualified taxon identifier."""

    if not isinstance(value, str):
        raise InputValidationError("Taxonomy identifiers must be text.")
    normalised = value.strip()
    if not normalised or len(normalised) > 512 or _SAFE_TEXT.fullmatch(normalised) is None:
        raise InputValidationError("Taxonomy identifier is empty, unsafe or overlong.")
    if normalised.isdigit() and int(normalised) <= 0:
        raise InputValidationError("Numeric taxonomy identifiers must be positive.")
    return normalised


def _normalise_identifiers(*, values: Sequence[str]) -> tuple[str, ...]:
    """Validate selector identifiers while preserving deliberate input order."""

    if isinstance(values, (str, bytes)):
        raise InputValidationError("Taxonomy selectors must be a sequence of identifiers.")
    identifiers = tuple(_normalise_taxon_id(value=value) for value in values)
    if len(identifiers) != len(set(identifiers)):
        raise InputValidationError("A taxonomy selector contains duplicate identifiers.")
    return identifiers


def _validate_selection_contradictions(
    *, graph: TaxonomyGraph, selection: SelectionSpec
) -> None:
    """Reject impossible AND-composed predicates before any group query."""

    required = set(selection.identifiers(selector_type="REQUIRED_EXACT"))
    include = set(selection.identifiers(selector_type="INCLUDE_CLADE"))
    only = set(selection.identifiers(selector_type="ONLY_IN_CLADE"))
    excluded_exact = set(selection.identifiers(selector_type="EXCLUDE_EXACT"))
    excluded_clades = set(selection.identifiers(selector_type="EXCLUDE_CLADE"))
    direct = required.intersection(excluded_exact)
    if direct:
        raise InputValidationError(
            "The same exact taxon is required and excluded: " + "; ".join(sorted(direct))
        )
    excluded_universe = set(excluded_exact)
    for taxon_id in excluded_clades:
        excluded_universe.update(graph.descendants(taxon_id=taxon_id, terminals_only=True))
    blocked_required = required.intersection(excluded_universe)
    if blocked_required:
        raise InputValidationError(
            "A required exact taxon is covered by an exclusion: "
            + "; ".join(sorted(blocked_required))
        )
    represented = set(graph.represented_taxon_ids)
    for taxon_id in sorted(include):
        eligible = set(graph.descendants(taxon_id=taxon_id, terminals_only=True))
        eligible.intersection_update(represented)
        if not eligible:
            raise InputValidationError(
                f"Included clade {taxon_id} has no represented reviewed terminal taxa."
            )
        if eligible.issubset(excluded_universe):
            raise InputValidationError(
                f"Included clade {taxon_id} is fully covered by active exclusions."
            )
    allowed = _only_intersection(graph=graph, only_taxon_ids=only)
    if only and not allowed.intersection(represented):
        raise InputValidationError(
            "Only-in clades have no represented terminal intersection."
        )
    outside_required = required.difference(allowed) if only else set()
    if outside_required:
        raise InputValidationError(
            "Required exact taxa fall outside the only-in intersection: "
            + "; ".join(sorted(outside_required))
        )
    for taxon_id in sorted(include):
        if only:
            contributors = set(graph.descendants(taxon_id=taxon_id, terminals_only=True))
            contributors.intersection_update(allowed)
            contributors.intersection_update(represented)
            contributors.difference_update(excluded_universe)
            if not contributors:
                raise InputValidationError(
                    f"Included clade {taxon_id} cannot contribute within the only-in scope."
                )
    if only and allowed.issubset(excluded_universe):
        raise InputValidationError("Active exclusions cover the complete only-in scope.")


def _only_intersection(*, graph: TaxonomyGraph, only_taxon_ids: Iterable[str]) -> set[str]:
    """Return terminal intersection for zero or more only-in clades."""

    identifiers = tuple(sorted(only_taxon_ids))
    if not identifiers:
        return set(graph.terminal_ids)
    descendant_sets = [
        set(graph.descendants(taxon_id=value, terminals_only=True))
        for value in identifiers
    ]
    allowed = descendant_sets[0]
    for descendants in descendant_sets[1:]:
        allowed.intersection_update(descendants)
    return allowed


def _predicate_matches_group(
    *,
    graph: TaxonomyGraph,
    predicate: SelectionPredicate,
    present_taxa: frozenset[str],
    has_unmapped: bool,
    selection: SelectionSpec,
) -> bool:
    """Evaluate one selector without changing its documented meaning."""

    if predicate.selector_type == "REQUIRED_EXACT":
        return predicate.taxon_id in present_taxa
    descendants = graph.descendants(taxon_id=predicate.taxon_id, terminals_only=True)
    if predicate.selector_type == "INCLUDE_CLADE":
        return bool(present_taxa.intersection(descendants))
    if predicate.selector_type == "ONLY_IN_CLADE":
        allowed = _only_intersection(
            graph=graph,
            only_taxon_ids=selection.identifiers(selector_type="ONLY_IN_CLADE"),
        )
        return bool(present_taxa) and present_taxa.issubset(allowed) and not has_unmapped
    if predicate.selector_type == "EXCLUDE_EXACT":
        return predicate.taxon_id not in present_taxa
    if predicate.selector_type == "EXCLUDE_CLADE":
        return not present_taxa.intersection(descendants)
    raise InputValidationError(f"Unknown selector type: {predicate.selector_type}")


def _predicate_failure(*, predicate: SelectionPredicate, has_unmapped: bool) -> str:
    """Return one deterministic, user-facing predicate failure reason."""

    descriptions = {
        "REQUIRED_EXACT": "required exact taxon is absent from the group",
        "INCLUDE_CLADE": "included clade contributes no mapped member",
        "ONLY_IN_CLADE": (
            "only-in scope cannot be proven because a member is outside or unmapped"
            if has_unmapped
            else "a mapped member lies outside the only-in scope"
        ),
        "EXCLUDE_EXACT": "excluded exact taxon occurs in the group",
        "EXCLUDE_CLADE": "an excluded clade descendant occurs in the group",
    }
    return (
        f"{predicate.selector_type}:{predicate.taxon_id}: "
        f"{descriptions[predicate.selector_type]}"
    )


def _outside_scope_taxa(
    *, graph: TaxonomyGraph, selection: SelectionSpec, present_taxa: frozenset[str]
) -> frozenset[str]:
    """Return mapped group taxa outside the active include/only display scope."""

    include_ids = selection.identifiers(selector_type="INCLUDE_CLADE")
    only_ids = selection.identifiers(selector_type="ONLY_IN_CLADE")
    if include_ids:
        scope: set[str] = set()
        for taxon_id in include_ids:
            scope.update(graph.descendants(taxon_id=taxon_id, terminals_only=True))
    elif only_ids:
        scope = _only_intersection(graph=graph, only_taxon_ids=only_ids)
    else:
        return frozenset()
    return frozenset(present_taxa.difference(scope))


def _positive_count(*, value: object, allow_zero: bool = False) -> int:
    """Parse one integer group count with explicit lower bound."""

    if not isinstance(value, int) or isinstance(value, bool):
        raise InputValidationError("Group species-member counts must be integers.")
    minimum = 0 if allow_zero else 1
    if value < minimum:
        raise InputValidationError(
            f"Group species-member counts must be at least {minimum}."
        )
    return value


def _parse_boolean(*, value: str, line_number: int) -> bool:
    """Parse one explicit expected-universe Boolean value."""

    normalised = value.casefold()
    if normalised in {"true", "yes", "1", "included"}:
        return True
    if normalised in {"false", "no", "0", "excluded"}:
        return False
    raise InputValidationError(
        f"Expected-taxon TSV line {line_number} has invalid included value {value!r}."
    )

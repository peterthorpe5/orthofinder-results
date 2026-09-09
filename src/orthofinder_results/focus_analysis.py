"""Deterministic focus-protein selection and machine-readable cluster exports."""

from __future__ import annotations

import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from orthofinder_interrogation_app.focus import FocusProtein, FocusProteinAuthority

from .distances import DISTANCE_STATISTIC_FIELDS
from .errors import InputValidationError
from .io_utils import read_tsv, write_tsv
from .statistics import GROUP_STATISTIC_FIELDS

FOCUS_MATCH_FIELDS = (
    "run_id",
    "group_type",
    "hierarchy_node",
    "group_id",
    "legacy_orthogroup_id",
    "gene_tree_parent_clade",
    "seed_id",
    "matched_member_id",
    "matched_species_label",
    "matched_internal_ids",
    "match_authority",
    "seed_protein_names",
    "seed_categories",
    "seed_review_statuses",
    "seed_ubiquitin_go_statuses",
    "seed_organisms",
    "seed_evidence_type",
    "seed_source",
    "seed_sequence_identifiers",
    "seed_sequence_species",
    "seed_protein_sequence_length",
    "annotation_scope",
    "focus_authority_name",
    "focus_authority_sha256",
)

FOCUS_AUDIT_FIELDS = (
    "seed_id",
    "enabled",
    "match_status",
    "matched_group_count",
    "matched_group_ids",
    "matched_member_count",
    "matched_member_ids",
    "matched_species_count",
    "matched_species_labels",
    "match_authorities",
    "seed_protein_names",
    "seed_categories",
    "seed_review_statuses",
    "seed_ubiquitin_go_statuses",
    "seed_organisms",
    "seed_evidence_type",
    "seed_source",
    "seed_sequence_identifiers",
    "seed_sequence_species",
    "seed_protein_sequence_length",
    "annotation_scope",
    "focus_authority_name",
    "focus_authority_sha256",
)

_GROUP_RESULT_FIELDS = tuple(
    field for field in GROUP_STATISTIC_FIELDS if field != "source_file"
)
_DISTANCE_RESULT_FIELDS = tuple(
    field
    for field in DISTANCE_STATISTIC_FIELDS
    if field
    not in {"run_id", "group_type", "hierarchy_node", "group_id", "source_file"}
)
FOCUS_CLUSTER_RESULT_FIELDS = (
    *_GROUP_RESULT_FIELDS,
    "group_source_file",
    "matched_seed_count",
    "matched_seed_ids",
    "matched_e3_member_count",
    "matched_e3_member_ids",
    "matched_e3_species_count",
    "matched_e3_species_labels",
    "match_authorities",
    "seed_protein_names",
    "seed_categories",
    "seed_review_statuses",
    "seed_ubiquitin_go_statuses",
    "seed_organisms",
    "focus_authority_name",
    "focus_authority_sha256",
    *_DISTANCE_RESULT_FIELDS,
    "distance_source_file",
    "distance_sampling_fraction",
    "expected_sample_pair_count",
    "resolved_pair_fraction",
    "full_group_distance_matrix",
    "distance_interquartile_range",
    "distance_range",
    "distance_coefficient_of_variation",
)

_UNIPROT_PATTERN = re.compile(r"^(?:sp|tr)\|([^|\s]+)\|([^|\s]+)$")


@dataclass(frozen=True)
class FocusSelection:
    """Complete deterministic match state for one group collection."""

    authority: FocusProteinAuthority
    group_type: str
    hierarchy_node: str
    match_rows: tuple[dict[str, Any], ...]
    audit_rows: tuple[dict[str, Any], ...]
    group_statistics: tuple[dict[str, str], ...]
    required_members_by_group: Mapping[str, tuple[str, ...]]

    @property
    def group_keys(self) -> frozenset[str]:
        """Return all selected collision-safe group keys."""

        return frozenset(_group_key(row=row) for row in self.group_statistics)

    @property
    def matched_seed_count(self) -> int:
        """Return the number of enabled seed identifiers matched at least once."""

        return sum(row["match_status"] == "MATCHED" for row in self.audit_rows)

    @property
    def selected_tree_ids(self) -> frozenset[str]:
        """Return possible tree identifiers for all selected groups."""

        identifiers: set[str] = set()
        for row in self.group_statistics:
            identifiers.update(
                value
                for value in (
                    row.get("legacy_orthogroup_id", ""),
                    row.get("group_id", ""),
                )
                if value
            )
        return frozenset(identifiers)


def select_focus_groups(
    *,
    tables_dir: Path,
    run_id: str,
    authority: FocusProteinAuthority,
    group_type: str,
    hierarchy_node: str,
) -> FocusSelection:
    """Match an exact protein authority to one published group collection.

    Args:
        tables_dir: Directory containing compressed TSV analytical relations.
        run_id: Immutable resource run identifier.
        authority: Validated focus-protein authority.
        group_type: ``HOG`` or ``LEGACY_ORTHOGROUP``.
        hierarchy_node: Exact HOG node, or an empty string for the flat authority.

    Returns:
        Complete match, audit, group and required-member state.

    Raises:
        InputValidationError: If controls or table relationships are inconsistent.
    """

    if group_type not in {"HOG", "LEGACY_ORTHOGROUP"}:
        raise InputValidationError(f"Unsupported focus group type: {group_type}")
    if group_type == "LEGACY_ORTHOGROUP" and hierarchy_node:
        raise InputValidationError(
            "Legacy orthogroup focus selection requires the ROOT hierarchy."
        )
    identifiers = frozenset(authority.identifiers)
    sequence_matches, internal_ids = _matched_sequence_aliases(
        path=tables_dir / "sequences.tsv.gz",
        run_id=run_id,
        focus_identifiers=identifiers,
    )
    relation = (
        "hog_memberships.tsv.gz"
        if group_type == "HOG"
        else "legacy_orthogroup_memberships.tsv.gz"
    )
    raw_matches: set[tuple[str, ...]] = set()
    membership_metadata: dict[str, dict[str, str]] = {}
    required: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(path=tables_dir / relation):
        if (
            row["run_id"] != run_id
            or row["group_type"] != group_type
            or row["hierarchy_node"] != hierarchy_node
        ):
            continue
        member_species = (row["member_id"], row["species_label"])
        matches = sequence_matches.get(member_species)
        if matches is None:
            matches = _matching_aliases(
                member_id=row["member_id"],
                internal_id="",
                focus_identifiers=identifiers,
            )
        if not matches:
            continue
        key = _group_key(row=row)
        membership_metadata[key] = row
        required[key].add(row["member_id"])
        for seed_id, match_authority in matches:
            raw_matches.add(
                (
                    key,
                    seed_id,
                    row["member_id"],
                    row["species_label"],
                    match_authority,
                )
            )
    metadata = {record.identifier: record for record in authority.records}
    match_rows = tuple(
        _match_record(
            group=membership_metadata[key],
            seed=metadata[seed_id],
            member_id=member_id,
            species_label=species_label,
            internal_ids=internal_ids.get((member_id, species_label), ()),
            match_authority=match_authority,
            authority=authority,
        )
        for key, seed_id, member_id, species_label, match_authority in sorted(raw_matches)
    )
    group_statistics = _selected_group_statistics(
        path=tables_dir / "group_statistics.tsv.gz",
        selected_keys=frozenset(membership_metadata),
    )
    audit_rows = _focus_audit_rows(authority=authority, matches=match_rows)
    return FocusSelection(
        authority=authority,
        group_type=group_type,
        hierarchy_node=hierarchy_node,
        match_rows=match_rows,
        audit_rows=audit_rows,
        group_statistics=group_statistics,
        required_members_by_group={
            key: tuple(sorted(members)) for key, members in sorted(required.items())
        },
    )


def publish_focus_selection(*, tables_dir: Path, selection: FocusSelection) -> None:
    """Publish complete focus match and seed-audit relations.

    Args:
        tables_dir: Analytical TSV output directory.
        selection: Completed deterministic focus selection.
    """

    write_tsv(
        path=tables_dir / "e3_seed_matches.tsv.gz",
        fieldnames=FOCUS_MATCH_FIELDS,
        records=selection.match_rows,
    )
    write_tsv(
        path=tables_dir / "e3_seed_catalogue_audit.tsv.gz",
        fieldnames=FOCUS_AUDIT_FIELDS,
        records=selection.audit_rows,
    )


def publish_focus_cluster_results(
    *,
    tables_dir: Path,
    selection: FocusSelection,
    distance_summaries: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, Any], ...]:
    """Publish one analysis-ready row for every selected focus cluster.

    Args:
        tables_dir: Analytical TSV output directory.
        selection: Complete focus selection.
        distance_summaries: Distance result for each selected group.

    Returns:
        Published rows in deterministic group order.

    Raises:
        InputValidationError: If selected groups lack one unambiguous summary.
    """

    summaries: dict[str, Mapping[str, Any]] = {}
    for summary in distance_summaries:
        key = _group_key(row=summary)
        if key in summaries:
            raise InputValidationError(
                f"Focus group has multiple distance summaries: {key}"
            )
        summaries[key] = summary
    selected_keys = selection.group_keys
    missing = sorted(selected_keys.difference(summaries))
    if missing:
        raise InputValidationError(
            "Focus groups lack explicit distance summaries: " + ";".join(missing[:20])
        )
    matches_by_group: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in selection.match_rows:
        matches_by_group[_group_key(row=row)].append(row)
    rows = tuple(
        _cluster_result_record(
            statistic=statistic,
            matches=matches_by_group[_group_key(row=statistic)],
            distance=summaries[_group_key(row=statistic)],
            authority=selection.authority,
        )
        for statistic in selection.group_statistics
    )
    write_tsv(
        path=tables_dir / "e3_cluster_results.tsv.gz",
        fieldnames=FOCUS_CLUSTER_RESULT_FIELDS,
        records=rows,
    )
    return rows


def _matched_sequence_aliases(
    *,
    path: Path,
    run_id: str,
    focus_identifiers: frozenset[str],
) -> tuple[
    dict[tuple[str, str], tuple[tuple[str, str], ...]],
    dict[tuple[str, str], tuple[str, ...]],
]:
    """Return exact matching aliases and internal IDs by member and species.

    Args:
        path: Published sequence-identifier relation.
        run_id: Immutable resource run identifier.
        focus_identifiers: Exact enabled focus identifiers.

    Returns:
        Match tuples and all associated internal IDs by canonical member/species.
    """

    matches: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)
    internal_ids: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in read_tsv(path=path):
        if row["run_id"] != run_id:
            continue
        member_species = (row["member_id"], row["species_label"])
        row_matches = _matching_aliases(
            member_id=row["member_id"],
            internal_id=row["internal_id"],
            focus_identifiers=focus_identifiers,
        )
        if row_matches:
            matches[member_species].update(row_matches)
            if row["internal_id"]:
                internal_ids[member_species].add(row["internal_id"])
    return (
        {key: tuple(sorted(values)) for key, values in matches.items()},
        {key: tuple(sorted(values)) for key, values in internal_ids.items()},
    )


def _matching_aliases(
    *, member_id: str, internal_id: str, focus_identifiers: frozenset[str]
) -> tuple[tuple[str, str], ...]:
    """Return focus identifiers matching controlled exact member aliases.

    Args:
        member_id: Canonical OrthoFinder membership identifier.
        internal_id: Optional OrthoFinder internal identifier.
        focus_identifiers: Exact identifiers from the focus authority.

    Returns:
        Unique ``(focus identifier, match authority)`` pairs.
    """

    aliases = [(member_id, "PROTEIN_ID")]
    if internal_id:
        aliases.append((internal_id, "ORTHOFINDER_INTERNAL_ID"))
    uniprot = _UNIPROT_PATTERN.fullmatch(member_id)
    if uniprot is not None:
        aliases.extend(
            (
                (uniprot.group(1), "UNIPROT_ACCESSION"),
                (uniprot.group(2), "UNIPROT_ENTRY"),
            )
        )
    return tuple(
        sorted(
            {
                (alias, authority)
                for alias, authority in aliases
                if alias in focus_identifiers
            }
        )
    )


def _selected_group_statistics(
    *, path: Path, selected_keys: frozenset[str]
) -> tuple[dict[str, str], ...]:
    """Load exactly one statistic row for every selected focus group.

    Args:
        path: Complete group-statistics relation.
        selected_keys: Collision-safe selected group keys.

    Returns:
        Selected statistic records sorted by group identity.

    Raises:
        InputValidationError: If a selected key is missing or duplicated.
    """

    selected: dict[str, dict[str, str]] = {}
    for row in read_tsv(path=path):
        key = _group_key(row=row)
        if key not in selected_keys:
            continue
        if key in selected:
            raise InputValidationError(f"Duplicate group-statistics row: {key}")
        selected[key] = row
    missing = sorted(selected_keys.difference(selected))
    if missing:
        raise InputValidationError(
            "Focus memberships lack group statistics: " + ";".join(missing[:20])
        )
    return tuple(selected[key] for key in sorted(selected))


def _focus_audit_rows(
    *, authority: FocusProteinAuthority, matches: Sequence[Mapping[str, Any]]
) -> tuple[dict[str, Any], ...]:
    """Return one coverage-audit row for every configured seed.

    Args:
        authority: Complete focus authority including disabled records.
        matches: Exact member/group matches.

    Returns:
        Audit rows in deterministic seed order.
    """

    by_seed: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in matches:
        by_seed[str(row["seed_id"])].append(row)
    rows: list[dict[str, Any]] = []
    for seed in authority.records:
        seed_matches = by_seed.get(seed.identifier, [])
        groups = sorted({str(row["group_id"]) for row in seed_matches})
        members = sorted({str(row["matched_member_id"]) for row in seed_matches})
        species = sorted({str(row["matched_species_label"]) for row in seed_matches})
        authorities = sorted({str(row["match_authority"]) for row in seed_matches})
        rows.append(
            {
                "seed_id": seed.identifier,
                "enabled": seed.enabled,
                "match_status": (
                    "DISABLED" if not seed.enabled else "MATCHED" if seed_matches else "UNMATCHED"
                ),
                "matched_group_count": len(groups),
                "matched_group_ids": ";".join(groups),
                "matched_member_count": len(members),
                "matched_member_ids": ";".join(members),
                "matched_species_count": len(species),
                "matched_species_labels": ";".join(species),
                "match_authorities": ";".join(authorities),
                **_seed_metadata(seed=seed),
                "focus_authority_name": authority.source_name,
                "focus_authority_sha256": authority.sha256,
            }
        )
    return tuple(rows)


def _match_record(
    *,
    group: Mapping[str, str],
    seed: FocusProtein,
    member_id: str,
    species_label: str,
    internal_ids: Sequence[str],
    match_authority: str,
    authority: FocusProteinAuthority,
) -> dict[str, Any]:
    """Return one provenance-rich seed-to-member-to-group match row.

    Args:
        group: Matched membership record.
        seed: Matched seed metadata.
        member_id: Canonical matched member identifier.
        species_label: Exact matched species label.
        internal_ids: Associated OrthoFinder internal identifiers.
        match_authority: Exact alias authority used for the match.
        authority: Complete input authority provenance.

    Returns:
        Schema-complete match record.
    """

    return {
        "run_id": group["run_id"],
        "group_type": group["group_type"],
        "hierarchy_node": group["hierarchy_node"],
        "group_id": group["group_id"],
        "legacy_orthogroup_id": group.get("legacy_orthogroup_id", ""),
        "gene_tree_parent_clade": group.get("gene_tree_parent_clade", ""),
        "seed_id": seed.identifier,
        "matched_member_id": member_id,
        "matched_species_label": species_label,
        "matched_internal_ids": ";".join(sorted(internal_ids)),
        "match_authority": match_authority,
        **_seed_metadata(seed=seed),
        "focus_authority_name": authority.source_name,
        "focus_authority_sha256": authority.sha256,
    }


def _seed_metadata(*, seed: FocusProtein) -> dict[str, str]:
    """Return stable analytical metadata for one focus seed.

    Args:
        seed: Validated focus record.

    Returns:
        Metadata fields shared by match and audit relations.
    """

    return {
        "seed_protein_names": seed.protein_name,
        "seed_categories": seed.category,
        "seed_review_statuses": seed.review_status,
        "seed_ubiquitin_go_statuses": seed.ubiquitin_go_status,
        "seed_organisms": seed.organism,
        "seed_evidence_type": seed.evidence_type,
        "seed_source": seed.source,
        "seed_sequence_identifiers": seed.sequence_identifiers,
        "seed_sequence_species": seed.sequence_species,
        "seed_protein_sequence_length": seed.protein_sequence_length,
        "annotation_scope": seed.annotation_scope,
    }


def _cluster_result_record(
    *,
    statistic: Mapping[str, Any],
    matches: Sequence[Mapping[str, Any]],
    distance: Mapping[str, Any],
    authority: FocusProteinAuthority,
) -> dict[str, Any]:
    """Join one group, its E3 matches and its distance distribution.

    Args:
        statistic: Complete group-statistics record.
        matches: All exact E3 seed matches for the group.
        distance: One explicit distance summary.
        authority: Focus authority provenance.

    Returns:
        Flat analysis-ready E3 cluster record.
    """

    seeds = sorted({str(row["seed_id"]) for row in matches})
    members = sorted({str(row["matched_member_id"]) for row in matches})
    species = sorted({str(row["matched_species_label"]) for row in matches})
    expected_pairs = _expected_pairs(value=distance.get("sampled_member_count"))
    resolved_pairs = _integer(value=distance.get("distance_pair_count"))
    unresolved_pairs = _integer(value=distance.get("unresolved_pair_count"))
    total_members = _integer(value=distance.get("total_member_count"))
    sampled_members = _integer(value=distance.get("sampled_member_count"))
    has_distances = resolved_pairs > 0
    q25 = _number(value=distance.get("q25_distance")) if has_distances else None
    q75 = _number(value=distance.get("q75_distance")) if has_distances else None
    minimum = _number(value=distance.get("minimum_distance")) if has_distances else None
    maximum = _number(value=distance.get("maximum_distance")) if has_distances else None
    mean = _number(value=distance.get("mean_distance")) if has_distances else None
    deviation = (
        _number(value=distance.get("population_stddev_distance"))
        if has_distances
        else None
    )
    distance_fields = {
        field: (
            ""
            if not has_distances and field in _distance_numeric_fields()
            else distance.get(field, "")
        )
        for field in _DISTANCE_RESULT_FIELDS
    }
    return {
        **{field: statistic.get(field, "") for field in _GROUP_RESULT_FIELDS},
        "group_source_file": statistic.get("source_file", ""),
        "matched_seed_count": len(seeds),
        "matched_seed_ids": ";".join(seeds),
        "matched_e3_member_count": len(members),
        "matched_e3_member_ids": ";".join(members),
        "matched_e3_species_count": len(species),
        "matched_e3_species_labels": ";".join(species),
        "match_authorities": _joined_values(matches=matches, field="match_authority"),
        "seed_protein_names": _joined_values(matches=matches, field="seed_protein_names"),
        "seed_categories": _joined_values(matches=matches, field="seed_categories"),
        "seed_review_statuses": _joined_values(
            matches=matches, field="seed_review_statuses"
        ),
        "seed_ubiquitin_go_statuses": _joined_values(
            matches=matches, field="seed_ubiquitin_go_statuses"
        ),
        "seed_organisms": _joined_values(matches=matches, field="seed_organisms"),
        "focus_authority_name": authority.source_name,
        "focus_authority_sha256": authority.sha256,
        **distance_fields,
        "distance_source_file": distance.get("source_file", ""),
        "distance_sampling_fraction": (
            sampled_members / total_members if total_members else ""
        ),
        "expected_sample_pair_count": expected_pairs,
        "resolved_pair_fraction": (
            resolved_pairs / expected_pairs if expected_pairs else ""
        ),
        "full_group_distance_matrix": bool(
            total_members >= 2
            and sampled_members == total_members
            and resolved_pairs + unresolved_pairs == expected_pairs
        ),
        "distance_interquartile_range": (
            q75 - q25 if q25 is not None and q75 is not None else ""
        ),
        "distance_range": (
            maximum - minimum
            if minimum is not None and maximum is not None
            else ""
        ),
        "distance_coefficient_of_variation": (
            deviation / mean
            if deviation is not None and mean is not None and mean > 0
            else ""
        ),
    }


def _joined_values(*, matches: Sequence[Mapping[str, Any]], field: str) -> str:
    """Return sorted distinct non-empty values from match records.

    Args:
        matches: Match records for one group.
        field: Exact match-record field.

    Returns:
        Semicolon-delimited deterministic values.
    """

    return ";".join(sorted({str(row.get(field, "")) for row in matches if row.get(field)}))


def _distance_numeric_fields() -> frozenset[str]:
    """Return summary fields whose unavailable values must remain blank."""

    return frozenset(
        {
            "minimum_distance",
            "q05_distance",
            "q25_distance",
            "median_distance",
            "mean_distance",
            "q75_distance",
            "q95_distance",
            "maximum_distance",
            "population_stddev_distance",
            "mean_comparable_sites",
        }
    )


def _expected_pairs(*, value: Any) -> int:
    """Return the complete unordered-pair count for a sampled member count.

    Args:
        value: Integer-like sampled member count.

    Returns:
        Non-negative expected pair count.
    """

    count = _integer(value=value)
    return count * (count - 1) // 2 if count >= 2 else 0


def _integer(*, value: Any) -> int:
    """Return an integer-like value, treating an empty value as zero.

    Args:
        value: Integer-like analytical value.

    Returns:
        Parsed integer or zero for missing values.
    """

    return 0 if value in {None, ""} else int(value)


def _number(*, value: Any) -> float | None:
    """Return a finite float-like analytical value or ``None`` when absent.

    Args:
        value: Numeric or empty analytical value.

    Returns:
        Parsed float or ``None``.
    """

    return None if value in {None, ""} else float(value)


def _group_key(*, row: Mapping[str, Any]) -> str:
    """Return a collision-safe run-scoped group key.

    Args:
        row: Group-like analytical record.

    Returns:
        Group type, hierarchy and identifier joined by a reserved separator.
    """

    return "|".join(
        (
            str(row.get("group_type", "")),
            str(row.get("hierarchy_node", "")),
            str(row.get("group_id", "")),
        )
    )

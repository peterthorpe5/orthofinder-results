"""Matched-background dispersion benchmarking at the cluster level."""

from __future__ import annotations

import hashlib
import logging
import math
import random
import statistics
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .benchmark_authority import BenchmarkAuthority
from .errors import InputValidationError
from .focus_analysis import FocusSelection
from .io_utils import read_tsv, write_tsv

_LOGGER = logging.getLogger("orthofinder_results.benchmark_analysis")
MINIMUM_BACKGROUND_GROUPS = 3
DISPERSION_METRICS = (
    "mean_distance",
    "median_distance",
    "population_stddev_distance",
    "distance_interquartile_range",
    "distance_coefficient_of_variation",
)
MATCHED_CONTROL_FIELDS = (
    "run_id",
    "anchor_group_type",
    "anchor_hierarchy_node",
    "anchor_group_id",
    "control_group_type",
    "control_hierarchy_node",
    "control_group_id",
    "control_rank",
    "matching_score",
    "anchor_member_count",
    "control_member_count",
    "anchor_species_count",
    "control_species_count",
    "anchor_mean_copies_per_species",
    "control_mean_copies_per_species",
    "anchor_single_copy_fraction",
    "control_single_copy_fraction",
    "control_reused",
    "matching_variables",
)
PROFILE_FIELDS = (
    "run_id",
    "group_type",
    "hierarchy_node",
    "group_id",
    "profile_id",
    "profile_class",
    "profile_subclass",
    "membership_role",
    "marker_count",
    "marker_ids",
    "protein_identifiers",
    "matched_member_count",
    "matched_member_ids",
    "matched_species_labels",
    "authority_names",
    "authority_sha256s",
)
MARKER_MATCH_FIELDS = (
    "run_id",
    "group_type",
    "hierarchy_node",
    "group_id",
    "legacy_orthogroup_id",
    "marker_id",
    "protein_identifier",
    "protein_entry",
    "marker_name",
    "benchmark_class",
    "benchmark_subclass",
    "domain_architecture",
    "matched_member_id",
    "matched_species_label",
    "match_authority",
    "evidence_type",
    "source_title",
    "source_doi",
    "source_table",
    "source_version",
    "benchmark_authority_name",
    "benchmark_authority_sha256",
)
MARKER_AUDIT_FIELDS = (
    "marker_id",
    "protein_identifier",
    "protein_entry",
    "marker_name",
    "benchmark_class",
    "benchmark_subclass",
    "enabled",
    "match_status",
    "matched_group_count",
    "matched_group_ids",
    "matched_member_count",
    "matched_member_ids",
    "evidence_type",
    "source_title",
    "source_doi",
    "source_table",
    "source_version",
    "benchmark_authority_name",
    "benchmark_authority_sha256",
)
CLUSTER_RESULT_FIELDS = (
    "run_id",
    "group_type",
    "hierarchy_node",
    "group_id",
    "legacy_orthogroup_id",
    "gene_tree_parent_clade",
    "member_count",
    "species_count",
    "single_copy_species_count",
    "max_copies_per_species",
    "mean_copies_per_species",
    "profile_ids",
    "profile_classes",
    "profile_subclasses",
    "membership_roles",
    "distance_method",
    "computation_status",
    "member_identifier_resolution",
    "total_member_count",
    "sampled_member_count",
    "distance_pair_count",
    "unresolved_pair_count",
    "minimum_distance",
    "q05_distance",
    "q25_distance",
    "median_distance",
    "mean_distance",
    "q75_distance",
    "q95_distance",
    "maximum_distance",
    "population_stddev_distance",
    "distance_interquartile_range",
    "distance_coefficient_of_variation",
    "sampling_fraction",
    "pairwise_rows_persisted",
    "failure_reason",
)
BACKGROUND_STATISTIC_FIELDS = (
    "run_id",
    "profile_id",
    "profile_class",
    "profile_subclass",
    "comparison_scale",
    "metric",
    "eligible_group_count",
    "minimum_value",
    "q25_value",
    "median_value",
    "mean_value",
    "q75_value",
    "maximum_value",
    "population_stddev_value",
)
CONTRAST_FIELDS = (
    "run_id",
    "contrast_id",
    "target_profile_id",
    "reference_profile_id",
    "comparison_scale",
    "metric",
    "target_group_count",
    "reference_group_count",
    "target_median",
    "reference_median",
    "median_difference",
    "median_difference_ci_low",
    "median_difference_ci_high",
    "cliffs_delta",
    "mann_whitney_u",
    "p_value_two_sided",
    "fdr_family",
    "fdr_q_value",
    "status",
    "interpretation",
)
INDIVIDUAL_COMPARISON_FIELDS = (
    "run_id",
    "group_type",
    "hierarchy_node",
    "group_id",
    "background_profile_id",
    "comparison_scale",
    "metric",
    "observed_value",
    "background_group_count",
    "background_median",
    "difference_from_background_median",
    "empirical_percentile",
    "lower_tail_p_value",
    "upper_tail_p_value",
    "two_sided_p_value",
    "fdr_family",
    "fdr_q_value",
    "leave_one_out",
    "status",
)
CLASSIFICATION_FIELDS = (
    "run_id",
    "group_type",
    "hierarchy_node",
    "group_id",
    "matched_control_count",
    "mean_distance",
    "mean_distance_control_median",
    "mean_distance_control_percentile",
    "central_divergence_class",
    "distance_sd",
    "distance_sd_control_median",
    "distance_sd_control_percentile",
    "heterogeneity_class",
    "classification_rule",
    "status",
)


@dataclass(frozen=True)
class BenchmarkPlan:
    """Target panels and deterministic non-focus matched controls."""

    run_id: str
    e3_selection: FocusSelection
    marker_selection: FocusSelection
    marker_authority: BenchmarkAuthority
    target_statistics: tuple[dict[str, str], ...]
    control_statistics: tuple[dict[str, str], ...]
    control_links: tuple[dict[str, Any], ...]

    @property
    def target_group_keys(self) -> frozenset[str]:
        """Return all unique biological-panel group keys."""

        return frozenset(_group_key(row=row) for row in self.target_statistics)

    @property
    def control_group_keys(self) -> frozenset[str]:
        """Return all selected non-focus control group keys."""

        return frozenset(_group_key(row=row) for row in self.control_statistics)

    @property
    def selected_group_keys(self) -> frozenset[str]:
        """Return the complete distance-analysis group set."""

        return self.target_group_keys | self.control_group_keys

    @property
    def required_members_by_group(self) -> dict[str, tuple[str, ...]]:
        """Return target proteins that bounded samples must retain."""

        required: dict[str, set[str]] = defaultdict(set)
        for selection in (self.e3_selection, self.marker_selection):
            for key, members in selection.required_members_by_group.items():
                required[key].update(members)
        return {key: tuple(sorted(values)) for key, values in required.items()}

    @property
    def selected_tree_ids(self) -> frozenset[str]:
        """Return possible tree identifiers for targets and controls."""

        identifiers = set()
        for row in (*self.target_statistics, *self.control_statistics):
            identifiers.update(
                value
                for value in (
                    row.get("legacy_orthogroup_id", ""),
                    row.get("group_id", ""),
                )
                if value
            )
        return frozenset(identifiers)


def build_benchmark_plan(
    *,
    tables_dir: Path,
    run_id: str,
    e3_selection: FocusSelection,
    marker_selection: FocusSelection,
    marker_authority: BenchmarkAuthority,
    controls_per_group: int,
) -> BenchmarkPlan:
    """Construct target union and covariate-matched non-focus controls.

    Matching never reads a distance value. It uses group size, represented
    species, mean copy count and single-copy fraction at the same HOG level.

    Args:
        tables_dir: Published analytical TSV staging directory.
        run_id: Immutable run identifier.
        e3_selection: Exact E3 seed-to-group selection.
        marker_selection: Exact housekeeping and R/NLR marker selection.
        marker_authority: Reviewed marker metadata.
        controls_per_group: Unique background controls requested per target.

    Returns:
        Complete deterministic benchmark plan.

    Raises:
        InputValidationError: If controls are invalid or cannot be selected.
    """

    if not 1 <= controls_per_group <= 50:
        raise InputValidationError("controls_per_group must be between 1 and 50.")
    target_by_key = {
        _group_key(row=row): row
        for row in (*e3_selection.group_statistics, *marker_selection.group_statistics)
    }
    target_statistics = tuple(target_by_key[key] for key in sorted(target_by_key))
    candidates = [
        row
        for row in read_tsv(path=tables_dir / "group_statistics.tsv.gz")
        if row["run_id"] == run_id
        and row["group_type"] == "HOG"
        and row["hierarchy_node"] == "N0"
        and int(row["member_count"]) >= 2
        and row.get("legacy_orthogroup_id", "")
        and _group_key(row=row) not in target_by_key
    ]
    controls, links = _match_controls(
        anchors=target_statistics,
        candidates=candidates,
        controls_per_group=controls_per_group,
    )
    if not controls:
        raise InputValidationError("No eligible non-focus matched controls were found.")
    _LOGGER.info(
        "Dispersion benchmark plan prepared: targets=%s, controls=%s, links=%s",
        len(target_statistics),
        len(controls),
        len(links),
    )
    return BenchmarkPlan(
        run_id=run_id,
        e3_selection=e3_selection,
        marker_selection=marker_selection,
        marker_authority=marker_authority,
        target_statistics=target_statistics,
        control_statistics=controls,
        control_links=links,
    )


def publish_benchmark_inputs(*, tables_dir: Path, plan: BenchmarkPlan) -> None:
    """Publish marker matches, audits, profiles and matched-control links."""

    marker_matches, marker_audit = _marker_rows(plan=plan)
    profiles = build_profile_rows(plan=plan)
    write_tsv(
        path=tables_dir / "benchmark_marker_matches.tsv.gz",
        fieldnames=MARKER_MATCH_FIELDS,
        records=marker_matches,
    )
    write_tsv(
        path=tables_dir / "benchmark_marker_audit.tsv.gz",
        fieldnames=MARKER_AUDIT_FIELDS,
        records=marker_audit,
    )
    write_tsv(
        path=tables_dir / "benchmark_group_profiles.tsv.gz",
        fieldnames=PROFILE_FIELDS,
        records=profiles,
    )
    write_tsv(
        path=tables_dir / "benchmark_matched_controls.tsv.gz",
        fieldnames=MATCHED_CONTROL_FIELDS,
        records=plan.control_links,
    )


def publish_benchmark_results(
    *,
    tables_dir: Path,
    plan: BenchmarkPlan,
    distance_summaries: Sequence[Mapping[str, Any]],
    bootstrap_resamples: int,
) -> dict[str, int]:
    """Publish cluster, background, contrast and individual benchmark tables.

    Args:
        tables_dir: Analytical staging directory.
        plan: Completed benchmark selection plan.
        distance_summaries: One explicit summary per attempted group.
        bootstrap_resamples: Deterministic percentile-bootstrap iterations.

    Returns:
        Relation row counts for manifests and quality control.
    """

    if not 100 <= bootstrap_resamples <= 100_000:
        raise InputValidationError(
            "bootstrap_resamples must be between 100 and 100,000."
        )
    profiles = build_profile_rows(plan=plan)
    cluster_rows = _cluster_result_rows(
        plan=plan,
        profiles=profiles,
        distance_summaries=distance_summaries,
    )
    values = _metric_values(cluster_rows=cluster_rows)
    residuals = _matched_residuals(plan=plan, values=values)
    background_rows = _background_statistics(
        run_id=plan.run_id,
        profiles=profiles,
        values=values,
        residuals=residuals,
    )
    contrast_rows = _profile_contrasts(
        run_id=plan.run_id,
        profiles=profiles,
        residuals=residuals,
        bootstrap_resamples=bootstrap_resamples,
    )
    individual_rows = _individual_comparisons(
        run_id=plan.run_id,
        plan=plan,
        profiles=profiles,
        values=values,
        residuals=residuals,
    )
    classifications = _classification_rows(
        run_id=plan.run_id,
        plan=plan,
        values=values,
    )
    outputs = (
        ("benchmark_cluster_results.tsv.gz", CLUSTER_RESULT_FIELDS, cluster_rows),
        (
            "benchmark_background_statistics.tsv.gz",
            BACKGROUND_STATISTIC_FIELDS,
            background_rows,
        ),
        ("benchmark_contrasts.tsv.gz", CONTRAST_FIELDS, contrast_rows),
        (
            "benchmark_individual_comparisons.tsv.gz",
            INDIVIDUAL_COMPARISON_FIELDS,
            individual_rows,
        ),
        (
            "benchmark_cluster_classifications.tsv.gz",
            CLASSIFICATION_FIELDS,
            classifications,
        ),
    )
    for filename, fields, rows in outputs:
        write_tsv(path=tables_dir / filename, fieldnames=fields, records=rows)
    counts = {
        "benchmark_cluster_count": len(cluster_rows),
        "benchmark_background_statistic_count": len(background_rows),
        "benchmark_contrast_count": len(contrast_rows),
        "benchmark_individual_comparison_count": len(individual_rows),
        "benchmark_classification_count": len(classifications),
    }
    _LOGGER.info(
        "Dispersion statistics published: clusters=%s, contrasts=%s, "
        "individual_comparisons=%s",
        len(cluster_rows),
        len(contrast_rows),
        len(individual_rows),
    )
    return counts


def build_profile_rows(*, plan: BenchmarkPlan) -> tuple[dict[str, Any], ...]:
    """Return one aggregated group row per biological background profile."""

    records: dict[tuple[str, str], dict[str, set[str] | str]] = {}

    def add(
        *,
        group: Mapping[str, Any],
        profile_id: str,
        profile_class: str,
        profile_subclass: str,
        role: str,
        marker_id: str = "",
        protein_identifier: str = "",
        member_id: str = "",
        species_label: str = "",
        authority_name: str = "",
        authority_sha256: str = "",
    ) -> None:
        key = (_group_key(row=group), profile_id)
        record = records.setdefault(
            key,
            {
                "run_id": str(group["run_id"]),
                "group_type": str(group["group_type"]),
                "hierarchy_node": str(group["hierarchy_node"]),
                "group_id": str(group["group_id"]),
                "profile_id": profile_id,
                "profile_class": profile_class,
                "profile_subclass": profile_subclass,
                "membership_role": role,
                "marker_ids": set(),
                "protein_identifiers": set(),
                "matched_member_ids": set(),
                "matched_species_labels": set(),
                "authority_names": set(),
                "authority_sha256s": set(),
            },
        )
        for field, value in (
            ("marker_ids", marker_id),
            ("protein_identifiers", protein_identifier),
            ("matched_member_ids", member_id),
            ("matched_species_labels", species_label),
            ("authority_names", authority_name),
            ("authority_sha256s", authority_sha256),
        ):
            if value:
                container = record[field]
                if isinstance(container, set):
                    container.add(value)

    e3_authority = plan.e3_selection.authority
    for match in plan.e3_selection.match_rows:
        add(
            group=match,
            profile_id="E3_ALL",
            profile_class="E3",
            profile_subclass="ALL",
            role="TARGET",
            marker_id=str(match["seed_id"]),
            protein_identifier=str(match["seed_id"]),
            member_id=str(match["matched_member_id"]),
            species_label=str(match["matched_species_label"]),
            authority_name=e3_authority.source_name,
            authority_sha256=e3_authority.sha256,
        )
        for category in _categories(value=str(match.get("seed_categories", ""))):
            add(
                group=match,
                profile_id=f"E3_CATEGORY::{category}",
                profile_class="E3",
                profile_subclass=category,
                role="TARGET",
                marker_id=str(match["seed_id"]),
                protein_identifier=str(match["seed_id"]),
                member_id=str(match["matched_member_id"]),
                species_label=str(match["matched_species_label"]),
                authority_name=e3_authority.source_name,
                authority_sha256=e3_authority.sha256,
            )
    markers = plan.marker_authority.by_identifier()
    marker_authority = plan.marker_authority
    for match in plan.marker_selection.match_rows:
        marker = markers[str(match["seed_id"])]
        profile_all = f"{marker.benchmark_class}_ALL"
        common = {
            "group": match,
            "profile_class": marker.benchmark_class,
            "role": "TARGET",
            "marker_id": marker.marker_id,
            "protein_identifier": marker.protein_identifier,
            "member_id": str(match["matched_member_id"]),
            "species_label": str(match["matched_species_label"]),
            "authority_name": marker_authority.source_name,
            "authority_sha256": marker_authority.sha256,
        }
        add(
            **common,
            profile_id=profile_all,
            profile_subclass="ALL",
        )
        add(
            **common,
            profile_id=(
                f"{marker.benchmark_class}_SUBCLASS::{marker.benchmark_subclass}"
            ),
            profile_subclass=marker.benchmark_subclass,
        )
    for group in plan.control_statistics:
        add(
            group=group,
            profile_id="MATCHED_NON_FOCUS",
            profile_class="NON_FOCUS_CONTROL",
            profile_subclass="MATCHED",
            role="MATCHED_CONTROL",
        )
    rows = []
    for key in sorted(records):
        record = records[key]
        marker_ids = record["marker_ids"]
        members = record["matched_member_ids"]
        assert isinstance(marker_ids, set) and isinstance(members, set)
        rows.append(
            {
                field: (
                    len(marker_ids)
                    if field == "marker_count"
                    else len(members)
                    if field == "matched_member_count"
                    else ";".join(sorted(record[field]))
                    if isinstance(record.get(field), set)
                    else record.get(field, "")
                )
                for field in PROFILE_FIELDS
            }
        )
    return tuple(rows)


def benjamini_hochberg(*, p_values: Sequence[float]) -> tuple[float, ...]:
    """Return monotonic Benjamini–Hochberg adjusted values."""

    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in p_values):
        raise ValueError("p_values must be finite values between zero and one.")
    count = len(p_values)
    if not count:
        return ()
    ordered = sorted(enumerate(p_values), key=lambda item: (item[1], item[0]))
    adjusted = [1.0] * count
    running = 1.0
    for reverse_rank, (index, value) in enumerate(reversed(ordered), start=1):
        rank = count - reverse_rank + 1
        running = min(running, value * count / rank)
        adjusted[index] = min(1.0, running)
    return tuple(adjusted)


def mann_whitney_test(
    *, target: Sequence[float], reference: Sequence[float]
) -> tuple[float, float, float]:
    """Return target U, two-sided asymptotic p and Cliff's delta.

    Tied observations receive average ranks and a tie-corrected variance. The
    implementation uses cluster summaries as units; protein pairs are never
    treated as independent observations.
    """

    left = tuple(float(value) for value in target)
    right = tuple(float(value) for value in reference)
    if not left or not right or any(
        not math.isfinite(value) for value in (*left, *right)
    ):
        raise ValueError("Both samples require finite observations.")
    labelled = sorted(
        [(value, 0) for value in left] + [(value, 1) for value in right]
    )
    rank_sum = 0.0
    tie_sum = 0.0
    index = 0
    while index < len(labelled):
        end = index + 1
        while end < len(labelled) and labelled[end][0] == labelled[index][0]:
            end += 1
        average_rank = (index + 1 + end) / 2
        rank_sum += average_rank * sum(
            sample == 0 for _, sample in labelled[index:end]
        )
        tie_size = end - index
        tie_sum += tie_size**3 - tie_size
        index = end
    n_left = len(left)
    n_right = len(right)
    u_value = rank_sum - n_left * (n_left + 1) / 2
    pair_count = n_left * n_right
    delta = (2 * u_value / pair_count) - 1
    total = n_left + n_right
    variance = n_left * n_right / 12 * (
        total + 1 - tie_sum / (total * (total - 1))
    )
    if variance <= 0:
        p_value = 1.0
    else:
        difference = abs(u_value - pair_count / 2)
        corrected = max(0.0, difference - 0.5)
        z_score = corrected / math.sqrt(variance)
        p_value = math.erfc(z_score / math.sqrt(2))
    return u_value, min(1.0, p_value), delta


def _match_controls(
    *,
    anchors: Sequence[Mapping[str, str]],
    candidates: Sequence[Mapping[str, str]],
    controls_per_group: int,
) -> tuple[tuple[dict[str, str], ...], tuple[dict[str, Any], ...]]:
    """Greedily select unique nearest covariate controls."""

    available = {_group_key(row=row): row for row in candidates}
    selected: dict[str, dict[str, str]] = {}
    links: list[dict[str, Any]] = []
    for anchor in sorted(anchors, key=lambda row: _group_key(row=row)):
        ranked = sorted(
            available.values(),
            key=lambda row: (_matching_score(anchor=anchor, control=row), _group_key(row=row)),
        )[:controls_per_group]
        if len(ranked) < controls_per_group:
            raise InputValidationError(
                f"Only {len(ranked)} unique matched controls remain for "
                f"{_group_key(row=anchor)}; requested {controls_per_group}."
            )
        for rank, control in enumerate(ranked, start=1):
            control_key = _group_key(row=control)
            score = _matching_score(anchor=anchor, control=control)
            selected[control_key] = dict(control)
            del available[control_key]
            links.append(
                _control_link(anchor=anchor, control=control, rank=rank, score=score)
            )
    return (
        tuple(selected[key] for key in sorted(selected)),
        tuple(links),
    )


def _matching_score(
    *, anchor: Mapping[str, str], control: Mapping[str, str]
) -> float:
    """Return a distance-blind structural matching score."""

    member_term = abs(
        math.log2((int(anchor["member_count"]) + 1) / (int(control["member_count"]) + 1))
    )
    anchor_species = int(anchor["species_count"])
    species_term = abs(anchor_species - int(control["species_count"])) / max(
        1, anchor_species
    )
    copy_term = abs(
        math.log2(
            (float(anchor["mean_copies_per_species"]) + 0.25)
            / (float(control["mean_copies_per_species"]) + 0.25)
        )
    )
    single_copy_term = abs(
        _single_copy_fraction(row=anchor) - _single_copy_fraction(row=control)
    )
    return member_term + (2 * species_term) + copy_term + single_copy_term


def _control_link(
    *,
    anchor: Mapping[str, str],
    control: Mapping[str, str],
    rank: int,
    score: float,
) -> dict[str, Any]:
    """Return one explicit target-to-control matching record."""

    return {
        "run_id": anchor["run_id"],
        "anchor_group_type": anchor["group_type"],
        "anchor_hierarchy_node": anchor["hierarchy_node"],
        "anchor_group_id": anchor["group_id"],
        "control_group_type": control["group_type"],
        "control_hierarchy_node": control["hierarchy_node"],
        "control_group_id": control["group_id"],
        "control_rank": rank,
        "matching_score": score,
        "anchor_member_count": int(anchor["member_count"]),
        "control_member_count": int(control["member_count"]),
        "anchor_species_count": int(anchor["species_count"]),
        "control_species_count": int(control["species_count"]),
        "anchor_mean_copies_per_species": float(
            anchor["mean_copies_per_species"]
        ),
        "control_mean_copies_per_species": float(
            control["mean_copies_per_species"]
        ),
        "anchor_single_copy_fraction": _single_copy_fraction(row=anchor),
        "control_single_copy_fraction": _single_copy_fraction(row=control),
        "control_reused": False,
        "matching_variables": (
            "member_count;species_count;mean_copies_per_species;"
            "single_copy_species_fraction"
        ),
    }


def _single_copy_fraction(*, row: Mapping[str, str]) -> float:
    """Return represented species with exactly one member as a fraction."""

    species = int(row["species_count"])
    return int(row["single_copy_species_count"]) / species if species else 0.0


def _marker_rows(
    *, plan: BenchmarkPlan
) -> tuple[tuple[dict[str, Any], ...], tuple[dict[str, Any], ...]]:
    """Return enriched marker match and complete audit records."""

    markers = plan.marker_authority.by_identifier()
    matches = []
    by_identifier: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for match in plan.marker_selection.match_rows:
        marker = markers[str(match["seed_id"])]
        row = {
            "run_id": match["run_id"],
            "group_type": match["group_type"],
            "hierarchy_node": match["hierarchy_node"],
            "group_id": match["group_id"],
            "legacy_orthogroup_id": match["legacy_orthogroup_id"],
            "marker_id": marker.marker_id,
            "protein_identifier": marker.protein_identifier,
            "protein_entry": marker.protein_entry,
            "marker_name": marker.marker_name,
            "benchmark_class": marker.benchmark_class,
            "benchmark_subclass": marker.benchmark_subclass,
            "domain_architecture": marker.domain_architecture,
            "matched_member_id": match["matched_member_id"],
            "matched_species_label": match["matched_species_label"],
            "match_authority": match["match_authority"],
            "evidence_type": marker.evidence_type,
            "source_title": marker.source_title,
            "source_doi": marker.source_doi,
            "source_table": marker.source_table,
            "source_version": marker.source_version,
            "benchmark_authority_name": plan.marker_authority.source_name,
            "benchmark_authority_sha256": plan.marker_authority.sha256,
        }
        matches.append(row)
        by_identifier[marker.protein_identifier].append(row)
    audit = []
    for marker in plan.marker_authority.records:
        marker_matches = by_identifier.get(marker.protein_identifier, [])
        groups = sorted({str(row["group_id"]) for row in marker_matches})
        members = sorted({str(row["matched_member_id"]) for row in marker_matches})
        audit.append(
            {
                "marker_id": marker.marker_id,
                "protein_identifier": marker.protein_identifier,
                "protein_entry": marker.protein_entry,
                "marker_name": marker.marker_name,
                "benchmark_class": marker.benchmark_class,
                "benchmark_subclass": marker.benchmark_subclass,
                "enabled": marker.enabled,
                "match_status": (
                    "DISABLED"
                    if not marker.enabled
                    else "MATCHED"
                    if marker_matches
                    else "UNMATCHED"
                ),
                "matched_group_count": len(groups),
                "matched_group_ids": ";".join(groups),
                "matched_member_count": len(members),
                "matched_member_ids": ";".join(members),
                "evidence_type": marker.evidence_type,
                "source_title": marker.source_title,
                "source_doi": marker.source_doi,
                "source_table": marker.source_table,
                "source_version": marker.source_version,
                "benchmark_authority_name": plan.marker_authority.source_name,
                "benchmark_authority_sha256": plan.marker_authority.sha256,
            }
        )
    return tuple(matches), tuple(audit)


def _cluster_result_rows(
    *,
    plan: BenchmarkPlan,
    profiles: Sequence[Mapping[str, Any]],
    distance_summaries: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, Any], ...]:
    """Join selected group structure, profiles and distance summaries."""

    statistics = {
        _group_key(row=row): row
        for row in (*plan.target_statistics, *plan.control_statistics)
    }
    summaries = {_group_key(row=row): row for row in distance_summaries}
    missing = sorted(set(statistics).difference(summaries))
    if missing:
        raise InputValidationError(
            "Benchmark groups lack explicit distance summaries: "
            + ";".join(missing[:20])
        )
    by_group: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for profile in profiles:
        by_group[_group_key(row=profile)].append(profile)
    rows = []
    for key in sorted(statistics):
        group = statistics[key]
        summary = summaries[key]
        group_profiles = by_group[key]
        mean = _optional_float(value=summary.get("mean_distance"))
        deviation = _optional_float(
            value=summary.get("population_stddev_distance")
        )
        q25 = _optional_float(value=summary.get("q25_distance"))
        q75 = _optional_float(value=summary.get("q75_distance"))
        total = int(summary.get("total_member_count") or 0)
        sampled = int(summary.get("sampled_member_count") or 0)
        row = {
            "run_id": plan.run_id,
            "group_type": group["group_type"],
            "hierarchy_node": group["hierarchy_node"],
            "group_id": group["group_id"],
            "legacy_orthogroup_id": group.get("legacy_orthogroup_id", ""),
            "gene_tree_parent_clade": group.get("gene_tree_parent_clade", ""),
            "member_count": int(group["member_count"]),
            "species_count": int(group["species_count"]),
            "single_copy_species_count": int(group["single_copy_species_count"]),
            "max_copies_per_species": int(group["max_copies_per_species"]),
            "mean_copies_per_species": float(group["mean_copies_per_species"]),
            "profile_ids": _profile_values(rows=group_profiles, field="profile_id"),
            "profile_classes": _profile_values(
                rows=group_profiles, field="profile_class"
            ),
            "profile_subclasses": _profile_values(
                rows=group_profiles, field="profile_subclass"
            ),
            "membership_roles": _profile_values(
                rows=group_profiles, field="membership_role"
            ),
            "distance_method": summary.get("distance_method", ""),
            "computation_status": summary.get("computation_status", ""),
            "member_identifier_resolution": summary.get(
                "member_identifier_resolution", ""
            ),
            "total_member_count": total,
            "sampled_member_count": sampled,
            "distance_pair_count": int(summary.get("distance_pair_count") or 0),
            "unresolved_pair_count": int(summary.get("unresolved_pair_count") or 0),
            "minimum_distance": summary.get("minimum_distance", ""),
            "q05_distance": summary.get("q05_distance", ""),
            "q25_distance": summary.get("q25_distance", ""),
            "median_distance": summary.get("median_distance", ""),
            "mean_distance": summary.get("mean_distance", ""),
            "q75_distance": summary.get("q75_distance", ""),
            "q95_distance": summary.get("q95_distance", ""),
            "maximum_distance": summary.get("maximum_distance", ""),
            "population_stddev_distance": summary.get(
                "population_stddev_distance", ""
            ),
            "distance_interquartile_range": (
                q75 - q25 if q25 is not None and q75 is not None else ""
            ),
            "distance_coefficient_of_variation": (
                deviation / mean
                if deviation is not None and mean is not None and mean > 0
                else ""
            ),
            "sampling_fraction": sampled / total if total else "",
            "pairwise_rows_persisted": key in plan.target_group_keys,
            "failure_reason": summary.get("failure_reason", ""),
        }
        rows.append(row)
    return tuple(rows)


def _metric_values(
    *, cluster_rows: Sequence[Mapping[str, Any]]
) -> dict[str, dict[str, float]]:
    """Return finite metric values by group key."""

    values: dict[str, dict[str, float]] = {}
    for row in cluster_rows:
        if int(row["distance_pair_count"]) <= 0:
            continue
        metrics = {}
        for metric in DISPERSION_METRICS:
            value = _optional_float(value=row.get(metric))
            if value is not None:
                metrics[metric] = value
        if metrics:
            values[_group_key(row=row)] = metrics
    return values


def _matched_residuals(
    *, plan: BenchmarkPlan, values: Mapping[str, Mapping[str, float]]
) -> dict[str, dict[str, float]]:
    """Return target value minus own matched-control median by metric."""

    controls = _controls_by_anchor(plan=plan)
    residuals = {}
    for target_key in plan.target_group_keys:
        target = values.get(target_key, {})
        target_residuals = {}
        for metric, observed in target.items():
            background = [
                values[key][metric]
                for key in controls.get(target_key, ())
                if metric in values.get(key, {})
            ]
            if len(background) >= MINIMUM_BACKGROUND_GROUPS:
                target_residuals[metric] = observed - statistics.median(background)
        if target_residuals:
            residuals[target_key] = target_residuals
    return residuals


def _background_statistics(
    *,
    run_id: str,
    profiles: Sequence[Mapping[str, Any]],
    values: Mapping[str, Mapping[str, float]],
    residuals: Mapping[str, Mapping[str, float]],
) -> tuple[dict[str, Any], ...]:
    """Summarise raw and matched-residual cluster distributions."""

    profile_groups, metadata = _profile_groups(profiles=profiles)
    rows = []
    for profile_id in sorted(profile_groups):
        for scale, source in (("RAW", values), ("MATCHED_RESIDUAL", residuals)):
            for metric in DISPERSION_METRICS:
                sample = sorted(
                    source[key][metric]
                    for key in profile_groups[profile_id]
                    if metric in source.get(key, {})
                )
                if not sample:
                    continue
                rows.append(
                    {
                        "run_id": run_id,
                        "profile_id": profile_id,
                        "profile_class": metadata[profile_id][0],
                        "profile_subclass": metadata[profile_id][1],
                        "comparison_scale": scale,
                        "metric": metric,
                        "eligible_group_count": len(sample),
                        "minimum_value": sample[0],
                        "q25_value": _quantile(values=sample, fraction=0.25),
                        "median_value": statistics.median(sample),
                        "mean_value": statistics.fmean(sample),
                        "q75_value": _quantile(values=sample, fraction=0.75),
                        "maximum_value": sample[-1],
                        "population_stddev_value": (
                            statistics.pstdev(sample) if len(sample) > 1 else 0.0
                        ),
                    }
                )
    return tuple(rows)


def _profile_contrasts(
    *,
    run_id: str,
    profiles: Sequence[Mapping[str, Any]],
    residuals: Mapping[str, Mapping[str, float]],
    bootstrap_resamples: int,
) -> tuple[dict[str, Any], ...]:
    """Test planned biological-profile contrasts on matched residuals."""

    profile_groups, _ = _profile_groups(profiles=profiles)
    pairs = _planned_profile_pairs(profile_groups=profile_groups)
    rows = []
    for target_id, reference_id in pairs:
        target_keys = profile_groups[target_id]
        if reference_id == "E3_OTHER":
            category_keys = target_keys
            reference_keys = profile_groups["E3_ALL"].difference(category_keys)
        else:
            reference_keys = profile_groups[reference_id]
        for metric in DISPERSION_METRICS:
            target = [
                residuals[key][metric]
                for key in sorted(target_keys)
                if metric in residuals.get(key, {})
            ]
            reference = [
                residuals[key][metric]
                for key in sorted(reference_keys)
                if metric in residuals.get(key, {})
            ]
            status = (
                "TESTED"
                if len(target) >= MINIMUM_BACKGROUND_GROUPS
                and len(reference) >= MINIMUM_BACKGROUND_GROUPS
                else "INSUFFICIENT_GROUPS"
            )
            row: dict[str, Any] = {
                "run_id": run_id,
                "contrast_id": f"{target_id}__VS__{reference_id}",
                "target_profile_id": target_id,
                "reference_profile_id": reference_id,
                "comparison_scale": "MATCHED_RESIDUAL",
                "metric": metric,
                "target_group_count": len(target),
                "reference_group_count": len(reference),
                "target_median": "" if not target else statistics.median(target),
                "reference_median": "" if not reference else statistics.median(reference),
                "median_difference": "",
                "median_difference_ci_low": "",
                "median_difference_ci_high": "",
                "cliffs_delta": "",
                "mann_whitney_u": "",
                "p_value_two_sided": "",
                "fdr_family": f"PROFILE_CONTRAST::{metric}",
                "fdr_q_value": "",
                "status": status,
                "interpretation": (
                    "Positive differences indicate greater dispersion after structural "
                    "matching; negative differences indicate greater compactness."
                ),
            }
            if status == "TESTED":
                u_value, p_value, delta = mann_whitney_test(
                    target=target, reference=reference
                )
                difference = statistics.median(target) - statistics.median(reference)
                low, high = _bootstrap_median_difference(
                    target=target,
                    reference=reference,
                    resamples=bootstrap_resamples,
                    seed_label=(target_id, reference_id, metric),
                )
                row.update(
                    {
                        "median_difference": difference,
                        "median_difference_ci_low": low,
                        "median_difference_ci_high": high,
                        "cliffs_delta": delta,
                        "mann_whitney_u": u_value,
                        "p_value_two_sided": p_value,
                    }
                )
            rows.append(row)
    _apply_fdr(rows=rows, p_field="p_value_two_sided")
    return tuple(rows)


def _individual_comparisons(
    *,
    run_id: str,
    plan: BenchmarkPlan,
    profiles: Sequence[Mapping[str, Any]],
    values: Mapping[str, Mapping[str, float]],
    residuals: Mapping[str, Mapping[str, float]],
) -> tuple[dict[str, Any], ...]:
    """Compare each target cluster against all eligible backgrounds."""

    profile_groups, _ = _profile_groups(profiles=profiles)
    comparison_profiles = [
        profile_id
        for profile_id in sorted(profile_groups)
        if profile_id == "E3_ALL"
        or profile_id in {"HOUSEKEEPING_ALL", "R_NLR_ALL"}
        or profile_id.startswith("E3_CATEGORY::")
        or profile_id.startswith("HOUSEKEEPING_SUBCLASS::")
        or profile_id.startswith("R_NLR_SUBCLASS::")
    ]
    group_lookup = {
        _group_key(row=row): row for row in plan.target_statistics
    }
    controls = _controls_by_anchor(plan=plan)
    rows = []
    for group_key in sorted(plan.target_group_keys):
        group = group_lookup[group_key]
        for metric in DISPERSION_METRICS:
            if metric not in values.get(group_key, {}):
                continue
            observed_raw = values[group_key][metric]
            own_controls = [
                values[key][metric]
                for key in controls.get(group_key, ())
                if metric in values.get(key, {})
            ]
            rows.append(
                _individual_row(
                    run_id=run_id,
                    group=group,
                    background_profile_id="MATCHED_NON_FOCUS",
                    scale="RAW_MATCHED_CONTROL",
                    metric=metric,
                    observed=observed_raw,
                    background=own_controls,
                    leave_one_out=False,
                )
            )
            if metric not in residuals.get(group_key, {}):
                continue
            observed_residual = residuals[group_key][metric]
            for profile_id in comparison_profiles:
                group_keys = set(profile_groups[profile_id])
                leave_one_out = group_key in group_keys
                group_keys.discard(group_key)
                background = [
                    residuals[key][metric]
                    for key in sorted(group_keys)
                    if metric in residuals.get(key, {})
                ]
                rows.append(
                    _individual_row(
                        run_id=run_id,
                        group=group,
                        background_profile_id=profile_id,
                        scale="MATCHED_RESIDUAL",
                        metric=metric,
                        observed=observed_residual,
                        background=background,
                        leave_one_out=leave_one_out,
                    )
                )
    _apply_fdr(rows=rows, p_field="two_sided_p_value")
    return tuple(rows)


def _individual_row(
    *,
    run_id: str,
    group: Mapping[str, Any],
    background_profile_id: str,
    scale: str,
    metric: str,
    observed: float,
    background: Sequence[float],
    leave_one_out: bool,
) -> dict[str, Any]:
    """Return one exact empirical cluster-to-background comparison."""

    status = (
        "TESTED"
        if len(background) >= MINIMUM_BACKGROUND_GROUPS
        else "INSUFFICIENT_GROUPS"
    )
    row = {
        "run_id": run_id,
        "group_type": group["group_type"],
        "hierarchy_node": group["hierarchy_node"],
        "group_id": group["group_id"],
        "background_profile_id": background_profile_id,
        "comparison_scale": scale,
        "metric": metric,
        "observed_value": observed,
        "background_group_count": len(background),
        "background_median": "",
        "difference_from_background_median": "",
        "empirical_percentile": "",
        "lower_tail_p_value": "",
        "upper_tail_p_value": "",
        "two_sided_p_value": "",
        "fdr_family": f"INDIVIDUAL::{background_profile_id}::{metric}",
        "fdr_q_value": "",
        "leave_one_out": leave_one_out,
        "status": status,
    }
    if status != "TESTED":
        return row
    lower = (1 + sum(value <= observed for value in background)) / (
        len(background) + 1
    )
    upper = (1 + sum(value >= observed for value in background)) / (
        len(background) + 1
    )
    percentile = (
        sum(value < observed for value in background)
        + 0.5 * sum(value == observed for value in background)
    ) / len(background)
    median = statistics.median(background)
    row.update(
        {
            "background_median": median,
            "difference_from_background_median": observed - median,
            "empirical_percentile": percentile,
            "lower_tail_p_value": lower,
            "upper_tail_p_value": upper,
            "two_sided_p_value": min(1.0, 2 * min(lower, upper)),
        }
    )
    return row


def _classification_rows(
    *,
    run_id: str,
    plan: BenchmarkPlan,
    values: Mapping[str, Mapping[str, float]],
) -> tuple[dict[str, Any], ...]:
    """Classify central divergence and heterogeneity against own controls."""

    controls = _controls_by_anchor(plan=plan)
    lookup = {_group_key(row=row): row for row in plan.target_statistics}
    rows = []
    for group_key in sorted(plan.target_group_keys):
        group = lookup[group_key]
        mean_result = _classification(
            observed=values.get(group_key, {}).get("mean_distance"),
            background=[
                values[key]["mean_distance"]
                for key in controls.get(group_key, ())
                if "mean_distance" in values.get(key, {})
            ],
        )
        spread_result = _classification(
            observed=values.get(group_key, {}).get(
                "population_stddev_distance"
            ),
            background=[
                values[key]["population_stddev_distance"]
                for key in controls.get(group_key, ())
                if "population_stddev_distance" in values.get(key, {})
            ],
        )
        status = (
            "CLASSIFIED"
            if mean_result[2] != "INSUFFICIENT"
            and spread_result[2] != "INSUFFICIENT"
            else "INSUFFICIENT_CONTROLS"
        )
        rows.append(
            {
                "run_id": run_id,
                "group_type": group["group_type"],
                "hierarchy_node": group["hierarchy_node"],
                "group_id": group["group_id"],
                "matched_control_count": len(controls.get(group_key, ())),
                "mean_distance": values.get(group_key, {}).get("mean_distance", ""),
                "mean_distance_control_median": mean_result[0],
                "mean_distance_control_percentile": mean_result[1],
                "central_divergence_class": mean_result[2],
                "distance_sd": values.get(group_key, {}).get(
                    "population_stddev_distance", ""
                ),
                "distance_sd_control_median": spread_result[0],
                "distance_sd_control_percentile": spread_result[1],
                "heterogeneity_class": spread_result[2],
                "classification_rule": (
                    "<=10th percentile COMPACT/LOW; >=90th percentile "
                    "DISPERSED/HIGH; otherwise TYPICAL"
                ),
                "status": status,
            }
        )
    return tuple(rows)


def _classification(
    *, observed: float | None, background: Sequence[float]
) -> tuple[float | str, float | str, str]:
    """Return background median, percentile and a calibrated class."""

    if observed is None or len(background) < MINIMUM_BACKGROUND_GROUPS:
        return "", "", "INSUFFICIENT"
    percentile = (
        sum(value < observed for value in background)
        + 0.5 * sum(value == observed for value in background)
    ) / len(background)
    if percentile <= 0.1:
        label = "COMPACT_OR_LOW"
    elif percentile >= 0.9:
        label = "DISPERSED_OR_HIGH"
    else:
        label = "TYPICAL"
    return statistics.median(background), percentile, label


def _planned_profile_pairs(
    *, profile_groups: Mapping[str, set[str]]
) -> tuple[tuple[str, str], ...]:
    """Return explicit biological profile contrasts, including E3 categories."""

    pairs = []
    for pair in (
        ("E3_ALL", "HOUSEKEEPING_ALL"),
        ("E3_ALL", "R_NLR_ALL"),
        ("R_NLR_ALL", "HOUSEKEEPING_ALL"),
    ):
        if pair[0] in profile_groups and pair[1] in profile_groups:
            pairs.append(pair)
    for reference in sorted(profile_groups):
        if reference.startswith(("HOUSEKEEPING_SUBCLASS::", "R_NLR_SUBCLASS::")):
            if "E3_ALL" in profile_groups:
                pairs.append(("E3_ALL", reference))
    for profile_id in sorted(profile_groups):
        if not profile_id.startswith("E3_CATEGORY::"):
            continue
        for reference in ("HOUSEKEEPING_ALL", "R_NLR_ALL"):
            if reference in profile_groups:
                pairs.append((profile_id, reference))
        if "E3_ALL" in profile_groups:
            pairs.append((profile_id, "E3_OTHER"))
    return tuple(pairs)


def _profile_groups(
    *, profiles: Sequence[Mapping[str, Any]]
) -> tuple[dict[str, set[str]], dict[str, tuple[str, str]]]:
    """Return profile group sets and class labels."""

    groups: dict[str, set[str]] = defaultdict(set)
    metadata = {}
    for row in profiles:
        profile_id = str(row["profile_id"])
        groups[profile_id].add(_group_key(row=row))
        metadata[profile_id] = (
            str(row["profile_class"]),
            str(row["profile_subclass"]),
        )
    return dict(groups), metadata


def _controls_by_anchor(*, plan: BenchmarkPlan) -> dict[str, tuple[str, ...]]:
    """Return deterministic control group keys by anchor key."""

    controls: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for row in plan.control_links:
        anchor = _group_key(
            row={
                "group_type": row["anchor_group_type"],
                "hierarchy_node": row["anchor_hierarchy_node"],
                "group_id": row["anchor_group_id"],
            }
        )
        control = _group_key(
            row={
                "group_type": row["control_group_type"],
                "hierarchy_node": row["control_hierarchy_node"],
                "group_id": row["control_group_id"],
            }
        )
        controls[anchor].append((int(row["control_rank"]), control))
    return {
        key: tuple(control for _, control in sorted(values))
        for key, values in controls.items()
    }


def _apply_fdr(*, rows: list[dict[str, Any]], p_field: str) -> None:
    """Apply BH correction independently within each declared family."""

    families: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        if row["status"] == "TESTED" and row[p_field] != "":
            families[str(row["fdr_family"])].append(index)
    for indices in families.values():
        adjusted = benjamini_hochberg(
            p_values=[float(rows[index][p_field]) for index in indices]
        )
        for index, q_value in zip(indices, adjusted, strict=True):
            rows[index]["fdr_q_value"] = q_value


def _bootstrap_median_difference(
    *,
    target: Sequence[float],
    reference: Sequence[float],
    resamples: int,
    seed_label: Sequence[str],
) -> tuple[float, float]:
    """Return deterministic percentile interval for a median difference."""

    digest = hashlib.sha256("\x1f".join(seed_label).encode("utf-8")).digest()
    generator = random.Random(int.from_bytes(digest[:8], "big"))
    differences = []
    for _ in range(resamples):
        left = [target[generator.randrange(len(target))] for _ in target]
        right = [reference[generator.randrange(len(reference))] for _ in reference]
        differences.append(statistics.median(left) - statistics.median(right))
    differences.sort()
    return (
        _quantile(values=differences, fraction=0.025),
        _quantile(values=differences, fraction=0.975),
    )


def _categories(*, value: str) -> tuple[str, ...]:
    """Return conservative distinct category labels from catalogue text."""

    labels = {
        label.strip()
        for fragment in value.split(";")
        for label in fragment.split(",")
        if label.strip()
    }
    return tuple(sorted(labels)) or ("UNCLASSIFIED",)


def _profile_values(*, rows: Iterable[Mapping[str, Any]], field: str) -> str:
    """Return sorted semicolon-delimited profile values."""

    return ";".join(sorted({str(row[field]) for row in rows if row.get(field)}))


def _optional_float(*, value: Any) -> float | None:
    """Return a finite float or ``None`` for an empty/non-finite value."""

    if value in {None, ""}:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _quantile(*, values: Sequence[float], fraction: float) -> float:
    """Return a linearly interpolated quantile from sorted finite values."""

    if not values:
        raise ValueError("values must not be empty.")
    position = (len(values) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return float(values[lower])
    weight = position - lower
    return float(values[lower] * (1 - weight) + values[upper] * weight)


def _group_key(*, row: Mapping[str, Any]) -> str:
    """Return one collision-safe group key."""

    return "|".join(
        (
            str(row.get("group_type", "")),
            str(row.get("hierarchy_node", "")),
            str(row.get("group_id", "")),
        )
    )

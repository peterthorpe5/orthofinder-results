"""Matched-background cluster-dispersion benchmark application page."""

from __future__ import annotations

from typing import Any, Mapping

import streamlit as st

from .exports import render_table_downloads
from .figures import (
    benchmark_contrast_figure,
    benchmark_distribution_figure,
    benchmark_individual_figure,
)
from .guidance import render_graph_guidance
from .models import GroupKey
from .queries import BENCHMARK_METRICS, OrthoFinderQueryService

_METRIC_LABELS = {
    "mean_distance": "Average pair distance",
    "median_distance": "Median pair distance",
    "population_stddev_distance": "Pair-distance spread (SD)",
    "distance_interquartile_range": "Pair-distance interquartile range",
    "distance_coefficient_of_variation": "Relative pair-distance spread",
}
_SCALE_LABELS = {
    "RAW": "Raw cluster statistic",
    "MATCHED_RESIDUAL": "Difference from own matched-control median",
}
_PROFILE_HELP = {
    "Biological profile": "Named E3 or reference panel represented by the cluster.",
    "Class": "Broad biological marker class.",
    "Subtype": "E3 category or marker subclass retained from the source authority.",
    "Role": "Biological target or structurally matched non-focus control.",
    "Clusters": "Distinct OrthoFinder clusters assigned to this profile.",
}
_CONTRAST_HELP = {
    "Target profile": "Profile whose residual is subtracted from the reference profile.",
    "Reference profile": "Comparison profile on the right side of the contrast.",
    "Metric": "Cluster-level dispersion statistic tested.",
    "Target clusters": "Eligible target clusters used as statistical replicates.",
    "Reference clusters": "Eligible reference clusters used as statistical replicates.",
    "Median difference": "Target median residual minus reference median residual.",
    "95% CI low": "Lower deterministic bootstrap confidence limit.",
    "95% CI high": "Upper deterministic bootstrap confidence limit.",
    "Cliff's delta": "Non-parametric effect size; sign follows target minus reference.",
    "Mann–Whitney p": "Tie-corrected two-sided cluster-level test p value.",
    "BH-FDR q": "Benjamini–Hochberg adjusted value within the declared metric family.",
    "Status": "TESTED or an explicit reason that inference was unavailable.",
}
_INDIVIDUAL_HELP = {
    "Background": "Empirical profile used to place the selected cluster.",
    "Scale": "Raw metric for own controls or matched residual for pooled profiles.",
    "Metric": "Cluster-level dispersion statistic.",
    "Selected cluster": "Observed statistic or matched-control residual.",
    "Background clusters": "Eligible clusters after leave-one-out when required.",
    "Background median": "Middle value across the eligible background clusters.",
    "Difference": "Selected-cluster value minus the background median.",
    "Percentile": "Empirical rank of the selected cluster within the background.",
    "Two-sided p": "Finite-sample empirical tail probability.",
    "BH-FDR q": "Benjamini–Hochberg value across tested clusters in this family.",
    "Leave one out": "Whether this cluster was removed from its own background.",
    "Status": "TESTED or an explicit insufficient-background state.",
}


def render_dispersion_benchmarks(*, service: OrthoFinderQueryService) -> None:
    """Render pooled and individual matched-background dispersion results.

    Args:
        service: Read-only query service for one validated completed resource.
    """

    st.header("Calibrated dispersion benchmarks")
    st.write(
        "This page tests—not assumes—whether E3, housekeeping-candidate and R/NLR-"
        "candidate clusters differ in evolutionary dispersion. Each OrthoFinder cluster "
        "is one statistical replicate; individual protein pairs are never treated as "
        "independent observations."
    )
    required = (
        "benchmark_group_profiles",
        "benchmark_cluster_results",
        "benchmark_contrasts",
        "benchmark_individual_comparisons",
        "benchmark_cluster_classifications",
    )
    missing = tuple(
        relation for relation in required if not service.has_relation(relation=relation)
    )
    if missing:
        st.info(
            "This resource predates matched-background dispersion benchmarking. Run the "
            "dispersion-benchmark cluster workflow and open its completed resource. Missing "
            f"relations: {', '.join(missing)}."
        )
        return
    with st.expander("Study design, hypotheses and safeguards", expanded=True):
        st.markdown(
            """
            - **Housekeeping candidates and R/NLR candidates are comparison panels, not
              predetermined answers.** Either panel may prove compact, dispersed or mixed.
            - Each target cluster receives unique non-focus controls matched without looking at
              distance: protein count, represented-species count, mean copies per species and
              single-copy-species fraction.
            - A **matched residual** is the target value minus the median of its own controls.
              Positive means more dispersed than those controls; negative means more compact.
            - Planned profile tests report effect size, deterministic bootstrap confidence
              intervals, Mann–Whitney p values and Benjamini–Hochberg FDR q values.
            - Individual-cluster profile tests use leave-one-out where necessary. Missing or
              insufficient results remain explicit and are never converted to biological zero.
            """
        )
    profiles = service.benchmark_profiles()
    profile_rows = tuple(_display_profile(row=row) for row in profiles)
    st.subheader("Available biological profiles")
    st.dataframe(profile_rows, width="stretch", hide_index=True)
    render_table_downloads(
        records=profile_rows,
        file_stem="orthofinder_dispersion_benchmark_profiles",
        key="benchmark_profiles_download",
        tsv_label="Download profile inventory as TSV",
        excel_label="Download profile inventory as formatted Excel",
        column_definitions=_PROFILE_HELP,
        workbook_title="OrthoFinder dispersion profiles",
    )

    st.subheader("Compare biological backgrounds")
    controls = st.columns(2)
    metric_label = controls[0].selectbox(
        "Dispersion statistic",
        tuple(_METRIC_LABELS.values()),
        help=(
            "Average and median describe central divergence; SD and interquartile range "
            "describe heterogeneity; relative spread divides SD by the mean."
        ),
    )
    metric = _metric_from_label(label=metric_label)
    scale_label = controls[1].selectbox(
        "Comparison scale",
        tuple(_SCALE_LABELS.values()),
        index=1,
        help=(
            "Matched residuals adjust every biological target for its own structurally "
            "matched controls. Raw values retain the original distance scale."
        ),
    )
    scale = next(key for key, value in _SCALE_LABELS.items() if value == scale_label)
    available_profiles = tuple(str(row["profile_id"]) for row in profiles)
    defaults = tuple(
        profile
        for profile in ("E3_ALL", "HOUSEKEEPING_ALL", "R_NLR_ALL")
        if profile in available_profiles
    )
    selected_profiles = tuple(
        st.multiselect(
            "Profiles shown",
            available_profiles,
            default=defaults,
            help="A cluster can appear in more than one profile when marker classes overlap.",
        )
    )
    distribution_rows = service.benchmark_distribution(
        profile_ids=selected_profiles,
        metric=metric,
        comparison_scale=scale,
    )
    if distribution_rows:
        render_graph_guidance(key="benchmark_distribution")
        st.plotly_chart(
            benchmark_distribution_figure(
                rows=distribution_rows,
                metric_label=metric_label,
                scale_label=scale_label,
            ),
            width="stretch",
            config={"displaylogo": False},
        )
    else:
        st.warning("No eligible cluster-level values match these profile controls.")

    contrast_rows = service.benchmark_contrasts(metric=metric)
    st.subheader("Planned profile contrasts")
    st.caption(
        "These tests always use matched-control residuals, regardless of the distribution "
        "scale selected above. Red forest-plot points have BH-FDR q≤0.05."
    )
    if any(row["status"] == "TESTED" for row in contrast_rows):
        render_graph_guidance(key="benchmark_contrast")
        st.plotly_chart(
            benchmark_contrast_figure(rows=contrast_rows, metric_label=metric_label),
            width="stretch",
            config={"displaylogo": False},
        )
    displayed_contrasts = tuple(
        _display_contrast(row=row) for row in contrast_rows
    )
    st.dataframe(displayed_contrasts, width="stretch", hide_index=True)
    render_table_downloads(
        records=displayed_contrasts,
        file_stem=f"orthofinder_dispersion_contrasts_{metric}",
        key=f"benchmark_contrasts_download_{metric}",
        tsv_label="Download these profile contrasts as TSV",
        excel_label="Download these profile contrasts as formatted Excel",
        column_definitions=_CONTRAST_HELP,
        workbook_title="OrthoFinder dispersion contrasts",
    )

    st.subheader("Test one cluster against every background")
    catalogue = service.benchmark_cluster_catalogue()
    if not catalogue:
        st.warning("No biological target clusters are available for individual tests.")
        return
    labels = tuple(_cluster_label(row=row) for row in catalogue)
    selected_label = st.selectbox(
        "Target cluster",
        labels,
        help=(
            "Includes every E3, housekeeping-candidate and R/NLR-candidate cluster. "
            "Use the exact composite group identity when exporting or citing it."
        ),
    )
    selected_index = labels.index(selected_label)
    selected = catalogue[selected_index]
    key = GroupKey(
        run_id=str(selected["run_id"]),
        group_type=str(selected["group_type"]),
        hierarchy_node=str(selected["hierarchy_node"]),
        group_id=str(selected["group_id"]),
    )
    individual_rows = service.benchmark_individual_comparisons(
        key=key,
        metric=metric,
    )
    tabs = st.tabs(("Own matched controls", "Biological profile backgrounds"))
    scales = (
        (tabs[0], "RAW_MATCHED_CONTROL", "Raw cluster statistic"),
        (
            tabs[1],
            "MATCHED_RESIDUAL",
            "Difference from own matched-control median",
        ),
    )
    for tab, row_scale, row_scale_label in scales:
        with tab:
            rows = tuple(
                row
                for row in individual_rows
                if row["comparison_scale"] == row_scale
            )
            tested = tuple(row for row in rows if row["status"] == "TESTED")
            if tested:
                render_graph_guidance(key="benchmark_individual")
                st.plotly_chart(
                    benchmark_individual_figure(
                        rows=tested,
                        metric_label=metric_label,
                        scale_label=row_scale_label,
                    ),
                    width="stretch",
                    config={"displaylogo": False},
                )
            else:
                st.info(
                    "No comparison at this scale has the minimum three eligible "
                    "background clusters. The explicit insufficient rows remain below."
                )
            st.dataframe(
                tuple(_display_individual(row=row) for row in rows),
                width="stretch",
                hide_index=True,
            )
    displayed_individual = tuple(
        _display_individual(row=row) for row in individual_rows
    )
    render_table_downloads(
        records=displayed_individual,
        file_stem=(
            f"orthofinder_individual_dispersion_{key.group_id}_{metric}"
        ),
        key=f"benchmark_individual_download_{key.group_id}_{metric}",
        tsv_label="Download all individual comparisons as TSV",
        excel_label="Download all individual comparisons as formatted Excel",
        column_definitions=_INDIVIDUAL_HELP,
        workbook_title="OrthoFinder individual dispersion tests",
    )


def _metric_from_label(*, label: str) -> str:
    """Return the exact metric key for one user-facing label."""

    for metric in BENCHMARK_METRICS:
        if _METRIC_LABELS[metric] == label:
            return metric
    raise ValueError(f"Unknown benchmark metric label: {label}")


def _display_profile(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one profile row with plain-language headings."""

    return {
        "Biological profile": row["profile_id"],
        "Class": row["profile_class"],
        "Subtype": row["profile_subclass"],
        "Role": row["membership_role"],
        "Clusters": row["group_count"],
    }


def _display_contrast(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one contrast row with plain-language headings."""

    return {
        "Target profile": row["target_profile_id"],
        "Reference profile": row["reference_profile_id"],
        "Metric": _METRIC_LABELS[str(row["metric"])],
        "Target clusters": row["target_group_count"],
        "Reference clusters": row["reference_group_count"],
        "Median difference": row["median_difference"],
        "95% CI low": row["median_difference_ci_low"],
        "95% CI high": row["median_difference_ci_high"],
        "Cliff's delta": row["cliffs_delta"],
        "Mann–Whitney p": row["p_value_two_sided"],
        "BH-FDR q": row["fdr_q_value"],
        "Status": row["status"],
    }


def _display_individual(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one individual comparison row with intuitive headings."""

    return {
        "Background": row["background_profile_id"],
        "Scale": row["comparison_scale"],
        "Metric": _METRIC_LABELS[str(row["metric"])],
        "Selected cluster": row["observed_value"],
        "Background clusters": row["background_group_count"],
        "Background median": row["background_median"],
        "Difference": row["difference_from_background_median"],
        "Percentile": row["empirical_percentile"],
        "Two-sided p": row["two_sided_p_value"],
        "BH-FDR q": row["fdr_q_value"],
        "Leave one out": row["leave_one_out"],
        "Status": row["status"],
    }


def _cluster_label(*, row: Mapping[str, Any]) -> str:
    """Return an informative stable selector label for one benchmark target."""

    return (
        f"{row['group_type']} | {row['hierarchy_node'] or 'ROOT'} | "
        f"{row['group_id']} — {row['profile_ids']}"
    )

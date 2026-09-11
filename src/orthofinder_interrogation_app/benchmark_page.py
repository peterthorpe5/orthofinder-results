"""Matched-background cluster-dispersion results and marker browser."""

from __future__ import annotations

import re
from collections import Counter
from typing import Any, Mapping, Sequence

import streamlit as st

from .benchmark_labels import benchmark_class_label, benchmark_profile_label
from .documentation_page import render_page_guidance
from .evolutionary_page import _store_active_group
from .exports import render_plotly_figure, render_table_downloads
from .figures import (
    benchmark_classification_figure,
    benchmark_contrast_figure,
    benchmark_distribution_figure,
    benchmark_individual_figure,
    benchmark_marker_coverage_figure,
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
_BROAD_PROFILES = ("E3_ALL", "HOUSEKEEPING_ALL", "R_NLR_ALL")
_HEADLINE_PAIRS = (
    ("E3_ALL", "HOUSEKEEPING_ALL"),
    ("E3_ALL", "R_NLR_ALL"),
    ("R_NLR_ALL", "HOUSEKEEPING_ALL"),
)
_PROFILE_HELP = {
    "Biological profile": "Readable name for the E3 or reference panel.",
    "Profile ID": "Exact machine-readable profile identifier retained in result files.",
    "Class": "Broad biological marker class.",
    "Subtype": "E3 category or marker subclass retained from the source authority.",
    "Role": "Biological target or structurally matched non-focus control.",
    "Clusters": "Distinct OrthoFinder clusters assigned to this profile.",
}
_MARKER_HELP = {
    "Biological profile": "Selected biological panel containing this marker match.",
    "Profile ID": "Exact machine-readable profile identifier.",
    "Class": "Broad E3, housekeeping-reference or R/NLR-reference class.",
    "Subtype": "E3 category or reviewed comparison-panel subclass.",
    "Marker gene / locus": "Authority marker ID, such as an Arabidopsis locus.",
    "Protein accession": "Reference protein accession supplied by the authority.",
    "Protein entry": "Reviewed protein entry name when available.",
    "Protein name": "Protein description supplied by the versioned authority.",
    "OrthoFinder group": "Exact group containing the matched protein.",
    "Matched protein": "Exact protein identifier found in this OrthoFinder run.",
    "Species": "Exact OrthoFinder species label for the matched protein.",
    "Domain architecture": "Reviewed R/NLR architecture when supplied.",
    "Evidence": "Evidence class recorded by the marker authority.",
    "Source": "Source title or E3 catalogue provenance.",
    "DOI": "Source DOI when one was supplied.",
    "Source version": "Version or release of the marker source.",
    "Authority": "Filename or authority name used by the completed run.",
}
_MEMBER_HELP = {
    "Species": "Exact OrthoFinder species label for this cluster member.",
    "Protein ID": "Exact protein or sequence identifier in the cluster.",
    "Defines selected profile": (
        "True when this protein is one of the matched markers defining the selected profile."
    ),
    "Parent legacy orthogroup": "Linked flat orthogroup when OrthoFinder supplied it.",
    "Gene-tree parent clade": "OrthoFinder gene-tree clade used to define this HOG.",
}
_CLASSIFICATION_HELP = {
    "Group": "Exact OrthoFinder group identifier.",
    "Biological class": "One or several biological panels represented by this cluster.",
    "Profiles": "Exact machine-readable profile memberships.",
    "Average pair distance": "Observed average distance in the analysed matrix.",
    "Average-distance percentile": (
        "Rank relative to this cluster's own structurally matched controls."
    ),
    "Central divergence class": "Compact, typical or dispersed matched-control label.",
    "Distance spread (SD)": "Population SD of pair distances in the analysed matrix.",
    "Spread percentile": "Distance-SD rank relative to this cluster's own controls.",
    "Heterogeneity class": "Low, typical or high matched-control spread label.",
    "Matched controls": "Number of distance-blind matched controls used.",
    "Status": "CLASSIFIED or an explicit insufficient-control state.",
}
_CLASSIFICATION_COUNT_HELP = {
    "Central divergence class": "Matched-control class based on average pair distance.",
    "Heterogeneity class": "Matched-control class based on pair-distance population SD.",
    "Clusters": "Number of biological target clusters in this class combination.",
}
_CONTRAST_HELP = {
    "Target profile": "Readable profile whose residual is compared with the reference.",
    "Target profile ID": "Exact machine-readable target profile.",
    "Reference profile": "Readable comparison profile.",
    "Reference profile ID": "Exact machine-readable reference profile.",
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
    "Background": "Readable empirical profile used to place the selected cluster.",
    "Background ID": "Exact machine-readable background profile.",
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
    """Render clear benchmark findings, markers and detailed statistics.

    Args:
        service: Read-only query service for one validated completed resource.
    """

    st.header("Calibrated dispersion results")
    render_page_guidance(key="dispersion_benchmarks")
    st.write(
        "Compare E3 clusters with housekeeping-reference, R/NLR-reference and "
        "structure-matched non-focus clusters. The page leads with the biological result; "
        "the complete statistics and provenance remain available underneath."
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

    profiles = service.benchmark_profiles()
    catalogue = service.benchmark_cluster_catalogue()
    classifications = service.benchmark_classifications()
    metric_label = st.selectbox(
        "Dispersion measure used throughout this page",
        tuple(_METRIC_LABELS.values()),
        help=(
            "Average and median describe central divergence. SD and interquartile range "
            "describe heterogeneity. Relative spread divides SD by the mean."
        ),
    )
    metric = _metric_from_label(label=metric_label)
    contrast_rows = service.benchmark_contrasts(metric=metric)

    summary_tab, genes_tab, details_tab = st.tabs(
        ("Results at a glance", "Genes and clusters", "Detailed statistics")
    )
    with summary_tab:
        _render_results_summary(
            service=service,
            profiles=profiles,
            catalogue=catalogue,
            classifications=classifications,
            contrast_rows=contrast_rows,
            metric=metric,
            metric_label=metric_label,
        )
    with genes_tab:
        _render_marker_browser(service=service, profiles=profiles)
    with details_tab:
        _render_detailed_results(
            service=service,
            profiles=profiles,
            catalogue=catalogue,
            contrast_rows=contrast_rows,
            metric=metric,
            metric_label=metric_label,
        )


def _render_results_summary(
    *,
    service: OrthoFinderQueryService,
    profiles: Sequence[Mapping[str, Any]],
    catalogue: Sequence[Mapping[str, Any]],
    classifications: Sequence[Mapping[str, Any]],
    contrast_rows: Sequence[Mapping[str, Any]],
    metric: str,
    metric_label: str,
) -> None:
    """Render plain-language findings and two complementary interactive graphs."""

    st.subheader("What was analysed?")
    columns = st.columns(5)
    columns[0].metric(
        "Biological clusters",
        f"{len(catalogue):,}",
        help="Distinct E3, housekeeping-reference or R/NLR-reference target clusters.",
    )
    for column, profile_id, label in zip(
        columns[1:4],
        _BROAD_PROFILES,
        ("E3 clusters", "Housekeeping-reference clusters", "R/NLR-reference clusters"),
        strict=True,
    ):
        column.metric(
            label,
            f"{_profile_group_count(profiles=profiles, profile_id=profile_id):,}",
            help="A cluster can occur in more than one profile if marker classes overlap.",
        )
    significant = sum(
        row.get("status") == "TESTED"
        and row.get("fdr_q_value") is not None
        and float(row["fdr_q_value"]) <= 0.05
        for row in contrast_rows
    )
    tested = sum(row.get("status") == "TESTED" for row in contrast_rows)
    columns[4].metric(
        "FDR-significant contrasts",
        f"{significant:,} of {tested:,}",
        help=f"Planned {metric_label.lower()} contrasts with BH-FDR q≤0.05.",
    )

    st.subheader("Headline comparisons")
    st.caption(
        "These statements use cluster-level matched-control residuals. Positive differences "
        "mean that the target profile is more dispersed than the reference after calibration."
    )
    headline_rows = _headline_contrasts(rows=contrast_rows)
    if not headline_rows:
        st.info("No broad E3, housekeeping and R/NLR contrasts were available for this measure.")
    for row in headline_rows:
        message, significant_result = _contrast_summary(
            row=row,
            metric_label=metric_label,
        )
        if significant_result:
            st.success(message)
        else:
            st.info(message)

    broad_profiles = tuple(
        profile_id
        for profile_id in _BROAD_PROFILES
        if any(str(row["profile_id"]) == profile_id for row in profiles)
    )
    st.subheader("How do the biological profiles compare?")
    distribution_rows = service.benchmark_distribution(
        profile_ids=broad_profiles,
        metric=metric,
        comparison_scale="MATCHED_RESIDUAL",
    )
    if distribution_rows:
        render_graph_guidance(key="benchmark_distribution")
        render_plotly_figure(
            figure=benchmark_distribution_figure(
                rows=distribution_rows,
                metric_label=metric_label,
                scale_label=_SCALE_LABELS["MATCHED_RESIDUAL"],
            ),
            file_stem=f"orthofinder_broad_profile_distribution_{metric}",
            key=f"benchmark_summary_distribution_{metric}",
        )
    else:
        st.warning("No matched-control residuals were available for the broad profiles.")

    st.subheader("Which individual clusters are compact, dispersed or heterogeneous?")
    if classifications:
        render_graph_guidance(key="benchmark_classification")
        render_plotly_figure(
            figure=benchmark_classification_figure(rows=classifications),
            file_stem="orthofinder_matched_control_cluster_map",
            key="benchmark_summary_classification",
        )
        displayed = tuple(_display_classification(row=row) for row in classifications)
        summary = Counter(
            (
                str(row.get("central_divergence_class", "UNAVAILABLE")),
                str(row.get("heterogeneity_class", "UNAVAILABLE")),
            )
            for row in classifications
        )
        with st.expander("Classification counts and complete cluster table"):
            count_rows = tuple(
                {
                    "Central divergence class": central.replace("_OR_", " / ").title(),
                    "Heterogeneity class": spread.replace("_OR_", " / ").title(),
                    "Clusters": count,
                }
                for (central, spread), count in sorted(summary.items())
            )
            st.dataframe(count_rows, width="stretch", hide_index=True)
            render_table_downloads(
                records=count_rows,
                file_stem="orthofinder_dispersion_classification_counts",
                key="benchmark_classification_count_download",
                tsv_label="Download classification counts as TSV",
                excel_label="Download classification counts as formatted Excel",
                column_definitions=_CLASSIFICATION_COUNT_HELP,
                workbook_title="OrthoFinder dispersion classification counts",
            )
            st.dataframe(
                displayed,
                width="stretch",
                hide_index=True,
                column_config=_column_config(descriptions=_CLASSIFICATION_HELP),
            )
        render_table_downloads(
            records=displayed,
            file_stem="orthofinder_dispersion_cluster_classifications",
            key="benchmark_classification_download",
            tsv_label="Download every cluster classification as TSV",
            excel_label="Download every cluster classification as formatted Excel",
            column_definitions=_CLASSIFICATION_HELP,
            workbook_title="OrthoFinder dispersion classifications",
        )
    else:
        st.info("No matched-control cluster classifications are available.")

    with st.expander("Study design, hypotheses and safeguards"):
        st.markdown(
            """
            - **Housekeeping and R/NLR panels are comparisons, not predetermined answers.**
              Either may prove compact, dispersed or mixed.
            - Each target cluster has unique non-focus controls matched without using distance:
              protein count, represented-species count, mean copies per species and
              single-copy-species fraction.
            - A **matched residual** is the target value minus the median of its own controls.
            - Profile tests use one cluster as one replicate and report bootstrap intervals,
              Mann–Whitney p values, Cliff's delta and Benjamini–Hochberg FDR q values.
            - Individual protein pairs are never treated as independent statistical replicates.
            """
        )


def _render_marker_browser(
    *,
    service: OrthoFinderQueryService,
    profiles: Sequence[Mapping[str, Any]],
) -> None:
    """Render profile-defining genes and complete selected-cluster membership."""

    st.subheader("Genes and proteins defining each biological profile")
    st.write(
        "Choose a profile to see the exact E3 seeds or Arabidopsis reference markers that "
        "matched this OrthoFinder run. This is the direct answer to, for example, “which "
        "housekeeping genes were used, and which clusters contain them?”"
    )
    available = tuple(
        str(row["profile_id"])
        for row in profiles
        if str(row["membership_role"]) == "TARGET"
    )
    if not available:
        st.info("No biological target profiles are stored in this resource.")
        return
    default_profile = "HOUSEKEEPING_ALL" if "HOUSEKEEPING_ALL" in available else available[0]
    selected_profile = st.selectbox(
        "Biological profile or gene set",
        available,
        index=available.index(default_profile),
        format_func=lambda value: benchmark_profile_label(profile_id=value),
        help=(
            "Pooled profiles provide the complete E3, housekeeping or R/NLR set. "
            "Subclass profiles restrict the table to that exact category."
        ),
    )
    marker_rows = service.benchmark_profile_markers(profile_ids=(selected_profile,))
    if not marker_rows:
        st.info(
            "This profile has no marker-level relation in the opened resource. Its cluster-level "
            "statistics remain available in Detailed statistics."
        )
        return
    marker_ids = {str(row["marker_id"]) for row in marker_rows}
    group_keys = {
        (str(row["group_type"]), str(row["hierarchy_node"]), str(row["group_id"]))
        for row in marker_rows
    }
    matched_proteins = {str(row["matched_member_id"]) for row in marker_rows}
    metrics = st.columns(3)
    metrics[0].metric("Marker genes / proteins", f"{len(marker_ids):,}")
    metrics[1].metric("Matched OrthoFinder clusters", f"{len(group_keys):,}")
    metrics[2].metric("Matched proteins in this run", f"{len(matched_proteins):,}")
    render_graph_guidance(key="benchmark_marker_coverage")
    render_plotly_figure(
        figure=benchmark_marker_coverage_figure(rows=marker_rows),
        file_stem=(
            "orthofinder_marker_coverage_"
            f"{_safe_file_token(value=selected_profile)}"
        ),
        key=f"benchmark_marker_coverage_{_safe_file_token(value=selected_profile)}",
    )
    if len(marker_ids) > 40:
        st.caption(
            "The graph shows the 40 markers matching the most distinct clusters; "
            "the table and downloads below retain every matched marker."
        )
    displayed_markers = tuple(_display_marker(row=row) for row in marker_rows)
    st.dataframe(
        displayed_markers,
        width="stretch",
        hide_index=True,
        column_config=_column_config(descriptions=_MARKER_HELP),
    )
    token = _safe_file_token(value=selected_profile)
    render_table_downloads(
        records=displayed_markers,
        file_stem=f"orthofinder_{token}_matched_genes_and_clusters",
        key=f"benchmark_markers_download_{token}",
        tsv_label="Download these genes and matched clusters as TSV",
        excel_label="Download these genes and matched clusters as formatted Excel",
        column_definitions=_MARKER_HELP,
        workbook_title=(
            "OrthoFinder markers: "
            f"{benchmark_profile_label(profile_id=selected_profile)}"
        ),
    )

    st.subheader("Inspect one marker's complete cluster")
    labels = tuple(_marker_match_label(row=row) for row in marker_rows)
    selected_label = st.selectbox(
        "Marker and matched cluster",
        labels,
        help="The same marker may appear more than once if it matched several proteins or groups.",
    )
    selected = marker_rows[labels.index(selected_label)]
    key = _group_key(row=selected, run_id=service.resource.run_id)
    group = service.get_group(key=key)
    group_metrics = st.columns(4)
    group_metrics[0].metric("Proteins in cluster", f"{int(group['member_count']):,}")
    group_metrics[1].metric("Species represented", f"{int(group['species_count']):,}")
    group_metrics[2].metric(
        "Stored average distance",
        _format_number(value=group.get("mean_distance")),
    )
    group_metrics[3].metric(
        "Stored distance spread (SD)",
        _format_number(value=group.get("population_stddev_distance")),
    )
    st.button(
        "Open this cluster's complete visual suite",
        on_click=_open_group_in_explorer,
        kwargs={"key": key},
    )
    related_marker_members = {
        str(row["matched_member_id"])
        for row in marker_rows
        if str(row["group_type"]) == key.group_type
        and str(row["hierarchy_node"]) == key.hierarchy_node
        and str(row["group_id"]) == key.group_id
    }
    members = service.get_group_members(key=key)
    displayed_members = tuple(
        _display_member(row=row, marker_members=related_marker_members) for row in members
    )
    st.dataframe(
        displayed_members,
        width="stretch",
        hide_index=True,
        column_config=_column_config(descriptions=_MEMBER_HELP),
    )
    if len(members) == int(group["member_count"]):
        render_table_downloads(
            records=displayed_members,
            file_stem=f"{key.group_id}_{token}_cluster_members",
            key=f"benchmark_members_download_{key.group_id}_{token}",
            tsv_label="Download every protein in this cluster as TSV",
            excel_label="Download every protein in this cluster as formatted Excel",
            column_definitions=_MEMBER_HELP,
            workbook_title=f"Cluster members: {key.group_id}",
        )
    else:
        st.warning(
            "The cluster exceeds the safe interactive membership limit. Use its DuckDB or "
            "Parquet membership relation for a complete programmatic export."
        )


def _render_detailed_results(
    *,
    service: OrthoFinderQueryService,
    profiles: Sequence[Mapping[str, Any]],
    catalogue: Sequence[Mapping[str, Any]],
    contrast_rows: Sequence[Mapping[str, Any]],
    metric: str,
    metric_label: str,
) -> None:
    """Render complete profile and individual-cluster statistical controls."""

    profile_rows = tuple(_display_profile(row=row) for row in profiles)
    with st.expander("Available biological profiles and exact profile IDs"):
        st.dataframe(
            profile_rows,
            width="stretch",
            hide_index=True,
            column_config=_column_config(descriptions=_PROFILE_HELP),
        )
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
    scale_label = controls[0].selectbox(
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
    defaults = tuple(profile for profile in _BROAD_PROFILES if profile in available_profiles)
    selected_profiles = tuple(
        controls[1].multiselect(
            "Profiles shown",
            available_profiles,
            default=defaults,
            format_func=lambda value: benchmark_profile_label(profile_id=value),
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
        render_plotly_figure(
            figure=benchmark_distribution_figure(
                rows=distribution_rows,
                metric_label=metric_label,
                scale_label=scale_label,
            ),
            file_stem=f"orthofinder_profile_distribution_{metric}_{scale}",
            key=f"benchmark_detailed_distribution_{metric}_{scale}",
        )
    else:
        st.warning("No eligible cluster-level values match these profile controls.")

    st.subheader("Planned profile contrasts")
    st.caption(
        "The forest plot now starts with the broad comparisons rather than all 420 rows. "
        "Choose additional targets or references to expand it. Red points have BH-FDR q≤0.05."
    )
    target_options = tuple(
        dict.fromkeys(str(row["target_profile_id"]) for row in contrast_rows)
    )
    reference_options = tuple(
        dict.fromkeys(str(row["reference_profile_id"]) for row in contrast_rows)
    )
    filter_columns = st.columns(3)
    target_default = ("E3_ALL",) if "E3_ALL" in target_options else target_options[:1]
    reference_default = tuple(
        value for value in ("HOUSEKEEPING_ALL", "R_NLR_ALL") if value in reference_options
    )
    selected_targets = tuple(
        filter_columns[0].multiselect(
            "Target profiles",
            target_options,
            default=target_default,
            format_func=lambda value: benchmark_profile_label(profile_id=value),
        )
    )
    selected_references = tuple(
        filter_columns[1].multiselect(
            "Reference profiles",
            reference_options,
            default=reference_default or reference_options[:2],
            format_func=lambda value: benchmark_profile_label(profile_id=value),
        )
    )
    significant_only = filter_columns[2].checkbox(
        "Show only BH-FDR q≤0.05",
        value=False,
        help="Insufficient and non-significant rows remain in the complete download below.",
    )
    filtered_contrasts = tuple(
        row
        for row in contrast_rows
        if str(row["target_profile_id"]) in selected_targets
        and str(row["reference_profile_id"]) in selected_references
        and (
            not significant_only
            or (
                row.get("fdr_q_value") is not None
                and float(row["fdr_q_value"]) <= 0.05
            )
        )
    )
    if any(row["status"] == "TESTED" for row in filtered_contrasts):
        render_graph_guidance(key="benchmark_contrast")
        render_plotly_figure(
            figure=benchmark_contrast_figure(
                rows=filtered_contrasts,
                metric_label=metric_label,
            ),
            file_stem=f"orthofinder_profile_contrasts_{metric}",
            key=f"benchmark_contrast_forest_{metric}",
        )
    else:
        st.info("No tested contrasts match the current target and reference filters.")
    displayed_filtered = tuple(_display_contrast(row=row) for row in filtered_contrasts)
    if displayed_filtered:
        st.dataframe(
            displayed_filtered,
            width="stretch",
            hide_index=True,
            column_config=_column_config(descriptions=_CONTRAST_HELP),
        )
        render_table_downloads(
            records=displayed_filtered,
            file_stem=f"orthofinder_dispersion_filtered_contrasts_{metric}",
            key=f"benchmark_filtered_contrasts_download_{metric}",
            tsv_label="Download displayed contrasts as TSV",
            excel_label="Download displayed contrasts as formatted Excel",
            column_definitions=_CONTRAST_HELP,
            workbook_title="Displayed OrthoFinder dispersion contrasts",
        )
    displayed_contrasts = tuple(_display_contrast(row=row) for row in contrast_rows)
    render_table_downloads(
        records=displayed_contrasts,
        file_stem=f"orthofinder_dispersion_all_contrasts_{metric}",
        key=f"benchmark_contrasts_download_{metric}",
        tsv_label="Download all contrasts for this measure as TSV",
        excel_label="Download all contrasts for this measure as formatted Excel",
        column_definitions=_CONTRAST_HELP,
        workbook_title="OrthoFinder dispersion contrasts",
    )

    st.subheader("Test one cluster against every background")
    if not catalogue:
        st.warning("No biological target clusters are available for individual tests.")
        return
    labels = tuple(_cluster_label(row=row) for row in catalogue)
    selected_label = st.selectbox(
        "Target cluster",
        labels,
        help=(
            "Includes every E3, housekeeping-reference and R/NLR-reference cluster. "
            "Use the exact composite group identity when exporting or citing it."
        ),
    )
    selected = catalogue[labels.index(selected_label)]
    key = _group_key(row=selected, run_id=str(selected["run_id"]))
    individual_rows = service.benchmark_individual_comparisons(key=key, metric=metric)
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
                row for row in individual_rows if row["comparison_scale"] == row_scale
            )
            tested_rows = tuple(row for row in rows if row["status"] == "TESTED")
            if tested_rows:
                render_graph_guidance(key="benchmark_individual")
                figure_key = (
                    "benchmark_individual_"
                    f"{key.group_type}_{key.hierarchy_node}_{key.group_id}_"
                    f"{metric}_{row_scale}"
                )
                render_plotly_figure(
                    figure=benchmark_individual_figure(
                        rows=tested_rows,
                        metric_label=metric_label,
                        scale_label=row_scale_label,
                    ),
                    file_stem=(
                        f"orthofinder_individual_{key.group_id}_{metric}_{row_scale}"
                    ),
                    key=figure_key,
                )
            else:
                st.info(
                    "No comparison at this scale has the minimum three eligible "
                    "background clusters. Explicit insufficient rows remain below."
                )
            displayed_scale = tuple(_display_individual(row=row) for row in rows)
            if displayed_scale:
                st.dataframe(
                    displayed_scale,
                    width="stretch",
                    hide_index=True,
                    column_config=_column_config(descriptions=_INDIVIDUAL_HELP),
                )
                render_table_downloads(
                    records=displayed_scale,
                    file_stem=(
                        f"orthofinder_individual_{key.group_id}_{metric}_{row_scale}"
                    ),
                    key=(
                        f"benchmark_individual_table_{key.group_id}_{metric}_{row_scale}"
                    ),
                    tsv_label="Download displayed comparisons as TSV",
                    excel_label="Download displayed comparisons as formatted Excel",
                    column_definitions=_INDIVIDUAL_HELP,
                    workbook_title="Displayed individual dispersion comparisons",
                )
    displayed_individual = tuple(
        _display_individual(row=row) for row in individual_rows
    )
    render_table_downloads(
        records=displayed_individual,
        file_stem=f"orthofinder_individual_dispersion_{key.group_id}_{metric}",
        key=f"benchmark_individual_download_{key.group_id}_{metric}",
        tsv_label="Download all individual comparisons as TSV",
        excel_label="Download all individual comparisons as formatted Excel",
        column_definitions=_INDIVIDUAL_HELP,
        workbook_title="OrthoFinder individual dispersion tests",
    )


def _headline_contrasts(
    *, rows: Sequence[Mapping[str, Any]]
) -> tuple[Mapping[str, Any], ...]:
    """Return broad biological contrasts in a stable interpretation order."""

    indexed = {
        (str(row["target_profile_id"]), str(row["reference_profile_id"])): row
        for row in rows
    }
    return tuple(indexed[pair] for pair in _HEADLINE_PAIRS if pair in indexed)


def _contrast_summary(
    *, row: Mapping[str, Any], metric_label: str
) -> tuple[str, bool]:
    """Return a plain-language result and whether it passed FDR control."""

    target = benchmark_profile_label(profile_id=str(row["target_profile_id"]))
    reference = benchmark_profile_label(profile_id=str(row["reference_profile_id"]))
    if row.get("status") != "TESTED" or row.get("median_difference") is None:
        return (
            f"{target} versus {reference}: insufficient eligible clusters for inference.",
            False,
        )
    difference = float(row["median_difference"])
    direction = (
        "more dispersed"
        if difference > 0
        else "more compact"
        if difference < 0
        else "equal"
    )
    q_value = row.get("fdr_q_value")
    significant = q_value is not None and float(q_value) <= 0.05
    interval = ""
    if row.get("median_difference_ci_low") is not None and row.get(
        "median_difference_ci_high"
    ) is not None:
        interval = (
            f"; 95% CI {_format_number(value=row['median_difference_ci_low'])} to "
            f"{_format_number(value=row['median_difference_ci_high'])}"
        )
    statistics = (
        f"median difference {_format_number(value=difference)}{interval}; "
        f"Cliff's delta {_format_number(value=row.get('cliffs_delta'))}; "
        f"BH-FDR q {_format_number(value=q_value)}"
    )
    if significant:
        return (
            f"FDR-significant: {target} were {direction} than {reference} for "
            f"{metric_label.lower()} ({statistics}).",
            True,
        )
    return (
        f"No FDR-significant difference was detected between {target} and {reference} "
        f"for {metric_label.lower()}. The observed target profile was {direction} "
        f"({statistics}).",
        False,
    )


def _profile_group_count(
    *, profiles: Sequence[Mapping[str, Any]], profile_id: str
) -> int:
    """Return a pooled profile's distinct cluster count or zero."""

    return sum(
        int(row["group_count"])
        for row in profiles
        if str(row["profile_id"]) == profile_id
    )


def _metric_from_label(*, label: str) -> str:
    """Return the exact metric key for one user-facing label."""

    for metric in BENCHMARK_METRICS:
        if _METRIC_LABELS[metric] == label:
            return metric
    raise ValueError(f"Unknown benchmark metric label: {label}")


def _display_profile(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one profile row with plain-language headings."""

    profile_id = str(row["profile_id"])
    return {
        "Biological profile": benchmark_profile_label(profile_id=profile_id),
        "Profile ID": profile_id,
        "Class": benchmark_class_label(profile_classes=str(row["profile_class"])),
        "Subtype": str(row["profile_subclass"]).replace("_", " "),
        "Role": str(row["membership_role"]).replace("_", " ").title(),
        "Clusters": row["group_count"],
    }


def _display_marker(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one matched marker row with explicit group and provenance fields."""

    profile_id = str(row["profile_id"])
    return {
        "Biological profile": benchmark_profile_label(profile_id=profile_id),
        "Profile ID": profile_id,
        "Class": benchmark_class_label(profile_classes=str(row["profile_class"])),
        "Subtype": str(row["profile_subclass"]).replace("_", " "),
        "Marker gene / locus": row["marker_id"],
        "Protein accession": row["protein_identifier"],
        "Protein entry": row["protein_entry"],
        "Protein name": row["marker_name"],
        "OrthoFinder group": row["group_id"],
        "Matched protein": row["matched_member_id"],
        "Species": row["matched_species_label"],
        "Domain architecture": row["domain_architecture"],
        "Evidence": row["evidence_type"],
        "Source": row["source_title"],
        "DOI": row["source_doi"],
        "Source version": row["source_version"],
        "Authority": row["authority_name"],
    }


def _display_member(
    *, row: Mapping[str, Any], marker_members: set[str]
) -> dict[str, Any]:
    """Return one cluster member and flag profile-defining matched proteins."""

    member_id = str(row["member_id"])
    return {
        "Species": row["species_label"],
        "Protein ID": member_id,
        "Defines selected profile": member_id in marker_members,
        "Parent legacy orthogroup": row["legacy_orthogroup_id"],
        "Gene-tree parent clade": row["gene_tree_parent_clade"],
    }


def _display_classification(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one matched-control classification with readable headings."""

    return {
        "Group": row["group_id"],
        "Biological class": benchmark_class_label(
            profile_classes=str(row.get("profile_classes", ""))
        ),
        "Profiles": row.get("profile_ids", ""),
        "Average pair distance": row.get("mean_distance"),
        "Average-distance percentile": row.get("mean_distance_control_percentile"),
        "Central divergence class": str(
            row.get("central_divergence_class", "")
        ).replace("_OR_", " / ").title(),
        "Distance spread (SD)": row.get("distance_sd"),
        "Spread percentile": row.get("distance_sd_control_percentile"),
        "Heterogeneity class": str(row.get("heterogeneity_class", "")).replace(
            "_OR_", " / "
        ).title(),
        "Matched controls": row.get("matched_control_count"),
        "Status": row.get("status"),
    }


def _display_contrast(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one contrast row with readable names and exact identifiers."""

    target_id = str(row["target_profile_id"])
    reference_id = str(row["reference_profile_id"])
    return {
        "Target profile": benchmark_profile_label(profile_id=target_id),
        "Target profile ID": target_id,
        "Reference profile": benchmark_profile_label(profile_id=reference_id),
        "Reference profile ID": reference_id,
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

    background_id = str(row["background_profile_id"])
    return {
        "Background": benchmark_profile_label(profile_id=background_id),
        "Background ID": background_id,
        "Scale": str(row["comparison_scale"]).replace("_", " ").title(),
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

    profiles = ", ".join(
        benchmark_profile_label(profile_id=value)
        for value in str(row["profile_ids"]).split(";")
        if value
    )
    return (
        f"{row['group_type']} | {row['hierarchy_node'] or 'ROOT'} | "
        f"{row['group_id']} — {profiles}"
    )


def _marker_match_label(*, row: Mapping[str, Any]) -> str:
    """Return a unique marker, group and matched-protein selector label."""

    name = str(row.get("marker_name", "")).strip()
    description = f" — {name[:80]}" if name else ""
    return (
        f"{row['marker_id']} | {row['group_id']} | {row['matched_member_id']}"
        f"{description}"
    )


def _group_key(*, row: Mapping[str, Any], run_id: str) -> GroupKey:
    """Return a collision-safe group key from a benchmark result row."""

    return GroupKey(
        run_id=run_id,
        group_type=str(row["group_type"]),
        hierarchy_node=str(row["hierarchy_node"]),
        group_id=str(row["group_id"]),
    )


def _open_group_in_explorer(*, key: GroupKey) -> None:
    """Store a benchmark cluster and navigate from a Streamlit callback."""

    _store_active_group(key=key)
    st.session_state["app_page"] = "Cluster explorer"


def _column_config(*, descriptions: Mapping[str, str]) -> dict[str, Any]:
    """Return hoverable Streamlit column descriptions."""

    return {
        label: st.column_config.Column(label=label, help=description)
        for label, description in descriptions.items()
    }


def _safe_file_token(*, value: str) -> str:
    """Return a filesystem-safe deterministic export filename token."""

    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    return token[:120] or "benchmark_profile"


def _format_number(*, value: object) -> str:
    """Format an optional numeric value without converting absence to zero."""

    if value is None or value == "":
        return "Unavailable"
    return f"{float(value):.4g}"

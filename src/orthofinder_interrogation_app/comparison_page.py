"""Multi-cluster distance and dispersion comparison workspace."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import streamlit as st

from orthofinder_results.errors import OrthoFinderResultsError

from .dispersion import PcoaGeometry, classical_pcoa
from .distance_data import DistanceAnalysisProvider, GroupAnalysis
from .evolutionary_page import (
    _comparison_keys,
    _load_catalog,
    _store_comparison_keys,
)
from .exports import render_table_downloads
from .figures import (
    comparison_distribution_figure,
    comparison_pcoa_figure,
    comparison_summary_figure,
)
from .guidance import render_graph_guidance
from .queries import OrthoFinderQueryService

_LOGGER = logging.getLogger("orthofinder_interrogation_app.comparison_page")
_COMPARISON_COLUMN_HELP = {
    "Cluster": "Composite group identity: group system, species-tree level and group ID.",
    "Proteins in full group": "Complete group membership before bounded sampling.",
    "Proteins analysed": "Proteins represented in this exact comparison matrix.",
    "Protein pairs": "Number of exact pair distances compared.",
    "Smallest distance": "Minimum exact pair distance in the displayed matrix.",
    "Lower quartile": "25th percentile of exact pair distances.",
    "Median distance": "Middle exact pair distance.",
    "Average distance": "Mean exact pair distance; smaller is generally more compact.",
    "Upper quartile": "75th percentile of exact pair distances.",
    "Largest distance": "Maximum exact pair distance.",
    "Distance spread (SD)": "Population standard deviation of exact pair distances.",
    "Distance method": "Calculation method used for this matrix.",
    "Calculation scope": "Whether the matrix is complete or a deterministic sample.",
    "Analysis source": "Stored resource, embedded pilot report or portable-tree sidecar.",
}


def render_cluster_comparison(
    *, resource: Any, service: OrthoFinderQueryService, cache_dir: Path
) -> None:
    """Render a two-to-twelve-group comparison using exact displayed distances.

    Args:
        resource: Validated immutable resource identity.
        service: Read-only DuckDB query service.
        cache_dir: Persistent sidecar cache outside the completed resource.
    """

    st.header("Compare clusters")
    st.caption(
        "Compare compactness and distance spread across 2–12 groups. Start with the summary "
        "and distributions; use the PCoA panels to compare within-group shape."
    )
    with st.expander("How to make a defensible comparison"):
        st.markdown(
            """
            - Check **Proteins analysed** and **Calculation scope** before comparing groups of
              very different sizes.
            - Compare the mean or median together with **Distance spread (SD)** and the full
              violin/CDF shape; no single statistic captures multi-scale dispersion.
            - Each PCoA panel has its own coordinate system. Compare internal shape and fit,
              never absolute position, angle or orientation between panels.
            - Method, sampling and analysis-source fields remain in the authority table and
              paired TSV/Excel downloads.
            """
        )
    keys = tuple(key for key in _comparison_keys() if key.run_id == resource.run_id)
    if not keys:
        st.info(
            "The comparison workspace is empty. Find a group first, then choose “Add selected "
            "group to comparison” in Find groups, Taxonomic search or Explore one cluster."
        )
        return
    labels = tuple(key.display_label() for key in keys)
    selected_labels = tuple(
        st.multiselect(
            "Groups retained in this comparison",
            labels,
            default=labels,
            help="Deselect a group to remove it from the persistent workspace.",
        )
    )
    selected = tuple(key for key in keys if key.display_label() in selected_labels)
    action_columns = st.columns(2)
    if action_columns[0].button("Apply group removals"):
        _store_comparison_keys(keys=selected)
        keys = selected
    if action_columns[1].button("Clear comparison workspace"):
        _store_comparison_keys(keys=())
        st.info("Comparison workspace cleared.")
        return
    else:
        keys = selected
    if len(keys) < 2:
        st.warning("Retain at least two groups for comparison.")
        return
    if len(keys) > 12:
        st.warning("Comparison is limited to twelve groups.")
        return
    controls = st.columns(2)
    max_members = int(
        controls[0].select_slider(
            "Maximum proteins for each new calculation",
            options=(50, 100, 250, 500),
            value=250,
            help=(
                "Controls on-demand tree calculations. Persisted matrices retain their "
                "exact published sample."
            ),
        )
    )
    nearest_neighbours = int(
        controls[1].slider(
            "Neighbours retained per protein",
            1,
            10,
            3,
            help="Stored in the reusable analysis payload for each selected cluster.",
        )
    )
    provider = DistanceAnalysisProvider(
        service=service,
        cache_dir=cache_dir,
        report_catalog=_load_catalog(resource=resource),
    )
    analyses: list[GroupAnalysis] = []
    failures: list[tuple[str, str]] = []
    with st.spinner("Loading exact distance analyses for the comparison…"):
        for key in keys:
            try:
                analyses.append(
                    provider.analyse(
                        key=key,
                        max_members=max_members,
                        nearest_neighbours=nearest_neighbours,
                    )
                )
            except OrthoFinderResultsError as error:
                _LOGGER.warning("Comparison group unavailable: %s: %s", key.display_label(), error)
                failures.append((key.display_label(), str(error)))
    for label, reason in failures:
        st.warning(f"{label}: {reason}")
    if len(analyses) < 2:
        st.error("Fewer than two selected groups have usable exact displayed distances.")
        return
    summaries = _comparison_summaries(analyses=tuple(analyses))
    render_graph_guidance(key="comparison_summary")
    st.plotly_chart(
        comparison_summary_figure(summaries=summaries),
        width="stretch",
        config={"displaylogo": False},
    )
    distributions = {
        analysis.key.display_label(): tuple(
            float(row["distance"]) for row in analysis.distances
        )
        for analysis in analyses
    }
    st.subheader("Compare complete displayed distance distributions")
    st.write(
        "Violin plots reveal distribution shape and central range. Empirical CDFs show the "
        "fraction of protein pairs at or below each distance."
    )
    distribution_tabs = st.tabs(("Violin distributions", "Cumulative distributions"))
    with distribution_tabs[0]:
        render_graph_guidance(key="comparison_violin")
        st.plotly_chart(
            comparison_distribution_figure(
                distance_groups=distributions,
                mode="VIOLIN",
            ),
            width="stretch",
            config={"displaylogo": False},
        )
    with distribution_tabs[1]:
        render_graph_guidance(key="comparison_ecdf")
        st.plotly_chart(
            comparison_distribution_figure(
                distance_groups=distributions,
                mode="ECDF",
            ),
            width="stretch",
            config={"displaylogo": False},
        )
    geometries: dict[str, PcoaGeometry] = {}
    for analysis in analyses:
        try:
            geometries[analysis.key.display_label()] = classical_pcoa(
                rows=analysis.distances,
                members=analysis.members,
            )
        except OrthoFinderResultsError as error:
            st.warning(f"PCoA unavailable for {analysis.key.display_label()}: {error}")
    if len(geometries) >= 2:
        render_graph_guidance(key="comparison_pcoa")
        st.plotly_chart(
            comparison_pcoa_figure(geometries=geometries),
            width="stretch",
            config={"displaylogo": False},
        )
        st.caption(
            "Panels use consistent species colours but independent coordinate systems. "
            "Arms or gaps are diagnostic patterns, not automatic subfamilies."
        )
    st.subheader("Comparison data and provenance")
    st.dataframe(
        tuple(_display_comparison_summary_row(row=row) for row in summaries),
        width="stretch",
        hide_index=True,
        column_config=_column_config(descriptions=_COMPARISON_COLUMN_HELP),
    )
    displayed_summaries = tuple(
        _display_comparison_summary_row(row=row) for row in summaries
    )
    render_table_downloads(
        records=displayed_summaries,
        file_stem="orthofinder_cluster_comparison",
        key="cluster_comparison_download",
        tsv_label="Download comparison summary as TSV",
        excel_label="Download comparison summary as formatted Excel",
        column_definitions=_COMPARISON_COLUMN_HELP,
        workbook_title="OrthoFinder cluster comparison",
    )


def _comparison_summaries(
    *, analyses: tuple[GroupAnalysis, ...]
) -> tuple[dict[str, Any], ...]:
    """Return explicit method- and sampling-labelled comparison rows.

    Args:
        analyses: Two or more exact displayed group analyses.

    Returns:
        Deterministically labelled comparison summary records.
    """

    rows = []
    for analysis in analyses:
        summary = analysis.summary
        rows.append(
            {
                "label": analysis.key.display_label(),
                "group_type": analysis.key.group_type,
                "hierarchy_node": analysis.key.hierarchy_node,
                "group_id": analysis.key.group_id,
                "analytical_members": int(summary["total_member_count"]),
                "displayed_members": int(summary["sampled_member_count"]),
                "pair_count": int(summary["distance_pair_count"]),
                "minimum_distance": float(summary["minimum_distance"]),
                "q25_distance": float(summary["q25_distance"]),
                "median_distance": float(summary["median_distance"]),
                "mean_distance": float(summary["mean_distance"]),
                "q75_distance": float(summary["q75_distance"]),
                "maximum_distance": float(summary["maximum_distance"]),
                "population_stddev_distance": float(
                    summary["population_stddev_distance"]
                ),
                "distance_method": str(summary["distance_method"]),
                "computation_status": str(summary["computation_status"]),
                "analysis_source": analysis.source,
            }
        )
    rows.sort(key=lambda row: str(row["label"]))
    return tuple(rows)


def _display_comparison_summary_row(*, row: dict[str, Any]) -> dict[str, Any]:
    """Return readable comparison fields for display and paired downloads."""

    return {
        "Cluster": row["label"],
        "Proteins in full group": row["analytical_members"],
        "Proteins analysed": row["displayed_members"],
        "Protein pairs": row["pair_count"],
        "Smallest distance": row["minimum_distance"],
        "Lower quartile": row["q25_distance"],
        "Median distance": row["median_distance"],
        "Average distance": row["mean_distance"],
        "Upper quartile": row["q75_distance"],
        "Largest distance": row["maximum_distance"],
        "Distance spread (SD)": row["population_stddev_distance"],
        "Distance method": row["distance_method"],
        "Calculation scope": row["computation_status"],
        "Analysis source": row["analysis_source"],
    }


def _column_config(*, descriptions: dict[str, str]) -> dict[str, Any]:
    """Return Streamlit column definitions with hoverable help descriptions."""

    return {
        label: st.column_config.Column(label=label, help=description)
        for label, description in descriptions.items()
    }

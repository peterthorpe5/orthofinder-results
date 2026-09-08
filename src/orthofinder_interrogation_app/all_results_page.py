"""Complete, selectable exports of persisted group-distance summaries."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import Any

import streamlit as st

from orthofinder_results.errors import InputValidationError

from .exports import render_table_downloads
from .models import DistanceResultFilters
from .queries import MAX_ALL_DISTANCE_RESULTS, OrthoFinderQueryService

_LOGGER = logging.getLogger("orthofinder_interrogation_app.all_results_page")
_RESULT_COLUMNS = {
    "Run": (
        "run_id",
        "Immutable identifier of the completed OrthoFinder resource.",
    ),
    "Group system": (
        "group_type",
        "HOG is hierarchical; LEGACY_ORTHOGROUP is the flat Orthogroups.tsv set.",
    ),
    "Species-tree level": (
        "hierarchy_node",
        "Species-tree node at which HOG membership is defined; ROOT is the flat collection.",
    ),
    "Group ID": (
        "group_id",
        "Exact OrthoFinder group identifier within the stated system and level.",
    ),
    "Parent legacy orthogroup": (
        "legacy_orthogroup_id",
        "Linked flat orthogroup identifier when OrthoFinder supplied one.",
    ),
    "Gene-tree parent clade": (
        "gene_tree_parent_clade",
        "Gene-tree clade used to define this HOG when available.",
    ),
    "Proteins in full group": (
        "member_count",
        "Complete protein membership before bounded distance sampling.",
    ),
    "Species represented": (
        "species_count",
        "Number of sampled species contributing at least one protein.",
    ),
    "Single-copy species": (
        "single_copy_species_count",
        "Represented species contributing exactly one protein.",
    ),
    "Highest copies in one species": (
        "max_copies_per_species",
        "Largest protein copy count observed for one represented species.",
    ),
    "Average copies per represented species": (
        "mean_copies_per_species",
        "Full-group protein count divided by represented species.",
    ),
    "Distance method": (
        "distance_method",
        "Exact calculation method used for the selected distance matrix.",
    ),
    "Calculation scope": (
        "computation_status",
        "Whether the stored matrix is complete or a deterministic bounded sample.",
    ),
    "Protein identifier resolution": (
        "member_identifier_resolution",
        "Rule used to reconcile group proteins with the distance authority.",
    ),
    "Distance authority proteins": (
        "total_member_count",
        "Full protein count declared by the distance summary.",
    ),
    "Proteins analysed": (
        "sampled_member_count",
        "Proteins represented in the stored pair-distance matrix.",
    ),
    "Sampling fraction": (
        "sampling_fraction",
        "Proteins analysed divided by proteins in the distance authority.",
    ),
    "Full-group matrix": (
        "full_group_matrix",
        "True when every protein in the distance authority was analysed.",
    ),
    "Protein pairs": (
        "distance_pair_count",
        "Number of stored protein-to-protein distance values.",
    ),
    "Unresolved pairs": (
        "unresolved_pair_count",
        "Candidate pairs that could not produce a valid distance.",
    ),
    "Smallest distance": (
        "minimum_distance",
        "Minimum stored pairwise distance.",
    ),
    "5th percentile distance": (
        "q05_distance",
        "Distance below which five per cent of stored protein pairs fall.",
    ),
    "Lower quartile distance": (
        "q25_distance",
        "Distance below which 25 per cent of stored protein pairs fall.",
    ),
    "Median distance": (
        "median_distance",
        "Middle stored pairwise distance.",
    ),
    "Average distance": (
        "mean_distance",
        "Mean stored pairwise distance; smaller generally means a more compact group.",
    ),
    "Upper quartile distance": (
        "q75_distance",
        "Distance below which 75 per cent of stored protein pairs fall.",
    ),
    "95th percentile distance": (
        "q95_distance",
        "Distance below which 95 per cent of stored protein pairs fall.",
    ),
    "Largest distance": (
        "maximum_distance",
        "Maximum stored pairwise distance.",
    ),
    "Distance spread (SD)": (
        "population_stddev_distance",
        "Population standard deviation of stored pair distances.",
    ),
    "Interquartile distance range": (
        "interquartile_range",
        "Upper quartile minus lower quartile; robust central distance spread.",
    ),
    "Relative distance spread": (
        "relative_distance_spread",
        "Population SD divided by mean distance; unavailable when the mean is zero.",
    ),
    "Average comparable sites": (
        "mean_comparable_sites",
        "Mean aligned sites used per pair when the distance method defines this quantity.",
    ),
    "Source file": (
        "source_file",
        "Resource-relative or provenance source recorded for the distance calculation.",
    ),
    "Failure reason": (
        "failure_reason",
        "Recorded calculation warning or failure detail; normally empty for successful rows.",
    ),
}
_DEFAULT_COLUMNS = (
    "Group system",
    "Species-tree level",
    "Group ID",
    "Proteins in full group",
    "Species represented",
    "Proteins analysed",
    "Sampling fraction",
    "Full-group matrix",
    "Protein pairs",
    "Smallest distance",
    "Lower quartile distance",
    "Median distance",
    "Average distance",
    "Upper quartile distance",
    "Largest distance",
    "Distance spread (SD)",
    "Interquartile distance range",
    "Relative distance spread",
    "Distance method",
    "Calculation scope",
)


def render_all_distance_results(*, service: OrthoFinderQueryService) -> None:
    """Render a selectable, complete persisted-distance results export.

    Args:
        service: Read-only query service for one validated resource.
    """

    st.header("All distance results")
    st.write(
        "Build one dataset-wide table containing every cluster with a successful persisted "
        "distance summary. Choose the biological and provenance columns needed for analysis, "
        "then download the complete selection as TSV or formatted Excel."
    )
    with st.expander("What is included—and what is not", expanded=False):
        st.markdown(
            """
            - The table contains **one preferred stored result per cluster**. A complete matrix is
              preferred to a sample; otherwise the largest stored sample is selected.
            - Mean, median, quantiles and SD describe only the proteins stated in **Proteins
              analysed**. Missing values are never treated as biological zero.
            - A method or scope filter is applied before choosing the preferred row.
            - On-demand sidecar caches are deliberately not mixed into this immutable,
              resource-level export. Publish them in a rebuilt resource before treating them as
              dataset-wide authority.
            """
        )
    facets = service.distance_result_facets()
    filter_columns = st.columns(4)
    group_type = filter_columns[0].selectbox(
        "Group system",
        ("All", *facets["group_types"]),
        help="Restrict the export to one group authority.",
    )
    hierarchy_options = tuple("ROOT" if not value else value for value in facets["hierarchy_nodes"])
    hierarchy = filter_columns[1].selectbox(
        "Species-tree level",
        ("All", *hierarchy_options),
        help="Restrict the export to one exact hierarchy node.",
    )
    distance_method = filter_columns[2].selectbox(
        "Distance method",
        ("All", *facets["distance_methods"]),
        help="Restrict the export to one exact stored calculation method.",
    )
    computation_status = filter_columns[3].selectbox(
        "Calculation scope",
        ("All", *facets["computation_statuses"]),
        help="Restrict the export to one exact complete or sampled status.",
    )
    selected_labels = tuple(
        st.multiselect(
            "Columns to include",
            options=tuple(_RESULT_COLUMNS),
            default=_DEFAULT_COLUMNS,
            help=(
                "The preview and both downloads use exactly this order. The Excel workbook "
                "adds a definitions sheet for the selected columns."
            ),
        )
    )
    if not selected_labels:
        st.warning("Select at least one output column.")
        return
    filters = DistanceResultFilters(
        group_type="" if group_type == "All" else group_type,
        hierarchy_node=(
            None if hierarchy == "All" else "" if hierarchy == "ROOT" else hierarchy
        ),
        distance_method="" if distance_method == "All" else distance_method,
        computation_status=(
            "" if computation_status == "All" else computation_status
        ),
    )
    raw_columns = tuple(_RESULT_COLUMNS[label][0] for label in selected_labels)
    try:
        result = service.distance_results(filters=filters, columns=raw_columns)
    except InputValidationError as error:
        _LOGGER.exception("Complete distance-result export failed")
        st.error(str(error))
        return
    if not result.rows:
        st.warning("No persisted distance results match the selected filters.")
        return
    displayed = _display_result_rows(rows=result.rows, selected_labels=selected_labels)
    metric_columns = st.columns(3)
    metric_columns[0].metric(
        "Clusters in complete export",
        f"{result.total_rows:,}",
        help="One preferred persisted distance summary per matching cluster.",
    )
    metric_columns[1].metric(
        "Selected columns",
        f"{len(selected_labels):,}",
        help="Columns retained in both downloads and the preview.",
    )
    metric_columns[2].metric(
        "Interactive export limit",
        f"{MAX_ALL_DISTANCE_RESULTS:,}",
        help="Narrow the filters if a future resource exceeds this defensive memory bound.",
    )
    preview = displayed[:500]
    st.subheader("Preview")
    st.caption(
        f"Showing {len(preview):,} of {len(displayed):,} rows. Both downloads contain all "
        "matching rows, not only this browser preview."
    )
    definitions = {label: _RESULT_COLUMNS[label][1] for label in selected_labels}
    st.dataframe(
        preview,
        width="stretch",
        hide_index=True,
        column_config={
            label: st.column_config.Column(label=label, help=definitions[label])
            for label in selected_labels
        },
    )
    render_table_downloads(
        records=displayed,
        file_stem=f"{service.resource.run_id}_all_distance_results",
        key="all_distance_results_download",
        tsv_label="Download all selected results as TSV",
        excel_label="Download all selected results as formatted Excel",
        column_definitions=definitions,
        workbook_title=f"OrthoFinder distance results: {service.resource.run_id}",
    )


def _display_result_rows(
    *,
    rows: tuple[dict[str, Any], ...],
    selected_labels: tuple[str, ...],
) -> tuple[dict[str, Any], ...]:
    """Map selected raw result fields to plain-language headings."""

    unsupported = tuple(label for label in selected_labels if label not in _RESULT_COLUMNS)
    if unsupported:
        raise InputValidationError(
            "Unsupported displayed result columns: " + "; ".join(unsupported)
        )
    displayed = []
    for row in rows:
        record: dict[str, Any] = {}
        for label in selected_labels:
            raw_name = _RESULT_COLUMNS[label][0]
            value = row.get(raw_name)
            if label == "Species-tree level" and value == "":
                value = "ROOT"
            record[label] = value
        displayed.append(record)
    return tuple(displayed)


def all_result_column_definitions() -> Mapping[str, str]:
    """Return the complete immutable plain-language column dictionary."""

    return {label: definition for label, (_, definition) in _RESULT_COLUMNS.items()}

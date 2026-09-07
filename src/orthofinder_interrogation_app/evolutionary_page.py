"""Streamlit page for bounded, method-labelled evolutionary views."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import streamlit as st

from orthofinder_results.errors import InputValidationError

from .figures import (
    distance_matrix_figure,
    linked_member_ids,
    nearest_neighbour_figure,
    pcoa_figure,
    phylogram_figure,
    shepard_figure,
)
from .report_data import VisualisationCatalog, load_visualisation_catalog
from .tsv import records_to_tsv

_LOGGER = logging.getLogger("orthofinder_interrogation_app.evolutionary_page")


def render_evolutionary_views(*, resource: Any) -> None:
    """Render all bounded evolutionary panels for one calculated group."""

    st.header("Evolutionary views")
    st.caption(
        "These views cover only groups whose bounded distance pilot was embedded in this "
        "completed resource. Group searching remains complete in DuckDB."
    )
    if resource.report_path is None:
        st.warning(
            "This schema-2 resource has no offline report, so its bounded phylogram and "
            "projection payload is unavailable."
        )
        return
    try:
        report_path = Path(resource.report_path)
        statistic = report_path.stat()
        catalog = _cached_catalog(
            report_path=str(report_path),
            expected_run_id=str(resource.run_id),
            size=statistic.st_size,
            modified_ns=statistic.st_mtime_ns,
        )
    except (OSError, InputValidationError) as error:
        _LOGGER.exception("Could not load bounded evolutionary views")
        st.error(f"Evolutionary views could not be loaded: {error}")
        return
    if not catalog.networks:
        st.info("No groups with bounded evolutionary views are present in this report.")
        return
    selected_label = st.selectbox(
        "Group with calculated distances",
        catalog.labels(),
        help="Only the explicitly calculated and embedded distance groups appear here.",
    )
    entry = catalog.get(label=selected_label)
    _render_distance_scope(entry=entry)
    members = _members(entry=entry)
    member_species = {str(row["member_id"]): str(row["species_label"]) for row in members}
    species = tuple(sorted(set(member_species.values())))
    controls = st.columns(2)
    selected_species = tuple(
        controls[0].multiselect(
            "Linked species selection",
            species,
            help="Highlights every displayed member of each selected species in all views.",
        )
    )
    selected_members = tuple(
        controls[1].multiselect(
            "Linked member selection",
            tuple(sorted(member_species)),
            help="Adds exact displayed proteins to the selection in all views.",
        )
    )
    linked = linked_member_ids(
        entry=entry,
        selected_members=selected_members,
        selected_species=selected_species,
    )
    if linked:
        st.info(
            f"Linked selection: {len(linked):,} of {len(members):,} displayed members "
            f"across {len({member_species[value] for value in linked}):,} species."
        )
    tabs = st.tabs(
        (
            "PCoA + diagnostics",
            "Shepard plot",
            "Branch-length phylogram",
            "Exact distance matrix",
            "Nearest-neighbour topology",
            "Selected members",
        )
    )
    with tabs[0]:
        _render_pcoa(entry=entry, linked=linked)
    with tabs[1]:
        _render_shepard(entry=entry)
    with tabs[2]:
        _render_phylogram(entry=entry, linked=linked)
    with tabs[3]:
        _render_matrix(entry=entry, linked=linked)
    with tabs[4]:
        _render_topology(
            entry=entry,
            linked=linked,
            nearest_neighbours=catalog.nearest_neighbours,
        )
    with tabs[5]:
        _render_selected_members(members=members, linked=linked)


@st.cache_data(show_spinner="Loading bounded evolutionary records…")
def _cached_catalog(
    *, report_path: str, expected_run_id: str, size: int, modified_ns: int
) -> VisualisationCatalog:
    """Cache a report catalog while including file identity in the cache key."""

    del size, modified_ns
    return load_visualisation_catalog(
        report_path=Path(report_path), expected_run_id=expected_run_id
    )


def _render_distance_scope(*, entry: dict[str, Any]) -> None:
    """Render sample size, distance method and medoid scope before the charts."""

    summary = _required_mapping(entry=entry, key="distanceSummary")
    columns = st.columns(5)
    columns[0].metric("Analytical members", f"{int(summary['total_member_count']):,}")
    columns[1].metric("Displayed sample", f"{int(summary['sampled_member_count']):,}")
    columns[2].metric("Exact pairs", f"{int(summary['distance_pair_count']):,}")
    columns[3].metric("Mean distance", _format_optional(summary.get("mean_distance")))
    columns[4].metric("Distance SD", _format_optional(summary.get("population_stddev_distance")))
    st.caption(
        f"Distance method: {summary.get('distance_method', 'unavailable')}; calculation "
        f"status: {summary.get('computation_status', 'unavailable')}. Every view below "
        "describes this displayed sample, not uncalculated full-group pair distances."
    )
    medoid = _required_mapping(entry=entry, key="medoid")
    if medoid.get("member_id"):
        st.caption(
            f"★ Sample medoid: {medoid['member_id']} (mean distance "
            f"{_format_optional(medoid.get('mean_distance'))} within the displayed "
            "sample). It is not a reconstructed ancestor or guaranteed full-group centre."
        )


def _render_pcoa(*, entry: dict[str, Any], linked: frozenset[str]) -> None:
    """Render PCoA coordinates and conservative quality diagnostics."""

    projection = _required_mapping(entry=entry, key="distanceProjection")
    if projection.get("status") != "COMPLETE_DISTANCE_PCOA_2D":
        st.warning(f"PCoA unavailable: {projection.get('reason', projection.get('status'))}")
        return
    diagnostics = st.columns(5)
    diagnostics[0].metric("Fit category", str(projection.get("quality_category", "")))
    diagnostics[1].metric(
        "Axes 1+2 positive inertia",
        _format_percent(projection.get("two_axis_positive_inertia_fraction")),
    )
    diagnostics[2].metric(
        "Distance correlation", _format_optional(projection.get("distance_correlation"))
    )
    diagnostics[3].metric(
        "Normalised stress", _format_optional(projection.get("normalised_stress"))
    )
    diagnostics[4].metric(
        "Negative inertia", _format_percent(projection.get("negative_inertia_fraction"))
    )
    category = str(projection.get("quality_category", ""))
    explanation = str(projection.get("quality_explanation", ""))
    if category == "POOR":
        st.warning(explanation)
    elif category == "MODERATE":
        st.info(explanation)
    else:
        st.success(explanation)
    st.plotly_chart(
        pcoa_figure(entry=entry, selected_members=linked),
        use_container_width=True,
        config={"displaylogo": False},
    )
    st.caption(
        "PCoA is a diagnostic approximation of exact pair distances. Apparent arms, gaps "
        "or angles are not subfamilies and require confirmation in the phylogram and exact "
        "distance matrix."
    )


def _render_shepard(*, entry: dict[str, Any]) -> None:
    """Render the stratified Shepard diagnostic and its complete-pair scope."""

    projection = _required_mapping(entry=entry, key="distanceProjection")
    try:
        figure = shepard_figure(entry=entry)
    except InputValidationError as error:
        st.warning(str(error))
        return
    st.plotly_chart(figure, use_container_width=True, config={"displaylogo": False})
    st.caption(
        f"Deterministic {int(projection.get('shepard_point_count', 0)):,}-point summary "
        f"stratified across {int(projection.get('shepard_total_pair_count', 0)):,} exact "
        "input pairs. The diagonal denotes perfect two-dimensional preservation."
    )


def _render_phylogram(*, entry: dict[str, Any], linked: frozenset[str]) -> None:
    """Render the pruned branch-length resolved-gene-tree view."""

    tree = _required_mapping(entry=entry, key="phylogram")
    show_labels = st.checkbox("Show every displayed leaf label", value=False)
    try:
        figure = phylogram_figure(
            entry=entry,
            selected_members=linked,
            show_all_labels=show_labels,
        )
    except InputValidationError as error:
        st.warning(str(error))
        return
    st.plotly_chart(figure, use_container_width=True, config={"displaylogo": False})
    st.caption(
        f"{tree.get('status', 'Unavailable')} · resolved tree "
        f"{tree.get('treeId', 'unavailable')} · {int(tree.get('displayedLeafCount', 0)):,} "
        f"of {int(tree.get('requestedMemberCount', 0)):,} requested displayed leaves. "
        "Horizontal lengths are quantitative branch lengths; vertical spacing is not."
    )
    unresolved = tree.get("unresolvedMembers")
    if isinstance(unresolved, list) and unresolved:
        st.warning(f"Unresolved displayed tree members: {len(unresolved):,}.")


def _render_matrix(*, entry: dict[str, Any], linked: frozenset[str]) -> None:
    """Render the exact matrix with optional linked-selection subsetting."""

    restrict = st.checkbox(
        "Restrict the matrix to the linked selection",
        value=bool(len(linked) >= 2),
        disabled=len(linked) < 2,
        help="Select at least two displayed members directly or through species.",
    )
    selected = tuple(sorted(linked)) if restrict else ()
    try:
        figure = distance_matrix_figure(entry=entry, selected_members=selected)
    except InputValidationError as error:
        st.warning(str(error))
        return
    st.plotly_chart(figure, use_container_width=True, config={"displaylogo": False})
    matrix = _required_mapping(entry=entry, key="distanceMatrix")
    st.caption(
        f"{matrix.get('status', 'Unavailable')} · order: "
        f"{matrix.get('orderMethod', 'unavailable')} · every displayed cell retains one "
        "exact supplied pairwise distance."
    )


def _render_topology(
    *,
    entry: dict[str, Any],
    linked: frozenset[str],
    nearest_neighbours: int,
) -> None:
    """Render a separate force-directed nearest-neighbour topology."""

    metrics = _required_mapping(entry=entry, key="networkMetrics")
    include_connectors = st.checkbox(
        "Show layout-only component connectors",
        value=False,
        help="Dashed connectors are not nearest-neighbour relationships.",
    )
    try:
        figure = nearest_neighbour_figure(
            entry=entry,
            selected_members=linked,
            include_connectors=include_connectors,
        )
    except InputValidationError as error:
        st.warning(str(error))
        return
    st.plotly_chart(figure, use_container_width=True, config={"displaylogo": False})
    st.caption(
        f"Solid edges retain up to {nearest_neighbours:,} nearest neighbours per displayed "
        f"protein. Raw graph: {int(metrics.get('rawComponentCount', 0)):,} components, "
        f"{int(metrics.get('rawIsolateCount', 0)):,} isolates and "
        f"{int(metrics.get('nearestNeighbourEdgeCount', 0)):,} neighbour edges. The "
        "force-layout spacing is non-quantitative; components are not OrthoFinder splits."
    )


def _render_selected_members(
    *, members: tuple[dict[str, str], ...], linked: frozenset[str]
) -> None:
    """Render and export the linked displayed-member set."""

    visible = tuple(row for row in members if not linked or str(row["member_id"]) in linked)
    st.dataframe(visible, use_container_width=True, hide_index=True)
    st.download_button(
        "Download displayed selection as TSV",
        data=records_to_tsv(records=visible),
        file_name="orthofinder_displayed_member_selection.tsv",
        mime="text/tab-separated-values",
    )


def _members(*, entry: dict[str, Any]) -> tuple[dict[str, str], ...]:
    """Return validated concise displayed memberships."""

    raw_members = entry.get("members")
    if not isinstance(raw_members, list):
        raise InputValidationError("Visual entry lacks displayed members.")
    rows = []
    for raw_row in raw_members:
        if not isinstance(raw_row, dict):
            raise InputValidationError("Visual entry contains a malformed member row.")
        member_id = str(raw_row.get("member_id", ""))
        species = str(raw_row.get("species_label", ""))
        if not member_id or not species:
            raise InputValidationError("Visual entry contains an unlabelled member row.")
        rows.append({"member_id": member_id, "species_label": species})
    return tuple(rows)


def _required_mapping(*, entry: dict[str, Any], key: str) -> dict[str, Any]:
    """Return a required nested visual record."""

    value = entry.get(key)
    if not isinstance(value, dict):
        raise InputValidationError(f"Visual entry lacks a valid {key} record.")
    return value


def _format_optional(value: object) -> str:
    """Format an optional finite numeric diagnostic."""

    if value is None or value == "":
        return "Unavailable"
    return f"{float(value):.4g}"


def _format_percent(value: object) -> str:
    """Format a fractional diagnostic as a percentage."""

    if value is None or value == "":
        return "Unavailable"
    return f"{100.0 * float(value):.1f}%"

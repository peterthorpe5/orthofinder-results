"""Streamlit cluster explorer with linked exact-distance views."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Mapping

import streamlit as st

from orthofinder_results.errors import InputValidationError, OrthoFinderResultsError

from .dispersion import classical_pcoa, member_dispersion_rows, species_dispersion_rows
from .distance_data import DistanceAnalysisProvider, GroupAnalysis
from .exports import render_table_downloads
from .figures import (
    distance_distribution_figure,
    distance_matrix_figure,
    linked_member_ids,
    medoid_distance_figure,
    member_dispersion_figure,
    nearest_neighbour_figure,
    pcoa_3d_figure,
    pcoa_axis_figure,
    pcoa_figure,
    phylogram_figure,
    shepard_figure,
    species_pair_heatmap_figure,
)
from .force_view import force_directed_html
from .guidance import render_graph_guidance
from .models import GroupKey
from .queries import OrthoFinderQueryService
from .report_data import VisualisationCatalog, load_visualisation_catalog

_LOGGER = logging.getLogger("orthofinder_interrogation_app.evolutionary_page")
ACTIVE_GROUP_STATE = "orthofinder_active_group"
COMPARISON_STATE = "orthofinder_comparison_groups"

_MEMBER_DISPERSION_HELP = {
    "Protein ID": "Exact identifier of one analysed protein.",
    "Species": "Exact species label for that protein.",
    "Average distance to all others": (
        "Mean distance from this protein to every other analysed protein."
    ),
    "Median distance to all others": (
        "Median distance from this protein to every other analysed protein."
    ),
    "Distance spread (SD)": (
        "Population standard deviation of this protein's distances to all others."
    ),
    "Nearest protein": "Analysed protein with the smallest pair distance from this protein.",
    "Nearest protein's species": "Species label of the nearest analysed protein.",
    "Nearest distance": "Exact distance to the nearest analysed protein.",
    "Sample medoid": (
        "Yes for the analysed protein with the smallest average distance to all others."
    ),
}
_SPECIES_DISPERSION_HELP = {
    "Species A": "First species in the deterministic species-pair label.",
    "Species B": "Second species in the deterministic species-pair label.",
    "Pair class": "Within species when both endpoints share a label; otherwise between species.",
    "Protein pairs": "Number of exact protein-pair distances in this species pairing.",
    "Average pair distance": "Mean exact distance for this species pairing.",
    "Median pair distance": "Median exact distance for this species pairing.",
    "Distance spread (SD)": "Population standard deviation for this species pairing.",
    "Smallest pair distance": "Minimum exact distance in this species pairing.",
    "Largest pair distance": "Maximum exact distance in this species pairing.",
}
_PAIR_DISTANCE_HELP = {
    "Protein A": "First endpoint protein identifier.",
    "Species A": "Species of the first endpoint.",
    "Protein B": "Second endpoint protein identifier.",
    "Species B": "Species of the second endpoint.",
    "Pair distance": "Exact supplied or calculated distance between the two proteins.",
    "Distance method": "Method used to calculate this pair distance.",
    "Calculation scope": "Whether the matrix is complete or a deterministic bounded sample.",
    "Comparable sites": "Aligned positions compared when sequence distances are used.",
    "Mismatch sites": "Differing aligned positions when sequence distances are used.",
}
_MEMBER_HELP = {
    "Protein ID": "Exact protein identifier represented in the displayed analysis.",
    "Species": "Exact species label assigned to this protein in the resource.",
}


def render_evolutionary_views(
    *,
    resource: Any,
    service: OrthoFinderQueryService | None = None,
    cache_dir: Path | None = None,
) -> None:
    """Render the linked cluster explorer or its legacy report fallback.

    Args:
        resource: Validated immutable resource identity.
        service: Optional read-only DuckDB service used by the full explorer.
        cache_dir: Optional persistent analysis sidecar outside the resource.
    """

    catalog = _load_catalog(resource=resource)
    if service is None or cache_dir is None:
        _render_report_only(resource=resource, catalog=catalog)
        return
    st.header("Explore one cluster")
    st.caption(
        "Inspect one group through linked views. Begin with Distance spread for compactness, "
        "then use the gene-tree phylogram and exact distance heatmap to confirm patterns."
    )
    key = _select_group(service=service, catalog=catalog)
    if key is None:
        st.info("Choose a precomputed group or enter one exact group identifier to begin.")
        return
    controls = st.columns((2, 2, 2, 3))
    max_members = int(
        controls[0].select_slider(
            "Maximum proteins for a new calculation",
            options=(50, 100, 250, 500),
            value=250,
            help=(
                "On-demand tree calculations above this limit use an order-independent "
                "deterministic sample. Existing persisted matrices retain their exact "
                "published sample."
            ),
        )
    )
    nearest_neighbours = int(
        controls[1].slider(
            "Neighbours retained per protein",
            min_value=1,
            max_value=10,
            value=3,
            help="Controls topology edges, not OrthoFinder group membership.",
        )
    )
    force_recompute = controls[2].checkbox(
        "Recalculate saved on-demand result",
        value=False,
        help=(
            "Rebuild only a schema-3 on-demand sidecar result. Stored schema-2 and DuckDB "
            "distance authorities are never recalculated or modified."
        ),
    )
    controls[3].caption(
        "New bounded results are saved in the configured user cache. The completed "
        "resource remains read-only."
    )
    provider = DistanceAnalysisProvider(
        service=service,
        cache_dir=cache_dir,
        report_catalog=catalog,
    )
    try:
        with st.spinner("Loading or calculating the bounded exact-distance analysis…"):
            analysis = provider.analyse(
                key=key,
                max_members=max_members,
                nearest_neighbours=nearest_neighbours,
                force_recompute=force_recompute,
            )
    except OrthoFinderResultsError as error:
        st.warning(str(error))
        return
    _comparison_control(key=key)
    _render_analysis(analysis=analysis, nearest_neighbours=nearest_neighbours)


@st.cache_data(show_spinner="Loading embedded pilot records…")
def _cached_catalog(
    *, report_path: str, expected_run_id: str, size: int, modified_ns: int
) -> VisualisationCatalog:
    """Cache a report catalogue while including immutable file identity in its key."""

    del size, modified_ns
    return load_visualisation_catalog(
        report_path=Path(report_path), expected_run_id=expected_run_id
    )


def _load_catalog(*, resource: Any) -> VisualisationCatalog | None:
    """Return a validated optional report catalogue with visible failures."""

    if resource.report_path is None:
        return None
    try:
        report_path = Path(resource.report_path)
        statistic = report_path.stat()
        return _cached_catalog(
            report_path=str(report_path),
            expected_run_id=str(resource.run_id),
            size=statistic.st_size,
            modified_ns=statistic.st_mtime_ns,
        )
    except (OSError, InputValidationError) as error:
        _LOGGER.exception("Could not load embedded evolutionary views")
        st.error(f"Embedded pilot views could not be loaded: {error}")
        return None


def _render_report_only(
    *, resource: Any, catalog: VisualisationCatalog | None
) -> None:
    """Retain the report-only API used by small integrations and fallback tests."""

    st.header("Evolutionary views")
    if resource.report_path is None:
        st.warning("This resource has no embedded offline visual report.")
        return
    if catalog is None:
        return
    if not catalog.networks:
        st.info("No groups with bounded evolutionary views are present in this report.")
        return
    selected_label = st.selectbox("Group with calculated distances", catalog.labels())
    entry = catalog.get(label=selected_label)
    _render_distance_scope(entry=entry)
    members = _members(entry=entry)
    linked = linked_member_ids(entry=entry)
    tabs = st.tabs(
        (
            "PCoA view & fit",
            "Distance-fit check",
            "Gene-tree phylogram",
            "Distance heatmap",
            "Nearest-neighbour topology",
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
    _render_selected_members(members=members, linked=linked)


def _select_group(
    *, service: OrthoFinderQueryService, catalog: VisualisationCatalog | None
) -> GroupKey | None:
    """Render compact precomputed/manual group controls and retain one active key."""

    stored = _state_group()
    precomputed = _catalog_keys(catalog=catalog, run_id=service.resource.run_id)
    with st.expander("Choose cluster", expanded=stored is None):
        if precomputed:
            labels = tuple(key.display_label() for key in precomputed)
            label_to_key = dict(zip(labels, precomputed, strict=True))
            default_index = 0
            if stored is not None and stored.display_label() in label_to_key:
                default_index = labels.index(stored.display_label())
            selected = st.selectbox(
                "Groups with stored pilot distances",
                labels,
                index=default_index,
                help="These open immediately from the existing exact report matrices.",
            )
            if st.button("Use selected stored group"):
                stored = label_to_key[selected]
                _store_active_group(key=stored)
        with st.form("exact_cluster_key"):
            columns = st.columns(3)
            group_types = service.list_group_types()
            initial_type = stored.group_type if stored is not None else group_types[0]
            group_type = columns[0].selectbox(
                "Group system",
                group_types,
                index=group_types.index(initial_type) if initial_type in group_types else 0,
                help=(
                    "HOG is a hierarchical orthogroup at one species-tree level; "
                    "LEGACY_ORTHOGROUP is the flat Orthogroups.tsv collection."
                ),
            )
            nodes = service.list_hierarchy_nodes(group_type=group_type)
            initial_node = stored.hierarchy_node if stored is not None else nodes[0]
            hierarchy_node = columns[1].selectbox(
                "Species-tree level",
                nodes,
                format_func=lambda value: value or "ROOT",
                index=nodes.index(initial_node) if initial_node in nodes else 0,
                help="Exact OrthoFinder node at which this HOG membership is defined.",
            )
            group_id = columns[2].text_input(
                "Exact group ID",
                value=stored.group_id if stored is not None else "",
                help="Enter the complete group identifier, for example N0.HOG0000001.",
            )
            submitted = st.form_submit_button("Open cluster", type="primary")
        if submitted:
            candidate = GroupKey(
                run_id=service.resource.run_id,
                group_type=group_type,
                hierarchy_node=hierarchy_node,
                group_id=group_id.strip(),
            )
            if not candidate.group_id:
                st.warning("Enter an exact group identifier.")
            else:
                service.get_group(key=candidate)
                _store_active_group(key=candidate)
                stored = candidate
    if stored is None and precomputed:
        stored = precomputed[0]
        _store_active_group(key=stored)
    if stored is not None:
        st.subheader(stored.display_label())
    return stored


def _catalog_keys(
    *, catalog: VisualisationCatalog | None, run_id: str
) -> tuple[GroupKey, ...]:
    """Convert valid report composite keys to group models."""

    if catalog is None:
        return ()
    keys = []
    for label in catalog.labels():
        parts = label.split("|", maxsplit=2)
        if len(parts) != 3 or not parts[0] or not parts[2]:
            continue
        keys.append(
            GroupKey(
                run_id=run_id,
                group_type=parts[0],
                hierarchy_node=parts[1],
                group_id=parts[2],
            )
        )
    return tuple(keys)


def _render_analysis(*, analysis: GroupAnalysis, nearest_neighbours: int) -> None:
    """Render every linked cluster view from one exact analysis record."""

    entry = analysis.visual_entry
    _render_distance_scope(entry=entry)
    with st.expander("Technical analysis details", expanded=False):
        st.caption(
            f"Analysis source: {analysis.source}; cache: {analysis.cache_status}; "
            f"tree authority: {analysis.tree_authority or 'not required/unavailable'}."
        )
        st.write(
            "These codes preserve exactly where the displayed matrix came from. They are "
            "also retained in machine-readable downloads and logs."
        )
    member_species = {
        str(row["member_id"]): str(row["species_label"]) for row in analysis.members
    }
    columns = st.columns(2)
    selected_species = tuple(
        columns[0].multiselect(
            "Highlight species across all views",
            tuple(sorted(set(member_species.values()))),
            key=f"linked_species_{_state_token(key=analysis.key)}",
            help=(
                "Highlights every analysed protein from the selected species in all linked "
                "plots; it does not remove other proteins."
            ),
        )
    )
    selected_members = tuple(
        columns[1].multiselect(
            "Highlight proteins across all views",
            tuple(sorted(member_species)),
            key=f"linked_members_{_state_token(key=analysis.key)}",
            help="Highlights exact protein identifiers in each linked plot.",
        )
    )
    linked = linked_member_ids(
        entry=entry,
        selected_members=selected_members,
        selected_species=selected_species,
    )
    if linked:
        st.info(
            f"Linked selection: {len(linked):,} of {len(member_species):,} displayed "
            f"members across {len({member_species[value] for value in linked}):,} species."
        )
    try:
        geometry = classical_pcoa(rows=analysis.distances, members=analysis.members)
        member_rows = member_dispersion_rows(
            rows=analysis.distances,
            members=analysis.members,
        )
    except InputValidationError as error:
        st.error(f"Exact dispersion analysis failed validation: {error}")
        return
    tabs = st.tabs(
        (
            "Interactive network",
            "Distance spread",
            "PCoA views",
            "Gene-tree phylogram",
            "Distance heatmap",
            "Protein-pair distances",
            "Proteins",
        )
    )
    with tabs[0]:
        _render_force_network(
            entry=entry,
            linked=linked,
            nearest_neighbours=nearest_neighbours,
        )
    with tabs[1]:
        _render_dispersion(
            analysis=analysis,
            member_rows=member_rows,
            linked=linked,
        )
    with tabs[2]:
        _render_enhanced_pcoa(entry=entry, geometry=geometry, linked=linked)
    with tabs[3]:
        _render_phylogram(entry=entry, linked=linked)
    with tabs[4]:
        _render_matrix(entry=entry, linked=linked)
    with tabs[5]:
        _render_pair_table(analysis=analysis)
    with tabs[6]:
        _render_selected_members(members=analysis.members, linked=linked)


def _render_force_network(
    *, entry: dict[str, Any], linked: frozenset[str], nearest_neighbours: int
) -> None:
    """Render the promoted draggable topology and a static alternative."""

    st.write(
        "Use this view to inspect local neighbourhoods and disconnected regions. Node "
        "positions are produced by the interactive layout and are not a distance scale."
    )
    controls = st.columns(3)
    connectors = controls[0].checkbox(
        "Show layout-only connectors",
        value=False,
        help="Dashed connectors join raw components only for layout continuity.",
    )
    labels = controls[1].checkbox(
        "Show every protein label",
        value=False,
        help="Large labels can obscure dense groups; hover remains available when disabled.",
    )
    physics = controls[2].checkbox(
        "Keep layout movement active",
        value=True,
        help="Disable after arranging nodes if a stable view is easier to inspect.",
    )
    try:
        document = force_directed_html(
            entry=entry,
            selected_members=linked,
            include_connectors=connectors,
            show_labels=labels,
            physics_enabled=physics,
        )
    except InputValidationError as error:
        st.warning(str(error))
        return
    render_graph_guidance(key="interactive_network")
    st.iframe(document, height=760)
    metrics = _required_mapping(entry=entry, key="networkMetrics")
    st.caption(
        f"Drag, zoom and pin nodes interactively. Solid edges retain up to "
        f"{nearest_neighbours:,} nearest neighbours per displayed protein. Raw graph: "
        f"{int(metrics.get('rawComponentCount', 0)):,} components and "
        f"{int(metrics.get('rawIsolateCount', 0)):,} isolates. Force-layout spacing is "
        "non-quantitative; components are not OrthoFinder splits."
    )
    with st.expander("Static nearest-neighbour topology (alternative view)"):
        _render_topology(
            entry=entry,
            linked=linked,
            nearest_neighbours=nearest_neighbours,
        )


def _render_dispersion(
    *,
    analysis: GroupAnalysis,
    member_rows: tuple[dict[str, Any], ...],
    linked: frozenset[str],
) -> None:
    """Render complementary exact-distance dispersion summaries."""

    st.write(
        "These views show whether exact pair distances form a tight, broad or multi-scale "
        "distribution. Compare within-species and between-species distances before treating "
        "a group as uniformly compact."
    )
    with st.expander("How to read the distance-spread views"):
        st.markdown(
            """
            - **Histogram:** where pair distances occur most often.
            - **Violin + box:** distribution shape, median and central range for within- and
              between-species pairs.
            - **Empirical CDF:** the fraction of pairs at or below each distance; a curve farther
              left represents generally smaller distances.
            - **Per-protein dispersion:** identifies central and peripheral proteins.
            - **Species-pair heatmap:** compares average distance among each pair of species.
            """
        )
    render_graph_guidance(key="distance_distribution")
    st.plotly_chart(
        distance_distribution_figure(rows=analysis.distances),
        width="stretch",
        config={"displaylogo": False},
    )
    render_graph_guidance(key="member_dispersion")
    st.plotly_chart(
        member_dispersion_figure(rows=member_rows, selected_members=linked),
        width="stretch",
        config={"displaylogo": False},
    )
    secondary = st.tabs(("Distance from sample medoid", "Species-pair heatmap", "Data tables"))
    with secondary[0]:
        render_graph_guidance(key="medoid_distance")
        st.plotly_chart(
            medoid_distance_figure(rows=analysis.distances, member_rows=member_rows),
            width="stretch",
            config={"displaylogo": False},
        )
        st.caption(
            "The medoid minimises mean distance only within the displayed exact or "
            "deterministically sampled members. It is not an ancestor."
        )
    with secondary[1]:
        render_graph_guidance(key="species_pair_heatmap")
        st.plotly_chart(
            species_pair_heatmap_figure(rows=analysis.distances),
            width="stretch",
            config={"displaylogo": False},
        )
    with secondary[2]:
        species_rows = species_dispersion_rows(rows=analysis.distances)
        st.subheader("Per-protein distance summary")
        st.dataframe(
            tuple(_display_member_dispersion_row(row=row) for row in member_rows),
            width="stretch",
            hide_index=True,
            column_config=_column_config(descriptions=_MEMBER_DISPERSION_HELP),
        )
        displayed_members = tuple(
            _display_member_dispersion_row(row=row) for row in member_rows
        )
        render_table_downloads(
            records=displayed_members,
            file_stem=f"{analysis.key.group_id}_member_dispersion",
            key=f"member_dispersion_{_state_token(key=analysis.key)}",
            tsv_label="Download per-protein dispersion as TSV",
            excel_label="Download per-protein dispersion as formatted Excel",
            column_definitions=_MEMBER_DISPERSION_HELP,
            workbook_title=f"Per-protein dispersion: {analysis.key.group_id}",
        )
        st.subheader("Species-pair distance summary")
        st.dataframe(
            tuple(_display_species_dispersion_row(row=row) for row in species_rows),
            width="stretch",
            hide_index=True,
            column_config=_column_config(descriptions=_SPECIES_DISPERSION_HELP),
        )
        displayed_species = tuple(
            _display_species_dispersion_row(row=row) for row in species_rows
        )
        render_table_downloads(
            records=displayed_species,
            file_stem=f"{analysis.key.group_id}_species_pair_dispersion",
            key=f"species_dispersion_{_state_token(key=analysis.key)}",
            tsv_label="Download species-pair dispersion as TSV",
            excel_label="Download species-pair dispersion as formatted Excel",
            column_definitions=_SPECIES_DISPERSION_HELP,
            workbook_title=f"Species-pair dispersion: {analysis.key.group_id}",
        )


def _render_enhanced_pcoa(
    *, entry: dict[str, Any], geometry: Any, linked: frozenset[str]
) -> None:
    """Render original, selectable-axis and rotatable 3D PCoA diagnostics."""

    st.write(
        "PCoA approximates every exact pair distance with two or three plotting axes. "
        "Rotate the 3D view and use the fit measures below to decide how much confidence "
        "to place in apparent screen geometry."
    )
    with st.expander("How to judge a PCoA view"):
        st.markdown(
            """
            - **Variation shown** is the positive coordinate-space inertia retained by the
              displayed axes; larger is better, but it is not a biological variance estimate.
            - **Distortion (stress)** compares plotted and exact distances; lower is better.
            - **Distance agreement** is their Pearson correlation; nearer 1 is better.
            - The **distance-fit check (Shepard plot)** shows exact distance against plotted
              distance directly. Curvature or broad scatter indicates distortion.
            """
        )
    metrics = st.columns(6)
    metrics[0].metric(
        "Variation shown in 2D",
        _format_percent(sum(geometry.axis_fractions[:2])),
        help="Fraction of positive PCoA inertia retained by axes 1 and 2.",
    )
    metrics[1].metric(
        "Variation shown in 3D",
        _format_percent(geometry.cumulative_fraction),
        help="Fraction of positive PCoA inertia retained by axes 1, 2 and 3.",
    )
    metrics[2].metric(
        "2D distortion (stress)",
        _format_optional(geometry.stress_2d),
        help="Normalised distance distortion in two axes; lower values are better.",
    )
    metrics[3].metric(
        "3D distortion (stress)",
        _format_optional(geometry.stress_3d),
        help="Normalised distance distortion in three axes; lower values are better.",
    )
    metrics[4].metric(
        "2D distance agreement",
        _format_optional(geometry.correlation_2d),
        help="Pearson correlation between exact and two-dimensional plotted distances.",
    )
    metrics[5].metric(
        "3D distance agreement",
        _format_optional(geometry.correlation_3d),
        help="Pearson correlation between exact and three-dimensional plotted distances.",
    )
    if geometry.stress_3d < geometry.stress_2d:
        st.info(
            "The third axis reduces distance distortion, but the exact matrix and "
            "branch-length phylogram remain the quantitative authorities."
        )
    views = st.tabs(
        (
            "Rotate 3D",
            "Choose two axes",
            "Original 2D diagnostic",
            "Distance-fit check (Shepard)",
        )
    )
    with views[0]:
        render_graph_guidance(key="pcoa_3d")
        st.plotly_chart(
            pcoa_3d_figure(geometry=geometry, selected_members=linked),
            width="stretch",
            config={"displaylogo": False},
        )
    with views[1]:
        axis_controls = st.columns(2)
        horizontal = int(
            axis_controls[0].selectbox(
                "Horizontal PCoA axis",
                (1, 2, 3),
                index=0,
                help="Choose one of the first three positive-coordinate axes.",
            )
        )
        vertical_options = tuple(value for value in (1, 2, 3) if value != horizontal)
        vertical = int(
            axis_controls[1].selectbox(
                "Vertical PCoA axis",
                vertical_options,
                index=0,
                help="Choose a different positive-coordinate axis.",
            )
        )
        render_graph_guidance(key="pcoa_axes")
        st.plotly_chart(
            pcoa_axis_figure(
                geometry=geometry,
                horizontal_axis=horizontal,
                vertical_axis=vertical,
                selected_members=linked,
            ),
            width="stretch",
            config={"displaylogo": False},
        )
    with views[2]:
        _render_pcoa(entry=entry, linked=linked)
    with views[3]:
        _render_shepard(entry=entry)
    st.caption(
        "PCoA is a diagnostic approximation. Apparent arms, gaps and angles are not "
        "subfamilies without confirmation in the exact matrix and phylogram."
    )


def _render_pair_table(*, analysis: GroupAnalysis) -> None:
    """Render, filter and export actual member-to-member distances."""

    controls = st.columns(3)
    st.write(
        "This is the quantitative authority behind the distribution plots and heatmap. "
        "Each row is one exact pair among the proteins analysed for this view."
    )
    member_text = controls[0].text_input(
        "Either protein ID contains",
        help="Literal, case-insensitive text matched against both pair endpoints.",
    )
    selected_species = tuple(
        controls[1].multiselect(
            "Either endpoint species",
            tuple(
                sorted(
                    {
                        str(row[field])
                        for row in analysis.distances
                        for field in ("species_a", "species_b")
                    }
                )
            ),
        )
    )
    pair_scope = controls[2].selectbox(
        "Species relationship",
        ("All pairs", "Within species", "Between species"),
        help="Choose whether pair endpoints must share or differ in their species label.",
    )
    filtered = _filter_distance_rows(
        rows=analysis.distances,
        member_text=member_text,
        species=selected_species,
        pair_scope=pair_scope,
    )
    st.caption(
        f"Showing {len(filtered):,} of {len(analysis.distances):,} exact displayed pairs."
    )
    display = tuple(_display_pair_distance_row(row=row) for row in filtered)
    st.dataframe(
        display,
        width="stretch",
        hide_index=True,
        column_config=_column_config(descriptions=_PAIR_DISTANCE_HELP),
    )
    render_table_downloads(
        records=display,
        file_stem=f"{analysis.key.group_id}_member_pair_distances",
        key=f"pair_distances_{_state_token(key=analysis.key)}",
        tsv_label="Download filtered protein-pair distances as TSV",
        excel_label="Download filtered protein-pair distances as formatted Excel",
        column_definitions=_PAIR_DISTANCE_HELP,
        workbook_title=f"Protein-pair distances: {analysis.key.group_id}",
    )


def _filter_distance_rows(
    *,
    rows: tuple[dict[str, Any], ...],
    member_text: str,
    species: tuple[str, ...],
    pair_scope: str,
) -> tuple[dict[str, Any], ...]:
    """Apply literal endpoint, species and within/between pair filters."""

    if pair_scope not in {"All pairs", "Within species", "Between species"}:
        raise InputValidationError(f"Unsupported pair scope: {pair_scope}")
    term = member_text.strip().casefold()
    selected_species = set(species)
    output = []
    for row in rows:
        member_a, member_b = str(row["member_a"]), str(row["member_b"])
        species_a, species_b = str(row["species_a"]), str(row["species_b"])
        if term and term not in member_a.casefold() and term not in member_b.casefold():
            continue
        if selected_species and not ({species_a, species_b} & selected_species):
            continue
        if pair_scope == "Within species" and species_a != species_b:
            continue
        if pair_scope == "Between species" and species_a == species_b:
            continue
        output.append(row)
    return tuple(output)


def _comparison_control(*, key: GroupKey) -> None:
    """Add one exact group to the two-to-twelve comparison workspace."""

    basket = _comparison_keys()
    if key in basket:
        st.success(f"This group is in the comparison workspace ({len(basket):,}/12).")
        return
    if st.button("Add this cluster to comparison", disabled=len(basket) >= 12):
        _store_comparison_keys(keys=(*basket, key))
        st.success(f"Added to comparison workspace ({len(basket) + 1:,}/12).")


def _state_group() -> GroupKey | None:
    """Return the active group from safe primitive Streamlit state."""

    raw = st.session_state.get(ACTIVE_GROUP_STATE)
    if not isinstance(raw, dict):
        return None
    fields = ("run_id", "group_type", "hierarchy_node", "group_id")
    if any(field not in raw for field in fields):
        return None
    return GroupKey(**{field: str(raw[field]) for field in fields})


def _store_active_group(*, key: GroupKey) -> None:
    """Store one active group using serialisable exact fields."""

    st.session_state[ACTIVE_GROUP_STATE] = _key_record(key=key)


def _comparison_keys() -> tuple[GroupKey, ...]:
    """Return valid unique comparison groups from Streamlit state."""

    raw = st.session_state.get(COMPARISON_STATE, [])
    if not isinstance(raw, list):
        return ()
    keys = []
    for record in raw:
        if not isinstance(record, dict):
            continue
        try:
            key = GroupKey(
                run_id=str(record["run_id"]),
                group_type=str(record["group_type"]),
                hierarchy_node=str(record["hierarchy_node"]),
                group_id=str(record["group_id"]),
            )
        except KeyError:
            continue
        if key not in keys:
            keys.append(key)
    return tuple(keys[:12])


def _store_comparison_keys(*, keys: tuple[GroupKey, ...]) -> None:
    """Persist a bounded deduplicated comparison workspace."""

    unique = tuple(dict.fromkeys(keys))
    if len(unique) > 12:
        raise InputValidationError("Comparison workspace cannot exceed 12 groups.")
    st.session_state[COMPARISON_STATE] = [_key_record(key=key) for key in unique]


def _key_record(*, key: GroupKey) -> dict[str, str]:
    """Return a serialisable exact group-key record."""

    return {
        "run_id": key.run_id,
        "group_type": key.group_type,
        "hierarchy_node": key.hierarchy_node,
        "group_id": key.group_id,
    }


def _state_token(*, key: GroupKey) -> str:
    """Return a stable Streamlit widget token for one group."""

    return "_".join((key.group_type, key.hierarchy_node or "ROOT", key.group_id))


def _render_distance_scope(*, entry: dict[str, Any]) -> None:
    """Render sample size, exact pair count and distance compactness."""

    summary = _required_mapping(entry=entry, key="distanceSummary")
    columns = st.columns(6)
    columns[0].metric(
        "Proteins in full group",
        f"{int(summary['total_member_count']):,}",
        help="Complete group membership before any bounded distance sampling.",
    )
    columns[1].metric(
        "Proteins analysed",
        f"{int(summary['sampled_member_count']):,}",
        help="Proteins included in every distance view below.",
    )
    columns[2].metric(
        "Protein pairs compared",
        f"{int(summary['distance_pair_count']):,}",
        help="Exact number of pair distances represented below.",
    )
    columns[3].metric(
        "Average pair distance",
        _format_optional(summary.get("mean_distance")),
        help="Mean distance across all displayed protein pairs; smaller is more compact.",
    )
    columns[4].metric(
        "Median pair distance",
        _format_optional(summary.get("median_distance")),
        help="Middle displayed pair distance, which is less affected by extremes than the mean.",
    )
    columns[5].metric(
        "Distance spread (SD)",
        _format_optional(summary.get("population_stddev_distance")),
        help=(
            "Population standard deviation across displayed pair distances; smaller means "
            "the distances are more consistent."
        ),
    )
    st.caption(
        "Every view below describes the exact displayed matrix or deterministic bounded "
        "sample—not uncalculated pairs from omitted proteins."
    )
    with st.expander("Distance calculation details", expanded=False):
        st.code(
            f"method={summary.get('distance_method', 'unavailable')}\n"
            f"status={summary.get('computation_status', 'unavailable')}\n"
            f"analysed_proteins={summary.get('sampled_member_count', 'unavailable')}\n"
            f"full_group_proteins={summary.get('total_member_count', 'unavailable')}",
            language=None,
        )
    medoid = _required_mapping(entry=entry, key="medoid")
    if medoid.get("member_id"):
        st.caption(
            f"★ Sample medoid: {medoid['member_id']} (mean distance "
            f"{_format_optional(medoid.get('mean_distance'))}). It is not a reconstructed "
            "ancestor or guaranteed full-group centre."
        )


def _render_pcoa(*, entry: dict[str, Any], linked: frozenset[str]) -> None:
    """Render the original PCoA with conservative quality diagnostics."""

    projection = _required_mapping(entry=entry, key="distanceProjection")
    if projection.get("status") != "COMPLETE_DISTANCE_PCOA_2D":
        st.warning(f"PCoA unavailable: {projection.get('reason', projection.get('status'))}")
        return
    diagnostics = st.columns(5)
    diagnostics[0].metric(
        "2D fit guidance",
        str(projection.get("quality_category", "")),
        help="Conservative display guidance, not a biological pass/fail classification.",
    )
    diagnostics[1].metric(
        "Variation shown in 2D",
        _format_percent(projection.get("two_axis_positive_inertia_fraction")),
        help="Fraction of positive PCoA inertia retained by axes 1 and 2.",
    )
    diagnostics[2].metric(
        "2D distance agreement",
        _format_optional(projection.get("distance_correlation")),
        help="Pearson correlation between exact and plotted pair distances; nearer 1 is better.",
    )
    diagnostics[3].metric(
        "2D distortion (stress)",
        _format_optional(projection.get("normalised_stress")),
        help="Normalised distortion between exact and plotted pair distances; lower is better.",
    )
    diagnostics[4].metric(
        "Non-Euclidean signal",
        _format_percent(projection.get("negative_inertia_fraction")),
        help=(
            "Absolute negative-inertia fraction. Larger values indicate that the exact "
            "distance matrix is less faithfully represented by Euclidean coordinates."
        ),
    )
    category = str(projection.get("quality_category", ""))
    explanation = str(projection.get("quality_explanation", ""))
    if category == "POOR":
        st.warning(explanation)
    elif category == "MODERATE":
        st.info(explanation)
    else:
        st.success(explanation)
    render_graph_guidance(key="pcoa_2d")
    st.plotly_chart(
        pcoa_figure(entry=entry, selected_members=linked),
        width="stretch",
        config={"displaylogo": False},
    )


def _render_shepard(*, entry: dict[str, Any]) -> None:
    """Render the deterministic Shepard diagnostic."""

    projection = _required_mapping(entry=entry, key="distanceProjection")
    try:
        figure = shepard_figure(entry=entry)
    except InputValidationError as error:
        st.warning(str(error))
        return
    st.write(
        "Points close to the diagonal have similar exact and plotted distances. Broad "
        "scatter or systematic curvature reveals where the PCoA view distorts them."
    )
    render_graph_guidance(key="shepard")
    st.plotly_chart(figure, width="stretch", config={"displaylogo": False})
    st.caption(
        f"Deterministic {int(projection.get('shepard_point_count', 0)):,}-point summary "
        f"across {int(projection.get('shepard_total_pair_count', 0)):,} exact input pairs."
    )


def _render_phylogram(*, entry: dict[str, Any], linked: frozenset[str]) -> None:
    """Render the pruned branch-length gene-tree view."""

    tree = _required_mapping(entry=entry, key="phylogram")
    st.write(
        "This pruned gene tree retains horizontal branch lengths for the analysed proteins. "
        "Use it to test whether patterns in PCoA or the network follow tree structure."
    )
    show_labels = st.checkbox(
        "Show every protein label",
        value=False,
        help="Hover labels remain available when dense permanent labels are hidden.",
    )
    try:
        figure = phylogram_figure(
            entry=entry,
            selected_members=linked,
            show_all_labels=show_labels,
        )
    except InputValidationError as error:
        st.warning(str(error))
        return
    render_graph_guidance(key="phylogram")
    st.plotly_chart(figure, width="stretch", config={"displaylogo": False})
    authority = tree.get("treeAuthority", "RESOLVED_GENE_TREE")
    st.caption(
        f"{tree.get('status', 'Unavailable')} · {authority} "
        f"{tree.get('treeId', 'unavailable')} · {int(tree.get('displayedLeafCount', 0)):,} "
        f"of {int(tree.get('requestedMemberCount', 0)):,} requested leaves. Horizontal "
        "lengths are quantitative branch lengths; vertical spacing is not."
    )
    unresolved = tree.get("unresolvedMembers")
    if isinstance(unresolved, list) and unresolved:
        st.warning(f"Unresolved displayed tree members: {len(unresolved):,}.")


def _render_matrix(*, entry: dict[str, Any], linked: frozenset[str]) -> None:
    """Render the exact matrix with optional linked-selection subsetting."""

    st.write(
        "Each heatmap cell is one exact pair distance. Use this quantitative view to "
        "confirm compact blocks, outliers or gaps suggested by the layouts."
    )
    restrict = st.checkbox(
        "Restrict the matrix to the linked selection",
        value=bool(len(linked) >= 2),
        disabled=len(linked) < 2,
        help="Select at least two linked proteins or species above to enable this subset.",
    )
    selected = tuple(sorted(linked)) if restrict else ()
    try:
        figure = distance_matrix_figure(entry=entry, selected_members=selected)
    except InputValidationError as error:
        st.warning(str(error))
        return
    render_graph_guidance(key="distance_heatmap")
    st.plotly_chart(figure, width="stretch", config={"displaylogo": False})
    matrix = _required_mapping(entry=entry, key="distanceMatrix")
    st.caption(
        f"{matrix.get('status', 'Unavailable')} · order: "
        f"{matrix.get('orderMethod', 'unavailable')}. Every displayed cell retains one "
        "exact supplied pairwise distance."
    )


def _render_topology(
    *, entry: dict[str, Any], linked: frozenset[str], nearest_neighbours: int
) -> None:
    """Render a static nearest-neighbour topology fallback."""

    st.write(
        "This static view preserves nearest-neighbour connections but not quantitative "
        "screen distance. It is an alternative to the draggable network above."
    )
    try:
        figure = nearest_neighbour_figure(
            entry=entry,
            selected_members=linked,
            include_connectors=False,
        )
    except InputValidationError as error:
        st.warning(str(error))
        return
    render_graph_guidance(key="nearest_neighbour")
    st.plotly_chart(figure, width="stretch", config={"displaylogo": False})
    st.caption(
        f"Solid edges retain up to {nearest_neighbours:,} neighbours. Static layout "
        "spacing remains non-quantitative."
    )


def _render_selected_members(
    *, members: tuple[dict[str, Any], ...], linked: frozenset[str]
) -> None:
    """Render and export displayed or linked member records."""

    visible = tuple(
        row for row in members if not linked or str(row["member_id"]) in linked
    )
    st.write(
        "These are the proteins represented in the displayed distance analysis. Linked "
        "selection above can restrict this table without changing the full group."
    )
    st.dataframe(
        tuple(_display_member_row(row=row) for row in visible),
        width="stretch",
        hide_index=True,
        column_config=_column_config(descriptions=_MEMBER_HELP),
    )
    displayed = tuple(_display_member_row(row=row) for row in visible)
    render_table_downloads(
        records=displayed,
        file_stem="orthofinder_displayed_protein_selection",
        key="displayed_member_selection",
        tsv_label="Download displayed protein selection as TSV",
        excel_label="Download displayed protein selection as formatted Excel",
        column_definitions=_MEMBER_HELP,
        workbook_title="OrthoFinder displayed protein selection",
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


def _display_member_dispersion_row(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return readable per-protein distance-summary headings."""

    return {
        "Protein ID": row["member_id"],
        "Species": row["species_label"],
        "Average distance to all others": row["mean_distance"],
        "Median distance to all others": row["median_distance"],
        "Distance spread (SD)": row["population_stddev_distance"],
        "Nearest protein": row["nearest_member_id"],
        "Nearest protein's species": row["nearest_species_label"],
        "Nearest distance": row["nearest_distance"],
        "Sample medoid": "Yes" if row.get("is_sample_medoid") else "No",
    }


def _display_species_dispersion_row(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return readable species-pair distance-summary headings."""

    return {
        "Species A": row["species_a"],
        "Species B": row["species_b"],
        "Pair class": row["pair_class"],
        "Protein pairs": row["pair_count"],
        "Average pair distance": row["mean_distance"],
        "Median pair distance": row["median_distance"],
        "Distance spread (SD)": row["population_stddev_distance"],
        "Smallest pair distance": row["minimum_distance"],
        "Largest pair distance": row["maximum_distance"],
    }


def _display_pair_distance_row(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return readable pair-distance headings while retaining optional site fields."""

    return {
        "Protein A": row["member_a"],
        "Species A": row["species_a"],
        "Protein B": row["member_b"],
        "Species B": row["species_b"],
        "Pair distance": row["distance"],
        "Distance method": row.get("distance_method") or "Unavailable",
        "Calculation scope": row.get("computation_status") or "Unavailable",
        "Comparable sites": row.get("comparable_sites"),
        "Mismatch sites": row.get("mismatch_sites"),
    }


def _display_member_row(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return concise readable fields for one analysed protein."""

    return {
        "Protein ID": row["member_id"],
        "Species": row["species_label"],
    }


def _column_config(*, descriptions: Mapping[str, str]) -> dict[str, Any]:
    """Return Streamlit column definitions with hoverable help descriptions."""

    return {
        label: st.column_config.Column(label=label, help=description)
        for label, description in descriptions.items()
    }


def _required_mapping(*, entry: Mapping[str, Any], key: str) -> dict[str, Any]:
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

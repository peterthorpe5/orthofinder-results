"""Streamlit cluster explorer with linked exact-distance views."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Mapping

import streamlit as st

from orthofinder_results.errors import InputValidationError, OrthoFinderResultsError

from .dispersion import classical_pcoa, member_dispersion_rows, species_dispersion_rows
from .distance_data import DistanceAnalysisProvider, GroupAnalysis
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
from .models import GroupKey
from .queries import OrthoFinderQueryService
from .report_data import VisualisationCatalog, load_visualisation_catalog
from .tsv import records_to_tsv

_LOGGER = logging.getLogger("orthofinder_interrogation_app.evolutionary_page")
ACTIVE_GROUP_STATE = "orthofinder_active_group"
COMPARISON_STATE = "orthofinder_comparison_groups"


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
    st.header("Cluster explorer")
    st.caption(
        "Inspect one group through several linked views. Exact pair distances and the "
        "branch-length phylogram retain quantitative scale; force-layout spacing and PCoA "
        "screen geometry require the stated interpretation checks."
    )
    key = _select_group(service=service, catalog=catalog)
    if key is None:
        st.info("Choose a precomputed group or enter one exact group identifier to begin.")
        return
    controls = st.columns((2, 2, 2, 3))
    max_members = int(
        controls[0].select_slider(
            "Maximum lazy-analysis members",
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
            "Nearest neighbours",
            min_value=1,
            max_value=10,
            value=3,
            help="Controls topology edges, not OrthoFinder group membership.",
        )
    )
    force_recompute = controls[2].checkbox(
        "Recalculate lazy cache",
        value=False,
        help="Persisted schema-2 and DuckDB distances are never recalculated.",
    )
    controls[3].caption(
        "On-demand results are cached in the configured user sidecar. The completed "
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
            "PCoA + diagnostics",
            "Shepard plot",
            "Branch-length phylogram",
            "Exact distance matrix",
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
                "Precomputed pilot groups",
                labels,
                index=default_index,
                help="These open immediately from the existing exact report matrices.",
            )
            if st.button("Use selected pilot group"):
                stored = label_to_key[selected]
                _store_active_group(key=stored)
        with st.form("exact_cluster_key"):
            columns = st.columns(3)
            group_types = service.list_group_types()
            initial_type = stored.group_type if stored is not None else group_types[0]
            group_type = columns[0].selectbox(
                "Exact group type",
                group_types,
                index=group_types.index(initial_type) if initial_type in group_types else 0,
            )
            nodes = service.list_hierarchy_nodes(group_type=group_type)
            initial_node = stored.hierarchy_node if stored is not None else nodes[0]
            hierarchy_node = columns[1].selectbox(
                "Exact hierarchy node",
                nodes,
                format_func=lambda value: value or "ROOT",
                index=nodes.index(initial_node) if initial_node in nodes else 0,
            )
            group_id = columns[2].text_input(
                "Exact group identifier",
                value=stored.group_id if stored is not None else "",
            )
            submitted = st.form_submit_button("Open exact cluster", type="primary")
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
    st.caption(
        f"Analysis source: {analysis.source}; cache: {analysis.cache_status}; tree authority: "
        f"{analysis.tree_authority or 'not required/unavailable'}."
    )
    member_species = {
        str(row["member_id"]): str(row["species_label"]) for row in analysis.members
    }
    columns = st.columns(2)
    selected_species = tuple(
        columns[0].multiselect(
            "Linked species selection",
            tuple(sorted(set(member_species.values()))),
            key=f"linked_species_{_state_token(key=analysis.key)}",
        )
    )
    selected_members = tuple(
        columns[1].multiselect(
            "Linked member selection",
            tuple(sorted(member_species)),
            key=f"linked_members_{_state_token(key=analysis.key)}",
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
            "Interactive force network",
            "Dispersion dashboard",
            "PCoA diagnostics",
            "Branch-length phylogram",
            "Exact matrix",
            "Pair distances",
            "Members",
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

    controls = st.columns(3)
    connectors = controls[0].checkbox(
        "Show layout-only connectors",
        value=False,
        help="Dashed connectors join raw components only for layout continuity.",
    )
    labels = controls[1].checkbox("Show every member label", value=False)
    physics = controls[2].checkbox("Keep node physics active", value=True)
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
    st.iframe(document, height=760)
    metrics = _required_mapping(entry=entry, key="networkMetrics")
    st.caption(
        f"Drag, zoom and pin nodes interactively. Solid edges retain up to "
        f"{nearest_neighbours:,} nearest neighbours per displayed protein. Raw graph: "
        f"{int(metrics.get('rawComponentCount', 0)):,} components and "
        f"{int(metrics.get('rawIsolateCount', 0)):,} isolates. Force-layout spacing is "
        "non-quantitative; components are not OrthoFinder splits."
    )
    with st.expander("Static nearest-neighbour topology"):
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

    st.plotly_chart(
        distance_distribution_figure(rows=analysis.distances),
        width="stretch",
        config={"displaylogo": False},
    )
    st.plotly_chart(
        member_dispersion_figure(rows=member_rows, selected_members=linked),
        width="stretch",
        config={"displaylogo": False},
    )
    secondary = st.tabs(("Distance from sample medoid", "Species-pair heatmap", "Tables"))
    with secondary[0]:
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
        st.plotly_chart(
            species_pair_heatmap_figure(rows=analysis.distances),
            width="stretch",
            config={"displaylogo": False},
        )
    with secondary[2]:
        species_rows = species_dispersion_rows(rows=analysis.distances)
        st.subheader("Per-member dispersion")
        st.dataframe(member_rows, width="stretch", hide_index=True)
        st.download_button(
            "Download per-member dispersion as TSV",
            data=records_to_tsv(records=member_rows),
            file_name=f"{analysis.key.group_id}_member_dispersion.tsv",
            mime="text/tab-separated-values",
        )
        st.subheader("Species-pair dispersion")
        st.dataframe(species_rows, width="stretch", hide_index=True)
        st.download_button(
            "Download species-pair dispersion as TSV",
            data=records_to_tsv(records=species_rows),
            file_name=f"{analysis.key.group_id}_species_pair_dispersion.tsv",
            mime="text/tab-separated-values",
        )


def _render_enhanced_pcoa(
    *, entry: dict[str, Any], geometry: Any, linked: frozenset[str]
) -> None:
    """Render original, selectable-axis and rotatable 3D PCoA diagnostics."""

    metrics = st.columns(6)
    metrics[0].metric("2D positive inertia", _format_percent(sum(geometry.axis_fractions[:2])))
    metrics[1].metric("3D positive inertia", _format_percent(geometry.cumulative_fraction))
    metrics[2].metric("2D stress", _format_optional(geometry.stress_2d))
    metrics[3].metric("3D stress", _format_optional(geometry.stress_3d))
    metrics[4].metric("2D correlation", _format_optional(geometry.correlation_2d))
    metrics[5].metric("3D correlation", _format_optional(geometry.correlation_3d))
    if geometry.stress_3d < geometry.stress_2d:
        st.info(
            "The third axis reduces distance distortion, but the exact matrix and "
            "branch-length phylogram remain the quantitative authorities."
        )
    views = st.tabs(("Rotatable 3D", "Choose two axes", "Original 2D", "Shepard plot"))
    with views[0]:
        st.plotly_chart(
            pcoa_3d_figure(geometry=geometry, selected_members=linked),
            width="stretch",
            config={"displaylogo": False},
        )
    with views[1]:
        axis_controls = st.columns(2)
        horizontal = int(axis_controls[0].selectbox("Horizontal axis", (1, 2, 3), index=0))
        vertical_options = tuple(value for value in (1, 2, 3) if value != horizontal)
        vertical = int(axis_controls[1].selectbox("Vertical axis", vertical_options, index=0))
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
    member_text = controls[0].text_input("Endpoint member contains")
    selected_species = tuple(
        controls[1].multiselect(
            "Endpoint species",
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
        "Pair scope", ("All pairs", "Within species", "Between species")
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
    display_fields = (
        "member_a",
        "species_a",
        "member_b",
        "species_b",
        "distance",
        "distance_method",
        "computation_status",
        "comparable_sites",
        "mismatch_sites",
    )
    display = tuple({field: row.get(field, "") for field in display_fields} for row in filtered)
    st.dataframe(display, width="stretch", hide_index=True)
    st.download_button(
        "Download filtered member-to-member distances as TSV",
        data=records_to_tsv(records=display),
        file_name=f"{analysis.key.group_id}_member_pair_distances.tsv",
        mime="text/tab-separated-values",
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
    columns[0].metric("Analytical members", f"{int(summary['total_member_count']):,}")
    columns[1].metric("Displayed sample", f"{int(summary['sampled_member_count']):,}")
    columns[2].metric("Exact pairs", f"{int(summary['distance_pair_count']):,}")
    columns[3].metric("Mean distance", _format_optional(summary.get("mean_distance")))
    columns[4].metric("Median distance", _format_optional(summary.get("median_distance")))
    columns[5].metric(
        "Distance SD", _format_optional(summary.get("population_stddev_distance"))
    )
    st.caption(
        f"Method: {summary.get('distance_method', 'unavailable')}; status: "
        f"{summary.get('computation_status', 'unavailable')}. Every view describes the "
        "displayed exact or deterministic sample, not uncalculated full-group pairs."
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
    st.plotly_chart(figure, width="stretch", config={"displaylogo": False})
    st.caption(
        f"Deterministic {int(projection.get('shepard_point_count', 0)):,}-point summary "
        f"across {int(projection.get('shepard_total_pair_count', 0)):,} exact input pairs."
    )


def _render_phylogram(*, entry: dict[str, Any], linked: frozenset[str]) -> None:
    """Render the pruned branch-length gene-tree view."""

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

    restrict = st.checkbox(
        "Restrict the matrix to the linked selection",
        value=bool(len(linked) >= 2),
        disabled=len(linked) < 2,
    )
    selected = tuple(sorted(linked)) if restrict else ()
    try:
        figure = distance_matrix_figure(entry=entry, selected_members=selected)
    except InputValidationError as error:
        st.warning(str(error))
        return
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

    try:
        figure = nearest_neighbour_figure(
            entry=entry,
            selected_members=linked,
            include_connectors=False,
        )
    except InputValidationError as error:
        st.warning(str(error))
        return
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
    st.dataframe(visible, width="stretch", hide_index=True)
    st.download_button(
        "Download displayed member selection as TSV",
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

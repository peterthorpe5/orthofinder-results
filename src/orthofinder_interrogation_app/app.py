"""Streamlit interface for generic OrthoFinder resource interrogation."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

import streamlit as st

from orthofinder_interrogation_app.comparison_page import render_cluster_comparison
from orthofinder_interrogation_app.distance_data import default_cache_directory
from orthofinder_interrogation_app.evolutionary_page import (
    _comparison_keys,
    _store_active_group,
    _store_comparison_keys,
    render_evolutionary_views,
)
from orthofinder_interrogation_app.launcher import (
    CACHE_ENVIRONMENT_VARIABLE,
    LOG_ENVIRONMENT_VARIABLE,
    RESOURCE_ENVIRONMENT_VARIABLE,
    TAXONOMY_ENVIRONMENT_VARIABLE,
)
from orthofinder_interrogation_app.models import GroupKey, GroupSearchFilters
from orthofinder_interrogation_app.queries import (
    MAX_GROUP_MEMBER_ROWS,
    OrthoFinderQueryService,
)
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_interrogation_app.taxonomy_page import render_taxonomy_search
from orthofinder_interrogation_app.tsv import records_to_tsv
from orthofinder_results import __version__
from orthofinder_results.errors import OrthoFinderResultsError
from orthofinder_results.io_utils import configure_logging

_LOGGER = logging.getLogger("orthofinder_interrogation_app.app")
_PAGES = (
    "Overview",
    "Find groups",
    "Cluster explorer",
    "Compare clusters",
    "Taxonomic search",
    "Offline report",
    "Help",
)
_INCLUDE_LABEL_TO_MODE = {
    "Any selected species": "ANY",
    "Every selected species": "ALL",
    "Exactly this species set": "EXACT_SET",
}
_DISTANCE_LABEL_TO_MODE = {
    "Any persisted distance state": "ANY",
    "Persisted distances available": "CALCULATED",
    "No persisted distances": "NOT_CALCULATED",
}
_SORT_LABEL_TO_MODE = {
    "Group identifier": "GROUP_ID_ASC",
    "Largest groups": "MEMBER_COUNT_DESC",
    "Broadest species representation": "SPECIES_COUNT_DESC",
    "Most compact persisted mean": "MEAN_DISTANCE_ASC",
    "Most divergent persisted mean": "MEAN_DISTANCE_DESC",
    "Lowest persisted distance variability": "DISTANCE_SD_ASC",
}


def main() -> None:
    """Render the complete local read-only application."""

    st.set_page_config(
        page_title="OrthoFinder Interrogation",
        page_icon="🧬",
        layout="wide",
    )
    _inject_style()
    _configure_application_logging()
    st.title("OrthoFinder Interrogation")
    st.caption(
        "Search the complete group authority, calculate bounded tree distances on demand, "
        "and compare evolutionary dispersion without modifying the completed resource."
    )
    default_resource = os.environ.get(RESOURCE_ENVIRONMENT_VARIABLE, "")
    resource_text = st.sidebar.text_input(
        "Resource directory or DuckDB",
        value=default_resource,
        help="Choose a completed resource directory or its DuckDB file.",
    )
    page_name = st.sidebar.radio("Page", _PAGES, key="app_page")
    default_cache = os.environ.get(
        CACHE_ENVIRONMENT_VARIABLE,
        str(default_cache_directory()),
    )
    with st.sidebar.expander("Local sidecars and logs"):
        cache_text = st.text_input(
            "Analysis cache directory",
            value=default_cache,
            help="Must remain outside the immutable completed resource.",
        )
        taxonomy_path_text = st.text_input(
            "Reviewed taxonomy TSV",
            value=os.environ.get(TAXONOMY_ENVIRONMENT_VARIABLE, ""),
            help="Optional versioned sidecar; mappings are never guessed.",
        )
    st.sidebar.caption(f"Application package {__version__}")
    if not resource_text.strip():
        st.info("Provide a completed resource path in the sidebar to begin.")
        return
    try:
        resource = open_resource(path=Path(resource_text))
        service = OrthoFinderQueryService(resource=resource)
        cache_dir = Path(cache_text)
        _render_resource_identity(resource=resource)
        if page_name == "Overview":
            _render_overview(service=service)
        elif page_name == "Find groups":
            _render_group_search(service=service)
        elif page_name == "Cluster explorer":
            render_evolutionary_views(
                resource=resource,
                service=service,
                cache_dir=cache_dir,
            )
        elif page_name == "Compare clusters":
            render_cluster_comparison(
                resource=resource,
                service=service,
                cache_dir=cache_dir,
            )
        elif page_name == "Taxonomic search":
            render_taxonomy_search(
                service=service,
                taxonomy_path_text=taxonomy_path_text,
            )
        elif page_name == "Offline report":
            _render_offline_report(resource=resource)
        else:
            _render_help()
    except OrthoFinderResultsError as error:
        _LOGGER.exception("Application request failed")
        st.error(f"The resource could not be opened or queried: {error}")


def _inject_style() -> None:
    """Apply a restrained readable width and compact metric styling."""

    st.markdown(
        """
        <style>
        .stMainBlockContainer {max-width: 1680px; padding-top: 2rem;}
        [data-testid="stMetric"] {
          border: 1px solid #e2e8f0;
          border-radius: 0.55rem;
          padding: 0.65rem 0.8rem;
          background: #fbfdff;
        }
        [data-testid="stSidebar"] {min-width: 19rem;}
        </style>
        """,
        unsafe_allow_html=True,
    )


def _configure_application_logging() -> None:
    """Configure console and optional persistent logging for this process."""

    log_text = os.environ.get(LOG_ENVIRONMENT_VARIABLE, "").strip()
    configure_logging(log_path=Path(log_text) if log_text else None, verbose=False)


def _render_resource_identity(*, resource: Any) -> None:
    """Render concise immutable identity above every application page."""

    with st.expander("Resource identity", expanded=False):
        columns = st.columns(5)
        columns[0].metric("Run", resource.run_id)
        columns[1].metric("OrthoFinder", resource.orthofinder_version)
        columns[2].metric("Adapter", resource.adapter_name)
        columns[3].metric("Resource schema", resource.schema_version)
        columns[4].metric("Resource package", resource.resource_package_version)
        st.code(str(resource.database_path), language=None)
        st.caption(
            "Resource package is the version that built this database; application package "
            "in the sidebar is the currently installed viewer."
        )


def _render_overview(*, service: OrthoFinderQueryService) -> None:
    """Render an actionable run summary instead of raw metadata."""

    st.header("Run overview")
    counts = service.overview_counts()
    columns = st.columns(5)
    columns[0].metric("Groups", f"{counts['group_count']:,}")
    columns[1].metric("Species", f"{counts['species_count']:,}")
    columns[2].metric("Group/species rows", f"{counts['group_species_statistic_count']:,}")
    columns[3].metric("Precomputed distance groups", f"{counts['distance_group_count']:,}")
    columns[4].metric("Portable gene trees", f"{counts['portable_tree_count']:,}")
    if counts["portable_tree_count"]:
        st.success(
            "This schema-3 resource can calculate bounded patristic distances for any group "
            "that maps unambiguously to a portable gene tree. Results use an external cache."
        )
    else:
        st.warning(
            "This schema-2 resource preserves its precomputed pilot groups, but other groups "
            "need a schema-3 rebuild before the app can calculate distances on demand."
        )
    st.subheader("Group authorities")
    authorities = tuple(_display_authority_row(row=row) for row in service.overview_authorities())
    st.dataframe(authorities, width="stretch", hide_index=True)
    actions = st.columns(3)
    actions[0].subheader("1. Find")
    actions[0].write(
        "Filter groups by exact identifiers, included species, inclusive rules, and rejected "
        "species."
    )
    actions[1].subheader("2. Explore")
    actions[1].write(
        "Open a cluster for force topology, 2D/3D PCoA, exact matrices and several dispersion "
        "views."
    )
    actions[2].subheader("3. Compare")
    actions[2].write(
        "Collect 2–12 groups and compare mean, SD, medians, distributions and independent "
        "PCoA panels."
    )


def _render_group_search(*, service: OrthoFinderQueryService) -> None:
    """Render progressive exact group filters and reusable selected-group actions."""

    st.header("Find groups")
    st.caption(
        "Species controls use exact labels from this run. Required species can use ANY, ALL "
        "or EXACT-SET semantics; rejected species exclude a group if any selected label occurs."
    )
    species = service.list_species()
    group_types = service.list_group_types()
    with st.form("group_filters"):
        first_row = st.columns(4)
        group_type = first_row[0].selectbox("Group type", ("Any", *group_types))
        nodes = service.list_hierarchy_nodes(
            group_type="" if group_type == "Any" else group_type
        )
        node_labels = tuple("ROOT" if not node else node for node in nodes)
        hierarchy_label = first_row[1].selectbox("Hierarchy node", ("Any", *node_labels))
        group_text = first_row[2].text_input("Group identifier contains")
        member_text = first_row[3].text_input("Member identifier contains")
        species_row = st.columns((3, 2, 3))
        included_species = tuple(
            species_row[0].multiselect("Species to include", options=species)
        )
        include_label = species_row[1].selectbox(
            "Inclusion rule",
            tuple(_INCLUDE_LABEL_TO_MODE),
            index=1,
        )
        excluded_species = tuple(
            species_row[2].multiselect("Reject groups containing", options=species)
        )
        with st.expander("Advanced size, distance, sorting and paging filters"):
            size_row = st.columns(4)
            minimum_members = int(
                size_row[0].number_input("Minimum members", min_value=1, value=1)
            )
            maximum_members_enabled = size_row[1].checkbox("Set maximum members")
            maximum_members = int(
                size_row[1].number_input(
                    "Maximum members",
                    min_value=1,
                    value=max(minimum_members, 1000),
                    disabled=not maximum_members_enabled,
                )
            )
            minimum_species = int(
                size_row[2].number_input("Minimum represented species", min_value=1, value=1)
            )
            maximum_species_enabled = size_row[3].checkbox("Set maximum species")
            maximum_species = int(
                size_row[3].number_input(
                    "Maximum represented species",
                    min_value=1,
                    value=max(minimum_species, len(species)),
                    disabled=not maximum_species_enabled,
                )
            )
            distance_row = st.columns(4)
            distance_label = distance_row[0].selectbox(
                "Persisted distance availability",
                tuple(_DISTANCE_LABEL_TO_MODE),
            )
            mean_filter_enabled = distance_row[1].checkbox("Limit persisted mean")
            maximum_mean_distance = float(
                distance_row[1].number_input(
                    "Maximum mean distance",
                    min_value=0.0,
                    value=1.0,
                    disabled=not mean_filter_enabled,
                )
            )
            sd_filter_enabled = distance_row[2].checkbox("Limit persisted SD")
            maximum_distance_sd = float(
                distance_row[2].number_input(
                    "Maximum distance SD",
                    min_value=0.0,
                    value=1.0,
                    disabled=not sd_filter_enabled,
                )
            )
            sort_label = distance_row[3].selectbox("Sort", tuple(_SORT_LABEL_TO_MODE))
            paging = st.columns(2)
            page_size = int(
                paging[0].selectbox("Rows per page", (25, 50, 100, 250, 500), index=2)
            )
            page_number = int(paging[1].number_input("Page", min_value=1, value=1))
        st.form_submit_button("Find groups", type="primary")
    hierarchy_node = (
        None
        if hierarchy_label == "Any"
        else ""
        if hierarchy_label == "ROOT"
        else hierarchy_label
    )
    filters = GroupSearchFilters(
        group_type="" if group_type == "Any" else group_type,
        hierarchy_node=hierarchy_node,
        group_id_contains=group_text,
        member_id_contains=member_text,
        included_species=included_species,
        include_mode=_INCLUDE_LABEL_TO_MODE[include_label],
        excluded_species=excluded_species,
        minimum_member_count=minimum_members,
        maximum_member_count=maximum_members if maximum_members_enabled else None,
        minimum_species_count=minimum_species,
        maximum_species_count=maximum_species if maximum_species_enabled else None,
        distance_availability=_DISTANCE_LABEL_TO_MODE[distance_label],
        maximum_mean_distance=maximum_mean_distance if mean_filter_enabled else None,
        maximum_distance_sd=maximum_distance_sd if sd_filter_enabled else None,
        sort_mode=_SORT_LABEL_TO_MODE[sort_label],
        page_size=page_size,
        page_number=page_number,
    )
    result = service.search_groups(filters=filters)
    start = filters.offset + 1 if result.rows else 0
    finish = filters.offset + len(result.rows)
    st.subheader("Matching groups")
    st.caption(f"Showing {start:,}–{finish:,} of {result.total_rows:,} matching groups.")
    if not result.rows:
        st.warning("No groups match the selected filters.")
        return
    st.dataframe(
        tuple(_display_group_row(row=row) for row in result.rows),
        width="stretch",
        hide_index=True,
    )
    st.download_button(
        "Download this result page as TSV",
        data=records_to_tsv(records=result.rows),
        file_name="orthofinder_group_search.tsv",
        mime="text/tab-separated-values",
    )
    labels_to_keys = {
        _row_group_key(row=row).display_label(): _row_group_key(row=row)
        for row in result.rows
    }
    selected_label = st.selectbox("Selected matching group", tuple(labels_to_keys))
    selected_key = labels_to_keys[selected_label]
    actions = st.columns(2)
    if actions[0].button("Open selected group in Cluster explorer", type="primary"):
        _store_active_group(key=selected_key)
        st.session_state["app_page"] = "Cluster explorer"
        st.rerun()
    basket = _comparison_keys()
    if actions[1].button(
        "Add selected group to comparison",
        disabled=selected_key in basket or len(basket) >= 12,
    ):
        _store_comparison_keys(keys=(*basket, selected_key))
        st.success(f"Added to comparison workspace ({len(basket) + 1:,}/12).")
    with st.expander("Preview selected group membership and copy counts"):
        _render_group_detail(service=service, key=selected_key)


def _render_group_detail(*, service: OrthoFinderQueryService, key: GroupKey) -> None:
    """Render one group's exact copy counts, members and persisted distance summary."""

    group = service.get_group(key=key)
    st.subheader(key.display_label())
    columns = st.columns(4)
    columns[0].metric("Members", f"{group['member_count']:,}")
    columns[1].metric("Species", f"{group['species_count']:,}")
    columns[2].metric("Maximum copies/species", f"{group['max_copies_per_species']:,}")
    columns[3].metric("Mean copies/species", _format_number(group["mean_copies_per_species"]))
    distance_count = group.get("distance_pair_count")
    if distance_count:
        distance_columns = st.columns(4)
        distance_columns[0].metric("Persisted mean", _format_number(group["mean_distance"]))
        distance_columns[1].metric(
            "Persisted SD", _format_number(group["population_stddev_distance"])
        )
        distance_columns[2].metric("Persisted median", _format_number(group["median_distance"]))
        distance_columns[3].metric("Persisted pairs", f"{distance_count:,}")
        st.caption(
            f"Method: {group['distance_method']}; status: {group['computation_status']}; "
            f"sampled members: {group['sampled_member_count']:,} of "
            f"{group['total_member_count']:,}."
        )
    elif service.has_relation(relation="tree_payloads"):
        st.info(
            "No distance pairs were precomputed for this group. Cluster explorer can try a "
            "bounded on-demand calculation from its portable gene tree."
        )
    else:
        st.info(
            "No distances were calculated for this schema-2 group. A schema-3 resource "
            "rebuild is required for on-demand calculation."
        )
    species_rows = service.get_group_species(key=key)
    member_rows = service.get_group_members(key=key)
    species_tab, member_tab = st.tabs(("Species copy counts", "Members"))
    with species_tab:
        st.dataframe(species_rows, width="stretch", hide_index=True)
        st.download_button(
            "Download species copy counts as TSV",
            data=records_to_tsv(records=species_rows),
            file_name=f"{key.group_id}_species_copy_counts.tsv",
            mime="text/tab-separated-values",
        )
    with member_tab:
        st.dataframe(member_rows, width="stretch", hide_index=True)
        if len(member_rows) == group["member_count"]:
            st.download_button(
                "Download complete membership as TSV",
                data=records_to_tsv(records=member_rows),
                file_name=f"{key.group_id}_members.tsv",
                mime="text/tab-separated-values",
            )
        else:
            st.warning(
                f"The browser preview is limited to {MAX_GROUP_MEMBER_ROWS:,} members; "
                f"this group contains {group['member_count']:,}. Use DuckDB/Parquet for "
                "a complete export."
            )


def _render_offline_report(*, resource: Any) -> None:
    """Expose the immutable offline report as a download."""

    st.header("Offline report")
    if resource.report_path is None:
        st.warning("This resource does not contain an offline HTML report.")
        return
    st.write(
        "The report remains a self-contained immutable export. The application now adds "
        "lazy schema-3 analysis, actual pair tables, 3D PCoA and cross-group comparison."
    )
    st.download_button(
        "Download self-contained HTML report",
        data=resource.report_path.read_bytes(),
        file_name=resource.report_path.name,
        mime="text/html",
    )


def _render_help() -> None:
    """Explain search semantics and conservative scientific interpretation."""

    st.header("Help and interpretation")
    st.markdown(
        """
        - **ANY selected species** retains groups containing at least one selected label.
          **ALL** requires every selected label but allows other species. **EXACT SET**
          requires precisely the selected sampled species.
        - **Reject groups containing** removes any group containing at least one rejected
          exact species label.
        - Mean, median and SD describe only the stated distance method and displayed exact
          or deterministic member sample. Missing values are not biological zero.
        - A schema-3 resource stores checksum-bound portable gene trees, not all possible
          pair matrices. The app calculates a selected bounded matrix and caches it outside
          the immutable resource. Schema 2 continues to expose its original pilot matrices.
        - Taxonomic ancestry uses a separately reviewed mapping. Unmapped, pending and
          ambiguous labels never become inferred outsiders. **Sampled exclusive** means
          exclusive only among the species in this analysis, not universal absence.
        - PCoA and force layouts are diagnostic. PCoA arms are not automatic subfamilies,
          force-layout spacing is non-quantitative, and nearest-neighbour components are
          not OrthoFinder splits. Confirm detailed geometry in the exact matrix and
          branch-length phylogram.
        - The sample medoid is the displayed member with the lowest mean distance to other
          displayed members. It is not a reconstructed ancestor or guaranteed full-group
          centre.
        """
    )


def _display_authority_row(*, row: dict[str, Any]) -> dict[str, Any]:
    """Return concise headings for one overview authority row."""

    return {
        "Group type": row["group_type"],
        "Hierarchy": row["hierarchy_node"] or "ROOT",
        "Groups": row["group_count"],
        "Minimum members": row["minimum_members"],
        "Median members": row["median_members"],
        "Maximum members": row["maximum_members"],
        "Mean represented species": row["mean_species"],
        "Maximum represented species": row["maximum_species"],
    }


def _display_group_row(*, row: dict[str, Any]) -> dict[str, Any]:
    """Return concise user-facing headings for a group-search row."""

    return {
        "Group type": row["group_type"],
        "Hierarchy": row["hierarchy_node"] or "ROOT",
        "Group": row["group_id"],
        "Members": row["member_count"],
        "Species": row["species_count"],
        "Maximum copies/species": row["max_copies_per_species"],
        "Mean copies/species": row["mean_copies_per_species"],
        "Persisted distance status": row["computation_status"] or "Not calculated",
        "Persisted mean distance": row["mean_distance"],
        "Persisted distance SD": row["population_stddev_distance"],
    }


def _row_group_key(*, row: dict[str, Any]) -> GroupKey:
    """Return one composite group key from a search result."""

    return GroupKey(
        run_id=str(row["run_id"]),
        group_type=str(row["group_type"]),
        hierarchy_node=str(row["hierarchy_node"]),
        group_id=str(row["group_id"]),
    )


def _format_number(value: object) -> str:
    """Format an optional numeric value without converting absence to zero."""

    if value is None:
        return "Unavailable"
    return f"{float(value):.4g}"


if __name__ == "__main__":
    main()

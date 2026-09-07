"""Streamlit user interface for generic OrthoFinder resource interrogation."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

import streamlit as st

from orthofinder_interrogation_app.evolutionary_page import render_evolutionary_views
from orthofinder_interrogation_app.launcher import (
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
    "Group search",
    "Evolutionary views",
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
    "Any": "ANY",
    "Calculated": "CALCULATED",
    "Not calculated": "NOT_CALCULATED",
}
_SORT_LABEL_TO_MODE = {
    "Largest groups": "MEMBER_COUNT_DESC",
    "Broadest species representation": "SPECIES_COUNT_DESC",
    "Group identifier": "GROUP_ID_ASC",
    "Most compact by mean distance": "MEAN_DISTANCE_ASC",
    "Most divergent by mean distance": "MEAN_DISTANCE_DESC",
    "Lowest distance variability": "DISTANCE_SD_ASC",
}


def main() -> None:
    """Render the complete Streamlit application."""

    st.set_page_config(
        page_title="OrthoFinder Interrogation",
        page_icon="🧬",
        layout="wide",
    )
    _configure_application_logging()
    st.title("OrthoFinder Interrogation")
    st.caption(
        "Read-only exploration of a versioned OrthoFinder resource. "
        "Completed resources are never modified."
    )
    default_resource = os.environ.get(RESOURCE_ENVIRONMENT_VARIABLE, "")
    resource_text = st.sidebar.text_input(
        "Resource directory or DuckDB",
        value=default_resource,
        help="Choose a completed resource directory or its DuckDB file.",
    )
    page_name = st.sidebar.radio("Page", _PAGES)
    taxonomy_path_text = ""
    if page_name == "Taxonomic search":
        taxonomy_path_text = st.sidebar.text_input(
            "Reviewed taxonomy TSV",
            value=os.environ.get(TAXONOMY_ENVIRONMENT_VARIABLE, ""),
            help="Optional versioned sidecar mapping; the completed resource stays read-only.",
        )
    st.sidebar.caption(f"Application package {__version__}")
    if not resource_text.strip():
        st.info("Provide a completed resource path in the sidebar to begin.")
        return
    try:
        resource = open_resource(path=Path(resource_text))
        service = OrthoFinderQueryService(resource=resource)
        _render_resource_identity(resource=resource)
        if page_name == "Overview":
            _render_overview(service=service)
        elif page_name == "Group search":
            _render_group_search(service=service)
        elif page_name == "Evolutionary views":
            render_evolutionary_views(resource=resource)
        elif page_name == "Taxonomic search":
            render_taxonomy_search(service=service, taxonomy_path_text=taxonomy_path_text)
        elif page_name == "Offline report":
            _render_offline_report(resource=resource)
        else:
            _render_help()
    except OrthoFinderResultsError as error:
        _LOGGER.exception("Application request failed")
        st.error(f"The resource could not be opened or queried: {error}")


def _configure_application_logging() -> None:
    """Configure console and optional persistent logging for this process."""

    log_text = os.environ.get(LOG_ENVIRONMENT_VARIABLE, "").strip()
    configure_logging(log_path=Path(log_text) if log_text else None, verbose=False)


def _render_resource_identity(*, resource: Any) -> None:
    """Render immutable run identity above every application page."""

    with st.expander("Resource identity", expanded=False):
        columns = st.columns(5)
        columns[0].metric("Run", resource.run_id)
        columns[1].metric("OrthoFinder", resource.orthofinder_version)
        columns[2].metric("Adapter", resource.adapter_name)
        columns[3].metric("Schema", resource.schema_version)
        columns[4].metric("Created by package", resource.resource_package_version)
        st.code(str(resource.database_path), language=None)


def _render_overview(*, service: OrthoFinderQueryService) -> None:
    """Render complete-authority resource counts and scope."""

    st.header("Run overview")
    counts = service.overview_counts()
    columns = st.columns(4)
    columns[0].metric("Groups", f"{counts['group_count']:,}")
    columns[1].metric("Group/species rows", f"{counts['group_species_statistic_count']:,}")
    columns[2].metric("Species", f"{counts['species_count']:,}")
    columns[3].metric("Groups with distances", f"{counts['distance_group_count']:,}")
    st.info(
        "Group and species counts use the complete DuckDB authority. Detailed distance "
        "statistics are available only for groups calculated in this resource."
    )
    st.subheader("Available authorities")
    st.write(
        {
            "Group types": ", ".join(service.list_group_types()),
            "Hierarchy levels": ", ".join(
                node or "ROOT" for node in service.list_hierarchy_nodes()
            ),
            "Exact species labels": len(service.list_species()),
        }
    )


def _render_group_search(*, service: OrthoFinderQueryService) -> None:
    """Render lazy group filters, one result page and selected-group details."""

    st.header("Group search")
    st.caption("Species include/reject filters use exact labels stored in this OrthoFinder run.")
    species = service.list_species()
    group_types = service.list_group_types()
    with st.form("group_filters"):
        first_row = st.columns(4)
        group_type = first_row[0].selectbox("Group type", ("Any", *group_types))
        nodes = service.list_hierarchy_nodes(group_type="" if group_type == "Any" else group_type)
        node_labels = tuple("ROOT" if not node else node for node in nodes)
        hierarchy_label = first_row[1].selectbox("Hierarchy node", ("Any", *node_labels))
        group_text = first_row[2].text_input("Group identifier contains")
        member_text = first_row[3].text_input("Member identifier contains")

        second_row = st.columns(3)
        included_species = tuple(second_row[0].multiselect("Species to include", options=species))
        include_label = second_row[1].selectbox("Include rule", tuple(_INCLUDE_LABEL_TO_MODE))
        excluded_species = tuple(
            second_row[2].multiselect("Reject groups containing", options=species)
        )

        third_row = st.columns(4)
        minimum_members = int(third_row[0].number_input("Minimum members", min_value=1, value=1))
        maximum_members_enabled = third_row[1].checkbox("Set maximum members")
        maximum_members = int(
            third_row[1].number_input(
                "Maximum members",
                min_value=1,
                value=max(minimum_members, 1000),
                disabled=not maximum_members_enabled,
            )
        )
        minimum_species = int(
            third_row[2].number_input("Minimum represented species", min_value=1, value=1)
        )
        maximum_species_enabled = third_row[3].checkbox("Set maximum species")
        maximum_species = int(
            third_row[3].number_input(
                "Maximum represented species",
                min_value=1,
                value=max(minimum_species, len(species)),
                disabled=not maximum_species_enabled,
            )
        )

        fourth_row = st.columns(4)
        distance_label = fourth_row[0].selectbox(
            "Distance availability", tuple(_DISTANCE_LABEL_TO_MODE)
        )
        mean_filter_enabled = fourth_row[1].checkbox("Limit mean distance")
        maximum_mean_distance = float(
            fourth_row[1].number_input(
                "Maximum mean distance",
                min_value=0.0,
                value=1.0,
                disabled=not mean_filter_enabled,
            )
        )
        sd_filter_enabled = fourth_row[2].checkbox("Limit distance SD")
        maximum_distance_sd = float(
            fourth_row[2].number_input(
                "Maximum distance SD",
                min_value=0.0,
                value=1.0,
                disabled=not sd_filter_enabled,
            )
        )
        sort_label = fourth_row[3].selectbox("Sort", tuple(_SORT_LABEL_TO_MODE))
        final_row = st.columns(2)
        page_size = int(final_row[0].selectbox("Rows per page", (25, 50, 100, 250, 500), index=2))
        page_number = int(final_row[1].number_input("Page", min_value=1, value=1))
        st.form_submit_button("Search groups", type="primary")

    hierarchy_node = (
        None if hierarchy_label == "Any" else "" if hierarchy_label == "ROOT" else hierarchy_label
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
    table_rows = [_display_group_row(row=row) for row in result.rows]
    st.dataframe(table_rows, use_container_width=True, hide_index=True)
    st.download_button(
        "Download this result page as TSV",
        data=records_to_tsv(records=result.rows),
        file_name="orthofinder_group_search.tsv",
        mime="text/tab-separated-values",
    )
    labels_to_keys = {
        GroupKey(
            run_id=str(row["run_id"]),
            group_type=str(row["group_type"]),
            hierarchy_node=str(row["hierarchy_node"]),
            group_id=str(row["group_id"]),
        ).display_label(): GroupKey(
            run_id=str(row["run_id"]),
            group_type=str(row["group_type"]),
            hierarchy_node=str(row["hierarchy_node"]),
            group_id=str(row["group_id"]),
        )
        for row in result.rows
    }
    selected_label = st.selectbox("Inspect one matching group", tuple(labels_to_keys))
    _render_group_detail(service=service, key=labels_to_keys[selected_label])


def _render_group_detail(*, service: OrthoFinderQueryService, key: GroupKey) -> None:
    """Render one group's exact copy counts, members and distance summary."""

    group = service.get_group(key=key)
    st.divider()
    st.header(key.display_label())
    columns = st.columns(4)
    columns[0].metric("Members", f"{group['member_count']:,}")
    columns[1].metric("Species", f"{group['species_count']:,}")
    columns[2].metric("Maximum copies/species", f"{group['max_copies_per_species']:,}")
    columns[3].metric("Mean copies/species", _format_number(group["mean_copies_per_species"]))
    distance_count = group.get("distance_pair_count")
    if distance_count:
        st.subheader("Within-group compactness")
        distance_columns = st.columns(4)
        distance_columns[0].metric("Mean distance", _format_number(group["mean_distance"]))
        distance_columns[1].metric(
            "Distance SD", _format_number(group["population_stddev_distance"])
        )
        distance_columns[2].metric("Median distance", _format_number(group["median_distance"]))
        distance_columns[3].metric("Pair distances", f"{distance_count:,}")
        st.caption(
            f"Method: {group['distance_method']}; status: {group['computation_status']}; "
            f"sampled members: {group['sampled_member_count']:,} of "
            f"{group['total_member_count']:,}. Lower mean and SD indicate a more compact "
            "calculated distance distribution, but interpretation remains method- and "
            "sampling-dependent."
        )
    else:
        st.info("Distance compactness was not calculated for this group in this resource.")
    species_rows = service.get_group_species(key=key)
    member_rows = service.get_group_members(key=key)
    species_tab, member_tab = st.tabs(("Species copy counts", "Members"))
    with species_tab:
        st.dataframe(species_rows, use_container_width=True, hide_index=True)
        st.download_button(
            "Download species copy counts as TSV",
            data=records_to_tsv(records=species_rows),
            file_name=f"{key.group_id}_species_copy_counts.tsv",
            mime="text/tab-separated-values",
        )
    with member_tab:
        st.dataframe(member_rows, use_container_width=True, hide_index=True)
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
                f"this group contains {group['member_count']:,}. Use the DuckDB/Parquet "
                "authority for a complete export."
            )


def _render_offline_report(*, resource: Any) -> None:
    """Expose the existing immutable offline report as a download."""

    st.header("Offline report")
    if resource.report_path is None:
        st.warning("This resource does not contain an offline HTML report.")
        return
    st.write(
        "The report remains an immutable, self-contained export. Its bounded visual "
        "payload also supplies schema-2 phylograms to the application; DuckDB remains "
        "the complete search authority."
    )
    st.download_button(
        "Download self-contained HTML report",
        data=resource.report_path.read_bytes(),
        file_name=resource.report_path.name,
        mime="text/html",
    )


def _render_help() -> None:
    """Explain current search and scientific interpretation contracts."""

    st.header("Help and interpretation")
    st.markdown(
        """
        - **Any selected species** retains groups containing at least one selected label.
        - **Every selected species** retains groups containing all selected labels; other
          species may also occur.
        - **Exactly this species set** requires all selected labels and no other represented
          species in this run.
        - **Reject groups containing** removes any group containing at least one rejected
          species.
        - Mean and SD distances describe only the stated distance method and calculated
          member sample. Missing distance values are not biological zero.
        - These exact-label filters do not imply taxonomic ancestry. Descendant-aware
          filtering uses a separately reviewed NCBI taxonomy mapping. Unmapped and
          ambiguous labels remain visible and never become inferred outsiders.
        - **Contains** requires at least one reviewed sampled descendant. **Exclusive**
          and **near-exclusive** are scoped only to this analysis's sampled species;
          neither is a universal absence claim.
        - **Enriched** uses one-sided Fisher exact tests on species presence and
          Benjamini–Hochberg correction across the complete selected group authority.
        - PCoA and force layouts are diagnostic screen views. PCoA arms are not
          subfamilies, force-layout spacing is non-quantitative, and graph components
          are not OrthoFinder splits.
        """
    )


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
        "Distance status": row["computation_status"] or "Not calculated",
        "Mean distance": row["mean_distance"],
        "Distance SD": row["population_stddev_distance"],
    }


def _format_number(value: object) -> str:
    """Format an optional numeric value without converting absence to zero."""

    if value is None:
        return "Unavailable"
    return f"{float(value):.4g}"


if __name__ == "__main__":
    main()

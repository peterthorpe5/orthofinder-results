"""Streamlit interface for generic OrthoFinder resource interrogation."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Mapping

import streamlit as st

from orthofinder_interrogation_app.all_results_page import render_all_distance_results
from orthofinder_interrogation_app.benchmark_page import render_dispersion_benchmarks
from orthofinder_interrogation_app.comparison_page import render_cluster_comparison
from orthofinder_interrogation_app.coverage_page import render_selection_coverage_tree
from orthofinder_interrogation_app.distance_data import default_cache_directory
from orthofinder_interrogation_app.documentation_page import (
    render_glossary_page,
    render_methods_page,
    render_page_guidance,
)
from orthofinder_interrogation_app.evolutionary_page import (
    _comparison_keys,
    _store_active_group,
    _store_comparison_keys,
    render_evolutionary_views,
)
from orthofinder_interrogation_app.exports import render_table_downloads
from orthofinder_interrogation_app.focus_page import render_focus_clusters
from orthofinder_interrogation_app.launcher import (
    CACHE_ENVIRONMENT_VARIABLE,
    EXPECTED_TAXA_ENVIRONMENT_VARIABLE,
    FOCUS_ENVIRONMENT_VARIABLE,
    LOG_ENVIRONMENT_VARIABLE,
    RESOURCE_ENVIRONMENT_VARIABLE,
    TAXONOMY_ENVIRONMENT_VARIABLE,
)
from orthofinder_interrogation_app.models import GroupKey, GroupSearchFilters
from orthofinder_interrogation_app.protein_page import render_protein_search
from orthofinder_interrogation_app.queries import (
    MAX_GROUP_MEMBER_ROWS,
    OrthoFinderQueryService,
)
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_interrogation_app.taxonomy_page import render_taxonomy_search
from orthofinder_results import __version__
from orthofinder_results.errors import OrthoFinderResultsError
from orthofinder_results.io_utils import configure_logging

_LOGGER = logging.getLogger("orthofinder_interrogation_app.app")
_PAGES = (
    "Overview",
    "Focus protein clusters",
    "Find a protein",
    "Find groups",
    "Cluster explorer",
    "Compare clusters",
    "Dispersion benchmarks",
    "All distance results",
    "Taxonomic search",
    "Selection coverage tree",
    "Offline report",
    "Methods",
    "Glossary",
    "Help",
)
_PAGE_LABELS = {
    "Overview": "Summary",
    "Focus protein clusters": "Focus protein clusters",
    "Find a protein": "Find a gene / protein",
    "Find groups": "Find groups",
    "Cluster explorer": "Explore one cluster",
    "Compare clusters": "Compare clusters",
    "Dispersion benchmarks": "Calibrated dispersion",
    "All distance results": "All distance results",
    "Taxonomic search": "Taxonomic search",
    "Selection coverage tree": "Selection coverage tree",
    "Offline report": "Download report",
    "Methods": "Methods & provenance",
    "Glossary": "Glossary",
    "Help": "Using the app",
}
_INCLUDE_LABEL_TO_MODE = {
    "At least one selected species": "ANY",
    "Every selected species": "ALL",
    "Exactly the selected species set": "EXACT_SET",
}
_DISTANCE_LABEL_TO_MODE = {
    "Any stored-distance state": "ANY",
    "Stored distances available": "CALCULATED",
    "No stored distances": "NOT_CALCULATED",
}
_SORT_LABEL_TO_MODE = {
    "Group identifier": "GROUP_ID_ASC",
    "Largest groups": "MEMBER_COUNT_DESC",
    "Broadest species representation": "SPECIES_COUNT_DESC",
    "Smallest stored mean distance": "MEAN_DISTANCE_ASC",
    "Largest stored mean distance": "MEAN_DISTANCE_DESC",
    "Smallest stored distance spread": "DISTANCE_SD_ASC",
}
_AUTHORITY_COLUMN_HELP = {
    "Group system": (
        "HOG means a hierarchical orthogroup at one species-tree level; "
        "LEGACY_ORTHOGROUP is the flat Orthogroups.tsv collection."
    ),
    "Species-tree level": (
        "The OrthoFinder node at which HOG membership is defined. ROOT is used for "
        "a flat legacy orthogroup collection."
    ),
    "Number of groups": "Number of group records in this exact system and hierarchy level.",
    "Smallest group (proteins)": "Fewest protein members in any group in this collection.",
    "Typical group (proteins)": "Median protein-member count across this collection.",
    "Largest group (proteins)": "Greatest protein-member count in this collection.",
    "Average species per group": "Mean number of represented species per group.",
    "Most species in one group": "Largest represented-species count in one group.",
}
_GROUP_COLUMN_HELP = {
    "Group system": _AUTHORITY_COLUMN_HELP["Group system"],
    "Species-tree level": _AUTHORITY_COLUMN_HELP["Species-tree level"],
    "Group ID": "Exact OrthoFinder group identifier within the stated system and level.",
    "Proteins in group": "Complete membership count for this group record.",
    "Species represented": "Number of sampled species contributing at least one protein.",
    "Highest copies in one species": "Largest protein copy count observed for one species.",
    "Average copies per represented species": (
        "Protein count divided by represented species; species absent from the group are "
        "not included in this mean."
    ),
    "Stored distance coverage": (
        "Whether a complete or deterministic bounded pair-distance matrix is already "
        "stored for this group. Not calculated does not mean zero distance."
    ),
    "Stored mean pair distance": "Mean across the stored pair-distance matrix, when present.",
    "Stored distance spread (SD)": (
        "Population standard deviation across stored pair distances; a smaller value "
        "means the displayed proteins have more similar pairwise distances."
    ),
}
_SPECIES_COLUMN_HELP = {
    "Species": "Exact species label recorded by OrthoFinder for this run.",
    "Proteins from species": "Number of group members contributed by this species.",
    "Share of group": "Fraction of all proteins in this group contributed by this species.",
}
_MEMBER_COLUMN_HELP = {
    "Species": "Exact OrthoFinder species label for this protein.",
    "Protein ID": "Exact protein or sequence identifier stored in the group membership.",
    "Parent legacy orthogroup": "Linked flat orthogroup identifier when OrthoFinder supplied it.",
    "Gene-tree parent clade": "OrthoFinder gene-tree clade used to define the HOG when present.",
}
_DISTANCE_STATUS_LABELS = {
    "COMPLETE": "Complete distances stored",
    "EXACT": "Complete distances stored",
    "DETERMINISTIC_MEMBER_SAMPLE": "Sampled distances stored",
    "PORTABLE_TREE_LAZY": "Calculated on demand from portable tree",
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
        "Data resource",
        value=default_resource,
        help=(
            "Choose a completed orthofinder-results directory or its "
            "duckdb/orthofinder_results.duckdb file. The app opens it read-only."
        ),
        placeholder="Completed resource directory or DuckDB file",
    )
    page_name = st.sidebar.radio(
        "Go to",
        _PAGES,
        key="app_page",
        format_func=lambda value: _PAGE_LABELS[value],
        help="Start at Summary, then move from finding groups to exploring and comparing them.",
    )
    default_cache = os.environ.get(
        CACHE_ENVIRONMENT_VARIABLE,
        str(default_cache_directory()),
    )
    with st.sidebar.expander("Advanced: cache, taxonomy & logs"):
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
        expected_path_text = st.text_input(
            "Expected-taxon TSV",
            value=os.environ.get(EXPECTED_TAXA_ENVIRONMENT_VARIABLE, ""),
            help="Optional reviewed expected universe for Selection coverage tree.",
        )
        focus_path_text = st.text_input(
            "Focus-protein TSV",
            value=os.environ.get(FOCUS_ENVIRONMENT_VARIABLE, ""),
            help=(
                "Optional replacement for the packaged E3 seed authority. Plain and "
                "gzip-compressed TSV are accepted by the launcher."
            ),
        )
    st.sidebar.caption(f"Viewer version {__version__}")
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
        elif page_name == "Focus protein clusters":
            render_focus_clusters(
                resource=resource,
                service=service,
                cache_dir=cache_dir,
                focus_path_text=focus_path_text,
            )
        elif page_name == "Find a protein":
            render_protein_search(
                resource=resource,
                service=service,
                cache_dir=cache_dir,
            )
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
        elif page_name == "Dispersion benchmarks":
            render_dispersion_benchmarks(service=service)
        elif page_name == "All distance results":
            render_all_distance_results(service=service)
        elif page_name == "Taxonomic search":
            render_taxonomy_search(
                service=service,
                taxonomy_path_text=taxonomy_path_text,
            )
        elif page_name == "Selection coverage tree":
            render_selection_coverage_tree(
                service=service,
                taxonomy_path_text=taxonomy_path_text,
                expected_path_text=expected_path_text,
                focus_path_text=focus_path_text,
            )
        elif page_name == "Offline report":
            _render_offline_report(resource=resource)
        elif page_name == "Methods":
            render_methods_page(resource=resource)
        elif page_name == "Glossary":
            render_glossary_page()
        else:
            _render_help()
    except OrthoFinderResultsError as error:
        _LOGGER.exception("Application request failed")
        st.error(f"The resource could not be opened or queried: {error}")


def _inject_style() -> None:
    """Apply readable widths, metric styling and accessible navigation text."""

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
        [data-testid="stSidebar"] [data-testid="stWidgetLabel"] p {
          font-size: 1rem !important;
          font-weight: 650 !important;
        }
        [data-testid="stSidebar"] div[role="radiogroup"] label p {
          font-size: 1.075rem !important;
          line-height: 1.45rem !important;
        }
        [data-testid="stSidebar"] div[role="radiogroup"] label {
          padding-top: 0.12rem;
          padding-bottom: 0.12rem;
        }
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
        columns[0].metric(
            "Run",
            resource.run_id,
            help="The immutable identifier for this completed OrthoFinder analysis.",
        )
        columns[1].metric(
            "OrthoFinder version",
            resource.orthofinder_version,
            help="The OrthoFinder version recorded when this analysis was created.",
        )
        columns[2].metric(
            "Input layout",
            resource.adapter_name,
            help="The version-specific parser used to read the OrthoFinder directory.",
        )
        columns[3].metric(
            "Data schema",
            resource.schema_version,
            help="The version of the portable orthofinder-results data contract.",
        )
        columns[4].metric(
            "Built by package",
            resource.resource_package_version,
            help="The orthofinder-results version that built this immutable resource.",
        )
        st.code(str(resource.database_path), language=None)
        st.caption(
            "Resource package is the version that built this database; application package "
            "in the sidebar is the currently installed viewer."
        )


def _render_overview(*, service: OrthoFinderQueryService) -> None:
    """Render a guided landing page before exposing technical detail."""

    st.header("Dataset summary")
    render_page_guidance(key="overview")
    st.write(
        "This page is the starting point. It describes what is in the completed run, "
        "what can be analysed now, and which page to use for each biological question."
    )
    counts = service.overview_counts()
    authorities = tuple(
        _display_authority_row(row=row) for row in service.overview_authorities()
    )
    columns = st.columns(5)
    columns[0].metric(
        "Group records",
        f"{counts['group_count']:,}",
        help=(
            "All retained HOG records across every hierarchy level, plus legacy "
            "orthogroups when present. This is not a count of unique root families."
        ),
    )
    columns[1].metric(
        "Species in this run",
        f"{counts['species_count']:,}",
        help="Exact species labels present in this OrthoFinder analysis.",
    )
    columns[2].metric(
        "Group collections",
        f"{len(authorities):,}",
        help="Distinct combinations of group system and species-tree hierarchy level.",
    )
    columns[3].metric(
        "Groups with stored distances",
        f"{counts['distance_group_count']:,}",
        help="Groups that already contain at least one persisted pairwise-distance result.",
    )
    columns[4].metric(
        "Portable gene trees",
        f"{counts['portable_tree_count']:,}",
        help=(
            "Checksum-bound gene trees available for bounded on-demand distance "
            "calculation. Schema-2 resources normally show zero here."
        ),
    )
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
    st.subheader("Choose your question")
    st.markdown("#### Which clusters contain our E3 proteins of interest?")
    st.write(
        "Use the packaged, versioned E3 seed authority by default, or substitute a reviewed "
        "custom protein list. Matching clusters retain the full distance and visual suite."
    )
    st.button(
        "Find E3 focus clusters",
        key="summary_focus_clusters",
        type="primary",
        on_click=_navigate_to,
        kwargs={"page": "Focus protein clusters"},
    )
    actions = st.columns(3)
    actions[0].markdown("#### Which cluster contains my protein?")
    actions[0].write(
        "Find an exact protein or internal ID, then inspect its cluster and distances."
    )
    actions[0].button(
        "Find a gene or protein",
        key="summary_find_protein",
        on_click=_navigate_to,
        kwargs={"page": "Find a protein"},
    )
    actions[1].markdown("#### Which groups contain my species?")
    actions[1].write(
        "Require one, every, or exactly a set of sampled species, and reject unwanted species."
    )
    actions[1].button(
        "Find groups",
        key="summary_find_groups",
        on_click=_navigate_to,
        kwargs={"page": "Find groups"},
    )
    actions[2].markdown("#### How compact or dispersed is one group?")
    actions[2].write(
        "Inspect exact distances, distributions, PCoA, a gene-tree phylogram and networks."
    )
    actions[2].button(
        "Explore one cluster",
        key="summary_explore",
        on_click=_navigate_to,
        kwargs={"page": "Cluster explorer"},
    )
    more_actions = st.columns(3)
    more_actions[0].markdown("#### Which groups focus on a lineage?")
    more_actions[0].write(
        "Search reviewed taxonomic descendants using contains, enrichment or exclusivity."
    )
    more_actions[0].button(
        "Search taxonomy",
        key="summary_taxonomy",
        on_click=_navigate_to,
        kwargs={"page": "Taxonomic search"},
    )
    more_actions[1].markdown("#### How do candidate groups differ?")
    more_actions[1].write(
        "Collect 2–12 groups and compare distance spread, compactness and PCoA shape."
    )
    more_actions[1].button(
        "Compare clusters",
        key="summary_compare",
        on_click=_navigate_to,
        kwargs={"page": "Compare clusters"},
    )
    more_actions[2].markdown("#### Which clusters have distance results?")
    more_actions[2].write(
        "Select columns and export every persisted cluster-distance summary in one table."
    )
    more_actions[2].button(
        "All distance results",
        key="summary_all_results",
        on_click=_navigate_to,
        kwargs={"page": "All distance results"},
    )
    coverage_action = st.columns(3)
    coverage_action[0].markdown("#### How does a combined taxonomy selection behave?")
    coverage_action[0].write(
        "Build a reviewed selection/coverage tree and audit exact, include, only and "
        "exclusion predicates across E3-focus clusters by default."
    )
    coverage_action[0].button(
        "Build selection coverage tree",
        key="summary_coverage_tree",
        on_click=_navigate_to,
        kwargs={"page": "Selection coverage tree"},
    )
    coverage_action[1].markdown("#### Are E3 clusters unusually dispersed?")
    coverage_action[1].write(
        "Test E3s and individual E3 subtypes against housekeeping, R/NLR and "
        "structure-matched non-focus cluster backgrounds."
    )
    coverage_action[1].button(
        "Open calibrated dispersion",
        key="summary_dispersion_benchmarks",
        on_click=_navigate_to,
        kwargs={"page": "Dispersion benchmarks"},
        disabled=not service.has_relation(relation="benchmark_contrasts"),
        help=(
            "Available after the dispersion-benchmark cluster workflow has built "
            "the required background relations."
        ),
    )

    with st.expander("Detailed group collections and hierarchy levels", expanded=False):
        st.write(
            "A HOG can appear at several species-tree levels. Rows below describe each "
            "retained group collection separately; they should not be summed as unique "
            "biological families."
        )
        st.dataframe(
            authorities,
            width="stretch",
            hide_index=True,
            column_config=_column_config(descriptions=_AUTHORITY_COLUMN_HELP),
        )
        render_table_downloads(
            records=authorities,
            file_stem="orthofinder_group_authority_summary",
            key="overview_group_authority_download",
            tsv_label="Download group authority summary as TSV",
            excel_label="Download group authority summary as formatted Excel",
            column_definitions=_AUTHORITY_COLUMN_HELP,
            workbook_title="OrthoFinder group authority summary",
        )
        st.caption(
            f"The database also contains {counts['group_species_statistic_count']:,} "
            "group-by-species copy-count rows used by species and taxonomy filters."
        )

    with st.expander("Scope of this standalone app and the wider E3 workflow"):
        st.markdown(
            """
            **Available here:** generic OrthoFinder group discovery, a replaceable protein-focus
            authority (the packaged E3 seed list is the default), species and copy-number
            profiles, reviewed taxonomy selection/coverage trees, bounded tree distances,
            compactness and dispersion views, gene-tree inspection, matched-background
            cluster statistics, and comparison within one run.

            **Planned generic extensions:** explicit nested-HOG and split/merge interrogation,
            followed by comparisons between independently versioned OrthoFinder runs.

            **Kept in the separate E3 application:** downstream E3-ligase prioritisation,
            expression, experimental evidence, structures, ligandable-pocket conservation
            and chemistry starting points. The two applications can later exchange versioned
            group links
            without putting E3-specific assumptions into this reusable backend.
            """
        )


def _render_group_search(*, service: OrthoFinderQueryService) -> None:
    """Render progressive exact group filters and reusable selected-group actions."""

    st.header("Find groups")
    render_page_guidance(key="group_search")
    st.caption(
        "Start with species or an identifier; open Advanced filters only when needed. "
        "Required species and rejected species use the exact labels in this dataset."
    )
    with st.expander("How species inclusion works", expanded=False):
        st.markdown(
            """
            - **At least one selected species** keeps a group containing any selected label.
            - **Every selected species** requires all selected labels but allows other species.
            - **Exactly the selected species set** allows no additional sampled species.
            - **Reject groups containing** removes a group if any rejected label is present.
            """
        )
    species = service.list_species()
    group_types = service.list_group_types()
    with st.form("group_filters"):
        first_row = st.columns(4)
        group_type = first_row[0].selectbox(
            "Group system",
            ("Any", *group_types),
            help=_AUTHORITY_COLUMN_HELP["Group system"],
        )
        nodes = service.list_hierarchy_nodes(
            group_type="" if group_type == "Any" else group_type
        )
        node_labels = tuple("ROOT" if not node else node for node in nodes)
        hierarchy_label = first_row[1].selectbox(
            "Species-tree level",
            ("Any", *node_labels),
            help=_AUTHORITY_COLUMN_HELP["Species-tree level"],
        )
        group_text = first_row[2].text_input(
            "Group ID contains",
            help="Literal, case-insensitive text within the OrthoFinder group identifier.",
        )
        member_text = first_row[3].text_input(
            "Protein ID contains",
            help="Literal, case-insensitive text within any member identifier.",
        )
        species_row = st.columns((3, 2, 3))
        included_species = tuple(
            species_row[0].multiselect(
                "Required species",
                options=species,
                help="Select sampled species that a matching group must contain.",
            )
        )
        include_label = species_row[1].selectbox(
            "Inclusion rule",
            tuple(_INCLUDE_LABEL_TO_MODE),
            index=1,
            help="Controls how the required-species selection is interpreted.",
        )
        excluded_species = tuple(
            species_row[2].multiselect(
                "Reject groups containing",
                options=species,
                help="Any occurrence of a selected species excludes the whole group.",
            )
        )
        with st.expander("Advanced size, distance, sorting and paging filters"):
            size_row = st.columns(4)
            minimum_members = int(
                size_row[0].number_input(
                    "Minimum proteins in group",
                    min_value=1,
                    value=1,
                    help="Complete group-member count, not the displayed distance sample.",
                )
            )
            maximum_members_enabled = size_row[1].checkbox("Set maximum group size")
            maximum_members = int(
                size_row[1].number_input(
                    "Maximum proteins in group",
                    min_value=1,
                    value=max(minimum_members, 1000),
                    disabled=not maximum_members_enabled,
                )
            )
            minimum_species = int(
                size_row[2].number_input(
                    "Minimum species represented",
                    min_value=1,
                    value=1,
                    help="Species must contribute at least one protein to count as represented.",
                )
            )
            maximum_species_enabled = size_row[3].checkbox("Set maximum species breadth")
            maximum_species = int(
                size_row[3].number_input(
                    "Maximum species represented",
                    min_value=1,
                    value=max(minimum_species, len(species)),
                    disabled=not maximum_species_enabled,
                )
            )
            distance_row = st.columns(4)
            distance_label = distance_row[0].selectbox(
                "Stored distance availability",
                tuple(_DISTANCE_LABEL_TO_MODE),
                help=(
                    "Filters only matrices already stored in this resource. Schema-3 "
                    "on-demand capability is not precomputed for every group."
                ),
            )
            mean_filter_enabled = distance_row[1].checkbox("Limit stored mean distance")
            maximum_mean_distance = float(
                distance_row[1].number_input(
                    "Maximum stored mean distance",
                    min_value=0.0,
                    value=1.0,
                    disabled=not mean_filter_enabled,
                )
            )
            sd_filter_enabled = distance_row[2].checkbox("Limit stored distance spread")
            maximum_distance_sd = float(
                distance_row[2].number_input(
                    "Maximum stored distance SD",
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
    displayed_results = tuple(_display_group_row(row=row) for row in result.rows)
    st.dataframe(
        displayed_results,
        width="stretch",
        hide_index=True,
        column_config=_column_config(descriptions=_GROUP_COLUMN_HELP),
    )
    render_table_downloads(
        records=displayed_results,
        file_stem="orthofinder_group_search",
        key="group_search_download",
        tsv_label="Download this result page as TSV",
        excel_label="Download this result page as formatted Excel",
        column_definitions=_GROUP_COLUMN_HELP,
        workbook_title="OrthoFinder group search",
    )
    labels_to_keys = {
        _row_group_key(row=row).display_label(): _row_group_key(row=row)
        for row in result.rows
    }
    selected_label = st.selectbox("Selected matching group", tuple(labels_to_keys))
    selected_key = labels_to_keys[selected_label]
    actions = st.columns(2)
    actions[0].button(
        "Explore selected cluster",
        type="primary",
        on_click=_open_group_in_explorer,
        kwargs={"key": selected_key},
    )
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
    columns[0].metric(
        "Proteins in group",
        f"{group['member_count']:,}",
        help=_GROUP_COLUMN_HELP["Proteins in group"],
    )
    columns[1].metric(
        "Species represented",
        f"{group['species_count']:,}",
        help=_GROUP_COLUMN_HELP["Species represented"],
    )
    columns[2].metric(
        "Highest copies in one species",
        f"{group['max_copies_per_species']:,}",
        help=_GROUP_COLUMN_HELP["Highest copies in one species"],
    )
    columns[3].metric(
        "Average copies per represented species",
        _format_number(group["mean_copies_per_species"]),
        help=_GROUP_COLUMN_HELP["Average copies per represented species"],
    )
    distance_count = group.get("distance_pair_count")
    if distance_count:
        distance_columns = st.columns(4)
        distance_columns[0].metric(
            "Stored mean pair distance",
            _format_number(group["mean_distance"]),
            help=_GROUP_COLUMN_HELP["Stored mean pair distance"],
        )
        distance_columns[1].metric(
            "Stored distance spread (SD)",
            _format_number(group["population_stddev_distance"]),
            help=_GROUP_COLUMN_HELP["Stored distance spread (SD)"],
        )
        distance_columns[2].metric(
            "Stored median pair distance",
            _format_number(group["median_distance"]),
            help="Middle pair distance in the stored matrix.",
        )
        distance_columns[3].metric(
            "Protein pairs compared",
            f"{distance_count:,}",
            help="Number of pairwise comparisons in the stored matrix.",
        )
        with st.expander("Technical distance provenance"):
            st.code(
                f"method={group['distance_method']}\n"
                f"status={group['computation_status']}\n"
                f"sampled_members={group['sampled_member_count']}\n"
                f"total_members={group['total_member_count']}",
                language=None,
            )
    elif service.has_relation(relation="tree_payloads"):
        st.info(
            "No distance pairs were precomputed for this group. Explore one cluster can try a "
            "bounded on-demand calculation from its portable gene tree."
        )
    else:
        st.info(
            "No distances were calculated for this schema-2 group. A schema-3 resource "
            "rebuild is required for on-demand calculation."
        )
    species_rows = service.get_group_species(key=key)
    member_rows = service.get_group_members(key=key)
    species_tab, member_tab = st.tabs(("Species copy counts", "Protein members"))
    with species_tab:
        displayed_species = tuple(_display_species_row(row=row) for row in species_rows)
        st.dataframe(
            displayed_species,
            width="stretch",
            hide_index=True,
            column_config=_column_config(descriptions=_SPECIES_COLUMN_HELP),
        )
        render_table_downloads(
            records=displayed_species,
            file_stem=f"{key.group_id}_species_copy_counts",
            key=f"group_species_{key.display_label()}",
            tsv_label="Download species copy counts as TSV",
            excel_label="Download species copy counts as formatted Excel",
            column_definitions=_SPECIES_COLUMN_HELP,
            workbook_title=f"Species copy counts: {key.group_id}",
        )
    with member_tab:
        displayed_members = tuple(_display_member_row(row=row) for row in member_rows)
        st.dataframe(
            displayed_members,
            width="stretch",
            hide_index=True,
            column_config=_column_config(descriptions=_MEMBER_COLUMN_HELP),
        )
        if len(member_rows) == group["member_count"]:
            render_table_downloads(
                records=displayed_members,
                file_stem=f"{key.group_id}_members",
                key=f"group_members_{key.display_label()}",
                tsv_label="Download complete membership as TSV",
                excel_label="Download complete membership as formatted Excel",
                column_definitions=_MEMBER_COLUMN_HELP,
                workbook_title=f"Group membership: {key.group_id}",
            )
        else:
            st.warning(
                f"The browser preview is limited to {MAX_GROUP_MEMBER_ROWS:,} members; "
                f"this group contains {group['member_count']:,}. Use DuckDB/Parquet for "
                "a complete export."
            )


def _render_offline_report(*, resource: Any) -> None:
    """Expose the immutable offline report as a download."""

    st.header("Download the offline report")
    render_page_guidance(key="offline_report")
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
    """Render task-oriented help for navigating the application."""

    st.header("Using the app")
    st.write(
        "Look for the **?** beside controls and column headings for help in context. "
        "Every analytical page and figure also has a closed interpretation panel. Use "
        "**Methods & provenance** for the workflow and **Glossary** for exact terms."
    )
    with st.expander("Getting started", expanded=True):
        st.markdown(
            """
            1. Open **Summary** to understand the run and its distance capability.
            2. Use **Find a gene / protein** to locate a stored identifier, inspect every
               cluster membership and highlight it throughout a selected cluster.
            3. Use **Find groups** for exact species, identifier, size and copy-number questions.
            4. Send an interesting result to **Explore one cluster** for distances and trees.
            5. Add 2–12 groups to **Compare clusters** when the same distance method and
               sampling scope are scientifically comparable.
            6. Use **All distance results** to select columns and export every cluster with a
               persisted distance summary as formatted Excel or TSV.
            7. Open **Methods & provenance** for the complete analytical route and the exact
               capabilities of the resource currently open.
            8. Search **Glossary** whenever a biological, statistical, tree or data-contract
               term is unfamiliar.

            **Taxonomic search** is a separate, stricter workflow because descendant claims
            require a reviewed species-to-taxonomy mapping. The packaged Results_Feb26
            authority is used only when all 60 exact labels match; every other dataset gets
            its own downloadable review template.
            """
        )
    with st.expander("Groups, HOGs and species-tree levels"):
        st.markdown(
            """
            - A **legacy orthogroup** is a flat group from `Orthogroups.tsv`.
            - A **hierarchical orthogroup (HOG)** is defined at a named node of the species
              tree. The same broad family can therefore have related records at several levels.
            - **Group records** counts all retained group-and-level records. It is not the
              number of unique root families.
            - **Species represented** counts species with at least one protein in a group.
            - **Average copies per represented species** excludes species absent from that group.
            """
        )
    with st.expander("Finding one gene or protein"):
        st.markdown(
            """
            The dedicated search accepts canonical membership IDs and, when available,
            OrthoFinder internal IDs from `SequenceIDs.txt`. Exact search is case-sensitive;
            the optional contains search is literal and case-insensitive. Descriptive gene
            symbols or annotations are searchable only when they occur in the stored identifier.

            Selecting a result opens the complete cluster visual suite and highlights the
            canonical protein. A separate nearest-to-farthest table reports every displayed
            distance involving it. If a schema-3 calculation must sample a large group, the
            requested protein is forcibly retained. A schema-2 matrix cannot be expanded when
            the protein was not included in its original persisted sample.
            """
        )
    with st.expander("Focus proteins and the packaged E3 authority"):
        st.markdown(
            """
            **Focus protein clusters** matches a versioned protein list to canonical IDs,
            OrthoFinder internal IDs and controlled UniProt accession/entry aliases. The packaged
            E3 seed-evidence table is selected by default for this project, but `--focus-proteins`
            or the page upload accepts a replacement authority without code changes.

            A match means that a cluster contains at least one configured focus protein. It does
            not automatically assign E3 function to every other cluster member. Open a match to
            inspect its exact distances, dispersion plots, PCoA, gene-tree phylogram and networks.
            """
        )
    with st.expander("Required species, exact sets and rejected species"):
        st.markdown(
            """
            - **At least one selected species** retains groups containing any selected label.
            - **Every selected species** requires all selected labels and allows other species.
            - **Exactly the selected species set** requires precisely those sampled species.
            - **Reject groups containing** removes a group if any rejected exact label occurs.

            These filters use literal labels from this run; spelling alone is never treated as
            evidence of taxonomic ancestry.
            """
        )
    with st.expander("Distances, compactness and sampling"):
        st.markdown(
            """
            - **Pair distance** is the stated distance between two displayed proteins. For
              tree-based analyses it is the sum of branch lengths joining the two leaves.
            - **Mean** and **median** describe the centre of all displayed pair distances.
              Smaller values generally indicate a more compact displayed group under that method.
            - **Distance spread (population SD)** describes how variable those pair distances
              are. A group may have a small mean but a broad spread, so inspect both.
            - **Proteins in full group** is the complete membership count. **Proteins analysed**
              is the exact or deterministic bounded sample used for the displayed matrix.
            - Missing distance values mean **not calculated or unavailable**, never zero.
            - The **sample medoid** is the displayed protein with the lowest mean distance to
              other displayed proteins. It is not an ancestor or guaranteed full-group centre.
            """
        )
    with st.expander("Networks, PCoA, phylograms and exact matrices"):
        st.markdown(
            """
            - The **interactive network** keeps a few nearest-neighbour edges. Its layout helps
              exploration, but screen spacing is not a distance scale and components are not
              automatically OrthoFinder splits.
            - **PCoA** places proteins in two or three dimensions to approximate all pair
              distances. Read it with its stress, distance agreement and Shepard diagnostic;
              apparent arms or gaps are hypotheses, not automatic subfamilies.
            - The **branch-length phylogram** preserves horizontal gene-tree branch length;
              vertical spacing is only layout.
            - The **distance heatmap** and protein-pair table retain the exact displayed
              quantitative values and should be used to confirm patterns seen in layouts.
            """
        )
    with st.expander("Reviewed taxonomy and exclusivity claims"):
        st.markdown(
            """
            Taxonomic ancestry is accepted only from a separately reviewed mapping generated
            for the current dataset. Unmapped, pending and ambiguous labels remain explicit and
            never become inferred descendants or reviewed outsiders.

            **Contains** requires target descendants. **Enriched** tests species presence and
            reports multiple-testing-corrected q-values. **Exclusive within sampled analysis**
            means exclusive only among species included in this OrthoFinder run—not universal
            biological absence. **Near-exclusive** permits only the configured sampled outsiders
            and unresolved labels, all of which remain listed.
            """
        )
    with st.expander("Selection coverage tree"):
        st.markdown(
            """
            This page draws a **reviewed taxonomy coverage tree**, not the OrthoFinder species
            tree or a HOG gene tree. Required exact, include-clade, only-in-clade, exact-exclude
            and clade-exclude selectors all compose with AND. Several only-in clades use their
            terminal intersection, and any unmapped group member makes that restriction fail
            closed.

            Selection state and dataset coverage are deliberately separate. **Expected no data**
            means a reviewed expected taxon is not represented in this imported dataset; it is
            never a biological absence claim. The complete download pairs Newick with a style
            table and includes TSV audits, SVG, PDF, JSON provenance and SHA-256 checksums.
            """
        )
    with st.expander("Resource schemas, caches and provenance"):
        st.markdown(
            """
            Schemas 3 and 4 store checksum-bound portable gene trees rather than billions of pair
            matrices. The app can calculate one bounded matrix when requested and caches it in
            a user sidecar outside the immutable resource. Schema 4 also adds the calibrated-
            dispersion profile, control, contrast, individual-test and classification relations.
            Schema 2 remains readable but normally exposes distances only for its original pilot
            groups. Method, status, source and cache details remain available under **Technical
            analysis details**, **Methods & provenance** and relevant downloadable tables.
            """
        )


def _display_authority_row(*, row: dict[str, Any]) -> dict[str, Any]:
    """Return plain-language headings for one overview authority row."""

    return {
        "Group system": row["group_type"],
        "Species-tree level": row["hierarchy_node"] or "ROOT",
        "Number of groups": row["group_count"],
        "Smallest group (proteins)": row["minimum_members"],
        "Typical group (proteins)": row["median_members"],
        "Largest group (proteins)": row["maximum_members"],
        "Average species per group": row["mean_species"],
        "Most species in one group": row["maximum_species"],
    }


def _display_group_row(*, row: dict[str, Any]) -> dict[str, Any]:
    """Return plain-language headings for a group-search row."""

    return {
        "Group system": row["group_type"],
        "Species-tree level": row["hierarchy_node"] or "ROOT",
        "Group ID": row["group_id"],
        "Proteins in group": row["member_count"],
        "Species represented": row["species_count"],
        "Highest copies in one species": row["max_copies_per_species"],
        "Average copies per represented species": row["mean_copies_per_species"],
        "Stored distance coverage": _distance_status_label(
            value=row.get("computation_status")
        ),
        "Stored mean pair distance": row["mean_distance"],
        "Stored distance spread (SD)": row["population_stddev_distance"],
    }


def _display_species_row(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return readable per-species copy-count fields for browser display."""

    return {
        "Species": row["species_label"],
        "Proteins from species": row["species_member_count"],
        "Share of group": row["member_fraction"],
    }


def _display_member_row(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return readable membership fields while hiding repeated composite keys."""

    return {
        "Species": row["species_label"],
        "Protein ID": row["member_id"],
        "Parent legacy orthogroup": row.get("legacy_orthogroup_id") or "Unavailable",
        "Gene-tree parent clade": row.get("gene_tree_parent_clade") or "Unavailable",
    }


def _distance_status_label(*, value: object) -> str:
    """Translate a stored computation code without discarding unknown provenance."""

    if value is None or str(value).strip() == "":
        return "Not calculated"
    code = str(value)
    return _DISTANCE_STATUS_LABELS.get(code, code.replace("_", " ").title())


def _column_config(*, descriptions: Mapping[str, str]) -> dict[str, Any]:
    """Return Streamlit column definitions with hoverable help descriptions."""

    return {
        label: st.column_config.Column(label=label, help=description)
        for label, description in descriptions.items()
    }


def _navigate_to(*, page: str) -> None:
    """Set one validated page from a Streamlit widget callback.

    Args:
        page: Exact internal page identifier.

    Raises:
        ValueError: If the page is not registered by the application.
    """

    if page not in _PAGES:
        raise ValueError(f"Unknown application page: {page}")
    st.session_state["app_page"] = page


def _open_group_in_explorer(*, key: GroupKey) -> None:
    """Store one selected group and navigate from a widget callback.

    Args:
        key: Collision-safe selected group key.
    """

    _store_active_group(key=key)
    _navigate_to(page="Cluster explorer")


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

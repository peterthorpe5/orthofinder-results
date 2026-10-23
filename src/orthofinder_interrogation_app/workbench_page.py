"""OrthoFinder navigation and complete input handoff into protein research."""

from __future__ import annotations

import logging
import os
from functools import partial
from pathlib import Path

import streamlit as st

from orthofinder_results.errors import OrthoFinderResultsError

from .evolutionary_page import ACTIVE_GROUP_STATE
from .exports import render_table_downloads
from .launcher import TERMINAL_MOTIF_ENVIRONMENT_VARIABLE
from .models import GroupKey
from .queries import OrthoFinderQueryService
from .research.common import ResearchConfig, ResearchError, load_config
from .research.ui import queue_action, render_research_workbench, split_identifiers
from .workbench_adapter import select_group_inputs, selected_group_fasta

LOGGER = logging.getLogger(__name__)


def open_group_workbench(*, key: GroupKey) -> None:
    """Open research with the exact current group retained in serialisable state.

    Args:
        key: Exact run, group authority, hierarchy node and group ID.
    """
    st.session_state[ACTIVE_GROUP_STATE] = {
        name: getattr(key, name) for name in ("run_id", "group_type", "hierarchy_node", "group_id")
    }
    st.session_state["app_page"] = "Protein research workbench"
    st.session_state["research_section"] = "OrthoFinder groups"


def render_workbench_page(*, service: OrthoFinderQueryService) -> None:
    """Render the integrated workbench with a run-aware local input section.

    Args:
        service: Read-only service for the current OrthoFinder run.
    """
    try:
        config = load_config()
        root = service.resource.resource_path.resolve()
        if config.workspace.resolve().is_relative_to(root):
            raise ResearchError("Choose a research workspace outside the completed resource")
        render_research_workbench(input_selector=partial(render_group_inputs, service=service))
    except (ResearchError, OrthoFinderResultsError, OSError) as error:
        LOGGER.exception("Could not open integrated protein research")
        st.error(str(error))


def render_group_inputs(*, config: ResearchConfig, service: OrthoFinderQueryService) -> None:
    """Select full local group memberships for exact global mapping or FASTA analysis.

    Args:
        config: Shared research workspace and tool profile.
        service: Current immutable OrthoFinder run.
    """
    st.subheader("Start from OrthoFinder groups")
    st.write(
        "Select exact groups or the largest groups at one hierarchy level. The complete "
        "member list is retained with its run, species and group context. Start with one "
        "group, collect its proteins, then compare its sequences or find the corresponding "
        "global DeepClust families."
    )
    with st.expander("Aim, steps and outputs", expanded=True):
        st.markdown(
            "1. Choose the group authority and species-tree level.\n"
            "2. Enter exact group IDs or rank local groups by size, then prepare the inputs.\n"
            "3. Import their local sequences or map their IDs to the server's DeepClust "
            "database. Track the job and use its next-step button.\n\n"
            "**Outputs:** complete membership with group/species context, downloadable "
            "TSV/Excel tables and a queued family import or global mapping. Local size "
            "ranking describes this OrthoFinder run. It is separate from global cluster size."
        )
    active = st.session_state.get(ACTIVE_GROUP_STATE, {})
    if not isinstance(active, dict) or active.get("run_id") != service.resource.run_id:
        active = {}
    types = list(service.list_group_types())
    if not types:
        st.info("This resource contains no selectable groups.")
        return
    group_type = st.selectbox(
        "OrthoFinder group authority",
        options=types,
        index=types.index(active["group_type"]) if active.get("group_type") in types else 0,
        help="HOGs are defined at a species-tree node; flat legacy orthogroups have no node.",
    )
    nodes = list(service.list_hierarchy_nodes(group_type=group_type))
    if not nodes:
        st.info("No hierarchy levels are available for this group authority.")
        return
    node = st.selectbox(
        "Species-tree hierarchy level",
        options=nodes,
        index=nodes.index(active["hierarchy_node"]) if active.get("hierarchy_node") in nodes else 0,
        format_func=lambda value: value or "ROOT (flat orthogroups)",
        help="Choose the biological ancestor at which the HOG is defined. IDs at different "
        "levels describe different memberships, even if their text resembles one another.",
    )
    with st.form("orthofinder_research_inputs"):
        exact_text = st.text_area(
            "Exact local group IDs (optional; overrides size ranking)",
            value=active.get("group_id", "")
            if active.get("group_type") == group_type and active.get("hierarchy_node") == node
            else "",
            help="Paste IDs from the selected run, authority and hierarchy level. Separate "
            "with whitespace or commas. Exact IDs override the size-ranking selection.",
        )
        maximum = st.number_input(
            "Largest local groups",
            min_value=1,
            max_value=500,
            value=1,
            help="Number of complete local groups, ordered by membership size. Begin with "
            "one; increase after reviewing the results. Exact IDs override this ranking "
            "but must fit within this limit.",
        )
        species = st.multiselect(
            "Groups containing species (optional)",
            options=service.list_species(),
            help="Filter which groups qualify. All members of each qualifying group remain "
            "in the output, including proteins from other species.",
        )
        all_species = st.checkbox(
            "Require every selected OrthoFinder species",
            help="When selected, each group must contain all chosen species; otherwise "
            "any chosen species is sufficient. No effect when no species are selected.",
        )
        query = st.text_input(
            "Groups containing protein ID text (optional)",
            help="Match text in a member ID to find its groups. The output retains their "
            "complete memberships, rather than just the matching protein.",
        )
        prepared = st.form_submit_button("Prepare complete group inputs", type="primary")
    signature = (
        str(service.resource.database_path),
        service.resource.run_id,
        group_type,
        node,
        exact_text,
        maximum,
        tuple(species),
        all_species,
        query,
    )
    if prepared:
        try:
            frame = select_group_inputs(
                service=service,
                group_type=group_type,
                hierarchy_node=node,
                group_ids=split_identifiers(value=exact_text),
                maximum_groups=maximum,
                included_species=species,
                require_all_species=all_species,
                member_id_contains=query,
            )
        except (ResearchError, OrthoFinderResultsError) as error:
            st.error(str(error))
            return
        st.session_state["workbench_inputs"] = frame
        st.session_state["workbench_input_signature"] = signature
    if st.session_state.get("workbench_input_signature") != signature:
        st.info("Prepare the current selection to see its complete membership.")
        return
    frame = st.session_state["workbench_inputs"]
    if frame.empty:
        st.info("No groups match these filters.")
        return
    st.dataframe(frame, hide_index=True)
    render_table_downloads(
        records=frame.to_dict(orient="records"),
        file_stem="orthofinder_research_inputs",
        key="research_local_context",
        workbook_title="OrthoFinder protein research inputs",
    )
    identifiers = frame["input_id"].drop_duplicates().tolist()
    st.caption(f"{len(identifiers):,} distinct protein IDs; every selected group is complete.")
    mapping_ready = bool(config.sequence_map and config.index)
    if st.button(
        "Queue selected OrthoFinder proteins for global mapping",
        disabled=not mapping_ready,
        help="Resolve exact IDs and declared aliases using the server's DeepClust mapping "
        "and cluster index. This queues a local worker job; it does not download the database.",
    ):
        queue_action(
            config=config,
            operation="resolve",
            parameters={
                "identifiers": identifiers,
                "context": frame.to_dict(orient="records"),
            },
        )
    st.caption(
        "Exact global mapping requires the installed DeepClust sequence map and cluster index. "
        "For absent IDs, export these local sequences and use representative homology search."
    )
    if not mapping_ready:
        st.warning("Global mapping needs sequence_map and index in the research profile.")
    if st.button(
        "Import selected local protein sequences",
        help="Collect the complete selected proteins from embedded sequences or the "
        "configured sequence Parquet sidecar. Missing or ambiguous sequences stop the "
        "import. The completed family retains its OrthoFinder context.",
    ):
        try:
            text = os.environ.get(TERMINAL_MOTIF_ENVIRONMENT_VARIABLE, "").strip()
            fasta = selected_group_fasta(
                service=service,
                inputs=frame,
                sidecar=Path(text) if text else None,
            )
            queue_action(
                config=config,
                operation="import_fasta",
                fasta=fasta,
                parameters={
                    "source_label": f"OrthoFinder run {service.resource.run_id}",
                    "source_context": frame.to_dict(orient="records"),
                },
            )
            st.info("Track this job. When it completes, use Analyse this family to continue.")
        except (ResearchError, OrthoFinderResultsError) as error:
            st.error(str(error))

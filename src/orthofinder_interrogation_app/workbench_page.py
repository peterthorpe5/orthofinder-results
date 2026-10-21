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
        "member list is retained with its run, species and group context. Local size rank "
        "describes this OrthoFinder selection; E3 candidate scores belong to the E3 app."
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
    )
    with st.form("orthofinder_research_inputs"):
        exact_text = st.text_area(
            "Exact local group IDs (optional; overrides size ranking)",
            value=active.get("group_id", "")
            if active.get("group_type") == group_type and active.get("hierarchy_node") == node
            else "",
        )
        maximum = st.number_input("Largest local groups", min_value=1, max_value=500, value=200)
        species = st.multiselect(
            "Groups containing species (optional)", options=service.list_species()
        )
        all_species = st.checkbox("Require every selected OrthoFinder species")
        query = st.text_input("Groups containing protein ID text (optional)")
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
    if st.button("Queue selected OrthoFinder proteins for global mapping"):
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
    if st.button("Import selected local protein sequences"):
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
            st.info("After the worker finishes, select family.fasta under Family analysis.")
        except (ResearchError, OrthoFinderResultsError) as error:
            st.error(str(error))

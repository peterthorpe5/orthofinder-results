"""Streamlit workflow for default or user-replaceable focus protein clusters."""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import streamlit as st

from orthofinder_results.errors import InputValidationError, OrthoFinderResultsError

from .evolutionary_page import render_selected_group_visualisations
from .exports import render_table_downloads
from .focus import (
    FocusProteinAuthority,
    bundled_focus_path,
    focus_template,
    parse_focus_proteins,
    read_focus_proteins,
)
from .models import FocusClusterFilters, GroupKey
from .queries import OrthoFinderQueryService

_LOGGER = logging.getLogger("orthofinder_interrogation_app.focus_page")
_FOCUS_RESULT_STATE = "orthofinder_focus_cluster_result"
_FOCUS_SIGNATURE_STATE = "orthofinder_focus_cluster_signature"
_FOCUS_OPEN_STATE = "orthofinder_focus_open_group"
_RESULT_HELP = {
    "Group system": "HOG is hierarchical; LEGACY_ORTHOGROUP is the flat group set.",
    "Species-tree level": "Exact OrthoFinder hierarchy node; ROOT denotes a flat set.",
    "Group ID": "Exact OrthoFinder group identifier in the stated collection.",
    "Matched focus IDs": "Exact configured focus identifiers found in this cluster.",
    "Matched focus proteins": "Distinct cluster proteins matching at least one focus alias.",
    "Matching identifiers": "Configured identifiers matched through exact controlled aliases.",
    "Matched protein IDs": "Canonical protein identifiers stored by OrthoFinder.",
    "Matched species": "Exact resource species labels contributing matched proteins.",
    "Matched through": "Canonical, internal, UniProt accession or UniProt entry authority.",
    "Proteins in group": "Complete membership count for this group record.",
    "Species represented": "Sampled species contributing at least one group member.",
    "Stored distance status": "Preferred stored matrix state; missing is not zero distance.",
    "Stored mean distance": "Mean of stored pair distances when available.",
    "Stored distance spread (SD)": "Population SD across stored pair distances.",
}


def render_focus_clusters(
    *,
    resource: Any,
    service: OrthoFinderQueryService,
    cache_dir: Path,
    focus_path_text: str,
) -> None:
    """Render a default E3 or custom protein-authority cluster workflow.

    Args:
        resource: Validated immutable resource identity.
        service: Read-only query service for the resource.
        cache_dir: Persistent on-demand analysis sidecar.
        focus_path_text: Optional launcher-provided custom authority path.
    """

    st.header("Focus protein clusters")
    st.write(
        "Start with the packaged E3 seed-evidence authority, or replace it with a reviewed "
        "project protein list. Matching is exact and identifies clusters containing at least "
        "one configured protein; it does not assign E3 function to every cluster member."
    )
    with st.expander("What is in the default E3 focus authority?", expanded=False):
        st.markdown(
            """
            The bundled authority contains the versioned protein-level seed evidence inherited
            from the E3 project. Its accession, category, evidence type, organism and source
            provenance are retained. Users can provide a smaller or different TSV without code
            changes. Exact UniProt accessions and entry names are recognised only when they can
            be parsed unambiguously from a stored `sp|ACCESSION|ENTRY` or `tr|...` identifier.

            A matched seed means **this cluster contains a configured focus protein**. It is a
            prioritisation flag, not proof that all paralogues or orthologues share the same
            molecular function.
            """
        )
    template = parse_focus_proteins(data=focus_template(), source_name="custom_focus.tsv")
    render_table_downloads(
        records=tuple(_focus_record(row=row) for row in template.records),
        file_stem="custom_focus_proteins_template",
        key="focus_protein_template",
        tsv_label="Download custom focus template as TSV",
        excel_label="Download custom focus template as formatted Excel",
        workbook_title="Custom OrthoFinder focus proteins",
    )
    uploaded = st.file_uploader(
        "Optional custom focus-protein TSV",
        type=("tsv", "txt"),
        help=(
            "Upload takes precedence over --focus-proteins. Leave empty to use the "
            "packaged E3 seed authority. Gzip-compressed TSV is supported by the launcher path."
        ),
        key="focus_authority_upload",
    )
    try:
        authority = load_focus_authority(
            focus_path_text=focus_path_text,
            uploaded_data=uploaded.getvalue() if uploaded is not None else None,
            uploaded_name=uploaded.name if uploaded is not None else "",
        )
    except InputValidationError as error:
        _LOGGER.exception("Focus protein authority validation failed")
        st.error(f"Focus protein authority could not be used: {error}")
        return
    summary = authority.summary()
    metrics = st.columns(5)
    metrics[0].metric("Configured proteins", f"{summary['records']:,}")
    metrics[1].metric("Enabled proteins", f"{summary['enabled']:,}")
    metrics[2].metric("E3/categories", f"{summary['categories']:,}")
    metrics[3].metric("Evidence types", f"{summary['evidence_types']:,}")
    metrics[4].metric("Annotated organisms", f"{summary['organisms']:,}")
    st.caption(
        f"Authority: {authority.source_name}; SHA-256: {authority.sha256}. "
        "The checksum is retained in downstream selection manifests."
    )
    group_types = service.list_group_types()
    default_type = "HOG" if "HOG" in group_types else group_types[0]
    with st.form("focus_cluster_controls"):
        controls = st.columns(3)
        group_type = controls[0].selectbox(
            "Group system",
            group_types,
            index=group_types.index(default_type),
            help="Choose one exact group collection; HOG N0 is the default E3 view.",
        )
        nodes = service.list_hierarchy_nodes(group_type=group_type)
        default_node = "N0" if "N0" in nodes else nodes[0]
        node = controls[1].selectbox(
            "Species-tree level",
            nodes,
            index=nodes.index(default_node),
            format_func=lambda value: value or "ROOT",
        )
        maximum = int(
            controls[2].selectbox(
                "Maximum matching clusters",
                (500, 2_000, 5_000, 10_000),
                index=1,
                help="A declared browser and export bound; narrowing remains reproducible.",
            )
        )
        submitted = st.form_submit_button("Find focus clusters", type="primary")
    signature = (authority.sha256, group_type, node, maximum, service.resource.run_id)
    should_query = submitted or st.session_state.get(_FOCUS_SIGNATURE_STATE) != signature
    if should_query:
        try:
            with st.spinner("Matching the reviewed focus authority to exact protein aliases…"):
                result = service.search_focus_clusters(
                    filters=FocusClusterFilters(
                        protein_identifiers=authority.identifiers,
                        group_type=group_type,
                        hierarchy_node=node,
                        maximum_rows=maximum,
                    )
                )
        except (InputValidationError, OrthoFinderResultsError) as error:
            _LOGGER.exception("Focus cluster search failed")
            st.error(f"Focus clusters could not be queried: {error}")
            return
        st.session_state[_FOCUS_RESULT_STATE] = result
        st.session_state[_FOCUS_SIGNATURE_STATE] = signature
    result = st.session_state.get(_FOCUS_RESULT_STATE)
    if result is None:
        return
    result_metrics = st.columns(4)
    result_metrics[0].metric("Matching clusters", f"{result.total_rows:,}")
    result_metrics[1].metric(
        "Focus IDs found",
        f"{result.matched_focus_identifiers:,}",
        help="Distinct submitted focus identifiers matched anywhere in this group collection.",
    )
    result_metrics[2].metric(
        "Focus IDs not found",
        f"{result.submitted_focus_identifiers - result.matched_focus_identifiers:,}",
    )
    result_metrics[3].metric("Clusters displayed", f"{len(result.rows):,}")
    if not result.rows:
        st.warning(
            "No configured focus protein matched this exact group system and hierarchy. "
            "This reports identifier coverage only; it is not evidence of biological absence."
        )
        return
    if result.truncated:
        st.warning(
            f"{result.total_rows:,} clusters matched but {len(result.rows):,} are retained by "
            "the declared limit. Increase the bound or use a smaller focus authority before "
            "treating this table as complete."
        )
    metadata = {row.identifier: row for row in authority.enabled_records}
    displayed = tuple(_display_focus_result(row=row, metadata=metadata) for row in result.rows)
    st.subheader("Clusters containing configured focus proteins")
    st.dataframe(
        displayed,
        width="stretch",
        hide_index=True,
        column_config={
            label: st.column_config.Column(label=label, help=description)
            for label, description in _RESULT_HELP.items()
        },
    )
    render_table_downloads(
        records=displayed,
        file_stem=f"{service.resource.run_id}_{group_type}_{node or 'ROOT'}_focus_clusters",
        key="focus_cluster_results",
        tsv_label="Download focus clusters as TSV",
        excel_label="Download focus clusters as formatted Excel",
        column_definitions=_RESULT_HELP,
        workbook_title="OrthoFinder focus protein clusters",
    )
    labels = {
        _focus_option(row=row): index for index, row in enumerate(result.rows)
    }
    selected_label = st.selectbox("Focus cluster to inspect", tuple(labels))
    selected = result.rows[labels[selected_label]]
    selected_key = GroupKey(
        run_id=str(selected["run_id"]),
        group_type=str(selected["group_type"]),
        hierarchy_node=str(selected["hierarchy_node"]),
        group_id=str(selected["group_id"]),
    )
    if st.button("Show this focus cluster's complete visual suite", type="primary"):
        st.session_state[_FOCUS_OPEN_STATE] = selected_key.display_label()
    if st.session_state.get(_FOCUS_OPEN_STATE) == selected_key.display_label():
        members = str(selected["matched_member_ids"]).split(";")
        focus_member = members[0]
        st.info(
            f"Highlighting configured focus protein **{focus_member}** in "
            f"**{selected_key.display_label()}**."
        )
        render_selected_group_visualisations(
            resource=resource,
            service=service,
            cache_dir=cache_dir,
            key=selected_key,
            focus_member=focus_member,
        )


def load_focus_authority(
    *, focus_path_text: str, uploaded_data: bytes | None, uploaded_name: str = ""
) -> FocusProteinAuthority:
    """Load upload, launcher path or bundled E3 authority in clear precedence order."""

    if uploaded_data is not None:
        return parse_focus_proteins(
            data=uploaded_data,
            source_name=uploaded_name or "uploaded_focus_proteins.tsv",
            source_sha256=hashlib.sha256(uploaded_data).hexdigest(),
        )
    if focus_path_text.strip():
        path = Path(focus_path_text).expanduser().resolve()
    else:
        path = bundled_focus_path()
    try:
        stat = path.stat()
    except OSError as error:
        raise InputValidationError(f"Focus protein TSV is unavailable: {path}") from error
    return _read_focus_authority_cached(
        path_text=str(path),
        modified_ns=stat.st_mtime_ns,
        size_bytes=stat.st_size,
    )


@st.cache_data(show_spinner=False)
def _read_focus_authority_cached(
    *, path_text: str, modified_ns: int, size_bytes: int
) -> FocusProteinAuthority:
    """Cache a path authority by exact path, modification time and physical size."""

    del modified_ns, size_bytes
    return read_focus_proteins(path=Path(path_text))


def _focus_record(*, row: Any) -> dict[str, Any]:
    """Return one focus authority row with editable canonical headings."""

    return {
        "protein_identifier": row.identifier,
        "protein_name": row.protein_name,
        "category": row.category,
        "evidence_type": row.evidence_type,
        "organism": row.organism,
        "source": row.source,
        "enabled": row.enabled,
        "note": row.note,
    }


def _display_focus_result(
    *, row: Mapping[str, Any], metadata: Mapping[str, Any]
) -> dict[str, Any]:
    """Return one focus-cluster result with provenance-rich readable headings."""

    identifiers = tuple(str(row["matched_focus_identifiers"]).split(";"))
    categories = sorted(
        {
            str(metadata[identifier].category)
            for identifier in identifiers
            if identifier in metadata and metadata[identifier].category
        }
    )
    evidence = sorted(
        {
            str(metadata[identifier].evidence_type)
            for identifier in identifiers
            if identifier in metadata and metadata[identifier].evidence_type
        }
    )
    return {
        "Group system": row["group_type"],
        "Species-tree level": row["hierarchy_node"] or "ROOT",
        "Group ID": row["group_id"],
        "Matched focus IDs": row["matched_focus_count"],
        "Matched focus proteins": row["matched_protein_count"],
        "Matching identifiers": row["matched_focus_identifiers"],
        "Matched protein IDs": row["matched_member_ids"],
        "Matched species": row["matched_species_labels"],
        "Matched through": row["match_authorities"],
        "Focus categories": ";".join(categories),
        "Evidence types": ";".join(evidence),
        "Proteins in group": row["member_count"],
        "Species represented": row["species_count"],
        "Stored distance status": row.get("computation_status") or "Not calculated",
        "Stored mean distance": row.get("mean_distance"),
        "Stored distance spread (SD)": row.get("population_stddev_distance"),
    }


def _focus_option(*, row: Mapping[str, Any]) -> str:
    """Return one unambiguous focus-cluster selector label."""

    return (
        f"{row['group_type']} | {row['hierarchy_node'] or 'ROOT'} | {row['group_id']} "
        f"({row['matched_focus_count']} focus IDs)"
    )

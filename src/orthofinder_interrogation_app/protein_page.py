"""Protein-centred discovery, cluster inspection and distance reporting."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import streamlit as st

from orthofinder_results.errors import InputValidationError, OrthoFinderResultsError

from .evolutionary_page import render_selected_group_visualisations
from .exports import render_table_downloads
from .models import GroupKey, ProteinSearchFilters
from .queries import OrthoFinderQueryService

_LOGGER = logging.getLogger("orthofinder_interrogation_app.protein_page")
_MATCH_MODE_LABELS = {
    "Exact identifier": "EXACT",
    "Identifier contains text": "CONTAINS",
}
_MATCH_SOURCE_LABELS = {
    "PROTEIN_ID": "Protein ID",
    "ORTHOFINDER_INTERNAL_ID": "OrthoFinder internal ID",
    "PROTEIN_AND_INTERNAL_ID": "Protein and internal ID",
}
_RESULT_COLUMN_HELP = {
    "Protein ID": "Canonical protein identifier stored in the group membership.",
    "OrthoFinder internal ID": (
        "Internal sequence identifier from SequenceIDs.txt, when available."
    ),
    "Species": "Exact species label assigned to this protein in the resource.",
    "Matched through": "Identifier authority that matched the submitted search text.",
    "Matched identifier": "Exact stored identifier that satisfied the search.",
    "Group system": (
        "HOG is hierarchical; LEGACY_ORTHOGROUP is the flat Orthogroups.tsv set."
    ),
    "Species-tree level": "Species-tree node defining this HOG; ROOT denotes the flat set.",
    "Group ID": "Exact OrthoFinder group identifier at the stated level.",
    "Parent legacy orthogroup": "Linked flat orthogroup identifier when supplied.",
    "Gene-tree parent clade": "Gene-tree clade used by OrthoFinder to define the HOG.",
    "Proteins in group": "Complete protein membership of this group record.",
    "Species represented": "Number of species contributing at least one protein.",
    "Highest copies in one species": "Largest protein copy count for one species.",
    "Average copies per represented species": (
        "Protein count divided by the represented-species count."
    ),
    "Stored distance coverage": (
        "Whether a complete or bounded pair-distance matrix is already stored."
    ),
    "Proteins in stored matrix": "Proteins represented in the preferred stored matrix.",
    "Stored protein pairs": "Pair distances in the preferred stored matrix.",
    "Stored average distance": "Mean stored pair distance when already calculated.",
    "Stored distance spread (SD)": "Population SD across stored pair distances.",
}
_DISTANCE_STATUS_LABELS = {
    "COMPLETE": "Complete distances stored",
    "EXACT": "Complete distances stored",
    "DETERMINISTIC_MEMBER_SAMPLE": "Sampled distances stored",
    "PORTABLE_TREE_LAZY": "Calculated on demand from portable tree",
}


def render_protein_search(
    *,
    resource: Any,
    service: OrthoFinderQueryService,
    cache_dir: Path,
) -> None:
    """Render protein lookup, matching clusters and focused visual analysis.

    Args:
        resource: Validated immutable resource identity.
        service: Read-only query service for the resource.
        cache_dir: Persistent on-demand analysis cache outside the resource.
    """

    st.header("Find a gene or protein")
    st.write(
        "Enter a protein identifier to find every HOG and legacy orthogroup containing it. "
        "Choose one result to highlight that protein across the full cluster visualisation "
        "suite and obtain its direct distances to the other analysed proteins."
    )
    with st.expander("Which identifiers can be searched?", expanded=False):
        st.markdown(
            """
            - **Protein ID** searches the exact identifier stored in OrthoFinder membership
              tables. This is normally the first space-delimited identifier from the FASTA
              header.
            - **OrthoFinder internal ID** searches values such as `0_123` from
              `SequenceIDs.txt`, where that authority is available.
            - Descriptive gene names, symbols and functional annotations can only be found if
              they are part of the stored protein identifier; this resource does not guess
              aliases from external databases.
            - One protein can legitimately occur in a flat orthogroup and in several HOGs at
              different species-tree levels. These are returned as separate cluster records.
            """
        )
    with st.form("protein_identifier_search"):
        controls = st.columns((4, 2, 2, 1))
        query = controls[0].text_input(
            "Protein or OrthoFinder internal ID",
            placeholder="For example Q9SA03 or 0_123",
            help=(
                "Exact matching is case-sensitive. Contains matching is literal "
                "and case-insensitive."
            ),
        )
        match_label = controls[1].selectbox(
            "How to match",
            tuple(_MATCH_MODE_LABELS),
            help="Use Exact identifier whenever the complete ID is known.",
        )
        group_type_label = controls[2].selectbox(
            "Group system",
            ("All", *service.list_group_types()),
            help="Search every authority or restrict the result to HOGs or flat orthogroups.",
        )
        maximum_rows = int(
            controls[3].selectbox(
                "Maximum matches",
                (100, 250, 500, 1_000),
                index=2,
                help="A defensive browser limit; exact searches normally return far fewer rows.",
            )
        )
        st.form_submit_button("Find this protein", type="primary")
    if not query.strip():
        st.info("Enter an identifier above to begin the protein-centred search.")
        return
    try:
        result = service.search_proteins(
            filters=ProteinSearchFilters(
                query=query,
                match_mode=_MATCH_MODE_LABELS[match_label],
                group_type="" if group_type_label == "All" else group_type_label,
                maximum_rows=maximum_rows,
            )
        )
    except (InputValidationError, OrthoFinderResultsError) as error:
        _LOGGER.exception("Protein search failed")
        st.error(str(error))
        return
    if not result.rows:
        st.warning(
            "No stored protein or OrthoFinder internal identifier matches this search. "
            "Check the exact identifier and selected group system."
        )
        return
    st.subheader("Matching cluster memberships")
    summary_columns = st.columns(3)
    summary_columns[0].metric(
        "Matching cluster records",
        f"{result.total_rows:,}",
        help="Separate HOG hierarchy levels and the legacy orthogroup are counted separately.",
    )
    summary_columns[1].metric(
        "Matching proteins",
        f"{len({str(row['member_id']) for row in result.rows}):,}",
        help="Distinct canonical protein identifiers in the returned bounded rows.",
    )
    summary_columns[2].metric(
        "Species represented",
        f"{len({str(row['species_label']) for row in result.rows}):,}",
        help="Distinct exact species labels among the returned matches.",
    )
    if result.truncated:
        st.warning(
            f"The search matched {result.total_rows:,} cluster records; only the first "
            f"{len(result.rows):,} are shown. Use an exact identifier or a narrower group "
            "system before interpreting the result as complete."
        )
    displayed = tuple(_display_protein_result(row=row) for row in result.rows)
    st.dataframe(
        displayed,
        width="stretch",
        hide_index=True,
        column_config={
            label: st.column_config.Column(label=label, help=description)
            for label, description in _RESULT_COLUMN_HELP.items()
        },
    )
    render_table_downloads(
        records=displayed,
        file_stem=f"{query}_protein_cluster_matches",
        key="protein_cluster_matches_download",
        tsv_label="Download matching clusters as TSV",
        excel_label="Download matching clusters as formatted Excel",
        column_definitions=_RESULT_COLUMN_HELP,
        workbook_title=f"OrthoFinder protein search: {query}",
    )
    options = _protein_result_options(rows=result.rows)
    selected_label = st.selectbox(
        "Cluster to inspect",
        tuple(options),
        help=(
            "Each HOG level is a distinct record. Select the biologically relevant level "
            "before interpreting distances and topology."
        ),
    )
    selected = result.rows[options[selected_label]]
    key = _result_group_key(row=selected)
    st.info(
        f"Focusing on **{selected['member_id']}** from **{selected['species_label']}** "
        f"in **{key.display_label()}**."
    )
    render_selected_group_visualisations(
        resource=resource,
        service=service,
        cache_dir=cache_dir,
        key=key,
        focus_member=str(selected["member_id"]),
    )


def _display_protein_result(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one protein search result with readable headings.

    Args:
        row: Raw query result.

    Returns:
        Ordered display and export record.
    """

    return {
        "Protein ID": row["member_id"],
        "OrthoFinder internal ID": row.get("internal_id") or "",
        "Species": row["species_label"],
        "Matched through": _MATCH_SOURCE_LABELS.get(
            str(row["match_source"]), str(row["match_source"])
        ),
        "Matched identifier": row["matched_identifier"],
        "Group system": row["group_type"],
        "Species-tree level": row["hierarchy_node"] or "ROOT",
        "Group ID": row["group_id"],
        "Parent legacy orthogroup": row["legacy_orthogroup_id"],
        "Gene-tree parent clade": row["gene_tree_parent_clade"],
        "Proteins in group": row["member_count"],
        "Species represented": row["species_count"],
        "Highest copies in one species": row["max_copies_per_species"],
        "Average copies per represented species": row["mean_copies_per_species"],
        "Stored distance coverage": _distance_status_label(
            value=row.get("computation_status")
        ),
        "Proteins in stored matrix": row.get("sampled_member_count"),
        "Stored protein pairs": row.get("distance_pair_count"),
        "Stored average distance": row.get("mean_distance"),
        "Stored distance spread (SD)": row.get("population_stddev_distance"),
    }


def _protein_result_options(
    *, rows: tuple[dict[str, Any], ...]
) -> dict[str, int]:
    """Return unique readable selector labels mapped to result positions.

    Args:
        rows: Bounded raw protein search results.

    Returns:
        Insertion-ordered label-to-row-index mapping.

    Raises:
        InputValidationError: If no result rows are supplied.
    """

    if not rows:
        raise InputValidationError("Protein result selection requires at least one row.")
    options = {}
    for index, row in enumerate(rows):
        node = str(row["hierarchy_node"]) or "ROOT"
        label = (
            f"{index + 1}. {row['member_id']} [{row['species_label']}] → "
            f"{row['group_type']} | {node} | {row['group_id']}"
        )
        options[label] = index
    return options


def _result_group_key(*, row: Mapping[str, Any]) -> GroupKey:
    """Return the exact group key for one protein search result.

    Args:
        row: Raw result containing the four composite identity fields.

    Returns:
        Valid group key.
    """

    return GroupKey(
        run_id=str(row["run_id"]),
        group_type=str(row["group_type"]),
        hierarchy_node=str(row["hierarchy_node"]),
        group_id=str(row["group_id"]),
    )


def _distance_status_label(*, value: object) -> str:
    """Translate one stored-distance state without hiding unfamiliar codes."""

    if value is None or str(value).strip() == "":
        return "Not calculated"
    code = str(value)
    return _DISTANCE_STATUS_LABELS.get(code, code.replace("_", " ").title())

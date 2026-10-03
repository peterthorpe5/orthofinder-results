"""Streamlit page for conserved C-terminal protein motifs across HOGs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import plotly.graph_objects as go
import streamlit as st

from orthofinder_results.errors import InputValidationError

from .documentation_page import render_page_guidance
from .exports import render_plotly_figure, render_table_downloads
from .models import ResourceIdentity
from .terminal_motif import (
    DEFAULT_MOTIF,
    DEFAULT_THRESHOLD,
    motif_group_members,
    motif_group_summary,
    motif_species,
    validate_motif,
    validate_sequence_sidecar,
)

_SUMMARY_HELP = {
    "HOG": "Root or selected hierarchy-level HOG identifier.",
    "Proteins with sequences": "Distinct sequence-bearing proteins used as the denominator.",
    "Species represented": "Species contributing at least one sequence-bearing protein.",
    "Matching proteins": "Proteins ending with the exact selected motif.",
    "Matching fraction": "Matching proteins divided by all sequence-bearing proteins in the HOG.",
    "Matching species": "Species with at least one protein ending with the motif.",
    "Matching species list": "Exact OrthoFinder labels for matching species.",
}
_MEMBER_HELP = {
    "Species": "Exact species label from the OrthoFinder run.",
    "Protein ID": "Original member identifier mapped by SequenceIDs.txt.",
    "Internal ID": "OrthoFinder internal sequence identifier.",
    "Sequence length": "Protein length after removal of a terminal FASTA stop marker.",
    "Observed terminus": "Last residues, using the same length as the selected motif.",
    "Matches motif": "True only when the complete exact C-terminal motif matches.",
    "Sequence": "Complete protein sequence used for the call.",
}


def render_terminal_motif_page(
    *, resource: ResourceIdentity, sidecar_path_text: str
) -> None:
    """Render flexible terminal-motif conservation discovery and exports."""

    st.header("C-terminal motif conservation")
    render_page_guidance(key="terminal_motif")
    st.write(
        "Find HOGs in which an exact protein C-terminal motif is conserved across "
        "proteins and species. The grant-facing default is terminal asparagine (N), "
        "but any canonical amino-acid suffix up to 100 residues can be tested."
    )
    if not sidecar_path_text.strip():
        st.info(
            "This analysis needs the complete-proteome sequence sidecar. Build it once "
            "from the same OrthoFinder Results_* directory, then launch the app with "
            "--terminal-motif-parquet. The completed resource itself is not modified."
        )
        st.code(
            "orthofinder-terminal-motif-build \\\n"
            "  --orthofinder-results-dir /path/to/Results_Feb26 \\\n"
            "  --output-parquet /path/to/terminal_motif_sequences.parquet",
            language="bash",
        )
        return
    try:
        sidecar = validate_sequence_sidecar(path=Path(sidecar_path_text))
        species = motif_species(sidecar_path=sidecar)
    except InputValidationError as error:
        st.error(str(error))
        return

    controls = st.columns((1.0, 1.4, 1.0, 1.2))
    motif_text = controls[0].text_input(
        "Exact C-terminal motif", value=DEFAULT_MOTIF, help="Canonical one-letter codes only."
    )
    threshold_percent = controls[1].slider(
        "Minimum matching proteins", min_value=0, max_value=100,
        value=int(DEFAULT_THRESHOLD * 100), step=1,
        help="Percentage of sequence-bearing proteins in a HOG ending with the exact motif.",
    )
    minimum_species = controls[2].number_input(
        "Minimum species", min_value=1, max_value=max(1, len(species)),
        value=min(3, max(1, len(species))), step=1,
    )
    maximum_rows = controls[3].number_input(
        "Maximum HOGs", min_value=10, max_value=20_000, value=2_000, step=10
    )
    required_species = st.multiselect(
        "Species that must have at least one matching protein",
        options=species,
        help=(
            "Use this to require Arabidopsis, human or other focal species. Leaving it "
            "empty does not privilege any species."
        ),
    )
    hierarchy_node = st.text_input(
        "HOG hierarchy node", value="N0",
        help="N0 is the root HOG level in the current analysis; exact labels are retained.",
    )
    try:
        motif = validate_motif(motif=motif_text)
        rows = motif_group_summary(
            resource=resource,
            sidecar_path=sidecar,
            motif=motif,
            threshold=float(threshold_percent) / 100.0,
            minimum_species=int(minimum_species),
            required_species=required_species,
            hierarchy_node=hierarchy_node,
            maximum_rows=int(maximum_rows),
        )
    except InputValidationError as error:
        st.error(str(error))
        return
    st.caption(
        f"{len(rows):,} HOGs pass the current filters. Denominators include only proteins "
        "with successfully reconciled sequences."
    )
    if not rows:
        st.warning("No HOG passes the current motif, conservation and species filters.")
        return
    display_rows = tuple(_display_summary(row=row) for row in rows)
    metrics = st.columns(4)
    metrics[0].metric("Passing HOGs", f"{len(rows):,}")
    metrics[1].metric("Motif", motif)
    metrics[2].metric("Threshold", f"{threshold_percent}%")
    metrics[3].metric("Sequence species", f"{len(species):,}")

    st.subheader("Conservation landscape")
    with st.expander("What this graph shows and how to interpret it"):
        st.write(
            "Each point is one HOG. The horizontal axis is the fraction of its "
            "sequence-bearing proteins ending with the exact motif; the vertical axis "
            "is the number of species with at least one match. Larger points contain "
            "more proteins. Strong pilot candidates lie towards the upper right, but "
            "large paralogue expansions and uneven species sampling must be checked in "
            "the protein table. This is not evidence of Cereblon binding or degradation."
        )
    figure = _motif_figure(rows=rows, motif=motif, threshold=threshold_percent / 100.0)
    render_plotly_figure(
        figure=figure,
        file_stem=f"terminal_motif_{motif}_conservation",
        key="terminal_motif_conservation_figure",
        pdf_width=1600,
        pdf_height=1000,
    )
    st.dataframe(display_rows, width="stretch", hide_index=True)
    render_table_downloads(
        records=display_rows,
        file_stem=f"terminal_motif_{motif}_hog_summary",
        key="terminal_motif_summary_tsv",
        column_definitions=_SUMMARY_HELP,
        workbook_title=f"C-terminal motif {motif} HOG summary",
    )

    st.subheader("Inspect one HOG")
    selected = st.selectbox(
        "HOG", options=[str(row["group_id"]) for row in rows], key="terminal_motif_hog"
    )
    members = motif_group_members(
        resource=resource,
        sidecar_path=sidecar,
        motif=motif,
        group_id=selected,
        hierarchy_node=hierarchy_node,
    )
    member_rows = tuple(_display_member(row=row) for row in members)
    st.dataframe(member_rows, width="stretch", hide_index=True)
    render_table_downloads(
        records=member_rows,
        file_stem=f"{selected}_{motif}_terminal_calls",
        key="terminal_motif_member_tsv",
        column_definitions=_MEMBER_HELP,
        workbook_title=f"{selected} C-terminal motif calls",
    )
    st.download_button(
        "Download selected HOG as FASTA",
        data=_members_to_fasta(rows=members),
        file_name=f"{selected}_{motif}_terminal_candidates.fasta",
        mime="text/x-fasta",
        key="terminal_motif_fasta",
    )


def _display_summary(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one readable HOG summary row."""

    return {
        "HOG": row["group_id"],
        "Proteins with sequences": row["sequence_count"],
        "Species represented": row["species_count"],
        "Matching proteins": row["matching_sequence_count"],
        "Matching fraction": row["matching_fraction"],
        "Matching species": row["matching_species_count"],
        "Matching species list": row["matching_species"] or "",
    }


def _display_member(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one readable protein-level motif call."""

    return {
        "Species": row["species_label"],
        "Protein ID": row["member_id"],
        "Internal ID": row["internal_id"],
        "Sequence length": row["sequence_length"],
        "Observed terminus": row["observed_terminus"],
        "Matches motif": row["motif_match"],
        "Sequence": row["sequence"],
    }


def _motif_figure(
    *, rows: Sequence[Mapping[str, Any]], motif: str, threshold: float
) -> go.Figure:
    """Build the interactive protein-fraction versus species-breadth plot."""

    figure = go.Figure(
        go.Scatter(
            x=[float(row["matching_fraction"]) for row in rows],
            y=[int(row["matching_species_count"]) for row in rows],
            mode="markers",
            customdata=[
                [row["group_id"], row["sequence_count"], row["species_count"]]
                for row in rows
            ],
            marker={
                "size": [max(7, min(30, 5 + int(row["sequence_count"]) ** 0.5)) for row in rows],
                "color": [int(row["matching_species_count"]) for row in rows],
                "colorscale": "Viridis",
                "showscale": True,
                "colorbar": {"title": "Matching species"},
                "opacity": 0.75,
            },
            hovertemplate=(
                "HOG=%{customdata[0]}<br>Matching fraction=%{x:.1%}<br>"
                "Matching species=%{y}<br>Sequences=%{customdata[1]}<br>"
                "Represented species=%{customdata[2]}<extra></extra>"
            ),
        )
    )
    figure.add_vline(x=threshold, line_dash="dash", line_color="#d95f02")
    figure.update_layout(
        title=f"HOG conservation of exact C-terminal motif {motif}",
        xaxis_title="Fraction of sequence-bearing proteins matching motif",
        yaxis_title="Species with at least one matching protein",
        xaxis={"tickformat": ".0%", "range": [0, 1.01]},
        template="plotly_white",
        height=700,
    )
    return figure


def _members_to_fasta(*, rows: Sequence[Mapping[str, Any]]) -> str:
    """Serialise selected HOG members as wrapped FASTA."""

    chunks: list[str] = []
    for row in rows:
        header = (
            f">{row['member_id']} species={row['species_label']} "
            f"motif_match={row['motif_match']}"
        )
        sequence = str(row["sequence"])
        wrapped = (sequence[index:index + 80] for index in range(0, len(sequence), 80))
        chunks.extend((header, *wrapped))
    return "\n".join(chunks) + "\n"

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
from .taxonomy import (
    TaxonomyAuthority,
    read_matching_bundled_taxonomy,
    read_taxonomy_mapping,
)
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
    "Focal lineage represented": "Reviewed descendant species present in this HOG.",
    "Focal lineage matching": "Reviewed descendant species with at least one motif match.",
    "Focal lineage coverage": "Matching descendants divided by represented descendants.",
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
    *,
    resource: ResourceIdentity,
    sidecar_path_text: str,
    taxonomy_path_text: str = "",
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
        st.warning(
            "No HOG motif analysis has run in this view. The open DuckDB contains HOG "
            "memberships, but not the complete amino-acid sequences needed to inspect "
            "their terminal residues."
        )
        st.subheader("One-time setup")
        st.write(
            "Build the sequence sidecar from the same OrthoFinder Results_* directory, "
            "copy it beside the completed resource if necessary, and restart the viewer "
            "with `--terminal-motif-parquet`. This does not rerun OrthoFinder or alter "
            "the completed resource. HOG and taxonomic controls appear after it loads."
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

    authority = _load_taxonomy(species=species, taxonomy_path_text=taxonomy_path_text)

    controls = st.columns((1.0, 1.4, 1.0, 1.2))
    motif_text = controls[0].text_input(
        "Exact C-terminal motif", value=DEFAULT_MOTIF, help="Canonical one-letter codes only."
    )
    threshold_percent = controls[1].slider(
        "Minimum matching proteins",
        min_value=0,
        max_value=100,
        value=int(DEFAULT_THRESHOLD * 100),
        step=1,
        help="Percentage of sequence-bearing proteins in a HOG ending with the exact motif.",
    )
    minimum_species = controls[2].number_input(
        "Minimum species",
        min_value=1,
        max_value=max(1, len(species)),
        value=min(3, max(1, len(species))),
        step=1,
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
        "HOG hierarchy node",
        value="N0",
        help="N0 is the root HOG level in the current analysis; exact labels are retained.",
    )
    focal_taxon_id: int | None = None
    required_taxon_ids: tuple[int, ...] = ()
    excluded_taxon_ids: tuple[int, ...] = ()
    lineage_threshold = 0.0
    if authority is None:
        st.warning(
            "No reviewed taxonomy authority exactly matches these sequence labels. "
            "Exact-species filtering remains available, but lineage filtering is "
            "disabled. Supply a reviewed taxonomy TSV under Advanced settings."
        )
    else:
        st.subheader("Taxonomic conservation filters")
        options = authority.taxon_options()
        option_by_label = {option.display_label(): option for option in options}
        labels = tuple(option_by_label)
        taxonomy_controls = st.columns((1.5, 1.0))
        focal_label = taxonomy_controls[0].selectbox(
            "Focal lineage",
            options=("All reviewed sampled species", *labels),
            help=(
                "Restrict conservation calculations to reviewed sampled descendants "
                "of this taxon, for example flowering plants, eudicots or humans."
            ),
        )
        lineage_threshold = (
            taxonomy_controls[1].slider(
                "Minimum matching descendants",
                min_value=0,
                max_value=100,
                value=80,
                step=1,
                help=(
                    "Percentage of represented descendant species that must contain at "
                    "least one protein ending with the motif."
                ),
            )
            / 100.0
        )
        if focal_label != "All reviewed sampled species":
            focal_taxon_id = option_by_label[focal_label].taxon_id
        required_labels = st.multiselect(
            "Lineages that must contain a matching descendant",
            options=labels,
            help=(
                "Optional cross-clade requirement. For example, require both a plant "
                "lineage and Homo to contain at least one matching protein."
            ),
        )
        excluded_labels = st.multiselect(
            "Lineages that must not contain a matching descendant",
            options=labels,
            help="Optional negative control; leave empty for the usual discovery analysis.",
        )
        required_taxon_ids = tuple(option_by_label[label].taxon_id for label in required_labels)
        excluded_taxon_ids = tuple(option_by_label[label].taxon_id for label in excluded_labels)
    try:
        motif = validate_motif(motif=motif_text)
        query_rows = motif_group_summary(
            resource=resource,
            sidecar_path=sidecar,
            motif=motif,
            threshold=float(threshold_percent) / 100.0,
            minimum_species=int(minimum_species),
            required_species=required_species,
            hierarchy_node=hierarchy_node,
            maximum_rows=20_000 if authority is not None else int(maximum_rows),
        )
        rows = _apply_taxonomy_filters(
            rows=query_rows,
            authority=authority,
            focal_taxon_id=focal_taxon_id,
            minimum_lineage_fraction=lineage_threshold,
            required_taxon_ids=required_taxon_ids,
            excluded_taxon_ids=excluded_taxon_ids,
        )[: int(maximum_rows)]
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
    if authority is not None:
        st.subheader("Taxonomic distribution of the selected HOG")
        taxon_rows = _species_distribution(rows=members, authority=authority)
        st.dataframe(taxon_rows, width="stretch", hide_index=True)
        render_table_downloads(
            records=taxon_rows,
            file_stem=f"{selected}_{motif}_taxonomic_distribution",
            key="terminal_motif_taxonomy_tsv",
            column_definitions={},
            workbook_title=f"{selected} taxonomic motif distribution",
        )
    st.subheader("Protein-level calls")
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
        "Focal lineage represented": row.get("focal_represented_species", ""),
        "Focal lineage matching": row.get("focal_matching_species", ""),
        "Focal lineage coverage": row.get("focal_matching_fraction", ""),
    }


def _load_taxonomy(
    *, species: tuple[str, ...], taxonomy_path_text: str
) -> TaxonomyAuthority | None:
    """Load an exact reviewed taxonomy authority without guessing labels."""

    if taxonomy_path_text.strip():
        return read_taxonomy_mapping(path=Path(taxonomy_path_text), expected_species=species)
    return read_matching_bundled_taxonomy(expected_species=species)


def _species_set(value: Any) -> set[str]:
    """Parse a semicolon-delimited species aggregation."""

    return {item for item in str(value or "").split("; ") if item}


def _apply_taxonomy_filters(
    *,
    rows: Sequence[Mapping[str, Any]],
    authority: TaxonomyAuthority | None,
    focal_taxon_id: int | None,
    minimum_lineage_fraction: float,
    required_taxon_ids: Sequence[int],
    excluded_taxon_ids: Sequence[int],
) -> tuple[dict[str, Any], ...]:
    """Annotate and filter HOG summaries with reviewed descendant sets."""

    if authority is None:
        return tuple(dict(row) for row in rows)
    focal = set(
        authority.reviewed_species
        if focal_taxon_id is None
        else authority.target_species(taxon_id=focal_taxon_id)
    )
    required = [set(authority.target_species(taxon_id=value)) for value in required_taxon_ids]
    excluded = [set(authority.target_species(taxon_id=value)) for value in excluded_taxon_ids]
    result: list[dict[str, Any]] = []
    for source in rows:
        represented = _species_set(source.get("represented_species"))
        matching = _species_set(source.get("matching_species"))
        focal_represented = represented & focal
        focal_matching = matching & focal
        fraction = len(focal_matching) / len(focal_represented) if focal_represented else 0.0
        if not focal_represented or fraction < minimum_lineage_fraction:
            continue
        if any(not matching.intersection(target) for target in required):
            continue
        if any(matching.intersection(target) for target in excluded):
            continue
        row = dict(source)
        row["focal_represented_species"] = len(focal_represented)
        row["focal_matching_species"] = len(focal_matching)
        row["focal_matching_fraction"] = fraction
        result.append(row)
    return tuple(result)


def _species_distribution(
    *, rows: Sequence[Mapping[str, Any]], authority: TaxonomyAuthority
) -> tuple[dict[str, Any], ...]:
    """Summarise motif calls by reviewed species for one HOG."""

    records = {row.workflow_species_label: row for row in authority.reviewed_records}
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row["species_label"]), []).append(row)
    output = []
    for species, members in sorted(grouped.items()):
        record = records.get(species)
        matches = sum(bool(row["motif_match"]) for row in members)
        output.append(
            {
                "OrthoFinder species": species,
                "Accepted species": record.accepted_species_name if record else "Unresolved",
                "NCBI taxon ID": record.ncbi_taxon_id if record else "",
                "Proteins": len(members),
                "Matching proteins": matches,
                "Species has motif": matches > 0,
            }
        )
    return tuple(output)


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


def _motif_figure(*, rows: Sequence[Mapping[str, Any]], motif: str, threshold: float) -> go.Figure:
    """Build the interactive protein-fraction versus species-breadth plot."""

    figure = go.Figure(
        go.Scatter(
            x=[float(row["matching_fraction"]) for row in rows],
            y=[int(row["matching_species_count"]) for row in rows],
            mode="markers",
            customdata=[
                [row["group_id"], row["sequence_count"], row["species_count"]] for row in rows
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
            f">{row['member_id']} species={row['species_label']} motif_match={row['motif_match']}"
        )
        sequence = str(row["sequence"])
        wrapped = (sequence[index : index + 80] for index in range(0, len(sequence), 80))
        chunks.extend((header, *wrapped))
    return "\n".join(chunks) + "\n"

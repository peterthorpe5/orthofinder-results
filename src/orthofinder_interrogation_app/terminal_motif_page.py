"""Streamlit page for taxonomic protein-sequence motif conservation."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from orthofinder_results.errors import InputValidationError

from .documentation_page import render_page_guidance
from .exports import render_plotly_figure, render_table_downloads
from .expression import (
    CONTEXT_COLUMNS,
    expression_available,
    expression_context_records,
    expression_dimensions,
    expression_group_species,
    expression_group_summaries,
    expression_heatmap_cells,
    expression_member_evidence,
    expression_selected_member_evidence,
)
from .models import ResourceIdentity
from .taxonomy import (
    TaxonomyAuthority,
    read_matching_bundled_taxonomy,
    read_taxonomy_mapping,
)
from .terminal_motif import (
    DEFAULT_MOTIF,
    DEFAULT_THRESHOLD,
    EXACT_C_TERMINAL,
    REGEX_ANYWHERE,
    REGEX_C_TERMINAL,
    SequenceSearch,
    motif_group_members,
    motif_group_summary,
    motif_hierarchy_nodes,
    motif_species,
    sequence_authority_available,
    validate_sequence_search,
    validate_sequence_sidecar,
)

_SUMMARY_HELP = {
    "Group": "Exact HOG or original OrthoFinder orthogroup identifier.",
    "Members": "All distinct proteins assigned to the group.",
    "Sequences assessed": "Members with a reconciled non-empty protein sequence.",
    "Sequence coverage": "Sequences assessed divided by all group members.",
    "Analysis proteins": "Assessed proteins in the selected primary species set.",
    "Matching analysis proteins": "Primary-set proteins matching the search.",
    "Analysis protein match": "Matching divided by assessed primary-set proteins.",
    "Analysis species assessed": "Primary species represented by an assessed sequence.",
    "Analysis species matching": "Primary species with at least one matching protein.",
    "Analysis species match": "Matching divided by assessed primary species.",
    "All matching species": "All sampled species with at least one matching protein.",
}
_SPECIES_HELP = {
    "OrthoFinder species": "Exact species label in the completed OrthoFinder run.",
    "Accepted species": "Reviewed scientific name from the taxonomy authority.",
    "NCBI taxon ID": "Reviewed NCBI Taxonomy identifier.",
    "Analysis role": "Primary denominator, required comparison, excluded or other sampled.",
    "Published proteins": "All group members for this species.",
    "Assessed": "Members with a reconciled sequence.",
    "Matching": "Assessed proteins matching the active search.",
    "Unavailable": "Members without an assessable sequence.",
    "Match fraction": "Matching divided by assessed proteins for this species.",
    "Species result": "Human-readable evidence state; unavailable is not measured zero.",
}
_MEMBER_HELP = {
    "Accepted species": "Reviewed species name when available.",
    "OrthoFinder species": "Exact label used by the source OrthoFinder run.",
    "NCBI taxon ID": "Reviewed NCBI Taxonomy identifier.",
    "Analysis role": "How this species participates in the current comparison.",
    "Protein ID": "Original member identifier mapped by SequenceIDs.txt.",
    "Internal ID": "OrthoFinder internal sequence identifier.",
    "Source FASTA": "Original proteome FASTA recorded by OrthoFinder.",
    "Protein description": "Description retained from the original FASTA header when available.",
    "Raw FASTA header": "Complete source header retained for identifier and provenance audit.",
    "Sequence available": "Whether a non-empty reconciled sequence was assessable.",
    "Sequence length": "Protein length after removal of a terminal FASTA stop marker.",
    "Matches search": "Whether this protein matches the active exact or regex search.",
    "Matched sequence": "First matching sequence fragment, or the exact suffix.",
    "Sequence": "Complete protein sequence; hidden in the on-screen table by default.",
}
_CANDIDATE_VIEW = "Motif candidates"
_TAXONOMY_VIEW = "Taxonomic conservation"
_GROUP_VIEW = "Protein & group explorer"
_EXPRESSION_VIEW = "RNA-seq explorer"


def render_terminal_motif_page(
    *,
    resource: ResourceIdentity,
    sidecar_path_text: str,
    taxonomy_path_text: str = "",
    expression_only: bool = False,
) -> None:
    """Render flexible sequence conservation or its RNA-seq explorer.

    Args:
        resource: Validated immutable OrthoFinder resource.
        sidecar_path_text: Optional external sequence-sidecar path for older resources.
        taxonomy_path_text: Optional reviewed taxonomy authority path.
        expression_only: Open directly into the RNA-seq results workspace.
    """

    st.header("RNA-seq explorer" if expression_only else "Protein motif conservation")
    render_page_guidance(key="terminal_motif")
    if expression_only:
        st.write(
            "Explore transcript evidence for motif-qualified OrthoFinder groups across the "
            "reviewed RNA-seq species panel. First define the sequence, group and taxonomic "
            "filters; the heatmap, UpSet intersections and downloadable evidence tables then "
            "update for that exact candidate set."
        )
        if not expression_available(resource=resource):
            _render_expression_setup_state()
            return
        st.success(
            "RNA-seq evidence is available for the reviewed 12-plant panel plus Homo sapiens. "
            "All other OrthoFinder species remain available to the motif and taxonomic filters "
            "but are not assessed for expression."
        )
    else:
        st.write(
            "Find OrthoFinder groups whose proteins share an exact C-terminal motif or, "
            "when explicitly enabled, a regular-expression pattern. Qualification is "
            "calculated over a transparent primary species set; other lineages remain "
            "available as independent comparison evidence."
        )
    sidecar = Path(sidecar_path_text) if sidecar_path_text.strip() else None
    if not sequence_authority_available(resource=resource, sidecar_path=sidecar):
        _render_setup_state()
        return
    try:
        if sidecar is not None:
            sidecar = validate_sequence_sidecar(path=sidecar)
        species = motif_species(resource=resource, sidecar_path=sidecar)
        authority = _load_taxonomy(
            species=species,
            taxonomy_path_text=taxonomy_path_text,
        )
    except InputValidationError as error:
        st.error(str(error))
        return

    selected_view = (
        _EXPRESSION_VIEW
        if expression_only
        else _render_result_workspace_selector(resource=resource)
    )
    st.subheader("Define candidate groups" if expression_only else "Search and filters")
    try:
        search = _render_search_controls()
        group_type, hierarchy_node = _render_group_controls(resource=resource)
        filters = _render_conservation_controls(species_count=len(species))
        taxonomy = _render_taxonomy_controls(
            species=species,
            authority=authority,
        )
        query_rows = motif_group_summary(
            resource=resource,
            sidecar_path=sidecar,
            motif=search.expression,
            search_mode=search.mode,
            threshold=filters["threshold"],
            minimum_species=filters["minimum_species"],
            minimum_proteins=filters["minimum_proteins"],
            analysis_species=taxonomy["analysis_species"],
            required_species=taxonomy["required_species"],
            group_type=group_type,
            hierarchy_node=hierarchy_node,
            maximum_rows=20_000 if authority is not None else filters["maximum_rows"],
        )
        rows = _apply_taxonomy_filters(
            rows=query_rows,
            authority=authority,
            minimum_lineage_fraction=filters["lineage_threshold"],
            required_taxon_ids=taxonomy["required_taxon_ids"],
            excluded_taxon_ids=taxonomy["excluded_taxon_ids"],
        )[: filters["maximum_rows"]]
    except InputValidationError as error:
        st.error(str(error))
        return

    st.caption(
        f"{len(rows):,} groups pass the current filters. Matching fractions use only "
        "successfully reconciled sequences; sequence coverage is reported separately."
    )
    if not rows:
        st.warning("No group passes the current sequence, conservation and taxonomy filters.")
        return
    _render_metrics(
        rows=rows,
        search=search,
        threshold=filters["threshold"],
        species_count=len(taxonomy["analysis_species"]),
    )
    roles = _species_roles(
        species=species,
        analysis_species=taxonomy["analysis_species"],
        authority=authority,
        required_taxon_ids=taxonomy["required_taxon_ids"],
        excluded_taxon_ids=taxonomy["excluded_taxon_ids"],
    )
    if expression_only:
        st.subheader("Explore RNA-seq evidence")
    if selected_view == _CANDIDATE_VIEW:
        _render_candidate_overview(
            rows=rows,
            search=search,
            threshold=filters["threshold"],
            group_type=group_type,
        )
    elif selected_view == _TAXONOMY_VIEW:
        _render_taxonomic_matrix(
            rows=rows,
            species=taxonomy["analysis_species"],
            search=search,
        )
    elif selected_view == _GROUP_VIEW:
        _render_selected_group(
            rows=rows,
            resource=resource,
            sidecar=sidecar,
            search=search,
            group_type=group_type,
            hierarchy_node=hierarchy_node,
            authority=authority,
            roles=roles,
        )
    else:
        _render_expression_evidence(
            rows=rows,
            resource=resource,
            group_type=group_type,
            hierarchy_node=hierarchy_node,
        )


def render_rna_seq_explorer_page(
    *,
    resource: ResourceIdentity,
    sidecar_path_text: str,
    taxonomy_path_text: str = "",
) -> None:
    """Render a direct, discoverable entry point to RNA-seq evidence.

    Args:
        resource: Validated immutable OrthoFinder resource.
        sidecar_path_text: Optional external sequence-sidecar path for older resources.
        taxonomy_path_text: Optional reviewed taxonomy authority path.
    """

    render_terminal_motif_page(
        resource=resource,
        sidecar_path_text=sidecar_path_text,
        taxonomy_path_text=taxonomy_path_text,
        expression_only=True,
    )


def _render_result_workspace_selector(*, resource: ResourceIdentity) -> str:
    """Render a prominent lazy workspace switch and return the selected view."""

    st.subheader("Choose a results workspace")
    st.caption(
        "Only the selected workspace is queried and rendered. Changing workspace does not "
        "alter the sequence, orthology or taxonomic filters below."
    )
    views = [_CANDIDATE_VIEW, _TAXONOMY_VIEW, _GROUP_VIEW]
    if expression_available(resource=resource):
        views.append(_EXPRESSION_VIEW)
        st.success(
            "RNA-seq explorer available — 12 focal plant species plus Homo sapiens."
        )
    else:
        st.info(
            "RNA-seq explorer is unavailable in this resource. Motif discovery remains fully "
            "available; open the dedicated RNA-seq page for the rebuild requirements."
        )
    selected = st.segmented_control(
        "Results workspace",
        options=views,
        default=_CANDIDATE_VIEW,
        selection_mode="single",
        key="terminal_motif_result_workspace",
        help=(
            "Choose the result type before adjusting the shared filters. The RNA-seq workspace "
            "contains the heatmap, UpSet intersections and complete evidence tables."
        ),
    )
    return str(selected or _CANDIDATE_VIEW)


def _render_expression_setup_state() -> None:
    """Explain why an older resource cannot expose the RNA-seq explorer."""

    st.warning(
        "The opened resource does not contain the complete RNA-seq evidence relations, so the "
        "explorer cannot run from this dataset."
    )
    st.markdown(
        "A schema-5 resource must contain `expression_context`, "
        "`expression_member_mapping`, `expression_member_summary` and "
        "`expression_group_summary`. This is an allowed missing resource: every non-expression "
        "page continues to work, and no missing expression value is interpreted as zero."
    )


def _render_setup_state() -> None:
    """Explain the one-time sequence-authority setup."""

    st.warning(
        "No group motif analysis has run in this view. The open DuckDB contains group "
        "memberships, but not the complete amino-acid sequences needed for motif search."
    )
    st.subheader("One-time setup")
    st.write(
        "Build the sequence sidecar from the same OrthoFinder Results_* directory and "
        "restart the viewer with `--terminal-motif-parquet`. This does not rerun "
        "OrthoFinder. Exact, regex and taxonomic controls appear only after it loads."
    )
    st.code(
        "orthofinder-terminal-motif-build \\\n"
        "  --orthofinder-results-dir /path/to/Results_Feb26 \\\n"
        "  --output-parquet /path/to/orthofinder_sequences.parquet",
        language="bash",
    )


def _render_search_controls() -> SequenceSearch:
    """Render exact-versus-regex controls and return a validated search."""

    st.subheader("Sequence search")
    use_regex = st.toggle(
        "Enable regular-expression search",
        value=False,
        help=(
            "Off by default. Exact C-terminal matching is used unless you explicitly "
            "enable regex mode."
        ),
    )
    if not use_regex:
        expression = st.text_input(
            "Exact C-terminal sequence",
            value=DEFAULT_MOTIF,
            help="Canonical one-letter amino-acid codes; no regex interpretation.",
        )
        return validate_sequence_search(expression=expression, mode=EXACT_C_TERMINAL)
    mode_label = st.radio(
        "Regular-expression location",
        ("Anywhere in the protein", "C-terminus only"),
        horizontal=True,
        help=(
            "Anywhere finds the first occurrence at any position. C-terminus only "
            "anchors the expression automatically; a trailing $ is optional."
        ),
    )
    expression = st.text_input(
        "Protein regular expression",
        value="",
        placeholder="Example: N[^P][ST]",
        help=(
            "Applied to upper-case amino-acid sequences with DuckDB's RE2 engine. "
            "Regular expressions are never enabled implicitly."
        ),
    )
    mode = REGEX_ANYWHERE if mode_label == "Anywhere in the protein" else REGEX_C_TERMINAL
    return validate_sequence_search(expression=expression, mode=mode)


def _render_group_controls(*, resource: ResourceIdentity) -> tuple[str, str]:
    """Render OrthoFinder grouping authority controls."""

    st.subheader("Orthology grouping")
    label = st.radio(
        "Grouping authority",
        ("Hierarchical orthogroups (HOGs)", "Original OrthoFinder orthogroups"),
        horizontal=True,
    )
    if label == "Original OrthoFinder orthogroups":
        return "LEGACY_ORTHOGROUP", ""
    nodes = motif_hierarchy_nodes(resource=resource)
    default_index = nodes.index("N0") if "N0" in nodes else 0
    hierarchy_node = st.selectbox(
        "HOG hierarchy node",
        options=nodes,
        index=default_index,
        help="N0 is the root HOG authority in this study; other available levels remain explicit.",
    )
    return "HOG", hierarchy_node


def _render_conservation_controls(*, species_count: int) -> dict[str, Any]:
    """Render bounded conservation and result-size controls."""

    controls = st.columns((1.4, 1.0, 1.0, 1.0))
    threshold = controls[0].slider(
        "Minimum matching analysis proteins",
        min_value=0,
        max_value=100,
        value=int(DEFAULT_THRESHOLD * 100),
        step=1,
    ) / 100.0
    minimum_proteins = int(
        controls[1].number_input(
            "Minimum analysis proteins", min_value=1, max_value=1_000_000, value=2, step=1
        )
    )
    minimum_species = int(
        controls[2].number_input(
            "Minimum analysis species",
            min_value=1,
            max_value=max(1, species_count),
            value=min(2, max(1, species_count)),
            step=1,
        )
    )
    maximum_rows = int(
        controls[3].number_input(
            "Maximum groups", min_value=10, max_value=20_000, value=2_000, step=10
        )
    )
    lineage_threshold = st.slider(
        "Minimum analysis species containing a match",
        min_value=0,
        max_value=100,
        value=80,
        step=1,
        help=(
            "Percentage of assessed primary species with at least one matching protein. "
            "This is separate from the protein-level threshold."
        ),
    ) / 100.0
    return {
        "threshold": threshold,
        "minimum_proteins": minimum_proteins,
        "minimum_species": minimum_species,
        "maximum_rows": maximum_rows,
        "lineage_threshold": lineage_threshold,
    }


def _render_taxonomy_controls(
    *, species: tuple[str, ...], authority: TaxonomyAuthority | None
) -> dict[str, Any]:
    """Render primary and comparison taxonomic controls."""

    st.subheader("Taxonomic conservation filters")
    required_species = tuple(
        st.multiselect(
            "Exact species that must contain a matching protein",
            options=species,
            help="Optional exact-label requirement, independent of broader lineage filters.",
        )
    )
    if authority is None:
        st.warning(
            "No reviewed taxonomy authority exactly matches these sequence labels. "
            "Lineage controls are disabled; exact species controls remain available."
        )
        selected = tuple(
            st.multiselect(
                "Primary species included in conservation calculations",
                options=species,
                default=species,
            )
        )
        return {
            "analysis_species": selected,
            "required_species": required_species,
            "required_taxon_ids": (),
            "excluded_taxon_ids": (),
        }
    options = authority.taxon_options()
    option_by_label = {option.display_label(): option for option in options}
    labels = tuple(option_by_label)
    focal_label = st.selectbox(
        "Primary lineage used for conservation calculations",
        options=("All reviewed sampled species", *labels),
        help=(
            "Only assessed descendants of this lineage enter the primary protein and "
            "species denominators. Other taxa remain comparison evidence."
        ),
    )
    if focal_label == "All reviewed sampled species":
        analysis_species = authority.reviewed_species
    else:
        analysis_species = authority.target_species(
            taxon_id=option_by_label[focal_label].taxon_id
        )
    with st.expander("Optional exact primary-species subset"):
        restrict_species = st.toggle("Restrict the primary lineage to selected species")
        if restrict_species:
            analysis_species = tuple(
                st.multiselect(
                    "Primary species",
                    options=analysis_species,
                    default=analysis_species,
                )
            )
    required_labels = st.multiselect(
        "Comparison lineages that must contain a matching descendant",
        options=labels,
        help="For example, require evidence in both a plant lineage and Homo.",
    )
    excluded_labels = st.multiselect(
        "Lineages that must not contain a matching descendant",
        options=labels,
        help="Optional negative-comparison filter; leave empty for discovery.",
    )
    return {
        "analysis_species": tuple(analysis_species),
        "required_species": required_species,
        "required_taxon_ids": tuple(
            option_by_label[label].taxon_id for label in required_labels
        ),
        "excluded_taxon_ids": tuple(
            option_by_label[label].taxon_id for label in excluded_labels
        ),
    }


def _render_metrics(
    *,
    rows: Sequence[Mapping[str, Any]],
    search: SequenceSearch,
    threshold: float,
    species_count: int,
) -> None:
    """Render compact result metrics."""

    metrics = st.columns(4)
    metrics[0].metric("Qualifying groups", f"{len(rows):,}")
    metrics[1].metric("Active search", search.expression)
    metrics[2].metric("Protein threshold", f"{threshold:.0%}")
    metrics[3].metric("Primary species", f"{species_count:,}")


def _render_candidate_overview(
    *,
    rows: Sequence[Mapping[str, Any]],
    search: SequenceSearch,
    threshold: float,
    group_type: str,
) -> None:
    """Render the candidate landscape and complete group table."""

    st.subheader("Conservation landscape")
    with st.expander("What this graph shows and how to interpret it"):
        st.write(
            "Each point is one group. The horizontal axis is the matching fraction among "
            "assessed proteins in the primary species set. The vertical axis is the number "
            "of primary species with a match. Point size reflects assessed proteins and "
            "colour shows sequence coverage. Strong candidates lie towards the upper right "
            "with high coverage. This is sequence evidence, not evidence of binding or degradation."
        )
    figure = _motif_figure(rows=rows, search=search, threshold=threshold)
    stem = _search_stem(search=search)
    render_plotly_figure(
        figure=figure,
        file_stem=f"{stem}_conservation_landscape",
        key="terminal_motif_conservation_figure",
        pdf_width=1600,
        pdf_height=1000,
    )
    display_rows = tuple(_display_summary(row=row) for row in rows)
    st.dataframe(display_rows, width="stretch", hide_index=True)
    render_table_downloads(
        records=display_rows,
        file_stem=f"{stem}_{group_type.lower()}_summary",
        key="terminal_motif_summary_tsv",
        column_definitions=_SUMMARY_HELP,
        workbook_title=f"Protein search {search.expression} group summary",
    )


def _render_taxonomic_matrix(
    *,
    rows: Sequence[Mapping[str, Any]],
    species: Sequence[str],
    search: SequenceSearch,
) -> None:
    """Render a bounded group-by-species evidence-state heatmap."""

    st.subheader("Group-by-species conservation matrix")
    maximum_groups = st.slider(
        "Groups shown in matrix", min_value=5, max_value=100, value=40, step=5
    )
    with st.expander("What this heatmap shows and how to interpret it"):
        st.write(
            "Rows are the highest-ranked passing groups and columns are primary species. "
            "Gold means at least one matching protein; blue means an assessed sequence but "
            "no match; grey means the species is represented but unassessed; white means it "
            "is not represented in that group. Look for broad gold blocks rather than a "
            "signal driven by one expanded species."
        )
    figure = _taxonomic_heatmap(
        rows=rows[:maximum_groups],
        species=species,
        search=search,
    )
    render_plotly_figure(
        figure=figure,
        file_stem=f"{_search_stem(search=search)}_taxonomic_matrix",
        key="terminal_motif_taxonomic_matrix",
        pdf_width=1900,
        pdf_height=max(900, 28 * min(maximum_groups, len(rows))),
    )


def _render_selected_group(
    *,
    rows: Sequence[Mapping[str, Any]],
    resource: ResourceIdentity,
    sidecar: Path | None,
    search: SequenceSearch,
    group_type: str,
    hierarchy_node: str,
    authority: TaxonomyAuthority | None,
    roles: Mapping[str, str],
) -> None:
    """Render species and protein evidence for one selected group."""

    st.subheader("Inspect one group")
    selected = st.selectbox(
        "Group",
        options=[str(row["group_id"]) for row in rows],
        key="terminal_motif_group",
    )
    members = motif_group_members(
        resource=resource,
        sidecar_path=sidecar,
        motif=search.expression,
        search_mode=search.mode,
        group_id=selected,
        group_type=group_type,
        hierarchy_node=hierarchy_node,
    )
    species_rows = _species_distribution(
        rows=members,
        authority=authority,
        roles=roles,
    )
    st.markdown("#### Species evidence")
    species_figure = _species_evidence_figure(rows=species_rows, group_id=selected)
    render_plotly_figure(
        figure=species_figure,
        file_stem=f"{selected}_{_search_stem(search=search)}_species_evidence",
        key="terminal_motif_species_evidence_figure",
        pdf_width=1800,
        pdf_height=1000,
    )
    st.dataframe(species_rows, width="stretch", hide_index=True)
    render_table_downloads(
        records=species_rows,
        file_stem=f"{selected}_{_search_stem(search=search)}_species_evidence",
        key="terminal_motif_taxonomy_tsv",
        column_definitions=_SPECIES_HELP,
        workbook_title=f"{selected} species motif evidence",
    )
    st.markdown("#### Protein-level calls")
    filtered_members, show_sequences = _filter_member_rows(rows=members, roles=roles)
    member_rows = tuple(
        _display_member(
            row=row,
            authority=authority,
            roles=roles,
            show_sequence=show_sequences,
        )
        for row in filtered_members
    )
    st.dataframe(member_rows, width="stretch", hide_index=True)
    export_rows = tuple(
        _display_member(
            row=row,
            authority=authority,
            roles=roles,
            show_sequence=True,
        )
        for row in filtered_members
    )
    render_table_downloads(
        records=export_rows,
        file_stem=f"{selected}_{_search_stem(search=search)}_protein_calls",
        key="terminal_motif_member_tsv",
        column_definitions=_MEMBER_HELP,
        workbook_title=f"{selected} protein motif calls",
    )
    st.download_button(
        "Download filtered protein sequences as FASTA",
        data=_members_to_fasta(rows=filtered_members),
        file_name=f"{selected}_{_search_stem(search=search)}_proteins.fasta",
        mime="text/x-fasta",
        key="terminal_motif_fasta",
    )
    if expression_available(resource=resource):
        st.markdown("#### RNA-seq mapping and expression status")
        with st.expander("How to interpret the RNA-seq evidence states"):
            st.write(
                "Mappings are exact and species scoped. MAPPED_UNIQUE means one Atlas "
                "gene matched at the best identifier tier. AMBIGUOUS and NOT_MAPPED are "
                "unavailable evidence, not negative expression. NO_EXPRESSION_RECORDS "
                "means a gene mapped but no compatible Atlas context was published. TPM "
                "and FPKM are never combined within an experiment. The production RNA-seq "
                "scope is the reviewed 12-plant panel plus Homo sapiens. Proteins from other "
                "OrthoFinder species are not assessed here; they are not negative mappings."
            )
        expression_rows = expression_member_evidence(
            resource=resource,
            group_type=group_type,
            hierarchy_node=hierarchy_node,
            group_id=selected,
        )
        st.dataframe(expression_rows, width="stretch", hide_index=True)
        render_table_downloads(
            records=expression_rows,
            file_stem=f"{selected}_rna_seq_member_evidence",
            key="terminal_motif_expression_members",
            workbook_title=f"{selected} RNA-seq member evidence",
        )


def _render_expression_evidence(
    *,
    rows: Sequence[Mapping[str, Any]],
    resource: ResourceIdentity,
    group_type: str,
    hierarchy_node: str,
) -> None:
    """Render lazy, unit-safe heatmap and cross-species UpSet evidence."""

    st.subheader("RNA-seq expression evidence")
    with st.expander("What these analyses show and their limitations"):
        st.write(
            "This section asks whether proteins in the motif-qualified groups have "
            "species-scoped Expression Atlas evidence, and in which tissues or contexts. "
            "The production evidence panel contains 12 focal plant species plus Homo sapiens; "
            "the remaining OrthoFinder species stay available in every non-expression view. "
            "It does not test whether the motif causes expression or whether Cereblon "
            "regulates transcript abundance. Protein accumulation is post-transcriptional, "
            "so RNA-seq is supporting biological context rather than a substitute for "
            "protein-level validation. Missing, ambiguous and unmapped evidence stays blank."
        )
    dimensions = expression_dimensions(resource=resource)
    units = dimensions["units"]
    if not units:
        st.info("The resource contains the expression schema but no mapped context rows.")
        return
    available_groups = [str(row["group_id"]) for row in rows]
    default_groups = available_groups[: min(12, len(available_groups))]
    selected_groups = tuple(
        st.multiselect(
            "Groups compared (maximum 25)",
            options=available_groups,
            default=default_groups,
            help="The heatmap is deliberately bounded so labels and PDF exports remain readable.",
        )
    )
    if not selected_groups:
        st.info("Select at least one group to display RNA-seq evidence.")
        return
    if len(selected_groups) > 25:
        st.error("Select no more than 25 groups for the expression heatmap.")
        return
    controls = st.columns((1.2, 1.0, 1.4, 0.8))
    context_label = controls[0].selectbox(
        "Heatmap context",
        options=tuple(CONTEXT_COLUMNS),
        index=1 if "Organism part / tissue" in CONTEXT_COLUMNS else 0,
    )
    default_unit = units.index("TPM") if "TPM" in units else 0
    unit = controls[1].selectbox("Expression unit", options=units, index=default_unit)
    selected_species = tuple(
        controls[2].multiselect(
            "Expression species",
            options=dimensions["species"],
            default=dimensions["species"],
            help="An empty selection means all expression-bearing species.",
        )
    )
    effective_species = selected_species or tuple(dimensions["species"])
    log_transform = controls[3].toggle("log2(1 + value)", value=True)
    species_rows = expression_group_species(
        resource=resource,
        group_type=group_type,
        hierarchy_node=hierarchy_node,
        group_ids=selected_groups,
        species=effective_species,
    )
    _render_expression_metrics(rows=species_rows, group_count=len(selected_groups))
    st.caption(
        f"Every selected-evidence table below is restricted to {len(effective_species):,} "
        "chosen expression species and the selected expression unit."
    )
    heatmap_cells = expression_heatmap_cells(
        resource=resource,
        group_type=group_type,
        hierarchy_node=hierarchy_node,
        group_ids=selected_groups,
        context_column=CONTEXT_COLUMNS[context_label],
        expression_unit=unit,
        species=effective_species,
    )
    st.markdown("#### Cross-species expression heatmap")
    if not heatmap_cells:
        st.info("No mapped expression contexts match the current unit and species selection.")
    else:
        heatmap = _expression_heatmap_figure(
            cells=heatmap_cells,
            selected_groups=selected_groups,
            log_transform=log_transform,
        )
        render_plotly_figure(
            figure=heatmap,
            file_stem="motif_group_rna_seq_expression_heatmap",
            key="terminal_motif_expression_heatmap",
            pdf_width=2100,
            pdf_height=max(1000, 52 * len(selected_groups)),
        )

    st.markdown("#### Cross-species expression-evidence intersections")
    ranked_species = _rank_expression_species(rows=species_rows)
    upset_species = tuple(
        st.multiselect(
            "Species included in UpSet intersections (maximum 8)",
            options=ranked_species,
            default=ranked_species[: min(6, len(ranked_species))],
            help=(
                "A group belongs to a species set when at least one member has one or more "
                "compatible expression contexts."
            ),
        )
    )
    if len(upset_species) > 8:
        st.error("Select no more than eight species for a readable UpSet plot.")
    elif upset_species:
        intersections = _expression_intersections(
            rows=species_rows,
            group_ids=selected_groups,
            species=upset_species,
        )
        upset = _expression_upset_figure(
            intersections=intersections,
            species=upset_species,
        )
        render_plotly_figure(
            figure=upset,
            file_stem="motif_group_rna_seq_species_upset",
            key="terminal_motif_expression_upset",
            pdf_width=1900,
            pdf_height=1200,
        )
        st.dataframe(intersections, width="stretch", hide_index=True)
        render_table_downloads(
            records=intersections,
            file_stem="motif_group_rna_seq_species_intersections",
            key="terminal_motif_expression_intersections",
            workbook_title="RNA-seq species evidence intersections",
        )

    st.markdown("#### Selected RNA-seq evidence tables")
    st.write(
        "These exports follow the selected groups, expression species and expression unit. "
        "The aggregated table contains the exact values plotted in the heatmap; the member "
        "table retains uniquely mapped, ambiguous, unmapped and unavailable evidence states."
    )
    if heatmap_cells:
        st.markdown("##### Aggregated heatmap cells")
        st.dataframe(heatmap_cells, width="stretch", hide_index=True)
        render_table_downloads(
            records=heatmap_cells,
            file_stem="selected_rna_seq_heatmap_cells",
            key="terminal_motif_expression_heatmap_cells",
            workbook_title="Selected RNA-seq heatmap cells",
        )

    st.markdown("##### Group-by-species mapping and expression coverage")
    if species_rows:
        st.dataframe(species_rows, width="stretch", hide_index=True)
        render_table_downloads(
            records=species_rows,
            file_stem="selected_rna_seq_group_species_summary",
            key="terminal_motif_expression_species_summary",
            workbook_title="Selected RNA-seq group-by-species summary",
        )
    else:
        st.info("No group members occur in the selected expression species.")

    member_result = expression_selected_member_evidence(
        resource=resource,
        group_type=group_type,
        hierarchy_node=hierarchy_node,
        group_ids=selected_groups,
        species=effective_species,
    )
    st.markdown("##### Protein mapping and missingness states")
    if member_result.truncated:
        st.warning(
            f"The member evidence exceeded {member_result.maximum_rows:,} rows. Narrow the "
            "group or species selection before treating the export as complete."
        )
    if member_result.rows:
        preview_rows = member_result.rows[:2_000]
        st.caption(
            f"Showing {len(preview_rows):,} of {len(member_result.rows):,} retained member "
            "rows on screen; both downloads contain every retained row."
        )
        st.dataframe(preview_rows, width="stretch", hide_index=True)
        render_table_downloads(
            records=member_result.rows,
            file_stem="selected_rna_seq_member_mapping_and_missingness",
            key="terminal_motif_expression_selected_members",
            workbook_title="Selected RNA-seq member evidence",
        )
    else:
        st.info("No protein member evidence matches the selected groups and species.")

    st.markdown("##### Underlying context-level expression records")
    st.write(
        "Load the exact Expression Atlas rows underlying the current unit, group and species "
        "selection. This is optional because context tables can be much larger than the "
        "aggregated heatmap. Missing and unmapped proteins remain in the member table above; "
        "only proteins with observed context rows occur here."
    )
    raw_controls = st.columns((1.4, 1.0))
    load_context_rows = raw_controls[0].toggle(
        "Load underlying context records",
        value=False,
        key="terminal_motif_load_expression_context",
    )
    maximum_context_rows = int(
        raw_controls[1].selectbox(
            "Maximum context rows",
            options=(10_000, 25_000, 50_000, 100_000),
            index=1,
            disabled=not load_context_rows,
            help="Narrow the groups or species if this protective bound is reached.",
        )
    )
    if load_context_rows:
        context_result = expression_context_records(
            resource=resource,
            group_type=group_type,
            hierarchy_node=hierarchy_node,
            group_ids=selected_groups,
            expression_unit=unit,
            species=effective_species,
            maximum_rows=maximum_context_rows,
        )
        if context_result.truncated:
            st.warning(
                f"More than {context_result.maximum_rows:,} context rows matched. The table "
                "and downloads are truncated; narrow the groups/species or increase the "
                "protective row limit before treating the export as complete."
            )
        if context_result.rows:
            preview_rows = context_result.rows[:2_000]
            st.caption(
                f"Showing {len(preview_rows):,} of {len(context_result.rows):,} retained "
                "context rows on screen; both downloads contain every retained row."
            )
            st.dataframe(preview_rows, width="stretch", hide_index=True)
            render_table_downloads(
                records=context_result.rows,
                file_stem="selected_rna_seq_context_records",
                key="terminal_motif_expression_context_records",
                workbook_title="Selected RNA-seq context records",
            )
        else:
            st.info("No observed context records match the current selection.")

    st.markdown("#### Whole-panel group summary")
    st.caption(
        "This packaged group-level authority covers the complete assessed RNA-seq panel. "
        "Use the selected tables above for a species-restricted export."
    )
    summaries = expression_group_summaries(
        resource=resource,
        group_type=group_type,
        hierarchy_node=hierarchy_node,
        group_ids=selected_groups,
    )
    st.dataframe(summaries, width="stretch", hide_index=True)
    render_table_downloads(
        records=summaries,
        file_stem="motif_group_rna_seq_summary",
        key="terminal_motif_expression_group_summary",
        workbook_title="Motif-group RNA-seq summary",
    )


def _render_expression_metrics(
    *, rows: Sequence[Mapping[str, Any]], group_count: int
) -> None:
    """Render expression mapping and observation totals for selected groups."""

    member_count = sum(int(row["member_count"]) for row in rows)
    mapped = sum(int(row["mapped_member_count"]) for row in rows)
    observed = sum(int(row["expression_observed_member_count"]) for row in rows)
    metrics = st.columns(4)
    metrics[0].metric("Groups", f"{group_count:,}")
    metrics[1].metric("Group members", f"{member_count:,}")
    metrics[2].metric(
        "Uniquely mapped",
        f"{mapped:,}",
        help="Exact best-tier, species-scoped mappings to one Atlas gene.",
    )
    metrics[3].metric(
        "With expression contexts",
        f"{observed:,}",
        help="Mapped proteins with at least one compatible Atlas context row.",
    )


def _expression_heatmap_figure(
    *,
    cells: Sequence[Mapping[str, Any]],
    selected_groups: Sequence[str],
    log_transform: bool,
) -> go.Figure:
    """Build a group-by-species/context expression heatmap with blank missing cells."""

    contexts = tuple(
        dict.fromkeys(
            f"{row['species_label']} — {row['context_label']}" for row in cells
        )
    )
    lookup = {
        (str(row["group_id"]), f"{row['species_label']} — {row['context_label']}"): row
        for row in cells
    }
    values = []
    custom = []
    for group_id in selected_groups:
        value_row = []
        custom_row = []
        for context in contexts:
            row = lookup.get((group_id, context))
            if row is None:
                value_row.append(None)
                custom_row.append([None, 0, 0, None])
                continue
            raw = float(row["median_expression"])
            value_row.append(math.log2(1.0 + max(0.0, raw)) if log_transform else raw)
            custom_row.append(
                [
                    raw,
                    int(row["mapped_member_count"]),
                    int(row["experiment_count"]),
                    float(row["positive_context_fraction"]),
                ]
            )
        values.append(value_row)
        custom.append(custom_row)
    unit = str(cells[0]["expression_unit"])
    colour_title = f"log2(1 + {unit})" if log_transform else unit
    figure = go.Figure(
        go.Heatmap(
            z=values,
            x=contexts,
            y=list(selected_groups),
            customdata=custom,
            colorscale=((0.0, "#ffffff"), (0.50, "#fcae91"), (1.0, "#cb181d")),
            colorbar={"title": colour_title},
            hoverongaps=False,
            hovertemplate=(
                "Group=%{y}<br>Species / context=%{x}<br>"
                f"Median {unit}=%{{customdata[0]:.4g}}<br>"
                "Mapped proteins=%{customdata[1]}<br>Experiments=%{customdata[2]}<br>"
                "Positive-context fraction=%{customdata[3]:.1%}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title="Median RNA-seq expression by OrthoFinder group and biological context",
        xaxis_title="Species and biological context",
        yaxis_title="OrthoFinder group",
        yaxis={"autorange": "reversed"},
        xaxis={"tickangle": -45},
        template="plotly_white",
        height=max(700, 40 * len(selected_groups)),
    )
    return figure


def _rank_expression_species(*, rows: Sequence[Mapping[str, Any]]) -> tuple[str, ...]:
    """Rank species by the number of selected groups with observed expression."""

    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        if int(row["expression_observed_member_count"]) > 0:
            counts[str(row["species_label"])] += 1
    return tuple(sorted(counts, key=lambda value: (-counts[value], value)))


def _expression_intersections(
    *,
    rows: Sequence[Mapping[str, Any]],
    group_ids: Sequence[str],
    species: Sequence[str],
) -> tuple[dict[str, Any], ...]:
    """Return exact UpSet combinations of species with observed RNA-seq evidence."""

    observed = {
        (str(row["group_id"]), str(row["species_label"]))
        for row in rows
        if int(row["expression_observed_member_count"]) > 0
    }
    groups_by_combination: dict[tuple[str, ...], list[str]] = defaultdict(list)
    for group_id in group_ids:
        combination = tuple(label for label in species if (group_id, label) in observed)
        groups_by_combination[combination].append(group_id)
    ordered = sorted(
        groups_by_combination.items(),
        key=lambda item: (-len(item[1]), -len(item[0]), item[0]),
    )
    return tuple(
        {
            "Intersection": index,
            "Species with expression evidence": "; ".join(combination) or "None selected",
            "Species count": len(combination),
            "Group count": len(groups),
            "Groups": "; ".join(sorted(groups)),
        }
        for index, (combination, groups) in enumerate(ordered, start=1)
    )


def _expression_upset_figure(
    *, intersections: Sequence[Mapping[str, Any]], species: Sequence[str]
) -> go.Figure:
    """Build an UpSet-style bar-and-membership-matrix figure."""

    x_values = [int(row["Intersection"]) for row in intersections]
    figure = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=(0.46, 0.54),
        vertical_spacing=0.04,
    )
    figure.add_bar(
        x=x_values,
        y=[int(row["Group count"]) for row in intersections],
        marker_color="#147d78",
        hovertemplate="Intersection %{x}<br>Groups=%{y}<extra></extra>",
        showlegend=False,
        row=1,
        col=1,
    )
    selected_by_x = {
        int(row["Intersection"]): set(
            str(row["Species with expression evidence"]).split("; ")
        )
        for row in intersections
    }
    for x_value in x_values:
        selected = selected_by_x[x_value]
        active_indices = [index for index, label in enumerate(species) if label in selected]
        if len(active_indices) > 1:
            figure.add_trace(
                go.Scatter(
                    x=[x_value, x_value],
                    y=[min(active_indices), max(active_indices)],
                    mode="lines",
                    line={"color": "#374151", "width": 2},
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=2,
                col=1,
            )
        figure.add_trace(
            go.Scatter(
                x=[x_value] * len(species),
                y=list(range(len(species))),
                mode="markers",
                marker={
                    "size": 11,
                    "color": [
                        "#147d78" if label in selected else "#d1d5db" for label in species
                    ],
                },
                customdata=[[label, label in selected] for label in species],
                hovertemplate=(
                    "Intersection=%{x}<br>Species=%{customdata[0]}<br>"
                    "Expression evidence=%{customdata[1]}<extra></extra>"
                ),
                showlegend=False,
            ),
            row=2,
            col=1,
        )
    figure.update_yaxes(title_text="Groups", row=1, col=1)
    figure.update_yaxes(
        tickmode="array",
        tickvals=list(range(len(species))),
        ticktext=[label.replace("_", " ") for label in species],
        autorange="reversed",
        title_text="Species",
        row=2,
        col=1,
    )
    figure.update_xaxes(title_text="Exact expression-evidence intersection", row=2, col=1)
    figure.update_layout(
        title="UpSet intersections of observed RNA-seq evidence across species",
        template="plotly_white",
        height=max(720, 55 * len(species) + 420),
        bargap=0.25,
    )
    return figure


def _filter_member_rows(
    *, rows: Sequence[Mapping[str, Any]], roles: Mapping[str, str]
) -> tuple[tuple[Mapping[str, Any], ...], bool]:
    """Render member filters and return the selected rows."""

    controls = st.columns((1.0, 1.2, 1.4))
    state = controls[0].selectbox(
        "Match state",
        ("All proteins", "Matching", "Not matching", "Sequence unavailable"),
    )
    role_options = tuple(sorted(set(roles.values())))
    selected_roles = tuple(
        controls[1].multiselect("Analysis role", options=role_options, default=role_options)
    )
    species_options = tuple(sorted({str(row["species_label"]) for row in rows}))
    selected_species = tuple(
        controls[2].multiselect("Species", options=species_options, default=species_options)
    )
    identifier_text = st.text_input(
        "Protein identifier contains",
        value="",
        help="Literal case-insensitive filter over the original protein identifier.",
    ).strip().casefold()
    show_sequences = st.toggle("Show complete sequences in the on-screen table", value=False)
    filtered = []
    for row in rows:
        available = bool(row["sequence_available"])
        matched = bool(row["motif_match"])
        species = str(row["species_label"])
        if state == "Matching" and not matched:
            continue
        if state == "Not matching" and (matched or not available):
            continue
        if state == "Sequence unavailable" and available:
            continue
        if roles.get(species, "UNRESOLVED") not in selected_roles:
            continue
        if species not in selected_species:
            continue
        if identifier_text and identifier_text not in str(row["member_id"]).casefold():
            continue
        filtered.append(row)
    st.caption(f"{len(filtered):,} of {len(rows):,} protein rows shown.")
    return tuple(filtered), show_sequences


def _load_taxonomy(
    *, species: tuple[str, ...], taxonomy_path_text: str
) -> TaxonomyAuthority | None:
    """Load an exact reviewed taxonomy authority without guessing labels."""

    if taxonomy_path_text.strip():
        return read_taxonomy_mapping(path=Path(taxonomy_path_text), expected_species=species)
    return read_matching_bundled_taxonomy(expected_species=species)


def _apply_taxonomy_filters(
    *,
    rows: Sequence[Mapping[str, Any]],
    authority: TaxonomyAuthority | None,
    minimum_lineage_fraction: float,
    required_taxon_ids: Sequence[int],
    excluded_taxon_ids: Sequence[int],
) -> tuple[dict[str, Any], ...]:
    """Annotate and filter summaries using assessed primary species."""

    required = (
        [set(authority.target_species(taxon_id=value)) for value in required_taxon_ids]
        if authority is not None
        else []
    )
    excluded = (
        [set(authority.target_species(taxon_id=value)) for value in excluded_taxon_ids]
        if authority is not None
        else []
    )
    result: list[dict[str, Any]] = []
    for source in rows:
        assessed = _species_set(source.get("analysis_assessed_species"))
        analysis_matching = _species_set(source.get("analysis_matching_species"))
        all_matching = _species_set(source.get("matching_species"))
        fraction = len(analysis_matching) / len(assessed) if assessed else 0.0
        if not assessed or fraction < minimum_lineage_fraction:
            continue
        if any(not all_matching.intersection(target) for target in required):
            continue
        if any(all_matching.intersection(target) for target in excluded):
            continue
        row = dict(source)
        row["analysis_species_match_fraction"] = fraction
        result.append(row)
    return tuple(result)


def _species_roles(
    *,
    species: Sequence[str],
    analysis_species: Sequence[str],
    authority: TaxonomyAuthority | None,
    required_taxon_ids: Sequence[int],
    excluded_taxon_ids: Sequence[int],
) -> dict[str, str]:
    """Assign mutually exclusive roles for result display."""

    primary = set(analysis_species)
    required = set()
    excluded = set()
    if authority is not None:
        for taxon_id in required_taxon_ids:
            required.update(authority.target_species(taxon_id=taxon_id))
        for taxon_id in excluded_taxon_ids:
            excluded.update(authority.target_species(taxon_id=taxon_id))
    roles = {}
    for label in species:
        if label in primary:
            role = "PRIMARY_ANALYSIS"
        elif label in required:
            role = "REQUIRED_COMPARISON"
        elif label in excluded:
            role = "EXCLUDED_COMPARISON"
        else:
            role = "OTHER_SAMPLED"
        roles[label] = role
    return roles


def _species_distribution(
    *,
    rows: Sequence[Mapping[str, Any]],
    authority: TaxonomyAuthority | None,
    roles: Mapping[str, str],
) -> tuple[dict[str, Any], ...]:
    """Summarise motif calls and missingness by species for one group."""

    records = (
        {row.workflow_species_label: row for row in authority.reviewed_records}
        if authority is not None
        else {}
    )
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row["species_label"]), []).append(row)
    output = []
    for species, members in sorted(grouped.items()):
        record = records.get(species)
        assessed = sum(bool(row["sequence_available"]) for row in members)
        matches = sum(bool(row["motif_match"]) for row in members)
        unavailable = len(members) - assessed
        if assessed == 0:
            result = "Sequence unavailable"
        elif matches == 0:
            result = "Assessed; no match"
        elif matches == assessed:
            result = "All assessed proteins match"
        else:
            result = "Some assessed proteins match"
        output.append(
            {
                "OrthoFinder species": species,
                "Accepted species": record.accepted_species_name if record else "Unresolved",
                "NCBI taxon ID": record.ncbi_taxon_id if record else "",
                "Analysis role": roles.get(species, "UNRESOLVED"),
                "Published proteins": len(members),
                "Assessed": assessed,
                "Matching": matches,
                "Unavailable": unavailable,
                "Match fraction": matches / assessed if assessed else None,
                "Species result": result,
            }
        )
    return tuple(output)


def _display_summary(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return one readable group summary row."""

    return {
        "Group": row["group_id"],
        "Members": row["total_member_count"],
        "Sequences assessed": row["sequence_count"],
        "Sequence coverage": row["sequence_coverage"],
        "Analysis proteins": row["analysis_sequence_count"],
        "Matching analysis proteins": row["analysis_matching_sequence_count"],
        "Analysis protein match": row["analysis_matching_fraction"],
        "Analysis species assessed": row["analysis_assessed_species_count"],
        "Analysis species matching": row["analysis_matching_species_count"],
        "Analysis species match": row["analysis_species_match_fraction"],
        "All matching species": row["matching_species"] or "",
    }


def _display_member(
    *,
    row: Mapping[str, Any],
    authority: TaxonomyAuthority | None,
    roles: Mapping[str, str],
    show_sequence: bool,
) -> dict[str, Any]:
    """Return one readable protein-level sequence call."""

    species = str(row["species_label"])
    records = (
        {item.workflow_species_label: item for item in authority.reviewed_records}
        if authority is not None
        else {}
    )
    record = records.get(species)
    result = {
        "Accepted species": record.accepted_species_name if record else "Unresolved",
        "OrthoFinder species": species,
        "NCBI taxon ID": record.ncbi_taxon_id if record else "",
        "Analysis role": roles.get(species, "UNRESOLVED"),
        "Protein ID": row["member_id"],
        "Internal ID": row["internal_id"] or "",
        "Source FASTA": row.get("source_fasta") or "",
        "Protein description": _protein_description(
            raw_header=str(row.get("raw_header") or ""),
            member_id=str(row["member_id"]),
        ),
        "Raw FASTA header": row.get("raw_header") or "",
        "Sequence available": row["sequence_available"],
        "Sequence length": row["sequence_length"] or "",
        "Matches search": row["motif_match"],
        "Matched sequence": row["matched_sequence"] or "",
    }
    if show_sequence:
        result["Sequence"] = row["sequence"] or ""
    return result


def _protein_description(*, raw_header: str, member_id: str) -> str:
    """Return the source description after the first FASTA identifier token."""

    text = raw_header.strip()
    if not text:
        return ""
    first, separator, description = text.partition(" ")
    if separator and first in {member_id, member_id.split(maxsplit=1)[0]}:
        return description.strip()
    return description.strip() if separator else ""


def _motif_figure(
    *, rows: Sequence[Mapping[str, Any]], search: SequenceSearch, threshold: float
) -> go.Figure:
    """Build the protein-fraction versus species-breadth landscape."""

    figure = go.Figure(
        go.Scatter(
            x=[float(row["analysis_matching_fraction"]) for row in rows],
            y=[int(row["analysis_matching_species_count"]) for row in rows],
            mode="markers",
            customdata=[
                [
                    row["group_id"],
                    row["analysis_sequence_count"],
                    row["analysis_assessed_species_count"],
                    row["sequence_coverage"],
                ]
                for row in rows
            ],
            marker={
                "size": [
                    max(7, min(30, 5 + int(row["analysis_sequence_count"]) ** 0.5))
                    for row in rows
                ],
                "color": [float(row["sequence_coverage"]) for row in rows],
                "colorscale": "Tealgrn",
                "cmin": 0,
                "cmax": 1,
                "showscale": True,
                "colorbar": {"title": "Sequence coverage", "tickformat": ".0%"},
                "opacity": 0.78,
            },
            hovertemplate=(
                "Group=%{customdata[0]}<br>Analysis protein match=%{x:.1%}<br>"
                "Matching analysis species=%{y}<br>Analysis sequences=%{customdata[1]}<br>"
                "Assessed analysis species=%{customdata[2]}<br>"
                "Overall sequence coverage=%{customdata[3]:.1%}<extra></extra>"
            ),
        )
    )
    figure.add_vline(x=threshold, line_dash="dash", line_color="#d95f02")
    figure.update_layout(
        title=f"Group conservation of {search.display_label}",
        xaxis_title="Fraction of assessed primary proteins matching",
        yaxis_title="Primary species with at least one matching protein",
        xaxis={"tickformat": ".0%", "range": [0, 1.01]},
        template="plotly_white",
        height=700,
    )
    return figure


def _taxonomic_heatmap(
    *,
    rows: Sequence[Mapping[str, Any]],
    species: Sequence[str],
    search: SequenceSearch,
) -> go.Figure:
    """Build a group-by-species categorical evidence heatmap."""

    selected_species = tuple(species)
    matrix = []
    hover = []
    for row in rows:
        represented = _species_set(row.get("represented_species"))
        assessed = _species_set(row.get("assessed_species"))
        matching = _species_set(row.get("matching_species"))
        values = []
        labels = []
        for label in selected_species:
            if label in matching:
                values.append(3)
                labels.append("Matching")
            elif label in assessed:
                values.append(2)
                labels.append("Assessed; no match")
            elif label in represented:
                values.append(1)
                labels.append("Sequence unavailable")
            else:
                values.append(0)
                labels.append("Not represented")
        matrix.append(values)
        hover.append(labels)
    figure = go.Figure(
        go.Heatmap(
            z=matrix,
            x=selected_species,
            y=[str(row["group_id"]) for row in rows],
            customdata=hover,
            zmin=0,
            zmax=3,
            colorscale=[
                [0.00, "#ffffff"],
                [0.24, "#ffffff"],
                [0.25, "#b8b8b8"],
                [0.49, "#b8b8b8"],
                [0.50, "#2677a8"],
                [0.74, "#2677a8"],
                [0.75, "#f2b134"],
                [1.00, "#f2b134"],
            ],
            colorbar={
                "title": "Evidence",
                "tickvals": [0, 1, 2, 3],
                "ticktext": ["Absent", "Unavailable", "No match", "Match"],
            },
            hovertemplate=(
                "Group=%{y}<br>Species=%{x}<br>Status=%{customdata}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title=f"Taxonomic distribution of {search.display_label}",
        xaxis_title="Primary species",
        yaxis_title="OrthoFinder group",
        template="plotly_white",
        height=max(650, 24 * len(rows)),
        xaxis={"tickangle": -45},
    )
    return figure


def _species_evidence_figure(
    *, rows: Sequence[Mapping[str, Any]], group_id: str
) -> go.Figure:
    """Build a stacked assessed/matching/unavailable chart for one group."""

    labels = [str(row["Accepted species"]) for row in rows]
    matching = [int(row["Matching"]) for row in rows]
    assessed_no_match = [int(row["Assessed"]) - int(row["Matching"]) for row in rows]
    unavailable = [int(row["Unavailable"]) for row in rows]
    figure = go.Figure()
    figure.add_bar(name="Matching", x=labels, y=matching, marker_color="#d9a21b")
    figure.add_bar(
        name="Assessed; no match", x=labels, y=assessed_no_match, marker_color="#2677a8"
    )
    figure.add_bar(name="Unavailable", x=labels, y=unavailable, marker_color="#a8a8a8")
    figure.update_layout(
        title=f"Species-level protein evidence for {group_id}",
        xaxis_title="Species",
        yaxis_title="Proteins",
        barmode="stack",
        template="plotly_white",
        height=650,
        xaxis={"tickangle": -45},
    )
    return figure


def _search_stem(*, search: SequenceSearch) -> str:
    """Return a stable portable search label for filenames."""

    mode = search.mode.lower()
    expression = "".join(
        character if character.isalnum() else "_"
        for character in search.expression
    )
    return f"{mode}_{expression[:50].strip('_') or 'pattern'}"


def _species_set(value: Any) -> set[str]:
    """Parse a semicolon-delimited species aggregation."""

    return {item for item in str(value or "").split("; ") if item}


def _members_to_fasta(*, rows: Sequence[Mapping[str, Any]]) -> str:
    """Serialise selected group members with available sequences as wrapped FASTA."""

    chunks: list[str] = []
    for row in rows:
        sequence = str(row.get("sequence") or "")
        if not sequence:
            continue
        header = (
            f">{row['member_id']} species={row['species_label']} "
            f"sequence_match={row['motif_match']}"
        )
        wrapped = (sequence[index : index + 80] for index in range(0, len(sequence), 80))
        chunks.extend((header, *wrapped))
    return "\n".join(chunks) + ("\n" if chunks else "")

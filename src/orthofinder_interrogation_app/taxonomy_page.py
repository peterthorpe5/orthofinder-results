"""Streamlit page for reviewed, descendant-aware taxonomic group searches."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import streamlit as st

from orthofinder_results.errors import InputValidationError

from .evolutionary_page import (
    _comparison_keys,
    _store_active_group,
    _store_comparison_keys,
)
from .exports import render_table_downloads
from .models import GroupKey, TaxonomySearchFilters
from .queries import OrthoFinderQueryService
from .taxonomy import (
    TAXONOMY_COLUMNS,
    TaxonomyAuthority,
    parse_taxonomy_mapping,
    read_taxonomy_mapping,
    taxonomy_audit_rows,
    taxonomy_template_rows,
)

_LOGGER = logging.getLogger("orthofinder_interrogation_app.taxonomy_page")
_MODE_LABELS = {
    "Contains target descendants": "CONTAINS",
    "Enriched in target descendants": "ENRICHED",
    "Exclusive within sampled analysis": "SAMPLED_EXCLUSIVE",
    "Near-exclusive within sampled analysis": "NEAR_EXCLUSIVE",
}
_TAXONOMY_COLUMN_HELP = {
    "Group ID": "Exact OrthoFinder group identifier at the selected hierarchy level.",
    "Proteins in group": "Complete protein-member count for this group record.",
    "Species represented": "All sampled species contributing at least one protein.",
    "Target descendants represented": "Reviewed target-descendant species found in the group.",
    "Target descendant coverage": (
        "Represented reviewed target descendants divided by all reviewed target descendants "
        "sampled in this run."
    ),
    "Mapped target purity": (
        "Represented target descendants divided by all represented REVIEWED species."
    ),
    "Reviewed outsiders": "Represented REVIEWED species outside the selected target lineage.",
    "Outsider species": "Exact labels of represented reviewed outsiders.",
    "Unresolved species": "Represented labels that are unmapped, ambiguous or awaiting review.",
    "Stored distance coverage": "Whether this group already has a stored distance matrix.",
    "Stored mean pair distance": "Mean of stored pair distances when available.",
    "Stored median pair distance": "Median of stored pair distances when available.",
    "Stored distance spread (SD)": "Population standard deviation of stored pair distances.",
    "Enrichment odds ratio": "One-sided species-presence enrichment effect estimate.",
    "Fisher p-value": "Unadjusted one-sided Fisher exact-test p-value.",
    "Adjusted q-value (BH)": (
        "Benjamini–Hochberg-adjusted p-value across every tested group in the selected authority."
    ),
}
_SPECIES_RESULT_HELP = {
    "Species": "Exact species label represented in the selected group.",
    "Proteins from species": "Number of group proteins contributed by this species.",
    "Share of group": "Fraction of the group's proteins contributed by this species.",
}
_TAXONOMY_TEMPLATE_HELP = {
    "workflow_species_label": "Exact species label used by this OrthoFinder resource.",
    "source_species_name": "Species text inferred only for manual review, not accepted taxonomy.",
    "accepted_species_name": "Human-reviewed accepted scientific name.",
    "ncbi_taxon_id": "Human-reviewed positive NCBI taxonomy identifier for the species.",
    "parent_taxon_id": "NCBI taxonomy identifier of the accepted immediate parent.",
    "parent_taxon_name": "Accepted name of the immediate parent taxon.",
    "lineage_taxon_ids": "Semicolon-separated ordered lineage of NCBI taxonomy identifiers.",
    "lineage_names": "Semicolon-separated lineage names in the same order as lineage IDs.",
    "mapping_status": "REVIEWED, PENDING_REVIEW, UNMAPPED or AMBIGUOUS.",
    "mapping_method": "Method used to propose or confirm this mapping.",
    "mapping_source": "Authoritative taxonomy source or database.",
    "source_date": "Date on which the mapping source was accessed.",
    "source_version": "Version or release identifier of the taxonomy source.",
    "reviewed_by": "Person who accepted the mapping.",
    "reviewed_at_utc": "UTC date and time at which the mapping was accepted.",
    "review_note": "Free-text rationale, ambiguity or review action.",
}


def render_taxonomy_search(*, service: OrthoFinderQueryService, taxonomy_path_text: str) -> None:
    """Render mapping audit, target selection and four taxonomic searches."""

    st.header("Taxonomic search")
    st.caption(
        "Find groups that contain, favour or are restricted to descendants of a reviewed "
        "taxon. The mapping is generated from the current dataset and is never fixed to one "
        "study's species list."
    )
    with st.expander("What the four taxonomic searches mean", expanded=False):
        st.markdown(
            """
            - **Contains target descendants:** at least the requested number of reviewed target
              species is represented; outsiders are allowed.
            - **Enriched in target descendants:** target species occur more often than expected
              in the selected group collection, with multiple-testing correction.
            - **Exclusive within sampled analysis:** no reviewed outsider or unresolved sampled
              label is represented. This is not a universal absence claim.
            - **Near-exclusive within sampled analysis:** allows only the configured number of
              outsiders or unresolved sampled labels, all of which remain visible.
            """
        )
    species = service.list_species()
    template_rows = taxonomy_template_rows(species=species)
    render_table_downloads(
        records=template_rows,
        fieldnames=TAXONOMY_COLUMNS,
        file_stem=f"{service.resource.run_id}_taxonomy_mapping_template",
        key="taxonomy_template_download",
        tsv_label="Download taxonomy review template as TSV",
        excel_label="Download taxonomy review template as formatted Excel",
        column_definitions=_TAXONOMY_TEMPLATE_HELP,
        workbook_title=f"Taxonomy review template: {service.resource.run_id}",
    )
    st.caption(
        "The formatted workbook is convenient for review; save the completed mapping as UTF-8 "
        "TSV before loading it into the app or command-line mapper."
    )
    uploaded = st.file_uploader(
        "Use a reviewed taxonomy TSV for this browser session",
        type=("tsv", "txt"),
        help="An upload takes precedence over the optional launcher --taxonomy-map path.",
    )
    try:
        authority = _mapping_authority(
            expected_species=species,
            taxonomy_path_text=taxonomy_path_text,
            uploaded_data=uploaded.getvalue() if uploaded is not None else None,
        )
    except InputValidationError as error:
        _LOGGER.exception("Taxonomy mapping validation failed")
        st.error(f"Taxonomy mapping could not be used: {error}")
        return
    if authority is None:
        st.info(
            "Load a reviewed mapping above or launch with --taxonomy-map. The template "
            "contains every exact species label and deliberately starts as UNMAPPED. The "
            "orthofinder-taxonomy-map command can add exact NCBI taxdump candidates, but "
            "they remain PENDING_REVIEW until a person approves them."
        )
        return
    _render_mapping_audit(authority=authority)
    options = authority.taxon_options()
    if not options:
        st.warning("The mapping contains no REVIEWED records, so descendant searches are disabled.")
        return
    labels_to_taxa = {option.display_label(): option for option in options}
    group_types = service.list_group_types()
    with st.form("taxonomy_search_filters"):
        first = st.columns(4)
        taxon_label = first[0].selectbox(
            "Target lineage",
            tuple(labels_to_taxa),
            help="Reviewed NCBI taxon whose sampled descendants define the target set.",
        )
        group_type = first[1].selectbox(
            "Group system",
            group_types,
            help="HOG is hierarchical; LEGACY_ORTHOGROUP is the flat Orthogroups.tsv set.",
        )
        nodes = service.list_hierarchy_nodes(group_type=group_type)
        hierarchy_node = first[2].selectbox(
            "Species-tree level",
            nodes,
            format_func=lambda value: value or "ROOT",
            help="Exact species-tree node at which HOG membership is defined.",
        )
        mode_label = first[3].selectbox(
            "Search question",
            tuple(_MODE_LABELS),
            help="Open the explanation above before interpreting exclusivity or enrichment.",
        )
        mode = _MODE_LABELS[mode_label]

        second = st.columns(3)
        minimum_target = int(
            second[0].number_input(
                "Minimum target descendants represented",
                min_value=1,
                value=1,
                help="Minimum number of reviewed sampled target species present in a group.",
            )
        )
        minimum_coverage = float(
            second[1].number_input(
                "Minimum target descendant coverage",
                min_value=0.0,
                max_value=1.0,
                value=0.0,
                step=0.05,
                help="Represented target descendants / reviewed target descendants in this run.",
            )
        )
        minimum_purity = float(
            second[2].number_input(
                "Minimum mapped target purity",
                min_value=0.0,
                max_value=1.0,
                value=0.0,
                step=0.05,
                help="Target / all represented REVIEWED species; unresolved labels are separate.",
            )
        )

        third = st.columns(4)
        maximum_outside = int(
            third[0].number_input(
                "Maximum reviewed outsiders",
                min_value=0,
                value=1,
                disabled=mode != "NEAR_EXCLUSIVE",
            )
        )
        maximum_unresolved = int(
            third[1].number_input(
                "Maximum unresolved labels",
                min_value=0,
                value=0,
                disabled=mode != "NEAR_EXCLUSIVE",
            )
        )
        maximum_q = float(
            third[2].number_input(
                "Maximum adjusted q-value (BH)",
                min_value=0.0,
                max_value=1.0,
                value=0.05,
                step=0.01,
                disabled=mode != "ENRICHED",
            )
        )
        minimum_odds = float(
            third[3].number_input(
                "Minimum enrichment odds ratio",
                min_value=0.0,
                value=1.0,
                step=0.5,
                disabled=mode != "ENRICHED",
            )
        )
        fourth = st.columns(2)
        page_size = int(fourth[0].selectbox("Rows per page", (25, 50, 100, 250, 500), index=2))
        page_number = int(fourth[1].number_input("Page", min_value=1, value=1))
        st.form_submit_button("Search taxonomy", type="primary")
    option = labels_to_taxa[taxon_label]
    filters = TaxonomySearchFilters(
        group_type=group_type,
        hierarchy_node=hierarchy_node,
        target_taxon_id=option.taxon_id,
        mode=mode,
        minimum_target_species_count=minimum_target,
        minimum_target_coverage=minimum_coverage,
        minimum_mapped_purity=minimum_purity,
        maximum_outside_species_count=maximum_outside,
        maximum_unresolved_species_count=maximum_unresolved,
        maximum_q_value=maximum_q,
        minimum_odds_ratio=minimum_odds,
        page_size=page_size,
        page_number=page_number,
    )
    target_species = authority.target_species(taxon_id=option.taxon_id)
    with st.expander(
        f"Reviewed sampled descendants used for {option.name} ({len(target_species):,})",
        expanded=False,
    ):
        st.write(target_species)
    try:
        result = service.search_taxonomy_groups(authority=authority, filters=filters)
    except InputValidationError as error:
        _LOGGER.exception("Taxonomic group search failed")
        st.error(f"Taxonomic search could not be completed: {error}")
        return
    _render_search_scope(mode=mode, option=option, result=result)
    if not result.rows:
        st.warning("No groups match the selected taxonomic criteria.")
        return
    displayed = tuple(_display_taxonomy_row(row=row) for row in result.rows)
    st.dataframe(
        displayed,
        width="stretch",
        hide_index=True,
        column_config=_column_config(descriptions=_TAXONOMY_COLUMN_HELP),
    )
    render_table_downloads(
        records=displayed,
        file_stem=f"orthofinder_taxonomy_{mode.lower()}",
        key=f"taxonomy_result_{mode.lower()}",
        tsv_label="Download this taxonomic result page as TSV",
        excel_label="Download this taxonomic result page as formatted Excel",
        column_definitions=_TAXONOMY_COLUMN_HELP,
        workbook_title=f"OrthoFinder taxonomic search: {mode.lower()}",
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
    selected_label = st.selectbox("Inspect one taxonomic match", tuple(labels_to_keys))
    key = labels_to_keys[selected_label]
    actions = st.columns(2)
    if actions[0].button("Explore selected taxonomic match", type="primary"):
        _store_active_group(key=key)
        st.session_state["app_page"] = "Cluster explorer"
        st.rerun()
    basket = _comparison_keys()
    if actions[1].button(
        "Add taxonomic match to comparison",
        disabled=key in basket or len(basket) >= 12,
    ):
        _store_comparison_keys(keys=(*basket, key))
        st.success(f"Added to comparison workspace ({len(basket) + 1:,}/12).")
    st.subheader(key.display_label())
    species_result = tuple(
        _display_species_result_row(row=row) for row in service.get_group_species(key=key)
    )
    st.dataframe(
        species_result,
        width="stretch",
        hide_index=True,
        column_config=_column_config(descriptions=_SPECIES_RESULT_HELP),
    )
    render_table_downloads(
        records=species_result,
        file_stem=f"{key.group_id}_taxonomic_match_species",
        key=f"taxonomy_species_{key.display_label()}",
        tsv_label="Download selected group species as TSV",
        excel_label="Download selected group species as formatted Excel",
        column_definitions=_SPECIES_RESULT_HELP,
        workbook_title=f"Taxonomic match species: {key.group_id}",
    )


def _mapping_authority(
    *,
    expected_species: tuple[str, ...],
    taxonomy_path_text: str,
    uploaded_data: bytes | None,
) -> TaxonomyAuthority | None:
    """Load an uploaded mapping or optional launcher path with clear precedence."""

    if uploaded_data is not None:
        return parse_taxonomy_mapping(data=uploaded_data, expected_species=expected_species)
    if taxonomy_path_text.strip():
        return read_taxonomy_mapping(
            path=Path(taxonomy_path_text), expected_species=expected_species
        )
    return None


def _render_mapping_audit(*, authority: TaxonomyAuthority) -> None:
    """Render mapping status without hiding ambiguous or missing labels."""

    counts = authority.summary()
    columns = st.columns(5)
    for column, status in zip(
        columns,
        ("REVIEWED", "PENDING_REVIEW", "UNMAPPED", "AMBIGUOUS", "MISSING"),
        strict=True,
    ):
        column.metric(
            status.replace("_", " ").title(),
            f"{counts[status]:,}",
            help=_mapping_status_help(status=status),
        )
    audit = taxonomy_audit_rows(authority=authority)
    unresolved = tuple(row for row in audit if row["mapping_status"] != "REVIEWED")
    if unresolved:
        st.warning(
            f"{len(unresolved):,} sampled labels are unresolved or pending review. They "
            "never count as reviewed outsiders or descendants and remain visible in every "
            "matching group."
        )
        with st.expander(
            f"Unresolved and pending mapping rows ({len(unresolved):,})",
            expanded=True,
        ):
            st.dataframe(unresolved, width="stretch", hide_index=True)
    with st.expander(f"Complete taxonomy mapping audit ({len(audit):,})", expanded=False):
        st.dataframe(audit, width="stretch", hide_index=True)
        render_table_downloads(
            records=audit,
            file_stem="orthofinder_taxonomy_mapping_audit",
            key="taxonomy_mapping_audit",
            tsv_label="Download mapping audit as TSV",
            excel_label="Download mapping audit as formatted Excel",
            column_definitions=_TAXONOMY_TEMPLATE_HELP,
            workbook_title="OrthoFinder taxonomy mapping audit",
        )


def _render_search_scope(*, mode: str, option: Any, result: Any) -> None:
    """State sampled-universe limits and matching-page range."""

    start = (result.page_number - 1) * result.page_size + 1 if result.rows else 0
    finish = (result.page_number - 1) * result.page_size + len(result.rows)
    st.subheader("Matching groups")
    st.caption(
        f"Showing {start:,}–{finish:,} of {result.total_rows:,} matches for "
        f"{option.name} (NCBI taxon {option.taxon_id}). Sampled mapping universe: "
        f"{result.target_species_count:,} reviewed target, "
        f"{result.outside_species_count:,} reviewed outside and "
        f"{result.unresolved_species_count:,} unresolved species labels."
    )
    if mode == "ENRICHED":
        st.info(
            f"One-sided Fisher exact tests used species presence across "
            f"{result.tested_group_count:,} groups at the selected group authority and "
            "hierarchy. Q-values use Benjamini–Hochberg correction across that complete "
            "tested set. Unresolved species are excluded from the contingency universe."
        )
    elif mode == "SAMPLED_EXCLUSIVE":
        st.info(
            "Exclusive means no reviewed outside or unresolved species is represented in "
            "the group among the species sampled in this OrthoFinder analysis. It is not a "
            "universal biological absence claim."
        )
    elif mode == "NEAR_EXCLUSIVE":
        st.info(
            "Near-exclusive permits only the configured sampled outsiders/unresolved labels; "
            "every retained outsider is listed."
        )


def _display_taxonomy_row(*, row: dict[str, Any]) -> dict[str, Any]:
    """Return plain-language headings while retaining outsider identities."""

    displayed = {
        "Group ID": row["group_id"],
        "Proteins in group": row["member_count"],
        "Species represented": row["species_count"],
        "Target descendants represented": row["target_species_count"],
        "Target descendant coverage": row["target_coverage"],
        "Mapped target purity": row["mapped_target_fraction"],
        "Reviewed outsiders": row["outside_species_count"],
        "Outsider species": row["outsider_species_labels"],
        "Unresolved species": row["unresolved_species_labels"],
        "Stored distance coverage": _distance_status_label(
            value=row.get("computation_status")
        ),
        "Stored mean pair distance": row.get("mean_distance"),
        "Stored median pair distance": row.get("median_distance"),
        "Stored distance spread (SD)": row.get("population_stddev_distance"),
    }
    if "enrichment_q_value" in row:
        displayed.update(
            {
                "Enrichment odds ratio": row["enrichment_odds_ratio"],
                "Fisher p-value": row["enrichment_p_value"],
                "Adjusted q-value (BH)": row["enrichment_q_value"],
            }
        )
    return displayed


def _display_species_result_row(*, row: dict[str, Any]) -> dict[str, Any]:
    """Return readable copy-count fields for one represented species."""

    return {
        "Species": row["species_label"],
        "Proteins from species": row["species_member_count"],
        "Share of group": row["member_fraction"],
    }


def _distance_status_label(*, value: object) -> str:
    """Return a readable stored-distance status without hiding unknown codes."""

    if value is None or str(value).strip() == "":
        return "Not calculated"
    labels = {
        "COMPLETE": "Complete distances stored",
        "EXACT": "Complete distances stored",
        "DETERMINISTIC_MEMBER_SAMPLE": "Sampled distances stored",
    }
    code = str(value)
    return labels.get(code, code.replace("_", " ").title())


def _mapping_status_help(*, status: str) -> str:
    """Return one conservative explanation for a taxonomy mapping status."""

    descriptions = {
        "REVIEWED": "A person approved the species identity, taxon and lineage.",
        "PENDING_REVIEW": "A candidate exists but cannot support descendant claims yet.",
        "UNMAPPED": "No accepted taxon has been assigned.",
        "AMBIGUOUS": "More than one plausible mapping remains unresolved.",
        "MISSING": "The reviewed sidecar has no row for an expected dataset label.",
    }
    return descriptions.get(status, "Unrecognised status retained from the mapping authority.")


def _column_config(*, descriptions: dict[str, str]) -> dict[str, Any]:
    """Return Streamlit column definitions with hoverable help descriptions."""

    return {
        label: st.column_config.Column(label=label, help=description)
        for label, description in descriptions.items()
    }

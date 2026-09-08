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
from .models import GroupKey, TaxonomySearchFilters
from .queries import OrthoFinderQueryService
from .taxonomy import (
    TaxonomyAuthority,
    parse_taxonomy_mapping,
    read_taxonomy_mapping,
    taxonomy_audit_rows,
    taxonomy_template,
)
from .tsv import records_to_tsv

_LOGGER = logging.getLogger("orthofinder_interrogation_app.taxonomy_page")
_MODE_LABELS = {
    "Contains target descendants": "CONTAINS",
    "Enriched in target descendants": "ENRICHED",
    "Exclusive within sampled analysis": "SAMPLED_EXCLUSIVE",
    "Near-exclusive within sampled analysis": "NEAR_EXCLUSIVE",
}


def render_taxonomy_search(*, service: OrthoFinderQueryService, taxonomy_path_text: str) -> None:
    """Render mapping audit, target selection and four taxonomic searches."""

    st.header("Taxonomic search")
    st.caption(
        "Descendant logic uses only an explicit reviewed TSV mapping. Workflow labels are "
        "never interpreted as taxonomy from their spelling. The workflow is generated from "
        "the current resource species set, so it is not tied to this study's 60 species."
    )
    species = service.list_species()
    st.download_button(
        "Download taxonomy review template",
        data=taxonomy_template(species=species),
        file_name=f"{service.resource.run_id}_taxonomy_mapping_template.tsv",
        mime="text/tab-separated-values",
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
        taxon_label = first[0].selectbox("Target taxon", tuple(labels_to_taxa))
        group_type = first[1].selectbox("Group type", group_types)
        nodes = service.list_hierarchy_nodes(group_type=group_type)
        hierarchy_node = first[2].selectbox(
            "Hierarchy node", nodes, format_func=lambda value: value or "ROOT"
        )
        mode_label = first[3].selectbox("Search mode", tuple(_MODE_LABELS))
        mode = _MODE_LABELS[mode_label]

        second = st.columns(3)
        minimum_target = int(
            second[0].number_input("Minimum represented target species", min_value=1, value=1)
        )
        minimum_coverage = float(
            second[1].number_input(
                "Minimum target coverage",
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
                "Maximum reviewed outsider species",
                min_value=0,
                value=1,
                disabled=mode != "NEAR_EXCLUSIVE",
            )
        )
        maximum_unresolved = int(
            third[1].number_input(
                "Maximum unresolved species",
                min_value=0,
                value=0,
                disabled=mode != "NEAR_EXCLUSIVE",
            )
        )
        maximum_q = float(
            third[2].number_input(
                "Maximum BH q-value",
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
    st.dataframe(displayed, width="stretch", hide_index=True)
    st.download_button(
        "Download this taxonomic result page as TSV",
        data=records_to_tsv(records=result.rows),
        file_name=f"orthofinder_taxonomy_{mode.lower()}.tsv",
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
    selected_label = st.selectbox("Inspect one taxonomic match", tuple(labels_to_keys))
    key = labels_to_keys[selected_label]
    actions = st.columns(2)
    if actions[0].button("Open taxonomic match in Cluster explorer", type="primary"):
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
    st.dataframe(service.get_group_species(key=key), width="stretch", hide_index=True)


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
        column.metric(status.replace("_", " ").title(), f"{counts[status]:,}")
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
        st.download_button(
            "Download mapping audit as TSV",
            data=records_to_tsv(records=audit),
            file_name="orthofinder_taxonomy_mapping_audit.tsv",
            mime="text/tab-separated-values",
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
    """Return concise headings while retaining outsider identities."""

    displayed = {
        "Group": row["group_id"],
        "Members": row["member_count"],
        "Species": row["species_count"],
        "Target species": row["target_species_count"],
        "Target coverage": row["target_coverage"],
        "Mapped target purity": row["mapped_target_fraction"],
        "Reviewed outsiders": row["outside_species_count"],
        "Outsider species": row["outsider_species_labels"],
        "Unresolved species": row["unresolved_species_labels"],
        "Persisted distance status": row.get("computation_status") or "Not calculated",
        "Persisted mean distance": row.get("mean_distance"),
        "Persisted median distance": row.get("median_distance"),
        "Persisted distance SD": row.get("population_stddev_distance"),
    }
    if "enrichment_q_value" in row:
        displayed.update(
            {
                "Enrichment odds ratio": row["enrichment_odds_ratio"],
                "Fisher p-value": row["enrichment_p_value"],
                "BH q-value": row["enrichment_q_value"],
            }
        )
    return displayed

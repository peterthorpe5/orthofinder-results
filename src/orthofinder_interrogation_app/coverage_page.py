"""Streamlit page for reproducible taxonomy selection and dataset coverage trees."""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Any, Mapping

import streamlit as st

from orthofinder_results import __version__
from orthofinder_results.errors import InputValidationError, OrthoFinderResultsError

from .coverage_exports import (
    SelectionCoverageRun,
    build_selection_coverage_run,
    coverage_export_files,
    coverage_export_zip,
)
from .coverage_tree_render import coverage_tree_figure
from .exports import render_table_downloads
from .focus import FocusProteinAuthority
from .focus_page import load_focus_authority
from .models import FocusClusterFilters, GroupKey
from .queries import OrthoFinderQueryService
from .taxonomy import (
    DEFAULT_TAXONOMY_FILENAME,
    TAXONOMY_COLUMNS,
    TAXONOMY_TEMPLATE_HELP,
    TaxonomyAuthority,
    parse_taxonomy_mapping,
    read_matching_bundled_taxonomy,
    read_taxonomy_mapping,
    taxonomy_template_rows,
)
from .taxonomy_selection import (
    ExpectedTaxaAuthority,
    SelectionSpec,
    TaxonomyGraph,
    build_taxonomy_graph,
    default_expected_taxa,
    make_selection,
    parse_expected_taxa,
    read_expected_taxa,
    selection_summary,
)

_LOGGER = logging.getLogger("orthofinder_interrogation_app.coverage_page")
_LAST_RUN_STATE = "orthofinder_last_valid_coverage_run"
_REBUILD_STATE = "orthofinder_coverage_rebuild"
_ACTIVE_NODE_STATE = "orthofinder_coverage_active_node"
_SELECTOR_KEYS = {
    "REQUIRED_EXACT": "coverage_required_exact",
    "INCLUDE_CLADE": "coverage_include_clade",
    "ONLY_IN_CLADE": "coverage_only_in_clade",
    "EXCLUDE_EXACT": "coverage_exclude_exact",
    "EXCLUDE_CLADE": "coverage_exclude_clade",
}
_COVERAGE_HELP = {
    "Taxon ID": "Stable identifier under the recorded taxonomy authority.",
    "Taxon": "Reviewed accepted taxon name.",
    "Rank": "Reviewed rank from the pinned mapping lineage.",
    "Selection state": "Predicate role, kept independent from dataset coverage.",
    "Coverage state": "Represented, expected-no-data, ancestor, outsider or not expected.",
    "Represented labels": "Exact OrthoFinder input labels reviewed to this terminal taxon.",
    "Represented descendants": "Displayed bounded terminal descendants represented in input.",
    "Expected-no-data descendants": (
        "Expected displayed descendants not represented in this dataset; not biological absence."
    ),
}
_EVALUATION_HELP = {
    "Group ID": "Exact OrthoFinder group identifier in the selected collection.",
    "Passes all predicates": "True only when every active selector passes under AND semantics.",
    "Mapped proteins": "Group members whose exact species label has a reviewed mapping.",
    "Unmapped proteins": "Members whose lineage cannot be asserted; only-in fails closed.",
    "Represented taxon IDs": "Reviewed group terminal taxa.",
    "Outside-scope taxon IDs": "Mapped hits outside include scope; outsiders remain allowed.",
    "Unmapped species labels": "Exact unresolved input labels, never placed on the tree.",
    "Failure reasons": "Deterministic reasons for every failed predicate.",
    "Predicate audit": "Pass/fail state for each normalised selector.",
}
_EXPECTED_HELP = {
    "taxon_id": "Reviewed terminal identifier present in the taxonomy mapping.",
    "reason": "Why this taxon belongs in the declared dataset-coverage universe.",
    "source": "Reviewed manifest, sampling plan or other bounded authority.",
    "included": "True includes the taxon; false retains an audited inactive row.",
}


def render_selection_coverage_tree(
    *,
    service: OrthoFinderQueryService,
    taxonomy_path_text: str,
    expected_path_text: str,
    focus_path_text: str,
) -> None:
    """Render validated selectors, a taxonomy coverage tree and group audit.

    Args:
        service: Read-only query service for one completed resource.
        taxonomy_path_text: Optional reviewed mapping supplied at launch.
        expected_path_text: Optional expected-taxon universe supplied at launch.
        focus_path_text: Optional custom protein focus supplied at launch.
    """

    st.header("Selection coverage tree")
    st.write(
        "Construct exact, clade, only-in and exclusion predicates against a pinned reviewed "
        "taxonomy. This is a **taxonomy/species coverage tree**—not the OrthoFinder species "
        "tree and not a HOG or gene-tree phylogram."
    )
    st.info(
        "Expected no data always means ‘not represented in this dataset’. It is never a "
        "claim that a gene, HOG or biological function is absent from that organism."
    )
    species = service.list_species()
    st.subheader("Prepare the species-to-taxonomy mapping")
    st.write(
        "This upload is a completed **taxonomy review table with one row per exact "
        "OrthoFinder species label**. It is not an E3 protein list, a gene table, a species "
        "tree or an NCBI taxdump. Download the dataset-specific template below, populate "
        "reviewed names, taxon IDs and lineages, then upload the completed UTF-8 TSV."
    )
    with st.expander("What must this file contain?", expanded=True):
        st.markdown(
            """
            - Keep every `workflow_species_label` exactly as supplied; it links the mapping to
              this resource.
            - A usable row has `mapping_status=REVIEWED`, an accepted species name, a stable
              terminal taxon ID, and ordered lineage IDs/names. Pending, ambiguous and
              unmapped rows remain visible but cannot support clade predicates.
            - Record the mapping source/version and reviewer fields so the decision is
              auditable. Candidate mappings produced by `orthofinder-taxonomy-map` still need
              human approval before their status becomes `REVIEWED`.
            """
        )
        template_rows = taxonomy_template_rows(species=species)
        render_table_downloads(
            records=template_rows,
            fieldnames=TAXONOMY_COLUMNS,
            file_stem=f"{service.resource.run_id}_taxonomy_mapping_template",
            key="coverage_taxonomy_template_download",
            tsv_label="Download this dataset's taxonomy review template",
            excel_label="Download the formatted review workbook",
            column_definitions=TAXONOMY_TEMPLATE_HELP,
            workbook_title=f"Taxonomy review template: {service.resource.run_id}",
        )
    taxonomy_upload = st.file_uploader(
        "Upload the completed reviewed taxonomy mapping TSV",
        type=("tsv", "txt"),
        help=(
            "Use the dataset-specific template provided immediately above. Upload takes "
            "precedence over --taxonomy-map. Only REVIEWED rows enter the tree; pending, "
            "ambiguous and unmapped labels remain in a separate inventory."
        ),
        key="coverage_taxonomy_upload",
    )
    try:
        authority = _load_taxonomy_authority(
            expected_species=species,
            taxonomy_path_text=taxonomy_path_text,
            uploaded_data=(
                taxonomy_upload.getvalue() if taxonomy_upload is not None else None
            ),
        )
        if authority is None:
            st.warning(
                "Complete and upload the species-to-taxonomy review template above to enable "
                "the selection coverage tree. Rows still marked UNMAPPED or PENDING_REVIEW "
                "cannot be used for clade selection."
            )
            return
        if taxonomy_upload is None and not taxonomy_path_text.strip():
            st.success(
                "Using the packaged Results_Feb26 mapping because all 60 exact species "
                "labels match this resource. The file is "
                f"`data/{DEFAULT_TAXONOMY_FILENAME}`: 59 mappings are REVIEWED and "
                "`Leismania_major` remains visibly UNMAPPED pending confirmation of the "
                "apparent spelling error."
            )
        graph = build_taxonomy_graph(authority=authority)
    except InputValidationError as error:
        _LOGGER.exception("Coverage taxonomy authority failed validation")
        st.error(f"Reviewed taxonomy mapping could not be used: {error}")
        return
    expected_upload = st.file_uploader(
        "Optional expected-taxon universe TSV",
        type=("tsv", "txt"),
        help=(
            "Columns: taxon_id, reason, source, included. Leave empty to use the exact "
            "reviewed species supplied to this OrthoFinder run."
        ),
        key="coverage_expected_upload",
    )
    try:
        expected = _load_expected_authority(
            graph=graph,
            expected_path_text=expected_path_text,
            uploaded_data=expected_upload.getvalue() if expected_upload is not None else None,
        )
    except InputValidationError as error:
        _LOGGER.exception("Expected-taxon authority failed validation")
        st.error(f"Expected-taxon universe could not be used: {error}")
        return
    with st.expander("Expected-taxon universe and editable template", expanded=False):
        st.write(
            "Download the current bounded universe, edit it as UTF-8 TSV, then upload it "
            "above. A taxon expected without input data must first have a reviewed terminal "
            "row in the taxonomy mapping, normally with `role=expected`."
        )
        expected_records = tuple(
            {
                "taxon_id": row.taxon_id,
                "reason": row.reason,
                "source": row.source,
                "included": row.included,
            }
            for row in expected.records
        )
        render_table_downloads(
            records=expected_records,
            file_stem="expected_taxa",
            key="coverage_expected_taxa_template",
            tsv_label="Download expected taxa as editable TSV",
            excel_label="Download expected taxa as formatted Excel",
            fieldnames=("taxon_id", "reason", "source", "included"),
            column_definitions=_EXPECTED_HELP,
            workbook_title="OrthoFinder expected-taxon universe",
        )
    mapping_counts = authority.summary()
    unresolved_count = sum(
        mapping_counts[key]
        for key in ("PENDING_REVIEW", "UNMAPPED", "AMBIGUOUS", "MISSING")
    )
    top = st.columns(6)
    top[0].metric("Taxonomy authority", graph.authority or "Unspecified")
    top[1].metric("Pinned release", graph.release or "Unspecified")
    top[2].metric("Reviewed input taxa", f"{len(graph.represented_taxon_ids):,}")
    top[3].metric(
        "Unmapped input labels",
        f"{unresolved_count:,}",
    )
    top[4].metric("Expected taxa", f"{len(expected.included_taxon_ids):,}")
    top[5].metric(
        "Expected without input data",
        f"{len(set(expected.included_taxon_ids).difference(graph.represented_taxon_ids)):,}",
    )
    st.caption(
        f"Mapping SHA-256: {graph.mapping_sha256}; expected-universe SHA-256: "
        f"{expected.sha256}. Rendering and evaluation are offline."
    )
    node_options = _node_options(graph=graph)
    terminal_options = {
        label: identifier
        for label, identifier in node_options.items()
        if identifier in graph.terminal_ids
    }
    with st.form("selection_coverage_controls"):
        st.subheader("Selection predicates")
        st.caption(
            "Every active selector composes with logical AND. Several only-in clades use "
            "the intersection of their mapped terminal sets."
        )
        row_one = st.columns(3)
        row_one[0].multiselect(
            "Required exact taxa",
            tuple(terminal_options),
            key=_SELECTOR_KEYS["REQUIRED_EXACT"],
            help="Every selected exact terminal taxon must occur in a passing group.",
        )
        row_one[1].multiselect(
            "Include clades",
            tuple(node_options),
            key=_SELECTOR_KEYS["INCLUDE_CLADE"],
            help="Every selected clade must contribute at least one member; outsiders are allowed.",
        )
        row_one[2].multiselect(
            "Only in clades",
            tuple(node_options),
            key=_SELECTOR_KEYS["ONLY_IN_CLADE"],
            help="No mapped member may lie outside the intersection; unmapped members fail closed.",
        )
        row_two = st.columns(3)
        row_two[0].multiselect(
            "Exclude exact taxa",
            tuple(terminal_options),
            key=_SELECTOR_KEYS["EXCLUDE_EXACT"],
            help="No selected exact terminal may occur in a passing group.",
        )
        row_two[1].multiselect(
            "Exclude clades",
            tuple(node_options),
            key=_SELECTOR_KEYS["EXCLUDE_CLADE"],
            help="No descendant of any selected clade may occur.",
        )
        compact = row_two[2].checkbox(
            "Compact unary neutral lineages",
            value=False,
            help="Selected and excluded nodes are never collapsed.",
        )
        st.subheader("Group evaluation scope")
        scope = st.columns(4)
        focus_limited = scope[0].checkbox(
            "Only clusters containing focus proteins",
            value=True,
            help=(
                "Enabled by default for the current E3 goals. Clear it to evaluate every "
                "group in the selected collection."
            ),
        )
        group_types = service.list_group_types()
        default_type = "HOG" if "HOG" in group_types else group_types[0]
        group_type = scope[1].selectbox(
            "Group system",
            group_types,
            index=group_types.index(default_type),
            key="coverage_group_type",
        )
        hierarchy_nodes = service.list_hierarchy_nodes(group_type=group_type)
        default_node = "N0" if "N0" in hierarchy_nodes else hierarchy_nodes[0]
        hierarchy_node = scope[2].selectbox(
            "Species-tree level",
            hierarchy_nodes,
            index=hierarchy_nodes.index(default_node),
            format_func=lambda value: value or "ROOT",
            key="coverage_hierarchy_node",
        )
        maximum_groups = int(
            scope[3].selectbox(
                "Maximum groups",
                (2_000, 10_000, 50_000, 100_000, 250_000),
                index=1,
                help="Fail closed when the chosen scope exceeds this declared audit bound.",
            )
        )
        submitted = st.form_submit_button(
            "Build selection coverage tree",
            type="primary",
        )
    rebuild = submitted or bool(st.session_state.pop(_REBUILD_STATE, False))
    if rebuild:
        try:
            selection = _selection_from_widgets(
                graph=graph,
                node_options=node_options,
            )
            focus_authority: FocusProteinAuthority | None = None
            group_ids: tuple[str, ...] = ()
            if focus_limited:
                focus_authority = load_focus_authority(
                    focus_path_text=focus_path_text,
                    uploaded_data=None,
                )
                focus_result = service.search_focus_clusters(
                    filters=FocusClusterFilters(
                        protein_identifiers=focus_authority.identifiers,
                        group_type=group_type,
                        hierarchy_node=hierarchy_node,
                        maximum_rows=maximum_groups,
                    )
                )
                if focus_result.truncated:
                    raise InputValidationError(
                        f"Focus authority matches {focus_result.total_rows:,} clusters, but the "
                        "current focus-query safety bound retained only "
                        f"{len(focus_result.rows):,}. Narrow the focus list before building a "
                        "complete selection audit."
                    )
                group_ids = tuple(str(row["group_id"]) for row in focus_result.rows)
            group_rows = (
                service.get_group_species_collection(
                    group_type=group_type,
                    hierarchy_node=hierarchy_node,
                    group_ids=group_ids,
                )
                if group_ids or not focus_limited
                else ()
            )
            run = build_selection_coverage_run(
                graph=graph,
                expected=expected,
                selection=selection,
                group_species_rows=group_rows,
                run_id=service.resource.run_id,
                resource_identity=_resource_identity(service=service),
                package_version=__version__,
                compact=compact,
                maximum_groups=maximum_groups,
                focus_authority_name=(
                    focus_authority.source_name if focus_authority is not None else ""
                ),
                focus_authority_sha256=(
                    focus_authority.sha256 if focus_authority is not None else ""
                ),
                focus_limited=focus_limited,
            )
            st.session_state[_LAST_RUN_STATE] = run
        except (InputValidationError, OrthoFinderResultsError) as error:
            _LOGGER.exception("Selection coverage build failed")
            st.error(f"Selection is invalid and was not applied: {error}")
            if st.session_state.get(_LAST_RUN_STATE) is not None:
                st.warning("The last valid, synchronised tree and group result remain below.")
    run = st.session_state.get(_LAST_RUN_STATE)
    if run is None:
        st.info(
            "Choose predicates and build the tree. The E3 focus authority is the default group "
            "scope, but every selector and taxonomy mapping remains generic and replaceable."
        )
        return
    _render_coverage_run(run=run, service=service)


def _render_coverage_run(
    *, run: SelectionCoverageRun, service: OrthoFinderQueryService
) -> None:
    """Render one internally synchronised last-valid selection result."""

    st.subheader("Current valid selection")
    for statement in selection_summary(selection=run.selection):
        st.write(f"- {statement}")
    summary = run.summary()
    metrics = st.columns(6)
    labels = (
        ("Displayed taxa", "displayed_nodes"),
        ("Represented taxa", "represented_terminal_taxa"),
        ("Expected, no input data", "expected_no_data_taxa"),
        ("Outside-scope hits", "outside_scope_taxa_with_hits"),
        ("Groups evaluated", "evaluated_groups"),
        ("Groups passing", "passing_groups"),
    )
    for column, (label, key) in zip(metrics, labels, strict=True):
        column.metric(label, f"{summary[key]:,}")
    event = st.plotly_chart(
        coverage_tree_figure(tree=run.tree),
        width="stretch",
        config={"displaylogo": False, "scrollZoom": True},
        key="selection_coverage_plot",
        on_select="rerun",
        selection_mode="points",
    )
    selected_points = event.selection.points if event is not None else ()
    if selected_points and selected_points[0].get("customdata"):
        custom = selected_points[0]["customdata"]
        selected_id = str(custom[0] if isinstance(custom, (list, tuple)) else custom)
        if selected_id in run.tree.node_by_id:
            st.session_state[_ACTIVE_NODE_STATE] = selected_id
    _render_node_actions(run=run)
    with st.expander("Legend and interpretation", expanded=True):
        st.markdown(
            """
            - **Filled circle:** represented by at least one exact reviewed input label.
            - **Hollow circle:** reviewed expected taxon not represented in this dataset.
            - **Amber diamond:** a represented hit outside the included-clade scope.
            - **Blue outline/branch:** an active required, include or only-in selector.
            - **Red outline/dashed branch:** an exact or clade exclusion.
            - **Grey square:** an ancestor needed to connect displayed terminals.

            Colours are paired with shapes, outlines and dash styles. Unmapped input labels are
            listed separately and never placed speculatively on the tree.
            """
        )
    coverage_rows = tuple(_coverage_display(row=row.as_record()) for row in run.tree.nodes)
    st.subheader("Taxon coverage audit")
    st.dataframe(
        coverage_rows,
        width="stretch",
        hide_index=True,
        column_config={
            label: st.column_config.Column(label=label, help=description)
            for label, description in _COVERAGE_HELP.items()
        },
    )
    render_table_downloads(
        records=coverage_rows,
        file_stem="taxon_coverage",
        key="coverage_table_download",
        tsv_label="Download coverage table as TSV",
        excel_label="Download coverage table as formatted Excel",
        column_definitions=_COVERAGE_HELP,
        workbook_title="OrthoFinder selection taxon coverage",
    )
    st.subheader("Filtered groups and predicate audit")
    evaluations = tuple(_evaluation_display(row=row.as_record()) for row in run.evaluations)
    if evaluations:
        st.dataframe(
            evaluations,
            width="stretch",
            hide_index=True,
            column_config={
                label: st.column_config.Column(label=label, help=description)
                for label, description in _EVALUATION_HELP.items()
            },
        )
        render_table_downloads(
            records=evaluations,
            file_stem="group_taxon_evaluation",
            key="coverage_group_evaluation_download",
            tsv_label="Download group audit as TSV",
            excel_label="Download group audit as formatted Excel",
            column_definitions=_EVALUATION_HELP,
            workbook_title="OrthoFinder group taxonomy evaluation",
        )
        passing = tuple(row for row in run.evaluations if row.predicate_pass)
        if passing:
            labels_to_groups = {
                f"{row.group_type} | {row.hierarchy_node or 'ROOT'} | {row.group_id}": row
                for row in passing
            }
            selected_label = st.selectbox(
                "Passing group to open",
                tuple(labels_to_groups),
            )
            selected = labels_to_groups[selected_label]
            st.button(
                "Open passing group in cluster explorer",
                on_click=_open_passing_group,
                kwargs={
                    "run_id": selected.run_id,
                    "group_type": selected.group_type,
                    "hierarchy_node": selected.hierarchy_node,
                    "group_id": selected.group_id,
                },
            )
    else:
        st.info("No groups occurred in the selected focus authority and collection.")
    if run.unmapped_rows:
        st.subheader("Unmapped input labels")
        st.warning(
            "These exact source labels have no reviewed placement. They are retained for "
            "audit and cause only-in predicates to fail closed when represented in a group."
        )
        st.dataframe(run.unmapped_rows, width="stretch", hide_index=True)
    files = coverage_export_files(run=run)
    with st.expander("Reproducible selection package", expanded=False):
        st.write(
            "The ZIP contains TSV, Newick plus styles, SVG, PDF, JSON provenance and "
            "SHA-256 checksums generated from this exact synchronised selection."
        )
        st.download_button(
            "Download complete selection package",
            data=coverage_export_zip(run=run),
            file_name=f"{service.resource.run_id}_selection_coverage.zip",
            mime="application/zip",
            key="coverage_complete_zip",
        )
        download_columns = st.columns(4)
        for index, name in enumerate(
            (
                "selection_coverage_tree.svg",
                "selection_coverage_tree.pdf",
                "selection_manifest.json",
                "selection_coverage_tree.newick",
            )
        ):
            download_columns[index].download_button(
                f"Download {name.rsplit('.', maxsplit=1)[-1].upper()}",
                data=files[name],
                file_name=name,
                key=f"coverage_{name}",
            )


def _render_node_actions(*, run: SelectionCoverageRun) -> None:
    """Offer keyboard-accessible valid actions for a clicked or selected tree node."""

    options = {
        f"{node.taxon_name} | {node.taxon_rank} | {node.taxon_id}": node.taxon_id
        for node in sorted(
            run.tree.nodes,
            key=lambda row: (row.taxon_name.casefold(), row.taxon_name, row.taxon_id),
        )
    }
    active_id = st.session_state.get(_ACTIVE_NODE_STATE, run.tree.root_taxon_id)
    labels = tuple(options)
    selected_index = next(
        (index for index, label in enumerate(labels) if options[label] == active_id),
        0,
    )
    selected_label = st.selectbox(
        "Tree node actions (also updated by clicking a node)",
        labels,
        index=selected_index,
        key="coverage_node_action_selector",
        help="Keyboard users can choose the same taxon here without using the chart.",
    )
    taxon_id = options[selected_label]
    st.session_state[_ACTIVE_NODE_STATE] = taxon_id
    actions = st.columns(6)
    actions[0].button(
        "Require exact",
        disabled=taxon_id not in run.graph.terminal_ids,
        on_click=_add_node_selector,
        kwargs={"selector_type": "REQUIRED_EXACT", "taxon_id": taxon_id, "run": run},
    )
    actions[1].button(
        "Include clade",
        on_click=_add_node_selector,
        kwargs={"selector_type": "INCLUDE_CLADE", "taxon_id": taxon_id, "run": run},
    )
    actions[2].button(
        "Only in clade",
        on_click=_add_node_selector,
        kwargs={"selector_type": "ONLY_IN_CLADE", "taxon_id": taxon_id, "run": run},
    )
    actions[3].button(
        "Exclude exact",
        disabled=taxon_id not in run.graph.terminal_ids,
        on_click=_add_node_selector,
        kwargs={"selector_type": "EXCLUDE_EXACT", "taxon_id": taxon_id, "run": run},
    )
    actions[4].button(
        "Exclude clade",
        on_click=_add_node_selector,
        kwargs={"selector_type": "EXCLUDE_CLADE", "taxon_id": taxon_id, "run": run},
    )
    actions[5].button(
        "Clear node selection",
        on_click=_clear_node_selectors,
        kwargs={"taxon_id": taxon_id, "run": run},
    )


def _add_node_selector(
    *, selector_type: str, taxon_id: str, run: SelectionCoverageRun
) -> None:
    """Add one tree-node predicate through the same multiselect widget state."""

    key = _SELECTOR_KEYS[selector_type]
    label = next(
        label
        for label, identifier in _node_options(graph=run.graph).items()
        if identifier == taxon_id
    )
    values = list(st.session_state.get(key, ()))
    if label not in values:
        values.append(label)
    st.session_state[key] = values
    st.session_state[_REBUILD_STATE] = True


def _clear_node_selectors(*, taxon_id: str, run: SelectionCoverageRun) -> None:
    """Remove one node from every selector widget and rebuild synchronously."""

    labels = _node_options(graph=run.graph)
    target_labels = {label for label, identifier in labels.items() if identifier == taxon_id}
    for key in _SELECTOR_KEYS.values():
        st.session_state[key] = [
            label
            for label in st.session_state.get(key, ())
            if label not in target_labels
        ]
    st.session_state[_REBUILD_STATE] = True


def _open_passing_group(
    *, run_id: str, group_type: str, hierarchy_node: str, group_id: str
) -> None:
    """Store one validated passing group and move to the cluster explorer."""

    from .evolutionary_page import _store_active_group

    _store_active_group(
        key=GroupKey(
            run_id=run_id,
            group_type=group_type,
            hierarchy_node=hierarchy_node,
            group_id=group_id,
        )
    )
    st.session_state["app_page"] = "Cluster explorer"


def _selection_from_widgets(
    *,
    graph: TaxonomyGraph,
    node_options: Mapping[str, str],
) -> SelectionSpec:
    """Translate visible selector labels into one validated domain selection."""

    def values(selector_type: str) -> tuple[str, ...]:
        """Resolve one widget's labels without accepting stale or invented values."""

        labels = tuple(st.session_state.get(_SELECTOR_KEYS[selector_type], ()))
        unknown = tuple(label for label in labels if label not in node_options)
        if unknown:
            raise InputValidationError(
                "A selector contains a stale taxonomy label: " + "; ".join(unknown)
            )
        return tuple(node_options[label] for label in labels)

    return make_selection(
        graph=graph,
        required_exact=values("REQUIRED_EXACT"),
        include_clade=values("INCLUDE_CLADE"),
        only_in_clade=values("ONLY_IN_CLADE"),
        exclude_exact=values("EXCLUDE_EXACT"),
        exclude_clade=values("EXCLUDE_CLADE"),
    )


def _load_taxonomy_authority(
    *,
    expected_species: tuple[str, ...],
    taxonomy_path_text: str,
    uploaded_data: bytes | None,
) -> TaxonomyAuthority | None:
    """Load upload, launcher path or exact-label packaged study default."""

    if uploaded_data is not None:
        return parse_taxonomy_mapping(data=uploaded_data, expected_species=expected_species)
    if taxonomy_path_text.strip():
        return read_taxonomy_mapping(
            path=Path(taxonomy_path_text),
            expected_species=expected_species,
        )
    return read_matching_bundled_taxonomy(expected_species=expected_species)


def _load_expected_authority(
    *, graph: TaxonomyGraph, expected_path_text: str, uploaded_data: bytes | None
) -> ExpectedTaxaAuthority:
    """Load an explicit expected universe or use reviewed input taxa by default."""

    if uploaded_data is not None:
        return parse_expected_taxa(
            data=uploaded_data,
            source_name="uploaded_expected_taxa.tsv",
            graph=graph,
        )
    if expected_path_text.strip():
        return read_expected_taxa(path=Path(expected_path_text), graph=graph)
    return default_expected_taxa(graph=graph)


def _node_options(*, graph: TaxonomyGraph) -> dict[str, str]:
    """Return unique, deterministic selector labels for every reviewed graph node."""

    return {
        f"{node.name} | {node.rank} | ID {node.taxon_id}": node.taxon_id
        for node in sorted(
            graph.nodes,
            key=lambda row: (row.name.casefold(), row.name, row.taxon_id),
        )
    }


def _resource_identity(*, service: OrthoFinderQueryService) -> str:
    """Return a small immutable resource-manifest checksum or database identity."""

    root = service.resource.resource_path
    manifest = (
        root / "run_manifest.json"
        if root.is_dir()
        else root.parent.parent / "run_manifest.json"
    )
    if manifest.is_file():
        return "run_manifest_sha256:" + hashlib.sha256(manifest.read_bytes()).hexdigest()
    stat = service.resource.database_path.stat()
    return (
        f"duckdb:{service.resource.run_id}:size={stat.st_size}:"
        f"schema={service.resource.schema_version}"
    )


def _coverage_display(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return readable headings for one taxon coverage record."""

    return {
        "Taxon ID": row["taxon_id"],
        "Taxon": row["taxon_name"],
        "Rank": row["taxon_rank"],
        "Selection state": row["selection_state"],
        "Coverage state": row["coverage_state"],
        "Represented labels": row["represented_input_labels"],
        "Represented descendants": row["represented_terminal_count"],
        "Expected-no-data descendants": row["expected_no_data_terminal_count"],
        "Excluded descendants": row["excluded_terminal_count"],
    }


def _evaluation_display(*, row: Mapping[str, Any]) -> dict[str, Any]:
    """Return readable headings for one group predicate evaluation."""

    return {
        "Group system": row["group_type"],
        "Species-tree level": row["hierarchy_node"] or "ROOT",
        "Group ID": row["group_id"],
        "Passes all predicates": row["predicate_pass"],
        "Mapped proteins": row["mapped_member_count"],
        "Unmapped proteins": row["unmapped_member_count"],
        "Represented taxon IDs": row["represented_taxon_ids"],
        "Outside-scope taxon IDs": row["outside_selected_scope_taxon_ids"],
        "Unmapped species labels": row["unmapped_species_labels"],
        "Failure reasons": row["predicate_failure_reasons"],
        "Predicate audit": row["predicate_audit"],
        "Selection manifest ID": row["selection_manifest_id"],
    }

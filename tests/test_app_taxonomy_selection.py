"""Tests for generic reviewed-taxonomy selection and coverage semantics."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from orthofinder_interrogation_app import taxonomy_selection
from orthofinder_interrogation_app.coverage_exports import (
    EXPORT_FILENAMES,
    build_selection_coverage_run,
    coverage_export_files,
    coverage_export_zip,
    publish_selection_coverage_run,
)
from orthofinder_interrogation_app.coverage_tree import (
    _least_common_root,
    _reconcile_tree,
    build_coverage_tree,
    coverage_records,
    tree_edge_records,
    tree_to_newick,
)
from orthofinder_interrogation_app.coverage_tree_render import (
    coverage_tree_figure,
    coverage_tree_layout,
    style_records,
    tree_to_pdf,
    tree_to_svg,
)
from orthofinder_interrogation_app.taxonomy import (
    TaxonomyAuthority,
    TaxonomyRecord,
    parse_taxonomy_mapping,
)
from orthofinder_interrogation_app.taxonomy_selection import (
    ExpectedTaxaAuthority,
    ExpectedTaxon,
    build_taxonomy_graph,
    default_expected_taxa,
    evaluate_group_rows,
    evaluate_group_taxonomy,
    expected_taxa_to_tsv,
    input_taxon_mapping_records,
    make_selection,
    parse_expected_taxa,
    read_expected_taxa,
    selection_predicate_records,
    selection_summary,
    unmapped_input_records,
)
from orthofinder_interrogation_app.tsv import records_to_tsv
from orthofinder_results.errors import InputValidationError


def _record(
    *,
    label: str,
    name: str,
    taxon_id: int | None,
    lineage: tuple[tuple[int, str, str], ...],
    role: str = "input",
    authority: str = "NCBI Taxonomy",
    authority_taxon_id: str = "",
    status: str = "REVIEWED",
) -> TaxonomyRecord:
    """Build one concise reviewed mapping record for domain tests."""

    parent = lineage[-1] if lineage else None
    return TaxonomyRecord(
        workflow_species_label=label,
        source_species_name=name,
        accepted_species_name=name,
        ncbi_taxon_id=taxon_id,
        parent_taxon_id=parent[0] if parent is not None else None,
        parent_taxon_name=parent[1] if parent is not None else "",
        lineage_taxon_ids=tuple(item[0] for item in lineage),
        lineage_names=tuple(item[1] for item in lineage),
        mapping_status=status,
        mapping_method="manual exact review" if status == "REVIEWED" else "",
        mapping_source="pinned taxonomy fixture" if status == "REVIEWED" else "",
        source_date="2026-09-01" if status == "REVIEWED" else "",
        source_version="2026-08" if status == "REVIEWED" else "",
        reviewed_by="Tester" if status == "REVIEWED" else "",
        reviewed_at_utc="2026-09-01T00:00:00Z" if status == "REVIEWED" else "",
        review_note="Reviewed exact fixture." if status == "REVIEWED" else "Unresolved.",
        lineage_ranks=tuple(item[2] for item in lineage),
        taxon_rank="cultivar" if authority_taxon_id else "species",
        taxonomy_authority=authority if status == "REVIEWED" else "",
        taxonomy_release="2026-08" if status == "REVIEWED" else "",
        source_name_original=label,
        authority_taxon_id=authority_taxon_id,
        role=role,
    )


def _reviewed_graph():
    """Return a graph spanning plants, an infraspecific taxon and one unresolved input."""

    root = ((1, "root", "no rank"), (33090, "Viridiplantae", "kingdom"))
    flowering = (*root, (3398, "Magnoliopsida", "class"))
    grass = (*flowering, (4479, "Poaceae", "family"))
    solanaceae = (*flowering, (4081, "Solanaceae", "family"))
    records = (
        _record(label="Rice", name="Oryza sativa", taxon_id=4530, lineage=grass),
        _record(label="Maize", name="Zea mays", taxon_id=4577, lineage=grass),
        _record(
            label="Potato",
            name="Solanum tuberosum",
            taxon_id=4113,
            lineage=solanaceae,
        ),
        _record(
            label="Arabidopsis",
            name="Arabidopsis thaliana",
            taxon_id=3702,
            lineage=flowering,
        ),
        _record(
            label="Wheat_expected",
            name="Triticum aestivum",
            taxon_id=4565,
            lineage=grass,
            role="expected",
        ),
        _record(
            label="Potato_red_cultivar",
            name="Solanum tuberosum Red Example",
            taxon_id=None,
            lineage=(*solanaceae, (4113, "Solanum tuberosum", "species")),
            authority="Project Cultivar Registry",
            authority_taxon_id="CULT-RED",
        ),
        _record(
            label="Unknown_input",
            name="",
            taxon_id=None,
            lineage=(),
            status="UNMAPPED",
        ),
    )
    authority = TaxonomyAuthority(
        records=records,
        expected_species=(
            "Arabidopsis",
            "Maize",
            "Potato",
            "Potato_red_cultivar",
            "Rice",
            "Unknown_input",
        ),
        mapping_sha256="a" * 64,
    )
    return build_taxonomy_graph(authority=authority)


def _group_rows() -> tuple[dict[str, object], ...]:
    """Return four deliberately different group compositions."""

    compositions = {
        "G1": {"Rice": 1, "Maize": 1, "Arabidopsis": 1},
        "G2": {"Rice": 1, "Potato": 1},
        "G3": {"Arabidopsis": 2},
        "G4": {"Rice": 1, "Unknown_input": 1},
    }
    return tuple(
        {
            "run_id": "synthetic",
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": group_id,
            "species_label": species,
            "species_member_count": count,
        }
        for group_id, species_counts in compositions.items()
        for species, count in species_counts.items()
    )


def _expected(graph) -> ExpectedTaxaAuthority:
    """Return input taxa plus one reviewed expected-no-data grass."""

    records = tuple(
        ExpectedTaxon(
            taxon_id=taxon_id,
            reason="Reviewed sampling manifest",
            source="fixture manifest",
        )
        for taxon_id in sorted(set(graph.represented_taxon_ids).union({"4565"}))
    )
    return ExpectedTaxaAuthority(records=records, source_name="expected.tsv", sha256="b" * 64)


@pytest.mark.parametrize(
    ("selector", "identifier", "passing"),
    [
        ("required_exact", "4530", {"G1", "G2", "G4"}),
        ("include_clade", "4479", {"G1", "G2", "G4"}),
        ("only_in_clade", "4479", set()),
        ("exclude_exact", "4113", {"G1", "G3", "G4"}),
        ("exclude_clade", "4081", {"G1", "G3", "G4"}),
    ],
)
def test_individual_selector_semantics(selector: str, identifier: str, passing: set[str]) -> None:
    """Each exact, include, only and exclude selector retains its named meaning."""

    graph = _reviewed_graph()
    selection = make_selection(graph=graph, **{selector: (identifier,)})
    evaluations = evaluate_group_rows(
        graph=graph,
        selection=selection,
        group_species_rows=_group_rows(),
    )
    assert {row.group_id for row in evaluations if row.predicate_pass} == passing


def test_combined_include_and_exact_exclusion_reports_outsiders() -> None:
    """Poaceae plus exact potato exclusion allows and reports other plant hits."""

    graph = _reviewed_graph()
    selection = make_selection(
        graph=graph,
        include_clade=("4479",),
        exclude_exact=("4113",),
    )
    evaluations = evaluate_group_rows(
        graph=graph,
        selection=selection,
        group_species_rows=_group_rows(),
    )
    by_group = {row.group_id: row for row in evaluations}
    assert {row.group_id for row in evaluations if row.predicate_pass} == {"G1", "G4"}
    assert by_group["G1"].outside_selected_scope_taxon_ids == ("3702",)
    assert by_group["G2"].predicate_failure_reasons == (
        "EXCLUDE_EXACT:4113: excluded exact taxon occurs in the group",
    )
    assert by_group["G4"].unmapped_member_count == 1
    assert by_group["G4"].predicate_pass


def test_multiple_include_clades_use_and_and_only_clades_use_intersection() -> None:
    """Several includes need separate contributors while only-in uses intersection."""

    graph = _reviewed_graph()
    includes = make_selection(graph=graph, include_clade=("4479", "3398"))
    evaluation = evaluate_group_taxonomy(
        graph=graph,
        selection=includes,
        run_id="run",
        group_type="HOG",
        hierarchy_node="N0",
        group_id="group",
        species_counts={"Rice": 1},
    )
    assert evaluation.predicate_pass
    only = make_selection(graph=graph, only_in_clade=("33090", "4479"))
    assert evaluate_group_taxonomy(
        graph=graph,
        selection=only,
        run_id="run",
        group_type="HOG",
        hierarchy_node="N0",
        group_id="group",
        species_counts={"Rice": 1},
    ).predicate_pass
    assert not evaluate_group_taxonomy(
        graph=graph,
        selection=only,
        run_id="run",
        group_type="HOG",
        hierarchy_node="N0",
        group_id="group",
        species_counts={"Rice": 1, "Arabidopsis": 1},
    ).predicate_pass
    assert not evaluate_group_taxonomy(
        graph=graph,
        selection=only,
        run_id="run",
        group_type="HOG",
        hierarchy_node="N0",
        group_id="group",
        species_counts={"Rice": 1, "Unknown_input": 1},
    ).predicate_pass


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"required_exact": ("4113",), "exclude_exact": ("4113",)}, "required and excluded"),
        ({"include_clade": ("4479",), "exclude_clade": ("4479",)}, "fully covered"),
        ({"required_exact": ("4113",), "only_in_clade": ("4479",)}, "outside"),
        ({"only_in_clade": ("4479", "4081")}, "no represented terminal intersection"),
        ({"include_clade": ("4081",), "only_in_clade": ("4479",)}, "cannot contribute"),
        ({"include_clade": ("999999",)}, "Unknown"),
    ],
)
def test_selection_contradictions_fail_before_group_evaluation(
    kwargs: dict[str, tuple[str, ...]], message: str
) -> None:
    """Impossible predicates fail closed with an actionable pre-query error."""

    with pytest.raises(InputValidationError, match=message):
        make_selection(graph=_reviewed_graph(), **kwargs)


def test_graph_rejects_inconsistent_metadata_and_mixed_releases() -> None:
    """The same identifier cannot silently change name, rank, parent or release."""

    graph = _reviewed_graph()
    assert "Project Cultivar Registry:CULT-RED" in graph.terminal_ids
    base = _record(
        label="Rice",
        name="Oryza sativa",
        taxon_id=4530,
        lineage=((1, "root", "no rank"), (4479, "Poaceae", "family")),
    )
    conflict = _record(
        label="Maize",
        name="Zea mays",
        taxon_id=4577,
        lineage=((1, "different root", "no rank"), (4479, "Poaceae", "family")),
    )
    with pytest.raises(InputValidationError, match="inconsistent name"):
        build_taxonomy_graph(
            authority=TaxonomyAuthority(
                records=(base, conflict),
                expected_species=("Maize", "Rice"),
            )
        )
    release_conflict = replace(
        conflict,
        lineage_names=("root", "Poaceae"),
        taxonomy_release="2026-09",
    )
    with pytest.raises(InputValidationError, match="mixes releases"):
        build_taxonomy_graph(
            authority=TaxonomyAuthority(
                records=(base, release_conflict),
                expected_species=("Maize", "Rice"),
            )
        )


def test_expected_taxa_contract_and_default_round_trip(tmp_path: Path) -> None:
    """Explicit expected taxa are bounded, terminal-only and checksum-versioned."""

    graph = _reviewed_graph()
    default = default_expected_taxa(graph=graph)
    assert set(default.included_taxon_ids) == set(graph.represented_taxon_ids)
    payload = (
        b"taxon_id\treason\tsource\tincluded\n"
        b"4530\tinput list\tstudy manifest\ttrue\n"
        b"4565\tplanned sample\tstudy manifest\tyes\n"
    )
    parsed = parse_expected_taxa(data=payload, source_name="expected.tsv", graph=graph)
    assert parsed.included_taxon_ids == ("4530", "4565")
    assert b"planned sample" in expected_taxa_to_tsv(authority=parsed)
    path = tmp_path / "expected.tsv"
    path.write_bytes(payload)
    assert read_expected_taxa(path=path, graph=graph) == parsed
    for invalid, message in (
        (payload.replace(b"4565", b"999999"), "Unknown"),
        (payload.replace(b"4565", b"4479"), "not a reviewed terminal"),
        (payload.replace(b"yes", b"maybe"), "invalid included"),
    ):
        with pytest.raises(InputValidationError, match=message):
            parse_expected_taxa(data=invalid, source_name="bad.tsv", graph=graph)


def test_tree_keeps_selection_and_coverage_states_independent() -> None:
    """Expected-no-data, represented outsiders and selectors remain orthogonal."""

    graph = _reviewed_graph()
    selection = make_selection(
        graph=graph,
        include_clade=("4479",),
        exclude_exact=("4113",),
    )
    tree = build_coverage_tree(
        graph=graph,
        selection=selection,
        expected=_expected(graph),
        outside_hit_taxon_ids=("3702",),
    )
    nodes = tree.node_by_id
    assert nodes["4565"].coverage_state == "EXPECTED_NO_DATA"
    assert not nodes["4565"].is_represented
    assert nodes["3702"].coverage_state == "OUTSIDE_SELECTED_SCOPE_WITH_HITS"
    assert nodes["4479"].selection_state == "INCLUDE_CLADE"
    assert nodes["4113"].selection_state == "EXCLUDED_EXACT"
    assert len(tree.edges) == len(tree.nodes) - 1
    assert coverage_records(tree=tree)
    assert tree_edge_records(tree=tree)


def test_compact_tree_never_collapses_selected_or_excluded_nodes() -> None:
    """Compact lineage mode removes only unary neutral internal nodes."""

    graph = _reviewed_graph()
    selection = make_selection(
        graph=graph,
        include_clade=("4479",),
        exclude_clade=("4081",),
    )
    full = build_coverage_tree(
        graph=graph,
        selection=selection,
        expected=_expected(graph),
    )
    compact = build_coverage_tree(
        graph=graph,
        selection=selection,
        expected=_expected(graph),
        compact=True,
    )
    assert len(compact.nodes) <= len(full.nodes)
    assert {"4479", "4081"}.issubset(compact.node_by_id)
    assert len(compact.edges) == len(compact.nodes) - 1


def test_newick_and_renderers_escape_labels_and_are_deterministic() -> None:
    """Apostrophes and markup remain data in Newick, HTML/SVG and figure output."""

    graph = _reviewed_graph()
    altered_nodes = tuple(
        replace(node, name="O'ryza <unsafe>") if node.taxon_id == "4530" else node
        for node in graph.nodes
    )
    graph = replace(graph, nodes=altered_nodes)
    tree = build_coverage_tree(
        graph=graph,
        selection=make_selection(graph=graph),
        expected=_expected(graph),
    )
    newick = tree_to_newick(tree=tree)
    assert "O''ryza <unsafe>" in newick
    assert newick == tree_to_newick(tree=tree)
    svg = tree_to_svg(tree=tree)
    assert b"O&#x27;ryza &lt;unsafe&gt;" in svg
    assert svg == tree_to_svg(tree=tree)
    pdf = tree_to_pdf(tree=tree)
    assert pdf.startswith(b"%PDF-1.4") and pdf.endswith(b"%%EOF\n")
    assert pdf == tree_to_pdf(tree=tree)
    assert coverage_tree_layout(tree=tree) == coverage_tree_layout(tree=tree)
    assert len(coverage_tree_figure(tree=tree).data) == len(tree.edges) + 1
    assert len(style_records(tree=tree)) == len(tree.nodes)


def test_mapping_predicate_summary_and_unmapped_derived_tables() -> None:
    """All required derived audit tables retain exact source labels and meanings."""

    graph = _reviewed_graph()
    selection = make_selection(graph=graph, required_exact=("4530",))
    assert "Require exact taxon Oryza sativa" in selection_summary(selection=selection)[0]
    assert selection_predicate_records(selection=selection)[0]["normalised_logical_mode"] == "AND"
    mappings = input_taxon_mapping_records(graph=graph)
    assert any(row["mapping_status"] == "UNMAPPED" for row in mappings)
    required_bridge_columns = {
        "workflow_species_label",
        "accepted_species_name",
        "ncbi_taxon_id",
        "lineage_taxon_ids",
        "lineage_names",
        "lineage_ranks",
        "taxon_rank",
        "mapping_status",
        "role",
    }
    assert required_bridge_columns.issubset(mappings[0])
    cultivar = next(
        row for row in mappings if row["workflow_species_label"] == "Potato_red_cultivar"
    )
    assert cultivar["ncbi_taxon_id"] == ""
    assert cultivar["authority_taxon_id"] == "CULT-RED"
    assert cultivar["reviewed_terminal_taxon_id"].endswith(":CULT-RED")
    reloaded_bridge = parse_taxonomy_mapping(
        data=records_to_tsv(records=mappings),
        expected_species=graph.represented_labels,
    )
    assert len(reloaded_bridge.all_reviewed_records) == 6
    missing_graph = replace(
        graph,
        input_mappings=(
            replace(
                graph.input_mappings[0],
                mapping_status="MISSING",
                mapping_reason="No reviewed mapping row was supplied.",
            ),
            *graph.input_mappings[1:],
        ),
    )
    missing_bridge = input_taxon_mapping_records(graph=missing_graph)
    assert missing_bridge[0]["mapping_status"] == "UNMAPPED"
    assert missing_bridge[0]["derived_mapping_status"] == "MISSING"
    assert parse_taxonomy_mapping(
        data=records_to_tsv(records=missing_bridge),
        expected_species=graph.represented_labels,
    ).all_reviewed_records
    unresolved = unmapped_input_records(graph=graph, group_species_rows=_group_rows())
    assert unresolved[0]["workflow_species_label"] == "Unknown_input"
    assert unresolved[0]["evaluated_group_member_occurrences"] == 1


def test_group_row_validation_rejects_duplicates_counts_and_limits() -> None:
    """Duplicate membership summaries and unsafe count/limit states fail closed."""

    graph = _reviewed_graph()
    rows = _group_rows()
    with pytest.raises(InputValidationError, match="duplicate row"):
        evaluate_group_rows(
            graph=graph,
            selection=make_selection(graph=graph),
            group_species_rows=(*rows, rows[0]),
        )
    bad = (dict(rows[0], species_member_count=0),)
    with pytest.raises(InputValidationError, match="at least"):
        evaluate_group_rows(
            graph=graph,
            selection=make_selection(graph=graph),
            group_species_rows=bad,
        )
    with pytest.raises(InputValidationError, match="limit"):
        evaluate_group_rows(
            graph=graph,
            selection=make_selection(graph=graph),
            group_species_rows=rows,
            maximum_groups=1,
        )


def test_complete_coverage_exports_reconcile_and_publish_atomically(
    tmp_path: Path,
) -> None:
    """TSV, Newick, JSON, SVG, PDF and checksums share one manifest identity."""

    graph = _reviewed_graph()
    run = build_selection_coverage_run(
        graph=graph,
        expected=_expected(graph),
        selection=make_selection(
            graph=graph,
            include_clade=("4479",),
            exclude_exact=("4113",),
        ),
        group_species_rows=_group_rows(),
        run_id="synthetic",
        resource_identity="run-manifest-sha256:1234",
        package_version="0.7.0",
        focus_authority_name="e3.tsv.gz",
        focus_authority_sha256="c" * 64,
        focus_limited=True,
        created_at_utc="2026-09-09T12:00:00Z",
    )
    assert {row.group_id for row in run.passing_evaluations} == {"G1", "G4"}
    files = coverage_export_files(run=run)
    assert set(files) == set(EXPORT_FILENAMES)
    assert files["selection_coverage_tree.pdf"].startswith(b"%PDF")
    assert b"EXPECTED_NO_DATA" in files["taxon_coverage.tsv"]
    assert b"biological absence" in files["selection_manifest.json"]
    assert coverage_export_zip(run=run) == coverage_export_zip(run=run)
    target = tmp_path / "published"
    manifest = publish_selection_coverage_run(run=run, output_dir=target)
    assert manifest["selection_manifest_id"] == run.selection_manifest_id
    assert tuple(sorted(path.name for path in target.iterdir())) == tuple(sorted(EXPORT_FILENAMES))
    dry_target = tmp_path / "dry"
    assert (
        publish_selection_coverage_run(
            run=run,
            output_dir=dry_target,
            dry_run=True,
        )["summary"]["passing_groups"]
        == 2
    )
    assert not dry_target.exists()
    with pytest.raises(InputValidationError, match="already exists"):
        publish_selection_coverage_run(run=run, output_dir=target)


def test_coverage_run_rejects_bad_provenance_and_timestamp() -> None:
    """Incomplete provenance and non-UTC timestamps fail before publication."""

    graph = _reviewed_graph()
    arguments = {
        "graph": graph,
        "expected": _expected(graph),
        "selection": make_selection(graph=graph),
        "group_species_rows": _group_rows(),
        "run_id": "synthetic",
        "resource_identity": "identity",
        "package_version": "0.7.0",
    }
    with pytest.raises(InputValidationError, match="run_id"):
        build_selection_coverage_run(**dict(arguments, run_id=""))
    with pytest.raises(InputValidationError, match="must be UTC"):
        build_selection_coverage_run(
            **arguments,
            created_at_utc="2026-09-09T12:00:00+01:00",
        )


def test_taxonomy_graph_accessors_fail_closed_for_unknowns_and_cycles() -> None:
    """Unknown node lookups and cyclic parent authorities never return partial paths."""

    graph = _reviewed_graph()
    with pytest.raises(InputValidationError, match="Unknown"):
        graph.descendants(taxon_id="missing")
    with pytest.raises(InputValidationError, match="Unknown"):
        graph.ancestors(taxon_id="missing")
    altered = tuple(
        replace(node, parent_taxon_id="33090") if node.taxon_id == "1" else node
        for node in graph.nodes
    )
    with pytest.raises(InputValidationError, match="cycle"):
        replace(graph, nodes=altered).ancestors(taxon_id="33090")
    with pytest.raises(InputValidationError, match="Unknown selector"):
        make_selection(graph=graph).identifiers(selector_type="OTHER")


def test_taxonomy_graph_rejects_empty_terminal_and_disjoint_roots() -> None:
    """Reviewed rows require stable terminals and one shared displayed root."""

    missing_terminal = _record(
        label="Missing",
        name="Missing terminal",
        taxon_id=None,
        lineage=((1, "root", "no rank"),),
    )
    with pytest.raises(InputValidationError, match="terminal ID"):
        build_taxonomy_graph(
            authority=TaxonomyAuthority(
                records=(missing_terminal,),
                expected_species=("Missing",),
            )
        )
    first = _record(
        label="First",
        name="First species",
        taxon_id=101,
        lineage=((1, "root one", "no rank"),),
    )
    second = _record(
        label="Second",
        name="Second species",
        taxon_id=102,
        lineage=((2, "root two", "no rank"),),
    )
    with pytest.raises(InputValidationError, match="share one displayed root"):
        build_taxonomy_graph(
            authority=TaxonomyAuthority(
                records=(first, second),
                expected_species=("First", "Second"),
            )
        )
    unresolved = _record(
        label="Unknown",
        name="",
        taxon_id=None,
        lineage=(),
        status="UNMAPPED",
    )
    with pytest.raises(InputValidationError, match="requires reviewed"):
        build_taxonomy_graph(
            authority=TaxonomyAuthority(
                records=(unresolved,),
                expected_species=("Unknown",),
            )
        )


def test_missing_input_mapping_and_expected_file_io_are_auditable(tmp_path: Path) -> None:
    """Missing input mapping rows and unavailable expected files remain explicit."""

    record = _record(
        label="Rice",
        name="Oryza sativa",
        taxon_id=4530,
        lineage=((1, "root", "no rank"),),
    )
    graph = build_taxonomy_graph(
        authority=TaxonomyAuthority(
            records=(record,),
            expected_species=("Missing_input", "Rice"),
        )
    )
    missing = graph.mapping_by_label["Missing_input"]
    assert missing.mapping_status == "MISSING"
    with pytest.raises(InputValidationError, match="unavailable"):
        read_expected_taxa(path=tmp_path / "missing.tsv", graph=graph)


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        (b"", "empty"),
        (b"taxon_id\treason\tsource\tincluded\n4530\twhy\tsrc\ttrue\x00\n", "unsafe"),
        (b"taxon_id\treason\tsource\tincluded\n4530\t\xff\tsrc\ttrue\n", "UTF-8"),
        (b"taxon_id\treason\tsource\n4530\twhy\tsrc\n", "unique columns"),
        (b"taxon_id\treason\tsource\tincluded\n4530\twhy\tsrc\ttrue\textra\n", "extra"),
        (b"taxon_id\treason\tsource\tincluded\n4530\t\tsrc\ttrue\n", "provenance"),
        (b"taxon_id\treason\tsource\tincluded\n", "1–"),
        (b"taxon_id\treason\tsource\tincluded\n4530\twhy\tsrc\tfalse\n", "no included"),
        (
            b"taxon_id\treason\tsource\tincluded\n4530\twhy\tsrc\ttrue\n4530\tagain\tsrc\ttrue\n",
            "duplicate",
        ),
    ],
)
def test_expected_taxa_parser_rejects_malformed_universes(payload: bytes, message: str) -> None:
    """Expected-taxon manifests are complete, bounded and unambiguous."""

    with pytest.raises(InputValidationError, match=message):
        parse_expected_taxa(
            data=payload,
            source_name="bad.tsv",
            graph=_reviewed_graph(),
        )


def test_expected_taxa_parser_enforces_record_bound(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The expected universe fails closed before an excessive tree is built."""

    monkeypatch.setattr(taxonomy_selection, "MAX_EXPECTED_TAXA", 1)
    with pytest.raises(InputValidationError, match="1–1"):
        parse_expected_taxa(
            data=(
                b"taxon_id\treason\tsource\tincluded\n4530\twhy\tsrc\ttrue\n4565\twhy\tsrc\ttrue\n"
            ),
            source_name="large.tsv",
            graph=_reviewed_graph(),
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"required_exact": ("4479",)}, "terminal"),
        ({"required_exact": ("4530", "4530")}, "duplicate"),
        ({"required_exact": "4530"}, "sequence"),
        ({"required_exact": ("",)}, "empty"),
        ({"required_exact": ("0",)}, "positive"),
        ({"required_exact": (1,)}, "must be text"),
        (
            {"required_exact": ("4530",), "exclude_clade": ("4479",)},
            "covered by an exclusion",
        ),
        (
            {
                "include_clade": ("4479",),
                "exclude_exact": ("4530", "4577", "4565"),
            },
            "fully covered",
        ),
        (
            {"only_in_clade": ("4479",), "exclude_clade": ("4479",)},
            "complete only-in scope",
        ),
    ],
)
def test_additional_selector_validation_is_fail_closed(
    kwargs: dict[str, object], message: str
) -> None:
    """Terminal roles, identifier shapes and indirect conflicts are pre-query errors."""

    with pytest.raises(InputValidationError, match=message):
        make_selection(graph=_reviewed_graph(), **kwargs)  # type: ignore[arg-type]


def test_group_evaluation_validates_identity_species_counts_and_bounds() -> None:
    """Malformed group summary rows cannot alter predicate outcomes."""

    graph = _reviewed_graph()
    selection = make_selection(graph=graph)
    base = {
        "graph": graph,
        "selection": selection,
        "run_id": "run",
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "group_id": "G",
        "species_counts": {"Rice": 1},
    }
    with pytest.raises(InputValidationError, match="identity"):
        evaluate_group_taxonomy(**dict(base, run_id=1))
    with pytest.raises(InputValidationError, match="requires group"):
        evaluate_group_taxonomy(**dict(base, group_id=""))
    with pytest.raises(InputValidationError, match="species labels"):
        evaluate_group_taxonomy(**dict(base, species_counts={"": 1}))
    with pytest.raises(InputValidationError, match="must be integers"):
        evaluate_group_taxonomy(**dict(base, species_counts={"Rice": 1.5}))
    for maximum in (True, 0, 250_001):
        with pytest.raises(InputValidationError, match="maximum_groups"):
            evaluate_group_rows(
                graph=graph,
                selection=selection,
                group_species_rows=(),
                maximum_groups=maximum,  # type: ignore[arg-type]
            )


def test_tree_bounds_empty_universe_and_outside_hits_are_validated() -> None:
    """Unsafe tree limits and unreviewed outsider claims fail before rendering."""

    graph = _reviewed_graph()
    empty_expected = ExpectedTaxaAuthority(records=(), source_name="none", sha256="x")
    for maximum in (True, 0, 5_001):
        with pytest.raises(InputValidationError, match="maximum_nodes"):
            build_coverage_tree(
                graph=graph,
                selection=make_selection(graph=graph),
                expected=_expected(graph),
                maximum_nodes=maximum,  # type: ignore[arg-type]
            )
    with pytest.raises(InputValidationError, match="unknown or non-terminal"):
        build_coverage_tree(
            graph=graph,
            selection=make_selection(graph=graph),
            expected=ExpectedTaxaAuthority(
                records=(ExpectedTaxon("4479", "why", "src"),),
                source_name="bad",
                sha256="x",
            ),
        )
    with pytest.raises(InputValidationError, match="Outside-scope"):
        build_coverage_tree(
            graph=graph,
            selection=make_selection(graph=graph),
            expected=_expected(graph),
            outside_hit_taxon_ids=("4565",),
        )
    with pytest.raises(InputValidationError, match="empty"):
        build_coverage_tree(
            graph=replace(graph, represented_labels=()),
            selection=make_selection(graph=graph),
            expected=empty_expected,
        )
    with pytest.raises(InputValidationError, match="requires"):
        build_coverage_tree(
            graph=graph,
            selection=make_selection(graph=graph),
            expected=_expected(graph),
            maximum_nodes=1,
        )


def test_least_common_root_and_tree_reconciliation_detect_corruption() -> None:
    """Disconnected, contradictory and miscounted derived trees cannot be exported."""

    with pytest.raises(InputValidationError, match="empty paths"):
        _least_common_root(paths=())
    with pytest.raises(InputValidationError, match="no common root"):
        _least_common_root(paths=(("a",), ("b",)))
    graph = _reviewed_graph()
    expected = _expected(graph)
    tree = build_coverage_tree(
        graph=graph,
        selection=make_selection(graph=graph),
        expected=expected,
    )
    terminals = set(graph.represented_taxon_ids).union(expected.included_taxon_ids)
    with pytest.raises(InputValidationError, match="root is missing"):
        _reconcile_tree(tree=replace(tree, root_taxon_id="missing"), displayed_terminals=terminals)
    with pytest.raises(InputValidationError, match="disconnected or contains duplicate"):
        _reconcile_tree(tree=replace(tree, edges=tree.edges[:-1]), displayed_terminals=terminals)
    duplicate_edges = (tree.edges[0], tree.edges[0], *tree.edges[2:])
    with pytest.raises(InputValidationError, match="duplicate edges"):
        _reconcile_tree(
            tree=replace(tree, edges=duplicate_edges),
            displayed_terminals=terminals,
        )
    bad_edge = replace(tree.edges[0], child_taxon_id="missing")
    with pytest.raises(InputValidationError, match="undisplayed"):
        _reconcile_tree(
            tree=replace(tree, edges=(bad_edge, *tree.edges[1:])),
            displayed_terminals=terminals,
        )
    child_id = tree.edges[0].child_taxon_id
    wrong_parent = next(
        node.taxon_id for node in tree.nodes if node.taxon_id != tree.edges[0].parent_taxon_id
    )
    altered_nodes = tuple(
        replace(node, parent_taxon_id=wrong_parent) if node.taxon_id == child_id else node
        for node in tree.nodes
    )
    with pytest.raises(InputValidationError, match="disagree"):
        _reconcile_tree(
            tree=replace(tree, nodes=altered_nodes),
            displayed_terminals=terminals,
        )
    with pytest.raises(InputValidationError, match="omitted"):
        _reconcile_tree(tree=tree, displayed_terminals=terminals.union({"missing"}))
    wrong_count = replace(
        tree.nodes[0],
        represented_terminal_count=tree.nodes[0].represented_terminal_count + 1,
    )
    with pytest.raises(InputValidationError, match="counts do not reconcile"):
        _reconcile_tree(
            tree=replace(tree, nodes=(wrong_count, *tree.nodes[1:])),
            displayed_terminals=terminals,
        )


def test_renderer_styles_cover_only_scope_expected_and_invalid_colours() -> None:
    """Accessible styles remain consistent for only-in, no-data and invalid colours."""

    from orthofinder_interrogation_app import coverage_tree_render

    graph = _reviewed_graph()
    tree = build_coverage_tree(
        graph=graph,
        selection=make_selection(graph=graph, only_in_clade=("4479",)),
        expected=_expected(graph),
    )
    styles = {row["taxon_id"]: row for row in style_records(tree=tree)}
    assert styles["4565"]["node_outline_style"] == "dashed"
    assert styles["4530"]["incoming_edge_colour"] == "#0072B2"
    pdf = tree_to_pdf(tree=tree)
    assert b"[3 2] 0 d" in pdf
    for colour in ("bad", "#GG0000"):
        with pytest.raises(InputValidationError, match="colour"):
            coverage_tree_render._hex_rgb(value=colour)


def test_export_contract_io_failures_and_timestamp_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Publication cleans staging and reports format, checksum and filesystem failures."""

    from orthofinder_interrogation_app import coverage_exports

    graph = _reviewed_graph()
    run = build_selection_coverage_run(
        graph=graph,
        expected=_expected(graph),
        selection=make_selection(graph=graph),
        group_species_rows=_group_rows(),
        run_id="synthetic",
        resource_identity="identity",
        package_version="0.7.0",
        created_at_utc="2026-09-09T12:00:00Z",
    )
    monkeypatch.setattr(coverage_exports, "EXPORT_FILENAMES", ())
    with pytest.raises(InputValidationError, match="contract"):
        coverage_export_files(run=run)
    monkeypatch.undo()
    blocked = tmp_path / "blocked"
    blocked.write_text("file", encoding="utf-8")
    with pytest.raises(InputValidationError, match="parent cannot be created"):
        publish_selection_coverage_run(run=run, output_dir=blocked / "result")
    with pytest.raises(InputValidationError, match="not valid ISO"):
        build_selection_coverage_run(
            graph=graph,
            expected=_expected(graph),
            selection=make_selection(graph=graph),
            group_species_rows=(),
            run_id="synthetic",
            resource_identity="identity",
            package_version="0.7.0",
            created_at_utc="not-a-date",
        )


def test_atomic_publication_rejects_checksum_mismatch_and_write_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Staging is removed after either corrupted reads or filesystem write errors."""

    graph = _reviewed_graph()
    run = build_selection_coverage_run(
        graph=graph,
        expected=_expected(graph),
        selection=make_selection(graph=graph),
        group_species_rows=(),
        run_id="synthetic",
        resource_identity="identity",
        package_version="0.7.0",
        created_at_utc="2026-09-09T12:00:00Z",
    )
    checksum_target = tmp_path / "checksum_failure"
    original_read_bytes = Path.read_bytes

    def corrupt_staging_read(path: Path) -> bytes:
        """Return altered bytes only while publication verifies staging."""

        payload = original_read_bytes(path)
        return payload + b"corrupt" if ".incoming." in str(path) else payload

    with monkeypatch.context() as patch:
        patch.setattr(Path, "read_bytes", corrupt_staging_read)
        with pytest.raises(InputValidationError, match="Checksum verification failed"):
            publish_selection_coverage_run(run=run, output_dir=checksum_target)
    assert not checksum_target.exists()
    assert not tuple(tmp_path.glob(".checksum_failure.incoming.*"))

    write_target = tmp_path / "write_failure"
    original_write_bytes = Path.write_bytes

    def fail_staging_write(path: Path, data: bytes) -> int:
        """Raise a representative filesystem error only inside staging."""

        if ".incoming." in str(path):
            raise OSError("synthetic write failure")
        return original_write_bytes(path, data)

    with monkeypatch.context() as patch:
        patch.setattr(Path, "write_bytes", fail_staging_write)
        with pytest.raises(InputValidationError, match="publication failed"):
            publish_selection_coverage_run(run=run, output_dir=write_target)
    assert not write_target.exists()
    assert not tuple(tmp_path.glob(".write_failure.incoming.*"))

"""Tests for reviewed taxonomy parsing, statistics and DuckDB searches."""

from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

from orthofinder_interrogation_app import taxonomy
from orthofinder_interrogation_app.models import TaxonomySearchFilters
from orthofinder_interrogation_app.queries import (
    OrthoFinderQueryService,
    _taxonomy_filter_condition,
    _taxonomy_summary_query,
)
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_interrogation_app.taxonomy import (
    TAXONOMY_COLUMNS,
    add_fisher_enrichment,
    benjamini_hochberg,
    parse_taxonomy_mapping,
    read_taxonomy_mapping,
    taxonomy_audit_rows,
    taxonomy_template,
)
from orthofinder_results.errors import InputValidationError


def _authority(*, path: Path):
    """Return the fixture mapping validated against all four run labels."""

    return read_taxonomy_mapping(
        path=path,
        expected_species=("Species_A", "Species_B", "Species_C", "Species_D"),
    )


def _service(*, resource: Path) -> OrthoFinderQueryService:
    """Return the real read-only query service for the compact application resource."""

    return OrthoFinderQueryService(resource=open_resource(path=resource))


def _rewrite_mapping(*, source: Path, change) -> bytes:
    """Return a fixture mapping after applying a row-level test mutation."""

    rows = list(csv.DictReader(io.StringIO(source.read_text()), delimiter="\t"))
    change(rows)
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output, fieldnames=TAXONOMY_COLUMNS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode()


def test_reviewed_mapping_descendants_options_audit_and_template(
    taxonomy_mapping_file: Path,
) -> None:
    """Reviewed decisions alone drive descendant targets and coverage reporting."""

    authority = _authority(path=taxonomy_mapping_file)
    assert authority.summary() == {
        "REVIEWED": 3,
        "PENDING_REVIEW": 0,
        "UNMAPPED": 0,
        "AMBIGUOUS": 1,
        "MISSING": 0,
    }
    assert authority.reviewed_species == ("Species_A", "Species_B", "Species_C")
    assert authority.unresolved_species == ("Species_D",)
    assert authority.target_species(taxon_id=10) == ("Species_A", "Species_B")
    assert authority.target_species(taxon_id=101) == ("Species_A",)
    with pytest.raises(InputValidationError, match="positive integer"):
        authority.target_species(taxon_id=0)
    options = {row.taxon_id: row for row in authority.taxon_options()}
    assert options[10].display_label() == "Target clade | NCBI taxon 10"
    audit = taxonomy_audit_rows(authority=authority)
    assert audit[-1]["mapping_status"] == "AMBIGUOUS"

    partial = parse_taxonomy_mapping(
        data=b"\t".join(column.encode() for column in TAXONOMY_COLUMNS)
        + b"\nSpecies_A\tSpecies alpha\t\t\t\t\t\t\tUNMAPPED\t\t\t\t\t\t\t\n",
        expected_species=("Species_A", "Species_B"),
    )
    assert partial.summary()["MISSING"] == 1
    assert taxonomy_audit_rows(authority=partial)[1]["mapping_status"] == "MISSING"
    template = taxonomy_template(species=("Species_B", "Species_A"))
    template_authority = parse_taxonomy_mapping(
        data=template, expected_species=("Species_A", "Species_B")
    )
    assert template_authority.summary()["UNMAPPED"] == 2


def test_contains_exclusive_near_exclusive_and_enrichment_queries(
    application_resource: Path, taxonomy_mapping_file: Path
) -> None:
    """All four modes retain exact, sampled-universe semantics."""

    service = _service(resource=application_resource)
    authority = _authority(path=taxonomy_mapping_file)
    common = {
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "target_taxon_id": 10,
        "sort_mode": "unused",
    }
    del common["sort_mode"]
    contains = service.search_taxonomy_groups(
        authority=authority,
        filters=TaxonomySearchFilters(mode="CONTAINS", **common),
    )
    assert [row["group_id"] for row in contains.rows] == [
        "N0.HOG1",
        "N0.HOG3",
        "N0.HOG2",
    ]
    assert contains.target_species_count == 2
    assert contains.outside_species_count == 1
    assert contains.unresolved_species_count == 1

    exclusive = service.search_taxonomy_groups(
        authority=authority,
        filters=TaxonomySearchFilters(mode="SAMPLED_EXCLUSIVE", **common),
    )
    assert [row["group_id"] for row in exclusive.rows] == ["N0.HOG1"]
    near = service.search_taxonomy_groups(
        authority=authority,
        filters=TaxonomySearchFilters(
            mode="NEAR_EXCLUSIVE",
            maximum_outside_species_count=1,
            maximum_unresolved_species_count=0,
            **common,
        ),
    )
    assert near.total_rows == 3
    hog_two = next(row for row in near.rows if row["group_id"] == "N0.HOG2")
    assert hog_two["outsider_species_labels"] == "Species_C"
    assert hog_two["unresolved_species_labels"] == ""

    enriched = service.search_taxonomy_groups(
        authority=authority,
        filters=TaxonomySearchFilters(
            mode="ENRICHED",
            maximum_q_value=1.0,
            minimum_odds_ratio=0.0,
            page_size=2,
            **common,
        ),
    )
    assert enriched.tested_group_count == 3
    assert enriched.total_rows == 3
    assert len(enriched.rows) == 2
    assert all("enrichment_q_value" in row for row in enriched.rows)


def test_taxonomy_query_rejects_wrong_authority_empty_target_and_excessive_universe(
    application_resource: Path,
    taxonomy_mapping_file: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Mapping/resource mismatches and unsafe enrichment materialisation fail visibly."""

    service = _service(resource=application_resource)
    authority = _authority(path=taxonomy_mapping_file)
    filters = TaxonomySearchFilters(group_type="HOG", hierarchy_node="N0", target_taxon_id=10)
    wrong = parse_taxonomy_mapping(
        data=taxonomy_template(species=("Other",)), expected_species=("Other",)
    )
    with pytest.raises(InputValidationError, match="exact species set"):
        service.search_taxonomy_groups(authority=wrong, filters=filters)
    with pytest.raises(InputValidationError, match="No reviewed sampled species"):
        service.search_taxonomy_groups(
            authority=authority,
            filters=TaxonomySearchFilters(
                group_type="HOG", hierarchy_node="N0", target_taxon_id=999
            ),
        )
    monkeypatch.setattr("orthofinder_interrogation_app.queries.MAX_TAXONOMY_TEST_GROUPS", 2)
    with pytest.raises(InputValidationError, match="interactive limit"):
        service.search_taxonomy_groups(
            authority=authority,
            filters=TaxonomySearchFilters(
                group_type="HOG",
                hierarchy_node="N0",
                target_taxon_id=10,
                mode="ENRICHED",
            ),
        )


def test_fisher_enrichment_and_bh_validation() -> None:
    """Species-presence tests are corrected stably and reject invalid universes."""

    rows = (
        {"group_id": "a", "target_species_count": 2, "outside_species_count": 0},
        {"group_id": "b", "target_species_count": 1, "outside_species_count": 1},
    )
    tested = add_fisher_enrichment(rows=rows, total_target_species=2, total_outside_species=2)
    assert tested[0]["enrichment_p_value"] < tested[1]["enrichment_p_value"]
    assert tested[0]["enrichment_q_value"] <= tested[1]["enrichment_q_value"]
    assert benjamini_hochberg(p_values=()) == ()
    assert benjamini_hochberg(p_values=(0.01, 0.04, 0.03)) == pytest.approx((0.03, 0.04, 0.04))
    with pytest.raises(InputValidationError, match="one reviewed target"):
        add_fisher_enrichment(rows=rows, total_target_species=0, total_outside_species=2)
    with pytest.raises(InputValidationError, match="p-values"):
        benjamini_hochberg(p_values=(float("nan"),))
    with pytest.raises(InputValidationError, match="exceeds the universe"):
        add_fisher_enrichment(
            rows=({"target_species_count": 3, "outside_species_count": 0},),
            total_target_species=2,
            total_outside_species=2,
        )
    with pytest.raises(InputValidationError, match="outside-species"):
        add_fisher_enrichment(
            rows=({"target_species_count": 1, "outside_species_count": 3},),
            total_target_species=2,
            total_outside_species=2,
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda rows: rows.append(dict(rows[0])), "duplicate workflow"),
        (lambda rows: rows[0].update(mapping_status="GUESSED"), "unsupported"),
        (lambda rows: rows[0].update(ncbi_taxon_id="bad"), "non-integer"),
        (lambda rows: rows[0].update(ncbi_taxon_id="0"), "non-positive"),
        (lambda rows: rows[0].update(lineage_names="root"), "unequal lineage"),
        (
            lambda rows: rows[0].update(
                lineage_taxon_ids="1;10;101",
                lineage_names="cellular organisms;Target clade;Species alpha",
            ),
            "own taxon",
        ),
        (lambda rows: rows[0].update(parent_taxon_id="1"), "final lineage"),
        (lambda rows: rows[0].update(parent_taxon_name=""), "together"),
        (lambda rows: rows[0].update(source_date="yesterday"), "provenance date"),
        (lambda rows: rows[0].update(reviewed_at_utc="2026-09-07"), "timezone"),
        (lambda rows: rows[0].update(lineage_taxon_ids="1;1"), "duplicate lineage"),
        (lambda rows: rows[0].update(lineage_taxon_ids="1;;10"), "empty lineage"),
        (lambda rows: rows[0].update(accepted_species_name=""), "lacks"),
    ],
)
def test_mapping_row_validation(taxonomy_mapping_file: Path, mutation, message: str) -> None:
    """Malformed review decisions fail with controlled field-specific errors."""

    data = _rewrite_mapping(source=taxonomy_mapping_file, change=mutation)
    with pytest.raises(InputValidationError, match=message):
        parse_taxonomy_mapping(
            data=data,
            expected_species=("Species_A", "Species_B", "Species_C", "Species_D"),
        )


def test_mapping_file_and_payload_boundaries(tmp_path: Path, taxonomy_mapping_file: Path) -> None:
    """Absent files, invalid encodings/headings and foreign labels are rejected."""

    with pytest.raises(InputValidationError, match="unavailable"):
        read_taxonomy_mapping(path=tmp_path / "missing.tsv", expected_species=())
    with pytest.raises(InputValidationError, match="1–"):
        parse_taxonomy_mapping(data=b"", expected_species=())
    with pytest.raises(InputValidationError, match="NUL"):
        parse_taxonomy_mapping(data=b"a\x00b", expected_species=())
    with pytest.raises(InputValidationError, match="UTF-8"):
        parse_taxonomy_mapping(data=b"\xff", expected_species=())
    with pytest.raises(InputValidationError, match="lacks required"):
        parse_taxonomy_mapping(data=b"wrong\nvalue\n", expected_species=())
    duplicate_heading = ("\t".join((*TAXONOMY_COLUMNS, "workflow_species_label")) + "\n").encode()
    with pytest.raises(InputValidationError, match="duplicate column"):
        parse_taxonomy_mapping(data=duplicate_heading, expected_species=())
    extra_field = (
        "\t".join(TAXONOMY_COLUMNS)
        + "\n"
        + "\t".join(("Species_A", *("" for _ in TAXONOMY_COLUMNS), "extra"))
        + "\n"
    ).encode()
    with pytest.raises(InputValidationError, match="extra tab"):
        parse_taxonomy_mapping(data=extra_field, expected_species=("Species_A",))
    foreign = _rewrite_mapping(
        source=taxonomy_mapping_file,
        change=lambda rows: rows[0].update(workflow_species_label="Other"),
    )
    with pytest.raises(InputValidationError, match="absent from this resource"):
        parse_taxonomy_mapping(
            data=foreign,
            expected_species=("Species_A", "Species_B", "Species_C", "Species_D"),
        )
    with pytest.raises(InputValidationError, match="unique and non-empty"):
        parse_taxonomy_mapping(
            data=taxonomy_mapping_file.read_bytes(),
            expected_species=("Species_A", "Species_A"),
        )
    empty_file = tmp_path / "empty.tsv"
    empty_file.touch()
    with pytest.raises(InputValidationError, match="size must"):
        read_taxonomy_mapping(path=empty_file, expected_species=())
    old_limit = taxonomy.MAX_TAXONOMY_BYTES
    taxonomy.MAX_TAXONOMY_BYTES = 1
    try:
        with pytest.raises(InputValidationError, match="size must"):
            read_taxonomy_mapping(
                path=taxonomy_mapping_file,
                expected_species=(
                    "Species_A",
                    "Species_B",
                    "Species_C",
                    "Species_D",
                ),
            )
    finally:
        taxonomy.MAX_TAXONOMY_BYTES = old_limit


def test_mapping_conflicting_taxon_names_are_rejected(
    taxonomy_mapping_file: Path,
) -> None:
    """One NCBI identifier cannot silently acquire multiple lineage names."""

    data = _rewrite_mapping(
        source=taxonomy_mapping_file,
        change=lambda rows: rows[1].update(
            lineage_names="cellular organisms;Different target name"
        ),
    )
    authority = parse_taxonomy_mapping(
        data=data,
        expected_species=("Species_A", "Species_B", "Species_C", "Species_D"),
    )
    with pytest.raises(InputValidationError, match="conflicting names"):
        authority.taxon_options()


def test_internal_taxonomy_sql_builders_reject_empty_or_statistical_modes() -> None:
    """Dynamic SQL builders accept only mapped labels and non-statistical modes."""

    enriched = TaxonomySearchFilters(
        group_type="HOG",
        hierarchy_node="N0",
        target_taxon_id=10,
        mode="ENRICHED",
    )
    with pytest.raises(InputValidationError, match="requires mapped"):
        _taxonomy_summary_query(
            run_id="run",
            filters=enriched,
            categories=(),
            target_species_count=0,
        )
    with pytest.raises(InputValidationError, match="does not support"):
        _taxonomy_filter_condition(filters=enriched)


def test_results_feb26_example_mapping_is_complete_and_auditable() -> None:
    """The shipped 60-species mapping preserves its one unresolved workflow typo."""

    path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "results_feb26_ncbi_taxonomy_mapping_20260907.tsv"
    )
    with path.open(encoding="utf-8", newline="") as handle:
        rows = tuple(csv.DictReader(handle, delimiter="\t"))
    expected = tuple(sorted(row["workflow_species_label"] for row in rows))
    authority = read_taxonomy_mapping(path=path, expected_species=expected)
    assert len(expected) == 60
    assert authority.summary() == {
        "REVIEWED": 59,
        "PENDING_REVIEW": 0,
        "UNMAPPED": 1,
        "AMBIGUOUS": 0,
        "MISSING": 0,
    }
    assert authority.unresolved_species == ("Leismania_major",)
    assert len(authority.target_species(taxon_id=33208)) == 25

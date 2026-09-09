"""Tests for replaceable focus authorities and aggregated cluster searches."""

from __future__ import annotations

from pathlib import Path

import duckdb
import pytest

from orthofinder_interrogation_app import focus, focus_page
from orthofinder_interrogation_app.focus import (
    bundled_focus_path,
    focus_template,
    parse_focus_proteins,
    read_focus_proteins,
)
from orthofinder_interrogation_app.models import FocusClusterFilters
from orthofinder_interrogation_app.queries import OrthoFinderQueryService
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_results.errors import InputValidationError


def test_bundled_e3_focus_authority_is_versioned_and_complete() -> None:
    """The reviewed 1,000-seed E3 catalogue is the exact packaged default."""

    authority = read_focus_proteins(path=bundled_focus_path())
    assert len(authority.records) == 1_000
    assert len(authority.identifiers) == 1_000
    assert authority.sha256 == (
        "10945b7cae2212f3bc425bd80a37481c96f3194742724bdabe41c0148e73db52"
    )
    assert "A0A060D0U3" in authority.identifiers
    first = next(row for row in authority.records if row.identifier == "A0A060D0U3")
    assert first.sequence_identifiers == "tr|A0A060D0U3|A0A060D0U3_MAIZE"
    assert first.review_status == "unreviewed"
    assert first.category == "Ring finger"
    assert authority.summary()["organisms"] > 1


def test_focus_template_and_controlled_aliases_are_user_editable() -> None:
    """Canonical and inherited headings parse without silent alias precedence."""

    template = parse_focus_proteins(data=focus_template(), source_name="custom.tsv")
    assert template.identifiers == ("Q9SA03",)
    inherited = parse_focus_proteins(
        data=(
            b"accession\tevidence_type\tsource\n"
            b"A1\tRING\treviewed\nB2\tHECT\treviewed\n"
        ),
        source_name="known.tsv",
    )
    assert inherited.identifiers == ("A1", "B2")
    with pytest.raises(InputValidationError, match="multiple aliases"):
        parse_focus_proteins(
            data=b"accession\tprotein_identifier\nA1\tA1\n",
            source_name="ambiguous.tsv",
        )


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        (b"name\nenabled_only\n", "identifier"),
        (b"protein_identifier\tenabled\nA1\tmaybe\n", "enabled"),
        (b"protein_identifier\nA1\nA1\n", "duplicate"),
        (b"protein_identifier\nA1\na1\n", "case-colliding"),
    ],
)
def test_focus_authority_rejects_ambiguous_records(
    payload: bytes, message: str
) -> None:
    """Unsafe Boolean, duplicate and case-collision states fail closed."""

    with pytest.raises(InputValidationError, match=message):
        parse_focus_proteins(data=payload, source_name="bad.tsv")


def test_focus_cluster_query_aggregates_canonical_and_internal_matches(
    application_resource: Path,
) -> None:
    """Several exact identifiers map to one cluster without duplicate inflation."""

    service = OrthoFinderQueryService(
        resource=open_resource(path=application_resource)
    )
    result = service.search_focus_clusters(
        filters=FocusClusterFilters(
            protein_identifiers=("alpha_1", "0_1", "absent"),
            group_type="HOG",
            hierarchy_node="N0",
        )
    )
    assert result.total_rows == 1
    assert result.matched_focus_identifiers == 2
    assert result.submitted_focus_identifiers == 3
    assert result.rows[0]["group_id"] == "N0.HOG1"
    assert result.rows[0]["matched_focus_count"] == 2
    assert result.rows[0]["matched_protein_count"] == 2
    assert result.rows[0]["matched_focus_identifiers"] == "0_1;alpha_1"
    assert result.rows[0]["mean_distance"] == pytest.approx(0.1)


def test_focus_cluster_query_resolves_uniprot_accession_and_entry(
    application_resource: Path,
) -> None:
    """Controlled pipe identifiers expose exact accession and entry aliases."""

    database = application_resource / "duckdb" / "orthofinder_results.duckdb"
    connection = duckdb.connect(str(database))
    try:
        connection.execute(
            "UPDATE hog_memberships SET member_id = 'sp|Q9SA03|FB27_ARATH' "
            "WHERE member_id = 'alpha_1'"
        )
        connection.execute(
            "UPDATE sequences SET member_id = 'sp|Q9SA03|FB27_ARATH' "
            "WHERE member_id = 'alpha_1'"
        )
    finally:
        connection.close()
    service = OrthoFinderQueryService(resource=open_resource(path=application_resource))
    result = service.search_focus_clusters(
        filters=FocusClusterFilters(
            protein_identifiers=("Q9SA03", "FB27_ARATH"),
            group_type="HOG",
            hierarchy_node="N0",
        )
    )
    assert result.matched_focus_identifiers == 2
    assert result.rows[0]["matched_protein_count"] == 1
    assert result.rows[0]["matched_member_ids"] == "sp|Q9SA03|FB27_ARATH"
    assert result.rows[0]["match_authorities"] == "UNIPROT_ACCESSION;UNIPROT_ENTRY"


def test_focus_cluster_controls_reject_invalid_authorities(
    application_resource: Path,
) -> None:
    """Bounds, group systems and legacy hierarchy semantics are defensive."""

    with pytest.raises(InputValidationError, match="unique"):
        FocusClusterFilters(
            protein_identifiers=("A", "A"),
            group_type="HOG",
            hierarchy_node="N0",
        )
    service = OrthoFinderQueryService(
        resource=open_resource(path=application_resource)
    )
    with pytest.raises(InputValidationError, match="Unsupported focus"):
        service.search_focus_clusters(
            filters=FocusClusterFilters(
                protein_identifiers=("A",),
                group_type="OTHER",
                hierarchy_node="N0",
            )
        )
    with pytest.raises(InputValidationError, match="ROOT"):
        service.search_focus_clusters(
            filters=FocusClusterFilters(
                protein_identifiers=("A",),
                group_type="LEGACY_ORTHOGROUP",
                hierarchy_node="N0",
            )
        )


def test_group_species_collection_is_exact_bounded_and_deterministic(
    application_resource: Path,
) -> None:
    """Coverage queries can retrieve all or exact focus groups without mutation."""

    service = OrthoFinderQueryService(
        resource=open_resource(path=application_resource)
    )
    focused = service.get_group_species_collection(
        group_type="HOG",
        hierarchy_node="N0",
        group_ids=("N0.HOG2", "N0.HOG1"),
    )
    assert {row["group_id"] for row in focused} == {"N0.HOG1", "N0.HOG2"}
    assert len(focused) == 4
    complete = service.get_group_species_collection(
        group_type="HOG",
        hierarchy_node="N0",
    )
    assert len(complete) == 7
    with pytest.raises(InputValidationError, match="contains 7"):
        service.get_group_species_collection(
            group_type="HOG",
            hierarchy_node="N0",
            maximum_rows=1,
        )
    with pytest.raises(InputValidationError, match="unique"):
        service.get_group_species_collection(
            group_type="HOG",
            hierarchy_node="N0",
            group_ids=("N0.HOG1", "N0.HOG1"),
        )


def test_focus_page_authority_precedence_and_display_helpers(tmp_path: Path) -> None:
    """Upload, path and bundled focus authorities retain exact provenance and labels."""

    path = tmp_path / "focus.tsv"
    path.write_bytes(b"protein_identifier\tcategory\nenabled_A\tRING\n")
    from_path = focus_page.load_focus_authority(
        focus_path_text=str(path),
        uploaded_data=None,
    )
    assert from_path.identifiers == ("enabled_A",)
    uploaded = focus_page.load_focus_authority(
        focus_path_text=str(tmp_path / "missing.tsv"),
        uploaded_data=b"protein_identifier\nUPLOAD\n",
        uploaded_name="uploaded.tsv",
    )
    assert uploaded.identifiers == ("UPLOAD",)
    assert focus_page.load_focus_authority(
        focus_path_text="",
        uploaded_data=None,
    ).identifiers[0]
    record = focus_page._focus_record(row=uploaded.records[0])
    assert record["protein_identifier"] == "UPLOAD"
    displayed = focus_page._display_focus_result(
        row={
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": "N0.HOG1",
            "matched_focus_count": 1,
            "matched_protein_count": 1,
            "matched_focus_identifiers": "UPLOAD",
            "matched_member_ids": "sp|UPLOAD|ENTRY",
            "matched_species_labels": "Species_A",
            "match_authorities": "UNIPROT_ACCESSION",
            "member_count": 3,
            "species_count": 2,
            "computation_status": "EXACT",
            "mean_distance": 0.1,
            "population_stddev_distance": 0.02,
        },
        metadata={"UPLOAD": uploaded.records[0]},
    )
    assert displayed["Matched focus IDs"] == 1
    assert "HOG | N0 | N0.HOG1" in focus_page._focus_option(
        row={
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": "N0.HOG1",
            "matched_focus_count": 1,
        }
    )


def test_focus_file_reader_handles_plain_gzip_and_io_failures(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Focus authorities fail clearly for missing, empty, corrupt and oversized files."""

    plain = tmp_path / "focus.tsv"
    plain.write_bytes(b"protein_identifier\nB2\n")
    assert read_focus_proteins(path=plain).identifiers == ("B2",)
    compressed = tmp_path / "focus.tsv.gz"
    compressed.write_bytes(__import__("gzip").compress(b"protein_identifier\nA1\n"))
    assert read_focus_proteins(path=compressed).identifiers == ("A1",)
    with pytest.raises(InputValidationError, match="unavailable"):
        read_focus_proteins(path=tmp_path / "missing.tsv")
    empty = tmp_path / "empty.tsv"
    empty.write_bytes(b"")
    with pytest.raises(InputValidationError, match="must contain"):
        read_focus_proteins(path=empty)
    corrupt = tmp_path / "corrupt.tsv.gz"
    corrupt.write_bytes(b"not-gzip")
    with pytest.raises(InputValidationError, match="could not be read"):
        read_focus_proteins(path=corrupt)
    expanded = tmp_path / "expanded.tsv.gz"
    expanded.write_bytes(__import__("gzip").compress(b"A" * 2_000))
    monkeypatch.setattr(focus, "MAX_FOCUS_BYTES", 100)
    with pytest.raises(InputValidationError, match="Expanded"):
        read_focus_proteins(path=expanded)
    monkeypatch.setattr(focus, "files", lambda _package: tmp_path)
    with pytest.raises(InputValidationError, match="Bundled"):
        bundled_focus_path()


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        (b"", "empty"),
        (b"protein_identifier\nA1\x00\n", "NUL"),
        (b"protein_identifier\n\xff\n", "UTF-8"),
        (b"protein_identifier\tprotein_identifier\nA1\tA1\n", "headings"),
        (b"protein_identifier\nA1\textra\n", "extra fields"),
        (b"protein_identifier\n\n", "no data"),
        (b"protein_identifier\tenabled\nA1\tfalse\n", "no enabled"),
        (b"protein_identifier\n\"unterminated\n", "could not be parsed"),
    ],
)
def test_focus_parser_rejects_malformed_payloads(payload: bytes, message: str) -> None:
    """Malformed focus files never become a partial identifier authority."""

    with pytest.raises(InputValidationError, match=message):
        parse_focus_proteins(data=payload, source_name="bad.tsv")


def test_focus_parser_bounds_text_rows_and_boolean_aliases(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Record bounds, unsafe values and explicit disabled rows are validated."""

    original_maximum = focus.MAX_FOCUS_RECORDS
    with pytest.raises(InputValidationError, match="empty identifier"):
        parse_focus_proteins(
            data=b"protein_identifier\tprotein_name\n\tmissing\n",
            source_name="bad.tsv",
        )
    with pytest.raises(InputValidationError, match="unsafe text"):
        parse_focus_proteins(
            data=(
                "protein_identifier\tprotein_name\nA1\t"
                + "x" * (focus.MAX_FOCUS_FIELD_CHARACTERS + 1)
                + "\n"
            ).encode(),
            source_name="bad.tsv",
        )
    monkeypatch.setattr(focus, "MAX_FOCUS_RECORDS", 1)
    with pytest.raises(InputValidationError, match="exceeds"):
        parse_focus_proteins(
            data=b"protein_identifier\nA1\nB2\n",
            source_name="bad.tsv",
        )
    monkeypatch.setattr(focus, "MAX_FOCUS_RECORDS", original_maximum)
    authority = parse_focus_proteins(
        data=b"protein_identifier\tenabled\nA1\tyes\nB2\tdisabled\n",
        source_name="",
        source_sha256="checksum",
    )
    assert authority.identifiers == ("A1",)
    assert authority.source_name == "uploaded_focus_proteins.tsv"
    assert authority.sha256 == "checksum"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"protein_identifiers": "A1"}, "sequence"),
        ({"protein_identifiers": ()}, "non-empty"),
        ({"protein_identifiers": ("",)}, "non-empty"),
        ({"protein_identifiers": ("A\n",)}, "unsafe"),
        ({"protein_identifiers": ("A",), "group_type": ""}, "group_type"),
        ({"protein_identifiers": ("A",), "hierarchy_node": 1}, "hierarchy_node"),
        ({"protein_identifiers": ("A",), "maximum_rows": True}, "integer"),
        ({"protein_identifiers": ("A",), "maximum_rows": 0}, "between"),
    ],
)
def test_focus_filter_model_rejects_unsafe_controls(
    kwargs: dict[str, object], message: str
) -> None:
    """Every focus-query input is normalised before DuckDB execution."""

    defaults: dict[str, object] = {
        "protein_identifiers": ("A",),
        "group_type": "HOG",
        "hierarchy_node": "N0",
    }
    defaults.update(kwargs)
    with pytest.raises(InputValidationError, match=message):
        FocusClusterFilters(**defaults)  # type: ignore[arg-type]


def test_focus_path_loader_reports_missing_file(tmp_path: Path) -> None:
    """A missing launcher-configured focus path produces an actionable error."""

    with pytest.raises(InputValidationError, match="unavailable"):
        focus_page.load_focus_authority(
            focus_path_text=str(tmp_path / "missing.tsv"),
            uploaded_data=None,
        )

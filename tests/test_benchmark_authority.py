"""Tests for the reviewed Arabidopsis dispersion-marker authority."""

from __future__ import annotations

from pathlib import Path

import pytest

from orthofinder_results.benchmark_authority import (
    BENCHMARK_FIELDS,
    MAX_BENCHMARK_BYTES,
    bundled_benchmark_path,
    read_benchmark_authority,
)
from orthofinder_results.errors import InputValidationError


def test_bundled_benchmark_authority_is_complete_and_checksum_bound() -> None:
    """Packaged housekeeping and pan-NLRome markers retain exact provenance."""

    authority = read_benchmark_authority(path=bundled_benchmark_path())
    classes = {record.benchmark_class for record in authority.records}
    assert len(authority.records) == 112
    assert classes == {"HOUSEKEEPING", "R_NLR"}
    assert sum(record.benchmark_class == "HOUSEKEEPING" for record in authority.records) == 13
    assert sum(record.benchmark_class == "R_NLR" for record in authority.records) == 99
    assert authority.sha256 == (
        "33f30fcb31966dfde4136a8140bb9db00ada078b7e0ed01da6be73301b9a8b6f"
    )
    focus = authority.to_focus_authority()
    assert focus.identifiers[0]
    assert all(record.taxon_id == "3702" for record in focus.records)
    assert any("doi:10.1016/j.cell.2019.07.038" in record.source for record in focus.records)


def test_benchmark_authority_rejects_schema_taxon_and_duplicate_errors(
    tmp_path: Path,
) -> None:
    """Malformed or biologically conflicting marker files fail closed."""

    heading = "\t".join(BENCHMARK_FIELDS)
    values = {
        "marker_id": "AT1G00001",
        "protein_identifier": "P1",
        "protein_entry": "ONE_ARATH",
        "marker_name": "Marker one",
        "benchmark_class": "HOUSEKEEPING",
        "benchmark_subclass": "REFERENCE_CANDIDATE",
        "domain_architecture": "",
        "evidence_type": "REVIEWED",
        "organism": "Arabidopsis thaliana",
        "taxon_id": "3702",
        "source_title": "Source",
        "source_doi": "",
        "source_table": "Table",
        "source_version": "1",
        "enabled": "true",
        "note": "Evaluate empirically.",
    }

    def row() -> str:
        return "\t".join(values[field] for field in BENCHMARK_FIELDS)

    wrong_heading = tmp_path / "wrong_heading.tsv"
    wrong_heading.write_text("marker_id\nA\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="headings"):
        read_benchmark_authority(path=wrong_heading)

    wrong_taxon = tmp_path / "wrong_taxon.tsv"
    values["taxon_id"] = "9606"
    wrong_taxon.write_text(f"{heading}\n{row()}\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="taxon 3702"):
        read_benchmark_authority(path=wrong_taxon)

    duplicate = tmp_path / "duplicate.tsv"
    values["taxon_id"] = "3702"
    duplicate.write_text(f"{heading}\n{row()}\n{row()}\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="duplicate protein"):
        read_benchmark_authority(path=duplicate)

    invalid_boolean = tmp_path / "invalid_boolean.tsv"
    values["enabled"] = "perhaps"
    invalid_boolean.write_text(f"{heading}\n{row()}\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="invalid enabled"):
        read_benchmark_authority(path=invalid_boolean)


def test_benchmark_authority_rejects_unsafe_and_incomplete_files(
    tmp_path: Path,
) -> None:
    """I/O, row-shape, required-field and class errors remain explicit."""

    heading = "\t".join(BENCHMARK_FIELDS)
    values = {
        "marker_id": "AT1G00001",
        "protein_identifier": "P1",
        "protein_entry": "ONE_ARATH",
        "marker_name": "Marker one",
        "benchmark_class": "HOUSEKEEPING",
        "benchmark_subclass": "REFERENCE_CANDIDATE",
        "domain_architecture": "",
        "evidence_type": "REVIEWED",
        "organism": "Arabidopsis thaliana",
        "taxon_id": "3702",
        "source_title": "Source",
        "source_doi": "",
        "source_table": "",
        "source_version": "1",
        "enabled": "false",
        "note": "Evaluate empirically.",
    }

    def row() -> str:
        return "\t".join(values[field] for field in BENCHMARK_FIELDS)

    with pytest.raises(InputValidationError, match="unavailable"):
        read_benchmark_authority(path=tmp_path / "missing.tsv")

    for name, content in (
        ("empty.tsv", b""),
        ("oversized.tsv", b"x" * (MAX_BENCHMARK_BYTES + 1)),
    ):
        path = tmp_path / name
        path.write_bytes(content)
        with pytest.raises(InputValidationError, match="must contain"):
            read_benchmark_authority(path=path)

    invalid_utf8 = tmp_path / "invalid_utf8.tsv"
    invalid_utf8.write_bytes(b"\xff\xfe\xfd")
    with pytest.raises(InputValidationError, match="readable UTF-8"):
        read_benchmark_authority(path=invalid_utf8)

    nul = tmp_path / "nul.tsv"
    nul.write_bytes(f"{heading}\n{row()}\n".encode() + b"\x00")
    with pytest.raises(InputValidationError, match="NUL byte"):
        read_benchmark_authority(path=nul)

    extra = tmp_path / "extra.tsv"
    extra.write_text(f"{heading}\n{row()}\textra\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="extra fields"):
        read_benchmark_authority(path=extra)

    no_records = tmp_path / "no_records.tsv"
    no_records.write_text(f"{heading}\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="must contain"):
        read_benchmark_authority(path=no_records)

    missing_name = tmp_path / "missing_name.tsv"
    values["marker_name"] = ""
    missing_name.write_text(f"{heading}\n{row()}\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="lacks: marker_name"):
        read_benchmark_authority(path=missing_name)

    unsupported = tmp_path / "unsupported.tsv"
    values["marker_name"] = "Marker one"
    values["benchmark_class"] = "OTHER"
    unsupported.write_text(f"{heading}\n{row()}\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="unsupported class"):
        read_benchmark_authority(path=unsupported)

    unsafe = tmp_path / "unsafe.tsv"
    values["benchmark_class"] = "HOUSEKEEPING"
    values["note"] = "unsafe\x1ftext"
    unsafe.write_text(f"{heading}\n{row()}\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="unsafe text"):
        read_benchmark_authority(path=unsafe)


def test_benchmark_authority_rejects_conflicting_marker_classification(
    tmp_path: Path,
) -> None:
    """One locus cannot silently acquire two benchmark classifications."""

    heading = "\t".join(BENCHMARK_FIELDS)
    common = {
        "marker_id": "AT1G00001",
        "protein_entry": "ONE_ARATH",
        "marker_name": "Marker one",
        "domain_architecture": "",
        "evidence_type": "REVIEWED",
        "organism": "Arabidopsis thaliana",
        "taxon_id": "3702",
        "source_title": "Source",
        "source_doi": "",
        "source_table": "",
        "source_version": "1",
        "enabled": "yes",
        "note": "Evaluate empirically.",
    }

    def row(*, protein: str, marker_class: str, subclass: str) -> str:
        values = {
            **common,
            "protein_identifier": protein,
            "benchmark_class": marker_class,
            "benchmark_subclass": subclass,
        }
        return "\t".join(values[field] for field in BENCHMARK_FIELDS)

    path = tmp_path / "conflicting.tsv"
    path.write_text(
        f"{heading}\n"
        f"{row(protein='P1', marker_class='HOUSEKEEPING', subclass='REFERENCE')}\n"
        f"{row(protein='P2', marker_class='R_NLR', subclass='TNL')}\n",
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="conflicting classes"):
        read_benchmark_authority(path=path)

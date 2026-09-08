"""Tests for arbitrary-dataset NCBI taxdump candidate generation."""

from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

from orthofinder_interrogation_app import taxonomy_cli, taxonomy_reference
from orthofinder_interrogation_app.taxonomy import (
    TAXONOMY_COLUMNS,
    parse_taxonomy_mapping,
)
from orthofinder_interrogation_app.taxonomy_reference import build_taxonomy_candidates
from orthofinder_results.errors import InputValidationError


def _write_taxdump(*, root: Path) -> Path:
    """Write a tiny valid taxdump with unique, ambiguous and synonym matches."""

    root.mkdir(parents=True, exist_ok=True)
    (root / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "10\t|\tExample clade\t|\t\t|\tscientific name\t|\n"
        "101\t|\tSpecies A\t|\t\t|\tscientific name\t|\n"
        "101\t|\tSpecies alpha alias\t|\t\t|\tsynonym\t|\n"
        "102\t|\tSpecies B accepted\t|\t\t|\tscientific name\t|\n"
        "102\t|\tSpecies B\t|\t\t|\tsynonym\t|\n"
        "103\t|\tSpecies C\t|\t\t|\tscientific name\t|\n"
        "104\t|\tAlternative C\t|\t\t|\tscientific name\t|\n"
        "104\t|\tSpecies C\t|\t\t|\tsynonym\t|\n",
        encoding="utf-8",
    )
    (root / "nodes.dmp").write_text(
        "1\t|\t1\t|\n"
        "10\t|\t1\t|\n"
        "101\t|\t10\t|\n"
        "102\t|\t10\t|\n"
        "103\t|\t10\t|\n"
        "104\t|\t10\t|\n",
        encoding="utf-8",
    )
    return root


def _rows(data: bytes) -> list[dict[str, str]]:
    """Decode generated TSV bytes into records."""

    return list(csv.DictReader(io.StringIO(data.decode("utf-8")), delimiter="\t"))


def test_generic_candidates_cover_every_species_without_auto_review(
    tmp_path: Path,
) -> None:
    """Unique names remain pending while ambiguity and missing names stay unresolved."""

    taxdump = _write_taxdump(root=tmp_path / "taxdump")
    content = build_taxonomy_candidates(
        species=("Species_D", "Species_C", "Species_B", "Species_A"),
        taxdump_dir=taxdump,
        source_date="2026-09-07",
        source_version="taxdump-fixture-sha256",
    )
    rows = {row["workflow_species_label"]: row for row in _rows(content)}
    assert tuple(rows) == ("Species_A", "Species_B", "Species_C", "Species_D")
    assert rows["Species_A"]["mapping_status"] == "PENDING_REVIEW"
    assert rows["Species_A"]["mapping_method"] == "EXACT_SCIENTIFIC_NAME_CANDIDATE"
    assert rows["Species_A"]["lineage_taxon_ids"] == "1;10"
    assert rows["Species_A"]["parent_taxon_id"] == "10"
    assert rows["Species_B"]["mapping_method"] == "EXACT_NCBI_NAME_CANDIDATE"
    assert rows["Species_B"]["accepted_species_name"] == "Species B accepted"
    assert rows["Species_C"]["mapping_status"] == "AMBIGUOUS"
    assert "103:Species C" in rows["Species_C"]["review_note"]
    assert "104:Alternative C" in rows["Species_C"]["review_note"]
    assert rows["Species_D"]["mapping_status"] == "UNMAPPED"
    authority = parse_taxonomy_mapping(
        data=content,
        expected_species=("Species_A", "Species_B", "Species_C", "Species_D"),
    )
    assert authority.summary() == {
        "REVIEWED": 0,
        "PENDING_REVIEW": 2,
        "UNMAPPED": 1,
        "AMBIGUOUS": 1,
        "MISSING": 0,
    }
    assert authority.target_species(taxon_id=10) == ()


def test_taxonomy_candidate_cli_writes_only_external_persistent_sidecar(
    application_resource: Path,
    persistent_test_root: Path,
    tmp_path: Path,
) -> None:
    """The CLI creates a complete sidecar and refuses replacement by default."""

    taxdump = _write_taxdump(root=tmp_path / "taxdump")
    output = persistent_test_root / "taxonomy/review_candidates.tsv"
    arguments = [
        "--resource-dir",
        str(application_resource),
        "--taxdump-dir",
        str(taxdump),
        "--source-date",
        "2026-09-07",
        "--source-version",
        "fixture-v1",
        "--output-tsv",
        str(output),
        "--verbose",
    ]
    assert taxonomy_cli.main(arguments) == 0
    assert output.is_file()
    assert len(_rows(output.read_bytes())) == 4
    with pytest.raises(SystemExit) as error:
        taxonomy_cli.main(arguments)
    assert error.value.code == 2
    assert taxonomy_cli.main([*arguments, "--force"]) == 0
    inside = application_resource / "review.tsv"
    assert taxonomy_cli.main(
        [
            *arguments[:-3],
            "--output-tsv",
            str(inside),
        ]
    ) == 2


@pytest.mark.parametrize(
    ("species", "source_date", "source_version", "message"),
    [
        ((), "2026-09-07", "v1", "unique and non-empty"),
        (("Species_A", "Species_A"), "2026-09-07", "v1", "unique and non-empty"),
        ((" ",), "2026-09-07", "v1", "unique and non-empty"),
        (("Species_A",), "07-09-2026", "v1", "ISO"),
        (("Species_A",), "2026-09-07", " ", "must not be empty"),
    ],
)
def test_candidate_metadata_validation(
    tmp_path: Path,
    species: tuple[str, ...],
    source_date: str,
    source_version: str,
    message: str,
) -> None:
    """Species identities and taxdump provenance must be explicit and valid."""

    taxdump = _write_taxdump(root=tmp_path / "taxdump")
    with pytest.raises(InputValidationError, match=message):
        build_taxonomy_candidates(
            species=species,
            taxdump_dir=taxdump,
            source_date=source_date,
            source_version=source_version,
        )


@pytest.mark.parametrize(
    ("filename", "content", "message"),
    [
        ("names.dmp", "", "non-empty names.dmp"),
        ("names.dmp", "bad\n", "expected four fields"),
        ("names.dmp", "bad\t|\tName\t|\t\t|\tscientific name\t|\n", "taxon ID"),
        ("names.dmp", "0\t|\tName\t|\t\t|\tscientific name\t|\n", "record"),
        ("nodes.dmp", "bad\n", "expected two fields"),
        ("nodes.dmp", "bad\t|\t1\t|\n", "taxon ID"),
        ("nodes.dmp", "0\t|\t1\t|\n", "Non-positive"),
        ("nodes.dmp", "1\t|\t1\t|\n1\t|\t1\t|\n", "Duplicate"),
    ],
)
def test_malformed_taxdump_records_fail_with_line_context(
    tmp_path: Path,
    filename: str,
    content: str,
    message: str,
) -> None:
    """Malformed local authorities fail deterministically before candidate publication."""

    taxdump = _write_taxdump(root=tmp_path / "taxdump")
    (taxdump / filename).write_text(content, encoding="utf-8")
    with pytest.raises(InputValidationError, match=message):
        build_taxonomy_candidates(
            species=("Species_A",),
            taxdump_dir=taxdump,
            source_date="2026-09-07",
            source_version="v1",
        )


def test_lineage_and_scientific_name_integrity_checks(tmp_path: Path) -> None:
    """Missing nodes, cycles and absent/duplicate scientific names are rejected."""

    taxdump = _write_taxdump(root=tmp_path / "missing_node")
    (taxdump / "nodes.dmp").write_text("1\t|\t1\t|\n10\t|\t1\t|\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="lacks candidate taxon 101"):
        build_taxonomy_candidates(
            species=("Species_A",),
            taxdump_dir=taxdump,
            source_date="2026-09-07",
            source_version="v1",
        )

    taxdump = _write_taxdump(root=tmp_path / "cycle")
    (taxdump / "nodes.dmp").write_text(
        "1\t|\t1\t|\n10\t|\t101\t|\n101\t|\t10\t|\n",
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="cycle"):
        build_taxonomy_candidates(
            species=("Species_A",),
            taxdump_dir=taxdump,
            source_date="2026-09-07",
            source_version="v1",
        )

    taxdump = _write_taxdump(root=tmp_path / "missing_scientific")
    text = (taxdump / "names.dmp").read_text(encoding="utf-8")
    (taxdump / "names.dmp").write_text(
        text.replace(
            "10\t|\tExample clade\t|\t\t|\tscientific name",
            "10\t|\tExample clade\t|\t\t|\tsynonym",
        ),
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="lacks scientific names"):
        build_taxonomy_candidates(
            species=("Species_A",),
            taxdump_dir=taxdump,
            source_date="2026-09-07",
            source_version="v1",
        )

    taxdump = _write_taxdump(root=tmp_path / "duplicate_scientific")
    with (taxdump / "names.dmp").open("a", encoding="utf-8") as handle:
        handle.write("101\t|\tSpecies A duplicate\t|\t\t|\tscientific name\t|\n")
    with pytest.raises(InputValidationError, match="multiple scientific names"):
        build_taxonomy_candidates(
            species=("Species_A",),
            taxdump_dir=taxdump,
            source_date="2026-09-07",
            source_version="v1",
        )


def test_taxdump_helpers_cover_root_and_empty_inputs(tmp_path: Path) -> None:
    """Low-level parsing handles root taxa, empty requirements and source labels."""

    taxdump = _write_taxdump(root=tmp_path / "taxdump")
    assert taxonomy_reference._source_species_name(label="  Species__A  ") == "Species A"
    assert taxonomy_reference._normalise_name(value="  SpEcIeS   A ") == "species a"
    assert taxonomy_reference._dmp_fields(line="1\t|\troot\t|\n") == ("1", "root")
    assert taxonomy_reference._read_scientific_names(
        path=taxdump / "names.dmp",
        required_taxon_ids=set(),
    ) == {}
    parents = taxonomy_reference._read_parent_map(path=taxdump / "nodes.dmp")
    assert taxonomy_reference._lineage_ids(taxon_id=1, parent_by_id=parents) == ()
    with pytest.raises(InputValidationError, match="missing node"):
        taxonomy_reference._lineage_ids(
            taxon_id=101,
            parent_by_id={101: 10},
        )


def test_sidecar_location_validator_and_parser_contract(
    application_resource: Path,
    persistent_test_root: Path,
) -> None:
    """CLI parser remains all-named and resource-internal output is prohibited."""

    parser = taxonomy_cli.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([])
    taxonomy_cli._validate_sidecar_location(
        output=persistent_test_root / "external.tsv",
        resource_path=application_resource,
    )
    with pytest.raises(InputValidationError, match="outside"):
        taxonomy_cli._validate_sidecar_location(
            output=application_resource / "internal.tsv",
            resource_path=application_resource,
        )
    assert TAXONOMY_COLUMNS[0] == "workflow_species_label"

"""Tests for generic checksum-bound RNA-seq evidence integration."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from orthofinder_results.errors import InputValidationError, PublicationError
from orthofinder_results.expression_analysis import (
    CONTEXT_FIELDS,
    _create_source_views,
    _header_aliases,
    _iter_aliases,
    _parquet_query,
    _remove_empty_directory,
    _validate_atlas_rows,
    _write_query,
    _write_query_parquet,
    publish_expression_evidence,
    read_expression_manifest,
)
from orthofinder_results.io_utils import read_tsv, sha256_file, write_tsv
from orthofinder_results.parsers import MEMBERSHIP_FIELDS, SEQUENCE_FIELDS
from orthofinder_results.pipeline import run_pipeline


def _write_manifest(*, root: Path) -> Path:
    """Create one expression and metadata partition plus its verified manifest."""

    expression = root / "expression.parquet"
    metadata = root / "metadata.parquet"
    expression_rows = []
    for unit, values in (("TPM", (3.0, 0.1)), ("FPKM", (30.0, 10.0))):
        for context, value in zip(("leaf", "root"), values, strict=True):
            expression_rows.append(
                {
                    "source_database": "Expression Atlas",
                    "experiment_accession": "E-TEST-1",
                    "species_column": "Species_A",
                    "gene_id": "AT1G00010",
                    "gene_name": "GENEA",
                    "sample_or_condition": context,
                    "expression_value": value,
                    "expression_minimum": None,
                    "expression_lower_quartile": None,
                    "expression_median": None,
                    "expression_upper_quartile": None,
                    "expression_maximum": None,
                    "expression_value_statistic": "value",
                    "expression_summary_type": "single_value",
                    "expression_unit": unit,
                    "source_file": "atlas.tsv",
                    "source_file_sha256": "a" * 64,
                }
            )
    for gene_id in ("AT2G00020", "AT2G00021"):
        expression_rows.append(
            {
                **expression_rows[0],
                "gene_id": gene_id,
                "gene_name": "SHARED",
                "sample_or_condition": "seedling",
                "expression_value": 1.0,
                "expression_unit": "TPM",
            }
        )
    pq.write_table(pa.Table.from_pylist(expression_rows), expression)
    metadata_rows = [
        {
            "experiment_accession": "E-TEST-1",
            "species_column": "Species_A",
            "sample_or_condition": context,
            "atlas_group_label": context,
            "assay_ids": f"assay_{context}",
            "assay_count": 1,
            "organism_part": context,
            "developmental_stage": "adult",
            "genotype": "wild type",
            "cultivar": "",
            "treatment": "none",
            "condition": "control",
            "source_file": "metadata.tsv",
            "source_file_sha256": "b" * 64,
            "expression_file_sha256": "a" * 64,
        }
        for context in ("leaf", "root", "seedling")
    ]
    pq.write_table(pa.Table.from_pylist(metadata_rows), metadata)
    manifest = root / "e3_workflow_expression_resources.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "resource_id",
                "resource_type",
                "species_column",
                "dataset",
                "path",
                "sha256",
                "include",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for resource_id, resource_type, path in (
            ("expression:a", "atlas_expression_long", expression),
            ("metadata:a", "atlas_sample_metadata_wide", metadata),
        ):
            writer.writerow(
                {
                    "resource_id": resource_id,
                    "resource_type": resource_type,
                    "species_column": "Species_A",
                    "dataset": path.stem,
                    "path": path,
                    "sha256": sha256_file(path=path),
                    "include": "true",
                }
            )
    return manifest


def _rewrite_manifest(*, path: Path, rows: list[dict[str, str]]) -> None:
    """Rewrite a synthetic resource manifest after one deliberate mutation."""

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "resource_id",
                "resource_type",
                "species_column",
                "dataset",
                "path",
                "sha256",
                "include",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _manifest_rows(*, path: Path) -> list[dict[str, str]]:
    """Return mutable rows from one synthetic expression manifest."""

    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_orthofinder_tables(*, root: Path) -> Path:
    """Create the four pipeline tables needed by expression publication."""

    tables = root / "tables"
    tables.mkdir()
    sequences = (
        {
            "run_id": "run",
            "internal_id": "0_0",
            "species_index": "0",
            "species_label": "Species_A",
            "source_fasta": "Species_A.fa",
            "raw_header": "sp|P11111|PROTA_SPECIES GN=GENEA",
            "member_id": "sp|P11111|PROTA_SPECIES",
            "source_file": "SequenceIDs.txt",
            "source_line": 1,
        },
        {
            "run_id": "run",
            "internal_id": "0_1",
            "species_index": "0",
            "species_label": "Species_A",
            "source_fasta": "Species_A.fa",
            "raw_header": "protein_two GN=SHARED",
            "member_id": "protein_two",
            "source_file": "SequenceIDs.txt",
            "source_line": 2,
        },
        {
            "run_id": "run",
            "internal_id": "1_0",
            "species_index": "1",
            "species_label": "Species_B",
            "source_fasta": "Species_B.fa",
            "raw_header": "protein_three",
            "member_id": "protein_three",
            "source_file": "SequenceIDs.txt",
            "source_line": 3,
        },
    )
    write_tsv(
        path=tables / "sequences.tsv.gz",
        fieldnames=SEQUENCE_FIELDS,
        records=sequences,
    )
    memberships = tuple(
        {
            "run_id": "run",
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": "N0.HOG1",
            "legacy_orthogroup_id": "OG1",
            "gene_tree_parent_clade": "N0",
            "species_label": row["species_label"],
            "member_id": row["member_id"],
            "source_file": "N0.tsv",
            "source_row": 2,
        }
        for row in sequences
    )
    write_tsv(
        path=tables / "hog_memberships.tsv.gz",
        fieldnames=MEMBERSHIP_FIELDS,
        records=memberships,
    )
    write_tsv(
        path=tables / "legacy_orthogroup_memberships.tsv.gz",
        fieldnames=MEMBERSHIP_FIELDS,
        records=(),
    )
    return tables


def test_expression_manifest_is_checksum_bound(tmp_path: Path) -> None:
    """A changed included partition is rejected before analytical work."""

    manifest = _write_manifest(root=tmp_path)
    authority = read_expression_manifest(path=manifest)
    assert len(authority.records) == 2
    authority.records[0].path.write_bytes(b"changed")
    with pytest.raises(InputValidationError, match="checksum differs"):
        read_expression_manifest(path=manifest)


def test_expression_manifest_rejects_incomplete_or_unsafe_authorities(
    tmp_path: Path,
) -> None:
    """Every manifest identity, schema and inclusion decision is validated."""

    with pytest.raises(InputValidationError, match="missing or empty"):
        read_expression_manifest(path=tmp_path / "missing.tsv")
    empty = tmp_path / "empty.tsv"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(InputValidationError, match="missing or empty"):
        read_expression_manifest(path=empty)
    headerless = tmp_path / "headerless.tsv"
    headerless.write_text("\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="lacks fields"):
        read_expression_manifest(path=headerless)
    incomplete = tmp_path / "incomplete.tsv"
    incomplete.write_text("resource_id\tinclude\nonly\ttrue\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="lacks fields"):
        read_expression_manifest(path=incomplete)

    cases = (
        ("include", "perhaps", "Invalid include"),
        ("resource_id", "", "Empty or duplicate"),
        ("resource_type", "unsafe", "Unsupported expression resource_type"),
        ("species_column", "", "lacks species_column"),
        ("path", str(tmp_path / "missing.parquet"), "resource is missing or empty"),
        ("sha256", "invalid", "invalid SHA-256"),
    )
    for index, (field, value, message) in enumerate(cases):
        root = tmp_path / f"case_{index}"
        root.mkdir()
        manifest = _write_manifest(root=root)
        rows = _manifest_rows(path=manifest)
        rows[0][field] = value
        _rewrite_manifest(path=manifest, rows=rows)
        with pytest.raises(InputValidationError, match=message):
            read_expression_manifest(path=manifest, verify_checksums=False)

    duplicate_root = tmp_path / "duplicate"
    duplicate_root.mkdir()
    duplicate_manifest = _write_manifest(root=duplicate_root)
    duplicate_rows = _manifest_rows(path=duplicate_manifest)
    duplicate_rows[1]["resource_id"] = duplicate_rows[0]["resource_id"]
    _rewrite_manifest(path=duplicate_manifest, rows=duplicate_rows)
    with pytest.raises(InputValidationError, match="Empty or duplicate"):
        read_expression_manifest(path=duplicate_manifest, verify_checksums=False)

    excluded_root = tmp_path / "excluded"
    excluded_root.mkdir()
    excluded_manifest = _write_manifest(root=excluded_root)
    excluded_rows = _manifest_rows(path=excluded_manifest)
    excluded_rows[0]["include"] = "false"
    _rewrite_manifest(path=excluded_manifest, rows=excluded_rows)
    with pytest.raises(InputValidationError, match="no included atlas_expression_long"):
        read_expression_manifest(path=excluded_manifest, verify_checksums=False)


def test_expression_controls_aliases_and_queries_fail_closed(tmp_path: Path) -> None:
    """Bounds, reviewed aliases and controlled Parquet scans reject unsafe input."""

    with pytest.raises(InputValidationError, match="At least one"):
        _parquet_query(paths=())
    aliases = _header_aliases(
        member_id="sp|P12345.2|ENTRY_ARATH",
        raw_header="sp|P12345.2|ENTRY_ARATH GN=AT1G01010.1 gene=EXAMPLE",
    )
    assert ("uniprot_accession", "P12345.2", 1, "SequenceIDs.txt") in aliases
    assert any(alias_type.endswith("_versionless") for alias_type, *_rest in aliases)

    source = tmp_path / "source"
    source.mkdir()
    authority = read_expression_manifest(path=_write_manifest(root=source))
    empty_tables = tmp_path / "empty_tables"
    empty_tables.mkdir()
    write_tsv(
        path=empty_tables / "sequences.tsv.gz",
        fieldnames=SEQUENCE_FIELDS,
        records=(),
    )
    with pytest.raises(InputValidationError, match="No exact member identifiers"):
        publish_expression_evidence(
            tables_dir=empty_tables,
            work_dir=tmp_path / "empty_work",
            run_id="run",
            authority=authority,
            threads=1,
            memory_limit_mb=256,
        )
    tables = _write_orthofinder_tables(root=tmp_path)
    invalid_controls = (
        ({"minimum_expression_value": float("nan")}, "finite and non-negative"),
        ({"broad_positive_fraction": 2.0}, "between zero and one"),
        ({"threads": 0}, "between 1 and 256"),
        ({"memory_limit_mb": 255}, "between 256"),
    )
    for overrides, message in invalid_controls:
        arguments = {
            "tables_dir": tables,
            "work_dir": tmp_path / "work",
            "run_id": "run",
            "authority": authority,
        }
        arguments.update(overrides)
        with pytest.raises(InputValidationError, match=message):
            publish_expression_evidence(**arguments)

    bad_aliases = tmp_path / "bad_aliases.tsv"
    write_tsv(
        path=bad_aliases,
        fieldnames=("species_label", "member_id"),
        records=({"species_label": "Species_A", "member_id": "protein_two"},),
    )
    with pytest.raises(InputValidationError, match="alias TSV lacks fields"):
        publish_expression_evidence(
            tables_dir=tables,
            work_dir=tmp_path / "work_aliases",
            run_id="run",
            authority=authority,
            additional_aliases_path=bad_aliases,
            threads=1,
            memory_limit_mb=256,
        )

    invalid_tier = tmp_path / "invalid_tier.tsv"
    write_tsv(
        path=invalid_tier,
        fieldnames=(
            "species_label",
            "member_id",
            "alias_type",
            "alias_value",
            "mapping_tier",
            "source",
        ),
        records=(
            {
                "species_label": "Species_A",
                "member_id": "protein_two",
                "alias_type": "gene",
                "alias_value": "AT2G00020",
                "mapping_tier": "not-an-integer",
                "source": "test",
            },
        ),
    )
    with pytest.raises(InputValidationError, match="mapping_tier is not an integer"):
        publish_expression_evidence(
            tables_dir=tables,
            work_dir=tmp_path / "work_tier",
            run_id="run",
            authority=authority,
            additional_aliases_path=invalid_tier,
            threads=1,
            memory_limit_mb=256,
        )


def test_expression_source_views_and_rows_validate_before_mapping(tmp_path: Path) -> None:
    """Atlas schemas, units, values and metadata keys fail before member mapping."""

    connection = duckdb.connect()
    try:
        with pytest.raises(InputValidationError, match="no expression partitions"):
            _create_source_views(
                connection=connection,
                expression_paths=(),
                all_expression_paths=(),
                metadata_paths=(),
            )
    finally:
        connection.close()

    bad_expression = tmp_path / "bad_expression.parquet"
    pq.write_table(pa.table({"wrong": ["value"]}), bad_expression)
    connection = duckdb.connect()
    try:
        with pytest.raises(InputValidationError, match="Expression Parquet lacks fields"):
            _create_source_views(
                connection=connection,
                expression_paths=(bad_expression,),
                all_expression_paths=(bad_expression,),
                metadata_paths=(),
            )
    finally:
        connection.close()

    source = tmp_path / "valid_source"
    source.mkdir()
    authority = read_expression_manifest(path=_write_manifest(root=source))
    expression_path = next(
        record.path
        for record in authority.records
        if record.resource_type == "atlas_expression_long"
    )
    bad_metadata = tmp_path / "bad_metadata.parquet"
    pq.write_table(pa.table({"wrong": ["value"]}), bad_metadata)
    connection = duckdb.connect()
    try:
        with pytest.raises(InputValidationError, match="metadata Parquet lacks fields"):
            _create_source_views(
                connection=connection,
                expression_paths=(expression_path,),
                all_expression_paths=(expression_path,),
                metadata_paths=(bad_metadata,),
            )
    finally:
        connection.close()

    connection = duckdb.connect()
    try:
        _create_source_views(
            connection=connection,
            expression_paths=(),
            all_expression_paths=(expression_path,),
            metadata_paths=(),
        )
        _validate_atlas_rows(connection=connection)
    finally:
        connection.close()

    for unit, value, duplicate_metadata, message in (
        ("COUNTS", 1.0, False, "outside TPM/FPKM"),
        ("TPM", -1.0, False, "negative or non-finite"),
        ("TPM", 1.0, True, "duplicated context keys"),
    ):
        connection = duckdb.connect()
        try:
            connection.execute(
                "CREATE TABLE atlas_expression(expression_unit VARCHAR, "
                "expression_value DOUBLE)"
            )
            connection.execute("INSERT INTO atlas_expression VALUES (?, ?)", [unit, value])
            connection.execute(
                "CREATE TABLE atlas_metadata(species_column VARCHAR, "
                "experiment_accession VARCHAR, sample_or_condition VARCHAR)"
            )
            connection.execute("INSERT INTO atlas_metadata VALUES ('Species_A','E-1','leaf')")
            if duplicate_metadata:
                connection.execute(
                    "INSERT INTO atlas_metadata VALUES ('Species_A','E-1','leaf')"
                )
            with pytest.raises(InputValidationError, match=message):
                _validate_atlas_rows(connection=connection)
        finally:
            connection.close()


def test_expression_alias_and_output_edge_cases_are_explicit(tmp_path: Path) -> None:
    """Empty aliases, NUL aliases and output-column drift cannot silently publish."""

    sequence_path = tmp_path / "sequences.tsv.gz"
    write_tsv(
        path=sequence_path,
        fieldnames=SEQUENCE_FIELDS,
        records=(
            {
                "run_id": "run",
                "internal_id": "0_0",
                "species_index": "0",
                "species_label": "Species_A",
                "source_fasta": "Species_A.fa",
                "raw_header": "raw_token description",
                "member_id": "member",
                "source_file": "SequenceIDs.txt",
                "source_line": 1,
            },
        ),
    )
    empty_aliases = tmp_path / "empty_aliases.tsv.gz"
    write_tsv(
        path=empty_aliases,
        fieldnames=(
            "species_label",
            "member_id",
            "alias_type",
            "alias_value",
            "mapping_tier",
            "source",
        ),
        records=(),
    )
    assert tuple(
        _iter_aliases(
            sequence_path=sequence_path,
            run_id="run",
            additional_aliases_path=empty_aliases,
        )
    )
    assert ("raw_primary_token", "raw_token", 1, "SequenceIDs.txt") in _header_aliases(
        member_id="member",
        raw_header="raw_token description",
    )

    nul_aliases = tmp_path / "nul_aliases.tsv.gz"
    write_tsv(
        path=nul_aliases,
        fieldnames=(
            "species_label",
            "member_id",
            "alias_type",
            "alias_value",
            "mapping_tier",
            "source",
        ),
        records=(
            {
                "species_label": "Species_A",
                "member_id": "member",
                "alias_type": "unsafe",
                "alias_value": "bad\x00alias",
                "mapping_tier": 1,
                "source": "test",
            },
            {
                "species_label": "Species_A",
                "member_id": "member",
                "alias_type": "aaa_preferred",
                "alias_value": "member",
                "mapping_tier": 1,
                "source": "reviewed",
            },
        ),
    )
    aliases = tuple(
        _iter_aliases(
            sequence_path=sequence_path,
            run_id="run",
            additional_aliases_path=nul_aliases,
        )
    )
    assert all("\x00" not in row["alias_value"] for row in aliases)
    member_alias = next(row for row in aliases if row["alias_value"] == "member")
    assert member_alias["alias_type"] == "aaa_preferred"

    invalid_aliases = tmp_path / "invalid_aliases.tsv.gz"
    write_tsv(
        path=invalid_aliases,
        fieldnames=(
            "species_label",
            "member_id",
            "alias_type",
            "alias_value",
            "mapping_tier",
            "source",
        ),
        records=(
            {
                "species_label": "Species_A",
                "member_id": "member",
                "alias_type": "reviewed",
                "alias_value": "",
                "mapping_tier": 101,
                "source": "test",
            },
        ),
    )
    with pytest.raises(InputValidationError, match="mapping_tier between 1 and 100"):
        tuple(
            _iter_aliases(
                sequence_path=sequence_path,
                run_id="run",
                additional_aliases_path=invalid_aliases,
            )
        )

    connection = duckdb.connect()
    try:
        connection.execute("CREATE TABLE one_column(value VARCHAR)")
        with pytest.raises(PublicationError, match="query columns differ"):
            _write_query(
                connection=connection,
                query="SELECT value FROM one_column",
                path=tmp_path / "wrong.tsv.gz",
                fieldnames=CONTEXT_FIELDS,
            )
        embedded_text = 'annotation with\ttab, "quote" and\nnewline'
        connection.execute("INSERT INTO one_column VALUES (?)", [embedded_text])
        context_parquet = tmp_path / "context.parquet"
        assert (
            _write_query_parquet(
                connection=connection,
                query="SELECT value FROM one_column",
                path=context_parquet,
                fieldnames=("value",),
            )
            == 1
        )
        assert (
            connection.execute(
                "SELECT value FROM read_parquet(?)", [str(context_parquet)]
            ).fetchone()[0]
            == embedded_text
        )
        assert (
            _write_query_parquet(
                connection=connection,
                query="SELECT value FROM one_column",
                path=context_parquet,
                fieldnames=("value",),
            )
            == 1
        )
        with pytest.raises(PublicationError, match="query columns differ"):
            _write_query_parquet(
                connection=connection,
                query="SELECT value FROM one_column",
                path=tmp_path / "wrong.parquet",
                fieldnames=("different",),
            )
    finally:
        connection.close()

    nonempty = tmp_path / "nonempty_spill"
    nonempty.mkdir()
    (nonempty / "spill.tmp").write_text("retained", encoding="utf-8")
    _remove_empty_directory(path=nonempty)
    assert nonempty.is_dir()


def test_expression_publication_keeps_mapping_states_and_preferred_units(
    tmp_path: Path,
) -> None:
    """Unique, ambiguous and unavailable mappings remain explicit and TPM wins."""

    source = tmp_path / "source"
    source.mkdir()
    manifest = _write_manifest(root=source)
    tables = _write_orthofinder_tables(root=tmp_path)
    publication = publish_expression_evidence(
        tables_dir=tables,
        work_dir=tmp_path / "work",
        run_id="run",
        authority=read_expression_manifest(path=manifest),
        minimum_expression_value=0.5,
        broad_positive_fraction=0.5,
        threads=1,
        memory_limit_mb=256,
    )
    assert publication.counts["expression_member_count"] == 3
    mappings = tuple(read_tsv(path=tables / "expression_member_mapping.tsv.gz"))
    statuses = {row["member_id"]: row["mapping_status"] for row in mappings}
    assert statuses == {
        "protein_three": "NOT_MAPPED",
        "protein_two": "AMBIGUOUS",
        "sp|P11111|PROTA_SPECIES": "MAPPED_UNIQUE",
    }
    context_connection = duckdb.connect()
    try:
        context_rows = context_connection.execute(
            "SELECT expression_unit, expression_context "
            "FROM read_parquet(?)",
            [str(tables / "expression_context.parquet")],
        ).fetchall()
    finally:
        context_connection.close()
    contexts = tuple(
        {"expression_unit": row[0], "expression_context": row[1]}
        for row in context_rows
    )
    assert {row["expression_unit"] for row in contexts} == {"TPM"}
    assert {row["expression_context"] for row in contexts} == {"leaf", "root"}
    summaries = tuple(read_tsv(path=tables / "expression_member_summary.tsv.gz"))
    mapped = next(row for row in summaries if row["mapping_status"] == "MAPPED_UNIQUE")
    assert mapped["positive_context_count"] == "1"
    assert mapped["broad_expression_supported"] == "true"
    groups = tuple(read_tsv(path=tables / "expression_group_summary.tsv.gz"))
    assert groups[0]["unique_mapped_member_count"] == "1"
    assert groups[0]["ambiguous_member_count"] == "1"


def test_pipeline_publishes_schema5_expression_and_sequence_resource(
    orthofinder2_results: Path,
    persistent_test_root: Path,
    tmp_path: Path,
) -> None:
    """A rebuilt resource embeds sequences, expression tables and provenance."""

    working = orthofinder2_results / "WorkingDirectory"
    (working / "Species0.fa").write_text(
        ">0_0\nMAAN\n>0_1\nMQQN\n",
        encoding="utf-8",
    )
    (working / "Species1.fa").write_text(">1_0\nMTTN\n", encoding="utf-8")
    manifest = _write_manifest(root=tmp_path)
    aliases = tmp_path / "expression_aliases.tsv"
    write_tsv(
        path=aliases,
        fieldnames=(
            "species_label",
            "member_id",
            "alias_type",
            "alias_value",
            "mapping_tier",
            "source",
        ),
        records=(
            {
                "species_label": "Species_A",
                "member_id": "protA",
                "alias_type": "reviewed_gene_id",
                "alias_value": "AT1G00010",
                "mapping_tier": 1,
                "source": "synthetic reviewed authority",
            },
        ),
    )
    output = persistent_test_root / "expression_resource"
    manifest_record = run_pipeline(
        results_dir=orthofinder2_results,
        output_dir=output,
        run_id="expression-run",
        work_dir=persistent_test_root / "expression_work",
        alignment_dir=None,
        distance_source="NONE",
        distance_group_type="AUTO",
        distance_hierarchy_node="N0",
        distance_max_groups=0,
        distance_max_members=10,
        parse_gene_trees=False,
        report_max_statistic_rows=100,
        report_max_groups=10,
        report_max_members=10,
        report_nearest_neighbours=2,
        resume=False,
        force=False,
        verbose=False,
        expression_manifest_path=manifest,
        expression_aliases_path=aliases,
        expression_threads=1,
        expression_memory_mb=256,
        include_protein_sequences=True,
    )
    assert manifest_record["schema_version"] == 5
    assert manifest_record["counts"]["protein_sequence_count"] == 3
    assert manifest_record["counts"]["expression_unique_mapping_count"] == 1
    assert not (output / "evidence/protein_sequences.parquet").exists()
    assert not (output / "tables/expression_context.tsv.gz").exists()
    assert not (output / "tables/expression_context.parquet").exists()
    published = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    assert published["rna_seq_expression"]["unit_selection_policy"].startswith("TPM")
    connection = duckdb.connect(
        str(output / "duckdb/orthofinder_results.duckdb"), read_only=True
    )
    try:
        relations = {
            row[0]
            for row in connection.execute(
                "SELECT table_name FROM information_schema.tables"
            ).fetchall()
        }
        assert {
            "expression_member_mapping",
            "expression_member_summary",
            "expression_context",
            "expression_group_summary",
            "expression_import_audit",
            "protein_sequences",
        }.issubset(relations)
        assert connection.execute("SELECT count(*) FROM protein_sequences").fetchone()[0] == 3
    finally:
        connection.close()

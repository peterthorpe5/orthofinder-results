"""Tests for complete-proteome C-terminal motif discovery."""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from orthofinder_interrogation_app.resource import open_resource
from orthofinder_interrogation_app.terminal_motif import (
    motif_group_members,
    motif_group_summary,
    motif_species,
    validate_motif,
    validate_sequence_sidecar,
)
from orthofinder_interrogation_app.terminal_motif_cli import build_sequence_sidecar, main
from orthofinder_results.errors import InputValidationError


def _write_sidecar(*, path: Path) -> Path:
    """Write sequence rows matching the shared application resource."""

    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "internal_id": "0_0", "species_label": "Species_A",
                    "member_id": "alpha_1", "sequence": "MAAN",
                },
                {
                    "internal_id": "0_1", "species_label": "Species_A",
                    "member_id": "alpha_2", "sequence": "MQQN",
                },
                {
                    "internal_id": "1_0", "species_label": "Species_B",
                    "member_id": "beta_1", "sequence": "MTTN",
                },
            ]
        ),
        path,
        compression="zstd",
    )
    return path


def test_motif_validation_and_sidecar_schema(tmp_path: Path) -> None:
    """Motifs and external sequence contracts fail closed."""

    assert validate_motif(motif=" n ") == "N"
    assert validate_motif(motif="gg") == "GG"
    for invalid in ("", "NX", "N-", "N" * 101):
        with pytest.raises(InputValidationError):
            validate_motif(motif=invalid)
    with pytest.raises(InputValidationError, match="must be text"):
        validate_motif(motif=1)  # type: ignore[arg-type]
    with pytest.raises(InputValidationError, match="missing or empty"):
        validate_sequence_sidecar(path=tmp_path / "missing.parquet")
    unreadable = tmp_path / "unreadable.parquet"
    unreadable.write_text("not parquet", encoding="utf-8")
    with pytest.raises(InputValidationError, match="not readable Parquet"):
        validate_sequence_sidecar(path=unreadable)
    bad = tmp_path / "bad.parquet"
    pq.write_table(pa.table({"wrong": ["value"]}), bad)
    with pytest.raises(InputValidationError, match="lacks required"):
        validate_sequence_sidecar(path=bad)


def test_motif_queries_return_exact_hog_and_members(
    application_resource: Path, tmp_path: Path
) -> None:
    """Protein fractions, species breadth and exact member calls agree."""

    sidecar = _write_sidecar(path=tmp_path / "sequences.parquet")
    resource = open_resource(path=application_resource)
    assert motif_species(sidecar_path=sidecar) == ("Species_A", "Species_B")
    rows = motif_group_summary(
        resource=resource,
        sidecar_path=sidecar,
        motif="N",
        threshold=0.8,
        minimum_species=2,
        required_species=("Species_A", "Species_B"),
        hierarchy_node="N0",
    )
    assert rows[0]["group_id"] == "N0.HOG1"
    assert rows[0]["matching_fraction"] == 1.0
    members = motif_group_members(
        resource=resource,
        sidecar_path=sidecar,
        motif="N",
        group_id="N0.HOG1",
    )
    assert len(members) == 3
    assert all(row["motif_match"] for row in members)
    assert {row["observed_terminus"] for row in members} == {"N"}
    for kwargs, message in (
        ({"threshold": 2.0}, "threshold"),
        ({"minimum_species": 0}, "Minimum species"),
        ({"maximum_rows": 0}, "Maximum result rows"),
    ):
        with pytest.raises(InputValidationError, match=message):
            motif_group_summary(
                resource=resource,
                sidecar_path=sidecar,
                motif="N",
                **kwargs,
            )


def test_sequence_sidecar_builder_and_cli(orthofinder2_results: Path, tmp_path: Path) -> None:
    """The builder reconciles complete internal FASTA identifiers atomically."""

    working = orthofinder2_results / "WorkingDirectory"
    (working / "Species0.fa").write_text(
        ">0_0\nMAAN\n>0_1\nMQQN\n", encoding="utf-8"
    )
    (working / "Species1.fa").write_text(">1_0\nMTTA\n", encoding="utf-8")
    output = tmp_path / "terminal.parquet"
    assert build_sequence_sidecar(
        results_dir=orthofinder2_results, output_path=output
    ) == 3
    assert motif_species(sidecar_path=output) == ("Species_A", "Species_B")
    assert main(
        [
            "--orthofinder-results-dir", str(orthofinder2_results),
            "--output-parquet", str(tmp_path / "cli.parquet"),
        ]
    ) == 0
    with pytest.raises(InputValidationError, match="must end in .parquet"):
        build_sequence_sidecar(
            results_dir=orthofinder2_results,
            output_path=tmp_path / "wrong.tsv",
        )


def test_builder_rejects_partial_reconciliation(
    orthofinder2_results: Path, tmp_path: Path
) -> None:
    """Missing source sequences cannot produce a misleading partial sidecar."""

    working = orthofinder2_results / "WorkingDirectory"
    (working / "Species0.fa").write_text(">0_0\nMAAN\n", encoding="utf-8")
    (working / "Species1.fa").write_text(">1_0\nMTTA\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="reconciliation failed"):
        build_sequence_sidecar(
            results_dir=orthofinder2_results,
            output_path=tmp_path / "partial.parquet",
        )
    assert main(
        [
            "--orthofinder-results-dir", str(tmp_path / "missing_results"),
            "--output-parquet", str(tmp_path / "never.parquet"),
        ]
    ) == 2


def test_builder_accepts_original_fasta_names(
    orthofinder2_results: Path, tmp_path: Path
) -> None:
    """Original FASTA names are a supported fallback to Species-index files."""

    working = orthofinder2_results / "WorkingDirectory"
    (working / "Species_A.fa").write_text(
        ">0_0\nMAAN\n>0_1\nMQQN\n", encoding="utf-8"
    )
    (working / "Species_B.faa").write_text(">1_0\nMTTA\n", encoding="utf-8")
    output = tmp_path / "fallback.parquet"
    assert build_sequence_sidecar(
        results_dir=working,
        output_path=output,
    ) == 3


def test_builder_rejects_duplicate_internal_identifiers(
    orthofinder2_results: Path, tmp_path: Path
) -> None:
    """Duplicate SequenceIDs cannot be silently collapsed."""

    working = orthofinder2_results / "WorkingDirectory"
    (working / "SequenceIDs.txt").write_text(
        "0_0: protA\n0_0: protA_duplicate\n1_0: protB\n",
        encoding="utf-8",
    )
    with pytest.raises(InputValidationError, match="duplicate internal"):
        build_sequence_sidecar(
            results_dir=orthofinder2_results,
            output_path=tmp_path / "duplicates.parquet",
        )

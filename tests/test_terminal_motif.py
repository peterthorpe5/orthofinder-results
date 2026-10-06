"""Tests for complete-proteome C-terminal motif discovery."""

from __future__ import annotations

import shutil
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from orthofinder_interrogation_app import terminal_motif_page
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_interrogation_app.taxonomy import TaxonomyAuthority, TaxonomyRecord
from orthofinder_interrogation_app.terminal_motif import (
    EXACT_C_TERMINAL,
    REGEX_ANYWHERE,
    REGEX_C_TERMINAL,
    motif_group_members,
    motif_group_summary,
    motif_species,
    sequence_authority_available,
    validate_motif,
    validate_sequence_search,
    validate_sequence_sidecar,
)
from orthofinder_interrogation_app.terminal_motif_cli import build_sequence_sidecar, main
from orthofinder_interrogation_app.terminal_motif_page import (
    _apply_taxonomy_filters,
    _display_member,
    _display_summary,
    _filter_member_rows,
    _load_taxonomy,
    _members_to_fasta,
    _motif_figure,
    _protein_description,
    _render_group_controls,
    _render_search_controls,
    _search_stem,
    _species_distribution,
    _species_evidence_figure,
    _species_roles,
    _taxonomic_heatmap,
)
from orthofinder_results.errors import InputValidationError


def _write_sidecar(*, path: Path) -> Path:
    """Write sequence rows matching the shared application resource."""

    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "internal_id": "0_0",
                    "species_label": "Species_A",
                    "member_id": "alpha_1",
                    "sequence": "MAAN",
                },
                {
                    "internal_id": "0_1",
                    "species_label": "Species_A",
                    "member_id": "alpha_2",
                    "sequence": "MQQN",
                },
                {
                    "internal_id": "1_0",
                    "species_label": "Species_B",
                    "member_id": "beta_1",
                    "sequence": "MTTN",
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


def test_sequence_search_requires_explicit_valid_regex() -> None:
    """Exact suffixes and the two opt-in regex locations remain distinct."""

    exact = validate_sequence_search(expression=" n ", mode=EXACT_C_TERMINAL)
    assert exact.expression == "N"
    assert exact.query_expression == "N"
    anywhere = validate_sequence_search(expression=r"N[^P][ST]", mode=REGEX_ANYWHERE)
    assert anywhere.query_expression == r"N[^P][ST]"
    terminal = validate_sequence_search(expression=r"N.{2}", mode=REGEX_C_TERMINAL)
    assert terminal.query_expression == r"(?:N.{2})$"
    anchored = validate_sequence_search(expression=r"N$", mode=REGEX_C_TERMINAL)
    assert anchored.query_expression == r"N$"
    for expression in ("", "(", "N\nP"):
        with pytest.raises(InputValidationError):
            validate_sequence_search(expression=expression, mode=REGEX_ANYWHERE)
    with pytest.raises(InputValidationError, match="must be text"):
        validate_sequence_search(expression=1, mode=REGEX_ANYWHERE)  # type: ignore[arg-type]
    with pytest.raises(InputValidationError, match="Unsupported sequence-search mode"):
        validate_sequence_search(expression="N", mode="UNSAFE")
    with pytest.raises(InputValidationError, match="500 characters"):
        validate_sequence_search(expression="N" * 501, mode=REGEX_ANYWHERE)


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
    assert rows[0]["represented_species"] == "Species_A; Species_B"
    members = motif_group_members(
        resource=resource,
        sidecar_path=sidecar,
        motif="N",
        group_id="N0.HOG1",
    )
    assert len(members) == 3
    assert all(row["motif_match"] for row in members)
    assert {row["matched_sequence"] for row in members} == {"N"}
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


def test_motif_queries_use_embedded_sequences_without_a_sidecar(
    application_resource: Path, tmp_path: Path
) -> None:
    """A standalone DuckDB retains complete motif functionality."""

    resource_without_sequences = open_resource(path=application_resource)
    assert not sequence_authority_available(resource=resource_without_sequences)
    with pytest.raises(InputValidationError, match="No embedded protein-sequence"):
        motif_species(resource=resource_without_sequences)
    with pytest.raises(InputValidationError, match="No embedded protein-sequence"):
        motif_group_summary(resource=resource_without_sequences, motif="N")

    sidecar = _write_sidecar(path=tmp_path / "sequences.parquet")
    database = tmp_path / "standalone.duckdb"
    shutil.copy2(
        application_resource / "duckdb/orthofinder_results.duckdb",
        database,
    )
    connection = duckdb.connect(str(database))
    try:
        connection.execute(
            "CREATE TABLE protein_sequences AS SELECT * FROM read_parquet(?)",
            [str(sidecar)],
        )
    finally:
        connection.close()
    resource = open_resource(path=database)
    assert sequence_authority_available(resource=resource)
    assert motif_species(resource=resource) == ("Species_A", "Species_B")
    rows = motif_group_summary(
        resource=resource,
        motif="N",
        threshold=0.8,
        minimum_species=2,
        required_species=("Species_A", "Species_B"),
        hierarchy_node="N0",
    )
    assert rows[0]["group_id"] == "N0.HOG1"
    members = motif_group_members(
        resource=resource,
        motif="N",
        group_id="N0.HOG1",
    )
    assert len(members) == 3
    assert all(row["motif_match"] for row in members)


def test_sequence_sidecar_builder_and_cli(orthofinder2_results: Path, tmp_path: Path) -> None:
    """The builder reconciles complete internal FASTA identifiers atomically."""

    working = orthofinder2_results / "WorkingDirectory"
    (working / "Species0.fa").write_text(">0_0\nMAAN\n>0_1\nMQQN\n", encoding="utf-8")
    (working / "Species1.fa").write_text(">1_0\nMTTA\n", encoding="utf-8")
    output = tmp_path / "terminal.parquet"
    assert build_sequence_sidecar(results_dir=orthofinder2_results, output_path=output) == 3
    assert motif_species(sidecar_path=output) == ("Species_A", "Species_B")
    assert (
        main(
            [
                "--orthofinder-results-dir",
                str(orthofinder2_results),
                "--output-parquet",
                str(tmp_path / "cli.parquet"),
            ]
        )
        == 0
    )
    with pytest.raises(InputValidationError, match="must end in .parquet"):
        build_sequence_sidecar(
            results_dir=orthofinder2_results,
            output_path=tmp_path / "wrong.tsv",
        )


def test_taxonomic_hog_filter_uses_reviewed_descendants() -> None:
    """Lineage coverage and required/excluded clades operate on reviewed labels."""

    def record(label: str, taxon_id: int, lineage: tuple[int, ...]) -> TaxonomyRecord:
        return TaxonomyRecord(
            workflow_species_label=label,
            source_species_name=label,
            accepted_species_name=label,
            ncbi_taxon_id=taxon_id,
            parent_taxon_id=lineage[-1],
            parent_taxon_name="parent",
            lineage_taxon_ids=lineage,
            lineage_names=tuple(f"Taxon {value}" for value in lineage),
            mapping_status="REVIEWED",
            mapping_method="manual",
            mapping_source="test",
            source_date="2026-10-03",
            source_version="test",
            reviewed_by="tester",
            reviewed_at_utc="2026-10-03T00:00:00Z",
            review_note="test",
        )

    authority = TaxonomyAuthority(
        records=(
            record("Arabidopsis", 3702, (1, 33090, 71240)),
            record("Rice", 4530, (1, 33090, 4447)),
            record("Human", 9606, (1, 33208, 9605)),
        ),
        expected_species=("Arabidopsis", "Rice", "Human"),
    )
    rows = (
        {
            "group_id": "N0.HOG1",
            "represented_species": "Arabidopsis; Human; Rice",
            "matching_species": "Arabidopsis; Human",
            "analysis_assessed_species": "Arabidopsis; Rice",
            "analysis_matching_species": "Arabidopsis",
        },
    )
    filtered = _apply_taxonomy_filters(
        rows=rows,
        authority=authority,
        minimum_lineage_fraction=0.5,
        required_taxon_ids=(9605,),
        excluded_taxon_ids=(),
    )
    assert filtered[0]["analysis_species_match_fraction"] == 0.5
    assert not _apply_taxonomy_filters(
        rows=rows,
        authority=authority,
        minimum_lineage_fraction=0.8,
        required_taxon_ids=(),
        excluded_taxon_ids=(),
    )
    assert not _apply_taxonomy_filters(
        rows=rows,
        authority=authority,
        minimum_lineage_fraction=0.0,
        required_taxon_ids=(),
        excluded_taxon_ids=(9605,),
    )
    assert (
        _apply_taxonomy_filters(
            rows=rows,
            authority=None,
            minimum_lineage_fraction=0.0,
            required_taxon_ids=(),
            excluded_taxon_ids=(),
        )
        == ({**rows[0], "analysis_species_match_fraction": 0.5},)
    )
    assert not _apply_taxonomy_filters(
        rows=rows,
        authority=authority,
        minimum_lineage_fraction=0.0,
        required_taxon_ids=(4447,),
        excluded_taxon_ids=(),
    )

    distribution = _species_distribution(
        rows=(
            {
                "species_label": "Unknown",
                "motif_match": False,
                "sequence_available": False,
            },
        ),
        authority=authority,
        roles={"Unknown": "OTHER_SAMPLED"},
    )
    assert distribution[0]["Accepted species"] == "Unresolved"


def test_taxonomic_visual_helpers_keep_all_evidence_states() -> None:
    """Visual and download helpers preserve match, no-match and missing states."""

    def record(label: str, taxon_id: int, lineage: tuple[int, ...]) -> TaxonomyRecord:
        return TaxonomyRecord(
            workflow_species_label=label,
            source_species_name=label,
            accepted_species_name=f"Accepted {label}",
            ncbi_taxon_id=taxon_id,
            parent_taxon_id=lineage[-1],
            parent_taxon_name="parent",
            lineage_taxon_ids=lineage,
            lineage_names=tuple(f"Taxon {value}" for value in lineage),
            mapping_status="REVIEWED",
            mapping_method="manual",
            mapping_source="test",
            source_date="2026-10-05",
            source_version="test",
            reviewed_by="tester",
            reviewed_at_utc="2026-10-05T00:00:00Z",
            review_note="test",
        )

    authority = TaxonomyAuthority(
        records=(
            record("Primary", 11, (1, 10)),
            record("Required", 21, (1, 20)),
            record("Excluded", 31, (1, 30)),
            record("Other", 41, (1, 40)),
        ),
        expected_species=("Primary", "Required", "Excluded", "Other"),
    )
    roles = _species_roles(
        species=authority.expected_species,
        analysis_species=("Primary",),
        authority=authority,
        required_taxon_ids=(20,),
        excluded_taxon_ids=(30,),
    )
    assert roles == {
        "Primary": "PRIMARY_ANALYSIS",
        "Required": "REQUIRED_COMPARISON",
        "Excluded": "EXCLUDED_COMPARISON",
        "Other": "OTHER_SAMPLED",
    }
    member_rows = (
        {
            "species_label": "Primary",
            "member_id": "p1",
            "internal_id": "0_0",
            "source_fasta": "Primary.fa",
            "raw_header": "p1 useful protein",
            "sequence_available": True,
            "sequence_length": 4,
            "motif_match": True,
            "matched_sequence": "N",
            "sequence": "MAAN",
        },
        {
            "species_label": "Primary",
            "member_id": "p2",
            "internal_id": "0_1",
            "source_fasta": "Primary.fa",
            "raw_header": "p2",
            "sequence_available": True,
            "sequence_length": 4,
            "motif_match": False,
            "matched_sequence": "",
            "sequence": "MAAA",
        },
        {
            "species_label": "Required",
            "member_id": "r1",
            "internal_id": "1_0",
            "source_fasta": "Required.fa",
            "raw_header": "",
            "sequence_available": True,
            "sequence_length": 4,
            "motif_match": False,
            "matched_sequence": "",
            "sequence": "MQQQ",
        },
        {
            "species_label": "Excluded",
            "member_id": "e1",
            "internal_id": "2_0",
            "source_fasta": "Excluded.fa",
            "raw_header": "e1 unavailable",
            "sequence_available": False,
            "sequence_length": None,
            "motif_match": False,
            "matched_sequence": "",
            "sequence": None,
        },
    )
    distribution = _species_distribution(
        rows=member_rows,
        authority=authority,
        roles=roles,
    )
    assert {row["Species result"] for row in distribution} == {
        "Some assessed proteins match",
        "Assessed; no match",
        "Sequence unavailable",
    }
    all_match = _species_distribution(
        rows=(member_rows[0],),
        authority=authority,
        roles=roles,
    )
    assert all_match[0]["Species result"] == "All assessed proteins match"
    display = _display_member(
        row=member_rows[0],
        authority=authority,
        roles=roles,
        show_sequence=True,
    )
    assert display["Protein description"] == "useful protein"
    assert display["Sequence"] == "MAAN"
    assert _protein_description(raw_header="different description", member_id="p1") == (
        "description"
    )
    assert _protein_description(raw_header="", member_id="p1") == ""
    assert _protein_description(raw_header="p1", member_id="p1") == ""
    fasta = _members_to_fasta(rows=member_rows)
    assert ">p1 species=Primary sequence_match=True\nMAAN" in fasta
    assert ">e1" not in fasta
    assert _members_to_fasta(rows=()) == ""
    evidence_figure = _species_evidence_figure(rows=distribution, group_id="N0.HOG1")
    assert len(evidence_figure.data) == 3

    search = validate_sequence_search(expression="N", mode=EXACT_C_TERMINAL)
    summaries = (
        {
            "group_id": "N0.HOG1",
            "total_member_count": 4,
            "sequence_count": 3,
            "sequence_coverage": 0.75,
            "analysis_sequence_count": 2,
            "analysis_matching_sequence_count": 1,
            "analysis_matching_fraction": 0.5,
            "analysis_assessed_species_count": 1,
            "analysis_matching_species_count": 1,
            "analysis_species_match_fraction": 1.0,
            "represented_species": "Primary; Required; Excluded",
            "assessed_species": "Primary; Required",
            "matching_species": "Primary",
        },
    )
    assert _display_summary(row=summaries[0])["Sequences assessed"] == 3
    assert len(_motif_figure(rows=summaries, search=search, threshold=0.8).data) == 1
    matrix = _taxonomic_heatmap(
        rows=summaries,
        species=("Primary", "Required", "Excluded", "Other"),
        search=search,
    )
    assert list(matrix.data[0].z[0]) == [3, 2, 1, 0]
    assert _search_stem(search=search) == "exact_c_terminal_N"
    regex = validate_sequence_search(expression="N.*$", mode=REGEX_ANYWHERE)
    assert _search_stem(search=regex) == "regex_anywhere_N"


def test_regex_group_and_member_controls_cover_explicit_modes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """UI control helpers keep regex opt-in and protein filters deterministic."""

    class SearchControls:
        """Minimal Streamlit control surface for sequence-search tests."""

        def __init__(self, *, location: str) -> None:
            self.location = location

        def subheader(self, _label: str) -> None:
            return None

        def toggle(self, _label: str, **_kwargs: object) -> bool:
            return True

        def radio(self, label: str, _options: object, **_kwargs: object) -> str:
            if label == "Grouping authority":
                return "Original OrthoFinder orthogroups"
            return self.location

        def text_input(self, _label: str, **_kwargs: object) -> str:
            return r"N[^P][ST]"

    monkeypatch.setattr(
        terminal_motif_page,
        "st",
        SearchControls(location="Anywhere in the protein"),
    )
    anywhere = _render_search_controls()
    assert anywhere.mode == REGEX_ANYWHERE
    assert _render_group_controls(resource=object()) == ("LEGACY_ORTHOGROUP", "")

    monkeypatch.setattr(
        terminal_motif_page,
        "st",
        SearchControls(location="C-terminus only"),
    )
    terminal = _render_search_controls()
    assert terminal.mode == REGEX_C_TERMINAL

    rows = (
        {
            "species_label": "Species_A",
            "member_id": "alpha_match",
            "sequence_available": True,
            "motif_match": True,
        },
        {
            "species_label": "Species_A",
            "member_id": "alpha_no_match",
            "sequence_available": True,
            "motif_match": False,
        },
        {
            "species_label": "Species_B",
            "member_id": "beta_missing",
            "sequence_available": False,
            "motif_match": False,
        },
    )
    roles = {"Species_A": "PRIMARY_ANALYSIS", "Species_B": "OTHER_SAMPLED"}

    class MemberControls:
        """Minimal Streamlit control surface for protein-row filters."""

        def __init__(
            self,
            *,
            state: str,
            selected_roles: tuple[str, ...] | None = None,
            selected_species: tuple[str, ...] | None = None,
            identifier: str = "",
        ) -> None:
            self.state = state
            self.selected_roles = selected_roles
            self.selected_species = selected_species
            self.identifier = identifier

        def columns(self, _widths: object) -> tuple[MemberControls, ...]:
            return (self, self, self)

        def selectbox(self, _label: str, _options: object) -> str:
            return self.state

        def multiselect(
            self, label: str, *, options: object, default: object
        ) -> object:
            del options
            if label == "Analysis role" and self.selected_roles is not None:
                return self.selected_roles
            if label == "Species" and self.selected_species is not None:
                return self.selected_species
            return default

        def text_input(self, _label: str, **_kwargs: object) -> str:
            return self.identifier

        def toggle(self, _label: str, **_kwargs: object) -> bool:
            return True

        def caption(self, _value: str) -> None:
            return None

    expected = {
        "Matching": ("alpha_match",),
        "Not matching": ("alpha_no_match",),
        "Sequence unavailable": ("beta_missing",),
    }
    for state, member_ids in expected.items():
        monkeypatch.setattr(
            terminal_motif_page,
            "st",
            MemberControls(state=state),
        )
        filtered, show_sequences = _filter_member_rows(rows=rows, roles=roles)
        assert tuple(row["member_id"] for row in filtered) == member_ids
        assert show_sequences

    monkeypatch.setattr(
        terminal_motif_page,
        "st",
        MemberControls(
            state="All proteins",
            selected_roles=("PRIMARY_ANALYSIS",),
            selected_species=("Species_A",),
            identifier="not-present",
        ),
    )
    filtered, _show_sequences = _filter_member_rows(rows=rows, roles=roles)
    assert filtered == ()


def test_custom_taxonomy_path_is_forwarded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A selected reviewed mapping is loaded with the exact sidecar labels."""

    expected = TaxonomyAuthority(records=(), expected_species=("Species_A",))
    mapping = tmp_path / "taxonomy.tsv"
    monkeypatch.setattr(
        "orthofinder_interrogation_app.terminal_motif_page.read_taxonomy_mapping",
        lambda **kwargs: (
            expected if kwargs == {"path": mapping, "expected_species": ("Species_A",)} else None
        ),
    )
    assert _load_taxonomy(species=("Species_A",), taxonomy_path_text=str(mapping)) is expected
    monkeypatch.setattr(
        "orthofinder_interrogation_app.terminal_motif_page.read_matching_bundled_taxonomy",
        lambda **kwargs: (
            expected if kwargs == {"expected_species": ("Species_A",)} else None
        ),
    )
    assert _load_taxonomy(species=("Species_A",), taxonomy_path_text="") is expected


def test_builder_rejects_partial_reconciliation(orthofinder2_results: Path, tmp_path: Path) -> None:
    """Missing source sequences cannot produce a misleading partial sidecar."""

    working = orthofinder2_results / "WorkingDirectory"
    (working / "Species0.fa").write_text(">0_0\nMAAN\n", encoding="utf-8")
    (working / "Species1.fa").write_text(">1_0\nMTTA\n", encoding="utf-8")
    with pytest.raises(InputValidationError, match="reconciliation failed"):
        build_sequence_sidecar(
            results_dir=orthofinder2_results,
            output_path=tmp_path / "partial.parquet",
        )
    assert (
        main(
            [
                "--orthofinder-results-dir",
                str(tmp_path / "missing_results"),
                "--output-parquet",
                str(tmp_path / "never.parquet"),
            ]
        )
        == 2
    )


def test_builder_accepts_original_fasta_names(orthofinder2_results: Path, tmp_path: Path) -> None:
    """Original FASTA names are a supported fallback to Species-index files."""

    working = orthofinder2_results / "WorkingDirectory"
    (working / "Species_A.fa").write_text(">0_0\nMAAN\n>0_1\nMQQN\n", encoding="utf-8")
    (working / "Species_B.faa").write_text(">1_0\nMTTA\n", encoding="utf-8")
    output = tmp_path / "fallback.parquet"
    assert (
        build_sequence_sidecar(
            results_dir=working,
            output_path=output,
        )
        == 3
    )


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

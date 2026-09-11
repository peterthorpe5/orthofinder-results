"""Tests for page guidance, methods and the scientific glossary."""

from __future__ import annotations

from pathlib import Path

import pytest

from orthofinder_interrogation_app import documentation_page
from orthofinder_interrogation_app.documentation_page import (
    GlossaryEntry,
    MethodStep,
)
from orthofinder_interrogation_app.models import ResourceIdentity
from orthofinder_results.errors import InputValidationError


class _FakeExpander:
    """Context manager for one captured guidance expander."""

    def __enter__(self) -> "_FakeExpander":
        """Enter the fake expander."""

        return self

    def __exit__(self, *args: object) -> None:
        """Leave without suppressing exceptions."""


class _FakeStreamlit:
    """Capture page-guidance headings and Markdown."""

    def __init__(self) -> None:
        """Initialise captured values."""

        self.labels: list[str] = []
        self.markdown_values: list[str] = []

    def expander(self, label: str, *, expanded: bool) -> _FakeExpander:
        """Capture the closed expander label."""

        assert not expanded
        self.labels.append(label)
        return _FakeExpander()

    def markdown(self, value: str) -> None:
        """Capture rendered Markdown."""

        self.markdown_values.append(value)


def test_page_guidance_registry_is_complete_and_result_led() -> None:
    """Every analytical and documentation page has substantive interpretation help."""

    assert len(documentation_page.PAGE_GUIDANCE) == 13
    for key, record in documentation_page.PAGE_GUIDANCE.items():
        assert key
        assert len(record.title) >= 5
        assert len(record.purpose) >= 50
        assert len(record.results) >= 50
        assert len(record.read_order) >= 50
        assert len(record.caution) >= 50
        assert documentation_page.page_guidance(key=key) is record
    with pytest.raises(InputValidationError, match="Unknown page guidance"):
        documentation_page.page_guidance(key="missing")


def test_render_page_guidance_exposes_results_and_qualification(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The closed page dropdown contains all required interpretation sections."""

    fake = _FakeStreamlit()
    monkeypatch.setattr(documentation_page, "st", fake)
    documentation_page.render_page_guidance(key="dispersion_benchmarks")
    assert fake.labels == ["How to read this page: Calibrated dispersion"]
    rendered = fake.markdown_values[0]
    assert "Purpose" in rendered
    assert "What the results show" in rendered
    assert "Recommended reading order" in rendered
    assert "Important qualification" in rendered
    monkeypatch.setattr(documentation_page, "st", None)
    with pytest.raises(RuntimeError, match="Streamlit"):
        documentation_page.render_page_guidance(key="overview")


def test_glossary_is_unique_searchable_and_export_ready() -> None:
    """The glossary covers core terms and supports literal topic-aware filtering."""

    entries = documentation_page.GLOSSARY_ENTRIES
    assert len(entries) >= 80
    terms = [entry.term for entry in entries]
    assert len(terms) == len(set(terms))
    for required in (
        "Hierarchical orthogroup (HOG)",
        "Patristic distance",
        "Principal coordinates analysis (PCoA)",
        "Benjamini–Hochberg FDR",
        "Matched residual",
        "Selection coverage tree",
    ):
        assert required in terms
    pcoa = documentation_page.filter_glossary_entries(
        entries=entries,
        query="pcoa",
        topic="All topics",
    )
    assert any(entry.term == "Principal coordinates analysis (PCoA)" for entry in pcoa)
    statistics = documentation_page.filter_glossary_entries(
        entries=entries,
        query="",
        topic="Statistics",
    )
    assert statistics and all(entry.topic == "Statistics" for entry in statistics)
    assert list(statistics) == sorted(statistics, key=lambda entry: entry.term.casefold())
    records = documentation_page.glossary_records(entries=pcoa)
    assert records[0]["Term"]
    assert "How to interpret it here" in records[0]
    with pytest.raises(InputValidationError, match="Unknown glossary topic"):
        documentation_page.filter_glossary_entries(
            entries=entries,
            query="",
            topic="Unknown",
        )


def test_method_steps_are_ordered_complete_and_defensive() -> None:
    """The method table follows the full workflow and rejects ambiguous ordering."""

    records = documentation_page.method_records()
    assert [record["Step"] for record in records] == list(range(1, 14))
    assert records[0]["Stage"] == "Completed OrthoFinder analysis"
    assert records[-1]["Stage"] == "Immutable publication and read-only interrogation"
    assert all(record["Important qualification"] for record in records)
    repeated = (
        MethodStep(1, "First", "method", "reason", "output", "caution"),
        MethodStep(1, "Again", "method", "reason", "output", "caution"),
    )
    with pytest.raises(InputValidationError, match="unique and increasing"):
        documentation_page.method_records(steps=repeated)


def test_resource_method_records_report_capabilities_without_mutation(
    tmp_path: Path,
) -> None:
    """Methods provenance distinguishes viewer, schema and stored capabilities."""

    benchmark_relations = frozenset(
        {
            "tree_payloads",
            "benchmark_group_profiles",
            "benchmark_cluster_results",
            "benchmark_contrasts",
            "benchmark_individual_comparisons",
            "benchmark_cluster_classifications",
        }
    )
    resource = ResourceIdentity(
        resource_path=tmp_path,
        database_path=tmp_path / "result.duckdb",
        report_path=None,
        run_id="fixture_run",
        schema_version=4,
        resource_package_version="0.9.0",
        orthofinder_version="2.5.5",
        adapter_name="orthofinder_2",
        primary_group_authority="HOG",
        relations=benchmark_relations,
    )
    records = documentation_page.resource_method_records(resource=resource)
    indexed = {record["Property"]: record["Value"] for record in records}
    assert indexed["Run ID"] == "fixture_run"
    assert indexed["Resource schema"] == "4"
    assert indexed["Portable tree payloads"] == "Available"
    assert indexed["Calibrated dispersion relations"] == "Available"
    assert indexed["Application access"] == "Read-only"


def test_documentation_pages_fail_explicitly_without_streamlit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Information pages never fail later with an opaque attribute error."""

    monkeypatch.setattr(documentation_page, "st", None)
    with pytest.raises(RuntimeError, match="Streamlit"):
        documentation_page.render_glossary_page()
    resource = ResourceIdentity(
        resource_path=Path("resource"),
        database_path=Path("resource.duckdb"),
        report_path=None,
        run_id="run",
        schema_version=4,
        resource_package_version="0.9.0",
        orthofinder_version="2.5.5",
        adapter_name="orthofinder_2",
        primary_group_authority="HOG",
    )
    with pytest.raises(RuntimeError, match="Streamlit"):
        documentation_page.render_methods_page(resource=resource)


def test_glossary_records_accept_a_small_external_subset() -> None:
    """Glossary table transformation remains generic and deterministic."""

    entry = GlossaryEntry(
        term="Term",
        topic="Topic",
        definition="Definition",
        interpretation="Interpretation",
        related_terms="Related",
    )
    assert documentation_page.glossary_records(entries=(entry,)) == (
        {
            "Term": "Term",
            "Topic": "Topic",
            "Definition": "Definition",
            "How to interpret it here": "Interpretation",
            "Related terms": "Related",
        },
    )

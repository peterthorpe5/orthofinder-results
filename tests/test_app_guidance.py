"""Tests for complete, consistent graph interpretation guidance."""

from __future__ import annotations

from pathlib import Path

import pytest

from orthofinder_interrogation_app import guidance
from orthofinder_results.errors import InputValidationError


class _FakeExpander:
    """Minimal context manager for one fake help section."""

    def __enter__(self) -> "_FakeExpander":
        """Enter the fake help panel."""

        return self

    def __exit__(self, *args: object) -> None:
        """Leave without suppressing exceptions."""


class _FakeStreamlit:
    """Capture one graph-help section."""

    def __init__(self) -> None:
        """Initialise captured labels and Markdown."""

        self.labels: list[str] = []
        self.markdown_values: list[str] = []

    def expander(self, label: str, *, expanded: bool) -> _FakeExpander:
        """Capture the closed help-panel label."""

        assert not expanded
        self.labels.append(label)
        return _FakeExpander()

    def markdown(self, value: str) -> None:
        """Capture rendered guidance text."""

        self.markdown_values.append(value)


def test_graph_registry_is_complete_and_plain_language() -> None:
    """Every supported graph has all three required interpretation fields."""

    assert len(guidance.GRAPH_GUIDANCE) == 21
    for key, record in guidance.GRAPH_GUIDANCE.items():
        assert key
        assert record.title
        assert len(record.shows) >= 40
        assert len(record.interpretation) >= 40
        assert len(record.limitation) >= 40
        assert guidance.graph_guidance(key=key) is record
    with pytest.raises(InputValidationError, match="Unknown graph"):
        guidance.graph_guidance(key="missing")


def test_render_graph_guidance_has_required_sections(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The UI visibly separates description, interpretation and limitation."""

    fake = _FakeStreamlit()
    monkeypatch.setattr(guidance, "st", fake)
    guidance.render_graph_guidance(key="shepard")
    assert fake.labels == ["How to read this graph: PCoA distance-fit (Shepard) plot"]
    rendered = fake.markdown_values[0]
    assert "What this graph shows" in rendered
    assert "How to interpret it" in rendered
    assert "Important limitation" in rendered
    monkeypatch.setattr(guidance, "st", None)
    with pytest.raises(RuntimeError, match="Streamlit"):
        guidance.render_graph_guidance(key="shepard")


def test_every_rendered_graph_has_a_guidance_panel() -> None:
    """Static source checks prevent new graphs from silently omitting help."""

    source_root = Path(__file__).resolve().parents[1] / "src"
    files = (
        source_root / "orthofinder_interrogation_app" / "evolutionary_page.py",
        source_root / "orthofinder_interrogation_app" / "comparison_page.py",
        source_root / "orthofinder_interrogation_app" / "benchmark_page.py",
    )
    source = "\n".join(path.read_text(encoding="utf-8") for path in files)
    graph_count = source.count("st.plotly_chart(") + source.count("st.iframe(")
    guidance_count = source.count("render_graph_guidance(key=")
    assert graph_count == 22
    assert guidance_count == graph_count

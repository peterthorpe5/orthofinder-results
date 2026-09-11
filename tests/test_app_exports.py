"""Unit tests for paired TSV and formatted Excel application exports."""

from __future__ import annotations

import ast
from datetime import date, datetime
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import pytest

from orthofinder_interrogation_app import exports
from orthofinder_results.errors import InputValidationError


def test_safe_file_stem_and_scalar_normalisation_are_defensive() -> None:
    """Filenames and uncommon scientific values remain portable and loss-aware."""

    assert exports.safe_file_stem(value="N0.HOG1: selected rows") == ("N0.HOG1_selected_rows")
    assert exports.normalise_excel_scalar(value=None) is None
    assert exports.normalise_excel_scalar(value={"b", "a"}) == '["a", "b"]'
    assert exports.normalise_excel_scalar(value=1234567890123456) == ("1234567890123456")
    assert exports.normalise_excel_scalar(value=float("inf")) == "inf"
    assert exports.normalise_excel_scalar(value=object()).startswith("<object object at")
    with pytest.raises(TypeError, match="must be text"):
        exports.safe_file_stem(value=1)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="cannot be empty"):
        exports.safe_file_stem(value=" /// ")


@pytest.mark.parametrize(
    ("column_name", "values", "expected"),
    [
        ("Group ID", (123,), "text"),
        ("Protein pairs", (1, 2), "integer"),
        ("Average distance", (0.1, 0.2), "decimal"),
        ("Adjusted q value", (1.0e-8, 0.05), "scientific"),
        ("Sampling fraction", (0.25, 1.0), "percentage"),
        ("Full-group matrix", (True, False), "logical"),
        ("Source date", (date(2026, 9, 8),), "date"),
        ("Created", (datetime(2026, 9, 8, 12, 30),), "datetime"),
        ("Notes", ("review",), "text"),
        ("Optional distance", (None,), "text"),
    ],
)
def test_excel_format_kind_uses_scientific_semantics(
    column_name: str,
    values: tuple[object, ...],
    expected: str,
) -> None:
    """Identifiers, counts, fractions and measurements receive stable formats."""

    assert (
        exports.excel_format_kind(
            column_name=column_name,
            values=values,
        )
        == expected
    )
    with pytest.raises(ValueError, match="non-empty"):
        exports.excel_format_kind(column_name="", values=values)


def test_excel_column_width_is_readable_and_bounded() -> None:
    """Column widths expand for content without becoming unusably narrow or wide."""

    assert exports.excel_column_width(column_name="Rank", values=(1, 2)) == 12.0
    assert exports.excel_column_width(column_name="Description", values=("x" * 200,)) == 50.0
    with pytest.raises(ValueError, match="non-empty"):
        exports.excel_column_width(column_name="", values=(1,))


def test_records_to_excel_bytes_has_filters_freeze_formats_and_dictionary() -> None:
    """The workbook is a typed, filterable table with a definitions worksheet."""

    records = (
        {
            "Group ID": "=2+2",
            "Protein pairs": 3,
            "Average distance": 0.123456,
            "Sampling fraction": 0.5,
            "Full-group matrix": True,
            "Source date": date(2026, 9, 8),
        },
        {
            "Group ID": "N0.HOG2",
            "Protein pairs": 6,
            "Average distance": 1.5,
            "Sampling fraction": 1.0,
            "Full-group matrix": False,
            "Source date": date(2026, 9, 9),
        },
    )
    payload = exports.records_to_excel_bytes(
        records=records,
        column_definitions={"Average distance": "Mean exact pair distance."},
        workbook_title="Fixture results",
    )
    assert payload.startswith(b"PK")

    with ZipFile(BytesIO(payload)) as archive:
        sheet_xml = archive.read("xl/worksheets/sheet1.xml").decode("utf-8")
        dictionary_xml = archive.read("xl/worksheets/sheet2.xml").decode("utf-8")
        table_xml = archive.read("xl/tables/table1.xml").decode("utf-8")
        dictionary_table = archive.read("xl/tables/table2.xml").decode("utf-8")
        styles_xml = archive.read("xl/styles.xml").decode("utf-8")
        strings_xml = archive.read("xl/sharedStrings.xml").decode("utf-8")
        properties_xml = archive.read("docProps/core.xml").decode("utf-8")

    assert 'state="frozen"' in sheet_xml
    assert 'ySplit="1"' in sheet_xml
    assert 'showGridLines="0"' in sheet_xml
    assert "<cols>" in sheet_xml
    assert "<f>" not in sheet_xml
    assert '<autoFilter ref="A1:F3"' in table_xml
    assert 'name="TableStyleMedium2"' in table_xml
    assert "0.0000" in styles_xml
    assert "0.0%" in styles_xml
    assert "yyyy-mm-dd" in styles_xml
    assert "=2+2" in strings_xml
    assert "Mean exact pair distance." in strings_xml
    assert "Plain-language definition" in strings_xml
    assert 'state="frozen"' in dictionary_xml
    assert 'name="OrthoFinderColumnDefinitions"' in dictionary_table
    assert "Fixture results" in properties_xml


def test_records_to_excel_bytes_rejects_unsafe_shapes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Invalid headings, undeclared values and excessive dimensions fail clearly."""

    with pytest.raises(InputValidationError, match="non-empty field"):
        exports.records_to_excel_bytes(records=())
    with pytest.raises(InputValidationError, match="must be unique"):
        exports.records_to_excel_bytes(
            records=({"one": 1},),
            fieldnames=("one", "one"),
        )
    with pytest.raises(InputValidationError, match="undeclared fields"):
        exports.records_to_excel_bytes(
            records=({"one": 1, "two": 2},),
            fieldnames=("one",),
        )
    monkeypatch.setattr(exports, "MAX_EXCEL_DATA_ROWS", 1)
    with pytest.raises(InputValidationError, match="at most 1"):
        exports.records_to_excel_bytes(records=({"one": 1}, {"one": 2}))
    monkeypatch.setattr(exports, "MAX_EXCEL_COLUMNS", 1)
    with pytest.raises(InputValidationError, match="at most 1"):
        exports.records_to_excel_bytes(records=({"one": 1, "two": 2},))


def test_empty_excel_export_keeps_headers_and_filters() -> None:
    """An explicit empty selection remains a useful, documented workbook."""

    payload = exports.records_to_excel_bytes(
        records=(),
        fieldnames=("Group ID", "Average distance"),
    )
    with ZipFile(BytesIO(payload)) as archive:
        sheet_xml = archive.read("xl/worksheets/sheet1.xml").decode("utf-8")
        strings_xml = archive.read("xl/sharedStrings.xml").decode("utf-8")
    assert '<autoFilter ref="A1:B1"' in sheet_xml
    assert "Group ID" in strings_xml


class _FakeColumn:
    """Context manager returned by the fake Streamlit column layout."""

    def __enter__(self) -> "_FakeColumn":
        """Enter the fake column context."""

        return self

    def __exit__(self, *args: object) -> None:
        """Leave the fake context without suppressing errors."""


class _FakeStreamlit:
    """Capture paired download calls without starting a server."""

    def __init__(self) -> None:
        """Initialise an empty call list."""

        self.calls: list[dict[str, object]] = []

    def columns(self, *, spec: int) -> tuple[_FakeColumn, _FakeColumn]:
        """Return exactly two fake download columns."""

        assert spec == 2
        return _FakeColumn(), _FakeColumn()

    def download_button(self, **kwargs: object) -> None:
        """Record one rendered download control."""

        self.calls.append(kwargs)


def test_render_table_downloads_preserves_tsv_and_adds_excel(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Each table receives neighbouring TSV and formatted Excel controls."""

    fake = _FakeStreamlit()
    monkeypatch.setattr(exports, "st", fake)
    exports.render_table_downloads(
        records=({"Group ID": "N0.HOG1", "Average distance": 0.2},),
        file_stem="N0.HOG1 results",
        key="download",
        column_definitions={"Average distance": "Mean pair distance."},
    )
    assert [call["file_name"] for call in fake.calls] == [
        "N0.HOG1_results.tsv",
        "N0.HOG1_results.xlsx",
    ]
    assert fake.calls[0]["mime"] == "text/tab-separated-values"
    assert fake.calls[1]["mime"] == exports.EXCEL_MIME_TYPE
    assert bytes(fake.calls[1]["data"]).startswith(b"PK")
    with pytest.raises(ValueError, match="non-empty key"):
        exports.render_table_downloads(
            records=({"Group ID": "g"},),
            file_stem="result",
            key="",
        )
    monkeypatch.setattr(exports, "st", None)
    with pytest.raises(RuntimeError, match="Streamlit"):
        exports.render_table_downloads(
            records=({"Group ID": "g"},),
            file_stem="result",
            key="result",
        )


class _FakeFigure:
    """Return deterministic PDF bytes from a Plotly-compatible interface."""

    def __init__(self, *, payload: object = b"%PDF-1.7\nfixture") -> None:
        """Store the renderer payload and captured keyword arguments."""

        self.payload = payload
        self.calls: list[dict[str, object]] = []

    def to_image(self, **kwargs: object) -> object:
        """Capture export arguments and return the configured payload."""

        self.calls.append(kwargs)
        return self.payload


class _BrokenFigure:
    """Raise an environmental error from the static renderer."""

    def to_image(self, **kwargs: object) -> bytes:
        """Raise a renderer failure for defensive error testing."""

        del kwargs
        raise RuntimeError("missing browser")


class _FakePlotStreamlit:
    """Capture one Plotly render and its deferred download control."""

    def __init__(self) -> None:
        """Initialise captured chart and button calls."""

        self.charts: list[dict[str, object]] = []
        self.calls: list[dict[str, object]] = []

    def plotly_chart(self, figure: object, **kwargs: object) -> str:
        """Capture the figure and return a stable event sentinel."""

        self.charts.append({"figure": figure, **kwargs})
        return "event"

    def download_button(self, **kwargs: object) -> None:
        """Capture one download control."""

        self.calls.append(kwargs)


def test_plotly_figure_pdf_bytes_validate_renderer_and_dimensions() -> None:
    """PDF generation accepts valid bytes and explains renderer failures."""

    figure = _FakeFigure()
    payload = exports.plotly_figure_to_pdf_bytes(
        figure=figure,
        width=1200,
        height=800,
    )
    assert payload.startswith(b"%PDF")
    assert figure.calls == [{"format": "pdf", "width": 1200, "height": 800, "scale": 1}]
    with pytest.raises(InputValidationError, match="width"):
        exports.plotly_figure_to_pdf_bytes(figure=figure, width=100)
    with pytest.raises(InputValidationError, match="height"):
        exports.plotly_figure_to_pdf_bytes(figure=figure, height=100)
    with pytest.raises(InputValidationError, match="Plotly-compatible"):
        exports.plotly_figure_to_pdf_bytes(figure=object())
    with pytest.raises(InputValidationError, match="Kaleido"):
        exports.plotly_figure_to_pdf_bytes(figure=_BrokenFigure())
    with pytest.raises(InputValidationError, match="valid PDF"):
        exports.plotly_figure_to_pdf_bytes(figure=_FakeFigure(payload=b"not-pdf"))


def test_render_plotly_figure_adds_deferred_pdf_download(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Every Plotly result receives an on-demand PDF generated from that figure."""

    fake = _FakePlotStreamlit()
    figure = _FakeFigure()
    monkeypatch.setattr(exports, "st", fake)
    event = exports.render_plotly_figure(
        figure=figure,
        file_stem="N0.HOG1 result figure",
        key="figure_key",
        on_select="rerun",
        selection_mode="points",
    )
    assert event == "event"
    assert fake.charts[0]["key"] == "figure_key"
    assert fake.charts[0]["on_select"] == "rerun"
    assert fake.calls[0]["file_name"] == "N0.HOG1_result_figure.pdf"
    assert fake.calls[0]["mime"] == exports.PDF_MIME_TYPE
    deferred = fake.calls[0]["data"]
    assert callable(deferred)
    assert deferred().startswith(b"%PDF")
    with pytest.raises(ValueError, match="non-empty key"):
        exports.render_plotly_figure(
            figure=figure,
            file_stem="figure",
            key="",
        )
    monkeypatch.setattr(exports, "st", None)
    with pytest.raises(RuntimeError, match="Streamlit"):
        exports.render_plotly_figure(
            figure=figure,
            file_stem="figure",
            key="figure",
        )


def test_interactive_html_download_is_validated_and_preserves_document(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Interactive networks retain a self-contained HTML download beside the PDF view."""

    fake = _FakePlotStreamlit()
    monkeypatch.setattr(exports, "st", fake)
    exports.render_html_download(
        document="<!doctype html><html><body>network</body></html>",
        file_stem="interactive network",
        key="network_html",
    )
    assert fake.calls[0]["file_name"] == "interactive_network.html"
    assert fake.calls[0]["mime"] == exports.HTML_MIME_TYPE
    assert "network" in str(fake.calls[0]["data"])
    with pytest.raises(InputValidationError, match="non-empty document"):
        exports.render_html_download(document=" ", file_stem="network", key="key")
    with pytest.raises(ValueError, match="non-empty key"):
        exports.render_html_download(document="<html />", file_stem="network", key="")
    monkeypatch.setattr(exports, "st", None)
    with pytest.raises(RuntimeError, match="Streamlit"):
        exports.render_html_download(
            document="<html />",
            file_stem="network",
            key="key",
        )


def _named_call_count(*, node: ast.AST, function_name: str) -> int:
    """Count calls to one imported function name below an AST node."""

    return sum(
        isinstance(candidate, ast.Call)
        and isinstance(candidate.func, ast.Name)
        and candidate.func.id == function_name
        for candidate in ast.walk(node)
    )


def _streamlit_call_count(*, node: ast.AST, method_name: str) -> int:
    """Count calls to one ``st`` method below an AST node."""

    return sum(
        isinstance(candidate, ast.Call)
        and isinstance(candidate.func, ast.Attribute)
        and candidate.func.attr == method_name
        and isinstance(candidate.func.value, ast.Name)
        and candidate.func.value.id == "st"
        for candidate in ast.walk(node)
    )


def test_every_result_table_and_figure_has_its_required_export_controls() -> None:
    """Guard the app-wide TSV, Excel, PDF and interactive-export contract."""

    package = Path(__file__).resolve().parents[1] / "src/orthofinder_interrogation_app"
    page_paths = tuple(sorted(package.glob("*_page.py"))) + (package / "app.py",)
    audited_tables = 0
    plotted_figures = 0
    for path in page_paths:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for function in (
            node
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        ):
            table_count = _streamlit_call_count(node=function, method_name="dataframe")
            download_count = _named_call_count(
                node=function,
                function_name="render_table_downloads",
            )
            assert download_count >= table_count, (
                f"{path.name}:{function.name} renders {table_count} table(s) but only "
                f"{download_count} paired TSV/Excel control(s)."
            )
            audited_tables += table_count
        plotted_figures += _named_call_count(
            node=tree,
            function_name="render_plotly_figure",
        )
    assert audited_tables >= 25
    assert plotted_figures >= 20

    coverage_source = (package / "coverage_page.py").read_text(encoding="utf-8")
    assert "st.plotly_chart(" in coverage_source
    assert 'files["selection_coverage_tree.pdf"]' in coverage_source

    evolutionary_source = (package / "evolutionary_page.py").read_text(
        encoding="utf-8"
    )
    assert "st.iframe(" in evolutionary_source
    assert "render_html_download(" in evolutionary_source
    assert "_render_topology(" in evolutionary_source

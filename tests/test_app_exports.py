"""Unit tests for paired TSV and formatted Excel application exports."""

from __future__ import annotations

from datetime import date, datetime
from io import BytesIO
from zipfile import ZipFile

import pytest

from orthofinder_interrogation_app import exports
from orthofinder_results.errors import InputValidationError


def test_safe_file_stem_and_scalar_normalisation_are_defensive() -> None:
    """Filenames and uncommon scientific values remain portable and loss-aware."""

    assert exports.safe_file_stem(value="N0.HOG1: selected rows") == (
        "N0.HOG1_selected_rows"
    )
    assert exports.normalise_excel_scalar(value=None) is None
    assert exports.normalise_excel_scalar(value={"b", "a"}) == '["a", "b"]'
    assert exports.normalise_excel_scalar(value=1234567890123456) == (
        "1234567890123456"
    )
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

    assert exports.excel_format_kind(
        column_name=column_name,
        values=values,
    ) == expected
    with pytest.raises(ValueError, match="non-empty"):
        exports.excel_format_kind(column_name="", values=values)


def test_excel_column_width_is_readable_and_bounded() -> None:
    """Column widths expand for content without becoming unusably narrow or wide."""

    assert exports.excel_column_width(column_name="Rank", values=(1, 2)) == 12.0
    assert exports.excel_column_width(
        column_name="Description", values=("x" * 200,)
    ) == 50.0
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

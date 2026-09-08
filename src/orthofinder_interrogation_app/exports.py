"""Safe paired TSV and formatted Excel exports for application tables."""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from datetime import date, datetime
from io import BytesIO
from numbers import Integral, Real
from typing import Any

import xlsxwriter

from orthofinder_results.errors import InputValidationError

from .tsv import records_to_tsv

try:
    import streamlit as st
except ModuleNotFoundError:  # pragma: no cover - checked by dependency smoke tests.
    st = None  # type: ignore[assignment]

EXCEL_MIME_TYPE = (
    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)
MAX_EXCEL_DATA_ROWS = 1_048_575
MAX_EXCEL_COLUMNS = 16_384
_INVALID_FILE_STEM = re.compile(r"[^A-Za-z0-9_.-]+")
_IDENTIFIER_COLUMN = re.compile(
    r"(^|[ _])(accession|checksum|digest|identifier|id)([ _]|$)",
    flags=re.IGNORECASE,
)
_INTEGER_COLUMN = re.compile(
    r"(^|[ _])(count|index|number|rank|members?|proteins?|pairs?|species)([ _]|$)",
    flags=re.IGNORECASE,
)
_SCIENTIFIC_COLUMN = re.compile(
    r"(^|[ _])(e[ _]?value|fdr|p[ _]?value|q[ _]?value)([ _]|$)",
    flags=re.IGNORECASE,
)
_PERCENTAGE_COLUMN = re.compile(
    r"fraction|coverage|purity|proportion|share|sampling percentage",
    flags=re.IGNORECASE,
)
_NARRATIVE_COLUMN = re.compile(
    r"definition|description|caution|reason|note|source file|limitation",
    flags=re.IGNORECASE,
)


def safe_file_stem(*, value: str) -> str:
    """Return a portable, non-empty download filename stem.

    Args:
        value: Proposed filename without an extension.

    Returns:
        Filename-safe stem containing letters, numbers, dots, dashes and
        underscores.

    Raises:
        TypeError: If ``value`` is not text.
        ValueError: If normalisation removes the complete value.
    """

    if not isinstance(value, str):
        raise TypeError("The export filename stem must be text.")
    normalised = _INVALID_FILE_STEM.sub("_", value.strip()).strip("._")
    if not normalised:
        raise ValueError("The export filename stem cannot be empty.")
    return normalised


def normalise_excel_scalar(*, value: Any) -> Any:
    """Return a safe scalar accepted by XlsxWriter.

    Args:
        value: Raw table value.

    Returns:
        A string, finite number, Boolean, date/time or ``None``. Collections
        become deterministic JSON text and integers longer than Excel's
        precise 15-digit range become text.
    """

    if value is None:
        return None
    if isinstance(value, (list, tuple, dict, set)):
        serialisable = sorted(value) if isinstance(value, set) else value
        return json.dumps(serialisable, ensure_ascii=False, default=str)
    if isinstance(value, Real) and not isinstance(value, (bool, Integral)):
        numeric = float(value)
        return numeric if math.isfinite(numeric) else str(numeric)
    if isinstance(value, Integral) and not isinstance(value, bool):
        integer = int(value)
        return integer if len(str(abs(integer))) <= 15 else str(integer)
    if isinstance(value, (str, bool, date, datetime)):
        return value
    return str(value)


def excel_format_kind(*, column_name: str, values: Sequence[Any]) -> str:
    """Choose a conservative Excel format for one scientific column.

    Args:
        column_name: User-facing column heading.
        values: Raw values from the column.

    Returns:
        One of ``text``, ``integer``, ``decimal``, ``scientific``,
        ``percentage``, ``date``, ``datetime`` or ``logical``.

    Raises:
        ValueError: If the column heading is empty.
    """

    if not isinstance(column_name, str) or not column_name.strip():
        raise ValueError("Excel columns must have non-empty text headings.")
    if _IDENTIFIER_COLUMN.search(column_name.replace("_", " ")):
        return "text"
    sample = tuple(
        normalise_excel_scalar(value=value)
        for value in values[:1000]
        if value is not None
    )
    sample = tuple(value for value in sample if value is not None)
    if sample and all(isinstance(value, bool) for value in sample):
        return "logical"
    if sample and all(isinstance(value, datetime) for value in sample):
        return "datetime"
    if sample and all(isinstance(value, date) for value in sample):
        return "date"
    numeric = sample and all(
        isinstance(value, Real) and not isinstance(value, bool) for value in sample
    )
    if numeric:
        if _SCIENTIFIC_COLUMN.search(column_name.replace("_", " ")):
            return "scientific"
        if _PERCENTAGE_COLUMN.search(column_name.replace("_", " ")):
            return "percentage"
        integer_like = all(float(value).is_integer() for value in sample)
        if integer_like or _INTEGER_COLUMN.search(column_name.replace("_", " ")):
            return "integer"
        return "decimal"
    return "text"


def excel_column_width(*, column_name: str, values: Sequence[Any]) -> float:
    """Return a readable Excel column width bounded to 12–50 characters.

    Args:
        column_name: User-facing column heading.
        values: Raw values from the column.

    Returns:
        Width in Excel character units.

    Raises:
        ValueError: If the heading is empty.
    """

    if not isinstance(column_name, str) or not column_name.strip():
        raise ValueError("Excel columns must have non-empty text headings.")
    lengths = [len(column_name)]
    for value in values[:500]:
        normalised = normalise_excel_scalar(value=value)
        if normalised is not None:
            lengths.append(len(str(normalised)))
    content_width = max(lengths, default=len(column_name)) + 2
    if _NARRATIVE_COLUMN.search(column_name):
        content_width = max(content_width, 28)
    return float(min(50, max(12, content_width)))


def records_to_excel_bytes(
    *,
    records: Sequence[Mapping[str, Any]],
    fieldnames: Sequence[str] | None = None,
    column_definitions: Mapping[str, str] | None = None,
    workbook_title: str = "OrthoFinder results",
) -> bytes:
    """Create a filterable, formatted scientific workbook from records.

    The first sheet contains the exact supplied rows. The second sheet defines
    every selected column. Text is always written as text, so identifiers that
    begin with ``=``, ``+``, ``-`` or ``@`` cannot become Excel formulas.

    Args:
        records: Ordered table records.
        fieldnames: Optional explicit output headings and order.
        column_definitions: Optional plain-language descriptions keyed by heading.
        workbook_title: Short title stored in workbook metadata.

    Returns:
        Complete XLSX workbook bytes.

    Raises:
        InputValidationError: If the table shape, headings or records are unsafe.
    """

    headings = _validated_headings(records=records, fieldnames=fieldnames)
    if len(records) > MAX_EXCEL_DATA_ROWS:
        raise InputValidationError(
            f"Excel export contains {len(records):,} data rows; the format allows "
            f"at most {MAX_EXCEL_DATA_ROWS:,}. Narrow the result first."
        )
    if len(headings) > MAX_EXCEL_COLUMNS:
        raise InputValidationError(
            f"Excel export contains {len(headings):,} columns; the format allows "
            f"at most {MAX_EXCEL_COLUMNS:,}."
        )
    expected = set(headings)
    for index, record in enumerate(records, start=1):
        unknown = sorted(set(record).difference(expected))
        if unknown:
            raise InputValidationError(
                f"Excel row {index} contains undeclared fields: {'; '.join(unknown)}"
            )

    columns = {
        heading: tuple(record.get(heading) for record in records) for heading in headings
    }
    output = BytesIO()
    workbook = xlsxwriter.Workbook(
        output,
        {
            "in_memory": True,
            "strings_to_formulas": False,
            "strings_to_urls": False,
        },
    )
    workbook.set_properties(
        {
            "title": workbook_title,
            "subject": "Filterable export from OrthoFinder Interrogation",
            "author": "Peter Thorpe and collaborators",
            "comments": (
                "The Results sheet contains the exact bounded rows selected in the app."
            ),
        }
    )
    header_format = workbook.add_format(
        {
            "bold": True,
            "font_color": "#FFFFFF",
            "bg_color": "#1F4E78",
            "border": 1,
            "border_color": "#A6A6A6",
            "align": "left",
            "valign": "vcenter",
        }
    )
    cell_formats = _workbook_cell_formats(workbook=workbook)
    worksheet = workbook.add_worksheet("Results")
    worksheet.hide_gridlines(option=2)
    worksheet.freeze_panes(row=1, col=0)
    worksheet.set_zoom(zoom=90)
    worksheet.set_row(row=0, height=26)
    table_columns = []
    for column_index, heading in enumerate(headings):
        values = columns[heading]
        kind = excel_format_kind(column_name=heading, values=values)
        worksheet.set_column(
            first_col=column_index,
            last_col=column_index,
            width=excel_column_width(column_name=heading, values=values),
            cell_format=cell_formats[kind],
        )
        worksheet.write_string(0, column_index, heading, header_format)
        table_columns.append({"header": heading, "header_format": header_format})
        identifier_column = bool(_IDENTIFIER_COLUMN.search(heading.replace("_", " ")))
        for row_index, raw_value in enumerate(values, start=1):
            _write_excel_cell(
                worksheet=worksheet,
                row_index=row_index,
                column_index=column_index,
                value=normalise_excel_scalar(value=raw_value),
                cell_format=cell_formats[kind],
                identifier_column=identifier_column,
            )
    if records:
        worksheet.add_table(
            first_row=0,
            first_col=0,
            last_row=len(records),
            last_col=len(headings) - 1,
            options={
                "name": "OrthoFinderResultsTable",
                "style": "Table Style Medium 2",
                "columns": table_columns,
                "autofilter": True,
                "banded_rows": True,
            },
        )
    else:
        worksheet.autofilter(0, 0, 0, len(headings) - 1)
    _write_definition_sheet(
        workbook=workbook,
        headings=headings,
        columns=columns,
        definitions=column_definitions or {},
        header_format=header_format,
    )
    workbook.close()
    return output.getvalue()


def render_table_downloads(
    *,
    records: Sequence[Mapping[str, Any]],
    file_stem: str,
    key: str,
    tsv_label: str = "Download as TSV",
    excel_label: str = "Download as formatted Excel",
    fieldnames: Sequence[str] | None = None,
    column_definitions: Mapping[str, str] | None = None,
    workbook_title: str = "OrthoFinder results",
) -> None:
    """Render neighbouring TSV and formatted Excel download controls.

    Args:
        records: Exact records exported by both controls.
        file_stem: Shared safe filename without an extension.
        key: Stable Streamlit key for the TSV button.
        tsv_label: User-facing TSV button label.
        excel_label: User-facing Excel button label.
        fieldnames: Optional explicit output headings and order.
        column_definitions: Optional definitions included in the workbook.
        workbook_title: Workbook metadata title.

    Raises:
        ValueError: If ``key`` is empty.
        RuntimeError: If Streamlit is unavailable.
    """

    if not isinstance(key, str) or not key.strip():
        raise ValueError("Download controls require a non-empty key.")
    if st is None:
        raise RuntimeError("Streamlit is required to render download controls.")
    stem = safe_file_stem(value=file_stem)
    headings = _validated_headings(records=records, fieldnames=fieldnames)
    tsv_column, excel_column = st.columns(spec=2)
    with tsv_column:
        st.download_button(
            label=tsv_label,
            data=records_to_tsv(records=records, fieldnames=headings),
            file_name=f"{stem}.tsv",
            mime="text/tab-separated-values",
            key=key,
        )
    with excel_column:
        st.download_button(
            label=excel_label,
            data=records_to_excel_bytes(
                records=records,
                fieldnames=headings,
                column_definitions=column_definitions,
                workbook_title=workbook_title,
            ),
            file_name=f"{stem}.xlsx",
            mime=EXCEL_MIME_TYPE,
            key=f"{key}_excel",
        )


def _validated_headings(
    *,
    records: Sequence[Mapping[str, Any]],
    fieldnames: Sequence[str] | None,
) -> tuple[str, ...]:
    """Return validated output headings in deterministic order."""

    headings = tuple(str(value) for value in (fieldnames or (tuple(records[0]) if records else ())))
    if not headings or any(not heading.strip() for heading in headings):
        raise InputValidationError("Table export requires non-empty field names.")
    if len(set(headings)) != len(headings):
        raise InputValidationError("Table export field names must be unique.")
    return headings


def _workbook_cell_formats(*, workbook: Any) -> dict[str, Any]:
    """Return the complete semantic Excel format registry."""

    base = {
        "align": "left",
        "valign": "vcenter",
        "border": 1,
        "border_color": "#D9E2F3",
    }
    numeric = {**base, "align": "right"}
    return {
        "text": workbook.add_format({**base, "text_wrap": False}),
        "integer": workbook.add_format({**numeric, "num_format": "#,##0"}),
        "decimal": workbook.add_format({**numeric, "num_format": "0.0000"}),
        "scientific": workbook.add_format({**numeric, "num_format": "0.00E+00"}),
        "percentage": workbook.add_format({**numeric, "num_format": "0.0%"}),
        "date": workbook.add_format({**base, "num_format": "yyyy-mm-dd"}),
        "datetime": workbook.add_format(
            {**base, "num_format": "yyyy-mm-dd hh:mm:ss"}
        ),
        "logical": workbook.add_format({**base, "align": "centre"}),
    }


def _write_excel_cell(
    *,
    worksheet: Any,
    row_index: int,
    column_index: int,
    value: Any,
    cell_format: Any,
    identifier_column: bool,
) -> None:
    """Write one safe, typed Excel cell without formula interpretation."""

    if value is None:
        worksheet.write_blank(row_index, column_index, None, cell_format)
    elif isinstance(value, str) or identifier_column:
        worksheet.write_string(row_index, column_index, str(value), cell_format)
    elif isinstance(value, bool):
        worksheet.write_boolean(row_index, column_index, value, cell_format)
    elif isinstance(value, (date, datetime)):
        worksheet.write_datetime(row_index, column_index, value, cell_format)
    else:
        worksheet.write_number(row_index, column_index, float(value), cell_format)


def _write_definition_sheet(
    *,
    workbook: Any,
    headings: tuple[str, ...],
    columns: Mapping[str, Sequence[Any]],
    definitions: Mapping[str, str],
    header_format: Any,
) -> None:
    """Add a filterable plain-language column dictionary worksheet."""

    worksheet = workbook.add_worksheet("Column definitions")
    worksheet.hide_gridlines(option=2)
    worksheet.freeze_panes(row=1, col=0)
    worksheet.set_zoom(zoom=90)
    dictionary_headings = ("Column", "Plain-language definition", "Excel data type")
    widths = (34, 72, 20)
    for index, (heading, width) in enumerate(zip(dictionary_headings, widths, strict=True)):
        worksheet.set_column(index, index, width)
        worksheet.write_string(0, index, heading, header_format)
    cell_format = workbook.add_format(
        {
            "align": "left",
            "valign": "top",
            "border": 1,
            "border_color": "#D9E2F3",
            "text_wrap": True,
        }
    )
    for row_index, heading in enumerate(headings, start=1):
        values = columns[heading]
        row = (
            heading,
            definitions.get(heading, "No additional definition supplied."),
            excel_format_kind(column_name=heading, values=values),
        )
        for column_index, value in enumerate(row):
            worksheet.write_string(row_index, column_index, str(value), cell_format)
        worksheet.set_row(row_index, 36)
    worksheet.add_table(
        first_row=0,
        first_col=0,
        last_row=len(headings),
        last_col=len(dictionary_headings) - 1,
        options={
            "name": "OrthoFinderColumnDefinitions",
            "style": "Table Style Medium 2",
            "columns": [
                {"header": heading, "header_format": header_format}
                for heading in dictionary_headings
            ],
            "autofilter": True,
            "banded_rows": True,
        },
    )

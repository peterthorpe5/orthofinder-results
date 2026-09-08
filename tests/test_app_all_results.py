"""Tests for the complete persisted-distance results page."""

from __future__ import annotations

import pytest

from orthofinder_interrogation_app import all_results_page
from orthofinder_results.errors import InputValidationError


def test_display_result_rows_preserves_selection_order_and_root_label() -> None:
    """Raw authority fields map to selected readable columns without reordering."""

    rows = (
        {
            "group_id": "OG1",
            "hierarchy_node": "",
            "mean_distance": 0.25,
        },
    )
    displayed = all_results_page._display_result_rows(
        rows=rows,
        selected_labels=("Group ID", "Species-tree level", "Average distance"),
    )
    assert displayed == (
        {
            "Group ID": "OG1",
            "Species-tree level": "ROOT",
            "Average distance": 0.25,
        },
    )
    with pytest.raises(InputValidationError, match="Unsupported displayed"):
        all_results_page._display_result_rows(
            rows=rows,
            selected_labels=("Unknown",),
        )


def test_all_result_column_dictionary_is_complete_and_distance_focused() -> None:
    """Every selectable field has useful workbook and hover help."""

    definitions = all_results_page.all_result_column_definitions()
    assert tuple(definitions) == tuple(all_results_page._RESULT_COLUMNS)
    assert "Average distance" in definitions
    assert "Distance spread (SD)" in definitions
    assert "Protein pairs" in definitions
    assert all(len(value) >= 20 for value in definitions.values())
    assert set(all_results_page._DEFAULT_COLUMNS).issubset(definitions)

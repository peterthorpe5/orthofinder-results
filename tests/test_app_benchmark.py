"""Tests for calibrated dispersion page helpers and figures."""

from __future__ import annotations

import pytest

from orthofinder_interrogation_app.benchmark_page import (
    _cluster_label,
    _display_contrast,
    _display_individual,
    _display_profile,
    _metric_from_label,
)
from orthofinder_interrogation_app.figures import (
    benchmark_contrast_figure,
    benchmark_distribution_figure,
    benchmark_individual_figure,
)


def test_benchmark_distribution_and_contrast_figures_retain_cluster_units() -> None:
    """Profile points and effect intervals remain labelled as cluster summaries."""

    distributions = (
        {
            "profile_id": "E3_ALL",
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": "N0.HOG1",
            "metric_value": 0.4,
        },
        {
            "profile_id": "HOUSEKEEPING_ALL",
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": "N0.HOG2",
            "metric_value": -0.1,
        },
    )
    figure = benchmark_distribution_figure(
        rows=distributions,
        metric_label="Average pair distance",
        scale_label="Matched residual",
    )
    assert len(figure.data) == 2
    assert figure.layout.yaxis.title.text.endswith("(Matched residual)")

    contrast = {
        "status": "TESTED",
        "target_profile_id": "E3_ALL",
        "reference_profile_id": "HOUSEKEEPING_ALL",
        "median_difference": 0.5,
        "median_difference_ci_low": 0.2,
        "median_difference_ci_high": 0.8,
        "cliffs_delta": 0.7,
        "p_value_two_sided": 0.01,
        "fdr_q_value": 0.02,
    }
    forest = benchmark_contrast_figure(
        rows=(contrast,),
        metric_label="Average pair distance",
    )
    assert forest.data[0].x[0] == pytest.approx(0.5)
    assert forest.data[0].marker.color[0] == "#b91c1c"


def test_individual_background_figure_and_display_helpers_are_intuitive() -> None:
    """Individual tests expose percentile, FDR and leave-one-out semantics."""

    row = {
        "background_profile_id": "HOUSEKEEPING_ALL",
        "comparison_scale": "MATCHED_RESIDUAL",
        "metric": "mean_distance",
        "observed_value": 0.4,
        "background_group_count": 12,
        "background_median": -0.1,
        "difference_from_background_median": 0.5,
        "empirical_percentile": 0.92,
        "two_sided_p_value": 0.04,
        "fdr_q_value": 0.08,
        "leave_one_out": True,
        "status": "TESTED",
    }
    figure = benchmark_individual_figure(
        rows=(row,),
        metric_label="Average pair distance",
        scale_label="Matched residual",
    )
    assert len(figure.data) == 2
    assert len(figure.layout.shapes) == 1
    displayed = _display_individual(row=row)
    assert displayed["BH-FDR q"] == pytest.approx(0.08)
    assert displayed["Leave one out"] is True

    profile = _display_profile(
        row={
            "profile_id": "E3_ALL",
            "profile_class": "E3",
            "profile_subclass": "ALL",
            "membership_role": "TARGET",
            "group_count": 10,
        }
    )
    assert profile["Clusters"] == 10
    contrast_row = {
        **row,
        "target_profile_id": "E3_ALL",
        "reference_profile_id": "HOUSEKEEPING_ALL",
        "target_group_count": 10,
        "reference_group_count": 12,
        "median_difference": 0.5,
        "median_difference_ci_low": 0.2,
        "median_difference_ci_high": 0.8,
        "cliffs_delta": 0.7,
        "p_value_two_sided": 0.01,
    }
    assert _display_contrast(row=contrast_row)["Metric"] == "Average pair distance"
    assert _metric_from_label(label="Average pair distance") == "mean_distance"
    with pytest.raises(ValueError, match="Unknown benchmark"):
        _metric_from_label(label="wrong")
    assert _cluster_label(
        row={
            "group_type": "HOG",
            "hierarchy_node": "N0",
            "group_id": "N0.HOG1",
            "profile_ids": "E3_ALL",
        }
    ) == "HOG | N0 | N0.HOG1 — E3_ALL"

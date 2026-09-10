"""Tests for calibrated dispersion page helpers and figures."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

import orthofinder_interrogation_app.benchmark_page as benchmark_page
from orthofinder_interrogation_app.benchmark_labels import (
    benchmark_class_label,
    benchmark_profile_label,
)
from orthofinder_interrogation_app.benchmark_page import (
    _cluster_label,
    _contrast_summary,
    _display_classification,
    _display_contrast,
    _display_individual,
    _display_marker,
    _display_member,
    _display_profile,
    _format_number,
    _group_key,
    _headline_contrasts,
    _marker_match_label,
    _metric_from_label,
    _open_group_in_explorer,
    _safe_file_token,
)
from orthofinder_interrogation_app.figures import (
    benchmark_classification_figure,
    benchmark_contrast_figure,
    benchmark_distribution_figure,
    benchmark_individual_figure,
    benchmark_marker_coverage_figure,
)
from orthofinder_results.errors import InputValidationError


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
    ) == "HOG | N0 | N0.HOG1 — All E3 clusters"


def test_benchmark_labels_and_plain_language_findings_preserve_authority() -> None:
    """Readable labels and findings retain exact IDs and cautious inference."""

    assert benchmark_profile_label(profile_id="HOUSEKEEPING_ALL").startswith("All")
    assert benchmark_profile_label(
        profile_id="R_NLR_SUBCLASS::TNL"
    ) == "R/NLR subclass — TNL"
    assert benchmark_profile_label(profile_id="NEW_CODE") == "NEW CODE"
    assert benchmark_profile_label(profile_id="") == "Unlabelled profile"
    assert benchmark_class_label(profile_classes="") == "Unlabelled profile"
    assert benchmark_class_label(profile_classes="E3;R_NLR") == "Mixed: E3 + R/NLR reference"
    row = {
        "target_profile_id": "E3_ALL",
        "reference_profile_id": "HOUSEKEEPING_ALL",
        "status": "TESTED",
        "median_difference": 0.5,
        "median_difference_ci_low": 0.2,
        "median_difference_ci_high": 0.8,
        "cliffs_delta": 0.7,
        "fdr_q_value": 0.02,
    }
    message, significant = _contrast_summary(
        row=row,
        metric_label="Average pair distance",
    )
    assert significant is True
    assert "more dispersed" in message and "BH-FDR q 0.02" in message
    assert _headline_contrasts(rows=(row,)) == (row,)
    unavailable = {**row, "status": "INSUFFICIENT_GROUPS", "median_difference": None}
    message, significant = _contrast_summary(
        row=unavailable,
        metric_label="Average pair distance",
    )
    assert significant is False and "insufficient" in message
    non_significant = {
        **row,
        "median_difference": -0.2,
        "median_difference_ci_low": None,
        "median_difference_ci_high": None,
        "fdr_q_value": 0.5,
    }
    message, significant = _contrast_summary(
        row=non_significant,
        metric_label="Average pair distance",
    )
    assert significant is False and "more compact" in message


def test_marker_membership_and_classification_views_are_explicit() -> None:
    """Marker provenance, complete membership and matched-control classes stay distinct."""

    marker = _display_marker(
        row={
            "profile_id": "HOUSEKEEPING_ALL",
            "profile_class": "HOUSEKEEPING",
            "profile_subclass": "ALL",
            "marker_id": "AT1G01010",
            "protein_identifier": "QHK1",
            "protein_entry": "HK1_ARATH",
            "marker_name": "Fixture housekeeping protein",
            "group_id": "N0.HOG3",
            "matched_member_id": "a3",
            "matched_species_label": "Species_A",
            "domain_architecture": "",
            "evidence_type": "reviewed",
            "source_title": "Fixture source",
            "source_doi": "10.0000/fixture",
            "source_version": "2026-09",
            "authority_name": "benchmarks.tsv",
        }
    )
    assert marker["Marker gene / locus"] == "AT1G01010"
    assert marker["OrthoFinder group"] == "N0.HOG3"
    marker_figure = benchmark_marker_coverage_figure(
        rows=(
            {
                "marker_id": "AT1G01010",
                "marker_name": "Fixture housekeeping protein",
                "group_type": "HOG",
                "hierarchy_node": "N0",
                "group_id": "N0.HOG3",
                "matched_member_id": "a3",
                "matched_species_label": "Species_A",
            },
        )
    )
    assert marker_figure.data[0].x[0] == 1
    with pytest.raises(InputValidationError, match="must be an integer"):
        benchmark_marker_coverage_figure(rows=(), maximum=True)
    with pytest.raises(InputValidationError, match="between 1 and 100"):
        benchmark_marker_coverage_figure(rows=(), maximum=101)
    member = _display_member(
        row={
            "species_label": "Species_A",
            "member_id": "a3",
            "legacy_orthogroup_id": "OG3",
            "gene_tree_parent_clade": "N0",
        },
        marker_members={"a3"},
    )
    assert member["Defines selected profile"] is True
    classification = {
        "group_type": "HOG",
        "hierarchy_node": "N0",
        "group_id": "N0.HOG1",
        "profile_ids": "E3_ALL",
        "profile_classes": "E3",
        "mean_distance": 0.2,
        "mean_distance_control_percentile": 0.9,
        "central_divergence_class": "DISPERSED",
        "distance_sd": 0.05,
        "distance_sd_control_percentile": 0.1,
        "heterogeneity_class": "LOW",
        "matched_control_count": 3,
        "status": "CLASSIFIED",
    }
    displayed = _display_classification(row=classification)
    assert displayed["Biological class"] == "E3"
    figure = benchmark_classification_figure(rows=(classification,))
    assert len(figure.data) == 1
    assert len(figure.layout.shapes) == 4
    assert _safe_file_token(value="HOUSEKEEPING::ALL / unsafe") == (
        "HOUSEKEEPING_ALL_unsafe"
    )
    assert _safe_file_token(value=" /// ") == "benchmark_profile"
    assert _format_number(value=None) == "Unavailable"
    assert _format_number(value=1.23456) == "1.235"


def test_cluster_labels_keys_and_navigation_are_collision_safe(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Cluster selectors and callbacks retain the complete composite identity."""

    row = {
        "group_type": "LEGACY_ORTHOGROUP",
        "hierarchy_node": "",
        "group_id": "OG1",
        "profile_ids": "E3_ALL;R_NLR_ALL",
    }
    assert "ROOT" in _cluster_label(row=row)
    assert _marker_match_label(
        row={"marker_id": "AT1G1", "group_id": "OG1", "matched_member_id": "p1"}
    ) == "AT1G1 | OG1 | p1"
    key = _group_key(row=row, run_id="run")
    captured = []
    monkeypatch.setattr(
        benchmark_page,
        "_store_active_group",
        lambda *, key: captured.append(key),
    )
    monkeypatch.setattr(
        benchmark_page,
        "st",
        SimpleNamespace(session_state={}),
    )
    _open_group_in_explorer(key=key)
    assert captured == [key]
    assert benchmark_page.st.session_state["app_page"] == "Cluster explorer"

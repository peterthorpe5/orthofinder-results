"""Tests for exact dispersion summaries, comparison figures and force views."""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import pytest

from orthofinder_interrogation_app import dispersion, figures, force_view
from orthofinder_interrogation_app.distance_data import DistanceAnalysisProvider
from orthofinder_interrogation_app.models import GroupKey
from orthofinder_interrogation_app.queries import OrthoFinderQueryService
from orthofinder_interrogation_app.report_data import load_visualisation_catalog
from orthofinder_interrogation_app.resource import open_resource
from orthofinder_results.errors import InputValidationError


@pytest.fixture
def exact_analysis(application_resource: Path, persistent_test_root: Path):
    """Recover the fixture's exact report-backed analysis."""

    resource = open_resource(path=application_resource)
    service = OrthoFinderQueryService(resource=resource)
    assert resource.report_path is not None
    catalog = load_visualisation_catalog(
        report_path=resource.report_path,
        expected_run_id=resource.run_id,
    )
    return DistanceAnalysisProvider(
        service=service,
        cache_dir=persistent_test_root / "dispersion_cache",
        report_catalog=catalog,
    ).analyse(key=GroupKey("test_run", "HOG", "N0", "N0.HOG1"), max_members=3)


def test_exact_geometry_member_and_species_dispersion(exact_analysis) -> None:
    """Every dispersion authority derives from the same complete exact matrix."""

    member_ids, species, matrix = dispersion.exact_distance_matrix(
        rows=exact_analysis.distances,
        members=exact_analysis.members,
    )
    assert member_ids == ("alpha_1", "alpha_2", "beta_1")
    assert species == ("Species_A", "Species_A", "Species_B")
    assert matrix.tolist() == [[0.0, 0.1, 0.2], [0.1, 0.0, 0.3], [0.2, 0.3, 0.0]]
    geometry = dispersion.classical_pcoa(
        rows=exact_analysis.distances,
        members=exact_analysis.members,
    )
    assert geometry.member_ids == member_ids
    assert len(geometry.coordinates) == 3
    assert geometry.axis_fractions[0] > 0
    assert 0 <= geometry.stress_3d <= geometry.stress_2d + 1e-12
    assert geometry.correlation_2d is not None
    member_rows = dispersion.member_dispersion_rows(
        rows=exact_analysis.distances,
        members=exact_analysis.members,
    )
    assert member_rows[0]["member_id"] == "alpha_1"
    assert member_rows[0]["is_sample_medoid"] is True
    assert member_rows[-1]["nearest_member_id"] in member_ids
    classes = dispersion.distance_class_rows(rows=exact_analysis.distances)
    assert [row["pair_class"] for row in classes].count("Within species") == 1
    species_rows = dispersion.species_dispersion_rows(rows=exact_analysis.distances)
    assert {(row["species_a"], row["species_b"]) for row in species_rows} == {
        ("Species_A", "Species_A"),
        ("Species_A", "Species_B"),
    }
    between = next(row for row in species_rows if row["pair_class"] == "Between species")
    assert between["pair_count"] == 2
    assert between["mean_distance"] == pytest.approx(0.25)


@pytest.mark.parametrize(
    ("members", "rows", "message"),
    [
        (({"member_id": "", "species_label": "A"},), (), "identifier and species"),
        (
            (
                {"member_id": "a", "species_label": "A"},
                {"member_id": "a", "species_label": "B"},
            ),
            (),
            "multiple species",
        ),
        (({"member_id": "a", "species_label": "A"},), (), "At least two"),
        (
            (
                {"member_id": "a", "species_label": "A"},
                {"member_id": "b", "species_label": "B"},
            ),
            ({"member_a": "a", "member_b": "a", "distance": 1},),
            "invalid displayed endpoint",
        ),
        (
            (
                {"member_id": "a", "species_label": "A"},
                {"member_id": "b", "species_label": "B"},
            ),
            ({"member_a": "a", "member_b": "b", "distance": "bad"},),
            "not numeric",
        ),
        (
            (
                {"member_id": "a", "species_label": "A"},
                {"member_id": "b", "species_label": "B"},
            ),
            ({"member_a": "a", "member_b": "b", "distance": -1},),
            "non-finite or negative",
        ),
        (
            (
                {"member_id": "a", "species_label": "A"},
                {"member_id": "b", "species_label": "B"},
            ),
            (
                {"member_a": "a", "member_b": "b", "distance": 1},
                {"member_a": "b", "member_b": "a", "distance": 1},
            ),
            "Duplicate",
        ),
        (
            (
                {"member_id": "a", "species_label": "A"},
                {"member_id": "b", "species_label": "B"},
                {"member_id": "c", "species_label": "C"},
            ),
            ({"member_a": "a", "member_b": "b", "distance": 1},),
            "requires 3 pairs",
        ),
    ],
)
def test_exact_matrix_rejects_incomplete_or_invalid_authorities(
    members, rows, message: str
) -> None:
    """Bad identities, distances, duplicates and incomplete pair sets fail visibly."""

    with pytest.raises(InputValidationError, match=message):
        dispersion.exact_distance_matrix(rows=rows, members=members)


def test_pcoa_and_classification_edge_cases(
    exact_analysis,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Degenerate geometry and missing endpoint taxonomy have controlled outcomes."""

    two_members = (
        {"member_id": "a", "species_label": "A"},
        {"member_id": "b", "species_label": "B"},
    )
    zero_row = ({"member_a": "a", "member_b": "b", "distance": 0.0},)
    with pytest.raises(InputValidationError, match="no positive"):
        dispersion.classical_pcoa(rows=zero_row, members=two_members)
    monkeypatch.setattr(
        dispersion.np.linalg,
        "eigh",
        lambda value: (_ for _ in ()).throw(np.linalg.LinAlgError("fixture")),
    )
    with pytest.raises(InputValidationError, match="eigendecomposition"):
        dispersion.classical_pcoa(
            rows=exact_analysis.distances,
            members=exact_analysis.members,
        )
    assert dispersion._normalised_stress(
        exact=np.asarray([0.0]),
        projected=np.asarray([1.0]),
    ) == 0.0
    assert dispersion._distance_correlation(
        exact=np.asarray([1.0]),
        projected=np.asarray([1.0]),
    ) is None
    assert dispersion._distance_correlation(
        exact=np.asarray([1.0, 1.0]),
        projected=np.asarray([1.0, 2.0]),
    ) is None
    with pytest.raises(InputValidationError, match="endpoint species"):
        dispersion.distance_class_rows(
            rows=({"member_a": "a", "member_b": "b", "distance": 1},)
        )
    with pytest.raises(InputValidationError, match="endpoint species"):
        dispersion.species_dispersion_rows(
            rows=({"member_a": "a", "member_b": "b", "distance": 1},)
        )


def test_all_new_dispersion_and_comparison_figures(exact_analysis) -> None:
    """The enhanced single- and multi-group visual repertoire renders real traces."""

    geometry = dispersion.classical_pcoa(
        rows=exact_analysis.distances,
        members=exact_analysis.members,
    )
    selected = ("alpha_1",)
    three_dimensional = figures.pcoa_3d_figure(
        geometry=geometry,
        selected_members=selected,
    )
    assert len(three_dimensional.data) == 2
    assert three_dimensional.layout.scene.zaxis.title.text.startswith("Axis 3")
    axes = figures.pcoa_axis_figure(
        geometry=geometry,
        horizontal_axis=1,
        vertical_axis=3,
        selected_members=selected,
    )
    assert len(axes.data) == 2
    assert "axes 1 and 3" in axes.layout.title.text
    with pytest.raises(InputValidationError, match="between one and three"):
        figures.pcoa_axis_figure(
            geometry=geometry,
            horizontal_axis=0,
            vertical_axis=2,
        )
    with pytest.raises(InputValidationError, match="must differ"):
        figures.pcoa_axis_figure(
            geometry=geometry,
            horizontal_axis=2,
            vertical_axis=2,
        )

    member_rows = dispersion.member_dispersion_rows(
        rows=exact_analysis.distances,
        members=exact_analysis.members,
    )
    assert len(figures.distance_distribution_figure(rows=exact_analysis.distances).data) == 6
    assert len(
        figures.member_dispersion_figure(rows=member_rows, selected_members=selected).data
    ) == 2
    assert len(
        figures.medoid_distance_figure(
            rows=exact_analysis.distances,
            member_rows=member_rows,
        ).data
    ) == 1
    assert len(figures.species_pair_heatmap_figure(rows=exact_analysis.distances).data) == 1

    summary = {
        "label": exact_analysis.key.display_label(),
        "mean_distance": 0.2,
        "median_distance": 0.2,
        "population_stddev_distance": 0.08,
    }
    assert len(figures.comparison_summary_figure(summaries=(summary,)).data) == 2
    groups = {"first": (0.1, 0.2), "second": (0.3, 0.4)}
    assert len(
        figures.comparison_distribution_figure(
            distance_groups=groups,
            mode="VIOLIN",
        ).data
    ) == 2
    assert len(
        figures.comparison_distribution_figure(
            distance_groups=groups,
            mode="ECDF",
        ).data
    ) == 2
    with pytest.raises(InputValidationError, match="Unsupported"):
        figures.comparison_distribution_figure(distance_groups=groups, mode="BOX")
    comparisons = {"group one": geometry, "group two": geometry}
    assert len(figures.comparison_pcoa_figure(geometries=comparisons).data) == 2
    with pytest.raises(InputValidationError, match="between 2 and 12"):
        figures.comparison_pcoa_figure(geometries={"only": geometry})
    assert figures._stable_species_colour(species="Species_A").startswith("hsl(")
    assert figures._series_colour(index=12) == figures._series_colour(index=0)


def test_medoid_figure_requires_one_medoid(exact_analysis) -> None:
    """A medoid distance profile cannot silently invent or choose a centre."""

    with pytest.raises(InputValidationError, match="Exactly one"):
        figures.medoid_distance_figure(
            rows=exact_analysis.distances,
            member_rows=exact_analysis.members,
        )
    two = [dict(row, is_sample_medoid=True) for row in exact_analysis.members[:2]]
    with pytest.raises(InputValidationError, match="Exactly one"):
        figures.medoid_distance_figure(
            rows=exact_analysis.distances,
            member_rows=two,
        )


def test_force_directed_view_is_inline_draggable_and_filterable(exact_analysis) -> None:
    """The restored force view is self-contained and exposes explicit layout controls."""

    entry = exact_analysis.visual_entry
    without_connectors = force_view.force_directed_html(
        entry=entry,
        selected_members=("alpha_1",),
        include_connectors=False,
        show_labels=True,
        physics_enabled=False,
    )
    assert "vis.Network" in without_connectors
    assert "dragNodes" in without_connectors
    assert '"enabled": false' in without_connectors
    with_connectors = force_view.force_directed_html(
        entry=entry,
        include_connectors=True,
        physics_enabled=True,
    )
    assert "COMPONENT_CONNECTOR" in with_connectors
    assert len(with_connectors) >= len(without_connectors) - 2_000
    assert force_view._node_colour(
        raw_node={"color": {"background": "#abcdef"}},
        selected=True,
        medoid=False,
    )["background"] == "#abcdef"
    assert '"enabled": true' in force_view._network_options(physics_enabled=True)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda entry: entry.update(nodes="wrong"), "node or edge arrays"),
        (lambda entry: entry.update(nodes=[]), "requires 1"),
        (lambda entry: entry["nodes"].append("wrong"), "malformed node"),
        (
            lambda entry: entry["nodes"][1].update(id=entry["nodes"][0]["id"]),
            "unique identifiers",
        ),
        (lambda entry: entry["edges"].append("wrong"), "malformed edge"),
        (
            lambda entry: entry["edges"].append(
                {"from": "missing", "to": "alpha_1", "edgeType": "NEAREST_NEIGHBOUR"}
            ),
            "invalid endpoint",
        ),
        (lambda entry: entry["nodes"][0].update(title="</script>"), "unsafe text"),
    ],
)
def test_force_view_rejects_malformed_or_unsafe_graphs(
    exact_analysis,
    mutation,
    message: str,
) -> None:
    """Graph embedding rejects malformed identities, endpoints and HTML delimiters."""

    entry = copy.deepcopy(exact_analysis.visual_entry)
    mutation(entry)
    with pytest.raises(InputValidationError, match=message):
        force_view.force_directed_html(entry=entry)


def test_force_edge_bound_and_safe_text_helpers(
    exact_analysis,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Large edge sets and control characters are rejected before HTML generation."""

    monkeypatch.setattr(force_view, "MAX_FORCE_EDGES", 1)
    with pytest.raises(InputValidationError, match="edge safety"):
        force_view.force_directed_html(
            entry=exact_analysis.visual_entry,
            include_connectors=True,
        )
    with pytest.raises(InputValidationError, match="unsafe text"):
        force_view._safe_graph_text(value="bad\x00text", role="fixture")

"""Plotly figures for bounded evolutionary visualisation records."""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import networkx as nx
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from orthofinder_results.errors import InputValidationError

from .benchmark_labels import benchmark_class_label, benchmark_profile_label
from .dispersion import PcoaGeometry, distance_class_rows, species_dispersion_rows

_DEFAULT_COLOUR = "#5b6475"
_HIGHLIGHT_BORDER = "#111827"


def pcoa_3d_figure(
    *, geometry: PcoaGeometry, selected_members: Iterable[str] = ()
) -> go.Figure:
    """Build an interactive three-axis PCoA diagnostic.

    Args:
        geometry: Exact-matrix PCoA coordinates and diagnostics.
        selected_members: Linked canonical members to emphasise.

    Returns:
        Species-coloured rotatable three-dimensional Plotly figure.
    """

    selected = {str(value) for value in selected_members}
    by_species: dict[str, list[int]] = defaultdict(list)
    for index, species in enumerate(geometry.species_labels):
        by_species[species].append(index)
    figure = go.Figure()
    for species in sorted(by_species):
        indices = by_species[species]
        ids = [geometry.member_ids[index] for index in indices]
        chosen = [not selected or member_id in selected for member_id in ids]
        figure.add_trace(
            go.Scatter3d(
                x=[geometry.coordinates[index][0] for index in indices],
                y=[geometry.coordinates[index][1] for index in indices],
                z=[geometry.coordinates[index][2] for index in indices],
                mode="markers",
                name=species,
                customdata=[[member_id, species] for member_id in ids],
                marker={
                    "color": _stable_species_colour(species=species),
                    "size": [8 if active else 4 for active in chosen],
                    "opacity": 0.9 if not selected else 0.75,
                    "line": {
                        "color": [
                            _HIGHLIGHT_BORDER if active and selected else "#ffffff"
                            for active in chosen
                        ],
                        "width": 2 if selected else 0.3,
                    },
                },
                hovertemplate=(
                    "Member=%{customdata[0]}<br>Species=%{customdata[1]}<br>"
                    "Axis 1=%{x:.6g}<br>Axis 2=%{y:.6g}<br>Axis 3=%{z:.6g}"
                    "<extra></extra>"
                ),
            )
        )
    fractions = geometry.axis_fractions
    figure.update_layout(
        title="Three-axis complete-distance PCoA diagnostic",
        scene={
            "xaxis_title": f"Axis 1 ({100 * fractions[0]:.1f}%)",
            "yaxis_title": f"Axis 2 ({100 * fractions[1]:.1f}%)",
            "zaxis_title": f"Axis 3 ({100 * fractions[2]:.1f}%)",
            "aspectmode": "data",
        },
        height=760,
        template="plotly_white",
        legend_title="Species",
    )
    return figure


def pcoa_axis_figure(
    *,
    geometry: PcoaGeometry,
    horizontal_axis: int,
    vertical_axis: int,
    selected_members: Iterable[str] = (),
) -> go.Figure:
    """Build a selectable two-axis view of the three-axis PCoA solution.

    Args:
        geometry: Exact-matrix PCoA coordinates and diagnostics.
        horizontal_axis: One-based horizontal axis in the range one to three.
        vertical_axis: One-based vertical axis in the range one to three.
        selected_members: Linked canonical members to emphasise.

    Returns:
        Interactive two-dimensional diagnostic scatter.

    Raises:
        InputValidationError: If axes are invalid or identical.
    """

    if horizontal_axis not in {1, 2, 3} or vertical_axis not in {1, 2, 3}:
        raise InputValidationError("PCoA axes must be between one and three.")
    if horizontal_axis == vertical_axis:
        raise InputValidationError("Horizontal and vertical PCoA axes must differ.")
    selected = {str(value) for value in selected_members}
    by_species: dict[str, list[int]] = defaultdict(list)
    for index, species in enumerate(geometry.species_labels):
        by_species[species].append(index)
    figure = go.Figure()
    x_index, y_index = horizontal_axis - 1, vertical_axis - 1
    for species in sorted(by_species):
        indices = by_species[species]
        ids = [geometry.member_ids[index] for index in indices]
        active = [not selected or member_id in selected for member_id in ids]
        figure.add_trace(
            go.Scattergl(
                x=[geometry.coordinates[index][x_index] for index in indices],
                y=[geometry.coordinates[index][y_index] for index in indices],
                mode="markers",
                name=species,
                customdata=[[member_id, species] for member_id in ids],
                marker={
                    "color": _stable_species_colour(species=species),
                    "size": [13 if chosen and selected else 8 for chosen in active],
                    "opacity": [1.0 if chosen else 0.12 for chosen in active],
                    "line": {
                        "color": [
                            _HIGHLIGHT_BORDER if chosen and selected else "#ffffff"
                            for chosen in active
                        ],
                        "width": [2 if chosen and selected else 0.5 for chosen in active],
                    },
                },
                hovertemplate=(
                    "Member=%{customdata[0]}<br>Species=%{customdata[1]}<br>"
                    f"Axis {horizontal_axis}=%{{x:.6g}}<br>"
                    f"Axis {vertical_axis}=%{{y:.6g}}<extra></extra>"
                ),
            )
        )
    fractions = geometry.axis_fractions
    figure.update_layout(
        title=f"Complete-distance PCoA: axes {horizontal_axis} and {vertical_axis}",
        xaxis_title=(
            f"PCoA axis {horizontal_axis} ({100 * fractions[x_index]:.1f}% positive inertia)"
        ),
        yaxis_title=(
            f"PCoA axis {vertical_axis} ({100 * fractions[y_index]:.1f}% positive inertia)"
        ),
        height=680,
        template="plotly_white",
        dragmode="lasso",
        legend_title="Species",
    )
    figure.update_yaxes(scaleanchor="x", scaleratio=1)
    return figure


def distance_distribution_figure(*, rows: Sequence[Mapping[str, Any]]) -> go.Figure:
    """Build histogram, violin and ECDF views of exact pair distances.

    Args:
        rows: Species-decorated pairwise distance records.

    Returns:
        Three-panel distribution figure separating within- and between-species pairs.
    """

    classified = distance_class_rows(rows=rows)
    figure = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=("Histogram", "Violin + box", "Empirical CDF"),
    )
    colours = {"Within species": "#2c7fb8", "Between species": "#d95f0e"}
    for pair_class in ("Within species", "Between species"):
        values = [
            float(row["distance"])
            for row in classified
            if row["pair_class"] == pair_class
        ]
        if not values:
            continue
        colour = colours[pair_class]
        figure.add_trace(
            go.Histogram(
                x=values,
                name=pair_class,
                marker_color=colour,
                opacity=0.60,
                histnorm="probability",
                legendgroup=pair_class,
            ),
            row=1,
            col=1,
        )
        figure.add_trace(
            go.Violin(
                x=[pair_class] * len(values),
                y=values,
                name=pair_class,
                box_visible=True,
                meanline_visible=True,
                line_color=colour,
                fillcolor=colour,
                opacity=0.60,
                legendgroup=pair_class,
                showlegend=False,
            ),
            row=1,
            col=2,
        )
        ordered = sorted(values)
        figure.add_trace(
            go.Scattergl(
                x=ordered,
                y=[(index + 1) / len(ordered) for index in range(len(ordered))],
                mode="lines",
                name=pair_class,
                line={"color": colour, "width": 2.5},
                legendgroup=pair_class,
                showlegend=False,
            ),
            row=1,
            col=3,
        )
    figure.update_xaxes(title_text="Exact pairwise distance", row=1, col=1)
    figure.update_yaxes(title_text="Fraction of pairs", row=1, col=1)
    figure.update_yaxes(title_text="Exact pairwise distance", row=1, col=2)
    figure.update_xaxes(title_text="Exact pairwise distance", row=1, col=3)
    figure.update_yaxes(title_text="Cumulative fraction", range=[0, 1], row=1, col=3)
    figure.update_layout(
        title="Within-group distance dispersion",
        barmode="overlay",
        height=560,
        template="plotly_white",
        legend_title="Endpoint class",
    )
    return figure


def member_dispersion_figure(
    *, rows: Sequence[Mapping[str, Any]], selected_members: Iterable[str] = ()
) -> go.Figure:
    """Plot each member's mean distance and variability to all peers.

    Args:
        rows: Per-member dispersion summary rows.
        selected_members: Linked members to emphasise.

    Returns:
        Ranked mean-distance scatter with population-SD error bars.
    """

    selected = {str(value) for value in selected_members}
    ordered = sorted(
        rows,
        key=lambda row: (float(row["mean_distance"]), str(row["member_id"])),
    )
    figure = go.Figure()
    for species in sorted({str(row["species_label"]) for row in ordered}):
        indices = [
            index for index, row in enumerate(ordered) if str(row["species_label"]) == species
        ]
        species_rows = [ordered[index] for index in indices]
        figure.add_trace(
            go.Scattergl(
                x=[index + 1 for index in indices],
                y=[float(row["mean_distance"]) for row in species_rows],
                mode="markers",
                name=species,
                customdata=[
                    [
                        row["member_id"],
                        row["nearest_member_id"],
                        row["nearest_distance"],
                        row["maximum_distance"],
                    ]
                    for row in species_rows
                ],
                error_y={
                    "type": "data",
                    "array": [
                        float(row["population_stddev_distance"]) for row in species_rows
                    ],
                    "visible": True,
                    "thickness": 0.8,
                },
                marker={
                    "color": _stable_species_colour(species=species),
                    "size": [
                        15
                        if bool(row.get("is_sample_medoid"))
                        or str(row["member_id"]) in selected
                        else 8
                        for row in species_rows
                    ],
                    "symbol": [
                        "star" if bool(row.get("is_sample_medoid")) else "circle"
                        for row in species_rows
                    ],
                    "line": {"color": "#ffffff", "width": 0.5},
                },
                hovertemplate=(
                    "Member=%{customdata[0]}<br>Mean distance=%{y:.6g}<br>"
                    "Nearest=%{customdata[1]} (%{customdata[2]:.6g})<br>"
                    "Maximum=%{customdata[3]:.6g}<extra></extra>"
                ),
            )
        )
    figure.update_layout(
        title="Member centrality and peripherality",
        xaxis_title="Member rank from compact centre to periphery",
        yaxis_title="Mean distance to every other displayed member (± population SD)",
        height=600,
        template="plotly_white",
        legend_title="Species",
    )
    return figure


def medoid_distance_figure(
    *,
    rows: Sequence[Mapping[str, Any]],
    member_rows: Sequence[Mapping[str, Any]],
) -> go.Figure:
    """Plot every displayed member's exact distance from the sample medoid.

    Args:
        rows: Species-decorated exact pairwise distances.
        member_rows: Per-member summaries containing the sample-medoid flag.

    Returns:
        Sorted sample-medoid distance profile.
    """

    medoids = [str(row["member_id"]) for row in member_rows if row.get("is_sample_medoid")]
    if len(medoids) != 1:
        raise InputValidationError("Exactly one sample medoid is required.")
    medoid = medoids[0]
    species_by_member = {
        str(row["member_id"]): str(row["species_label"]) for row in member_rows
    }
    values = [(medoid, 0.0)]
    for row in rows:
        left, right = str(row["member_a"]), str(row["member_b"])
        if left == medoid:
            values.append((right, float(row["distance"])))
        elif right == medoid:
            values.append((left, float(row["distance"])))
    values.sort(key=lambda item: (item[1], item[0]))
    figure = go.Figure(
        go.Bar(
            x=[member for member, _ in values],
            y=[distance for _, distance in values],
            marker_color=[
                _stable_species_colour(species=species_by_member[member])
                for member, _ in values
            ],
            customdata=[[species_by_member[member]] for member, _ in values],
            hovertemplate=(
                "Member=%{x}<br>Species=%{customdata[0]}<br>"
                "Distance from sample medoid=%{y:.6g}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title=f"Exact distance from sample medoid {medoid}",
        xaxis_title="Displayed members, sorted by medoid distance",
        yaxis_title="Exact pairwise distance",
        height=560,
        template="plotly_white",
    )
    figure.update_xaxes(showticklabels=len(values) <= 80)
    return figure


def species_pair_heatmap_figure(*, rows: Sequence[Mapping[str, Any]]) -> go.Figure:
    """Build a mean-distance heatmap across represented species pairs.

    Args:
        rows: Species-decorated pairwise distance records.

    Returns:
        Symmetric species-pair mean-distance heatmap with missing cells explicit.
    """

    summaries = species_dispersion_rows(rows=rows)
    species = sorted(
        {str(row[field]) for row in summaries for field in ("species_a", "species_b")}
    )
    index = {label: position for position, label in enumerate(species)}
    matrix = np.full((len(species), len(species)), np.nan, dtype=float)
    counts = np.zeros((len(species), len(species)), dtype=int)
    for row in summaries:
        left, right = index[str(row["species_a"])], index[str(row["species_b"])]
        matrix[left, right] = matrix[right, left] = float(row["mean_distance"])
        counts[left, right] = counts[right, left] = int(row["pair_count"])
    figure = go.Figure(
        go.Heatmap(
            z=matrix,
            x=species,
            y=species,
            customdata=counts,
            colorscale="Viridis",
            colorbar={"title": "Mean distance"},
            hovertemplate=(
                "Species A=%{y}<br>Species B=%{x}<br>Mean distance=%{z:.6g}<br>"
                "Pairs=%{customdata:,}<extra></extra>"
            ),
            hoverongaps=False,
        )
    )
    figure.update_layout(
        title="Mean pairwise distance by represented species pair",
        xaxis_title="Species",
        yaxis_title="Species",
        height=max(620, min(1_300, 240 + 14 * len(species))),
        template="plotly_white",
    )
    figure.update_yaxes(autorange="reversed")
    return figure


def comparison_summary_figure(*, summaries: Sequence[Mapping[str, Any]]) -> go.Figure:
    """Compare group mean distance, SD and median on a shared scale.

    Args:
        summaries: Labelled distance-summary records.

    Returns:
        Horizontal mean-distance plot with population-SD error bars and medians.
    """

    ordered = sorted(
        summaries,
        key=lambda row: (float(row["mean_distance"]), str(row["label"])),
    )
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=[float(row["mean_distance"]) for row in ordered],
            y=[str(row["label"]) for row in ordered],
            mode="markers",
            name="Mean ± population SD",
            error_x={
                "type": "data",
                "array": [float(row["population_stddev_distance"]) for row in ordered],
                "visible": True,
            },
            marker={"size": 11, "color": "#2c7fb8"},
        )
    )
    figure.add_trace(
        go.Scatter(
            x=[float(row["median_distance"]) for row in ordered],
            y=[str(row["label"]) for row in ordered],
            mode="markers",
            name="Median",
            marker={"size": 10, "symbol": "diamond", "color": "#d95f0e"},
        )
    )
    figure.update_layout(
        title="Between-group compactness comparison",
        xaxis_title="Exact displayed pairwise distance",
        yaxis_title="Group",
        height=max(500, 140 + 55 * len(ordered)),
        template="plotly_white",
    )
    return figure


def comparison_distribution_figure(
    *, distance_groups: Mapping[str, Sequence[float]], mode: str
) -> go.Figure:
    """Compare exact group distance distributions as violins or ECDFs.

    Args:
        distance_groups: Group label to complete displayed distance vector.
        mode: ``VIOLIN`` or ``ECDF``.

    Returns:
        Shared-axis comparison figure.

    Raises:
        InputValidationError: If the requested display mode is unsupported.
    """

    if mode not in {"VIOLIN", "ECDF"}:
        raise InputValidationError(f"Unsupported comparison distribution mode: {mode}")
    figure = go.Figure()
    for index, (label, raw_values) in enumerate(sorted(distance_groups.items())):
        values = sorted(float(value) for value in raw_values)
        colour = _series_colour(index=index)
        if mode == "VIOLIN":
            figure.add_trace(
                go.Violin(
                    x=[label] * len(values),
                    y=values,
                    name=label,
                    box_visible=True,
                    meanline_visible=True,
                    line_color=colour,
                    fillcolor=colour,
                    opacity=0.65,
                    showlegend=False,
                )
            )
        else:
            figure.add_trace(
                go.Scattergl(
                    x=values,
                    y=[(position + 1) / len(values) for position in range(len(values))],
                    mode="lines",
                    name=label,
                    line={"color": colour, "width": 2.2},
                )
            )
    if mode == "VIOLIN":
        x_title, y_title = "Group", "Exact displayed pairwise distance"
        title = "Exact distance distributions across groups"
    else:
        x_title, y_title = "Exact displayed pairwise distance", "Cumulative fraction"
        title = "Empirical distance distributions across groups"
    figure.update_layout(
        title=title,
        xaxis_title=x_title,
        yaxis_title=y_title,
        height=620,
        template="plotly_white",
        legend_title="Group",
    )
    return figure


def comparison_pcoa_figure(*, geometries: Mapping[str, PcoaGeometry]) -> go.Figure:
    """Build consistent-colour small-multiple PCoA panels for several groups.

    Each group is projected independently. Coordinates and orientation must not be
    compared between panels; the panels compare within-group dispersion patterns.

    Args:
        geometries: Two to twelve labelled group geometries.

    Returns:
        Small-multiple axes-one-and-two PCoA figure.

    Raises:
        InputValidationError: If the comparison size is outside two to twelve groups.
    """

    if not 2 <= len(geometries) <= 12:
        raise InputValidationError("PCoA comparison requires between 2 and 12 groups.")
    labels = sorted(geometries)
    columns = min(3, len(labels))
    rows = math.ceil(len(labels) / columns)
    figure = make_subplots(rows=rows, cols=columns, subplot_titles=labels)
    for panel, label in enumerate(labels):
        geometry = geometries[label]
        row, column = panel // columns + 1, panel % columns + 1
        figure.add_trace(
            go.Scattergl(
                x=[coordinate[0] for coordinate in geometry.coordinates],
                y=[coordinate[1] for coordinate in geometry.coordinates],
                mode="markers",
                name=label,
                customdata=[
                    [member, species]
                    for member, species in zip(
                        geometry.member_ids,
                        geometry.species_labels,
                        strict=True,
                    )
                ],
                marker={
                    "color": [
                        _stable_species_colour(species=species)
                        for species in geometry.species_labels
                    ],
                    "size": 6,
                    "opacity": 0.80,
                    "line": {"color": "#ffffff", "width": 0.3},
                },
                showlegend=False,
                hovertemplate=(
                    "Member=%{customdata[0]}<br>Species=%{customdata[1]}<br>"
                    "Axis 1=%{x:.6g}<br>Axis 2=%{y:.6g}<extra></extra>"
                ),
            ),
            row=row,
            col=column,
        )
        figure.update_xaxes(title_text="Axis 1", row=row, col=column)
        figure.update_yaxes(
            title_text="Axis 2",
            scaleanchor=f"x{panel + 1}" if panel else "x",
            scaleratio=1,
            row=row,
            col=column,
        )
    figure.update_layout(
        title="Independent within-group PCoA small multiples",
        height=max(620, 480 * rows),
        template="plotly_white",
    )
    return figure


def benchmark_distribution_figure(
    *, rows: Sequence[Mapping[str, Any]], metric_label: str, scale_label: str
) -> go.Figure:
    """Build cluster-level box-and-point distributions for benchmark profiles.

    Args:
        rows: Profile-labelled cluster observations.
        metric_label: Plain-language selected dispersion metric.
        scale_label: Plain-language raw or matched-residual scale.

    Returns:
        Interactive profile distributions whose points are cluster replicates.
    """

    figure = go.Figure()
    profiles = sorted(
        {str(row["profile_id"]) for row in rows},
        key=lambda value: benchmark_profile_label(profile_id=value),
    )
    for index, profile_id in enumerate(profiles):
        selected = [row for row in rows if str(row["profile_id"]) == profile_id]
        display_label = benchmark_profile_label(profile_id=profile_id)
        figure.add_trace(
            go.Box(
                x=[display_label] * len(selected),
                y=[float(row["metric_value"]) for row in selected],
                name=display_label,
                boxpoints="all",
                jitter=0.35,
                pointpos=0,
                marker={"size": 6, "opacity": 0.68},
                line={"color": _series_colour(index=index)},
                fillcolor=_series_colour(index=index),
                opacity=0.62,
                customdata=[
                    [
                        profile_id,
                        row["group_type"],
                        row["hierarchy_node"],
                        row["group_id"],
                    ]
                    for row in selected
                ],
                hovertemplate=(
                    "Profile=%{x}<br>Profile ID=%{customdata[0]}<br>"
                    "Value=%{y:.6g}<br>Group=%{customdata[1]} | "
                    "%{customdata[2]} | %{customdata[3]}<extra></extra>"
                ),
            )
        )
    figure.update_layout(
        title=f"Cluster-level {metric_label} by biological profile",
        xaxis_title="Biological profile",
        yaxis_title=f"{metric_label} ({scale_label})",
        height=650,
        template="plotly_white",
        showlegend=False,
    )
    return figure


def benchmark_marker_coverage_figure(
    *, rows: Sequence[Mapping[str, Any]], maximum: int = 40
) -> go.Figure:
    """Show how profile-defining markers map to clusters and proteins.

    Args:
        rows: Marker-to-cluster matches from the benchmark authority.
        maximum: Maximum markers displayed after deterministic ranking.

    Returns:
        Horizontal cluster-count bars with matched-protein hover details.

    Raises:
        InputValidationError: If the display bound is outside one to 100.
    """

    if not isinstance(maximum, int) or isinstance(maximum, bool):
        raise InputValidationError("maximum must be an integer.")
    if not 1 <= maximum <= 100:
        raise InputValidationError("maximum must be between 1 and 100.")
    markers: dict[str, dict[str, Any]] = {}
    for row in rows:
        marker_id = str(row["marker_id"])
        marker = markers.setdefault(
            marker_id,
            {
                "marker_id": marker_id,
                "marker_name": str(row.get("marker_name", "")),
                "clusters": set(),
                "members": set(),
                "species": set(),
            },
        )
        marker["clusters"].add(
            (
                str(row["group_type"]),
                str(row["hierarchy_node"]),
                str(row["group_id"]),
            )
        )
        marker["members"].add(str(row["matched_member_id"]))
        marker["species"].add(str(row["matched_species_label"]))
    ranked = sorted(
        markers.values(),
        key=lambda marker: (
            len(marker["clusters"]),
            len(marker["members"]),
            str(marker["marker_id"]),
        ),
        reverse=True,
    )[:maximum]
    ranked.reverse()
    labels = [
        (
            f"{marker['marker_id']} — {str(marker['marker_name'])[:55]}"
            if marker["marker_name"]
            else str(marker["marker_id"])
        )
        for marker in ranked
    ]
    figure = go.Figure(
        go.Bar(
            x=[len(marker["clusters"]) for marker in ranked],
            y=labels,
            orientation="h",
            marker_color="#2c7fb8",
            customdata=[
                [
                    marker["marker_id"],
                    len(marker["members"]),
                    len(marker["species"]),
                ]
                for marker in ranked
            ],
            hovertemplate=(
                "Marker=%{customdata[0]}<br>Matched clusters=%{x:,}"
                "<br>Matched proteins=%{customdata[1]:,}"
                "<br>Matched species=%{customdata[2]:,}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title=f"Marker-to-cluster coverage (up to {maximum} markers shown)",
        xaxis_title="Distinct matched OrthoFinder clusters",
        yaxis_title="Profile-defining marker",
        height=max(480, min(1_500, 180 + 28 * len(ranked))),
        template="plotly_white",
    )
    figure.update_xaxes(dtick=1)
    return figure


def benchmark_classification_figure(*, rows: Sequence[Mapping[str, Any]]) -> go.Figure:
    """Build an interactive map of each target against its own matched controls.

    Args:
        rows: Cluster classifications enriched with biological profile fields.

    Returns:
        Mean-distance and distance-spread percentiles for classified clusters.
    """

    classified = [
        row
        for row in rows
        if row.get("status") == "CLASSIFIED"
        and row.get("mean_distance_control_percentile") is not None
        and row.get("distance_sd_control_percentile") is not None
    ]
    figure = go.Figure()
    class_labels = sorted(
        {
            benchmark_class_label(profile_classes=str(row.get("profile_classes", "")))
            for row in classified
        }
    )
    for index, class_label in enumerate(class_labels):
        selected = [
            row
            for row in classified
            if benchmark_class_label(
                profile_classes=str(row.get("profile_classes", ""))
            )
            == class_label
        ]
        figure.add_trace(
            go.Scattergl(
                x=[float(row["mean_distance_control_percentile"]) for row in selected],
                y=[float(row["distance_sd_control_percentile"]) for row in selected],
                mode="markers",
                name=class_label,
                marker={
                    "size": 9,
                    "opacity": 0.72,
                    "color": _series_colour(index=index),
                    "line": {"color": "#ffffff", "width": 0.6},
                },
                customdata=[
                    [
                        row["group_type"],
                        row["hierarchy_node"],
                        row["group_id"],
                        row.get("profile_ids", ""),
                        row.get("mean_distance"),
                        row.get("distance_sd"),
                        row.get("central_divergence_class", ""),
                        row.get("heterogeneity_class", ""),
                    ]
                    for row in selected
                ],
                hovertemplate=(
                    "Group=%{customdata[0]} | %{customdata[1]} | %{customdata[2]}"
                    "<br>Profiles=%{customdata[3]}<br>Average pair distance="
                    "%{customdata[4]:.5g}<br>Distance SD=%{customdata[5]:.5g}"
                    "<br>Central divergence=%{customdata[6]}"
                    "<br>Heterogeneity=%{customdata[7]}<extra></extra>"
                ),
            )
        )
    for threshold in (0.1, 0.9):
        figure.add_vline(x=threshold, line_dash="dot", line_color="#9ca3af")
        figure.add_hline(y=threshold, line_dash="dot", line_color="#9ca3af")
    figure.update_layout(
        title="Each biological cluster relative to its own matched controls",
        xaxis_title="Average-distance percentile (left = compact; right = dispersed)",
        yaxis_title="Distance-spread percentile (low = uniform; high = heterogeneous)",
        height=700,
        template="plotly_white",
        legend_title="Biological class",
        hovermode="closest",
    )
    figure.update_xaxes(range=[-0.03, 1.03], tickformat=".0%")
    figure.update_yaxes(range=[-0.03, 1.03], tickformat=".0%")
    return figure


def benchmark_contrast_figure(
    *, rows: Sequence[Mapping[str, Any]], metric_label: str
) -> go.Figure:
    """Build a confidence-interval forest plot for tested profile contrasts.

    Args:
        rows: Precomputed contrast records for one metric.
        metric_label: Plain-language selected metric.

    Returns:
        Median matched-residual differences with 95% bootstrap intervals.
    """

    tested = [
        row
        for row in rows
        if row.get("status") == "TESTED"
        and row.get("median_difference") is not None
        and row.get("median_difference_ci_low") is not None
        and row.get("median_difference_ci_high") is not None
    ]
    tested.sort(key=lambda row: float(row["median_difference"]))
    labels = [
        f"{benchmark_profile_label(profile_id=str(row['target_profile_id']))} vs "
        f"{benchmark_profile_label(profile_id=str(row['reference_profile_id']))}"
        for row in tested
    ]
    differences = [float(row["median_difference"]) for row in tested]
    q_values = [
        None if row.get("fdr_q_value") is None else float(row["fdr_q_value"])
        for row in tested
    ]
    figure = go.Figure(
        go.Scatter(
            x=differences,
            y=labels,
            mode="markers",
            marker={
                "size": 10,
                "color": [
                    "#b91c1c" if value is not None and value <= 0.05 else "#355c8a"
                    for value in q_values
                ],
            },
            error_x={
                "type": "data",
                "symmetric": False,
                "array": [
                    float(row["median_difference_ci_high"]) - difference
                    for row, difference in zip(tested, differences, strict=True)
                ],
                "arrayminus": [
                    difference - float(row["median_difference_ci_low"])
                    for row, difference in zip(tested, differences, strict=True)
                ],
                "visible": True,
            },
            customdata=[
                [
                    row["target_profile_id"],
                    row["reference_profile_id"],
                    row["cliffs_delta"],
                    row["p_value_two_sided"],
                    row["fdr_q_value"],
                ]
                for row in tested
            ],
            hovertemplate=(
                "Target ID=%{customdata[0]}<br>Reference ID=%{customdata[1]}"
                "<br>Median difference=%{x:.6g}<br>Cliff's delta="
                "%{customdata[2]:.4g}<br>p=%{customdata[3]:.4g}"
                "<br>FDR q=%{customdata[4]:.4g}"
                "<extra></extra>"
            ),
        )
    )
    figure.add_vline(x=0.0, line_dash="dash", line_color="#6b7280")
    figure.update_layout(
        title=f"Matched-background contrasts: {metric_label}",
        xaxis_title=(
            "Median target minus reference matched-control residual "
            "(positive = greater dispersion)"
        ),
        yaxis_title="Planned contrast",
        height=max(520, 170 + 38 * len(tested)),
        template="plotly_white",
    )
    return figure


def benchmark_individual_figure(
    *, rows: Sequence[Mapping[str, Any]], metric_label: str, scale_label: str
) -> go.Figure:
    """Compare one cluster value with several empirical background medians.

    Args:
        rows: Individual comparison records sharing a metric and scale.
        metric_label: Plain-language dispersion metric.
        scale_label: Plain-language comparison scale.

    Returns:
        Horizontal observed-versus-background-median comparison plot.
    """

    tested = [
        row
        for row in rows
        if row.get("status") == "TESTED" and row.get("background_median") is not None
    ]
    tested.sort(key=lambda row: str(row["background_profile_id"]))
    backgrounds = [
        benchmark_profile_label(profile_id=str(row["background_profile_id"]))
        for row in tested
    ]
    figure = go.Figure()
    for position, row in enumerate(tested):
        observed = float(row["observed_value"])
        median = float(row["background_median"])
        figure.add_shape(
            type="line",
            x0=min(observed, median),
            x1=max(observed, median),
            y0=position,
            y1=position,
            line={"color": "#9ca3af", "width": 2},
        )
    figure.add_trace(
        go.Scatter(
            x=[float(row["background_median"]) for row in tested],
            y=backgrounds,
            mode="markers",
            name="Background median",
            marker={"symbol": "diamond", "size": 11, "color": "#6b7280"},
        )
    )
    figure.add_trace(
        go.Scatter(
            x=[float(row["observed_value"]) for row in tested],
            y=backgrounds,
            mode="markers",
            name="Selected cluster",
            marker={"size": 12, "color": "#d95f0e"},
            customdata=[
                [
                    row["empirical_percentile"],
                    row["two_sided_p_value"],
                    row["fdr_q_value"],
                    row["background_group_count"],
                ]
                for row in tested
            ],
            hovertemplate=(
                "Observed=%{x:.6g}<br>Percentile=%{customdata[0]:.3f}"
                "<br>p=%{customdata[1]:.4g}<br>FDR q=%{customdata[2]:.4g}"
                "<br>Background clusters=%{customdata[3]}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title=f"Selected cluster versus empirical backgrounds: {metric_label}",
        xaxis_title=f"{metric_label} ({scale_label})",
        yaxis_title="Comparison background",
        height=max(480, 170 + 48 * len(tested)),
        template="plotly_white",
    )
    return figure


def _stable_species_colour(*, species: str) -> str:
    """Return a deterministic HSL colour derived from an exact species label."""

    digest = hashlib.sha256(species.encode("utf-8")).digest()
    hue = int.from_bytes(digest[:4], "big") / (2**32) * 360
    saturation = 58 + digest[4] % 19
    lightness = 42 + digest[5] % 17
    return f"hsl({hue:.1f},{saturation}%,{lightness}%)"


def _series_colour(*, index: int) -> str:
    """Return a high-contrast deterministic colour for a comparison series."""

    palette = (
        "#2c7fb8",
        "#d95f0e",
        "#31a354",
        "#756bb1",
        "#e7298a",
        "#636363",
        "#1b9e77",
        "#e6ab02",
        "#a6761d",
        "#66a61e",
        "#7570b3",
        "#e41a1c",
    )
    return palette[index % len(palette)]


def linked_member_ids(
    *,
    entry: Mapping[str, Any],
    selected_members: Iterable[str] = (),
    selected_species: Iterable[str] = (),
) -> frozenset[str]:
    """Resolve linked member and species selections to exact member identifiers.

    Args:
        entry: One validated report visual entry.
        selected_members: Individually selected member identifiers.
        selected_species: Species whose displayed members should all be selected.

    Returns:
        Exact displayed members selected by either control.
    """

    members = _member_species(entry=entry)
    requested_members = {str(value) for value in selected_members}
    requested_species = {str(value) for value in selected_species}
    selected = requested_members.intersection(members)
    selected.update(
        member_id for member_id, species in members.items() if species in requested_species
    )
    return frozenset(selected)


def pcoa_figure(*, entry: Mapping[str, Any], selected_members: Iterable[str] = ()) -> go.Figure:
    """Build a distance-PCoA scatter plot with linked selection highlighting."""

    nodes = _visual_nodes(entry=entry)
    selected = {str(value) for value in selected_members}
    by_species: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for node in nodes:
        by_species[str(node.get("species", "Unlabelled"))].append(node)
    figure = go.Figure()
    for species in sorted(by_species):
        records = by_species[species]
        ids = [str(row["id"]) for row in records]
        active = [_is_active(member_id=value, selected=selected) for value in ids]
        figure.add_trace(
            go.Scattergl(
                x=[float(row["projectionX"]) for row in records],
                y=[float(row["projectionY"]) for row in records],
                mode="markers",
                name=species,
                customdata=[[member_id, species] for member_id in ids],
                marker={
                    "color": [_safe_colour(row.get("speciesColour")) for row in records],
                    "size": [14 if chosen else 8 for chosen in active],
                    "opacity": [1.0 if chosen else 0.14 for chosen in active],
                    "line": {
                        "color": [_HIGHLIGHT_BORDER if chosen else "#ffffff" for chosen in active],
                        "width": [2 if chosen else 0.5 for chosen in active],
                    },
                    "symbol": [
                        "star" if bool(row.get("isMedoid")) else "circle" for row in records
                    ],
                },
                hovertemplate=(
                    "Member=%{customdata[0]}<br>Species=%{customdata[1]}<br>"
                    "Axis 1=%{x:.6g}<br>Axis 2=%{y:.6g}<extra></extra>"
                ),
            )
        )
    projection = _mapping(entry.get("distanceProjection"), label="distance projection")
    axis_one = 100.0 * float(projection.get("axis_1_positive_inertia_fraction", 0.0))
    axis_two = 100.0 * float(projection.get("axis_2_positive_inertia_fraction", 0.0))
    figure.update_layout(
        title="Complete-distance PCoA (two-dimensional diagnostic)",
        xaxis_title=f"PCoA axis 1 ({axis_one:.1f}% positive inertia)",
        yaxis_title=f"PCoA axis 2 ({axis_two:.1f}% positive inertia)",
        legend_title="Species",
        height=680,
        template="plotly_white",
        dragmode="lasso",
    )
    figure.update_yaxes(scaleanchor="x", scaleratio=1)
    return figure


def shepard_figure(*, entry: Mapping[str, Any]) -> go.Figure:
    """Build a Shepard plot comparing exact and projected pair distances."""

    projection = _mapping(entry.get("distanceProjection"), label="distance projection")
    raw_points = projection.get("shepard_points")
    if not isinstance(raw_points, list) or not raw_points:
        raise InputValidationError("The selected group has no Shepard-plot points.")
    points: list[tuple[float, float]] = []
    for raw_point in raw_points:
        if not isinstance(raw_point, list) or len(raw_point) != 2:
            raise InputValidationError("The Shepard-plot payload contains a malformed point.")
        actual, projected = float(raw_point[0]), float(raw_point[1])
        if not math.isfinite(actual) or not math.isfinite(projected):
            raise InputValidationError("The Shepard-plot payload contains a non-finite point.")
        points.append((actual, projected))
    maximum = max(max(point) for point in points)
    figure = go.Figure()
    figure.add_trace(
        go.Scattergl(
            x=[point[0] for point in points],
            y=[point[1] for point in points],
            mode="markers",
            name="Deterministic stratified pairs",
            marker={"size": 5, "opacity": 0.45, "color": "#355c8a"},
            hovertemplate=(
                "Input distance=%{x:.6g}<br>Projected 2D distance=%{y:.6g}<extra></extra>"
            ),
        )
    )
    figure.add_trace(
        go.Scatter(
            x=[0.0, maximum],
            y=[0.0, maximum],
            mode="lines",
            name="Perfect preservation",
            line={"color": "#b7791f", "dash": "dash"},
            hoverinfo="skip",
        )
    )
    figure.update_layout(
        title="Shepard plot",
        xaxis_title="Exact input pairwise distance",
        yaxis_title="Distance in the two-dimensional PCoA",
        height=560,
        template="plotly_white",
    )
    figure.update_yaxes(scaleanchor="x", scaleratio=1)
    return figure


def phylogram_figure(
    *,
    entry: Mapping[str, Any],
    selected_members: Iterable[str] = (),
    show_all_labels: bool = False,
) -> go.Figure:
    """Build a rectangular branch-length phylogram from a pruned resolved tree."""

    tree = _mapping(entry.get("phylogram"), label="phylogram")
    status = str(tree.get("status", ""))
    if "PRUNED_PHYLOGRAM" not in status:
        raise InputValidationError(
            f"The selected group has no usable pruned phylogram: {tree.get('reason', status)}"
        )
    raw_nodes = tree.get("nodes")
    raw_edges = tree.get("edges")
    if not isinstance(raw_nodes, list) or not isinstance(raw_edges, list):
        raise InputValidationError("The phylogram payload lacks node or edge arrays.")
    nodes = {
        str(row.get("id", "")): row for row in raw_nodes if isinstance(row, dict) and row.get("id")
    }
    children: dict[str, list[str]] = defaultdict(list)
    horizontal_x: list[float | None] = []
    horizontal_y: list[float | None] = []
    for edge in raw_edges:
        if not isinstance(edge, dict):
            continue
        parent_id, child_id = str(edge.get("parentId", "")), str(edge.get("childId", ""))
        parent, child = nodes.get(parent_id), nodes.get(child_id)
        if parent is None or child is None:
            continue
        children[parent_id].append(child_id)
        horizontal_x.extend([float(parent["x"]), float(child["x"]), None])
        horizontal_y.extend([float(child["y"]), float(child["y"]), None])
    vertical_x: list[float | None] = []
    vertical_y: list[float | None] = []
    for parent_id, child_ids in children.items():
        child_y = [float(nodes[child_id]["y"]) for child_id in child_ids]
        parent_x = float(nodes[parent_id]["x"])
        vertical_x.extend([parent_x, parent_x, None])
        vertical_y.extend([min(child_y), max(child_y), None])
    figure = go.Figure()
    figure.add_trace(_tree_line_trace(x=vertical_x, y=vertical_y))
    figure.add_trace(_tree_line_trace(x=horizontal_x, y=horizontal_y))
    leaves = [row for row in nodes.values() if bool(row.get("isLeaf"))]
    selected = {str(value) for value in selected_members}
    ids = [str(row.get("memberId", "")) for row in leaves]
    active = [_is_active(member_id=value, selected=selected) for value in ids]
    labels = [
        member_id if show_all_labels or chosen or bool(row.get("isMedoid")) else ""
        for member_id, chosen, row in zip(ids, active, leaves, strict=True)
    ]
    figure.add_trace(
        go.Scattergl(
            x=[float(row["x"]) for row in leaves],
            y=[float(row["y"]) for row in leaves],
            mode="markers+text",
            text=labels,
            textposition="middle right",
            customdata=[
                [ids[index], str(row.get("species", ""))] for index, row in enumerate(leaves)
            ],
            marker={
                "color": [_safe_colour(row.get("colour")) for row in leaves],
                "size": [13 if chosen else 7 for chosen in active],
                "opacity": [1.0 if chosen else 0.14 for chosen in active],
                "line": {
                    "color": [_HIGHLIGHT_BORDER if chosen else "#ffffff" for chosen in active],
                    "width": [2 if chosen else 0.5 for chosen in active],
                },
                "symbol": ["star" if bool(row.get("isMedoid")) else "circle" for row in leaves],
            },
            name="Displayed tree leaves",
            hovertemplate=(
                "Member=%{customdata[0]}<br>Species=%{customdata[1]}<br>"
                "Root distance=%{x:.6g}<extra></extra>"
            ),
        )
    )
    height = max(520, min(1_800, 180 + 7 * len(leaves)))
    figure.update_layout(
        title="Pruned resolved-gene-tree phylogram",
        xaxis_title="Cumulative branch length from displayed root",
        yaxis={"visible": False},
        height=height,
        template="plotly_white",
        showlegend=False,
        hovermode="closest",
    )
    return figure


def distance_matrix_figure(
    *, entry: Mapping[str, Any], selected_members: Sequence[str] = ()
) -> go.Figure:
    """Build an exact symmetric heatmap, optionally restricted to linked members."""

    labels, matrix = distance_matrix_array(entry=entry)
    selected = tuple(dict.fromkeys(str(value) for value in selected_members))
    if selected:
        missing = sorted(set(selected).difference(labels))
        if missing:
            raise InputValidationError(
                "Selected matrix members are absent from the exact matrix: "
                + "; ".join(missing[:10])
            )
        indices = [labels.index(member_id) for member_id in selected]
        labels = selected
        matrix = matrix[np.ix_(indices, indices)]
    figure = go.Figure(
        data=go.Heatmap(
            z=matrix,
            x=labels,
            y=labels,
            colorscale="YlOrRd",
            colorbar={"title": "Exact distance"},
            hovertemplate=("Row=%{y}<br>Column=%{x}<br>Exact distance=%{z:.8g}<extra></extra>"),
        )
    )
    figure.update_layout(
        title=f"Exact complete displayed distance matrix (n={len(labels):,})",
        xaxis_title="Displayed member",
        yaxis_title="Displayed member",
        height=max(650, min(1_300, 250 + 4 * len(labels))),
        template="plotly_white",
    )
    figure.update_yaxes(autorange="reversed")
    return figure


def distance_matrix_array(*, entry: Mapping[str, Any]) -> tuple[tuple[str, ...], np.ndarray]:
    """Reconstruct the exact symmetric matrix from its upper-triangle payload."""

    payload = _mapping(entry.get("distanceMatrix"), label="distance matrix")
    if payload.get("status") != "EXACT_COMPLETE_DISPLAYED_MATRIX":
        raise InputValidationError(
            f"The selected group has no exact complete displayed matrix: "
            f"{payload.get('reason', payload.get('status', 'unknown status'))}"
        )
    raw_order, raw_values = payload.get("memberOrder"), payload.get("upperTriangle")
    if not isinstance(raw_order, list) or not isinstance(raw_values, list):
        raise InputValidationError("The exact distance matrix payload is malformed.")
    labels = tuple(str(value) for value in raw_order)
    if len(labels) < 2 or len(set(labels)) != len(labels):
        raise InputValidationError("The exact distance matrix member order is invalid.")
    expected = len(labels) * (len(labels) - 1) // 2
    if len(raw_values) != expected:
        raise InputValidationError(
            f"The exact matrix requires {expected:,} values; observed {len(raw_values):,}."
        )
    matrix = np.zeros((len(labels), len(labels)), dtype=float)
    value_index = 0
    for left in range(len(labels)):
        for right in range(left + 1, len(labels)):
            value = float(raw_values[value_index])
            if not math.isfinite(value) or value < 0:
                raise InputValidationError("The exact matrix contains an invalid distance.")
            matrix[left, right] = value
            matrix[right, left] = value
            value_index += 1
    return labels, matrix


def nearest_neighbour_figure(
    *,
    entry: Mapping[str, Any],
    selected_members: Iterable[str] = (),
    include_connectors: bool = False,
) -> go.Figure:
    """Build a deterministic force layout for nearest-neighbour topology."""

    nodes = _visual_nodes(entry=entry)
    raw_edges = entry.get("edges")
    if not isinstance(raw_edges, list):
        raise InputValidationError("The nearest-neighbour payload lacks an edge array.")
    graph = nx.Graph()
    for node in nodes:
        graph.add_node(str(node["id"]))
    retained_edges: list[dict[str, Any]] = []
    for raw_edge in raw_edges:
        if not isinstance(raw_edge, dict):
            continue
        edge_type = str(raw_edge.get("edgeType", ""))
        if edge_type == "COMPONENT_CONNECTOR" and not include_connectors:
            continue
        if edge_type not in {"NEAREST_NEIGHBOUR", "COMPONENT_CONNECTOR"}:
            continue
        left, right = str(raw_edge.get("from", "")), str(raw_edge.get("to", ""))
        if left not in graph or right not in graph:
            continue
        distance = max(float(raw_edge.get("distance", 0.0)), 0.0)
        graph.add_edge(left, right, spring_strength=1.0 / (distance + 1e-9))
        retained_edges.append(raw_edge)
    if not graph:
        raise InputValidationError("The nearest-neighbour graph contains no nodes.")
    positions = nx.spring_layout(
        graph,
        seed=1729,
        iterations=150,
        weight="spring_strength",
    )
    figure = go.Figure()
    for edge_type, colour, dash, name in (
        ("NEAREST_NEIGHBOUR", "#7786a8", "solid", "Nearest-neighbour edges"),
        ("COMPONENT_CONNECTOR", "#b7791f", "dash", "Layout-only connectors"),
    ):
        x_values: list[float | None] = []
        y_values: list[float | None] = []
        for edge in retained_edges:
            if str(edge.get("edgeType")) != edge_type:
                continue
            left, right = str(edge["from"]), str(edge["to"])
            x_values.extend([float(positions[left][0]), float(positions[right][0]), None])
            y_values.extend([float(positions[left][1]), float(positions[right][1]), None])
        if x_values:
            figure.add_trace(
                go.Scattergl(
                    x=x_values,
                    y=y_values,
                    mode="lines",
                    line={"color": colour, "width": 1.2, "dash": dash},
                    name=name,
                    hoverinfo="skip",
                )
            )
    selected = {str(value) for value in selected_members}
    ids = [str(node["id"]) for node in nodes]
    active = [_is_active(member_id=value, selected=selected) for value in ids]
    figure.add_trace(
        go.Scattergl(
            x=[float(positions[member_id][0]) for member_id in ids],
            y=[float(positions[member_id][1]) for member_id in ids],
            mode="markers",
            customdata=[
                [ids[index], str(node.get("species", ""))] for index, node in enumerate(nodes)
            ],
            marker={
                "color": [_safe_colour(node.get("speciesColour")) for node in nodes],
                "size": [14 if chosen else 8 for chosen in active],
                "opacity": [1.0 if chosen else 0.14 for chosen in active],
                "line": {
                    "color": [_HIGHLIGHT_BORDER if chosen else "#ffffff" for chosen in active],
                    "width": [2 if chosen else 0.5 for chosen in active],
                },
                "symbol": ["star" if bool(node.get("isMedoid")) else "circle" for node in nodes],
            },
            name="Displayed members",
            hovertemplate="Member=%{customdata[0]}<br>Species=%{customdata[1]}<extra></extra>",
        )
    )
    figure.update_layout(
        title="Nearest-neighbour topology (force layout)",
        xaxis={"visible": False},
        yaxis={"visible": False},
        height=700,
        template="plotly_white",
        hovermode="closest",
    )
    return figure


def _member_species(*, entry: Mapping[str, Any]) -> dict[str, str]:
    """Return exact displayed member-to-species mappings."""

    raw_members = entry.get("members")
    if not isinstance(raw_members, list):
        raise InputValidationError("The visual entry lacks a member array.")
    result: dict[str, str] = {}
    for row in raw_members:
        if not isinstance(row, dict):
            raise InputValidationError("The visual entry contains a malformed member.")
        member_id = str(row.get("member_id", ""))
        if not member_id or member_id in result:
            raise InputValidationError("The visual entry contains duplicate/empty members.")
        result[member_id] = str(row.get("species_label", ""))
    return result


def _visual_nodes(*, entry: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return validated displayed nodes with finite PCoA coordinates."""

    raw_nodes = entry.get("nodes")
    if not isinstance(raw_nodes, list):
        raise InputValidationError("The visual entry lacks a node array.")
    nodes: list[dict[str, Any]] = []
    for raw_node in raw_nodes:
        if not isinstance(raw_node, dict) or not raw_node.get("id"):
            raise InputValidationError("The visual entry contains a malformed node.")
        try:
            x_value = float(raw_node["projectionX"])
            y_value = float(raw_node["projectionY"])
        except (KeyError, TypeError, ValueError) as error:
            raise InputValidationError(
                "The visual entry contains a node without valid PCoA coordinates."
            ) from error
        if not math.isfinite(x_value) or not math.isfinite(y_value):
            raise InputValidationError("The visual entry contains non-finite coordinates.")
        nodes.append(raw_node)
    return nodes


def _mapping(value: object, *, label: str) -> Mapping[str, Any]:
    """Return a mapping or raise a user-facing payload error."""

    if not isinstance(value, Mapping):
        raise InputValidationError(f"The visual entry lacks a valid {label} record.")
    return value


def _tree_line_trace(*, x: Sequence[float | None], y: Sequence[float | None]) -> go.Scattergl:
    """Return one non-interactive phylogram line trace."""

    return go.Scattergl(
        x=x,
        y=y,
        mode="lines",
        line={"color": "#657084", "width": 1},
        hoverinfo="skip",
        showlegend=False,
    )


def _safe_colour(value: object) -> str:
    """Return a bounded CSS colour token from the trusted report palette."""

    colour = str(value or _DEFAULT_COLOUR)
    if len(colour) > 64 or any(character in colour for character in "<>;{}"):
        return _DEFAULT_COLOUR
    return colour


def _is_active(*, member_id: str, selected: set[str]) -> bool:
    """Treat every member as active until a linked selection is present."""

    return not selected or member_id in selected

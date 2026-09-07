"""Plotly figures for bounded evolutionary visualisation records."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import networkx as nx
import numpy as np
import plotly.graph_objects as go

from orthofinder_results.errors import InputValidationError

_DEFAULT_COLOUR = "#5b6475"
_HIGHLIGHT_BORDER = "#111827"


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

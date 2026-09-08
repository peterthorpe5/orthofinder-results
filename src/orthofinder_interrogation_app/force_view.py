"""Self-contained draggable force-directed nearest-neighbour view."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

from pyvis.network import Network

from orthofinder_results.errors import InputValidationError

MAX_FORCE_NODES = 500
MAX_FORCE_EDGES = 25_000


def force_directed_html(
    *,
    entry: Mapping[str, Any],
    selected_members: Iterable[str] = (),
    include_connectors: bool = False,
    show_labels: bool = False,
    physics_enabled: bool = True,
) -> str:
    """Return a self-contained draggable vis-network document.

    Args:
        entry: Validated group visualisation record.
        selected_members: Linked canonical members to emphasise.
        include_connectors: Include explicit layout-only component connectors.
        show_labels: Display every member identifier beside its node.
        physics_enabled: Keep force-layout physics active after stabilisation.

    Returns:
        Complete HTML document with JavaScript and CSS embedded locally.

    Raises:
        InputValidationError: If the graph payload is malformed or unsafe.
    """

    raw_nodes, raw_edges = entry.get("nodes"), entry.get("edges")
    if not isinstance(raw_nodes, list) or not isinstance(raw_edges, list):
        raise InputValidationError("Force-directed view lacks node or edge arrays.")
    if not 1 <= len(raw_nodes) <= MAX_FORCE_NODES:
        raise InputValidationError(
            f"Force-directed view requires 1–{MAX_FORCE_NODES:,} nodes; "
            f"observed {len(raw_nodes):,}."
        )
    if len(raw_edges) > MAX_FORCE_EDGES:
        raise InputValidationError(
            f"Force-directed view exceeds the {MAX_FORCE_EDGES:,}-edge safety limit."
        )
    selected = {str(value) for value in selected_members}
    network = Network(
        height="720px",
        width="100%",
        bgcolor="#ffffff",
        font_color="#172033",
        directed=False,
        notebook=False,
        cdn_resources="in_line",
    )
    node_ids: set[str] = set()
    for raw_node in raw_nodes:
        if not isinstance(raw_node, dict):
            raise InputValidationError("Force-directed view contains a malformed node.")
        node_id = _safe_graph_text(value=raw_node.get("id", ""), role="node identifier")
        if not node_id or node_id in node_ids:
            raise InputValidationError("Force-directed nodes require unique identifiers.")
        node_ids.add(node_id)
        is_medoid = bool(raw_node.get("isMedoid", False))
        is_selected = node_id in selected
        species = _safe_graph_text(
            value=raw_node.get("species", "Unlabelled"), role="species label"
        )
        member_label = _safe_graph_text(
            value=raw_node.get("memberLabel", node_id), role="member label"
        )
        label = member_label if show_labels or is_medoid or is_selected else ""
        colour = _node_colour(raw_node=raw_node, selected=is_selected, medoid=is_medoid)
        network.add_node(
            node_id,
            label=label,
            title=_safe_graph_text(
                value=raw_node.get("title", f"{member_label} | {species}"),
                role="node title",
            ),
            color=colour,
            shape="star" if is_medoid else "dot",
            size=30 if is_medoid or is_selected else 16,
            borderWidth=5 if is_medoid or is_selected else 1,
            species=species,
        )
    for raw_edge in raw_edges:
        if not isinstance(raw_edge, dict):
            raise InputValidationError("Force-directed view contains a malformed edge.")
        edge_type = str(raw_edge.get("edgeType", ""))
        if edge_type == "COMPONENT_CONNECTOR" and not include_connectors:
            continue
        if edge_type not in {"NEAREST_NEIGHBOUR", "COMPONENT_CONNECTOR"}:
            continue
        left = _safe_graph_text(value=raw_edge.get("from", ""), role="edge endpoint")
        right = _safe_graph_text(value=raw_edge.get("to", ""), role="edge endpoint")
        if left not in node_ids or right not in node_ids or left == right:
            raise InputValidationError("Force-directed edge has an invalid endpoint.")
        connector = edge_type == "COMPONENT_CONNECTOR"
        network.add_edge(
            left,
            right,
            title=_safe_graph_text(
                value=raw_edge.get("title", edge_type), role="edge title"
            ),
            color="#b7791f" if connector else "#718096",
            width=2.0 if connector else 1.1,
            dashes=connector,
            edgeType=edge_type,
        )
    network.set_options(_network_options(physics_enabled=physics_enabled))
    document = network.generate_html(notebook=False)
    if not isinstance(document, str) or "vis.Network" not in document:
        raise InputValidationError("Force-directed HTML generation failed.")
    return document


def _node_colour(
    *, raw_node: Mapping[str, Any], selected: bool, medoid: bool
) -> dict[str, Any]:
    """Return selected, medoid and default node colours."""

    raw_colour = raw_node.get("color")
    if isinstance(raw_colour, dict):
        background = str(raw_colour.get("background", "#5b6475"))
    else:
        background = str(raw_node.get("speciesColour", "#5b6475"))
    border = "#111827" if selected else "#7a4e00" if medoid else "#4b5563"
    return {
        "background": background,
        "border": border,
        "highlight": {"background": background, "border": "#111827"},
        "hover": {"background": background, "border": "#111827"},
    }


def _network_options(*, physics_enabled: bool) -> str:
    """Return deterministic local interaction and physics options as JSON text."""

    enabled = "true" if physics_enabled else "false"
    return """
    {
      "interaction": {
        "dragNodes": true,
        "dragView": true,
        "hover": true,
        "keyboard": {"enabled": true},
        "multiselect": true,
        "navigationButtons": true,
        "tooltipDelay": 120,
        "zoomView": true
      },
      "nodes": {"font": {"size": 13, "face": "Arial"}},
      "edges": {"smooth": false},
      "physics": {
        "enabled": %s,
        "barnesHut": {
          "gravitationalConstant": -5200,
          "centralGravity": 0.18,
          "springLength": 105,
          "springConstant": 0.035,
          "damping": 0.16,
          "avoidOverlap": 0.35
        },
        "stabilization": {"enabled": true, "iterations": 450, "fit": true}
      }
    }
    """ % enabled


def _safe_graph_text(*, value: object, role: str) -> str:
    """Return graph text that cannot terminate an embedded script element.

    Args:
        value: Candidate identifier or display text.
        role: User-facing field role for controlled errors.

    Returns:
        Original text when it contains no HTML delimiters or control characters.

    Raises:
        InputValidationError: If unsafe embedded-HTML characters are present.
    """

    text = str(value)
    if not text or "<" in text or ">" in text or any(ord(character) < 32 for character in text):
        raise InputValidationError(f"Force-directed {role} contains unsafe text.")
    return text

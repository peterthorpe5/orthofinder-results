"""Accessible interactive and static rendering for taxonomy coverage trees."""

from __future__ import annotations

import html
import io
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from orthofinder_results.errors import InputValidationError

from .coverage_tree import CoverageTree, TaxonCoverage

if TYPE_CHECKING:
    import plotly.graph_objects as go

_BLUE = "#0072B2"
_LIGHT_BLUE = "#56B4E9"
_RED = "#D55E00"
_AMBER = "#E69F00"
_GREEN = "#009E73"
_GREY = "#6B7280"
_DARK = "#1F2937"


@dataclass(frozen=True)
class TreePosition:
    """One deterministic rectangular-tree plot coordinate."""

    taxon_id: str
    x: float
    y: float


@dataclass(frozen=True)
class NodeVisual:
    """Accessible colour-plus-shape styling for one coverage node."""

    fill_colour: str
    outline_colour: str
    shape: str
    outline_width: float
    outline_style: str
    legend_label: str


def coverage_tree_layout(*, tree: CoverageTree) -> tuple[TreePosition, ...]:
    """Return deterministic rooted rectangular-tree coordinates."""

    children = tree.children_by_id
    nodes = tree.node_by_id
    leaf_order: list[str] = []

    def visit(taxon_id: str) -> None:
        """Append displayed leaves in stable depth-first order."""

        if not children.get(taxon_id):
            leaf_order.append(taxon_id)
            return
        for child in children[taxon_id]:
            visit(child)

    visit(tree.root_taxon_id)
    if not leaf_order:
        raise InputValidationError("Coverage tree contains no displayed leaves.")
    y_by_id = {taxon_id: float(index) for index, taxon_id in enumerate(leaf_order)}

    def assign_internal(taxon_id: str) -> float:
        """Assign internal nodes to the mean of direct-child positions."""

        child_ids = children.get(taxon_id, ())
        if not child_ids:
            return y_by_id[taxon_id]
        child_y = [assign_internal(child) for child in child_ids]
        y_by_id[taxon_id] = sum(child_y) / len(child_y)
        return y_by_id[taxon_id]

    assign_internal(tree.root_taxon_id)
    root_depth = nodes[tree.root_taxon_id].lineage_depth
    return tuple(
        TreePosition(
            taxon_id=taxon_id,
            x=float(nodes[taxon_id].lineage_depth - root_depth),
            y=y_by_id[taxon_id],
        )
        for taxon_id in sorted(
            nodes,
            key=lambda value: (
                nodes[value].lineage_depth,
                nodes[value].taxon_name.casefold(),
                value,
            ),
        )
    )


def coverage_tree_figure(*, tree: CoverageTree) -> "go.Figure":
    """Build an interactive taxonomy tree with accessible orthogonal states."""

    try:
        import plotly.graph_objects as go
    except ModuleNotFoundError as error:  # pragma: no cover - dependency gate.
        raise InputValidationError(
            "Interactive coverage trees require the optional application dependencies."
        ) from error

    positions = {row.taxon_id: row for row in coverage_tree_layout(tree=tree)}
    nodes = tree.node_by_id
    figure = go.Figure()
    for edge in tree.edges:
        parent = positions[edge.parent_taxon_id]
        child = positions[edge.child_taxon_id]
        colour, width, dash = _edge_visual(node=nodes[edge.child_taxon_id])
        figure.add_trace(
            go.Scatter(
                x=[parent.x, child.x, child.x],
                y=[parent.y, parent.y, child.y],
                mode="lines",
                line={"color": colour, "width": width, "dash": dash},
                hoverinfo="skip",
                showlegend=False,
            )
        )
    ordered = [nodes[position.taxon_id] for position in positions.values()]
    visuals = [_node_visual(node=node) for node in ordered]
    hover = [
        "<br>".join(
            (
                f"Taxon: {html.escape(node.taxon_name)}",
                f"Rank: {html.escape(node.taxon_rank)}",
                f"ID: {html.escape(node.taxon_id)}",
                f"Selection: {html.escape(node.selection_state)}",
                f"Coverage: {html.escape(node.coverage_state)}",
                f"Represented descendants: {node.represented_terminal_count}",
                f"Expected-no-data descendants: {node.expected_no_data_terminal_count}",
            )
        )
        for node in ordered
    ]
    figure.add_trace(
        go.Scatter(
            x=[positions[node.taxon_id].x for node in ordered],
            y=[positions[node.taxon_id].y for node in ordered],
            mode="markers+text",
            text=[
                f"{node.taxon_name} ({node.taxon_rank}; {node.taxon_id})"
                for node in ordered
            ],
            textposition="middle right",
            customdata=[[node.taxon_id] for node in ordered],
            hovertext=hover,
            hoverinfo="text",
            marker={
                "color": [visual.fill_colour for visual in visuals],
                "symbol": [_plotly_symbol(visual=visual) for visual in visuals],
                "size": [15 if node.is_selected or node.is_excluded else 11 for node in ordered],
                "line": {
                    "color": [visual.outline_colour for visual in visuals],
                    "width": [visual.outline_width for visual in visuals],
                },
            },
            name="Taxonomy nodes",
            showlegend=False,
        )
    )
    height = min(2_400, max(560, 28 * len(_leaf_ids(tree=tree))))
    figure.update_layout(
        title="Selection coverage tree — reviewed taxonomy authority",
        template="plotly_white",
        height=height,
        dragmode="pan",
        hovermode="closest",
        margin={"l": 30, "r": 360, "t": 75, "b": 40},
        xaxis={"title": "Reviewed lineage depth", "showgrid": False, "zeroline": False},
        yaxis={"showticklabels": False, "showgrid": False, "zeroline": False},
    )
    return figure


def style_records(*, tree: CoverageTree) -> tuple[dict[str, Any], ...]:
    """Return one visual-style audit row per displayed taxon."""

    records = []
    for node in tree.nodes:
        visual = _node_visual(node=node)
        edge_colour, edge_width, edge_style = _edge_visual(node=node)
        records.append(
            {
                "taxon_id": node.taxon_id,
                "selection_state": node.selection_state,
                "coverage_state": node.coverage_state,
                "node_fill_colour": visual.fill_colour,
                "node_outline_colour": visual.outline_colour,
                "node_shape": visual.shape,
                "node_outline_width": visual.outline_width,
                "node_outline_style": visual.outline_style,
                "incoming_edge_colour": edge_colour,
                "incoming_edge_width": edge_width,
                "incoming_edge_style": edge_style,
                "legend_label": visual.legend_label,
            }
        )
    return tuple(records)


def tree_to_svg(*, tree: CoverageTree) -> bytes:
    """Render a self-contained deterministic SVG with an accessible legend."""

    positions = {row.taxon_id: row for row in coverage_tree_layout(tree=tree)}
    nodes = tree.node_by_id
    leaves = _leaf_ids(tree=tree)
    x_scale = 145.0
    y_scale = 25.0
    left = 55.0
    top = 65.0
    legend_height = 165.0
    maximum_depth = max(position.x for position in positions.values())
    width = max(1_400.0, left + maximum_depth * x_scale + 620.0)
    height = max(500.0, top + max(1, len(leaves) - 1) * y_scale + legend_height)
    output = io.StringIO()
    output.write(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" '
        f'height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}" '
        'role="img" aria-labelledby="title desc">\n'
    )
    output.write('<title id="title">Selection coverage tree</title>\n')
    output.write(
        '<desc id="desc">Reviewed taxonomy selection and dataset coverage; '
        'selection and coverage are encoded independently.</desc>\n'
    )
    output.write('<rect width="100%" height="100%" fill="#ffffff"/>\n')
    output.write(
        '<text x="30" y="32" font-family="Arial,sans-serif" font-size="20" '
        'font-weight="bold" fill="#1F2937">Selection coverage tree</text>\n'
    )
    for edge in tree.edges:
        parent = positions[edge.parent_taxon_id]
        child = positions[edge.child_taxon_id]
        colour, edge_width, dash = _edge_visual(node=nodes[edge.child_taxon_id])
        dash_attribute = ' stroke-dasharray="7 5"' if dash == "dash" else ""
        parent_x, parent_y = left + parent.x * x_scale, top + parent.y * y_scale
        child_x, child_y = left + child.x * x_scale, top + child.y * y_scale
        output.write(
            f'<path d="M {parent_x:.2f} {parent_y:.2f} H {child_x:.2f} '
            f'V {child_y:.2f}" fill="none" stroke="{colour}" '
            f'stroke-width="{edge_width:.1f}"{dash_attribute}/>\n'
        )
    for taxon_id in sorted(
        positions,
        key=lambda value: (positions[value].y, positions[value].x, value),
    ):
        node = nodes[taxon_id]
        visual = _node_visual(node=node)
        x = left + positions[taxon_id].x * x_scale
        y = top + positions[taxon_id].y * y_scale
        output.write(_svg_marker(x=x, y=y, visual=visual, size=7.0))
        label = html.escape(f"{node.taxon_name} ({node.taxon_rank}; ID {node.taxon_id})")
        output.write(
            f'<text x="{x + 12:.2f}" y="{y + 4:.2f}" '
            'font-family="Arial,sans-serif" font-size="12" fill="#1F2937">'
            f"{label}</text>\n"
        )
    legend_y = height - 128.0
    output.write(
        f'<text x="30" y="{legend_y - 24:.2f}" font-family="Arial,sans-serif" '
        'font-size="15" font-weight="bold" fill="#1F2937">Legend</text>\n'
    )
    for index, (label, visual) in enumerate(_legend_visuals()):
        column, row = index % 3, index // 3
        x = 45.0 + column * 430.0
        y = legend_y + row * 34.0
        output.write(_svg_marker(x=x, y=y, visual=visual, size=7.0))
        output.write(
            f'<text x="{x + 14:.2f}" y="{y + 4:.2f}" '
            'font-family="Arial,sans-serif" font-size="12" fill="#1F2937">'
            f"{html.escape(label)}</text>\n"
        )
    output.write(
        f'<text x="30" y="{height - 18:.2f}" font-family="Arial,sans-serif" '
        'font-size="11" fill="#6B7280">Expected no data means not represented in '
        'this dataset; it is not a claim of biological absence.</text>\n'
    )
    output.write("</svg>\n")
    return output.getvalue().encode("utf-8")


def tree_to_pdf(*, tree: CoverageTree) -> bytes:
    """Render a dependency-free deterministic one-page vector PDF."""

    positions = {row.taxon_id: row for row in coverage_tree_layout(tree=tree)}
    nodes = tree.node_by_id
    maximum_x = max((row.x for row in positions.values()), default=1.0) or 1.0
    maximum_y = max((row.y for row in positions.values()), default=1.0) or 1.0
    commands = [
        "BT /F1 16 Tf 40 560 Td (Selection coverage tree) Tj ET",
        "BT /F1 8 Tf 40 545 Td "
        "(Reviewed taxonomy; coverage is not evidence of biological absence.) Tj ET",
    ]
    xy: dict[str, tuple[float, float]] = {}
    for taxon_id, position in positions.items():
        xy[taxon_id] = (
            55.0 + 470.0 * position.x / maximum_x,
            70.0 + 445.0 * position.y / maximum_y,
        )
    for edge in tree.edges:
        parent_x, parent_y = xy[edge.parent_taxon_id]
        child_x, child_y = xy[edge.child_taxon_id]
        colour, width, dash = _edge_visual(node=nodes[edge.child_taxon_id])
        red, green, blue = _hex_rgb(value=colour)
        commands.append(f"{red:.3f} {green:.3f} {blue:.3f} RG {width:.2f} w")
        commands.append("[5 4] 0 d" if dash == "dash" else "[] 0 d")
        commands.append(
            f"{parent_x:.2f} {parent_y:.2f} m {child_x:.2f} {parent_y:.2f} l "
            f"{child_x:.2f} {child_y:.2f} l S"
        )
    for taxon_id in sorted(xy, key=lambda value: (xy[value][1], xy[value][0], value)):
        node = nodes[taxon_id]
        visual = _node_visual(node=node)
        x, y = xy[taxon_id]
        commands.extend(_pdf_marker(x=x, y=y, visual=visual))
        label = _pdf_text(
            value=f"{node.taxon_name} [{node.taxon_rank}; {node.taxon_id}]"[:62]
        )
        commands.append(f"BT /F1 5 Tf {x + 5:.2f} {y - 1.5:.2f} Td ({label}) Tj ET")
    commands.extend(
        (
            "0 0 0 RG [] 0 d",
            "BT /F1 7 Tf 560 515 Td (Legend:) Tj ET",
            "BT /F1 6 Tf 560 500 Td (filled circle = represented input) Tj ET",
            "BT /F1 6 Tf 560 488 Td (hollow circle = expected no data) Tj ET",
            "BT /F1 6 Tf 560 476 Td (amber diamond = outside-scope hit) Tj ET",
            "BT /F1 6 Tf 560 464 Td (blue outline = selected) Tj ET",
            "BT /F1 6 Tf 560 452 Td (red outline = excluded) Tj ET",
        )
    )
    return _minimal_pdf(content="\n".join(commands).encode("latin-1", errors="replace"))


def _leaf_ids(*, tree: CoverageTree) -> tuple[str, ...]:
    """Return displayed leaves in deterministic depth-first order."""

    children = tree.children_by_id
    leaves: list[str] = []

    def visit(taxon_id: str) -> None:
        """Append terminal displayed nodes recursively."""

        if not children.get(taxon_id):
            leaves.append(taxon_id)
            return
        for child in children[taxon_id]:
            visit(child)

    visit(tree.root_taxon_id)
    return tuple(leaves)


def _node_visual(*, node: TaxonCoverage) -> NodeVisual:
    """Return colour-plus-shape semantics for one node."""

    coverage = {
        "REPRESENTED_IN_INPUT": (_GREEN, "circle", "Represented in input"),
        "EXPECTED_NO_DATA": ("#FFFFFF", "circle-open", "Expected but not represented"),
        "ANCESTOR_OF_REPRESENTED": ("#D1D5DB", "square", "Ancestor of represented"),
        "OUTSIDE_SELECTED_SCOPE_WITH_HITS": (_AMBER, "diamond", "Outside-scope hit"),
        "NOT_IN_EXPECTED_UNIVERSE": ("#FFFFFF", "square-open", "Not expected"),
    }
    fill, shape, label = coverage.get(node.coverage_state, (_GREY, "square", "Other"))
    if node.is_excluded:
        outline, width, style = _RED, 3.5, "dashed"
    elif node.is_selected:
        outline, width, style = _BLUE, 3.5, "solid"
    elif node.coverage_state == "EXPECTED_NO_DATA":
        outline, width, style = _DARK, 1.5, "dashed"
    else:
        outline, width, style = _DARK, 1.2, "solid"
    return NodeVisual(fill, outline, shape, width, style, label)


def _edge_visual(*, node: TaxonCoverage) -> tuple[str, float, str]:
    """Return incoming branch colour, width and dash style."""

    if node.is_excluded or node.selection_state == "EXCLUDED_BY_CLADE":
        return _RED, 2.8, "dash"
    if "ONLY_IN_SCOPE" in node.selection_state:
        return _BLUE, 3.4, "solid"
    if "INCLUDE_CLADE" in node.selection_state or "REQUIRED_EXACT" in node.selection_state:
        return _LIGHT_BLUE, 3.2, "solid"
    return "#9CA3AF", 1.3, "solid"


def _plotly_symbol(*, visual: NodeVisual) -> str:
    """Translate renderer-neutral node shapes to Plotly symbols."""

    return {
        "circle": "circle",
        "circle-open": "circle-open",
        "diamond": "diamond",
        "square": "square",
        "square-open": "square-open",
    }[visual.shape]


def _svg_marker(*, x: float, y: float, visual: NodeVisual, size: float) -> str:
    """Return one escaped SVG marker carrying colour and shape semantics."""

    dash = ' stroke-dasharray="3 2"' if visual.outline_style == "dashed" else ""
    common = (
        f'fill="{visual.fill_colour}" stroke="{visual.outline_colour}" '
        f'stroke-width="{visual.outline_width:.1f}"{dash}'
    )
    if visual.shape.startswith("circle"):
        return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{size:.2f}" {common}/>\n'
    if visual.shape == "diamond":
        points = (
            f"{x:.2f},{y - size:.2f} {x + size:.2f},{y:.2f} "
            f"{x:.2f},{y + size:.2f} {x - size:.2f},{y:.2f}"
        )
        return f'<polygon points="{points}" {common}/>\n'
    return (
        f'<rect x="{x - size:.2f}" y="{y - size:.2f}" width="{2 * size:.2f}" '
        f'height="{2 * size:.2f}" {common}/>\n'
    )


def _legend_visuals() -> tuple[tuple[str, NodeVisual], ...]:
    """Return canonical exported legend entries."""

    return (
        ("Represented input", NodeVisual(_GREEN, _DARK, "circle", 1.2, "solid", "")),
        (
            "Expected, not represented in this dataset",
            NodeVisual("#FFFFFF", _DARK, "circle-open", 1.5, "dashed", ""),
        ),
        (
            "Outside selected scope with hits",
            NodeVisual(_AMBER, _DARK, "diamond", 1.2, "solid", ""),
        ),
        ("Selected predicate", NodeVisual("#FFFFFF", _BLUE, "square-open", 3.5, "solid", "")),
        ("Excluded exact/clade", NodeVisual("#FFFFFF", _RED, "square-open", 3.5, "dashed", "")),
        ("Ancestor of represented", NodeVisual("#D1D5DB", _DARK, "square", 1.2, "solid", "")),
    )


def _hex_rgb(*, value: str) -> tuple[float, float, float]:
    """Convert a validated hexadecimal colour to PDF RGB fractions."""

    if len(value) != 7 or not value.startswith("#"):
        raise InputValidationError(f"Unsupported renderer colour: {value}")
    try:
        channels = tuple(int(value[index:index + 2], 16) / 255 for index in (1, 3, 5))
    except ValueError as error:
        raise InputValidationError(f"Unsupported renderer colour: {value}") from error
    return channels  # type: ignore[return-value]


def _pdf_marker(*, x: float, y: float, visual: NodeVisual) -> tuple[str, ...]:
    """Return vector PDF drawing commands for one small accessible marker."""

    fill = _hex_rgb(value=visual.fill_colour)
    line = _hex_rgb(value=visual.outline_colour)
    prefix = (
        f"{fill[0]:.3f} {fill[1]:.3f} {fill[2]:.3f} rg "
        f"{line[0]:.3f} {line[1]:.3f} {line[2]:.3f} RG "
        f"{visual.outline_width:.2f} w "
        f"{'[3 2] 0 d' if visual.outline_style == 'dashed' else '[] 0 d'}"
    )
    size = 3.4
    if visual.shape == "diamond":
        path = (
            f"{x:.2f} {y + size:.2f} m {x + size:.2f} {y:.2f} l "
            f"{x:.2f} {y - size:.2f} l {x - size:.2f} {y:.2f} l h B"
        )
    elif visual.shape.startswith("circle"):
        control = size * 0.55228475
        path = (
            f"{x + size:.2f} {y:.2f} m "
            f"{x + size:.2f} {y + control:.2f} {x + control:.2f} "
            f"{y + size:.2f} {x:.2f} {y + size:.2f} c "
            f"{x - control:.2f} {y + size:.2f} {x - size:.2f} "
            f"{y + control:.2f} {x - size:.2f} {y:.2f} c "
            f"{x - size:.2f} {y - control:.2f} {x - control:.2f} "
            f"{y - size:.2f} {x:.2f} {y - size:.2f} c "
            f"{x + control:.2f} {y - size:.2f} {x + size:.2f} "
            f"{y - control:.2f} {x + size:.2f} {y:.2f} c h B"
        )
    else:
        path = f"{x - size:.2f} {y - size:.2f} {2 * size:.2f} {2 * size:.2f} re B"
    return prefix, path


def _pdf_text(*, value: str) -> str:
    """Escape one bounded PDF literal string."""

    return value.encode("latin-1", errors="replace").decode("latin-1").replace(
        "\\", "\\\\"
    ).replace("(", "\\(").replace(")", "\\)")


def _minimal_pdf(*, content: bytes) -> bytes:
    """Wrap a vector content stream in a deterministic standards-valid PDF."""

    objects = (
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 842 595] "
        b"/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        b"<< /Length " + str(len(content)).encode("ascii") + b" >>\nstream\n"
        + content
        + b"\nendstream",
    )
    output = io.BytesIO()
    output.write(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = [0]
    for index, payload in enumerate(objects, start=1):
        offsets.append(output.tell())
        output.write(f"{index} 0 obj\n".encode("ascii"))
        output.write(payload)
        output.write(b"\nendobj\n")
    xref = output.tell()
    output.write(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    output.write(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        output.write(f"{offset:010d} 00000 n \n".encode("ascii"))
    output.write(
        f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
        f"startxref\n{xref}\n%%EOF\n".encode("ascii")
    )
    return output.getvalue()

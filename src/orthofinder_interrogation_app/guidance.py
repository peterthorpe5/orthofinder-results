"""Consistent biological interpretation guidance for every application graph."""

from __future__ import annotations

from dataclasses import dataclass

from orthofinder_results.errors import InputValidationError

try:
    import streamlit as st
except ModuleNotFoundError:  # pragma: no cover - checked by dependency smoke tests.
    st = None  # type: ignore[assignment]


@dataclass(frozen=True)
class GraphGuidance:
    """Plain-language explanation attached to one graph type."""

    title: str
    shows: str
    interpretation: str
    limitation: str


GRAPH_GUIDANCE = {
    "interactive_network": GraphGuidance(
        title="Interactive nearest-neighbour network",
        shows=(
            "Proteins as draggable nodes, coloured by species, with a small number of nearest "
            "distance neighbours retained as solid edges."
        ),
        interpretation=(
            "Look for proteins that repeatedly connect locally, peripheral nodes and disconnected "
            "regions. Use linked selection to follow the same proteins in quantitative views."
        ),
        limitation=(
            "Force-layout screen spacing is not a branch-length or distance scale, and connected "
            "components are not automatically biological subfamilies."
        ),
    ),
    "nearest_neighbour": GraphGuidance(
        title="Static nearest-neighbour topology",
        shows=(
            "The same retained local neighbour edges in a reproducible static force layout."
        ),
        interpretation=(
            "Use it to compare connectivity, isolates and local neighbourhoods when the draggable "
            "network is inconvenient."
        ),
        limitation=(
            "Only retained neighbour edges are shown and screen distance remains non-quantitative."
        ),
    ),
    "distance_distribution": GraphGuidance(
        title="Within-cluster distance distributions",
        shows=(
            "A histogram, violin-and-box summary and empirical cumulative distribution of every "
            "displayed pair distance, split into within- and between-species pairs."
        ),
        interpretation=(
            "A narrow left-shifted distribution suggests a compact displayed sample. Multiple "
            "peaks, long tails or separated within/between curves indicate heterogeneous spread."
        ),
        limitation=(
            "The distribution describes the analysed proteins and stated distance method only; "
            "large groups may be represented by a deterministic sample."
        ),
    ),
    "member_dispersion": GraphGuidance(
        title="Per-protein distance dispersion",
        shows=(
            "Each analysed protein's average distance to all other displayed proteins, grouped "
            "by species and ordered from central to peripheral."
        ),
        interpretation=(
            "Lower values identify relatively central proteins; high values flag peripheral or "
            "divergent proteins worth checking in the matrix and phylogram."
        ),
        limitation=(
            "Centrality is relative to the displayed sample and does not identify an ancestor or "
            "prove that a high-distance protein is erroneous."
        ),
    ),
    "medoid_distance": GraphGuidance(
        title="Distance from the sample medoid",
        shows=(
            "The exact distance from every analysed protein to the displayed protein with the "
            "smallest average distance to all others."
        ),
        interpretation=(
            "Short bars indicate proteins close to this practical sample centre; separated high "
            "bars can reveal peripheral branches or distinct distance scales."
        ),
        limitation=(
            "The medoid is an observed protein in the displayed sample, not a reconstructed "
            "ancestor or a guaranteed centre of the complete group."
        ),
    ),
    "species_pair_heatmap": GraphGuidance(
        title="Species-pair mean-distance heatmap",
        shows=(
            "Mean pairwise distance for every represented species pair among the analysed proteins."
        ),
        interpretation=(
            "Compare cells to find species pairs with unusually small or large mean distances and "
            "inspect the underlying pair counts before drawing conclusions."
        ),
        limitation=(
            "A mean can hide broad or multi-modal variation, and cells based on few proteins are "
            "less stable than densely represented cells."
        ),
    ),
    "pcoa_3d": GraphGuidance(
        title="Rotatable three-axis PCoA",
        shows=(
            "A three-dimensional approximation of all exact distances among analysed proteins."
        ),
        interpretation=(
            "Rotate and zoom to inspect whether a third axis separates points that overlap in 2D; "
            "read the 3D stress and distance agreement alongside the shape."
        ),
        limitation=(
            "PCoA is an approximation. Apparent arms, gaps and angles are not automatic "
            "subfamilies and require confirmation in exact distances and the phylogram."
        ),
    ),
    "pcoa_axes": GraphGuidance(
        title="Selectable two-axis PCoA",
        shows=(
            "Any chosen pair of the first three positive PCoA coordinate axes."
        ),
        interpretation=(
            "Compare alternative axis pairs to reveal structure compressed in axes 1 and 2, "
            "while checking how much positive inertia those axes retain."
        ),
        limitation=(
            "Choosing visually attractive axes does not improve the underlying fit and can "
            "overemphasise low-information dimensions."
        ),
    ),
    "pcoa_2d": GraphGuidance(
        title="Original two-axis PCoA diagnostic",
        shows=(
            "The first two principal-coordinate axes used by the original offline report."
        ),
        interpretation=(
            "Use species colours and linked highlights to locate broad distance patterns, then "
            "judge reliability from stress, distance agreement and retained inertia."
        ),
        limitation=(
            "Two-dimensional proximity can be distorted, especially when stress is high or the "
            "first two axes retain little positive inertia."
        ),
    ),
    "shepard": GraphGuidance(
        title="PCoA distance-fit (Shepard) plot",
        shows=(
            "Exact pair distance on one axis and the corresponding plotted PCoA distance on the "
            "other, with the ideal diagonal as reference."
        ),
        interpretation=(
            "Tight points near the diagonal indicate faithful distance representation. Broad "
            "scatter, curvature or systematic offsets identify distortion."
        ),
        limitation=(
            "The plotted points may be a deterministic diagnostic sample when the exact matrix "
            "contains too many pairs for a responsive graph."
        ),
    ),
    "phylogram": GraphGuidance(
        title="Branch-length gene-tree phylogram",
        shows=(
            "The gene-tree paths connecting analysed proteins, with horizontal length retaining "
            "the supplied branch-length scale."
        ),
        interpretation=(
            "Confirm whether network or PCoA separation follows long branches, shallow clades or "
            "single peripheral leaves."
        ),
        limitation=(
            "Vertical ordering is layout only, the tree may be pruned to analysed proteins, and "
            "gene-tree clades are not automatically equivalent to OrthoFinder HOG boundaries."
        ),
    ),
    "distance_heatmap": GraphGuidance(
        title="Exact pair-distance heatmap",
        shows=(
            "Every exact distance among the displayed proteins, ordered consistently along "
            "both axes."
        ),
        interpretation=(
            "Compact blocks of low values suggest locally similar proteins; warm bands or isolated "
            "cells reveal divergent subsets and outliers. Hover for exact values."
        ),
        limitation=(
            "Colour perception depends on the displayed scale and ordering; use the pair table for "
            "exact values and do not compare colours across independently scaled heatmaps."
        ),
    ),
    "comparison_summary": GraphGuidance(
        title="Cross-cluster compactness summary",
        shows=(
            "Average pair distance with distance spread for every selected cluster."
        ),
        interpretation=(
            "Clusters farther left are more compact on average under the same method; shorter "
            "spread indicators suggest more consistent distances."
        ),
        limitation=(
            "Means and SDs are not fully comparable when methods, sampling fractions or biological "
            "scope differ; inspect the provenance table."
        ),
    ),
    "comparison_violin": GraphGuidance(
        title="Cross-cluster violin distributions",
        shows=(
            "The complete displayed pair-distance distribution for each selected cluster."
        ),
        interpretation=(
            "Compare median, central width, tails and multiple modes rather than relying on one "
            "mean. Narrow, left-shifted violins indicate more compact displayed distances."
        ),
        limitation=(
            "Visual width represents density rather than protein count, and sampled matrices can "
            "contain different numbers of proteins and pairs."
        ),
    ),
    "comparison_ecdf": GraphGuidance(
        title="Cross-cluster cumulative distance distributions",
        shows=(
            "For each cluster, the fraction of displayed protein pairs at or below each distance."
        ),
        interpretation=(
            "A curve that rises earlier and remains farther left generally represents smaller pair "
            "distances. Crossings indicate that compactness depends on the chosen quantile."
        ),
        limitation=(
            "Curves should be compared only with attention to distance method, calculation scope "
            "and sample size."
        ),
    ),
    "comparison_pcoa": GraphGuidance(
        title="Small-multiple cluster PCoA comparison",
        shows=(
            "A separate two-axis PCoA diagnostic for each selected cluster using consistent "
            "species colours."
        ),
        interpretation=(
            "Compare within-panel shape, overlap and fit diagnostics to identify clusters with "
            "different internal distance structure."
        ),
        limitation=(
            "Panels have independent coordinate systems: absolute positions, rotations, axis "
            "ranges and angles must not be compared between clusters."
        ),
    ),
    "benchmark_distribution": GraphGuidance(
        title="Biological-profile dispersion distributions",
        shows=(
            "One point per OrthoFinder cluster, grouped into E3, E3 subtype, "
            "housekeeping-candidate, R/NLR-candidate or matched non-focus profiles."
        ),
        interpretation=(
            "Compare medians, overlap and cluster-to-cluster variability. On the matched-"
            "residual scale, positive values mean more dispersed than that cluster's own "
            "structurally matched controls; negative values mean more compact."
        ),
        limitation=(
            "The marker panels are hypotheses and reference candidates, not known compact "
            "or dispersed truths. A cluster can belong to more than one biological profile."
        ),
    ),
    "benchmark_contrast": GraphGuidance(
        title="FDR-controlled profile contrast forest plot",
        shows=(
            "Median differences between planned biological profiles after each cluster is "
            "centred on its own matched controls, with deterministic 95% bootstrap intervals."
        ),
        interpretation=(
            "Positive estimates indicate greater target dispersion than the reference; "
            "negative estimates indicate greater compactness. Red points have BH-FDR q≤0.05."
        ),
        limitation=(
            "Intervals describe cluster-level sampling variation. Matching reduces measured "
            "size/copy/species confounding but cannot remove unmeasured confounding."
        ),
    ),
    "benchmark_individual": GraphGuidance(
        title="Individual cluster versus empirical backgrounds",
        shows=(
            "The selected cluster and the median of each eligible background, joined by a "
            "line. Hover text supplies percentile, empirical p value, FDR q value and sample size."
        ),
        interpretation=(
            "Read the selected cluster relative to its own matched controls first, then compare "
            "its matched residual with pooled E3, housekeeping, R/NLR and E3 subtype backgrounds."
        ),
        limitation=(
            "Individual empirical tests can be coarse for small backgrounds. Leave-one-out is "
            "used when the selected cluster belongs to the comparison profile."
        ),
    ),
}


def graph_guidance(*, key: str) -> GraphGuidance:
    """Return validated guidance for one registered graph key.

    Args:
        key: Stable graph registry key.

    Returns:
        Immutable plain-language guidance.

    Raises:
        InputValidationError: If the graph key is unknown.
    """

    try:
        return GRAPH_GUIDANCE[key]
    except KeyError as error:
        raise InputValidationError(f"Unknown graph guidance key: {key}") from error


def render_graph_guidance(*, key: str) -> None:
    """Render a consistent expandable interpretation section for one graph.

    Args:
        key: Stable graph registry key.

    Raises:
        RuntimeError: If Streamlit is unavailable.
        InputValidationError: If the graph key is unknown.
    """

    if st is None:
        raise RuntimeError("Streamlit is required to render graph guidance.")
    guidance = graph_guidance(key=key)
    with st.expander(f"How to read this graph: {guidance.title}", expanded=False):
        st.markdown(
            f"**What this graph shows**  \n{guidance.shows}\n\n"
            f"**How to interpret it**  \n{guidance.interpretation}\n\n"
            f"**Important limitation**  \n{guidance.limitation}"
        )

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
    patterns: str
    confirmation: str
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
        patterns=(
            "Dense local neighbourhoods suggest short retained neighbour distances; isolates or "
            "small components indicate that their nearest retained links do not join the "
            "main network."
        ),
        confirmation=(
            "Confirm any apparent separation in the exact distance heatmap and branch-length "
            "phylogram; use the PCoA only after checking its fit diagnostics."
        ),
        limitation=(
            "Force-layout screen spacing is not a branch-length or distance scale, and connected "
            "components are not automatically biological subfamilies."
        ),
    ),
    "nearest_neighbour": GraphGuidance(
        title="Static nearest-neighbour topology",
        shows=("The same retained local neighbour edges in a reproducible static force layout."),
        interpretation=(
            "Use it to compare connectivity, isolates and local neighbourhoods when the draggable "
            "network is inconvenient."
        ),
        patterns=(
            "A connected core with peripheral leaves suggests a shared local neighbourhood plus "
            "more distant members; several components indicate missing retained cross-links."
        ),
        confirmation=(
            "Use the pair-distance table to inspect the edges and the phylogram to determine "
            "whether components follow gene-tree branches."
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
        patterns=(
            "A right-shift indicates larger typical distances; broad or multi-modal shapes "
            "indicate "
            "that one centre value does not describe every pair well."
        ),
        confirmation=(
            "Compare mean with median and SD, then locate the responsible proteins or species in "
            "the per-protein table, heatmap and phylogram."
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
        patterns=(
            "A smooth range suggests graded divergence; a small high-distance set can indicate a "
            "long branch, divergent paralogues, alignment problems or identifier issues."
        ),
        confirmation=(
            "Select the high-distance proteins and inspect their exact pair rows, species identity "
            "and placement on the branch-length tree."
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
        patterns=(
            "Many similar bars indicate a broadly even radius around the medoid; a long tail "
            "identifies proteins much farther from that observed centre."
        ),
        confirmation=(
            "Check whether the same long-tail proteins are peripheral in the member-dispersion "
            "view and lie on long branches in the phylogram."
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
        patterns=(
            "Cool diagonal cells can reflect similar within-species copies; warmer off-diagonal "
            "blocks can reflect deeper species or paralogue structure."
        ),
        confirmation=(
            "Use the downloadable species-pair table for counts and exact means, then inspect the "
            "protein-level heatmap when a cell is based on mixed distances."
        ),
        limitation=(
            "A mean can hide broad or multi-modal variation, and cells based on few proteins are "
            "less stable than densely represented cells."
        ),
    ),
    "pcoa_3d": GraphGuidance(
        title="Rotatable three-axis PCoA",
        shows=("A three-dimensional approximation of all exact distances among analysed proteins."),
        interpretation=(
            "Rotate and zoom to inspect whether a third axis separates points that overlap in 2D; "
            "read the 3D stress and distance agreement alongside the shape."
        ),
        patterns=(
            "Stable separation across several viewing angles is more credible than a gap visible "
            "from one angle, but neither is a formal cluster assignment."
        ),
        confirmation=(
            "Check the 3D fit metrics and Shepard plot, then confirm the proteins' exact distances "
            "and tree positions."
        ),
        limitation=(
            "PCoA is an approximation. Apparent arms, gaps and angles are not automatic "
            "subfamilies and require confirmation in exact distances and the phylogram."
        ),
    ),
    "pcoa_axes": GraphGuidance(
        title="Selectable two-axis PCoA",
        shows=("Any chosen pair of the first three positive PCoA coordinate axes."),
        interpretation=(
            "Compare alternative axis pairs to reveal structure compressed in axes 1 and 2, "
            "while checking how much positive inertia those axes retain."
        ),
        patterns=(
            "Structure appearing only on a low-inertia axis may be real but explains a small part "
            "of the positive coordinate signal."
        ),
        confirmation=(
            "Report the selected axes and their inertia, and verify any apparent group in the "
            "exact heatmap and phylogram."
        ),
        limitation=(
            "Choosing visually attractive axes does not improve the underlying fit and can "
            "overemphasise low-information dimensions."
        ),
    ),
    "pcoa_2d": GraphGuidance(
        title="Original two-axis PCoA diagnostic",
        shows=("The first two principal-coordinate axes used by the original offline report."),
        interpretation=(
            "Use species colours and linked highlights to locate broad distance patterns, then "
            "judge reliability from stress, distance agreement and retained inertia."
        ),
        patterns=(
            "Tight clouds indicate similar plotted distances; arms, gaps or overlap may reflect "
            "real distance structure or projection distortion."
        ),
        confirmation=(
            "Read the quality badge and Shepard plot before using geometry, then return to the "
            "exact matrix and branch-length tree."
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
        patterns=(
            "Points above the diagonal are expanded in the projection; points below it are "
            "compressed. A fan shape indicates error that changes with distance scale."
        ),
        confirmation=(
            "Use this plot with the numerical stress, correlation and retained-inertia values; "
            "do not judge PCoA fidelity from a single diagnostic alone."
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
        patterns=(
            "Long terminal branches identify divergent leaves; long internal branches separate "
            "sets of leaves; shallow branching indicates comparatively small tree-path distances."
        ),
        confirmation=(
            "Hover or use the heatmap and pair table for exact distances, and retain the tree "
            "authority and pruning status when reporting the pattern."
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
        patterns=(
            "Symmetric blocks suggest subsets with small internal distances; a uniformly warm row "
            "or column identifies one protein distant from most others."
        ),
        confirmation=(
            "Use the downloadable pair table for exact values and the phylogram to determine "
            "whether a block follows a coherent branch."
        ),
        limitation=(
            "Colour perception depends on the displayed scale and ordering; use the pair table for "
            "exact values and do not compare colours across independently scaled heatmaps."
        ),
    ),
    "comparison_summary": GraphGuidance(
        title="Cross-cluster compactness summary",
        shows=("Average pair distance with distance spread for every selected cluster."),
        interpretation=(
            "Clusters farther left are more compact on average under the same method; shorter "
            "spread indicators suggest more consistent distances."
        ),
        patterns=(
            "Separated centres indicate different typical divergence; overlapping spread bars "
            "show that the within-cluster distance ranges remain similar or broad."
        ),
        confirmation=(
            "Inspect the provenance table first, then compare full violin and ECDF distributions "
            "rather than ranking clusters from means alone."
        ),
        limitation=(
            "Means and SDs are not fully comparable when methods, sampling fractions or biological "
            "scope differ; inspect the provenance table."
        ),
    ),
    "comparison_violin": GraphGuidance(
        title="Cross-cluster violin distributions",
        shows=("The complete displayed pair-distance distribution for each selected cluster."),
        interpretation=(
            "Compare median, central width, tails and multiple modes rather than relying on one "
            "mean. Narrow, left-shifted violins indicate more compact displayed distances."
        ),
        patterns=(
            "Wide sections mark dense distance ranges; long tails flag uncommon large values; "
            "several bulges suggest multiple distance scales."
        ),
        confirmation=(
            "Use the ECDF for quantile-wise comparison and inspect each cluster's heatmap and tree "
            "when distributions overlap or have several modes."
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
        patterns=(
            "Parallel separated curves support a broad shift; strong crossings show that one "
            "cluster has more small distances but also a heavier tail."
        ),
        confirmation=(
            "Compare medians, IQRs and tails in the table and violin plot; do not reduce crossing "
            "curves to a single compactness ordering."
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
        patterns=(
            "Different apparent shapes can reflect genuine internal structure or unequal PCoA fit; "
            "similar shapes do not establish comparable branch-length scales."
        ),
        confirmation=(
            "Check each panel's fit measures and compare its full distance distribution and "
            "phylogram before describing a structural difference."
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
            "composition-matched controls; negative values mean more compact."
        ),
        patterns=(
            "Separated medians suggest a profile-level shift, whereas broad overlap indicates "
            "substantial cluster-to-cluster variation even when centres differ."
        ),
        confirmation=(
            "Use the planned forest plot for effect size, interval and BH-FDR inference; profile "
            "boxplots alone are descriptive."
        ),
        limitation=(
            "The marker panels are hypotheses and reference candidates, not known compact "
            "or dispersed truths. A cluster can belong to more than one biological profile."
        ),
    ),
    "benchmark_marker_coverage": GraphGuidance(
        title="Profile-marker cluster coverage",
        shows=(
            "Each bar is one profile-defining E3 seed or Arabidopsis reference marker. "
            "Bar length is the number of distinct OrthoFinder clusters that contained an "
            "exact matched protein."
        ),
        interpretation=(
            "Longer bars identify markers represented in several retained clusters. Hover to "
            "see matched-protein and sampled-species counts, then use the table to inspect the "
            "exact group and protein identifiers."
        ),
        patterns=(
            "A marker in one cluster has a one-to-one coverage path in this authority; several "
            "clusters can reflect hierarchy, duplication, identifier mapping or profile overlap."
        ),
        confirmation=(
            "Inspect the complete marker table and selected cluster membership before treating "
            "multiple bars or matches as biological expansion."
        ),
        limitation=(
            "This is mapping coverage, not expression, copy number or evolutionary dispersion. "
            "Only enabled authority markers matched in this exact dataset appear."
        ),
    ),
    "benchmark_classification": GraphGuidance(
        title="Matched-control cluster map",
        shows=(
            "One point per biological target cluster. Its horizontal position is the "
            "average-distance percentile and its vertical position is the distance-spread "
            "percentile relative to that cluster's own matched non-focus controls."
        ),
        interpretation=(
            "Points on the left are comparatively compact and points on the right are "
            "comparatively dispersed. Low points have comparatively uniform distances; high "
            "points are more heterogeneous. Hover to identify the cluster and its profiles."
        ),
        patterns=(
            "Upper-right points combine high central divergence and high spread; lower-right "
            "points are uniformly divergent; upper-left points mix compact centres with long tails."
        ),
        confirmation=(
            "Open a point in the cluster explorer and inspect its raw distribution, "
            "sampling scope, "
            "heatmap and phylogram before assigning a biological explanation."
        ),
        limitation=(
            "Percentiles are calibrated separately against each cluster's small matched-control "
            "set. With three controls they are deliberately coarse and are descriptive labels, "
            "not FDR-controlled profile tests."
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
        patterns=(
            "An interval wholly on one side of zero supports a directional effect estimate; an "
            "interval crossing zero shows that the compatible range includes no median difference."
        ),
        confirmation=(
            "Read target and reference cluster counts, Cliff's delta and the BH-FDR q "
            "value together; "
            "then inspect the underlying profile distributions for overlap and outliers."
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
        patterns=(
            "A long connecting line indicates a large difference from the background median; a "
            "high percentile places the cluster near the dispersed end of that background."
        ),
        confirmation=(
            "Check background size, leave-one-out status, empirical p value and FDR q value, then "
            "open the cluster's exact distance and tree views."
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
            f"**What common result patterns mean**  \n{guidance.patterns}\n\n"
            f"**What to check next**  \n{guidance.confirmation}\n\n"
            f"**Important limitation**  \n{guidance.limitation}"
        )

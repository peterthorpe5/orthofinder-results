"""In-application methods, glossary and result-interpretation guidance."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

from orthofinder_results.errors import InputValidationError

from .exports import render_table_downloads
from .models import ResourceIdentity

try:
    import streamlit as st
except ModuleNotFoundError:  # pragma: no cover - checked by dependency smoke tests.
    st = None  # type: ignore[assignment]


@dataclass(frozen=True)
class PageGuidance:
    """Plain-language orientation for one application page."""

    title: str
    purpose: str
    results: str
    read_order: str
    caution: str


@dataclass(frozen=True)
class GlossaryEntry:
    """One scientific or application term and its interpretation."""

    term: str
    topic: str
    definition: str
    interpretation: str
    related_terms: str


@dataclass(frozen=True)
class MethodStep:
    """One reproducible stage in the analysis and viewer workflow."""

    number: int
    title: str
    method: str
    rationale: str
    output: str
    caution: str


PAGE_GUIDANCE = {
    "overview": PageGuidance(
        title="Dataset summary",
        purpose=(
            "Summarises the imported OrthoFinder run, the retained group authorities and "
            "which analytical capabilities are stored in the opened resource."
        ),
        results=(
            "Counts describe records in this resource: groups, sampled species, hierarchy "
            "collections, persisted distance results and portable trees."
        ),
        read_order=(
            "Confirm the run identity first, then inspect the headline counts and choose the "
            "biological question that best matches the analysis you want to perform."
        ),
        caution=(
            "A HOG retained at several hierarchy levels contributes several group records. "
            "Do not sum those rows and interpret the result as unique biological families."
        ),
    ),
    "focus_clusters": PageGuidance(
        title="Focus protein clusters",
        purpose=(
            "Maps a reviewed focus-protein authority to exact proteins and groups in this "
            "OrthoFinder run; the packaged E3 catalogue is the project default."
        ),
        results=(
            "Each row is a cluster containing at least one matched focus protein, with the "
            "matched identifiers and retained provenance shown explicitly."
        ),
        read_order=(
            "Review the authority and matching summary, inspect exact matched proteins, then "
            "open a cluster to examine its complete membership and evolutionary distances."
        ),
        caution=(
            "A focus match prioritises a cluster for inspection. It does not transfer E3 or "
            "other focus function automatically to every paralogue and orthologue in the group."
        ),
    ),
    "protein_search": PageGuidance(
        title="Gene and protein search",
        purpose=(
            "Finds an exact or literal partial identifier and reports every retained HOG or "
            "legacy-orthogroup membership containing that protein."
        ),
        results=(
            "Rows are membership records rather than unique proteins: the same protein can occur "
            "at several HOG hierarchy levels and in a linked legacy orthogroup."
        ),
        read_order=(
            "Prefer exact search, confirm the canonical protein and species, then choose the "
            "specific group-and-level record required for tree or distance analysis."
        ),
        caution=(
            "Contains search is literal and bounded. Gene symbols or descriptions are found "
            "only when they occur in the identifiers stored by the source run."
        ),
    ),
    "group_search": PageGuidance(
        title="Group search",
        purpose=(
            "Filters the complete group authority by identifier, membership, species presence, "
            "copy number, size and stored-distance state."
        ),
        results=(
            "Each row is one unambiguous run, group-system, hierarchy-node and group-ID record; "
            "the table reports complete membership statistics and any stored distance summary."
        ),
        read_order=(
            "Choose the group system and level, apply biological filters, review the number of "
            "matches, then inspect one group before comparing or exporting it."
        ),
        caution=(
            "Not calculated means that no stored distance result exists for that record; it "
            "never means a biological distance of zero. Paginated downloads contain that page."
        ),
    ),
    "cluster_explorer": PageGuidance(
        title="Single-cluster evolutionary views",
        purpose=(
            "Combines membership, exact or bounded patristic distances, dispersion summaries, "
            "PCoA diagnostics, a branch-length phylogram and sparse topology for one group."
        ),
        results=(
            "The complementary views describe central divergence, heterogeneity, peripheral "
            "proteins, species structure and the gene-tree paths underlying the distance matrix."
        ),
        read_order=(
            "Start with calculation provenance, inspect distance distributions, judge PCoA fit, "
            "then confirm apparent patterns in the phylogram, exact heatmap and pair table."
        ),
        caution=(
            "Force-layout spacing is non-quantitative and PCoA is an approximation. Large groups "
            "may use a declared deterministic member sample rather than every protein."
        ),
    ),
    "cluster_comparison": PageGuidance(
        title="Cross-cluster comparison",
        purpose=(
            "Places two to twelve selected clusters beside one another using their stored or "
            "on-demand distance summaries, distributions and PCoA diagnostics."
        ),
        results=(
            "Tables and figures show which displayed clusters are more compact, more variable or "
            "differently shaped under their stated distance calculations."
        ),
        read_order=(
            "Check method and sampling comparability, compare mean and spread, inspect complete "
            "distributions, then use each cluster's tree and exact matrix for confirmation."
        ),
        caution=(
            "The selected clusters are not independent experimental replicates by default. "
            "This page is descriptive and does not apply matched-control or FDR inference."
        ),
    ),
    "dispersion_benchmarks": PageGuidance(
        title="Calibrated dispersion",
        purpose=(
            "Tests whether E3 and reference-profile clusters differ in evolutionary dispersion "
            "after each target is centred on composition-matched non-focus controls."
        ),
        results=(
            "The page reports cluster-level distributions, profile contrasts, effect sizes, "
            "confidence intervals, FDR values and descriptive individual-cluster classifications."
        ),
        read_order=(
            "Read the headline profile contrasts first, check the BH-FDR q value and interval, "
            "then inspect profile overlap and individual clusters without treating trends as tests."
        ),
        caution=(
            "Housekeeping and R/NLR panels are candidate comparators, not known truths. Protein "
            "pairs are not inferential replicates, and three-control percentiles are "
            "deliberately coarse."
        ),
    ),
    "all_distance_results": PageGuidance(
        title="Dataset-wide stored distance results",
        purpose=(
            "Collects one preferred persisted distance summary per eligible group for filtering, "
            "review and complete table export."
        ),
        results=(
            "Rows report the group identity, calculation method and scope, analysed proteins, "
            "pair count and descriptive distance statistics."
        ),
        read_order=(
            "Filter by authority, hierarchy level, method and status; choose informative columns; "
            "then export the complete filtered result rather than relying on the bounded preview."
        ),
        caution=(
            "Absence from this table means no preferred persisted result was available. Direct "
            "comparison still requires compatible distance methods and sampling scopes."
        ),
    ),
    "taxonomic_search": PageGuidance(
        title="Reviewed taxonomic search",
        purpose=(
            "Finds groups containing, enriched for, sampled-exclusive to or near-exclusive to a "
            "selected taxon using a reviewed mapping for this exact species-label set."
        ),
        results=(
            "Rows report sampled descendant coverage, mapped purity, outsiders, unresolved labels "
            "and, for enrichment, odds ratios and BH-adjusted q values."
        ),
        read_order=(
            "Audit the mapping first, choose the biological claim, inspect its sampled universe, "
            "then review outsiders and unresolved labels for every candidate group."
        ),
        caution=(
            "Sampled exclusivity applies only to this imported species set. Unmapped labels never "
            "become guessed descendants or reviewed outsiders."
        ),
    ),
    "selection_coverage": PageGuidance(
        title="Selection coverage tree",
        purpose=(
            "Builds a reviewed taxonomy tree showing how exact, clade-inclusion, only-in and "
            "exclusion predicates affect a bounded group-selection authority."
        ),
        results=(
            "The tree distinguishes represented, expected-but-unrepresented, ancestor, outsider "
            "and unresolved states; audit tables explain why every group passed or failed."
        ),
        read_order=(
            "Confirm the mapping and expected universe, compose selectors, inspect the visible "
            "tree state, then use the predicate and unmapped-label audits before accepting groups."
        ),
        caution=(
            "This is a taxonomy coverage tree, not the OrthoFinder species tree or a gene tree. "
            "Expected no data is not evidence of biological absence."
        ),
    ),
    "offline_report": PageGuidance(
        title="Offline report",
        purpose=(
            "Downloads the self-contained HTML snapshot published with the immutable resource."
        ),
        results=(
            "The report retains the bounded summaries and visual payload prepared when the "
            "resource "
            "was built, allowing review without running Streamlit."
        ),
        read_order=(
            "Confirm the resource identity, download the HTML, and open it locally in a modern "
            "browser; use the live app for newer on-demand or schema-aware views."
        ),
        caution=(
            "The offline report is a fixed publication snapshot. It does not acquire later cached "
            "distance calculations or viewer-only interface additions."
        ),
    ),
    "methods": PageGuidance(
        title="Methods and provenance",
        purpose=(
            "Explains the complete route from an already completed OrthoFinder run to the portable "
            "resource, calibrated statistics, app views and downloads."
        ),
        results=(
            "The page separates upstream OrthoFinder inference, resource construction, distance "
            "analysis, calibration, visualisation and read-only interrogation."
        ),
        read_order=(
            "Start with the opened-resource identity, follow the numbered workflow, then use the "
            "downloadable methods table when preparing a manuscript or analysis record."
        ),
        caution=(
            "The viewer does not rerun OrthoFinder or alter the formal resource. Exact upstream "
            "alignment, tree and reconciliation settings remain properties of the source run."
        ),
    ),
    "glossary": PageGuidance(
        title="Scientific and application glossary",
        purpose=(
            "Defines the biological, phylogenetic, statistical, visualisation and data-provenance "
            "terms used throughout the app."
        ),
        results=(
            "Each entry provides a concise definition, an interpretation for this app and related "
            "terms that help connect results across pages."
        ),
        read_order=(
            "Search a word or filter by topic; retain All topics when exploring unfamiliar terms; "
            "export the complete dictionary for offline use."
        ),
        caution=(
            "Definitions describe this package's data contract and displays. They do not replace "
            "the original OrthoFinder method papers or dataset-specific biological validation."
        ),
    ),
}


GLOSSARY_ENTRIES = (
    GlossaryEntry(
        "Analysis cache",
        "Data and provenance",
        "A user-writable sidecar holding bounded distance analyses calculated on demand.",
        "It improves repeat viewing but is never written inside the completed resource.",
        "Immutable resource; checksum; computation status",
    ),
    GlossaryEntry(
        "Authority",
        "Data and provenance",
        "A versioned input or table treated as the named source for identifiers or "
        "classifications.",
        "Authority does not imply biological truth; it means the source and version are explicit.",
        "Focus authority; taxonomy authority; provenance",
    ),
    GlossaryEntry(
        "Benjamini–Hochberg FDR",
        "Statistics",
        "A procedure controlling the expected false-discovery proportion across a test family.",
        "Use the adjusted q value—not only the unadjusted p value—when assessing planned "
        "contrasts.",
        "FDR q value; p value; test family",
    ),
    GlossaryEntry(
        "Biological profile",
        "Calibrated dispersion",
        "A named collection of clusters defined by E3 seeds or Arabidopsis comparison markers.",
        "Profiles support cluster-level comparisons; one cluster may occur in several profiles.",
        "E3 profile; housekeeping candidate; R/NLR candidate",
    ),
    GlossaryEntry(
        "Bootstrap confidence interval",
        "Statistics",
        "An interval obtained by repeatedly resampling clusters with replacement.",
        "The app reports deterministic 95% intervals around median profile differences.",
        "Median difference; statistical replicate",
    ),
    GlossaryEntry(
        "Bounded analysis",
        "Distances and trees",
        "An analysis limited to a declared maximum number of deterministically selected proteins.",
        "It keeps visualisation and pair counts tractable while retaining the stated "
        "sampling scope.",
        "Deterministic member sample; sampling fraction",
    ),
    GlossaryEntry(
        "Canonical protein identifier",
        "Groups and proteins",
        "The exact sequence or protein label retained in OrthoFinder membership tables.",
        "Use it for reproducible searches and exports; aliases are accepted only when unambiguous.",
        "Internal sequence ID; focus protein",
    ),
    GlossaryEntry(
        "Central divergence class",
        "Calibrated dispersion",
        "A descriptive compact, typical or dispersed label based on matched-control percentiles.",
        "It summarises average pair distance and is not an FDR-controlled profile result.",
        "Average pair distance; matched controls; percentile",
    ),
    GlossaryEntry(
        "Checksum",
        "Data and provenance",
        "A cryptographic digest used to detect any change to a source or published payload.",
        "Matching checksums support integrity; they do not assess biological correctness.",
        "Manifest; immutable resource; provenance",
    ),
    GlossaryEntry(
        "Cliff's delta",
        "Statistics",
        "A non-parametric effect size comparing the probability that values from one "
        "group exceed another.",
        "Its sign follows target minus reference; magnitude describes separation, not "
        "significance.",
        "Mann–Whitney test; effect size; median difference",
    ),
    GlossaryEntry(
        "Coefficient of variation",
        "Distances and trees",
        "Population standard deviation divided by the mean when the mean is positive.",
        "It expresses relative spread but becomes unstable near a mean of zero.",
        "Population SD; mean pair distance; heterogeneity",
    ),
    GlossaryEntry(
        "Comparison scale",
        "Calibrated dispersion",
        "The raw cluster statistic or its difference from that cluster's matched-control median.",
        "Matched residuals improve comparability; raw values retain the original "
        "branch-length scale.",
        "Raw statistic; matched residual",
    ),
    GlossaryEntry(
        "Complete distance matrix",
        "Distances and trees",
        "All pair distances among every protein included in the stated analysis set.",
        "Complete refers to analysed members; a bounded sample can have a complete matrix "
        "of its sample.",
        "Protein pair; bounded analysis; exact matrix",
    ),
    GlossaryEntry(
        "Computation status",
        "Data and provenance",
        "A code stating whether distances are exact, sampled, calculated lazily or unavailable.",
        "Always read it with analysed-member counts and method before comparing results.",
        "Distance method; sampling scope; unavailable",
    ),
    GlossaryEntry(
        "Copy number",
        "Groups and proteins",
        "The number of proteins from one sampled species assigned to a specified group record.",
        "It is group- and hierarchy-level-specific and is not a genome-wide gene-family estimate.",
        "Represented species; single-copy species fraction",
    ),
    GlossaryEntry(
        "Deterministic member sample",
        "Distances and trees",
        "A reproducible subset selected when a group exceeds the declared distance-member limit.",
        "Repeated runs with the same authority and settings retain the same members; "
        "required markers are forced in.",
        "Bounded analysis; sampling fraction",
    ),
    GlossaryEntry(
        "Distance agreement",
        "Visualisation",
        "Pearson correlation between exact pair distances and distances in a PCoA display.",
        "Values nearer one indicate better agreement, but must be read with stress and "
        "retained inertia.",
        "PCoA; stress; Shepard plot",
    ),
    GlossaryEntry(
        "Distance interquartile range",
        "Distances and trees",
        "The difference between the 75th and 25th percentiles of pair distances.",
        "It describes central spread and is less influenced by extreme pairs than "
        "standard deviation.",
        "Population SD; heterogeneity",
    ),
    GlossaryEntry(
        "Distance method",
        "Distances and trees",
        "The declared procedure used to obtain pairwise values, such as resolved-tree "
        "patristic distance.",
        "Only compare numerical values directly when methods and calculation scopes are "
        "compatible.",
        "Patristic distance; computation status",
    ),
    GlossaryEntry(
        "DuckDB",
        "Data and provenance",
        "The portable analytical database containing typed versions of the published "
        "result tables.",
        "The app opens it read-only and uses bounded queries rather than loading every "
        "row at once.",
        "Parquet; TSV; schema",
    ),
    GlossaryEntry(
        "E3 profile",
        "Calibrated dispersion",
        "A pooled or subclassed set of HOGs containing proteins from the reviewed E3 seed "
        "catalogue.",
        "Membership prioritises relevant clusters; it does not annotate every member as "
        "an E3 ligase.",
        "Biological profile; focus authority",
    ),
    GlossaryEntry(
        "Effect size",
        "Statistics",
        "A measure of the direction and magnitude of a difference, independent of "
        "sample-size significance.",
        "Interpret Cliff's delta and median difference alongside uncertainty and FDR q values.",
        "Cliff's delta; median difference; confidence interval",
    ),
    GlossaryEntry(
        "Empirical cumulative distribution (ECDF)",
        "Visualisation",
        "The fraction of observed pair distances at or below each distance value.",
        "A curve that rises earlier generally represents smaller distances; crossings "
        "indicate distributional differences.",
        "Distance distribution; pair distance",
    ),
    GlossaryEntry(
        "Empirical p value",
        "Statistics",
        "A finite-sample tail probability calculated from the observed background ranks.",
        "With small backgrounds it has few possible values and should not be over-interpreted.",
        "Percentile; leave-one-out; FDR q value",
    ),
    GlossaryEntry(
        "Enrichment odds ratio",
        "Taxonomy",
        "A presence-based association between group membership and a reviewed target taxon.",
        "Values above one indicate enrichment in the sampled universe; use the adjusted q "
        "value for inference.",
        "Taxonomic enrichment; FDR q value",
    ),
    GlossaryEntry(
        "Exact species set",
        "Taxonomy",
        "A filter requiring precisely the selected sampled species and no additional "
        "species labels.",
        "It uses exact labels from this run and does not infer ancestry from spelling.",
        "Represented species; taxonomy authority",
    ),
    GlossaryEntry(
        "FDR q value",
        "Statistics",
        "A p value adjusted within a declared multiple-testing family using Benjamini–Hochberg.",
        "A threshold such as q≤0.05 must be chosen before interpreting a result as "
        "FDR-significant.",
        "Benjamini–Hochberg FDR; p value",
    ),
    GlossaryEntry(
        "Focus authority",
        "Groups and proteins",
        "A versioned list of proteins used to find and prioritise clusters of interest.",
        "The packaged authority contains E3 seeds, but a reviewed replacement TSV can be supplied.",
        "Authority; focus protein; E3 profile",
    ),
    GlossaryEntry(
        "Force-directed layout",
        "Visualisation",
        "A display that moves connected nodes to make a sparse network easier to inspect.",
        "Connectivity is informative, but screen spacing and angles are not evolutionary "
        "distances.",
        "Nearest-neighbour edge; sparse topology",
    ),
    GlossaryEntry(
        "Gene tree",
        "Distances and trees",
        "A phylogenetic tree describing relationships among homologous gene or protein sequences.",
        "The displayed tree may be pruned to analysed proteins while retaining connecting "
        "branch lengths.",
        "Resolved gene tree; phylogram; species tree",
    ),
    GlossaryEntry(
        "Gene-tree parent clade",
        "Groups and proteins",
        "The OrthoFinder gene-tree clade recorded as the source of a hierarchical "
        "orthogroup membership.",
        "It provides provenance for the HOG record and is not a descriptive gene annotation.",
        "HOG; hierarchy node; gene tree",
    ),
    GlossaryEntry(
        "Group record",
        "Groups and proteins",
        "One group identified by run ID, group system, hierarchy node and group ID.",
        "Use the complete composite identity because the same textual ID can occur in "
        "different contexts.",
        "HOG; legacy orthogroup; run ID",
    ),
    GlossaryEntry(
        "Heterogeneity class",
        "Calibrated dispersion",
        "A descriptive low, typical or high label based on the distance-SD "
        "matched-control percentile.",
        "It concerns variation among pair distances, not the central divergence of the cluster.",
        "Population SD; central divergence class",
    ),
    GlossaryEntry(
        "Hierarchical orthogroup (HOG)",
        "Groups and proteins",
        "A set of genes descended from one ancestral gene at a stated species-tree node.",
        "Related HOGs can be retained at several hierarchy levels, so level is part of "
        "the identity.",
        "Hierarchy node; orthogroup; legacy orthogroup",
    ),
    GlossaryEntry(
        "Hierarchy node",
        "Groups and proteins",
        "The named species-tree level at which an OrthoFinder HOG is defined.",
        "N0 and other node labels are identifiers from the source run, not rank names.",
        "HOG; species tree",
    ),
    GlossaryEntry(
        "Housekeeping-reference candidate",
        "Calibrated dispersion",
        "An Arabidopsis marker selected as a plausible stable or core-cellular comparison panel.",
        "It is a hypothesis-bearing comparator, not a guarantee that its HOG is "
        "evolutionarily compact.",
        "Biological profile; R/NLR-reference candidate",
    ),
    GlossaryEntry(
        "Immutable resource",
        "Data and provenance",
        "A completed, versioned output directory that the viewer is not permitted to modify.",
        "New analyses use a sidecar cache or a new run ID; formal published tables remain "
        "unchanged.",
        "Analysis cache; manifest; run ID",
    ),
    GlossaryEntry(
        "Individual cluster comparison",
        "Calibrated dispersion",
        "Placement of one selected cluster against its controls and eligible biological "
        "backgrounds.",
        "Read its own matched controls first, then compare matched residuals with pooled profiles.",
        "Empirical p value; leave-one-out; matched residual",
    ),
    GlossaryEntry(
        "Internal sequence ID",
        "Groups and proteins",
        "An OrthoFinder-specific identifier mapped to a canonical sequence label when "
        "SequenceIDs.txt is available.",
        "The app resolves it for searching but retains the canonical authority in results.",
        "Canonical protein identifier",
    ),
    GlossaryEntry(
        "Leave-one-out comparison",
        "Statistics",
        "A comparison that removes the selected cluster from a background to which it belongs.",
        "It prevents the test item from contributing to its own reference distribution.",
        "Individual cluster comparison; empirical p value",
    ),
    GlossaryEntry(
        "Legacy orthogroup",
        "Groups and proteins",
        "A flat group imported from OrthoFinder's Orthogroups.tsv authority.",
        "It has no HOG hierarchy node and is labelled ROOT in composite display identities.",
        "Orthogroup; HOG; group record",
    ),
    GlossaryEntry(
        "Mann–Whitney U test",
        "Statistics",
        "A rank-based test comparing two independent collections of cluster-level values.",
        "The benchmark uses a tie-corrected two-sided form and adjusts planned tests by FDR.",
        "Cliff's delta; FDR q value; statistical replicate",
    ),
    GlossaryEntry(
        "Manifest",
        "Data and provenance",
        "The completed run record containing identity, schema, versions, counts and "
        "publication status.",
        "The app checks it against DuckDB before presenting the resource as valid.",
        "Checksum; immutable resource; schema",
    ),
    GlossaryEntry(
        "Marker",
        "Calibrated dispersion",
        "A named reference gene or protein used to assign a biological profile to a matched HOG.",
        "Marker-to-cluster coverage is mapping provenance, not an expression or dispersion result.",
        "Biological profile; marker coverage",
    ),
    GlossaryEntry(
        "Matched control",
        "Calibrated dispersion",
        "A non-focus cluster selected without using distance and matched on group "
        "composition variables.",
        "Matching uses protein count, represented species, mean copies and "
        "single-copy-species fraction.",
        "Matched residual; composition-matched non-focus cluster",
    ),
    GlossaryEntry(
        "Matched residual",
        "Calibrated dispersion",
        "A target cluster statistic minus the median statistic of its own matched controls.",
        "Positive values mean more dispersed than controls; negative values mean more compact.",
        "Matched control; comparison scale",
    ),
    GlossaryEntry(
        "Mean pair distance",
        "Distances and trees",
        "The arithmetic mean of all pair distances in the analysed matrix.",
        "Smaller values usually indicate central compactness under the same method and scope.",
        "Median pair distance; population SD",
    ),
    GlossaryEntry(
        "Median difference",
        "Statistics",
        "The target profile's median cluster residual minus the reference profile's "
        "median residual.",
        "Positive values favour greater target dispersion; uncertainty is shown by the "
        "bootstrap interval.",
        "Matched residual; bootstrap confidence interval",
    ),
    GlossaryEntry(
        "Median pair distance",
        "Distances and trees",
        "The middle pair-distance value after sorting all analysed pairs.",
        "It describes central divergence with less sensitivity to extreme distances than the mean.",
        "Mean pair distance; distance distribution",
    ),
    GlossaryEntry(
        "Medoid",
        "Distances and trees",
        "The analysed protein with the smallest mean distance to all other analysed proteins.",
        "It is a practical sample centre, not a reconstructed ancestor or necessarily a "
        "consensus sequence.",
        "Mean pair distance; central protein",
    ),
    GlossaryEntry(
        "Nearest-neighbour edge",
        "Visualisation",
        "A retained connection from one protein to one of its closest partners by exact distance.",
        "Sparse edges make local topology inspectable but omit most pair relationships.",
        "Force-directed layout; exact distance matrix",
    ),
    GlossaryEntry(
        "Newick",
        "Distances and trees",
        "A compact text format encoding tree topology and, when present, branch lengths.",
        "The resource stores compressed, checksum-bound Newick payloads for portable tree "
        "analysis.",
        "Gene tree; phylogram; checksum",
    ),
    GlossaryEntry(
        "Non-Euclidean signal",
        "Visualisation",
        "The negative-inertia contribution arising when a distance matrix is not "
        "perfectly Euclidean.",
        "Larger values warn that low-dimensional PCoA coordinates cannot faithfully "
        "represent all distances.",
        "PCoA; stress; positive inertia",
    ),
    GlossaryEntry(
        "OrthoFinder",
        "Groups and proteins",
        "Upstream software that infers orthogroups, hierarchical orthogroups, gene trees "
        "and orthology relationships.",
        "This package interrogates a completed run; it does not reproduce the upstream "
        "inference interactively.",
        "Orthogroup; HOG; resolved gene tree",
    ),
    GlossaryEntry(
        "Orthogroup",
        "Groups and proteins",
        "Genes descended from a single gene in the last common ancestor of the sampled species.",
        "The app distinguishes flat legacy orthogroups from node-specific hierarchical "
        "orthogroups.",
        "Legacy orthogroup; HOG",
    ),
    GlossaryEntry(
        "Pair distance",
        "Distances and trees",
        "The declared evolutionary or sequence distance between two analysed proteins.",
        "For resolved-tree patristic analysis it is the sum of branch lengths joining the leaves.",
        "Patristic distance; protein pair",
    ),
    GlossaryEntry(
        "Parquet",
        "Data and provenance",
        "A typed columnar file format published for efficient programmatic analysis.",
        "Use it or DuckDB for large complete tables; use TSV for transparent exchange.",
        "DuckDB; TSV",
    ),
    GlossaryEntry(
        "Patristic distance",
        "Distances and trees",
        "The sum of branch lengths along the tree path connecting two leaves.",
        "Values inherit the branch-length units and assumptions of the source resolved gene tree.",
        "Resolved gene tree; pair distance",
    ),
    GlossaryEntry(
        "Percentile",
        "Statistics",
        "The empirical rank of a selected value within a stated background distribution.",
        "High percentiles indicate relatively large values; three controls yield only "
        "coarse ranks.",
        "Empirical p value; matched control",
    ),
    GlossaryEntry(
        "Phylogram",
        "Visualisation",
        "A tree drawing in which horizontal branch length retains the source evolutionary scale.",
        "Use it to confirm whether apparent distance patterns follow long branches or clades.",
        "Gene tree; patristic distance; cladogram",
    ),
    GlossaryEntry(
        "Principal coordinates analysis (PCoA)",
        "Visualisation",
        "A low-dimensional coordinate approximation of all pair distances in a matrix.",
        "Treat apparent clusters as hypotheses and judge fit with inertia, stress, "
        "agreement and the Shepard plot.",
        "Positive inertia; stress; distance agreement",
    ),
    GlossaryEntry(
        "Population standard deviation",
        "Distances and trees",
        "The root-mean-square deviation of all analysed pair distances from their mean.",
        "It measures distance heterogeneity; it is not a standard error of independent "
        "protein pairs.",
        "Heterogeneity; coefficient of variation",
    ),
    GlossaryEntry(
        "Positive inertia",
        "Visualisation",
        "The positive coordinate-space contribution used to describe how much PCoA "
        "structure axes retain.",
        "Variation shown is diagnostic geometry, not a biological variance estimate.",
        "PCoA; non-Euclidean signal",
    ),
    GlossaryEntry(
        "Protein pair",
        "Distances and trees",
        "One unordered pair of distinct analysed proteins and its associated distance.",
        "Pairs within a cluster are correlated and are not treated as inferential replicates.",
        "Pair distance; statistical replicate",
    ),
    GlossaryEntry(
        "Provenance",
        "Data and provenance",
        "Information identifying the source file, method, version, sampling scope and "
        "transformation of a result.",
        "Retain it when exporting or citing a result so its analytical meaning remains "
        "recoverable.",
        "Authority; manifest; checksum",
    ),
    GlossaryEntry(
        "R/NLR-reference candidate",
        "Calibrated dispersion",
        "An Arabidopsis resistance or NLR marker used as a plausible comparison panel.",
        "It tests a biological hypothesis and is not assumed in advance to be dispersed.",
        "Biological profile; housekeeping-reference candidate",
    ),
    GlossaryEntry(
        "Raw statistic",
        "Calibrated dispersion",
        "The cluster's uncentred distance summary on the source branch-length scale.",
        "Use it for direct description; use matched residuals for calibrated profile comparison.",
        "Matched residual; comparison scale",
    ),
    GlossaryEntry(
        "Represented species",
        "Groups and proteins",
        "A sampled species contributing at least one protein to the selected group.",
        "Absent species are not included in mean copies per represented species.",
        "Copy number; sampled universe",
    ),
    GlossaryEntry(
        "Resolved gene tree",
        "Distances and trees",
        "An OrthoFinder gene tree after reconciliation or resolution used as the "
        "preferred tree authority.",
        "Patristic values depend on its topology and branch lengths; the app does not "
        "re-estimate them.",
        "Gene tree; patristic distance",
    ),
    GlossaryEntry(
        "Resource schema",
        "Data and provenance",
        "The versioned contract defining required tables, columns and capabilities of a "
        "completed resource.",
        "Viewer version and resource schema are distinct; compatible viewers can open "
        "older schemas read-only.",
        "Immutable resource; viewer version",
    ),
    GlossaryEntry(
        "Run ID",
        "Data and provenance",
        "The immutable identifier for one imported OrthoFinder analysis and its published "
        "resource.",
        "Use a new run ID when species or source results change; never merge identities silently.",
        "Manifest; group record",
    ),
    GlossaryEntry(
        "Sampled exclusivity",
        "Taxonomy",
        "A claim that no reviewed represented outsider occurs within this imported "
        "species universe.",
        "It is not evidence that the group is universally absent from all unsampled taxa.",
        "Near-exclusive; sampled universe",
    ),
    GlossaryEntry(
        "Sampling fraction",
        "Distances and trees",
        "Analysed protein count divided by the complete group membership count.",
        "A value below one signals that the displayed distance analysis is a bounded sample.",
        "Deterministic member sample; bounded analysis",
    ),
    GlossaryEntry(
        "Selection coverage tree",
        "Taxonomy",
        "A reviewed taxonomy display of selector state and dataset coverage.",
        "It is neither the OrthoFinder species tree nor a cluster gene tree.",
        "Taxonomy authority; expected-no-data",
    ),
    GlossaryEntry(
        "Shepard plot",
        "Visualisation",
        "A scatter plot of exact input distances against distances in a PCoA display.",
        "Points close to the diagonal indicate faithful representation; scatter or "
        "curvature shows distortion.",
        "PCoA; distance agreement; stress",
    ),
    GlossaryEntry(
        "Single-copy species fraction",
        "Groups and proteins",
        "The fraction of represented species contributing exactly one protein to a group.",
        "It is one composition variable used to match non-focus controls without "
        "examining distance.",
        "Copy number; matched control",
    ),
    GlossaryEntry(
        "Species-pair heatmap",
        "Visualisation",
        "A matrix of mean distances for every represented pair of species.",
        "Check underlying pair counts because a mean can conceal broad variation or sparse cells.",
        "Pair distance; represented species",
    ),
    GlossaryEntry(
        "Species tree",
        "Distances and trees",
        "The phylogeny of sampled species used by OrthoFinder to define hierarchy nodes "
        "and reconcile gene trees.",
        "Do not confuse it with a gene tree or the separate taxonomy coverage tree.",
        "Hierarchy node; gene tree; selection coverage tree",
    ),
    GlossaryEntry(
        "Statistical replicate",
        "Statistics",
        "The independent observational unit used by an inferential test.",
        "Dispersion profile tests use clusters as replicates; correlated protein pairs "
        "remain within-cluster observations.",
        "Protein pair; Mann–Whitney test",
    ),
    GlossaryEntry(
        "Stress",
        "Visualisation",
        "A normalised measure of distortion between exact pair distances and plotted PCoA "
        "distances.",
        "Lower is better, but it must be considered with retained inertia, agreement and "
        "the Shepard plot.",
        "PCoA; distance agreement; Shepard plot",
    ),
    GlossaryEntry(
        "Composition-matched non-focus cluster",
        "Calibrated dispersion",
        "A control cluster matched on group size, species breadth and copy-number "
        "structure, not 3D protein structure.",
        "It provides a local empirical baseline selected without using the outcome "
        "distance values.",
        "Matched control; matched residual",
    ),
    GlossaryEntry(
        "Taxonomy authority",
        "Taxonomy",
        "A reviewed mapping from exact source species labels to taxon identifiers and ancestry.",
        "Pending, ambiguous and unmapped labels remain explicit and are never guessed from names.",
        "Sampled exclusivity; selection coverage tree",
    ),
    GlossaryEntry(
        "Test family",
        "Statistics",
        "The declared collection of related hypotheses adjusted together for multiple testing.",
        "FDR q values are interpretable only with the family in which adjustment was performed.",
        "Benjamini–Hochberg FDR; FDR q value",
    ),
    GlossaryEntry(
        "TSV",
        "Data and provenance",
        "A tab-separated plain-text table with one record per line.",
        "It is the transparent exchange format used by downloads; identifiers are not "
        "comma-separated.",
        "Excel; Parquet; DuckDB",
    ),
    GlossaryEntry(
        "Unavailable",
        "Data and provenance",
        "An explicit state indicating that a value could not be calculated or was not retained.",
        "It is never converted to zero and should be inspected with its status or reason field.",
        "Computation status; missing value",
    ),
    GlossaryEntry(
        "Analysed proteins",
        "Distances and trees",
        "The exact full or bounded set of proteins included in a displayed distance calculation.",
        "Use this denominator, rather than complete group size, when interpreting pairs and plots.",
        "Complete group; sampling fraction; bounded analysis",
    ),
    GlossaryEntry(
        "Complete group",
        "Groups and proteins",
        "Every protein assigned to the selected OrthoFinder group record before distance sampling.",
        "Its size can exceed the analysed set shown in quadratic distance visualisations.",
        "Analysed proteins; group record; sampling fraction",
    ),
    GlossaryEntry(
        "Expected no data",
        "Taxonomy",
        "A reviewed expected taxon that is not represented in the imported dataset.",
        "It records dataset coverage and is never evidence that a gene or function is "
        "biologically absent.",
        "Selection coverage tree; sampled universe",
    ),
    GlossaryEntry(
        "Formatted Excel",
        "Data and provenance",
        "An XLSX export containing exact displayed rows, filters, typed formats and a "
        "definitions sheet.",
        "Use it for review and reporting; TSV remains the transparent plain-text exchange "
        "authority.",
        "TSV; Parquet; DuckDB",
    ),
    GlossaryEntry(
        "Interactive network",
        "Visualisation",
        "A draggable sparse nearest-neighbour view with hover labels and linked selections.",
        "Download HTML to retain interaction or use the static topology PDF for reporting.",
        "Force-directed layout; nearest-neighbour edge",
    ),
    GlossaryEntry(
        "Mapped purity",
        "Taxonomy",
        "The fraction of reviewed represented species that belong to the selected target taxon.",
        "Unresolved labels are reported separately and do not become reviewed target or "
        "outside species.",
        "Target descendant coverage; sampled exclusivity",
    ),
    GlossaryEntry(
        "Marker coverage",
        "Calibrated dispersion",
        "The number of exact matched OrthoFinder clusters or proteins associated with a marker.",
        "It describes mapping success and multiplicity, not expression or evolutionary dispersion.",
        "Marker; biological profile",
    ),
    GlossaryEntry(
        "Near-exclusive",
        "Taxonomy",
        "A sampled restriction allowing only a declared number of reviewed outsiders or "
        "unresolved labels.",
        "Every permitted exception remains visible; the claim is not universal across "
        "unsampled taxa.",
        "Sampled exclusivity; mapped purity",
    ),
    GlossaryEntry(
        "PDF figure export",
        "Data and provenance",
        "A static PDF generated on demand from the plotted figure for reporting or manuscript use.",
        "Interactive selections are not embedded; WebGL traces can be rasterised inside the PDF.",
        "Interactive network; provenance",
    ),
    GlossaryEntry(
        "ROOT",
        "Groups and proteins",
        "The display label used when a flat legacy orthogroup has no HOG hierarchy node.",
        "It prevents an empty hierarchy field being mistaken for missing group identity.",
        "Legacy orthogroup; hierarchy node",
    ),
    GlossaryEntry(
        "Sampled universe",
        "Taxonomy",
        "The exact species labels imported into the current OrthoFinder resource and "
        "reviewed mapping.",
        "Coverage, enrichment and exclusivity conclusions are bounded to this declared universe.",
        "Represented species; sampled exclusivity",
    ),
    GlossaryEntry(
        "Sampling scope",
        "Distances and trees",
        "The declaration of whether a result covers every group member or a deterministic "
        "bounded subset.",
        "Compare cluster statistics only after checking that their methods and sampling "
        "scopes are compatible.",
        "Computation status; sampling fraction",
    ),
    GlossaryEntry(
        "Sparse topology",
        "Visualisation",
        "A network retaining only a small number of nearest-neighbour edges from the full "
        "pair matrix.",
        "It reveals local connectivity while deliberately omitting most quantitative pair "
        "relationships.",
        "Nearest-neighbour edge; exact distance matrix",
    ),
    GlossaryEntry(
        "Target descendant coverage",
        "Taxonomy",
        "Represented reviewed target descendants divided by all reviewed target "
        "descendants sampled in the run.",
        "It describes breadth within the sampled target universe, not abundance or "
        "universal presence.",
        "Mapped purity; sampled universe",
    ),
    GlossaryEntry(
        "Tree authority",
        "Distances and trees",
        "The exact source tree identifier and type from which topology or patristic "
        "values were obtained.",
        "Retain it with exported results because alternative trees can yield different distances.",
        "Resolved gene tree; provenance",
    ),
    GlossaryEntry(
        "Viewer version",
        "Data and provenance",
        "The installed application/package release presenting an immutable resource.",
        "It can add explanations and compatible views without changing the resource "
        "schema or results.",
        "Resource schema; immutable resource",
    ),
)


METHOD_STEPS = (
    MethodStep(
        1,
        "Completed OrthoFinder analysis",
        (
            "The workflow begins with a completed OrthoFinder result directory. OrthoFinder "
            "performed the upstream orthogroup, HOG, species-tree, gene-tree and orthology "
            "inference."
        ),
        "Using completed outputs preserves the original evolutionary authority and avoids "
        "hidden re-analysis.",
        "Versioned OrthoFinder files plus the source run metadata.",
        "Alignment, tree-inference and reconciliation settings belong to the source "
        "OrthoFinder run.",
    ),
    MethodStep(
        2,
        "Version-aware discovery and validation",
        (
            "orthofinder-results detected the supported OrthoFinder layout, selected the relevant "
            "adapter, inventoried inputs and rejected missing or contradictory authorities."
        ),
        "Failing before publication prevents a partially understood directory becoming an "
        "analytical resource.",
        "Validated layout, source inventory, adapter identity and checksums.",
        "No filename or species relationship is guessed when an authority is absent or ambiguous.",
    ),
    MethodStep(
        3,
        "Normalised group and membership authority",
        (
            "Legacy orthogroups and every available HOG hierarchy level were converted to "
            "long-form "
            "membership, group-statistic and group-by-species copy-number tables."
        ),
        "A composite key prevents group IDs from being silently equated across runs, "
        "systems or hierarchy levels.",
        "TSV, Parquet and DuckDB relations keyed by run, group system, hierarchy node and "
        "group ID.",
        "Rows from related HOG levels are related records, not automatically unique families.",
    ),
    MethodStep(
        4,
        "Portable tree authority",
        (
            "Preferred resolved gene trees were normalised and stored as checksum-bound compressed "
            "Newick payloads, with tree inventory and node records retained where available."
        ),
        "Portable trees permit later bounded analysis without retaining access to the "
        "original large results directory.",
        "Verified tree payloads, inventories and optional tree-node tables.",
        "A tree remains unavailable when no valid leaf mapping or usable branch-length "
        "authority exists.",
    ),
    MethodStep(
        5,
        "Patristic distances and bounded sampling",
        (
            "Pairwise evolutionary distances were calculated as the summed branch lengths "
            "connecting "
            "each pair of analysed leaves. Groups above the declared limit used a "
            "deterministic sample "
            "of at most 250 proteins, with requested focus or marker proteins retained."
        ),
        "The bound controls quadratic pair growth while preserving reproducibility and "
        "target representation.",
        "Pair-distance rows, calculation status, analysed-member count and sampling fraction.",
        "Results describe the analysed members; sampled groups are not presented as "
        "full-group matrices.",
    ),
    MethodStep(
        6,
        "Cluster-level dispersion summaries",
        (
            "Each eligible matrix was summarised using mean and median pair distance, population "
            "standard deviation, interquartile range and coefficient of variation when the "
            "mean was positive."
        ),
        "Centre and spread answer different questions and should not be collapsed into one score.",
        "One descriptive record per cluster, metric, method and sampling scope.",
        "Correlated protein pairs are within-cluster observations, not independent "
        "inferential replicates.",
    ),
    MethodStep(
        7,
        "E3 and reference-profile mapping",
        (
            "Versioned E3 seeds, housekeeping-reference candidates and R/NLR-reference candidates "
            "were matched to exact proteins and HOGs, retaining marker identity and source "
            "provenance."
        ),
        "Named profiles make biological hypotheses testable while keeping their evidence "
        "auditable.",
        "Marker audit, marker matches and cluster-to-profile membership tables.",
        "A marker match does not transfer its function automatically to every cluster member.",
    ),
    MethodStep(
        8,
        "Distance-blind matched controls",
        (
            "Each biological target received three unique non-focus control clusters "
            "selected without "
            "examining distance. Matching considered protein count, represented-species "
            "count, mean "
            "copies per species and single-copy-species fraction."
        ),
        "Local composition matching reduces obvious group-size, breadth and "
        "copy-structure confounding.",
        "Auditable target-to-control assignments and matching scores.",
        "Matching cannot remove unmeasured confounding and does not match "
        "three-dimensional protein structure.",
    ),
    MethodStep(
        9,
        "Matched residuals and planned profile contrasts",
        (
            "For each metric, the target statistic was centred by subtracting the median "
            "of its controls. "
            "Profiles were compared using median residual differences, deterministic "
            "1,000-resample "
            "bootstrap 95% intervals, Cliff's delta and tie-corrected two-sided Mann–Whitney tests."
        ),
        "Clusters, rather than the many dependent protein pairs, form the inferential units.",
        "Profile contrast records containing sample sizes, effects, intervals, p values "
        "and status.",
        "Reference profiles are hypotheses to test, not labels that force an expected direction.",
    ),
    MethodStep(
        10,
        "Multiple testing and individual-cluster placement",
        (
            "Benjamini–Hochberg correction was applied within declared metric-specific "
            "test families. "
            "Each target was also ranked against its controls and eligible pooled profiles using "
            "finite-sample empirical tails and leave-one-out where it belonged to the background."
        ),
        "FDR controls planned family-wise interpretation; leave-one-out prevents self-comparison.",
        "Adjusted contrast tables, individual comparisons and transparent percentile "
        "classifications.",
        "Three-control ranks and empirical p values are coarse; classifications are descriptive.",
    ),
    MethodStep(
        11,
        "Complementary visual diagnostics",
        (
            "The viewer combines exact tables and heatmaps with distributions, sparse "
            "nearest-neighbour "
            "topology, PCoA plus Shepard diagnostics, and a branch-length phylogram."
        ),
        "No single display captures central divergence, heterogeneity, topology and tree "
        "structure.",
        "Interactive figures, linked selections and manuscript-ready PDF exports.",
        "Force layouts are non-quantitative and PCoA geometry must be judged by its fit "
        "diagnostics.",
    ),
    MethodStep(
        12,
        "Reviewed taxonomy and coverage",
        (
            "Taxonomic searches and coverage trees use a separately reviewed mapping from "
            "exact run "
            "labels to taxon IDs. Unmapped, pending and ambiguous labels remain explicit "
            "and fail closed "
            "where a strict only-in claim requires complete placement."
        ),
        "Explicit review prevents spelling or naming conventions being mistaken for ancestry.",
        "Taxonomy audits, sampled-universe searches and reproducible coverage-tree packages.",
        "Exclusivity is bounded to the sampled and reviewed dataset, never universal absence.",
    ),
    MethodStep(
        13,
        "Immutable publication and read-only interrogation",
        (
            "The pipeline published compressed TSV, typed Parquet, DuckDB, manifest, QC "
            "and offline "
            "report outputs through checksum validation and guarded atomic publication. "
            "The Streamlit "
            "application validates and opens the completed DuckDB read-only."
        ),
        "Separating construction from interrogation protects the formal analytical record.",
        "One versioned resource plus optional user-writable analysis cache and viewer downloads.",
        "Viewer-only releases can add explanation or compatible exports but cannot change "
        "stored statistics.",
    ),
)


_GLOSSARY_HELP = {
    "Term": "Term exactly as used in the application.",
    "Topic": "Broad scientific or technical area.",
    "Definition": "Concise meaning within the orthofinder-results data contract.",
    "How to interpret it here": "Practical guidance for reading this application's results.",
    "Related terms": "Other glossary entries that provide useful context.",
}
_METHOD_HELP = {
    "Step": "Stable order in the resource and analysis workflow.",
    "Stage": "Name of the analytical or publication stage.",
    "What was done": "Method implemented by OrthoFinder or orthofinder-results.",
    "Why it was done": "Scientific or reproducibility purpose of the stage.",
    "Published result": "Primary output or authority produced by the stage.",
    "Important qualification": "Boundary that must be retained when interpreting the output.",
}
_RESOURCE_HELP = {
    "Property": "Resource identity or capability field.",
    "Value": "Value validated for the resource currently open in the app.",
    "Meaning": "Why the field matters for reproducible interpretation.",
}


def page_guidance(*, key: str) -> PageGuidance:
    """Return validated guidance for one application page.

    Args:
        key: Stable page-guidance key.

    Returns:
        Immutable page guidance.

    Raises:
        InputValidationError: If the key is unknown.
    """

    try:
        return PAGE_GUIDANCE[key]
    except KeyError as error:
        raise InputValidationError(f"Unknown page guidance key: {key}") from error


def render_page_guidance(*, key: str) -> None:
    """Render an expandable result-orientation panel.

    Args:
        key: Stable page-guidance key.

    Raises:
        RuntimeError: If Streamlit is unavailable.
        InputValidationError: If the key is unknown.
    """

    if st is None:
        raise RuntimeError("Streamlit is required to render page guidance.")
    guidance = page_guidance(key=key)
    with st.expander(f"How to read this page: {guidance.title}", expanded=False):
        st.markdown(
            f"**Purpose**  \n{guidance.purpose}\n\n"
            f"**What the results show**  \n{guidance.results}\n\n"
            f"**Recommended reading order**  \n{guidance.read_order}\n\n"
            f"**Important qualification**  \n{guidance.caution}"
        )


def filter_glossary_entries(
    *,
    entries: Sequence[GlossaryEntry],
    query: str,
    topic: str,
) -> tuple[GlossaryEntry, ...]:
    """Filter glossary entries by literal text and exact topic.

    Args:
        entries: Candidate glossary entries.
        query: Case-insensitive literal search text.
        topic: Exact topic or ``All topics``.

    Returns:
        Alphabetically ordered matching entries.

    Raises:
        InputValidationError: If a requested topic is unavailable.
    """

    available = {entry.topic for entry in entries}
    if topic != "All topics" and topic not in available:
        raise InputValidationError(f"Unknown glossary topic: {topic}")
    needle = query.strip().casefold()
    matches = []
    for entry in entries:
        if topic != "All topics" and entry.topic != topic:
            continue
        searchable = " ".join(
            (
                entry.term,
                entry.topic,
                entry.definition,
                entry.interpretation,
                entry.related_terms,
            )
        ).casefold()
        if not needle or needle in searchable:
            matches.append(entry)
    return tuple(sorted(matches, key=lambda entry: entry.term.casefold()))


def glossary_records(*, entries: Sequence[GlossaryEntry]) -> tuple[dict[str, Any], ...]:
    """Return user-facing records for glossary display and export.

    Args:
        entries: Glossary entries to transform.

    Returns:
        Ordered table records.
    """

    return tuple(
        {
            "Term": entry.term,
            "Topic": entry.topic,
            "Definition": entry.definition,
            "How to interpret it here": entry.interpretation,
            "Related terms": entry.related_terms,
        }
        for entry in entries
    )


def method_records(*, steps: Sequence[MethodStep] = METHOD_STEPS) -> tuple[dict[str, Any], ...]:
    """Return user-facing records for method display and export.

    Args:
        steps: Ordered method stages.

    Returns:
        Ordered table records.

    Raises:
        InputValidationError: If step numbers are not unique and increasing.
    """

    numbers = tuple(step.number for step in steps)
    if numbers != tuple(sorted(set(numbers))):
        raise InputValidationError("Method step numbers must be unique and increasing.")
    return tuple(
        {
            "Step": step.number,
            "Stage": step.title,
            "What was done": step.method,
            "Why it was done": step.rationale,
            "Published result": step.output,
            "Important qualification": step.caution,
        }
        for step in steps
    )


def resource_method_records(*, resource: ResourceIdentity) -> tuple[dict[str, str], ...]:
    """Describe the validated identity and capabilities of an opened resource.

    Args:
        resource: Validated immutable resource identity.

    Returns:
        Records suitable for display and paired export.
    """

    benchmark_relations = {
        "benchmark_group_profiles",
        "benchmark_cluster_results",
        "benchmark_contrasts",
        "benchmark_individual_comparisons",
        "benchmark_cluster_classifications",
    }
    values = (
        (
            "Run ID",
            resource.run_id,
            "Immutable identity of the imported OrthoFinder analysis.",
        ),
        (
            "Source OrthoFinder version",
            resource.orthofinder_version,
            "Version reported by the upstream completed run.",
        ),
        (
            "Parser adapter",
            resource.adapter_name,
            "Version-aware layout reader selected during resource construction.",
        ),
        (
            "Resource schema",
            str(resource.schema_version),
            "Portable data-contract version validated by the viewer.",
        ),
        (
            "Resource builder version",
            resource.resource_package_version,
            "orthofinder-results release that published the immutable resource.",
        ),
        (
            "Primary group authority",
            resource.primary_group_authority,
            "Preferred HOG or legacy group collection recorded by the run.",
        ),
        (
            "Portable tree payloads",
            "Available" if "tree_payloads" in resource.relations else "Unavailable",
            "Determines whether bounded tree distances can be calculated lazily.",
        ),
        (
            "Calibrated dispersion relations",
            "Available" if benchmark_relations.issubset(resource.relations) else "Unavailable",
            "Determines whether profile contrasts and individual calibration can be shown.",
        ),
        (
            "Application access",
            "Read-only",
            "The formal DuckDB and completed resource are never modified by the viewer.",
        ),
    )
    return tuple(
        {"Property": property_name, "Value": value, "Meaning": meaning}
        for property_name, value, meaning in values
    )


def render_methods_page(*, resource: ResourceIdentity) -> None:
    """Render the complete analysis methods and opened-resource provenance.

    Args:
        resource: Validated immutable resource identity.

    Raises:
        RuntimeError: If Streamlit is unavailable.
    """

    if st is None:
        raise RuntimeError("Streamlit is required to render the methods page.")
    st.header("Methods and provenance")
    st.write(
        "This page separates what OrthoFinder inferred upstream from what "
        "orthofinder-results calculated, published and displays. It describes the full "
        "scientific route without modifying the resource currently open in the app."
    )
    render_page_guidance(key="methods")

    st.subheader("Opened resource")
    resource_rows = resource_method_records(resource=resource)
    st.dataframe(
        resource_rows,
        width="stretch",
        hide_index=True,
        column_config={
            label: st.column_config.Column(label=label, help=description)
            for label, description in _RESOURCE_HELP.items()
        },
    )
    render_table_downloads(
        records=resource_rows,
        file_stem=f"{resource.run_id}_methods_resource_identity",
        key="methods_resource_identity_download",
        tsv_label="Download resource identity as TSV",
        excel_label="Download resource identity as formatted Excel",
        column_definitions=_RESOURCE_HELP,
        workbook_title="OrthoFinder resource methods identity",
    )

    st.subheader("Workflow from OrthoFinder outputs to the app")
    st.caption(
        "Open each numbered stage for the method, rationale, published output and the "
        "qualification that should accompany interpretation."
    )
    for step in METHOD_STEPS:
        with st.expander(f"{step.number}. {step.title}", expanded=step.number == 1):
            st.markdown(
                f"**What was done**  \n{step.method}\n\n"
                f"**Why it was done**  \n{step.rationale}\n\n"
                f"**Published result**  \n{step.output}\n\n"
                f"**Important qualification**  \n{step.caution}"
            )

    method_rows = method_records()
    with st.expander("Complete methods table and downloads", expanded=False):
        st.dataframe(
            method_rows,
            width="stretch",
            hide_index=True,
            column_config={
                label: st.column_config.Column(label=label, help=description)
                for label, description in _METHOD_HELP.items()
            },
        )
        render_table_downloads(
            records=method_rows,
            file_stem="orthofinder_interrogation_methods",
            key="methods_complete_download",
            tsv_label="Download methods as TSV",
            excel_label="Download methods as formatted Excel",
            column_definitions=_METHOD_HELP,
            workbook_title="OrthoFinder interrogation methods",
        )

    st.subheader("Interpretive boundary")
    st.info(
        "The completed resource is the analytical authority. The app performs read-only queries "
        "and, where supported, bounded calculations in a separate user cache. Viewer-only "
        "updates can improve explanations and exports but cannot change stored biological results."
    )


def render_glossary_page() -> None:
    """Render a searchable and fully downloadable application glossary.

    Raises:
        RuntimeError: If Streamlit is unavailable.
    """

    if st is None:
        raise RuntimeError("Streamlit is required to render the glossary page.")
    st.header("Glossary")
    st.write(
        "Search the terms used in controls, tables, figures and methods. Definitions are "
        "specific to how orthofinder-results stores and presents the analysis."
    )
    render_page_guidance(key="glossary")
    controls = st.columns(2)
    query = controls[0].text_input(
        "Search glossary",
        value="",
        placeholder="For example: HOG, FDR, residual or PCoA",
        help="Literal case-insensitive search across terms, definitions and related terms.",
    )
    topics = ("All topics",) + tuple(sorted({entry.topic for entry in GLOSSARY_ENTRIES}))
    topic = controls[1].selectbox(
        "Topic",
        topics,
        help="Restrict results to one scientific or technical area.",
    )
    entries = filter_glossary_entries(
        entries=GLOSSARY_ENTRIES,
        query=query,
        topic=topic,
    )
    st.caption(f"{len(entries):,} of {len(GLOSSARY_ENTRIES):,} terms shown.")
    if not entries:
        st.info("No glossary terms match the current search and topic.")
        return
    records = glossary_records(entries=entries)
    st.dataframe(
        records,
        width="stretch",
        hide_index=True,
        column_config={
            label: st.column_config.Column(label=label, help=description)
            for label, description in _GLOSSARY_HELP.items()
        },
    )
    render_table_downloads(
        records=records,
        file_stem="orthofinder_interrogation_glossary",
        key="glossary_download",
        tsv_label="Download displayed glossary as TSV",
        excel_label="Download displayed glossary as formatted Excel",
        column_definitions=_GLOSSARY_HELP,
        workbook_title="OrthoFinder interrogation glossary",
    )

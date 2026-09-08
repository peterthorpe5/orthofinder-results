# Changelog

## Unreleased

## 0.4.1 - 2026-09-08

- Make Summary the guided opening page, place biological questions and routes
  before technical group-collection tables, and explain schema-2 versus
  schema-3 distance capability in plain language.
- Increase sidebar navigation type and replace navigation, control, metric and
  tab labels that assumed knowledge of internal implementation terms.
- Add contextual `?` help to principal controls, summary measures and displayed
  table headings while preserving stable raw field names in TSV downloads.
- Translate stored-distance states for display, distinguish complete group size
  from analysed sample size, and explain mean, median, population SD and sample
  medoid scope beside the relevant views.
- Add task-specific guidance to force networks, multi-view distance dispersion,
  two- and three-dimensional PCoA, Shepard diagnostics, gene-tree phylograms,
  exact heatmaps, protein-pair tables and cross-cluster comparisons.
- Replace the dense help page with an expandable getting-started guide and
  scientific glossary covering HOG hierarchy, species rules, distance sampling,
  visual interpretation, reviewed taxonomy, exclusivity claims and provenance.
- State the boundary between this generic standalone app and the separate E3
  prioritisation workflow, while retaining nested-HOG, split/merge and cross-run
  comparison as explicit generic development goals.
- Extend unit and Streamlit end-to-end tests for the new landing page, readable
  display adapters, navigation, help content and unchanged raw authorities.

## 0.4.0 - 2026-09-07

- Bump the physical resource contract to schema 3 and publish one compressed,
  checksum-bound portable Newick payload per preferred resolved or fallback
  original gene tree.
- Calculate exact bounded patristic matrices lazily for selected groups without
  modifying the completed resource, including exact HOG membership restriction,
  canonical/internal/species-prefixed leaf resolution and explicit failures.
- Add content-addressed, gzip-compressed, atomically published sidecar caches
  keyed by run, composite group identity, tree checksum and calculation controls;
  default to persistent macOS or XDG user cache locations rather than `/tmp`.
- Preserve schema-2 pilot matrices and schema-3 persisted DuckDB pairs as
  higher-priority distance authorities for backward compatibility.
- Restore a self-contained draggable force-directed view with separately
  toggleable layout-only connectors and retain the static nearest-neighbour view.
- Add rotatable 3D and selectable-axis PCoA alongside the original 2D diagnostic,
  with separate 2D/3D inertia, stress and distance-correlation reporting.
- Add within/between histograms, violins and empirical CDFs; per-member
  mean/median/SD and nearest-neighbour summaries; medoid-distance profiles; and
  species-pair mean-distance heatmaps and TSV exports.
- Add filtered exact member-to-member distance tables and linked member/species
  selection across the cluster explorer.
- Add a 2–12-cluster comparison workspace for mean, population SD, median,
  complete displayed distributions and independent stable-colour PCoA panels.
- Redesign the overview and group-search flow around actionable find, explore and
  compare tasks while retaining `ANY`, `ALL`, `EXACT_SET` and rejected-species
  filters.
- Add the dataset-generic `orthofinder-taxonomy-map` command, which resolves the
  current resource's exact species labels against a local NCBI taxdump and emits
  only `PENDING_REVIEW`, `AMBIGUOUS` or `UNMAPPED` candidates until human review.
- Expand unit, integration and Streamlit end-to-end tests to cover the new
  portable backend, cache corruption, taxonomy authority, dispersion mathematics,
  visual repertoire, force-network defences and comparison workflow.

## 0.3.0 - 2026-09-07

- Port the bounded distance pilot into separate application pages for PCoA and
  conservative fit diagnostics, Shepard plots, branch-length phylograms, exact
  distance matrices and nearest-neighbour topology.
- Validate the schema-2 offline report's run-bound JSON payload as the bounded
  compatibility authority for pruned phylograms not persisted separately in
  older DuckDB resources.
- Link exact member and species selections across every evolutionary view and
  retain sampled-medoid, distance-method and displayed/full-group scope.
- Add a versioned reviewed-taxonomy TSV contract, audit table, review template,
  upload/sidecar workflows and explicit `REVIEWED`, `UNMAPPED`, `AMBIGUOUS` and
  missing states.
- Add NCBI taxon-ID/descendant searches for `contains`, `enriched`, `exclusive
  within sampled analysis` and `near-exclusive` semantics.
- Calculate enrichment with one-sided Fisher exact species-presence tests and
  Benjamini–Hochberg correction across the complete selected group authority;
  retain sampled outsiders and unresolved labels in results.
- Select the first available local Streamlit port from 8501 upwards when no
  exact `--server-port` is supplied.
- Add explicit Plotly, NetworkX and SciPy application dependencies and expand
  unit, integration and Streamlit end-to-end coverage.

## 0.2.0 - 2026-09-07

- Add the first independent Streamlit application package and the fully named
  `orthofinder-interrogation-app` launcher.
- Validate completed resource directories and direct DuckDB files before opening
  them, require schema 2, and keep every application query read-only.
- Add complete-authority overview counts and bounded lazy group searches.
- Add exact species filters for any selected species, every selected species or
  exactly the selected species set, plus explicit rejection of groups containing
  any selected species.
- Add group/member identifier searches, group-size and species-breadth bounds,
  distance-availability controls and deterministic sorting.
- Join existing mean, median and population-SD distance summaries to group search
  results while retaining calculation method and sampled/full status.
- Add selected-group copy-count and complete-membership views with TSV downloads.
- Retain the v0.1.5 offline HTML report as a downloadable companion view.
- Add unit, integration and Streamlit end-to-end coverage for the application.

- Move the generic package into its own repository while retaining its focused
  five-commit development history and exact v0.1.5 source tree.
- Add explicit migration provenance, MIT licence text, canonical repository
  metadata and a pinned GitHub Actions quality gate.
- Correct the installation path for the standalone repository and allow pip to
  verify declared dependencies rather than recommending `--no-deps`.
- Make no scientific calculation, schema or report-rendering changes.

## 0.1.5 - 2026-09-01

- Keep the PCoA view first but render points without edges by default, using
  equal x/y scaling, real ticks and axis titles containing per-group retained
  positive inertia.
- Add conservative `POOR`, `MODERATE` and `BETTER` two-dimensional fit guidance
  from retained inertia, input-versus-map distance correlation and normalised
  stress; these categories are explicitly display guidance rather than
  biological gates.
- Add a deterministic bounded Shepard plot so input and projected distances can
  be compared directly for every displayed group.
- Add a branch-length-scaled rectangular phylogram pruned to the displayed
  proteins while retaining the resolved gene tree's horizontal branch lengths,
  species colours and sampled-medoid role.
- Prefer normalised tree tables and permit a report-only compatibility fallback
  to an original tree only after its SHA-256 matches the immutable tree
  inventory; changed or missing sources remain explicitly unavailable.
- Add an exact bounded distance-matrix heatmap containing every supplied
  displayed pair distance and using complete phylogram leaf order when
  available.
- Retain the force-directed sparse nearest-neighbour topology as a separate,
  explicitly non-quantitative fourth detailed view.
- Replace the uninformative mean-distance-versus-fixed-sample-size plot with
  mean distance versus analytical group size while retaining displayed-sample
  provenance in point details.
- Calculate compact exact cluster-size, species-breadth and copy-number bins
  from the complete group-statistics authority; label size/breadth and distance
  panels separately as embedded-sample and selected-pilot views.
- Add a ranked projection/topology diagnostics table and raw, `log1p` and
  presence/absence shading for exact group-by-species copy counts.
- Preserve report-only regeneration from older immutable resources without
  repeating OrthoFinder parsing, distance calculations, Parquet or DuckDB work.

## 0.1.4 - 2026-09-01

- Clarify that every selectable protein network represents one HOG and expose
  the displayed, supplied and analytical group sizes separately.
- Make node fill colour species-specific by default, with deterministic colours
  that remain stable across expanded runs, exact collision avoidance and a
  searchable, clickable species legend.
- Identify the sampled medoid only from a complete displayed pairwise-distance
  matrix; retain its species fill while marking it with a gold star and border.
- Report the raw nearest-neighbour component and isolate counts, and add the
  minimum number of explicitly dashed, toggleable component connectors needed
  to keep a disconnected HOG visible as one layout.
- Use a two-dimensional classical-MDS/PCoA projection of the complete displayed
  distance matrix as the default node geometry, with positive inertia,
  pairwise-distance correlation and normalised stress reported per HOG.
- Present the PCoA map first and retain the force-directed neighbour topology as
  a separate, simultaneously available, explicitly non-quantitative view; both
  share HOG, species, search, medoid, label and connector controls.
- Declare NumPy as a direct numerical dependency.
- Log report-network progress and, per HOG, the displayed size, species count,
  raw components, connectors, projection status and projection stress.
- Hide dense member labels by default while retaining label, member and species
  search and selection controls.
- Add explicit x- and y-axis titles to every summary histogram and scatter plot,
  and state the row, column and cell semantics of the copy-count heatmap.
- State that force-directed positions are exploratory and that edge tooltips and
  distance summaries, rather than geometric spacing, are the quantitative
  authorities.

## 0.1.3 - 2026-09-01

- Fix the JavaScript function-name collision that stopped every run-wide chart
  after the static HTML headings appeared.
- Add a visible in-page rendering error so browser failures cannot remain silent.
- Reduce standalone report size by retaining only browser-used fields and by
  embedding fixed-width distance histogram bins instead of every pair value.
- Lower the default embedded group-summary bound to 20,000 and enforce a
  browser-safety ceiling of 50,000 rows; deterministic stratified sampling retains
  representation across group types and hierarchy levels.
- Ensure distance-backed networks use the exact distance-sampled member identifiers
  (or a deterministic subset when the report member bound is smaller).
- Add `--action report` to regenerate a separate HTML from a completed resource
  without repeating source parsing, distance calculations, Parquet or DuckDB work.
- Use scheduler-supplied node-local work storage for report input scans and remove
  temporary copies after success or failure.
- Add start, finish, elapsed-time and row/size logging for major stages, membership
  sources, every distance group, Parquet conversion and verified publication.
- Publish machine-readable `logs/stage_metrics.tsv` with every completed resource.

## 0.1.2 - 2026-08-31

- Declare every TSV column type during Parquet streaming so a late text value,
  such as HOG node `N0` after blank legacy rows, cannot conflict with null inference.
- Resolve exact, species-prefixed and OrthoFinder-internal tree leaves to canonical
  membership identifiers with explicit ambiguity and collision checks.
- Record the tree member-identifier resolution method in distance summaries.
- Add a queryable per-group/per-species copy-number relation for taxonomic and
  cluster-splitting analyses.
- Stream analytical TSV authorities directly to gzip and read compressed TSVs
  transparently during Parquet, DuckDB and HTML publication.
- Remove failed staging and incomplete copy directories by default, with an
  explicit `--keep-failed-work` diagnostic opt-in.
- Add self-contained cluster-size, species-breadth, copy-number, size/breadth,
  distance-coverage and group/species heatmap visualisations to the offline HTML.
- Bump the physical resource contract to schema version 2.

## 0.1.1 - 2026-08-31

- Stage Slurm analyses in job-specific node-local temporary storage by default.
- Add checksum-verified cross-filesystem copying followed by atomic persistent
  publication.
- Retain explicit ``--work-dir`` overrides for other schedulers and local runs.
- Default the bundled Dundee Slurm launcher to the `barton` partition.
- Add pycodestyle and pydocstyle to the reproducible Conda environment used by
  the complete release gate.

## 0.1.0 - 2026-08-31

- Add version-aware OrthoFinder 2 and 3 result discovery.
- Add long-form species, sequence, orthogroup and all-level HOG publication.
- Add tree inventories and optional Newick node/edge publication.
- Add explicit pairwise alignment and patristic-distance calculations.
- Add automatic alignment-to-resolved-tree distance fallback with per-cluster
  source provenance and explicit unavailable reasons.
- Add versioned Parquet, DuckDB, TSV, QC and provenance resources.
- Add an offline interactive HTML report with a Cytoscape-style network view.

# Changelog

## Unreleased

## 0.9.0 - 2026-09-10

- Add a checksum-bound Arabidopsis comparison authority containing 13
  housekeeping-reference candidates and 99 reviewed R/NLR candidates. These
  panels encode hypotheses to test and never predetermine compact or dispersed
  outcomes.
- Add `orthofinder-results --action dispersion-benchmark` for every E3,
  housekeeping-candidate and R/NLR-candidate HOG at `N0`, plus unique non-focus
  controls selected without examining distances and matched on protein count,
  represented-species count, mean copies per species and single-copy fraction.
- Measure cluster-level mean, median, population SD, interquartile range and
  coefficient of variation; retain protein pairs as within-cluster observations
  rather than incorrectly treating them as independent inferential replicates.
- Test planned E3, E3-category, housekeeping and R/NLR profile contrasts on
  matched-control residuals using tie-corrected Mann–Whitney tests, Cliff's delta,
  deterministic bootstrap confidence intervals and Benjamini–Hochberg FDR.
- Compare every individual biological target cluster with its own matched controls,
  pooled E3, housekeeping and R/NLR backgrounds, and every eligible subclass or
  E3 category using leave-one-out empirical tests and family-specific FDR.
- Publish complete compressed TSV, typed Parquet and DuckDB relations for marker
  matching, control matching, cluster statistics, background summaries, planned
  contrasts, individual tests and explicit compact/typical/dispersed classifications.
- Add a **Calibrated dispersion** application page with selectable metrics and
  scales, cluster-level distribution plots, confidence-interval forest plots,
  individual-cluster comparison views, interpretation help and paired TSV or
  formatted-Excel downloads.
- Add a portable Slurm wrapper that stages OrthoFinder inputs and authorities under
  scheduler `TMPDIR`, delegates checksum-verified publication to the pipeline and
  validates every required compressed result before reporting success.
- Bump the physical resource contract to schema 4 and retain read-only support for
  schema-2 and schema-3 resources in the application.

## 0.8.1 - 2026-09-10

- Attach matching checksum-verified portable tree payloads to persisted DuckDB
  distance analyses, restoring the branch-length phylogram for precomputed E3
  clusters without recalculating their distance matrices.
- Fix the E3 Slurm wrapper to pass the persistent formal output to the pipeline
  while retaining all resource construction beneath scheduler `TMPDIR`; this
  removes a conflict with the persistent-output safety policy.
- Delegate cross-filesystem checksum verification and atomic publication to the
  tested pipeline publication layer, while retaining the wrapper's destination
  lock, required-output validation and failure cleanup.

## 0.8.0 - 2026-09-09

- Replace the implicit default with the exact supplied 1,000-record
  `e3_seed_catalogue.tsv`, retain the previous broad authority only for
  reproducibility, and accept the rich seed catalogue through explicit
  primary-to-associated metadata fallbacks.
- Add `orthofinder-results --action e3-precursor` for every focus-matched HOG at
  `N0`, fixed to resolved-gene-tree patristic distances and prohibited from
  silently applying a largest-group limit.
- Force every matched E3 protein into a bounded deterministic member sample so
  the primary seeds remain represented in downstream distances and graphics.
- Publish `e3_cluster_results.tsv.gz` with one row per selected cluster,
  complete group/copy statistics, seed provenance, min/quantile/median/mean/SD/
  max distances, derived dispersion measures and explicit sampling coverage.
- Publish exact seed-to-member matches, a complete matched/unmatched seed audit,
  selected-cluster pairwise distances, matching typed Parquet relations and
  DuckDB tables; unavailable distances retain blank measures and a reason.
- Restrict portable and normalised gene-tree expansion to the focus-related
  tree identifiers for precursor runs while retaining a checksum inventory of
  every discovered tree.
- Add a scheduler wrapper that requires cluster `TMPDIR`, stages heavy inputs
  and first output on node-local storage, validates gzip outputs, checksum-
  compares the persistent copy and uses a guarded atomic publication rename.
- Fix Summary and result-table navigation with Streamlit callbacks so clicking
  **Find E3 focus clusters** no longer mutates an already-instantiated page
  widget; cover the real button click in a headless end-to-end test.
- Package the audited 60-label `Results_Feb26` taxonomy mapping as a conditional
  application default, selecting it only for an exact resource label-set match;
  retain 59 reviewed mappings and keep the apparent `Leismania_major` spelling
  error explicitly unmapped pending source confirmation.
- Explain directly beside the selection-tree uploader that it expects a
  species-to-taxonomy review table, not a protein list, tree or taxdump, and
  provide the current dataset's editable TSV and formatted workbook there.
- Extend unit, CLI, pipeline, Slurm and Streamlit end-to-end coverage to 421
  passing tests while retaining the enforced 95% branch-coverage quality gate.

## 0.7.0 - 2026-09-09

- Add a generic **Selection coverage tree** page and command-line action with
  exact-terminal, include-clade, only-in-clade, exact-exclusion and
  clade-exclusion predicates composed under explicit AND semantics.
- Reject direct and indirect selector contradictions before querying groups;
  several only-in clades use their terminal-set intersection and unresolved
  group members make only-in tests fail closed.
- Keep selection state independent from dataset-coverage state, distinguish
  represented, expected-no-data, ancestor, outsider and not-expected taxa, and
  state that expected-no-data is never evidence of biological absence.
- Build a deterministic minimal-ancestor taxonomy tree from a pinned reviewed
  mapping, retain stable authority or NCBI taxon identifiers, and optionally
  compact only unary neutral lineages.
- Add keyboard-equivalent node actions, synchronised last-valid controls and
  group results, accessible colour-plus-shape styling, and a safe transition
  from a passing group to the complete cluster explorer.
- Export complete reproducibility packages containing taxonomy, mapping,
  expected-universe, predicate, coverage, group and unmapped-label TSVs plus
  Newick, style data, SVG, PDF, JSON provenance and SHA-256 checksums.
- Add `orthofinder-results --action coverage-tree` with repeated named selector
  options, validate-only and dry-run modes, bounded group/node controls and
  verified atomic publication outside the immutable resource.
- Package a versioned 43,066-record E3 seed-protein authority as the default
  focus while keeping the query engine generic and allowing a plain or
  gzip-compressed replacement TSV in the app and CLI.
- Add a **Focus protein clusters** page that reports exact matching identifiers,
  proteins, species and alias authorities, then opens any result with its focus
  protein highlighted throughout the existing distance and visual suite.
- Parse controlled UniProt `sp|ACCESSION|ENTRY` and `tr|ACCESSION|ENTRY` aliases
  for exact accession or entry-name searches without heuristic identifier
  stripping.
- Preserve header-only group-audit TSVs for valid zero-match focus runs and
  strengthen control-character, taxonomy topology, checksum and filesystem
  failure validation.
- Keep the scheduler-TMP publication wrapper testable on macOS/BSD and Linux
  through portable UTC timestamps and an atomic publication lock around the
  same-filesystem directory rename.
- Extend unit, command-line and headless Streamlit end-to-end coverage to 401
  passing tests while retaining the enforced 95% branch-coverage quality gate.

## 0.6.0 - 2026-09-08

- Add a dedicated **Find a gene / protein** page that searches canonical
  OrthoFinder membership identifiers and, when available, internal identifiers
  from `SequenceIDs.txt`.
- Provide case-sensitive exact matching by default and an explicitly bounded,
  literal, case-insensitive contains mode whose wildcard characters are safely
  escaped.
- Return every matching HOG hierarchy record and legacy orthogroup separately,
  together with species, copy-count, group-size and preferred persisted-distance
  summaries.
- Open any matching membership directly in the complete cluster visualisation
  suite, with the requested protein highlighted across linked views.
- Add a protein-centred nearest-to-farthest exact-distance table, summary measures
  and paired TSV/formatted-Excel downloads.
- Guarantee that a requested protein is retained when a schema-3 portable-tree
  analysis must deterministically sample a large cluster, and include that
  requirement in the content-addressed cache identity.
- Report explicitly when an older schema-2 persisted sample omitted the requested
  protein and therefore cannot be expanded without a portable-tree rebuild.
- Extend model, SQL, sampling, cache, page-helper and Streamlit end-to-end tests
  for canonical IDs, aliases, literal matching, defensive bounds and focused
  distance analysis.

## 0.5.0 - 2026-09-08

- Add an **All distance results** page that exports one preferred successful
  persisted distance summary per cluster, with exact group-system, hierarchy,
  method and calculation-scope filters and a selectable 34-column schema.
- Prefer complete matrices over sampled matrices for the same cluster and,
  when only samples exist, select the largest stored sample deterministically.
- Add distance-centred result fields including sampling fraction, full-matrix
  status, five-number summaries, population SD, interquartile range and relative
  distance spread, while keeping source and calculation provenance selectable.
- Add formatted Excel workbooks beside every application TSV table download.
  Workbooks use frozen headings, filterable banded tables, readable bounded
  column widths, typed scientific number formats and a column-definitions sheet.
- Write identifiers and all other text explicitly as text to prevent spreadsheet
  formula interpretation, and enforce Excel row and column limits before export.
- Add a consistent expandable interpretation panel to every application graph,
  explaining what it shows, how to interpret it and its principal limitation.
- Keep on-demand sidecar calculations separate from the immutable dataset-wide
  export until they are published in a rebuilt resource, avoiding mixed
  authorities in reproducible results.
- Extend unit and Streamlit end-to-end coverage for result selection, query
  ranking, workbook internals, safe cell handling, graph-guidance completeness
  and the new application route.

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

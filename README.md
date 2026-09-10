# orthofinder-results

[![CI](https://github.com/peterthorpe5/orthofinder-results/actions/workflows/ci.yml/badge.svg)](https://github.com/peterthorpe5/orthofinder-results/actions/workflows/ci.yml)

`orthofinder-results` turns a completed OrthoFinder result directory into a
versioned, queryable and portable analytical resource. Its parsers, schema and
selection engine remain dataset-generic. For this project, the viewer ships a
replaceable E3 seed-protein authority as its default focus; a match is a
prioritisation flag, not an inferred functional annotation.

The package supports completed OrthoFinder 2 and OrthoFinder 3 layouts. For
OrthoFinder 3, hierarchical orthogroups are the primary authority (recorded as
`HOG` in the data contract) and a legacy `Orthogroups.tsv` file is optional.
Every hierarchy level (`N*.tsv`) is kept, so later analyses can work at the
root, at an internal species-tree node, or across several levels without
rebuilding the package.

This project is independent of, and is not endorsed by, the OrthoFinder
authors. Cite OrthoFinder itself when using its results.

The canonical repository is
[`peterthorpe5/orthofinder-results`](https://github.com/peterthorpe5/orthofinder-results).
The generic package was extracted from its original development location in
`E3_project_draft` without bringing E3-specific code or ranking assumptions
with it. See the
[migration provenance](https://github.com/peterthorpe5/orthofinder-results/blob/main/MIGRATION_PROVENANCE.md)
for the exact source commit, subtree hash and commit mapping.

## What one run publishes

Each successful run creates one immutable output directory containing:

- long-form legacy orthogroup and all-level hierarchical-orthogroup memberships;
- per-group and per-group/per-species copy-number statistics;
- species and sequence identifier mappings when OrthoFinder retained them;
- a checksum inventory of every discovered input, gene tree and resolved tree;
- one checksum-bound, compressed portable Newick payload per preferred gene tree,
  allowing later bounded distance calculations without access to the original run;
- a normalised species tree, and optionally every gene tree, as node/edge tables;
- optional aligned-sequence pairwise distances and per-cluster distributions;
- optional exact E3/focus selection, seed audit and one-row-per-cluster
  compressed result authority;
- matching gzip-compressed TSV and typed Parquet tables;
- a physical, portable DuckDB containing the same analytical relations;
- explicit QC checks, a complete run manifest and a persistent run log; and
- a self-contained offline HTML report with linked distance, phylogeny,
  sparse-topology and comparative visual summaries.

The HTML is the final publication stage. It embeds the JavaScript and CSS it
needs, so the report opens without internet access. Every network selector entry
is one OrthoFinder group at the stated authority and hierarchy level. Node fill
colour denotes species by default, a searchable legend can select all displayed
members from one species, and the sampled medoid is marked with a gold star. The
raw nearest-neighbour component count remains explicit;
dashed minimum-distance component connectors can be shown for a unified layout
or hidden to recover the raw graph. The first detailed view is a points-only
two-dimensional classical-MDS/PCoA projection of every displayed pairwise
distance. It has equal geometric axis scaling, explicit ticks, per-axis positive
inertia, a conservative fit-guidance badge and a bounded Shepard plot. Edges are
hidden by default. Correlation between projected and input distances, normalised
stress and retained inertia remain visible because no 2D map can preserve every
high-dimensional distance exactly.

Fit guidance is deliberately conservative. `POOR` is shown when axes 1+2 retain
less than 25% of positive inertia, distance correlation is below 0.60 or
undefined, or normalised stress exceeds 0.60. `MODERATE` is shown when the fit is
not poor but retained inertia is below 40%, correlation is below 0.80 or stress
exceeds 0.45. Other projections are labelled `BETTER`, not “validated”. These
thresholds guide display interpretation and do not define biological clusters.

A branch-length phylogram follows the PCoA. It prunes the checksum-verified
resolved gene tree to the displayed proteins while retaining horizontal branch
lengths. An exact displayed distance-matrix heatmap retains every supplied pair
distance and uses the phylogram leaf order when every displayed member resolves.
The separate force-directed nearest-neighbour topology remains explicitly
non-quantitative. Member, species, sampled-medoid and label interactions are
linked across the views.

Opening histograms use compact exact bins calculated from the complete group
authority, rather than the bounded embedded rows. The size-versus-breadth panel
is still identified as a deterministic sample, while distance and projection
panels are labelled as the selected network pilot. Mean distance is plotted
against analytical group size, not the fixed displayed sample size. The exact
group-by-species copy heatmap offers raw, `log1p` and presence/absence display
scales while retaining raw counts on hover. Every quantitative chart has
explicit axes, population and denominator wording.

Browser visualisation is deliberately bounded. Large groups are selected and
sampled deterministically for rendering, while the full compressed TSV, Parquet
and DuckDB tables remain the analytical authorities. Only fields used by the
browser are embedded. Exact pairwise distances are retained only for the bounded
displayed matrices, alongside compact histograms and sparse neighbour edges. The
default is 20,000 group summaries and the enforced ceiling is 50,000. Those
summaries are deterministically stratified across group type and hierarchy node,
rather than taking only the first rows from one level. The requested
network-group/member bounds must also remain within a
one-million-pair browser payload budget; the default 25 groups by 250 displayed
members is below that ceiling. The report states these limits rather than
pretending that a browser can safely render millions of nodes.

## Data contract designed for later questions

All group records use this composite identity:

```text
(run_id, group_type, hierarchy_node, group_id)
```

This prevents an identifier such as `OG0001686` or `N0.HOG0002084` from being
silently equated across independent OrthoFinder runs. Expanded runs with new
species must receive a new `run_id`. Future split, merge and overlap analyses
can then map member sets across runs without overwriting either source.

The schema preserves:

- the exact source species label and member identifier;
- one authoritative copy-count row per represented group/species pair;
- hierarchical-orthogroup node, parent clade and legacy orthogroup link when
  supplied;
- source filename and row number;
- tree type, tree identifier, branch lengths and node relationships;
- the exact distance method and `EXACT` or sampled status; and
- checksums and OrthoFinder/package versions.

These fields support later phylogeny, cluster-size, taxonomic distribution,
cluster splitting and expanded-species comparisons without redesigning the
raw ingestion layer.

## Installation

Conda is recommended on the cluster:

```bash
cd /path/to/orthofinder-results
conda env create --file environment.yml
conda activate orthofinder_results
python -m pip install --editable .
./run_tests.sh
```

For an existing Python 3.11+ environment:

```bash
python -m pip install --editable '.[app,dev]'
./run_tests.sh
```

The `app` extra installs the independent Streamlit viewer. It is not required
when the package is used only to build resources on a cluster.

## Complete E3 precursor resource

Version 0.8.0 adds a dedicated cluster action for the downstream E3 and motif
workflow. It matches the packaged or user-supplied focus authority exactly to
`HOG` groups at `N0`, calculates resolved-gene-tree patristic distances for
every matching cluster, and publishes an analysis-ready compressed table:

```bash
orthofinder-results \
  --action e3-precursor \
  --results-dir /path/to/completed/OrthoFinder/results \
  --output-dir /path/to/new/e3_resource \
  --run-id my_orthofinder_e3_run \
  --work-dir "${TMPDIR}/orthofinder_e3_work" \
  --distance-max-members 250
```

`tables/e3_cluster_results.tsv.gz` has one row for every matched cluster and
includes matched seed identities plus minimum, quartiles, median, mean,
population SD and maximum distances, calculation scope and provenance. The
companion `tables/pairwise_distances.tsv.gz` retains every calculated pair;
`e3_seed_catalogue_audit.tsv.gz` retains matched and unmatched seeds. Missing
distance calculations remain explicit `UNAVAILABLE` rows with blank statistics.

For Slurm, [`slurm/e3_precursor.sbatch`](slurm/e3_precursor.sbatch) requires a
scheduler-provided `TMPDIR` and stages source reads and resource construction
there. The pipeline checksum-verifies its hidden persistent copy before atomic
publication while the wrapper holds a destination lock. See
the [complete E3 precursor guide](docs/e3_precursor.md) for the exact contract,
outputs and downstream parsing guidance.

## Calibrated dispersion benchmarks

Version 0.9.0 adds a matched-background analysis that tests, rather than assumes,
whether E3, housekeeping-reference candidates and R/NLR candidates occupy clusters
with different evolutionary dispersion. The bundled Arabidopsis marker authority is
replaceable. Its classes and subclasses are biological comparison panels; they are
not expected-result labels.

```bash
orthofinder-results \
  --action dispersion-benchmark \
  --results-dir /path/to/completed/OrthoFinder/results \
  --output-dir /persistent/path/to/new/dispersion_resource \
  --run-id my_dispersion_benchmark \
  --work-dir "${TMPDIR}/orthofinder_dispersion_work" \
  --distance-max-members 250 \
  --benchmark-controls-per-group 3 \
  --benchmark-bootstrap-resamples 1000
```

Every biological target cluster receives unique non-focus controls matched without
using distance values. The matching variables are protein count, represented-species
count, mean copies per species and single-copy-species fraction. Statistical
replicates are clusters, never the many correlated protein pairs within a cluster.
The pipeline reports five complementary measures: mean and median distance,
population SD, interquartile range and coefficient of variation.

Profile contrasts use a matched residual: the target cluster statistic minus the
median statistic of its own controls. They report the median difference, a
deterministic bootstrap 95% confidence interval, Cliff's delta, a tie-corrected
two-sided Mann–Whitney p value and Benjamini–Hochberg FDR q value. Individual
clusters are also tested against their own controls, pooled E3, housekeeping and
R/NLR backgrounds, every E3 category and every eligible benchmark subclass.
Leave-one-out prevents a selected cluster from becoming part of its own reference.

The primary machine-readable outputs are:

- `tables/benchmark_marker_matches.tsv.gz`: the exact housekeeping and R/NLR
  genes/proteins, matched run proteins, species, clusters and source provenance;
- `tables/benchmark_group_profiles.tsv.gz`: the E3 or reference-marker genes that
  define each biological profile and cluster;
- `tables/benchmark_cluster_results.tsv.gz`: one row per biological target or
  matched-control cluster, with structure, distance summaries and sampling scope;
- `tables/benchmark_contrasts.tsv.gz`: planned profile-level tests;
- `tables/benchmark_individual_comparisons.tsv.gz`: each target cluster against
  every eligible background;
- `tables/benchmark_cluster_classifications.tsv.gz`: transparent empirical
  compact/typical/dispersed and low/typical/high-heterogeneity labels;
- `tables/benchmark_background_statistics.tsv.gz`: raw and matched-residual profile
  summaries; and
- `tables/benchmark_matched_controls.tsv.gz`: the auditable control assignment and
  matching score.

See the [dispersion benchmark guide](docs/dispersion_benchmark.md) for the full
scientific contract, authority provenance, output dictionary, cluster commands and
interpretation safeguards.

Version 0.9.1 is a schema-compatible viewer update for the same completed resource.
Its calibrated-dispersion page leads with plain-language broad comparisons, adds an
interactive matched-control cluster map, and provides a **Genes and clusters** tab.
That tab shows the exact E3, housekeeping or R/NLR markers used, their matched
proteins and clusters, every member of a selected cluster, and paired TSV or
formatted-Excel downloads. It does not require a new cluster run or DuckDB transfer.

## Interactive application

Version 0.8.0 provides the read-only standalone application. It opens either a
completed resource directory or its `duckdb/orthofinder_results.duckdb` file:

Version 0.8.1 also attaches each matching checksum-verified portable tree to
an existing persisted distance matrix, so precomputed E3 clusters retain their
branch-length phylogram without rebuilding the completed resource.

```bash
orthofinder-interrogation-app \
  --resource-dir /path/to/completed/resource
```

The launcher validates the manifest, schema, run identity and required DuckDB
relations before starting a local Streamlit server. Every database query is
read-only. If port 8501 is occupied, the launcher selects the first available
port from 8501 upwards; `--server-port` still requests one exact port. Optional
persistent cache and log paths remain outside the immutable resource:

```bash
orthofinder-interrogation-app \
  --resource-dir /path/to/completed/resource \
  --cache-dir "$HOME/Library/Caches/orthofinder-results" \
  --log-file "$HOME/orthofinder_results_app_logs/app.log"
```

The default cache is `$HOME/Library/Caches/orthofinder-results` on macOS and the
XDG user cache (or `$HOME/.cache/orthofinder-results`) on Linux. It never assumes
that a Mac has `/tmp`. Cache files are content-addressed by run, group, tree
checksum and analysis controls, gzip-compressed, atomically published with
user-only permissions and reproducible. The app refuses to put its cache inside
a completed resource.

The packaged focus authority lives in
[`src/orthofinder_interrogation_app/data`](src/orthofinder_interrogation_app/data).
The default is the exact 1,000-record `e3_seed_catalogue.tsv` supplied for the
current project, with seed identifiers, associated category/organism metadata,
review state and sequence provenance. Its checksum is documented beside the
file. The previous broader 43,066-record authority remains packaged for
reproducibility but is not selected implicitly. Copy and edit
`custom_focus_proteins.template.tsv`, then use the replacement without changing
code:

```bash
orthofinder-interrogation-app \
  --resource-dir /path/to/completed/resource \
  --focus-proteins /path/to/my_focus_proteins.tsv
```

Matching is exact against canonical member IDs, OrthoFinder internal IDs and
unambiguous accession/entry fields in canonical UniProt pipe identifiers. The
authority does not infer new E3 annotations. The default file and every custom
file are checksum-bound in downstream selection manifests.

The application provides:

- a guided Summary landing page organised around biological questions, with
  contextual help, readable table headings and an expandable scientific glossary;
- a **Focus protein clusters** page that starts from the packaged 1,000-record
  E3 seed catalogue, reports exact matched accessions and cluster membership,
  and accepts a replacement TSV without code changes;
- a dedicated gene/protein search across canonical membership identifiers and
  available OrthoFinder internal IDs, returning every matching HOG level and flat
  orthogroup before opening a protein-focused cluster view;
- bounded group/member searches with exact group type and hierarchy;
- exact stored-species filters with `ANY`, `ALL` and `EXACT_SET` semantics;
- rejection of every group containing any selected excluded species;
- member-count, species-count and persisted-distance bounds and sorting;
- per-species copy counts, member tables and filtered pair-distance TSV exports;
- a draggable inline force network, with optional layout-only component
  connectors, plus a separate static nearest-neighbour topology;
- original 2D, selectable-axis 2D and rotatable 3D PCoA, with separate 2D/3D
  retained-inertia, stress and distance-correlation diagnostics;
- a Shepard plot, branch-length phylogram and exact displayed distance matrix;
- a nearest-to-farthest distance table for a searched protein, with that protein
  highlighted throughout the linked cluster views and forcibly retained in any
  newly calculated bounded sample;
- histogram, violin, empirical-CDF, medoid-distance, member-centrality and
  species-pair heatmap views of within-group dispersion; and
- a 2–12-group workspace comparing means, population SDs, medians, exact
  displayed distributions and independent PCoA small multiples;
- an **All distance results** page for selecting, previewing and exporting the
  required biological, statistical and provenance columns across every cluster
  with a successful persisted distance calculation;
- a generic **Selection coverage tree** that combines exact, clade, only-in and
  exclusion predicates, evaluates E3-focus clusters by default, keeps selection
  state separate from dataset coverage, and exports a reconciled audit package;
- and a three-part expandable guide beside every graph describing what it shows,
  how to interpret it and its most important limitation.

Every downloadable application table is offered as both UTF-8 TSV and formatted
Excel. The `.xlsx` workbook freezes the top row, supplies filter dropdowns on a
banded Excel table, applies readable bounded column widths and scientific number
formats, and includes a **Column definitions** worksheet. Identifiers are written
as text so values beginning with spreadsheet formula characters remain inert.

Member and species selections are linked across one cluster's panels. PCoA and
force layouts retain their diagnostic or non-quantitative warnings; the exact
matrix, pair table and branch-length phylogram remain quantitative. Every view
states the exact method, full analytical membership, displayed exact or
deterministic sample, analysis source and cache state.

Schema-2 resources remain supported: their validated run-bound offline-report
matrices keep the original pilot groups working. Schema 3 adds a `tree_payloads`
relation containing one preferred resolved (or fallback original) gene tree for
portable on-demand analysis. For a selected group without persisted pairs, the
app resolves its parent OG tree, restricts it to exact HOG membership, calculates
at most 500 members and 124,750 exact pairs, then stores only the reproducible
sidecar result. It does not precompute billions of mostly unused pairs and never
modifies the completed DuckDB.

Exact species labels do not establish taxonomic ancestry. The taxonomic-search
contract is dataset-generic: it reads the exact species authority from the
selected resource and contains no fixed list of the 60 pilot species. A reviewed
row records the workflow/source/accepted names, NCBI taxon and parent IDs/names,
aligned lineage IDs/names, mapping method/source/version, reviewer, review time
and note. `PENDING_REVIEW`, `UNMAPPED`, `AMBIGUOUS` and missing labels remain
visible and never silently become descendants or reviewed outsiders.

The page provides an empty review template in both TSV and formatted Excel. The
reviewed mapping supplied back to the application remains versionable UTF-8 TSV.
A local extracted NCBI `taxdump` can also generate exact-name candidates for
every species in any resource:

```bash
orthofinder-taxonomy-map \
  --resource-dir /path/to/completed/resource \
  --taxdump-dir /path/to/extracted_ncbi_taxdump \
  --source-date 2026-09-07 \
  --source-version taxdump-20260907-sha256-REPLACE_ME \
  --output-tsv /persistent/path/taxonomy_candidates.tsv
```

Unique exact scientific-name or synonym matches are deliberately marked
`PENDING_REVIEW`, never `REVIEWED`; multiple matches remain `AMBIGUOUS` and
missing matches remain `UNMAPPED`. Confirm each taxon and lineage, then change
only accepted decisions to `REVIEWED` and complete `reviewed_by`,
`reviewed_at_utc` and `review_note` before descendant searches.

The mapping can be uploaded for one browser session or supplied as a sidecar:

```bash
orthofinder-interrogation-app \
  --resource-dir /path/to/completed/resource \
  --taxonomy-map /path/to/reviewed_taxonomy_mapping.tsv \
  --log-file /path/to/logs/orthofinder_interrogation_app.log
```

[`examples/results_feb26_ncbi_taxonomy_mapping_20260907.tsv`](examples/results_feb26_ncbi_taxonomy_mapping_20260907.tsv)
is the audited 60-label mapping for the `Results_Feb26` resource. An identical
packaged copy is the application's conditional default: it is selected only
when all 60 exact resource labels match, so it can never leak into another
dataset. It was resolved against the official NCBI `taxdump.tar.gz` snapshot
published on 7 September 2026. It contains 59 reviewed unique exact-name
matches and deliberately retains `Leismania_major` as `UNMAPPED`; confirm and
document that apparent workflow misspelling against the original FASTA
provenance before changing its status. Other datasets receive their own
label-derived review template and must supply a reviewed mapping.

Taxon-ID/descendant searches use only `REVIEWED` rows and provide four explicit
semantics:

- `contains` requires represented reviewed descendants;
- `enriched` applies a one-sided Fisher exact species-presence test against the
  reviewed sampled target/outside universe and Benjamini–Hochberg correction
  across every group at the selected exact authority and hierarchy;
- `exclusive within sampled analysis` rejects reviewed outsiders and unresolved
  represented labels; and
- `near-exclusive` applies declared outsider/unresolved limits and lists every
  retained outsider.

Exclusive results are claims only about species sampled in this OrthoFinder run,
not universal biological absence claims.

### Selection coverage tree

This workflow uses the reviewed mapping as a versioned offline taxonomy
authority. It does not draw the OrthoFinder species tree or a group gene tree.
Every selected predicate is joined by logical AND:

| Selector | Passing-group requirement |
|---|---|
| Required exact taxon | Every selected reviewed terminal occurs. |
| Include clade | Each selected clade contributes at least one reviewed terminal; outsiders are retained and reported. |
| Only in clade | Every mapped terminal lies within the selected clade; unresolved members fail closed. Multiple clades use their intersection. |
| Exclude exact taxon | No selected reviewed terminal occurs. |
| Exclude clade | No reviewed descendant of the selected clade occurs. |

The expected-taxon universe is independent of those selectors. By default it
contains every reviewed species supplied to the OrthoFinder run. A replacement
TSV can include reviewed terminal taxa expected by a sampling plan even when
they have no input data. **Expected no data** means only “not represented in
this imported dataset”; it is not evidence that a gene or biological function
is absent.

The application evaluates clusters containing configured focus proteins by
default, which centres the current project on E3-containing groups. Clear the
focus checkbox to evaluate the complete selected group authority. Downloadable
packages contain reconciled TSV audits, Newick plus a separate style table,
SVG, PDF, a JSON selection manifest and SHA-256 checksums.

The equivalent command-line workflow supports repeated named selectors:

```bash
orthofinder-results \
  --action coverage-tree \
  --resource-dir /path/to/completed/resource \
  --taxonomy-map /path/to/reviewed_taxonomy_mapping.tsv \
  --expected-taxa /path/to/expected_taxa.tsv \
  --focus-proteins /path/to/my_focus_proteins.tsv \
  --require-exact-tax-id 3702 \
  --include-clade-tax-id 33090 \
  --exclude-clade-tax-id 4751 \
  --coverage-group-type HOG \
  --coverage-hierarchy-node N0 \
  --coverage-max-groups 100000 \
  --output-dir /persistent/project/selection_coverage/e3_plants_v1
```

Omit `--focus-proteins` to use the packaged E3 authority. Add
`--coverage-all-groups` for a complete non-focus scope. `--validate-only`
checks authorities and predicates without group querying or output;
`--dry-run` performs the bounded query and export construction without writing.
The output must be a new path outside the completed resource. See
[`docs/selection_coverage_tree.md`](docs/selection_coverage_tree.md) for the
export contract and scheduler-TMP workflow.

## Inspect before running

Inspection is read-only and writes its complete output to a named persistent
file:

```bash
orthofinder-results \
  --action inspect \
  --results-dir /path/to/OrthoFinder/Results_Feb26 \
  --inspection-output "$HOME/orthofinder_results_feb26_inspection.json"
```

The inspection records the detected OrthoFinder version, chosen adapter,
primary group authority and available result capabilities.

## Run locally

All CLI controls are named. No source file beneath `--results-dir` is modified.

```bash
orthofinder-results \
  --action run \
  --results-dir /path/to/OrthoFinder/Results_Feb26 \
  --output-dir /persistent/project/orthofinder_results/results_feb26_v0_4_0 \
  --run-id results_feb26 \
  --work-dir /persistent/project/orthofinder_results/work \
  --report-max-groups 25 \
  --report-max-members 250 \
  --report-nearest-neighbours 3
```

If OrthoFinder produced `MultipleSequenceAlignments`, they are discovered. An
external aligned-FASTA directory can instead be declared:

```bash
orthofinder-results \
  --action run \
  --results-dir /path/to/OrthoFinder/Results_Feb26 \
  --output-dir /persistent/project/orthofinder_results/results_feb26_with_distances \
  --run-id results_feb26 \
  --work-dir /persistent/project/orthofinder_results/work \
  --alignment-dir /persistent/project/alignments \
  --distance-group-type HOG \
  --distance-hierarchy-node N0 \
  --distance-max-members 250
```

`--distance-source AUTO` is the default. It prefers aligned-sequence distances
when recognised alignments exist and otherwise uses resolved gene-tree branch
lengths. Tree-backed hierarchical-orthogroup calculations match the group's
legacy orthogroup link when available and restrict the tree to that group's
members. OrthoFinder 3 layouts without a legacy link can match a tree by the
hierarchical-orthogroup identifier itself.
Use `--distance-source NONE` to disable distances, or name
`ALIGNED_SEQUENCE`/`RESOLVED_GENE_TREE` to require one authority and fail when
it is unavailable.

Alignment filenames are treated as exact group identifiers. Protein distances
are amino-acid p-distances with pairwise deletion of gaps and ambiguous
residues. Unequal sequence lengths fail explicitly. Groups above the member
limit use a deterministic hash sample, recorded as
`DETERMINISTIC_MEMBER_SAMPLE`; smaller groups are recorded as `EXACT`.

For tree-backed calculations, clusters are ordered by decreasing member count
and `--distance-max-groups` applies to that deterministic order. A value of
zero requests every eligible group; choose an explicit bound for a first
cluster pilot. Missing trees or member-name mismatches produce a per-cluster
`UNAVAILABLE` summary with `failure_reason`, not a fabricated zero distance.
Pair and summary records carry the exact alignment or tree `source_file`.
Tree calculations preserve canonical membership identifiers in outputs while
explicitly resolving exact, species-prefixed and OrthoFinder-internal tree-leaf
aliases. The recorded `member_identifier_resolution` reports which mapping was
used. Missing, duplicate and ambiguous mappings fail that cluster explicitly;
the package never strips prefixes heuristically.

For an ordinary complete run, schema 3 or newer publishes checksum-verified preferred
gene trees as compressed portable payloads regardless of the precomputed
distance-group bound. Keeping `--distance-max-groups 25` therefore preserves a
fast opening pilot while the app can calculate other selected groups lazily.
The dedicated E3 precursor is different: it publishes the portable trees and
distance matrices for every focus-matched group, while omitting unrelated tree
payloads. Do not set the ordinary run bound to zero merely to support the app:
on a large run that would attempt every eligible pair matrix and can be
computationally and spatially prohibitive.

Use `--parse-gene-trees` only when normalised nodes and edges for every gene
tree are required. Tree files are checksum-inventoried even when their nodes
are not expanded, so a resumed run cannot miss a changed tree.

## Resume and replacement

`--resume` reuses an output only when the completed manifest, run identifier,
package version and full source digest match. `--force` never deletes an old
result: it moves it to a timestamped `.superseded.*` path before publication.
The two options are mutually exclusive.

The formal `--output-dir` must be persistent and paths beneath `/tmp` or
`/private/tmp` are rejected for that role. `--work-dir` may use node-local
temporary storage. Same-filesystem staging is published by atomic rename.
Cross-filesystem staging is copied into a hidden incoming directory beside the
formal output; every manifested file size and SHA-256 checksum is verified
before that directory is atomically renamed into place. A formal output is
therefore never exposed as complete while copying is still in progress.

Analytical table authorities are streamed directly to `.tsv.gz`; a multi-GB
uncompressed intermediate is not created. The gzip tables remain stream-readable
and are converted to typed, Zstandard-compressed Parquet before DuckDB is built.

On failure, the full traceback is written to the persistent Slurm error log and
partial staging/copy directories are removed by default. Use
`--keep-failed-work` only for a diagnostic rerun when the partial files themselves
are needed. Completed formal outputs and superseded outputs are never removed by
this cleanup policy.

## Regenerate only the standalone HTML

The report action can build a new compact report from an existing completed resource.
It does not mutate that resource or repeat OrthoFinder parsing, distance
calculation, Parquet conversion or DuckDB construction. The HTML and log must be
outside the immutable resource directory:

```bash
orthofinder-results \
  --action report \
  --resource-dir /persistent/project/orthofinder_results/results_feb26_v0_1_2 \
  --report-output /persistent/project/orthofinder_results/reports/results_feb26_report_v0_1_5.html \
  --log-output /persistent/project/orthofinder_results/reports/results_feb26_report_v0_1_5.log \
  --report-max-statistic-rows 20000 \
  --report-max-groups 25 \
  --report-max-members 250 \
  --report-nearest-neighbours 3
```

When the Slurm wrapper supplies its job-specific `--work-dir`, the compressed
relations needed by this action are copied to node-local storage, scanned there,
and removed after success or failure. This avoids high-volume repeated reads from
shared project storage while preserving the completed resource.

Normalised resolved-tree nodes and edges are preferred for phylograms. If an
older resource checksum-inventoried the tree but did not expand its nodes, the
report action may parse the exact distance-summary source only after its
SHA-256 matches the immutable inventory. Missing or changed sources produce an
explicit unavailable phylogram; they are never accepted silently.

## Slurm wrapper

The submission wrapper writes both standard output and error to persistent log
files, which is suitable for `mosh` sessions:

```bash
./submit_orthofinder_results_slurm.sh \
  --slurm-log-dir "$HOME/orthofinder_results_logs" \
  -- \
  --action run \
  --results-dir /path/to/OrthoFinder/Results_Feb26 \
  --output-dir /persistent/project/orthofinder_results/results_feb26_v0_4_0 \
  --run-id results_feb26
```

The wrapper prints the job identifier, exact output/error log paths and the
`squeue` command. Unless `--work-dir` is explicitly supplied, each Slurm job
uses a private directory below `${TMPDIR}`. The dedicated E3 and selection-
coverage and dispersion-benchmark wrappers fail closed if the scheduler did not
provide this variable;
they never substitute a Mac or cluster `/tmp` path. Only a completed, checksum-
verified result is copied to `--output-dir`. The Dundee launcher defaults to the `barton`
account and partition, requests ordinary resources, and does not select a
long-duration QoS. An explicit `--work-dir` remains available for clusters
with a different scratch policy.

The persistent `logs/run.log` now records stage boundaries and elapsed times,
million-row membership progress, every selected distance group, compressed TSV
and Parquet sizes, and report payload size. Scheduler stdout/stderr additionally
records cross-filesystem copy and checksum-validation timing. The resource's
`logs/stage_metrics.tsv` provides major computation-stage timings in a queryable
tab-separated form.

## Query examples

```sql
-- Largest root-level hierarchical orthogroups (HOG rows).
SELECT group_id, member_count, species_count, max_copies_per_species
FROM group_statistics
WHERE group_type = 'HOG' AND hierarchy_node = 'N0'
ORDER BY member_count DESC
LIMIT 25;

-- Distance coverage and central tendency for one group.
SELECT group_id, distance_method, computation_status,
       sampled_member_count, distance_pair_count,
       mean_distance, median_distance
FROM distance_statistics
WHERE run_id = 'results_feb26'
  AND hierarchy_node = 'N0'
  AND group_id = 'N0.HOG0002084';

-- Candidate splits: root hierarchical orthogroups represented at child levels.
SELECT legacy_orthogroup_id, hierarchy_node,
       count(DISTINCT group_id) AS child_hog_count
FROM hog_memberships
WHERE legacy_orthogroup_id <> ''
GROUP BY legacy_orthogroup_id, hierarchy_node
HAVING count(DISTINCT group_id) > 1
ORDER BY child_hog_count DESC;
```

Open the database with:

```bash
duckdb /path/to/output/duckdb/orthofinder_results.duckdb
```

## Scope of version 0.9

Version 0.9.0 adds calibrated cluster-dispersion benchmarking to the v0.8.1
complete E3 precursor. The exact 1,000-record E3 seed catalogue and Arabidopsis
comparison-marker authority are current project defaults, but users can replace
both. The underlying identifier, membership, distance, tree, matching and
statistics engines remain data driven. The standalone app owns OrthoFinder group
membership, copy number, species breadth, reviewed taxonomy, distances,
compactness, trees, within-run comparison and matched-background inference.
Version 0.9.1 improves the presentation and marker/membership interrogation of the
same schema-4 authority without changing its scientific calculations.
Explicit nested-HOG interrogation and cross-run cluster lineage (stable overlap
scores plus split/merge classification) remain later, separately tested generic
layers.

The dataset-wide page deliberately reports persisted resource results only. It
does not silently mix calculations from a user's mutable on-demand sidecar into
an immutable run-level export. Rebuild or publish those additional calculations
before treating them as dataset-wide authority.

Protein lookup searches identifiers stored by OrthoFinder plus controlled
accession and entry-name fields parsed from canonical UniProt pipe identifiers;
it does not infer gene symbols or descriptive aliases from external annotation
databases. Exact search is case-sensitive, while the optional contains mode is
literal, case-insensitive and browser-bounded. A selected protein is always
retained in a new schema-3 bounded tree calculation. If a schema-2 pilot matrix
omitted that protein, the app reports the limitation because the immutable
sample cannot be enlarged without a portable gene tree.

E3-ligase ranking, expression, experimental evidence, structures, conserved
ligandable pockets and chemistry starting points remain in the separate E3
application. A later versioned link contract can connect an exact composite
OrthoFinder group identity to that evidence without building E3-specific
assumptions into this reusable package.

## Development quality gate

```bash
./run_tests.sh
```

The gate runs unit and integration tests with branch coverage, pycodestyle,
pydocstyle, Ruff, Python compilation and shell syntax checks. Coverage must be
at least 95%.

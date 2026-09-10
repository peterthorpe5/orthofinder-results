# Calibrated cluster-dispersion benchmarking

## Purpose

The dispersion workflow asks whether selected biological classes occupy unusually
compact or dispersed OrthoFinder clusters. It does not encode an expected ordering.
Housekeeping-reference candidates and R/NLR candidates are empirical comparison
panels, and either may be compact, dispersed, mixed or unsupported in a particular
dataset.

The default biological profiles are:

- all E3 focus clusters and each E3 catalogue category;
- all housekeeping-reference candidates and their subclasses;
- all R/NLR candidates and their domain-architecture subclasses; and
- unique, distance-blind non-focus controls matched separately to every target.

An individual E3, housekeeping or R/NLR cluster is compared with its own matched
controls and with every eligible pooled profile or subclass. When the selected
cluster belongs to the reference profile, that cluster is removed before comparison.

## Authority provenance

The packaged `arabidopsis_dispersion_benchmarks.tsv` contains 112 records:

| Panel | Records | Provenance role |
|---|---:|---|
| Stable-expression references | 5 | Candidate reference genes from Czechowski et al. (2005), DOI `10.1104/pp.105.063743` |
| Core-cellular references | 8 | Reviewed Arabidopsis proteins chosen a priori as conventional cellular-reference candidates |
| R/NLR candidates | 99 | Col-0/Araport11 loci in Van de Weyer et al. (2019), supplementary table S3a, intersected with reviewed UniProt release 2026_03 proteins |

The authority checksum is
`33f30fcb31966dfde4136a8140bb9db00ada078b7e0ed01da6be73301b9a8b6f`.
Every row retains its exact protein identifier, Arabidopsis locus, profile class,
subclass, evidence label and source version. Marker matches and non-matches are
published so a partially represented panel cannot silently masquerade as complete.

## Statistical unit and estimands

One OrthoFinder cluster is one statistical replicate. Pairwise protein distances are
correlated observations used to summarise that cluster; they are never treated as
independent replicates for profile inference.

For every successfully calculated cluster, the workflow reports:

- mean pairwise distance: overall divergence;
- median pairwise distance: robust central divergence;
- population SD: absolute heterogeneity;
- interquartile range: robust heterogeneity; and
- coefficient of variation: SD divided by the mean, when the mean is positive.

Large clusters use the package's deterministic bounded member sample. Exact versus
sampled scope, sample fraction, requested marker retention, member resolution and
any failure remain explicit in every result.

## Distance-blind matching

Each E3 or benchmark-marker target is assigned a requested number of unique
non-focus control clusters. Controls are ranked using only protein count,
represented-species count, mean copies per represented species and the fraction of
represented species with one copy. Distance values and outcome classifications are
not examined during matching. Controls are not reused across targets.

The matched residual for target cluster \(i\) and metric \(m\) is

\[
r_{im} = y_{im} - \operatorname{median}(y_{c_1m}, \ldots, y_{c_km}).
\]

A positive residual means more dispersion than that cluster's matched controls; a
negative residual means more compactness. It is an adjustment within this dataset,
not a universal biological scale.

## Profile and individual inference

Planned profile contrasts compare cluster residuals. Each test reports target and
reference cluster counts, group medians, their difference, a deterministic
bootstrap 95% confidence interval, Cliff's delta, a tie-corrected two-sided
Mann–Whitney p value and Benjamini–Hochberg q value. FDR families are declared by
metric.

Individual-cluster comparisons use finite-sample empirical tail probabilities.
They cover the cluster's own raw matched controls and the matched-residual
distributions for pooled E3s, every E3 category, pooled housekeeping candidates,
each housekeeping subclass, pooled R/NLR candidates and each R/NLR subclass.
Within-profile tests are leave-one-out. Individual FDR families are defined by
background and metric across tested target clusters.

At least three eligible background clusters are required. Smaller panels are
reported as `INSUFFICIENT_GROUPS`; missing measurements are not converted to zero.
With only three own controls, empirical individual p values are necessarily coarse.
Use more controls when node-local time and space allow and individual-tail
resolution is important.

## Cluster execution

The wrapper requires seven positional arguments. A dash selects each packaged
authority:

```bash
PACKAGE_ROOT=/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/orthofinder-results
CONDA_ENV=orthofinder_results
RESULTS_DIR=/home/pthorpe001/data/2026_E3_protac/SSD_back_up_July_2026/Erin_Butterfield_data/Main_folder/OrthoFinder/Results_Feb26
RUN_ID=results_feb26_dispersion_benchmark_v0_9_0
OUTPUT_ROOT=/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/orthofinder_results_resources

sbatch \
  --job-name=of_dispersion \
  --partition=barton \
  --cpus-per-task=9 \
  --mem=96G \
  --time=2-00:00:00 \
  --output="${OUTPUT_ROOT}/slurm_logs/${RUN_ID}/orthofinder_results_%j.out" \
  --error="${OUTPUT_ROOT}/slurm_logs/${RUN_ID}/orthofinder_results_%j.err" \
  "${PACKAGE_ROOT}/slurm/dispersion_benchmark.sbatch" \
  "${PACKAGE_ROOT}" \
  "${CONDA_ENV}" \
  "${RESULTS_DIR}" \
  - \
  - \
  "${RUN_ID}" \
  "${OUTPUT_ROOT}/${RUN_ID}" \
  --distance-max-members 250 \
  --benchmark-controls-per-group 3 \
  --benchmark-bootstrap-resamples 1000
```

Create the persistent log directory before submission. The cluster scheduler must
provide an absolute writable `TMPDIR`. The wrapper stages the completed OrthoFinder
source and both authorities there. The pipeline uses node-local work space and
publishes the formal result through its checksum-verified atomic publication layer.
The wrapper refuses an existing formal output.

## Machine-readable outputs

| File | Unit | Purpose |
|---|---|---|
| `benchmark_marker_matches.tsv.gz` | marker-to-cluster match | Exact resolved marker evidence |
| `benchmark_marker_audit.tsv.gz` | authority marker | Complete matched/unmatched audit |
| `benchmark_group_profiles.tsv.gz` | cluster-profile membership | E3, benchmark-subclass and control labels |
| `benchmark_matched_controls.tsv.gz` | target-control pair | Matching variables, rank and score |
| `benchmark_cluster_results.tsv.gz` | selected cluster | Structure, sampling and five distance metrics |
| `benchmark_background_statistics.tsv.gz` | profile-scale-metric | Raw and matched-residual descriptive summaries |
| `benchmark_contrasts.tsv.gz` | planned profile contrast | Effect, confidence interval, p and FDR q |
| `benchmark_individual_comparisons.tsv.gz` | target-background-metric | Leave-one-out empirical placement and FDR q |
| `benchmark_cluster_classifications.tsv.gz` | target cluster | Empirical classification against own controls |

The same relations are typed in Parquet and DuckDB. Pairwise rows are retained for
biological target clusters so the application can render their full distance and
tree visual suite. Matched-control pair matrices are summarised but not persisted,
which prevents the database from expanding merely to support calibration.

## Application interpretation and gene downloads

Viewer version 0.9.1 opens an existing schema-4 benchmark resource directly; it
does not require resource reconstruction. The calibrated-dispersion page is divided
into three result-led tabs:

- **Results at a glance** reports the broad preplanned comparisons in plain language,
  then shows profile distributions and each target cluster's average-distance and
  distance-spread percentiles relative to its own matched controls.
- **Genes and clusters** lists the exact E3 seeds or Arabidopsis markers defining a
  selected profile. It retains the marker locus, protein accession and name, matched
  OrthoFinder member, species, cluster, evidence and source provenance. A second
  table contains every protein in a selected cluster and flags the matched proteins
  that define the selected profile.
- **Detailed statistics** retains selectable raw or matched-residual distributions,
  filtered confidence-interval forest plots, complete contrast downloads and each
  individual cluster's empirical placement against every eligible background.

Tables on the gene page have both TSV and formatted-Excel downloads. The Excel
workbooks freeze and filter the heading row, set readable column widths and include
the column definitions used by the interface. The marker coverage graph is bounded
to the 40 most widely represented markers for readability, but the table and both
downloads always retain every selected marker-to-cluster match.

## Interpretation safeguards

- A significant result is evidence for a distributional difference in this run,
  conditional on the selected markers, matching variables and distance sampling.
- A non-significant result is not proof of equivalence or absence of biological
  differences.
- A marker match identifies a cluster containing the listed protein. It does not
  transfer that protein's function to every cluster member.
- Effect direction, interval width, eligible cluster counts and sampling scope must
  be interpreted with the FDR q value.
- Apparent subclass results with fewer than three eligible clusters remain
  descriptive only.

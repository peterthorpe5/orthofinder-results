# Selection coverage tree

This workflow answers a bounded selection question against a reviewed taxonomy
and one exact OrthoFinder group authority. It does not alter the completed
resource. The interactive page and `--action coverage-tree` use the same domain
models, contradiction checks, group evaluation and export builder.

## Required authorities

- `--resource-dir` identifies one validated completed resource or standalone
  DuckDB.
- `--taxonomy-map` is a reviewed, versioned taxonomy mapping. Only `REVIEWED`
  rows are placed on the tree; unresolved labels remain in the audit.
- `--expected-taxa` is optional. Its columns are `taxon_id`, `reason`, `source`
  and `included`. When omitted, reviewed OrthoFinder input taxa define the
  expected universe.
- `--focus-proteins` is optional. When omitted, the versioned packaged E3 seed
  authority is used. `--coverage-all-groups` disables focus restriction.

The expected universe can name only reviewed terminal taxa from the selected
mapping. To represent a taxon expected by the sampling plan but absent from the
OrthoFinder input, add a reviewed mapping row with `role=expected`, then include
that stable taxon ID in the expected-taxa file.

## Validation and publication

Use `--validate-only` first to validate the resource, mapping, expected universe
and selector compatibility. Use `--dry-run` to perform bounded group evaluation
and construct every export in memory without writing. Normal publication writes
to a hidden sibling staging directory, verifies every file after writing and
atomically renames it to a new final output path.

The package contains:

- `taxonomy_nodes.tsv`, `input_taxon_mapping.tsv` and
  `unmapped_input_taxa.tsv`;
- `expected_taxa.tsv`, `selection_predicates.tsv`, `taxon_coverage.tsv`,
  `group_taxon_evaluation.tsv` and `tree_edges.tsv`;
- `selection_coverage_tree.newick` and
  `selection_coverage_tree_styles.tsv`;
- `selection_coverage_tree.svg` and `selection_coverage_tree.pdf`;
- `selection_manifest.json`; and
- `checksums.sha256`.

`input_taxon_mapping.tsv` is also a directly reloadable E3 taxonomy bridge. It
retains the required workflow label, accepted name, NCBI ID, aligned lineage
vectors, terminal rank, mapping status and role. Authority-specific cultivar
IDs remain in their own field and are never written into `ncbi_taxon_id`.
`group_taxon_evaluation.tsv` carries the optional group-level E3 bridge fields,
including mapped/unmapped counts, represented/outside taxon IDs, predicate
outcome, failure reasons and selection-manifest ID.

Newick carries topology and stable labelled taxon IDs. Selection and coverage
styles remain in the companion TSV so a plain phylogenetic viewer cannot silently
reinterpret display state as branch length or evolutionary inference.

## Cluster execution through scheduler-managed TMP

The supplied `slurm/selection_coverage_tree.sbatch` requires `${TMPDIR}` and
copies the complete read-only resource plus sidecar authorities into a unique
job directory there. It builds and verifies the selection package on node-local
storage, then uses a hidden incoming directory on persistent storage and a
same-parent atomic rename. The final path is therefore never exposed as complete
while copying.

Example submission for the Dundee Barton partition:

```bash
REPO_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/orthofinder-results"
RESOURCE_DIR="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/orthofinder_results_resources/results_feb26_of2_5_5_standalone_pilot25_v0_1_5"
TAXONOMY_MAP="${REPO_ROOT}/examples/results_feb26_ncbi_taxonomy_mapping_20260907.tsv"
OUTPUT_DIR="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/orthofinder_results_resources/selection_coverage_e3_v0_7_0"
LOG_DIR="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/orthofinder_results_resources/slurm_logs/selection_coverage_e3_v0_7_0"

mkdir -p "${LOG_DIR}"
sbatch \
  --account=barton \
  --partition=barton \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task=4 \
  --mem=32G \
  --time=08:00:00 \
  --job-name=of_selection \
  --output="${LOG_DIR}/selection_coverage_%j.out" \
  --error="${LOG_DIR}/selection_coverage_%j.err" \
  "${REPO_ROOT}/slurm/selection_coverage_tree.sbatch" \
  "${REPO_ROOT}" \
  orthofinder_results \
  "${RESOURCE_DIR}" \
  "${TAXONOMY_MAP}" \
  - \
  - \
  "${OUTPUT_DIR}" \
  --include-clade-tax-id 33090 \
  --coverage-group-type HOG \
  --coverage-hierarchy-node N0 \
  --coverage-max-groups 100000
```

The two `-` values mean “use the default expected universe” and “use the
packaged E3 focus authority”. Replace them with persistent expected-taxa and
custom-focus paths when required. Do not pass `--action`, `--resource-dir`,
`--taxonomy-map`, `--expected-taxa`, `--focus-proteins` or `--output-dir` again
in the trailing selector options.

# Complete E3 precursor analysis

The `e3-precursor` action builds a complete schema-3 resource for every `HOG`
at hierarchy node `N0` containing at least one exact identifier from a validated
focus authority. It is the machine-readable precursor for the E3 motif workflow;
the Streamlit app consumes the same DuckDB and tables for visual review.

## Fixed scientific contract

The dedicated action deliberately fixes the following choices:

- group collection: `HOG` at `N0`;
- distance authority: resolved-gene-tree patristic branch length;
- selected groups: every focus-matched group, with no largest-group truncation;
- tree portability: preferred resolved Newick payloads for selected groups;
- matching: exact canonical protein ID, OrthoFinder internal ID, or an exact
  accession/entry field from a syntactically valid `sp|ACCESSION|ENTRY` or
  `tr|ACCESSION|ENTRY` identifier; and
- sampling: at most `--distance-max-members` members per group, while forcing
  every matched focus protein into a deterministic sample.

The default focus authority is the packaged `data/e3_seed_catalogue.tsv` with
1,000 records and SHA-256
`10945b7cae2212f3bc425bd80a37481c96f3194742724bdabe41c0148e73db52`.
Pass `--focus-proteins` to substitute a reviewed custom TSV without changing
code. The exact input is copied into resource provenance and its checksum is
repeated in each E3 analytical relation.

`--distance-max-groups` must remain zero. This means every matched cluster is
attempted; it does not mean an unbounded all-pairs matrix inside a very large
cluster. A group larger than the member bound receives an explicitly labelled
`DETERMINISTIC_MEMBER_SAMPLE`. Increase the bound only after estimating the
quadratic pair count, `n × (n - 1) / 2`.

## Machine-readable outputs

The primary output is:

```text
tables/e3_cluster_results.tsv.gz
```

It contains exactly one row per selected E3 cluster. Fields include complete
group size/copy statistics, matched seed/member/species identities, distance
method and status, analysed member/pair counts, minimum, 5th percentile, first
quartile, median, mean, third quartile, 95th percentile, maximum, population
standard deviation, interquartile range, range, coefficient of variation,
sampling fraction, resolved-pair fraction and provenance.

Distance values remain empty when no pairs were calculated. A singleton,
missing tree or identifier-resolution failure remains in the file with
`computation_status=UNAVAILABLE` and `failure_reason`; it is never converted to
a biological zero.

Companion authorities are:

| Path | Grain and purpose |
|---|---|
| `tables/pairwise_distances.tsv.gz` | One exact calculated member pair for selected E3 clusters only |
| `tables/e3_seed_matches.tsv.gz` | One seed/member/species/group/match-authority relationship |
| `tables/e3_seed_catalogue_audit.tsv.gz` | One row per configured seed, including unmatched and disabled seeds |
| `tables/e3_cluster_results.parquet` | Typed columnar copy of the primary cluster table |
| `duckdb/orthofinder_results.duckdb` | Complete generic and E3-specific query backend |
| `provenance/focus_protein_authority.tsv` | Byte-for-byte focus input for a plain TSV run |
| `qc/validation_checks.tsv` | Generic and E3 reconciliation checks |
| `run_manifest.json` | Checksums, counts, versions, inputs and output inventory |

## Direct command

Run in new scheduler-managed storage and publish to a new output directory:

```bash
orthofinder-results \
  --action e3-precursor \
  --results-dir /path/to/completed/OrthoFinder/results \
  --output-dir /path/to/new/e3_resource \
  --run-id results_feb26_e3_precursor_v0_8_0 \
  --work-dir "${TMPDIR}/orthofinder_e3_work" \
  --distance-max-members 250
```

The Slurm wrapper is preferred because it stages all source reads and initial
outputs below scheduler-provided `TMPDIR`, validates the compressed outputs,
copies to a hidden incoming directory on persistent storage, checksum-compares
the complete copy and performs a locked same-filesystem rename:

```text
slurm/e3_precursor.sbatch
```

Pass `-` as its focus-authority argument to use the packaged catalogue.

## Lightweight downstream parsing

Inspect the primary table without copying the DuckDB:

```bash
gzip -cd tables/e3_cluster_results.tsv.gz | sed -n '1,6p'
```

Keep the complete composite identity `(run_id, group_type, hierarchy_node,
group_id)` in downstream joins. `group_id` alone must not be equated between
independent OrthoFinder runs.

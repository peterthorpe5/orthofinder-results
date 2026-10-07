# Protein-motif conservation and RNA-seq evidence

## Scientific question

This workflow identifies OrthoFinder groups containing proteins that match a user-defined
sequence pattern and asks how broadly that pattern and transcript evidence are distributed
across the sampled phylogeny. The default pilot question is whether proteins ending in
asparagine (`N`) occur in broadly conserved plant groups that could justify experimental
testing as candidate Cereblon substrates. The implementation is deliberately generic:

- the exact C-terminal suffix can contain any canonical amino-acid sequence;
- regular expressions can be enabled explicitly and searched anywhere in a protein or
  anchored at the C-terminus;
- either root/child hierarchical orthogroups or original flat orthogroups can be queried;
- primary, required-comparison and excluded lineages are selected from reviewed taxonomy;
- RNA-seq evidence is attached without labelling any protein as a Cereblon substrate.

A matching motif, orthology and RNA-seq evidence can prioritise pilot experiments. None of
them demonstrates Cereblon binding, ubiquitination, protein accumulation or degradation.

## Sequence authority and search denominator

The DuckDB `protein_sequences` relation is constructed from the same completed OrthoFinder
run as the group tables. `WorkingDirectory/SequenceIDs.txt` supplies the exact internal,
species and published member identifiers; `WorkingDirectory/Species*.fa` supplies the
sequences. Publication fails unless every identifier reconciles to exactly one sequence.
Each protein is stored once, rather than duplicated across every HOG hierarchy level.
The viewer remains compatible with the external Parquet sequence sidecars produced for
older schema-4 resources.

Resource construction streams the FASTA authority directly into typed Parquet batches of
at most 25,000 proteins. It does not retain the complete 1.4-million-protein sequence set,
a duplicate list of Python records and one monolithic Arrow table at the same time. RNA-seq
integration is deliberately completed before this sequence publication stage so native
Arrow allocator memory cannot raise the expression analysis's starting resident set.

The app distinguishes four quantities:

1. all published members of the selected group;
2. members with an assessable sequence;
3. assessed members belonging to the selected primary species set; and
4. primary-set members matching the active search.

The protein conservation fraction is item 4 divided by item 3. The species conservation
fraction is primary species with at least one match divided by primary species with at
least one assessed sequence. Missing sequences remain unavailable and never count as
non-matching sequence evidence. Sequence coverage is displayed separately.

Exact mode is the default. Regex mode is off until the user enables it. Regex patterns use
DuckDB's RE2-compatible engine over upper-case amino-acid sequences. C-terminal regex mode
adds an end anchor unless the expression already ends in `$` or `\Z`.

## RNA-seq input authority

The production rebuild uses the corrected Expression Atlas authority:

```text
/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/
analysis/expression_atlas_rebuild_v0_5_1_20260804/
manifests/e3_workflow_expression_resources.tsv
```

The manifest inventories included expression and metadata Parquet partitions, exact species
labels and SHA-256 digests. The production wrapper filters it to the 12 focal plants in the
project species authority plus *Homo sapiens*, and every selected file is verified before
resource construction. The other 47 OrthoFinder species remain fully available for orthology,
taxonomy, phylogeny and motif analysis but are explicitly `NOT_ASSESSED` for RNA-seq. They are
not counted as zero expression, absent genes or mapping failures. The filtered authority and
selected partition digests remain in the input inventory.

## Identifier mapping

Aliases are derived generically from the original `SequenceIDs.txt` member and FASTA header:

- exact member identifier;
- first raw FASTA-header token;
- UniProt accession and entry name for canonical `sp|...|...` or `tr|...|...` tokens;
- `GN=`, `gene=`, `gene_id=`, `locus=` and `locus_tag=` values; and
- a versionless form when an identifier ends in a dot and digits.

An optional reviewed alias TSV can add exact species/member-specific aliases. Matching is
case-insensitive but exact and is always restricted to the same species. The lowest mapping
tier wins. One gene at that tier is `MAPPED_UNIQUE`; several genes are `AMBIGUOUS`; none is
`NOT_MAPPED`. Ambiguous and not-mapped proteins remain unavailable rather than being guessed.

## Unit selection and missingness

For each species and Expression Atlas experiment, TPM is selected when any TPM records are
available. FPKM is used only when TPM is absent for that exact species/experiment. TPM and
FPKM are never pooled into one numerical distribution or heatmap. The selected unit is
retained in every context and summary.

A context is positive at the configured threshold (default `0.5`) and a member is called
broadly expressed when its positive-context fraction reaches the configured value (default
`0.5`). These are transparent screening thresholds, not significance tests.

The following states are distinct:

- a valid expression value of zero;
- a uniquely mapped gene with no included context records;
- an ambiguous identifier mapping;
- an unmapped identifier; and
- absent or unavailable metadata.

Only the first is measured zero. The others remain missing or unavailable.

## Heatmap

The application exposes RNA-seq in two equivalent ways: choose **RNA-seq explorer** directly
in the sidebar, or select the **RNA-seq explorer** workspace at the top of **Protein motif
conservation**. Both routes use the same sequence, orthology and taxonomy filters. The direct
page is intended for users whose question begins with expression; the workspace switcher is
useful when moving between candidate, taxonomic, protein and expression evidence for the same
search. Only the selected workspace runs its larger queries and figures.

The heatmap rows are selected motif-qualified groups. Columns are exact species plus the
selected biological context: Expression Atlas sample label, organism part/tissue,
developmental stage or condition. Each populated cell is the median expression value across
mapped group members and matching source rows for that species/context/unit. Hover text gives
the untransformed median, mapped-protein count, experiment count and positive-context fraction.

`log2(1 + value)` is available to reduce visual domination by very high values. It is a
display transformation only; the exported cell table retains the original unit and values.
Blank cells mean that compatible evidence was unavailable, not that expression was zero.

Interpret broad coloured bands as transcript evidence spanning several species or contexts.
Interpret isolated cells cautiously: they can reflect lineage specificity, incomplete
mapping, differences in sampled tissues or experiment availability. RNA-seq cannot determine
whether a protein accumulates in a Cereblon mutant.

## UpSet plot

The UpSet view summarises exact intersections of species with observed expression evidence.
A group belongs to a species set when at least one group member in that species has one or
more compatible expression contexts. The upper bars count groups in each exact combination;
the lower matrix shows which species define that combination. The display is bounded to eight
selected species for readability, while the downloadable intersection table lists every group
in each displayed combination.

The plot answers a coverage question: “for which species do we have observed transcript
evidence for these groups?” It is not an expression-level comparison and does not make an
absence claim for species lacking evidence.

## Selected RNA-seq evidence downloads

The **Expression species** control governs the heatmap query, the species available to the
UpSet plot, the group-by-species coverage table and the detailed exports. Four complementary
tables are available as TSV and formatted Excel:

1. **Aggregated heatmap cells** contain the exact untransformed values behind the displayed
   heatmap, including group, species, selected context label, expression unit, median, range,
   contributing member/experiment counts and positive-context fraction.
2. **Group-by-species coverage** reports members, mapped members, members with observed
   expression contexts and broad-expression calls for each selected group/species pair.
3. **Protein mapping and missingness** retains every selected group member, including exact
   mapping status, matched identifiers, selected units, observed-context counts and evidence
   status. This is the authoritative table for distinguishing an observed zero from an
   ambiguous, unmapped or otherwise unavailable protein.
4. **Underlying context records** contain the observed Expression Atlas rows used by the
   aggregation: experiment, sample/condition, tissue, stage, genotype, treatment, value,
   unit, positivity call and checksum-bound source provenance. Loading this table is explicit
   because it can be much larger than the summaries.

TPM and FPKM are queried and exported separately. The detailed context query has a selectable
protective bound of 10,000–100,000 rows, and the member export is bounded at 100,000 rows.
The app previews at most 2,000 detailed records while both downloads contain all retained
rows. If a bound is reached, the app labels the export as truncated and instructs the user to
narrow the group/species selection or increase the context bound; it never implies that a
partial table is complete. Proteins without observed contexts occur in the mapping/missingness
table but cannot occur in the observed context table.

## Published schema-5 outputs

| Output | Unit | Purpose |
|---|---|---|
| DuckDB `protein_sequences` | protein | Exact complete-proteome sequence authority |
| `expression_identifier_aliases.tsv.gz` | protein alias | Auditable identifiers used for mapping |
| `expression_member_mapping.tsv.gz` | protein | Unique, ambiguous or not-mapped state |
| `expression_member_summary.tsv.gz` | protein | Experiments, contexts, range and broad-expression call |
| DuckDB `expression_context` | protein/context | Unit-safe expression and biological metadata |
| `expression_group_summary.tsv.gz` | group | Mapping and observed-expression coverage |
| `expression_import_audit.tsv.gz` | resource build | Source counts, thresholds and policies |

Small mapping, summary and audit tables are retained as compressed TSV and embedded as typed
DuckDB relations. Expression integration runs in a fresh process against a file-backed work
database. Protein aliases first join a deduplicated species/gene catalogue, so one alias is
not multiplied by every tissue and condition before its identity is resolved. The largest
member-context relation is streamed directly to typed Parquet rather than materialised in
RAM. High-volume expression-context and protein-sequence authorities are then materialised
once in the final DuckDB and their construction Parquet files are removed after checkpoint.
This prevents duplicate multi-gigabyte completed outputs without changing any app query.
Visible app tables have exact TSV and formatted-Excel downloads;
every quantitative Plotly figure has a deferred PDF download.

Final DuckDB publication uses explicit thread and memory controls. If a relation cannot be
materialised within that buffer-manager limit, DuckDB may spill temporary construction data
to a dedicated directory beneath scheduler-local work storage. The expression worker logs
each substage's elapsed time, process RSS, work-database size and spill size and returns a
small validated metadata record to the parent. A native crash therefore cannot corrupt the
formal resource and can be localised to one named substage. Construction directories are
removed only after a successful checkpoint. These controls bound database construction
without changing stored values or identifiers. Expression group denominators include only
the reviewed species that were actually assessed.

## Performance and integrity boundary

Full SHA-256 verification is a construction-time operation. The viewer validates the local
manifest, schema, run identity and required DuckDB relations once per process, opens DuckDB
read-only, and does not repeatedly checksum the external Atlas partitions. Result views are
lazy: selecting the RNA-seq view triggers its bounded queries, while other motif views do not.
Low-cardinality selectors are cached. TSV, Excel and PDF files are generated only after their
download control is activated.

If the four required expression relations are absent, the dedicated explorer reports an
allowed-missing resource rather than showing an empty result. Motif, taxonomy, phylogeny and
dispersion pages remain usable, and unavailable expression evidence is never converted to
zero. Versions 0.11.5 and 0.11.6 are read-compatible with a completed version-0.11.4
schema-5 resource; these viewer updates do not require rebuilding or copying that resource
again.

## Safe resource replacement

Build and validate a new run ID beside the current resource. Do not delete the old directory
before the new manifest is complete, every compressed table passes `gzip -t`, DuckDB opens
read-only, the embedded `protein_sequences` and `expression_context` relations are non-empty,
and the app renders the motif, heatmap
and UpSet views. After Mac transfer and checksum validation, the old resource can first be
moved to a dated retirement directory. Permanent deletion should be a later explicit step.

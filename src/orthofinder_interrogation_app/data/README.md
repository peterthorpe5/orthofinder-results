# Focus-protein authorities

`e3_seed_catalogue.tsv` is the default, replaceable focus profile. It is the
exact 1,000-record reviewed project handover supplied on 2026-09-09, including
its seed/sequence identifiers, associated names, categories, organisms,
review states and source metadata. Its SHA-256 checksum is
`10945b7cae2212f3bc425bd80a37481c96f3194742724bdabe41c0148e73db52`.

`default_e3_focus_proteins.tsv.gz` is retained as the previous broader
43,066-record authority for reproducibility, but it is no longer selected by
default.

The generic OrthoFinder engine does not classify proteins as E3 ligases. It only
maps exact identifiers from the selected focus authority to the immutable
OrthoFinder membership tables. A match is evidence that a cluster contains a
listed seed, not proof that every cluster member has E3 function.

To use another project list, copy `custom_focus_proteins.template.tsv`, edit it
as UTF-8 TSV, and select it in the application or pass it with the named
`--focus-proteins` option. The required heading may be `seed_id`,
`protein_identifier` or `accession`. The optional `enabled` field accepts
true/false, yes/no or 1/0. Rich catalogue rows use explicit seed fields first
and associated metadata only as a declared fallback; this is not heuristic
annotation transfer.

## Conditional Results_Feb26 taxonomy mapping

`results_feb26_ncbi_taxonomy_mapping_20260907.tsv` is the packaged 60-label
species-to-taxonomy authority for the current study. It is byte-identical to
the provenance example under `examples/` and has SHA-256 checksum
`7d576f19355f4c44842a812a387060279bc830025b482e45d7f8942e98de7bdc`.

The application selects this mapping automatically only when the opened
resource has exactly the same 60 workflow species labels. It contains 59
`REVIEWED` decisions. `Leismania_major` remains `UNMAPPED` because its apparent
misspelling of *Leishmania major* must be confirmed against the original FASTA
provenance before it can support descendant or exclusivity claims. Any other
dataset receives a fresh label-derived review template; this file is never a
global species list.

## Arabidopsis dispersion comparison authority

`arabidopsis_dispersion_benchmarks.tsv` is the default, replaceable marker
authority for calibrated dispersion analyses. It contains 112 exact reviewed
Arabidopsis protein identifiers: 13 housekeeping-reference candidates and 99
R/NLR candidates. Its SHA-256 checksum is
`33f30fcb31966dfde4136a8140bb9db00ada078b7e0ed01da6be73301b9a8b6f`.

The housekeeping panel contains five stable-expression references drawn from
Czechowski et al. (2005) plus eight conventional core-cellular reference
candidates. The R/NLR panel is the exact intersection of the Col-0/Araport11 loci
listed in Van de Weyer et al. (2019), supplementary table S3a, with reviewed
Arabidopsis proteins in UniProt release 2026_03. Source title, DOI, table, release,
class and subclass remain attached to every marker row.

These labels define comparison panels and hypotheses; they do not assert that a
matched cluster is compact or dispersed. The pipeline tests that proposition from
the selected OrthoFinder dataset and preserves unmatched markers in a complete
audit. To use another species or panel, copy the TSV, retain its exact headings,
replace the records with reviewed protein identifiers and pass it through
`--benchmark-proteins`. The package never transfers Arabidopsis classifications to
another dataset by identifier spelling or taxonomic inference.

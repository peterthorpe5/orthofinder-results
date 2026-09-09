# Focus-protein authorities

`default_e3_focus_proteins.tsv.gz` is the default, replaceable focus profile.
It contains 43,066 known E3 seed accessions inherited from the E3 discovery
authority, together with their evidence category, source organism and provenance.
Its SHA-256 checksum is
`9e5ba99e751651be37e9abbaee445ef931e414c0c0ea3f70d7553eeecebcf8f7`.

The generic OrthoFinder engine does not classify proteins as E3 ligases. It only
maps exact identifiers from the selected focus authority to the immutable
OrthoFinder membership tables. A match is evidence that a cluster contains a
listed seed, not proof that every cluster member has E3 function.

To use another project list, copy `custom_focus_proteins.template.tsv`, edit it
as UTF-8 TSV, and select it in the application or pass it with the named
`--focus-proteins` option. The required heading may be `protein_identifier` or
`accession`. The optional `enabled` field accepts true/false, yes/no or 1/0.

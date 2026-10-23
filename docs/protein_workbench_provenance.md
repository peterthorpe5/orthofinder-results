# Protein workbench source provenance

## OrthoFinder baseline

- Repository: https://github.com/peterthorpe5/orthofinder-results
- Branch inspected: `main`.
- Exact base commit: `df71821c4abf37e8c0c99c868bc36ed3379b90da`.
- Base package version: 0.11.6; updated package version: 0.12.0.
- Published resource schema remains 5; resources supported by the existing viewer
  remain supported. No resource migration is required for the workbench.

## Supplied E3 sources

The user supplied two archives. They were inspected directly from their local
attachment copies. The v0.19.1 source overlay supplied the final generic engine,
viewer bundles, frontend source and exchange correction. The older archive was
used to understand the preceding architecture and documentation.

| Supplied archive | SHA-256 |
| --- | --- |
| e3_protein_research_v0_19_0.zip | b5de9a8c6029be7ca7b1f7d0000916ca76a17ddd71d0cf9a8622a45f518e0d97 |
| e3_deepclust_exchange_v0_19_1.zip | 422e5cc457118ca48e4348ba1372e5ebf037084581fe7ee9c2641fd0c19085ea |

The E3 resource adapter and E3 application were not copied into this independent
package. Imports, package-data lookup, command names, application version records
and the research configuration environment variable were adapted to OrthoFinder.
The default workspaces and commands remain distinct so both applications can be
installed in the same environment.

## Added OrthoFinder integration

- Exact run/group-system/hierarchy/group selection and current-group navigation.
- Complete-membership and sequence-authority checks before local FASTA transfer.
- Reusable local family ingestion with exact source context.
- Standard host-app TSV, Excel and deferred figure export controls.
- Source metadata inheritance and topology-preserving Auspice streamtree labels.
- Portable upstream extraction wrapper and real scientific smoke script.
- Mac installation, selected-export and operation guides, plus integration tests.

## Upstream extraction compatibility

The separately cloned upstream repository was read and executed without source
modifications at commit `b406fec8460d8044697c057bcb8c34f63549d98e`:
https://github.com/drostlab/deepclust_dataretrieval

A synthetic official-schema cluster spanning two physical row groups was
extracted through this project's shell wrapper. The test checks unchanged exact
IDs, all three requested members and exclusion of an unrelated cluster. This
establishes interface compatibility for the small fixture, rather than performance
or availability of the full global database. The upstream code is not bundled.

## Scientific and provider validation scope

A real local workflow used DIAMOND 2.2.8, MAFFT 7.505 and FastTree 2.1.11, followed
by exact membership auditing, alignment statistics and Auspice JSON generation.
These development versions are not claims that other tool versions produce
identical clusters or inferred trees. Tool commands and versions are recorded
for each job.

Model and annotation adapters use deterministic fixtures. No live model endpoint,
GPU checkpoint, remote sequence submission or full external annotation database
was run during this update. The three local viewer bundles retain the supplied
asset checksums and notices. The final delivery's validation record reports the
available automated checks and any browser-validation limitation.

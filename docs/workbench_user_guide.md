# Protein research workbench user guide

Aim: collect a relevant protein family, compare its sequences and interpret conservation, relationships and optional structure/annotation evidence.

Start with one family. Find/collect → analyse → inspect results. Models and Annotations can also start directly from an unaligned protein FASTA.

## Starting from OrthoFinder

Open OrthoFinder groups, choose one HOG or orthogroup at the correct hierarchy level, and Prepare complete group inputs. Inspect its full membership and species. Import selected local protein sequences if available, or queue global ID mapping when the DeepClust indices are installed. After import completes, use the next-step button in Jobs to choose family.fasta for analysis. Local group selection does not select a medoid or change the original OrthoFinder resource.

## Working on de.NBI

Launch with launch_orthofinder_denbi.sh on your Mac. The OrthoFinder app and its worker run beside the installed DeepClust database on the server. Reuse E3's data/tool settings in a separate OrthoFinder profile and job workspace. Ranking, ID mapping and complete selected-cluster retrieval use server-local files; do not download the database to your Mac. Closing the tunnel leaves server jobs running. Reconnect and open Jobs to inspect results. Optional remote model/annotation providers still require their explicit submission choice. See docs/denbi.md for setup, logs and recovery.

## Step 1 · Find and collect a protein family

Find global DeepClust clusters containing your proteins, then retrieve their sequences.

1. Choose your starting point: project groups, sequence identifiers, a known global centroid, or an exported cluster FASTA.
2. If you start with identifiers, queue mapping and inspect its completed result. Mapping finds membership; it does not retrieve the sequences.
3. Select mapped clusters and queue complete retrieval. Wait for COMPLETE in Jobs.

**Output:** one unaligned protein FASTA per cluster, membership and length summaries. **Next:** open Family analysis and select a retrieved FASTA. Retrieval uses the configured database on the analysis server; it does not download the entire DeepClust database to your laptop.

## Step 2 · Compare the sequences in your family

Split a broad family into more similar subclusters, or align an existing family directly.

1. Choose a task and an unaligned protein FASTA from a completed job or an upload.
2. For a broad family, start with Stricter subclusters. For a small family you can choose Alignment and tree immediately.
3. To align one subcluster, select that completed subclustering job as your input and choose a subcluster from its audited membership.

**Output:** audited subclusters, or a multiple sequence alignment and optional tree. **Next:** inspect a completed alignment in Alignments and trees. Subclustering and alignment are separate jobs; creating subclusters does not automatically align them.

## Step 3 · Inspect conservation and relationships

Explore a completed multiple sequence alignment, its conservation statistics and tree.

1. Select a completed alignment job. This page displays results; it does not run an alignment.
2. Start with JBrowse MSA for the tree and alignment together. Use its search and zoom controls to locate sequences and residues.
3. Try Alignment Viewer 2 for sorting by similarity to the reference, or Statistics for per-position and per-sequence evidence.

**Output:** interactive views and downloadable standalone alignment HTML. Complete FASTA, Newick and TSV files are in Jobs. Tree branch lengths represent divergence, not elapsed time; alignments and trees do not demonstrate experimental function.

## Optional · Predict and inspect protein structures

Request ESMFold models for selected proteins and inspect the returned structure evidence.

1. Select an unaligned protein FASTA, then paste exact sequence IDs to restrict the request. Start with one protein; a model job accepts at most 100 proteins.
2. Choose a configured route. External routes require your explicit sequence submission selection. Atlas accepts proteins up to 400 amino acids in this app.
3. Queue the job, then select its completed result below to inspect confidence and download PDB models.

**Output:** PDB models, residue evidence and a per-protein success/failure table. An optional annotation overlay must describe exactly the same sequence. A failed provider request remains visible; a model is never substituted for a different sequence.

## Optional · Predict domains and functional annotations

Obtain computational evidence for proteins with missing or incomplete annotation.

1. Select an unaligned protein FASTA and optionally restrict it to exact sequence IDs.
2. Choose InterProScan for domains/signatures, or installed eggNOG-mapper for orthology-based annotation. Only eggNOG uses the taxonomic scope controls.
3. For remote InterProScan, explicitly allow sequence submission. Queue the job, then inspect the completed annotation audit and feature tracks below.

**Output:** annotations, per-input no-hit records and positional Nightingale tracks where domain coordinates exist. These are predictions, not experimentally validated functions. Eukaryota/2759 is the default eggNOG scope, not a prokaryote-only workflow.

## Track progress and download results

Follow submitted work, inspect its outcome and obtain the complete result files.

1. Refresh the job list and choose a job by its task, timestamp and identifier.
2. QUEUED means waiting for a worker; RUNNING means a worker has claimed it. A running job may be waiting for disk, a mounted file or an external service.
3. COMPLETE means results were published; choose the next analysis or prepare a verified download. FAILED means inspect the error and diagnostics before retrying.

Only queued jobs can be cancelled here. The browser submits jobs and a separate worker executes them using the same research profile/workspace. Closing the browser does not cancel queued or running work. The list shows the 200 most recent jobs.

## Glossary

**Local group / HOG:** A family/group in the source project. HOG means hierarchical orthogroup. Its local ID is not automatically a global DeepClust cluster ID.

**Sequence identifier:** The exact ID of one protein; in FASTA it is the first whitespace-delimited part of the header after >.

**Global cluster / centroid:** A published DeepClust family, named by its centroid sequence ID. An existing member can map you to that family.

**Representative / medoid:** A protein chosen to represent a family; a medoid is selected using distances. Exact identifier mapping does not calculate a medoid. The separate homology search compares all supplied query proteins against global centroids.

**Subcluster:** A local, stricter grouping inside your selected sequences. It is specific to the source job and thresholds, and does not change global DeepClust.

**Mutual coverage:** The aligned region must cover the required fraction of both the member sequence and its assigned centroid.

**Reference / query:** The selected protein used to compare alignment members. It can be chosen explicitly and is not necessarily a biological representative.

**Job / worker:** A job is a saved request; a separate worker executes it. Refresh Jobs to see its status and access the outputs.

**Artifact / manifest:** An artifact is an output file. The manifest records its checksum, parameters and provenance for reproducibility.

**No annotation hit:** No supported annotation was found under the chosen method and database. This does not prove a protein has no function.

## Troubleshooting

- **Queued for a long time:** a worker may be busy, stopped or using a different profile. Refresh Jobs; ask the operator to check the worker and workspace. Low CPU use alone does not establish that a job stopped.
- **Unmatched ID:** check the exact accession/version and source prefix. Unmatched means no exact/declared-alias mapping, not no biological homologues.
- **No indexed cluster:** some mapped clusters lack a published row-group index. Keep their mapping audit; do not assume they can be retrieved here.
- **Missing tool or database:** an operator must configure that optional route. Opening this page does not install tools or download databases.
- **Alignment too large:** select a smaller audited subcluster or explicitly enable sampling; the application does not silently truncate a family.
- **Failed job:** read the error in Jobs, then expand technical details for the worker log and request. Correct the input/setup before submitting again.

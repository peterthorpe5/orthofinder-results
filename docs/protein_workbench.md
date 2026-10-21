# Protein research workbench — OrthoFinder results v0.12.0

The workbench starts from **your actual OrthoFinder memberships**, selected
DeepClust exports, or an ordinary protein FASTA. It ports the generic research
engine from the supplied E3 app v0.19.1 and adds an OrthoFinder-specific input
adapter. No E3 app installation or E3 resource is needed.

## What appears where

| Research section | Work available |
| --- | --- |
| OrthoFinder groups | Exact HOGs/orthogroups, largest local groups, species filters, member-ID text filters, complete membership context and sequence transfer. |
| Discover | Size-ranked global DeepClust catalogue; exact sequence-to-cluster mapping; selected cluster import; exact centroid lists; indexed extraction where the database is installed. |
| Family analysis | DIAMOND DeepClust with stricter identity/mutual coverage; MAFFT alignment and optional FastTree tree; optional search against global representative sequences. |
| Alignments and trees | JBrowse MSA, Alignment Viewer 2, alignment statistics, Auspice dataset download and optional embedded Auspice service. |
| Models | Explicit ESMFold requests, PDB downloads, residue evidence, C-alpha trace, confidence profile and sequence-matched feature overlays. |
| Annotations | Local InterProScan, local eggNOG-mapper or explicitly enabled remote EMBL-EBI InterProScan; Nightingale positional features where available. |
| Jobs | Durable status, complete verified artifacts, source records, tool commands, logs and failures. |
| Setup | Current workspace, tool availability and configured data sources; no database or model downloads. |

All tabular research previews include TSV and formatted Excel controls. Large
complete artifacts remain available through Jobs or the worker filesystem.
Charts use the application's existing deferred PDF export controls.

## Install and launch on your Mac

Apply the supplied update from `APPLY_AND_RUN.md`, then use the existing
OrthoFinder environment:

```bash
conda activate orthofinder_results
python -m pip install --editable '.[app]'
conda install -c conda-forge -c bioconda 'diamond>=2.2.1' mafft fasttree
```

Copy `config/protein_research.example.toml` outside the repository. A minimal
profile is enough for imported FASTAs and normal protein analysis:

```toml
workspace = "/Users/yourname/orthofinder_protein_research"
threads = 4
memory_gb = 4
```

Use an absolute path or `~` for workspace paths. TOML does not expand `$HOME`.
The app and its worker must use the **same profile and workspace**:

```bash
orthofinder-interrogation-app \
    --resource-dir /absolute/path/to/completed_resource \
    --research-config /absolute/path/research.toml
```

In a second terminal:

```bash
conda activate orthofinder_results
orthofinder-protein-workbench --config /absolute/path/research.toml worker
```

Open **Protein research workbench** in the sidebar. Its **OrthoFinder groups**
section is the natural starting point. The existing resource schema remains 5;
there is no need to rerun the resource-building Slurm job.

For protein-only work without an OrthoFinder resource:

```bash
orthofinder-protein-workbench --config /absolute/path/research.toml serve --port 8503
```

## Start from an OrthoFinder group

1. In **Explore one cluster** or a gene/protein result, choose **Open this group
   in protein research workbench**. The exact run, group system and hierarchy
   level are retained.
2. Choose **Prepare complete group inputs**. Alternatively, clear the exact-ID
   field and select the largest local groups, initially 200, at one hierarchy level.
3. Species filters and protein-ID text filters constrain the size-ranked search.
   Explicit group IDs override ranking and search filters.
4. Download the full input context as TSV or Excel. A local size-selection rank
   is not an E3 candidate-priority score; those scores belong to the E3 app.
5. Choose **Queue selected OrthoFinder proteins for global mapping** if the
   DeepClust sequence map and cluster index are installed on the worker.
6. Or choose **Import selected local protein sequences** to create a reusable
   local family input. After the worker finishes, select `family.fasta` under
   **Family analysis**.

The transfer uses embedded `protein_sequences` if present, or the existing
`--terminal-motif-parquet` sequence sidecar. Membership-only resources do not
contain protein sequences merely because they contain protein IDs. Supply the
sequence sidecar or upload a matching FASTA when needed.

Sequence joins use both exact member ID and species label. Missing sequences,
duplicate sequence records and protein IDs reused across species stop FASTA
transfer. Group selection stops at 10,000 distinct protein IDs and never silently
cuts a group to fit. For larger selections, export a suitable FASTA near the data
and use the CLI or split the selection deliberately.

For a complete ordinary FASTA exceeding the browser limit, stream a private
snapshot through the CLI:

```bash
orthofinder-protein-workbench --config /absolute/path/research.toml import-fasta \
    --fasta /absolute/path/local_family.fasta \
    --source-label "Complete local protein family export"
```

This produces a reusable `family.fasta` after the worker completes. It does not
claim that an ordinary local FASTA has independently verified global membership.

## DeepClust access and selected exports

See [deepclust_exchange.md](deepclust_exchange.md). The practical laptop route is:

1. Determine exact global centroids using installed mapping/index files or a
   database holder. An OrthoFinder HOG ID is not a DeepClust centroid.
2. Download `centroids.id.txt` from Discover: one unchanged centroid per line.
3. Have the extraction machine run the official tool in per-cluster output mode.
4. Transfer only the selected FASTAs to the Mac.
5. Choose **Discover → Import extracted cluster**, declare its source and centroid,
   optionally supply the independently obtained member count, then queue import.
6. Select the resulting `cluster.fasta` under **Family analysis**.

Import validates protein symbols, unique identifiers and the presence of the exact
centroid. A supplied member count must match. Without an independent count, the
audit records `COUNT_NOT_SUPPLIED`: it does not claim global completeness from
FASTA content alone. Original cluster IDs and accession versions are retained.

The public data release is a large multipart archive. This update neither
downloads it automatically nor supplies access to a database holder. Global
size ranking requires the cluster index; exact member mapping requires the
sequence map. Imported FASTAs do not supply a global catalogue.

## Reclustering, alignments and visual interpretation

**Stricter subclusters** runs native `diamond deepclust`, then exact
representative-to-member realignment. Initial exploratory settings are 70%
identity, 80% mutual coverage and E-value at most 1e-10. Both coverages and a
bitscore greater than 20 are audited. Settings should be chosen for the family.
Assignments failing the exact audit are retained as singletons.

These are representative-based subclusters. They do not guarantee that every
pair of members passes the thresholds, and are not new inferred orthogroups.
The original OrthoFinder and global DeepClust memberships remain separate.

For **Alignment and tree**, select the completed family's FASTA. An audited
`subcluster_id` restricts alignment to that subcluster. Select an exact reference
ID, align the whole selection or explicitly request deterministic sampling, and
enable tree inference if needed. Sampled selections are labelled and exported;
oversized selections otherwise fail instead of being silently truncated.

MAFFT results include per-column gaps, consensus and entropy; reference-to-member
identity, overlap and raw BLOSUM62 scores; and selection/alignment summaries.
JBrowse MSA and Alignment Viewer 2 are bundled locally. Standalone alignment HTML
downloads work without a CDN. Alignment Viewer 2 provides BLOSUM62-to-query sorting.
Nightingale shows exact sequence coordinates and available positional domains.

FastTree uses a protein WAG model with gamma rate correction. For substantially
divergent mixed families, align biologically justified subclusters before drawing
evolutionary conclusions. Divergence branch lengths are not dates, and the
displayed root is not asserted to be an ancestral protein.

## Auspice, focus and streamtrees

The worker creates `tree.nwk` and Auspice v2 `auspice.json`, and stages the latter
as `WORKSPACE/datasets/JOB_ID.json`. Species labels from imported local inputs and
the latest audited subcluster memberships are traced through checksummed source
jobs. Configured reviewed taxonomy can add exact taxon IDs.

Maximal multi-tip clades whose tips share a subcluster receive `subcluster`
branch labels. A non-monophyletic sequence-similarity subcluster can yield several
distinctly labelled parts. Mixed or unknown clades remain ordinary branches.
The topology is unchanged. If no supported multi-tip clades exist, streamtree
labels are omitted. For a tree of one subcluster, one stream can be expected.

Auspice is a separate service. Its recent focus-on-selected and streamtree UI is
available when using a suitable current service. For the verified 3.1.0 package,
use Node 24 or 26 and npm 11:

```bash
npm install --prefix "$HOME/auspice-service" auspice@3.1.0
"$HOME/auspice-service/node_modules/.bin/auspice" view \
    --datasetDir /absolute/path/orthofinder_protein_research/datasets \
    --host 127.0.0.1 --port 4000
```

Set `auspice_url = "http://localhost:4000"` in the profile and restart/reload.
The browser must be able to reach that URL. An app accessed remotely needs a
corresponding reachable service/proxy address. Choose **Auspice tree** under
**Alignments and trees**; use the Auspice sidebar to focus selected tips or enable
streamtrees using the `subcluster` branch-label key. Without a service, download
the JSON for your own Auspice installation. Auspice itself is not bundled.

## Missing structures and functional annotation

Prediction is an explicit action on selected proteins. Available routes are the
public Atlas endpoint, an operator-deployed Hugging Face/custom endpoint, and
optional local ESMFold. Endpoint availability is not guaranteed. A Hugging Face
model repository by itself does not establish a working inference API.

The 400-residue bound belongs to the public endpoint route and the default custom
endpoint setting, rather than being a universal ESMFold algorithm limit. Longer
proteins require a configured supported endpoint limit or sufficient local
resources. Proteins are not automatically truncated to meet a service limit.

Local prediction is optional:

```bash
python -m pip install --editable '.[folding]'
```

The first explicit local prediction downloads the pinned checkpoint. GPU/CPU
memory requirements depend on sequence length. Do not install model weights just
to browse completed results. Successful jobs provide PDB, residue evidence,
confidence summaries and interactive C-alpha geometry; these are computational
models, not experimental structures or a claim of validated ligandability.

For missing annotation, install/configure local InterProScan and its databases,
or local eggNOG-mapper with its own compatible database. The default eggNOG scope
is **Eukaryota (2759)** and can be changed. Remote EMBL-EBI InterProScan requires
a contact email and explicit remote submission. Protein data are sent only when
that route is requested. No external predictions were submitted during development.

Nightingale displays InterProScan's positional evidence. Whole-protein eggNOG
assignments are not turned into invented residue coordinates. No-hit proteins and
individual provider failures remain explicit. Model/annotation overlays require
the same exact sequence SHA-256, preventing accidental isoform overlays.

## Workers, limits and provenance

Scientific work runs outside Streamlit through a SQLite-backed queue. Jobs keep
validated request snapshots, commands, tool versions, source paths, timestamps,
failure diagnostics and SHA-256 hashes of output artifacts. Next stages accept
only completed source artifacts whose hashes still match.

Use `doctor`, `jobs`, `cancel --job-id ID`, `worker --once`, or
`worker --max-jobs 1` through `orthofinder-protein-workbench --config PROFILE`.
Cancellation applies to queued work. A crashed worker can leave RUNNING work;
inspect its logs and partial files, stop that worker, then submit a new job.
No automatic restart of ambiguous partial scientific jobs is claimed.

Slurm example:

```bash
export RESEARCH_ENV_PREFIX=/absolute/path/to/conda/envs/orthofinder_results
export ORTHOFINDER_RESEARCH_CONFIG=/absolute/path/research.toml
sbatch scripts/research_worker.slurm
```

The supplied script uses `--account=barton --partition=general`. Its 4 CPUs and
16 GB allocation must match the profile and actual tools. It expects cluster
`TMPDIR`; edit it for your site if necessary. Queue storage shared across web and
HPC processes must support reliable POSIX locking. On separate machines, copy
completed cluster FASTAs explicitly; this application does not submit SSH or
Slurm jobs remotely. Local laptop and HPC workspaces can be separate.

Default bounds: browser FASTA 20 MB; one request JSON 2 MB; at most 1,000 exact
global clusters per extraction job; 8 million extracted sequences; 5,000 aligned
sequences and 25 million alignment cells. Large artifact downloads are not loaded
wholesale into the browser. These bounds are configurable where indicated in
the example profile. Bounds are safeguards, not hardware capacity promises.

## Validation and development

```bash
python -m pip install --editable '.[dev]'
python -m coverage run --module pytest -q
python -m coverage report
python -m ruff check src tests scripts
python -m pycodestyle src/orthofinder_interrogation_app/research
python -m pydocstyle src/orthofinder_interrogation_app/research
```

For a real small external-tool check:

```bash
python scripts/smoke_protein_workbench.py \
    --config /absolute/path/research.toml --output /new/empty/smoke_workspace
```

The smoke uses four synthetic sequence variants and requests no network provider.
Optional integration tests use `WORKBENCH_SCIENTIFIC_CONFIG` and
`DEEPCLUST_RETRIEVAL_REPO`. Viewer sources, exact npm versions, checksum inventory
and third-party notices are included under `frontend/` and packaged viewer data.

## Primary references

- [DeepClust data release](https://objectstore.hpccloud.mpcdf.mpg.de/deepclust/index.html)
- [Official extraction tool](https://github.com/drostlab/deepclust_dataretrieval)
- [DIAMOND documentation](https://github.com/bbuchfink/diamond/wiki)
- [JBrowse MSA](https://gmod.org/JBrowseMSA/)
- [Nightingale](https://github.com/ebi-webcomponents/nightingale)
- [Alignment Viewer 2](https://fast.alignmentviewer.org/?sort-by=blosum-score-to-query)
- [Auspice streamtree requirements](https://docs.nextstrain.org/projects/auspice/en/latest/advanced-functionality/streamtrees.html)
- [Auspice March 2026 update](https://nextstrain.org/blog/2026-03-12-auspice-updates)
- [ESMFold](https://github.com/facebookresearch/esm)
- [Hugging Face ESM](https://huggingface.co/docs/transformers/model_doc/esm)
- [InterProScan](https://interproscan-docs.readthedocs.io/)
- [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)

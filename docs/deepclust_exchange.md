# Selected DeepClust exports

For your **de.NBI/Sabrina deployment**, DeepClust is already on the server.
Use the [server launcher and local-retrieval workflow](denbi.md): configure the
installed index, sequence map and sequences there, and retrieve selected
clusters with the app's server worker. This exchange guide describes the
alternative export/import route when transferring selected FASTAs between
machines; a different database holder is not required for Sabrina.

The official `drostlab/deepclust_dataretrieval` tool extracts clusters efficiently
from an **already extracted DeepClust database on the extraction machine**. Its
Git repository contains code. It does not provide a public remote cluster API.
The public download listing currently supplies the Parquet release in large
multipart archives. Installing the tool on your Mac does not bypass that data
requirement.

The recommended workflow is to obtain selected per-cluster exports from someone
who already hosts the database, then transfer those FASTAs. The OrthoFinder
workbench supports the entire import and analysis side of that exchange.

## 1. Resolve the cluster IDs

An exact global centroid can look like `NR_UGC79169.1`. Keep the `NR_`, `UNIPROT_`
or other source namespace, case and accession version. An OrthoFinder ID such as
`N0.HOG1` is a different identifier and must not be substituted as a centroid.

The official sequence map has columns **SEQID** and **CLUSTERID**. The workbench
accepts these names and the earlier documented normalised aliases. Submitted
UniProt/NR aliases are tested explicitly; no arbitrary version stripping, fuzzy
equivalence or inferred orthology is used for exact global membership.

With the index and sequence-map files installed, queue mapping under either
**OrthoFinder groups** or **Discover → Sequence identifiers**. Multiple global
matches and unmatched inputs remain in the audit. A homology hit against a global
representative is separate evidence; it does not prove exact cluster membership.

Without those installed files, ask the database holder to resolve the exact
source accessions using their sequence map, or obtain their mapping table. The
app cannot infer a correct global cluster merely from a local HOG number.

## 2. Export a centroid request list

Under **Discover**, enter the chosen global centroids in **Global cluster IDs to
retrieve in full**, then choose **Download centroid list for official extraction
tool**. `centroids.id.txt` is plain text with **no heading**, one centroid per
line, for example:

```text
NR_UGC79169.1
UNIPROT_P12345
```

The second line illustrates naming only; use actual resolved centroids for your
own proteins. Downloading a list does not establish that an arbitrary entered ID
exists in the release.

## 3. Extraction on the database-holder machine

Clone the tool once and record its revision. The wrapper was compatibility-tested
against upstream commit `b406fec8460d8044697c057bcb8c34f63549d98e`.

```bash
git clone https://github.com/drostlab/deepclust_dataretrieval.git \
    "$HOME/software/deepclust_dataretrieval"
git -C "$HOME/software/deepclust_dataretrieval" rev-parse HEAD
```

Install DuckDB and PyArrow in the extraction environment. The global data files
must already exist. From the updated OrthoFinder repository:

```bash
bash scripts/extract_deepclust_selected.sh \
    --retrieval-repo "$HOME/software/deepclust_dataretrieval" \
    --centroid-file /absolute/path/centroids.id.txt \
    --database /shared/deepclust/joined_with_index_RowGroupFinal.parquet \
    --index /shared/deepclust/clust_index_RowGroup.parquet \
    --output /absolute/path/new_selected_exports \
    --threads 1 \
    --python /absolute/path/to/environment/bin/python
```

The wrapper preserves upstream's `-path_to_DCD`, `-path_to_output` and
`-path_to_index` spelling: those three flags use **one leading dash**. It adds the
trailing output-directory separator required by upstream's filename construction.
It uses `--per-clust-output 1`, preserving original sequence IDs, and processes
one requested cluster at a time. It refuses non-empty output directories.

macOS requires one upstream worker; the default is `--threads 1`. Linux can use
more, but even selected row groups may decompress to substantial memory.
Upstream builds/reuses a `persistent` DuckDB index beside the configured Parquet
index. That directory must be writable on its first run; keep this index tied to
the same release and a compatible DuckDB version. Request a new index location
or use the app's native read-only Parquet extraction if those constraints matter.

The output contains files such as `NR_UGC79169.1.fa`. Ask the provider to supply:

- the exact centroid and release/source label;
- the per-cluster FASTA with unchanged headers;
- the independently obtained published member count, when possible;
- file checksums and the extractor revision.

This project does not vendor the upstream extractor. It includes a wrapper and
a deterministic optional integration test. Review upstream's own licence and
citation information when installing it.

If the database is already mounted on the same machine as the worker, configure
its `index`, `sequences` and `sequence_map` paths instead. Native extraction reads
only the physical row groups listed for the exact selected centroids, checks
every extracted cluster against `NMEMBER`, and avoids creating upstream's
persistent index. It still requires installed global data.

## 4. Transfer only selected FASTAs

For example, from your Mac, using your actual approved host/path:

```bash
mkdir -p "$HOME/data/deepclust_selected"
rsync --archive --verbose \
    user@database-holder:/absolute/path/new_selected_exports/ \
    "$HOME/data/deepclust_selected/"
```

No automatic SSH connection, remote job submission or database download is made
by the workbench. Transfer only files the provider is authorised to share.

## 5. Import in the app

Choose **Discover → Import extracted cluster**. Upload the **one-cluster** FASTA,
enter its exact global centroid, source label and optionally its published
member count. Queue the import, then run the worker with the same profile.

Checks include valid protein symbols, unique sequence IDs, the exact centroid
among the headers and agreement with a supplied expected count. Combined
upstream output mode deliberately prefixes/renames IDs; use per-cluster mode
instead. Files with a missing centroid, duplicates or count mismatches fail
visibly and are not published as successful imports.

The result provides `cluster.fasta`, `members.tsv`, `clusters.tsv` and
`cluster_summary.tsv`, with sequence/input hashes and provenance. When no count
was supplied, it records `COUNT_NOT_SUPPLIED`. A count match checks the supplied
number; it is not independent verification that every member belongs to that
global cluster. The provider's export and release provenance remain necessary.

## 6. Large imports through the CLI

Browser uploads are limited to 20 MB. Larger exports can be snapshotted directly
without loading the whole file into the web process:

```bash
orthofinder-protein-workbench --config /absolute/path/research.toml import-deepclust \
    --fasta "$HOME/data/deepclust_selected/NR_UGC79169.1.fa" \
    --cluster-id NR_UGC79169.1 \
    --source-label "Provider name; DeepClust release; extraction date" \
    --expected-members 1960408
orthofinder-protein-workbench --config /absolute/path/research.toml worker --once
```

The example count is from upstream's documented example; replace it with the
provider's count for your export or omit `--expected-members` when unknown.
The CLI prints a job ID. A worker consumes the checksummed private input
snapshot, so subsequently changing the original FASTA cannot alter the queued
analysis. There must be room for the input snapshot, validated outputs and
temporary indexing; this is not a zero-copy import.

## 7. Continue analysis

After import completes, select `cluster.fasta` under **Family analysis**. Create
stricter subclusters, inspect their exact membership and sizes, and align a
chosen subcluster or an explicitly chosen sample. Sequence limits still apply
after import. Importing a million proteins does not imply that aligning them all
on a laptop is feasible.

Reclustering is a new analysis with its own settings. It does not rewrite the
global release or the source OrthoFinder groups. Keep the global centroid and
original source mapping when reporting local subcluster results.

Primary sources: [upstream extraction documentation](https://github.com/drostlab/deepclust_dataretrieval)
and [the DeepClust data release](https://objectstore.hpccloud.mpcdf.mpg.de/deepclust/index.html).

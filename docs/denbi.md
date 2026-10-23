# OrthoFinder on de.NBI / Sabrina

The app and analysis worker run on Sabrina, beside its installed DeepClust
database. Your Mac displays the web interface through an SSH tunnel. Selected
cluster extraction, DIAMOND reclustering, MAFFT alignment and FastTree analysis
run on the server. There is no need to transfer DeepClust to your Mac or obtain
exports from a different database holder for this deployment.

This launcher adapts the supplied E3 Sabrina launcher v1.0.2. It retains SSH
alias authentication, persistent tmux services, readiness checks, retained logs
and reconnect behaviour. It has separate OrthoFinder services and defaults:

| Setting | OrthoFinder default |
| --- | --- |
| SSH alias | `denbi` |
| Mac address | `http://127.0.0.1:8510` |
| Server bind | `127.0.0.1:8506` |
| tmux session | `orthofinder-research` |
| tmux windows | `app`, `worker` |
| Conda environment | `/home/ubuntu/data/pete/software/miniforge3/envs/orthofinder_results` |
| Server checkout | `/home/ubuntu/data/pete/software/orthofinder-results` |
| Research profile | `/home/ubuntu/data/orthofinder_research/config/research.toml` |
| Job workspace | Set by that profile; recommended `/home/ubuntu/data/orthofinder_research/workspace` |
| Completed resource | `/home/ubuntu/data/orthofinder_research/resources/current` |
| Logs, temporary files and caches | Under `/home/ubuntu/data/orthofinder_research` |

The environment and resource defaults are **new installation destinations**,
not a claim that those paths already exist. The Miniforge and existing E3
profile paths below come from your supplied launcher. The exact completed
OrthoFinder resource path must be supplied from your server or linked at
`resources/current`. A completed resource contains `orthofinder_results.duckdb`;
an original OrthoFinder `Results_*` directory is a different input.

## First installation: transfer the updated code from your Mac

Apply the download to your local Git checkout using `APPLY_AND_RUN.md` first.
No server-side GitHub access is required. From the root of the updated checkout:

```bash
ssh -T -o ClearAllForwardings=yes denbi \
    'mkdir -p /home/ubuntu/data/pete/software/orthofinder-results'
rsync -av --exclude .git --exclude .venv --exclude __pycache__ \
    --exclude .pytest_cache --exclude .pytest_work --exclude .coverage \
    --exclude build --exclude dist --exclude '*.egg-info' \
    -e 'ssh -o ClearAllForwardings=yes' ./ \
    denbi:/home/ubuntu/data/pete/software/orthofinder-results/
```

Rsync does not use `--delete`. Keep datasets and job workspaces outside the
checkout. Before an update to an existing installation, let the OrthoFinder
queue finish and stop its services as described under **Update the app**.

## First installation: create the environment on Sabrina

Connect from your Mac:

```bash
ssh -o ClearAllForwardings=yes denbi
```

Then run these commands **on Sabrina**:

```bash
OF_CONDA=/home/ubuntu/data/pete/software/miniforge3/bin/conda
OF_ENV=/home/ubuntu/data/pete/software/miniforge3/envs/orthofinder_results
OF_REPO=/home/ubuntu/data/pete/software/orthofinder-results
mkdir -p /home/ubuntu/data/orthofinder_research/tmp
mkdir -p /home/ubuntu/data/orthofinder_research/cache/pip
export TMPDIR=/home/ubuntu/data/orthofinder_research/tmp
export PIP_CACHE_DIR=/home/ubuntu/data/orthofinder_research/cache/pip
"$OF_CONDA" create --prefix "$OF_ENV" -c conda-forge -c bioconda \
    python=3.11 tmux 'diamond>=2.2.1' mafft fasttree
cd "$OF_REPO"
"$OF_ENV/bin/python" -m pip install --editable '.[app]'
```

If this exact environment already exists, use `conda install --prefix "$OF_ENV"`
with the package list instead of `conda create`. The optional local ESMFold,
InterProScan and eggNOG installations remain separate scientific resources;
the normal app does not download weights or annotation databases during launch.

## Reuse the working DeepClust database settings

The attached E3 launcher uses this existing profile:

`/home/ubuntu/data/protein_research/config/research.toml`

Create a separate OrthoFinder profile from it, **on Sabrina**:

```bash
OF_ENV=/home/ubuntu/data/pete/software/miniforge3/envs/orthofinder_results
"$OF_ENV/bin/orthofinder-protein-workbench" copy-profile \
    --source /home/ubuntu/data/protein_research/config/research.toml \
    --output /home/ubuntu/data/orthofinder_research/config/research.toml \
    --workspace /home/ubuntu/data/orthofinder_research/workspace
"$OF_ENV/bin/orthofinder-protein-workbench" \
    --config /home/ubuntu/data/orthofinder_research/config/research.toml doctor
```

`copy-profile` preserves every data, executable, provider and resource-limit
setting, changing only the workspace. It does not copy databases, modify E3's
profile, copy jobs or overwrite an existing profile. If a profile already exists,
inspect and keep it. Use an actual alternate E3 profile path if the other chat's
working setup differs from the attachment.

Check the `doctor` result: `index`, `sequences` and `sequence_map` should address
the server's installed DeepClust files. Retain any reviewed taxonomy and global
centroid DIAMOND database paths already configured. Default executable names
are resolved with the OrthoFinder environment first on `PATH`; explicit E3 tool
paths are preserved. Dataset path globs can be valid even when `doctor`'s simple
file-existence field is false; verify the actual mounted dataset.

Both apps may read the same database. Their profiles must use different job
workspace roots so the correct worker handles each queue. Set each profile's
threads and memory allowance to suit concurrent use of the server.

## Point to the completed OrthoFinder resource

Keep the resource on the server and use its actual absolute path with the
launcher. For a one-command daily launch, create the default link once, on
Sabrina, replacing the example source path with the real completed resource:

```bash
mkdir -p /home/ubuntu/data/orthofinder_research/resources
ln -s /absolute/server/path/to/completed_resource \
    /home/ubuntu/data/orthofinder_research/resources/current
```

`ln -s` does not replace an existing link or directory. If the link already
exists, inspect it rather than overwriting it. The launcher accepts a completed
resource directory or its `orthofinder_results.duckdb` file. No resource rebuild
is required for this application update; the resource schema remains 5.

If you only have the completed resource on your Mac, copy that resource to the
server using rsync. DeepClust itself is already there:

```bash
rsync -av -e 'ssh -o ClearAllForwardings=yes' \
    /absolute/mac/path/to/completed_resource/ \
    denbi:/home/ubuntu/data/orthofinder_research/resources/my_run/
```

A local sequence Parquet may be embedded at
`evidence/protein_sequences.parquet`, in which case the app detects it. Otherwise
retain your existing sidecar with the launcher's `--sequence-parquet` option.
This sidecar supplies **local OrthoFinder protein sequences**; it is distinct
from the global DeepClust `sequences` setting.

## Daily launch from your Mac

Exit the server shell. Run this from the local updated repository:

```bash
bash scripts/launch_orthofinder_denbi.sh
```

Or use explicit server paths without creating the resource link:

```bash
bash scripts/launch_orthofinder_denbi.sh \
    --resource-dir /absolute/server/path/to/completed_resource \
    --sequence-parquet /absolute/server/path/to/protein_sequences.parquet
```

Omit `--sequence-parquet` when it is not needed. All path arguments refer to
Sabrina. Other named options are `--ssh-host`, `--env-prefix`, `--research-config`,
`--state-dir`, `--taxonomy-map`, `--local-port` and `--remote-port`. To make
custom paths persistent, set `ORTHOFINDER_DENBI_RESOURCE`,
`ORTHOFINDER_DENBI_SEQUENCES`, `ORTHOFINDER_DENBI_TAXONOMY`,
`ORTHOFINDER_DENBI_ENV`, `ORTHOFINDER_DENBI_CONFIG`, `ORTHOFINDER_DENBI_STATE`
and/or `ORTHOFINDER_DENBI_HOST` in your Mac shell configuration.
Use the same custom settings for status/restart commands.

The launcher validates the resource/profile/sidecars, starts or reuses the two
owned tmux windows, waits for the app and opens your Mac browser. The app binds
to server localhost; no public firewall port needs opening. Existing E3 services
and their tunnel have independent names and ports.

## Guided scientific workflow

1. In the sidebar open **Protein research workbench → Start here**. Use the
   recommended **OrthoFinder groups** route or start from IDs, a global centroid
   or a protein FASTA.
2. For local groups, start with one HOG/orthogroup at the intended tree level.
   **Prepare complete group inputs**, review membership, then import local
   sequences or queue exact global mapping. **Track this job** opens its status.
3. For DeepClust, use **Discover** to rank global clusters or review the mapping.
   Select mapped centroids and use **Use selected mapped clusters for retrieval**.
   **Queue complete cluster retrieval** extracts the selected families from the
   installed database on Sabrina. The official extraction/export tool remains
   an alternative, rather than a requirement for this server-local route.
4. In **Jobs**, open a complete family with **Analyse this family**. Apply
   stricter identity and mutual coverage when appropriate. Use the completed
   subclustering result's next-step button and select an audited subcluster for
   alignment/tree analysis. These are explicit queued steps.
5. Inspect **Alignments and trees** in JBrowse MSA, Alignment Viewer 2, statistics
   or Auspice. Use **Models** and **Annotations** when needed. Optional external
   providers still require an explicit submission choice; local providers run
   on Sabrina if installed and configured.
6. Download selected results through **Jobs** or rsync large artifacts from the
   server workspace. For example, use the actual path displayed in Jobs:

```bash
rsync -av -e 'ssh -o ClearAllForwardings=yes' \
    denbi:/absolute/server/workspace/path/to/result.fasta ./
```

Every page has its aim, steps, expected outputs, next step and question-mark
field help. **Help & glossary** contains a downloadable user guide and common
failure explanations. Selecting a next-step button fills the next form; it does
not submit another scientific job.

## Status, reconnect and logs

From the Mac:

```bash
bash scripts/launch_orthofinder_denbi.sh --status
bash scripts/launch_orthofinder_denbi.sh --stop-tunnel
bash scripts/launch_orthofinder_denbi.sh
```

Stopping the tunnel or closing the Mac terminal leaves the app and server
worker running. Re-running the launcher reconnects without duplicating either
service. `--no-browser` prints the URL for manual use.

For app-only settings changes, preserve the running worker:

```bash
bash scripts/launch_orthofinder_denbi.sh --restart-app
```

A changed resource or remote port requires `--restart-app`. A changed research
profile/environment is refused while the old worker is live: finish its jobs
and replace the services using the update procedure below. This prevents a GUI
and worker from silently operating on different queues/settings.

On Sabrina, diagnostics live here:

```bash
tmux list-windows -t orthofinder-research
tmux capture-pane -p -S -80 -t orthofinder-research:=app
tail -n 80 /home/ubuntu/data/orthofinder_research/logs/orthofinder_server.log
tail -n 80 /home/ubuntu/data/orthofinder_research/logs/research_worker.log
```

Use the environment's `bin/tmux` full path if tmux is not in your interactive
server `PATH`. The structured application log is `logs/application.log`.
The launcher refuses to replace windows that lack its ownership marker.

- **Missing executable/profile/resource:** finish the first-install steps or
  pass the actual absolute server paths. The launcher does not guess datasets.
- **Mac port occupied:** choose `--local-port 8511`; the remote port can stay 8506.
- **Remote port occupied:** choose a free `--remote-port` and `--restart-app`;
  keep using that port on subsequent launches.
- **Tunnel failure:** run `--stop-tunnel` with the same host/ports, then reconnect.
- **Queued jobs:** inspect the worker window/log and check app and worker use
  the same OrthoFinder profile. The launcher starts one persistent worker.
- **Missing IDs:** exact mapping needs matching accessions/versions or declared
  aliases; absence of an exact map is not evidence of absent homology.
- **Database I/O or provider failure:** inspect the failed job's technical
  details. Reconnects do not resubmit failed jobs or automatically retry analyses.

## Update the app

First inspect **Jobs** and let running/queued OrthoFinder work finish. Then, on
Sabrina, stop only this session:

```bash
/home/ubuntu/data/pete/software/miniforge3/envs/orthofinder_results/bin/tmux \
    kill-session -t orthofinder-research
```

This command stops the OrthoFinder worker as well as its app. Do it only after
its jobs have finished. It does not affect the separate `e3-research` session.
Transfer the updated code with the earlier rsync command and reinstall on
Sabrina:

```bash
cd /home/ubuntu/data/pete/software/orthofinder-results
/home/ubuntu/data/pete/software/miniforge3/envs/orthofinder_results/bin/python \
    -m pip install --editable '.[app]'
```

Preserve the research profile, job workspace, database and resource. Run the
Mac launcher again. A worker interrupted during a job may leave that job
marked RUNNING; inspect its diagnostics and outputs before arranging recovery.
The launcher never silently requeues an interrupted scientific analysis.

The release is tested locally with controlled SSH/tmux boundaries. It has not
connected to, installed on, or verified your live Sabrina server.

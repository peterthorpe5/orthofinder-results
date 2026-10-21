# Protein workbench scripts

- `research_worker.slurm`: one queued task using account barton and partition general.
- `extract_deepclust_selected.sh`: portable wrapper for the separately installed
  official extractor and existing local global database.
- `smoke_protein_workbench.py`: real small DIAMOND/MAFFT/FastTree workflow.
- `build_viewer_smoke.py`: produce synthetic offline viewer HTML for browser checks.

All data/profile controls are named options. Read the workbench and DeepClust
exchange guides before using global data or HPC jobs.

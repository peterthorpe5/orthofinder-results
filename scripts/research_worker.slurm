#!/usr/bin/env bash
#SBATCH --job-name=protein-research
#SBATCH --account=barton
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=protein-research-%j.log

# Submit with environment variables already exported; no shell-embedded Python.
set -euo pipefail
: "${RESEARCH_ENV_PREFIX:?Set RESEARCH_ENV_PREFIX to the installed OrthoFinder environment}"
: "${ORTHOFINDER_RESEARCH_CONFIG:?Set ORTHOFINDER_RESEARCH_CONFIG to an absolute TOML profile path}"
: "${TMPDIR:?The cluster must supply per-job local TMPDIR}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
exec "$RESEARCH_ENV_PREFIX/bin/orthofinder-protein-workbench" \
    --config "$ORTHOFINDER_RESEARCH_CONFIG" worker --max-jobs 1

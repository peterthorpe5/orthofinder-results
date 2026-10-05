#!/bin/bash

# first build the resource

    PACKAGE_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/orthofinder-results"
    CONDA_ENV="orthofinder_results"
    RESULTS_DIR="/home/pthorpe001/data/2026_E3_protac/SSD_back_up_July_2026/Erin_Butterfield_data/Main_folder/OrthoFinder/Results_Feb26"
    EXPRESSION_MANIFEST="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/expression_atlas_rebuild_v0_5_1_20260804/manifests/e3_workflow_expression_resources.tsv"
    RUN_ID="results_feb26_motif_expression_v0_11_0"
    OUTPUT_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/orthofinder_results_resources"
    LOG_DIR="${OUTPUT_ROOT}/slurm_logs/${RUN_ID}"

    test ! -e "${OUTPUT_ROOT}/${RUN_ID}" || {
        echo "STOP: output already exists: ${OUTPUT_ROOT}/${RUN_ID}" >&2
        exit 2
    }
    mkdir -p "$LOG_DIR"
    export ORTHOFINDER_EXPRESSION_MEMORY_MB=32768

    JOB_ID=$(sbatch --parsable \
        --job-name=of_motif_rna \
        --partition=barton \
        --cpus-per-task=13 \
        --mem=96G \
        --time=2-00:00:00 \
        --output="${LOG_DIR}/orthofinder_results_%j.out" \
        --error="${LOG_DIR}/orthofinder_results_%j.err" \
        "${PACKAGE_ROOT}/slurm/dispersion_benchmark.sbatch" \
        "$PACKAGE_ROOT" \
        "$CONDA_ENV" \
        "$RESULTS_DIR" \
        - \
        - \
        "$EXPRESSION_MANIFEST" \
        "$RUN_ID" \
        "${OUTPUT_ROOT}/${RUN_ID}" \
        --distance-max-members 250 \
        --benchmark-controls-per-group 3 \
        --benchmark-bootstrap-resamples 1000)

    echo "Submitted job: $JOB_ID"
    echo "Slurm stdout: ${LOG_DIR}/orthofinder_results_${JOB_ID}.out"
    echo "Slurm stderr: ${LOG_DIR}/orthofinder_results_${JOB_ID}.err"
    squeue -j "$JOB_ID"

# NEW VERSION for the app. 

    RESOURCE_DIR="$HOME/Downloads/results_feb26_dispersion_benchmark_v0_9_0"
    SIDECAR="$HOME/Downloads/orthofinder_results_sequence_sidecars/results_feb26_terminal_motif_sequences_v0_10_1.parquet"
    APP_LOG_DIR="$HOME/orthofinder_results_app_logs"

    mkdir -p "$APP_LOG_DIR"

    orthofinder-interrogation-app \
        --resource-dir "$RESOURCE_DIR" \
        --terminal-motif-parquet "$SIDECAR" \
        --server-port 8502 \
        --log-file "$APP_LOG_DIR/terminal_motif_v0_10_1.log"

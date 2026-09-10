

    PROJECT_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac"
    REPO_ROOT="$PROJECT_ROOT/orthofinder-results"
    RESULTS_DIR="$PROJECT_ROOT/SSD_back_up_July_2026/Erin_Butterfield_data/Main_folder/OrthoFinder/Results_Feb26"
    RESOURCE_ROOT="$PROJECT_ROOT/orthofinder_results_resources"
    RUN_ID="results_feb26_e3_precursor_v0_8_0"
    OUTPUT_DIR="$RESOURCE_ROOT/$RUN_ID"
    LOG_DIR="$RESOURCE_ROOT/slurm_logs/$RUN_ID"

    test -d "$RESULTS_DIR" || echo "STOP: OrthoFinder results directory is missing"
    test -x "$REPO_ROOT/slurm/e3_precursor.sbatch" || echo "STOP: wrapper is missing"
    test -f "$REPO_ROOT/src/orthofinder_interrogation_app/data/e3_seed_catalogue.tsv" || echo "STOP: E3 catalogue is missing"
    test ! -e "$OUTPUT_DIR" || echo "STOP: choose a new RUN_ID; output already exists"
    mkdir -p "$LOG_DIR"


    JOB_ID=$(sbatch --parsable \
        --account=barton \
        --partition=barton \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=9 \
        --mem=64G \
        --time=48:00:00 \
        --job-name=of_e3_all \
        --export=ALL,TMPDIR=/tmp \
        --output="$LOG_DIR/orthofinder_results_%j.out" \
        --error="$LOG_DIR/orthofinder_results_%j.err" \
        "$REPO_ROOT/slurm/e3_precursor.sbatch" \
        "$REPO_ROOT" \
        orthofinder_results \
        "$RESULTS_DIR" \
        - \
        "$RUN_ID" \
        "$OUTPUT_DIR" \
        --distance-max-members 250)

    echo "Submitted job: $JOB_ID"
    echo "Slurm stdout: $LOG_DIR/orthofinder_results_${JOB_ID}.out"
    echo "Slurm stderr: $LOG_DIR/orthofinder_results_${JOB_ID}.err"
    squeue -j "$JOB_ID"

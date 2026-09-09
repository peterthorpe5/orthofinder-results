#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    echo "Usage: $0 PACKAGE_ROOT CONDA_ENV RESULTS_DIR FOCUS_OR_DASH RUN_ID PERSISTENT_OUTPUT [safe analysis options]" >&2
}

[[ $# -ge 6 ]] || { usage; exit 2; }
PACKAGE_ROOT=$1
CONDA_ENV=$2
RESULTS_SOURCE=$3
FOCUS_SOURCE=$4
RUN_ID=$5
PERSISTENT_OUTPUT=$6
shift 6

[[ -x "$PACKAGE_ROOT/run_orthofinder_results.sh" ]] || {
    echo "ERROR: package runner is not executable: $PACKAGE_ROOT/run_orthofinder_results.sh" >&2
    exit 2
}
[[ -d "$RESULTS_SOURCE" ]] || {
    echo "ERROR: completed OrthoFinder result directory is unavailable: $RESULTS_SOURCE" >&2
    exit 2
}
[[ "$FOCUS_SOURCE" == "-" || -f "$FOCUS_SOURCE" ]] || {
    echo "ERROR: focus-protein authority is unavailable: $FOCUS_SOURCE" >&2
    exit 2
}
[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || {
    echo "ERROR: run ID contains unsupported characters: $RUN_ID" >&2
    exit 2
}

for argument in "$@"; do
    case "$argument" in
        --action|--action=*|--results-dir|--results-dir=*|--output-dir|--output-dir=*|\
        --run-id|--run-id=*|--work-dir|--work-dir=*|--focus-proteins|\
        --focus-proteins=*|--alignment-dir|--alignment-dir=*|--distance-source|\
        --distance-source=*|--distance-group-type|--distance-group-type=*|\
        --distance-hierarchy-node|--distance-hierarchy-node=*|--distance-max-groups|\
        --distance-max-groups=*|--parse-gene-trees)
            echo "ERROR: selection, path and authority option is managed by this wrapper: $argument" >&2
            exit 2
            ;;
    esac
done

[[ -n "${TMPDIR:-}" && "$TMPDIR" == /* && -d "$TMPDIR" && -w "$TMPDIR" ]] || {
    echo "ERROR: scheduler-managed TMPDIR must be an absolute writable directory." >&2
    exit 2
}
[[ "$PERSISTENT_OUTPUT" == /* ]] || {
    echo "ERROR: persistent output must be an absolute path." >&2
    exit 2
}
OUTPUT_PARENT=$(dirname -- "$PERSISTENT_OUTPUT")
OUTPUT_NAME=$(basename -- "$PERSISTENT_OUTPUT")
[[ "$OUTPUT_PARENT" != "/" && "$OUTPUT_NAME" != "." && "$OUTPUT_NAME" != "/" ]] || {
    echo "ERROR: persistent output target is too broad." >&2
    exit 2
}
[[ ! -e "$PERSISTENT_OUTPUT" ]] || {
    echo "ERROR: persistent output already exists: $PERSISTENT_OUTPUT" >&2
    exit 2
}

umask 077
JOB_ID=${SLURM_JOB_ID:-manual_$$}
JOB_WORK_DIR=$(mktemp -d "${TMPDIR%/}/orthofinder_e3_${USER}_${JOB_ID}_XXXXXX")
PERSISTENT_INCOMING=""
PUBLISH_LOCK=""

cleanup() {
    local expected_incoming="$OUTPUT_PARENT/.${OUTPUT_NAME}.incoming.${JOB_ID}"
    local expected_lock="$OUTPUT_PARENT/.${OUTPUT_NAME}.publish.lock"
    if [[ -n "$PERSISTENT_INCOMING" && "$PERSISTENT_INCOMING" == "$expected_incoming" ]]; then
        rm -rf -- "$PERSISTENT_INCOMING"
    fi
    if [[ -n "$PUBLISH_LOCK" && "$PUBLISH_LOCK" == "$expected_lock" ]]; then
        rmdir "$PUBLISH_LOCK" 2>/dev/null || true
    fi
    case "$JOB_WORK_DIR" in
        "${TMPDIR%/}"/orthofinder_e3_*) rm -rf -- "$JOB_WORK_DIR" ;;
    esac
}
trap cleanup EXIT

echo "Host: $(hostname)"
echo "Started UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "Slurm job: $JOB_ID"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK:-unavailable}"
echo "Node-local work directory: $JOB_WORK_DIR"

SOURCE_KB=$(du -sk -- "$RESULTS_SOURCE" | awk '{print $1}')
AVAILABLE_KB=$(df -Pk "$TMPDIR" | awk 'NR == 2 {print $4}')
RESERVE_KB=$((10 * 1024 * 1024))
REQUIRED_KB=$((SOURCE_KB + RESERVE_KB))
[[ "$AVAILABLE_KB" -ge "$REQUIRED_KB" ]] || {
    echo "ERROR: scheduler TMPDIR has ${AVAILABLE_KB} KiB free; source staging plus the required analysis reserve needs at least ${REQUIRED_KB} KiB." >&2
    exit 2
}
echo "TMPDIR capacity check: source=${SOURCE_KB} KiB, free=${AVAILABLE_KB} KiB, analysis reserve=${RESERVE_KB} KiB"

RESULTS_STAGE="$JOB_WORK_DIR/orthofinder_results_source"
mkdir "$RESULTS_STAGE"
echo "Staging completed OrthoFinder inputs to scheduler TMPDIR."
rsync -a -- "$RESULTS_SOURCE/" "$RESULTS_STAGE/"

FOCUS_STAGE="$JOB_WORK_DIR/e3_seed_catalogue.tsv"
if [[ "$FOCUS_SOURCE" == "-" ]]; then
    BUNDLED_FOCUS="$PACKAGE_ROOT/src/orthofinder_interrogation_app/data/e3_seed_catalogue.tsv"
    [[ -f "$BUNDLED_FOCUS" ]] || {
        echo "ERROR: packaged E3 seed catalogue is unavailable: $BUNDLED_FOCUS" >&2
        exit 2
    }
    cp -- "$BUNDLED_FOCUS" "$FOCUS_STAGE"
else
    cp -- "$FOCUS_SOURCE" "$FOCUS_STAGE"
fi

LOCAL_OUTPUT="$JOB_WORK_DIR/completed_resource"
PIPELINE_WORK="$JOB_WORK_DIR/pipeline_work"
"$PACKAGE_ROOT/run_orthofinder_results.sh" \
    --conda-env "$CONDA_ENV" \
    -- \
    --action e3-precursor \
    --results-dir "$RESULTS_STAGE" \
    --focus-proteins "$FOCUS_STAGE" \
    --run-id "$RUN_ID" \
    --output-dir "$LOCAL_OUTPUT" \
    --work-dir "$PIPELINE_WORK" \
    "$@"

[[ -f "$LOCAL_OUTPUT/run_manifest.json" ]] || {
    echo "ERROR: completed local resource lacks run_manifest.json." >&2
    exit 2
}
grep -Eq '"status"[[:space:]]*:[[:space:]]*"complete"' \
    "$LOCAL_OUTPUT/run_manifest.json" || {
        echo "ERROR: local run manifest is not marked complete." >&2
        exit 2
    }
for required_file in \
    tables/e3_cluster_results.tsv.gz \
    tables/e3_seed_catalogue_audit.tsv.gz \
    tables/e3_seed_matches.tsv.gz \
    tables/pairwise_distances.tsv.gz \
    duckdb/orthofinder_results.duckdb \
    qc/validation_checks.tsv; do
    [[ -s "$LOCAL_OUTPUT/$required_file" ]] || {
        echo "ERROR: completed local resource lacks required output: $required_file" >&2
        exit 2
    }
done
gzip -t \
    "$LOCAL_OUTPUT/tables/e3_cluster_results.tsv.gz" \
    "$LOCAL_OUTPUT/tables/e3_seed_catalogue_audit.tsv.gz" \
    "$LOCAL_OUTPUT/tables/e3_seed_matches.tsv.gz" \
    "$LOCAL_OUTPUT/tables/pairwise_distances.tsv.gz"

mkdir -p -- "$OUTPUT_PARENT"
PERSISTENT_INCOMING="$OUTPUT_PARENT/.${OUTPUT_NAME}.incoming.${JOB_ID}"
[[ ! -e "$PERSISTENT_INCOMING" ]] || {
    echo "ERROR: persistent incoming path already exists: $PERSISTENT_INCOMING" >&2
    exit 2
}
mkdir -m 700 -- "$PERSISTENT_INCOMING"
echo "Copying completed resource to persistent incoming path: $PERSISTENT_INCOMING"
rsync -a -- "$LOCAL_OUTPUT/" "$PERSISTENT_INCOMING/"
VERIFY_OUTPUT=$(rsync -a --checksum --dry-run --itemize-changes -- \
    "$LOCAL_OUTPUT/" "$PERSISTENT_INCOMING/")
[[ -z "$VERIFY_OUTPUT" ]] || {
    echo "ERROR: persistent copy differs from the completed TMPDIR resource." >&2
    echo "$VERIFY_OUTPUT" >&2
    exit 2
}

PUBLISH_LOCK_CANDIDATE="$OUTPUT_PARENT/.${OUTPUT_NAME}.publish.lock"
mkdir -- "$PUBLISH_LOCK_CANDIDATE" || {
    echo "ERROR: another publication holds the output lock: $PUBLISH_LOCK_CANDIDATE" >&2
    exit 2
}
PUBLISH_LOCK="$PUBLISH_LOCK_CANDIDATE"
[[ ! -e "$PERSISTENT_OUTPUT" ]] || {
    echo "ERROR: persistent output appeared during publication: $PERSISTENT_OUTPUT" >&2
    exit 2
}
mv "$PERSISTENT_INCOMING" "$PERSISTENT_OUTPUT"
PERSISTENT_INCOMING=""
rmdir "$PUBLISH_LOCK"
PUBLISH_LOCK=""

echo "Published completed E3 precursor resource: $PERSISTENT_OUTPUT"
echo "Primary cluster table: $PERSISTENT_OUTPUT/tables/e3_cluster_results.tsv.gz"
echo "Pairwise distance table: $PERSISTENT_OUTPUT/tables/pairwise_distances.tsv.gz"
sha256sum \
    "$PERSISTENT_OUTPUT/tables/e3_cluster_results.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/pairwise_distances.tsv.gz" \
    "$PERSISTENT_OUTPUT/run_manifest.json"
echo "Finished UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"

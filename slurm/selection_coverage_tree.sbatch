#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    echo "Usage: $0 PACKAGE_ROOT CONDA_ENV RESOURCE_DIR TAXONOMY_MAP EXPECTED_OR_DASH FOCUS_OR_DASH PERSISTENT_OUTPUT [coverage options]" >&2
}

[[ $# -ge 7 ]] || { usage; exit 2; }
PACKAGE_ROOT=$1
CONDA_ENV=$2
RESOURCE_SOURCE=$3
TAXONOMY_SOURCE=$4
EXPECTED_SOURCE=$5
FOCUS_SOURCE=$6
PERSISTENT_OUTPUT=$7
shift 7

[[ -x "$PACKAGE_ROOT/run_orthofinder_results.sh" ]] || {
    echo "ERROR: package runner is not executable: $PACKAGE_ROOT/run_orthofinder_results.sh" >&2
    exit 2
}
[[ -d "$RESOURCE_SOURCE" ]] || {
    echo "ERROR: completed resource directory is unavailable: $RESOURCE_SOURCE" >&2
    exit 2
}
[[ -f "$TAXONOMY_SOURCE" ]] || {
    echo "ERROR: reviewed taxonomy mapping is unavailable: $TAXONOMY_SOURCE" >&2
    exit 2
}
[[ "$EXPECTED_SOURCE" == "-" || -f "$EXPECTED_SOURCE" ]] || {
    echo "ERROR: expected-taxon authority is unavailable: $EXPECTED_SOURCE" >&2
    exit 2
}
[[ "$FOCUS_SOURCE" == "-" || -f "$FOCUS_SOURCE" ]] || {
    echo "ERROR: focus-protein authority is unavailable: $FOCUS_SOURCE" >&2
    exit 2
}

for argument in "$@"; do
    case "$argument" in
        --action|--action=*|--resource-dir|--resource-dir=*|--taxonomy-map|--taxonomy-map=*|\
        --expected-taxa|--expected-taxa=*|--focus-proteins|--focus-proteins=*|\
        --output-dir|--output-dir=*)
            echo "ERROR: path/action option is managed by this Slurm wrapper: $argument" >&2
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
JOB_WORK_DIR=$(mktemp -d "${TMPDIR%/}/orthofinder_selection_${USER}_${JOB_ID}_XXXXXX")
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
        "${TMPDIR%/}"/orthofinder_selection_*) rm -rf -- "$JOB_WORK_DIR" ;;
    esac
}
trap cleanup EXIT

echo "Host: $(hostname)"
echo "Started UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "Slurm job: $JOB_ID"
echo "Node-local work directory: $JOB_WORK_DIR"

RESOURCE_STAGE="$JOB_WORK_DIR/resource"
mkdir "$RESOURCE_STAGE"
echo "Staging immutable resource to scheduler TMPDIR."
rsync -a -- "$RESOURCE_SOURCE/" "$RESOURCE_STAGE/"

TAXONOMY_STAGE="$JOB_WORK_DIR/taxonomy_mapping.tsv"
cp -- "$TAXONOMY_SOURCE" "$TAXONOMY_STAGE"
COMMAND=(
    --action coverage-tree
    --resource-dir "$RESOURCE_STAGE"
    --taxonomy-map "$TAXONOMY_STAGE"
    --output-dir "$JOB_WORK_DIR/selection_coverage_output"
)
if [[ "$EXPECTED_SOURCE" != "-" ]]; then
    EXPECTED_STAGE="$JOB_WORK_DIR/expected_taxa.tsv"
    cp -- "$EXPECTED_SOURCE" "$EXPECTED_STAGE"
    COMMAND+=(--expected-taxa "$EXPECTED_STAGE")
fi
if [[ "$FOCUS_SOURCE" != "-" ]]; then
    FOCUS_STAGE="$JOB_WORK_DIR/$(basename -- "$FOCUS_SOURCE")"
    cp -- "$FOCUS_SOURCE" "$FOCUS_STAGE"
    COMMAND+=(--focus-proteins "$FOCUS_STAGE")
fi
COMMAND+=("$@")

"$PACKAGE_ROOT/run_orthofinder_results.sh" \
    --conda-env "$CONDA_ENV" \
    -- \
    "${COMMAND[@]}"

LOCAL_OUTPUT="$JOB_WORK_DIR/selection_coverage_output"
[[ -f "$LOCAL_OUTPUT/checksums.sha256" ]] || {
    echo "ERROR: completed local output lacks checksums.sha256." >&2
    exit 2
}
(cd "$LOCAL_OUTPUT" && sha256sum --check checksums.sha256)

mkdir -p -- "$OUTPUT_PARENT"
PERSISTENT_INCOMING="$OUTPUT_PARENT/.${OUTPUT_NAME}.incoming.${JOB_ID}"
[[ ! -e "$PERSISTENT_INCOMING" ]] || {
    echo "ERROR: persistent incoming path already exists: $PERSISTENT_INCOMING" >&2
    exit 2
}
mkdir -m 700 -- "$PERSISTENT_INCOMING"
echo "Copying verified result to persistent incoming path: $PERSISTENT_INCOMING"
rsync -a -- "$LOCAL_OUTPUT/" "$PERSISTENT_INCOMING/"
(cd "$PERSISTENT_INCOMING" && sha256sum --check checksums.sha256)
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

echo "Published completed selection coverage package: $PERSISTENT_OUTPUT"
echo "Finished UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"

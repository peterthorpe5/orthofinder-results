#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    echo "Usage: $0 PACKAGE_ROOT CONDA_ENV RESULTS_DIR FOCUS_OR_DASH BENCHMARK_OR_DASH EXPRESSION_MANIFEST_OR_DASH RUN_ID PERSISTENT_OUTPUT [safe analysis options]" >&2
}

[[ $# -ge 8 ]] || { usage; exit 2; }
PACKAGE_ROOT=$1
CONDA_ENV=$2
RESULTS_SOURCE=$3
FOCUS_SOURCE=$4
BENCHMARK_SOURCE=$5
EXPRESSION_MANIFEST_SOURCE=$6
RUN_ID=$7
PERSISTENT_OUTPUT=$8
shift 8

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
[[ "$BENCHMARK_SOURCE" == "-" || -f "$BENCHMARK_SOURCE" ]] || {
    echo "ERROR: dispersion-marker authority is unavailable: $BENCHMARK_SOURCE" >&2
    exit 2
}
[[ "$EXPRESSION_MANIFEST_SOURCE" == "-" || -f "$EXPRESSION_MANIFEST_SOURCE" ]] || {
    echo "ERROR: expression-resource manifest is unavailable: $EXPRESSION_MANIFEST_SOURCE" >&2
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
        --focus-proteins=*|--benchmark-proteins|--benchmark-proteins=*|\
        --expression-manifest|--expression-manifest=*|--expression-aliases|\
        --expression-aliases=*|--expression-threads|--expression-threads=*|\
        --expression-memory-mb|--expression-memory-mb=*|--include-protein-sequences|\
        --alignment-dir|--alignment-dir=*|--distance-source|--distance-source=*|\
        --distance-group-type|--distance-group-type=*|--distance-hierarchy-node|\
        --distance-hierarchy-node=*|--distance-max-groups|--distance-max-groups=*|\
        --parse-gene-trees)
            echo "ERROR: path, selection and authority option is managed by this wrapper: $argument" >&2
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
JOB_WORK_DIR=$(mktemp -d "${TMPDIR%/}/orthofinder_benchmark_${USER}_${JOB_ID}_XXXXXX")
PUBLISH_LOCK=""

cleanup() {
    local expected_lock="$OUTPUT_PARENT/.${OUTPUT_NAME}.publish.lock"
    if [[ -n "$PUBLISH_LOCK" && "$PUBLISH_LOCK" == "$expected_lock" ]]; then
        rmdir "$PUBLISH_LOCK" 2>/dev/null || true
    fi
    case "$JOB_WORK_DIR" in
        "${TMPDIR%/}"/orthofinder_benchmark_*) rm -rf -- "$JOB_WORK_DIR" ;;
    esac
}
trap cleanup EXIT

checksum_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum -- "$1"
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 -- "$1"
    else
        echo "ERROR: sha256sum or shasum is required for final checksums." >&2
        return 2
    fi
}

echo "Host: $(hostname)"
echo "Started UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "Slurm job: $JOB_ID"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK:-unavailable}"
echo "Node-local work directory: $JOB_WORK_DIR"

SOURCE_KB=$(du -sk -- "$RESULTS_SOURCE" | awk '{print $1}')
AVAILABLE_KB=$(df -Pk "$TMPDIR" | awk 'NR == 2 {print $4}')
SCRATCH_RESERVE_GIB=${ORTHOFINDER_SCRATCH_RESERVE_GIB:-200}
[[ "$SCRATCH_RESERVE_GIB" =~ ^[0-9]+$ ]] && \
    [[ "$SCRATCH_RESERVE_GIB" -ge 1 ]] || {
    echo "ERROR: ORTHOFINDER_SCRATCH_RESERVE_GIB must be a positive integer." >&2
    exit 2
}
RESERVE_KB=$((SCRATCH_RESERVE_GIB * 1024 * 1024))
REQUIRED_KB=$((SOURCE_KB + RESERVE_KB))
[[ "$AVAILABLE_KB" -ge "$REQUIRED_KB" ]] || {
    echo "ERROR: scheduler TMPDIR has ${AVAILABLE_KB} KiB free; source staging plus the required analysis reserve needs at least ${REQUIRED_KB} KiB." >&2
    exit 2
}
echo "TMPDIR capacity check: source=${SOURCE_KB} KiB, free=${AVAILABLE_KB} KiB, analysis reserve=${RESERVE_KB} KiB (${SCRATCH_RESERVE_GIB} GiB)"

mkdir -p -- "$OUTPUT_PARENT"
PUBLISH_LOCK_CANDIDATE="$OUTPUT_PARENT/.${OUTPUT_NAME}.publish.lock"
mkdir -- "$PUBLISH_LOCK_CANDIDATE" || {
    echo "ERROR: another publication holds the output lock: $PUBLISH_LOCK_CANDIDATE" >&2
    exit 2
}
PUBLISH_LOCK="$PUBLISH_LOCK_CANDIDATE"
[[ ! -e "$PERSISTENT_OUTPUT" ]] || {
    echo "ERROR: persistent output appeared during setup: $PERSISTENT_OUTPUT" >&2
    exit 2
}

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
    cp "$BUNDLED_FOCUS" "$FOCUS_STAGE"
else
    cp "$FOCUS_SOURCE" "$FOCUS_STAGE"
fi

BENCHMARK_STAGE="$JOB_WORK_DIR/arabidopsis_dispersion_benchmarks.tsv"
if [[ "$BENCHMARK_SOURCE" == "-" ]]; then
    BUNDLED_BENCHMARK="$PACKAGE_ROOT/src/orthofinder_interrogation_app/data/arabidopsis_dispersion_benchmarks.tsv"
    [[ -f "$BUNDLED_BENCHMARK" ]] || {
        echo "ERROR: packaged dispersion-marker authority is unavailable: $BUNDLED_BENCHMARK" >&2
        exit 2
    }
    cp "$BUNDLED_BENCHMARK" "$BENCHMARK_STAGE"
else
    cp "$BENCHMARK_SOURCE" "$BENCHMARK_STAGE"
fi

EXPRESSION_ARGUMENTS=(--include-protein-sequences)
if [[ "$EXPRESSION_MANIFEST_SOURCE" != "-" ]]; then
    BUNDLED_EXPRESSION_SCOPE="$PACKAGE_ROOT/src/orthofinder_interrogation_app/data/expression_species_scope.tsv"
    EXPRESSION_SCOPE_SOURCE=${ORTHOFINDER_EXPRESSION_SPECIES_SCOPE:-$BUNDLED_EXPRESSION_SCOPE}
    [[ -s "$EXPRESSION_SCOPE_SOURCE" ]] || {
        echo "ERROR: expression species scope is unavailable: $EXPRESSION_SCOPE_SOURCE" >&2
        exit 2
    }
    EXPRESSION_MANIFEST_STAGE="$JOB_WORK_DIR/expression_resources.tsv"
    awk -F '\t' -v OFS='\t' '
        NR == FNR {
            if (FNR == 1) {
                for (column = 1; column <= NF; column++) {
                    if ($column == "species_label") scope_column = column
                }
                if (!scope_column) exit 41
                next
            }
            if ($scope_column != "") allowed[$scope_column] = 1
            next
        }
        FNR == 1 {
            for (column = 1; column <= NF; column++) {
                if ($column == "species_column") species_column = column
            }
            if (!species_column) exit 42
            print
            next
        }
        $species_column in allowed {
            print
            selected++
        }
        END {
            if (!selected) exit 43
        }
    ' "$EXPRESSION_SCOPE_SOURCE" "$EXPRESSION_MANIFEST_SOURCE" \
        > "$EXPRESSION_MANIFEST_STAGE" || {
        echo "ERROR: could not restrict the expression manifest to the reviewed focal-species scope." >&2
        exit 2
    }
    EXPRESSION_RESOURCE_COUNT=$(( $(wc -l < "$EXPRESSION_MANIFEST_STAGE") - 1 ))
    echo "RNA-seq scope: 12 target plant species plus Homo_sapiens; manifest resources=${EXPRESSION_RESOURCE_COUNT}"
    EXPRESSION_MEMORY_MB=${ORTHOFINDER_EXPRESSION_MEMORY_MB:-49152}
    [[ "$EXPRESSION_MEMORY_MB" =~ ^[0-9]+$ ]] && \
        [[ "$EXPRESSION_MEMORY_MB" -ge 256 ]] || {
        echo "ERROR: ORTHOFINDER_EXPRESSION_MEMORY_MB must be at least 256 MiB." >&2
        exit 2
    }
    EXPRESSION_THREADS=${ORTHOFINDER_EXPRESSION_THREADS:-${SLURM_CPUS_PER_TASK:-4}}
    [[ "$EXPRESSION_THREADS" =~ ^[0-9]+$ ]] && \
        [[ "$EXPRESSION_THREADS" -ge 1 ]] || {
        echo "ERROR: ORTHOFINDER_EXPRESSION_THREADS must be a positive integer." >&2
        exit 2
    }
    if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]] && \
        [[ "$EXPRESSION_THREADS" -gt "$SLURM_CPUS_PER_TASK" ]]; then
        echo "ERROR: expression threads exceed the allocated Slurm CPUs." >&2
        exit 2
    fi
    echo "RNA-seq DuckDB controls: threads=${EXPRESSION_THREADS}, memory_limit_mib=${EXPRESSION_MEMORY_MB}"
    EXPRESSION_ARGUMENTS+=(
        --expression-manifest "$EXPRESSION_MANIFEST_STAGE"
        --expression-threads "$EXPRESSION_THREADS"
        --expression-memory-mb "$EXPRESSION_MEMORY_MB"
    )
fi

PIPELINE_WORK="$JOB_WORK_DIR/pipeline_work"
"$PACKAGE_ROOT/run_orthofinder_results.sh" \
    --conda-env "$CONDA_ENV" \
    -- \
    --action dispersion-benchmark \
    --results-dir "$RESULTS_STAGE" \
    --focus-proteins "$FOCUS_STAGE" \
    --benchmark-proteins "$BENCHMARK_STAGE" \
    --run-id "$RUN_ID" \
    --output-dir "$PERSISTENT_OUTPUT" \
    --work-dir "$PIPELINE_WORK" \
    "${EXPRESSION_ARGUMENTS[@]}" \
    "$@"

[[ -f "$PERSISTENT_OUTPUT/run_manifest.json" ]] || {
    echo "ERROR: published resource lacks run_manifest.json." >&2
    exit 2
}
grep -Eq '"status"[[:space:]]*:[[:space:]]*"complete"' \
    "$PERSISTENT_OUTPUT/run_manifest.json" || {
        echo "ERROR: published run manifest is not marked complete." >&2
        exit 2
    }
for required_file in \
    tables/e3_cluster_results.tsv.gz \
    tables/benchmark_marker_audit.tsv.gz \
    tables/benchmark_group_profiles.tsv.gz \
    tables/benchmark_matched_controls.tsv.gz \
    tables/benchmark_cluster_results.tsv.gz \
    tables/benchmark_background_statistics.tsv.gz \
    tables/benchmark_contrasts.tsv.gz \
    tables/benchmark_individual_comparisons.tsv.gz \
    tables/benchmark_cluster_classifications.tsv.gz \
    tables/pairwise_distances.tsv.gz \
    duckdb/orthofinder_results.duckdb \
    qc/validation_checks.tsv; do
    [[ -s "$PERSISTENT_OUTPUT/$required_file" ]] || {
        echo "ERROR: published resource lacks required output: $required_file" >&2
        exit 2
    }
done
if [[ "$EXPRESSION_MANIFEST_SOURCE" != "-" ]]; then
    for required_file in \
        tables/expression_identifier_aliases.tsv.gz \
        tables/expression_member_mapping.tsv.gz \
        tables/expression_member_summary.tsv.gz \
        tables/expression_group_summary.tsv.gz \
        tables/expression_import_audit.tsv.gz; do
        [[ -s "$PERSISTENT_OUTPUT/$required_file" ]] || {
            echo "ERROR: published resource lacks required RNA-seq output: $required_file" >&2
            exit 2
        }
    done
fi
gzip -t \
    "$PERSISTENT_OUTPUT/tables/e3_cluster_results.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/benchmark_marker_audit.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/benchmark_group_profiles.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/benchmark_matched_controls.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/benchmark_cluster_results.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/benchmark_background_statistics.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/benchmark_contrasts.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/benchmark_individual_comparisons.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/benchmark_cluster_classifications.tsv.gz" \
    "$PERSISTENT_OUTPUT/tables/pairwise_distances.tsv.gz"
if [[ "$EXPRESSION_MANIFEST_SOURCE" != "-" ]]; then
    gzip -t \
        "$PERSISTENT_OUTPUT/tables/expression_identifier_aliases.tsv.gz" \
        "$PERSISTENT_OUTPUT/tables/expression_member_mapping.tsv.gz" \
        "$PERSISTENT_OUTPUT/tables/expression_member_summary.tsv.gz" \
        "$PERSISTENT_OUTPUT/tables/expression_group_summary.tsv.gz" \
        "$PERSISTENT_OUTPUT/tables/expression_import_audit.tsv.gz"
fi

rmdir "$PUBLISH_LOCK"
PUBLISH_LOCK=""

echo "Validated completed matched-background dispersion resource: $PERSISTENT_OUTPUT"
echo "Cluster statistics: $PERSISTENT_OUTPUT/tables/benchmark_cluster_results.tsv.gz"
echo "Profile contrasts: $PERSISTENT_OUTPUT/tables/benchmark_contrasts.tsv.gz"
echo "Individual comparisons: $PERSISTENT_OUTPUT/tables/benchmark_individual_comparisons.tsv.gz"
echo "Protein sequence evidence: DuckDB relation protein_sequences"
if [[ "$EXPRESSION_MANIFEST_SOURCE" != "-" ]]; then
    echo "RNA-seq contexts: DuckDB relation expression_context"
    echo "RNA-seq group summary: $PERSISTENT_OUTPUT/tables/expression_group_summary.tsv.gz"
fi
checksum_file "$PERSISTENT_OUTPUT/tables/benchmark_cluster_results.tsv.gz"
checksum_file "$PERSISTENT_OUTPUT/tables/benchmark_contrasts.tsv.gz"
checksum_file "$PERSISTENT_OUTPUT/run_manifest.json"
echo "Finished UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"

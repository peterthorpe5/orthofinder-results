#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=""
PYTHON_EXECUTABLE="python"
PACKAGE_ARGS=()

while (($#)); do
    case "$1" in
        --conda-env)
            [[ $# -ge 2 ]] || { echo "ERROR: --conda-env requires a value." >&2; exit 2; }
            CONDA_ENV=$2
            shift 2
            ;;
        --python-executable)
            [[ $# -ge 2 ]] || { echo "ERROR: --python-executable requires a value." >&2; exit 2; }
            PYTHON_EXECUTABLE=$2
            shift 2
            ;;
        --)
            shift
            PACKAGE_ARGS+=("$@")
            break
            ;;
        *)
            PACKAGE_ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ -n "$CONDA_ENV" ]]; then
    command -v conda >/dev/null 2>&1 || {
        echo "ERROR: conda is required for --conda-env ${CONDA_ENV}." >&2
        exit 2
    }
    CONDA_PREFIX_PATH=$(conda run --name "$CONDA_ENV" env | awk -F= '
        $1 == "CONDA_PREFIX" { print substr($0, index($0, "=") + 1); exit }
    ')
    [[ -n "$CONDA_PREFIX_PATH" && -d "$CONDA_PREFIX_PATH/lib" ]] || {
        echo "ERROR: could not resolve the library directory for Conda environment ${CONDA_ENV}." >&2
        exit 2
    }
    exec conda run --no-capture-output --name "$CONDA_ENV" \
        env \
        LD_LIBRARY_PATH="$CONDA_PREFIX_PATH/lib" \
        MALLOC_ARENA_MAX=2 \
        "$PYTHON_EXECUTABLE" -m orthofinder_results "${PACKAGE_ARGS[@]}"
fi

exec env MALLOC_ARENA_MAX=2 \
    "$PYTHON_EXECUTABLE" -m orthofinder_results "${PACKAGE_ARGS[@]}"

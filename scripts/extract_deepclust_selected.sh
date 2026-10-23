#!/usr/bin/env bash
# Run the upstream extractor beside an installed, extracted DeepClust release.
set -euo pipefail

retrieval_repo=""
centroid_file=""
database=""
index=""
output=""
threads="1"
python_executable="python3"

usage() {
    cat <<'TEXT'
Usage: extract_deepclust_selected.sh --retrieval-repo PATH --centroid-file PATH
       --database PATH --index PATH --output PATH [--threads N] [--python PATH]

The database and index must already exist on the extraction machine.
Output is one FASTA per exact centroid, with original sequence identifiers.
TEXT
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help) usage; exit 0 ;;
        --retrieval-repo|--centroid-file|--database|--index|--output|--threads|--python)
            [[ $# -ge 2 && -n "$2" ]] || { usage >&2; exit 2; }
            case "$1" in
                --retrieval-repo) retrieval_repo="$2" ;;
                --centroid-file) centroid_file="$2" ;;
                --database) database="$2" ;;
                --index) index="$2" ;;
                --output) output="$2" ;;
                --threads) threads="$2" ;;
                --python) python_executable="$2" ;;
            esac
            shift 2 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -f "$retrieval_repo/cluster_index.py" && -s "$centroid_file" && -f "$database" \
    && -e "$index" && -n "$output" ]] || { usage >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ && ${#threads} -le 3 && "$threads" -le 256 ]] || {
    printf 'Threads must be an integer from 1 to 256.\n' >&2; exit 2;
}
case "$(uname -s)" in
    Darwin|MINGW*|MSYS*)
        [[ "$threads" == "1" ]] || {
            printf 'The upstream extractor requires --threads 1 on this platform.\n' >&2
            exit 2
        } ;;
esac
while IFS= read -r centroid || [[ -n "$centroid" ]]; do
    [[ -n "$centroid" && ! "$centroid" =~ [[:space:][:cntrl:]] ]] || {
        printf 'Centroid files require one exact non-empty ID per line, without a header.\n' >&2
        exit 2
    }
    [[ "$centroid" =~ ^[A-Za-z0-9_.:+\|=-]+$ && "$centroid" != "." && "$centroid" != ".." ]] || {
        printf 'Unsupported centroid characters for the upstream SQL and filename interface.\n' >&2
        exit 2
    }
done < "$centroid_file"
if [[ -e "$output" && ( ! -d "$output" || -n "$(ls -A "$output")" ) ]]; then
    printf 'Use a new or empty output directory; existing exports will not be overwritten.\n' >&2
    exit 2
fi
mkdir -p "$output"
exec "$python_executable" "$retrieval_repo/cluster_index.py" \
    --centroid_file "$centroid_file" --per-clust-output 1 --threads "$threads" \
    --max_num_of_cluster_at_once 1 --verbose 1 \
    -path_to_DCD "$database" -path_to_output "${output%/}/" -path_to_index "$index"

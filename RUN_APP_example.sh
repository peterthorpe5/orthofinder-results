ENV_PREFIX="/Users/PThorpe001/miniforge3/envs/orthofinder_results"
RESOURCE_DIR="$HOME/Downloads/results_feb26_dispersion_benchmark_v0_9_0"
APP_LOG_DIR="$HOME/orthofinder_results_app_logs"

mkdir -p "$APP_LOG_DIR"

"$ENV_PREFIX/bin/orthofinder-interrogation-app" \
    --resource-dir "$RESOURCE_DIR" \
    --server-port 8502 \
    --log-file "$APP_LOG_DIR/orthofinder_app_$(date -u +%Y%m%dT%H%M%SZ).log"

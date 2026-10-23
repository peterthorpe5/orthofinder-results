#!/usr/bin/env bash
# Launch the OrthoFinder app on Sabrina and open it through SSH on this Mac.
# Usage: bash launch_orthofinder_denbi.sh --help
# Bash 3.2+; no Python is embedded in this script. Version 1.0.0.

set -Eeuo pipefail
trap 'printf "[OrthoFinder] Unexpected failure at line %s (exit %s).\n" "$LINENO" "$?" >&2' ERR

OF_MODE="start"
OF_OPEN_BROWSER="yes"
OF_SSH_HOST="${ORTHOFINDER_DENBI_HOST:-denbi}"
OF_LOCAL_PORT="8510"
OF_REMOTE_PORT="8506"
OF_RESTART_APP="no"
OF_ENV="${ORTHOFINDER_DENBI_ENV:-/home/ubuntu/data/pete/software/miniforge3/envs/orthofinder_results}"
OF_CONFIG="${ORTHOFINDER_DENBI_CONFIG:-/home/ubuntu/data/orthofinder_research/config/research.toml}"
OF_RESOURCES="${ORTHOFINDER_DENBI_RESOURCE:-/home/ubuntu/data/orthofinder_research/resources/current}"
OF_STATE="${ORTHOFINDER_DENBI_STATE:-/home/ubuntu/data/orthofinder_research}"
OF_SEQUENCE="${ORTHOFINDER_DENBI_SEQUENCES:-}"
OF_TAXONOMY="${ORTHOFINDER_DENBI_TAXONOMY:-}"

while (($#)); do
    case "$1" in
        --help|-h)
            cat <<'USAGE'
Run on your Mac:
  bash launch_orthofinder_denbi.sh                 Start/reuse OrthoFinder and open the browser.
  bash launch_orthofinder_denbi.sh --no-browser    Start/reuse OrthoFinder and print its URL.
  bash launch_orthofinder_denbi.sh --status        Inspect services, jobs and the tunnel.
  bash launch_orthofinder_denbi.sh --stop-tunnel   Disconnect this launcher's tunnel.
  bash launch_orthofinder_denbi.sh --local-port 8510
                                        Choose the Mac port (default 8510).
  bash launch_orthofinder_denbi.sh --remote-port 8506
                                        Choose the Sabrina port (default 8506).
  bash launch_orthofinder_denbi.sh --restart-app   Restart only the owned OrthoFinder GUI.
  bash launch_orthofinder_denbi.sh --resource-dir /absolute/server/resource
  bash launch_orthofinder_denbi.sh --research-config /absolute/server/research.toml

Other named options: --ssh-host, --env-prefix, --state-dir,
--sequence-parquet, --taxonomy-map. Paths refer to Sabrina, not your Mac.
Defaults use a separate orthofinder_results environment and orthofinder_research
workspace under /home/ubuntu/data. See docs/denbi.md for first-time installation.
The server reads its installed DeepClust database and runs the analysis worker.
Only the browser interface crosses the SSH tunnel; database files stay on Sabrina.

The OrthoFinder app and analysis worker run on Sabrina through your existing SSH alias
"denbi". Their tmux session survives closing the Mac terminal. Stopping the
tunnel leaves remote analyses running. Run the launcher again to reconnect.
USAGE
            exit 0
            ;;
        --version)
            printf 'OrthoFinder Sabrina launcher 1.0.0\n'
            exit 0
            ;;
        --no-browser) OF_OPEN_BROWSER="no" ;;
        --restart-app) OF_RESTART_APP="yes" ;;
        --ssh-host|--env-prefix|--research-config|--resource-dir|--state-dir|--sequence-parquet|--taxonomy-map)
            if (($# < 2)) || [[ -z "$2" || "$2" == *$'\n'* || "$2" == *$'\r'* ]]; then
                printf '[OrthoFinder] %s requires a non-empty value without newlines.\n' "$1" >&2
                exit 2
            fi
            if [[ "$1" != "--ssh-host" && "$2" != /* ]]; then
                printf '[OrthoFinder] %s requires an absolute server path.\n' "$1" >&2
                exit 2
            fi
            case "$1" in
                --ssh-host)
                    if [[ ! "$2" =~ ^[A-Za-z0-9_][A-Za-z0-9_.@-]*$ ]]; then
                        printf '[OrthoFinder] Invalid SSH alias or user@host.\n' >&2
                        exit 2
                    fi
                    OF_SSH_HOST="$2" ;;
                --env-prefix) OF_ENV="$2" ;;
                --research-config) OF_CONFIG="$2" ;;
                --resource-dir) OF_RESOURCES="$2" ;;
                --state-dir) OF_STATE="$2" ;;
                --sequence-parquet) OF_SEQUENCE="$2" ;;
                --taxonomy-map) OF_TAXONOMY="$2" ;;
            esac
            shift
            ;;
        --local-port|--remote-port)
            if (($# < 2)) || [[ ! "$2" =~ ^[1-9][0-9]{0,4}$ ]] \
                || (( $2 < 1024 || $2 > 65535 )); then
                printf '[OrthoFinder] %s requires a port between 1024 and 65535.\n' "$1" >&2
                exit 2
            fi
            if [[ "$1" == "--local-port" ]]; then
                OF_LOCAL_PORT="$2"
            else
                OF_REMOTE_PORT="$2"
            fi
            shift
            ;;
        --status|--stop-tunnel)
            if [[ "$OF_MODE" != "start" ]]; then
                printf '[OrthoFinder] Choose one action: --status or --stop-tunnel.\n' >&2
                exit 2
            fi
            OF_MODE="${1#--}"
            ;;
        *)
            printf '[OrthoFinder] Unknown option: %s. Use --help.\n' "$1" >&2
            exit 2
            ;;
    esac
    shift
done

if [[ "$OF_MODE" != "start" && "$OF_RESTART_APP" == "yes" ]]; then
    printf '[OrthoFinder] --restart-app applies only to a normal launch.\n' >&2
    exit 2
fi
for OF_PATH in "$OF_ENV" "$OF_CONFIG" "$OF_RESOURCES" "$OF_STATE"; do
    if [[ "$OF_PATH" != /* || "$OF_PATH" == *$'\n'* || "$OF_PATH" == *$'\r'* ]]; then
        printf '[OrthoFinder] Remote paths must be absolute and contain no newlines.\n' >&2
        exit 2
    fi
done
if [[ ! "$OF_SSH_HOST" =~ ^[A-Za-z0-9_][A-Za-z0-9_.@-]*$ ]]; then
    printf '[OrthoFinder] Invalid SSH alias or user@host.\n' >&2
    exit 2
fi
OF_URL="http://127.0.0.1:$OF_LOCAL_PORT"
# Short, host-specific socket names also fit OpenSSH's Unix socket limit on macOS.
OF_SOCKET_KEY="$(printf '%s' "$OF_SSH_HOST:$OF_LOCAL_PORT:$OF_REMOTE_PORT" | cksum | cut -d ' ' -f 1)"
OF_CONTROL_PATH="$HOME/.ssh/orthofinder_${OF_SOCKET_KEY}.sock"

if ! command -v ssh >/dev/null 2>&1; then
    printf '[OrthoFinder] OpenSSH (ssh) is required on this computer.\n' >&2
    exit 1
fi

if [[ "$OF_MODE" == "stop-tunnel" ]]; then
    if ssh -F /dev/null -S "$OF_CONTROL_PATH" -O check "$OF_SSH_HOST" >/dev/null 2>&1; then
        ssh -F /dev/null -S "$OF_CONTROL_PATH" -O exit "$OF_SSH_HOST"
        printf '[OrthoFinder] Tunnel disconnected. Remote OrthoFinder and analysis jobs remain running.\n'
    else
        printf '[OrthoFinder] This launcher has no active tunnel.\n'
    fi
    exit 0
fi

if ! command -v curl >/dev/null 2>&1; then
    printf '[OrthoFinder] curl is required for the app readiness check.\n' >&2
    exit 1
fi

of_remote_quote() {
    # Quote one argument for the remote shell without evaluating its contents.
    printf "'%s'" "${1//\'/\'\\\'\'}"
}
OF_REMOTE_COMMAND="bash -s --"
for OF_VALUE in "$OF_MODE" "$OF_REMOTE_PORT" "$OF_RESTART_APP" "$OF_ENV" "$OF_CONFIG" \
    "$OF_RESOURCES" "$OF_STATE" "$OF_SEQUENCE" "$OF_TAXONOMY"; do
    OF_REMOTE_COMMAND+=" $(of_remote_quote "$OF_VALUE")"
done

printf '[OrthoFinder] Launcher v1.0.0: connecting to Sabrina using SSH alias %s...\n' "$OF_SSH_HOST"
if ! ssh -T -o ClearAllForwardings=yes -o ConnectTimeout=15 \
    -o ServerAliveInterval=15 -o ServerAliveCountMax=3 \
    "$OF_SSH_HOST" "$OF_REMOTE_COMMAND" <<'OF_REMOTE'
set -Eeuo pipefail
trap 'printf "[OrthoFinder remote] Failure at line %s (exit %s).\n" "$LINENO" "$?" >&2' ERR

OF_ACTION="${1:-start}"
OF_PORT="${2:-8506}"
OF_RESTART="${3:-no}"
OF_ENV="${4:?Missing environment prefix}"
OF_CONFIG="${5:?Missing research profile}"
OF_RESOURCES="${6:?Missing completed resource}"
OF_STATE="${7:?Missing state directory}"
OF_SEQUENCE="${8:-}"
OF_TAXONOMY="${9:-}"
OF_SESSION="orthofinder-research"
OF_APP_WINDOW="$OF_SESSION:=app"
OF_WORKER_WINDOW="$OF_SESSION:=worker"
OF_TMUX="$OF_ENV/bin/tmux"
OF_WORKBENCH="$OF_ENV/bin/orthofinder-protein-workbench"
OF_APP="$OF_ENV/bin/orthofinder-interrogation-app"
OF_HEALTH="http://127.0.0.1:$OF_PORT/_stcore/health"
OF_LOG_ROOT="$OF_STATE/logs"
OF_APP_LOG="$OF_LOG_ROOT/orthofinder_server.log"
OF_WORKER_LOG="$OF_LOG_ROOT/research_worker.log"

of_window_dead() {
    # Read a service window's state using strict tmux target resolution.
    #
    # Args:
    #     $1: Exact service window target.
    #
    # Returns:
    #     Prints 0 (live), 1 (dead) or absent. Fails for multiple panes.
    local of_state
    if ! of_state="$("$OF_TMUX" list-panes -t "$1" -F '#{pane_dead}' 2>/dev/null)"; then
        printf 'absent\n'
        return 0
    fi
    case "$of_state" in
        0|1) printf '%s\n' "$of_state" ;;
        *)
            printf '[OrthoFinder remote] Expected one service pane in %s; inspect that window.\n' "$1" >&2
            return 1
            ;;
    esac
}

of_diagnostics() {
    # Show recent console output and the retained service log.
    #
    # Args:
    #     $1: Exact service window target.
    #     $2: Retained log file.
    #
    # Returns:
    #     Zero after printing available diagnostics to standard error.
    "$OF_TMUX" capture-pane -p -S -80 -t "$1" >&2 || true
    if [[ -f "$2" ]]; then
        printf '\n[OrthoFinder remote] Recent output from %s:\n' "$2" >&2
        tail -n 80 "$2" >&2
    fi
    return 0
}

of_wait_for_app() {
    # Wait for HTTP readiness while checking that the OF pane is still alive.
    #
    # Returns:
    #     Zero if ready; one on process exit or after the bounded wait.
    local of_attempt of_dead
    for ((of_attempt = 0; of_attempt < 30; of_attempt++)); do
        of_dead="$(of_window_dead "$OF_APP_WINDOW")" || return 1
        if [[ "$of_dead" != "0" ]]; then
            return 1
        fi
        if curl --fail --silent --max-time 1 "$OF_HEALTH" >/dev/null 2>&1; then
            return 0
        fi
        sleep 1
    done
    return 1
}

for OF_BINARY in "$OF_TMUX" "$OF_WORKBENCH" "$OF_APP"; do
    if [[ ! -x "$OF_BINARY" ]]; then
        printf '[OrthoFinder remote] Required executable is missing: %s\n' "$OF_BINARY" >&2
        exit 1
    fi
done
if ! command -v curl >/dev/null 2>&1; then
    printf '[OrthoFinder remote] curl is required on Sabrina.\n' >&2
    exit 1
fi

if [[ "$OF_ACTION" == "status" ]]; then
    if ! "$OF_TMUX" list-windows -t "$OF_SESSION"; then
        printf '[OrthoFinder remote] The OrthoFinder tmux session is not running.\n' >&2
        exit 1
    fi
    "$OF_WORKBENCH" --config "$OF_CONFIG" jobs
    OF_APP_DEAD="$(of_window_dead "$OF_APP_WINDOW")"
    OF_WORKER_DEAD="$(of_window_dead "$OF_WORKER_WINDOW")"
    if [[ "$OF_APP_DEAD" != "0" ]]; then
        printf '[OrthoFinder remote] The OrthoFinder service window is not running.\n' >&2
        of_diagnostics "$OF_APP_WINDOW" "$OF_APP_LOG"
        exit 1
    fi
    if ! curl --fail --silent --show-error --max-time 5 "$OF_HEALTH"; then
        printf '\n[OrthoFinder remote] OrthoFinder app is not responding on port %s.\n' "$OF_PORT" >&2
        of_diagnostics "$OF_APP_WINDOW" "$OF_APP_LOG"
        exit 1
    fi
    printf '\n'
    if [[ "$OF_WORKER_DEAD" != "0" ]]; then
        printf '[OrthoFinder remote] The research worker is not running.\n' >&2
        exit 1
    fi
    exit 0
fi

for OF_REQUIRED in "$OF_CONFIG" "$OF_RESOURCES"; do
    if [[ ! -r "$OF_REQUIRED" ]]; then
        printf '[OrthoFinder remote] Missing resource/profile: %s\n' "$OF_REQUIRED" >&2
        printf '[OrthoFinder remote] Follow docs/denbi.md or pass --resource-dir /absolute/server/resource.\n' >&2
        exit 1
    fi
done
mkdir -p "$OF_LOG_ROOT"
OF_ARGUMENTS=(
    --resource-dir "$OF_RESOURCES" --research-config "$OF_CONFIG"
    --cache-dir "$OF_STATE/cache/application" --log-file "$OF_LOG_ROOT/application.log"
    --server-address 127.0.0.1 --server-port "$OF_PORT" --headless
)
if [[ -n "$OF_SEQUENCE" ]]; then
    OF_ARGUMENTS+=(--terminal-motif-parquet "$OF_SEQUENCE")
fi
if [[ -n "$OF_TAXONOMY" ]]; then
    OF_ARGUMENTS+=(--taxonomy-map "$OF_TAXONOMY")
fi
OF_CONFIG_HASH="$(sha256sum "$OF_CONFIG" | cut -d ' ' -f 1)"
OF_WORKER_SIGNATURE="$OF_ENV|$OF_CONFIG|$OF_CONFIG_HASH"
OF_APP_SIGNATURE="$OF_WORKER_SIGNATURE|$OF_RESOURCES|$OF_PORT|$OF_SEQUENCE|$OF_TAXONOMY"
# Every managed window has an ownership marker and its exact launch settings.
for OF_TARGET in "$OF_APP_WINDOW" "$OF_WORKER_WINDOW"; do
    OF_DEAD="$(of_window_dead "$OF_TARGET")"
    if [[ "$OF_DEAD" != "absent" ]]; then
        OF_OWNER="$("$OF_TMUX" show-options -w -v -t "$OF_TARGET" @orthofinder_owner 2>/dev/null || true)"
        if [[ "$OF_OWNER" != "orthofinder-launcher-1" ]]; then
            printf '[OrthoFinder remote] Refusing to modify an unrecognised window: %s\n' "$OF_TARGET" >&2
            exit 1
        fi
    fi
done
OF_WORKER_DEAD="$(of_window_dead "$OF_WORKER_WINDOW")"
if [[ "$OF_WORKER_DEAD" == "0" ]]; then
    OF_EXISTING="$("$OF_TMUX" show-options -w -v -t "$OF_WORKER_WINDOW" @orthofinder_signature)"
    if [[ "$OF_EXISTING" != "$OF_WORKER_SIGNATURE" ]]; then
        printf '[OrthoFinder remote] The live worker uses different environment/profile settings.\n' >&2
        printf '[OrthoFinder remote] Let its jobs finish; see docs/denbi.md before replacing the worker.\n' >&2
        exit 1
    fi
fi
OF_APP_DEAD="$(of_window_dead "$OF_APP_WINDOW")"
if [[ "$OF_APP_DEAD" != "0" || "$OF_RESTART" == "yes" ]]; then
    printf '[OrthoFinder remote] Validating the completed OrthoFinder resource, profile and sidecars...\n'
    if ! "$OF_APP" "${OF_ARGUMENTS[@]}" --validate-only; then
        printf '[OrthoFinder remote] Resource validation failed; resolve the error above and retry.\n' >&2
        exit 1
    fi
fi

if [[ "$OF_APP_DEAD" == "0" && "$OF_RESTART" != "yes" ]]; then
    OF_EXISTING="$("$OF_TMUX" show-options -w -v -t "$OF_APP_WINDOW" @orthofinder_signature)"
    if [[ "$OF_EXISTING" != "$OF_APP_SIGNATURE" ]]; then
        printf '[OrthoFinder remote] App settings changed; use --restart-app to apply them.\n' >&2
        exit 1
    fi
fi
if [[ "$OF_APP_DEAD" == "0" && "$OF_RESTART" == "yes" ]]; then
    OF_STARTED_COMMAND="$("$OF_TMUX" list-panes -t "$OF_APP_WINDOW" -F '#{pane_start_command}')"
    if [[ "$OF_STARTED_COMMAND" != *"$OF_APP"* \
        || "$OF_STARTED_COMMAND" != *"$OF_CONFIG"* ]]; then
        printf '[OrthoFinder remote] Refusing to restart an unrecognised command in app.\n' >&2
        exit 1
    fi
    printf '[OrthoFinder remote] Restarting the OrthoFinder GUI; the research worker continues running.\n'
    "$OF_TMUX" kill-pane -t "$OF_APP_WINDOW"
    OF_APP_DEAD="absent"
fi

# Explicit cache/temp paths keep newly started services on the large data disk.
OF_TEMP="$OF_STATE/tmp"
OF_CACHE="$OF_STATE/cache/xdg"
mkdir -p "$OF_TEMP" "$OF_CACHE" "$OF_LOG_ROOT"
OF_SERVICE_ENV=(env "PATH=$OF_ENV/bin:$PATH" "TMPDIR=$OF_TEMP" "XDG_CACHE_HOME=$OF_CACHE" PYTHONUNBUFFERED=1)
OF_APP_COMMAND=(
    "${OF_SERVICE_ENV[@]}" "OF_LAUNCHER_LOG=$OF_APP_LOG"
    bash -c 'exec "$@" >>"$OF_LAUNCHER_LOG" 2>&1' of-app
    "$OF_APP" "${OF_ARGUMENTS[@]}"
)
OF_WORKER_COMMAND=(
    "${OF_SERVICE_ENV[@]}" "OF_LAUNCHER_LOG=$OF_WORKER_LOG"
    bash -c 'exec "$@" >>"$OF_LAUNCHER_LOG" 2>&1' of-worker
    "$OF_WORKBENCH" --config "$OF_CONFIG" worker
)
OF_WORKER_DEAD="$(of_window_dead "$OF_WORKER_WINDOW")"
if [[ "$OF_WORKER_DEAD" == "1" ]]; then
    "$OF_TMUX" respawn-pane -t "$OF_WORKER_WINDOW" \
        "${OF_WORKER_COMMAND[@]}"
elif [[ "$OF_WORKER_DEAD" != "0" ]]; then
    if "$OF_TMUX" has-session -t "$OF_SESSION" 2>/dev/null; then
        "$OF_TMUX" new-window -d -t "$OF_SESSION:" -n worker \
            "${OF_WORKER_COMMAND[@]}"
    else
        "$OF_TMUX" new-session -d -s "$OF_SESSION" -n worker \
            "${OF_WORKER_COMMAND[@]}"
    fi
fi
if [[ "$OF_WORKER_DEAD" != "0" ]]; then
    "$OF_TMUX" set-option -w -t "$OF_WORKER_WINDOW" @orthofinder_owner orthofinder-launcher-1
    "$OF_TMUX" set-option -w -t "$OF_WORKER_WINDOW" @orthofinder_signature "$OF_WORKER_SIGNATURE"
    "$OF_TMUX" set-option -w -t "$OF_WORKER_WINDOW" remain-on-exit on
fi
if [[ "$OF_APP_DEAD" == "1" ]]; then
    "$OF_TMUX" respawn-pane -t "$OF_APP_WINDOW" \
        "${OF_APP_COMMAND[@]}"
elif [[ "$OF_APP_DEAD" != "0" ]]; then
    "$OF_TMUX" new-window -d -t "$OF_SESSION:" -n app \
        "${OF_APP_COMMAND[@]}"
fi

if [[ "$OF_APP_DEAD" != "0" ]]; then
    "$OF_TMUX" set-option -w -t "$OF_APP_WINDOW" @orthofinder_owner orthofinder-launcher-1
    "$OF_TMUX" set-option -w -t "$OF_APP_WINDOW" @orthofinder_signature "$OF_APP_SIGNATURE"
    "$OF_TMUX" set-option -w -t "$OF_APP_WINDOW" remain-on-exit on
fi

printf '[OrthoFinder remote] Waiting for the OrthoFinder interface...\n'
if ! of_wait_for_app; then
    of_diagnostics "$OF_APP_WINDOW" "$OF_APP_LOG"
    printf '[OrthoFinder remote] OrthoFinder did not become ready on Sabrina port %s.\n' "$OF_PORT" >&2
    printf '[OrthoFinder remote] Use --restart-app to replace an earlier OrthoFinder GUI started on another port.\n' >&2
    exit 1
fi
OF_WORKER_DEAD="$(of_window_dead "$OF_WORKER_WINDOW")"
if [[ "$OF_WORKER_DEAD" != "0" ]]; then
    of_diagnostics "$OF_WORKER_WINDOW" "$OF_WORKER_LOG"
    printf '[OrthoFinder remote] The research worker stopped during startup.\n' >&2
    exit 1
fi
printf '[OrthoFinder remote] OrthoFinder app and research worker are running.\n'
OF_REMOTE
then
    printf '[OrthoFinder] Could not prepare the remote app. Check the diagnostic above.\n' >&2
    exit 1
fi

if [[ "$OF_MODE" == "status" ]]; then
    if ! ssh -F /dev/null -S "$OF_CONTROL_PATH" -O check "$OF_SSH_HOST"; then
        printf '[OrthoFinder] No tunnel owned by this launcher. Run the launcher to connect.\n' >&2
        exit 1
    fi
else
    OF_NEW_MASTER="no"
    if ssh -F /dev/null -S "$OF_CONTROL_PATH" -O check "$OF_SSH_HOST" >/dev/null 2>&1; then
        printf '[OrthoFinder] Reusing this launcher\047s SSH connection.\n'
    else
        # Keep alias, identity, proxy and host-key settings for authentication.
        # Clear inherited forwards before creating this dedicated connection.
        if ! ssh -fNT -M -S "$OF_CONTROL_PATH" \
            -o ControlPersist=yes -o ClearAllForwardings=yes \
            -o ConnectTimeout=15 \
            -o ServerAliveInterval=15 -o ServerAliveCountMax=3 "$OF_SSH_HOST"; then
            printf '[OrthoFinder] Could not establish the launcher\047s SSH connection.\n' >&2
            printf '[OrthoFinder] The remote OrthoFinder services remain running.\n' >&2
            exit 1
        fi
        OF_NEW_MASTER="yes"
    fi
    # Control requests use the already authenticated socket. Reading no config
    # here prevents inherited forwards from being added back. Setting
    # ClearAllForwardings=yes with -L would also erase the requested OF forward.
    # OpenSSH acknowledges an identical existing forward, so this is repeatable
    # and also repairs a master left without a forward by an interrupted run.
    printf '[OrthoFinder] Connecting Mac port %s to Sabrina port %s...\n' "$OF_LOCAL_PORT" "$OF_REMOTE_PORT"
    if ! ssh -F /dev/null -S "$OF_CONTROL_PATH" -O forward \
        -o ExitOnForwardFailure=yes \
        -L "127.0.0.1:$OF_LOCAL_PORT:127.0.0.1:$OF_REMOTE_PORT" "$OF_SSH_HOST"; then
        printf '[OrthoFinder] Tunnel failed for Mac port %s. Choose another --local-port if it is occupied.\n' "$OF_LOCAL_PORT" >&2
        if [[ "$OF_NEW_MASTER" == "yes" ]]; then
            if ! ssh -F /dev/null -S "$OF_CONTROL_PATH" -O exit "$OF_SSH_HOST" >/dev/null 2>&1; then
                printf '[OrthoFinder] Could not close the new SSH connection; use --stop-tunnel to retry.\n' >&2
            fi
        fi
        printf '[OrthoFinder] The remote OrthoFinder services remain running.\n' >&2
        exit 1
    fi
fi

if ! curl --fail --silent --show-error --max-time 3 \
    --retry 5 --retry-connrefused --retry-delay 1 --retry-max-time 15 \
    "$OF_URL/_stcore/health" >/dev/null; then
    printf '[OrthoFinder] The tunnel did not reach OF. Try --stop-tunnel, then run the launcher again.\n' >&2
    exit 1
fi

printf '\n[OrthoFinder] OrthoFinder app: %s\n' "$OF_URL"
printf '[OrthoFinder] Open Protein research workbench > Start here > OrthoFinder groups.\n'
if [[ "$OF_MODE" == "start" && "$OF_OPEN_BROWSER" == "yes" ]]; then
    if [[ "$(uname -s)" == "Darwin" ]] && command -v open >/dev/null 2>&1; then
        if ! open "$OF_URL"; then
            printf '[OrthoFinder] Browser could not open automatically; open the URL above.\n' >&2
        fi
    fi
fi

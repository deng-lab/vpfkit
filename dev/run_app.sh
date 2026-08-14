#!/usr/bin/env bash
# dev/run_app.sh — start the ViroProfiler-viewer Shiny app from the working tree.
#
#   bash dev/run_app.sh                      # localhost:7474, foreground, tails the log
#   bash dev/run_app.sh --lan                # 0.0.0.0:7474, reachable from the LAN
#   bash dev/run_app.sh --lan --detach       # same, but returns and keeps running
#   bash dev/run_app.sh --port 8000 --detach
#   bash dev/run_app.sh --stop               # stop whatever this script started
#   bash dev/run_app.sh --status
#
# --lan binds every interface. The app can then also read any `.rds` the account
# running it can read, through the "Path on this server" input, to anyone who can
# reach the port. That is the point on a trusted lab network and the risk anywhere
# else; pass --no-server-path to turn that input off while keeping upload working.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOG="$REPO_ROOT/.shiny_dev.log"
PID_FILE="$REPO_ROOT/.shiny_dev.pid"

PORT=7474
HOST=127.0.0.1
DETACH=0
BROWSER=1
SERVER_PATH=TRUE
ACTION=start

while [[ $# -gt 0 ]]; do
  case "$1" in
    --lan)             HOST=0.0.0.0; BROWSER=0; shift ;;
    --host)            HOST="$2"; shift 2 ;;
    --port|-p)         PORT="$2"; shift 2 ;;
    --detach|-d)       DETACH=1; BROWSER=0; shift ;;
    --no-browser)      BROWSER=0; shift ;;
    --no-server-path)  SERVER_PATH=FALSE; shift ;;
    --stop)            ACTION=stop; shift ;;
    --status)          ACTION=status; shift ;;
    -h|--help)         sed -n '2,16p' "$0"; exit 0 ;;
    *)                 echo "[ERROR] unknown argument: $1" >&2; exit 2 ;;
  esac
done

die()  { echo "[ERROR] $*" >&2; exit 1; }
info() { echo "[vpfkit] $*"; }

running_pid() {
  [[ -f "$PID_FILE" ]] || return 1
  local p; p="$(cat "$PID_FILE")"
  kill -0 "$p" 2>/dev/null || return 1
  echo "$p"
}

stop_app() {
  local p
  if p="$(running_pid)"; then
    info "stopping PID $p"
    # The R process spawns nothing, but kill the group in case that changes.
    kill "$p" 2>/dev/null || true
    sleep 1
    kill -0 "$p" 2>/dev/null && kill -9 "$p" 2>/dev/null || true
  else
    info "not running"
  fi
  rm -f "$PID_FILE"
}

case "$ACTION" in
  stop)   stop_app; exit 0 ;;
  status)
    if p="$(running_pid)"; then
      info "running, PID $p"
      ss -ltnp 2>/dev/null | grep -F "pid=$p" || true
    else
      info "not running"
    fi
    exit 0 ;;
esac

# A previous instance holds the port; replace it rather than failing.
if running_pid >/dev/null; then
  info "replacing the instance already started by this script"
  stop_app
fi
if command -v ss >/dev/null && ss -ltn "sport = :$PORT" | grep -q LISTEN; then
  die "port $PORT is in use by something else. Pass --port <n>."
fi

info "starting on ${HOST}:${PORT} (server-path input: ${SERVER_PATH})"
info "log: $LOG"

# setsid so the app outlives the shell that launched it, which is what --detach
# is for; without it a background job dies when the terminal or the agent session
# that started it goes away.
setsid R --no-save --no-restore -q -e "
  options(shiny.port = ${PORT}, shiny.host = '${HOST}', golem.app.prod = FALSE)
  pkgload::load_all('${REPO_ROOT}', quiet = TRUE)
  run_app(allow_server_path = ${SERVER_PATH})
" >"$LOG" 2>&1 < /dev/null &

APP_PID=$!
echo "$APP_PID" > "$PID_FILE"

PROBE_HOST="$HOST"
[[ "$HOST" == "0.0.0.0" ]] && PROBE_HOST=127.0.0.1
for _ in $(seq 1 120); do
  curl -sf "http://${PROBE_HOST}:${PORT}" -o /dev/null 2>/dev/null && break
  kill -0 "$APP_PID" 2>/dev/null || { cat "$LOG"; die "R exited during startup"; }
  sleep 0.5
done
curl -sf "http://${PROBE_HOST}:${PORT}" -o /dev/null 2>/dev/null \
  || { tail -30 "$LOG"; die "app did not answer within 60s"; }

info "ready, PID $APP_PID"
if [[ "$HOST" == "0.0.0.0" ]]; then
  ip -4 -o addr show scope global 2>/dev/null |
    awk -v p="$PORT" '{split($4,a,"/"); print "[vpfkit] http://" a[1] ":" p}'
else
  info "http://${HOST}:${PORT}"
fi

if (( BROWSER )); then
  for opener in xdg-open open start; do
    command -v "$opener" >/dev/null && { "$opener" "http://${PROBE_HOST}:${PORT}"; break; }
  done
fi

if (( DETACH )); then
  info "detached. stop with: bash dev/run_app.sh --stop"
  exit 0
fi

info "Ctrl-C to stop"
trap 'stop_app' EXIT INT TERM
tail -f "$LOG"

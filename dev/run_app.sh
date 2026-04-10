#!/usr/bin/env bash
# dev/run_app.sh — start vpfkit Shiny app and open it in the browser
# Usage: bash dev/run_app.sh [port] [rds_file]
#   port      optional port number (default: 7474)
#   rds_file  optional path to a .rds TSE file to pre-load (not yet wired — placeholder)

set -euo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
PORT="${1:-7474}"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOG="$REPO_ROOT/.shiny_dev.log"
PID_FILE="$REPO_ROOT/.shiny_dev.pid"

# ── Helpers ───────────────────────────────────────────────────────────────────
die()  { echo "[ERROR] $*" >&2; exit 1; }
info() { echo "[vpfkit] $*"; }

# ── Kill any previous instance ────────────────────────────────────────────────
if [[ -f "$PID_FILE" ]]; then
  OLD_PID=$(cat "$PID_FILE")
  if kill -0 "$OLD_PID" 2>/dev/null; then
    info "Stopping previous app instance (PID $OLD_PID)…"
    kill "$OLD_PID" 2>/dev/null || true
    sleep 1
  fi
  rm -f "$PID_FILE"
fi

# ── Check port is free ────────────────────────────────────────────────────────
if lsof -i :"$PORT" -sTCP:LISTEN -t &>/dev/null; then
  die "Port $PORT is already in use. Pass a different port: bash dev/run_app.sh <port>"
fi

# ── Start app ─────────────────────────────────────────────────────────────────
info "Starting vpfkit on http://localhost:$PORT …"
info "Log: $LOG"

R --no-save --no-restore -q -e "
  options(shiny.port = $PORT, shiny.host = '127.0.0.1', golem.app.prod = FALSE)
  devtools::load_all('$REPO_ROOT', quiet = TRUE)
  run_app()
" >"$LOG" 2>&1 &

APP_PID=$!
echo "$APP_PID" > "$PID_FILE"

# ── Wait for the app to be ready ──────────────────────────────────────────────
MAX_WAIT=30
WAITED=0
until curl -sf "http://localhost:$PORT" -o /dev/null 2>/dev/null; do
  sleep 0.5
  WAITED=$(( WAITED + 1 ))
  if (( WAITED >= MAX_WAIT * 2 )); then
    die "App did not start within ${MAX_WAIT}s. Check $LOG for details."
  fi
  # Bail early if R exited with an error
  if ! kill -0 "$APP_PID" 2>/dev/null; then
    die "R process exited unexpectedly. Check $LOG for details."
  fi
done

info "App is ready (PID $APP_PID)"

# ── Open browser ──────────────────────────────────────────────────────────────
URL="http://localhost:$PORT"
if command -v open &>/dev/null; then          # macOS
  open "$URL"
elif command -v xdg-open &>/dev/null; then    # Linux
  xdg-open "$URL"
elif command -v start &>/dev/null; then       # Windows Git Bash
  start "$URL"
fi
info "Opened $URL"

# ── Tail log until Ctrl-C ─────────────────────────────────────────────────────
info "Press Ctrl-C to stop the app."
cleanup() {
  info "Stopping app (PID $APP_PID)…"
  kill "$APP_PID" 2>/dev/null || true
  rm -f "$PID_FILE"
}
trap cleanup EXIT INT TERM

tail -f "$LOG"

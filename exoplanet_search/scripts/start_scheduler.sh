#!/usr/bin/env bash
# Standalone Dask scheduler on the coordinator machine.
# SCHED_HOST defaults to localhost. Set it to the LAN or Tailscale IP ONLY when you want other
# machines to join; the scheduler has no auth, so never bind it to a public interface.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${SCHED_HOST:=127.0.0.1}" "${SCHED_PORT:=8786}" "${DASH_PORT:=8797}"
exec nice -n 19 .venv/bin/dask scheduler --host "$SCHED_HOST" --port "$SCHED_PORT" \
  --dashboard-address "$SCHED_HOST:$DASH_PORT"

#!/usr/bin/env bash
# Apple Silicon worker (macOS). Advertises MLX=1 so --backend auto uses the mlx.core BLS.
# Add CPU=... style extra processes with start_worker_cpu.sh if you also want CPU stars here.
# UNTESTED until the Mac is online.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${SCHEDULER:?set SCHEDULER=tcp://<coordinator-ip>:8786}" "${NWORKERS:=2}"
export EXO_DATA_ROOT="${EXO_DATA_ROOT:-$HOME/Library/Caches/exosearch}"
exec .venv/bin/dask worker "$SCHEDULER" --nworkers "$NWORKERS" --nthreads 1 --resources "MLX=1" \
  --memory-limit 8GiB --name "mlx-$(hostname -s)" --local-directory /tmp/dask-exosearch

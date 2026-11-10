#!/usr/bin/env bash
# GPU worker (e.g. GTX 1050 Ti, 4 GB). One process owns the GPU and advertises GPU=1, so
# --backend auto runs its stars with the CuPy BLS. UNTESTED until the GPU host is online.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${SCHEDULER:?set SCHEDULER=tcp://<coordinator-ip>:8786}"
export EXO_DATA_ROOT="${EXO_DATA_ROOT:-$HOME/.cache/exosearch}"   # local cache if no shared disk
exec .venv/bin/dask worker "$SCHEDULER" --nworkers 1 --nthreads 1 --resources "GPU=1" \
  --memory-limit 6GiB --name "cuda-$(hostname -s)" --local-directory /tmp/dask-exosearch

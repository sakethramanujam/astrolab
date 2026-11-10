#!/usr/bin/env bash
# CPU workers. Default 8 single-threaded processes = leaves most of a 20-core host free.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${SCHEDULER:=tcp://127.0.0.1:8786}" "${NWORKERS:=8}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
exec nice -n 19 .venv/bin/dask worker "$SCHEDULER" --nworkers "$NWORKERS" --nthreads 1 \
  --memory-limit 2GiB --name "cpu-$(hostname -s)" --local-directory /tmp/dask-exosearch

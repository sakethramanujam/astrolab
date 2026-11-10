#!/usr/bin/env bash
# One-command NVIDIA GPU worker bootstrap (tested target: GTX 1050 Ti, Pascal sm_61, 4 GB).
#   SCHEDULER=tcp://<scheduler-host>:8786 scripts/worker_cuda.sh
# Checks the driver and CUDA/CuPy, builds .venv if needed, runs a self-test benchmark on real
# TESS stars (CUDA BLS vs astropy), and only then joins the scheduler with --resources GPU=1.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${SCHEDULER:?set SCHEDULER=tcp://<scheduler-host>:8786}"
: "${PYTHON:=python3.12}" "${DASK_VERSION:=2026.8.0}" "${MIN_DRIVER:=525}"
export EXO_DATA_ROOT="${EXO_DATA_ROOT:-$HOME/.cache/exosearch}"
die(){ echo "worker_cuda: $*" >&2; exit 1; }

command -v nvidia-smi >/dev/null || die "nvidia-smi not found: install the NVIDIA driver (>= $MIN_DRIVER)"
IFS=, read -r NAME DRIVER CAP MEM < <(nvidia-smi --query-gpu=name,driver_version,compute_cap,memory.total --format=csv,noheader | head -1)
echo "GPU:$NAME  driver:$DRIVER  compute capability:$CAP  memory:$MEM"
(( ${DRIVER%%.*} >= MIN_DRIVER )) || die "driver $DRIVER < $MIN_DRIVER; CUDA 12 wheels need a newer driver"
awk -v c="$CAP" 'BEGIN{exit !(c+0 >= 5.0)}' || die "compute capability $CAP < 5.0 is not supported by CUDA 12"
awk -v c="$CAP" 'BEGIN{exit !(c+0 < 7.0)}' && echo "note: Pascal-class GPU (sm_${CAP/./}); cupy-cuda12x still ships kernels for it"

if [[ ! -x .venv/bin/python ]]; then "$PYTHON" -m venv .venv; fi
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q -r requirements.txt -r requirements-cuda.txt "dask[distributed]==$DASK_VERSION" "distributed==$DASK_VERSION"
.venv/bin/python - <<'PY'
import cupy as cp
p = cp.cuda.runtime.getDeviceProperties(0)
x = cp.arange(1_000_000, dtype=cp.float64)
assert float(cp.sum(x)) == 499999500000.0
print(f"CuPy {cp.__version__}: {p['name'].decode()} sm_{p['major']}{p['minor']}, CUDA runtime {cp.cuda.runtime.runtimeGetVersion()}")
PY
.venv/bin/python -m exosearch.selftest --backend cuda || die "self-test failed: not joining"
exec env SCHEDULER="$SCHEDULER" scripts/start_worker_cuda.sh

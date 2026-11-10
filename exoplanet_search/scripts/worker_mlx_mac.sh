#!/usr/bin/env bash
# One-command Apple Silicon (MLX) worker bootstrap for macOS.
#   SCHEDULER=tcp://<scheduler-host>:8786 scripts/worker_mlx_mac.sh
# Uses Homebrew Python 3.12, builds .venv, checks MLX sees the GPU, runs a self-test benchmark on
# real TESS stars (MLX BLS vs astropy), then joins the scheduler with --resources MLX=1.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${SCHEDULER:?set SCHEDULER=tcp://<scheduler-host>:8786}"
: "${DASK_VERSION:=2026.8.0}" "${NWORKERS:=2}"
export EXO_DATA_ROOT="${EXO_DATA_ROOT:-$HOME/Library/Caches/exosearch}"
die(){ echo "worker_mlx_mac: $*" >&2; exit 1; }
[[ "$(uname -s)" == Darwin && "$(uname -m)" == arm64 ]] || die "needs macOS on Apple Silicon"
PY="${PYTHON:-$(command -v python3.12 || true)}"
if [[ -z "$PY" ]]; then
  command -v brew >/dev/null || die "install Homebrew first (https://brew.sh)"
  brew install python@3.12; PY="$(brew --prefix)/bin/python3.12"
fi
if [[ ! -x .venv/bin/python ]]; then "$PY" -m venv .venv; fi
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q -r requirements.txt -r requirements-mlx.txt "dask[distributed]==$DASK_VERSION" "distributed==$DASK_VERSION"
.venv/bin/python -c "import mlx.core as mx; d=mx.default_device(); print('MLX', mx.__version__, d); assert 'gpu' in str(d).lower(), 'MLX is not using the GPU'"
.venv/bin/python -m exosearch.selftest --backend mlx || die "self-test failed: not joining"
exec env SCHEDULER="$SCHEDULER" NWORKERS="$NWORKERS" scripts/start_worker_mlx.sh

#!/usr/bin/env bash
# Build the isolated worker venv from the committed lockfile (user-level only; no sudo, no system Python edits).
#   scripts/worker_setup.sh [python3.12 interpreter]      (run from the exoplanet_search/ directory of a clone)
# The lock (requirements.lock) is universal (linux x86_64 / macOS arm64, CPython 3.12) and hash-pinned; it is
# generated on the coordinator with: uv pip compile requirements.txt -c <coordinator freeze> --universal
#   --python-version 3.12 --generate-hashes -o requirements.lock
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${1:-python3}"
"$PY" -c 'import sys; assert sys.version_info[:2] == (3, 12), sys.version' 
"$PY" -m venv .venv-worker
.venv-worker/bin/python -m pip install -q --upgrade pip
.venv-worker/bin/python -m pip install -q --require-hashes --no-deps -r requirements.lock
.venv-worker/bin/python -m pip check
.venv-worker/bin/python -c "import numpy, scipy, astropy, numba, wotan, lightkurve, batman; print('worker venv ok', numpy.__version__, scipy.__version__)"

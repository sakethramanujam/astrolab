#!/usr/bin/env bash
# Run one host's share of the v1.1 production manifest (clean checkout at the pinned commit, worker venv from
# requirements.lock, nice 19). Production refuses to start while configs/v1.1.json has no beta option.
#   scripts/v11_work.sh HOST_ID MANIFEST WORKERS OUTDIR
# Remote hosts set EXO_DATA_ROOT to their own data directory. On the coordinator (h1) run scripts/v11_backoff.sh
# next to it.
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=. EXO_HOST_ID="$1" PYTHONUNBUFFERED=1
exec nice -n 19 .venv-worker/bin/python -m exosearch.v11.run work --manifest "$2" --workers "$3" --out "$4"

#!/usr/bin/env bash
# Incremental nightly batch: already-searched stars are skipped, so this only does new work.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${SECTORS:=70}" "${N_TARGETS:=2000}" "${WORKERS:=8}"
for s in $SECTORS; do
  if [[ -n "${SCHEDULER:-}" ]]; then
    nice -n 19 .venv/bin/python run_search.py --sector "$s" --n-targets "$N_TARGETS" --scheduler "$SCHEDULER"
  else
    nice -n 19 .venv/bin/python run_search.py --sector "$s" --n-targets "$N_TARGETS" --workers "$WORKERS"
  fi
done

#!/usr/bin/env bash
# Nightly batch: work through the sector queue (resumes automatically), then publish the site.
#   MAX_SECTORS  sectors per night (1-3)          WORKERS   Dask worker processes (<= 8 here)
#   DEADLINE     don't start a new sector after HH (24h local)
#   SCHEDULER    tcp://<scheduler-host>:8786 to use an external scheduler + remote workers
set -euo pipefail
cd "$(dirname "$0")/.."
: "${MAX_SECTORS:=2}" "${WORKERS:=8}" "${DEADLINE:=6}"
exec 9>/tmp/exosearch-nightly.lock
flock -n 9 || { echo "another run is active"; exit 0; }
PY=(nice -n 19 ionice -c3 .venv/bin/python)
for s in $(.venv/bin/python -m exosearch.sectorqueue next --count "$MAX_SECTORS"); do
  h=$(date +%-H)
  if (( h >= DEADLINE && h < 18 )); then echo "past deadline, stopping"; break; fi
  echo "=== sector $s ($(date))"
  if [[ -n "${SCHEDULER:-}" ]]; then
    "${PY[@]}" run_search.py --sector "$s" --n-targets 0 --scheduler "$SCHEDULER" --prune-fits --plots 15
  else
    "${PY[@]}" run_search.py --sector "$s" --n-targets 0 --workers "$WORKERS" --prune-fits --plots 15
  fi
done
nice -n 19 scripts/publish_site.sh

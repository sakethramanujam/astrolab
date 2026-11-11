#!/usr/bin/env bash
# Run reliability batch steps one after another (each step is resumable, so a killed lane just restarts).
#   scripts/reliability_lane.sh [--after PID] "cmd 1" "cmd 2" ...
# Meant to be launched as a transient systemd user unit so it outlives the login session, e.g.
#   systemd-run --user --unit=exo-lane-a -p Nice=19 -p CPUSchedulingPolicy=idle -p IOSchedulingClass=idle \
#     --working-directory="$PWD" scripts/reliability_lane.sh "python -m exosearch.nulltest --mode scramble --realization 0 --workers 6"
# Every step registers its worker claim (exosearch/batch.py) so scripts/continuous.sh stays within COMBINED_CAP.
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=. OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
: "${MAX_LOAD:=18}"
if [[ "${1:-}" == "--after" ]]; then
  while kill -0 "$2" 2>/dev/null; do sleep 60; done
  shift 2
fi
for cmd in "$@"; do
  while awk -v l="$(cut -d' ' -f1 /proc/loadavg)" -v m="$MAX_LOAD" 'BEGIN{exit !(l>m)}'; do
    echo "$(date +%T) load above $MAX_LOAD, waiting"; sleep 120
  done
  echo "$(date +%T) START $cmd"
  for attempt in 1 2 3; do
    .venv/bin/$cmd && break
    echo "$(date +%T) step failed (attempt $attempt), retrying in 5 min"; sleep 300
  done
  echo "$(date +%T) DONE  $cmd"
done

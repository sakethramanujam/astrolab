#!/usr/bin/env bash
# Run steps once batch worker claims leave room under the batch share of the core cap.
#   scripts/queue_when_free.sh WORKERS "cmd 1" "cmd 2" ...
# Starts only after (claimed batch workers + WORKERS) <= BATCH_CAP (12 of the 16-core cap; the runner keeps 4) and the
# 1-min load <= MAX_LOAD, both holding for SUSTAIN seconds (so a lane's brief gap between steps does not count).
# Each step registers its own claim (exosearch/batch.py), so later lane steps and the runner see it.
# Launch as a transient user unit at idle priority, e.g.
#   systemd-run --user --unit=exo-subsample -p Nice=19 -p CPUSchedulingPolicy=idle -p IOSchedulingClass=idle \
#     --working-directory="$PWD" scripts/queue_when_free.sh 3 "python -m exosearch.subsample ..."
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=. OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
: "${CLAIM_DIR:=/tmp/exosearch-claims}" "${BATCH_CAP:=12}" "${MAX_LOAD:=18}" "${SUSTAIN:=600}" "${POLL:=60}"
W="$1"; shift
claimed(){ local n=0 f p; for f in "$CLAIM_DIR"/*.workers; do [[ -e "$f" ]] || continue; p=$(basename "$f" .workers)
  kill -0 "$p" 2>/dev/null && n=$((n + $(cat "$f"))); done; echo "$n"; }
room(){ local c; c=$(claimed); awk -v l="$(cut -d' ' -f1 /proc/loadavg)" -v m="$MAX_LOAD" -v c="$c" -v w="$W" -v cap="$BATCH_CAP" \
  'BEGIN{exit !((c + w) <= cap && l <= m)}'; }
for cmd in "$@"; do
  ok=0; echo "$(date +%T) waiting for room: $W workers (batch cap $BATCH_CAP, sustain ${SUSTAIN}s) for: $cmd"
  while (( ok < SUSTAIN )); do
    if room; then ok=$((ok + POLL)); else ok=0; fi
    sleep "$POLL"
  done
  echo "$(date +%T) START $cmd (claimed before start: $(claimed))"
  for attempt in 1 2 3; do
    .venv/bin/$cmd && break
    echo "$(date +%T) step failed (attempt $attempt), retrying in 5 min"; sleep 300
  done
  echo "$(date +%T) DONE  $cmd"
done

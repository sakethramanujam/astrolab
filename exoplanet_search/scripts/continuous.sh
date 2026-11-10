#!/usr/bin/env bash
# Work through the TESS sector queue back to back, around the clock, politely.
#  - one instance only (flock); everything at nice 19 / idle I/O (also set by the systemd unit)
#  - workers: DAY_WORKERS (07:00-23:00 local) or NIGHT_WORKERS, NIGHT only if load allows
#  - backs off while the 1-min load > MAX_LOAD or free space on the data disk < MIN_FREE_GB
#  - prefetches the next sector's light curves while the current one is searched
#  - prunes raw FITS after each sector and publishes the site
# Remote workers (opt-in): set EXO_LISTEN=<lan-or-tailscale-ip> to run a persistent scheduler on
# EXO_LISTEN:${EXO_PORT:-8786} that GPU/MLX/CPU workers elsewhere can join. Default: localhost only.
set -uo pipefail
cd "$(dirname "$0")/.."
: "${DAY_WORKERS:=8}" "${NIGHT_WORKERS:=12}" "${MAX_LOAD:=14}" "${MIN_FREE_GB:=100}"
: "${DATA_DISK:=${EXO_DATA_ROOT:-/media/nfs01/astro_data/tess}}" "${PUBLISH:=1}" "${IDLE_SLEEP:=1800}"
exec 9>/tmp/exosearch-run.lock
flock -n 9 || { echo "another exosearch run holds the lock"; exit 0; }
PY=(nice -n 19 ionice -c3 .venv/bin/python)
log(){ echo "$(date '+%F %T') $*"; }
load1(){ cut -d' ' -f1 /proc/loadavg; }
free_gb(){ df -BG --output=avail "$DATA_DISK" | tail -1 | tr -dc 0-9; }
wait_ok(){  # block until load and disk are fine
  while :; do
    local l f; l=$(load1); f=$(free_gb)
    if awk -v l="$l" -v m="$MAX_LOAD" 'BEGIN{exit !(l>m)}'; then log "load $l > $MAX_LOAD, backing off 5 min"; sleep 300; continue; fi
    if (( f < MIN_FREE_GB )); then log "only ${f} GB free < ${MIN_FREE_GB} GB, backing off 15 min"; sleep 900; continue; fi
    return 0
  done
}
workers(){
  local h; h=$(date +%-H)
  if (( h >= 23 || h < 7 )) && awk -v l="$(load1)" 'BEGIN{exit !(l<6)}'; then echo "$NIGHT_WORKERS"; else echo "$DAY_WORKERS"; fi
}

SCHED_ARGS=()
if [[ -n "${EXO_LISTEN:-}" ]]; then   # persistent scheduler + local CPU workers that remote workers can join
  PORT="${EXO_PORT:-8786}"
  nice -n 19 .venv/bin/dask scheduler --host "$EXO_LISTEN" --port "$PORT" --dashboard-address "$EXO_LISTEN:$((PORT+11))" &
  SCHED_PID=$!; sleep 5
  SCHEDULER="tcp://$EXO_LISTEN:$PORT" NWORKERS="$DAY_WORKERS" nice -n 19 scripts/start_worker_cpu.sh &
  WORK_PID=$!
  trap 'kill $WORK_PID $SCHED_PID 2>/dev/null' EXIT
  SCHED_ARGS=(--scheduler "tcp://$EXO_LISTEN:$PORT")
fi

PREFETCH_PID=""; PREFETCH_SECTOR=""
while :; do
  read -r -a NEXT <<<"$(.venv/bin/python -m exosearch.sectorqueue next --count 2)"
  if (( ${#NEXT[@]} == 0 )); then log "queue empty; sleeping"; sleep "$IDLE_SLEEP"; continue; fi
  S="${NEXT[0]}"; N2="${NEXT[1]:-}"
  if [[ -n "$PREFETCH_PID" ]]; then wait "$PREFETCH_PID" 2>/dev/null; PREFETCH_PID=""; fi
  wait_ok
  W=$(workers)
  log "=== sector $S with ${W} workers (prefetching ${N2:-none})"
  if [[ -n "$N2" && "$N2" != "$PREFETCH_SECTOR" ]]; then
    "${PY[@]}" run_search.py --sector "$N2" --n-targets 0 --download-only --download-threads 6 \
      > "/tmp/exosearch-prefetch-s$N2.log" 2>&1 &
    PREFETCH_PID=$!; PREFETCH_SECTOR="$N2"
  fi
  if ((${#SCHED_ARGS[@]})); then
    "${PY[@]}" run_search.py --sector "$S" --n-targets 0 "${SCHED_ARGS[@]}" --prune-fits --plots 15 --download-threads 8
  else
    "${PY[@]}" run_search.py --sector "$S" --n-targets 0 --workers "$W" --prune-fits --plots 15 --download-threads 8
  fi
  rc=$?
  if (( rc != 0 )); then log "sector $S exited $rc; pausing 10 min"; sleep 600; continue; fi
  "${PY[@]}" -m exosearch.deepvet --sector "$S" --only-review || log "deep vetting failed for $S (continuing)"
  if [[ "$PUBLISH" == 1 ]]; then nice -n 19 scripts/publish_site.sh || log "publish failed (continuing)"; fi
done

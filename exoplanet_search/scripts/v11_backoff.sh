#!/usr/bin/env bash
# Coordinator back-off for v1.1 production: pause (SIGSTOP) the v1.1 worker process tree while the 5-minute load
# average stays above MAX_LOAD, and resume (SIGCONT) once it falls below RESUME_LOAD. Other lanes are never touched.
#   scripts/v11_backoff.sh PID [MAX_LOAD=18] [RESUME_LOAD=15]
set -uo pipefail
PID="$1"; MAX="${2:-18}"; RES="${3:-15}"; paused=0
tree() { echo "$PID"; pgrep -P "$PID"; }
while kill -0 "$PID" 2>/dev/null; do
  l5=$(cut -d' ' -f2 /proc/loadavg)
  if (( !paused )) && awk -v l="$l5" -v m="$MAX" 'BEGIN{exit !(l>m)}'; then
    kill -STOP $(tree) 2>/dev/null; paused=1; echo "$(date +%T) load5 $l5 > $MAX: paused"
  elif (( paused )) && awk -v l="$l5" -v m="$RES" 'BEGIN{exit !(l<m)}'; then
    kill -CONT $(tree) 2>/dev/null; paused=0; echo "$(date +%T) load5 $l5 < $RES: resumed"
  fi
  sleep 60
done
(( paused )) && kill -CONT $(tree) 2>/dev/null
echo "$(date +%T) worker $PID finished"

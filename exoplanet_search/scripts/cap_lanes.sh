#!/usr/bin/env bash
# Cap the v1.0 reliability lanes at what the interim layer needs (agreed 2026-11-10; v1.0 cannot produce the frozen S*).
# Touches no running step: it waits for the last wanted step of each lane to log DONE, then stops that lane before
# (or within seconds of) its next, unwanted step, and for lane b starts the reduced replacement step.
#   lane a: flipped null S70+S106 (running) + scrambled r0      -> drop scrambled r1, r2
#   lane b: injections run 0 (running, S70+S106, 8,000)        -> replace runs 1-7 by run 1, S70 only (4,000)
#   lane c: stellar + labels                                   -> drop injections runs 8-14
# Run as: systemd-run --user --unit=exo-cap -p Nice=19 -p CPUSchedulingPolicy=idle -p IOSchedulingClass=idle \
#           --working-directory="$PWD" scripts/cap_lanes.sh
set -uo pipefail
cd "$(dirname "$0")/.."
SINCE="$(date '+%Y-%m-%d %H:%M:%S')"
done_in() { journalctl --user -u "$1" --since "$SINCE" -o cat 2>/dev/null | grep -qF "DONE  $2"; }
a=0; b=0; c=0
while (( a + b + c < 3 )); do
  if (( !a )) && done_in exo-lane-a "python -m exosearch.nulltest --mode scramble --realization 0 "; then
    systemctl --user stop exo-lane-a; echo "$(date +%T) lane a capped after scrambled r0"; a=1; fi
  if (( !b )) && done_in exo-lane-b "python -m exosearch.inject --config configs/injection_v1.json --runs 0 "; then
    systemctl --user stop exo-lane-b
    systemd-run --user --unit=exo-lane-b2 -p Nice=19 -p CPUSchedulingPolicy=idle -p IOSchedulingClass=idle \
      --working-directory="$PWD" scripts/reliability_lane.sh \
      "python -m exosearch.inject --config configs/injection_v1.json --runs 1 --sectors 70 --workers 3"
    echo "$(date +%T) lane b capped: run 1 S70 only (exo-lane-b2)"; b=1; fi
  if (( !c )) && done_in exo-lane-c "python -m exosearch.labels --v1"; then
    systemctl --user stop exo-lane-c; echo "$(date +%T) lane c capped after labels"; c=1; fi
  sleep 3
done

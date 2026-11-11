#!/usr/bin/env bash
# Remote-host side of the unattended v1.1 launch (run every 5 min by a user timer / launchd agent).
# Starts this host's production share once, and only once, after the coordinator has pushed
# exoplanet_search/runs/v1.1/go.json (written only when the beta TOST gate passed and both v1.1 checks were committed
# as passing). Otherwise does nothing. State goes to $EXO_DATA_ROOT/../v11_agent_status.
#   scripts/v11_host_agent.sh HOST_ID
set -uo pipefail
HOST="$1"; SUB=exoplanet_search
cd "$(dirname "$0")/../.."                                  # repo root
DATA="${EXO_DATA_ROOT:-$PWD/$SUB/data/tess}"
STATE="$PWD/$SUB/data/v11_agent_status"; OUT="$PWD/$SUB/data/v11prod"
say() { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $*" > "$STATE"; echo "$*"; }
[ -f "$OUT/started" ] && { say "already started: $(cat "$OUT/started")"; exit 0; }
git fetch -q origin exoplanet-search || { say "fetch failed"; exit 0; }
GO=$(git show origin/exoplanet-search:$SUB/runs/v1.1/go.json 2>/dev/null) || { say "waiting: no go.json"; exit 0; }
read -r R W M < <(python3 -c "import json,sys; g=json.loads(sys.argv[1]); h=sys.argv[2]; print(g['pinned_commit'], g['workers'][h], g['manifests'][h])" "$GO" "$HOST")
[ -n "$(git status --porcelain --untracked-files=no)" ] && { say "refusing: tracked files modified"; exit 0; }
mkdir -p "$OUT"
git show "origin/exoplanet-search:$M" > "$OUT/manifest_$HOST.jsonl.gz" || { say "manifest missing"; exit 0; }
git checkout -q --detach "$R" || { say "checkout $R failed"; exit 0; }
echo "$R $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/started"
EXO_DATA_ROOT="$DATA" nohup "$SUB/scripts/v11_work.sh" "$HOST" "$OUT/manifest_$HOST.jsonl.gz" "$W" "$OUT/rows" > "$OUT/work.log" 2>&1 &
say "started $HOST at ${R:0:8} with $W workers (pid $!)"

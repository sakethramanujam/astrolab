#!/usr/bin/env bash
# Publish process-only v1.2 status to GitHub Pages (exoplanet-search branch).
# Runs on h1. Does not touch the local git checkout or the running search.
set -euo pipefail

ENV_SH="${V12_ENV_SH:-/tmp/v12-search-h1/env.sh}"
if [[ -f "$ENV_SH" ]]; then
  # shellcheck disable=SC1090
  source "$ENV_SH"
fi
: "${EXO_DATA_ROOT:=/media/nfs01/astro_data/tess}"
PY="${PY:-python3}"
GH="${GH:-gh}"
REPO="${GITHUB_REPOSITORY:-sakethramanujam/astrolab}"
BRANCH="${STATUS_BRANCH:-exoplanet-search}"
SCRIPT_PATH="exoplanet_search/scripts/publish_status_v12.py"
CACHE="${STATUS_SCRIPT_CACHE:-${EXO_DATA_ROOT}/site_publish/publish_status_v12.py}"

mkdir -p "$(dirname "$CACHE")"
if [[ "${STATUS_OFFLINE:-0}" != "1" ]]; then
  "$GH" api -H "Accept: application/vnd.github.raw" \
    "repos/${REPO}/contents/${SCRIPT_PATH}?ref=${BRANCH}" \
    > "${CACHE}.tmp"
  mv "${CACHE}.tmp" "$CACHE"
elif [[ ! -f "$CACHE" ]]; then
  echo "no cached publisher at $CACHE and STATUS_OFFLINE=1" >&2
  exit 1
fi
exec "$PY" "$CACHE" --publish "$@"

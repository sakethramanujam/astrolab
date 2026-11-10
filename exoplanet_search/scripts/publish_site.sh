#!/usr/bin/env bash
# Export the site summary and push ONLY site/data + site/plots to the exoplanet-search branch.
# The GitHub Actions Pages workflow then redeploys the site. Never pushes raw data.
set -euo pipefail
cd "$(dirname "$0")/.."
BRANCH="${BRANCH:-exoplanet-search}"
.venv/bin/python -m exosearch.export

# privacy + size guards: abort rather than publish anything machine-specific.
# Generic patterns live here; machine/user names live in an untracked file next to the data
# (one regex per line), so the guard itself never publishes them.
GENERIC='192\.168\.|(^|[^0-9.])10\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}|(^|[^0-9.])100\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}|/home/|/Users/'
PRIV="${EXO_PRIVATE_PATTERNS:-${EXO_DATA_ROOT:-/media/nfs01/astro_data/tess}/private_patterns.txt}"
PAT="$GENERIC"
if [[ -s "$PRIV" ]]; then PAT="$GENERIC|$(grep -v '^\s*$' "$PRIV" | paste -sd'|')"; fi
if grep -rIl -i -E "$PAT" site; then echo "privacy guard tripped; not publishing" >&2; exit 1; fi
find site -name '*.png' -size +300k -print -delete

REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"
[ "$(git branch --show-current)" = "$BRANCH" ] || { echo "not on $BRANCH; not publishing" >&2; exit 1; }
git add -A -- exoplanet_search/site/data exoplanet_search/site/plots
if git diff --cached --quiet -- exoplanet_search/site/data exoplanet_search/site/plots; then echo "nothing new to publish"; exit 0; fi
if git diff --cached | grep -i -E "$PAT" >/dev/null; then echo "privacy guard tripped in diff" >&2; git reset -q; exit 1; fi
# commit ONLY the site data paths, even if other work is staged in this checkout
git -c commit.gpgsign=false commit -q -m "data: site summary $(date -u +%Y-%m-%dT%H:%MZ)" -- exoplanet_search/site/data exoplanet_search/site/plots

# HTTPS + gh credentials work from timers without an ssh-agent
URL="$(git remote get-url origin | sed -E 's#^git@github.com:#https://github.com/#')"
CRED=(-c credential.helper= -c "credential.helper=!gh auth git-credential")
git "${CRED[@]}" pull -q --rebase --autostash "$URL" "$BRANCH"
git "${CRED[@]}" push -q "$URL" "HEAD:$BRANCH"
echo "published $(git rev-parse --short HEAD)"

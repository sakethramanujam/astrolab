"""Search-path hash for the S70 v1.1 beta-registration nulls. Usage: python searchpath_hash.py <repo> <commit>
sha256 over sorted 'path\0sha256(content)\n' lines for the search-path files read from <commit> via git show,
plus the canonical JSON (sort_keys, compact) of configs/v1.1.json without the non-search sections."""
import hashlib, json, subprocess, sys
FILES = sorted("""exoplanet_search/exosearch/__init__.py exoplanet_search/exosearch/batch.py exoplanet_search/exosearch/bls.py
exoplanet_search/exosearch/config.py exoplanet_search/exosearch/deepvet.py exoplanet_search/exosearch/density.py
exoplanet_search/exosearch/holdout.py exoplanet_search/exosearch/inject.py exoplanet_search/exosearch/interim.py
exoplanet_search/exosearch/labels.py exoplanet_search/exosearch/lightcurve.py exoplanet_search/exosearch/nulltest.py
exoplanet_search/exosearch/quality.py exoplanet_search/exosearch/stellar.py exoplanet_search/exosearch/store.py
exoplanet_search/exosearch/targets.py exoplanet_search/exosearch/task.py exoplanet_search/exosearch/version.py
exoplanet_search/exosearch/vet.py exoplanet_search/exosearch/worker_guard.py exoplanet_search/exosearch/v11/__init__.py
exoplanet_search/exosearch/v11/grid.py exoplanet_search/exosearch/v11/kernels.py exoplanet_search/exosearch/v11/pipeline.py
exoplanet_search/exosearch/v11/prewhiten.py exoplanet_search/exosearch/v11/run.py exoplanet_search/exosearch/v11/scramble.py
exoplanet_search/exosearch/v11/search.py exoplanet_search/exosearch/v11/settings.py exoplanet_search/exosearch/v11/vet11.py
exoplanet_search/configs/v1.1_beta_prior.json exoplanet_search/configs/scramble_exclude_s0070_2026-11-11.csv""".split())
NON_SEARCH_KEYS = ("beta_option", "beta_option_note", "validation_pw2", "known_planets")
repo, rev = sys.argv[1], sys.argv[2]
show = lambda p: subprocess.run(["git", "-C", repo, "show", f"{rev}:{p}"], check=True, capture_output=True).stdout
h = hashlib.sha256()
for p in FILES:
    h.update(f"{p}\0{hashlib.sha256(show(p)).hexdigest()}\n".encode())
cfg = json.loads(show("exoplanet_search/configs/v1.1.json"))
for k in NON_SEARCH_KEYS: cfg.pop(k, None)
h.update(b"configs/v1.1.json#search_keys\0" + json.dumps(cfg, sort_keys=True, separators=(",", ":")).encode())
print(h.hexdigest())

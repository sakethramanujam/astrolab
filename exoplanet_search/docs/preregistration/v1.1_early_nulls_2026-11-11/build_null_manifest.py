"""Null-only manifest builder for the S70 v1.1 beta-registration nulls (launch plan stage B1; run.build has no task
filter). Calls exosearch.v11.run.build unchanged (write=False, production=False: no injection hosts), keeps only the
flipped null (invert, seed 0) and scrambled realization 0 (already omitted for the registered excluded hosts), drops
stars left with no task, and writes manifests with run.write_manifests. No 'real' or injection tasks are emitted."""
import sys
from exosearch.v11 import run, settings
import json
out, weights = sys.argv[1], [(h, float(x)) for h, x in (p.split(":") for p in sys.argv[2].split(","))]
jobs = run.build([70], weights, out, production=False, write=False)
keep = []
for j in jobs:
    j["tasks"] = [t for t in j["tasks"] if t["tag"] in ("invert", "scramble_r0")]
    if j["tasks"]:
        keep.append(j)
cfg = settings.cfg(); cfg_inj = json.loads((settings.ROOT / "configs" / "injection_v1.1.json").read_text())
run.write_manifests(keep, out, cfg, cfg_inj, weights, name="nulls_s70")
n_inv = sum(t["tag"] == "invert" for j in keep for t in j["tasks"]); n_scr = sum(t["tag"] == "scramble_r0" for j in keep for t in j["tasks"])
print(json.dumps(dict(stars=len(keep), invert=n_inv, scramble_r0=n_scr, excluded_from_scramble=n_inv - n_scr)))

"""Re-run stored v1.0.0 rows with a pinned, clean checkout of the search code, into a scratch file (never the DB),
and diff against the stored rows. Used for the 378 S105 rows of the pool that loaded uncommitted code
(docs/v1.0_disclosure_ties_and_numba_cache.md).

  git worktree add --detach /tmp/wt-83eac4dd 83eac4dd
  python -m exosearch.rerun_rows run --code-root /tmp/wt-83eac4dd/exoplanet_search --sector 105 \
      --tics configs/rerun_queue_v1.0_s105_dirty.txt --workers 2 --out /tmp/rerun378.csv.gz
  python -m exosearch.rerun_rows diff /tmp/rerun378.csv.gz --sector 105 --out docs/host_gate/rerun_s105_dirty_diff.json
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CMP = ("period", "duration", "t0", "S", "depth", "sde", "snr", "n_transits", "flags", "passed", "git_hash")


def _init(code_root):
    sys.path.insert(0, code_root)
    for m in [m for m in sys.modules if m == "exosearch" or m.startswith("exosearch.")]:
        del sys.modules[m]


def _one(args):
    tic, sector, fname, url = args
    import exosearch.task as task                      # resolved from code_root (see _init)
    r = task.process_star(tic, sector, fname, url, "auto")
    r.pop("host", None)
    r["code_ok"] = str(Path(task.__file__).resolve()).startswith(sys.path[0])
    return r


def run(code_root, sector, tics_file, workers, out):
    from multiprocessing import get_context
    from . import targets
    code_root = str(Path(code_root).resolve())
    if subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=code_root,
                      capture_output=True, text=True).stdout.strip():
        raise SystemExit("code root is not a clean checkout")
    head = subprocess.run(["git", "rev-parse", "--short=10", "HEAD"], cwd=code_root, capture_output=True, text=True).stdout.strip()
    tics = [int(l.split()[0]) for l in Path(tics_file).read_text().splitlines() if l.strip() and not l.startswith("#")]
    man = targets.fetch_manifest(sector).set_index("tic")
    jobs = [(t, sector, man.loc[t].fname, man.loc[t].url) for t in tics]
    with get_context("spawn").Pool(workers, initializer=_init, initargs=(code_root,)) as pool:
        rows = list(pool.imap_unordered(_one, jobs))
    d = pd.DataFrame(rows)
    d["rerun_code"] = head
    if not d.code_ok.all():
        raise SystemExit("some workers imported the search code from outside the pinned checkout")
    d.to_csv(out, index=False)
    print(f"{len(d)} rows ({(d.status == 'ok').sum()} ok) from code {head} -> {out}")


def diff(new_file, sector, out):
    from . import store
    new = pd.read_csv(new_file).set_index("tic")
    old = store.frame(store.connect(), sector).set_index("tic").loc[new.index]
    ok = (new.status == "ok") & (old.status == "ok")
    n, o = new[ok], old[ok]
    same_pd = (n.period == o.period) & (n.duration == o.duration)
    flags_n, flags_o = n.flags.fillna(""), o.flags.fillna("")
    res = dict(sector=sector, n=int(len(new)), n_ok_both=int(ok.sum()), rerun_code=str(new.rerun_code.iloc[0]),
               stored_git_hashes=sorted(map(str, old.git_hash.unique())),
               same_P_D=int(same_pd.sum()), same_P_D_frac=float(same_pd.mean()),
               P_rel_diff_max=float((n.period / o.period - 1).abs().max()),
               max_abs_dS=float((n.S - o.S).abs().max()), n_dS_gt_1e6=int(((n.S - o.S).abs() > 1e-6).sum()),
               max_abs_dt0=float((n.t0 - o.t0).abs().max()),
               flags_changed=int((flags_n != flags_o).sum()), passed_changed=int((n.passed != o.passed).sum()),
               bit_identical_all=int((same_pd & (n.S == o.S) & (flags_n == flags_o) & (n.passed == o.passed)
                                      & (n.t0 == o.t0) & (n.depth == o.depth)).sum()))
    ch = n.index[~same_pd | ((n.S - o.S).abs() > 1e-6) | (flags_n != flags_o) | (n.passed != o.passed)]
    res["changed"] = [dict(tic=int(t), P=[float(o.period[t]), float(n.period[t])], D=[float(o.duration[t]), float(n.duration[t])],
                           S=[float(o.S[t]), float(n.S[t])], flags=[flags_o[t], flags_n[t]],
                           passed=[int(o.passed[t]), int(n.passed[t])]) for t in ch]
    Path(out).write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "changed"}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run"); r.add_argument("--code-root", required=True); r.add_argument("--sector", type=int, required=True)
    r.add_argument("--tics", required=True); r.add_argument("--workers", type=int, default=2); r.add_argument("--out", required=True)
    d = sp.add_parser("diff"); d.add_argument("file"); d.add_argument("--sector", type=int, required=True); d.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "run":
        run(a.code_root, a.sector, a.tics, a.workers, a.out)
    else:
        diff(a.file, a.sector, a.out)

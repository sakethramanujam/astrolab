"""Reconstruct the uncommitted edits behind each '-dirty' v1.0.0 git stamp in a sector (disclosure input).

v1.0.0 stamps git_hash per row at write time, '-dirty' when `git diff -- exoplanet_search/exosearch` is non-empty.
The code that produced a row is the code its search process loaded at start (the per-star pool persists), so each
row is attributed to a process start (from the service log, passed in as --process-starts) and the tree at that
moment is reconstructed from the commit sequence: the edits uncommitted at time t are taken to be the exosearch/
files changed between HEAD(t) and the next commit touching exosearch/ (they are what got committed next).
Edits are classed 'pipeline' if they touch the per-star search/vetting import closure of run_search.py, else
'other' (analysis/deep-vetting/reliability code, docs, site).
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
# per-star search/vetting path of run_search.py (import closure; deepvet only contributes THRESH/REJECT to the hash)
PIPELINE = {"run_search.py", "exosearch/catalogs.py", "exosearch/config.py", "exosearch/plots.py",
            "exosearch/store.py", "exosearch/targets.py", "exosearch/vet.py", "exosearch/task.py",
            "exosearch/bls.py", "exosearch/lightcurve.py", "exosearch/version.py"}
HASH_ONLY = {"exosearch/deepvet.py"}


def git(*a):
    return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def commits():
    rows = [l.split(" ", 1) for l in git("log", "--format=%H %cI", "--first-parent", "-n", "400").splitlines()]
    c = pd.DataFrame(rows, columns=["sha", "t"])
    c["t"] = pd.to_datetime(c.t, utc=True)
    return c.sort_values("t").reset_index(drop=True)


def files_between(a, b):
    return sorted(f.removeprefix("exoplanet_search/") for f in
                  git("diff", "--name-only", a, b, "--", ".").splitlines())


def next_exo_commit(c, sha):
    i = c.index[c.sha.str.startswith(sha)][0]
    for j in range(i + 1, len(c)):
        if any(f.startswith("exosearch/") for f in files_between(c.sha[j - 1], c.sha[j])):
            return c.sha[j], c.t[j]
    return None, None


def classify(files):
    exo = [f for f in files if f.startswith("exosearch/") or f == "run_search.py"]
    pipe = sorted(set(exo) & PIPELINE)
    return pipe, sorted(set(exo) - PIPELINE)


def main(sector, process_starts, out, queue):
    from . import store
    d = store.frame(store.connect(), sector)
    d = d[d.pipeline_version == "1.0.0"].copy()
    d["pa"] = pd.to_datetime(d.processed_at, utc=True)
    starts = sorted(pd.Timestamp(s).tz_convert("UTC") for s in process_starts)
    d["proc_start"] = pd.Series(pd.NaT, index=d.index, dtype="datetime64[ns, UTC]")
    for s in starts:
        d.loc[d.pa >= s, "proc_start"] = s
    c = commits()
    procs = []
    for s, g in d.groupby("proc_start"):
        head = c[c.t <= s].iloc[-1]
        # tree state when the pool loaded its code: HEAD(t) plus the dirty flag of the first row stamped after start
        # if that row still carries HEAD(t) (a clean stamp there means nothing in exosearch/ was uncommitted)
        first = g.sort_values("pa").git_hash.iloc[0]
        if str(first).removesuffix("-dirty") != head.sha[:10]:
            raise SystemExit(f"first row after {s} is stamped {first}, not HEAD {head.sha[:10]}; cannot reconstruct")
        dirty_at_load = str(first).endswith("-dirty")
        nxt, nt = next_exo_commit(c, head.sha)
        files = files_between(head.sha, nxt) if (dirty_at_load and nxt) else []
        pipe, other = classify(files)
        procs.append(dict(process_start=s.isoformat(), head_at_start=head.sha[:10], rows=int(len(g)),
                          dirty_rows_stamped=int(g.git_hash.str.endswith("-dirty").sum()),
                          first_stamp=str(first), tree_dirty_at_load=bool(dirty_at_load), uncommitted_became=nxt[:10] if files else None,
                          uncommitted_committed_at=nt.isoformat() if files else None,
                          pipeline_files=pipe, other_files=other,
                          rows_flagged=int(len(g)) if pipe else 0))
        d.loc[g.index, "flag"] = bool(pipe)
    per_hash = []
    for h, g in d[d.git_hash.str.endswith("-dirty", na=False)].groupby("git_hash"):
        base = h.removesuffix("-dirty")
        nxt, nt = next_exo_commit(c, base)
        pipe, other = classify(files_between(base, nxt) if nxt else [])
        per_hash.append(dict(git_hash=h, rows=int(len(g)), first=g.pa.min().isoformat(), last=g.pa.max().isoformat(),
                             uncommitted_became=nxt[:10] if nxt else None, pipeline_files=pipe, other_files=other,
                             loaded_by_process=sorted({p.isoformat() for p in g.proc_start}),
                             rows_flagged=int(g.flag.sum())))
    per_hash.sort(key=lambda r: r["first"])
    flagged = d[d.flag == True]
    summary = dict(sector=sector, rows_v100=int(len(d)), rows_dirty_stamped=int(d.git_hash.str.endswith("-dirty").sum()),
                   n_dirty_hashes=len(per_hash), processes=procs, per_hash=per_hash, rows_flagged=int(len(flagged)),
                   rerun_queue=Path(queue).name if queue else None, rerun_launched=False)
    Path(out).write_text(json.dumps(summary, indent=1) + "\n")
    if queue:
        Path(queue).write_text("# v1.0.0 S%d rows produced by a process that loaded uncommitted search/vetting edits;\n"
                               "# queued for re-run with the committed code, NOT launched (see %s)\n" % (sector, Path(out).name)
                               + "\n".join(str(int(t)) for t in sorted(flagged.tic)) + "\n")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=105)
    ap.add_argument("--process-starts", nargs="+", required=True, help="ISO times the search pool (re)started")
    ap.add_argument("--out", required=True)
    ap.add_argument("--queue", default=None)
    a = ap.parse_args()
    s = main(a.sector, a.process_starts, a.out, a.queue)
    print(json.dumps({k: v for k, v in s.items() if k != "per_hash"}, indent=1))

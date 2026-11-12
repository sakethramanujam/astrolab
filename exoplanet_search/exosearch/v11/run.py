"""v1.1 runs: build per-host manifests on the coordinator, run them on each host, compare hosts for the gate.

  # coordinator (h1): manifests for the production populations, split by TIC across hosts
  python -m exosearch.v11.run build --sectors 70,106 --weights h1:0.2,h2:0.7,h3:0.1 --out DIR
  # coordinator: the 3-host gate manifest (configs/host_gate_tics.txt, real search)
  python -m exosearch.v11.run build-gate --out DIR
  # every host (clean checkout at the pinned commit, worker venv from requirements.lock, nice 19)
  python -m exosearch.v11.run work --manifest DIR/manifest_h2.jsonl.gz --host-id h2 --workers 12 --out OUTDIR
  # gate comparison
  python -m exosearch.v11.run compare gate_h1 gate_h2 gate_h3 --summary docs/host_gate/summary_v1.1.json

Host assignment: u = sha256("astrolab-v1.1-shard:<tic>") / 2^256 against the cumulative weights; a TIC's real, null and
injection tasks are in one job, so they always run on one host. `work` refuses to start on a dirty tree, on a commit
other than the manifest's pinned commit, or (production) while configs/v1.1.json has no beta option.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"          # before numpy is imported: one thread per worker (reproducible sums, no oversubscription)

import argparse
import gzip
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import settings

SHARD_SALT = "astrolab-v1.1-shard"
SCRAMBLE_REALIZATIONS = (0, 1, 2)


def shard_u(tic):
    return int.from_bytes(hashlib.sha256(f"{SHARD_SALT}:{int(tic)}".encode()).digest(), "big") / 2 ** 256


def assign(tic, weights):
    u, acc = shard_u(tic), 0.0
    tot = sum(w for _, w in weights)
    for h, w in weights:
        acc += w / tot
        if u < acc:
            return h
    return weights[-1][0]


def _nan(x):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else x


def star_params(sector, tics, log=print):
    from .. import stellar
    tp = stellar.tic_params(list(tics), log=log)
    gd = stellar.gaia_dr3(tp.GAIA.tolist(), log=log)
    out = {}
    for tic in tics:
        tic = int(tic)
        if tic not in tp.index:
            out[tic] = {}
            continue
        x = tp.loc[tic]
        rho, erho = stellar.tic_rho(x.to_dict())
        g = str(x.GAIA)
        gr = gd.loc[int(g)] if g.isdigit() and int(g) in gd.index else None
        R = float(x.rad) if x.rad == x.rad else None
        M = float(x.mass) if (x.mass == x.mass and x.mass > 0) else (
            float(10 ** (x.logg - 4.438) * R ** 2) if (x.logg == x.logg and R) else None)
        out[tic] = dict(rho=_nan(rho), rho_err=_nan(erho), teff=_nan(float(x.Teff)) if x.Teff == x.Teff else None,
                        logg=float(x.logg) if x.logg == x.logg else None, rad_tic=R, mass_tic=M,
                        ruwe=None if gr is None or pd.isna(gr.ruwe) else float(gr.ruwe),
                        nss=None if gr is None or pd.isna(gr.non_single_star) else int(gr.non_single_star))
    return out


def host_rm(j):
    """(R*, M*) for an injection host, as exosearch.inject.build_jobs: TIC radius (else the FITS header radius),
    TIC mass (else from logg and R*). None when unusable."""
    R = j.get("rad_tic") or j.get("radius")
    M = j.get("mass_tic")
    if not M and j.get("logg") is not None and R:
        M = float(10 ** (j["logg"] - 4.438) * R ** 2)
    if not (R and R > 0 and M and 0.05 < M < 10):
        return None
    return float(R), float(M)


def base_job(sector, tic, fname, radius, tmag, sp):
    from .. import holdout
    return dict(key=f"v11:{sector}:{tic}", tic=int(tic), sector=int(sector), fname=fname, holdout=bool(holdout.is_holdout(tic)),
                radius=_nan(float(radius)) if radius is not None else None, tmag=_nan(float(tmag)) if tmag is not None else None,
                **sp)


def build(sectors, weights, out, production=True, log=print, tics_filter=None, tasks_override=None, write=True):
    """Production manifests: per searched star real + flipped + 3 scrambled; injection hosts get their runs."""
    from .. import inject, nulltest, store, targets
    cfg = settings.cfg()
    cfg_inj = json.loads((settings.ROOT / "configs" / "injection_v1.1.json").read_text())
    prior_all = settings.beta_prior(cfg) or {}
    con = store.connect()
    jobs = []
    for s in sectors:
        df = store.frame(con, s)
        df = df[df.status == "ok"].drop_duplicates("tic")
        if tics_filter is not None:
            df = df[df.tic.isin(tics_filter)]
        man = targets.fetch_manifest(s).set_index("tic")
        df = df[df.tic.isin(man.index)]
        sp = star_params(s, df.tic.astype(int).tolist(), log=log)
        inj_hosts = set()
        if production and tasks_override is None:
            inj_hosts = set(inject.select_hosts(s, cfg_inj, con).tic.astype(int))
        runs = cfg_inj["runs"].get(str(s), [])

        scramble_exclude = settings.exclusion_tics(s, cfg)      # hard error if registered but missing / no tic column

        for r in df.itertuples():
            tic = int(r.tic)
            p = sp.get(tic, {})
            j = base_job(s, tic, man.loc[tic].fname, r.radius, r.tmag, p)
            if tasks_override is not None:
                tasks = [dict(t) for t in tasks_override]
            else:
                tasks = [dict(tag="real", mode="real"), dict(tag="invert", mode="invert", seed=0)]
                if tic not in scramble_exclude:
                    tasks += [dict(tag=f"scramble_r{k}", mode="scramble", realization=k, seed=nulltest.seed_for("scramble", s, tic, k))
                              for k in SCRAMBLE_REALIZATIONS]
                rm = host_rm(j) if tic in inj_hosts else None
                if rm:
                    j.update(rstar=rm[0], mstar=rm[1])
                    tasks += [dict(tag=f"inject_run{k}", mode="inject", run=k, seed=inject.seed_for(cfg_inj["seed"], s, k, tic))
                              for k in runs]
            j["tasks"] = tasks
            j["beta_prior"] = prior_all.get(str(s))
            j["host"] = assign(tic, weights)
            jobs.append(j)
    if write:
        write_manifests(jobs, out, cfg, cfg_inj, weights)
    return jobs


def write_manifests(jobs, out, cfg, cfg_inj, weights, name="manifest"):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    head = dict(pinned_commit=settings.commit(), cfg=cfg, cfg_inj=cfg_inj, weights=weights, created=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    summ = {}
    for h in sorted({j["host"] for j in jobs}):
        mine = [j for j in jobs if j["host"] == h]
        with gzip.open(out / f"{name}_{h}.jsonl.gz", "wt") as fh:
            fh.write(json.dumps(dict(header=head)) + "\n")
            for j in mine:
                fh.write(json.dumps(j, default=float) + "\n")
        summ[h] = dict(stars=len(mine), tasks=sum(len(j["tasks"]) for j in mine))
    (out / f"{name}_summary.json").write_text(json.dumps(dict(head, hosts=summ), indent=1))
    print(json.dumps(summ))


def read_manifest(path):
    with gzip.open(path, "rt") as fh:
        head = json.loads(fh.readline())["header"]
        jobs = [json.loads(l) for l in fh if l.strip()]
    return head, jobs


def build_gate(out, tics_file=None, sector=70):
    tics = [int(l.split()[0]) for l in Path(tics_file or settings.ROOT / "configs" / "host_gate_tics.txt").read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    jobs = build([sector], [("all", 1.0)], out, production=False, tics_filter=set(tics),
                 tasks_override=[dict(tag="real", mode="real")])
    return jobs


def primary_option(c, beta_option=None, dry_run=False):
    """Beta option used as the statistic. Production: only the configured one, and never unset."""
    if dry_run:
        return beta_option or "a"
    if c.get("beta_option") not in settings.BETA_OPTIONS:
        raise SystemExit("refusing to start production: configs/v1.1.json has no beta option "
                         "(set only by the dated registration applying docs/preregistration/beta_treatment.md)")
    if beta_option and beta_option != c["beta_option"]:
        raise SystemExit("production uses the configured beta option only")
    return c["beta_option"]


def work(manifest, host_id, workers, out, beta_option=None, dry_run=False, allow_other_host=False, limit=0):
    from .. import batch, worker_guard
    from . import pipeline
    hid = settings.host_id(host_id)
    worker_guard.require_clean(False)
    head, jobs = read_manifest(manifest)
    c = settings.cfg()
    if settings.commit() != head["pinned_commit"]:
        raise SystemExit("refusing to start: checkout is not at the manifest's pinned commit")
    if json.dumps(c, sort_keys=True) != json.dumps(head["cfg"], sort_keys=True):
        raise SystemExit("refusing to start: configs/v1.1.json differs from the manifest")
    primary = primary_option(c, beta_option, dry_run)
    mine = [j for j in jobs if j["host"] in (hid, "all")] if not allow_other_host else jobs
    if limit:
        mine = mine[:limit]
    prov = settings.provenance(hid, primary, production=not dry_run)
    for j in mine:
        j.update(cfg=head["cfg"], cfg_inj=head["cfg_inj"], primary=primary, prov=prov)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"provenance_{hid}.json").write_text(json.dumps(dict(prov, manifest=Path(manifest).name, n_jobs=len(mine),
                                                                 workers=workers, dry_run=dry_run), indent=1))
    t0 = time.time()
    with batch.claim(workers):
        df = batch.pmap(pipeline.run_job, mine, out / "parts", workers, key="job", task_key="key", part_size=200)
    wall = time.time() - t0
    ok = df[df.status == "ok"] if len(df) else df
    meta = dict(prov, n_rows=int(len(df)), n_ok=int(len(ok)), wall_s=wall, workers=workers,
                runtime_s_median=float(ok.runtime_s.median()) if len(ok) else None,
                tasks_per_hour=float(3600 * workers / ok.runtime_s.mean()) if len(ok) else None)
    (out / f"done_{hid}.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta))


GATE_COLS = ["S", "period", "P_grid", "duration", "P_grid_halfwidth"]


def load_rows(d):
    from .. import interim
    return interim.read_parts(Path(d) / "parts" / "part-*.parquet")


def compare(dirs, summary=None, tol=1e-6):
    ds = [load_rows(d) for d in dirs]
    ref = ds[0][ds[0].status == "ok"].set_index("key")
    hosts = []
    cols = GATE_COLS + [f"{k}_opt_{o}" for o in settings.BETA_OPTIONS for k in ("S", "period", "duration", "P_grid_halfwidth")]
    for d in ds:
        x = d[d.status == "ok"].set_index("key")
        common = ref.index.intersection(x.index)
        a, b = ref.loc[common], x.loc[common]
        fails = set()
        rep = {}
        for ccol in cols:
            if ccol not in a:
                continue
            va, vb = a[ccol].astype(float), b[ccol].astype(float)
            if ccol.startswith("S"):
                diff = (va - vb).abs()
                bad = common[(diff > tol) | (va.isna() != vb.isna())]
                rep[ccol] = dict(max_abs_diff=float(diff.max()) if len(diff) else 0.0, n_fail=int(len(bad)))
            else:
                same = (va == vb) | (va.isna() & vb.isna())
                bad = common[~same]
                rep[ccol] = dict(n_differ=int(len(bad)), max_abs_diff=float((va - vb).abs().max()) if len(va) else 0.0)
            fails |= set(bad)
        h = d.host_id.dropna().iloc[0] if "host_id" in d and d.host_id.notna().any() else "?"
        hosts.append(dict(host_id=h, arch=d.arch.dropna().iloc[0] if "arch" in d else None,
                          lib_versions=json.loads(d.lib_versions.dropna().iloc[0]) if "lib_versions" in d else None,
                          commit=d.commit.dropna().iloc[0] if "commit" in d else None,
                          tree=d.tree.dropna().iloc[0] if "tree" in d and d.tree.notna().any() else None,
                          n_rows=int(len(d)), n_ok=int(len(x)), n_err=int((d.status != "ok").sum()), n_compared=int(len(common)),
                          columns=rep, n_fail=len(fails), fail_keys=sorted(fails)[:20],
                          tie_band_gt1=int((x.n_tie_band > 1).sum()) if "n_tie_band" in x else None,
                          plateau_gt1=int((x.plateau_width > 1).sum()) if "plateau_width" in x else None,
                          prewhitened_stars=int((x.pw_prewhiten_n > 0).sum()) if "pw_prewhiten_n" in x else None,
                          runtime_s_median=float(x.runtime_s.median()) if len(x) else None,
                          passes=bool(len(fails) == 0 and (d.status != "ok").sum() == 0 and len(common) == len(ref))))
    trees = {h["tree"] for h in hosts}
    commits = {h["commit"] for h in hosts}
    tree_match = len(trees) == 1
    res = dict(gate="3-host gate on the v1.1 code (exosearch.v11, real search, sector 70)", tolerance_S=tol,
               exact=["period", "P_grid", "duration", "P_grid_halfwidth"], reference=hosts[0]["host_id"],
               n_stars=int(len(ref)), hosts=hosts, commits=sorted(map(str, commits)), trees=sorted(map(str, trees)),
               tree_match=tree_match, all_pass=bool(tree_match and all(h["passes"] for h in hosts)))
    if summary:
        Path(summary).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))
    return res


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    b = sp.add_parser("build"); b.add_argument("--sectors", default="70,106"); b.add_argument("--weights", required=True)
    b.add_argument("--out", required=True)
    g = sp.add_parser("build-gate"); g.add_argument("--out", required=True); g.add_argument("--tics-file")
    w = sp.add_parser("work"); w.add_argument("--manifest", required=True); w.add_argument("--host-id")
    w.add_argument("--workers", type=int, default=4); w.add_argument("--out", required=True)
    w.add_argument("--beta-option", choices=settings.BETA_OPTIONS)
    w.add_argument("--dry-run", action="store_true", help="gate/validation: not production rows; beta option from --beta-option")
    w.add_argument("--limit", type=int, default=0)
    c = sp.add_parser("compare"); c.add_argument("dirs", nargs="+"); c.add_argument("--summary")
    a = ap.parse_args()
    if a.cmd == "build":
        wts = [(h, float(x)) for h, x in (p.split(":") for p in a.weights.split(","))]
        build([int(s) for s in a.sectors.split(",")], wts, a.out)
    elif a.cmd == "build-gate":
        build_gate(a.out, a.tics_file)
    elif a.cmd == "work":
        work(a.manifest, a.host_id, a.workers, a.out, a.beta_option, a.dry_run, limit=a.limit)
    else:
        compare(a.dirs, a.summary)


if __name__ == "__main__":
    main()

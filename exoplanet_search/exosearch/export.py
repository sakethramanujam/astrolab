"""Export the public summary the GitHub Pages site reads: site/data/progress.json (+ CSV)
and small diagnostic PNGs for listed candidates.  Contains no hostnames or local paths.

    python -m exosearch.export [--site ../site] [--max-candidates 40]
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__, catalogs, config, plots, sectorqueue, store, targets

SITE = Path(__file__).resolve().parents[1] / "site"
SUBMISSIONS = Path(__file__).resolve().parents[1] / "submissions.json"   # manual: {"<tic>": {"where": "...", "ref": "...", "date": "..."}}
MAX_PNG = 300 * 1024


def _known(r):
    return bool(r.toi_period_match or r.ctoi_period_match or r.planet_period_match)


def _label(r):
    if r.toi_period_match and r.known_toi:
        return f"TOI-{r.known_toi}"
    if r.planet_period_match and r.known_planet:
        return r.known_planet
    if r.ctoi_period_match and r.known_ctoi:
        return f"CTOI {r.known_ctoi}"
    return None


def _status(r, subs):
    if str(r.tic) in subs:
        return "submitted"
    if _known(r):
        return "known"
    if r.passed == 1:
        return "needs_review"
    return "rejected"


def _f(x, nd=6):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return None
    return round(float(x), nd)


def _plot(r, site, dpi=80):
    rel = Path("plots") / f"s{int(r.sector):04d}" / f"TIC{int(r.tic)}_s{int(r.sector):04d}.png"
    out = site / rel
    if out.exists() and out.stat().st_size <= MAX_PNG:
        return str(rel)
    try:
        man = targets.fetch_manifest(int(r.sector)).set_index("tic")
        m = man.loc[int(r.tic)]
        path = targets.download_one(m.fname, m.url, config.lc_dir(int(r.sector)))  # re-fetch if pruned
        for d in (dpi, 64, 50):
            plots.diagnostic(path, r._asdict(), out, dpi=d)
            if out.stat().st_size <= MAX_PNG:
                return str(rel)
        out.unlink(missing_ok=True)
    except Exception as exc:  # a missing plot must not break the export
        print(f"plot failed for TIC {r.tic}: {exc}")
    return None


def build(site=SITE, max_candidates=80, max_rejected=12, plots_on=True):
    con = store.connect()
    df = store.frame(con)
    runs = store.runs(con)
    subs = json.loads(SUBMISSIONS.read_text()) if SUBMISSIONS.exists() else {}
    ok = df[df.status == "ok"].copy()
    try:
        ps_ok = bool(catalogs.load_catalogs().get("ps_ok"))
    except Exception:
        ps_ok = False

    q = {r["sector"]: r for r in sectorqueue.status()}
    sectors = []
    for s, g in ok.groupby("sector"):
        toi = g[g.known_toi.notna()]
        sectors.append(dict(
            sector=int(s), n_targets=q.get(s, {}).get("n_targets"), n_searched=int(len(g)),
            n_errors=int((df[(df.sector == s)].status != "ok").sum()), status=q.get(s, {}).get("status"),
            n_passed=int(g.passed.sum()), toi_hosts=int(len(toi)),
            toi_recovered=int(((toi.toi_period_match == 1) & (toi.passed == 1)).sum()),
            toi_period_found=int((toi.toi_period_match == 1).sum()),
            needs_review=int(((g.passed == 1) & ~g.apply(_known, axis=1)).sum()),
            median_star_s=_f(g.runtime_s.median(), 3)))

    ok["state"] = ok.apply(lambda r: _status(r, subs), axis=1)
    listed = pd.concat([
        ok[ok.state.isin(["submitted", "needs_review"])].sort_values("sde", ascending=False),   # always all of these
        ok[ok.state == "known"].sort_values("sde", ascending=False).head(max_candidates),
        # strongest rejects that still look physically plausible (skip broken photometry)
        ok[(ok.state == "rejected") & (ok.depth < 0.05) & (ok.depth > 0)]
        .sort_values("sde", ascending=False).head(max_rejected)])
    cands = []
    for r in listed.itertuples():
        cands.append(dict(
            tic=int(r.tic), sector=int(r.sector), status=r.state, known=_label(r),
            toi_host=r.known_toi if r.known_toi and not r.toi_period_match else None,
            tmag=_f(r.tmag, 2), period=_f(r.period, 5), duration_h=_f(r.duration * 24, 2),
            depth_ppm=_f(r.depth * 1e6, 0), sde=_f(r.sde, 2), snr=_f(r.snr, 1),
            odd_even_sigma=_f(r.odd_even_sigma, 2), secondary_sigma=_f(r.secondary_sigma, 2),
            n_transits=int(r.n_transits or 0), rp_earth=_f(r.rp_earth, 2), flags=r.flags or "",
            submission=subs.get(str(r.tic)),
            plot=_plot(r, site) if plots_on else None,
            exofop=f"https://exofop.ipac.caltech.edu/tess/target.php?id={int(r.tic)}"))

    toi = ok[ok.known_toi.notna()]
    sde_hist, sde_edges = np.histogram(ok.sde.clip(0, 40), bins=40, range=(0, 40))
    p_edges = np.geomspace(0.5, 14, 41)
    p_hist, _ = np.histogram(ok.period, bins=p_edges)
    rng = np.random.default_rng(0)
    sky = ok[["ra", "dec"]].dropna()
    sky = sky.iloc[rng.permutation(len(sky))[:4000]] if len(sky) > 4000 else sky
    flags = {}
    for f in ok["flags"].fillna(""):
        for x in filter(None, f.split(",")):
            flags[x] = flags.get(x, 0) + 1

    out = dict(
        generated_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        pipeline_version=__version__,
        totals=dict(sectors_with_results=len(sectors), sectors_available=len(q),
                    sectors_done=sum(1 for r in q.values() if r["status"] == "done"),
                    stars_searched=int(len(ok)), errors=int((df.status != "ok").sum()),
                    passed=int(ok.passed.sum()), toi_hosts=int(len(toi)),
                    toi_recovered=int(((toi.toi_period_match == 1) & (toi.passed == 1)).sum()),
                    toi_period_found=int((toi.toi_period_match == 1).sum()),
                    needs_review=int((ok.state == "needs_review").sum()),
                    submitted=int((ok.state == "submitted").sum()),
                    planet_archive_crossmatch=ps_ok),
        queue_next=sectorqueue.next_sectors(3),
        catalogs=catalogs.catalog_status(),
        sectors=sectors,
        runs=[dict(sector=int(r.sector), started_at=r.started_at, finished_at=r.finished_at,
                   n_searched=int(r.n_searched or 0), n_errors=int(r.n_errors or 0),
                   search_wall_s=_f(r.search_wall_s, 1), median_star_s=_f(r.median_star_s, 3),
                   workers=None if pd.isna(r.workers) else int(r.workers), backend=r.backend)
              for r in runs.itertuples()],
        candidates=cands,
        charts=dict(sde_hist=dict(edges=sde_edges.tolist(), counts=sde_hist.tolist()),
                    period_hist=dict(edges=[round(x, 3) for x in p_edges], counts=p_hist.tolist()),
                    flags=flags,
                    sky=[[round(a, 2), round(b, 2)] for a, b in sky.itertuples(index=False)]),
    )
    keep = {c["plot"] for c in cands if c["plot"]}
    for png in (site / "plots").glob("s*/*.png"):        # drop plots no longer listed
        if str(png.relative_to(site)) not in keep:
            png.unlink()
    (site / "data").mkdir(parents=True, exist_ok=True)
    (site / "data" / "progress.json").write_text(json.dumps(out, indent=1))
    pd.DataFrame(cands).drop(columns=["submission"]).to_csv(site / "data" / "candidates.csv", index=False)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=Path, default=SITE)
    ap.add_argument("--max-candidates", type=int, default=80)
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args()
    out = build(a.site, a.max_candidates, plots_on=not a.no_plots)
    t = out["totals"]
    print(f"exported: {t['stars_searched']} stars, {len(out['candidates'])} listed, "
          f"{t['toi_recovered']}/{t['toi_hosts']} TOI hosts recovered -> {a.site / 'data' / 'progress.json'}")


if __name__ == "__main__":
    main()

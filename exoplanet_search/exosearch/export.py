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

from . import __version__, catalogs, config, deepvet, plots, sectorqueue, store, targets

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


def _deep(r):
    f = getattr(r, "deep_flags", None)
    return set(filter(None, f.split(","))) if isinstance(f, str) else set()


def _status(r, subs):
    if str(r.tic) in subs:
        return "submitted"
    if _known(r):
        return "known"
    if r.passed == 1:
        return "rejected" if _deep(r) & deepvet.REJECT else "needs_review"
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


LC_VERSION = 2


def _bin(x, y, edges):
    idx = np.digitize(x, edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    sums = np.bincount(idx[ok], weights=y[ok], minlength=len(edges) - 1)
    cnt = np.bincount(idx[ok], minlength=len(edges) - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt > 0, sums / cnt, np.nan)


def _lcjson(r, site):
    """Small light-curve JSON for the interactive candidate view (folded windows + a coarse time series)."""
    rel = Path("data") / "lc" / f"s{int(r.sector):04d}" / f"TIC{int(r.tic)}_s{int(r.sector):04d}.json"
    out = site / rel
    if out.exists():
        try:
            if json.loads(out.read_text()).get("v") == LC_VERSION:
                return str(rel)
        except ValueError:
            pass
    try:
        from .task import prepare
        man = targets.fetch_manifest(int(r.sector)).set_index("tic")
        m = man.loc[int(r.tic)]
        path = targets.download_one(m.fname, m.url, config.lc_dir(int(r.sector)))
        t, f, _ = prepare(path)["binned"]
        good = np.isfinite(t) & np.isfinite(f)
        t, f = t[good], f[good]
        P, t0, dur = float(r.period), float(r.t0), float(r.duration)
        ph = ((t - t0 + 0.5 * P) % P) - 0.5 * P                      # days from mid-transit
        ph2 = ((t - t0) % P) - 0.5 * P                               # days from phase 0.5
        epoch = np.round((t - t0) / P).astype(int)
        win = min(max(3 * dur, 0.125), 0.25 * P)
        rng = np.random.default_rng(int(r.tic) % 2**32)

        def cap(mask, n=700):                                        # keep the JSON small (~15 kB)
            i = np.flatnonzero(mask)
            return np.sort(rng.choice(i, n, replace=False)) if len(i) > n else i
        a, b = cap(np.abs(ph) < win), cap(np.abs(ph2) < win, 400)
        ppm = (f - 1.0) * 1e6
        edges = np.linspace(-0.5, 0.5, 201)
        tb_edges = np.arange(t.min(), t.max() + 1 / 24, 1 / 24)        # 1-hour bins
        ts = _bin(t, ppm, tb_edges)
        tc = 0.5 * (tb_edges[1:] + tb_edges[:-1])
        keep = np.isfinite(ts)
        r2 = lambda x, n: [None if not np.isfinite(v) else round(float(v), n) for v in x]
        d = dict(v=LC_VERSION, tic=int(r.tic), sector=int(r.sector), period=P, duration_h=dur * 24,
                 depth_ppm=float(r.depth) * 1e6,
                 t0_offset=round(float(t0 - t.min()), 5),
                 near=dict(h=r2(ph[a] * 24, 3), ppm=r2(ppm[a], 0), e=epoch[a].tolist()),
                 sec=dict(h=r2(ph2[b] * 24, 3), ppm=r2(ppm[b], 0)),
                 fold=dict(phase=r2(0.5 * (edges[1:] + edges[:-1]), 4), ppm=r2(_bin(ph / P, ppm, edges), 1)),
                 series=dict(t=r2(tc[keep] - t.min(), 4), ppm=r2(ts[keep], 0)))
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(d, separators=(",", ":")))
        return str(rel)
    except Exception as exc:
        print(f"lc json failed for TIC {r.tic}: {exc}")
    return None


def build(site=SITE, max_candidates=80, max_rejected=12, plots_on=True):
    con = store.connect()
    df = store.frame(con)
    runs = store.runs(con)
    subs = json.loads(SUBMISSIONS.read_text()) if SUBMISSIONS.exists() else {}
    dv = deepvet.frame(deepvet.connect())
    if len(dv):
        df = df.merge(dv.drop(columns=["updated_at"]), on=["tic", "sector"], how="left")
    else:
        for k in deepvet.COLS:
            if k not in ("tic", "sector", "updated_at"):
                df[k] = None
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
        ok[(ok.state == "rejected") & (ok.passed == 1)].sort_values("sde", ascending=False),  # failed deep vetting
        ok[(ok.state == "rejected") & (ok.passed == 0) & (ok.depth < 0.05) & (ok.depth > 0)]
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
            deep_vetted=bool(isinstance(r.deep_flags, str)), deep_flags=r.deep_flags if isinstance(r.deep_flags, str) else None,
            centroid_offset_px=_f(r.centroid_offset_px, 2), centroid_offset_sigma=_f(r.centroid_offset_sigma, 1),
            gaia_n_can_mimic=None if r.gaia_n_can_mimic != r.gaia_n_can_mimic or r.gaia_n_can_mimic is None else int(r.gaia_n_can_mimic),
            gaia_mimic_min_sep=_f(r.gaia_mimic_min_sep, 1), eb_match=r.eb_match if isinstance(r.eb_match, str) else None,
            v_ratio=_f(r.v_ratio, 2), sec_scan_sigma=_f(r.sec_scan_sigma, 1), rp_rjup=_f(r.rp_rjup, 2),
            plot=_plot(r, site) if plots_on else None,
            lc=_lcjson(r, site) if plots_on else None,
            ra=_f(r.ra, 5), dec=_f(r.dec, 5), radius_rsun=_f(r.radius, 3), teff=_f(r.teff, 0),
            t0_btjd=_f(r.t0, 5), sec_scan_phase=_f(r.sec_scan_phase, 3), diff_peak_dist_px=_f(r.diff_peak_dist_px, 2),
            gaia_n_neighbours=None if r.gaia_n_neighbours != r.gaia_n_neighbours or r.gaia_n_neighbours is None else int(r.gaia_n_neighbours),
            gaia_mimic_brightest_dm=_f(r.gaia_mimic_brightest_dm, 2),
            exofop=f"https://exofop.ipac.caltech.edu/tess/target.php?id={int(r.tic)}"))

    # deep-vetting funnel over signals that passed the per-star screens with no catalogue match
    rev = ok[(ok.passed == 1) & ~ok.apply(_known, axis=1) & ~ok.tic.astype(str).isin(subs)]
    vet_ = rev[rev.deep_flags.apply(lambda x: isinstance(x, str))]
    fl = vet_.apply(_deep, axis=1) if len(vet_) else pd.Series(dtype=object)
    funnel = dict(
        unmatched_passed=int(len(rev)), deep_vetted=int(len(vet_)),
        survive_blend=int(sum(not (x & deepvet.BLEND) for x in fl)),
        survive_eb=int(sum(not (x & deepvet.EB) for x in fl)),
        survive_reject_checks=int(sum(not (x & deepvet.REJECT) for x in fl)),
        survive_all=int(sum(not (x & (deepvet.REJECT | deepvet.BLEND)) for x in fl)),
        by_flag={k: int(sum(k in x for x in fl)) for k in sorted(deepvet.REJECT | deepvet.BLEND |
                                                                {"neighbour_could_mimic", "centroid_inconclusive", "centroid_shift_minor"})})

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
        queue=[dict(sector=s_, status=r_["status"], n_targets=r_["n_targets"], n_searched=r_["n_searched"])
               for s_, r_ in sorted(q.items())],
        deep_funnel=funnel,
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
    keep = {c["plot"] for c in cands if c["plot"]} | {c["lc"] for c in cands if c.get("lc")}
    for png in [*(site / "plots").glob("s*/*.png"), *(site / "data" / "lc").glob("s*/*.json")]:   # drop files no longer listed
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

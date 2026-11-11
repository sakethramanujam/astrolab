"""Why did the search miss known TOIs? One bucket per missed host, first matching rule wins:

1 few_transits   : < 2 transits of the TOI ephemeris with >= half their cadences in the raw sector data
2 long_duration  : biweight window (0.75 d) < 3 x TOI duration -> the detrending eats the transit
3 masked         : >= 2 transits in the raw cadences but < 2 survive the quality mask + gap-edge trim
                   (momentum dumps, scattered light, downlink edges)
4 below_threshold: S measured at the TOI ephemeris < the detection threshold S*
5 outcompeted    : S at the TOI ephemeris >= S* but another peak (other planet, variability,
                   systematics) had higher S and the search keeps one peak per star
6 other          : e.g. period beyond half the baseline, ephemeris drift
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import batch, bls, config, labels, lightcurve, store, targets, task, vet

WINDOW = 0.75


def n_covered(t, P, t0, dur, cad):
    ph = ((t - t0 + 0.5 * P) % P) - 0.5 * P
    inn = np.abs(ph) < dur / 2
    if not inn.any():
        return 0
    ep = np.round((t[inn] - t0) / P).astype(int)
    _, c = np.unique(ep, return_counts=True)
    return int(np.sum(c >= 0.5 * dur / cad))


def s_at(tf, ff, P, t0, dur):
    ph = ((tf - t0 + 0.5 * P) % P) - 0.5 * P
    inn, out = np.abs(ph) < dur / 2, np.abs(ph) > dur
    if inn.sum() < 3 or out.sum() < 50:
        return 0.0, 0.0
    depth = float(np.mean(ff[out]) - np.mean(ff[inn]))
    cd = bls.noise_profile(tf, ff, [dur])[2][float(dur)][0]
    ntr = bls.transit_count(tf, P, t0, dur, float(np.median(np.diff(tf))))
    return (depth * np.sqrt(ntr) / cd if cd and cd > 0 else 0.0), depth


def diagnose(sector=70, matched_col="toi_period_match"):
    con = store.connect()
    df = store.frame(con, sector)
    miss = df[df.known_toi.notna() & (df[matched_col] == 0)]
    toi = labels.toi_table()
    man = targets.fetch_manifest(sector).set_index("tic")
    rows = []
    for r in miss.itertuples():
        ts = toi[(toi.tic == r.tic) & toi.period.notna() & (toi.period > 0)].copy()
        p = batch.lc_path(sector, man.loc[r.tic].fname)
        t, f, e, q, meta = lightcurve.load_spoc_raw(p)
        good = lightcurve.good_mask(f, e, q)
        med = np.nanmedian(f[good])
        cad = float(np.median(np.diff(t)))
        d = task.prepare_arrays(t[good], f[good] / med, e[good] / med, meta)
        tf, ff, _ = d["flat"]
        s_star, prov = vet.s_star(sector, meta.get("camera"))
        best = None
        for x in ts.itertuples():
            P, dur = float(x.period), float(x.dur_h) / 24 if x.dur_h == x.dur_h else 0.1
            t0 = float(x.epoch_bjd) - 2457000.0
            n_raw = n_covered(t, P, t0, dur, cad)
            n_good = n_covered(t[good], P, t0, dur, cad)
            n_kept = n_covered(tf, P, t0, dur, cad)
            S, dep = s_at(tf, ff, P, t0, dur)
            if n_raw < 2:
                b = "few_transits"
            elif WINDOW < 3 * dur:
                b = "long_duration"
            elif n_kept < 2:
                b = "masked"
            elif S < s_star:
                b = "below_threshold"
            elif (getattr(r, "S", None) or 0) > S or not r.period or abs(r.period / P - 1) > 0.01:
                b = "outcompeted"
            else:
                b = "other"
            cand = dict(tic=int(r.tic), toi=x.toi, disp=x.disp, period=P, dur_h=dur * 24, depth_ppm=x.depth_ppm,
                        n_raw=n_raw, n_quality=n_good, n_kept=n_kept, S_at_ephemeris=S, s_star=s_star,
                        best_peak_P=r.period, best_peak_sde=r.sde, bucket=b)
            # a host is explained by its most detectable TOI
            if best is None or S > best["S_at_ephemeris"]:
                best = cand
        if best is None:      # single-transit TOIs carry no period
            best = dict(tic=int(r.tic), toi=";".join(str(x) for x in toi[toi.tic == r.tic].toi), period=np.nan,
                        best_peak_P=r.period, best_peak_sde=r.sde, bucket="few_transits", note="TOI has no period")
        rows.append(best)
    return pd.DataFrame(rows)


def main():
    import argparse
    from .version import stamp
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=70)
    a = ap.parse_args()
    out = config.DATA_ROOT / "products" / f"v{stamp()['pipeline_version']}" / "missed_tois"
    out.mkdir(parents=True, exist_ok=True)
    d = diagnose(a.sector)
    d.to_csv(out / f"missed_tois_s{a.sector:04d}.csv", index=False)
    counts = d.bucket.value_counts().to_dict()
    (out / f"missed_tois_s{a.sector:04d}.json").write_text(json.dumps(dict(stamp(), n=len(d), buckets=counts), indent=1))
    print(len(d), counts)


if __name__ == "__main__":
    main()

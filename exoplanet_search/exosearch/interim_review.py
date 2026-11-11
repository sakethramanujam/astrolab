"""Interim-layer review diagnostics (analysis layer over frozen v1.0.0; nothing here changes the search).

  python -m exosearch.interim_review            # -> products/v<ver>/analysis_layer/review/{review.json, *.png}

* GPD starting point u: mean-residual-life (mean excess) and parameter-stability (xi, modified scale sigma - xi*u)
  plots of the vetted inverted-null max S; u chosen by a fixed rule (below) and the S* spread over the stable range
  quoted as a systematic.
* Per-camera S*: each camera solved at its own alpha = 1/N_cam and, for comparison, at the sector alpha = 1/N_sector.
* Star-quality cut sensitivity: S* and vetted completeness for beta_max cuts 2.5, 3, 4 and no cut.
* Known planets (TOI hosts) among dropped vs kept stars: share and recovery rate.
* Completeness inside the radius/depth range the too_deep and size caps allow; what too_deep costs at high SNR.
* Tail-fit FAP of a given S (per star, and expected post-vetting false alarms per sector).
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from . import interim, labels, quality, store, vet

OUT = interim.PROD / "analysis_layer" / "review"
U_GRID = np.round(np.arange(5.5, 10.01, 0.25), 2)
MIN_K = 50                     # tail points needed at the chosen u
SE_BOOT = 300
CUTS = (2.5, 3.0, 4.0, None)
MAX_DEPTH = vet.THRESH["max_depth"]
MAX_RP = vet.THRESH["max_rp_earth"]


def kept(null, q, cut):
    m = null.merge(q[["tic", "sector", "beta_max", "status"]].rename(columns={"status": "q_status"}), on=["tic", "sector"], how="left")
    ok = m.q_status.eq("ok")
    m["drop"] = ~ok if cut is None else ~(ok & (m.beta_max <= cut))
    return m


def pop_kept(q, sector, cut):
    x = q[(q.sector == sector) & (q.status == "ok")]
    return int(len(x)) if cut is None else int((x.beta_max <= cut).sum())


def fit_u(Sv_all, n, u):
    y = Sv_all[Sv_all > u] - u
    if len(y) < 5:
        return None
    xi, _, sig = stats.genpareto.fit(y, floc=0)
    return dict(u=float(u), k=int(len(y)), xi=float(xi), sigma=float(sig), sigma_star=float(sig - xi * u), mean_excess=float(y.mean()))


def s_star_at(Sv_all, n, p, u):
    r = interim.gpd_threshold(Sv_all, n, p, u)
    return r["s_star"]


def tail_rate(Sv_all, n, u, s):
    """Per-null-star probability of a vetted peak >= s (GPD above u, empirical below)."""
    k = int((Sv_all > u).sum())
    if s <= u:
        return float((Sv_all >= s).sum() / n)
    xi, _, sig = stats.genpareto.fit(Sv_all[Sv_all > u] - u, floc=0)
    return float(k / n * stats.genpareto.sf(s - u, xi, 0, sig))


def choose_u(diag):
    """Lowest u with >= MIN_K tail points whose xi lies inside +-1 bootstrap SE of every higher-u xi that still has
    >= 30 points (parameter stability), and whose mean-excess curve beyond it is linear within errors."""
    ok = [d for d in diag if d["k"] >= 30]
    for i, d in enumerate(ok):
        if d["k"] < MIN_K:
            continue
        if all(abs(e["xi"] - d["xi"]) <= max(d["xi_se"], e["xi_se"]) for e in ok[i + 1:]):
            return d["u"]
    return ok[0]["u"] if ok else interim.U_DEFAULT


def gpd_diagnostics(x, n_stars, sector):
    Sv = x.S.values[x.vetted.values]
    n = len(x)
    p = 1.0 / n_stars
    rng = np.random.default_rng(interim.SEED)
    diag = []
    for u in U_GRID:
        f = fit_u(Sv, n, u)
        if f is None:
            continue
        xs = []
        for _ in range(SE_BOOT):
            b = rng.choice(Sv[Sv > u], size=int((Sv > u).sum()), replace=True)
            xs.append(stats.genpareto.fit(b - u, floc=0)[0])
        y = Sv[Sv > u] - u
        f.update(xi_se=float(np.std(xs)), me_se=float(y.std(ddof=1) / np.sqrt(len(y))) if len(y) > 1 else None,
                 s_star=s_star_at(Sv, n, p, u))
        diag.append(f)
    u0 = choose_u(diag)
    stable = [d for d in diag if d["u"] >= u0 and d["k"] >= 30]
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), dpi=100)
    U = [d["u"] for d in diag]
    ax[0].errorbar(U, [d["mean_excess"] for d in diag], [2 * (d["me_se"] or 0) for d in diag], fmt="o", ms=3)
    ax[0].set_xlabel("threshold u"); ax[0].set_ylabel("mean excess  E[S - u | S > u]"); ax[0].set_title("mean residual life", fontsize=9)
    ax[1].errorbar(U, [d["xi"] for d in diag], [2 * d["xi_se"] for d in diag], fmt="o", ms=3)
    ax[1].set_xlabel("threshold u"); ax[1].set_ylabel("shape xi (2 SE)"); ax[1].set_title("parameter stability: shape", fontsize=9)
    ax[2].plot(U, [d["s_star"] for d in diag], "o-", ms=3)
    ax[2].set_xlabel("threshold u"); ax[2].set_ylabel("S* at alpha = 1/N_stars"); ax[2].set_title("S* vs starting point", fontsize=9)
    for a in ax:
        a.axvline(u0, color="#e45756", ls=":", lw=1)
        a2 = a.twiny(); a2.set_xlim(a.get_xlim()); a2.set_xticks(U[::4]); a2.set_xticklabels([d["k"] for d in diag][::4], fontsize=7)
        a2.set_xlabel("tail points k", fontsize=7)
    fig.suptitle(f"Sector {sector}: vetted inverted-null max S, GPD threshold choice (chosen u = {u0:g})", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / f"gpd_threshold_choice_s{sector:04d}.png"); plt.close(fig)
    return diag, u0, stable


def completeness(inj, q, cut, s_star, sector=70):
    d = inj[(inj.status == "ok") & (inj.sector == sector)].merge(q[["tic", "sector", "beta_max", "status"]].rename(columns={"status": "q_status"}),
                                                                 on=["tic", "sector"], how="left")
    ok = d.q_status.eq("ok")
    d = d[ok if cut is None else ok & (d.beta_max <= cut)].copy()
    failed = d.failed_checks.fillna("").str.split(",").map(lambda x: set(x) - {"", "low_s"})
    d["det"] = d.matched.astype(bool) & (d.S >= s_star)
    d["vet"] = d.det & (failed.map(len) == 0)
    d["failed"] = failed
    d["allowed"] = (d.depth_inj < MAX_DEPTH) & (d.rp_inj < MAX_RP)
    return d


def comp_table(d, edges=(0, 10, 15, 25, 50, 100, 1e9)):
    b = pd.cut(d.snr_exp, list(edges), right=False)
    g = d.groupby(b, observed=True).agg(n=("det", "size"), det=("det", "mean"), vet=("vet", "mean"))
    return {str(k): dict(n=int(v.n), detected=round(float(v.det), 3), vetted=round(float(v.vet), 3)) for k, v in g.iterrows()}


def toi_split(q, sectors=(70, 106)):
    t = labels.toi_table()
    t = t[t.disp.isin(["CP", "KP", "PC", "APC"])]
    con = store.connect()
    out = {}
    for s in sectors:
        df = store.frame(con, s)
        df = df[df.status == "ok"].merge(q[q.sector == s][["tic", "beta_max", "status"]].rename(columns={"status": "q_status"}), on="tic", how="left")
        df["dropped"] = ~(df.q_status.eq("ok") & (df.beta_max <= quality.rule()["threshold"]))
        df["toi_host"] = df.tic.isin(t.tic)
        r = {}
        for lab, g in (("dropped", df[df.dropped]), ("kept", df[~df.dropped])):
            h = g[g.toi_host]
            rec = int((h.toi_period_match == 1).sum()) if "toi_period_match" in h else None
            r[lab] = dict(stars=int(len(g)), toi_hosts=int(len(h)), toi_host_rate=round(len(h) / max(len(g), 1), 4),
                          recovered=rec, recovery_rate=round(rec / len(h), 3) if len(h) and rec is not None else None,
                          beta_max_median_hosts=float(h.beta_max.median()) if len(h) else None)
        a, b = r["dropped"], r["kept"]
        tab = [[a["toi_hosts"], a["stars"] - a["toi_hosts"]], [b["toi_hosts"], b["stars"] - b["toi_hosts"]]]
        r["fisher_overrep"] = dict(odds_ratio=float(stats.fisher_exact(tab)[0]), p=float(stats.fisher_exact(tab)[1]))
        out[str(s)] = r
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rule = quality.rule()
    q = interim.read_parts(interim.PROD / "analysis_layer" / "star_quality_s*_parts" / "part-*.parquet")
    null = interim.load_null("invert")
    inj = interim.read_parts(interim.PROD / "injections" / "*" / "part-*.parquet")
    res = dict(quality_rule=rule["rule"], n_null_rows=int(len(null)))
    s = 70
    # --- GPD starting point at the registered cut
    m = kept(null, q, rule["threshold"])
    x = m[(m.sector == s) & ~m["drop"] & ~m.holdout.astype(bool)]
    N = pop_kept(q, s, rule["threshold"])
    diag, u0, stable = gpd_diagnostics(x, N, s)
    Sv = x.S.values[x.vetted.values]
    n = len(x)
    p = 1.0 / N
    point = s_star_at(Sv, n, p, u0)
    rng = np.random.default_rng(interim.SEED)
    vv = np.where(x.vetted.values, x.S.values, -np.inf)
    boots = []
    for _ in range(2000):
        b = vv[rng.integers(0, n, n)]
        r_ = s_star_at(b[np.isfinite(b)], n, p, u0)
        if r_ is not None and np.isfinite(r_):
            boots.append(r_)
    sys_vals = [d["s_star"] for d in stable if d["s_star"]]
    res["gpd"] = dict(sector=s, n_null=n, n_stars=N, alpha=p, u_rule=f"lowest u with k>={MIN_K} whose xi is within 1 bootstrap SE of xi at every higher u with k>=30",
                      u=u0, k=int((Sv > u0).sum()), s_star=point, ci95_stat=[float(np.quantile(boots, .025)), float(np.quantile(boots, .975))],
                      sys_range=[float(min(sys_vals)), float(max(sys_vals))] if sys_vals else None,
                      sys_halfwidth=float((max(sys_vals) - min(sys_vals)) / 2) if sys_vals else None,
                      stable_us=[d["u"] for d in stable], diag=diag)
    # --- per camera: own alpha vs sector alpha
    cams = {}
    for cam, xc in x.groupby("camera"):
        Svc = xc.S.values[xc.vetted.values]
        frac = len(xc) / len(x)
        Ncam = N * frac
        k = int((Svc > u0).sum())
        cams[str(int(cam))] = dict(n_null=int(len(xc)), k_tail=k, N_cam_est=round(Ncam),
                                   s_star_own_alpha=s_star_at(Svc, len(xc), 1 / Ncam, u0) if k >= 15 else None,
                                   s_star_sector_alpha=s_star_at(Svc, len(xc), 1 / N, u0) if k >= 15 else None,
                                   vetted_rate_above_u=round(k / len(xc), 4))
    res["per_camera"] = cams
    # --- quality-cut sensitivity
    sens = {}
    for cut in CUTS:
        mc = kept(null, q, cut)
        xc = mc[(mc.sector == s) & ~mc["drop"] & ~mc.holdout.astype(bool)]
        Nc = pop_kept(q, s, cut)
        Svc = xc.S.values[xc.vetted.values]
        sst = s_star_at(Svc, len(xc), 1 / Nc, u0)
        dcomp = completeness(inj[~inj.holdout.astype(bool)], q, cut, sst, s)
        al = dcomp[dcomp.allowed]
        sens[str(cut)] = dict(stars_kept=Nc, stars_dropped=int((q.sector == s).sum() - Nc), n_null=int(len(xc)), k_tail=int((Svc > u0).sum()),
                              s_star=sst, inj_kept=int(len(dcomp)),
                              vetted_completeness_allowed_snr_ge_25=round(float(al[al.snr_exp >= 25].vet.mean()), 3),
                              detected_allowed_snr_ge_25=round(float(al[al.snr_exp >= 25].det.mean()), 3),
                              vetted_completeness_allowed_snr_15_25=round(float(al[(al.snr_exp >= 15) & (al.snr_exp < 25)].vet.mean()), 3))
    res["quality_sensitivity"] = sens
    # --- completeness inside the allowed radius/depth range at the registered cut
    d = completeness(inj, q, rule["threshold"], point, s)
    res["completeness_allowed"] = {sp: comp_table(d[d.allowed & (d.holdout.astype(bool) == (sp == "holdout"))]) for sp in ("tuning", "holdout")}
    res["completeness_all_radii"] = {sp: comp_table(d[d.holdout.astype(bool) == (sp == "holdout")]) for sp in ("tuning", "holdout")}
    res["allowed_fraction"] = round(float(d.allowed.mean()), 3)
    hi = d[d.det & ~d.vet & (d.snr_exp >= 50)]
    cnt = hi.failed.explode().value_counts()
    res["high_snr_failures"] = dict(n=int(len(hi)), by_check={k: int(v) for k, v in cnt.items()},
                                    too_deep_any=round(float(hi.failed.map(lambda f: "too_deep" in f).mean()), 3),
                                    too_deep_only=round(float(hi.failed.map(lambda f: f == {"too_deep"}).mean()), 3),
                                    outside_allowed=round(float((~hi.allowed).mean()), 3),
                                    odd_even_any=round(float(hi.failed.map(lambda f: "odd_even" in f).mean()), 3))
    # --- known planets in dropped vs kept
    res["toi_dropped_vs_kept"] = toi_split(q)
    # --- tail FAP for a given S
    res["tail_fap"] = {str(S): dict(per_star=tail_rate(Sv, n, u0, S), expected_false_alarms_per_sector=tail_rate(Sv, n, u0, S) * N)
                       for S in (14.8, 16.6, 24.1)}
    (OUT / "review.json").write_text(json.dumps(res, indent=1, default=float))
    small = {k: v for k, v in res.items() if k != "gpd"}
    small["gpd"] = {k: v for k, v in res["gpd"].items() if k != "diag"}
    print(json.dumps(small, indent=1, default=float))


if __name__ == "__main__":
    main()

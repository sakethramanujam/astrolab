"""Provisional reliability summaries for the site (site/data/reliability.json) and the archive.

Only binned counts and empirical quantiles live here. The statistical modelling (gamma-CDF completeness
with bootstrap, GPD/Gumbel tail fits for S*, ROC, vetting effectiveness E, DR25-style reliability) belongs
to analysis/; if analysis/make_reliability_json.py exists it is run instead and its output wins.

  python -m exosearch.reliability            # -> site/data/reliability.json (+ products/v*/summary/)
"""
from __future__ import annotations

import glob
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import config, deepvet, holdout, nulltest, store, vet
from .version import PIPELINE_VERSION, config_hash, stamp


def _clean(o):
    """JSON-safe copy: NaN/inf -> null, numpy scalars -> Python (browsers reject bare NaN)."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.generic):
        o = o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o

ROOT = Path(__file__).resolve().parents[1]
SITE_JSON = ROOT / "site" / "data" / "reliability.json"
PROD = config.DATA_ROOT / "products" / f"v{PIPELINE_VERSION}"
P_EDGES = np.geomspace(0.5, 13.5, 9)
R_EDGES = np.geomspace(0.5, 20.0, 9)
TMAG_EDGES = [0, 9, 10.5, 12, 13.5, 20]
SNR_EDGES = np.r_[np.arange(0, 30, 1.5), 30, 40, 60, 100]
S_EDGES = np.round(np.arange(0, 40.25, 0.25), 2)
GROUPS = {"s70": [70], "s106_107": [106, 107]}
ALL_CHECKS = (["few_transits", "odd_even", "secondary", "too_deep", "too_long", "variable", "rp_too_large", "no_dip",
               "common_period"] + ["centroid_offset", "dip_off_target", "known_eb", "v_shape", "secondary_scan",
               "rp_gt_2rjup", "density_mismatch", "unresolved_neighbour", "centroid_shift_minor", "ruwe_high", "gaia_nss"])


def _layer():
    """Interim analysis layer (quality cut + vetted-null S*), or None before exosearch.interim has run."""
    from . import interim, quality
    if not interim.THR_FILE.exists():
        return None
    thr = json.loads(interim.THR_FILE.read_text())
    ss = interim.s_star_by_sector()
    ref = ss.get(70) or (next(iter(ss.values())) if ss else None)
    return dict(thr=thr, q=quality.load(), s_star=ss, ref=ref, interim=interim)


MARGINAL_EFA = 0.1      # expected post-vetting false alarms per sector at the signal's S above which it is called marginal


def s_star_for(layer, sector):
    """(S*, borrowed_from) for a sector; sectors without their own null borrow the sector-70 value."""
    v = layer["s_star"].get(int(sector))
    return (v, None) if v else (layer["ref"], 70)


def _read(pattern):
    fs = sorted(glob.glob(str(pattern)))
    return pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def wilson(k, n, z=1.0):
    if n == 0:
        return [None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 4), round(c + h, 4)]


def group_of(sector):
    return next((g for g, ss in GROUPS.items() if int(sector) in ss), f"s{int(sector)}")


# ---------------------------------------------------------------- injections
def completeness(inj):
    if inj.empty:
        return None
    d = inj[inj.status == "ok"].copy()
    d["group"] = d.sector.map(group_of)
    d["split"] = np.where(d.holdout, "holdout", "tuning")
    d["pb"] = np.digitize(d.P_inj, P_EDGES) - 1
    d["rb"] = np.digitize(d.rp_inj, R_EDGES) - 1
    d["tb"] = np.digitize(d.tmag, TMAG_EDGES[1:-1])
    d["sb"] = np.digitize(d.snr_exp.fillna(0), SNR_EDGES) - 1
    out = {"p_edges": P_EDGES.round(3).tolist(), "r_edges": R_EDGES.round(3).tolist(), "tmag_edges": TMAG_EDGES,
           "snr_edges": SNR_EDGES.tolist(), "n": int(len(d)), "maps": {}, "snr_curves": {}, "summary": {}}
    for (g, sp), x in list(d.groupby(["group", "split"])) + [((g, "all"), x) for g, x in d.groupby("group")]:
        for tb in ["all"] + list(range(len(TMAG_EDGES) - 1)):
            y = x if tb == "all" else x[x.tb == tb]
            n = np.zeros((8, 8), int); det = n.copy(); pas = n.copy()
            ok = (y.pb >= 0) & (y.pb < 8) & (y.rb >= 0) & (y.rb < 8)
            for (pb, rb), z in y[ok].groupby(["pb", "rb"]):
                n[rb, pb], det[rb, pb], pas[rb, pb] = len(z), int(z.detected.sum()), int(z.passed_all.sum())
            out["maps"][f"{g}|{sp}|{tb}"] = dict(n=n.tolist(), detected=det.tolist(), passed=pas.tolist())
        c = x[(x.sb >= 0) & (x.sb < len(SNR_EDGES) - 1)].groupby("sb").agg(n=("detected", "size"), det=("detected", "sum"),
                                                                         pas=("passed_all", "sum"))
        out["snr_curves"][f"{g}|{sp}"] = {int(k): [int(v.n), int(v.det), int(v.pas)] for k, v in c.iterrows()}
        out["summary"][f"{g}|{sp}"] = dict(n=int(len(x)), detected=int(x.detected.sum()), passed=int(x.passed_all.sum()),
                                           matched=int(x.matched.sum()),
                                           aliases={str(k): int(v) for k, v in x.alias.value_counts().items()},
                                           fail_counts={c: int(x[x.detected][f"fail_{c}"].fillna(False).sum())
                                                        for c in nulltest.SCREEN_CHECKS + ["common_period"] + nulltest.DEEP_CHECKS
                                                        if f"fail_{c}" in x})
    return out


# ---------------------------------------------------------------- null tests
def nulls(n_df, layer=None):
    if n_df.empty:
        return None
    d = n_df[n_df.status == "ok"].copy()
    if layer is not None:                  # same star-quality cut as the search population and injections
        d = layer["interim"].with_quality(d, layer["q"])
        d = d[~d["drop"]].copy()
    fl = d["flags"].fillna("").str.split(",").map(lambda x: set(x) - {"", "low_s"})
    dfl = d.deep_flags.fillna("").str.split(",").map(set)
    # passes every check except the S threshold (deep checks are logged for S >= S_FLOOR)
    d["vet_ok"] = (fl.map(len) == 0) & ~d.get("fail_common_period", pd.Series(False, index=d.index)).fillna(False) \
        & dfl.map(lambda x: not (x & deepvet.REJECT)) & d.checks_logged.fillna(False)
    out = {"s_edges": S_EDGES.tolist(), "s_floor": nulltest.S_FLOOR, "hist": {}, "n_stars": {}, "threshold": {}}
    for (s, m), x in d.groupby(["sector", "mode"]):
        key = f"{int(s)}|{m}"
        h_all, _ = np.histogram(x.S.clip(0, 39.99), S_EDGES)
        h_vet, _ = np.histogram(x[x.vet_ok].S.clip(0, 39.99), S_EDGES)
        out["hist"][key] = dict(all=h_all.tolist(), vetted=h_vet.tolist(), n_real=int(x.realization.nunique()))
        out["n_stars"][key] = int(len(x) / max(x.realization.nunique(), 1))
    if layer is not None:
        for sec, r in layer["thr"]["sectors"].items():
            out["threshold"][sec] = dict(n=r["n_null"], n_stars=r["n_stars"], s_star=r["s_star"], ci68=r["ci68"], ci95=r["ci95"],
                                         k_tail=r["k_tail"], u=r["u"], xi=r["xi"], max_vetted=r["max_vetted"],
                                         sys_halfwidth=r.get("sys_halfwidth"), cameras="pooled")
        out["threshold_kind"] = "interim GPD fit, vetted inverted null, tuning stars passing the quality cut (provisional)"
        # flipped vs scrambled side by side (same quality cut, vetting, u rule); scrambled pending until its run lands
        sc = d[d["mode"] == "scramble"] if "mode" in d else d.iloc[:0]
        side = {}
        for sec, r in layer["thr"]["sectors"].items():
            side[sec] = dict(invert=dict(s_star=r["s_star"], ci95=r["ci95"], n=r["n_null"]), scramble=None)
            xs = sc[sc.sector == int(sec)]
            if len(xs):
                try:
                    from . import interim_review as rv
                    xq = layer["interim"].with_quality(xs.assign(vetted=_vetted(xs)), layer["q"])
                    xq = xq[~xq["drop"] & ~xq.holdout.astype(bool)]
                    f = layer["interim"].fit_sector(xq, r["n_stars"], u=r["u"], n_boot=500)
                    side[sec]["scramble"] = dict(s_star=f["s_star"], ci95=f["ci95"], n=f["n_null"], realizations=int(xs.realization.nunique()))
                except Exception as exc:
                    side[sec]["scramble"] = dict(error=type(exc).__name__)
        out["threshold_by_null"] = dict(side, note="S* is provisional until both nulls are complete; scramble = null when not yet run")
        return out
    inv = d[d["mode"] == "invert"]
    for (s, cam), x in inv.groupby(["sector", "camera"]):
        n = len(x)
        out["threshold"][f"{int(s)}:{int(cam)}"] = dict(n=n, s_star=round(float(np.quantile(x.S, 1 - 1 / n)), 3),
                                                       max=round(float(x.S.max()), 3))
    for s, x in inv.groupby("sector"):
        n = len(x)
        out["threshold"][f"{int(s)}"] = dict(n=n, s_star=round(float(np.quantile(x.S, 1 - 1 / n)), 3),
                                            max=round(float(x.S.max()), 3))
    return out


def _vetted(x):
    fl = x["flags"].fillna("").str.split(",").map(lambda f: set(f) - {"", "low_s"})
    dfl = x.deep_flags.fillna("").str.split(",").map(set)
    cp = x.fail_common_period.fillna(False).astype(bool) if "fail_common_period" in x else False
    return (fl.map(len) == 0) & ~cp & dfl.map(lambda f: not (f & deepvet.REJECT)) & x.checks_logged.fillna(False).astype(bool)


def fap(S, sector, n_df, layer=None):
    """Empirical FAP. With the analysis layer: post-vetting, on quality-passing tuning stars of the inverted null
    (k = null stars whose vetted max S reaches S). Without it: raw (pre-vetting) fraction."""
    if S is None:
        return None
    if layer is not None:
        x = layer["interim"].with_quality(layer["interim"].load_null("invert"), layer["q"])
        x = x[(x.sector == sector) & ~x["drop"] & ~x.holdout.astype(bool)]
        if x.empty:
            return None
        k = int((x.vetted & (x.S >= S)).sum())
        return dict(k=k, n=int(len(x)), fap=round((k + 1) / (len(x) + 1), 6), post_vetting=True)
    x = n_df[(n_df["mode"] == "invert") & (n_df.sector == sector) & (n_df.status == "ok")]
    if x.empty:
        return None
    k = int((x.S >= S).sum())
    return dict(k=k, n=int(len(x)), fap=round((k + 1) / (len(x) + 1), 6), post_vetting=False)


# ---------------------------------------------------------------- per-check rates (labelled signals)
def check_rates(lab):
    if lab.empty:
        return None
    d = lab[lab.status == "ok"].copy()
    d["group"] = d.sector.map(group_of)
    d["split"] = np.where(d.holdout, "holdout", "tuning")
    fl = d["flags"].fillna("").str.split(",").map(set) | d.deep_flags.fillna("").str.split(",").map(set)
    d["allflags"] = fl
    out = {}
    for (g, sp), x in d.groupby(["group", "split"]):
        rows = {}
        for c in ALL_CHECKS:
            r = {}
            for lab_ in ("planet", "fp"):
                y = x[x.label == lab_]
                if c == "known_eb" and lab_ == "fp":
                    y = y[~y.source.fillna("").str.startswith("TESS EB")]       # circular for catalogue EBs
                k = int(y.allflags.map(lambda s: c in s).sum())
                r[lab_] = dict(k=k, n=int(len(y)), ci68=wilson(k, len(y)))
            r["reject"] = c in deepvet.REJECT or c in nulltest.SCREEN_CHECKS + ["common_period"]
            rows[c] = r
        y = x.copy()
        rej = y.allflags.map(lambda s: bool(s & (deepvet.REJECT | set(nulltest.SCREEN_CHECKS) | {"common_period"})))
        rows["_any_reject"] = {l: dict(k=int(rej[y.label == l].sum()), n=int((y.label == l).sum()),
                                       ci68=wilson(int(rej[y.label == l].sum()), int((y.label == l).sum())))
                               for l in ("planet", "fp")}
        out[f"{g}|{sp}"] = rows
    return out


def build():
    analysis = ROOT / "analysis" / "make_reliability_json.py"
    if analysis.exists():                                   # the analysis PR owns the numbers once it lands
        subprocess.run([sys.executable, str(analysis)], check=True)
        return json.loads(SITE_JSON.read_text())
    inj = _read(PROD / "injections" / "*" / "part-*.parquet")
    nl = _read(PROD / "nulltests" / "*" / "null_*_maxS.parquet")
    if nl.empty:
        nl = _read(PROD / "nulltests" / "*" / "parts" / "part-*.parquet")
        if not nl.empty and "realization" not in nl:
            nl["realization"] = 0
    lab = _read(PROD / "labelled_v1" / "parts" / "part-*.parquet")
    missed = {}
    for f in glob.glob(str(PROD / "missed_tois" / "*.json")):
        missed[Path(f).stem.split("_")[-1]] = json.loads(Path(f).read_text())
    layer = _layer()
    if layer is not None:
        thr = dict(statistic="S = depth_box / (CDPP_D / sqrt(N_tr))", interim=True, provisional=True,
                   status="provisional until both the flipped (inverted) and the scrambled nulls are complete",
                   nulls_complete=dict(invert=(PROD / "nulltests" / "invert" / "null_invert_maxS.parquet").exists(),
                                       scramble=bool(glob.glob(str(PROD / "nulltests" / "scramble_r*" / "null_*_maxS.parquet")))),
                   kind="interim: GPD tail fit to the vetted inverted null (tuning stars passing the quality cut), "
                        "post-vetting false-alarm rate 1/N_stars per sector; bootstrap CI",
                   default=layer["ref"], by_sector={k: dict(s_star=v["s_star"], ci68=v["ci68"], ci95=v["ci95"], n_null=v["n_null"],
                                                            n_stars=v["n_stars"], k_tail=v["k_tail"], u=v.get("u"), u_stable=v.get("u_stable"),
                                                            sys_range=v.get("sys_range"), sys_halfwidth=v.get("sys_halfwidth"),
                                                            cameras="pooled: per-camera tails are too thin to fit (not published)")
                                              for k, v in layer["thr"]["sectors"].items()},
                   borrowed_note="sectors without their own inverted null yet use the sector-70 value",
                   pipeline_floor=json.loads(vet.THRESHOLD_FILE.read_text())["default"],
                   quality_cut=layer["thr"]["quality_cut"], population=layer["thr"]["population"])
        if not inj.empty:
            inj = inj[inj.status == "ok"]
            ss = {int(sec): s_star_for(layer, sec)[0] for sec in inj.sector.unique()}
            inj = layer["interim"].rescore_injections(inj, layer["q"], ss)
    else:
        thr = json.loads(vet.THRESHOLD_FILE.read_text())
    out = dict(stamp(), generated_by="exosearch.reliability (interim analysis layer)" if layer else "exosearch.reliability (provisional binned summaries)",
               holdout=dict(seed=holdout.SEED, fraction=holdout.FRACTION),
               detection_threshold=thr, completeness=completeness(inj), nulls=nulls(nl, layer),
               check_rates=check_rates(lab), missed_tois=missed)
    extra = {"review": PROD / "analysis_layer" / "review" / "review.json",
             "whitenoise_neff": PROD / "calibration" / "whitenoise_neff_s0070.json",
             "ffi_check_318707987": PROD / "analysis_layer" / "ffi_check" / "TIC318707987" / "ffi_check.json"}
    for k, f in extra.items():
        if f.exists():
            j = json.loads(f.read_text())
            if k == "review":
                j["gpd"] = {a: b for a, b in j["gpd"].items() if a != "diag"}
                j.pop("per_camera", None)          # per-camera S* not published: cameras pooled until the samples grow
            if k == "ffi_check_318707987":
                j = {a: j.get(a) for a in ("tess_point_on_silicon", "not_on_silicon", "verdict", "control_at_ref", "control_orbits_at_ref")}
                j["control_blind"] = (json.loads(f.read_text())["sectors"].get("70") or {}).get("blind_top")
            out[k] = j
    cand = PROD / "summary" / "survivors.json"
    if cand.exists():
        out["survivors"] = json.loads(cand.read_text())
    SITE_JSON.write_text(json.dumps(_clean(out), separators=(",", ":"), allow_nan=False, default=str))
    (PROD / "summary").mkdir(parents=True, exist_ok=True)
    (PROD / "summary" / "reliability.json").write_text(SITE_JSON.read_text())
    return out



def survivors(sector=70):
    """Re-vet the v0 deep-vetting survivors: v1 S at their own light curve, density/RUWE flags, empirical FAP."""
    from . import batch, targets, task
    dcon = deepvet.connect()
    df = store.frame(store.connect(), sector)
    dv = deepvet.frame(dcon, sector)
    c = deepvet.candidates(df, only_review=True).merge(dv, on=["tic", "sector"])
    v0_reject = {"centroid_offset", "dip_off_target", "known_eb", "v_shape", "secondary_scan", "rp_gt_2rjup"}
    fl = c.deep_flags.fillna("").str.split(",").map(lambda x: set(x) - {""})
    c = c[fl.map(lambda x: not (x & v0_reject))].copy()
    man = targets.fetch_manifest(sector).set_index("tic")
    nl = _read(PROD / "nulltests" / "invert" / "parts" / "part-*.parquet")
    layer = _layer()
    sstar = s_star_for(layer, sector)[0] if layer else None
    qs = layer["q"][layer["q"].sector == sector].set_index("tic") if layer else None
    tail = None
    if layer:
        from . import interim_review as rv
        tr = layer["thr"]["sectors"].get(str(sector)) or layer["thr"]["sectors"]["70"]
        x = layer["interim"].with_quality(layer["interim"].load_null("invert"), layer["q"])
        x = x[(x.sector == (sector if str(sector) in layer["thr"]["sectors"] else 70)) & ~x["drop"] & ~x.holdout.astype(bool)]
        Sv = x.S.values[x.vetted.values]
        tail = lambda S_: rv.tail_rate(Sv, len(x), tr["u"], S_) * tr["n_stars"]   # expected post-vetting false alarms / sector
    rows = []
    for r in c.itertuples():
        flags = sorted(set((r.deep_flags or "").split(",")) - {""})
        res = task.search_prepared(task.prepare(batch.lc_path(sector, man.loc[r.tic].fname)))
        same = abs(res["period"] / r.period - 1) < 0.01
        rej = sorted(set(flags) & deepvet.REJECT)
        efa = float(tail(res["S"])) if tail and same else None
        word = "marginal" if efa is not None and efa >= MARGINAL_EFA else "survives vetting"
        status = "rejected" if rej else (word if not flags else f"{word}, warnings")
        f = fap(res["S"], sector, nl, layer) if same else None
        qd = bool(qs.loc[r.tic, "drop"]) if qs is not None and r.tic in qs.index else None
        bmax = float(qs.loc[r.tic, "beta_max"]) if qs is not None and r.tic in qs.index else None
        if not rej and qd:
            status = "dropped (star-quality cut)"
        elif not rej and sstar is not None and same and res["S"] < sstar:
            status = "below S*"
        rows.append(dict(tic=int(r.tic), sector=sector, period=float(r.period), rp_earth=None if r.rp_earth != r.rp_earth else float(r.rp_earth),
                         S=float(res["S"]) if same else None, v1_period=float(res["period"]), v1_same_peak=bool(same),
                         beta_d=float(res["beta_d"]), fap=None if f is None else f["fap"], fap_k=None if f is None else f["k"],
                         fap_n=None if f is None else f["n"], rho_transit=r.rho_transit, rho_tic=r.rho_tic, rho_z=r.rho_z,
                         ruwe=r.ruwe, gaia_nss=r.gaia_nss, flags=flags, status=status, s_star=sstar,
                         above_s_star=None if sstar is None or not same else bool(res["S"] >= sstar),
                         beta_max=bmax, quality_drop=qd, expected_false_alarms_per_sector=efa))
    n_null = int(len(nl[(nl.sector == sector)])) if not nl.empty else 0
    out = dict(stamp(), sector=sector, n=len(rows), n_null=n_null,
               note=(f"The {len(rows)} sector-{sector} signals that survived the v0 deep checks, re-vetted with the v1 density and Gaia checks. "
                     f"S is the v1 statistic at the same peak. " + (
                         f"Interim S* = {sstar:.1f} (vetted inverted null). FAP = (k+1)/(N+1), k of N quality-passing tuning flipped light "
                         f"curves whose max S passes every check and reaches that S (post-vetting, partial null). "
                         f"Expected false alarms per sector at S from the GPD tail fit; 'marginal' when >= {MARGINAL_EFA}. "
                         f"No row is a detection or a confirmed planet." if sstar else
                         f"FAP = (k+1)/(N+1) with k of N flipped light curves reaching that S (N={n_null}, before vetting).")),
               rows=sorted(rows, key=lambda x: (not x["status"].startswith(("marginal", "survives")), -(x["S"] or 0))))
    (PROD / "summary").mkdir(parents=True, exist_ok=True)
    (PROD / "summary" / "survivors.json").write_text(json.dumps(_clean(out), indent=1, allow_nan=False, default=str))
    return out


if __name__ == "__main__":
    if "--survivors" in sys.argv:
        o = survivors(70)
        print({x["status"]: sum(1 for y in o["rows"] if y["status"] == x["status"]) for x in o["rows"]})
    o = build()
    print({k: (None if v is None else "ok") for k, v in o.items() if k in ("completeness", "nulls", "check_rates")},
          SITE_JSON.stat().st_size // 1024, "KB")

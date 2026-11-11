"""sigma_P for the candidate tables: half-width of the in-band plateau that contains the stored period
(registered tie-break, docs/preregistration/tie_break.md, amendment 1).

Each candidate star is re-searched with the interim fixed search (`host_gate run --fixed`, x86_64 host, 1 thread
per worker); the stored v1.0.0 period is located on that search's plateau list. If it lies on an in-band plateau
(|P - centre| <= sigma_P), sigma_P is that plateau's half-width; otherwise the stored peak is not in the tie band of
the fixed-grid search and sigma_P is half the local grid step, P^2 * df / 2 (flag sigma_P_in_band = False).

  python -m exosearch.sigma_p lists --outdir /tmp/sigp          # one TIC list per sector from site/data/candidates.csv
  (host_gate run --fixed --sector S --tics-file /tmp/sigp/tics_sS.txt --out /tmp/sigp/plat_sS.csv.gz, per sector)
  python -m exosearch.sigma_p match /tmp/sigp/plat_s*.csv.gz     # -> products/interim/sigma_p_v1.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import config

ROOT = Path(__file__).resolve().parents[1]
OUT = config.DATA_ROOT / "products" / "interim" / "sigma_p_v1.csv"


def lists(outdir, cands=ROOT / "site" / "data" / "candidates.csv"):
    c = pd.read_csv(cands)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    for s, g in c.groupby("sector"):
        Path(outdir, f"tics_s{int(s)}.txt").write_text("\n".join(map(str, sorted(set(g.tic.astype(int))))) + "\n")
        print(f"S{int(s)}: {g.tic.nunique()} stars")


def locate(P, plateaus, grid_df):
    for lo, hi, w, sp in plateaus:
        if abs(P - 0.5 * (lo + hi)) <= sp * (1 + 1e-9):
            return dict(sigma_P=sp, sigma_P_in_band=True, plateau_width=int(w), plateau_centre=0.5 * (lo + hi))
    return dict(sigma_P=0.5 * P * P * grid_df, sigma_P_in_band=False, plateau_width=1, plateau_centre=np.nan)


def match(files, out=OUT, stored=None):
    from . import store
    con = store.connect()
    rows = []
    for f in files:
        g = pd.read_csv(f)
        meta = json.loads(Path(str(f) + ".meta.json").read_text()) if Path(str(f) + ".meta.json").exists() else {}
        sector = int(meta.get("sector") or g.get("sector", pd.Series([0])).iloc[0])
        st = store.frame(con, sector).set_index("tic")
        for r in g[g.status == "ok"].itertuples():
            P = float(st.loc[int(r.tic)].period)
            loc = locate(P, json.loads(r.plateaus), float(r.grid_df))
            rows.append(dict(tic=int(r.tic), sector=sector, period_v10=P, **loc,
                             n_plateaus=int(r.n_plateaus), disjoint=bool(r.disjoint)))
    d = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(out, index=False)
    print(f"{len(d)} rows -> {out}; in band {d.sigma_P_in_band.mean():.3f}; plateau>1 {(d.plateau_width > 1).mean():.3f}; "
          f"median sigma_P {d.sigma_P.median():.5f} d")
    return d


def load(path=OUT):
    return pd.read_csv(path) if Path(path).exists() else None


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    l = sp.add_parser("lists"); l.add_argument("--outdir", required=True)
    m = sp.add_parser("match"); m.add_argument("files", nargs="+"); m.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    if a.cmd == "lists":
        lists(a.outdir)
    else:
        match(a.files, a.out)

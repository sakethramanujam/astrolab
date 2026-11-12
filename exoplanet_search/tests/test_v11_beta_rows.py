import gzip
import json

import numpy as np
import pandas as pd
import pytest

from exosearch import interim, quality, stat_options
from exosearch.v11 import beta_rows as br
from exosearch.v11 import beta_select as bs
from exosearch.v11 import registration_checks

TREE = "t" * 40
SHA = "e" * 64


def _rows(tics, task, seed=0, **kw):
    rng = np.random.default_rng(seed)
    n = len(tics)
    S = rng.gumbel(6.0, 0.8, n)
    d = dict(tic=list(tics), sector=70, task=task, status="ok", tree=TREE, commit="c" * 40, config_hash="f9f10beec4f1",
             pipeline_version="1.1.0", production=False, scramble_exclude_sha256_70=SHA, holdout=False,
             S=S, period=rng.uniform(0.5, 13, n), flags="", deep_flags=None, checks_logged=True)
    for o, f in zip("abcd", (1.0, 0.97, 0.99, 0.98)):
        d[f"S_opt_{o}"] = S * f
        d[f"P_grid_opt_{o}"] = d["period"]
    d["P_grid"] = d["period"]
    d.update(kw)
    return pd.DataFrame(d)


def _set(n=40, excluded=(1, 2, 3)):
    tics = list(range(1, n + 1))
    flip = _rows(tics, "invert")
    scr = _rows([t for t in tics if t not in excluded], "scramble_r0", seed=1)
    return flip, scr, {"invert": n, "scramble_r0": n - len(excluded)}


def test_dedupe_identical_duplicates_and_ok_supersedes_error():
    flip, scr, exp = _set()
    dup = flip.iloc[:5].copy()
    err = scr.iloc[:2].copy().assign(status="error", S=np.nan)
    df = pd.concat([flip.assign(_src=0), scr.assign(_src=1), dup.assign(_src=2), err.assign(_src=2)], ignore_index=True)
    x, rep = br.merge_rows(df, exp, tree=TREE)
    assert rep["n_per_task"] == exp and rep["n_duplicate_rows"] == 14 and len(x) == sum(exp.values())
    assert (x.status == "ok").all()


def test_conflicting_duplicates_rejected():
    flip, scr, exp = _set()
    bad = flip.iloc[:1].copy()
    bad["S_opt_c"] += 0.1
    df = pd.concat([flip.assign(_src=0), scr.assign(_src=1), bad.assign(_src=2)], ignore_index=True)
    with pytest.raises(br.RowsError, match="different results"):
        br.merge_rows(df, exp)
    bad = flip.iloc[:1].copy()
    bad["flags"] = "odd_even"                          # vetting inputs must agree too
    df = pd.concat([flip.assign(_src=0), scr.assign(_src=1), bad.assign(_src=2)], ignore_index=True)
    with pytest.raises(br.RowsError, match="different results"):
        br.merge_rows(df, exp)


def test_tree_mismatch_rejected():
    flip, scr, exp = _set()
    scr.loc[scr.index[0], "tree"] = "u" * 40
    df = pd.concat([flip.assign(_src=0), scr.assign(_src=1)], ignore_index=True)
    with pytest.raises(br.RowsError, match="tree"):
        br.merge_rows(df, exp)
    flip, scr, exp = _set()
    df = pd.concat([flip.assign(_src=0), scr.assign(_src=1)], ignore_index=True)
    with pytest.raises(br.RowsError, match="expected"):
        br.merge_rows(df, exp, tree="v" * 40)


def test_missing_and_error_rows_rejected():
    flip, scr, exp = _set()
    df = pd.concat([flip.iloc[1:].assign(_src=0), scr.assign(_src=1)], ignore_index=True)
    with pytest.raises(br.RowsError, match="row counts"):
        br.merge_rows(df, exp)
    df = pd.concat([flip.assign(_src=0), scr.assign(_src=1)], ignore_index=True)
    df.loc[0, "status"] = "error"
    with pytest.raises(br.RowsError, match="without an ok row"):
        br.merge_rows(df, exp)
    # manifest set: one task listed in the manifest has no row (counts alone would not see a swap)
    df = pd.concat([flip.assign(_src=0), scr.assign(_src=1)], ignore_index=True)
    keys = set(zip(df.tic.astype(int), df.task)) - {(40, "invert")} | {(99, "invert")}
    with pytest.raises(br.RowsError, match="1 missing"):
        br.merge_rows(df, exp, manifest_keys=keys)


def test_production_counts_are_registered_values():
    assert br.EXPECTED == {"invert": 12797, "scramble_r0": 12512} and sum(br.EXPECTED.values()) == 25309


def test_manifest_reader(tmp_path):
    p = tmp_path / "m.jsonl.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(json.dumps(dict(header={})) + "\n")
        fh.write(json.dumps(dict(tic=5, tasks=[dict(tag="invert"), dict(tag="scramble_r0")])) + "\n")
        fh.write(json.dumps(dict(tic=6, tasks=[dict(tag="invert")])) + "\n")
    assert br.manifest_tasks([p]) == {(5, "invert"), (5, "scramble_r0"), (6, "invert")}


def test_guard_excluded_scrambled_row_rejected_and_paired_set_clean():
    flip, scr, exp = _set()
    leak = _rows([2], "scramble_r0", seed=3)
    df = pd.concat([flip.assign(_src=0), scr.assign(_src=1), leak.assign(_src=1)], ignore_index=True)
    with pytest.raises(br.RowsError, match="excluded host"):
        br.merge_rows(df, {"invert": 40, "scramble_r0": 38}, exclusion_tics={1, 2, 3})
    with pytest.raises(br.RowsError, match="scramble_exclude_sha256"):
        br.merge_rows(pd.concat([flip.assign(_src=0), scr.assign(_src=1)]), exp, exclusion_sha256="0" * 64)
    # same population: excluded stars drop from the flipped side through the pairing; holdout and quality cut apply
    x, _ = br.merge_rows(pd.concat([flip.assign(_src=0), scr.assign(_src=1)]), exp, exclusion_tics={1, 2, 3},
                         exclusion_sha256=SHA)
    x.loc[x.tic == 10, "holdout"] = True
    q = pd.DataFrame(dict(tic=range(1, 41), sector=70, drop=[t == 11 for t in range(1, 41)], beta_max=1.0))
    ex = pd.DataFrame(dict(tic=[1, 2, 3]))
    xf, xs, rep = br.beta_frames(br.with_vetting(x), q, ex)
    assert list(xf.index) == list(xs.index)
    assert not ({1, 2, 3, 10, 11} & set(xf.index))
    assert rep["paired"] == 35 and rep["excluded_removed_from_flipped"] == 3 and rep["n_alpha"] == 39
    assert rep["guard"]["status"] == "ok"
    # an excluded star that reached both sides stops the registration
    xx = pd.concat([x, _rows([2], "scramble_r0", seed=4)], ignore_index=True)
    with pytest.raises(registration_checks.PopulationGuardError):
        br.beta_frames(br.with_vetting(xx), q, ex)


def test_vetting_matches_interim_definition():
    x = _rows(range(1, 201), "invert")
    x.loc[:59, "period"] = 1.2345                      # one period shared by 60 stars: common_period
    x.loc[60, "flags"] = "odd_even"
    x.loc[61, "flags"] = "low_s"                       # low_s does not count against vetting
    x.loc[62, "deep_flags"] = "v_shape"
    x.loc[63, "checks_logged"] = False
    v = br.with_vetting(x).set_index("tic").vetted
    assert not v.loc[1:60].any() and not v.loc[61] and v.loc[62] and not v.loc[63] and not v.loc[64]
    assert v.loc[65:].all()


def _old_path_inputs(n=1500, seed=7):
    rng = np.random.default_rng(seed)
    tics = np.arange(1000, 1000 + n)

    def frame(s):
        r = np.random.default_rng(s)
        return pd.DataFrame(dict(S_exact=r.gumbel(6.0, 0.9, n), beta_d=r.uniform(0.7, 1.6, n),
                                 duration=r.choice([0.04, 0.08, 0.13, 0.25], n), n_cad=r.integers(9000, 13000, n).astype(float),
                                 vetted=r.uniform(size=n) < 0.6), index=pd.Index(tics, name="tic"))
    pri = {f"{D:.2f}": (1.1, 0.04) for D in (0.04, 0.06, 0.08, 0.10, 0.13, 0.17, 0.21, 0.25)}
    return frame(seed), frame(seed + 1), pri, int(n * 1.2)


@pytest.mark.parametrize("opt", ["a", "c"])
def test_tost_path_unchanged_same_inputs_same_verdict(monkeypatch, opt):
    xf, xs, pri, N = _old_path_inputs()
    monkeypatch.setattr(bs, "N_BOOT", 150)
    monkeypatch.setattr(quality, "load", lambda *a, **k: pd.DataFrame(dict(sector=[70], n_cad=[11000.0])))
    monkeypatch.setattr(stat_options, "prior", lambda q: pri)
    monkeypatch.setattr(interim, "search_population", lambda s, q: {70: dict(kept=N)})
    monkeypatch.setattr(bs, "null_frame", lambda tag, q, s: xf if tag == "F" else xs)
    old = bs.sector_option((70, opt, "F", "S"))
    # the same statistics as v1.1 rows carry them: S_opt_<o> per option
    def v11(x):
        v, _ = stat_options.variants(x.S_exact.values, x.beta_d.values, x.duration.values, x.n_cad.values, pri)
        return x.assign(**{f"S_opt_{o}": v[o] for o in bs.OPTIONS})
    new = br.sector_gate(v11(xf), v11(xs), opt, N)
    assert new == old


def test_decide_is_choose_and_rule_constants_unchanged():
    assert (bs.TOST_MARGIN, bs.N_BOOT, bs.SEED, bs.ORDER, bs.SNR_BIN) == (0.5, 2000, 20261110, ("a", "d", "c", "b"), (15, 25))
    rng = np.random.default_rng(3)
    res = {}
    for s in bs.SECTORS:
        for o, ci in zip(bs.OPTIONS, ([-0.3, 0.4], [-0.6, 0.2], [-0.2, 0.45], [0.1, 0.3])):
            res[(s, o)] = dict(gate_pass=bs.tost(ci), ci95_contains_0=True, null_peaks_moved_gt1pct=int(rng.integers(0, 9)))
    base = rng.uniform(size=400) < 0.5
    ind = {"a": base, "b": base, "c": base | (rng.uniform(size=400) < 0.3), "d": base}
    comp = {o: {"vetted": float(v.mean())} for o, v in ind.items()}
    assert br.decide(res, comp, ind, sectors=bs.SECTORS) == bs.choose(res, comp, ind)
    # completeness is only required when two or more options pass
    assert br.decide(res, None, None, sectors=bs.SECTORS) is None
    none_pass = {k: dict(v, gate_pass=False) for k, v in res.items()}
    d = br.decide(none_pass, None, None, sectors=bs.SECTORS)
    assert d["chosen"] == "b" and d["gate_failed"] and "not launched" in d["reason"]
    one = {k: dict(v, gate_pass=(k[1] == "d")) for k, v in res.items()}
    d = br.decide(one, None, None, sectors=bs.SECTORS)
    assert d["chosen"] == "d" and d["passing"] == ["d"] and not d["gate_failed"]

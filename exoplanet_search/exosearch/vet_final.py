"""One final verdict for every row, whatever its mode (real, null or inject).

One verdict function for every row, whatever its mode (real, null or inject), matching the v1.1 real-row
recipe exactly:

    vet_amend.amended() + common_period flag + checks_logged

finish.py used to score real and null rows with the amended vetter (exosearch.v11.vet_amend.revet) plus
the common-period flag and the checks-logged requirement, but injections with the search-time column
passed_lc_checks (the OLD deepvet.REJECT set, still containing rp_gt_2rjup and density_mismatch) and
without the common-period term.  This module applies exactly the real-row recipe to every row.

No rule, threshold, beta rule or S* is changed here: deepvet.py and vet_amend.py are untouched and the
verdict is assembled only from their outputs and the two row-level terms above.

Light curves: never rebuilt here.  When the deep checks raised 'secondary_scan' the amended verdict needs
the rebuilt-light-curve scan (vet_amend.scan_depth + loo_scan_min) to decide rules (a)/(d).  When the caller
supplies that scan (finish.py rebuilds it for real and null rows exactly as vet_amend.revet_one does),
the rules are decided on it and the rebuilt sigma must reproduce the row's logged sec_scan_sigma
(revet_one's consistency check).  When no scan is supplied for such a row, status is 'needs_lc'.  Callers
bracket those rows: count them as rejected in the lower bound (needs_lc="reject", the default) and exempt
them in the upper bound (needs_lc="exempt", which drops only the 'secondary_scan' rejection; every other
rule still applies).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from . import vet
from .v11 import vet11, vet_amend

__all__ = [
    "Verdict",
    "final_verdict",
    "final_verdict_row",
    "check_scan_sigma",
    "COMMON_PERIOD",
    "CHECKS_NOT_LOGGED",
    "ODD_EVEN",
    "ODD_EVEN_TOO_FEW",
    "ODD_EVEN_NEEDS_LC",
    "STATUS_OK",
    "STATUS_NEEDS_LC",
    "VERDICT_KEYS",
    "RP_CAP_FLAG",
    "MISSING_RSTAR",
    "RP_UNDETERMINED",
    "RP_ERR_INCOMPLETE",
    "DEFAULT_RP_CAP",
    "rp_cap_cfg",
    "rp_cap_enabled",
    "rp_with_error",
    "rp_cap_flags",
    "odd_even_mask_cfg",
    "odd_even_mask_enabled",
    "gap_edges",
    "odd_even_transit_mask",
    "masked_odd_even",
]

COMMON_PERIOD = "common_period"          # vet.flag_common_periods group flag
CHECKS_NOT_LOGGED = "checks_not_logged"  # row['checks_logged'] false: deep checks were never run
ODD_EVEN = vet11.ODD_EVEN
# Distinct flag when < 2 odd or < 2 even transits remain after the odd-even transit mask.
# Does not pass light-curve checks, and is not counted as an `odd_even` depth-comparison failure.
ODD_EVEN_TOO_FEW = vet11.ODD_EVEN_TOO_FEW
ODD_EVEN_NEEDS_LC = "odd_even_needs_lc"  # mask is enabled but no light curve was supplied
STATUS_OK = "ok"
STATUS_NEEDS_LC = "needs_lc"
VERDICT_KEYS = ("passed_final", "first_failing_rule", "all_failing", "status_final",
                "fail_common_period_final", "checks_logged_final")

# --- v1.2 step 3: radius cap on Rp - 2 sigma (settings: vet.rp_cap, default OFF = v1.1 point-estimate cap) ---
RP_CAP_FLAG = "rp_lo2s_gt_cap"          # (Rp - 2 sigma_Rp) > cap; replaces the point-estimate 'rp_gt_30re' when on
MISSING_RSTAR = "missing_rstar"         # R* missing / non-finite / <= 0: Rp cannot be computed; never a cap pass
RP_UNDETERMINED = "rp_undetermined"     # R* valid but depth is non-finite or <= 0: cap cannot be evaluated; not a pass
# v1.2 step 5a: R* valid and depth > 0 (Rp computable) but depth_err or e_rad (TIC R* error, R_Sun) is missing,
# None, non-finite or < 0: sigma_Rp cannot be formed, so the cap is NOT applied to the point estimate (no sigma = 0
# fallback).  Not a pass; distinct from RP_CAP_FLAG, MISSING_RSTAR and RP_UNDETERMINED.
RP_ERR_INCOMPLETE = "rp_err_incomplete"
# cap value = vet_amend.RP_CAP_REARTH (v1.1).  cap_rearth_candidates: values calibration may choose from
# later (a setting only; nothing here selects among them); cap_rearth stays the value in force.
DEFAULT_RP_CAP = dict(enabled=False, cap_rearth=30.0, cap_rearth_candidates=[22.0, 25.0, 30.0])
_RP_CAP_EPS = 1e-9                      # R_Earth; boundary (Rp - 2 sigma == cap) is kept despite float noise

odd_even_mask_cfg = vet11.odd_even_mask_cfg
odd_even_mask_enabled = vet11.odd_even_mask_enabled
gap_edges = vet11.gap_edges
odd_even_transit_mask = vet11.odd_even_transit_mask
masked_odd_even = vet11.masked_odd_even


class Verdict(tuple):
    """final_verdict() result: (passed, first_failing_rule, all_failing, status).

    status is STATUS_OK, or STATUS_NEEDS_LC when the amended verdict needs a rebuilt light curve
    (deep 'secondary_scan', no scan given); passed then follows the needs_lc bound ("reject" or
    "exempt").  Iterable as a plain 4-tuple.
    """

    @property
    def passed(self) -> bool:
        return self[0]

    @property
    def first_failing_rule(self) -> Optional[str]:
        return self[1]

    @property
    def all_failing(self) -> List[str]:
        return list(self[2])

    @property
    def status(self) -> str:
        return self[3]


def _set(s) -> set:
    return vet_amend._set(s)


def _bool(x) -> bool:
    """fillna(False).astype(bool) semantics for a scalar row value (None/NaN -> False)."""
    if x is None:
        return False
    if isinstance(x, float) and x != x:
        return False
    if isinstance(x, str):
        return x.strip().lower() not in ("", "nan", "none", "false", "0")
    return bool(x)


def fail_common_period_of(row: Dict[str, Any], common_period: Optional[Any] = None) -> bool:
    """The common-period term: explicit argument, else the row's own flag/column."""
    if common_period is not None:
        return _bool(common_period)
    if "fail_common_period" in row:
        return _bool(row.get("fail_common_period"))
    return COMMON_PERIOD in _set(row.get("flags"))


def checks_logged_of(row: Dict[str, Any], checks_logged: Optional[Any] = None) -> bool:
    """The checks-logged term: explicit argument, else the row's own value (None/NaN -> False)."""
    if checks_logged is not None:
        return _bool(checks_logged)
    return _bool(row.get("checks_logged"))


def rp_cap_cfg(cfg: Optional[Dict[str, Any]] = None, rp_cap: Optional[Any] = None) -> Dict[str, Any]:
    """Resolve the Rp - 2 sigma radius-cap settings from `cfg["vet"]["rp_cap"]` (and an optional override).

    `rp_cap` may be a bool (enabled) or a dict with `enabled` / `cap_rearth` / `cap_rearth_candidates`.
    Configs without the keys resolve to `DEFAULT_RP_CAP` (rule off: the v1.1 point-estimate cap
    `rp_gt_30re` in vet_amend.amended).  `cap_rearth_candidates` (default [22, 25, 30]) is only exposed for
    calibration; it never changes `cap_rearth`.
    """
    out = dict(DEFAULT_RP_CAP)
    out["cap_rearth_candidates"] = list(DEFAULT_RP_CAP["cap_rearth_candidates"])
    v = cfg.get("vet") if isinstance(cfg, dict) else None
    for src in (v.get("rp_cap") if isinstance(v, dict) else None, rp_cap):
        if isinstance(src, bool):
            out["enabled"] = src
        elif isinstance(src, dict):
            if src.get("enabled") is not None:
                out["enabled"] = bool(src["enabled"])
            if src.get("cap_rearth") is not None:
                out["cap_rearth"] = float(src["cap_rearth"])
            if src.get("cap_rearth_candidates") is not None:
                out["cap_rearth_candidates"] = [float(x) for x in src["cap_rearth_candidates"]]
    return out


def rp_cap_enabled(cfg: Optional[Dict[str, Any]] = None, rp_cap: Optional[Any] = None) -> bool:
    return bool(rp_cap_cfg(cfg, rp_cap)["enabled"])


def _valid_err(x) -> float:
    """x as a float when explicitly present, finite and >= 0 (a literal 0 is valid), else NaN."""
    v = vet_amend._f(x)
    return v if (np.isfinite(v) and v >= 0) else float("nan")


def rp_with_error(row: Dict[str, Any], rstar: Optional[float] = None, e_rad: Optional[float] = None) -> Dict[str, Any]:
    """Rp (R_Earth) and its 1-sigma error for the Rp - 2 sigma radius cap.  The ONE implementation.

    Formulae (R_EARTH_PER_RSUN = exosearch.vet.R_EARTH_PER_RSUN = 109.1):

        Rp              = R* * sqrt(max(depth, 0)) * R_EARTH_PER_RSUN        (R* in R_Sun, depth a fraction)
        (sigma_Rp/Rp)^2 = (sigma_R*/R*)^2 + (depth_err / (2 * depth))^2      depth > 0
        lower2s         = Rp - 2 * sigma_Rp

    sigma_R* is `e_rad` (TIC R* error, R_Sun, the unit of R*).  Both terms must be valid: depth_err and
    e_rad each explicitly present, finite and >= 0.  A literal 0 is valid (that term contributes nothing,
    e.g. e_rad=0 gives the depth-only sigma_Rp = Rp * depth_err / (2 * depth)).  Missing / None / NaN / Inf /
    negative is incomplete: nothing is invented (no 10% R* error, no sigma = 0 fallback).

    R* is the row's `radius` (pipeline KEEP_META) and, only when that is missing/invalid, the job's `radius`
    (`rstar` argument).  e_rad is the row's `e_rad` when valid, else the job's (`e_rad` argument).  A stored
    `rp_earth` / `rp_rjup` is never used: it must not hide a missing R*, and Rp is always computed from the
    same `depth` as sigma_Rp.

    -> dict(status, rp, sigma, lower2s, rstar, e_rad) with status one of
       "ok"               Rp, sigma_Rp, lower2s are finite;
       MISSING_RSTAR      R* missing, non-finite or <= 0 (rp/sigma/lower2s NaN);
       RP_UNDETERMINED    R* valid but depth is non-finite or <= 0 (rp/sigma/lower2s NaN);
       RP_ERR_INCOMPLETE  Rp is computed but depth_err or e_rad is not valid (sigma/lower2s NaN).
    """
    nan = float("nan")
    r = vet_amend._f(row.get("radius"))
    if not (np.isfinite(r) and r > 0):
        r = vet_amend._f(rstar)
    er = _valid_err(row.get("e_rad"))
    if not np.isfinite(er):
        er = _valid_err(e_rad)
    out = dict(status="ok", rp=nan, sigma=nan, lower2s=nan, rstar=r if np.isfinite(r) else nan, e_rad=er)
    if not (np.isfinite(r) and r > 0):
        out["status"] = MISSING_RSTAR
        return out
    depth = vet_amend._f(row.get("depth"))
    if not (np.isfinite(depth) and depth > 0):
        out["status"] = RP_UNDETERMINED
        return out
    rp = float(r * np.sqrt(depth) * vet.R_EARTH_PER_RSUN)
    out["rp"] = rp
    derr = _valid_err(row.get("depth_err"))
    if not (np.isfinite(derr) and np.isfinite(er)):
        out["status"] = RP_ERR_INCOMPLETE
        return out
    sigma = float(rp * np.sqrt((er / r) ** 2 + (derr / (2.0 * depth)) ** 2))
    out.update(sigma=sigma, lower2s=rp - 2.0 * sigma)
    return out


def rp_cap_flags(row: Dict[str, Any], cfg: Optional[Dict[str, Any]] = None, rstar: Optional[float] = None,
                 rp_cap: Optional[Any] = None, e_rad: Optional[float] = None):
    """(flags, info) of the Rp - 2 sigma radius cap for one row; ({}, info) when the rule is off.

    Rejects only when (Rp - 2 sigma_Rp) > cap (strictly; equality, within 1e-9 R_Earth, is kept).
    `missing_rstar` / `rp_undetermined` / `rp_err_incomplete` are distinct from the cap flag and never pass;
    an incomplete error (depth_err or e_rad not valid) never falls back to the point estimate.
    """
    c = rp_cap_cfg(cfg, rp_cap)
    if not c["enabled"]:
        return set(), None
    info = rp_with_error(row, rstar, e_rad)
    info["cap_rearth"] = c["cap_rearth"]
    if info["status"] != "ok":
        return {info["status"]}, info
    if info["lower2s"] > c["cap_rearth"] + _RP_CAP_EPS:
        return {RP_CAP_FLAG}, info
    return set(), info


def _unpack_lc(lc: Any, dumps: Optional[Any] = None):
    if isinstance(lc, dict):
        if "binned" in lc:
            t, f, e = lc["binned"]
        elif "t" in lc:
            t, f, e = lc["t"], lc["f"], lc["e"]
        elif "flat" in lc:
            t, f, e = lc["flat"]
        else:
            raise ValueError("lc dict must contain ('t', 'f', 'e') or 'binned' or 'flat'")
        d = dumps if dumps is not None else lc.get("dumps")
        ge = lc.get("gap_edges")
        return t, f, e, d, ge
    if isinstance(lc, (tuple, list)):
        if len(lc) == 3:
            t, f, e = lc
            return t, f, e, dumps, None
        if len(lc) == 4:
            t, f, e, d_lc = lc
            return t, f, e, (dumps if dumps is not None else d_lc), None
        if len(lc) >= 5:
            t, f, e, d_lc, ge = lc[:5]
            return t, f, e, (dumps if dumps is not None else d_lc), ge
    raise ValueError("lc must be (t, f, e[, dumps[, gap_edges]]) or a dict")


def needs_lc_row(
    row: Dict[str, Any],
    scan: Optional[Dict[str, float]] = None,
    *,
    cfg: Optional[Dict[str, Any]] = None,
    lc: Optional[Any] = None,
    odd_even: Optional[Dict[str, Any]] = None,
    odd_even_mask: Optional[Any] = None,
) -> bool:
    """True when the verdict needs a light curve that was not supplied (deep 'secondary_scan' without
    `scan`, or the v1.2 odd-even transit mask enabled without `lc` or `odd_even`)."""
    sec_nlc = scan is None and "secondary_scan" in _set(row.get("deep_flags"))
    oe_nlc = odd_even_mask_enabled(cfg, odd_even_mask) and lc is None and odd_even is None
    return bool(sec_nlc or oe_nlc)


def check_scan_sigma(row: Dict[str, Any], scan: Dict[str, float]) -> None:
    """vet_amend.revet_one's consistency check for a supplied scan: the rebuilt secondary scan must
    reproduce the row's logged sec_scan_sigma (the rebuild runs at the row's own period and epoch), so
    the rules (a)/(d) decided on it are the ones the v1.1 revet decided.  Raises RuntimeError on a
    mismatch, exactly as revet_one does before calling amended()."""
    got = vet_amend._f((scan or {}).get("sec_scan_sigma"))
    want = vet_amend._f(row.get("sec_scan_sigma"))
    if not np.isclose(got, want, rtol=1e-6, atol=1e-6):
        raise RuntimeError(f"TIC {row.get('tic')}: rebuilt secondary scan {got} != row "
                           f"{row.get('sec_scan_sigma')}")


def final_verdict_row(
    row: Dict[str, Any],
    *,
    common_period: Optional[Any] = None,
    checks_logged: Optional[Any] = None,
    rho: Optional[float] = None,
    scan: Optional[Dict[str, float]] = None,
    needs_lc: str = "reject",
    cfg: Optional[Dict[str, Any]] = None,
    lc: Optional[Any] = None,
    dumps: Optional[Any] = None,
    odd_even: Optional[Dict[str, Any]] = None,
    odd_even_mask: Optional[Any] = None,
    rstar: Optional[float] = None,
    rp_cap: Optional[Any] = None,
    e_rad: Optional[float] = None,
) -> Dict[str, Any]:
    """Row dict + amended fields + the final verdict fields.

    common_period: the vet.flag_common_periods verdict for the row's group (bool); default reads the
    row's 'common_period' flag / 'fail_common_period' column.  rho: TIC density in solar units (as the
    manifest job's 'rho'); None falls back to a/R* from the duration, exactly as vet_amend.amended does.
    scan: vet_amend.scan_depth() output (plus loo_scan_min) when the caller did rebuild the light curve;
    for a deep 'secondary_scan' row it must reproduce the row's logged sec_scan_sigma (revet_one's check,
    RuntimeError otherwise).  needs_lc: policy for rows that require a light curve that was not supplied
    ("reject" = lower bound, "exempt" = upper bound, "raise" = raise ValueError).

    rstar / e_rad: the job's R* and TIC R* error (both R_Sun), used by the Rp - 2 sigma cap (rule on) only when
    the row's own `radius` / `e_rad` is missing or invalid; see `rp_with_error`.

    When the odd-even transit mask is enabled (via `cfg` or `odd_even_mask` or explicit `odd_even`),
    odd/even depths are recomputed on the light curve `lc` using only transits that pass the mask, and
    the stored unmasked `odd_even` flag is never used as a fallback. If `lc` (and `odd_even`) is missing
    when the mask is on, `status_final` is `STATUS_NEEDS_LC` (`needs_lc="reject"` fails with
    `ODD_EVEN_NEEDS_LC`; `needs_lc="raise"` raises `ValueError`).
    """
    if needs_lc not in ("reject", "exempt", "raise"):
        raise ValueError(f"needs_lc must be 'reject', 'exempt' or 'raise', got {needs_lc!r}")
    flags = _set(row.get("flags")) - {COMMON_PERIOD}
    mask_on = odd_even_mask_enabled(cfg, odd_even_mask) or (odd_even is not None)
    oe_res = odd_even
    oe_nlc = False

    if mask_on:
        # Never fall back to the stored unmasked odd_even flag when the mask is enabled.
        flags = flags - {ODD_EVEN, ODD_EVEN_TOO_FEW}
        if oe_res is None and lc is not None:
            t_lc, f_lc, e_lc, d_lc, ge_lc = _unpack_lc(lc, dumps)
            P = vet_amend._f(row.get("P_grid"))
            if not np.isfinite(P):
                P = vet_amend._f(row.get("period"))
            t0 = vet_amend._f(row.get("t0"))
            dur = vet_amend._f(row.get("duration"))
            oe_res = masked_odd_even(
                t_lc, f_lc, e_lc, P, t0, dur, cfg=cfg, dumps=d_lc, gap_edges_t=ge_lc, odd_even_mask=odd_even_mask
            )
        if oe_res is not None:
            if oe_res.get("flag"):
                flags.add(str(oe_res["flag"]))
        else:
            oe_nlc = True
            if needs_lc == "raise":
                raise ValueError(
                    f"TIC {row.get('tic')}: odd_even_mask is enabled but no light curve (lc) or odd_even was supplied"
                )

    row_eff = dict(row, flags=",".join(sorted(flags)))
    if scan is not None and "secondary_scan" in _set(row.get("deep_flags")):
        check_scan_sigma(row, scan)             # revet_one's check: the scan must be the row's own
    sec_nlc = scan is None and "secondary_scan" in _set(row.get("deep_flags"))
    if sec_nlc and needs_lc == "raise":
        raise ValueError(f"TIC {row.get('tic')}: secondary_scan requires a rebuilt scan")
    a = vet_amend.amended(row_eff, scan, rho)
    cap_flags, cap_info = rp_cap_flags(row, cfg, rstar, rp_cap, e_rad)
    if cap_info is not None:
        # rule on: drop the point-estimate cap, apply Rp - 2 sigma (the secondary-exemption 'small' test inside
        # amended() is a separate v1.1 rule and is unchanged)
        a["reject_v2"] = ",".join(sorted((_set(a["reject_v2"]) - {"rp_gt_30re"}) | cap_flags))
        a["passed_lc_checks_v2"] = not _set(a["reject_v2"])
        a.update(rp_cap_rp_v2=cap_info["rp"], rp_cap_sigma_v2=cap_info["sigma"],
                 rp_cap_lower2s_v2=cap_info["lower2s"], rp_cap_rstar_v2=cap_info["rstar"], rp_cap_e_rad_v2=cap_info["e_rad"],
                 rp_cap_rearth_v2=cap_info["cap_rearth"])
    rej = _set(a["reject_v2"])
    status = STATUS_NEEDS_LC if (sec_nlc or oe_nlc) else STATUS_OK
    if sec_nlc and needs_lc == "exempt":
        rej = rej - {"secondary_scan"}          # upper bound: only the undecided scan rule is waived
    if oe_nlc and needs_lc == "reject":
        rej = rej | {ODD_EVEN_NEEDS_LC}
        a["reject_v2"] = ",".join(sorted(rej))
        a["passed_lc_checks_v2"] = False
    cp = fail_common_period_of(row, common_period)
    cl = checks_logged_of(row, checks_logged)
    all_failing: List[str] = sorted(rej)
    if cp:
        all_failing.append(COMMON_PERIOD)
    if not cl:
        all_failing.append(CHECKS_NOT_LOGGED)
    out = dict(row, **a)
    if oe_res is not None:
        out.update(
            depth_odd_v2=oe_res.get("depth_odd"),
            depth_even_v2=oe_res.get("depth_even"),
            odd_even_sigma_v2=oe_res.get("odd_even_sigma"),
            odd_even_rel_v2=oe_res.get("odd_even_rel"),
            n_odd_kept_v2=oe_res.get("n_odd"),
            n_even_kept_v2=oe_res.get("n_even"),
            odd_even_flag_v2=oe_res.get("flag"),
            odd_even_loto_min_sigma_v2=oe_res.get("loto_min_sigma"),
            odd_even_loto_max_sigma_v2=oe_res.get("loto_max_sigma"),
        )
    out.update(passed_final=bool(not all_failing),
               first_failing_rule=all_failing[0] if all_failing else None,
               all_failing=",".join(all_failing),
               status_final=status,
               fail_common_period_final=cp,
               checks_logged_final=cl)
    return out


def final_verdict(
    row: Dict[str, Any],
    *,
    common_period: Optional[Any] = None,
    checks_logged: Optional[Any] = None,
    rho: Optional[float] = None,
    scan: Optional[Dict[str, float]] = None,
    needs_lc: str = "reject",
    cfg: Optional[Dict[str, Any]] = None,
    lc: Optional[Any] = None,
    dumps: Optional[Any] = None,
    odd_even: Optional[Dict[str, Any]] = None,
    odd_even_mask: Optional[Any] = None,
    rstar: Optional[float] = None,
    rp_cap: Optional[Any] = None,
    e_rad: Optional[float] = None,
) -> Verdict:
    """(passed, first_failing_rule, all_failing, status) for one row of any mode.

    Exactly the recipe the v1.1 finish.py applied to real rows: vet_amend.amended() + common_period +
    checks_logged, plus the v1.2 masked odd-even comparison when enabled in `cfg` / `odd_even_mask`.
    See `final_verdict_row` for the arguments.
    """
    d = final_verdict_row(
        row,
        common_period=common_period,
        checks_logged=checks_logged,
        rho=rho,
        scan=scan,
        needs_lc=needs_lc,
        cfg=cfg,
        lc=lc,
        dumps=dumps,
        odd_even=odd_even,
        odd_even_mask=odd_even_mask,
        rstar=rstar,
        rp_cap=rp_cap,
        e_rad=e_rad,
    )
    return Verdict((d["passed_final"], d["first_failing_rule"],
                    [x for x in d["all_failing"].split(",") if x], d["status_final"]))


def bracket(verdicts_lower: List[Verdict], verdicts_upper: List[Verdict]) -> Dict[str, int]:
    """Counts for a set of rows scored under both needs_lc bounds: survivors of each bound and the
    rows whose verdict the two bounds separate (the needs_lc rows that pass everything else)."""
    if len(verdicts_lower) != len(verdicts_upper):
        raise ValueError("lower/upper verdict lists must have the same length")
    lo = sum(1 for v in verdicts_lower if v.passed)
    hi = sum(1 for v in verdicts_upper if v.passed)
    return dict(survived_lower=lo, survived_upper=hi, needs_lc_gap=hi - lo)

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
from .v11 import vet_amend

__all__ = ["Verdict", "final_verdict", "final_verdict_row", "check_scan_sigma", "COMMON_PERIOD",
           "CHECKS_NOT_LOGGED", "STATUS_OK", "STATUS_NEEDS_LC", "VERDICT_KEYS"]

COMMON_PERIOD = "common_period"          # vet.flag_common_periods group flag
CHECKS_NOT_LOGGED = "checks_not_logged"  # row['checks_logged'] false: deep checks were never run
STATUS_OK = "ok"
STATUS_NEEDS_LC = "needs_lc"
VERDICT_KEYS = ("passed_final", "first_failing_rule", "all_failing", "status_final",
                "fail_common_period_final", "checks_logged_final")


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


def needs_lc_row(row: Dict[str, Any], scan: Optional[Dict[str, float]] = None) -> bool:
    """True when amended() would need the rebuilt light curve (deep 'secondary_scan', scan not given)."""
    return scan is None and "secondary_scan" in _set(row.get("deep_flags"))


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


def final_verdict_row(row: Dict[str, Any], *, common_period: Optional[Any] = None,
                      checks_logged: Optional[Any] = None, rho: Optional[float] = None,
                      scan: Optional[Dict[str, float]] = None, needs_lc: str = "reject") -> Dict[str, Any]:
    """Row dict + amended fields + the final verdict fields.

    common_period: the vet.flag_common_periods verdict for the row's group (bool); default reads the
    row's 'common_period' flag / 'fail_common_period' column.  rho: TIC density in solar units (as the
    manifest job's 'rho'); None falls back to a/R* from the duration, exactly as vet_amend.amended does.
    scan: vet_amend.scan_depth() output (plus loo_scan_min) when the caller did rebuild the light curve;
    for a deep 'secondary_scan' row it must reproduce the row's logged sec_scan_sigma (revet_one's check,
    RuntimeError otherwise).  needs_lc: bracket policy for rows that would need a rebuild ("reject" =
    lower bound, "exempt" = upper bound).
    Output keys: every amended() key plus passed_final, first_failing_rule, all_failing, status_final,
    fail_common_period_final, checks_logged_final.
    """
    if needs_lc not in ("reject", "exempt"):
        raise ValueError(f"needs_lc must be 'reject' or 'exempt', got {needs_lc!r}")
    flags = _set(row.get("flags"))
    # common_period is applied as its own term (as in finish.py: it flags the row after the amended
    # verdict), so it must not also enter the amended rejection set.
    row_eff = dict(row, flags=",".join(sorted(flags - {COMMON_PERIOD})))
    if scan is not None and "secondary_scan" in _set(row.get("deep_flags")):
        check_scan_sigma(row, scan)             # revet_one's check: the scan must be the row's own
    a = vet_amend.amended(row_eff, scan, rho)
    rej = _set(a["reject_v2"])
    status = STATUS_NEEDS_LC if needs_lc_row(row, scan) else STATUS_OK
    if status == STATUS_NEEDS_LC and needs_lc == "exempt":
        rej = rej - {"secondary_scan"}          # upper bound: only the undecided scan rule is waived
    cp = fail_common_period_of(row, common_period)
    cl = checks_logged_of(row, checks_logged)
    all_failing: List[str] = sorted(rej)
    if cp:
        all_failing.append(COMMON_PERIOD)
    if not cl:
        all_failing.append(CHECKS_NOT_LOGGED)
    out = dict(row, **a)
    out.update(passed_final=bool(not all_failing),
               first_failing_rule=all_failing[0] if all_failing else None,
               all_failing=",".join(all_failing),
               status_final=status,
               fail_common_period_final=cp,
               checks_logged_final=cl)
    return out


def final_verdict(row: Dict[str, Any], *, common_period: Optional[Any] = None,
                  checks_logged: Optional[Any] = None, rho: Optional[float] = None,
                  scan: Optional[Dict[str, float]] = None, needs_lc: str = "reject") -> Verdict:
    """(passed, first_failing_rule, all_failing, status) for one row of any mode.

    Exactly the recipe the v1.1 finish.py applied to real rows: vet_amend.amended() + common_period +
    checks_logged.  For a deep 'secondary_scan' row the caller supplies the rebuilt-light-curve scan
    (vet_amend.rebuild_flat -> scan_depth -> loo_scan_min, as revet_one does) and rules (a)/(d) are
    decided on it exactly as revet_one decides them; without a scan the status is 'needs_lc'.  See
    final_verdict_row for the arguments.
    """
    d = final_verdict_row(row, common_period=common_period, checks_logged=checks_logged, rho=rho,
                          scan=scan, needs_lc=needs_lc)
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

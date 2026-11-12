#!/usr/bin/env python3
"""Unattended v1.1 beta registration and production launch (coordinator, run by the exo-v11-launch timer).

Each run is idempotent and ends in exactly one state, written to STATUS (outside the repo):
  waiting_r0       the flipped null and scrambled r0 are not finished (no completed max-S
                   products, and no lane-a journal DONE lines)
  gate_failed      beta_select ran; no option passed the TOST gate. Result committed and pushed; nothing launched
  registered_blocked_checks
                   beta registered and committed (configs/v1.1.json beta_option set), but the 3-host gate or the
                   validation is not committed as passing on the current search code; re-checked every run
  launched         manifests + go.json pushed at the current origin head, coordinator share started. Remote hosts
                   start their own share when their agent sees go.json
  error            a step failed; see message (nothing is launched after a failed step)
Terminal states (gate_failed, launched) are never redone; the beta registration runs once. --dry-run evaluates the checks and prints the plan only.

  scripts/v11_launch.py [--dry-run]
"""
import argparse, datetime as dt, json, os, re, subprocess as sp, sys
from pathlib import Path

BRANCH = "exoplanet-search"
REPO = Path(__file__).resolve().parents[2]          # repo root (main checkout; only docs/configs/runs are committed)
SUB = "exoplanet_search"
DATA = Path(os.environ.get("EXO_DATA_ROOT", "/media/nfs01/astro_data/tess"))
OUTROOT = DATA / "products" / "v1.1.0"
STATUS = OUTROOT / "launch_status.json"
WT = Path("/tmp/exo-v11-launch")                    # scratch worktree for registration commits
PROD = Path("/tmp/exo-v11-prod")                    # coordinator production worktree, detached at the pinned commit
# About 50% of cores per host, nice 19; h1 also pauses while its 5-min load is above 18.
WORKERS = {"h1": 8, "h2": 9, "h3": 6, "h4": 6}
# TIC shares proportional to measured throughput (tasks/h per worker x workers), so all hosts finish together.
# Measured 2026-11-10 on the v1.1 gate / validation runs: h1 1251 (with the v1.0 lanes running), h2 3008
# (mixed-task validation), h3 573, h4 771. Expected capacity: h1 10.0k, h2 27.1k, h3 3.4k, h4 4.6k tasks/h.
WEIGHTS = "h1:0.222,h2:0.600,h3:0.076,h4:0.102"
DONE = ["DONE  python -m exosearch.nulltest --mode invert --sectors 70,106",
        "DONE  python -m exosearch.nulltest --mode scramble --realization 0"]
PRIVATE = re.compile(r"/(?:ho)me/|/(?:Us)ers/|192\.168\.|\.local\b")


def sh(*cmd, cwd=None, check=True, **kw):
    r = sp.run(cmd, cwd=cwd, text=True, capture_output=True, **kw)
    if check and r.returncode:
        raise RuntimeError(f"{' '.join(cmd)}: {r.stderr.strip()[-400:]}")
    return r.stdout


def status(state, **kw):
    OUTROOT.mkdir(parents=True, exist_ok=True)
    rec = dict(state=state, at_utc=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), **kw)
    th = OUTROOT / "thermal_h1.json"                          # written every minute by scripts/thermal_guard.sh
    if th.exists():
        try:
            rec["thermal_h1"] = json.loads(th.read_text())
        except ValueError:
            pass
    STATUS.write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps(rec, default=str))
    return rec


def r0_products_done():
    """True when both null max-S tables exist and were written after that run's provenance (start).

    These files are the end-of-run outputs of `exosearch.nulltest` (v1.0 interim layer). They remain
    after the transient exo-lane-a unit is stopped, so detection does not depend on the unit existing.
    """
    root = DATA / "products" / "v1.0.0" / "nulltests"
    for tag in ("invert", "scramble_r0"):
        parquet, prov = root / tag / f"null_{tag}_maxS.parquet", root / tag / "provenance.json"
        if not parquet.is_file() or parquet.stat().st_size <= 0:
            return False
        if prov.is_file() and parquet.stat().st_mtime < prov.stat().st_mtime:
            return False
    return True


def r0_journal_done():
    """DONE lines in the exo-lane-a journal, including after the transient unit is gone.

    If the unit still exists, only lines since its start are used, so an older finished run cannot
    satisfy a new one. Otherwise the retained journal is used in full.
    """
    cmd = ["journalctl", "--user", "-u", "exo-lane-a", "-o", "cat", "--no-pager"]
    since = sh("systemctl", "--user", "show", "exo-lane-a", "-p", "ExecMainStartTimestamp", "--value", check=False).strip()
    if since:
        cmd.extend(["--since", " ".join(since.split()[1:3])])  # drop weekday and zone: "YYYY-MM-DD HH:MM:SS"
    log = sh(*cmd, check=False)
    return all(any(" DONE " in l and d.split("  ", 1)[1] in l for l in log.splitlines()) for d in DONE)


def r0_done():
    """True once the flipped null and scrambled r0 have finished."""
    return r0_products_done() or r0_journal_done()


def origin_json(path):
    r = sp.run(["git", "show", f"origin/{BRANCH}:{SUB}/{path}"], cwd=REPO, text=True, capture_output=True)
    return json.loads(r.stdout) if r.returncode == 0 else None


def checks(registered=None):
    g, v = origin_json("docs/v1.1_checks/host_gate.json"), origin_json("docs/v1.1_checks/validation.json")
    cfg = origin_json("configs/v1.1.json")
    why = []
    if not (g and g.get("all_pass") is True):
        why.append("3-host gate summary not committed as passing")
    if not (v and v.get("all_pass") is True):
        why.append("validation report not committed as passing")
    if cfg is None or (registered and cfg.get("beta_option") != registered):
        why.append("configs/v1.1.json beta_option does not match the registration")
    # Both checks must have run on the search code that would be launched (settings.py only adds host ids;
    # validate.py and beta_select.py are report/selection code, not the search).
    commits = {h.get("commit") for h in (g or {}).get("hosts", [])} | set((v or {}).get("commit", []))
    spec = [f"{SUB}/exosearch", f"{SUB}/configs", f":(exclude){SUB}/exosearch/v11/settings.py",
            f":(exclude){SUB}/exosearch/v11/validate.py", f":(exclude){SUB}/exosearch/v11/beta_select.py",
            f":(exclude){SUB}/configs/v1.1.json"]
    def cfg_wo_beta(ref):
        r = sp.run(["git", "show", f"{ref}:{SUB}/configs/v1.1.json"], cwd=REPO, text=True, capture_output=True)
        return {k: x for k, x in json.loads(r.stdout).items() if k != "beta_option"} if r.returncode == 0 else None
    head_cfg = cfg_wo_beta(f"origin/{BRANCH}")
    for c in sorted(x for x in commits if x):
        if (sp.run(["git", "diff", "--quiet", c, f"origin/{BRANCH}", "--", *spec], cwd=REPO).returncode != 0
                or cfg_wo_beta(c) != head_cfg):
            why.append(f"search code or configs changed since check commit {c[:8]}; rerun the checks")
    return why


def push(wt, paths, msg):
    sh("git", "add", "--", *paths, cwd=wt)
    diff = sh("git", "diff", "--cached", cwd=wt)
    if PRIVATE.search(diff):
        raise RuntimeError("privacy grep hit; not committing")
    sh("git", "-c", "user.name=sakethramanujam", "commit", "-q", "-m", msg, cwd=wt)
    sh("git", "fetch", "-q", "origin", BRANCH, cwd=wt)
    sh("git", "rebase", "-q", f"origin/{BRANCH}", cwd=wt)
    sh("git", "push", "-q", "origin", f"HEAD:{BRANCH}", cwd=wt)            # never forced
    return sh("git", "rev-parse", "HEAD", cwd=wt).strip()


def worktree(path, ref):
    if path.exists():
        sh("git", "worktree", "remove", "--force", str(path), cwd=REPO, check=False)
    sh("git", "worktree", "add", "-q", "--detach", str(path), ref, cwd=REPO)


def register(env):
    """Run beta_select exactly as registered on the flipped null and scrambled r0; commit the result (and the
    beta option when the TOST gate passed). Returns the registration record."""
    worktree(WT, f"origin/{BRANCH}")
    (WT / SUB / ".venv").symlink_to(REPO / SUB / ".venv")
    day = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d")
    reg = Path("docs/preregistration") / f"beta_registration_{day}"
    # Same-population guard (amendment section 5): a failure raises here and registration stops before beta_select.
    sh(".venv/bin/python", "-m", "exosearch.v11.registration_checks", cwd=WT / SUB, env=env)
    sh(".venv/bin/python", "-m", "exosearch.v11.beta_select", "--out", str(reg), cwd=WT / SUB, env=env)
    dec = json.loads((WT / SUB / reg / "beta_selection.json").read_text())["decision"]
    ts = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [f"# Beta-treatment registration, {ts} (UTC)", "",
             "Applied unattended by scripts/v11_launch.py: docs/preregistration/beta_treatment.md with amendment 1,",
             "exactly as written, on the flipped null and scrambled r0 (nulls and injections only).", "",
             f"* TOST gate (90% CI of dS* inside [-0.5, +0.5], every sector): {dec['gate_tost']}",
             f"* Passing: {dec['passing']}", f"* Decision: **{dec['chosen']}**. {dec['reason']}",
             f"* Gate failed: {dec['gate_failed']}", "", "Full numbers: beta_selection.json in this directory."]
    (WT / SUB / reg / "README.md").write_text("\n".join(lines) + "\n")
    paths = [f"{SUB}/{reg}/README.md", f"{SUB}/{reg}/beta_selection.json"]
    if dec["gate_failed"]:
        R = push(WT, paths, f"beta registration {ts}: TOST gate FAILED for every option; v1.1 production not launched")
        return dict(commit=R, at_utc=ts, decision=dec, option=None)
    cfgp = WT / SUB / "configs/v1.1.json"
    cfg = json.loads(cfgp.read_text()); cfg["beta_option"] = dec["chosen"]
    cfgp.write_text(json.dumps(cfg, indent=2) + "\n")
    R = push(WT, paths + [f"{SUB}/configs/v1.1.json"],
             f"beta registration {ts}: option {dec['chosen']} (TOST gate passed: {dec['passing']}); configs/v1.1.json beta_option set")
    return dict(commit=R, at_utc=ts, decision=dec, option=dec["chosen"])


def launch(reg, env):
    """Pin the current origin head (search code verified by checks()), build manifests, push go.json, start h1."""
    worktree(WT, f"origin/{BRANCH}")
    (WT / SUB / ".venv").symlink_to(REPO / SUB / ".venv")
    P = sh("git", "rev-parse", "HEAD", cwd=WT).strip()
    man = Path("runs/v1.1")
    sh(".venv/bin/python", "-m", "exosearch.v11.run", "build", "--sectors", "70,106", "--weights", WEIGHTS,
       "--out", str(man), cwd=WT / SUB, env=env)
    go = dict(pinned_commit=P, beta_option=reg["option"], registration_commit=reg["commit"], workers=WORKERS,
              weights=WEIGHTS, manifests={h: f"{SUB}/{man}/manifest_{h}.jsonl.gz" for h in WORKERS})
    (WT / SUB / man / "go.json").write_text(json.dumps(go, indent=1) + "\n")
    M = push(WT, [f"{SUB}/{man}/{f}" for f in ("go.json", "manifest_summary.json")] +
             [go["manifests"][h] for h in WORKERS], f"v1.1 production manifests at {P[:8]} (S70 + S106); go")
    worktree(PROD, P)
    (PROD / SUB / ".venv-worker").symlink_to(REPO / SUB / ".venv-worker")
    sh("git", "fetch", "-q", "origin", BRANCH, cwd=REPO)
    mf = OUTROOT / "manifest_h1.jsonl.gz"
    mf.write_bytes(sp.run(["git", "show", f"{M}:{go['manifests']['h1']}"], cwd=REPO, capture_output=True, check=True).stdout)
    cmd = (f"{PROD}/{SUB}/scripts/v11_work.sh h1 {mf} {WORKERS['h1']} {OUTROOT}/h1 & w=$!; "
           f"{PROD}/{SUB}/scripts/v11_backoff.sh $w 18 15; wait $w")
    sh("systemd-run", "--user", "--unit", "exo-v11-h1", "--setenv", f"EXO_DATA_ROOT={DATA}", "bash", "-c", cmd)
    return dict(pinned_commit=P, manifest_commit=M, unit="exo-v11-h1")


def main(dry):
    prev = json.loads(STATUS.read_text()) if STATUS.exists() else {}
    if prev.get("state") in ("gate_failed", "launched"):
        print(f"terminal state {prev['state']}; nothing to do"); return
    reg = prev.get("registration")
    sh("git", "fetch", "-q", "origin", BRANCH, cwd=REPO)
    done = r0_done()
    if dry:
        why = checks(reg and reg.get("option"))
        print(json.dumps(dict(dry_run=True, r0_done=done, registered=bool(reg), launch_blocked_by=why,
                              status_file=str(STATUS), weights=WEIGHTS, workers=WORKERS,
                              would=("wait for r0" if not done else
                                     ("register beta (beta_select, commit), then " if not reg else "") +
                                     ("launch" if not why else "write status; no launch"))), indent=1))
        return
    if not done:
        return status("waiting_r0")
    env = dict(os.environ, PYTHONPATH=".", EXO_DATA_ROOT=str(DATA))
    try:
        if not reg:
            reg = register(env)
            if reg["option"] is None:
                return status("gate_failed", registration=reg)
            sh("git", "fetch", "-q", "origin", BRANCH, cwd=REPO)
        why = checks(reg["option"])
        if why:
            return status("registered_blocked_checks", registration=reg, reasons=why)
        return status("launched", registration=reg, **launch(reg, env))
    except Exception as e:                                    # nothing further runs after a failed step
        return status("error", registration=reg, message=str(e)[-800:])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    main(ap.parse_args().dry_run)

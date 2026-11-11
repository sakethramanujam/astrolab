#!/usr/bin/env python3
"""Unattended v1.1 beta registration and production launch (coordinator, run by the exo-v11-launch timer).

Each run is idempotent and ends in exactly one state, written to STATUS (outside the repo):
  waiting_r0       lane-a has not logged DONE for the flipped null and scrambled r0
  blocked_checks   docs/v1.1_checks/host_gate.json or validation.json missing or not all_pass on origin
  gate_failed      beta_select ran; no option passed the TOST gate. Result committed and pushed; nothing launched
  launched         beta registration committed, configs/v1.1.json beta_option set, manifests + go.json pushed,
                   coordinator share started. Remote hosts start their own share when their agent sees go.json
  error            a step failed; see message (nothing is launched after a failed step)
Terminal states (gate_failed, launched) are never redone. --dry-run evaluates the checks and prints the plan only.

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
WEIGHTS = "h1:0.15,h2:0.78,h3:0.07"
WORKERS = {"h1": 6, "h2": 13, "h3": 6}              # <= 75% of cores; h1 also backs off above load 18
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
    STATUS.write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps(rec, default=str))
    return rec


def r0_done():
    """DONE lines for the flipped null and scrambled r0 in the current lane-a run (journal since its start)."""
    since = sh("systemctl", "--user", "show", "exo-lane-a", "-p", "ExecMainStartTimestamp", "--value", check=False).strip()
    if not since:
        return False
    since = " ".join(since.split()[1:3])                     # drop weekday and zone: "YYYY-MM-DD HH:MM:SS"
    log = sh("journalctl", "--user", "-u", "exo-lane-a", "--since", since, "-o", "cat", "--no-pager", check=False)
    return all(any(" DONE " in l and d.split("  ", 1)[1] in l for l in log.splitlines()) for d in DONE)


def origin_json(path):
    r = sp.run(["git", "show", f"origin/{BRANCH}:{SUB}/{path}"], cwd=REPO, text=True, capture_output=True)
    return json.loads(r.stdout) if r.returncode == 0 else None


def checks():
    g, v = origin_json("docs/v1.1_checks/host_gate.json"), origin_json("docs/v1.1_checks/validation.json")
    cfg = origin_json("configs/v1.1.json")
    why = []
    if not (g and g.get("all_pass") is True):
        why.append("3-host gate summary not committed as passing")
    if not (v and v.get("all_pass") is True):
        why.append("validation report not committed as passing")
    if cfg is None or cfg.get("beta_option") is not None:
        why.append("configs/v1.1.json beta_option already set or unreadable")
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


def main(dry):
    prev = json.loads(STATUS.read_text()) if STATUS.exists() else {}
    if prev.get("state") in ("gate_failed", "launched"):
        print(f"terminal state {prev['state']}; nothing to do"); return
    sh("git", "fetch", "-q", "origin", BRANCH, cwd=REPO)
    done, why = r0_done(), checks()
    if dry:
        print(json.dumps(dict(dry_run=True, r0_done=done, blocked=why, status_file=str(STATUS), weights=WEIGHTS, workers=WORKERS,
                              would=("run beta_select, commit registration, launch if TOST passes" if done and not why
                                     else "write status and exit")), indent=1))
        return
    if not done:
        return status("waiting_r0")
    if why:
        return status("blocked_checks", reasons=why)
    try:
        if WT.exists():
            sh("git", "worktree", "remove", "--force", str(WT), cwd=REPO, check=False)
        sh("git", "worktree", "add", "-q", "--detach", str(WT), f"origin/{BRANCH}", cwd=REPO)
        (WT / SUB / ".venv").symlink_to(REPO / SUB / ".venv")
        day = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d")
        reg = Path("docs/preregistration") / f"beta_registration_{day}"
        env = dict(os.environ, PYTHONPATH=".", EXO_DATA_ROOT=str(DATA))
        sh(".venv/bin/python", "-m", "exosearch.v11.beta_select", "--out", str(reg), cwd=WT / SUB, env=env)
        rep = json.loads((WT / SUB / reg / "beta_selection.json").read_text())
        dec = rep["decision"]
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
            return status("gate_failed", commit=R, decision=dec)
        cfgp = WT / SUB / "configs/v1.1.json"
        cfg = json.loads(cfgp.read_text()); cfg["beta_option"] = dec["chosen"]
        cfgp.write_text(json.dumps(cfg, indent=2) + "\n")
        R = push(WT, paths + [f"{SUB}/configs/v1.1.json"],
                 f"beta registration {ts}: option {dec['chosen']} (TOST gate passed: {dec['passing']}); configs/v1.1.json beta_option set")
        man = Path("runs/v1.1")
        sh(".venv/bin/python", "-m", "exosearch.v11.run", "build", "--sectors", "70,106", "--weights", WEIGHTS,
           "--out", str(man), cwd=WT / SUB, env=env)
        go = dict(pinned_commit=R, beta_option=dec["chosen"], registered_utc=ts, workers=WORKERS,
                  manifests={h: f"{SUB}/{man}/manifest_{h}.jsonl.gz" for h in WORKERS})
        (WT / SUB / man / "go.json").write_text(json.dumps(go, indent=1) + "\n")
        M = push(WT, [f"{SUB}/{man}/{f}" for f in ("go.json", "manifest_summary.json")] +
                 [go["manifests"][h] for h in WORKERS], f"v1.1 production manifests at {R[:8]} (S70 + S106); go")
        if PROD.exists():
            sh("git", "worktree", "remove", "--force", str(PROD), cwd=REPO, check=False)
        sh("git", "worktree", "add", "-q", "--detach", str(PROD), R, cwd=REPO)
        (PROD / SUB / ".venv-worker").symlink_to(REPO / SUB / ".venv-worker")
        mf = OUTROOT / "manifest_h1.jsonl.gz"
        sh("git", "fetch", "-q", "origin", BRANCH, cwd=REPO)
        mf.write_bytes(sp.run(["git", "show", f"{M}:{go['manifests']['h1']}"], cwd=REPO, capture_output=True, check=True).stdout)
        cmd = (f"{PROD}/{SUB}/scripts/v11_work.sh h1 {mf} {WORKERS['h1']} {OUTROOT}/h1 & w=$!; "
               f"{PROD}/{SUB}/scripts/v11_backoff.sh $w 18 15; wait $w")
        sh("systemd-run", "--user", "--unit", "exo-v11-h1", "--setenv", f"EXO_DATA_ROOT={DATA}", "bash", "-c", cmd)
        return status("launched", registration_commit=R, manifest_commit=M, decision=dec, unit="exo-v11-h1")
    except Exception as e:                                    # nothing further runs after a failed step
        return status("error", message=str(e)[-800:])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    main(ap.parse_args().dry_run)

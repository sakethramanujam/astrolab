#!/usr/bin/env python3
"""Publish process-only v1.2 search status to the public Pages site.

Reads pipeline outputs (status files, logs, parquet *counts*) and writes
exoplanet_search/site/data/status_v12.json on the exoplanet-search branch
via the GitHub Contents API. It does not touch any local git checkout,
does not recompute science, and until the S69 gate passes it emits only
whitelisted process fields: stars/rows searched, null rows done, checks,
rate/ETA, thermal pauses, and next steps. No candidate counts, no score
statistics, no TICs, no hostnames.

Recommended command on h1, from the 3x-daily progress routine:

  source /tmp/v12-search-h1/env.sh && \\
    gh api -H "Accept: application/vnd.github.raw" \\
      repos/sakethramanujam/astrolab/contents/exoplanet_search/scripts/publish_status_v12.py?ref=exoplanet-search \\
      > /tmp/v12-search-h1/publish_status_v12.py && \\
    "$PY" /tmp/v12-search-h1/publish_status_v12.py --publish

Or, after the wrapper has been installed on the NAS:

  /media/nfs01/astro_data/tess/site_publish/publish_status_v12.sh
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import re
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

PT = ZoneInfo("America/Los_Angeles")
SCHEMA_VERSION = "astrolab.v12_status.v1"
DEFAULT_ROOT = Path("/media/nfs01/astro_data/tess")
DEFAULT_PROD = DEFAULT_ROOT / "products" / "v1.2.0"
DEFAULT_REPO = "sakethramanujam/astrolab"
DEFAULT_BRANCH = "exoplanet-search"
DEFAULT_PATH = "exoplanet_search/site/data/status_v12.json"
REGISTRATION_TAG = "v1.2-registered"
REGISTRATION_AT = "2026-11-18T17:02:23-08:00"
CODE_PIN = "bc05650"
REGISTERED_HOSTS = ["h1", "h3"]
ALLOWED_HOSTS = {"h1", "h3"}
PROCESS_COLUMNS = ["mode", "sector", "status", "host_id"]

# Nested allow-list. Anything not listed is dropped, then the publisher
# refuses to emit if unexpected keys survive or if science/privacy tokens appear.
ALLOWED = {
    "schema_version": None,
    "process_only": None,
    "generated_at": None,
    "generated_at_pt": None,
    "campaign": None,
    "phase": None,
    "previous": {"campaign": None, "outcome": None},
    "registration": {
        "tag": None,
        "code_pin": None,
        "registered_at": None,
        "frozen_before_search_rows": None,
    },
    "search": {
        "state": None,
        "started_at": None,
        "contributing_hosts": None,
        "registered_hosts": None,
        "mode_rows_done": None,
        "mode_rows_total": None,
        "mode_rows_ok": None,
        "mode_rows_error": None,
        "star_jobs_done": None,
        "star_jobs_total": None,
        "by_mode": None,
        "by_sector": None,
        "rate_mode_rows_per_hour": None,
        "eta_at": None,
        "eta_pt": None,
        "scope": None,
        "n_parts": None,
        "n_parts_unreadable": None,
    },
    "thermal": {
        "paused": None,
        "pause_events": None,
        "resume_events": None,
        "paused_hours": None,
    },
    "checks": [{"id": None, "result": None, "detail": None}],
    "next_steps": None,
    "notes": None,
}

# Keys that must never appear, even nested. "candidates" as a step name in
# next_steps/notes is a string value, not a key.
FORBIDDEN_KEYS = {
    "tic", "tic_id", "tics", "known_tic", "candidate", "candidates",
    "n_candidates", "sde", "rho_z", "s_star", "s_exact", "s_oot", "s_scan",
    "s_v10", "score", "scores", "period", "depth", "signals", "detections",
    "detected", "toi", "ctoi", "toi_hosts", "toi_recovered",
}

FORBIDDEN_TEXT = re.compile(
    r"\b(detected|detection|detections|signals?\s+above|above\s+s\*)\b",
    re.I,
)

GENERIC_PRIVACY = [
    r"192\.168\.",
    r"(^|[^0-9.])10\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}([^0-9]|$)",
    r"(^|[^0-9.])100\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}([^0-9]|$)",
    r"/home/",
    r"/Users/",
    r"\barrakis\b",
    r"\bspaceman1729\b",
    r"\bcostlybatman\b",
    r"\bngc1952\b",
    r"\bspiceman\b",
    r"\bspaceman\b",
]

NEXT_STEPS = [
    "Finish the registered search",
    "Calibration on S59–64",
    "Padding (0A.2.1)",
    "Lock S* (0.005 grid)",
    "Freeze the eligibility list",
    "S69 real search",
    "S69 gate (zero known-planet losses)",
    "Only then: candidates",
]

NOTES = [
    "Until the S69 gate passes, this feed is process only: stars and rows searched, null rows done, checks, rate, ETA, and next steps. It does not include candidate counts, score statistics, or target identities.",
    "Only h1 and h3 produce registered rows. A Mac check under translation did not match bit-for-bit.",
]

SCOPE = (
    "S59–64 real, flipped, and scrambled nulls; S69 nulls only. "
    "S69 real waits for the frozen eligibility list."
)


def die(msg: str, code: int = 1) -> None:
    print(msg, file=sys.stderr)
    sys.exit(code)


def load_json(path: Path):
    if not path.is_file():
        return None
    return json.loads(path.read_text())


def parse_status(path: Path) -> tuple[str, str | None]:
    if not path.is_file():
        return "unknown", None
    text = path.read_text().strip()
    if not text:
        return "unknown", None
    # e.g. "running 2026-11-18 17:18:57 PST"
    parts = text.split(None, 1)
    state = parts[0].strip().lower()
    rest = parts[1].strip() if len(parts) > 1 else None
    return state, rest


def parse_dt(value: str | None) -> datetime | None:
    if not value:
        return None
    s = value.strip()
    s = s.replace("Z", "+00:00")
    s = re.sub(r" PST$", "-0800", s)
    s = re.sub(r" PDT$", "-0700", s)
    # "2026-11-18 17:18:57-0800" or iso
    s = s.replace(" ", "T", 1) if re.match(r"\d{4}-\d{2}-\d{2} ", s) and "T" not in s else s
    if re.search(r"[+-]\d{4}$", s):
        s = s[:-5] + s[-5:-2] + ":" + s[-2:]
    try:
        dt = datetime.fromisoformat(s)
    except ValueError:
        return None
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=PT)
    return dt


def fmt_pt(dt: datetime) -> str:
    local = dt.astimezone(PT)
    return local.strftime("%Y-%m-%d %-I:%M %p PT")


def iso_offset(dt: datetime) -> str:
    return dt.astimezone(PT).isoformat(timespec="seconds")


def privacy_patterns(data_root: Path) -> list[re.Pattern]:
    raw = list(GENERIC_PRIVACY)
    priv = data_root / "private_patterns.txt"
    if priv.is_file():
        for line in priv.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                raw.append(line)
    # de-dupe, compile case-insensitive
    seen = set()
    out = []
    for p in raw:
        if p in seen:
            continue
        seen.add(p)
        out.append(re.compile(p, re.I))
    return out


def load_parts(prod: Path) -> dict:
    parts_dir = prod / "rows" / "parts"
    files = sorted(parts_dir.glob("part-*.parquet")) if parts_dir.is_dir() else []
    mode: Counter = Counter()
    sector_mode: Counter = Counter()
    status: Counter = Counter()
    host: Counter = Counter()
    n = 0
    unreadable = 0
    reader = None
    try:
        import pyarrow.parquet as pq  # type: ignore

        def _read(path: Path):
            table = pq.read_table(path, columns=PROCESS_COLUMNS)
            return table.num_rows, table.to_pydict()

        reader = _read
    except Exception:
        try:
            import pandas as pd  # type: ignore

            def _read_pd(path: Path):
                df = pd.read_parquet(path, columns=PROCESS_COLUMNS)
                return int(len(df)), {c: df[c].tolist() for c in PROCESS_COLUMNS}

            reader = _read_pd
        except Exception as exc:
            die(f"need pyarrow or pandas to count parquet rows (process columns only): {exc}")

    for path in files:
        try:
            nrows, cols = reader(path)
        except Exception:
            unreadable += 1
            continue
        n += nrows
        # only the process columns; never touch science fields
        modes = cols.get("mode") or []
        sectors = cols.get("sector") or []
        statuses = cols.get("status") or [None] * nrows
        hosts = cols.get("host_id") or [None] * nrows
        for i in range(nrows):
            m = modes[i] if i < len(modes) else None
            s = sectors[i] if i < len(sectors) else None
            st = statuses[i] if i < len(statuses) else None
            h = hosts[i] if i < len(hosts) else None
            if m is not None:
                mode[str(m)] += 1
            if st is not None:
                status[str(st)] += 1
            if h is not None:
                label = str(h)
                host[label if label in ALLOWED_HOSTS else "other"] += 1
            if s is not None and m is not None:
                sector_mode[(str(int(s) if str(s).isdigit() else s), str(m))] += 1

    by_mode = {k: int(mode.get(k, 0)) for k in ("real", "invert", "scramble")}
    # include any unexpected mode names only if they are simple tokens
    for k, v in mode.items():
        if k not in by_mode and re.fullmatch(r"[a-z][a-z0-9_]*", k or ""):
            by_mode[k] = int(v)

    by_sector: dict[str, dict] = {}
    for (sec, m), v in sector_mode.items():
        slot = by_sector.setdefault(sec, {"rows_done": 0, "by_mode": {}})
        slot["rows_done"] += int(v)
        if m in ("real", "invert", "scramble") or re.fullmatch(r"[a-z][a-z0-9_]*", m or ""):
            slot["by_mode"][m] = int(v)

    contributing = [h for h in REGISTERED_HOSTS if host.get(h)]
    if host.get("other"):
        contributing.append("other")

    return {
        "n": n,
        "by_mode": by_mode,
        "by_sector": by_sector,
        "status": {k: int(v) for k, v in status.items()},
        "hosts": contributing,
        "n_parts": len(files),
        "n_parts_unreadable": unreadable,
        "ok": int(status.get("ok", 0)),
        "error": int(sum(v for k, v in status.items() if k not in {"ok", None, "None"})),
    }


def parse_work_log(path: Path) -> dict:
    out: dict = {}
    if not path.is_file():
        return out
    text = path.read_text(errors="replace")
    # "83130 tasks, 20 done, 83110 to run on 10 workers"
    m = re.search(r"(\d+)\s+tasks,\s+(\d+)\s+done,\s+(\d+)\s+to run", text)
    if m:
        out["star_jobs_total"] = int(m.group(1))
        out["star_jobs_predone"] = int(m.group(2))
        out["star_jobs_to_run"] = int(m.group(3))
    lines = re.findall(
        r"^\s+(\d+)/(\d+)\s+([\d.]+)\s+s/task/worker\s+eta\s+(\d+)\s+min",
        text,
        re.M,
    )
    if lines:
        done, total, secs, eta_min = lines[-1]
        out["star_jobs_this_run"] = int(done)
        out["star_jobs_this_run_total"] = int(total)
        out["s_per_task_worker"] = float(secs)
        out["eta_min_log"] = int(eta_min)
    if "full run finished" in text or re.search(r"^done rc=", text, re.M):
        out["log_finished"] = True
    return out


def parse_thermal(prod: Path, started: datetime | None, now: datetime) -> dict:
    snap = load_json(prod / "thermal.json") or {}
    paused = bool(snap.get("paused", False))
    log_path = prod / "thermal.json.log"
    n_pause = n_resume = 0
    paused_s = 0.0
    cur = None
    if log_path.is_file():
        for line in log_path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                ev = json.loads(line)
            except json.JSONDecodeError:
                continue
            action = ev.get("action")
            t = parse_dt(ev.get("at_utc") or ev.get("at"))
            if t is None:
                continue
            if started is not None and t < started:
                # still count the state so a pause that began before start is handled
                pass
            if action == "paused":
                n_pause += 1
                if cur is None:
                    cur = t
            elif action == "resumed":
                n_resume += 1
                if cur is not None:
                    a, b = cur, t
                    if started is not None and b < started:
                        cur = None
                        continue
                    if started is not None and a < started:
                        a = started
                    paused_s += max(0.0, (b - a).total_seconds())
                    cur = None
        if cur is not None:
            a = cur
            if started is not None and a < started:
                a = started
            paused_s += max(0.0, (now - a).total_seconds())
    return {
        "paused": paused,
        "pause_events": n_pause,
        "resume_events": n_resume,
        "paused_hours": round(paused_s / 3600.0, 3),
    }


def parse_smoke(path: Path) -> str:
    if not path.is_file():
        return "unknown"
    text = path.read_text(errors="replace")
    if re.search(r"smoke work finished rc=0", text):
        return "passed"
    if re.search(r"smoke work finished rc=", text) or "SMOKE FAIL" in text:
        return "failed"
    if "SMOKE PASS" in text:
        return "passed"
    return "unknown"


def sector_totals(launch: dict) -> dict[str, int]:
    out: dict[str, int] = {}
    for sec, info in (launch.get("by_sector") or {}).items():
        if isinstance(info, dict) and "tasks" in info:
            out[str(sec)] = int(info["tasks"])
    return out


def build_status(prod: Path, now: datetime) -> dict:
    launch = load_json(prod / "launch_status.json") or {}
    run_id = load_json(prod / "run_id.json") or {}
    state, _ = parse_status(prod / "status")
    started = parse_dt(launch.get("started_at")) or parse_dt(run_id.get("written_at"))
    counts = load_parts(prod)
    work = parse_work_log(prod / "logs" / "work.log")
    thermal = parse_thermal(prod, started, now)

    mode_total = int(launch.get("n_tasks") or 0)
    star_total = int(work.get("star_jobs_total") or launch.get("n_stars") or 0)
    star_done = None
    if "star_jobs_this_run" in work:
        star_done = int(work.get("star_jobs_predone") or 0) + int(work["star_jobs_this_run"])
        if "star_jobs_this_run_total" in work and "star_jobs_predone" in work:
            # prefer the log's implied total when present
            implied = int(work["star_jobs_predone"]) + int(work["star_jobs_this_run_total"])
            if implied:
                star_total = implied

    if state == "unknown" and work.get("log_finished"):
        state = "done"
    if launch.get("state") and state == "unknown":
        state = str(launch["state"])

    totals = sector_totals(launch)
    by_sector = {}
    # stable sector order
    keys = sorted(set(totals) | set(counts["by_sector"]), key=lambda x: int(x) if str(x).isdigit() else 0)
    for sec in keys:
        slot = counts["by_sector"].get(sec, {"rows_done": 0, "by_mode": {}})
        item = {
            "rows_done": int(slot.get("rows_done") or 0),
            "by_mode": {k: int(v) for k, v in (slot.get("by_mode") or {}).items()},
        }
        if sec in totals:
            item["rows_total"] = int(totals[sec])
        by_sector[sec] = item

    # wall-clock rate, including heat-guard pauses (pauses are expected to continue)
    rate = None
    eta_at = None
    eta_pt = None
    if started is not None:
        wall_h = max(0.0, (now - started).total_seconds() / 3600.0)
        done_for_rate = counts["n"]
        if wall_h > 0 and done_for_rate > 0:
            rate = round(done_for_rate / wall_h, 1)
            remain = max(0, mode_total - counts["n"]) if mode_total else None
            if remain is not None and rate > 0:
                eta = now.timestamp() + (remain / rate) * 3600.0
                eta_dt = datetime.fromtimestamp(eta, tz=timezone.utc)
                eta_at = iso_offset(eta_dt)
                eta_pt = fmt_pt(eta_dt)

    smoke = parse_smoke(prod / "logs" / "smoke.log")
    pin_ok = False
    commit = str(launch.get("commit") or run_id.get("commit") or "")
    if commit.startswith(CODE_PIN):
        pin_ok = True

    checks = [
        {
            "id": "registration",
            "result": "passed",
            "detail": f"tag {REGISTRATION_TAG}, code pin {CODE_PIN}, frozen before any S59–64 or S69 search rows existed",
        },
        {
            "id": "pre_run_pin",
            "result": "passed" if pin_ok else "unknown",
            "detail": "search is pinned to the registered commit" if pin_ok else "could not confirm the pinned commit from launch_status",
        },
        {
            "id": "smoke_s69_20_star",
            "result": smoke,
            "detail": "20-star S69 smoke (nulls only) recorded in the run log",
        },
        {
            "id": "mac_rosetta_bit_identity",
            "result": "failed",
            "detail": "did not match bit-for-bit; does not produce registered rows",
        },
    ]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "process_only": True,
        "generated_at": iso_offset(now),
        "generated_at_pt": fmt_pt(now),
        "campaign": "v1.2",
        "phase": "registered_search",
        "previous": {"campaign": "v1.1", "outcome": "gate_failed"},
        "registration": {
            "tag": REGISTRATION_TAG,
            "code_pin": CODE_PIN,
            "registered_at": REGISTRATION_AT,
            "frozen_before_search_rows": True,
        },
        "search": {
            "state": state,
            "started_at": iso_offset(started) if started else None,
            "contributing_hosts": counts["hosts"] or ["h1"],
            "registered_hosts": list(REGISTERED_HOSTS),
            "mode_rows_done": int(counts["n"]),
            "mode_rows_total": mode_total,
            "mode_rows_ok": int(counts["ok"]),
            "mode_rows_error": int(counts["error"]),
            "star_jobs_done": star_done,
            "star_jobs_total": star_total or None,
            "by_mode": counts["by_mode"],
            "by_sector": by_sector,
            "rate_mode_rows_per_hour": rate,
            "eta_at": eta_at,
            "eta_pt": eta_pt,
            "scope": SCOPE,
            "n_parts": int(counts["n_parts"]),
            "n_parts_unreadable": int(counts["n_parts_unreadable"]),
        },
        "thermal": thermal,
        "checks": checks,
        "next_steps": list(NEXT_STEPS),
        "notes": list(NOTES),
    }
    return whitelist(payload, ALLOWED)


def whitelist(obj, schema):
    if schema is None:
        return obj
    if isinstance(schema, list):
        proto = schema[0] if schema else None
        if not isinstance(obj, list):
            return []
        return [whitelist(item, proto) for item in obj]
    if isinstance(schema, dict):
        if not isinstance(obj, dict):
            return {}
        out = {}
        for k, sub in schema.items():
            if k in obj and obj[k] is not None:
                out[k] = whitelist(obj[k], sub)
        return out
    return obj


def walk_keys(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield k
            yield from walk_keys(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from walk_keys(v)


def enforce(obj: dict, data_root: Path) -> str:
    text = json.dumps(obj, indent=1, ensure_ascii=False) + "\n"
    bad_keys = [k for k in walk_keys(obj) if str(k).lower() in FORBIDDEN_KEYS]
    if bad_keys:
        die(f"process-only guard: forbidden keys {sorted(set(bad_keys))}")
    if FORBIDDEN_TEXT.search(text):
        die("process-only guard: forbidden science wording in the status JSON")
    def unexpected(node, schema, path="$"):
        bad = []
        if schema is None:
            # free-form process maps (by_mode, by_sector): keys must be safe tokens
            if isinstance(node, dict):
                for k, v in node.items():
                    if not re.fullmatch(r"[A-Za-z0-9_]+", str(k) or ""):
                        bad.append(f"{path}.{k}")
                    bad.extend(unexpected(v, None, f"{path}.{k}"))
            elif isinstance(node, list):
                for i, v in enumerate(node):
                    bad.extend(unexpected(v, None, f"{path}[{i}]"))
            return bad
        if isinstance(schema, list):
            proto = schema[0] if schema else None
            if not isinstance(node, list):
                return [path]
            for i, v in enumerate(node):
                bad.extend(unexpected(v, proto, f"{path}[{i}]"))
            return bad
        if isinstance(schema, dict):
            if not isinstance(node, dict):
                return [path]
            for k, v in node.items():
                if k not in schema:
                    bad.append(f"{path}.{k}")
                else:
                    bad.extend(unexpected(v, schema[k], f"{path}.{k}"))
            return bad
        return []

    extra = unexpected(obj, ALLOWED)
    if extra:
        die(f"process-only guard: unexpected keys {extra}")
    for pat in privacy_patterns(data_root):
        if pat.search(text):
            die(f"privacy guard tripped ({pat.pattern})")
    if not obj.get("process_only"):
        die("process-only guard: process_only is not true")
    return text


def gh(*args: str, input_text: str | None = None, check: bool = True) -> subprocess.CompletedProcess:
    cmd = ["gh", *args]
    return subprocess.run(
        cmd,
        input=input_text,
        text=True,
        capture_output=True,
        check=check,
    )


def current_file_sha(repo: str, path: str, branch: str) -> tuple[str | None, str | None]:
    proc = gh(
        "api",
        "-H",
        "Accept: application/vnd.github+json",
        f"repos/{repo}/contents/{path}?ref={branch}",
        check=False,
    )
    if proc.returncode != 0:
        err = (proc.stderr or "") + (proc.stdout or "")
        if "404" in err or "Not Found" in err:
            return None, None
        die(f"gh api read failed: {err.strip() or proc.returncode}")
    data = json.loads(proc.stdout)
    sha = data.get("sha")
    content = None
    if data.get("encoding") == "base64" and data.get("content"):
        content = base64.b64decode(data["content"]).decode()
    return sha, content


def publish(text: str, repo: str, path: str, branch: str, now: datetime) -> str:
    sha, existing = current_file_sha(repo, path, branch)
    if existing == text:
        print("unchanged; not committing")
        return sha or ""
    body = {
        "message": f"data: v1.2 process status {now.astimezone(timezone.utc).strftime('%Y-%m-%dT%H:%MZ')}",
        "content": base64.b64encode(text.encode()).decode(),
        "branch": branch,
    }
    if sha:
        body["sha"] = sha
    proc = gh(
        "api",
        "--method",
        "PUT",
        f"repos/{repo}/contents/{path}",
        "--input",
        "-",
        input_text=json.dumps(body),
        check=False,
    )
    if proc.returncode != 0:
        die(f"gh api put failed: {(proc.stderr or proc.stdout or '').strip() or proc.returncode}")
    data = json.loads(proc.stdout)
    commit = ((data.get("commit") or {}).get("sha")) or ""
    print(f"published {commit} {path} on {branch}")
    return commit


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prod", type=Path, default=Path(os.environ.get("EXO_DATA_ROOT", DEFAULT_ROOT)) / "products" / "v1.2.0",
                    help="v1.2 products directory (read-only)")
    ap.add_argument("--data-root", type=Path, default=None, help="data root for private_patterns.txt")
    ap.add_argument("--out", type=Path, default=None, help="write JSON here (default: stdout)")
    ap.add_argument("--publish", action="store_true", help="commit the JSON via the GitHub Contents API")
    ap.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY", DEFAULT_REPO))
    ap.add_argument("--branch", default=os.environ.get("STATUS_BRANCH", DEFAULT_BRANCH))
    ap.add_argument("--path", default=DEFAULT_PATH)
    args = ap.parse_args(argv)

    prod = args.prod.resolve() if args.prod.exists() else args.prod
    if not prod.is_dir():
        die(f"products directory not found: {prod}")
    data_root = args.data_root or prod.parent.parent  # .../tess/products/v1.2.0 -> .../tess
    now = datetime.now(tz=PT)
    payload = build_status(prod, now)
    text = enforce(payload, data_root)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(text)

    if args.publish:
        publish(text, args.repo, args.path, args.branch, now)
    return 0


if __name__ == "__main__":
    sys.exit(main())

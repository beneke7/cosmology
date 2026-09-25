#!/usr/bin/env python3
"""Repeat the saved PRIMAT grid and check deterministic abundance outputs."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import time


TASK_DIR = Path(__file__).resolve().parent
REPO_DIR = TASK_DIR / "repo"
BASELINE = json.loads((TASK_DIR / "forward_grid.json").read_text())
sys.path.insert(0, str(REPO_DIR))
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

from primat.backend import run_bbn


common = {
    "DeltaNeff": 0.0,
    "tau_n": 878.4,
    "network": "small",
    "thermal_corrections": True,
    "cache_dir": str(TASK_DIR / "cache"),
    "show_progress": False,
}
current_commit = subprocess.check_output(
    ["git", "-C", str(REPO_DIR), "rev-parse", "HEAD"], text=True
).strip()
if current_commit != BASELINE["source"]["commit"]:
    raise SystemExit(f"source commit changed: {current_commit}")

rows = []
for original in BASELINE["runs"]:
    started = time.monotonic()
    values = run_bbn(
        dict(common, Omegabh2=original["Omegabh2"]),
        force_backend="python",
        log_backend=True,
        progress=False,
    )
    repeated = {
        "Omegabh2": original["Omegabh2"],
        "DoH": float(values["DoH"]),
        "10^5_DoH": float(values["DoH"]) * 1.0e5,
        "YPBBN": float(values["YPBBN"]),
        "Neff": float(values["Neff"]),
        "seconds": time.monotonic() - started,
    }
    checks = {
        key: repeated[key] == original[key]
        for key in ("DoH", "10^5_DoH", "YPBBN", "Neff")
    }
    rows.append({"repeat": repeated, "exact_matches": checks})
    print(json.dumps(rows[-1], sort_keys=True), flush=True)

passed = all(all(row["exact_matches"].values()) for row in rows)
record = {
    "status": "passed" if passed else "failed",
    "source_commit": current_commit,
    "source_data_manifest_sha256": __import__("hashlib").sha256(
        (TASK_DIR / "source_data_hashes.json").read_bytes()
    ).hexdigest(),
    "runs": rows,
}
(TASK_DIR / "reproducibility_check.json").write_text(json.dumps(record, indent=2) + "\n")
if not passed:
    raise SystemExit("one or more abundance outputs did not match exactly")
print("exact abundance-output match for all repeated points")

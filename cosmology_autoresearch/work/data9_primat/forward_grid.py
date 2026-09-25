#!/usr/bin/env python3
"""Run a bounded, standard-background PRIMAT v0.3.2 forward grid."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time


TASK_DIR = Path(__file__).resolve().parent
REPO_DIR = TASK_DIR / "repo"
PACKAGE_DIR = REPO_DIR / "primat"
OUTPUT = TASK_DIR / "forward_grid.json"
HASH_OUTPUT = TASK_DIR / "source_data_hashes.json"
OMEGA_GRID = [0.02150, 0.02175, 0.02205, 0.02225, 0.02250]
PDG_DH_1E5 = 2.508
PDG_DH_SIGMA_1E5 = 0.029
COMMIT = "21ff8f39fa18e3937e9fdf386cfa982361bfdfce"


def file_digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def hash_package_tree() -> dict:
    rows = []
    for path in sorted(PACKAGE_DIR.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        rows.append({
            "path": path.relative_to(REPO_DIR).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": file_digest(path),
        })
    aggregate = hashlib.sha256()
    for row in rows:
        aggregate.update(row["path"].encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(row["sha256"].encode("ascii"))
        aggregate.update(b"\n")
    result = {
        "repository": "https://github.com/CyrilPitrou/primat",
        "tag": "v0.3.2",
        "commit": COMMIT,
        "git_tree": subprocess.check_output(
            ["git", "-C", str(REPO_DIR), "rev-parse", "HEAD^{tree}"], text=True
        ).strip(),
        "package_tree_sha256": aggregate.hexdigest(),
        "package_file_count": len(rows),
        "package_bytes": sum(row["bytes"] for row in rows),
        "files": rows,
    }
    HASH_OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    sys.path.insert(0, str(REPO_DIR))
    hash_record = hash_package_tree()

    import numpy as np
    import scipy
    from primat.backend import run_bbn

    cache_dir = TASK_DIR / "cache"
    cache_dir.mkdir(exist_ok=True)
    common = {
        "DeltaNeff": 0.0,
        "tau_n": 878.4,
        "network": "small",
        "thermal_corrections": True,
        "cache_dir": str(cache_dir),
        "show_progress": False,
    }
    result = {
        "status": "running",
        "purpose": "standard-SBBN forward check; not a posterior or paper-exact reproduction",
        "source": hash_record,
        "environment": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "primat_version": "0.3.2 (pinned Git tag and commit; source import, no install)",
        },
        "calculation": {
            "backend": "pure Python, explicitly forced",
            "omega_b_h2_grid": OMEGA_GRID,
            "omega_b_h2_units": "dimensionless",
            "DeltaNeff": 0.0,
            "Neff_SM_default": 3.044,
            "tau_n_s": 878.4,
            "network": "small: 8 light nuclides and 12 thermonuclear reactions plus n-p weak rates",
            "standard_H_T": "PRIMAT StandardBackground with no extra density/component; DeltaNeff=0",
            "weak_rate_corrections": "upstream defaults, including radiative, finite-mass, thermal and spectral-distortion corrections",
            "nuclear_rates": "central tabulated PRIMAT v0.3.2 rate inputs; no rate nuisance sampling",
            "abundance_uncertainties": "not sampled or propagated in this point grid; PDG observational uncertainty used only for a descriptive residual",
            "PDG_DH_observation": {
                "quantity": "10^5 D/H",
                "mean": PDG_DH_1E5,
                "sigma": PDG_DH_SIGMA_1E5,
                "source": "PDG 2025, Big Bang Nucleosynthesis review, Eq. 24.2",
            },
            "broad_SBBN_prior_center": 0.02205,
            "broad_SBBN_prior_sigma": 0.00043,
            "broad_SBBN_prior_source": "PDG 2025, Big Bang Nucleosynthesis review, Eq. 24.6",
            "parameters_overridden": common,
            "parameters_not_overridden": "all other v0.3.2 defaults",
        },
        "runs": [],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")

    for omega in OMEGA_GRID:
        params = dict(common, Omegabh2=omega)
        started = time.monotonic()
        values = run_bbn(
            params,
            force_backend="python",
            log_backend=True,
            progress=False,
        )
        elapsed = time.monotonic() - started
        row = {
            "Omegabh2": omega,
            "DoH": float(values["DoH"]),
            "10^5_DoH": float(values["DoH"]) * 1.0e5,
            "YPBBN": float(values["YPBBN"]),
            "Neff": float(values["Neff"]),
            "seconds": elapsed,
        }
        if abs(omega - 0.02205) < 1e-12:
            row["PDG_observation_residual_sigma_obs_only"] = (
                row["10^5_DoH"] - PDG_DH_1E5
            ) / PDG_DH_SIGMA_1E5
        result["runs"].append(row)
        result["status"] = "running"
        OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(row, sort_keys=True), flush=True)

    result["status"] = "completed"
    result["elapsed_seconds"] = sum(row["seconds"] for row in result["runs"])
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(f"completed {len(result['runs'])} standard points in {result['elapsed_seconds']:.3f}s")


if __name__ == "__main__":
    main()

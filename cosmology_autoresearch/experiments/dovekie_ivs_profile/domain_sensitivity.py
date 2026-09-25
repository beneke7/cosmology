#!/usr/bin/env python3
"""Wider-box/denser-grid stress screen for the Dovekie IVS profile.

This is a second search using the production scorer under a deliberately
expanded *numerical domain*, not a new prior, posterior, or physical claim.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PROFILE_PATH = HERE / "profile.py"
BASE_RESULT_PATH = HERE / "result.json"
OUT_PATH = HERE / "domain_sensitivity.json"
SURFACE_PATH = HERE / "domain_surface.npz"
EXPANDED_BOUNDS = ((0.005, 0.95), (-5.0, 5.0))
GRID_SIZE = 61
SCORE_MATCH_TOL = 1.0e-5


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_profile_module():
    spec = importlib.util.spec_from_file_location("dovekie_ivs_profile_module", PROFILE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load profile scorer")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    t0 = time.perf_counter()
    if not BASE_RESULT_PATH.is_file():
        raise FileNotFoundError(BASE_RESULT_PATH)
    base = json.loads(BASE_RESULT_PATH.read_text(encoding="utf-8"))
    profile = load_profile_module()
    profile._DATA = profile.sn.load_inputs()
    profile.GRID_N_OM = GRID_SIZE
    profile.GRID_N_G = GRID_SIZE

    surface, scale = profile.evaluate_grid(EXPANDED_BOUNDS)
    finite = np.isfinite(surface) & (surface < profile.INVALID / 2)
    if not np.any(finite):
        raise RuntimeError("expanded profile grid contains no valid physical point")
    flat = int(np.argmin(np.where(finite, surface, np.inf)))
    ig, iom = np.unravel_index(flat, surface.shape)
    grid_min = (
        float(np.linspace(*EXPANDED_BOUNDS[0], GRID_SIZE)[iom]),
        float(np.linspace(*EXPANDED_BOUNDS[1], GRID_SIZE)[ig]),
    )
    attempts = profile.optimize_profile(EXPANDED_BOUNDS, grid_min)
    valid = [a for a in attempts if a["endpoint_valid_physical"]]
    if not valid:
        raise RuntimeError("expanded-box search found no valid optimizer endpoint")
    best = min(valid, key=lambda a: a["independent_endpoint_recheck_chi2"])
    point = [float(x) for x in best["endpoint_parameters"]]
    hits = profile.boundary_hits(point, EXPANDED_BOUNDS)
    refined = profile.validate_point(*point)
    best_chi2 = float(refined["profile"]["chi2"])
    base_chi2 = float(base["results"]["ivs_sn_only_minimum"]["profile"]["chi2"])
    delta = float(base_chi2 - best_chi2)
    if abs(best_chi2 - best["independent_endpoint_recheck_chi2"]) > 1e-7:
        raise ArithmeticError("refined expanded-box score disagrees with optimizer endpoint")

    np.savez_compressed(
        SURFACE_PATH,
        omega_m=np.linspace(*EXPANDED_BOUNDS[0], GRID_SIZE),
        g=np.linspace(*EXPANDED_BOUNDS[1], GRID_SIZE),
        chi2_surface=surface,
    )
    runtime = time.perf_counter() - t0
    result = {
        "status": "wider_domain_denser_grid_stress_screen_complete",
        "interpretation": "Expanded finite numerical search domain only; its bounds are not a prior. This sensitivity screen reuses the reviewed production score implementation and does not independently verify its equations or prove global optimality.",
        "base_result": str(BASE_RESULT_PATH.relative_to(ROOT)),
        "base_result_sha256": sha256(BASE_RESULT_PATH),
        "profile_script_sha256": sha256(PROFILE_PATH),
        "domain_sensitivity_script_sha256": sha256(Path(__file__)),
        "input_hashes": {
            "hubble_diagram": profile._DATA["hashes"]["des_dovekie_hd.csv"],
            "packed_precision": profile._DATA["hashes"]["des_dovekie_stat_sys.npz"],
        },
        "base_domain": base["search"]["final_bounds"],
        "expanded_domain": {
            "Omega_m0": list(EXPANDED_BOUNDS[0]),
            "g": list(EXPANDED_BOUNDS[1]),
            "role": "search stress test only, not a paper prior",
        },
        "grid_shape": [GRID_SIZE, GRID_SIZE],
        "scale_benchmark": scale,
        "optimizer_method": "same L-BFGS-B objective and 16 reference/corner/seeded starts; expanded-grid minimum added as an explicit start",
        "optimizer_attempts": attempts,
        "best_point": {
            "parameters": {"Omega_m0": point[0], "g": point[1]},
            "chi2": best_chi2,
            "delta_chi2_base_box_minimum_minus_expanded_minimum": delta,
            "boundary_hits": hits,
            "rechecked_score_abs_difference": abs(best_chi2 - best["independent_endpoint_recheck_chi2"]),
            "refined_physicality": refined["physicality"],
            "distance_quadrature_check": refined["distance_quadrature_check"],
        },
        "base_box_best_point": base["results"]["ivs_sn_only_minimum"]["parameters"],
        "numerical_interpretation": (
            "The expanded-box minimum matches the base-box score within the declared tolerance and has no expanded-box boundary contact."
            if abs(delta) <= SCORE_MATCH_TOL and not hits
            else "The expanded-box screen found a lower score or boundary contact; inspect the optimizer/grid records and extend/review before interpreting the original-box optimum."
        ),
        "score_match_tolerance": SCORE_MATCH_TOL,
        "surface_npz": str(SURFACE_PATH.relative_to(ROOT)),
        "surface_sha256": sha256(SURFACE_PATH),
        "runtime_seconds": float(runtime),
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "cpu_count_visible": int(os.cpu_count() or 1),
            "gpu_used": False,
        },
        "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --seconds 1200 -- .venv/bin/python experiments/dovekie_ivs_profile/domain_sensitivity.py",
    }
    OUT_PATH.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "runtime_seconds": runtime,
                      "best_point": result["best_point"]["parameters"],
                      "chi2": best_chi2, "delta_from_base": delta,
                      "boundary_hits": hits, "parallel_grid": scale["parallel_grid_used"]},
                     allow_nan=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independently recheck non-converged BAO mock fits with profiled alpha + Powell."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import datetime as dt
import json
import os
from pathlib import Path
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np
from scipy.linalg import cho_solve, solve_triangular
from scipy.optimize import minimize

from background_bao import BAOData, FIT_BOUNDS, predict_bao


def parameter_vector(fit: dict, model: str) -> np.ndarray:
    names = ["alpha", "Omega_m"] + (["w0"] if model == "wcdm" else []) + (["w0", "wa"] if model == "cpl" else [])
    return np.asarray([fit[name] for name in names], dtype=np.float64)


def recheck_one(task: tuple) -> dict:
    index, budget, key, fit, model, upper, z, value, observable, covariance = task
    data = BAOData(z, value, observable, covariance)
    chol = np.linalg.cholesky(covariance)
    wy = cho_solve((chol, True), value, check_finite=False)
    base_bounds = [tuple(row) for row in FIT_BOUNDS[model]]
    if model == "cpl":
        base_bounds[-1] = (base_bounds[-1][0], float(upper))
    shape_bounds = base_bounds[1:]
    lows = np.asarray([lo for lo, _ in shape_bounds], dtype=np.float64)
    widths = np.asarray([hi - lo for lo, hi in shape_bounds], dtype=np.float64)
    original_pars = parameter_vector(fit, model)

    def profiled(u: np.ndarray, details: bool = False):
        shape = lows + np.asarray(u, dtype=np.float64) * widths
        if model == "lcdm":
            omega_m, w0, wa = float(shape[0]), -1.0, 0.0
        elif model == "wcdm":
            omega_m, w0, wa = float(shape[0]), float(shape[1]), 0.0
        else:
            omega_m, w0, wa = map(float, shape)
        try:
            unit_prediction = predict_bao(z, observable, 1.0, omega_m, w0, wa)
            wf = cho_solve((chol, True), unit_prediction, check_finite=False)
            denominator = float(unit_prediction @ wf)
            if denominator <= 0.0 or not np.isfinite(denominator):
                return (1e100, None, shape) if details else 1e100
            alpha = float(np.clip((unit_prediction @ wy) / denominator, *base_bounds[0]))
            residual = value - alpha * unit_prediction
            whitened = solve_triangular(chol, residual, lower=True, check_finite=False)
            score = float(whitened @ whitened)
            if details:
                names = ["Omega_m"] + (["w0"] if model == "wcdm" else []) + (["w0", "wa"] if model == "cpl" else [])
                return score, {"alpha": alpha, **dict(zip(names, map(float, shape)))}, shape
            return score
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return (1e100, None, shape) if details else 1e100

    starts = []
    original_unit = (original_pars[1:] - lows) / widths
    starts.append(np.clip(original_unit, 0.0, 1.0))
    starts.append(np.full(len(shape_bounds), 0.5, dtype=np.float64))
    runs = []
    begin = time.perf_counter()
    for start_no, x0 in enumerate(starts):
        start_score, start_parameters, _ = profiled(x0, details=True)
        result = minimize(
            profiled, x0, method="Powell", bounds=[(0.0, 1.0)] * len(shape_bounds),
            options={"maxiter": 1000, "xtol": 1e-9, "ftol": 1e-12},
        )
        score, params, shape = profiled(result.x, details=True)
        runs.append({"start": "saved_optimizer_point" if start_no == 0 else "box_center",
                     "initial_chi2": start_score, "initial_parameters": start_parameters,
                     "success": bool(result.success), "message": str(result.message),
                     "nfev": int(result.nfev), "chi2": score, "parameters": params,
                     "normalized_shape": np.asarray(result.x).tolist()})
    lbfgsb = minimize(
        profiled, starts[0], method="L-BFGS-B", bounds=[(0.0, 1.0)] * len(shape_bounds),
        options={"maxiter": 10000, "maxls": 100, "ftol": 1e-14, "gtol": 1e-10},
    )
    lbfgsb_score, lbfgsb_params, _ = profiled(lbfgsb.x, details=True)
    runs.append({"start": "saved_optimizer_point", "method": "L-BFGS-B with profiled alpha",
                 "initial_chi2": profiled(starts[0]), "initial_parameters": None,
                 "success": bool(lbfgsb.success), "message": str(lbfgsb.message),
                 "nfev": int(lbfgsb.nfev), "chi2": lbfgsb_score,
                 "parameters": lbfgsb_params, "normalized_shape": np.asarray(lbfgsb.x).tolist()})
    elapsed = time.perf_counter() - begin
    direct_prediction = predict_bao(z, observable, *original_pars)
    direct_residual = value - direct_prediction
    direct_white = solve_triangular(chol, direct_residual, lower=True, check_finite=False)
    direct_chi2 = float(direct_white @ direct_white)
    best = min(runs, key=lambda row: min(row["chi2"], row["initial_chi2"]))
    return {
        "mock_index": int(index), "start_budget": budget, "candidate": key, "model": model,
        "original_optimizer_success": bool(fit["success"]), "original_optimizer_message": fit["message"],
        "original_reported_chi2": float(fit["chi2"]),
        "original_direct_chi2": direct_chi2,
        "original_chi2_recalculation_abs_difference": abs(direct_chi2 - float(fit["chi2"])),
        "alternative_optimizer_runs": runs, "best_checked_chi2": min(best["chi2"], best["initial_chi2"], direct_chi2),
        "best_checked_parameters": best["initial_parameters"] if best["initial_chi2"] <= min(best["chi2"], direct_chi2) else best["parameters"],
        "best_checked_minus_original_direct_chi2": float(min(best["chi2"], best["initial_chi2"], direct_chi2) - direct_chi2),
        "runtime_seconds": elapsed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", default="experiments/bao_robustness/result.json")
    parser.add_argument("--mean", default="context/data/desi_dr2_mean.txt")
    parser.add_argument("--cov", default="context/data/desi_dr2_cov.txt")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", default="experiments/bao_robustness/optimizer_recheck.json")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    result_path = Path(args.result)
    result = json.loads(result_path.read_text(encoding="utf-8"))
    data = BAOData.from_files(args.mean, args.cov)
    boot = result["adaptive_search_null_mock_calibration"]
    baseline = result["baseline_profile_fits"]["lcdm"]
    mean = predict_bao(data.z, data.observable, baseline["alpha"], baseline["Omega_m"], -1.0, 0.0)
    rng = np.random.default_rng(int(boot["random_seed"]))
    normals = rng.standard_normal((int(boot["n_mocks"]), len(data.z)))
    mocks = mean[None, :] + normals @ np.linalg.cholesky(data.covariance).T

    tasks = []
    for row in boot["realizations"]:
        mock_values = mocks[int(row["index"])]
        for budget_key, budget_name in (("primary_16_starts", "16_starts"),
                                        ("expanded_64_starts", "64_starts")):
            fits = row[budget_key]["candidate_fits"]
            for key, fit in fits.items():
                if fit["success"]:
                    continue
                model = fit["model"]
                upper = float(key.rsplit("_", 1)[1]) if key.startswith("cpl_w0_upper_") else None
                tasks.append((row["index"], budget_name, key, fit, model, upper,
                              data.z, mock_values, data.observable, data.covariance))

    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        checked = list(pool.map(recheck_one, tasks, chunksize=1))
    runtime = time.perf_counter() - started
    output = {
        "status": "complete",
        "classification": "diagnostic independent derivative-free recheck of L-BFGS-B candidates marked non-converged; does not itself prove global optimality",
        "source_result": args.result,
        "command": f".venv/bin/python scripts/run_bounded.py --seconds 1800 -- .venv/bin/python scripts/recheck_bao_optimizer.py --workers {args.workers}",
        "optimizer": "SciPy Powell over unit-scaled shape parameters with alpha analytically profiled using the full covariance (saved point and box center), plus a reconditioned unit-scaled L-BFGS-B retry",
        "covariance_handling": "Cholesky solves with the released full 13x13 covariance; direct residual chi-square recomputed at the original reported candidate",
        "reconstructed_null_generator": "NumPy default_rng PCG64; original fixed seed and full-covariance Cholesky transform",
        "random_seed": int(boot["random_seed"]), "n_mocks": int(boot["n_mocks"]),
        "failed_candidates_rechecked": len(checked), "workers": args.workers,
        "runtime_seconds": runtime,
        "summary": {
            "original_chi2_recalculation_max_abs_difference": max((r["original_chi2_recalculation_abs_difference"] for r in checked), default=0.0),
            "best_checked_minus_original_min": min((r["best_checked_minus_original_direct_chi2"] for r in checked), default=0.0),
            "best_checked_minus_original_max": max((r["best_checked_minus_original_direct_chi2"] for r in checked), default=0.0),
            "alternative_candidates_improved_by_more_than_1e-7": sum(r["best_checked_minus_original_direct_chi2"] < -1e-7 for r in checked),
            "failed_candidates_with_all_retries_successful": sum(all(s["success"] for s in r["alternative_optimizer_runs"]) for r in checked),
            "alternative_optimizer_candidate_count": 3 * len(checked),
        },
        "rechecks": checked,
    }
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["status"], "rechecked": len(checked),
                      "summary": output["summary"], "runtime_seconds": runtime,
                      "output": args.output}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

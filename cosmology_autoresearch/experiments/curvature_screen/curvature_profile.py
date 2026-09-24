#!/usr/bin/env python3
"""Curved-LCDM profile fit and PCG64 null-search calibration for DESI DR2 BAO."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import beta as beta_distribution


HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import background_bao as bg  # noqa: E402


MEAN_PATH = ROOT / "context/data/desi_dr2_mean.txt"
COV_PATH = ROOT / "context/data/desi_dr2_cov.txt"
FLAT_SEED_PATH = ROOT / "experiments/campaign_seed_bao/result.json"
MOCK_RESULT_PATHS = {
    20260924: ROOT / "experiments/bao_robustness/result.json",
    20260925: ROOT / "experiments/bao_robustness/null_extension_1000_seed20260925.json",
}
MOCK_RECHECK_PATH = ROOT / "experiments/bao_robustness/null_extension_1000_optimizer_recheck.json"
PAPER_PATH = ROOT / "context/papers/desi_dr2_bao.pdf"
LYA_PAPER_PATH = ROOT / "context/papers/desi_lya_ap.pdf"
OUTPUT_DIR = HERE.parent
RESULT_PATH = OUTPUT_DIR / "result.json"
JSONL_PATH = OUTPUT_DIR / "mock_realizations.jsonl"
CHECKPOINT_PATH = OUTPUT_DIR / "result_checkpoint.json"
PLOT_PNG_PATH = OUTPUT_DIR / "curvature_null_calibration.png"
PLOT_SVG_PATH = OUTPUT_DIR / "curvature_null_calibration.svg"

ALPHA_BOUNDS = (1.0e-6, 1.0e4)
OMEGA_M_BOUNDS = (0.05, 0.60)
CURVATURE_DOMAINS = (0.20, 0.10, 0.05)
NULL_SEEDS = ((20260924, 200), (20260925, 1000))
N_STARTS = 16
FLOAT = np.float64
_WORKER_CONTEXT = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def omega_de(omega_m: float, omega_k: float) -> float:
    return float(1.0 - FLOAT(omega_m) - FLOAT(omega_k))


def e2_frw(z: np.ndarray | float, omega_m: float, omega_k: float) -> np.ndarray:
    z_arr = np.asarray(z, dtype=FLOAT)
    one_plus_z = FLOAT(1.0) + z_arr
    # Flat closure, with DESI's convention that positive Omega_k is open.
    result = FLOAT(omega_m) * one_plus_z**3 + FLOAT(omega_k) * one_plus_z**2 + FLOAT(omega_de(omega_m, omega_k))
    if not np.all(np.isfinite(result)) or np.any(result <= 0.0):
        raise ValueError("E(z)^2 must be finite and positive within the profile domain")
    return np.asarray(result, dtype=FLOAT)


def curvature_map(chi: np.ndarray | float, omega_k: float) -> np.ndarray:
    chi_arr = np.asarray(chi, dtype=FLOAT)
    if omega_k == 0.0:
        return chi_arr.copy()
    y = FLOAT(omega_k) * chi_arr**2
    series = chi_arr * (
        FLOAT(1.0)
        + y / FLOAT(6.0)
        + y**2 / FLOAT(120.0)
        + y**3 / FLOAT(5040.0)
        + y**4 / FLOAT(362880.0)
    )
    if omega_k > 0.0:
        root = np.sqrt(FLOAT(omega_k))
        direct = np.sinh(root * chi_arr) / root
    else:
        root = np.sqrt(-FLOAT(omega_k))
        direct = np.sin(root * chi_arr) / root
    return np.asarray(np.where(np.abs(y) < FLOAT(1.0e-4), series, direct), dtype=FLOAT)


class BAOProfile:
    """Full-covariance likelihood with analytical profiling over alpha."""

    def __init__(self, z: np.ndarray, observable: np.ndarray, covariance: np.ndarray, order: int = 96):
        self.z = np.asarray(z, dtype=FLOAT)
        self.observable = np.asarray(observable, dtype=str)
        self.covariance = np.asarray(covariance, dtype=FLOAT)
        self.chol = np.linalg.cholesky(self.covariance)
        self.whitener = solve_triangular(self.chol, np.eye(len(self.z), dtype=FLOAT), lower=True, check_finite=False)
        self.sigma = np.sqrt(np.diag(self.covariance))
        self.unique_z, self.row_z_index = np.unique(self.z, return_inverse=True)
        self.order = int(order)
        self.nodes, self.weights = np.polynomial.legendre.leggauss(self.order)
        self.nodes = np.asarray(self.nodes, dtype=FLOAT)
        self.weights = np.asarray(self.weights, dtype=FLOAT)
        self.integral_z = (self.unique_z[:, None] / FLOAT(2.0)) * (self.nodes[None, :] + FLOAT(1.0))

    @classmethod
    def from_bao_data(cls, data: bg.BAOData, order: int = 96) -> "BAOProfile":
        return cls(data.z, data.observable, data.covariance, order)

    def unit_prediction(self, omega_m: float, omega_k: float, order: int | None = None) -> np.ndarray:
        if order is None or int(order) == self.order:
            nodes, weights, integral_z = self.nodes, self.weights, self.integral_z
        else:
            nodes, weights = np.polynomial.legendre.leggauss(int(order))
            nodes, weights = np.asarray(nodes, dtype=FLOAT), np.asarray(weights, dtype=FLOAT)
            integral_z = (self.unique_z[:, None] / FLOAT(2.0)) * (nodes[None, :] + FLOAT(1.0))
        e2_integral = e2_frw(integral_z, omega_m, omega_k)
        chi = (self.unique_z / FLOAT(2.0)) * np.sum(weights[None, :] / np.sqrt(e2_integral), axis=1, dtype=FLOAT)
        dm_unique = curvature_map(chi, omega_k)
        dm = dm_unique[self.row_z_index]
        dh = FLOAT(1.0) / np.sqrt(e2_frw(self.z, omega_m, omega_k))
        dv = np.cbrt(self.z * dm**2 * dh)
        output = np.empty_like(self.z, dtype=FLOAT)
        for label, values in (("DM_over_rs", dm), ("DH_over_rs", dh), ("DV_over_rs", dv)):
            mask = self.observable == label
            output[mask] = values[mask]
        if not np.all(np.isfinite(output)) or np.any(output <= 0.0):
            raise ValueError("unit-amplitude BAO predictions must be finite and positive")
        return output

    def profile_alpha(self, values: np.ndarray, unit_prediction: np.ndarray) -> tuple[float, float, np.ndarray]:
        y = np.asarray(values, dtype=FLOAT)
        white_y = self.whitener @ y
        white_q = self.whitener @ np.asarray(unit_prediction, dtype=FLOAT)
        alpha_raw = FLOAT(np.dot(white_q, white_y) / np.dot(white_q, white_q))
        alpha = FLOAT(np.clip(alpha_raw, ALPHA_BOUNDS[0], ALPHA_BOUNDS[1]))
        residual = white_y - alpha * white_q
        chi2 = FLOAT(np.dot(residual, residual))
        return float(chi2), float(alpha), alpha * np.asarray(unit_prediction, dtype=FLOAT)

    def profile(self, values: np.ndarray, omega_m: float, omega_k: float, order: int | None = None):
        q = self.unit_prediction(omega_m, omega_k, order)
        chi2, alpha, prediction = self.profile_alpha(values, q)
        return chi2, alpha, prediction, q


def scalar_e2(z: float, omega_m: float, omega_k: float) -> float:
    x = 1.0 + float(z)
    value = float(omega_m) * x**3 + float(omega_k) * x**2 + omega_de(omega_m, omega_k)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("scalar E(z)^2 must be finite and positive")
    return value


def scalar_curvature_map(chi: float, omega_k: float) -> float:
    if omega_k == 0.0:
        return float(chi)
    y = float(omega_k) * float(chi) ** 2
    if abs(y) < 1.0e-4:
        return float(chi) * (1.0 + y / 6.0 + y**2 / 120.0 + y**3 / 5040.0 + y**4 / 362880.0)
    if omega_k > 0:
        root = math.sqrt(omega_k)
        return math.sinh(root * chi) / root
    root = math.sqrt(-omega_k)
    return math.sin(root * chi) / root


def scalar_unit_prediction(z: np.ndarray, observable: np.ndarray, omega_m: float, omega_k: float) -> np.ndarray:
    """Independent scalar adaptive-quadrature prediction for observed-fit checks."""
    z = np.asarray(z, dtype=FLOAT)
    observable = np.asarray(observable, dtype=str)
    unique_z, row_index = np.unique(z, return_inverse=True)
    chi = np.asarray([
        quad(lambda zp: 1.0 / math.sqrt(scalar_e2(zp, omega_m, omega_k)), 0.0, float(zi),
             epsabs=1.0e-13, epsrel=1.0e-13, limit=200)[0]
        for zi in unique_z
    ], dtype=FLOAT)
    dm = np.asarray([scalar_curvature_map(float(value), omega_k) for value in chi], dtype=FLOAT)[row_index]
    dh = np.asarray([1.0 / math.sqrt(scalar_e2(float(zi), omega_m, omega_k)) for zi in z], dtype=FLOAT)
    dv = np.cbrt(z * dm**2 * dh)
    output = np.empty_like(z, dtype=FLOAT)
    for label, values in (("DM_over_rs", dm), ("DH_over_rs", dh), ("DV_over_rs", dv)):
        output[observable == label] = values[observable == label]
    return output


def starts_for_fit(curvature_bound: float | None, count: int, reference_omega_m: float) -> list[np.ndarray]:
    if curvature_bound is None:
        values = np.linspace(OMEGA_M_BOUNDS[0], OMEGA_M_BOUNDS[1], count - 1, dtype=FLOAT)
        return [np.asarray([x], dtype=FLOAT) for x in values] + [np.asarray([reference_omega_m], dtype=FLOAT)]
    if count == 16:
        om_grid = (OMEGA_M_BOUNDS[0], (OMEGA_M_BOUNDS[0] + OMEGA_M_BOUNDS[1]) / 2.0, OMEGA_M_BOUNDS[1])
        ok_grid = (-curvature_bound, -curvature_bound / 2.0, 0.0, curvature_bound / 2.0, curvature_bound)
        grid = [np.asarray([om, ok], dtype=FLOAT) for om in om_grid for ok in ok_grid]
        grid = [vector for vector in grid if not (vector[0] == om_grid[1] and vector[1] == 0.0)]
        return [np.asarray([0.3, 0.0], dtype=FLOAT), np.asarray([reference_omega_m, 0.0], dtype=FLOAT)] + grid
    om_grid = np.linspace(OMEGA_M_BOUNDS[0], OMEGA_M_BOUNDS[1], int(round(math.sqrt(count))), dtype=FLOAT)
    ok_grid = np.linspace(-curvature_bound, curvature_bound, int(round(math.sqrt(count))), dtype=FLOAT)
    grid = [np.asarray([om, ok], dtype=FLOAT) for om in om_grid for ok in ok_grid]
    center_index = min(range(len(grid)), key=lambda j: abs(grid[j][0] - reference_omega_m) + abs(grid[j][1]))
    grid[center_index] = np.asarray([reference_omega_m, 0.0], dtype=FLOAT)
    return grid


def bound_hits(shape: np.ndarray, bounds: list[tuple[float, float]]) -> list[dict[str, object]]:
    names = ["Omega_m"] if len(bounds) == 1 else ["Omega_m", "Omega_k"]
    hits = []
    for name, value, (lower, upper) in zip(names, shape, bounds):
        tolerance = 1.0e-7 * (1.0 + abs(float(value)))
        if abs(float(value) - lower) <= tolerance:
            hits.append({"parameter": name, "side": "lower", "bound": float(lower)})
        if abs(float(value) - upper) <= tolerance:
            hits.append({"parameter": name, "side": "upper", "bound": float(upper)})
    return hits


def profile_fit(
    context: BAOProfile,
    values: np.ndarray,
    curvature_bound: float | None,
    reference_omega_m: float,
    starts: int = N_STARTS,
    extra_start: np.ndarray | None = None,
    do_failure_recheck: bool = True,
) -> dict[str, object]:
    bounds = [OMEGA_M_BOUNDS] if curvature_bound is None else [OMEGA_M_BOUNDS, (-curvature_bound, curvature_bound)]
    init_vectors = starts_for_fit(curvature_bound, starts, reference_omega_m)
    if extra_start is not None:
        init_vectors.append(np.asarray(extra_start, dtype=FLOAT))

    def objective(shape: np.ndarray) -> float:
        try:
            om = float(shape[0])
            ok = 0.0 if curvature_bound is None else float(shape[1])
            return context.profile(values, om, ok)[0]
        except (ValueError, FloatingPointError, OverflowError, ZeroDivisionError):
            return 1.0e100

    attempts = []
    best = None
    for start_index, x0 in enumerate(init_vectors):
        result = minimize(
            objective,
            x0,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 1000, "ftol": 1.0e-12, "gtol": 1.0e-8, "maxls": 50},
        )
        row = {
            "start_index": int(start_index),
            "x0": [float(x) for x in x0],
            "x": [float(x) for x in result.x],
            "chi2": float(result.fun),
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message),
            "nit": int(result.nit),
            "nfev": int(result.nfev),
        }
        attempts.append(row)
        if np.isfinite(result.fun) and (best is None or float(result.fun) < best[0]):
            best = (float(result.fun), np.asarray(result.x, dtype=FLOAT), bool(result.success), "L-BFGS-B", start_index)
    if best is None:
        raise RuntimeError("all local optimizer candidates were nonfinite")

    retry = None
    if do_failure_recheck and not best[2]:
        retry_result = minimize(
            objective,
            best[1],
            method="Powell",
            bounds=bounds,
            options={"maxiter": 1500, "xtol": 1.0e-9, "ftol": 1.0e-12},
        )
        retry = {
            "method": "Powell from best L-BFGS-B candidate",
            "x": [float(x) for x in retry_result.x],
            "chi2": float(retry_result.fun),
            "success": bool(retry_result.success),
            "status": int(retry_result.status),
            "message": str(retry_result.message),
            "nit": int(retry_result.nit),
            "nfev": int(retry_result.nfev),
            "delta_chi2_vs_original_candidate": float(retry_result.fun - best[0]),
        }
        if retry_result.success and np.isfinite(retry_result.fun) and retry_result.fun < best[0]:
            best = (float(retry_result.fun), np.asarray(retry_result.x, dtype=FLOAT), True, "Powell retry", None)

    selected_shape = np.asarray(best[1], dtype=FLOAT)
    omega_m = float(selected_shape[0])
    omega_k = 0.0 if curvature_bound is None else float(selected_shape[1])
    chi2, alpha, prediction, _ = context.profile(values, omega_m, omega_k)
    return {
        "model": "flat_LCDM" if curvature_bound is None else "curved_LCDM",
        "Omega_m": omega_m,
        "Omega_k": omega_k,
        "Omega_Lambda": float(omega_de(omega_m, omega_k)),
        "alpha": float(alpha),
        "alpha_bounds": [float(ALPHA_BOUNDS[0]), float(ALPHA_BOUNDS[1])],
        "alpha_bound_active": bool(alpha in ALPHA_BOUNDS),
        "chi2": float(chi2),
        "prediction": [float(x) for x in prediction],
        "parameter_bound_hits": bound_hits(selected_shape, bounds),
        "selected_method": best[3],
        "selected_start_index": best[4],
        "selected_success": bool(best[2]),
        "starts_requested": int(starts),
        "starts_executed": len(init_vectors),
        "failed_start_count": int(sum(not entry["success"] for entry in attempts)),
        "optimizer_attempts": attempts,
        "powell_failure_recheck": retry,
        "all_search_bounds": [[float(x) for x in bound] for bound in bounds],
    }


def init_worker(z, observable, covariance):
    global _WORKER_CONTEXT
    _WORKER_CONTEXT = BAOProfile(z, observable, covariance, order=96)


def fit_mock(task: tuple[int, int, list[float], float]) -> dict[str, object]:
    seed, index, values_list, null_omega_m = task
    values = np.asarray(values_list, dtype=FLOAT)
    flat = profile_fit(_WORKER_CONTEXT, values, None, reference_omega_m=null_omega_m, starts=N_STARTS)
    curved = profile_fit(
        _WORKER_CONTEXT,
        values,
        0.20,
        reference_omega_m=float(flat["Omega_m"]),
        starts=N_STARTS,
    )
    # Curvature is nested: explicitly compare against the flat point in case
    # every numerical local run misses it by roundoff or reports a failed status.
    if float(curved["chi2"]) > float(flat["chi2"]):
        nested_chi2, nested_alpha, nested_prediction, _ = _WORKER_CONTEXT.profile(values, float(flat["Omega_m"]), 0.0)
        if nested_chi2 < float(curved["chi2"]):
            curved = dict(curved)
            curved.update({
                "Omega_m": float(flat["Omega_m"]),
                "Omega_k": 0.0,
                "Omega_Lambda": float(1.0 - float(flat["Omega_m"])),
                "alpha": float(nested_alpha),
                "chi2": float(nested_chi2),
                "prediction": [float(x) for x in nested_prediction],
                "selected_method": "explicit nested flat point",
                "selected_success": True,
            })
    return {
        "seed": int(seed),
        "index": int(index),
        "flat": flat,
        "curved": curved,
        "delta_chi2_curvature_improvement": float(flat["chi2"] - curved["chi2"]),
    }


def adaptive_profile(context: BAOProfile, values: np.ndarray, omega_m: float, omega_k: float) -> tuple[float, float, np.ndarray, np.ndarray]:
    q = scalar_unit_prediction(context.z, context.observable, omega_m, omega_k)
    chi2, alpha, pred = context.profile_alpha(values, q)
    return chi2, alpha, pred, q


def independent_adaptive_fit(
    context: BAOProfile,
    values: np.ndarray,
    main_fit: dict[str, object],
    curvature_bound: float | None,
) -> dict[str, object]:
    bounds = [OMEGA_M_BOUNDS] if curvature_bound is None else [OMEGA_M_BOUNDS, (-curvature_bound, curvature_bound)]

    def objective(shape: np.ndarray) -> float:
        try:
            omega_m = float(shape[0])
            omega_k = 0.0 if curvature_bound is None else float(shape[1])
            return adaptive_profile(context, values, omega_m, omega_k)[0]
        except (ValueError, FloatingPointError, OverflowError):
            return 1.0e100

    x0 = np.asarray([float(main_fit["Omega_m"])] if curvature_bound is None else [float(main_fit["Omega_m"]), float(main_fit["Omega_k"])], dtype=FLOAT)
    result = minimize(
        objective,
        x0,
        method="Powell",
        bounds=bounds,
        options={"maxiter": 1000, "xtol": 1.0e-10, "ftol": 1.0e-12},
    )
    shape = np.asarray(result.x, dtype=FLOAT)
    omega_m = float(shape[0])
    omega_k = 0.0 if curvature_bound is None else float(shape[1])
    chi2, alpha, prediction, q = adaptive_profile(context, values, omega_m, omega_k)
    q_main = context.unit_prediction(float(main_fit["Omega_m"]), float(main_fit["Omega_k"]))
    pred_main = float(main_fit["alpha"]) * q_main
    return {
        "method": "bounded Powell with independent scalar scipy.integrate.quad distances",
        "Omega_m": omega_m,
        "Omega_k": omega_k,
        "alpha": float(alpha),
        "chi2": float(chi2),
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "nit": int(result.nit),
        "nfev": int(result.nfev),
        "chi2_delta_from_main_96_node_profile": float(chi2 - float(main_fit["chi2"])),
        "max_abs_prediction_delta_from_main_96_node": float(np.max(np.abs(prediction - pred_main))),
        "prediction": [float(x) for x in prediction],
        "unit_prediction": [float(x) for x in q],
    }


def alpha_scalar_check(context: BAOProfile, values: np.ndarray, omega_m: float, omega_k: float, q: np.ndarray) -> dict[str, object]:
    white_y = context.whitener @ np.asarray(values, dtype=FLOAT)
    white_q = context.whitener @ np.asarray(q, dtype=FLOAT)

    def objective(alpha: float) -> float:
        residual = white_y - float(alpha) * white_q
        return float(np.dot(residual, residual))

    analytic = float(np.clip(np.dot(white_q, white_y) / np.dot(white_q, white_q), *ALPHA_BOUNDS))
    scalar = minimize_scalar(
        objective,
        method="bounded",
        bounds=ALPHA_BOUNDS,
        options={"xatol": 1.0e-12, "maxiter": 2000},
    )
    return {
        "Omega_m": float(omega_m),
        "Omega_k": float(omega_k),
        "analytic_profile_alpha": float(analytic),
        "bounded_scalar_alpha": float(scalar.x),
        "abs_alpha_difference": float(abs(analytic - scalar.x)),
        "analytic_chi2": float(objective(analytic)),
        "scalar_minimized_chi2": float(scalar.fun),
        "delta_chi2": float(scalar.fun - objective(analytic)),
        "scalar_success": bool(scalar.success),
    }


def boundary_steps(fit: dict[str, object], model_bound: float | None, context: BAOProfile, values: np.ndarray) -> list[dict[str, object]]:
    omega_m = float(fit["Omega_m"])
    omega_k = float(fit["Omega_k"])
    baseline_chi2 = float(fit["chi2"])
    result = []
    if abs(omega_m - OMEGA_M_BOUNDS[0]) < 2.0e-6:
        moved = min(OMEGA_M_BOUNDS[1], omega_m + 1.0e-5)
        result.append({"parameter": "Omega_m", "side": "lower", "step_into_domain": moved - omega_m,
                       "delta_chi2_inward": float(context.profile(values, moved, omega_k)[0] - baseline_chi2)})
    if abs(omega_m - OMEGA_M_BOUNDS[1]) < 2.0e-6:
        moved = max(OMEGA_M_BOUNDS[0], omega_m - 1.0e-5)
        result.append({"parameter": "Omega_m", "side": "upper", "step_into_domain": omega_m - moved,
                       "delta_chi2_inward": float(context.profile(values, moved, omega_k)[0] - baseline_chi2)})
    if model_bound is not None and abs(omega_k + model_bound) < 2.0e-6:
        moved = min(model_bound, omega_k + 1.0e-5)
        result.append({"parameter": "Omega_k", "side": "lower", "step_into_domain": moved - omega_k,
                       "delta_chi2_inward": float(context.profile(values, omega_m, moved)[0] - baseline_chi2)})
    if model_bound is not None and abs(omega_k - model_bound) < 2.0e-6:
        moved = max(-model_bound, omega_k - 1.0e-5)
        result.append({"parameter": "Omega_k", "side": "upper", "step_into_domain": omega_k - moved,
                       "delta_chi2_inward": float(context.profile(values, omega_m, moved)[0] - baseline_chi2)})
    return result


def observed_fit_suite(data: bg.BAOData, context: BAOProfile, flat_seed: dict[str, object]) -> dict[str, object]:
    y = np.asarray(data.value, dtype=FLOAT)
    baseline_om = float(flat_seed["Omega_m"])
    flat = profile_fit(context, y, None, reference_omega_m=baseline_om, starts=16)
    # As a compatibility anchor, compare its predictions and score with the
    # previously recorded flat-LCDM seed fit.
    seed_prediction = np.asarray(bg.predict_bao(
        data.z, data.observable, float(flat_seed["alpha"]), float(flat_seed["Omega_m"]), -1.0, 0.0,
    ), dtype=FLOAT)
    flat["seed_comparison"] = {
        "seed_chi2": float(flat_seed["chi2"]),
        "reproduced_chi2": float(flat["chi2"]),
        "chi2_difference": float(flat["chi2"] - float(flat_seed["chi2"])),
        "max_abs_prediction_difference_at_seed_parameters": float(np.max(np.abs(seed_prediction - np.asarray(flat_seed["prediction"], dtype=FLOAT)))),
    }
    flat["boundary_inward_steps"] = boundary_steps(flat, None, context, y)
    flat["adaptive_quadrature_check"] = independent_adaptive_fit(context, y, flat, None)
    flat["GL192_check"] = profile_order_check(context, y, flat, None, 192)
    flat["alpha_profile_check"] = alpha_scalar_check(context, y, float(flat["Omega_m"]), 0.0, context.unit_prediction(float(flat["Omega_m"]), 0.0))

    domains = {}
    for curvature_bound in CURVATURE_DOMAINS:
        primary = profile_fit(context, y, curvature_bound, reference_omega_m=baseline_om, starts=16)
        primary["delta_chi2_improvement_over_flat"] = float(flat["chi2"] - primary["chi2"])
        primary["boundary_inward_steps"] = boundary_steps(primary, curvature_bound, context, y)
        primary["adaptive_quadrature_at_minimum"] = adaptive_at_minimum(context, y, primary)
        primary["GL192_check"] = profile_order_check(context, y, primary, curvature_bound, 192)
        primary["alpha_profile_check"] = alpha_scalar_check(
            context, y, float(primary["Omega_m"]), float(primary["Omega_k"]),
            context.unit_prediction(float(primary["Omega_m"]), float(primary["Omega_k"])),
        )
        domains[str(curvature_bound)] = primary

    # Larger 64-start plus dense-grid refinements are observed-data stability
    # checks only; the fixed mock search contract remains 16 starts per fit.
    full_curved = domains["0.2"]
    curved64 = profile_fit(context, y, 0.20, reference_omega_m=baseline_om, starts=64)
    curved64["delta_chi2_improvement_over_flat"] = float(flat["chi2"] - curved64["chi2"])
    curved64["reference_16_start_chi2"] = float(full_curved["chi2"])
    curved64["chi2_delta_64_minus_16"] = float(curved64["chi2"] - float(full_curved["chi2"]))

    # A 61x61 deterministic coarse grid is refined locally; it is a search
    # diagnostic, not a proof of global optimality.
    om_values = np.linspace(*OMEGA_M_BOUNDS, 61, dtype=FLOAT)
    ok_values = np.linspace(-0.20, 0.20, 61, dtype=FLOAT)
    surface = np.empty((len(om_values), len(ok_values)), dtype=FLOAT)
    for i, om in enumerate(om_values):
        for j, ok in enumerate(ok_values):
            surface[i, j] = context.profile(y, float(om), float(ok))[0]
    grid_index = np.unravel_index(int(np.argmin(surface)), surface.shape)
    grid_best = np.asarray([om_values[grid_index[0]], ok_values[grid_index[1]]], dtype=FLOAT)
    grid_refined = profile_fit(context, y, 0.20, reference_omega_m=baseline_om, starts=16, extra_start=grid_best)
    grid_refined["delta_chi2_improvement_over_flat"] = float(flat["chi2"] - grid_refined["chi2"])
    grid_record = {
        "grid_points_per_axis": 61,
        "total_points": int(surface.size),
        "grid_chi2_min": float(surface[grid_index]),
        "grid_best_Omega_m": float(grid_best[0]),
        "grid_best_Omega_k": float(grid_best[1]),
        "grid_best_refinement_chi2": float(grid_refined["chi2"]),
        "grid_refinement_delta_chi2_vs_16_start": float(grid_refined["chi2"] - float(full_curved["chi2"])),
        "refined_fit": grid_refined,
        "globality_claim": False,
    }
    return {
        "flat_LCDM": flat,
        "curved_LCDM_by_abs_Omega_k_bound": domains,
        "curved_LCDM_full_domain_64_start_stability": curved64,
        "curved_LCDM_full_domain_coarse_grid_check": grid_record,
    }


def profile_order_check(context: BAOProfile, values: np.ndarray, fit: dict[str, object], curvature_bound: float | None, order: int) -> dict[str, object]:
    omega_m, omega_k = float(fit["Omega_m"]), float(fit["Omega_k"])
    chi96, alpha96, pred96, _ = context.profile(values, omega_m, omega_k, order=96)
    chi_n, alpha_n, pred_n, _ = context.profile(values, omega_m, omega_k, order=order)
    return {
        "nodes": {"reference": 96, "comparison": int(order)},
        "chi2_96": float(chi96),
        "chi2_comparison": float(chi_n),
        "chi2_difference": float(chi_n - chi96),
        "alpha_96": float(alpha96),
        "alpha_comparison": float(alpha_n),
        "max_abs_prediction_difference": float(np.max(np.abs(pred_n - pred96))),
    }


def adaptive_at_minimum(context: BAOProfile, values: np.ndarray, fit: dict[str, object]) -> dict[str, object]:
    chi2, alpha, prediction, q = adaptive_profile(context, values, float(fit["Omega_m"]), float(fit["Omega_k"]))
    pred96 = np.asarray(fit["prediction"], dtype=FLOAT)
    return {
        "method": "scalar adaptive quadrature evaluated at the 96-node best-fit shape",
        "chi2": float(chi2),
        "alpha": float(alpha),
        "chi2_delta_from_96_node_fit": float(chi2 - float(fit["chi2"])),
        "max_abs_prediction_delta_from_96_node_fit": float(np.max(np.abs(prediction - pred96))),
        "prediction": [float(x) for x in prediction],
        "unit_prediction": [float(x) for x in q],
    }


def cp_interval(successes: int, trials: int, confidence: float = 0.95) -> list[float]:
    tail = (1.0 - confidence) / 2.0
    lower = 0.0 if successes == 0 else float(beta_distribution.ppf(tail, successes, trials - successes + 1))
    upper = 1.0 if successes == trials else float(beta_distribution.ppf(1.0 - tail, successes + 1, trials - successes))
    return [lower, upper]


def calibration_summary(records: list[dict[str, object]], observed_delta: float) -> dict[str, object]:
    statistic = np.asarray([float(row["delta_chi2_curvature_improvement"]) for row in records], dtype=FLOAT)
    k = int(np.count_nonzero(statistic >= float(observed_delta)))
    n = int(len(statistic))
    return {
        "n": n,
        "exceedances_at_least_observed": k,
        "tail_fraction": float(k / n),
        "clopper_pearson_95pct_interval": cp_interval(k, n),
        "observed_delta_chi2": float(observed_delta),
        "delta_chi2_quantiles_linear": {
            "q50": float(np.quantile(statistic, 0.50, method="linear")),
            "q90": float(np.quantile(statistic, 0.90, method="linear")),
            "q95": float(np.quantile(statistic, 0.95, method="linear")),
            "q99": float(np.quantile(statistic, 0.99, method="linear")),
            "maximum": float(np.max(statistic)),
        },
        "minimum_delta_chi2": float(np.min(statistic)),
    }


def read_jsonl(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    rows = []
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def write_checkpoint(started: float, rows: list[dict[str, object]], completed_by_seed: dict[str, int], status: str) -> None:
    checkpoint = {
        "status": status,
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": float(time.perf_counter() - started),
        "completed_records": int(len(rows)),
        "completed_by_seed": completed_by_seed,
        "jsonl_path": str(JSONL_PATH.relative_to(ROOT)),
    }
    CHECKPOINT_PATH.write_text(json.dumps(checkpoint, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def run_null_mocks(data: bg.BAOData, null_alpha: float, null_omega_m: float, workers: int, started: float, checkpoint_every: int = 25) -> list[dict[str, object]]:
    existing = read_jsonl(JSONL_PATH)
    have = {(int(row["seed"]), int(row["index"])) for row in existing}
    records_by_key = {(int(row["seed"]), int(row["index"])): row for row in existing}
    all_keys = {(seed, index) for seed, count in NULL_SEEDS for index in range(count)}
    unexpected = have - all_keys
    if unexpected:
        raise ValueError(f"checkpoint contains unrecognized seed/index keys: {sorted(unexpected)[:5]}")
    pending_tasks = []
    for seed, count in NULL_SEEDS:
        rng = np.random.default_rng(seed)
        normals = rng.standard_normal((count, len(data.z)))
        null_mean = np.asarray(bg.predict_bao(
            data.z, data.observable,
            float(null_alpha), null_omega_m, -1.0, 0.0,
        ), dtype=FLOAT)
        cholesky = np.linalg.cholesky(data.covariance)
        mock_values = null_mean[None, :] + normals @ cholesky.T
        for index, values in enumerate(mock_values):
            if (seed, index) not in have:
                pending_tasks.append((seed, index, values.tolist(), null_omega_m))
    worker_n = min(int(workers), len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else int(workers))
    completed_since_checkpoint = 0
    with ProcessPoolExecutor(
        max_workers=worker_n,
        initializer=init_worker,
        initargs=(data.z, data.observable, data.covariance),
    ) as pool:
        for row in pool.map(fit_mock, pending_tasks, chunksize=1):
            key = (int(row["seed"]), int(row["index"]))
            records_by_key[key] = row
            with JSONL_PATH.open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(row, allow_nan=False, separators=(",", ":")) + "\n")
                stream.flush()
            completed_since_checkpoint += 1
            if completed_since_checkpoint >= checkpoint_every:
                have.update((seed, index) for seed, index in records_by_key)
                done_counts = {str(seed): sum(1 for k in records_by_key if k[0] == seed) for seed, _ in NULL_SEEDS}
                write_checkpoint(started, list(records_by_key.values()), done_counts, "running_checkpoint")
                completed_since_checkpoint = 0
                print(f"curvature mocks completed: {len(records_by_key)}/{len(all_keys)}", flush=True)
    records = [records_by_key[key] for key in sorted(records_by_key)]
    done_counts = {str(seed): sum(1 for row in records if int(row["seed"]) == seed) for seed, _ in NULL_SEEDS}
    write_checkpoint(started, records, done_counts, "mock_searches_complete")
    return records


def optimizer_audit(records: list[dict[str, object]]) -> dict[str, object]:
    starts = {"flat": 0, "curved": 0}
    failures = {"flat": 0, "curved": 0}
    fits = {"flat": 0, "curved": 0}
    rechecks = {"flat": 0, "curved": 0}
    recheck_successes = {"flat": 0, "curved": 0}
    max_recheck_improvement = {"flat": 0.0, "curved": 0.0}
    for row in records:
        for key, label in (("flat", "flat"), ("curved", "curved")):
            fit = row[key]
            starts[label] += len(fit["optimizer_attempts"])
            failures[label] += int(fit["failed_start_count"])
            fits[label] += 1
            retry = fit["powell_failure_recheck"]
            if retry is not None:
                rechecks[label] += 1
                recheck_successes[label] += int(bool(retry["success"]))
                max_recheck_improvement[label] = min(
                    max_recheck_improvement[label], float(retry["delta_chi2_vs_original_candidate"]),
                )
    return {
        "fit_count_by_model": fits,
        "start_count_by_model": starts,
        "unsuccessful_start_count_by_model": failures,
        "fits_with_best_candidate_failure_and_Powell_recheck": rechecks,
        "successful_Powell_rechecks": recheck_successes,
        "best_retry_minus_original_chi2_min_by_model": max_recheck_improvement,
        "all_observed_search_starts_retained": True,
        "all_recheck_attempts_retained_in_jsonl": True,
    }


def create_plot(records: list[dict[str, object]], observed_delta: float) -> dict[str, str]:
    stats = np.asarray([float(row["delta_chi2_curvature_improvement"]) for row in records], dtype=FLOAT)
    fig, ax = plt.subplots(figsize=(8.0, 4.8), constrained_layout=True)
    bins = np.linspace(float(np.min(stats)), float(np.max(stats)), 46)
    ax.hist(stats, bins=bins, color="#4477AA", alpha=0.82, edgecolor="white", linewidth=0.35)
    q95 = float(np.quantile(stats, 0.95, method="linear"))
    ax.axvline(observed_delta, color="#BB3344", linewidth=2.0, label=f"observed Δχ² = {observed_delta:.4g}")
    ax.axvline(q95, color="#228833", linewidth=1.7, linestyle="--", label=f"null q95 = {q95:.4g}")
    ax.set_xlabel("Δχ² = χ²(flat LCDM) − χ²(curved LCDM)")
    ax.set_ylabel("Mock count")
    ax.set_title("Finite-search curvature improvement under fitted flat-LCDM null")
    ax.legend(frameon=False)
    fig.savefig(PLOT_PNG_PATH, dpi=180)
    fig.savefig(PLOT_SVG_PATH)
    plt.close(fig)
    return {"png": str(PLOT_PNG_PATH.relative_to(ROOT)), "svg": str(PLOT_SVG_PATH.relative_to(ROOT))}


def domain_checks(context: BAOProfile) -> dict[str, object]:
    x_star = FLOAT(8.0 / 3.0)
    z_star = x_star - FLOAT(1.0)
    min_e2_formula = FLOAT(73.0 / 108.0)
    min_e2_eval = float(e2_frw(z_star, 0.05, -0.20))
    z_check = np.linspace(0.0, 2.33, 101, dtype=FLOAT)
    flat_shape = (0.3, 0.0)
    e0 = float(e2_frw(np.asarray([0.0]), *flat_shape)[0])
    chi_probe = 1.2
    curvature_series_probe = {
        str(ok): float(curvature_map(np.asarray([chi_probe]), ok)[0])
        for ok in (-1.0e-8, 0.0, 1.0e-8)
    }
    pos = []
    for om in (0.05, 0.3, 0.6):
        for ok in (-0.2, 0.0, 0.2):
            values = e2_frw(z_check, om, ok)
            pos.append(float(np.min(values)))
    flat_data_test = context.unit_prediction(0.297, 0.0) * 29.0
    original_flat = bg.predict_bao(context.z, context.observable, 29.0, 0.297, -1.0, 0.0)
    near_ok = {str(ok): float(np.max(np.abs(context.unit_prediction(0.297, ok) - context.unit_prediction(0.297, 0.0)))) for ok in (-1.0e-12, 0.0, 1.0e-12)}
    return {
        "Omega_de_minimum_over_bounds": float(1.0 - OMEGA_M_BOUNDS[1] - 0.20),
        "E2_analytic_minimum_from_contract": float(min_e2_formula),
        "E2_at_analytic_minimum": min_e2_eval,
        "analytic_minimum_redshift": float(z_star),
        "sampled_minimum_E2_across_parameter_redshift_checks": float(min(pos)),
        "sampled_minimum_E2_positive": bool(min(pos) > 0.0),
        "E2_at_z0": e0,
        "E0_is_unity": bool(e0 == 1.0),
        "curvature_kernel_chi1p2": curvature_series_probe,
        "open_kernel_greater_than_flat": bool(curvature_series_probe["1e-08"] > curvature_series_probe["0.0"]),
        "closed_kernel_less_than_flat": bool(curvature_series_probe["-1e-08"] < curvature_series_probe["0.0"]),
        "max_abs_prediction_error_at_Omega_k_zero_vs_background_bao": float(np.max(np.abs(flat_data_test - original_flat))),
        "max_abs_unit_distance_change_for_Omega_k_pm_1e-12": near_ok,
        "alpha_distance_unit_contract": "predictions are dimensionless D_M/rd, D_H/rd, D_V/rd with alpha=c/(H0*rd); H0 and rd remain separately unidentified",
        "positivity_argument": "E2=1+Omega_m*(x^3-1)+Omega_k*(x^2-1) is minimized at Omega_m=0.05, Omega_k=-0.20; its continuous minimum for z<=2.33 is 73/108 at x=8/3, as derived in work/theory2/curvature_contract.md",
    }


def generate_stream_hashes() -> dict[str, object]:
    entries = {}
    for seed, count in NULL_SEEDS:
        first = np.random.default_rng(seed).standard_normal((count, 13))
        again = np.random.default_rng(seed).standard_normal((count, 13))
        entries[str(seed)] = {
            "n_mocks": count,
            "bit_generator": type(np.random.default_rng(seed).bit_generator).__name__,
            "shape": list(first.shape),
            "repeat_bitwise": bool(np.array_equal(first, again)),
            "normals_sha256_little_endian_float64": hashlib.sha256(np.ascontiguousarray(first, dtype="<f8").tobytes()).hexdigest(),
        }
    return {
        "random_generator": "np.random.default_rng(seed), NumPy PCG64",
        "seed_streams": entries,
        "seeds_distinct": NULL_SEEDS[0][0] != NULL_SEEDS[1][0],
        "stream_digests_distinct": entries[str(NULL_SEEDS[0][0])]["normals_sha256_little_endian_float64"] != entries[str(NULL_SEEDS[1][0])]["normals_sha256_little_endian_float64"],
        "mock_transform": "fitted LCDM mean + standard_normal((N,13)) @ cholesky(full covariance).T",
        "full_covariance_draw": True,
    }


def hash_inputs(data: bg.BAOData) -> dict[str, object]:
    stored_sources = {}
    for name, source_path in (("original_200_mock_result", MOCK_RESULT_PATHS[20260924]), ("extension_1000_mock_result", MOCK_RESULT_PATHS[20260925])):
        doc = json.loads(source_path.read_text(encoding="utf-8"))
        recorded = doc["source"]
        stored_sources[name] = {
            "path": str(source_path.relative_to(ROOT)),
            "sha256": sha256(source_path),
            "data_hashes": recorded["data_files"],
            "code_hashes": recorded["code_files"],
        }
    paths = [MEAN_PATH, COV_PATH, FLAT_SEED_PATH, MOCK_RESULT_PATHS[20260924], MOCK_RESULT_PATHS[20260925], MOCK_RECHECK_PATH, PAPER_PATH, LYA_PAPER_PATH, ROOT / "scripts/background_bao.py", ROOT / "scripts/run_bao_robustness.py", ROOT / "work/theory2/curvature_contract.md"]
    return {
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in paths},
        "prior_null_output_source_hashes": stored_sources,
        "recorded_null_data_and_code_hashes_match_current_files": all(
            digest == sha256(ROOT / name)
            for doc_path in MOCK_RESULT_PATHS.values()
            for doc in [json.loads(doc_path.read_text(encoding="utf-8"))]
            for group in doc["source"].values()
            for name, digest in group.items()
        ),
        "mean_row_order": [{"z": float(z), "observable": str(obs), "value": float(value)} for z, obs, value in zip(data.z, data.observable, data.value)],
        "covariance_shape": list(data.covariance.shape),
        "covariance_cholesky_success": True,
        "covariance_min_eigenvalue": float(np.linalg.eigvalsh(data.covariance)[0]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--checkpoint-every", type=int, default=25)
    args = parser.parse_args()
    if args.workers < 1 or args.checkpoint_every < 1:
        parser.error("--workers and --checkpoint-every must be positive")

    started = time.perf_counter()
    data = bg.BAOData.from_files(MEAN_PATH, COV_PATH)
    context = BAOProfile.from_bao_data(data, order=96)
    seed_doc = json.loads(FLAT_SEED_PATH.read_text(encoding="utf-8"))
    baseline_fits = {row["model"]: row for row in seed_doc["models"]}
    NULL_MEAN = baseline_fits["lcdm"]
    # The null extension was generated around the same fixed best-fit mean.
    extension_baseline = json.loads(MOCK_RESULT_PATHS[20260925].read_text(encoding="utf-8"))["adaptive_search_null_mock_calibration"]["observed_searches"]["16_starts"]["candidate_fits"]["lcdm"]
    if abs(float(extension_baseline["alpha"]) - float(NULL_MEAN["alpha"])) > 1.0e-13 or abs(float(extension_baseline["Omega_m"]) - float(NULL_MEAN["Omega_m"])) > 1.0e-13:
        raise ValueError("the two stored mock ensembles use different fitted flat-LCDM null means")
    observed = observed_fit_suite(data, context, NULL_MEAN)
    records = run_null_mocks(
        data, float(NULL_MEAN["alpha"]), float(NULL_MEAN["Omega_m"]), args.workers, started, args.checkpoint_every,
    )
    records = sorted(records, key=lambda row: (int(row["seed"]), int(row["index"])))
    observed_delta = float(observed["flat_LCDM"]["chi2"] - observed["curved_LCDM_by_abs_Omega_k_bound"]["0.2"]["chi2"])
    seed_records = {
        str(seed): [row for row in records if int(row["seed"]) == seed]
        for seed, _ in NULL_SEEDS
    }
    calibration_by_seed = {
        seed: calibration_summary(rows, observed_delta)
        for seed, rows in seed_records.items()
    }
    pooled = calibration_summary(records, observed_delta)
    optimizer = optimizer_audit(records)
    plot_files = create_plot(records, observed_delta)
    domain = domain_checks(context)
    source_hashes = hash_inputs(data)

    all_flat_stats = np.asarray([float(row["flat"]["chi2"]) for row in records], dtype=FLOAT)
    all_curved_stats = np.asarray([float(row["curved"]["chi2"]) for row in records], dtype=FLOAT)
    all_delta = np.asarray([float(row["delta_chi2_curvature_improvement"]) for row in records], dtype=FLOAT)
    if len(records) != sum(n for _, n in NULL_SEEDS):
        raise RuntimeError(f"expected {sum(n for _, n in NULL_SEEDS)} mock records, got {len(records)}")
    if np.any(all_delta < -1.0e-8):
        raise RuntimeError("curved profile score is worse than its nested flat model by more than numerical tolerance")

    # Compact summary per realization is embedded in result.json. Full per-start
    # statuses and retry records remain in mock_realizations.jsonl.
    compact_rows = []
    for row in records:
        compact_rows.append({
            "seed": int(row["seed"]),
            "index": int(row["index"]),
            "flat_chi2": float(row["flat"]["chi2"]),
            "curved_chi2": float(row["curved"]["chi2"]),
            "delta_chi2_curvature_improvement": float(row["delta_chi2_curvature_improvement"]),
            "flat_parameters": {k: float(row["flat"][k]) for k in ("alpha", "Omega_m", "Omega_k")},
            "curved_parameters": {k: float(row["curved"][k]) for k in ("alpha", "Omega_m", "Omega_k")},
            "flat_failed_start_count": int(row["flat"]["failed_start_count"]),
            "curved_failed_start_count": int(row["curved"]["failed_start_count"]),
            "flat_powell_recheck": row["flat"]["powell_failure_recheck"],
            "curved_powell_recheck": row["curved"]["powell_failure_recheck"],
        })
    result = {
        "status": "complete_profiled_curvature_screen",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "13-row compressed DESI DR2 BAO profile likelihood with unchanged mean vector and full covariance; no early-time calibration, perturbations or posterior sampling",
        "interpretation": "Exploratory finite profile-search/bootstrap screen under the fitted flat-LCDM null; no posterior p-value, evidence, or discovery significance claim.",
        "curvature_convention": {
            "source": "DESI DR2 Results II, Eqs. (3)-(6), pp. 4-5; local PDF context/papers/desi_dr2_bao.pdf",
            "omega_k": "Omega_k=-K*c^2/(a0^2*H0^2)=1-Omega_tot; positive is open and uses sinh, negative is closed and uses sin",
            "friedmann": "E2=Omega_m*(1+z)^3+Omega_k*(1+z)^2+Omega_Lambda",
            "transverse_distance": "DM/rd=alpha*S_k(integral_0^z dz/E), with Sk=sinh(sqrt(Omega_k)*chi)/sqrt(Omega_k) for Omega_k>0; chi for zero; sin(sqrt(|Omega_k|)*chi)/sqrt(|Omega_k|) for Omega_k<0",
            "line_of_sight_distance": "DH/rd=alpha/E(z)",
            "volume_distance": "DV/rd=[z*(DM/rd)^2*(DH/rd)]^(1/3)",
        },
        "parameterization_and_search": {
            "flat_closure": "Omega_Lambda=1-Omega_m-Omega_k",
            "shape_bounds": {"Omega_m": list(OMEGA_M_BOUNDS), "Omega_k_primary": [-0.20, 0.20], "Omega_k_nested": [[-0.10, 0.10], [-0.05, 0.05]]},
            "alpha_bounds": list(ALPHA_BOUNDS),
            "alpha_profile": "clip((q^T C^-1 y)/(q^T C^-1 q), alpha_bounds), using Cholesky-whitened residual norm",
            "observed_search": "16 deterministic L-BFGS-B starts per fit; 61x61 coarse grid independently refined; full-domain 64-start stability fit; Powell scalar adaptive-quadrature reoptimization",
            "mock_search": "same 16 deterministic L-BFGS-B starts per flat and curved fit; per-start status stored; selected nonconverged minima receive bounded Powell recheck",
            "mock_curvature_domain": [-0.20, 0.20],
            "confidence_interval": "exact two-sided Clopper-Pearson interval at 95% from beta quantiles",
            "arithmetic": "float64 throughout; CPU only; full 13x13 covariance retained",
        },
        "inputs": source_hashes,
        "null_definition": {
            "fit": "saved best-fit flat LCDM profile from experiments/campaign_seed_bao/result.json",
            "alpha": float(NULL_MEAN["alpha"]),
            "Omega_m": float(NULL_MEAN["Omega_m"]),
            "random_streams": generate_stream_hashes(),
            "seed_runs": [{"seed": int(seed), "n_mocks": int(count)} for seed, count in NULL_SEEDS],
            "reuse": "regenerated PCG64 normal arrays and full-covariance Cholesky mock vectors according to scripts/run_bao_robustness.py; no original flat/CPL fit ensemble was rerun",
        },
        "observed_fits": observed,
        "observed_statistic": {
            "definition": "Delta chi2 = chi2_flat_LCDM - chi2_curved_LCDM; positive values favor the lower curved profile score",
            "value": observed_delta,
            "curved_full_domain_bound_hits": observed["curved_LCDM_by_abs_Omega_k_bound"]["0.2"]["parameter_bound_hits"],
            "statistic_domain": "Omega_k in [-0.20,0.20]",
        },
        "null_calibration": {
            "fixed_search_contract": {"flat_starts_per_fit": N_STARTS, "curved_starts_per_fit": N_STARTS, "curvature_bounds": [-0.20, 0.20], "parameter_starts": "deterministic grids including each mock's flat optimum embedded at Omega_k=0"},
            "by_seed": calibration_by_seed,
            "pooled": pooled,
            "optimizer_audit": optimizer,
            "mock_realization_statistics": compact_rows,
            "full_optimizer_status_records": {
                "path": str(JSONL_PATH.relative_to(ROOT)),
                "sha256": sha256(JSONL_PATH),
                "records": len(records),
                "per_record_contains": "all optimizer starts, bounds, success/status/message/nit/nfev, fit parameters and any Powell retry",
            },
        },
        "domain_checks": domain,
        "lyalpha_overlap_limit": {
            "statement": "The current 13-row BAO vector includes the Ly-alpha BAO block at z=2.33. DESI DR2 Ly-alpha full-shape analysis reports the BAO measurement from the same data and combines its AP and BAO information; do not append it as an independent likelihood without the relevant joint covariance.",
            "source": "DESI DR2 Ly-alpha full-shape paper, abstract and Table II, local context/papers/desi_lya_ap.pdf; its abstract says the BAO constraint is from the same data, and Table II labels AP and isotropic BAO as constraints from the same dataset.",
        },
        "runtime": {
            "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --seconds 1450 -- .venv/bin/python experiments/curvature_screen/curvature_profile.py --workers 20",
            "seconds": None,
            "workers_requested": int(args.workers),
            "workers_used": int(min(int(args.workers), len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else int(args.workers))),
            "cpu_affinity_count": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "platform": platform.platform(),
            "blas_thread_environment": {name: os.environ.get(name) for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
            "gpu_used": False,
        },
    }
    result["runtime"]["seconds"] = float(time.perf_counter() - started)
    plot_files = create_plot(records, observed_delta)
    result["plot_artifacts"] = {
        "files": plot_files,
        "sha256": {str(PLOT_PNG_PATH.relative_to(ROOT)): sha256(PLOT_PNG_PATH), str(PLOT_SVG_PATH.relative_to(ROOT)): sha256(PLOT_SVG_PATH)},
    }
    result["code_sha256"] = {str(HERE.relative_to(ROOT)): sha256(HERE)}
    RESULT_PATH.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    done_counts = {str(seed): sum(1 for row in records if int(row["seed"]) == seed) for seed, _ in NULL_SEEDS}
    write_checkpoint(started, records, done_counts, "complete")
    print(json.dumps({
        "status": result["status"],
        "runtime_seconds": result["runtime"]["seconds"],
        "observed_flat_chi2": observed["flat_LCDM"]["chi2"],
        "observed_curved_chi2": observed["curved_LCDM_by_abs_Omega_k_bound"]["0.2"]["chi2"],
        "observed_delta_chi2": observed_delta,
        "curvature": observed["curved_LCDM_by_abs_Omega_k_bound"]["0.2"]["Omega_k"],
        "bootstrap_pooled": pooled,
        "optimizer_audit": optimizer,
        "result": str(RESULT_PATH.relative_to(ROOT)),
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())

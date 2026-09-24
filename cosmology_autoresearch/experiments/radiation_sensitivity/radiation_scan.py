#!/usr/bin/env python3
"""CPU float64 radiation-sensitivity profile fits for the DESI DR2 BAO vector."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.optimize import differential_evolution, minimize, minimize_scalar


HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import background_bao as bg  # noqa: E402


FLOAT = np.float64
H0_GRID = (50.0, 55.0, 60.0, 65.0, 67.4, 70.0, 75.0, 80.0, 85.0, 90.0)
ALPHA_BOUNDS = (1.0e-6, 1.0e4)
BASELINE_SEED = 20260923
FIT_SEED = 20260924
MEAN_PATH = ROOT / "context/data/desi_dr2_mean.txt"
COV_PATH = ROOT / "context/data/desi_dr2_cov.txt"
SEED_RESULT_PATH = ROOT / "experiments/campaign_seed_bao/result.json"
OUTPUT_PATH = HERE.parent / "result.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def physical_radiation() -> dict[str, float | dict[str, float]]:
    """Return omega_gamma and omega_r from SI blackbody and critical densities."""
    c = FLOAT(299792458.0)
    k_b = FLOAT(1.380649e-23)
    h_p = FLOAT(6.62607015e-34)
    g = FLOAT(6.67430e-11)
    au = FLOAT(149597870700.0)
    pc = FLOAT(648000.0 / math.pi) * au
    h100 = FLOAT(100000.0) / (FLOAT(1.0e6) * pc)
    temperature = FLOAT(2.7255)
    n_eff = FLOAT(3.046)
    a_rad = FLOAT(8.0 * math.pi**5) * k_b**4 / (FLOAT(15.0) * h_p**3 * c**3)
    photon_energy_density = a_rad * temperature**4
    rho_critical_mass = FLOAT(3.0) * h100**2 / (FLOAT(8.0 * math.pi) * g)
    omega_gamma_h2 = (photon_energy_density / c**2) / rho_critical_mass
    neutrino_factor = FLOAT(1.0) + (FLOAT(7.0) / FLOAT(8.0)) * (FLOAT(4.0) / FLOAT(11.0)) ** (FLOAT(4.0) / FLOAT(3.0)) * n_eff
    omega_r_h2 = omega_gamma_h2 * neutrino_factor
    return {
        "constants_si": {
            "c_m_s": float(c),
            "k_B_J_K": float(k_b),
            "h_P_J_s": float(h_p),
            "G_m3_kg_s2": float(g),
            "AU_m": float(au),
            "pc_m_from_648000_over_pi_AU": float(pc),
            "T_CMB_K": float(temperature),
            "N_eff_massless": float(n_eff),
            "H100_s-1": float(h100),
            "a_rad_J_m-3_K-4": float(a_rad),
            "photon_energy_density_J_m-3": float(photon_energy_density),
            "critical_mass_density_at_H100_kg_m-3": float(rho_critical_mass),
            "neutrino_to_photon_energy_factor": float(neutrino_factor),
        },
        "omega_gamma_h2": float(omega_gamma_h2),
        "omega_r_h2": float(omega_r_h2),
        "omega_r_by_h0": {str(h0): float(omega_r_h2 / (FLOAT(h0) / FLOAT(100.0)) ** 2) for h0 in H0_GRID},
        "equations": {
            "radiation_constant": "a_rad = 8*pi^5*k_B^4/(15*h_P^3*c^3)",
            "critical_mass_density": "rho_crit,100 = 3*H100^2/(8*pi*G)",
            "photons": "omega_gamma_h2 = (a_rad*T_CMB^4/c^2)/rho_crit,100",
            "massless_neutrinos": "omega_r = omega_gamma*(1+(7/8)*(4/11)^(4/3)*N_eff)",
            "H0_scaling": "Omega_r(H0) = omega_r_h2/(H0/100 km s^-1 Mpc^-1)^2",
        },
    }


class ProfileLikelihood:
    def __init__(self, data: bg.BAOData, quadrature_order: int = 96):
        self.data = data
        self.nodes, self.weights = np.polynomial.legendre.leggauss(quadrature_order)
        self.nodes = np.asarray(self.nodes, dtype=FLOAT)
        self.weights = np.asarray(self.weights, dtype=FLOAT)
        self.chol = np.linalg.cholesky(data.covariance)
        self.white_y = solve_triangular(self.chol, data.value, lower=True, check_finite=False)
        self.sigma = np.sqrt(np.diag(data.covariance))

    @staticmethod
    def decode(model: str, shape: np.ndarray) -> tuple[float, float, float]:
        omega_m = FLOAT(shape[0])
        if model == "lcdm":
            return float(omega_m), -1.0, 0.0
        if model == "wcdm":
            return float(omega_m), float(shape[1]), 0.0
        if model == "cpl":
            return float(omega_m), float(shape[1]), float(shape[2])
        raise ValueError(model)

    @staticmethod
    def e2(z: np.ndarray | float, omega_m: float, w0: float, wa: float, omega_r: float) -> np.ndarray:
        z = np.asarray(z, dtype=FLOAT)
        zp = FLOAT(1.0) + z
        omega_de = FLOAT(1.0) - FLOAT(omega_m) - FLOAT(omega_r)
        if omega_de <= 0.0 or not np.isfinite(omega_de):
            raise ValueError("flat closure requires positive dark-energy density")
        log_fde = FLOAT(3.0) * (FLOAT(1.0) + FLOAT(w0) + FLOAT(wa)) * np.log1p(z)
        log_fde = log_fde - FLOAT(3.0) * FLOAT(wa) * z / zp
        with np.errstate(over="raise", invalid="raise"):
            e2 = FLOAT(omega_m) * zp**3 + FLOAT(omega_r) * zp**4 + omega_de * np.exp(log_fde)
        if not np.all(np.isfinite(e2)) or np.any(e2 <= 0.0):
            raise ValueError("nonpositive or nonfinite E(z)^2")
        return np.asarray(e2, dtype=FLOAT)

    def unit_prediction(self, model: str, shape: np.ndarray, omega_r: float, order: int | None = None) -> np.ndarray:
        omega_m, w0, wa = self.decode(model, shape)
        if order is None or order == len(self.nodes):
            nodes, weights = self.nodes, self.weights
        else:
            nodes, weights = np.polynomial.legendre.leggauss(order)
            nodes, weights = np.asarray(nodes, dtype=FLOAT), np.asarray(weights, dtype=FLOAT)
        z = self.data.z
        integration_z = (z[:, None] / FLOAT(2.0)) * (nodes[None, :] + FLOAT(1.0))
        e2_i = self.e2(integration_z, omega_m, w0, wa, omega_r)
        dm_unit = (z / FLOAT(2.0)) * np.sum(weights[None, :] / np.sqrt(e2_i), axis=1, dtype=FLOAT)
        dh_unit = FLOAT(1.0) / np.sqrt(self.e2(z, omega_m, w0, wa, omega_r))
        dv_unit = np.cbrt(z * dm_unit**2 * dh_unit)
        shape_prediction = np.empty_like(z, dtype=FLOAT)
        for label, values in (("DM_over_rs", dm_unit), ("DH_over_rs", dh_unit), ("DV_over_rs", dv_unit)):
            mask = self.data.observable == label
            shape_prediction[mask] = values[mask]
        if not np.all(np.isfinite(shape_prediction)) or np.any(shape_prediction <= 0.0):
            raise ValueError("unit-amplitude BAO prediction must be finite and positive")
        return shape_prediction

    def profile(self, model: str, shape: np.ndarray, omega_r: float, order: int | None = None) -> tuple[float, float, np.ndarray, np.ndarray]:
        q = self.unit_prediction(model, shape, omega_r, order)
        white_q = solve_triangular(self.chol, q, lower=True, check_finite=False)
        alpha_unbounded = FLOAT(np.dot(white_q, self.white_y) / np.dot(white_q, white_q))
        alpha = FLOAT(np.clip(alpha_unbounded, ALPHA_BOUNDS[0], ALPHA_BOUNDS[1]))
        white_residual = self.white_y - alpha * white_q
        chi2 = FLOAT(np.dot(white_residual, white_residual))
        return float(chi2), float(alpha), alpha * q, q

    def objective(self, model: str, shape: np.ndarray, omega_r: float) -> float:
        try:
            return self.profile(model, shape, omega_r)[0]
        except (ValueError, FloatingPointError, OverflowError, ZeroDivisionError):
            return 1.0e100


def shape_bounds(model: str) -> list[tuple[float, float]]:
    return [(float(lo), float(hi)) for lo, hi in bg.FIT_BOUNDS[model][1:]]


def start_vectors(model: str, baseline_shape: np.ndarray, random_seed: int, count: int = 16) -> list[np.ndarray]:
    bounds = shape_bounds(model)
    midpoint = np.asarray([(lo + hi) / 2.0 for lo, hi in bounds], dtype=FLOAT)
    nominal = np.asarray(([0.3, -1.0, 0.0])[: len(bounds)], dtype=FLOAT)
    nominal = np.asarray([np.clip(v, lo, hi) for v, (lo, hi) in zip(nominal, bounds)], dtype=FLOAT)
    seeded = np.asarray([np.clip(v, lo, hi) for v, (lo, hi) in zip(baseline_shape, bounds)], dtype=FLOAT)
    starts = [midpoint, nominal, seeded]
    rng = np.random.default_rng(random_seed)
    while len(starts) < count:
        starts.append(np.asarray([rng.uniform(lo, hi) for lo, hi in bounds], dtype=FLOAT))
    return starts[:count]


def fit_profile(
    likelihood: ProfileLikelihood,
    model: str,
    omega_r: float,
    baseline_shape: np.ndarray,
    random_seed: int,
) -> dict[str, object]:
    bounds = shape_bounds(model)
    objective = lambda x: likelihood.objective(model, np.asarray(x, dtype=FLOAT), omega_r)
    initial = start_vectors(model, baseline_shape, random_seed)
    main_runs = []
    best = None
    for x0 in initial:
        result = minimize(
            objective,
            x0,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000, "ftol": 1.0e-12, "gtol": 1.0e-8},
        )
        main_runs.append({
            "x0": [float(x) for x in x0],
            "x": [float(x) for x in result.x],
            "chi2": float(result.fun),
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message),
            "nit": int(result.nit),
            "nfev": int(result.nfev),
        })
        if np.isfinite(result.fun) and (best is None or result.fun < best[0]):
            best = (float(result.fun), np.asarray(result.x, dtype=FLOAT), "L-BFGS-B-start")

    # Independent population-based derivative-free bounded check in shape space.
    de = differential_evolution(
        objective,
        bounds,
        seed=random_seed + 997,
        strategy="best1bin",
        maxiter=500,
        popsize=10,
        tol=1.0e-9,
        atol=1.0e-10,
        mutation=(0.5, 1.0),
        recombination=0.7,
        polish=False,
        updating="immediate",
        workers=1,
    )
    de_record = {
        "x": [float(x) for x in de.x],
        "chi2": float(de.fun),
        "success": bool(de.success),
        "status": getattr(de, "status", None),
        "message": str(de.message),
        "nit": int(de.nit),
        "nfev": int(de.nfev),
        "method": "scipy.optimize.differential_evolution; polish=False",
    }
    if np.isfinite(de.fun) and (best is None or de.fun < best[0]):
        polished = minimize(
            objective,
            de.x,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000, "ftol": 1.0e-12, "gtol": 1.0e-8},
        )
        polish_record = {
            "x0": [float(x) for x in de.x],
            "x": [float(x) for x in polished.x],
            "chi2": float(polished.fun),
            "success": bool(polished.success),
            "status": int(polished.status),
            "message": str(polished.message),
            "nit": int(polished.nit),
            "nfev": int(polished.nfev),
        }
        if np.isfinite(polished.fun) and polished.fun < best[0]:
            best = (float(polished.fun), np.asarray(polished.x, dtype=FLOAT), "DE-seeded-L-BFGS-B")
    else:
        polish_record = None

    # If a box face contains the best point, run a separate bounded Powell
    # optimization on that face and check the feasible inward direction.
    face_checks = []
    best_shape = np.asarray(best[1], dtype=FLOAT)
    for index, ((lower, upper), value) in enumerate(zip(bounds, best_shape)):
        tolerance = 1.0e-7 * (1.0 + abs(float(value)))
        side = "lower" if abs(float(value) - lower) <= tolerance else (
            "upper" if abs(float(value) - upper) <= tolerance else None
        )
        if side is None:
            continue
        face_value = lower if side == "lower" else upper
        free = [j for j in range(len(bounds)) if j != index]
        face_seed = best_shape.copy()
        face_seed[index] = face_value
        if free:
            face_bounds = [bounds[j] for j in free]

            def embed(free_x: np.ndarray) -> np.ndarray:
                vector = face_seed.copy()
                vector[free] = free_x
                return vector

            face_result = minimize(
                lambda free_x: objective(embed(np.asarray(free_x, dtype=FLOAT))),
                face_seed[free],
                method="Powell",
                bounds=face_bounds,
                options={"maxiter": 2000, "xtol": 1.0e-10, "ftol": 1.0e-12},
            )
            face_x = embed(np.asarray(face_result.x, dtype=FLOAT))
            face_record = {
                "parameter_index": index,
                "side": side,
                "bound": float(face_value),
                "method": "scipy.optimize.minimize(method='Powell') with this shape fixed",
                "x": [float(x) for x in face_x],
                "chi2": float(face_result.fun),
                "success": bool(face_result.success),
                "status": int(face_result.status),
                "message": str(face_result.message),
                "nit": int(face_result.nit),
                "nfev": int(face_result.nfev),
            }
            face_checks.append(face_record)
            if np.isfinite(face_result.fun) and face_result.fun < best[0]:
                face_refined = minimize(
                    objective,
                    face_x,
                    method="L-BFGS-B",
                    bounds=bounds,
                    options={"maxiter": 2000, "ftol": 1.0e-12, "gtol": 1.0e-8},
                )
                face_record["refinement"] = {
                    "x": [float(x) for x in face_refined.x],
                    "chi2": float(face_refined.fun),
                    "success": bool(face_refined.success),
                    "status": int(face_refined.status),
                    "message": str(face_refined.message),
                    "nit": int(face_refined.nit),
                    "nfev": int(face_refined.nfev),
                }
                if np.isfinite(face_refined.fun) and face_refined.fun < best[0]:
                    best = (float(face_refined.fun), np.asarray(face_refined.x, dtype=FLOAT), "Powell-face-refinement")
        inward_step = (upper - lower) * 1.0e-5
        inner_x = best_shape.copy()
        inner_x[index] += inward_step if side == "lower" else -inward_step
        face_checks[-1]["inward_step"] = float(inward_step)
        face_checks[-1]["inward_delta_chi2"] = float(objective(inner_x) - objective(best_shape))

    final_shape = np.asarray(best[1], dtype=FLOAT)
    chi2, alpha, prediction, unit_q = likelihood.profile(model, final_shape, omega_r)
    bound_hits = []
    for index, ((lower, upper), value) in enumerate(zip(bounds, final_shape)):
        tolerance = 1.0e-7 * (1.0 + abs(float(value)))
        if abs(float(value) - lower) <= tolerance:
            bound_hits.append({"shape_parameter": index, "side": "lower", "bound": float(lower)})
        if abs(float(value) - upper) <= tolerance:
            bound_hits.append({"shape_parameter": index, "side": "upper", "bound": float(upper)})
    return {
        "model": model,
        "shape": [float(x) for x in final_shape],
        "parameters": {
            "alpha": float(alpha),
            "Omega_m": ProfileLikelihood.decode(model, final_shape)[0],
            "w0": ProfileLikelihood.decode(model, final_shape)[1],
            "wa": ProfileLikelihood.decode(model, final_shape)[2],
        },
        "Omega_r": float(omega_r),
        "Omega_de_from_closure": float(1.0 - final_shape[0] - omega_r),
        "chi2": float(chi2),
        "prediction": [float(x) for x in prediction],
        "unit_amplitude_prediction": [float(x) for x in unit_q],
        "alpha_profile": {
            "unbounded_solution": float(np.dot(
                solve_triangular(likelihood.chol, unit_q, lower=True, check_finite=False),
                likelihood.white_y,
            ) / np.dot(
                solve_triangular(likelihood.chol, unit_q, lower=True, check_finite=False),
                solve_triangular(likelihood.chol, unit_q, lower=True, check_finite=False),
            )),
            "bounds": [float(ALPHA_BOUNDS[0]), float(ALPHA_BOUNDS[1])],
            "bound_active": bool(alpha in ALPHA_BOUNDS),
        },
        "selected_method": best[2],
        "parameter_bound_hits": bound_hits,
        "optimizers": {
            "bounded_LBFGSB_multistart": main_runs,
            "differential_evolution": de_record,
            "DE_seeded_LBFGSB_refinement": polish_record,
            "boundary_face_Powell_checks": face_checks,
        },
    }


def covariance_shift(likelihood: ProfileLikelihood, base: np.ndarray, changed: np.ndarray, base_chi2: float) -> dict[str, object]:
    delta = np.asarray(changed, dtype=FLOAT) - np.asarray(base, dtype=FLOAT)
    white_delta = solve_triangular(likelihood.chol, delta, lower=True, check_finite=False)
    residual = likelihood.data.value - changed
    white_residual = solve_triangular(likelihood.chol, residual, lower=True, check_finite=False)
    changed_chi2 = FLOAT(np.dot(white_residual, white_residual))
    return {
        "delta_prediction": [float(x) for x in delta],
        "fractional_delta_prediction": [float(x) for x in delta / base],
        "delta_over_marginal_sigma": [float(x) for x in delta / likelihood.sigma],
        "full_covariance_delta_prediction_quadratic": float(np.dot(white_delta, white_delta)),
        "chi2_at_changed_prediction": float(changed_chi2),
        "delta_chi2_at_changed_prediction": float(changed_chi2 - base_chi2),
    }


def param_values(fit: dict[str, object]) -> dict[str, float]:
    pars = fit["parameters"]
    assert isinstance(pars, dict)
    return {name: float(pars[name]) for name in ("alpha", "Omega_m", "w0", "wa")}


def parameter_deltas(new: dict[str, float], old: dict[str, float]) -> dict[str, float]:
    return {name: float(new[name] - old[name]) for name in new}


def convergence_check(likelihood: ProfileLikelihood, fit: dict[str, object]) -> dict[str, object]:
    model = str(fit["model"])
    shape = np.asarray(fit["shape"], dtype=FLOAT)
    omega_r = float(fit["Omega_r"])
    chi96, alpha96, pred96, _ = likelihood.profile(model, shape, omega_r, order=96)
    chi192, alpha192, pred192, _ = likelihood.profile(model, shape, omega_r, order=192)
    omega_m, w0, wa = ProfileLikelihood.decode(model, shape)

    def adaptive_integral(z_value: float) -> float:
        value = quad(
            lambda zp: float(1.0 / np.sqrt(ProfileLikelihood.e2(zp, omega_m, w0, wa, omega_r))),
            0.0,
            float(z_value),
            epsabs=2.0e-13,
            epsrel=2.0e-13,
            limit=200,
        )[0]
        return float(value)

    dm = np.asarray([adaptive_integral(z) for z in likelihood.data.z], dtype=FLOAT)
    dh = 1.0 / np.sqrt(ProfileLikelihood.e2(likelihood.data.z, omega_m, w0, wa, omega_r))
    dv = np.cbrt(likelihood.data.z * dm**2 * dh)
    q = np.empty_like(dm)
    for name, values in (("DM_over_rs", dm), ("DH_over_rs", dh), ("DV_over_rs", dv)):
        q[likelihood.data.observable == name] = values[likelihood.data.observable == name]
    _, alpha_profile192, _, _ = likelihood.profile(model, shape, omega_r, order=192)
    # Re-profile alpha against adaptive distances using the same Cholesky norm.
    white_q = solve_triangular(likelihood.chol, q, lower=True, check_finite=False)
    white_alpha = FLOAT(np.clip(np.dot(white_q, likelihood.white_y) / np.dot(white_q, white_q), *ALPHA_BOUNDS))
    pred_adaptive = white_alpha * q
    residual = likelihood.data.value - pred_adaptive
    white_residual = solve_triangular(likelihood.chol, residual, lower=True, check_finite=False)
    chi_adaptive = float(np.dot(white_residual, white_residual))
    return {
        "quadrature_nodes": {"GL96": 96, "GL192": 192},
        "max_abs_prediction_GL96_vs_GL192": float(np.max(np.abs(pred96 - pred192))),
        "max_abs_prediction_GL96_vs_adaptive_quad": float(np.max(np.abs(pred96 - pred_adaptive))),
        "max_relative_prediction_GL96_vs_adaptive_quad": float(np.max(np.abs((pred96 - pred_adaptive) / pred_adaptive))),
        "chi2_GL96": float(chi96),
        "chi2_GL192": float(chi192),
        "chi2_adaptive_quad": float(chi_adaptive),
        "chi2_GL96_minus_adaptive_quad": float(chi96 - chi_adaptive),
        "max_abs_alpha_GL96_vs_GL192": float(abs(alpha96 - alpha192)),
        "adaptive_profile_alpha": float(white_alpha),
        "profile_alpha_call_reference": float(alpha_profile192),
    }


def main() -> int:
    start_time = time.perf_counter()
    data = bg.BAOData.from_files(MEAN_PATH, COV_PATH)
    likelihood = ProfileLikelihood(data, 96)
    seed_result = json.loads(SEED_RESULT_PATH.read_text(encoding="utf-8"))
    seed_by_model = {entry["model"]: entry for entry in seed_result["models"]}
    radiation = physical_radiation()
    omega_r_by_h0 = {float(h0): float(value) for h0, value in radiation["omega_r_by_h0"].items()}

    baseline = {}
    for model_index, model in enumerate(bg.FIT_BOUNDS):
        baseline_fit = fit_profile(
            likelihood,
            model,
            omega_r=0.0,
            baseline_shape=np.asarray(([0.3, -1.0, 0.0])[: len(shape_bounds(model))], dtype=FLOAT),
            random_seed=FIT_SEED + model_index * 1000,
        )
        expected = seed_by_model[model]
        seed_vector = np.asarray(expected["prediction"], dtype=FLOAT)
        recomputed = bg.predict_bao(
            data.z,
            data.observable,
            baseline_fit["parameters"]["alpha"],
            baseline_fit["parameters"]["Omega_m"],
            baseline_fit["parameters"]["w0"],
            baseline_fit["parameters"]["wa"],
        )
        baseline_fit["seed_comparison"] = {
            "seed_result_sha256": sha256(SEED_RESULT_PATH),
            "seed_chi2": float(expected["chi2"]),
            "reproduced_chi2": float(baseline_fit["chi2"]),
            "chi2_difference": float(baseline_fit["chi2"] - expected["chi2"]),
            "max_abs_prediction_difference": float(np.max(np.abs(recomputed - seed_vector))),
            "baseline_match_tolerance": 1.0e-6,
            "matched": bool(abs(float(baseline_fit["chi2"]) - float(expected["chi2"])) <= 1.0e-6),
        }
        baseline[model] = baseline_fit

    radiation_fits: dict[str, dict[str, object]] = {model: {} for model in bg.FIT_BOUNDS}
    fixed_shifts: dict[str, dict[str, object]] = {model: {} for model in bg.FIT_BOUNDS}
    refitted_shifts: dict[str, dict[str, object]] = {model: {} for model in bg.FIT_BOUNDS}
    convergence = []
    alpha_profile_checks = []
    for h_index, h0 in enumerate(H0_GRID):
        omega_r = omega_r_by_h0[h0]
        for model_index, model in enumerate(bg.FIT_BOUNDS):
            base = baseline[model]
            base_shape = np.asarray(base["shape"], dtype=FLOAT)
            base_params = param_values(base)
            fixed_prediction = float(base_params["alpha"]) * likelihood.unit_prediction(model, base_shape, omega_r)
            fixed_shifts[model][str(h0)] = {
                "H0_km_s_Mpc": h0,
                "Omega_r": omega_r,
                "alpha_held_at_no_radiation_best_fit": float(base_params["alpha"]),
                "shape_held_at_no_radiation_best_fit": [float(x) for x in base_shape],
                **covariance_shift(likelihood, np.asarray(base["prediction"]), fixed_prediction, float(base["chi2"])),
            }
            fit = fit_profile(
                likelihood,
                model,
                omega_r,
                base_shape,
                random_seed=FIT_SEED + 10000 + h_index * 100 + model_index * 1000,
            )
            fit["H0_km_s_Mpc"] = float(h0)
            fit["baseline_shape"] = [float(x) for x in base_shape]
            fit["refitted_parameter_shift_from_no_radiation"] = parameter_deltas(param_values(fit), base_params)
            fit["refitted_prediction_shift_from_no_radiation"] = covariance_shift(
                likelihood,
                np.asarray(base["prediction"], dtype=FLOAT),
                np.asarray(fit["prediction"], dtype=FLOAT),
                float(base["chi2"]),
            )
            fit["refitted_delta_chi2_from_no_radiation"] = float(fit["chi2"] - base["chi2"])
            radiation_fits[model][str(h0)] = fit
            refitted_shifts[model][str(h0)] = {
                "H0_km_s_Mpc": h0,
                "parameter_delta": fit["refitted_parameter_shift_from_no_radiation"],
                "prediction_delta": fit["refitted_prediction_shift_from_no_radiation"],
                "delta_chi2": fit["refitted_delta_chi2_from_no_radiation"],
                "chi2": float(fit["chi2"]),
                "parameters": fit["parameters"],
                "bound_hits": fit["parameter_bound_hits"],
                "selected_method": fit["selected_method"],
            }
            convergence.append({"model": model, "H0_km_s_Mpc": h0, **convergence_check(likelihood, fit)})

            # Check the analytic profile amplitude against an independent scalar minimizer.
            q = np.asarray(fit["unit_amplitude_prediction"], dtype=FLOAT)
            scalar = minimize_scalar(
                lambda alpha: float(np.dot(
                    solve_triangular(likelihood.chol, data.value - alpha * q, lower=True, check_finite=False),
                    solve_triangular(likelihood.chol, data.value - alpha * q, lower=True, check_finite=False),
                )),
                bounds=ALPHA_BOUNDS,
                method="bounded",
                options={"xatol": 1.0e-12, "maxiter": 2000},
            )
            alpha_profile_checks.append({
                "model": model,
                "H0_km_s_Mpc": h0,
                "analytic_alpha": float(fit["parameters"]["alpha"]),
                "bounded_scalar_alpha": float(scalar.x),
                "abs_alpha_difference": float(abs(scalar.x - float(fit["parameters"]["alpha"]))),
                "profile_chi2": float(fit["chi2"]),
                "scalar_chi2": float(scalar.fun),
                "delta_chi2": float(scalar.fun - float(fit["chi2"])),
                "success": bool(scalar.success),
            })

    rankings = {}
    for h0 in H0_GRID:
        rows = [(model, float(radiation_fits[model][str(h0)]["chi2"])) for model in bg.FIT_BOUNDS]
        rankings[str(h0)] = {
            "ordered_lowest_chi2_first": [name for name, _ in sorted(rows, key=lambda row: row[1])],
            "chi2_by_model": {name: value for name, value in rows},
            "improvement_vs_lcdm": {
                model: float(radiation_fits["lcdm"][str(h0)]["chi2"] - radiation_fits[model][str(h0)]["chi2"])
                for model in bg.FIT_BOUNDS
            },
            "incremental_wcdm_improvement_over_lcdm": float(
                radiation_fits["lcdm"][str(h0)]["chi2"] - radiation_fits["wcdm"][str(h0)]["chi2"]
            ),
            "incremental_cpl_improvement_over_wcdm": float(
                radiation_fits["wcdm"][str(h0)]["chi2"] - radiation_fits["cpl"][str(h0)]["chi2"]
            ),
        }

    def envelope(values: list[tuple[float, float]]) -> dict[str, float]:
        return {
            "min": float(min(value for _, value in values)),
            "min_at_H0_km_s_Mpc": float(min(values, key=lambda pair: pair[1])[0]),
            "max": float(max(value for _, value in values)),
            "max_at_H0_km_s_Mpc": float(max(values, key=lambda pair: pair[1])[0]),
        }

    envelope_summary = {}
    for model in bg.FIT_BOUNDS:
        envelope_summary[model] = {
            "Omega_r": envelope([(h0, omega_r_by_h0[h0]) for h0 in H0_GRID]),
            "fixed_baseline_delta_chi2": envelope([(h0, fixed_shifts[model][str(h0)]["delta_chi2_at_changed_prediction"]) for h0 in H0_GRID]),
            "fixed_baseline_max_abs_marginal_sigma_shift": envelope([(h0, max(abs(x) for x in fixed_shifts[model][str(h0)]["delta_over_marginal_sigma"])) for h0 in H0_GRID]),
            "fixed_baseline_full_covariance_quadratic": envelope([(h0, fixed_shifts[model][str(h0)]["full_covariance_delta_prediction_quadratic"]) for h0 in H0_GRID]),
            "refitted_delta_chi2": envelope([(h0, float(radiation_fits[model][str(h0)]["refitted_delta_chi2_from_no_radiation"])) for h0 in H0_GRID]),
            "refitted_max_abs_marginal_sigma_shift": envelope([(h0, max(abs(x) for x in radiation_fits[model][str(h0)]["refitted_prediction_shift_from_no_radiation"]["delta_over_marginal_sigma"])) for h0 in H0_GRID]),
            "refitted_full_covariance_quadratic": envelope([(h0, radiation_fits[model][str(h0)]["refitted_prediction_shift_from_no_radiation"]["full_covariance_delta_prediction_quadratic"]) for h0 in H0_GRID]),
            "refitted_parameter_delta": {
                name: envelope([(h0, float(radiation_fits[model][str(h0)]["refitted_parameter_shift_from_no_radiation"][name])) for h0 in H0_GRID])
                for name in ("alpha", "Omega_m", "w0", "wa")
            },
        }

    zero_rad_checks = []
    for model, fit in baseline.items():
        shape = np.asarray(fit["shape"], dtype=FLOAT)
        omega_m, w0, wa = ProfileLikelihood.decode(model, shape)
        z_check = np.asarray([0.0, 0.05, 0.51, 0.934, 2.33, 3.0], dtype=FLOAT)
        e2_here = ProfileLikelihood.e2(z_check, omega_m, w0, wa, 0.0)
        e2_bg = bg.e2_cpl(z_check, omega_m, w0, wa)
        pred_rad0 = likelihood.profile(model, shape, 0.0)
        pred_bg = bg.predict_bao(data.z, data.observable, fit["parameters"]["alpha"], omega_m, w0, wa)
        zero_rad_checks.append({
            "model": model,
            "max_abs_E2_difference_from_background_bao": float(np.max(np.abs(e2_here - e2_bg))),
            "E0": float(np.sqrt(ProfileLikelihood.e2(np.asarray([0.0]), omega_m, w0, wa, 0.0)[0])),
            "max_abs_prediction_difference_from_background_bao": float(np.max(np.abs(pred_rad0[2] - pred_bg))),
            "profile_chi2": float(pred_rad0[0]),
            "seed_chi2": float(fit["seed_comparison"]["seed_chi2"]),
        })

    nested_baseline_ranking = sorted(bg.FIT_BOUNDS, key=lambda name: float(baseline[name]["chi2"]))
    seed_matches = all(bool(baseline[model]["seed_comparison"]["matched"]) for model in bg.FIT_BOUNDS)
    max_de_gap = max(
        abs(float(fit["chi2"]) - float(fit["optimizers"]["differential_evolution"]["chi2"]))
        for by_h0 in radiation_fits.values()
        for fit in by_h0.values()
    )
    result = {
        "status": "exploratory_reference_model_profile_scan",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {
            "data": "13-row DESI DR2 compressed BAO mean and full covariance in unchanged mean-file order",
            "radiation_model": "reference background with photons plus massless neutrinos at N_eff=3.046; this is not the exact DESI collaboration physical model (e.g. massive neutrinos are not represented)",
            "late_time_domain": "flat background distances over 0<z<=3; no early-universe likelihood or sound-horizon prediction",
            "alpha": "c/(H0*r_d) is an independent profiled BAO amplitude, bounded exactly by [1e-6,1e4]",
            "shape_domains": {model: [[float(lo), float(hi)] for lo, hi in bg.FIT_BOUNDS[model][1:]] for model in bg.FIT_BOUNDS},
            "radiation_envelope": "H0 grid is a sensitivity domain only, not a prior, posterior support or H0 inference",
            "envelope_limitation": "extrema cover the stated finite H0 grid and the reported baseline/refitted optima; the entire shape domain was not exhaustively bounded",
            "claims": "no posterior, evidence, significance or published-constraint reproduction claim",
        },
        "inputs": {
            "mean_path": str(MEAN_PATH.relative_to(ROOT)),
            "mean_sha256": sha256(MEAN_PATH),
            "covariance_path": str(COV_PATH.relative_to(ROOT)),
            "covariance_sha256": sha256(COV_PATH),
            "seed_result_path": str(SEED_RESULT_PATH.relative_to(ROOT)),
            "seed_result_sha256": sha256(SEED_RESULT_PATH),
            "background_bao_path": str((ROOT / "scripts/background_bao.py").relative_to(ROOT)),
            "background_bao_sha256": sha256(ROOT / "scripts/background_bao.py"),
            "mean_order": [{"z": float(z), "observable": str(obs)} for z, obs in zip(data.z, data.observable)],
            "covariance_dimension": int(data.covariance.shape[0]),
            "covariance_symmetric": bool(np.allclose(data.covariance, data.covariance.T, rtol=1e-10, atol=1e-12)),
            "covariance_cholesky_success": True,
        },
        "physical_radiation_derivation": radiation,
        "parameterization": {
            "E2": "Omega_m*(1+z)^3 + Omega_r*(1+z)^4 + Omega_de*f_DE(z)",
            "flat_closure": "Omega_de=1-Omega_m-Omega_r",
            "lcdm_f_DE": "1",
            "wcdm_f_DE": "(1+z)^(3*(1+w0))",
            "cpl_f_DE": "exp(3*(1+w0+wa)*ln(1+z)-3*wa*z/(1+z))",
            "quadrature": "Gauss-Legendre 96 nodes (matching background_bao.py), independently checked with 192 nodes and scipy.integrate.quad",
            "alpha_profile": "clip((q^T C^-1 y)/(q^T C^-1 q),1e-6,1e4); objective formed from whitened residual norm",
            "optimization": "16-start bounded L-BFGS-B over shape parameters; independent bounded differential_evolution with polish disabled; Powell checks on any active best-fit face; shape bounds imported unchanged from background_bao.FIT_BOUNDS",
            "float_format": "NumPy float64 throughout; no GPU",
        },
        "execution": {
            "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --seconds 1200 -- .venv/bin/python experiments/radiation_sensitivity/radiation_scan.py",
            "random_seed_base": FIT_SEED,
            "baseline_seed_reference": BASELINE_SEED,
            "optimizer_starts_per_fit": 16,
            "optimizer_statuses_retained": True,
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_affinity_count": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
            "thread_environment": {name: os.environ.get(name) for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
            "runtime_seconds": None,
        },
        "h0_grid_km_s_mpc": list(H0_GRID),
        "no_radiation_baseline": baseline,
        "baseline_model_order_lowest_chi2_first": nested_baseline_ranking,
        "fixed_baseline_radiation_shifts": fixed_shifts,
        "radiation_refits": radiation_fits,
        "refitted_shifts": refitted_shifts,
        "model_ordering_and_improvements_by_h0": rankings,
        "h0_grid_envelopes": envelope_summary,
        "checks": {
            "all_seed_baseline_chi2_within_1e-6": bool(seed_matches),
            "zero_radiation_reduction_to_background_bao": zero_rad_checks,
            "profiled_alpha_vs_bounded_scalar_minimizer": alpha_profile_checks,
            "max_abs_profile_vs_scalar_alpha_difference": float(max(x["abs_alpha_difference"] for x in alpha_profile_checks)),
            "max_abs_profile_vs_scalar_delta_chi2": float(max(abs(x["delta_chi2"]) for x in alpha_profile_checks)),
            "quadrature_convergence_at_refitted_optima": convergence,
            "max_abs_GL96_vs_adaptive_prediction": float(max(x["max_abs_prediction_GL96_vs_adaptive_quad"] for x in convergence)),
            "max_abs_GL96_vs_GL192_prediction": float(max(x["max_abs_prediction_GL96_vs_GL192"] for x in convergence)),
            "max_abs_GL96_vs_adaptive_delta_chi2": float(max(abs(x["chi2_GL96_minus_adaptive_quad"]) for x in convergence)),
            "max_profile_fit_vs_independent_DE_delta_chi2": float(max_de_gap),
            "all_best_fit_de_checks_within_1e-6": bool(max_de_gap <= 1.0e-6),
            "analytic_E0_is_one_for_baseline_closure": bool(all(abs(x["E0"] - 1.0) <= 2.0e-15 for x in zero_rad_checks)),
            "tested_domain_includes_H0_grid_only": True,
        },
    }
    result["execution"]["runtime_seconds"] = float(time.perf_counter() - start_time)
    result["code_sha256"] = {
        str(HERE.relative_to(ROOT)): sha256(HERE),
        "scripts/background_bao.py": sha256(ROOT / "scripts/background_bao.py"),
    }
    OUTPUT_PATH.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "runtime_seconds": result["execution"]["runtime_seconds"],
        "omega_r_h2": radiation["omega_r_h2"],
        "omega_gamma_h2": radiation["omega_gamma_h2"],
        "baseline_model_order": nested_baseline_ranking,
        "baseline_chi2": {model: baseline[model]["chi2"] for model in bg.FIT_BOUNDS},
        "baseline_seed_matches": seed_matches,
        "rad_fit_chi2_by_H0": {
            model: {h0: radiation_fits[model][str(h0)]["chi2"] for h0 in H0_GRID}
            for model in bg.FIT_BOUNDS
        },
        "envelopes": envelope_summary,
        "checks": {
            "max_abs_GL96_vs_adaptive_prediction": result["checks"]["max_abs_GL96_vs_adaptive_prediction"],
            "max_abs_GL96_vs_GL192_prediction": result["checks"]["max_abs_GL96_vs_GL192_prediction"],
            "max_profile_fit_vs_independent_DE_delta_chi2": max_de_gap,
        },
        "output": str(OUTPUT_PATH.relative_to(ROOT)),
    }, indent=2, allow_nan=False))
    return 0 if seed_matches and max_de_gap <= 1.0e-6 else 2


if __name__ == "__main__":
    sys.exit(main())

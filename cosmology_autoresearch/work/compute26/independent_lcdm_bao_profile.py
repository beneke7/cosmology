#!/usr/bin/env python3
"""Independent CAMB 2.0.4 flat-LambdaCDM / DESI DR2 BAO profile.

This is a BAO-only finite-box profile, not a posterior and not a CMB-likelihood
reproduction.  It fits H0 and physical baryon/CDM densities directly.  This
avoids nested theta_MC root failures during broad-box exploration; the winning
point is also mapped back through CAMB's theta_MC setter with an expanded H0
root interval to audit the original runner's [10, 100] km/s/Mpc bracket.

The density bounds are the broad Table-I limits in the local, hash-pinned
Yang et al. DR2 paper (p. 10, Table I); H0 has an explicit broad stress box.
All score calculations use float64 and the complete released 13x13 covariance.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
from datetime import datetime, timezone

for _thread_var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_thread_var] = "1"

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
PINNED_SITE = PROJECT / "work/compute24/site"
sys.path.insert(0, str(PINNED_SITE))

import camb  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402
from scipy.integrate import quad  # noqa: E402
from scipy.optimize import brentq, differential_evolution  # noqa: E402

MEAN_PATH = PROJECT / "context/data/desi_dr2_mean.txt"
COV_PATH = PROJECT / "context/data/desi_dr2_cov.txt"
PAPER_PATH = PROJECT / "context/papers/interacting_de_desi_dr2_2026.pdf"
WHEEL_PATH = PROJECT / "context/code/camb-2.0.4-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl"
SOURCE_PATH = PROJECT / "context/code/CAMB_2.0.4_a6de8cc59124c8bbe3924f1c546f96cef73dcabb.tar.gz"
OUTPUT_PATH = HERE / "independent_lcdm_bao_profile.json"
PLOT_PATH = HERE / "independent_lcdm_bao_profile.png"

# The omega_b h^2 and omega_c h^2 ranges reproduce the deliberately broad
# finite ranges printed in Table I; H0's [10,1000] box makes the direct-H0
# reparameterization much wider than ordinary late-time measurements.
PARAMETER_NAMES = ("H0_km_s_Mpc", "omega_b_h2", "omega_c_h2")
BOUNDS = ((10.0, 1000.0), (0.005, 0.1), (0.001, 0.99))
FIXED = {
    "mnu_eV": 0.06,
    "N_eff": 3.044,
    "T_CMB_K": 2.7255,
    "num_massive_neutrinos": 1,
    "neutrino_hierarchy": "degenerate",
    "YHe": "CAMB BBN consistency default (YHe=None)",
    "tau_reio": 0.054,
    "As": 2.1e-9,
    "ns": 0.965,
    "Omega_k": 0.0,
    "w": -1.0,
    "wa": 0.0,
}
C_KM_S = np.float64(299792.458)
INVALID_SCORE = np.float64(1.0e80)
SEEDS = (20260925, 20260926, 20260927, 20260928)
MAXITER = 600
POPULATION_SIZE = 48
DE_TOL = 1.0e-9
DE_ATOL = 1.0e-10
SAME_MINIMUM_SCALED_LINF_TOL = 2.0e-3
SAME_MINIMUM_CHI2_TOL = 2.0e-5
THETA100_PAPER_RANGE = (0.5, 10.0)
ORIGINAL_THETA_H0_RANGE = (10.0, 100.0)
EXPANDED_THETA_H0_RANGE = (1.0, 2000.0)

_INVALID_COUNTS: dict[str, int] = {}
_INVALID_EXAMPLES: list[dict[str, object]] = []


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_bao() -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    records: list[tuple[float, float, str]] = []
    for row in MEAN_PATH.read_text(encoding="utf-8").splitlines():
        row = row.strip()
        if row and not row.startswith("#"):
            z_s, mean_s, observable = row.split()
            records.append((float(z_s), float(mean_s), observable))
    z = np.asarray([row[0] for row in records], dtype=np.float64)
    mean = np.asarray([row[1] for row in records], dtype=np.float64)
    observable = [row[2] for row in records]
    covariance = np.asarray(np.loadtxt(COV_PATH, dtype=np.float64), dtype=np.float64)
    if z.shape != (13,) or mean.shape != (13,) or covariance.shape != (13, 13):
        raise ValueError(f"expected 13 rows and a 13x13 covariance; got {z.shape}, {mean.shape}, {covariance.shape}")
    if not np.isfinite(covariance).all() or not np.allclose(covariance, covariance.T, rtol=0.0, atol=1e-12):
        raise ValueError("covariance is non-finite or asymmetric")
    np.linalg.cholesky(covariance)
    return z, mean, covariance, observable


Z, Y, COV, OBSERVABLE = load_bao()
CHOLESKY = np.linalg.cholesky(COV)


def make_params(h0: float, omega_b: float, omega_c: float) -> camb.CAMBparams:
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
    pars.set_cosmology(
        H0=float(h0),
        ombh2=float(omega_b),
        omch2=float(omega_c),
        omk=0.0,
        mnu=FIXED["mnu_eV"],
        nnu=FIXED["N_eff"],
        num_massive_neutrinos=FIXED["num_massive_neutrinos"],
        neutrino_hierarchy=FIXED["neutrino_hierarchy"],
        TCMB=FIXED["T_CMB_K"],
        YHe=None,
        tau=FIXED["tau_reio"],
    )
    pars.InitPower.set_params(As=FIXED["As"], ns=FIXED["ns"], r=0.0)
    return pars


def invalid(reason: str, point: np.ndarray, detail: str | None = None) -> np.float64:
    _INVALID_COUNTS[reason] = _INVALID_COUNTS.get(reason, 0) + 1
    if len(_INVALID_EXAMPLES) < 24:
        example: dict[str, object] = {"parameters": np.asarray(point, dtype=np.float64).tolist()}
        if detail:
            example["detail"] = detail[:240]
        _INVALID_EXAMPLES.append({"reason": reason, **example})
    return INVALID_SCORE


def model(point: np.ndarray | tuple[float, float, float] | list[float]) -> dict[str, object]:
    h0, omega_b, omega_c = (float(v) for v in point)
    if not all(math.isfinite(v) for v in (h0, omega_b, omega_c)):
        raise ValueError("non-finite parameter")
    pars = make_params(h0, omega_b, omega_c)
    if not (0.0 < pars.omegam < 1.0):
        raise ValueError(f"non-positive flat-Lambda density implied by CAMB Omega_m={pars.omegam:.16g}")
    background = camb.get_background(pars)
    omega_de = float(background.get_Omega("de", 0.0))
    if not math.isfinite(omega_de) or omega_de <= 0.0:
        raise ValueError(f"non-positive CAMB dark-energy density Omega_de={omega_de:.16g}")
    rdrag = float(background.get_derived_params()["rdrag"])
    if not math.isfinite(rdrag) or rdrag <= 0.0:
        raise ValueError(f"non-positive CAMB r_drag={rdrag}")

    hz = np.asarray(background.hubble_parameter(Z), dtype=np.float64)
    dm = np.asarray(background.comoving_radial_distance(Z), dtype=np.float64)
    dh = C_KM_S / hz
    prediction = np.empty(13, dtype=np.float64)
    for i, name in enumerate(OBSERVABLE):
        if name == "DM_over_rs":
            prediction[i] = dm[i] / rdrag
        elif name == "DH_over_rs":
            prediction[i] = dh[i] / rdrag
        elif name == "DV_over_rs":
            prediction[i] = np.power(Z[i] * dm[i] * dm[i] * dh[i], 1.0 / 3.0) / rdrag
        else:
            raise ValueError(f"unknown observable {name}")
    residual = prediction - Y
    whitened = np.linalg.solve(CHOLESKY, residual)
    chi2_cholesky = float(whitened @ whitened)
    chi2_dense = float(residual @ np.linalg.solve(COV, residual))
    if not np.isfinite(chi2_cholesky) or not np.isfinite(chi2_dense):
        raise ArithmeticError("non-finite score")
    if abs(chi2_cholesky - chi2_dense) > 2e-10 * max(1.0, abs(chi2_dense)):
        raise ArithmeticError(f"dense/Cholesky score mismatch: {chi2_dense} vs {chi2_cholesky}")
    derived = background.get_derived_params()
    return {
        "parameters": dict(zip(PARAMETER_NAMES, [h0, omega_b, omega_c])),
        "derived": {
            "Omega_m0_including_massive_nu": float(pars.omegam),
            "Omega_de0_CAMB": omega_de,
            "Omega_b0": float(pars.omegab),
            "Omega_c0": float(pars.omegac),
            "Omega_nu0": float(pars.omeganu),
            "r_drag_Mpc": rdrag,
            "100theta_star_CAMB": float(derived["thetastar"]),
            "z_drag": float(derived["zdrag"]),
        },
        "bao": {
            "chi2_full_covariance_cholesky": chi2_cholesky,
            "chi2_full_covariance_dense_solve": chi2_dense,
            "chi2_difference_abs": abs(chi2_cholesky - chi2_dense),
            "predictions_in_input_order": prediction.tolist(),
            "residuals_in_input_order": residual.tolist(),
            "cholesky_whitened_residuals": whitened.tolist(),
            "H_km_s_Mpc_at_rows": hz.tolist(),
            "DM_comoving_Mpc_at_rows": dm.tolist(),
        },
        "_objects": {"params": pars, "background": background},
    }


def objective(point: np.ndarray) -> float:
    point64 = np.asarray(point, dtype=np.float64)
    try:
        # Bounds are enforced by SciPy; explicit guard also protects local-polish calls.
        if any(value < lo or value > hi for value, (lo, hi) in zip(point64, BOUNDS)):
            return float(invalid("outside_declared_bounds", point64))
        result = model(point64)
        return float(result["bao"]["chi2_full_covariance_cholesky"])
    except ArithmeticError as exc:
        return float(invalid("arithmetic_or_score", point64, str(exc)))
    except Exception as exc:  # CAMB and flat-closure failures are invalid trial points.
        return float(invalid(type(exc).__name__, point64, str(exc)))


def feasible_population(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    population: list[np.ndarray] = []
    anchors = (
        np.asarray([67.057, 0.0224, 0.1200], dtype=np.float64),
        np.asarray([99.159, 0.06069, 0.23096], dtype=np.float64),
    )
    for candidate in anchors:
        if BOUNDS[0][0] < candidate[0] < BOUNDS[0][1] and (candidate[1] + candidate[2] + 0.0007) / (candidate[0] / 100.0) ** 2 < 0.98:
            population.append(candidate)
    attempts = 0
    while len(population) < POPULATION_SIZE and attempts < 100_000:
        attempts += 1
        candidate = np.asarray([rng.uniform(lo, hi) for lo, hi in BOUNDS], dtype=np.float64)
        h0, omega_b, omega_c = candidate
        # This conservative prefilter only builds an initially physical
        # population; the objective checks CAMB's actual closure exactly.
        if (omega_b + omega_c + 0.0007) / (h0 / 100.0) ** 2 < 0.98:
            population.append(candidate)
    if len(population) != POPULATION_SIZE:
        raise RuntimeError(f"could construct only {len(population)} physically interior starting points")
    return np.asarray(population, dtype=np.float64)


def independent_distance_check(best_parameters: np.ndarray, best: dict[str, object]) -> dict[str, object]:
    bg = best["_objects"]["background"]
    errors = []
    integrals = []
    for z, dm in zip(Z, best["bao"]["DM_comoving_Mpc_at_rows"]):
        integral, _ = quad(
            lambda z_eval: C_KM_S / float(bg.hubble_parameter(z_eval)),
            0.0,
            float(z),
            epsabs=1e-9,
            epsrel=1e-11,
            limit=250,
        )
        integrals.append(float(integral))
        errors.append(abs(float(dm) - float(integral)))
    return {
        "method": "SciPy adaptive Gauss-Kronrod quadrature of c/H(z) independently at all 13 redshifts",
        "c_km_s": float(C_KM_S),
        "integral_D_M_Mpc_in_input_order": integrals,
        "abs_CAMB_minus_integral_Mpc_in_input_order": errors,
        "maximum_absolute_error_Mpc": float(max(errors)),
        "all_rows_within_1e-6_Mpc": bool(max(errors) < 1e-6),
        "dense_vs_cholesky_abs_score_difference": best["bao"]["chi2_difference_abs"],
    }


def theta_mapping_diagnostic(best_parameters: np.ndarray) -> dict[str, object]:
    target_h0, omega_b, omega_c = (float(v) for v in best_parameters)

    def h0_from_theta(theta100: float, h0_range: tuple[float, float]) -> float:
        pars = camb.CAMBparams()
        pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
        pars.set_cosmology(
            cosmomc_theta=float(theta100) / 100.0,
            ombh2=omega_b,
            omch2=omega_c,
            omk=0.0,
            mnu=FIXED["mnu_eV"],
            nnu=FIXED["N_eff"],
            num_massive_neutrinos=FIXED["num_massive_neutrinos"],
            neutrino_hierarchy=FIXED["neutrino_hierarchy"],
            TCMB=FIXED["T_CMB_K"],
            YHe=None,
            tau=FIXED["tau_reio"],
            theta_H0_range=h0_range,
        )
        return float(pars.H0)

    theta_nodes = np.linspace(0.6, 2.0, 71, dtype=np.float64)
    sampled: list[tuple[float, float]] = []
    old_range_solutions: list[tuple[float, float]] = []
    for theta100 in theta_nodes:
        try:
            h0_expanded = h0_from_theta(float(theta100), EXPANDED_THETA_H0_RANGE)
            sampled.append((float(theta100), h0_expanded))
        except Exception:
            continue
        try:
            old_h0 = h0_from_theta(float(theta100), ORIGINAL_THETA_H0_RANGE)
            old_range_solutions.append((float(theta100), old_h0))
        except Exception:
            pass
    bracket = None
    for left, right in zip(sampled, sampled[1:]):
        if (left[1] - target_h0) * (right[1] - target_h0) <= 0.0:
            bracket = (left[0], right[0])
            break
    inferred: dict[str, object]
    if bracket is None:
        inferred = {"status": "no bracket found on theta100 scan [0.6,2.0]", "theta100_scan_points": len(sampled)}
    else:
        theta100 = float(
            brentq(
                lambda theta: h0_from_theta(theta, EXPANDED_THETA_H0_RANGE) - target_h0,
                bracket[0],
                bracket[1],
                xtol=2e-12,
                rtol=1e-13,
            )
        )
        reproduced_h0 = h0_from_theta(theta100, EXPANDED_THETA_H0_RANGE)
        try:
            old_range_h0 = h0_from_theta(theta100, ORIGINAL_THETA_H0_RANGE)
            old_range_status = "root exists in original [10,100] bracket"
        except Exception as exc:
            old_range_h0 = None
            old_range_status = f"original [10,100] bracket fails: {type(exc).__name__}: {exc}"
        inferred = {
            "status": "rooted back to same theta_MC input via CAMB set_cosmology",
            "implied_100theta_MC": theta100,
            "implied_100theta_MC_inside_paper_Table_I_range_0p5_to_10": bool(THETA100_PAPER_RANGE[0] <= theta100 <= THETA100_PAPER_RANGE[1]),
            "theta_solver_H0_range_km_s_Mpc": list(EXPANDED_THETA_H0_RANGE),
            "H0_recovered_km_s_Mpc": reproduced_h0,
            "H0_recovery_abs_error_km_s_Mpc": abs(reproduced_h0 - target_h0),
            "original_theta_solver_range_km_s_Mpc": list(ORIGINAL_THETA_H0_RANGE),
            "original_range_status": old_range_status,
            "original_range_H0_if_rooted": old_range_h0,
            "distance_to_original_low_edge_km_s_Mpc": target_h0 - ORIGINAL_THETA_H0_RANGE[0],
            "distance_to_original_high_edge_km_s_Mpc": ORIGINAL_THETA_H0_RANGE[1] - target_h0,
            "expanded_root_range_edge_distances_km_s_Mpc": [target_h0 - EXPANDED_THETA_H0_RANGE[0], EXPANDED_THETA_H0_RANGE[1] - target_h0],
        }
    return {
        "interpretation": "The independent search varies H0 directly. This inverse check establishes the CAMB cosmomc_theta mapping and shows whether the old root bracket clips the fitted physical point.",
        "expanded_range_theta_samples_with_root": len(sampled),
        "original_range_theta_samples_with_root": len(old_range_solutions),
        "implied_theta_mapping": inferred,
    }


def boundary_diagnostic(parameters: np.ndarray) -> dict[str, object]:
    margins = {}
    for name, value, (lo, hi) in zip(PARAMETER_NAMES, parameters, BOUNDS):
        margins[name] = {
            "lower_bound": lo,
            "value": float(value),
            "upper_bound": hi,
            "absolute_distance_to_lower": float(value - lo),
            "absolute_distance_to_upper": float(hi - value),
            "normalized_nearest_edge_margin": float(min(value - lo, hi - value) / (hi - lo)),
        }
    return {
        "parameter_bounds": margins,
        "minimum_normalized_edge_margin": float(min(row["normalized_nearest_edge_margin"] for row in margins.values())),
        "candidate_is_at_least_1_percent_from_every_parameter_edge": bool(min(row["normalized_nearest_edge_margin"] for row in margins.values()) >= 0.01),
    }


def baryon_edge_stress(best_parameters: np.ndarray) -> dict[str, object]:
    h0, _omega_b, omega_c = (float(v) for v in best_parameters)
    trial_values = (0.06, 0.08, 0.09, 0.095, 0.099, 0.1, 0.1001, 0.102, 0.105, 0.11, 0.12, 0.15, 0.2)
    rows = []
    for omega_b in trial_values:
        trial = model(np.asarray([h0, omega_b, omega_c], dtype=np.float64))
        rows.append(
            {
                "omega_b_h2": omega_b,
                "inside_declared_Table_I_range": bool(BOUNDS[1][0] <= omega_b <= BOUNDS[1][1]),
                "chi2_fixed_H0_and_omega_c": trial["bao"]["chi2_full_covariance_cholesky"],
                "Omega_m0": trial["derived"]["Omega_m0_including_massive_nu"],
                "r_drag_Mpc": trial["derived"]["r_drag_Mpc"],
            }
        )
    return {
        "purpose": "One-dimensional diagnostic with H0 and omega_c h2 held at the selected box minimum; points above 0.1 are outside the declared Table-I box and are not candidate fits.",
        "rows": rows,
        "chi2_at_upper_edge": next(row["chi2_fixed_H0_and_omega_c"] for row in rows if row["omega_b_h2"] == 0.1),
        "chi2_0p0001_above_upper_edge": next(row["chi2_fixed_H0_and_omega_c"] for row in rows if row["omega_b_h2"] == 0.1001),
    }


def same_interior_minimum(first: dict[str, object], second: dict[str, object]) -> dict[str, object]:
    x1 = np.asarray(first["polished_parameters"], dtype=np.float64)
    x2 = np.asarray(second["polished_parameters"], dtype=np.float64)
    width = np.asarray([hi - lo for lo, hi in BOUNDS], dtype=np.float64)
    scaled_linf = float(np.max(np.abs(x1 - x2) / width))
    delta_chi2 = abs(float(first["chi2"]) - float(second["chi2"]))
    return {
        "seed_pair": [first["seed"], second["seed"]],
        "scaled_parameter_Linf_distance": scaled_linf,
        "absolute_chi2_difference": delta_chi2,
        "scaled_parameter_tolerance": SAME_MINIMUM_SCALED_LINF_TOL,
        "chi2_tolerance": SAME_MINIMUM_CHI2_TOL,
        "same_minimum_within_declared_tolerance": bool(scaled_linf <= SAME_MINIMUM_SCALED_LINF_TOL and delta_chi2 <= SAME_MINIMUM_CHI2_TOL),
    }


def save_plot(fiducial: dict[str, object], selected: dict[str, object]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.arange(len(Z))
    sigma = np.sqrt(np.diag(COV))
    labels = [f"{z:g}\n{kind.replace('_over_rs', '')}" for z, kind in zip(Z, OBSERVABLE)]
    fig, (ax, resid_ax) = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True, gridspec_kw={"height_ratios": [2, 1]})
    ax.errorbar(x, Y, yerr=sigma, fmt="ko", ms=4, capsize=2, label="DESI DR2 mean (marginal σ)")
    ax.plot(x, fiducial["bao"]["predictions_in_input_order"], "s--", color="#377eb8", label="CAMB flat ΛCDM fiducial")
    ax.plot(x, selected["bao"]["predictions_in_input_order"], "o-", color="#e41a1c", label="Independent BAO profile")
    ax.set_ylabel("BAO distance / r_drag")
    ax.set_xticks(x, labels)
    ax.legend(frameon=False, ncol=3)
    ax.grid(alpha=0.25)
    resid_ax.axhline(0.0, color="k", lw=0.8)
    resid_ax.bar(x, selected["bao"]["cholesky_whitened_residuals"], color="#984ea3")
    resid_ax.set_ylabel("Cholesky-whitened residual")
    resid_ax.set_xticks(x, labels)
    resid_ax.grid(axis="y", alpha=0.25)
    fig.suptitle("CAMB 2.0.4 flat ΛCDM × DESI DR2 BAO (single-probe profile)")
    fig.savefig(PLOT_PATH, dpi=170)
    plt.close(fig)


def main() -> int:
    started_utc = datetime.now(timezone.utc).isoformat()
    start_clock = time.perf_counter()
    fiducial_parameters = np.asarray([67.0571399255, 0.0224, 0.12], dtype=np.float64)
    fiducial = model(fiducial_parameters)
    runs: list[dict[str, object]] = []
    for seed in SEEDS:
        run_start = time.perf_counter()
        initial = feasible_population(seed)
        fit = differential_evolution(
            objective,
            bounds=BOUNDS,
            init=initial,
            strategy="best1bin",
            maxiter=MAXITER,
            popsize=POPULATION_SIZE // len(BOUNDS),
            tol=DE_TOL,
            atol=DE_ATOL,
            seed=seed,
            polish=True,
            updating="immediate",
            workers=1,
            disp=False,
        )
        polished = np.asarray(fit.x, dtype=np.float64)
        # Explicit reevaluation makes optimizer sentinel values insufficient
        # to masquerade as an accepted candidate.
        try:
            checked = model(polished)
            checked_chi2 = float(checked["bao"]["chi2_full_covariance_cholesky"])
            checked_valid = math.isfinite(checked_chi2) and checked_chi2 < float(INVALID_SCORE)
        except Exception as exc:
            checked_chi2 = float(INVALID_SCORE)
            checked_valid = False
            checked_error = f"{type(exc).__name__}: {exc}"
        runs.append(
            {
                "seed": seed,
                "differential_evolution_success": bool(fit.success),
                "differential_evolution_message": str(fit.message),
                "nfev": int(fit.nfev),
                "nit": int(fit.nit),
                "runtime_seconds": time.perf_counter() - run_start,
                "chi2": checked_chi2,
                "optimizer_reported_chi2": float(fit.fun),
                "optimizer_parameters": np.asarray(fit.x, dtype=np.float64).tolist(),
                "polished_parameters": polished.tolist(),
                "checked_candidate_valid": checked_valid,
                **({"checked_candidate_error": checked_error} if not checked_valid and "checked_error" in locals() else {}),
                "population_energy_min": float(np.min(fit.population_energies)),
                "population_energy_median": float(np.median(fit.population_energies)),
                "population_energy_max": float(np.max(fit.population_energies)),
                "accepted_as_converged_seed": bool(fit.success and checked_valid),
            }
        )
        if "checked_error" in locals():
            del checked_error

    valid_runs = [row for row in runs if row["checked_candidate_valid"] and row["differential_evolution_success"]]
    pair_checks = [same_interior_minimum(valid_runs[i], valid_runs[j]) for i in range(len(valid_runs)) for j in range(i + 1, len(valid_runs))]
    matching_pairs = [row for row in pair_checks if row["same_minimum_within_declared_tolerance"]]
    converged = bool(matching_pairs)
    candidate_run = min((row for row in runs if row["checked_candidate_valid"]), key=lambda row: float(row["chi2"]))
    best_parameters = np.asarray(candidate_run["polished_parameters"], dtype=np.float64)
    best = model(best_parameters)
    clean_best = {key: value for key, value in best.items() if key != "_objects"}
    distance = independent_distance_check(best_parameters, best)
    theta_diagnostic = theta_mapping_diagnostic(best_parameters)
    boundary = boundary_diagnostic(best_parameters)
    save_plot(fiducial, best)

    all_scores_finite = all(math.isfinite(float(row["chi2"])) and float(row["chi2"]) < float(INVALID_SCORE) for row in runs)
    success_criteria = {
        "at_least_two_valid_DE_converged_seeds": len(valid_runs) >= 2,
        "at_least_two_converged_seeds_same_minimum_within_tolerances": converged,
        "minimum_pairwise_coordinate_tolerance_scaled_Linf": SAME_MINIMUM_SCALED_LINF_TOL,
        "minimum_pairwise_chi2_tolerance": SAME_MINIMUM_CHI2_TOL,
        "candidate_interior_by_1_percent_margin": boundary["candidate_is_at_least_1_percent_from_every_parameter_edge"],
        "all_seed_candidate_scores_are_valid_finite_chi2": all_scores_finite,
        "independent_distance_integrals_pass_1e-6_Mpc": distance["all_rows_within_1e-6_Mpc"],
        "dense_and_Cholesky_scores_agree": distance["dense_vs_cholesky_abs_score_difference"] <= 1e-10,
    }
    accepted = all(success_criteria.values())
    failed_seeds = [
        {"seed": row["seed"], "success": row["differential_evolution_success"], "valid": row["checked_candidate_valid"], "message": row["differential_evolution_message"], "chi2": row["chi2"]}
        for row in runs
        if not row["accepted_as_converged_seed"]
    ]
    ended_utc = datetime.now(timezone.utc).isoformat()
    status = "reproduced_interior_profile_minimum" if accepted else (
        "incomplete_boundary_pinned_profile" if matching_pairs and len(valid_runs) >= 2 else "incomplete_nonconverged_profile"
    )
    acceptance_interpretation = (
        "The strict profile gate passed."
        if accepted
        else "The best valid seeds agree, but the selected point is on a declared parameter boundary; it is retained only as an exploratory boundary profile and is not accepted as an interior minimum."
        if matching_pairs and len(valid_runs) >= 2 and not boundary["candidate_is_at_least_1_percent_from_every_parameter_edge"]
        else "The strict profile gate did not pass; no candidate is accepted."
    )
    report = {
        "schema_version": 1,
        "status": status,
        "candidate_acceptance": {
            "accepted": accepted,
            "interpretation": acceptance_interpretation,
            "rule": "At least two distinct valid DE seeds must report convergence to an interior point agreeing within the declared scaled-coordinate and chi2 tolerances; all numerical checks must pass.",
            "successful_valid_seed_count": len(valid_runs),
            "matching_seed_pairs": matching_pairs,
            "all_pairwise_checks_among_successful_valid_seeds": pair_checks,
            "failed_seeds": failed_seeds,
            "success_criteria": success_criteria,
        },
        "scope": {
            "model": "flat standard ΛCDM / general relativity with CAMB 2.0.4 background",
            "observational_likelihood": "DESI DR2 released 13-row Gaussian BAO mean and full 13x13 covariance only",
            "free_parameters": list(PARAMETER_NAMES),
            "fixed_physics": FIXED,
            "not_included": ["Planck likelihood", "ACT likelihood", "supernova likelihood", "dark-sector interaction", "posterior or evidence", "cross-survey likelihood"],
            "result_interpretation": "An optimizer profile in an explicit finite box; not a posterior and not a full Planck/ACT or target-paper likelihood reproduction.",
            "H0_domain_note": "H0 is varied directly to avoid CAMB theta_MC root-range failures in the broad search. The candidate's implied 100theta_MC is recovered after fitting. H0 is a computational profile coordinate, not a prior on H0.",
        },
        "sources": {
            "DESI_DR2_mean": {
                "path": str(MEAN_PATH.relative_to(PROJECT)),
                "sha256": sha256(MEAN_PATH),
                "official_upstream": "CobayaSampler/bao_data commit bb0c1c9009dc76d1391300e169e8df38fd1096db, desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",
            },
            "DESI_DR2_covariance": {
                "path": str(COV_PATH.relative_to(PROJECT)),
                "sha256": sha256(COV_PATH),
                "official_upstream": "CobayaSampler/bao_data commit bb0c1c9009dc76d1391300e169e8df38fd1096db, desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt",
            },
            "CAMB_wheel": {
                "path": str(WHEEL_PATH.relative_to(PROJECT)),
                "sha256": sha256(WHEEL_PATH),
                "version": camb.__version__,
                "repository_commit": "a6de8cc59124c8bbe3924f1c546f96cef73dcabb",
            },
            "CAMB_source_archive": {"path": str(SOURCE_PATH.relative_to(PROJECT)), "sha256": sha256(SOURCE_PATH)},
            "target_paper_prior_source": {
                "path": str(PAPER_PATH.relative_to(PROJECT)),
                "sha256": sha256(PAPER_PATH),
                "locator": "Yang et al., arXiv:2609.05660 v1, p. 10, Table I: omega_b h^2=[0.005,0.1], omega_c h^2=[0.001,0.99], 100theta_MC=[0.5,10].",
            },
            "input_row_order": [
                {"z": float(z), "observable": name, "mean": float(y)}
                for z, name, y in zip(Z, OBSERVABLE, Y)
            ],
        },
        "optimization": {
            "bounds": {name: list(bound) for name, bound in zip(PARAMETER_NAMES, BOUNDS)},
            "bounds_basis": "Omega_b h2 and Omega_c h2 use the target paper's very broad Table-I limits; H0 in [10,1000] km/s/Mpc is an explicit broad search coordinate after removing the fragile nested theta root.",
            "finite_box_is_not_a_posterior_prior": True,
            "seeds": list(SEEDS),
            "differential_evolution": {
                "strategy": "best1bin",
                "maxiter": MAXITER,
                "initial_population_size": POPULATION_SIZE,
                "tol": DE_TOL,
                "atol": DE_ATOL,
                "updating": "immediate",
                "polish": True,
                "workers": 1,
                "thread_pools": 1,
                "float_dtype": "IEEE-754 binary64 / NumPy float64",
            },
            "same_minimum_tolerances": {"scaled_parameter_Linf": SAME_MINIMUM_SCALED_LINF_TOL, "chi2": SAME_MINIMUM_CHI2_TOL},
            "seed_results": runs,
            "failed_seeds": failed_seeds,
            "invalid_objective_points": {"score_sentinel": float(INVALID_SCORE), "counts_by_reason": _INVALID_COUNTS, "first_examples_capped_at_24": _INVALID_EXAMPLES},
            "candidate_boundary_diagnostic": boundary,
            "theta_mc_root_diagnostic": theta_diagnostic,
        },
        "fiducial": {key: value for key, value in fiducial.items() if key != "_objects"},
        "selected_candidate": clean_best,
        "independent_checks": {
            "c_over_h_integral": distance,
            "baryon_upper_edge_stress": baryon_edge_stress(best_parameters),
            "full_covariance_dense_vs_cholesky": {
                "cholesky_score": best["bao"]["chi2_full_covariance_cholesky"],
                "dense_score": best["bao"]["chi2_full_covariance_dense_solve"],
                "absolute_difference": best["bao"]["chi2_difference_abs"],
                "pass_tolerance": 1e-10,
            },
        },
        "execution": {
            "started_utc": started_utc,
            "finished_utc": ended_utc,
            "runtime_seconds": time.perf_counter() - start_clock,
            "command": "cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute26/task_state.json --seconds 540 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute26/independent_lcdm_bao_profile.py",
            "runtime_cap_seconds": 540,
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "CAMB": camb.__version__,
            "pinned_CAMB_path": str(PINNED_SITE.relative_to(PROJECT)),
            "code_sha256": sha256(Path(__file__)),
            "worker_pool_size": 1,
        },
        "plot": str(PLOT_PATH.relative_to(PROJECT)),
    }
    OUTPUT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "accepted": accepted,
                "selected_chi2": best["bao"]["chi2_full_covariance_cholesky"],
                "parameters": best["parameters"],
                "derived": best["derived"],
                "successful_valid_seeds": len(valid_runs),
                "matching_pairs": matching_pairs,
                "failed_seeds": failed_seeds,
                "invalid_counts": _INVALID_COUNTS,
                "runtime_seconds": report["execution"]["runtime_seconds"],
                "result": str(OUTPUT_PATH.relative_to(PROJECT)),
                "plot": str(PLOT_PATH.relative_to(PROJECT)),
            },
            sort_keys=True,
        )
    )
    return 0 if accepted else 2


if __name__ == "__main__":
    raise SystemExit(main())

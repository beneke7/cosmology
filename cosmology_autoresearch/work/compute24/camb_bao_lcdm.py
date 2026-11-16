#!/usr/bin/env python3
"""Pinned CAMB 2.0.4 flat-LCDM baseline for the released DESI DR2 BAO vector.

This is a single-probe, finite-box profile screen, not the target paper's full
Planck+BAO+SN likelihood or a reproduction of its interacting model. CAMB's
standard LambdaCDM background supplies the physical sound horizon. The three
search coordinates and bounds are copied from the target paper's Table I,
but only DESI's released 13-row BAO Gaussian likelihood is evaluated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import sys
import time

import camb
import numpy as np
import scipy
from camb.baseconfig import CAMBError, CAMBFortranError
from scipy.integrate import quad
from scipy.optimize import differential_evolution


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MEAN_PATH = ROOT / "context/data/desi_dr2_mean.txt"
COV_PATH = ROOT / "context/data/desi_dr2_cov.txt"
PAPER_PATH = ROOT / "context/papers/interacting_de_desi_dr2_2026.pdf"
WHEEL_PATH = ROOT / "context/code/camb-2.0.4-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl"
OUTPUT_PATH = HERE / "camb_bao_lcdm_result.json"
PROFILE_CHECKPOINT_PATH = HERE / "profile_checkpoint.json"
PLOT_PATH = HERE / "camb_bao_lcdm_baseline.png"

# Physical-density / acoustic-scale boxes from target paper Table I. These are
# treated as a declared search box, not as posterior priors for this BAO-only fit.
BOUNDS = ((0.5, 10.0), (0.005, 0.1), (0.001, 0.99))
PARAMETER_NAMES = ("100theta_MC", "omega_b_h2", "omega_c_h2")
FIXED = {
    "mnu_eV": 0.06,
    "N_eff": 3.044,
    "T_CMB_K": 2.7255,
    "num_massive_neutrinos": 1,
    "neutrino_hierarchy": "degenerate (CAMB default for one massive species)",
    "YHe": "BBN consistency (CAMB default)",
    "tau_reio": 0.054,
    "As": 2.1e-9,
    "ns": 0.965,
    "Omega_k": 0.0,
    "w": -1.0,
    "wa": 0.0,
}
C_KM_S = 299792.458
BAD_SCORE = 1.0e100


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_bao() -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    rows: list[tuple[float, float, str]] = []
    for line in MEAN_PATH.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        z_s, value_s, kind = line.split()
        rows.append((float(z_s), float(value_s), kind))
    z = np.asarray([row[0] for row in rows], dtype=float)
    observed = np.asarray([row[1] for row in rows], dtype=float)
    kinds = [row[2] for row in rows]
    covariance = np.loadtxt(COV_PATH, dtype=float)
    if covariance.shape != (len(rows), len(rows)):
        raise ValueError(f"mean/covariance shape mismatch: {len(rows)} vs {covariance.shape}")
    if not np.all(np.isfinite(covariance)) or not np.allclose(covariance, covariance.T, rtol=0, atol=1e-12):
        raise ValueError("BAO covariance is non-finite or asymmetric")
    np.linalg.cholesky(covariance)
    return z, observed, covariance, kinds


Z_DATA, Y_DATA, COV, OBS_KINDS = read_bao()
CHOLESKY = np.linalg.cholesky(COV)


def make_params(theta100: float, omega_b: float, omega_c: float) -> camb.CAMBparams:
    pars = camb.CAMBparams()
    # Set Lambda before CAMB solves H0 from theta_MC; no background edits follow.
    pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
    pars.set_cosmology(
        cosmomc_theta=float(theta100) / 100.0,
        ombh2=float(omega_b),
        omch2=float(omega_c),
        omk=0.0,
        mnu=FIXED["mnu_eV"],
        nnu=FIXED["N_eff"],
        num_massive_neutrinos=FIXED["num_massive_neutrinos"],
        neutrino_hierarchy="degenerate",
        TCMB=FIXED["T_CMB_K"],
        YHe=None,
        tau=FIXED["tau_reio"],
        theta_H0_range=(10.0, 100.0),
    )
    return pars


def calculate_model(
    parameters: tuple[float, float, float] | list[float] | np.ndarray,
    *,
    return_background: bool = False,
) -> dict[str, object]:
    theta100, omega_b, omega_c = (float(x) for x in parameters)
    pars = make_params(theta100, omega_b, omega_c)
    if not np.isfinite(pars.omegam) or pars.omegam >= 1.0:
        raise ValueError("CAMB point has nonpositive flat-Lambda density")
    background = camb.get_background(pars)
    derived = background.get_derived_params()
    rdrag = float(derived["rdrag"])
    if not np.isfinite(rdrag) or rdrag <= 0:
        raise ValueError("CAMB returned nonpositive/nonfinite drag sound horizon")

    h_z = np.asarray(background.hubble_parameter(Z_DATA), dtype=float)
    dm_z = np.asarray(background.comoving_radial_distance(Z_DATA), dtype=float)
    dh_z = C_KM_S / h_z
    predictions = np.empty(len(Z_DATA), dtype=float)
    for i, kind in enumerate(OBS_KINDS):
        if kind == "DM_over_rs":
            predictions[i] = dm_z[i] / rdrag
        elif kind == "DH_over_rs":
            predictions[i] = dh_z[i] / rdrag
        elif kind == "DV_over_rs":
            predictions[i] = (Z_DATA[i] * dm_z[i] ** 2 * dh_z[i]) ** (1.0 / 3.0) / rdrag
        else:
            raise ValueError(f"unrecognized BAO observable {kind!r}")

    residual = predictions - Y_DATA
    whitened = np.linalg.solve(CHOLESKY, residual)
    score_cholesky = float(whitened @ whitened)
    score_dense = float(residual @ np.linalg.solve(COV, residual))
    if abs(score_cholesky - score_dense) > 2e-9 * max(1.0, abs(score_dense)):
        raise ArithmeticError("Cholesky and dense covariance scores disagree")

    output: dict[str, object] = {
        "parameters": {
            "100theta_MC": theta100,
            "omega_b_h2": omega_b,
            "omega_c_h2": omega_c,
        },
        "derived": {
            "H0_km_s_Mpc": float(pars.H0),
            "Omega_m0_including_massive_nu": float(pars.omegam),
            "Omega_b0": float(pars.omegab),
            "Omega_c0": float(pars.omegac),
            "Omega_nu0": float(pars.omeganu),
            "r_drag_Mpc": rdrag,
            "100theta_star_CAMB": float(derived["thetastar"]),
            "z_drag": float(derived["zdrag"]),
        },
        "bao": {
            "chi2_full_covariance": score_cholesky,
            "chi2_dense_solve_check": score_dense,
            "predictions_in_input_order": predictions.tolist(),
            "residuals_in_input_order": residual.tolist(),
            "whitened_residuals_in_cholesky_order": whitened.tolist(),
            "H_km_s_Mpc_at_rows": h_z.tolist(),
            "DM_comoving_Mpc_at_rows": dm_z.tolist(),
        },
    }
    if return_background:
        output["_objects"] = {"params": pars, "background": background}
    return output


def bao_chi2(parameters: tuple[float, float, float] | list[float] | np.ndarray) -> float:
    try:
        result = calculate_model(parameters)
        return float(result["bao"]["chi2_full_covariance"])
    except (ArithmeticError, CAMBError, CAMBFortranError, FloatingPointError, RuntimeError, ValueError, OverflowError):
        return BAD_SCORE


def one_dimensional_boundary_stress(best: np.ndarray, points: int = 41) -> dict[str, object]:
    rows: dict[str, object] = {}
    for dimension, name in enumerate(PARAMETER_NAMES):
        grid = np.linspace(*BOUNDS[dimension], points)
        scores: list[float] = []
        for value in grid:
            trial = best.copy()
            trial[dimension] = value
            scores.append(bao_chi2(trial))
        rows[name] = {
            "grid": grid.tolist(),
            "chi2_at_fixed_other_parameters": scores,
            "minimum_grid_chi2": float(np.min(scores)),
            "minimum_grid_value": float(grid[int(np.argmin(scores))]),
            "edge_minima": [float(scores[0]), float(scores[-1])],
        }
    return rows


def independent_distance_checks(best: np.ndarray) -> dict[str, object]:
    result = calculate_model(best, return_background=True)
    objects = result.pop("_objects")
    background = objects["background"]
    quadrature_errors: list[float] = []
    for z, dm, h_z in zip(
        Z_DATA,
        result["bao"]["DM_comoving_Mpc_at_rows"],
        result["bao"]["H_km_s_Mpc_at_rows"],
    ):
        integral, _ = quad(
            lambda z_eval: C_KM_S / float(background.hubble_parameter(z_eval)),
            0.0,
            float(z),
            epsabs=1e-9,
            epsrel=1e-11,
            limit=250,
        )
        quadrature_errors.append(abs(float(dm) - integral))
    # Return a clean object-free summary: CAMB proxy objects are not serializable.
    result.pop("_objects", None)
    return {
        "distance_integral_max_abs_Mpc": float(max(quadrature_errors)),
        "distance_integral_all_rows_pass_1e-6_Mpc": bool(max(quadrature_errors) < 1e-6),
        "DH_definition": "c/H(z), c=299792.458 km/s",
        "full_covariance_dense_vs_cholesky_max_abs": abs(
            float(result["bao"]["chi2_full_covariance"])
            - float(result["bao"]["chi2_dense_solve_check"])
        ),
        "background_identity_note": "The radial-distance check integrates c/H(z) independently with SciPy quad against CAMB's comoving-distance output.",
    }


def run_cmb_accuracy_check(theta100: float, omega_b: float, omega_c: float) -> dict[str, object]:
    spectra: dict[int, np.ndarray] = {}
    timings: dict[int, float] = {}
    for boost in (1, 2):
        pars = make_params(theta100, omega_b, omega_c)
        pars.InitPower.set_params(As=FIXED["As"], ns=FIXED["ns"], r=0.0)
        pars.set_accuracy(AccuracyBoost=float(boost), lAccuracyBoost=float(boost), lSampleBoost=float(boost))
        pars.set_for_lmax(2500, lens_potential_accuracy=boost)
        start = time.perf_counter()
        results = camb.get_results(pars)
        powers = results.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=False)
        timings[boost] = time.perf_counter() - start
        spectra[boost] = np.asarray(powers["total"][:, 0], dtype=float)

    low = spectra[1]
    high = spectra[2]
    common_length = min(len(low), len(high), 2501)
    low_common = low[:common_length]
    high_common = high[:common_length]
    ell = np.arange(common_length)
    sel = (ell >= 30) & (ell <= 2500) & np.isfinite(low_common) & np.isfinite(high_common)
    scale = np.maximum(np.abs(high_common[sel]), 1.0)
    delta = high_common[sel] - low_common[sel]
    max_rel = float(np.max(np.abs(delta) / scale))
    rms_rel = float(np.sqrt(np.mean((delta / scale) ** 2)))
    sample_ells = [2, 30, 100, 500, 1000, 1500, 2000, 2500]
    return {
        "definition": "CAMB total lensed TT D_l = l(l+1)C_l/(2pi) in microK^2; As=2.1e-9, ns=0.965, tau=0.054, scalar adiabatic default ICs.",
        "lmax": 2500,
        "accuracy_boosts": [1, 2],
        "runtime_seconds_by_boost": {str(k): v for k, v in timings.items()},
        "max_abs_relative_delta_l30_2500_normalized_by_max_abs_reference_1uK2": max_rel,
        "rms_relative_delta_same_normalization": rms_rel,
        "convergence_gate": {
            "max_relative_below_0_02": bool(max_rel < 0.02),
            "rms_relative_below_0_005": bool(rms_rel < 0.005),
        },
        "sample_Dl_TT_microK2": {
            str(l): {"boost1": float(low_common[l]), "boost2": float(high_common[l])}
            for l in sample_ells
            if l < len(high)
        },
        "common_multipole_count": int(common_length),
        "plot_ell": ell[2:common_length].tolist(),
        "plot_Dl_boost1": low_common[2:common_length].tolist(),
        "plot_Dl_boost2": high_common[2:common_length].tolist(),
    }


def make_plot(fiducial: dict[str, object], best: dict[str, object], cmb: dict[str, object]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sigma = np.sqrt(np.diag(COV))
    x = np.arange(len(Z_DATA))
    labels = [f"{z:g}\n{kind.replace('_over_rs', '')}" for z, kind in zip(Z_DATA, OBS_KINDS)]
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True, gridspec_kw={"height_ratios": [2, 1]})
    axes[0].errorbar(x, Y_DATA, yerr=sigma, fmt="ko", ms=4, capsize=2, label="DESI DR2 mean (marginal σ)")
    axes[0].plot(x, fiducial["bao"]["predictions_in_input_order"], "s--", color="#377eb8", label="CAMB ΛCDM fiducial")
    axes[0].plot(x, best["bao"]["predictions_in_input_order"], "o-", color="#e41a1c", label="CAMB ΛCDM BAO-only profile")
    axes[0].set_ylabel("BAO distance / r_drag")
    axes[0].set_xticks(x, labels)
    axes[0].legend(frameon=False, ncol=3)
    axes[0].grid(alpha=0.25)

    whitened = np.asarray(best["bao"]["whitened_residuals_in_cholesky_order"], dtype=float)
    axes[1].axhline(0.0, color="k", lw=0.8)
    axes[1].bar(x, whitened, color="#984ea3")
    axes[1].set_ylabel("Cholesky-whitened residual")
    axes[1].set_xticks(x, labels)
    axes[1].grid(axis="y", alpha=0.25)
    fig.suptitle("CAMB 2.0.4 flat ΛCDM × DESI DR2 BAO (single-probe screen)")
    fig.savefig(PLOT_PATH, dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    ell = np.asarray(cmb["plot_ell"], dtype=int)
    ax.plot(ell, cmb["plot_Dl_boost1"], lw=0.9, label="AccuracyBoost=1")
    ax.plot(ell, cmb["plot_Dl_boost2"], lw=0.9, ls="--", label="AccuracyBoost=2")
    ax.set_yscale("log")
    ax.set_xlim(2, 2500)
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"Lensed $D_\ell^{TT}$ [$\mu$K$^2$]")
    ax.set_title("Standard adiabatic ΛCDM CAMB baseline; convergence comparison, not a data fit")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.savefig(HERE / "camb_lcdm_tt_convergence.png", dpi=170)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--maxiter", type=int, default=60)
    parser.add_argument("--popsize", type=int, default=12)
    parser.add_argument("--seeds", type=int, nargs="+", default=[20260925, 20260926])
    args = parser.parse_args()
    if args.workers < 1 or args.maxiter < 1 or args.popsize < 4:
        raise SystemExit("workers/maxiter must be positive and popsize >= 4")
    available_cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    if args.workers > available_cpus:
        raise SystemExit(f"requested {args.workers} workers but only {available_cpus} CPUs are available to this process")

    start_all = time.perf_counter()
    fiducial_parameters = np.asarray([1.04, 0.0224, 0.1200], dtype=float)
    fiducial = calculate_model(fiducial_parameters)

    fit_runs: list[dict[str, object]] = []
    for seed in args.seeds:
        start = time.perf_counter()
        fit = differential_evolution(
            bao_chi2,
            bounds=BOUNDS,
            strategy="best1bin",
            maxiter=args.maxiter,
            popsize=args.popsize,
            tol=1e-7,
            atol=1e-8,
            seed=int(seed),
            polish=True,
            init="latinhypercube",
            updating="deferred",
            workers=args.workers,
            disp=False,
        )
        fit_runs.append(
            {
                "seed": int(seed),
                "success": bool(fit.success),
                "message": str(fit.message),
                "nfev": int(fit.nfev),
                "nit": int(fit.nit),
                "runtime_seconds": time.perf_counter() - start,
                "chi2": float(fit.fun),
                "x": np.asarray(fit.x, dtype=float).tolist(),
                "population_energies_min": float(np.min(fit.population_energies)),
                "population_energies_max": float(np.max(fit.population_energies)),
            }
        )
    best_fit_run = min(fit_runs, key=lambda row: float(row["chi2"]))
    best_parameters = np.asarray(best_fit_run["x"], dtype=float)
    best = calculate_model(best_parameters)
    boundary_stress = one_dimensional_boundary_stress(best_parameters)
    distance_checks = independent_distance_checks(best_parameters)
    profile_checkpoint = {
        "status": "profile_and_background_checks_complete_CMB_spectrum_stage_pending",
        "profile_runs": fit_runs,
        "selected_profile": best,
        "boundary_stress": boundary_stress,
        "distance_checks": distance_checks,
        "worker_count": args.workers,
        "source_script_sha256": sha256(Path(__file__)),
        "generated_unix_utc": time.time(),
    }
    PROFILE_CHECKPOINT_PATH.write_text(json.dumps(profile_checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    cmb = run_cmb_accuracy_check(*fiducial_parameters)

    make_plot(fiducial, best, cmb)
    report = {
        "schema_version": 1,
        "status": "completed_exploratory_CAMB_LCDM_BAO_profile",
        "scope": {
            "model": "flat LambdaCDM, Q=0, w=-1, standard GR/CAMB 2.0.4 background; standard adiabatic CMB spectra only as numerical baseline",
            "observational_likelihood": "DESI DR2 13-row compressed BAO Gaussian mean and full 13x13 covariance only",
            "not_included": ["Planck/ACT likelihood", "supernovae", "dark-sector interaction", "posterior/evidence", "cross-survey covariance", "target-paper modified solver"],
            "prior_interpretation": "Table-I parameter bounds copied as a finite search box; they are not a posterior prior or reproduction of the target paper's full likelihood.",
        },
        "source_and_model": {
            "CAMB_version": camb.__version__,
            "CAMB_repository_commit": "a6de8cc59124c8bbe3924f1c546f96cef73dcabb",
            "CAMB_wheel_sha256": sha256(WHEEL_PATH),
            "CAMB_wheel_license": "LGPLv3-arxiv expression, from wheel METADATA / local source LICENCE.txt",
            "CAMB_source_archive_sha256": sha256(ROOT / "context/code/CAMB_2.0.4_a6de8cc59124c8bbe3924f1c546f96cef73dcabb.tar.gz"),
            "fixed_physics": FIXED,
            "free_search_parameters": list(PARAMETER_NAMES),
            "bounds": {name: list(bound) for name, bound in zip(PARAMETER_NAMES, BOUNDS)},
            "theta_H0_solver_range_km_s_Mpc": [10.0, 100.0],
            "fiducial_parameters": dict(zip(PARAMETER_NAMES, fiducial_parameters.tolist())),
        },
        "inputs": {
            "bao_mean_path": str(MEAN_PATH.relative_to(ROOT)),
            "bao_mean_sha256": sha256(MEAN_PATH),
            "bao_covariance_path": str(COV_PATH.relative_to(ROOT)),
            "bao_covariance_sha256": sha256(COV_PATH),
            "target_paper_path": str(PAPER_PATH.relative_to(ROOT)),
            "target_paper_sha256": sha256(PAPER_PATH),
            "row_count": int(len(Y_DATA)),
            "row_order": [
                {"z": float(z), "observable": kind, "mean": float(y)}
                for z, kind, y in zip(Z_DATA, OBS_KINDS, Y_DATA)
            ],
        },
        "fiducial_DESI_score": fiducial,
        "profile_runs": fit_runs,
        "selected_profile": best,
        "boundary_stress": boundary_stress,
        "checks": {
        "distance_and_covariance": distance_checks,
        "CMB_spectrum_accuracy": cmb,
        "all_required_checks_pass": bool(
            distance_checks["distance_integral_all_rows_pass_1e-6_Mpc"]
            and cmb["convergence_gate"]["max_relative_below_0_02"]
            and cmb["convergence_gate"]["rms_relative_below_0_005"]
            and all(np.isfinite(float(row["chi2"])) for row in fit_runs)
            and any(bool(row["success"]) and float(row["chi2"]) < BAD_SCORE for row in fit_runs)
            and float(best["bao"]["chi2_full_covariance"]) < BAD_SCORE
        ),
        },
        "runtime": {
            "total_seconds": time.perf_counter() - start_all,
            "parallel_workers": args.workers,
            "available_cpus_in_affinity": available_cpus,
            "differential_evolution": {
                "maxiter": args.maxiter,
                "popsize": args.popsize,
                "strategy": "best1bin",
                "updating": "deferred",
                "polish": True,
            },
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "interpretation": {
            "claim_limit": "The result is a finite-box standard-physics BAO baseline. It cannot establish the target interaction, its posterior, or the CMB-likelihood result.",
            "comparison_limit": "Any chi-square difference from the prior two-parameter free-alpha screen is descriptive: this CAMB parameterization predicts r_drag and has three finite-box parameters, with massive-neutrino and BBN choices stated above.",
            "new_research_question": "Does enforcing a physical sound horizon under the target paper's very broad physical-density/thetaMC box leave the BAO-only LambdaCDM fit well identified, or do BAO data alone drive parameters to the search edges?",
        },
    }
    report["script_sha256"] = sha256(Path(__file__))
    OUTPUT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "all_required_checks_pass": report["checks"]["all_required_checks_pass"],
                "fiducial_chi2": fiducial["bao"]["chi2_full_covariance"],
                "selected_chi2": best["bao"]["chi2_full_covariance"],
                "selected_parameters": best["parameters"],
                "selected_derived": best["derived"],
                "profile_seeds": [
                    {key: run[key] for key in ("seed", "success", "chi2", "x", "nfev", "nit", "runtime_seconds")}
                    for run in fit_runs
                ],
                "distance_checks": distance_checks,
                "CMB_convergence": cmb["convergence_gate"],
                "CMB_delta": {
                    "max_rel": cmb["max_abs_relative_delta_l30_2500_normalized_by_max_abs_reference_1uK2"],
                    "rms_rel": cmb["rms_relative_delta_same_normalization"],
                },
                "runtime_seconds": report["runtime"]["total_seconds"],
                "result": str(OUTPUT_PATH.relative_to(ROOT)),
                "plot": str(PLOT_PATH.relative_to(ROOT)),
            },
            sort_keys=True,
        )
    )
    return 0 if report["checks"]["all_required_checks_pass"] else 2


if __name__ == "__main__":
    if "spawn" not in mp.get_all_start_methods():
        raise RuntimeError("CAMB profile requires a multiprocessing spawn context")
    mp.set_start_method("spawn", force=True)
    raise SystemExit(main())

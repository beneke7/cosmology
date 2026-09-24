#!/usr/bin/env python3
"""Boundary, conditional-prediction, Jacobian, and null-mock checks for DESI DR2 BAO.

This is a profile-likelihood screening experiment. It does not compute posterior
constraints or Bayesian evidence. The null mocks repeat the stated finite model
search under the fitted flat LCDM screening model.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import time

# Each process uses one BLAS thread; the process pool owns parallelism.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.stats import beta

from background_bao import BAOData, FIT_BOUNDS, fit_model, predict_bao


CPL_W0_UPPER_SEARCH = [-0.3, -0.1, 0.0, 0.1]
MODEL_STARTS = 16


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scalar_quad_prediction(data: BAOData, pars: tuple[float, float, float, float]) -> np.ndarray:
    """Independent scalar distance implementation using adaptive quadrature."""
    alpha, omega_m, w0, wa = pars
    omega_de = 1.0 - omega_m

    def e_of_z(z: float) -> float:
        zp = 1.0 + z
        de = math.exp(3.0 * ((1.0 + w0 + wa) * math.log(zp) - wa * z / zp))
        return math.sqrt(omega_m * zp**3 + omega_de * de)

    values = []
    for z, obs in zip(data.z, data.observable):
        zf = float(z)
        dm = alpha * quad(lambda x: 1.0 / e_of_z(x), 0.0, zf,
                          epsabs=2e-12, epsrel=2e-12)[0]
        dh = alpha / e_of_z(zf)
        if obs == "DM_over_rs":
            values.append(dm)
        elif obs == "DH_over_rs":
            values.append(dh)
        elif obs == "DV_over_rs":
            values.append((zf * dm * dm * dh) ** (1.0 / 3.0))
        else:
            raise ValueError(f"unsupported observable {obs!r}")
    return np.asarray(values, dtype=np.float64)


def subset_data(data: BAOData, keep: np.ndarray) -> BAOData:
    idx = np.asarray(keep, dtype=int)
    return BAOData(data.z[idx], data.value[idx], data.observable[idx],
                   data.covariance[np.ix_(idx, idx)])


def fit_cpl_at_upper(data: BAOData, upper: float, starts: int = MODEL_STARTS) -> dict:
    original = FIT_BOUNDS["cpl"]
    FIT_BOUNDS["cpl"] = [original[0], original[1], (-2.0, float(upper)), original[3]]
    try:
        result = fit_model(data, "cpl", starts=starts)
        result["parameter_bounds"] = FIT_BOUNDS["cpl"]
        result["starts"] = starts
        result["optimizer_seed"] = 20260923
        return result
    finally:
        FIT_BOUNDS["cpl"] = original


def profile_scan(data: BAOData) -> list[dict]:
    return [fit_cpl_at_upper(data, upper) | {"w0_upper": upper}
            for upper in CPL_W0_UPPER_SEARCH]


def response_analysis(data: BAOData, lcdm_fit: dict) -> dict:
    pars = np.asarray([lcdm_fit["alpha"], lcdm_fit["Omega_m"], -1.0, 0.0],
                      dtype=np.float64)
    vectorized = predict_bao(data.z, data.observable, *pars)
    independent = scalar_quad_prediction(data, tuple(map(float, pars)))
    np.testing.assert_allclose(vectorized, independent, rtol=2e-11, atol=2e-11)

    # Scale columns to finite, interpretable changes: 1% in alpha, 0.05 in
    # Omega_m, 0.3 in w0, and 1 in wa. SVD conclusions are therefore explicit
    # about parameter scaling rather than silently unit-dependent.
    scales = np.asarray([0.01 * pars[0], 0.05, 0.3, 1.0], dtype=np.float64)
    steps = [1e-3, 5e-4, 2.5e-4]
    jacobians = []
    for relative_step in steps:
        columns = []
        for j in range(4):
            h = relative_step * scales[j]
            samples = []
            for multiplier in (-2.0, -1.0, 1.0, 2.0):
                shifted = pars.copy()
                shifted[j] += multiplier * h
                samples.append(scalar_quad_prediction(data, tuple(map(float, shifted))))
            derivative = (samples[0] - 8.0 * samples[1] + 8.0 * samples[2] - samples[3]) / (12.0 * h)
            columns.append(derivative * scales[j])
        jacobians.append(np.column_stack(columns))

    jacobian_delta = float(np.linalg.norm(jacobians[-1] - jacobians[-2])
                           / max(np.linalg.norm(jacobians[-1]), np.finfo(float).tiny))
    chol = np.linalg.cholesky(data.covariance)
    whitened = solve_triangular(chol, jacobians[-1], lower=True, check_finite=False)
    nuisance = whitened[:, :2]
    q, _ = np.linalg.qr(nuisance, mode="reduced")
    ext = whitened[:, 2:]
    projected = ext - q @ (q.T @ ext)
    ext_singular = np.linalg.svd(projected, compute_uv=False)
    full_singular = np.linalg.svd(whitened, compute_uv=False)
    thresholds = np.finfo(float).eps * max(whitened.shape) * full_singular[0]
    ext_threshold = np.finfo(float).eps * max(projected.shape) * max(ext_singular[0], 1.0)
    ext_rank = int(np.count_nonzero(ext_singular > ext_threshold))
    projected_norms = np.linalg.norm(projected, axis=0)
    cos_ext = float(np.dot(projected[:, 0], projected[:, 1])
                    / (projected_norms[0] * projected_norms[1]))
    raw_norms = np.linalg.norm(ext, axis=0)
    nuisance_projection_fraction = [
        float(max(0.0, 1.0 - (projected_norms[j] / raw_norms[j]) ** 2))
        for j in range(2)
    ]
    return {
        "reference": "flat LCDM maximum-likelihood screen",
        "parameter_order": ["alpha", "Omega_m", "w0", "wa"],
        "parameter_scales": {"alpha": float(scales[0]), "Omega_m": float(scales[1]),
                             "w0": float(scales[2]), "wa": float(scales[3])},
        "units": {"BAO observables": "dimensionless distance / sound-horizon ratios",
                  "alpha": "dimensionless c/(H0*rd)", "Omega_m": "dimensionless",
                  "w0": "dimensionless", "wa": "dimensionless"},
        "finite_difference_relative_steps": steps,
        "relative_jacobian_change_last_step_halving": jacobian_delta,
        "independent_scalar_quad_max_prediction_difference": float(np.max(np.abs(vectorized - independent))),
        "full_whitened_jacobian_singular_values": full_singular.tolist(),
        "full_jacobian_rank_machine_epsilon_threshold": int(np.count_nonzero(full_singular > thresholds)),
        "extension_after_projection_over_alpha_omega_m": {
            "column_order": ["w0_step_0.3", "wa_step_1"],
            "singular_values": ext_singular.tolist(),
            "rank_machine_epsilon_threshold": ext_rank,
            "fraction_of_each_raw_whitened_extension_norm_projected_onto_nuisance_span":
                nuisance_projection_fraction,
            "cosine_between_residualized_w0_and_wa_directions": cos_ext,
        },
        "method": "five-point finite differences of an independent adaptive-quadrature implementation; full covariance whitened by Cholesky; extension directions projected orthogonally to alpha and Omega_m nuisance directions",
    }


def conditional_block_cv(data: BAOData) -> list[dict]:
    """Leave each redshift block out, keeping its full Gaussian conditioning."""
    results = []
    for z_hold in sorted(set(map(float, data.z))):
        test_idx = np.flatnonzero(np.isclose(data.z, z_hold, rtol=0.0, atol=1e-12))
        train_idx = np.asarray([i for i in range(len(data.z)) if i not in set(test_idx)], dtype=int)
        train = subset_data(data, train_idx)
        c_aa = data.covariance[np.ix_(train_idx, train_idx)]
        c_ta = data.covariance[np.ix_(test_idx, train_idx)]
        c_at = c_ta.T
        c_tt = data.covariance[np.ix_(test_idx, test_idx)]
        c_cond = c_tt - c_ta @ np.linalg.solve(c_aa, c_at)
        c_cond = (c_cond + c_cond.T) / 2.0
        chol_cond = np.linalg.cholesky(c_cond)
        fold = {"held_out_redshift": z_hold,
                "held_out_rows": [{"index": int(i), "observable": str(data.observable[i]),
                                   "value": float(data.value[i])} for i in test_idx],
                "n_train": int(len(train_idx)), "n_test": int(len(test_idx)),
                "conditional_covariance_min_eigenvalue": float(np.linalg.eigvalsh(c_cond)[0]),
                "models": {}}
        for model in ("lcdm", "wcdm", "cpl"):
            fit = fit_model(train, model, starts=MODEL_STARTS)
            params = [fit[k] for k in ("alpha", "Omega_m", "w0", "wa")]
            mu_all = predict_bao(data.z, data.observable, *params)
            delta_train = data.value[train_idx] - mu_all[train_idx]
            conditional_mean = (mu_all[test_idx]
                                + c_ta @ np.linalg.solve(c_aa, delta_train))
            residual = data.value[test_idx] - conditional_mean
            whitened_residual = solve_triangular(chol_cond, residual, lower=True,
                                                 check_finite=False)
            chi2_cond = float(whitened_residual @ whitened_residual)
            logdet = float(2.0 * np.sum(np.log(np.diag(chol_cond))))
            log_density = float(-0.5 * (len(test_idx) * math.log(2.0 * math.pi)
                                        + logdet + chi2_cond))
            fold["models"][model] = {
                "parameters": {k: fit[k] for k in ("alpha", "Omega_m", "w0", "wa")},
                "training_chi2": fit["chi2"], "training_success": fit["success"],
                "training_boundary_hits": fit["parameter_bound_hits"],
                "conditional_mean": conditional_mean.tolist(),
                "conditional_residual": residual.tolist(),
                "conditional_chi2": chi2_cond,
                "conditional_log_predictive_density": log_density,
            }
        fold["interpretation_note"] = "Plug-in maximum-likelihood predictions; this is one overlapping conditional fold, not a posterior predictive distribution or an independent-probe product."
        results.append(fold)
    return results


def candidate_search(data: BAOData, starts: int) -> dict:
    """One fully specified model/bound search at a fixed multi-start budget."""
    fits = {"lcdm": fit_model(data, "lcdm", starts=starts),
            "wcdm": fit_model(data, "wcdm", starts=starts)}
    for upper in CPL_W0_UPPER_SEARCH:
        fits[f"cpl_w0_upper_{upper:g}"] = fit_cpl_at_upper(data, upper, starts=starts)
    selected = min((key for key in fits if key != "lcdm"), key=lambda key: fits[key]["chi2"])
    return {
        "starts_per_fit": starts,
        "selected_candidate": selected,
        "selected_candidate_chi2": float(fits[selected]["chi2"]),
        "search_statistic_delta_chi2": float(fits["lcdm"]["chi2"] - fits[selected]["chi2"]),
        "candidate_fits": fits,
    }


def fit_mock_search(args: tuple) -> dict:
    """Fit one null realization with both the primary and expanded start budgets."""
    index, values, z, observable, covariance = args
    mock = BAOData(z, values, observable, covariance)
    primary = candidate_search(mock, MODEL_STARTS)
    expanded = candidate_search(mock, 64)
    return {"index": int(index), "primary_16_starts": primary,
            "expanded_64_starts": expanded,
            "search_statistic_delta_chi2": primary["search_statistic_delta_chi2"],
            "expanded_search_statistic_delta_chi2": expanded["search_statistic_delta_chi2"],
            "primary_failure_count": sum(not fit["success"] for fit in primary["candidate_fits"].values()),
            "expanded_failure_count": sum(not fit["success"] for fit in expanded["candidate_fits"].values()),
            "boundary_solution_count": sum(len(fit["parameter_bound_hits"])
                                            for fit in primary["candidate_fits"].values())}


def null_mock_calibration(data: BAOData, lcdm_fit: dict, n_mocks: int, seed: int,
                          workers: int) -> dict:
    pars = [lcdm_fit["alpha"], lcdm_fit["Omega_m"], -1.0, 0.0]
    mean = predict_bao(data.z, data.observable, *pars)
    chol = np.linalg.cholesky(data.covariance)
    rng = np.random.default_rng(seed)
    normals = rng.standard_normal((n_mocks, len(data.z)))
    mocks = mean[None, :] + normals @ chol.T
    common = [(i, mocks[i], data.z, data.observable, data.covariance)
              for i in range(n_mocks)]
    results = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, row in enumerate(pool.map(fit_mock_search, common, chunksize=1), start=1):
            results.append(row)
            if i % 10 == 0 or i == n_mocks:
                print(f"completed primary and expanded-start fits for {i}/{n_mocks} null mocks", flush=True)

    observed_16 = {"lcdm": lcdm_fit, "wcdm": fit_model(data, "wcdm", starts=MODEL_STARTS)}
    for fit in profile_scan(data):
        observed_16[f"cpl_w0_upper_{fit['w0_upper']:g}"] = fit
    selected_16 = min((key for key in observed_16 if key != "lcdm"),
                      key=lambda key: observed_16[key]["chi2"])
    obs_primary = {"starts_per_fit": MODEL_STARTS, "selected_candidate": selected_16,
                   "selected_candidate_chi2": observed_16[selected_16]["chi2"],
                   "search_statistic_delta_chi2": float(lcdm_fit["chi2"] - observed_16[selected_16]["chi2"]),
                   "candidate_fits": observed_16}
    obs_expanded = candidate_search(data, 64)
    exceedances = {}
    intervals = {}
    summaries = {}
    for label, statistic_name, observed_stat in (
            ("16_starts", "search_statistic_delta_chi2", obs_primary["search_statistic_delta_chi2"]),
            ("64_starts", "expanded_search_statistic_delta_chi2", obs_expanded["search_statistic_delta_chi2"])):
        statistics = np.asarray([row[statistic_name] for row in results])
        k = int(np.count_nonzero(statistics >= observed_stat))
        n = int(len(statistics))
        lower = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
        upper = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
        exceedances[label] = k
        intervals[label] = [lower, upper]
        summaries[label] = {"observed_delta_chi2": float(observed_stat),
                            "mock_tail_fraction_at_least_observed": k / n,
                            "clopper_pearson_95pct_interval": [lower, upper],
                            "null_search_statistic_summary": {
                                "median": float(np.median(statistics)),
                                "q90": float(np.quantile(statistics, 0.9)),
                                "q95": float(np.quantile(statistics, 0.95)),
                                "maximum": float(np.max(statistics)),
                            }}
    primary_fits = sum(len(row["primary_16_starts"]["candidate_fits"]) for row in results)
    expanded_fits = sum(len(row["expanded_64_starts"]["candidate_fits"]) for row in results)
    return {
        "null": "data generated from the observed best-fit flat LCDM background screen with the released full 13x13 Gaussian covariance",
        "fit_search_repeated_per_realization": {
            "models": ["lcdm", "constant_w", "CPL"],
            "cpl_w0_upper_bounds_searched": CPL_W0_UPPER_SEARCH,
            "all other bounds": FIT_BOUNDS,
            "primary_deterministic_multistarts_per_fit": MODEL_STARTS,
            "expanded_convergence_check_starts_per_fit": 64,
            "optimizer_seed_each_fit": 20260923,
        },
        "n_mocks": n, "random_generator": "NumPy default_rng PCG64",
        "random_seed": seed, "workers": workers,
        "observed_searches": {"16_starts": obs_primary, "64_starts": obs_expanded},
        "search_calibration_by_start_budget": summaries,
        "optimizer_nonconverged_candidate_fits": {
            "16_starts": int(sum(row["primary_failure_count"] for row in results)),
            "64_starts": int(sum(row["expanded_failure_count"] for row in results)),
            "16_start_candidate_fits_total": primary_fits,
            "64_start_candidate_fits_total": expanded_fits,
        },
        "boundary_solution_count_total_primary_search": int(sum(row["boundary_solution_count"] for row in results)),
        "realizations": results,
        "interpretation": "Coarse parametric-bootstrap false-positive audit for this exact adaptive profile-fit search; not a discovery significance or posterior/evidence calculation.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mean", default="context/data/desi_dr2_mean.txt")
    parser.add_argument("--cov", default="context/data/desi_dr2_cov.txt")
    parser.add_argument("--output", default="experiments/bao_robustness/result.json")
    parser.add_argument("--n-mocks", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260924)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--skip-mocks", action="store_true")
    args = parser.parse_args()
    if args.n_mocks < 1 or args.workers < 1:
        parser.error("--n-mocks and --workers must be positive")

    started_utc = dt.datetime.now(dt.timezone.utc)
    started = time.perf_counter()
    data = BAOData.from_files(args.mean, args.cov)
    baseline = {model: fit_model(data, model, starts=MODEL_STARTS)
                for model in ("lcdm", "wcdm", "cpl")}
    for model in baseline.values():
        model["starts"] = MODEL_STARTS
        model["optimizer_seed"] = 20260923
        model["parameter_bounds"] = FIT_BOUNDS[model["model"]]
    cpl_profile = profile_scan(data)
    cpl_profile_scan_cv = conditional_block_cv(data)
    response = response_analysis(data, baseline["lcdm"])
    mocks = None if args.skip_mocks else null_mock_calibration(
        data, baseline["lcdm"], args.n_mocks, args.seed, args.workers)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "exploratory_robustness_screen",
        "experiment_id": "desi_dr2_bao_profile_robustness_20260924",
        "start_utc": started_utc.isoformat(),
        "scope": "13-row compressed DESI DR2 BAO only; flat late-time background with radiation omitted; free dimensionless alpha=c/(H0*rd); correlated Gaussian likelihood; profile chi-squared only",
        "source": {
            "data_files": {str(Path(args.mean)): sha256(Path(args.mean)),
                           str(Path(args.cov)): sha256(Path(args.cov))},
            "code_files": {str(Path(__file__)): sha256(Path(__file__)),
                           "scripts/background_bao.py": sha256(Path("scripts/background_bao.py"))},
        },
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": __import__("scipy").__version__,
                        "logical_cpus_visible": os.cpu_count()},
        "assumptions_and_units": ["All BAO observables are dimensionless distances divided by rd.",
                                  "alpha=c/(H0*rd) is a free dimensionless amplitude; H0 and rd are not separately identified.",
                                  "Flat late-time CPL background only; radiation, curvature, perturbations, and a sound-horizon model are omitted.",
                                  "Covariance is used in the supplied mean-file order through Cholesky/linear solves."],
        "parameter_bounds": {key: value for key, value in FIT_BOUNDS.items()},
        "starts_per_fit": MODEL_STARTS,
        "optimizer_seed": 20260923,
        "baseline_profile_fits": baseline,
        "cpl_w0_upper_bound_sensitivity": cpl_profile,
        "redshift_block_conditional_cross_validation": cpl_profile_scan_cv,
        "rank_aware_whitened_response": response,
        "adaptive_search_null_mock_calibration": mocks,
        "runtime_seconds": time.perf_counter() - started,
        "command": ".venv/bin/python scripts/run_bao_robustness.py --n-mocks 200 --seed 20260924 --workers 16 --output experiments/bao_robustness/result.json",
        "random_seed": args.seed,
        "interpretation_and_limits": ["Exploratory profile minima are not posterior constraints or Bayesian evidence.",
                                      "The redshift-block scores are plug-in conditional Gaussian predictions; folds overlap and are not summed into an independent likelihood.",
                                      "The mock calibration evaluates the observed candidate-selection/search budget and is a coarse false-positive audit, not a discovery claim."],
    }
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "status": result["status"],
                      "runtime_seconds": result["runtime_seconds"],
                      "baseline_chi2": {key: value["chi2"] for key, value in baseline.items()},
                      "cpl_profile": [{"w0_upper": row["w0_upper"], "chi2": row["chi2"],
                                       "w0": row["w0"], "hits": row["parameter_bound_hits"]}
                                      for row in cpl_profile],
                      "mock_count": 0 if mocks is None else mocks["n_mocks"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

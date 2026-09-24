#!/usr/bin/env python3
"""Review the saved profile scan with an independent scalar implementation."""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path
import time

import numpy as np

from scalar_reference import (ROOT, MEAN, COV, SEED, digest, displacement,
                              physical_radiation_density, prediction, score)


RESULT = ROOT / "experiments/radiation_sensitivity/result.json"
WORKER_CODE = ROOT / "experiments/radiation_sensitivity/radiation_scan.py"
THEORY_CODE = ROOT / "work/theory2/check_radiation.py"
THEORY_RESULT = ROOT / "work/theory2/radiation_check.json"


def main():
    started = time.perf_counter()
    raw = json.loads(RESULT.read_text())
    numeric = np.loadtxt(MEAN, usecols=(0, 1))
    labels = np.loadtxt(MEAN, usecols=(2,), dtype=str)
    zs, observed = numeric.T
    cov = np.loadtxt(COV)
    rad = physical_radiation_density()
    assert raw["inputs"]["mean_sha256"] == digest(MEAN)
    assert raw["inputs"]["covariance_sha256"] == digest(COV)
    assert raw["inputs"]["seed_result_sha256"] == digest(SEED)
    assert raw["code_sha256"][str(WORKER_CODE.relative_to(ROOT))] == digest(WORKER_CODE)
    baseline_code = ROOT / "scripts/background_bao.py"
    assert digest(baseline_code) == "34bdf1d81397b39b740faa3e2f4a79748d5fa9d3bfdd8303619d631dbf8aa2e2"
    assert np.allclose(raw["physical_radiation_derivation"]["omega_r_h2"], rad["omega_r"], rtol=1e-13, atol=0)
    assert np.allclose(raw["physical_radiation_derivation"]["omega_gamma_h2"], rad["omega_gamma"], rtol=1e-13, atol=0)

    fits = []
    fits.extend(raw["no_radiation_baseline"].values())
    fits.extend(fit for by_h0 in raw["radiation_refits"].values() for fit in by_h0.values())
    records = []
    fixed_records = []
    optimizer_status_summary = {"LBFGSB_runs": 0, "LBFGSB_failed_runs": 0,
                                "DE_runs": 0, "DE_failed_runs": 0,
                                "face_checks": 0, "face_check_failed_runs": 0}
    for fit in fits:
        pars = [fit["parameters"][key] for key in ["alpha", "Omega_m", "w0", "wa"]]
        h0 = fit.get("H0_km_s_Mpc")
        omega_r = rad["omega_r"] / (h0 / 100)**2 if h0 is not None else 0.0
        assert abs(omega_r - fit["Omega_r"]) < 1e-16
        assert abs(fit["Omega_de_from_closure"] - (1 - pars[1] - omega_r)) < 1e-15
        assert 1e-6 <= pars[0] <= 1e4 and .05 <= pars[1] <= .6
        if fit["model"] in ("wcdm", "cpl"):
            assert -2 <= pars[2] <= -.3
        if fit["model"] == "cpl":
            assert -3 <= pars[3] <= 3
        independent, q_error = prediction(zs, labels, pars, omega_r)
        chi2 = score(observed, cov, independent)
        q = independent / pars[0]
        alpha_star = float(q @ np.linalg.solve(cov, observed)
                           / (q @ np.linalg.solve(cov, q)))
        alpha_star = float(np.clip(alpha_star, 1e-6, 1e4))
        pred_error = float(np.max(np.abs(independent - fit["prediction"])))
        chi_error = abs(chi2 - fit["chi2"])
        alpha_error = abs(alpha_star - pars[0])
        assert pred_error < 1e-10 and chi_error < 1e-9 and alpha_error < 1e-10
        record = {"model": fit["model"], "H0_km_s_Mpc": h0,
                  "chi2_independent": chi2, "max_abs_prediction_error": pred_error,
                  "abs_chi2_error": chi_error, "abs_alpha_profile_error": alpha_error,
                  "quad_reported_max_abs_error": q_error}
        if h0 is not None:
            base_fit = raw["no_radiation_baseline"][fit["model"]]
            base_pars = [base_fit["parameters"][key] for key in ["alpha", "Omega_m", "w0", "wa"]]
            base_pred, _ = prediction(zs, labels, base_pars, 0)
            fixed_pred, _ = prediction(zs, labels, base_pars, omega_r)
            fixed_worker = raw["fixed_baseline_radiation_shifts"][fit["model"]][str(h0)]
            fixed_delta = fixed_pred - base_pred
            fixed_error = float(np.max(np.abs(fixed_delta - fixed_worker["delta_prediction"])))
            assert fixed_error < 1e-10, "Fixed-parameter radiation check must hold alpha and shape fixed."
            refitted_displacement = displacement(independent - base_pred, cov)
            fixed_displacement = displacement(fixed_delta, cov)
            record["refitted_displacement"] = refitted_displacement
            record["delta_chi2_from_baseline"] = chi2 - score(observed, cov, base_pred)
            fixed_records.append({"model": fit["model"], "H0_km_s_Mpc": h0,
                                  "max_abs_delta_prediction_error": fixed_error,
                                  "max_abs_fractional_shift": float(np.max(np.abs(fixed_delta / base_pred))),
                                  **fixed_displacement})
        records.append(record)
        optimizers = fit["optimizers"]
        lbfgs = optimizers["bounded_LBFGSB_multistart"]
        optimizer_status_summary["LBFGSB_runs"] += len(lbfgs)
        optimizer_status_summary["LBFGSB_failed_runs"] += sum(not run["success"] for run in lbfgs)
        optimizer_status_summary["DE_runs"] += 1
        optimizer_status_summary["DE_failed_runs"] += int(not optimizers["differential_evolution"]["success"])
        for face in optimizers["boundary_face_Powell_checks"]:
            optimizer_status_summary["face_checks"] += 1
            optimizer_status_summary["face_check_failed_runs"] += int(not face["success"])
            assert face["inward_delta_chi2"] >= -1e-7

    theory = json.loads(THEORY_RESULT.read_text())
    assert theory["script_sha256"] == digest(THEORY_CODE)
    theory_checks = []
    for item in theory["models"]:
        pars = [item["reference_parameters"][key] for key in ["alpha", "Omega_m", "w0", "wa"]]
        base_pred, _ = prediction(zs, labels, pars, 0)
        base_score = score(observed, cov, base_pred)
        grid = []
        for h0 in np.linspace(50, 90, 41):
            changed, _ = prediction(zs, labels, pars, rad["omega_r"] / (h0 / 100)**2)
            shift = displacement(changed - base_pred, cov)
            grid.append({"H0_km_s_Mpc": float(h0),
                         "max_abs_shift_over_marginal_row_sigma": shift["max_abs_shift_over_marginal_sigma"],
                         "full_covariance_shift_mahalanobis": shift["full_covariance_displacement_norm"],
                         "abs_delta_chi2_at_same_fit_parameters": abs(score(observed, cov, changed) - base_score)})
        summary = item["H0_screen_41_point_grid_summary"]
        assert summary["point_count"] == len(grid) == 41
        for summary_key, value_key in [
            ("max_marginal_sigma_shift", "max_abs_shift_over_marginal_row_sigma"),
            ("max_full_covariance_displacement", "full_covariance_shift_mahalanobis"),
            ("max_abs_delta_chi2_at_same_fit_parameters", "abs_delta_chi2_at_same_fit_parameters"),
        ]:
            independent_max = max(grid, key=lambda row: row[value_key])
            saved = summary[summary_key]
            assert independent_max["H0_km_s_Mpc"] == saved["H0_km_s_Mpc"]
            error = abs(independent_max[value_key] - saved[value_key])
            assert error < 1e-10
            theory_checks.append({"model": item["model"], "metric": value_key,
                                  "independent_max": independent_max[value_key],
                                  "H0_at_max": independent_max["H0_km_s_Mpc"],
                                  "abs_difference_from_saved_theory_result": error})

    final = {
        "status": "independently_checked_saved_profile_scan",
        "reviewer": "/root/orchestrator2_radiation",
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "method": "Independent scalar adaptive distances and raw full-covariance dense solves; no worker/campaign model function imported.",
        "command": "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ./.venv/bin/python scripts/run_bounded.py --seconds 120 -- ./.venv/bin/python work/orchestrator2_radiation/review_saved_fits.py",
        "hashes_verified": True,
        "sha256": {str(path.relative_to(ROOT)): digest(path) for path in
                   [MEAN, COV, SEED, RESULT, WORKER_CODE, THEORY_CODE, THEORY_RESULT,
                    Path(__file__), Path(__file__).with_name("scalar_reference.py")]},
        "fits_checked": len(records),
        "max_prediction_difference": max(row["max_abs_prediction_error"] for row in records),
        "max_chi2_difference": max(row["abs_chi2_error"] for row in records),
        "max_alpha_profile_difference": max(row["abs_alpha_profile_error"] for row in records),
        "max_fixed_baseline_difference": max(row["max_abs_delta_prediction_error"] for row in fixed_records),
        "maximum_fixed_baseline_fractional_shift": max(row["max_abs_fractional_shift"] for row in fixed_records),
        "maximum_fixed_baseline_shift_over_marginal_sigma": max(row["max_abs_shift_over_marginal_sigma"] for row in fixed_records),
        "maximum_fixed_baseline_full_covariance_norm": max(row["full_covariance_displacement_norm"] for row in fixed_records),
        "maximum_refitted_shift_over_marginal_sigma": max(row["refitted_displacement"]["max_abs_shift_over_marginal_sigma"] for row in records if row["H0_km_s_Mpc"] is not None),
        "maximum_refitted_full_covariance_norm": max(row["refitted_displacement"]["full_covariance_displacement_norm"] for row in records if row["H0_km_s_Mpc"] is not None),
        "maximum_absolute_refitted_delta_chi2": max(abs(row["delta_chi2_from_baseline"]) for row in records if row["H0_km_s_Mpc"] is not None),
        "optimizer_status_counts": optimizer_status_summary,
        "theory_41_point_grid_checks": theory_checks,
        "profiles": records,
        "fixed_parameter_shifts": fixed_records,
        "limitations": "Checks saved optima and the given grid; no global mathematical optimizer certificate, posterior/evidence/significance or H0 inference. Effective massless-radiation reference only.",
        "runtime_seconds": time.perf_counter() - started,
    }
    output = Path(__file__).with_name("saved_fit_review.json")
    output.write_text(json.dumps(final, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: value for key, value in final.items()
                      if key not in ("profiles", "fixed_parameter_shifts", "sha256")}, indent=2))


if __name__ == "__main__":
    main()

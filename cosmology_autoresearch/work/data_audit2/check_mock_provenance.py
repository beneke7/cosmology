#!/usr/bin/env python3
"""Rebuild selected DESI BAO null mocks and independently check their flat fits.

Run from any directory with the campaign .venv:
  .venv/bin/python work/data_audit2/check_mock_provenance.py

This checker writes only work/data_audit2/provenance_check.json.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[2]
SEEDS = {20260924: (200, (0, 99, 199)), 20260925: (1000, (0, 499, 999))}
MEAN_PATH = ROOT / "context/data/desi_dr2_mean.txt"
COV_PATH = ROOT / "context/data/desi_dr2_cov.txt"
SEED_RESULT = ROOT / "experiments/campaign_seed_bao/result.json"
RESULT_PATHS = {
    20260924: ROOT / "experiments/bao_robustness/result.json",
    20260925: ROOT / "experiments/bao_robustness/null_extension_1000_seed20260925.json",
}
CURVATURE_JSONL = ROOT / "experiments/curvature_screen/mock_realizations.jsonl"
CURVATURE_RESULT = ROOT / "experiments/curvature_screen/result.json"
CURVATURE_RECORD = ROOT / "experiments/curvature_screen/record.json"
OUTPUT_PATH = Path(__file__).resolve().parent / "provenance_check.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_inputs():
    rows = []
    for line_no, raw in enumerate(MEAN_PATH.read_text(encoding="utf-8").splitlines(), 1):
        fields = raw.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) != 3:
            raise ValueError(f"mean line {line_no}: expected 3 fields")
        rows.append((float(fields[0]), float(fields[1]), fields[2]))
    z = np.asarray([row[0] for row in rows], dtype=np.float64)
    value = np.asarray([row[1] for row in rows], dtype=np.float64)
    observable = np.asarray([row[2] for row in rows], dtype=str)
    covariance = np.loadtxt(COV_PATH, dtype=np.float64, ndmin=2)
    if covariance.shape != (len(z), len(z)):
        raise ValueError(f"covariance shape {covariance.shape} does not match {len(z)} rows")
    if not np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-12):
        raise ValueError("covariance is not symmetric")
    chol = np.linalg.cholesky(covariance)
    return z, value, observable, covariance, chol, rows


def unit_flat_prediction(z: np.ndarray, observable: np.ndarray, omega_m: float) -> np.ndarray:
    """Independent adaptive-quadrature flat-LCDM prediction at alpha=1."""
    omega_de = 1.0 - float(omega_m)

    def inv_e(zz: float) -> float:
        zp = 1.0 + zz
        e2 = float(omega_m) * zp**3 + omega_de
        if not np.isfinite(e2) or e2 <= 0.0:
            raise ValueError("nonpositive flat-LCDM E^2")
        return 1.0 / np.sqrt(e2)

    unique_z, row_index = np.unique(z, return_inverse=True)
    chi = np.asarray([
        quad(inv_e, 0.0, float(zz), epsabs=2e-13, epsrel=2e-13, limit=200)[0]
        for zz in unique_z
    ], dtype=np.float64)[row_index]
    dh = np.asarray([inv_e(float(zz)) for zz in z], dtype=np.float64)
    dv = np.cbrt(z * chi**2 * dh)
    out = np.empty_like(z, dtype=np.float64)
    for name, values in (("DM_over_rs", chi), ("DH_over_rs", dh), ("DV_over_rs", dv)):
        out[observable == name] = values[observable == name]
    return out


def profile_flat(y: np.ndarray, z: np.ndarray, observable: np.ndarray, chol: np.ndarray) -> dict:
    yw = solve_triangular(chol, y, lower=True, check_finite=False)

    def at_omega(omega_m: float):
        unit = unit_flat_prediction(z, observable, omega_m)
        qw = solve_triangular(chol, unit, lower=True, check_finite=False)
        alpha = float(np.clip(np.dot(qw, yw) / np.dot(qw, qw), 1e-6, 1e4))
        residual_white = yw - alpha * qw
        score = float(np.dot(residual_white, residual_white))
        return score, alpha, alpha * unit

    grid = np.linspace(0.05, 0.60, 33, dtype=np.float64)
    scores = np.asarray([at_omega(float(om))[0] for om in grid])
    best_grid = int(np.argmin(scores))
    lo = float(grid[max(best_grid - 1, 0)])
    hi = float(grid[min(best_grid + 1, len(grid) - 1)])
    if lo == hi:
        omega_m = lo
    else:
        opt = minimize_scalar(lambda om: at_omega(float(om))[0], bounds=(lo, hi),
                              method="bounded", options={"xatol": 5e-15, "maxiter": 1000})
        candidates = [(float(opt.fun), float(opt.x)), (float(scores[best_grid]), float(grid[best_grid]))]
        score, omega_m = min(candidates)
    score, alpha, prediction = at_omega(omega_m)
    return {"Omega_m": omega_m, "alpha": alpha, "chi2": score,
            "prediction": prediction.tolist(), "grid_step": float(grid[1] - grid[0])}


def main() -> int:
    z, published_mean, observable, covariance, chol, mean_rows = load_inputs()
    source_seed = json.loads(SEED_RESULT.read_text(encoding="utf-8"))
    flat_seed = next(row for row in source_seed["models"] if row["model"] == "lcdm")
    null_alpha = float(flat_seed["alpha"])
    null_omega = float(flat_seed["Omega_m"])
    saved_null_mean = np.asarray(flat_seed["prediction"], dtype=np.float64)
    independent_null_mean = null_alpha * unit_flat_prediction(z, observable, null_omega)

    records = {}
    for raw in CURVATURE_JSONL.read_text(encoding="utf-8").splitlines():
        row = json.loads(raw)
        key = (int(row["seed"]), int(row["index"]))
        if key in records:
            raise ValueError(f"duplicate curvature record {key}")
        records[key] = row

    selected = []
    stream_summaries = {}
    all_normal_rows: set[bytes] = set()
    all_rows_unique = True
    for seed, (count, indices) in SEEDS.items():
        doc = json.loads(RESULT_PATHS[seed].read_text(encoding="utf-8"))
        cal = doc["adaptive_search_null_mock_calibration"]
        if int(cal["random_seed"]) != seed or int(cal["n_mocks"]) != count:
            raise ValueError(f"seed/count metadata mismatch for {seed}")
        rng = np.random.default_rng(seed)
        if type(rng.bit_generator).__name__ != "PCG64":
            raise ValueError("expected NumPy default_rng to use PCG64")
        normals = rng.standard_normal((count, len(z)))
        normal_digest = hashlib.sha256(np.ascontiguousarray(normals, dtype="<f8").tobytes()).hexdigest()
        curvature_result = json.loads(CURVATURE_RESULT.read_text(encoding="utf-8"))
        recorded_stream = curvature_result["null_definition"]["random_streams"]["seed_streams"][str(seed)]
        if normal_digest != recorded_stream["normals_sha256_little_endian_float64"]:
            raise ValueError(f"full normal stream digest mismatch for seed {seed}")
        for normal in normals:
            token = np.ascontiguousarray(normal, dtype="<f8").tobytes()
            if token in all_normal_rows:
                all_rows_unique = False
            all_normal_rows.add(token)
        mocks = saved_null_mean[None, :] + normals @ chol.T
        seed_rows = []
        for index in indices:
            key = (seed, index)
            if key not in records:
                raise ValueError(f"missing stored curvature fit record {key}")
            vector = mocks[index]
            fit = profile_flat(vector, z, observable, chol)
            stored = records[key]["flat"]
            stored_pred = np.asarray(stored["prediction"], dtype=np.float64)
            pred_max = float(np.max(np.abs(np.asarray(fit["prediction"]) - stored_pred)))
            score_diff = float(fit["chi2"] - float(stored["chi2"]))
            om_diff = float(fit["Omega_m"] - float(stored["Omega_m"]))
            alpha_diff = float(fit["alpha"] - float(stored["alpha"]))
            row = {
                "seed": seed, "index": index,
                "regenerated_vector_sha256_little_endian_float64": hashlib.sha256(np.ascontiguousarray(vector, dtype="<f8").tobytes()).hexdigest(),
                "vector_max_abs_diff_from_reconstruction_via_independent_scalar_mean": float(np.max(np.abs(vector - (independent_null_mean + normals[index] @ chol.T)))),
                "stored_flat_fit": {"Omega_m": stored["Omega_m"], "alpha": stored["alpha"], "chi2": stored["chi2"]},
                "independent_flat_fit": fit,
                "differences_independent_minus_stored": {"Omega_m": om_diff, "alpha": alpha_diff, "chi2": score_diff,
                                                         "max_abs_prediction": pred_max},
                "stored_fit_success": stored["selected_success"],
                "stored_failed_start_count": stored["failed_start_count"],
            }
            selected.append(row)
            seed_rows.append(row)
        stream_summaries[str(seed)] = {
            "n_mocks": count, "bit_generator": "PCG64",
            "normal_shape": list(normals.shape),
            "normal_stream_sha256_little_endian_float64": normal_digest,
            "recorded_normal_digest_matches": True,
            "indices_checked_zero_based": list(indices),
            "sample_mock_indices_are_consecutive_13-normal_draw_blocks": True,
            "checks": seed_rows,
        }

    data_hashes = {str(path.relative_to(ROOT)): sha256(path) for path in (MEAN_PATH, COV_PATH)}
    code_paths = [ROOT / "scripts/run_bao_robustness.py", ROOT / "scripts/background_bao.py",
                  ROOT / "experiments/curvature_screen/curvature_profile.py"]
    code_hashes = {str(path.relative_to(ROOT)): sha256(path) for path in code_paths}
    source_hash_audit = {}
    for seed, path in RESULT_PATHS.items():
        doc = json.loads(path.read_text(encoding="utf-8"))
        saved = doc["source"]
        saved_data = saved["data_files"]
        saved_code = saved["code_files"]
        script_entry = next((v for k, v in saved_code.items() if k.endswith("run_bao_robustness.py")), None)
        source_hash_audit[str(seed)] = {
            "result_path": str(path.relative_to(ROOT)), "result_sha256": sha256(path),
            "data_hashes_match_current": saved_data == data_hashes,
            "data_hashes_recorded": saved_data,
            "run_bao_robustness_hash_matches_current": script_entry == code_hashes["scripts/run_bao_robustness.py"],
            "run_bao_robustness_hash_recorded": script_entry,
            "background_bao_hash_matches_current": saved_code.get("scripts/background_bao.py") == code_hashes["scripts/background_bao.py"],
            "background_bao_hash_recorded": saved_code.get("scripts/background_bao.py"),
        }

    curvature_record = json.loads(CURVATURE_RECORD.read_text(encoding="utf-8"))
    record_inputs = curvature_record["input_sha256"]
    record_outputs = curvature_record["code_and_output_sha256"]
    curvature_hash_audit = {
        "record_path": str(CURVATURE_RECORD.relative_to(ROOT)),
        "record_sha256": sha256(CURVATURE_RECORD),
        "mean_hash_matches_record": record_inputs.get("context/data/desi_dr2_mean.txt") == data_hashes["context/data/desi_dr2_mean.txt"],
        "covariance_hash_matches_record": record_inputs.get("context/data/desi_dr2_cov.txt") == data_hashes["context/data/desi_dr2_cov.txt"],
        "seed_20260924_result_hash_matches_record": record_inputs.get("experiments/bao_robustness/result.json") == sha256(RESULT_PATHS[20260924]),
        "seed_20260925_result_hash_matches_record": record_inputs.get("experiments/bao_robustness/null_extension_1000_seed20260925.json") == sha256(RESULT_PATHS[20260925]),
        "curvature_code_hash_matches_record": record_outputs.get("experiments/curvature_screen/curvature_profile.py") == code_hashes["experiments/curvature_screen/curvature_profile.py"],
        "mock_jsonl_hash_matches_record": record_outputs.get("experiments/curvature_screen/mock_realizations.jsonl") == sha256(CURVATURE_JSONL),
        "curvature_result_hash_matches_record": record_outputs.get("experiments/curvature_screen/result.json") == sha256(CURVATURE_RESULT),
    }

    eig = np.linalg.eigvalsh(covariance)
    matrix = {
        "shape": list(covariance.shape), "symmetric_max_abs_difference": float(np.max(np.abs(covariance - covariance.T))),
        "cholesky_success": True, "minimum_eigenvalue": float(eig[0]),
        "maximum_eigenvalue": float(eig[-1]), "condition_number_2": float(eig[-1] / eig[0]),
        "row_order_sha256_source_mean": data_hashes["context/data/desi_dr2_mean.txt"],
        "covariance_sha256": data_hashes["context/data/desi_dr2_cov.txt"],
    }
    result = {
        "status": "pending_numeric_thresholds",
        "scope": "CPU-only deterministic provenance audit; no interacting-vacuum mock fits",
        "environment": {"python": sys.version.split()[0], "numpy": np.__version__, "platform": platform.platform()},
        "mock_generator_contract": {
            "random_generator": "np.random.default_rng(seed), NumPy PCG64",
            "draw_call": "standard_normal((n_mocks, 13))",
            "zero_based_indexing": "mock index i is row i of the one complete draw array, equivalent to its consecutive 13-normal row block",
            "mean": {"source": "campaign_seed_bao/result.json fitted flat LCDM prediction", "alpha": null_alpha,
                     "Omega_m": null_omega, "independent_scalar_quadrature_max_abs_difference": float(np.max(np.abs(independent_null_mean - saved_null_mean)))},
            "covariance_transform": "y_i = mean + normals[i] @ L.T; L = lower Cholesky factor of the full ordered 13x13 covariance",
            "seed_streams_are_separately_initialized": True,
            "exact_duplicate_normal_rows_across_both_finite_streams": not all_rows_unique,
            "overlap_caveat": "Distinct seeds and distinct verified PCG64 stream digests establish independently initialized streams for these runs; no finite check proves mathematical non-overlap of infinite PRNG sequences.",
            "streams": stream_summaries,
        },
        "covariance": matrix,
        "source_hash_audit": source_hash_audit,
        "curvature_record_hash_audit": curvature_hash_audit,
        "current_code_hashes": code_hashes,
        "saved_mock_records": {"path": str(CURVATURE_JSONL.relative_to(ROOT)), "rows": len(records),
                                "row_keys_unique": len(records) == 1200,
                                "sha256": sha256(CURVATURE_JSONL),
                                "curvature_result_sha256": sha256(CURVATURE_RESULT)},
        "matched_interaction_search_contract": {
            "data": "same released DESI DR2 13-row mean order and full 13x13 covariance; no Ly-alpha Results IV addition",
            "null": "fixed best-fit flat LCDM mean from campaign_seed_bao/result.json; regenerate and reuse all 1,200 vectors from seeds 20260924 and 20260925",
            "mock_fit_statistic": "Delta chi2 = flat-LCDM profiled score minus one-rate interacting-vacuum profiled score; gamma=0 is nested",
            "interaction_model": "current screen Q=Gamma*rho_x, w_x=-1, baryons separately conserved, dimensionless Gamma/H0; background-only, radiation and perturbations omitted",
            "nuisances": "profile alpha=c/(H0*rd) on [1e-6,1e4] using full covariance; fit Omega_m0; no H0, rd, or baryon fraction inference",
            "physicality": "repeat exact ODE domain and baryon/CDM positivity feasibility rule through z=2.33; preserve every invalid evaluation and optimizer status",
            "primary_domain_status": "resolved_and_frozen_before_observed_fit",
            "primary_Gamma_over_H0_search_domain": [-3.0, 3.0],
            "primary_domain_basis": "source paper Table-I box; declared in the locked compute contract before the observed optimization; nested [-1,1] and [-0.5,0.5] profiles retain the same interior minimum",
            "compact_theory_domain_status": "[-0.2,0.2] is an exploratory alternative only; its independent profile hits the lower boundary and it does not replace the primary null-search domain",
            "selection_warning": "1,200-mock one-model calibration corrects this exact interacting search only; claims selecting among CPL, curvature, and interaction also need the full across-model selection rule repeated per null or clearly remain lane-specific screens",
            "retained_fit_protocol": "rerun flat and interacting fits with the same prespecified start grids, optimizer, physical masks, alpha profiling and retry/failure disposition used on observed data; keep full per-mock fit/status rows and checkpoint output",
            "process_pool": {"recommended_workers": 20, "basis": "20 CPUs are in affinity; realizations are independent CPU ODE/profile jobs; cap BLAS/OpenMP/MKL/NumExpr at one thread and use chunksize 1", "gpu": "not justified for adaptive low-dimensional ODE fits without an end-to-end benchmark"},
            "timing_evidence_seconds": {
                "interaction_observed_39_start_plus_41x41_grid_one_cpu": 12.096865177154541,
                "CPL_selection_1000_mocks_20_workers": 711.5716000488028,
                "flat_curvature_1200_mocks_20_workers": 6.3639099495485425,
                "interpretation": "these timings differ in model cost and fit protocol; they do not justify a precise interaction-mock runtime",
            },
            "bounded_runtime_plan": "run the first 40 exact final mock fits on the 20-worker pool, checkpoint every 25, project total runtime from measured throughput, then resume the remaining 1,160 from those saved rows with explicit per-command timeout/checkpointing",
            "mock_fits_run_by_this_audit": 0,
        },
        "overall_numeric_checks": {
            "all_selected_stored_flat_fits_successful": all(row["stored_fit_success"] for row in selected),
            "max_abs_selected_fit_score_difference": max(abs(row["differences_independent_minus_stored"]["chi2"]) for row in selected),
            "max_abs_selected_fit_prediction_difference": max(row["differences_independent_minus_stored"]["max_abs_prediction"] for row in selected),
            "max_abs_selected_fit_Omega_m_difference": max(abs(row["differences_independent_minus_stored"]["Omega_m"]) for row in selected),
            "max_abs_selected_fit_alpha_difference": max(abs(row["differences_independent_minus_stored"]["alpha"]) for row in selected),
        },
    }
    source_hashes_ok = all(v["data_hashes_match_current"] and v["run_bao_robustness_hash_matches_current"] and v["background_bao_hash_matches_current"] for v in source_hash_audit.values())
    curvature_hashes_ok = all(curvature_hash_audit[key] for key in (
        "mean_hash_matches_record", "covariance_hash_matches_record",
        "seed_20260924_result_hash_matches_record", "seed_20260925_result_hash_matches_record",
        "curvature_code_hash_matches_record", "mock_jsonl_hash_matches_record",
        "curvature_result_hash_matches_record"))
    numeric = result["overall_numeric_checks"]
    records_ok = len(records) == 1200 and len(records) == len(set(records))
    numeric_ok = (all_rows_unique and numeric["all_selected_stored_flat_fits_successful"]
                  and numeric["max_abs_selected_fit_score_difference"] < 5e-10
                  and numeric["max_abs_selected_fit_prediction_difference"] < 1e-7
                  and numeric["max_abs_selected_fit_Omega_m_difference"] < 1e-7
                  and numeric["max_abs_selected_fit_alpha_difference"] < 1e-6
                  and matrix["minimum_eigenvalue"] > 0.0)
    result["status"] = "pass" if source_hashes_ok and curvature_hashes_ok and records_ok and numeric_ok else "fail"
    result["pass_criteria"] = {
        "both_prior_result_source_hash_sets_match_current": source_hashes_ok,
        "curvature_recorded_input_and_output_hashes_match_current": curvature_hashes_ok,
        "all_1200_fit_records_present_with_unique_keys": records_ok,
        "finite_streams_have_no_exact_duplicate_normal_rows": all_rows_unique,
        "selected_fit_score_abs_difference_lt_5e-10": numeric["max_abs_selected_fit_score_difference"] < 5e-10,
        "selected_prediction_max_difference_lt_1e-7": numeric["max_abs_selected_fit_prediction_difference"] < 1e-7,
        "covariance_positive_definite": matrix["minimum_eigenvalue"] > 0.0,
    }
    OUTPUT_PATH.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "output": str(OUTPUT_PATH.relative_to(ROOT)),
                      "selected": [{"seed": r["seed"], "index": r["index"], "chi2_diff": r["differences_independent_minus_stored"]["chi2"],
                                   "max_prediction_diff": r["differences_independent_minus_stored"]["max_abs_prediction"]} for r in selected],
                      "source_hashes_match": all(v["data_hashes_match_current"] and v["run_bao_robustness_hash_matches_current"] and v["background_bao_hash_matches_current"] for v in source_hash_audit.values())}, indent=2))
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

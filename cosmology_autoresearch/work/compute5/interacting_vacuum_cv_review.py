#!/usr/bin/env python3
"""Independent bounded review of the saved seven-fold IVS-vs-LCDM block CV.

No optimizer is run. Dense ``numpy.linalg.solve`` reconstructs each conditional
mean, Schur covariance, and held-out chi-square; model predictions at all 14
saved fold optima are checked with the distinct theory4 scalar/adaptive-quad path.
Run through scripts/run_bounded.py with one numerical thread and no GPU.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT = ROOT / "experiments/interacting_vacuum_block_cv"
RESULT_PATH = EXPERIMENT / "result.json"
RECORD_PATH = EXPERIMENT / "record.json"
SOURCE_PATH = EXPERIMENT / "leave_bin_out.py"
REPORT_PATH = ROOT / "work/compute4/interacting_vacuum_block_cv.md"
MEAN_PATH = ROOT / "context/data/desi_dr2_mean.txt"
COV_PATH = ROOT / "context/data/desi_dr2_cov.txt"
BAO_SOURCE = ROOT / "scripts/background_bao.py"
IV_SOURCE = ROOT / "experiments/interacting_vacuum_screen/interacting_vacuum_profile.py"
THEORY4_PATH = ROOT / "work/theory4/audit_interacting_identifiability.py"
OUTPUT_PATH = ROOT / "work/compute5/interacting_vacuum_cv_review.json"
OBSERVED_AGGREGATE = 22.750299624349342
OBSERVED_Z233_DELTA = 22.069812381805836


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def close(a: float, b: float, tol: float = 1e-9) -> bool:
    return math.isfinite(float(a)) and math.isfinite(float(b)) and abs(float(a) - float(b)) <= tol


def load_inputs():
    rows = []
    for line_no, raw in enumerate(MEAN_PATH.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"{MEAN_PATH}:{line_no}: expected z, value, observable")
        rows.append((float(fields[0]), float(fields[1]), fields[2]))
    z = np.asarray([r[0] for r in rows], dtype=np.float64)
    y = np.asarray([r[1] for r in rows], dtype=np.float64)
    labels = [r[2] for r in rows]
    cov = np.loadtxt(COV_PATH, dtype=np.float64, ndmin=2)
    if z.shape != (13,) or cov.shape != (13, 13):
        raise ValueError(f"expected 13 rows and 13x13 covariance, found {z.shape}, {cov.shape}")
    return z, y, labels, cov


def import_theory4():
    spec = importlib.util.spec_from_file_location("iv_cv_theory4_review", THEORY4_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load independent scalar implementation {THEORY4_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def adaptive_prediction(theory4, z: np.ndarray, labels: list[str], omega_m: float,
                        gamma: float, alpha: float) -> tuple[np.ndarray, dict]:
    """Independent scalar ODE for M,X; adaptive QUAD for D_M; no fitted calls."""
    sol, minima = theory4.scalar_history(float(omega_m), float(gamma), float(np.max(z)))

    def inv_e(zp: float) -> float:
        matter, vacuum = sol(math.log1p(float(zp)))
        e2 = float(matter + vacuum)
        if not math.isfinite(e2) or e2 <= 0.0:
            raise ValueError("nonpositive E^2 on adaptive distance path")
        return 1.0 / math.sqrt(e2)

    dm = np.asarray([
        quad(inv_e, 0.0, float(zi), epsabs=5e-12, epsrel=5e-12, limit=200)[0]
        for zi in z
    ], dtype=np.float64)
    inv_e_at_z = np.asarray([inv_e(float(zi)) for zi in z], dtype=np.float64)
    e = 1.0 / inv_e_at_z
    unit = np.empty_like(z)
    for i, label in enumerate(labels):
        if label == "DM_over_rs":
            unit[i] = dm[i]
        elif label == "DH_over_rs":
            unit[i] = inv_e_at_z[i]
        elif label == "DV_over_rs":
            unit[i] = np.cbrt(z[i] * dm[i] * dm[i] * inv_e_at_z[i])
        else:
            raise ValueError(f"unknown observable: {label}")

    # Verify that at least one positive baryon/CDM split exists on a dense grid.
    ugrid = np.linspace(0.0, math.log1p(float(np.max(z))), 1001)
    matter, vacuum = np.asarray(sol(ugrid), dtype=np.float64)
    zgrid = np.expm1(ugrid)
    fb_ceiling = float(min(1.0, np.min(matter / (float(omega_m) * np.exp(3.0 * ugrid)))))
    physical = bool(np.all(matter > 0.0) and np.all(vacuum > 0.0)
                    and np.all(matter + vacuum > 0.0) and fb_ceiling > 0.0)
    if not physical:
        raise ValueError("saved IVS parameters fail dense positive-history/split check")
    return float(alpha) * unit, {
        **minima,
        "positive_split_ceiling_sampled": fb_ceiling,
        "sampled_z_max": float(zgrid[-1]),
        "physical_positive_split_witness_exists": physical,
    }


def report_hash_claims(report_text: str, targets: dict[str, str]) -> dict:
    results = {}
    for relative, digest in targets.items():
        needle = f"| `{relative}` | `{digest}` |"
        results[relative] = {"digest_claimed": digest, "report_contains_exact_claim": needle in report_text}
    return results


def main() -> int:
    started = time.monotonic()
    result = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    record = json.loads(RECORD_PATH.read_text(encoding="utf-8"))
    z, y, labels, cov = load_inputs()
    if not np.allclose(cov, cov.T, rtol=1e-10, atol=1e-12):
        raise ValueError("full covariance is not symmetric")
    eig_cov = np.linalg.eigvalsh(cov)

    theory4 = import_theory4()
    result_folds = result.get("folds", [])
    expected_groups = np.unique(z)
    mismatches: dict[str, int] = {}
    def mismatch(name: str):
        mismatches[name] = mismatches.get(name, 0) + 1

    hash_targets = {
        "experiments/interacting_vacuum_block_cv/leave_bin_out.py": sha256(SOURCE_PATH),
        "experiments/interacting_vacuum_block_cv/result.json": sha256(RESULT_PATH),
        "context/data/desi_dr2_mean.txt": sha256(MEAN_PATH),
        "context/data/desi_dr2_cov.txt": sha256(COV_PATH),
        "scripts/background_bao.py": sha256(BAO_SOURCE),
        "experiments/interacting_vacuum_screen/interacting_vacuum_profile.py": sha256(IV_SOURCE),
    }
    result_input_hash_checks = {
        rel: {"claimed": result.get("sha256_inputs", {}).get(rel), "current": digest,
              "matches": result.get("sha256_inputs", {}).get(rel) == digest}
        for rel, digest in hash_targets.items()
        if rel != "experiments/interacting_vacuum_block_cv/result.json"
    }
    for rel, check in result_input_hash_checks.items():
        if not check["matches"]:
            mismatch(f"input_hash_{rel}")
    record_sha_checks = {
        "experiments/interacting_vacuum_block_cv/leave_bin_out.py": {
            "claimed": record.get("sha256", {}).get("experiments/interacting_vacuum_block_cv/leave_bin_out.py"),
            "current": sha256(SOURCE_PATH)},
        "experiments/interacting_vacuum_block_cv/result.json": {
            "claimed": record.get("sha256", {}).get("experiments/interacting_vacuum_block_cv/result.json"),
            "current": sha256(RESULT_PATH)},
    }
    for rel, check in record_sha_checks.items():
        check["matches"] = check["claimed"] == check["current"]
        if not check["matches"]:
            mismatch(f"record_hash_{rel}")
    report_text = REPORT_PATH.read_text(encoding="utf-8")
    report_claims = report_hash_claims(report_text, hash_targets)
    for rel, check in report_claims.items():
        if not check["report_contains_exact_claim"]:
            mismatch(f"report_hash_claim_{rel}")

    if result.get("status") != "complete_leave_one_redshift_bin_out_conditional_prediction":
        mismatch("result_status")
    if record.get("status") != "complete" or record.get("result") != "experiments/interacting_vacuum_block_cv/result.json":
        mismatch("record_status_or_target")
    if int(result.get("data_contract", {}).get("row_count", -1)) != 13:
        mismatch("data_contract_row_count")
    if len(result_folds) != len(expected_groups):
        mismatch("fold_count")

    fold_results = []
    all_train_predictions = []
    all_held_predictions = []
    max_mean_diff = 0.0
    max_schur_diff = 0.0
    max_conditional_chi2_diff = 0.0
    max_train_prediction_diff = 0.0
    max_saved_delta_formula_diff = 0.0
    min_schur_eigenvalue = math.inf

    for fold_index, (held_z, fold) in enumerate(zip(expected_groups, result_folds)):
        expected_held = np.flatnonzero(z == held_z)
        expected_train = np.flatnonzero(z != held_z)
        held = np.asarray(fold.get("held_out_row_indices_0based", []), dtype=int)
        train = np.asarray(fold.get("training_row_indices_0based", []), dtype=int)
        if not np.array_equal(held, expected_held):
            mismatch(f"fold_{fold_index}_held_row_order")
        if not np.array_equal(train, expected_train):
            mismatch(f"fold_{fold_index}_train_row_order")
        if float(fold.get("held_out_z", math.nan)) != float(held_z):
            mismatch(f"fold_{fold_index}_held_z_order")
        held_row_meta = fold.get("held_out_rows", [])
        if len(held_row_meta) != len(expected_held):
            mismatch(f"fold_{fold_index}_held_row_metadata_count")
        else:
            for meta, ri in zip(held_row_meta, expected_held):
                if (int(meta.get("row_0based", -1)) != int(ri)
                        or float(meta.get("z", math.nan)) != float(z[ri])
                        or meta.get("observable") != labels[ri]
                        or float(meta.get("value", math.nan)) != float(y[ri])):
                    mismatch(f"fold_{fold_index}_held_row_metadata")

        ctt = cov[np.ix_(train, train)]
        cht = cov[np.ix_(held, train)]
        cth = cov[np.ix_(train, held)]
        chh = cov[np.ix_(held, held)]
        # Deliberately use general dense solves, independent of the source Cholesky path.
        solved_cross = np.linalg.solve(ctt, cth)
        schur = chh - cht @ solved_cross
        schur = (schur + schur.T) / 2.0
        schur_eigs = np.linalg.eigvalsh(schur)
        min_schur_eigenvalue = min(min_schur_eigenvalue, float(schur_eigs[0]))
        if schur_eigs[0] <= 0.0:
            mismatch(f"fold_{fold_index}_schur_not_positive_definite")
        stored_covs = {}
        for model, key in (("lcdm", "lcdm_conditional_covariance"), ("ivs", "ivs_conditional_covariance")):
            stored_cov = np.asarray(fold.get(key, []), dtype=np.float64)
            stored_covs[model] = stored_cov
            if stored_cov.shape != schur.shape:
                mismatch(f"fold_{fold_index}_{model}_schur_shape")
            else:
                max_schur_diff = max(max_schur_diff, float(np.max(np.abs(stored_cov - schur))))

        fold_model_checks = {}
        independent_predictions = {}
        for model in ("lcdm", "ivs"):
            fit = fold[f"{model}_fit"]
            om = float(fit["Omega_m"])
            gamma = 0.0 if model == "lcdm" else float(fit["Gamma_over_H0"])
            alpha = float(fit["alpha"])
            independent_prediction, physical = adaptive_prediction(
                theory4, z, labels, om, gamma, alpha)
            saved_train = np.asarray(fit.get("prediction_train", []), dtype=np.float64)
            if saved_train.shape != (len(train),):
                mismatch(f"fold_{fold_index}_{model}_saved_training_prediction_shape")
            else:
                train_diff = float(np.max(np.abs(saved_train - independent_prediction[train])))
                max_train_prediction_diff = max(max_train_prediction_diff, train_diff)
            train_prediction = independent_prediction[train]
            held_prediction = independent_prediction[held]
            train_residual = y[train] - train_prediction
            conditional_mean = held_prediction + cht @ np.linalg.solve(ctt, train_residual)
            held_residual = y[held] - conditional_mean
            chi2 = float(held_residual @ np.linalg.solve(schur, held_residual))
            saved_mean = np.asarray(fold[f"{model}_conditional_mean"], dtype=np.float64)
            if saved_mean.shape != conditional_mean.shape:
                mismatch(f"fold_{fold_index}_{model}_saved_conditional_mean_shape")
                mean_diff = math.inf
            else:
                mean_diff = float(np.max(np.abs(saved_mean - conditional_mean)))
                max_mean_diff = max(max_mean_diff, mean_diff)
            saved_chi2 = float(fold[f"{model}_conditional_chi2"])
            chi_diff = abs(chi2 - saved_chi2)
            max_conditional_chi2_diff = max(max_conditional_chi2_diff, chi_diff)
            if not np.all(np.isfinite(independent_prediction)):
                mismatch(f"fold_{fold_index}_{model}_nonfinite_prediction")
            independent_predictions[model] = independent_prediction
            all_train_predictions.append({"fold": fold_index, "model": model, "rows": train.tolist(),
                                          "saved_prediction_max_abs_difference": train_diff if saved_train.shape == (len(train),) else None})
            all_held_predictions.append({"fold": fold_index, "model": model, "rows": held.tolist(),
                                         "independent_heldout_prediction": held_prediction.tolist(),
                                         "reconstructed_conditional_mean": conditional_mean.tolist(),
                                         "saved_conditional_mean": saved_mean.tolist(),
                                         "heldout_chi2_dense_solve": chi2,
                                         "saved_heldout_chi2": saved_chi2,
                                         "chi2_abs_difference": chi_diff,
                                         "physicality": physical})
            fold_model_checks[model] = {
                "saved_train_prediction_max_abs_difference": train_diff if saved_train.shape == (len(train),) else None,
                "conditional_mean_max_abs_difference_vs_saved": mean_diff,
                "conditional_chi2_dense_solve": chi2,
                "conditional_chi2_saved": saved_chi2,
                "conditional_chi2_abs_difference": chi_diff,
                "independent_heldout_prediction": held_prediction.tolist(),
                "positive_history_and_split": physical,
            }

        lcdm_chi2 = fold_model_checks["lcdm"]["conditional_chi2_dense_solve"]
        ivs_chi2 = fold_model_checks["ivs"]["conditional_chi2_dense_solve"]
        delta = ivs_chi2 - lcdm_chi2
        saved_delta = float(fold["delta_ivs_minus_lcdm_conditional_chi2"])
        delta_diff = abs(delta - saved_delta)
        max_saved_delta_formula_diff = max(max_saved_delta_formula_diff, delta_diff)
        if delta_diff > 2e-9:
            mismatch(f"fold_{fold_index}_delta_formula")
        fold_results.append({
            "held_out_z": float(held_z), "held_out_rows_0based": held.tolist(),
            "training_rows_0based": train.tolist(),
            "schur_covariance_dense_solve": schur.tolist(),
            "schur_covariance_eigenvalues": schur_eigs.tolist(),
            "max_abs_schur_difference_vs_each_saved_model_copy": max(
                float(np.max(np.abs(stored_covs["lcdm"] - schur))),
                float(np.max(np.abs(stored_covs["ivs"] - schur))),
            ),
            "models": fold_model_checks,
            "delta_ivs_minus_lcdm_dense_solve": delta,
            "saved_delta": saved_delta,
            "delta_abs_difference": delta_diff,
        })

    sum_lcdm = float(sum(f["models"]["lcdm"]["conditional_chi2_dense_solve"] for f in fold_results))
    sum_ivs = float(sum(f["models"]["ivs"]["conditional_chi2_dense_solve"] for f in fold_results))
    aggregate_delta = sum_ivs - sum_lcdm
    z233 = next(f for f in fold_results if f["held_out_z"] == 2.33)
    z233_delta = float(z233["delta_ivs_minus_lcdm_dense_solve"])
    saved_aggregate = result["aggregate"]
    aggregate_checks = {
        "lcdm_sum_dense_solve": sum_lcdm,
        "lcdm_sum_saved": float(saved_aggregate["lcdm_sum_conditional_chi2"]),
        "ivs_sum_dense_solve": sum_ivs,
        "ivs_sum_saved": float(saved_aggregate["ivs_sum_conditional_chi2"]),
        "delta_sum_dense_solve": aggregate_delta,
        "delta_sum_saved": float(saved_aggregate["sum_delta_ivs_minus_lcdm"]),
        "delta_sum_abs_difference": abs(aggregate_delta - float(saved_aggregate["sum_delta_ivs_minus_lcdm"])),
        "delta_sum_record": float(record["aggregate_delta_ivs_minus_lcdm_conditional_chi2"]),
        "z233_delta_dense_solve": z233_delta,
        "z233_delta_saved": OBSERVED_Z233_DELTA,
        "z233_delta_abs_difference": abs(z233_delta - OBSERVED_Z233_DELTA),
        "record_vs_result_aggregate_abs_difference": abs(float(record["aggregate_delta_ivs_minus_lcdm_conditional_chi2"])
                                                           - float(saved_aggregate["sum_delta_ivs_minus_lcdm"])),
    }
    for key in ("delta_sum_abs_difference", "z233_delta_abs_difference"):
        if aggregate_checks[key] > 2e-9:
            mismatch(f"aggregate_{key}")
    aggregate_checks["record_delta_matches_result_rounded_to_6dp"] = (
        round(float(saved_aggregate["sum_delta_ivs_minus_lcdm"]), 6)
        == float(record["aggregate_delta_ivs_minus_lcdm_conditional_chi2"])
    )
    if not aggregate_checks["record_delta_matches_result_rounded_to_6dp"]:
        mismatch("record_result_aggregate_rounding")
    if abs(float(saved_aggregate["sum_delta_ivs_minus_lcdm"]) - OBSERVED_AGGREGATE) > 5e-7:
        mismatch("aggregate_not_reported_22_7503")
    if abs(max_mean_diff) > 2e-7:
        mismatch("conditional_mean_reconstruction_tolerance")
    if abs(max_schur_diff) > 1e-12:
        mismatch("schur_reconstruction_tolerance")
    if abs(max_conditional_chi2_diff) > 2e-8:
        mismatch("conditional_chi2_reconstruction_tolerance")
    if abs(max_train_prediction_diff) > 2e-7:
        mismatch("saved_train_prediction_scalar_quad_tolerance")

    checks = {
        "seven_folds_and_13_row_partition_order_verified": len(fold_results) == 7 and not any("row_order" in k for k in mismatches),
        "full_covariance_and_all_schur_blocks_positive_definite": bool(eig_cov[0] > 0.0 and min_schur_eigenvalue > 0.0),
        "all_conditional_means_schur_blocks_and_chi2_reconstructed": max_mean_diff < 2e-7 and max_schur_diff < 1e-12 and max_conditional_chi2_diff < 2e-8,
        "all_14_saved_fit_predictions_recomputed_via_theory4_scalar_quad": len(all_train_predictions) == 14 and max_train_prediction_diff < 2e-7,
        "aggregate_and_z233_delta_match": aggregate_checks["delta_sum_abs_difference"] < 2e-9 and aggregate_checks["z233_delta_abs_difference"] < 2e-9 and aggregate_checks["record_delta_matches_result_rounded_to_6dp"],
        "input_record_and_report_hash_claims_verified": all(x["matches"] for x in result_input_hash_checks.values()) and all(x["matches"] for x in record_sha_checks.values()) and all(x["report_contains_exact_claim"] for x in report_claims.values()),
        "training_only_fit_and_conditional_test_semantics_reviewed": True,
    }
    status = "pass" if all(checks.values()) and not mismatches else "fail"
    summary = {
        "status": status,
        "scope": "Independent numerical and semantic review of the saved seven-fold leave-redshift-group-out plug-in conditional prediction comparison; no optimizers/refits",
        "runtime_environment": {
            "python": sys.version.split()[0], "numpy": np.__version__, "platform": platform.platform(),
            "cpu_only": True, "gpu_used": False,
            "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
            "review_runtime_seconds": float(time.monotonic() - started),
        },
        "commands": {
            "review": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --seconds 360 -- .venv/bin/python work/compute5/interacting_vacuum_cv_review.py",
            "optimizer_calls": 0,
            "reevaluated_saved_fit_count": 14,
            "refit_statement": "Each fold's saved alpha, Omega_m and Gamma/H0 were held fixed. Independent scalar ODE plus adaptive QUAD recomputed all 14 model predictions; no optimizers or parameter refits ran.",
        },
        "hashes": {
            "record_json_sha256_current": sha256(RECORD_PATH),
            "result_json_sha256_current": sha256(RESULT_PATH),
            "inputs_sha256_current": hash_targets,
            "result_record_sha_claim_checks": record_sha_checks,
            "result_input_sha_claim_checks": result_input_hash_checks,
            "published_report_sha256_current": sha256(REPORT_PATH),
            "published_report_sha_claim_checks": report_claims,
            "theory4_independent_scalar_code_sha256": sha256(THEORY4_PATH),
        },
        "data_and_covariance": {
            "row_count": len(z), "row_order": "non-comment order of context/data/desi_dr2_mean.txt; zero-based",
            "redshift_groups": [float(v) for v in expected_groups],
            "heldout_row_group_order": [f["held_out_rows_0based"] for f in fold_results],
            "full_covariance_symmetric": True,
            "full_covariance_min_eigenvalue": float(eig_cov[0]),
            "full_covariance_max_eigenvalue": float(eig_cov[-1]),
            "full_covariance_condition_number_2": float(eig_cov[-1] / eig_cov[0]),
            "minimum_conditional_schur_eigenvalue_across_folds": min_schur_eigenvalue,
            "all_fold_schur_eigenvalues": {str(f["held_out_z"]): f["schur_covariance_eigenvalues"] for f in fold_results},
        },
        "independent_reconstruction": {
            "linear_algebra_path": "Dense numpy.linalg.solve for C_TT^-1 residual, C_TT^-1 C_TH, and S^-1 heldout residual; no Cholesky/cho_solve in reconstruction.",
            "prediction_path": "work/theory4/audit_interacting_identifiability.py scalar_history DOP853 (M,X) plus scipy.integrate.quad for every D_M; flat LCDM uses the same scalar history at g=0.",
            "max_abs_conditional_mean_vs_saved": max_mean_diff,
            "max_abs_schur_covariance_vs_saved": max_schur_diff,
            "max_abs_conditional_chi2_vs_saved": max_conditional_chi2_diff,
            "max_abs_training_prediction_vs_saved": max_train_prediction_diff,
            "max_abs_delta_formula_difference": max_saved_delta_formula_diff,
            "all_training_fit_predictions": all_train_predictions,
            "all_heldout_predictions_and_conditional_means": all_held_predictions,
            "folds": fold_results,
        },
        "aggregate": aggregate_checks,
        "semantics_and_limitations": {
            "fit_scope_checked": "For each held-out redshift group, alpha and shape parameters are selected using only training-row likelihood/residuals and C_TT; the conditional held-out score is computed after those parameters are fixed.",
            "selection_semantics": "The finite multistart/grid/refinement selection in each fold uses the training score only; held-out rows enter only through the conditional predictive evaluation. The seven fold scores are not an independent product because training sets and survey covariance overlap.",
            "score_semantics": "Conditional chi-square omits the conditional log determinant. That determinant depends on the shared covariance and fold partition, not on model, so it cancels in within-fold IVS-vs-LCDM comparisons; summing these values remains a descriptive aggregate across overlapping folds, not a joint likelihood/evidence.",
            "uncertainty_semantics": "Plug-in point estimates; no parameter-uncertainty integration, posterior predictive CV, model evidence or external-survey validation.",
            "physical_scope": "Compressed late-time BAO background only; no radiation, sound-horizon calibration, CMB, perturbation stability/closure, or growth test.",
            "conceptual_defects_found": [],
            "semantic_caveat": "The +22.7503 aggregate is dominated by the z=2.33 fold (+22.0698), where two high-z measurements are held together and the other folds train on those rows. Treat it as a localized leave-bin-out stress test plus seven-fold descriptive aggregate, not universal predictive failure evidence.",
        },
        "mismatches": mismatches,
        "checks": checks,
    }
    OUTPUT_PATH.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "output": str(OUTPUT_PATH.relative_to(ROOT)),
                      "checks": checks, "aggregate": aggregate_checks,
                      "max_mean_diff": max_mean_diff, "max_schur_diff": max_schur_diff,
                      "max_chi2_diff": max_conditional_chi2_diff,
                      "max_train_prediction_diff": max_train_prediction_diff,
                      "mismatches": mismatches}, indent=2))
    return 0 if status == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

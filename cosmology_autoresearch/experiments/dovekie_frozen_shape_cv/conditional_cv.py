#!/usr/bin/env python3
"""Frozen DESI-BAO cosmological shapes, conditionally scored on DES-Dovekie SNe.

No cosmological shape is fit on the SN sample. On each predeclared redshift fold,
only a common SN magnitude intercept is learned from training rows. The release
supplies precision P=C^{-1}, so the code uses exact precision-block identities
and independently reconstructs one fold through the explicit covariance C.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import cho_solve


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "work/inference"))
sys.path.insert(0, str(ROOT / "experiments/interacting_vacuum_screen"))
sys.path.insert(0, str(ROOT / "work/theory4"))

import fit_dovekie as sn  # noqa: E402
import interacting_vacuum_profile as iv  # noqa: E402
import audit_interacting_identifiability as scalar_iv  # noqa: E402


OUT = ROOT / "experiments/dovekie_frozen_shape_cv"
RESULT_PATH = OUT / "result.json"
BAO_RESULT = ROOT / "experiments/interacting_vacuum_screen/result.json"
N_FOLDS = 4
REFERENCE_FOLD = 1
FB_PHYSICAL_WITNESS = 0.16
LOCKED_SHAPE_TARGETS = {
    "BAO_flat_LCDM": {"Omega_m": 0.2974618150, "g": 0.0},
    "BAO_interacting_vacuum": {"Omega_m": 0.3869710077, "g": -0.4666647299},
}
EDGE_SHIFT_FRACTION = 0.01
DISTANCE_CHECK_TOLERANCE_MAG = 2e-8
EXPECTED_CID_ORDER_SHA256 = "b7d5c6ad8dfdadf1e12443852b14006c0d9d0ceb4bcb46efb378f058ce5aa5a3"


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_numpy_to_native(value: Any) -> Any:
    """Serialize only explicit NumPy containers/scalars, preserving their values."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"unsupported JSON value type: {type(value).__name__}")


def z_quantile_bins(
    z: np.ndarray,
    n_folds: int = N_FOLDS,
    target_shifts_rows: dict[int, float] | None = None,
) -> tuple[list[float], list[np.ndarray]]:
    """Make contiguous zHD bins; optional fixed row-target shifts never split ties."""
    unique, counts = np.unique(np.asarray(z, dtype=np.float64), return_counts=True)
    if unique.size < n_folds:
        raise ValueError("fewer unique zHD values than requested folds")
    target_shifts_rows = target_shifts_rows or {}
    if any(not 1 <= key < n_folds for key in target_shifts_rows):
        raise ValueError("target shift key must name an internal one-based fold boundary")
    cumulative = np.cumsum(counts)
    cut_positions: list[int] = []
    for fold_number in range(1, n_folds):
        low = cut_positions[-1] + 1 if cut_positions else 1
        high = unique.size - (n_folds - fold_number) + 1
        if low >= high:
            raise ValueError("cannot form nonempty redshift folds")
        target = z.size * fold_number / n_folds + target_shifts_rows.get(fold_number, 0.0)
        candidates = range(low, high)
        cut = min(candidates, key=lambda p: (abs(float(cumulative[p - 1]) - target), p))
        cut_positions.append(cut)
    edges = [
        float(0.5 * (unique[p - 1] + unique[p]))
        for p in cut_positions
    ]
    labels = np.searchsorted(np.asarray(edges), z, side="right")
    groups = [np.flatnonzero(labels == i) for i in range(n_folds)]
    if any(group.size == 0 for group in groups):
        raise ValueError("redshift partition contains an empty fold")
    return edges, groups


def load_frozen_bao_shapes() -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """Load exact selected BAO-only coordinates and verify the predeclared lock."""
    source = json.loads(BAO_RESULT.read_text(encoding="utf-8"))
    fits = source["fits"]
    flat = fits["Gamma0_flat_LCDM"]
    ivs = fits["interacting_vacuum"]
    for name, fit in (("BAO_flat_LCDM", flat), ("BAO_interacting_vacuum", ivs)):
        if fit.get("selected_status") != 0 or fit.get("selected_success") is not True:
            raise ValueError(f"BAO source point {name} is not a selected successful optimizer result")
    shapes = {
        "BAO_flat_LCDM": {
            "model": "LCDM",
            "Omega_m": float(flat["Omega_m"]),
            "g": float(flat["Gamma_over_H0"]),
        },
        "BAO_interacting_vacuum": {
            "model": "IVS",
            "Omega_m": float(ivs["Omega_m"]),
            "g": float(ivs["Gamma_over_H0"]),
        },
    }
    for name, target in LOCKED_SHAPE_TARGETS.items():
        for parameter in ("Omega_m", "g"):
            if abs(shapes[name][parameter] - target[parameter]) > 1e-8:
                raise ValueError(
                    f"exact BAO result moved outside the predeclared frozen shape: "
                    f"{name}.{parameter}={shapes[name][parameter]} vs {target[parameter]}"
                )
    if shapes["BAO_flat_LCDM"]["g"] != 0.0:
        raise ValueError("LCDM source point is not exactly the nested zero-interaction model")
    return shapes, {
        "source_path": str(BAO_RESULT.relative_to(ROOT)),
        "predeclared_rounded_targets": LOCKED_SHAPE_TARGETS,
        "exact_selected_points": shapes,
        "selected_fit_chi2": {
            "BAO_flat_LCDM": float(flat["chi2"]),
            "BAO_interacting_vacuum": float(ivs["chi2"]),
        },
        "selected_optimizer_status": {
            "BAO_flat_LCDM": int(flat["selected_status"]),
            "BAO_interacting_vacuum": int(ivs["selected_status"]),
        },
    }


def luminosity_modulus(z_hd: np.ndarray, z_hel: np.ndarray, dc_dimensionless: np.ndarray) -> np.ndarray:
    dl_mpc = (1.0 + z_hel) * (sn.C_KM_S / sn.H0_GAUGE_KM_S_MPC) * dc_dimensionless
    if np.any(~np.isfinite(dl_mpc)) or np.any(dl_mpc <= 0.0):
        raise ValueError("nonpositive or nonfinite luminosity distance")
    return 5.0 * np.log10(dl_mpc) + 25.0


def predictions(
    inputs: dict[str, Any], shapes: dict[str, dict[str, Any]]
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    z_hd = inputs["z_hd"]
    z_hel = inputs["z_hel"]
    base: dict[str, np.ndarray] = {}
    checks: dict[str, Any] = {}

    lcdm_pars = np.asarray([shapes["BAO_flat_LCDM"]["Omega_m"]])
    lcdm_main = sn.background_mu(z_hd, z_hel, lcdm_pars, "LCDM")

    def inverse_e_lcdm(z: float) -> float:
        return 1.0 / np.sqrt(
            shapes["BAO_flat_LCDM"]["Omega_m"] * (1.0 + z) ** 3
            + 1.0 - shapes["BAO_flat_LCDM"]["Omega_m"]
        )

    dc_lcdm_quad = np.asarray(
        [quad(inverse_e_lcdm, 0.0, float(z), epsabs=5e-12, epsrel=5e-12, limit=200)[0] for z in z_hd]
    )
    lcdm_independent = luminosity_modulus(z_hd, z_hel, dc_lcdm_quad)
    base["BAO_flat_LCDM"] = lcdm_main
    checks["LCDM_max_independent_mu_abs_difference_mag"] = float(
        np.max(np.abs(lcdm_main - lcdm_independent))
    )

    ivs = shapes["BAO_interacting_vacuum"]
    ivs_ode = iv.integrate_profile(float(ivs["Omega_m"]), float(ivs["g"]), z_hd)
    ivs_main = luminosity_modulus(z_hd, z_hel, ivs_ode["distance"])

    independent_solution, independent_minima = scalar_iv.scalar_history(
        float(ivs["Omega_m"]), float(ivs["g"]), float(np.max(z_hd))
    )

    def inverse_e_ivs(z: float) -> float:
        matter, vacuum = independent_solution(np.log1p(float(z)))
        e2 = float(matter + vacuum)
        if not np.isfinite(e2) or e2 <= 0.0:
            raise ValueError("independent scalar IVS implementation has E^2 <= 0")
        return 1.0 / np.sqrt(e2)

    dc_ivs_quad = np.asarray(
        [quad(inverse_e_ivs, 0.0, float(z), epsabs=5e-12, epsrel=5e-12, limit=200)[0] for z in z_hd]
    )
    ivs_independent = luminosity_modulus(z_hd, z_hel, dc_ivs_quad)
    base["BAO_interacting_vacuum"] = ivs_main
    checks["IVS_max_independent_mu_abs_difference_mag"] = float(
        np.max(np.abs(ivs_main - ivs_independent))
    )
    checks["IVS_independent_scalar_minima"] = independent_minima
    checks["IVS_fixed_fb_0p16_physical_through_z2p33"] = iv.fixed_split_check(
        float(ivs["Omega_m"]), float(ivs["g"]), FB_PHYSICAL_WITNESS
    )

    # Exact zero-rate limit against the same flat LCDM distance implementation.
    gamma0 = iv.integrate_profile(float(ivs["Omega_m"]), 0.0, z_hd)
    gamma0_mu = luminosity_modulus(z_hd, z_hel, gamma0["distance"])
    lcdm_same_omega = sn.background_mu(
        z_hd, z_hel, np.asarray([float(ivs["Omega_m"])]), "LCDM"
    )
    checks["Gamma0_max_LCDM_mu_abs_difference_mag"] = float(
        np.max(np.abs(gamma0_mu - lcdm_same_omega))
    )
    checks["all_model_means_finite"] = all(np.all(np.isfinite(v)) for v in base.values())
    if checks["LCDM_max_independent_mu_abs_difference_mag"] > DISTANCE_CHECK_TOLERANCE_MAG:
        raise ValueError("independent LCDM scalar-quadrature distance check failed")
    if checks["IVS_max_independent_mu_abs_difference_mag"] > DISTANCE_CHECK_TOLERANCE_MAG:
        raise ValueError("independent IVS scalar/adaptive distance check failed")
    if checks["Gamma0_max_LCDM_mu_abs_difference_mag"] > DISTANCE_CHECK_TOLERANCE_MAG:
        raise ValueError("Gamma=0 did not recover flat LCDM at equal Omega_m")
    if not checks["all_model_means_finite"]:
        raise ValueError("nonfinite fixed-shape SN distance prediction")
    witness = checks["IVS_fixed_fb_0p16_physical_through_z2p33"]
    if not witness.get("physical_on_0_to_zmax", False):
        raise ValueError("declared f_b=0.16 IVS physicality witness failed through z=2.33")
    return base, checks


def solve_precision_block(
    precision: np.ndarray,
    y: np.ndarray,
    model_base: np.ndarray,
    train: np.ndarray,
    held: np.ndarray,
) -> dict[str, Any]:
    """Exact conditional Gaussian score, integrating the training-fit SN intercept."""
    p_tt = precision[np.ix_(train, train)]
    p_th = precision[np.ix_(train, held)]
    p_ht = precision[np.ix_(held, train)]
    p_hh = precision[np.ix_(held, held)]
    l_hh = np.linalg.cholesky(p_hh)

    def solve_hh(rhs: np.ndarray) -> np.ndarray:
        return cho_solve((l_hh, True), rhs, check_finite=False)

    one_t = np.ones(train.size)
    one_h = np.ones(held.size)
    r0_t = y[train] - model_base[train]

    # The Schur complement is the inverse marginal covariance of training data.
    q_one = p_tt @ one_t - p_th @ solve_hh(p_ht @ one_t)
    q_r0 = p_tt @ r0_t - p_th @ solve_hh(p_ht @ r0_t)
    offset_information = float(one_t @ q_one)
    if not np.isfinite(offset_information) or offset_information <= 0.0:
        raise ValueError("training-only intercept has nonpositive GLS information")
    offset_variance = 1.0 / offset_information
    offset = float(one_t @ q_r0 / offset_information)
    normal_equation_residual = abs(float(one_t @ (q_r0 - offset * q_one)))

    residual_train = y[train] - model_base[train] - offset
    b_one = solve_hh(p_ht @ one_t)
    v = one_h + b_one
    conditional_mean = model_base[held] + offset - solve_hh(p_ht @ residual_train)
    residual_held = y[held] - conditional_mean

    # P_HH is the precision of S_(H|T); Woodbury adds training uncertainty in M.
    s_inv_r = p_hh @ residual_held
    s_inv_v = p_hh @ v
    plugin_chi2 = float(residual_held @ s_inv_r)
    v_sinv_v = float(v @ s_inv_v)
    v_sinv_r = float(v @ s_inv_r)
    woodbury_denominator = 1.0 + offset_variance * v_sinv_v
    if woodbury_denominator <= 0.0 or not np.isfinite(woodbury_denominator):
        raise ValueError("invalid rank-one predictive covariance update")
    integrated_chi2 = float(
        plugin_chi2 - offset_variance * v_sinv_r**2 / woodbury_denominator
    )
    logdet_s = float(-2.0 * np.sum(np.log(np.diag(l_hh))))
    logdet_predictive = float(logdet_s + np.log1p(offset_variance * v_sinv_v))
    heldout_neg2logp = integrated_chi2 + logdet_predictive + held.size * np.log(2.0 * np.pi)
    plugin_neg2logp = plugin_chi2 + logdet_s + held.size * np.log(2.0 * np.pi)
    return {
        "offset_mag": offset,
        "offset_posterior_variance_mag2": offset_variance,
        "offset_information": offset_information,
        "training_intercept_normal_equation_abs_residual": normal_equation_residual,
        "conditional_mean": conditional_mean,
        "conditional_residual": residual_held,
        "rank_one_vector": v,
        "plugin_conditional_chi2": plugin_chi2,
        "integrated_intercept_conditional_chi2": integrated_chi2,
        "integrated_intercept_conditional_neg2logpredictive": float(heldout_neg2logp),
        "plugin_conditional_neg2logpredictive": float(plugin_neg2logp),
        "conditional_logdet_covariance": logdet_s,
        "integrated_predictive_logdet_covariance": logdet_predictive,
        "heldout_precision_cholesky_min_diag": float(np.min(np.diag(l_hh))),
        "woodbury_denominator": woodbury_denominator,
        "heldout_rows": held.tolist(),
        "training_rows_count": int(train.size),
    }


def covariance_reference(
    covariance: np.ndarray,
    precision_result: dict[str, Any],
    y: np.ndarray,
    model_base: np.ndarray,
    train: np.ndarray,
    held: np.ndarray,
) -> dict[str, float]:
    """Direct covariance-domain reconstruction for the predeclared reference fold."""
    c_tt = covariance[np.ix_(train, train)]
    c_ht = covariance[np.ix_(held, train)]
    c_th = covariance[np.ix_(train, held)]
    c_hh = covariance[np.ix_(held, held)]
    l_tt = np.linalg.cholesky(c_tt)

    def solve_tt(rhs: np.ndarray) -> np.ndarray:
        return cho_solve((l_tt, True), rhs, check_finite=False)

    one_t = np.ones(train.size)
    one_h = np.ones(held.size)
    q_one = solve_tt(one_t)
    q_r0 = solve_tt(y[train] - model_base[train])
    information = float(one_t @ q_one)
    variance = 1.0 / information
    offset = float(one_t @ q_r0 / information)
    residual_t = y[train] - model_base[train] - offset
    conditional_mean = model_base[held] + offset + c_ht @ solve_tt(residual_t)
    s_raw = c_hh - c_ht @ solve_tt(c_th)
    conditional_covariance_asymmetry = float(np.max(np.abs(s_raw - s_raw.T)))
    s = 0.5 * (s_raw + s_raw.T)
    l_s = np.linalg.cholesky(s)
    rank_one = one_h - c_ht @ solve_tt(one_t)
    predictive_cov = s + variance * np.outer(rank_one, rank_one)
    l_predictive = np.linalg.cholesky(0.5 * (predictive_cov + predictive_cov.T))
    residual_h = y[held] - conditional_mean
    whitened = np.linalg.solve(l_predictive, residual_h)
    integrated_chi2 = float(whitened @ whitened)
    s_whitened = np.linalg.solve(l_s, residual_h)
    plugin_chi2 = float(s_whitened @ s_whitened)
    logdet_predictive = float(2.0 * np.sum(np.log(np.diag(l_predictive))))
    return {
        "max_intercept_abs_difference": abs(offset - precision_result["offset_mag"]),
        "max_conditional_mean_abs_difference": float(
            np.max(np.abs(conditional_mean - precision_result["conditional_mean"]))
        ),
        "max_rank_one_vector_abs_difference": float(
            np.max(np.abs(rank_one - precision_result["rank_one_vector"]))
        ),
        "integrated_intercept_conditional_chi2_abs_difference": abs(
            integrated_chi2 - precision_result["integrated_intercept_conditional_chi2"]
        ),
        "plugin_conditional_chi2_abs_difference": abs(
            plugin_chi2 - precision_result["plugin_conditional_chi2"]
        ),
        "integrated_predictive_logdet_abs_difference": abs(
            logdet_predictive - precision_result["integrated_predictive_logdet_covariance"]
        ),
        "full_conditional_covariance_max_asymmetry_before_symmetrization": conditional_covariance_asymmetry,
    }


def main() -> None:
    started = time.perf_counter()
    inputs = sn.load_inputs()
    y = inputs["mu_obs"]
    z_hd = inputs["z_hd"]
    n = int(y.size)
    shapes, shape_provenance = load_frozen_bao_shapes()
    model_means, theory_checks = predictions(inputs, shapes)
    cid_order_hash = hashlib.sha256(
        ("\n".join(row["CID"] for row in inputs["records"]) + "\n").encode("utf-8")
    ).hexdigest()
    if n != 1820 or cid_order_hash != EXPECTED_CID_ORDER_SHA256:
        raise ValueError(f"row count/order differs from audited release: n={n}, cid_hash={cid_order_hash}")
    edges, groups = z_quantile_bins(z_hd)

    p = inputs["precision"]
    lp = inputs["precision_cholesky"]
    if p.shape != (n, n) or not np.array_equal(p, p.T):
        raise ValueError("release precision shape/order/symmetry check failed")
    source_files = [
        ROOT / "context/data/des_dovekie_hd.csv",
        ROOT / "context/data/des_dovekie_stat_sys.npz",
        ROOT / "experiments/interacting_vacuum_screen/result.json",
        ROOT / "work/inference/fit_dovekie.py",
        ROOT / "experiments/interacting_vacuum_screen/interacting_vacuum_profile.py",
        ROOT / "work/theory4/audit_interacting_identifiability.py",
        ROOT / "work/data_audit4/dovekie_frozen_shape_contract.md",
        ROOT / "work/data_audit4/dovekie_frozen_shape_contract.json",
        Path(__file__).resolve(),
    ]
    source_hashes = {str(path.relative_to(ROOT)): file_sha256(path) for path in source_files}
    covariance = cho_solve((lp, True), np.eye(n), check_finite=False)
    covariance_asymmetry = float(np.max(np.abs(covariance - covariance.T)))
    covariance = 0.5 * (covariance + covariance.T)
    lc = np.linalg.cholesky(covariance)
    inverse_residual = float(np.max(np.abs(p @ covariance - np.eye(n))))
    if covariance_asymmetry > 1e-10 or inverse_residual > 1e-7:
        raise ValueError("explicit inverse covariance failed symmetry or P*C=I tolerance")

    reference_tolerances = {
        "max_intercept_abs_difference": 1e-6,
        "max_conditional_mean_abs_difference": 2e-5,
        "max_rank_one_vector_abs_difference": 2e-5,
        "integrated_intercept_conditional_chi2_abs_difference": 1e-4,
        "plugin_conditional_chi2_abs_difference": 1e-4,
        "integrated_predictive_logdet_abs_difference": 1e-4,
        "full_conditional_covariance_max_asymmetry_before_symmetrization": 1e-8,
    }
    reference_results: dict[str, Any] = {}
    reference_held = groups[REFERENCE_FOLD]
    reference_train = np.setdiff1d(np.arange(n), reference_held, assume_unique=True)
    for name, base in model_means.items():
        precision_result = solve_precision_block(p, y, base, reference_train, reference_held)
        reference_results[name] = covariance_reference(
            covariance, precision_result, y, base, reference_train, reference_held
        )
    reference_pass = all(
        checks[key] <= reference_tolerances[key]
        for checks in reference_results.values()
        for key in reference_tolerances
    )
    if not reference_pass:
        failure_record = {
            "status": "failed_covariance_domain_reconstruction_tolerance",
            "comparative_scoring_performed": False,
            "failure_stage": "pre-score_actual_matrix_covariance_reference",
            "row_count": n,
            "cid_order_sha256": cid_order_hash,
            "input_sha256": inputs["hashes"],
            "source_file_sha256": source_hashes,
            "fixed_BAO_shapes": shape_provenance,
            "reference_fold_index": REFERENCE_FOLD,
            "reference_tolerances": reference_tolerances,
            "reference_checks": reference_results,
            "P_times_C_identity_max_abs_error": inverse_residual,
            "covariance_asymmetry_before_symmetrization": covariance_asymmetry,
            "runtime_seconds": float(time.perf_counter() - started),
            "disposition": "stop without emitting comparative SN scores; diagnose algebra/conditioning",
        }
        RESULT_PATH.write_text(json.dumps(failure_record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(failure_record, indent=2, sort_keys=True))
        raise SystemExit("direct covariance-domain reconstruction failed before comparative scoring")

    fold_results: list[dict[str, Any]] = []
    for fold_index, held in enumerate(groups):
        train_mask = np.ones(n, dtype=bool)
        train_mask[held] = False
        train = np.flatnonzero(train_mask)
        fold: dict[str, Any] = {
            "fold": fold_index,
            "zHD_min": float(np.min(z_hd[held])),
            "zHD_max": float(np.max(z_hd[held])),
            "heldout_count": int(held.size),
            "training_count": int(train.size),
            "heldout_row_indices": held.tolist(),
            "models": {},
        }
        for name, base in model_means.items():
            result = solve_precision_block(p, y, base, train, held)
            result.pop("conditional_residual")
            result.pop("rank_one_vector")
            fold["models"][name] = result

        fold_results.append(fold)

    # Predeclared robustness: shift each of the three equal-count target cuts
    # independently by +/-1% of N.  No score is consulted when choosing cuts;
    # equal-z groups remain intact and every variant is recorded in full.
    edge_sensitivity: list[dict[str, Any]] = []
    shift_rows = EDGE_SHIFT_FRACTION * n
    for boundary in range(1, N_FOLDS):
        for direction, sign in (("left", -1.0), ("right", 1.0)):
            shifts = {boundary: sign * shift_rows}
            varied_edges, varied_groups = z_quantile_bins(
                z_hd, target_shifts_rows=shifts
            )
            varied_folds: list[dict[str, Any]] = []
            for fold_index, held in enumerate(varied_groups):
                train_mask = np.ones(n, dtype=bool)
                train_mask[held] = False
                train = np.flatnonzero(train_mask)
                fold: dict[str, Any] = {
                    "fold": fold_index,
                    "zHD_min": float(np.min(z_hd[held])),
                    "zHD_max": float(np.max(z_hd[held])),
                    "heldout_count": int(held.size),
                    "training_count": int(train.size),
                    "heldout_row_indices": held.tolist(),
                    "models": {},
                }
                for model_name, base in model_means.items():
                    score = solve_precision_block(p, y, base, train, held)
                    score.pop("conditional_residual")
                    score.pop("rank_one_vector")
                    fold["models"][model_name] = score
                varied_folds.append(fold)
            varied_aggregate = {
                model_name: {
                    "integrated_conditional_chi2_sum": float(
                        sum(
                            f["models"][model_name]["integrated_intercept_conditional_chi2"]
                            for f in varied_folds
                        )
                    ),
                    "integrated_conditional_neg2logpredictive_sum": float(
                        sum(
                            f["models"][model_name]["integrated_intercept_conditional_neg2logpredictive"]
                            for f in varied_folds
                        )
                    ),
                    "plugin_conditional_chi2_sum": float(
                        sum(f["models"][model_name]["plugin_conditional_chi2"] for f in varied_folds)
                    ),
                    "plugin_conditional_neg2logpredictive_sum": float(
                        sum(
                            f["models"][model_name]["plugin_conditional_neg2logpredictive"]
                            for f in varied_folds
                        )
                    ),
                }
                for model_name in model_means
            }
            edge_sensitivity.append(
                {
                    "variant": f"boundary_{boundary}_{direction}",
                    "target_shift_rows": shifts,
                    "redshift_edges_zHD": varied_edges,
                    "folds": varied_folds,
                    "aggregate": varied_aggregate,
                    "delta_IVS_minus_LCDM_integrated_conditional_neg2logpredictive": float(
                        varied_aggregate["BAO_interacting_vacuum"][
                            "integrated_conditional_neg2logpredictive_sum"
                        ]
                        - varied_aggregate["BAO_flat_LCDM"][
                            "integrated_conditional_neg2logpredictive_sum"
                        ]
                    ),
                }
            )

    aggregate: dict[str, Any] = {}
    for name in model_means:
        aggregate[name] = float(
            sum(f["models"][name]["integrated_intercept_conditional_chi2"] for f in fold_results)
        )
        aggregate[f"{name}_integrated_conditional_neg2logpredictive_sum"] = float(
            sum(
                f["models"][name]["integrated_intercept_conditional_neg2logpredictive"]
                for f in fold_results
            )
        )
    aggregate["delta_IVS_minus_LCDM_integrated"] = float(
        aggregate["BAO_interacting_vacuum"] - aggregate["BAO_flat_LCDM"]
    )
    aggregate["delta_IVS_minus_LCDM_integrated_conditional_neg2logpredictive"] = float(
        aggregate["BAO_interacting_vacuum_integrated_conditional_neg2logpredictive_sum"]
        - aggregate["BAO_flat_LCDM_integrated_conditional_neg2logpredictive_sum"]
    )
    aggregate["LCDM_plugin"] = float(
        sum(f["models"]["BAO_flat_LCDM"]["plugin_conditional_chi2"] for f in fold_results)
    )
    aggregate["IVS_plugin"] = float(
        sum(f["models"]["BAO_interacting_vacuum"]["plugin_conditional_chi2"] for f in fold_results)
    )
    aggregate["delta_IVS_minus_LCDM_plugin"] = float(
        aggregate["IVS_plugin"] - aggregate["LCDM_plugin"]
    )
    aggregate["LCDM_plugin_neg2logpredictive_sum"] = float(
        sum(f["models"]["BAO_flat_LCDM"]["plugin_conditional_neg2logpredictive"] for f in fold_results)
    )
    aggregate["IVS_plugin_neg2logpredictive_sum"] = float(
        sum(
            f["models"]["BAO_interacting_vacuum"]["plugin_conditional_neg2logpredictive"]
            for f in fold_results
        )
    )
    aggregate["delta_IVS_minus_LCDM_plugin_neg2logpredictive"] = float(
        aggregate["IVS_plugin_neg2logpredictive_sum"]
        - aggregate["LCDM_plugin_neg2logpredictive_sum"]
    )

    result: dict[str, Any] = {
        "status": (
            "complete_exploratory_cross_probe_prediction_screen"
            if reference_pass
            else "failed_covariance_domain_reconstruction_tolerance"
        ),
        "objective": "Conditionally score the frozen DESI BAO flat-LCDM and IVS shape points on separate DES-Dovekie SN data, fitting only the training SN magnitude intercept.",
        "scope": "Fixed BAO-selected cosmological shapes; no SN shape refit, no joint BAO+SN likelihood, no H0 inference, and no posterior/evidence claim.",
        "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --seconds 1200 -- .venv/bin/python experiments/dovekie_frozen_shape_cv/conditional_cv.py",
        "runtime_seconds": float(time.perf_counter() - started),
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "numerical_threads": 1,
        },
        "inputs": {
            "row_count": n,
            "source_order_preserved": True,
            "unique_CID_count": len({row["CID"] for row in inputs["records"]}),
            "cid_order_sha256": cid_order_hash,
            "mean_and_precision_sha256": inputs["hashes"],
            "source_file_sha256": source_hashes,
            "release_precision_shape": list(p.shape),
            "release_precision_is_float64_after_exact_float32_unpack": True,
            "precision_min_eigenvalue_from_pinned_record": 4.56550e-6,
            "precision_max_eigenvalue_from_pinned_record": 271.84070,
            "max_abs_P_times_C_minus_I": inverse_residual,
            "covariance_asymmetry_before_symmetrization": covariance_asymmetry,
            "covariance_cholesky_min_diagonal": float(np.min(np.diag(lc))),
        },
        "fixed_BAO_shapes": shape_provenance,
        "nuisance": {
            "magnitude_offset_prior": "flat over the common additive distance-modulus intercept; integrated analytically conditional on each training fold",
            "h0_gauge_km_s_mpc": sn.H0_GAUGE_KM_S_MPC,
            "h0_interpretation": "arbitrary scale absorbed by the intercept; H0 is not inferred",
            "distance_units": "c/H0 in Mpc; D_L=(1+zHEL)(c/H0) integral_0^zHD dz/E(z); modulus in magnitudes",
            "intercept_marginal_predictive_covariance": "S_(H|T) + sigma_M^2 v v^T, with sigma_M^2=(1_T^T Q_T 1_T)^-1, Q_T=C_TT^-1=P_TT-P_TH P_HH^-1 P_HT, and v=1_H+P_HH^-1 P_HT 1_T",
            "plug_in_sensitivity": "also report the fixed training-offset score using S_(H|T) only",
        },
        "fold_contract": {
            "scheme": "four contiguous zHD bins with nearest equal-row-count cutpoints over unique redshift groups; ties never split",
            "redshift_edges_zHD": edges,
            "folds": fold_results,
            "fold_sum_interpretation": "descriptive sum of overlapping conditional holdout scores, not one joint likelihood or independent four-sample statistic",
            "predeclared_edge_sensitivity": {
                "rule": "For each internal quartile boundary in turn, shift its target cumulative row count by +/-1% of N, hold the other two targets fixed at equal-count quartiles, then choose the nearest unique-z cumulative cut (ties lower index); do not split equal-z groups.",
                "target_shift_fraction_of_rows": EDGE_SHIFT_FRACTION,
                "variants": edge_sensitivity,
                "status": "complete_six_predeclared_neighboring_edge_perturbations",
                "delta_range_integrated_neg2logpredictive": [
                    float(
                        min(
                            v["delta_IVS_minus_LCDM_integrated_conditional_neg2logpredictive"]
                            for v in edge_sensitivity
                        )
                    ),
                    float(
                        max(
                            v["delta_IVS_minus_LCDM_integrated_conditional_neg2logpredictive"]
                            for v in edge_sensitivity
                        )
                    ),
                ],
            },
        },
        "aggregate_predictive_chi2": aggregate,
        "independent_covariance_reference": {
            "reference_fold_index": REFERENCE_FOLD,
            "method": "form C=P^-1 by Cholesky solves, then directly factor C_TT and S_(H|T); independently compute training intercept, conditional mean, rank-one offset uncertainty, log determinant and held-out scores in covariance space",
            "max_full_covariance_asymmetry_before_symmetrization": covariance_asymmetry,
            "max_precision_times_covariance_identity_residual": inverse_residual,
            "direct_covariance_cholesky_positive": True,
            "pass": reference_pass,
            "tolerances": reference_tolerances,
            "fold_model_checks": reference_results,
        },
        "numerical_checks": {
            "status": "pass_after_pre-score_covariance_reference",
            "checks": theory_checks,
            "tolerances": {
                "max_independent_distance_modulus_difference_mag": DISTANCE_CHECK_TOLERANCE_MAG,
                "max_precision_times_covariance_identity_abs_error": 1e-7,
                "max_inverse_covariance_asymmetry_before_symmetrization": 1e-10,
                "max_precision_vs_covariance_domain_reference_differences": reference_tolerances,
                "fixed_fb_0p16_matter_and_cdm_positive_through_z2p33": True,
            },
        },
        "uncertainty_and_limits": [
            "Only the common magnitude intercept is integrated; the BAO-selected Omega_m and interaction rate are fixed points, not SN posteriors.",
            "The flat intercept measure is common to both models; its arbitrary constant cancels in the within-fold comparison.",
            "The predictive log-determinant is identical between models within a fold; the reported comparison uses conditional quadratic deviance and records the determinant for audit.",
            "Four folds are overlapping conditional assessments on one SN release. Their sum is descriptive, not a joint score or calibrated significance.",
            "Dovekie shares a subset of historical low-redshift SNe with other SN compilations; no second SN compilation is used here.",
            "The two probes are reported separately; no combined likelihood or claim of exact cross-probe independence is made.",
            "The calculation is a late-time background prediction only; it does not validate perturbations, CMB, growth, or the authors' exact modified-CAMB run.",
        ],
    }
    RESULT_PATH.write_text(
        json.dumps(result, indent=2, sort_keys=True, default=json_numpy_to_native) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": result["status"],
        "runtime_seconds": result["runtime_seconds"],
        "redshift_edges_zHD": edges,
        "aggregate_predictive_chi2": aggregate,
        "reference_fold_checks": reference_results,
        "numerical_checks": theory_checks,
        "result_path": str(RESULT_PATH.relative_to(ROOT)),
    }, indent=2, sort_keys=True, default=json_numpy_to_native))
    if not reference_pass:
        raise SystemExit("direct covariance-domain reconstruction exceeded preregistered tolerance")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent CAMB × DESI DR2 BAO profile over a standard-BBN omega_b region.

The likelihood, released-vector parser, CAMB setup, and full-covariance scoring
are implemented locally here.  This reads compute30 only as an audit of setup
principles; it does not import an earlier scorer or modify shared files.
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

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"

TASK_DIR = Path(__file__).resolve().parent
PROJECT = TASK_DIR.parents[1]
sys.path.insert(0, str(PROJECT / "work/compute24/site"))

import camb  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402
from scipy.optimize import brentq, minimize  # noqa: E402

MEAN_PATH = PROJECT / "context/data/desi_dr2_mean.txt"
COV_PATH = PROJECT / "context/data/desi_dr2_cov.txt"
WHEEL_PATH = PROJECT / "context/code/camb-2.0.4-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl"
AUDIT_PATH = PROJECT / "work/compute30/boundary_gradient_audit.py"
RESULT_PATH = TASK_DIR / "bbn_bao_profile.json"
PLOT_PATH = TASK_DIR / "bbn_bao_profile.png"
REPORT_PATH = TASK_DIR / "REPORT.md"

C_KM_S = 299792.458
H0_BOUNDS = (10.0, 1000.0)
OBH2_BOUNDS = (0.0205, 0.0240)
OCH2_BOUNDS = (0.001, 0.99)
THETA100_CHECK = (0.5, 10.0)
PRIOR_DEFINITIONS = (
    {
        "id": "poulin_2026_arxiv_v1_BBN_LCDM",
        "mean": 0.022017,
        "sigma_minus": 0.000217,
        "sigma_plus": 0.000221,
        "provenance": "Poulin et al., arXiv:2607.20635v1, Appendix A / Table 2, preprint; BBN-only LambdaCDM result from updated D/H = 2.508 +/- 0.030, 100 omega_b = 2.2017 (+0.0221/-0.0217) at 68 percent; source details supplied by orchestrator's source review",
        "source_url": "https://arxiv.org/html/2607.20635",
        "penalty_kind": "piecewise Gaussianized 68-percent interval; this is not the paper's exact likelihood profile",
    },
    {
        "id": "PDG_2025_SBBN",
        "mean": 0.02205,
        "sigma_minus": 0.00043,
        "sigma_plus": 0.00043,
        "provenance": "PDG 2025 SBBN value 0.02205 +/- 0.00043 supplied by orchestrator's source review; kept as a distinct sensitivity prior",
        "source_url": None,
        "penalty_kind": "symmetric Gaussian sensitivity; not combined with the D/H-based preprint summary",
    },
    {
        "id": "provisional_conservative_BBN",
        "mean": 0.0222,
        "sigma_minus": 0.0005,
        "sigma_plus": 0.0005,
        "provenance": "provisional user-directed conservative-BBN sensitivity assumption; source review not completed",
        "source_url": None,
        "penalty_kind": "symmetric Gaussian sensitivity",
    },
)
INVALID_SCORE = 1.0e80
FIXED_CAMB = {
    "mnu_eV": 0.06,
    "nnu": 3.044,
    "num_massive_neutrinos": 1,
    "neutrino_hierarchy": "degenerate",
    "TCMB_K": 2.7255,
    "YHe": None,  # CAMB BBN-consistent helium prediction
    "w": -1.0,
    "wa": 0.0,
    "Omega_k": 0.0,
}
COARSE_GRID = np.unique(np.concatenate((
    np.linspace(OBH2_BOUNDS[0], OBH2_BOUNDS[1], 15, dtype=np.float64),
    np.asarray([prior["mean"] for prior in PRIOR_DEFINITIONS], dtype=np.float64),
)))
INVALID_COUNTS: dict[str, int] = {}
INVALID_EXAMPLES: list[dict[str, object]] = []


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_data() -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray]:
    rows: list[tuple[float, float, str]] = []
    for raw in MEAN_PATH.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            z_text, value_text, name = line.split()
            rows.append((float(z_text), float(value_text), name))
    z = np.asarray([row[0] for row in rows], dtype=np.float64)
    mean = np.asarray([row[1] for row in rows], dtype=np.float64)
    names = [row[2] for row in rows]
    covariance = np.asarray(np.loadtxt(COV_PATH, dtype=np.float64), dtype=np.float64)
    expected_order = [
        "DV_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs",
        "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs",
        "DH_over_rs", "DH_over_rs", "DM_over_rs",
    ]
    if z.shape != (13,) or mean.shape != (13,) or names != expected_order or covariance.shape != (13, 13):
        raise ValueError("DESI input must be the released 13-row vector and 13x13 covariance in its pinned order")
    if not np.isfinite(z).all() or not np.isfinite(mean).all() or not np.isfinite(covariance).all():
        raise ValueError("non-finite DESI input")
    if not np.allclose(covariance, covariance.T, rtol=0.0, atol=1.0e-12):
        raise ValueError("DESI covariance is not symmetric")
    np.linalg.cholesky(covariance)
    return z, mean, names, covariance


Z, Y, OBSERVABLES, COV = load_data()
CHOLESKY = np.linalg.cholesky(COV)


def make_parameters(h0: float, obh2: float, och2: float) -> camb.CAMBparams:
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=FIXED_CAMB["w"], wa=FIXED_CAMB["wa"], dark_energy_model="fluid")
    pars.set_cosmology(
        H0=float(h0), ombh2=float(obh2), omch2=float(och2), omk=FIXED_CAMB["Omega_k"],
        mnu=FIXED_CAMB["mnu_eV"], nnu=FIXED_CAMB["nnu"],
        num_massive_neutrinos=FIXED_CAMB["num_massive_neutrinos"],
        neutrino_hierarchy=FIXED_CAMB["neutrino_hierarchy"], TCMB=FIXED_CAMB["TCMB_K"],
        YHe=FIXED_CAMB["YHe"],
    )
    return pars


def evaluate(h0: float, obh2: float, och2: float) -> dict[str, object]:
    """Compute all 13 BAO predictions and independently score the full covariance."""
    pars = make_parameters(h0, obh2, och2)
    omega_m = float(pars.omegam)
    if not math.isfinite(omega_m) or not 0.0 < omega_m < 1.0:
        raise ValueError(f"flat positive-Lambda closure failed: Omega_m={omega_m}")
    background = camb.get_background(pars)
    omega_de = float(background.get_Omega("de", 0.0))
    if not math.isfinite(omega_de) or omega_de <= 0.0:
        raise ValueError(f"CAMB returned non-positive Omega_de={omega_de}")
    derived = background.get_derived_params()
    rdrag = float(derived["rdrag"])
    hz = np.asarray(background.hubble_parameter(Z), dtype=np.float64)
    dm = np.asarray(background.comoving_radial_distance(Z), dtype=np.float64)
    if not math.isfinite(rdrag) or rdrag <= 0.0 or not np.isfinite(hz).all() or not np.isfinite(dm).all():
        raise ValueError("CAMB returned invalid drag scale or background distances")
    if np.any(hz <= 0.0) or np.any(dm <= 0.0):
        raise ValueError("CAMB returned non-positive H(z) or D_M(z)")
    dh = C_KM_S / hz
    prediction = np.empty(13, dtype=np.float64)
    for index, name in enumerate(OBSERVABLES):
        if name == "DV_over_rs":
            prediction[index] = (Z[index] * dm[index] * dm[index] * dh[index]) ** (1.0 / 3.0) / rdrag
        elif name == "DM_over_rs":
            prediction[index] = dm[index] / rdrag
        elif name == "DH_over_rs":
            prediction[index] = dh[index] / rdrag
        else:
            raise ValueError(f"unsupported observable {name}")
    residual = prediction - Y
    whitened = np.linalg.solve(CHOLESKY, residual)
    chi2_cholesky = float(whitened @ whitened)
    chi2_dense = float(residual @ np.linalg.solve(COV, residual))
    difference = abs(chi2_dense - chi2_cholesky)
    if not math.isfinite(chi2_dense) or not math.isfinite(chi2_cholesky):
        raise ArithmeticError("non-finite full-covariance score")
    if difference > 2.0e-10 * max(1.0, abs(chi2_dense)):
        raise ArithmeticError(f"dense and Cholesky scores disagree by {difference}")
    return {
        "chi2_dense_solve": chi2_dense,
        "chi2_cholesky_whitening": chi2_cholesky,
        "dense_cholesky_abs_delta": difference,
        "predictions_in_released_order": prediction.tolist(),
        "residuals_in_released_order": residual.tolist(),
        "closure": {
            "Omega_m0_including_massive_neutrino": omega_m,
            "Omega_de0_CAMB": omega_de,
            "Omega_b0": float(pars.omegab),
            "Omega_c0": float(pars.omegac),
            "Omega_nu0": float(pars.omeganu),
            "r_drag_Mpc": rdrag,
            "100theta_star_CAMB": float(derived["thetastar"]),
        },
    }


def record_invalid(reason: str, point: list[float], detail: str | None = None) -> float:
    INVALID_COUNTS[reason] = INVALID_COUNTS.get(reason, 0) + 1
    if len(INVALID_EXAMPLES) < 40:
        sample: dict[str, object] = {"reason": reason, "parameters": point}
        if detail:
            sample["detail"] = detail[:240]
        INVALID_EXAMPLES.append(sample)
    return INVALID_SCORE


def from_unit2(unit: np.ndarray) -> tuple[float, float]:
    return (
        H0_BOUNDS[0] + float(unit[0]) * (H0_BOUNDS[1] - H0_BOUNDS[0]),
        OCH2_BOUNDS[0] + float(unit[1]) * (OCH2_BOUNDS[1] - OCH2_BOUNDS[0]),
    )


def to_unit2(h0: float, och2: float) -> np.ndarray:
    return np.clip(np.asarray([
        (h0 - H0_BOUNDS[0]) / (H0_BOUNDS[1] - H0_BOUNDS[0]),
        (och2 - OCH2_BOUNDS[0]) / (OCH2_BOUNDS[1] - OCH2_BOUNDS[0]),
    ], dtype=np.float64), 0.0, 1.0)


def profile_objective(unit: np.ndarray, fixed_ob: float, evaluation_counter: list[int]) -> float:
    evaluation_counter[0] += 1
    unit = np.asarray(unit, dtype=np.float64)
    h0, och2 = from_unit2(unit)
    if np.any(unit < 0.0) or np.any(unit > 1.0):
        return record_invalid("optimizer_outside_box", [h0, fixed_ob, och2])
    try:
        return float(evaluate(h0, fixed_ob, och2)["chi2_cholesky_whitening"])
    except Exception as exc:
        return record_invalid(type(exc).__name__, [h0, fixed_ob, och2], str(exc))


PROFILE_STARTS = (
    (67.0, 0.12),
    (82.0, 0.18),
    (133.6, 0.43),
    (300.0, 0.30),
)
POWELL_OPTIONS = {"xtol": 2.0e-9, "ftol": 2.0e-12, "maxiter": 500, "maxfev": 2200}


def optimize_profile(obh2: float, starts: list[tuple[float, float]], label: str) -> dict[str, object]:
    objective_evaluations = [0]
    trials: list[dict[str, object]] = []
    for start_h0, start_och2 in starts:
        result = minimize(
            profile_objective, to_unit2(start_h0, start_och2),
            args=(float(obh2), objective_evaluations), method="Powell",
            bounds=((0.0, 1.0), (0.0, 1.0)), options=POWELL_OPTIONS,
        )
        h0, och2 = from_unit2(result.x)
        try:
            final = evaluate(h0, obh2, och2)
            score = float(final["chi2_cholesky_whitening"])
            valid = (
                H0_BOUNDS[0] <= h0 <= H0_BOUNDS[1]
                and OCH2_BOUNDS[0] <= och2 <= OCH2_BOUNDS[1]
                and float(final["closure"]["Omega_de0_CAMB"]) > 0.0
            )
        except Exception as exc:
            final = None
            score = record_invalid(type(exc).__name__, [h0, obh2, och2], f"final point: {exc}")
            valid = False
        trials.append({
            "start": {"H0_km_s_Mpc": start_h0, "omega_c_h2": start_och2},
            "optimizer_success": bool(result.success),
            "optimizer_status": int(result.status),
            "optimizer_message": str(result.message),
            "function_evaluations": int(result.nfev),
            "H0_km_s_Mpc": h0,
            "omega_b_h2": float(obh2),
            "omega_c_h2": och2,
            "chi2": score,
            "valid_flat_positive_Lambda": bool(valid),
            "closure": None if final is None else final["closure"],
        })
    valid_trials = [trial for trial in trials if trial["valid_flat_positive_Lambda"]]
    if not valid_trials:
        raise RuntimeError(f"all {label} optimizations failed at omega_b h2={obh2}")
    selected = min(valid_trials, key=lambda trial: float(trial["chi2"]))
    spread = [float(trial["chi2"]) for trial in valid_trials]
    selected_eval = evaluate(
        float(selected["H0_km_s_Mpc"]), obh2, float(selected["omega_c_h2"]),
    )
    return {
        "omega_b_h2_fixed": float(obh2),
        "best_valid_trial": selected,
        "all_starts": trials,
        "valid_start_chi2_span": [min(spread), max(spread)],
        "selected_score_checks": selected_eval,
        "optimizer_function_evaluations": objective_evaluations[0],
        "all_reported_starts_valid": len(valid_trials) == len(trials),
    }


def h0_for_theta(theta100: float, obh2: float, och2: float) -> float:
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
    pars.set_cosmology(
        cosmomc_theta=float(theta100) / 100.0,
        ombh2=float(obh2), omch2=float(och2), omk=0.0,
        mnu=FIXED_CAMB["mnu_eV"], nnu=FIXED_CAMB["nnu"],
        num_massive_neutrinos=FIXED_CAMB["num_massive_neutrinos"],
        neutrino_hierarchy=FIXED_CAMB["neutrino_hierarchy"], TCMB=FIXED_CAMB["TCMB_K"],
        YHe=None, theta_H0_range=(1.0, 2000.0),
    )
    return float(pars.H0)


def map_h0_to_theta(h0_target: float, obh2: float, och2: float) -> float:
    """Invert CAMB's theta setter with a bracket, rather than identify theta_star."""
    brackets = ((0.6, 2.2, 33), (0.5, 10.0, 65))
    for lo, hi, ngrid in brackets:
        grid = np.linspace(lo, hi, ngrid, dtype=np.float64)
        mapped: list[tuple[float, float]] = []
        for theta in grid:
            try:
                mapped.append((float(theta), h0_for_theta(float(theta), obh2, och2)))
            except Exception:
                continue
        for (theta_a, h0_a), (theta_b, h0_b) in zip(mapped, mapped[1:]):
            if (h0_a - h0_target) * (h0_b - h0_target) <= 0.0:
                return float(brentq(
                    lambda theta: h0_for_theta(theta, obh2, och2) - h0_target,
                    theta_a, theta_b, xtol=2.0e-11, rtol=1.0e-12,
                ))
    raise ValueError(f"no CAMB theta inverse bracket for H0={h0_target}, omega_b h2={obh2}, omega_c h2={och2}")


def theta_check(fit: dict[str, object]) -> dict[str, object]:
    best = fit["best_valid_trial"]
    h0 = float(best["H0_km_s_Mpc"])
    obh2 = float(best["omega_b_h2"])
    och2 = float(best["omega_c_h2"])
    theta100 = map_h0_to_theta(h0, obh2, och2)
    recovered_h0 = h0_for_theta(theta100, obh2, och2)
    return {
        "implied_100theta_MC": theta100,
        "inside_declared_theta100_check_interval": THETA100_CHECK[0] <= theta100 <= THETA100_CHECK[1],
        "target_H0_km_s_Mpc": h0,
        "CAMB_theta_setter_recovered_H0_km_s_Mpc": recovered_h0,
        "absolute_H0_recovery_error_km_s_Mpc": abs(recovered_h0 - h0),
    }


def prior_delta_chi2(obh2: float, prior: dict[str, object]) -> float:
    mean = float(prior["mean"])
    sigma = float(prior["sigma_minus"] if obh2 < mean else prior["sigma_plus"])
    return ((obh2 - mean) / sigma) ** 2


def penalized_objective(unit: np.ndarray, counter: list[int], prior: dict[str, object]) -> float:
    counter[0] += 1
    unit = np.asarray(unit, dtype=np.float64)
    h0 = H0_BOUNDS[0] + float(unit[0]) * (H0_BOUNDS[1] - H0_BOUNDS[0])
    obh2 = OBH2_BOUNDS[0] + float(unit[1]) * (OBH2_BOUNDS[1] - OBH2_BOUNDS[0])
    och2 = OCH2_BOUNDS[0] + float(unit[2]) * (OCH2_BOUNDS[1] - OCH2_BOUNDS[0])
    if np.any(unit < 0.0) or np.any(unit > 1.0):
        return record_invalid("optimizer_outside_box", [h0, obh2, och2])
    try:
        data_score = float(evaluate(h0, obh2, och2)["chi2_cholesky_whitening"])
        prior_penalty = prior_delta_chi2(obh2, prior)
        return data_score + prior_penalty
    except Exception as exc:
        return record_invalid(type(exc).__name__, [h0, obh2, och2], str(exc))


def optimize_gaussian_prior(prior: dict[str, object], starts: list[tuple[float, float, float]]) -> dict[str, object]:
    counter = [0]
    trials: list[dict[str, object]] = []
    for start_h0, start_obh2, start_och2 in starts:
        unit0 = np.asarray([
            (start_h0 - H0_BOUNDS[0]) / (H0_BOUNDS[1] - H0_BOUNDS[0]),
            (start_obh2 - OBH2_BOUNDS[0]) / (OBH2_BOUNDS[1] - OBH2_BOUNDS[0]),
            (start_och2 - OCH2_BOUNDS[0]) / (OCH2_BOUNDS[1] - OCH2_BOUNDS[0]),
        ], dtype=np.float64)
        result = minimize(
            penalized_objective, np.clip(unit0, 0.0, 1.0), args=(counter, prior), method="Powell",
            bounds=((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)), options=POWELL_OPTIONS,
        )
        h0 = H0_BOUNDS[0] + float(result.x[0]) * (H0_BOUNDS[1] - H0_BOUNDS[0])
        obh2 = OBH2_BOUNDS[0] + float(result.x[1]) * (OBH2_BOUNDS[1] - OBH2_BOUNDS[0])
        och2 = OCH2_BOUNDS[0] + float(result.x[2]) * (OCH2_BOUNDS[1] - OCH2_BOUNDS[0])
        try:
            final = evaluate(h0, obh2, och2)
            chi2 = float(final["chi2_cholesky_whitening"])
            penalty = prior_delta_chi2(obh2, prior)
            total = chi2 + penalty
            valid = (
                H0_BOUNDS[0] <= h0 <= H0_BOUNDS[1]
                and OBH2_BOUNDS[0] <= obh2 <= OBH2_BOUNDS[1]
                and OCH2_BOUNDS[0] <= och2 <= OCH2_BOUNDS[1]
                and float(final["closure"]["Omega_de0_CAMB"]) > 0.0
            )
        except Exception as exc:
            final = None
            chi2 = record_invalid(type(exc).__name__, [h0, obh2, och2], f"prior final point: {exc}")
            penalty = prior_delta_chi2(obh2, prior)
            total = chi2 + penalty
            valid = False
        trials.append({
            "start": {"H0_km_s_Mpc": start_h0, "omega_b_h2": start_obh2, "omega_c_h2": start_och2},
            "optimizer_success": bool(result.success),
            "optimizer_status": int(result.status),
            "optimizer_message": str(result.message),
            "function_evaluations": int(result.nfev),
            "H0_km_s_Mpc": h0,
            "omega_b_h2": obh2,
            "omega_c_h2": och2,
            "data_chi2": chi2,
            "gaussian_prior_delta_chi2": penalty,
            "penalized_objective": total,
            "valid_flat_positive_Lambda": bool(valid),
            "closure": None if final is None else final["closure"],
        })
    valid_trials = [trial for trial in trials if trial["valid_flat_positive_Lambda"]]
    if not valid_trials:
        raise RuntimeError("all Gaussian-prior optimizations failed")
    selected = min(valid_trials, key=lambda trial: float(trial["penalized_objective"]))
    selected_eval = evaluate(
        float(selected["H0_km_s_Mpc"]), float(selected["omega_b_h2"]), float(selected["omega_c_h2"]),
    )
    return {
        "prior": prior,
        "prior_convention": "Piecewise/symmetric Gaussian -2 log prior penalty using the reported 68-percent uncertainty on each side; additive normalization omitted; truncated to declared BBN scan bounds",
        "all_starts": trials,
        "valid_start_penalized_objective_span": [
            min(float(t["penalized_objective"]) for t in valid_trials),
            max(float(t["penalized_objective"]) for t in valid_trials),
        ],
        "best_valid_trial": selected,
        "selected_score_checks": selected_eval,
        "optimizer_function_evaluations": counter[0],
    }


def main() -> None:
    started = time.perf_counter()
    if camb.__version__ != "2.0.4":
        raise RuntimeError(f"expected pinned CAMB 2.0.4, imported {camb.__version__}")

    profile: list[dict[str, object]] = []
    previous: tuple[float, float] | None = None
    for obh2 in COARSE_GRID:
        starts = list(PROFILE_STARTS)
        if previous is not None:
            starts.insert(0, previous)
        fit = optimize_profile(float(obh2), starts, "coarse profile")
        profile.append(fit)
        best = fit["best_valid_trial"]
        previous = (float(best["H0_km_s_Mpc"]), float(best["omega_c_h2"]))

    coarse_best = min(profile, key=lambda item: float(item["best_valid_trial"]["chi2"]))
    coarse_index = profile.index(coarse_best)
    coarse_b_values = [float(item["omega_b_h2_fixed"]) for item in profile]
    # Prior means are additional grid nodes, so use the nominal 15-point scan
    # spacing rather than the smaller irregular gaps they create.
    coarse_step = (OBH2_BOUNDS[1] - OBH2_BOUNDS[0]) / 14.0
    best_b = float(coarse_best["omega_b_h2_fixed"])
    refine_lo = max(OBH2_BOUNDS[0], best_b - coarse_step)
    refine_hi = min(OBH2_BOUNDS[1], best_b + coarse_step)
    refined_values = np.linspace(refine_lo, refine_hi, 9, dtype=np.float64)
    refined: list[dict[str, object]] = []
    nearby = profile[max(0, coarse_index - 1):min(len(profile), coarse_index + 2)]
    refine_starts = [
        (float(item["best_valid_trial"]["H0_km_s_Mpc"]), float(item["best_valid_trial"]["omega_c_h2"]))
        for item in nearby
    ]
    refine_starts.extend(PROFILE_STARTS[:2])
    unique_starts: list[tuple[float, float]] = []
    for start in refine_starts:
        if start not in unique_starts:
            unique_starts.append(start)
    for obh2 in refined_values:
        refined.append(optimize_profile(float(obh2), unique_starts, "refined profile"))

    prior_fits: list[dict[str, object]] = []
    prior_starts_by_id: dict[str, list[tuple[float, float, float]]] = {}
    for prior in PRIOR_DEFINITIONS:
        prior_mean = float(prior["mean"])
        prior_seed_profiles = sorted(profile + refined, key=lambda item: abs(float(item["omega_b_h2_fixed"]) - prior_mean))
        prior_start_from_profile = prior_seed_profiles[0]["best_valid_trial"]
        prior_starts = [
            (
                float(prior_start_from_profile["H0_km_s_Mpc"]),
                float(prior_start_from_profile["omega_b_h2"]),
                float(prior_start_from_profile["omega_c_h2"]),
            ),
            (67.0, prior_mean, 0.12),
            (82.0, prior_mean, 0.18),
            (133.6, prior_mean, 0.43),
        ]
        prior_starts_by_id[str(prior["id"])] = prior_starts
        prior_fit = optimize_gaussian_prior(prior, prior_starts)
        prior_fit["theta_support_check"] = theta_check(prior_fit)
        prior_fits.append(prior_fit)

    for fit in profile + refined:
        fit["theta_support_check"] = theta_check(fit)

    dense_all = sorted(profile + refined, key=lambda item: float(item["omega_b_h2_fixed"]))
    dense_by_b: dict[float, dict[str, object]] = {}
    for item in dense_all:
        dense_by_b[float(item["omega_b_h2_fixed"])] = item
    dense_all = list(dense_by_b.values())
    priorized_grids: dict[str, list[dict[str, object]]] = {}
    grid_prior_minima: dict[str, dict[str, object]] = {}
    for prior in PRIOR_DEFINITIONS:
        prior_rows = [
            {
                "omega_b_h2": float(item["omega_b_h2_fixed"]),
                "data_chi2": float(item["best_valid_trial"]["chi2"]),
                "gaussian_prior_delta_chi2": prior_delta_chi2(float(item["omega_b_h2_fixed"]), prior),
                "total_delta_chi2_raw": float(item["best_valid_trial"]["chi2"])
                    + prior_delta_chi2(float(item["omega_b_h2_fixed"]), prior),
            }
            for item in dense_all
        ]
        priorized_grids[str(prior["id"])] = prior_rows
        grid_prior_minima[str(prior["id"])] = min(prior_rows, key=lambda item: float(item["total_delta_chi2_raw"]))
    best_data_trial = min(dense_all, key=lambda item: float(item["best_valid_trial"]["chi2"]))["best_valid_trial"]
    prior_fit_summary = {str(fit["prior"]["id"]): fit for fit in prior_fits}
    refine_best = min(refined, key=lambda item: float(item["best_valid_trial"]["chi2"]))
    refine_vs_coarse = {
        "coarse_best_omega_b_h2": best_b,
        "coarse_best_data_chi2": float(coarse_best["best_valid_trial"]["chi2"]),
        "refined_best_omega_b_h2": float(refine_best["omega_b_h2_fixed"]),
        "refined_best_data_chi2": float(refine_best["best_valid_trial"]["chi2"]),
        "refined_minus_coarse_chi2": float(refine_best["best_valid_trial"]["chi2"]) - float(coarse_best["best_valid_trial"]["chi2"]),
        "grid_definition": "15-point inclusive linear scan plus each distinct prior mean; local 9-point refinement spanning one nominal scan spacing on each available side of the coarse best point",
    }

    distinct_local_solutions: list[dict[str, object]] = []
    profile_start_count = 0
    best_basin_start_count = 0
    for item in profile + refined:
        slice_best = float(item["best_valid_trial"]["chi2"])
        for trial in item["all_starts"]:
            profile_start_count += 1
            if not trial["valid_flat_positive_Lambda"]:
                continue
            delta = float(trial["chi2"]) - slice_best
            if delta <= 1.0e-7:
                best_basin_start_count += 1
            elif delta > 1.0:
                distinct_local_solutions.append({
                    "omega_b_h2_fixed": float(item["omega_b_h2_fixed"]),
                    "start": trial["start"],
                    "returned_H0_km_s_Mpc": float(trial["H0_km_s_Mpc"]),
                    "returned_omega_c_h2": float(trial["omega_c_h2"]),
                    "chi2": float(trial["chi2"]),
                    "delta_chi2_vs_best_start_on_same_slice": delta,
                    "optimizer_message": trial["optimizer_message"],
                })
    local_minima_audit = {
        "profile_start_total": profile_start_count,
        "starts_reaching_slice_best_within_delta_chi2_1e-7": best_basin_start_count,
        "valid_distinct_suboptimal_solutions_delta_chi2_over_1_count": len(distinct_local_solutions),
        "valid_distinct_suboptimal_solutions": distinct_local_solutions,
        "summary": "Several broad high-H0/high-omega_c starts terminate at the omega_c h2 upper bound with valid positive-Lambda closure and very poor chi2; they are retained as alternate bounded local solutions. The low-chi2 basin is recovered independently by standard-scale starts. This documents multimodality and does not establish a mathematical global-minimum proof.",
    }

    # Plot the unconstrained BAO profile and each distinct prior sensitivity.
    ordered_b = np.asarray([float(item["omega_b_h2_fixed"]) for item in dense_all], dtype=np.float64)
    ordered_chi2 = np.asarray([float(item["best_valid_trial"]["chi2"]) for item in dense_all], dtype=np.float64)
    delta_data = ordered_chi2 - float(np.min(ordered_chi2))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
    axes[0].plot(ordered_b, delta_data, "o-", color="#205493", ms=4)
    axes[0].axvline(0.0222, color="#b34b4b", ls="--", lw=1.1, label="provisional-prior mean")
    axes[0].set(xlabel=r"fixed $\omega_b h^2$", ylabel=r"$\Delta\chi^2_{\rm BAO}$", title="Full-covariance BAO profile")
    axes[0].legend(frameon=False, fontsize=8)
    prior_colors = {
        "poulin_2026_arxiv_v1_BBN_LCDM": "#7d3c98",
        "PDG_2025_SBBN": "#d17c00",
        "provisional_conservative_BBN": "#457a32",
    }
    for prior in PRIOR_DEFINITIONS:
        prior_id = str(prior["id"])
        curve = ordered_chi2 + np.asarray([prior_delta_chi2(float(b), prior) for b in ordered_b], dtype=np.float64)
        curve -= float(np.min(curve))
        axes[1].plot(ordered_b, curve, "o-", color=prior_colors[prior_id], ms=3, lw=1.25, label=prior_id.replace("_", " "))
        axes[1].axvline(float(prior["mean"]), color=prior_colors[prior_id], ls="--", lw=0.9, alpha=0.85)
    axes[1].set(xlabel=r"fixed $\omega_b h^2$", ylabel=r"relative penalized $\Delta\chi^2$", title="Prior sensitivity (distinct inputs)")
    axes[1].legend(frameon=False, fontsize=7)
    for ax in axes:
        ax.grid(alpha=0.22)
    fig.suptitle("Flat-ΛCDM background fit to the released DESI DR2 BAO vector")
    fig.savefig(PLOT_PATH, dpi=180)
    plt.close(fig)

    runtime = time.perf_counter() - started
    try:
        available_cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        available_cpus = os.cpu_count()
    best_data_fit = min(dense_all, key=lambda item: float(item["best_valid_trial"]["chi2"]))
    result: dict[str, object] = {
        "status": "completed_standard_BBN_region_BAO_profile_and_three_distinct_prior_sensitivities",
        "finding_scope": "exploratory BAO-only flat-Lambda profile; not a posterior, evidence calculation, or DESI collaboration parameter reproduction",
        "inputs": {
            "row_count": int(Z.size),
            "observable_order": OBSERVABLES,
            "redshifts": Z.tolist(),
            "covariance_shape": list(COV.shape),
            "covariance_symmetric_abs_max": float(np.max(np.abs(COV - COV.T))),
            "covariance_positive_definite": True,
            "full_13_by_13_covariance_used": True,
            "mean_path": str(MEAN_PATH.relative_to(PROJECT)),
            "covariance_path": str(COV_PATH.relative_to(PROJECT)),
        },
        "model": {
            "description": "CAMB flat fluid Lambda background; w=-1, wa=0, Omega_k=0; YHe=None retains CAMB BBN-consistency relation",
            "c_km_s": C_KM_S,
            "camb_setup": FIXED_CAMB,
            "likelihood": "Gaussian released DESI DR2 compressed BAO vector; score residual^T covariance^-1 residual",
            "observables": {
                "DV_over_rd": "(z * DM^2 * DH)^(1/3) / rdrag",
                "DM_over_rd": "DM / rdrag",
                "DH_over_rd": "(c/H) / rdrag",
            },
        },
        "search_domain_and_prior": {
            "H0_km_s_Mpc_bounds": list(H0_BOUNDS),
            "omega_b_h2_fixed_profile_scan_bounds": list(OBH2_BOUNDS),
            "omega_b_grid_coarse": [float(x) for x in COARSE_GRID],
            "omega_c_h2_bounds": list(OCH2_BOUNDS),
            "100theta_MC_support_check_interval": list(THETA100_CHECK),
            "H0_is_a_profile_coordinate_not_a_gaussian_prior": True,
            "baryon_grid_rationale": "user-directed broad standard-BBN-motivated scan interval; provenance intentionally provisional pending source review",
            "distinct_BBN_prior_sensitivities": [
                {
                    **prior,
                    "exact_form": "piecewise/symmetric Gaussianized -2 log prior penalty based on reported 68-percent interval; additive normalization omitted; truncated to [0.0205,0.0240]",
                }
                for prior in PRIOR_DEFINITIONS
            ],
        },
        "profile_grid": profile,
        "profile_refinement": refined,
        "profile_convergence": refine_vs_coarse,
        "local_minima_audit": local_minima_audit,
        "profile_minimum_boundary_status": {
            "best_omega_b_h2_is_upper_scan_edge": abs(float(best_data_fit["omega_b_h2_fixed"]) - OBH2_BOUNDS[1]) <= 1.0e-14,
            "upper_scan_edge_omega_b_h2": OBH2_BOUNDS[1],
            "interior_minimum_established": False,
            "reason": "lowest sampled pure-BAO profile score occurs at the upper edge of the declared BBN-motivated scan interval",
        },
        "pure_BAO_best_valid_fit": {
            "profile_record": best_data_fit,
            "theta_support_check": best_data_fit["theta_support_check"],
        },
        "BBN_prior_sensitivity_fits": prior_fits,
        "grid_prior_augmented_minima": grid_prior_minima,
        "data_vs_prior_sensitivity": {
            "pure_BAO_minimum_chi2": float(best_data_fit["best_valid_trial"]["chi2"]),
            "pure_BAO_best_omega_b_h2": float(best_data_fit["omega_b_h2_fixed"]),
            "prior_fit_summaries": {
                prior_id: {
                    "best_omega_b_h2": float(fit["best_valid_trial"]["omega_b_h2"]),
                    "data_chi2": float(fit["best_valid_trial"]["data_chi2"]),
                    "prior_delta_chi2": float(fit["best_valid_trial"]["gaussian_prior_delta_chi2"]),
                    "penalized_objective": float(fit["best_valid_trial"]["penalized_objective"]),
                    "theta_support_check": fit["theta_support_check"],
                }
                for prior_id, fit in prior_fit_summary.items()
            },
            "grid_profile_penalized_minima": grid_prior_minima,
        },
        "quality_checks": {
            "every_selected_profile_score_compares_dense_and_Cholesky": True,
            "selected_score_max_dense_cholesky_abs_delta": float(max(
                [item["selected_score_checks"]["dense_cholesky_abs_delta"] for item in profile + refined]
                + [fit["selected_score_checks"]["dense_cholesky_abs_delta"] for fit in prior_fits]
            )),
            "all_profile_theta_checks_inside_interval": all(
                item["theta_support_check"]["inside_declared_theta100_check_interval"] for item in profile + refined
            ),
            "all_prior_fit_theta_inside_interval": all(
                bool(fit["theta_support_check"]["inside_declared_theta100_check_interval"]) for fit in prior_fits
            ),
            "every_selected_profile_has_positive_flat_Lambda_closure": all(
                item["best_valid_trial"]["valid_flat_positive_Lambda"] for item in profile + refined
            ),
            "all_prior_fits_have_positive_flat_Lambda_closure": all(
                bool(fit["best_valid_trial"]["valid_flat_positive_Lambda"]) for fit in prior_fits
            ),
            "within_start_chi2_spans": {
                "coarse_profile_max": float(max(
                    item["valid_start_chi2_span"][1] - item["valid_start_chi2_span"][0] for item in profile
                )),
                "refined_profile_max": float(max(
                    item["valid_start_chi2_span"][1] - item["valid_start_chi2_span"][0] for item in refined
                )),
                "prior_penalized_objective_span_by_id": {
                    str(fit["prior"]["id"]): float(
                        fit["valid_start_penalized_objective_span"][1]
                        - fit["valid_start_penalized_objective_span"][0]
                    ) for fit in prior_fits
                },
            },
            "optimizer_non_success_count": int(sum(
                not bool(trial["optimizer_success"])
                for item in profile + refined for trial in item["all_starts"]
            ) + sum(
                not bool(trial["optimizer_success"]) for fit in prior_fits for trial in fit["all_starts"]
            )),
            "invalid_objective_evaluations_by_reason": INVALID_COUNTS,
            "invalid_objective_evaluation_examples_first_40": INVALID_EXAMPLES,
        },
        "execution": {
            "runtime_seconds": runtime,
            "numerical_workers": 1,
            "numerical_library_threads": {key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
            "available_affinity_cpus": available_cpus,
            "gpu_used": False,
            "optimization": "bounded Powell local optimization on scaled unit-square coordinates; four or five deterministic starts per fixed omega_b slice; four starts for 3D prior fit",
            "powell_options": POWELL_OPTIONS,
            "profile_starts_H0_omega_c": [list(x) for x in PROFILE_STARTS],
            "prior_fit_starts_H0_omega_b_omega_c_by_prior": {
                prior_id: [list(x) for x in starts] for prior_id, starts in prior_starts_by_id.items()
            },
            "profile_optimizer_function_evaluations": int(sum(
                item["optimizer_function_evaluations"] for item in profile + refined
            )),
            "prior_optimizer_function_evaluations_by_id": {
                str(fit["prior"]["id"]): int(fit["optimizer_function_evaluations"]) for fit in prior_fits
            },
            "reproduction_command": "cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute31/no_campaign_deadline.json --seconds 850 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute31/bbn_bao_profile.py",
        },
        "versions": {
            "CAMB": camb.__version__,
            "NumPy": np.__version__,
            "SciPy": scipy.__version__,
            "Matplotlib": matplotlib.__version__,
            "Python": sys.version,
            "platform": platform.platform(),
        },
        "hashes": {
            "script_sha256": sha256(Path(__file__)),
            "mean_sha256": sha256(MEAN_PATH),
            "covariance_sha256": sha256(COV_PATH),
            "CAMB_2_0_4_wheel_sha256": sha256(WHEEL_PATH),
            "compute30_audit_script_sha256": sha256(AUDIT_PATH),
        },
        "failed_attempts": {
            "invalid_objective_evaluations": INVALID_COUNTS,
            "prior_run_attempts": [
                {
                    "attempt": "001",
                    "status": "numerical profile completed; report-format QA failed and was corrected before final rerun",
                    "runtime_seconds": 136.5834279647097,
                    "script_sha256": "5f16d3c635f79c6d56e6560f0ddadb3da2a21492838b0a04bc2b21f85567d39b",
                    "result_json_sha256": "b8761d5d6e7cb92b87041c8fcd13981ae0a53587a70f3411ff344cde348d39ed",
                    "failure": "Python interpreted LaTeX backslash sequences in the generated f-string report, producing SyntaxWarnings and malformed theta/math markup; calculation values were retained and the report template was changed to a raw f-string.",
                    "command": "cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute31/no_campaign_deadline.json --seconds 850 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute31/bbn_bao_profile.py",
                }
            ],
            "optimizer_trials_reporting_non_success": [
                trial for item in profile + refined for trial in item["all_starts"] if not trial["optimizer_success"]
            ] + [
                trial for fit in prior_fits for trial in fit["all_starts"] if not trial["optimizer_success"]
            ],
            "note": "Invalid physical/CAMB probes are penalized and retained by reason/example; a Powell non-success result is retained even if its final point is valid.",
        },
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result_hash = sha256(RESULT_PATH)
    selected_chi2 = float(best_data_fit["best_valid_trial"]["chi2"])
    selected_b = float(best_data_fit["omega_b_h2_fixed"])
    source_prior_result = prior_fit_summary["poulin_2026_arxiv_v1_BBN_LCDM"]["best_valid_trial"]
    pdg_prior_result = prior_fit_summary["PDG_2025_SBBN"]["best_valid_trial"]
    provisional_prior_result = prior_fit_summary["provisional_conservative_BBN"]["best_valid_trial"]
    max_span = float(result["quality_checks"]["within_start_chi2_spans"]["coarse_profile_max"])
    report = rf"""# Compute31: standard-BBN-region DESI DR2 BAO profile

Status: **exploratory, independently implemented CAMB × DESI DR2 BAO profile**. The data-only minimum is pinned to the upper scan edge \(\omega_bh^2=0.0240\), so the scan does not establish an interior baryon-density minimum. Three separate BBN prior sensitivities are reported. The user-directed conservative prior remains provisional; the source-matched 2026 input is a preprint translated to a local piecewise Gaussianized penalty rather than importing its likelihood. This is not a posterior, evidence calculation, or DESI collaboration parameter reproduction.

The released 13-row BAO mean vector and full 13×13 covariance were parsed and scored in their pinned order. The prediction code independently evaluates (D_V/r_d), (D_M/r_d), and (D_H/r_d=(c/H)/r_d) with local CAMB 2.0.4. The fixed setup is flat fluid ΛCDM, (m_\nu=0.06\,\mathrm{{eV}}), (N_{{eff}}=3.044), one degenerate massive species, (T_{{CMB}}=2.7255\,\mathrm{{K}}), and CAMB's BBN-consistent helium setting (`YHe=None`).

## Fit and checks

For each fixed \(\omega_bh^2\) in [{OBH2_BOUNDS[0]:.4f}, {OBH2_BOUNDS[1]:.4f}], the code reoptimized \(H_0\) and \(\omega_ch^2\) in \([{H0_BOUNDS[0]:g},{H0_BOUNDS[1]:g}]\,\mathrm{{km\,s^{{-1}}\,Mpc^{{-1}}}}\) and [{OCH2_BOUNDS[0]:.3f}, {OCH2_BOUNDS[1]:.2f}] using four or five deterministic Powell starts. The coarse profile has {len(profile)} points; the local refinement has {len(refined)} points. At the lowest sampled data-only score, \(\chi^2={selected_chi2:.10f}\), \(\omega_bh^2={selected_b:.7f}\), \(H_0={float(best_data_fit['best_valid_trial']['H0_km_s_Mpc']):.6f}\), and \(\omega_ch^2={float(best_data_fit['best_valid_trial']['omega_c_h2']):.7f}\). The refined-versus-coarse profile change is \(\Delta\chi^2={refine_vs_coarse['refined_minus_coarse_chi2']:.4g}\) at the respective minima; it is a resolution check, not proof of globality.

The coarse and refined sampled minima agree at the upper boundary with \(\Delta\chi^2={refine_vs_coarse['refined_minus_coarse_chi2']:.3g}\); this is a local resolution check. The multistart audit found {local_minima_audit['valid_distinct_suboptimal_solutions_delta_chi2_over_1_count']} valid but much worse bounded solutions, generally ending at the \(\omega_ch^2=0.99\) upper bound from broad high-\(H_0\) starts. These alternatives remain in the JSON; only {local_minima_audit['starts_reaching_slice_best_within_delta_chi2_1e-7']} of {local_minima_audit['profile_start_total']} fixed-slice starts landed within \(10^{{-7}}\) of that slice's best start. The lowest-scoring basin was independently recovered from standard-scale starts, but globality over the full continuous box is not proved. Powell success flags and every start are retained, including non-success outcomes. Selected fixed-slice dense-precision and Cholesky-whitening scores agree to at most {float(result['quality_checks']['selected_score_max_dense_cholesky_abs_delta']):.3g}. Every selected profile point was re-mapped through CAMB's `cosmomc_theta` setter and checked against broad \(100\theta_{{MC}}\in[0.5,10]\) support; all selected minima have positive flat-Λ closure.

Three distinct prior sensitivities were evaluated as alternatives and were never multiplied together. Poulin et al., [arXiv:2607.20635v1](https://arxiv.org/html/2607.20635), Appendix A/Table 2, report a BBN-only LambdaCDM result based on updated D/H = 2.508 +/- 0.030, 100 omega_b = 2.2017 (+0.0221/-0.0217) at 68% (preprint). I use a piecewise Gaussianized penalty with omega_b h2 = 0.022017 (+0.000221/-0.000217); this is not the source paper's exact profile. The separate PDG 2025 SBBN sensitivity supplied by the orchestrator is 0.02205 +/- 0.00043. The provisional conservative sensitivity is 0.0222 +/- 0.0005; its source provenance remains unreviewed. The reoptimized omega_b h2 values are {float(source_prior_result['omega_b_h2']):.7f} (preprint), {float(pdg_prior_result['omega_b_h2']):.7f} (PDG 2025), and {float(provisional_prior_result['omega_b_h2']):.7f} (provisional), respectively. The PDG entry's page/URL was not present in this worker's local inputs and remains to be attached by the orchestrator.

## Reproduction

The isolated script ran in {runtime:.2f} seconds with one numerical worker, numerical-library threads set to one, and no GPU. No installs or downloads were used. The 850-second command cap is below the task's 900-second limit:

```sh
cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute31/no_campaign_deadline.json --seconds 850 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute31/bbn_bao_profile.py
```

The machine-readable record contains exact bounds and prior form, all profile slices and optimizer starts, invalid probes and non-success outcomes, software versions, run time, input and code hashes, and the score/closure/theta checks. Result SHA-256: `{result_hash}`. Software: CAMB {camb.__version__}, NumPy {np.__version__}, SciPy {scipy.__version__}, Python {platform.python_version()}.
"""
    REPORT_PATH.write_text(report, encoding="utf-8")
    print(json.dumps({
        "result": str(RESULT_PATH),
        "report": str(REPORT_PATH),
        "plot": str(PLOT_PATH),
        "runtime_seconds": runtime,
        "result_sha256": result_hash,
        "profile_points": len(profile),
        "refinement_points": len(refined),
        "prior_fit_omega_b_h2": {
            "poulin_2026_arxiv_v1_BBN_LCDM": float(source_prior_result["omega_b_h2"]),
            "PDG_2025_SBBN": float(pdg_prior_result["omega_b_h2"]),
            "provisional_conservative_BBN": float(provisional_prior_result["omega_b_h2"]),
        },
    }, indent=2))


if __name__ == "__main__":
    main()

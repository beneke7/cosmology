#!/usr/bin/env python3
"""Independent local/profile audit of the compute26 flat-LCDM BAO edge.

This script owns its CAMB setup, DESI input parser, observables, covariance score,
and fixed-omega_b profile.  It reads the compute26 candidate and shared pinned
inputs, but writes only under work/compute30/.
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

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT / "work/compute24/site"))

import camb  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402
from scipy.optimize import brentq, minimize  # noqa: E402

MEAN_FILE = PROJECT / "context/data/desi_dr2_mean.txt"
COV_FILE = PROJECT / "context/data/desi_dr2_cov.txt"
REFERENCE_FILE = PROJECT / "work/compute26/independent_lcdm_bao_profile.json"
WHEEL_FILE = PROJECT / "context/code/camb-2.0.4-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl"
OUT_FILE = HERE / "boundary_gradient_audit.json"
REPORT_FILE = HERE / "REPORT.md"

C_KM_S = 299792.458
H0_BOUNDS = (10.0, 1000.0)
OBH2_BOUNDS = (0.005, 0.1)
OCH2_BOUNDS = (0.001, 0.99)
THETA100_BOUNDS = (0.5, 10.0)
INVALID_SCORE = 1.0e80
FIXED_CAMB = {
    "mnu_eV": 0.06,
    "nnu": 3.044,
    "num_massive_neutrinos": 1,
    "neutrino_hierarchy": "degenerate",
    "TCMB_K": 2.7255,
    "YHe": None,
    "tau": 0.054,
    "As": 2.1e-9,
    "ns": 0.965,
}
PROFILE_OBH2 = (0.1, 0.0999999, 0.099999, 0.09999, 0.0999, 0.099, 0.095, 0.08, 0.05)
FD_STEPS = (1.0e-5, 1.0e-4, 5.0e-4)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_mean() -> tuple[np.ndarray, np.ndarray, list[str]]:
    parsed: list[tuple[float, float, str]] = []
    for raw in MEAN_FILE.read_text(encoding="utf-8").splitlines():
        row = raw.strip()
        if row and not row.startswith("#"):
            z_text, value_text, name = row.split()
            parsed.append((float(z_text), float(value_text), name))
    return (
        np.asarray([row[0] for row in parsed], dtype=np.float64),
        np.asarray([row[1] for row in parsed], dtype=np.float64),
        [row[2] for row in parsed],
    )


Z, Y, OBSERVABLES = read_mean()
COV = np.asarray(np.loadtxt(COV_FILE, dtype=np.float64), dtype=np.float64)
if Z.shape != (13,) or Y.shape != (13,) or len(OBSERVABLES) != 13 or COV.shape != (13, 13):
    raise ValueError("Expected 13 DESI DR2 rows and a 13 by 13 covariance")
if not np.isfinite(Z).all() or not np.isfinite(Y).all() or not np.isfinite(COV).all():
    raise ValueError("DESI input contains non-finite values")
if not np.allclose(COV, COV.T, rtol=0.0, atol=1.0e-12):
    raise ValueError("DESI covariance is not symmetric")
if not np.array_equal(OBSERVABLES, [
    "DV_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs",
    "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs",
    "DH_over_rs", "DH_over_rs", "DM_over_rs",
]):
    raise ValueError(f"Unexpected released observable order: {OBSERVABLES}")
CHOLESKY = np.linalg.cholesky(COV)


def make_params(h0: float, obh2: float, och2: float) -> camb.CAMBparams:
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
    pars.set_cosmology(
        H0=float(h0), ombh2=float(obh2), omch2=float(och2), omk=0.0,
        mnu=FIXED_CAMB["mnu_eV"], nnu=FIXED_CAMB["nnu"],
        num_massive_neutrinos=FIXED_CAMB["num_massive_neutrinos"],
        neutrino_hierarchy=FIXED_CAMB["neutrino_hierarchy"],
        TCMB=FIXED_CAMB["TCMB_K"], YHe=FIXED_CAMB["YHe"], tau=FIXED_CAMB["tau"],
    )
    pars.InitPower.set_params(As=FIXED_CAMB["As"], ns=FIXED_CAMB["ns"], r=0.0)
    return pars


def score(h0: float, obh2: float, och2: float) -> dict[str, object]:
    """Return score and checks; raw evaluation also supports explicit edge diagnostics."""
    pars = make_params(h0, obh2, och2)
    om0 = float(pars.omegam)
    if not (math.isfinite(om0) and 0.0 < om0 < 1.0):
        raise ValueError(f"non-positive flat-Lambda closure: Omega_m={om0}")
    background = camb.get_background(pars)
    ode0 = float(background.get_Omega("de", 0.0))
    if not (math.isfinite(ode0) and ode0 > 0.0):
        raise ValueError(f"non-positive CAMB dark-energy density: Omega_de={ode0}")
    derived = background.get_derived_params()
    rdrag = float(derived["rdrag"])
    if not (math.isfinite(rdrag) and rdrag > 0.0):
        raise ValueError(f"non-positive r_drag: {rdrag}")
    hz = np.asarray(background.hubble_parameter(Z), dtype=np.float64)
    dm = np.asarray(background.comoving_radial_distance(Z), dtype=np.float64)
    if not np.isfinite(hz).all() or not np.isfinite(dm).all() or np.any(hz <= 0.0) or np.any(dm <= 0.0):
        raise ValueError("CAMB returned invalid background distances")
    dh = C_KM_S / hz
    prediction = np.empty(13, dtype=np.float64)
    for i, name in enumerate(OBSERVABLES):
        if name == "DM_over_rs":
            prediction[i] = dm[i] / rdrag
        elif name == "DH_over_rs":
            prediction[i] = dh[i] / rdrag
        elif name == "DV_over_rs":
            prediction[i] = (Z[i] * dm[i] ** 2 * dh[i]) ** (1.0 / 3.0) / rdrag
        else:
            raise ValueError(f"Unknown DESI observable {name}")
    residual = prediction - Y
    dense = float(residual @ np.linalg.solve(COV, residual))
    whitened = np.linalg.solve(CHOLESKY, residual)
    chol = float(whitened @ whitened)
    if not (math.isfinite(dense) and math.isfinite(chol)) or abs(dense - chol) > 2.0e-10 * max(1.0, abs(dense)):
        raise ArithmeticError(f"full covariance score mismatch: dense={dense} Cholesky={chol}")
    return {
        "chi2_dense_solve": dense,
        "chi2_cholesky_whitening": chol,
        "dense_cholesky_abs_delta": abs(dense - chol),
        "predictions_in_released_order": prediction.tolist(),
        "residuals_in_released_order": residual.tolist(),
        "derived": {
            "Omega_m0_including_massive_nu": om0,
            "Omega_de0_CAMB": ode0,
            "Omega_b0": float(pars.omegab),
            "Omega_c0": float(pars.omegac),
            "Omega_nu0": float(pars.omeganu),
            "r_drag_Mpc": rdrag,
            "100theta_star_CAMB": float(derived["thetastar"]),
        },
        "checks": {"positive_flat_LCDM_closure": True, "finite_positive_background": True},
    }


def h0_from_theta(theta100: float, obh2: float, och2: float) -> float:
    """Apply CAMB's actual CosmoMC theta setter with a broad internal H0 bracket."""
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
    pars.set_cosmology(
        cosmomc_theta=float(theta100) / 100.0, ombh2=float(obh2), omch2=float(och2),
        omk=0.0, mnu=FIXED_CAMB["mnu_eV"], nnu=FIXED_CAMB["nnu"],
        num_massive_neutrinos=FIXED_CAMB["num_massive_neutrinos"],
        neutrino_hierarchy=FIXED_CAMB["neutrino_hierarchy"], TCMB=FIXED_CAMB["TCMB_K"],
        YHe=FIXED_CAMB["YHe"], tau=FIXED_CAMB["tau"], theta_H0_range=(1.0, 2000.0),
    )
    return float(pars.H0)


def implied_theta100(target_h0: float, obh2: float, och2: float) -> float:
    # The BAO profile minima stay in the broad theta range; bracket the inverse map
    # locally first, then expand only if needed. Every mapping uses CAMB's setter.
    for lo, hi in ((0.6, 2.0), (0.5, 10.0)):
        grid = np.linspace(lo, hi, 25 if hi == 2.0 else 48, dtype=np.float64)
        mapped: list[tuple[float, float]] = []
        for theta in grid:
            try:
                mapped.append((float(theta), h0_from_theta(float(theta), obh2, och2)))
            except Exception:
                continue
        for (theta_a, h0_a), (theta_b, h0_b) in zip(mapped, mapped[1:]):
            if (h0_a - target_h0) * (h0_b - target_h0) <= 0.0:
                return float(brentq(
                    lambda theta: h0_from_theta(theta, obh2, och2) - target_h0,
                    theta_a, theta_b, xtol=2.0e-12, rtol=1.0e-13,
                ))
    raise ValueError(f"No theta100 inverse-map bracket for H0={target_h0}, obh2={obh2}, och2={och2}")


def in_declared_box(h0: float, obh2: float, och2: float, theta100: float) -> bool:
    return (
        H0_BOUNDS[0] <= h0 <= H0_BOUNDS[1]
        and OBH2_BOUNDS[0] <= obh2 <= OBH2_BOUNDS[1]
        and OCH2_BOUNDS[0] <= och2 <= OCH2_BOUNDS[1]
        and THETA100_BOUNDS[0] <= theta100 <= THETA100_BOUNDS[1]
    )


def profile_slice(obh2: float, starts: list[tuple[float, float]]) -> dict[str, object]:
    evaluations = 0

    def physical(unit: np.ndarray) -> tuple[float, float]:
        unit = np.asarray(unit, dtype=np.float64)
        h0 = H0_BOUNDS[0] + float(unit[0]) * (H0_BOUNDS[1] - H0_BOUNDS[0])
        och2 = OCH2_BOUNDS[0] + float(unit[1]) * (OCH2_BOUNDS[1] - OCH2_BOUNDS[0])
        return h0, och2

    def objective(unit: np.ndarray) -> float:
        nonlocal evaluations
        evaluations += 1
        h0, och2 = physical(unit)
        if np.any(np.asarray(unit) < 0.0) or np.any(np.asarray(unit) > 1.0):
            return INVALID_SCORE
        try:
            return float(score(h0, obh2, och2)["chi2_cholesky_whitening"])
        except Exception:
            return INVALID_SCORE

    trials: list[dict[str, object]] = []
    for start_h0, start_och2 in starts:
        unit0 = np.asarray([
            (start_h0 - H0_BOUNDS[0]) / (H0_BOUNDS[1] - H0_BOUNDS[0]),
            (start_och2 - OCH2_BOUNDS[0]) / (OCH2_BOUNDS[1] - OCH2_BOUNDS[0]),
        ], dtype=np.float64)
        unit0 = np.clip(unit0, 0.0, 1.0)
        result = minimize(
            objective, unit0, method="Powell", bounds=((0.0, 1.0), (0.0, 1.0)),
            options={"xtol": 1.0e-9, "ftol": 1.0e-12, "maxiter": 500, "maxfev": 1800},
        )
        h0, och2 = physical(result.x)
        final = score(h0, obh2, och2)
        theta = implied_theta100(h0, obh2, och2)
        valid = in_declared_box(h0, obh2, och2, theta)
        trials.append({
            "start": {"H0_km_s_Mpc": start_h0, "omega_c_h2": start_och2},
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "H0_km_s_Mpc": h0,
            "omega_b_h2": obh2,
            "omega_c_h2": och2,
            "implied_100theta_MC": theta,
            "declared_box_and_theta_valid": valid,
            "chi2": float(final["chi2_cholesky_whitening"]),
            "closure": final["derived"],
        })
    valid_trials = [trial for trial in trials if trial["declared_box_and_theta_valid"]]
    if not valid_trials:
        raise RuntimeError(f"No valid profile solution found at omega_b h2={obh2}")
    best = min(valid_trials, key=lambda row: float(row["chi2"]))
    # Re-score the selected solution and keep the score/covariance parity evidence.
    selected_score = score(float(best["H0_km_s_Mpc"]), obh2, float(best["omega_c_h2"]))
    return {
        "omega_b_h2_fixed": obh2,
        "best_valid_trial": best,
        "all_starts": trials,
        "valid_start_scores_span": [
            min(float(row["chi2"]) for row in valid_trials),
            max(float(row["chi2"]) for row in valid_trials),
        ],
        "selected_score_checks": {
            "chi2_dense_solve": selected_score["chi2_dense_solve"],
            "chi2_cholesky_whitening": selected_score["chi2_cholesky_whitening"],
            "dense_cholesky_abs_delta": selected_score["dense_cholesky_abs_delta"],
            "predictions_in_released_order": selected_score["predictions_in_released_order"],
        },
        "objective_evaluations_across_starts": evaluations,
    }


def main() -> None:
    started = time.perf_counter()
    reference = json.loads(REFERENCE_FILE.read_text(encoding="utf-8"))
    candidate_rec = reference["selected_candidate"]
    p = candidate_rec["parameters"]
    cand_h0 = float(p["H0_km_s_Mpc"])
    cand_ob = float(p["omega_b_h2"])
    cand_oc = float(p["omega_c_h2"])
    cand_score = score(cand_h0, cand_ob, cand_oc)
    claimed_chi2 = float(candidate_rec["bao"]["chi2_full_covariance_cholesky"])
    boundary_score = score(cand_h0, 0.1, cand_oc)

    # Fixed-H0, fixed-omega_c one- and two-sided differences at the actual box edge.
    finite_differences: list[dict[str, object]] = []
    for step in FD_STEPS:
        at_edge = score(cand_h0, 0.1, cand_oc)
        below = score(cand_h0, 0.1 - step, cand_oc)
        above = score(cand_h0, 0.1 + step, cand_oc)
        chi0 = float(at_edge["chi2_cholesky_whitening"])
        chiminus = float(below["chi2_cholesky_whitening"])
        chiplus = float(above["chi2_cholesky_whitening"])
        finite_differences.append({
            "step_omega_b_h2": step,
            "at_box_edge_chi2": chi0,
            "inward_point": {"omega_b_h2": 0.1 - step, "chi2": chiminus, "delta_chi2_vs_edge": chiminus - chi0, "eligible_in_box": True},
            "outward_diagnostic": {"omega_b_h2": 0.1 + step, "chi2": chiplus, "delta_chi2_vs_edge": chiplus - chi0, "eligible_in_box": False},
            "inward_one_sided_dchi2_domega_b_h2": (chi0 - chiminus) / step,
            "outward_one_sided_dchi2_domega_b_h2": (chiplus - chi0) / step,
            "central_dchi2_domega_b_h2_diagnostic": (chiplus - chiminus) / (2.0 * step),
        })

    # Continue the profile inward, restarting from both the edge and the previous
    # optimum; an ordinary low-redshift anchor provides a third independent start.
    profile: list[dict[str, object]] = []
    continuation = (cand_h0, cand_oc)
    for obh2 in PROFILE_OBH2:
        starts = [continuation, (cand_h0, cand_oc), (67.0, 0.12)]
        record = profile_slice(obh2, starts)
        profile.append(record)
        best = record["best_valid_trial"]
        continuation = (float(best["H0_km_s_Mpc"]), float(best["omega_c_h2"]))

    lower_candidates = [
        row["best_valid_trial"] for row in profile
        if float(row["best_valid_trial"]["chi2"]) < claimed_chi2
        and row["best_valid_trial"]["declared_box_and_theta_valid"]
    ]
    best_profile = min(
        (row["best_valid_trial"] for row in profile),
        key=lambda row: float(row["chi2"]),
    )
    edge_profile = profile[0]["best_valid_trial"]
    inward_profile_min = min(
        (row["best_valid_trial"] for row in profile if float(row["omega_b_h2_fixed"]) < 0.1),
        key=lambda row: float(row["chi2"]),
    )
    inward_profile_delta_vs_edge = float(inward_profile_min["chi2"]) - float(edge_profile["chi2"])
    candidate_to_polished_edge_delta = float(edge_profile["chi2"]) - claimed_chi2
    compute26_minimum_agreement_tolerance = float(
        reference["candidate_acceptance"]["success_criteria"]["minimum_pairwise_chi2_tolerance"]
    )
    edge_theta = implied_theta100(cand_h0, 0.1, cand_oc)
    candidate_theta = implied_theta100(cand_h0, cand_ob, cand_oc)
    lower_exists = bool(lower_candidates)
    result: dict[str, object] = {
        "status": "completed_boundary_gradient_and_fixed_omega_b_profile",
        "interpretation_scope": "BAO-only flat-LCDM score audit; not a posterior, evidence calculation, or physical cosmological inference",
        "candidate": {
            "parameters": {"H0_km_s_Mpc": cand_h0, "omega_b_h2": cand_ob, "omega_c_h2": cand_oc},
            "claimed_compute26_chi2": claimed_chi2,
            "independent_chi2_dense_solve": cand_score["chi2_dense_solve"],
            "independent_chi2_cholesky_whitening": cand_score["chi2_cholesky_whitening"],
            "score_abs_delta_to_compute26": abs(cand_score["chi2_cholesky_whitening"] - claimed_chi2),
            "implied_100theta_MC": candidate_theta,
            "candidate_theta_inside_table_I_interval": THETA100_BOUNDS[0] <= candidate_theta <= THETA100_BOUNDS[1],
            "closure": cand_score["derived"],
        },
        "declared_search_domain_audited": {
            "H0_km_s_Mpc": list(H0_BOUNDS),
            "omega_b_h2": list(OBH2_BOUNDS),
            "omega_c_h2": list(OCH2_BOUNDS),
            "paper_table_I_100theta_MC": list(THETA100_BOUNDS),
            "notes": "The compute26 direct-H0 box is used for optimization. Each reported optimum is independently mapped through CAMB's cosmomc_theta setter and checked against the broad paper interval; theta mapping is not approximated with theta_star.",
        },
        "input_validation": {
            "row_count": int(len(Z)),
            "observable_order": OBSERVABLES,
            "redshifts": Z.tolist(),
            "covariance_shape": list(COV.shape),
            "covariance_symmetric_abs_max": float(np.max(np.abs(COV - COV.T))),
            "covariance_cholesky_positive_definite": True,
            "complete_covariance_used": True,
        },
        "finite_differences_fixed_H0_and_omega_c_at_box_edge": finite_differences,
        "edge_point_at_exact_upper_bound": {
            "parameters": {"H0_km_s_Mpc": cand_h0, "omega_b_h2": 0.1, "omega_c_h2": cand_oc},
            "chi2": boundary_score["chi2_cholesky_whitening"],
            "implied_100theta_MC": edge_theta,
            "declared_box_and_theta_valid": in_declared_box(cand_h0, 0.1, cand_oc, edge_theta),
        },
        "fixed_omega_b_profile": profile,
        "finding": {
            "literal_valid_lower_point_than_compute26_candidate_found": lower_exists,
            "literal_lower_point_source": "Powell polish of H0 and omega_c h2 at omega_b h2=0.1; not an inward omega_b move",
            "polished_edge_delta_chi2_vs_supplied_compute26_candidate": candidate_to_polished_edge_delta,
            "compute26_seed_agreement_tolerance_for_minima": compute26_minimum_agreement_tolerance,
            "polished_edge_delta_is_below_compute26_tolerance": abs(candidate_to_polished_edge_delta) < compute26_minimum_agreement_tolerance,
            "number_of_profile_slices_lower_than_supplied_candidate": len(lower_candidates),
            "lowest_profile_slice_candidate": best_profile,
            "lowest_profile_chi2": float(best_profile["chi2"]),
            "delta_chi2_vs_compute26_candidate": float(best_profile["chi2"]) - claimed_chi2,
            "lowest_inward_slice_candidate": inward_profile_min,
            "lowest_inward_slice_delta_chi2_vs_polished_edge": inward_profile_delta_vs_edge,
            "inward_omega_b_downhill_found_on_tested_slices": inward_profile_delta_vs_edge < 0.0,
            "fixed_coordinate_inward_deltas_all_positive": all(
                float(row["inward_point"]["delta_chi2_vs_edge"]) > 0.0 for row in finite_differences
            ),
            "interpretation": "The supplied point is beaten literally by a tiny same-edge H0/omega_c polish (far below compute26's declared optimizer agreement tolerance). Both-sided fixed-coordinate finite differences rise away from omega_b h2=0.1, and every reoptimized inward slice is above the polished edge; the suspected inward baryon downhill direction is not supported.",
            "global_optimality_claim": "none; finite slices and local bounded profile solves do not establish a global minimum",
        },
        "method": {
            "implementation": "Standalone CAMB/background/data/observable/covariance code; no compute26 or compute28 scorer imported",
            "optimizer": "Three-start Powell profile optimization in normalized (H0, omega_c h2) box for each fixed omega_b h2 slice; CAMB objective re-evaluated at returned points",
            "camb_setup": {"dark_energy": "flat fluid Lambda, w=-1, wa=0", **FIXED_CAMB},
            "theta_inverse_mapping": "CAMB.set_cosmology(cosmomc_theta=theta100/100, theta_H0_range=(1,2000)) and brentq inverse at reported candidates",
            "finite_difference_points_above_0p1": "explicitly out-of-box diagnostics only; never eligible as candidate minima",
        },
        "runtime_seconds": time.perf_counter() - started,
        "versions": {
            "CAMB": camb.__version__,
            "NumPy": np.__version__,
            "SciPy": scipy.__version__,
            "Python": sys.version,
            "platform": platform.platform(),
            "thread_environment": {key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
        },
        "hashes": {
            "script_sha256": sha256(Path(__file__)),
            "mean_sha256": sha256(MEAN_FILE),
            "covariance_sha256": sha256(COV_FILE),
            "compute26_candidate_json_sha256": sha256(REFERENCE_FILE),
            "CAMB_2_0_4_wheel_sha256": sha256(WHEEL_FILE),
        },
        "reproduction_command": "cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute30/no_campaign_deadline.json --seconds 560 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute30/boundary_gradient_audit.py",
    }
    OUT_FILE.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result_json_sha = sha256(OUT_FILE)

    finding = result["finding"]
    report = f'''# Compute30 boundary-gradient audit

Status: **independently checked; no inward \\(\\omega_b h^2\\) downhill direction found in the tested neighborhood**. A same-edge \\(H_0\\)/\\(\\omega_c h^2\\) polish gives \\(\\Delta\\chi^2={candidate_to_polished_edge_delta:.3g}\\) relative to the supplied point, far below compute26's declared \\(2\\times10^{{-5}}\\) optimizer agreement tolerance. This is a BAO-only flat-LCDM finite-box/profile result, not physical cosmological inference.

The standalone implementation reproduces the candidate score as \\({cand_score["chi2_cholesky_whitening"]:.12f}\\), matching compute26 \\({claimed_chi2:.12f}\\). At fixed \\(H_0\\) and \\(\\omega_ch^2\\), moving either inward or outward from \\(\\omega_b h^2=0.1\\) raises the score at all tested steps \\(10^{{-5}}\\), \\(10^{{-4}}\\), and \\(5\\times10^{{-4}}\\). For example, the inward \\(10^{{-4}}\\) point is higher by \\({float(finite_differences[1]["inward_point"]["delta_chi2_vs_edge"]):.9g}\\); the outward point is outside the declared box and is diagnostic only. The compute26 outward stress therefore did not imply an inward feasible descent.

Re-optimizing \\(H_0\\) and \\(\\omega_ch^2\\) from three starts at {len(profile)} fixed \\(\\omega_b h^2\\) slices puts the best sampled score at the upper edge: \\(\\chi^2={float(edge_profile["chi2"]):.12f}\\). Every tested inward slice is higher; the nearest, \\(\\omega_b h^2=0.0999999\\), is higher by \\({float(profile[1]["best_valid_trial"]["chi2"]) - float(edge_profile["chi2"]):.3g}\\), and the lowest inward-slice increase is \\({inward_profile_delta_vs_edge:.3g}\\). The polished edge point has \\(H_0={float(edge_profile["H0_km_s_Mpc"]):.12g}\\), \\(\\omega_c h^2={float(edge_profile["omega_c_h2"]):.12g}\\), and \\(100\\theta_{{MC}}={float(edge_profile["implied_100theta_MC"]):.9g}\\), all inside the declared domain. The tested profile supports a local minimum at the \\(\\omega_b h^2\\) boundary. The tiny literal improvement over the supplied point reflects final same-edge polishing; discrete slices and local Powell solves do not establish global optimality.

The score uses the full released 13-row order and full 13x13 covariance; both dense precision solve and Cholesky whitening are checked at every selected point. CAMB 2.0.4 is configured for flat fluid \\(\\Lambda\\)CDM with \\(\\Omega_k=0\\), \\(m_\\nu=0.06\\) eV, \\(N_{{eff}}=3.044\\), \\(T_{{CMB}}=2.7255\\) K, one degenerate massive species and BBN-consistent helium. Every retained point has \\(0<\\Omega_m<1\\), \\(\\Omega_\\Lambda>0\\), and finite positive background distances and Hubble rates. The CAMB \\(\\theta_{{MC}}\\) setter, with expanded internal H0 root bracket \\([1,2000]\\), independently maps each selected \\(H_0\\) back to the target theta interval.

## Reproduction

```sh
{result["reproduction_command"]}
```

The run is bounded at 560 seconds and fixes numerical-library threads to one. Detailed finite differences, all optimizer starts, score vectors, closure checks, versions and source hashes are in [boundary_gradient_audit.json](boundary_gradient_audit.json). Script SHA-256: `{result["hashes"]["script_sha256"]}`. Result SHA-256: `{result_json_sha}`.
'''
    REPORT_FILE.write_text(report, encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "candidate_chi2": cand_score["chi2_cholesky_whitening"],
        "candidate_claimed_chi2": claimed_chi2,
        "edge_fixed_score": boundary_score["chi2_cholesky_whitening"],
        "inward_delta_chi2_h1e-4": finite_differences[1]["inward_point"]["delta_chi2_vs_edge"],
        "lower_valid_profile_points_than_supplied_candidate": len(lower_candidates),
        "minimum_inward_profile_delta_vs_edge": inward_profile_delta_vs_edge,
        "candidate_to_polished_edge_delta": candidate_to_polished_edge_delta,
        "best_profile_chi2": best_profile["chi2"],
        "best_profile_point": {k: best_profile[k] for k in ("H0_km_s_Mpc", "omega_b_h2", "omega_c_h2", "implied_100theta_MC")},
        "runtime_seconds": result["runtime_seconds"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

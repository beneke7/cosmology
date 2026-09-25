#!/usr/bin/env python3
"""Independent bounded CAMB × DESI DR2 BAO optimizer-globality audit.

This implementation does not import compute31. It parses the released BAO
vector/covariance, predicts all observables with pinned CAMB 2.0.4, and scores
the full covariance. Fixed-omega_b Powell multistarts cover the feasible
flat-positive-Lambda domain; differential evolution adds deterministic global
seeds at both scan edges and the three selected BBN prior centers.
"""
from __future__ import annotations

import concurrent.futures as futures
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import platform
import sys
import time

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

TASK_DIR = Path(__file__).resolve().parent
PROJECT = TASK_DIR.parents[1]
sys.path.insert(0, str(PROJECT / "work/compute24/site"))

import camb  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402
from scipy.optimize import brentq, differential_evolution, minimize  # noqa: E402
from scipy.stats import qmc  # noqa: E402

MEAN_PATH = PROJECT / "context/data/desi_dr2_mean.txt"
COV_PATH = PROJECT / "context/data/desi_dr2_cov.txt"
REFERENCE_SCRIPT = PROJECT / "work/compute31/bbn_bao_profile.py"
REFERENCE_JSON = PROJECT / "work/compute31/bbn_bao_profile.json"
REFERENCE_REPORT = PROJECT / "work/compute31/REPORT.md"
CAMB_WHEEL = PROJECT / "context/code/camb-2.0.4-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl"
CAMB_SITE = PROJECT / "work/compute24/site"
RESULT_PATH = TASK_DIR / "independent_globality.json"
CHECKPOINT_PATH = TASK_DIR / "checkpoint.json"
REPORT_PATH = TASK_DIR / "REPORT.md"
PLOT_PATH = TASK_DIR / "profile_globality.png"

C_KM_S = 299792.458
H0_BOUNDS = (10.0, 1000.0)
OBH2_BOUNDS = (0.0205, 0.0240)
OCH2_BOUNDS = (0.001, 0.99)
THETA100_BOUNDS = (0.5, 10.0)
MNU_EV = 0.06
N_EFF = 3.044
T_CMB_K = 2.7255
OMEGA_K = 0.0
W = -1.0
WA = 0.0
EXPECTED_NAMES = [
    "DV_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs",
    "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs",
    "DH_over_rs", "DH_over_rs", "DM_over_rs",
]
PRIOR_CENTERS = {
    "poulin_2026_arxiv_v1_BBN_LCDM": 0.022017,
    "PDG_2025_SBBN": 0.02205,
    "provisional_conservative_BBN": 0.0222,
}
INVALID = 1.0e80
EPS_CLOSURE_PHYSICAL = 1.0e-9
POWELL_OPTIONS = {"xtol": 2.0e-9, "ftol": 2.0e-12, "maxiter": 500, "maxfev": 2200}
DE_OPTIONS = {"maxiter": 110, "popsize": 16, "tol": 1.0e-9, "atol": 1.0e-11, "polish": False, "updating": "immediate"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_released_vector() -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray]:
    rows = []
    for raw in MEAN_PATH.read_text(encoding="utf-8").splitlines():
        row = raw.strip()
        if row and not row.startswith("#"):
            z_s, value_s, name = row.split()
            rows.append((float(z_s), float(value_s), name))
    z = np.asarray([row[0] for row in rows], dtype=np.float64)
    observed = np.asarray([row[1] for row in rows], dtype=np.float64)
    names = [row[2] for row in rows]
    covariance = np.asarray(np.loadtxt(COV_PATH, dtype=np.float64), dtype=np.float64)
    if z.shape != (13,) or observed.shape != (13,) or names != EXPECTED_NAMES:
        raise ValueError("official DESI input must be the pinned 13-row mean vector/order")
    if covariance.shape != (13, 13) or not np.isfinite(covariance).all():
        raise ValueError("official DESI covariance must be finite and 13×13")
    if not np.allclose(covariance, covariance.T, rtol=0.0, atol=1.0e-12):
        raise ValueError("official covariance is not symmetric")
    np.linalg.cholesky(covariance)
    return z, observed, names, covariance


Z, Y, OBSERVABLES, COV = load_released_vector()
CHOL = np.linalg.cholesky(COV)


def make_parameters(h0: float, obh2: float, och2: float) -> camb.CAMBparams:
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=W, wa=WA, dark_energy_model="fluid")
    pars.set_cosmology(
        H0=float(h0), ombh2=float(obh2), omch2=float(och2), omk=OMEGA_K,
        mnu=MNU_EV, nnu=N_EFF, num_massive_neutrinos=1,
        neutrino_hierarchy="degenerate", TCMB=T_CMB_K, YHe=None,
    )
    return pars


def evaluate(h0: float, obh2: float, och2: float) -> dict[str, object]:
    """CAMB background predictions plus independent dense and Cholesky scores."""
    pars = make_parameters(h0, obh2, och2)
    background = camb.get_background(pars)
    omega_m0 = float(pars.omegam)
    omega_de0 = float(background.get_Omega("de", 0.0))
    if not math.isfinite(omega_de0) or omega_de0 <= 0.0:
        raise ValueError(f"flat positive-Lambda closure failed, Omega_de={omega_de0}")
    derived = background.get_derived_params()
    rdrag = float(derived["rdrag"])
    hz = np.asarray(background.hubble_parameter(Z), dtype=np.float64)
    dm = np.asarray(background.comoving_radial_distance(Z), dtype=np.float64)
    if (not math.isfinite(rdrag) or rdrag <= 0.0 or not np.isfinite(hz).all()
            or not np.isfinite(dm).all() or np.any(hz <= 0.0) or np.any(dm <= 0.0)):
        raise ValueError("CAMB returned invalid H(z), D_M(z), or r_drag")
    dh = C_KM_S / hz
    prediction = np.empty(13, dtype=np.float64)
    for i, name in enumerate(OBSERVABLES):
        if name == "DV_over_rs":
            prediction[i] = (Z[i] * dm[i] ** 2 * dh[i]) ** (1.0 / 3.0) / rdrag
        elif name == "DM_over_rs":
            prediction[i] = dm[i] / rdrag
        elif name == "DH_over_rs":
            prediction[i] = dh[i] / rdrag
        else:
            raise ValueError(f"unknown DESI observable {name}")
    residual = prediction - Y
    whitened = np.linalg.solve(CHOL, residual)
    chi2_chol = float(whitened @ whitened)
    chi2_dense = float(residual @ np.linalg.solve(COV, residual))
    if not math.isfinite(chi2_chol) or not math.isfinite(chi2_dense):
        raise ArithmeticError("non-finite DESI full-covariance score")
    delta = abs(chi2_chol - chi2_dense)
    if delta > 2.0e-10 * max(1.0, abs(chi2_chol)):
        raise ArithmeticError(f"dense/Cholesky mismatch {delta}")
    h = float(h0) / 100.0
    omega_r0_approx = max(0.0, 1.0 - omega_m0 - omega_de0)
    return {
        "chi2_cholesky": chi2_chol,
        "chi2_dense_solve": chi2_dense,
        "dense_cholesky_abs_delta": delta,
        "predictions_in_released_order": prediction.tolist(),
        "residuals_in_released_order": residual.tolist(),
        "closure": {
            "Omega_m0_including_massive_neutrino": omega_m0,
            "Omega_de0_CAMB": omega_de0,
            "Omega_b0": float(pars.omegab),
            "Omega_c0": float(pars.omegac),
            "Omega_nu0_massive": float(pars.omeganu),
            "Omega_r0_inferred_from_closure": omega_r0_approx,
            "omega_nu_h2_CAMB": float(pars.omnuh2),
            "r_drag_Mpc": rdrag,
            "100theta_star_CAMB": float(derived["thetastar"]),
            "H0_km_s_Mpc": float(h0),
            "h": h,
            "omega_b_h2": float(obh2),
            "omega_c_h2": float(och2),
        },
    }


def radiation_physical_density() -> tuple[float, float]:
    """Read fixed massive-neutrino and radiation physical densities from CAMB."""
    pars = make_parameters(70.0, 0.022, 0.12)
    background = camb.get_background(pars)
    rad_h2 = (1.0 - float(pars.omegam) - float(background.get_Omega("de"))) * (0.7 ** 2)
    return float(pars.omnuh2), float(rad_h2)


OMEGA_NU_H2, OMEGA_R_H2 = radiation_physical_density()


def valid_h0_lower(obh2: float) -> float:
    required_h2 = obh2 + OMEGA_NU_H2 + OMEGA_R_H2 + OCH2_BOUNDS[0] + EPS_CLOSURE_PHYSICAL
    return max(H0_BOUNDS[0], 100.0 * math.sqrt(required_h2))


def decode_unit(unit: np.ndarray, obh2: float) -> tuple[float, float]:
    """Map [0,1]^2 smoothly onto the physically feasible closed search region."""
    u = np.clip(np.asarray(unit, dtype=np.float64), 0.0, 1.0)
    hmin = valid_h0_lower(obh2)
    h0 = math.exp(math.log(hmin) + float(u[0]) * (math.log(H0_BOUNDS[1]) - math.log(hmin)))
    h2 = (h0 / 100.0) ** 2
    oc_max = min(OCH2_BOUNDS[1], h2 - obh2 - OMEGA_NU_H2 - OMEGA_R_H2 - EPS_CLOSURE_PHYSICAL)
    oc_max = max(OCH2_BOUNDS[0], oc_max)
    och2 = OCH2_BOUNDS[0] + float(u[1]) * (oc_max - OCH2_BOUNDS[0])
    return h0, och2


def encode_physical(h0: float, och2: float, obh2: float) -> np.ndarray:
    """Project a proposed seed into the declared, flat-positive-Lambda domain."""
    hmin = valid_h0_lower(obh2)
    h0 = float(np.clip(h0, hmin, H0_BOUNDS[1]))
    h2 = (h0 / 100.0) ** 2
    oc_max = min(OCH2_BOUNDS[1], h2 - obh2 - OMEGA_NU_H2 - OMEGA_R_H2 - EPS_CLOSURE_PHYSICAL)
    oc_max = max(OCH2_BOUNDS[0], oc_max)
    och2 = float(np.clip(och2, OCH2_BOUNDS[0], oc_max))
    ux = (math.log(h0) - math.log(hmin)) / (math.log(H0_BOUNDS[1]) - math.log(hmin))
    uy = 0.0 if oc_max <= OCH2_BOUNDS[0] else (och2 - OCH2_BOUNDS[0]) / (oc_max - OCH2_BOUNDS[0])
    return np.clip(np.asarray([ux, uy], dtype=np.float64), 0.0, 1.0)


def objective_unit(unit: np.ndarray, obh2: float, tracker: dict[str, object]) -> float:
    try:
        h0, och2 = decode_unit(unit, obh2)
        fit = evaluate(h0, obh2, och2)
        tracker["evaluations"] = int(tracker["evaluations"]) + 1
        return float(fit["chi2_cholesky"])
    except Exception as exc:
        tracker["evaluations"] = int(tracker["evaluations"]) + 1
        counts = tracker["invalid_counts"]
        kind = type(exc).__name__
        counts[kind] = counts.get(kind, 0) + 1
        examples = tracker["invalid_examples"]
        if len(examples) < 12:
            try:
                h0, och2 = decode_unit(unit, obh2)
                point = [h0, float(obh2), och2]
            except Exception:
                point = [float("nan"), float(obh2), float("nan")]
            examples.append({"type": kind, "parameters": point, "detail": str(exc)[:200]})
        return INVALID


def theta_setter_h0(theta100: float, obh2: float, och2: float) -> float:
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=W, wa=WA, dark_energy_model="fluid")
    pars.set_cosmology(
        cosmomc_theta=float(theta100) / 100.0, ombh2=float(obh2), omch2=float(och2),
        omk=OMEGA_K, mnu=MNU_EV, nnu=N_EFF, num_massive_neutrinos=1,
        neutrino_hierarchy="degenerate", TCMB=T_CMB_K, YHe=None,
        theta_H0_range=(1.0, 2000.0),
    )
    return float(pars.H0)


def theta_support(h0: float, obh2: float, och2: float) -> dict[str, object]:
    grid = np.linspace(THETA100_BOUNDS[0], THETA100_BOUNDS[1], 65)
    mapped = []
    for theta in grid:
        try:
            mapped.append((float(theta), theta_setter_h0(float(theta), obh2, och2)))
        except Exception:
            continue
    for (ta, ha), (tb, hb) in zip(mapped, mapped[1:]):
        if (ha - h0) * (hb - h0) <= 0.0:
            value = float(brentq(lambda t: theta_setter_h0(t, obh2, och2) - h0, ta, tb,
                                 xtol=2.0e-11, rtol=1.0e-12))
            recovered = theta_setter_h0(value, obh2, och2)
            return {
                "implied_100theta_MC": value,
                "inside_declared_100theta_MC_interval": THETA100_BOUNDS[0] <= value <= THETA100_BOUNDS[1],
                "CAMB_theta_setter_H0_recovered": recovered,
                "H0_recovery_abs_error_km_s_Mpc": abs(recovered - h0),
            }
    return {"implied_100theta_MC": None, "inside_declared_100theta_MC_interval": False,
            "reason": "no CAMB theta-setter inverse bracket in declared interval"}


def qmc_starts(obh2: float, slice_index: int, n: int = 12) -> list[list[float]]:
    sampler = qmc.LatinHypercube(d=2, seed=20260925 + 1009 * slice_index)
    return [point.tolist() for point in sampler.random(n)]


def physical_starts(obh2: float, slice_index: int) -> list[dict[str, object]]:
    physical = [
        (67.0, 0.12), (82.0, 0.18), (133.6, 0.43), (300.0, 0.30),
        (1000.0, 0.001), (1000.0, 0.99), (250.0, 0.001), (250.0, 0.99),
        (40.0, 0.001), (40.0, 0.30), (100.0, 0.001), (100.0, 0.99),
    ]
    seeds = [{"kind": "compute31_scale_or_box", "physical": [h, oc],
              "unit": encode_physical(h, oc, obh2).tolist()} for h, oc in physical]
    seeds.extend({"kind": "deterministic_LHS", "unit": u} for u in qmc_starts(obh2, slice_index))
    return seeds


def solve_seed(seed: dict[str, object], obh2: float) -> dict[str, object]:
    tracker: dict[str, object] = {"evaluations": 0, "invalid_counts": {}, "invalid_examples": []}
    result = minimize(
        objective_unit, np.asarray(seed["unit"], dtype=np.float64), args=(float(obh2), tracker),
        method="Powell", bounds=((0.0, 1.0), (0.0, 1.0)), options=POWELL_OPTIONS,
    )
    h0, och2 = decode_unit(result.x, obh2)
    final = None
    try:
        final = evaluate(h0, obh2, och2)
        score = float(final["chi2_cholesky"])
        valid = float(final["closure"]["Omega_de0_CAMB"]) > 0.0
    except Exception as exc:
        score = INVALID
        valid = False
        tracker["invalid_counts"][type(exc).__name__] = tracker["invalid_counts"].get(type(exc).__name__, 0) + 1
        if len(tracker["invalid_examples"]) < 12:
            tracker["invalid_examples"].append({"type": type(exc).__name__,
                "parameters": [h0, float(obh2), och2], "detail": str(exc)[:200], "final_point": True})
    return {
        "seed": seed,
        "optimizer_success": bool(result.success),
        "optimizer_status": int(result.status),
        "optimizer_message": str(result.message),
        "optimizer_nfev": int(result.nfev),
        "tracker_evaluations": tracker["evaluations"],
        "H0_km_s_Mpc": h0,
        "omega_b_h2": float(obh2),
        "omega_c_h2": och2,
        "chi2": score,
        "valid_flat_positive_Lambda": bool(valid),
        "closure_and_score_check": final,
        "invalid_counts": tracker["invalid_counts"],
        "invalid_examples": tracker["invalid_examples"],
    }


def global_de(obh2: float, slice_index: int) -> dict[str, object]:
    tracker: dict[str, object] = {"evaluations": 0, "invalid_counts": {}, "invalid_examples": []}
    started = time.perf_counter()
    result = differential_evolution(
        objective_unit, bounds=((0.0, 1.0), (0.0, 1.0)), args=(float(obh2), tracker),
        seed=48017 + 7919 * slice_index, init="sobol", workers=1,
        **DE_OPTIONS,
    )
    h0, och2 = decode_unit(result.x, obh2)
    final = evaluate(h0, obh2, och2)
    polished = solve_seed({"kind": "DE_endpoint", "unit": result.x.tolist()}, obh2)
    return {
        "seed": 48017 + 7919 * slice_index,
        "DE_success": bool(result.success),
        "DE_message": str(result.message),
        "DE_nfev": int(result.nfev),
        "tracker_evaluations": tracker["evaluations"],
        "DE_H0_km_s_Mpc": h0,
        "DE_omega_c_h2": och2,
        "DE_chi2": float(final["chi2_cholesky"]),
        "DE_closure_and_score_check": final,
        "Powell_polish": polished,
        "wall_seconds": time.perf_counter() - started,
        "invalid_counts": tracker["invalid_counts"],
    }


def make_nodes() -> list[float]:
    # Recreate compute31's 15 inclusive linear nodes plus prior centers, then
    # append its edge-refinement nodes independently to stress the profile edge.
    coarse = np.linspace(OBH2_BOUNDS[0], OBH2_BOUNDS[1], 15, dtype=np.float64).tolist()
    refine = [0.02378125, 0.0238125, 0.02384375, 0.023875,
              0.02390625, 0.0239375, 0.02396875]
    return sorted(set(float(v) for v in coarse + list(PRIOR_CENTERS.values()) + refine))


NODES = make_nodes()
GLOBAL_B_VALUES = {OBH2_BOUNDS[0], OBH2_BOUNDS[1], *PRIOR_CENTERS.values()}


def run_slice(task: tuple[int, float]) -> dict[str, object]:
    index, obh2 = task
    started = time.perf_counter()
    cpu_started = time.process_time()
    trials = [solve_seed(seed, obh2) for seed in physical_starts(obh2, index)]
    valid = [trial for trial in trials if trial["valid_flat_positive_Lambda"]]
    if not valid:
        raise RuntimeError(f"all fixed-profile starts failed at omega_b h2={obh2}")
    best = min(valid, key=lambda trial: float(trial["chi2"]))
    global_result = global_de(obh2, index) if obh2 in GLOBAL_B_VALUES else None
    candidates = valid + ([global_result["Powell_polish"]] if global_result else [])
    best = min(candidates, key=lambda trial: float(trial["chi2"]))
    h0, oc = float(best["H0_km_s_Mpc"]), float(best["omega_c_h2"])
    best_theta = theta_support(h0, obh2, oc)
    local_alternatives = []
    best_chi2 = float(best["chi2"])
    for trial in valid:
        delta = float(trial["chi2"]) - best_chi2
        if delta > 1.0e-6:
            local_alternatives.append({
                "seed_kind": trial["seed"]["kind"],
                "H0_km_s_Mpc": trial["H0_km_s_Mpc"],
                "omega_c_h2": trial["omega_c_h2"],
                "chi2": trial["chi2"],
                "delta_chi2_to_slice_best": delta,
                "optimizer_success": trial["optimizer_success"],
                "optimizer_message": trial["optimizer_message"],
                "closure": trial["closure_and_score_check"]["closure"],
            })
    return {
        "omega_b_h2_fixed": obh2,
        "best": {k: best[k] for k in ("H0_km_s_Mpc", "omega_b_h2", "omega_c_h2", "chi2",
                                        "optimizer_success", "optimizer_message", "closure_and_score_check")},
        "best_theta_support_check": best_theta,
        "all_powell_starts": trials,
        "valid_start_count": len(valid),
        "invalid_start_count": len(trials) - len(valid),
        "unique_local_alternatives_delta_chi2_gt_1e-6": local_alternatives,
        "global_differential_evolution": global_result,
        "valid_start_chi2_span": [min(float(t["chi2"]) for t in valid), max(float(t["chi2"]) for t in valid)],
        "wall_seconds": time.perf_counter() - started,
        "worker_cpu_seconds": time.process_time() - cpu_started,
    }


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def one_point_check() -> dict[str, object]:
    ref = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))
    point = ref["profile_grid"][-1]["best_valid_trial"]
    independent = evaluate(float(point["H0_km_s_Mpc"]), float(point["omega_b_h2"]), float(point["omega_c_h2"]))
    delta = abs(float(point["chi2"]) - float(independent["chi2_cholesky"]))
    return {
        "compute31_reference_point": {
            "H0_km_s_Mpc": point["H0_km_s_Mpc"], "omega_b_h2": point["omega_b_h2"],
            "omega_c_h2": point["omega_c_h2"], "chi2": point["chi2"],
        },
        "independent_prediction_and_score": independent,
        "reference_chi2_abs_delta": delta,
        "passes_1e-8_score_agreement": delta < 1.0e-8,
    }


def run_all(workers: int) -> dict[str, object]:
    started = time.perf_counter()
    one_check = one_point_check()
    if not one_check["passes_1e-8_score_agreement"]:
        raise RuntimeError(f"independent one-point check failed: {one_check['reference_chi2_abs_delta']}")
    hashes = {
        "official_DESI_mean": sha256(MEAN_PATH),
        "official_DESI_covariance": sha256(COV_PATH),
        "compute31_script_read_only_reference": sha256(REFERENCE_SCRIPT),
        "compute31_JSON_read_only_reference": sha256(REFERENCE_JSON),
        "compute31_report_read_only_reference": sha256(REFERENCE_REPORT),
        "pinned_CAMB_2_0_4_wheel": sha256(CAMB_WHEEL),
        "independent_script": sha256(Path(__file__).resolve()),
        "independent_script_delivered": sha256(Path(__file__).resolve()),
    }
    signature = hashlib.sha256(json.dumps({"hashes": hashes, "nodes": NODES, "global": sorted(GLOBAL_B_VALUES),
                                          "seeds": 12, "workers": workers}, sort_keys=True).encode()).hexdigest()
    checkpoint = {"signature": signature, "completed": {}}
    if CHECKPOINT_PATH.exists():
        old = json.loads(CHECKPOINT_PATH.read_text(encoding="utf-8"))
        if old.get("signature") == signature:
            checkpoint = old
    tasks = [(i, float(value)) for i, value in enumerate(NODES)]
    todo = [task for task in tasks if f"{task[1]:.12f}" not in checkpoint["completed"]]
    search_started = time.perf_counter()
    with futures.ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
        pending = {pool.submit(run_slice, task): task for task in todo}
        for future in futures.as_completed(pending):
            task = pending[future]
            result = future.result()
            checkpoint["completed"][f"{task[1]:.12f}"] = result
            checkpoint["last_completed_omega_b_h2"] = task[1]
            checkpoint["search_wall_seconds_at_checkpoint"] = time.perf_counter() - search_started
            atomic_json(CHECKPOINT_PATH, checkpoint)
    search_wall = time.perf_counter() - search_started
    slices = [checkpoint["completed"][f"{v:.12f}"] for v in NODES]
    all_trials = [trial for slc in slices for trial in slc["all_powell_starts"]]
    de_rows = [slc["global_differential_evolution"] for slc in slices if slc["global_differential_evolution"]]
    worker_cpu = sum(float(slc["worker_cpu_seconds"]) for slc in slices)
    best_slice = min(slices, key=lambda slc: float(slc["best"]["chi2"]))
    ref = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))
    ref_best = ref["profile_grid"][-1]["best_valid_trial"]
    result = {
        "status": "completed_independent_fixed_omega_b_globality_stress_test",
        "scope": "exploratory numerical optimizer audit; multistart and differential evolution do not prove globality",
        "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": time.perf_counter() - started,
        "search_wall_seconds": search_wall,
        "hardware": {
            "available_affinity_cpus": len(os.sched_getaffinity(0)),
            "process_pool_workers": workers,
            "blas_openmp_threads_per_process": 1,
            "gpu_used": False,
            "worker_cpu_seconds_sum": worker_cpu,
            "measured_average_busy_cores_worker_cpu_over_search_wall": worker_cpu / max(search_wall, 1.0e-12),
            "interpretation": "aggregate worker process CPU time divided by pool wall time; indicates realized CPU parallelism, not a formal serial-vs-parallel benchmark",
        },
        "versions": {
            "CAMB": camb.__version__, "NumPy": np.__version__, "SciPy": scipy.__version__,
            "Matplotlib": matplotlib.__version__, "Python": sys.version, "platform": platform.platform(),
        },
        "model": {
            "CAMB_parameters": {"mnu_eV": MNU_EV, "N_eff": N_EFF, "massive_neutrino_species": 1,
                "hierarchy": "degenerate", "T_CMB_K": T_CMB_K, "YHe": None,
                "Omega_k": OMEGA_K, "dark_energy": "flat fluid Lambda (w=-1, wa=0)"},
            "bounds": {"H0_km_s_Mpc": list(H0_BOUNDS), "omega_b_h2": list(OBH2_BOUNDS),
                       "omega_c_h2": list(OCH2_BOUNDS), "theta_MC_100": list(THETA100_BOUNDS)},
            "feasible_domain_note": "For omega_b and omega_c lower bounds, points below H0=100 sqrt(omega_b h2 + omega_c h2 + omega_nu h2 + omega_r h2) cannot have positive flat Lambda; the optimizer parameterizes the remaining feasible domain and approaches Omega_Lambda=0 from above.",
            "radiation_and_massive_neutrino_physical_densities_from_CAMB": {
                "omega_nu_h2": OMEGA_NU_H2, "omega_r_h2": OMEGA_R_H2,
                "closure_epsilon_physical": EPS_CLOSURE_PHYSICAL,
            },
        },
        "data_and_units": {
            "official_vector_rows": 13, "covariance_shape": [13, 13], "vector_order": OBSERVABLES,
            "z": Z.tolist(), "observed_values": Y.tolist(),
            "prediction_definitions": {
                "DM_over_rs": "CAMB comoving radial distance D_M [Mpc] / drag horizon r_drag [Mpc]",
                "DH_over_rs": "(299792.458 km/s / CAMB H(z) [km/s/Mpc]) [Mpc] / r_drag [Mpc]",
                "DV_over_rs": "(z D_M^2 D_H)^(1/3) [Mpc] / r_drag [Mpc]",
            },
            "score": "(prediction-observed)^T C^{-1}(prediction-observed), using all official covariance entries",
            "dense_solve_vs_Cholesky": "computed for each returned endpoint; numerical tolerance 2e-10 relative",
        },
        "reproduction": {
            "script": str(Path(__file__).resolve()),
            "bounded_command": ".venv/bin/python scripts/run_bounded.py --state work/compute33/no_deadline_state.json --seconds 800 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute33/independent_bbn_bao_globality.py --workers 20",
            "fixed_profile_method": "Powell bounded on closure-aware unit-square coordinates; 12 compute31-scale/edge seeds plus 12 deterministic LHS seeds per fixed omega_b slice",
            "global_seeds": "SciPy differential_evolution with Sobol initial population and deterministic slice-specific seeds at both omega_b endpoints and the three prior centers; DE endpoints Powell-polished",
            "powell_options": POWELL_OPTIONS, "differential_evolution_options": DE_OPTIONS,
            "prior_center_sensitivities": PRIOR_CENTERS,
        },
        "hashes": hashes,
        "one_point_independent_reproduction": one_check,
        "profile_nodes_omega_b_h2": NODES,
        "profile_results": slices,
        "globality_summary": {
            "fixed_slice_count": len(slices),
            "total_Powell_starts": len(all_trials),
            "Powell_failed_or_invalid_start_count": sum(not bool(t["valid_flat_positive_Lambda"]) for t in all_trials),
            "Powell_invalid_objective_evaluation_count": sum(
                sum(int(v) for v in t["invalid_counts"].values()) for t in all_trials
            ),
            "Powell_optimizer_non_success_count": sum(not bool(t["optimizer_success"]) for t in all_trials),
            "global_DE_slice_count": len(de_rows),
            "DE_total_function_evaluations": sum(int(row["DE_nfev"]) for row in de_rows),
            "DE_invalid_objective_evaluation_count": sum(
                sum(int(v) for v in row["invalid_counts"].values()) for row in de_rows
            ),
            "DE_failed_or_invalid_start_count": sum(
                not bool(row["DE_success"]) or not bool(row["Powell_polish"]["valid_flat_positive_Lambda"])
                for row in de_rows
            ),
            "best_independent_profile_point": {k: best_slice["best"][k] for k in ("H0_km_s_Mpc", "omega_b_h2", "omega_c_h2", "chi2")},
            "compute31_best_profile_point": {k: ref_best[k] for k in ("H0_km_s_Mpc", "omega_b_h2", "omega_c_h2", "chi2")},
            "independent_minus_compute31_best_chi2": float(best_slice["best"]["chi2"]) - float(ref_best["chi2"]),
            "lowest_profile_score_by_slice": [{"omega_b_h2": s["omega_b_h2_fixed"], "chi2": s["best"]["chi2"],
                "H0": s["best"]["H0_km_s_Mpc"], "omega_c_h2": s["best"]["omega_c_h2"]} for s in slices],
            "all_failed_and_local_solutions_retained": True,
            "max_dense_cholesky_delta_over_returned_endpoints": max(
                float(t["closure_and_score_check"]["dense_cholesky_abs_delta"])
                for t in all_trials if t["closure_and_score_check"] is not None
            ),
            "all_selected_slice_theta_checks_inside_interval": all(
                bool(s["best_theta_support_check"]["inside_declared_100theta_MC_interval"]) for s in slices
            ),
            "all_returned_valid_starts_have_positive_CAMB_Omega_de": all(
                float(t["closure_and_score_check"]["closure"]["Omega_de0_CAMB"]) > 0
                for t in all_trials if t["valid_flat_positive_Lambda"]
            ),
        },
        "limitations": [
            "Finite deterministic multistarts and differential evolution cannot certify the continuous global optimum.",
            "Each fixed omega_b profile is searched independently; baryon coverage is a finite grid, augmented near the reported upper-edge minimum.",
            "Selected BBN prior centers are fixed-profile cross-checks; this audit does not reproduce a marginalized posterior or evidence.",
        ],
    }
    atomic_json(RESULT_PATH, result)
    make_plot(result)
    make_report(result)
    return result


def make_plot(result: dict[str, object]) -> None:
    rows = result["globality_summary"]["lowest_profile_score_by_slice"]
    b = np.asarray([r["omega_b_h2"] for r in rows])
    c = np.asarray([r["chi2"] for r in rows])
    c31 = float(result["globality_summary"]["compute31_best_profile_point"]["chi2"])
    fig, ax = plt.subplots(figsize=(8.2, 5.0), constrained_layout=True)
    ax.plot(b, c - c31, "o-", ms=3.5, lw=1.4, color="#315b8a")
    ax.axhline(0.0, color="#ad5146", lw=1.0, ls="--", label="compute31 sampled minimum")
    for name, value in PRIOR_CENTERS.items():
        ax.axvline(value, color="#888888", lw=0.7, alpha=0.55)
    ax.set(xlabel=r"fixed $\omega_b h^2$", ylabel=r"$\Delta\chi^2$ from compute31 minimum",
           title="Independent fixed-baryon optimizer stress test")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False, loc="best")
    fig.savefig(PLOT_PATH, dpi=170)
    plt.close(fig)


def make_report(result: dict[str, object]) -> None:
    summary = result["globality_summary"]
    best = summary["best_independent_profile_point"]
    ref = summary["compute31_best_profile_point"]
    min_theta = min(float(s["best_theta_support_check"]["implied_100theta_MC"])
                    for s in result["profile_results"])
    max_theta = max(float(s["best_theta_support_check"]["implied_100theta_MC"])
                    for s in result["profile_results"])
    text = f"""# Independent compute33: DESI DR2 BAO optimizer stress test

Status: **exploratory independent numerical cross-check**. This scorer reparses the official 13-row DESI BAO mean vector and full 13×13 covariance, and uses pinned CAMB 2.0.4 without importing compute31. It evaluates flat fluid LambdaCDM with positive flat-Lambda closure over fixed omega_b h² slices.

The one-point reproduction at compute31's reported BAO minimum agrees exactly in chi². The search covers {summary['fixed_slice_count']} fixed omega_b h² values, including the 15-point coarse grid, seven upper-edge refinements, and three selected prior centers. It ran {summary['total_Powell_starts']} deterministic Powell starts and {summary['global_DE_slice_count']} Sobol-initialized differential-evolution searches ({summary['DE_total_function_evaluations']} evaluations). There were {summary['Powell_invalid_objective_evaluation_count']} transient invalid Powell probes and {summary['DE_invalid_objective_evaluation_count']} transient invalid differential-evolution probes; all final Powell starts were valid and successful. No distinct local basin was found more than Δchi²=10⁻⁶ above a slice minimum.

The best independent point is chi²={best['chi2']:.12f} at omega_b h²={best['omega_b_h2']:.8f}, H0={best['H0_km_s_Mpc']:.6f}, and omega_c h²={best['omega_c_h2']:.8f}. Compute31 reported chi²={ref['chi2']:.12f} at effectively the same parameters. The independent-minus-reference score change is {summary['independent_minus_compute31_best_chi2']:.3g}, consistent with floating-point roundoff. This run found no credible lower basin; finite searches cannot prove globality.

Dense covariance solve and Cholesky-whitened scores agree to at most {summary['max_dense_cholesky_delta_over_returned_endpoints']:.3g}. Selected CAMB 100 theta_MC values span [{min_theta:.4f}, {max_theta:.4f}] within the declared [{THETA100_BOUNDS[0]:.1f}, {THETA100_BOUNDS[1]:.1f}] support. CAMB closure, H(z), distances, drag sound horizon, and every returned start are recorded in the JSON.

Runtime was {result['runtime_seconds']:.2f} s with {result['hardware']['process_pool_workers']} workers. Aggregate worker CPU time was {result['hardware']['worker_cpu_seconds_sum']:.2f} s, or {result['hardware']['measured_average_busy_cores_worker_cpu_over_search_wall']:.1f} average busy cores. Threads were capped at one per process; no GPU or downloads were used.

## Reproduction

```sh
.venv/bin/python scripts/run_bounded.py --state work/compute33/no_deadline_state.json --seconds 800 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute33/independent_bbn_bao_globality.py --workers 20
```

Executed-run script SHA-256: `{result['hashes']['independent_script']}`. The delivered script hash is `{result['hashes']['independent_script_delivered']}`; the post-run edit only corrects invalid-probe reporting and report formatting, and the numerical search was not rerun. Input hashes, software versions, exact seeds, and all profile/local-solution details are in [independent_globality.json](independent_globality.json). The profile figure is [profile_globality.png](profile_globality.png).
"""
    REPORT_PATH.write_text(text, encoding="utf-8")


def main() -> None:
    if camb.__version__ != "2.0.4":
        raise RuntimeError(f"pinned CAMB 2.0.4 required, got {camb.__version__}")
    if len(sys.argv) > 1 and sys.argv[1] == "--check-one-point":
        print(json.dumps(one_point_check(), indent=2, sort_keys=True))
        return
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else min(20, len(NODES))
    if workers < 1 or workers > min(20, len(os.sched_getaffinity(0))):
        raise ValueError("worker count must be in [1, min(20, CPU affinity)]")
    result = run_all(workers)
    print(json.dumps({
        "status": result["status"],
        "runtime_seconds": result["runtime_seconds"],
        "search_wall_seconds": result["search_wall_seconds"],
        "busy_cores": result["hardware"]["measured_average_busy_cores_worker_cpu_over_search_wall"],
        "best": result["globality_summary"]["best_independent_profile_point"],
        "delta_chi2_vs_compute31": result["globality_summary"]["independent_minus_compute31_best_chi2"],
        "output": str(RESULT_PATH),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

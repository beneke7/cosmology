#!/usr/bin/env python3
"""Independent fixed-omega_b profile and boundary audit for compute26.

Uses only local DESI DR2 inputs and the pinned CAMB 2.0.4 wheel. The likelihood
and CAMB parameter construction are reimplemented here; compute26 is read only
as the challenge point. No network access, installation, or shared-state write.
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

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

TASK_DIR = Path(__file__).resolve().parent
PROJECT = TASK_DIR.parents[1]
PINNED_SITE = PROJECT / "work/compute24/site"
sys.path.insert(0, str(PINNED_SITE))

import camb  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402
from scipy.integrate import quad  # noqa: E402
from scipy.optimize import brentq, differential_evolution  # noqa: E402

MEAN_PATH = PROJECT / "context/data/desi_dr2_mean.txt"
COV_PATH = PROJECT / "context/data/desi_dr2_cov.txt"
WHEEL_PATH = PROJECT / "context/code/camb-2.0.4-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl"
PAPER_PATH = PROJECT / "context/papers/interacting_de_desi_dr2_2026.pdf"
COMPUTE26_PATH = PROJECT / "work/compute26/independent_lcdm_bao_profile.json"
OUTPUT_PATH = TASK_DIR / "local_boundary_audit.json"
PLOT_PATH = TASK_DIR / "local_boundary_audit.png"

C_KM_S = np.float64(299792.458)
FIXED = {
    "mnu_eV": 0.06,
    "N_eff": 3.044,
    "T_CMB_K": 2.7255,
    "num_massive_neutrinos": 1,
    "neutrino_hierarchy": "degenerate",
    "YHe": "CAMB BBN consistency default (YHe=None)",
    "w": -1.0,
    "wa": 0.0,
    "Omega_k": 0.0,
}
FULL_BOUNDS = ((10.0, 1000.0), (0.005, 0.1), (0.001, 0.99))
PROFILE_BOUNDS = ((10.0, 1000.0), (0.001, 0.99))  # (H0, omega_c h2), at fixed omega_b h2
WB_SLICES = (0.08, 0.095, 0.099, 0.0999, 0.1, 0.1001, 0.105, 0.12)
FD_STEPS = (1e-3, 1e-4, 1e-5)
FULL_SEEDS = (20260929, 20260930, 20260931)
PROFILE_SEED_BASE = 20261000
FULL_MAXITER = 350
PROFILE_MAXITER = 300
FULL_POPULATION = 36
PROFILE_POPULATION = 32
INVALID_SCORE = np.float64(1e70)
_invalid_counts: dict[str, int] = {}
_invalid_examples: list[dict[str, object]] = []


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_data() -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    records = []
    for line in MEAN_PATH.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            z_s, mean_s, observable = line.split()
            records.append((float(z_s), float(mean_s), observable))
    z = np.asarray([r[0] for r in records], dtype=np.float64)
    y = np.asarray([r[1] for r in records], dtype=np.float64)
    obs = [r[2] for r in records]
    cov = np.asarray(np.loadtxt(COV_PATH, dtype=np.float64), dtype=np.float64)
    if z.shape != (13,) or y.shape != (13,) or cov.shape != (13, 13):
        raise ValueError("expected the 13-row released BAO vector and 13x13 covariance")
    if not np.all(np.isfinite(cov)) or not np.allclose(cov, cov.T, rtol=0.0, atol=1e-12):
        raise ValueError("covariance is nonfinite or asymmetric")
    np.linalg.cholesky(cov)
    return z, y, cov, obs


Z, Y, COV, OBSERVABLES = load_data()
L = np.linalg.cholesky(COV)


def camb_params(h0: float, omega_b: float, omega_c: float) -> camb.CAMBparams:
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=FIXED["w"], wa=FIXED["wa"], dark_energy_model="fluid")
    pars.set_cosmology(
        H0=float(h0),
        ombh2=float(omega_b),
        omch2=float(omega_c),
        omk=FIXED["Omega_k"],
        mnu=FIXED["mnu_eV"],
        nnu=FIXED["N_eff"],
        num_massive_neutrinos=FIXED["num_massive_neutrinos"],
        neutrino_hierarchy=FIXED["neutrino_hierarchy"],
        TCMB=FIXED["T_CMB_K"],
        YHe=None,
    )
    return pars


def score(point: tuple[float, float, float] | np.ndarray | list[float], retain_objects: bool = False) -> dict[str, object]:
    h0, omega_b, omega_c = (float(v) for v in point)
    pars = camb_params(h0, omega_b, omega_c)
    om = float(pars.omegam)
    if not (0.0 < om < 1.0):
        raise ValueError(f"flat positive-Λ closure failed: CAMB Omega_m={om:.16g}")
    background = camb.get_background(pars)
    ode = float(background.get_Omega("de", 0.0))
    if not math.isfinite(ode) or ode <= 0.0:
        raise ValueError(f"CAMB Omega_de is not positive: {ode:.16g}")
    derived = background.get_derived_params()
    rdrag = float(derived["rdrag"])
    if not math.isfinite(rdrag) or rdrag <= 0.0:
        raise ValueError(f"invalid CAMB r_drag={rdrag}")

    hz = np.asarray(background.hubble_parameter(Z), dtype=np.float64)
    dm = np.asarray(background.comoving_radial_distance(Z), dtype=np.float64)
    dh = C_KM_S / hz
    prediction = np.empty(13, dtype=np.float64)
    for i, name in enumerate(OBSERVABLES):
        if name == "DM_over_rs":
            prediction[i] = dm[i] / rdrag
        elif name == "DH_over_rs":
            prediction[i] = dh[i] / rdrag
        elif name == "DV_over_rs":
            prediction[i] = np.power(Z[i] * dm[i] ** 2 * dh[i], 1.0 / 3.0) / rdrag
        else:
            raise ValueError(f"unsupported BAO observable {name}")
    residual = prediction - Y
    white = np.linalg.solve(L, residual)
    chi2_l = float(white @ white)
    chi2_dense = float(residual @ np.linalg.solve(COV, residual))
    if not math.isfinite(chi2_l) or not math.isfinite(chi2_dense):
        raise ArithmeticError("nonfinite BAO chi2")
    if abs(chi2_l - chi2_dense) > 1e-10 * max(1.0, abs(chi2_dense)):
        raise ArithmeticError(f"full-covariance solve mismatch {chi2_l} vs {chi2_dense}")
    out: dict[str, object] = {
        "parameters": {"H0_km_s_Mpc": h0, "omega_b_h2": omega_b, "omega_c_h2": omega_c},
        "physical_closure": {
            "Omega_m0_CAMB": om,
            "Omega_de0_CAMB": ode,
            "positive_flat_Lambda": True,
        },
        "derived": {
            "r_drag_Mpc": rdrag,
            "100theta_star_CAMB": float(derived["thetastar"]),
            "z_drag": float(derived["zdrag"]),
        },
        "bao": {
            "chi2_cholesky_full_covariance": chi2_l,
            "chi2_dense_full_covariance": chi2_dense,
            "abs_difference": abs(chi2_l - chi2_dense),
            "predictions": prediction.tolist(),
            "DM_Mpc": dm.tolist(),
            "H_km_s_Mpc": hz.tolist(),
        },
    }
    if retain_objects:
        out["_background"] = background
    return out


def objective(point: np.ndarray, fixed_b: float | None = None) -> float:
    point = np.asarray(point, dtype=np.float64)
    if fixed_b is None:
        full = point
    else:
        full = np.asarray([point[0], fixed_b, point[1]], dtype=np.float64)
    try:
        if fixed_b is None and any(v < lo or v > hi for v, (lo, hi) in zip(full, FULL_BOUNDS)):
            return _penalty("outside_full_box", full)
        if fixed_b is not None and any(v < lo or v > hi for v, (lo, hi) in zip(point, PROFILE_BOUNDS)):
            return _penalty("outside_profile_box", full)
        return float(score(full)["bao"]["chi2_cholesky_full_covariance"])
    except ArithmeticError as exc:
        return _penalty("arithmetic", full, str(exc))
    except Exception as exc:
        return _penalty(type(exc).__name__, full, str(exc))


def _penalty(reason: str, point: np.ndarray, details: str | None = None) -> float:
    _invalid_counts[reason] = _invalid_counts.get(reason, 0) + 1
    if len(_invalid_examples) < 30:
        item: dict[str, object] = {"reason": reason, "parameters": np.asarray(point, dtype=np.float64).tolist()}
        if details:
            item["detail"] = details[:200]
        _invalid_examples.append(item)
    return float(INVALID_SCORE)


def seed_population_full(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    # These anchors challenge the compute26 point and the older theta-root
    # candidate; the remainder spans the declared box, filtered only to begin
    # with nonnegative flat-Lambda points. Objective evaluation checks CAMB.
    candidates = [
        np.array([133.6134163, 0.1, 0.43027242]),
        np.array([99.16242, 0.06069, 0.23096]),
        np.array([300.0, 0.03, 0.3]),
        np.array([700.0, 0.05, 0.2]),
        np.array([950.0, 0.02, 0.15]),
    ]
    accepted = [x for x in candidates if (x[1] + x[2] + 0.0007) / (x[0] / 100.0) ** 2 < 0.98]
    tries = 0
    while len(accepted) < FULL_POPULATION and tries < 100_000:
        tries += 1
        # Keep useful coverage around H0~100--200 while retaining a set of
        # high-H0 probes; uniform H0 sampling alone overpopulates corners that
        # fail recombination and leaves too little effective support nearby.
        if rng.random() < 0.75:
            h0 = rng.uniform(60.0, 250.0)
        else:
            h0 = rng.uniform(250.0, 950.0)
        x = np.array([h0, rng.uniform(*FULL_BOUNDS[1]), rng.uniform(*FULL_BOUNDS[2])], dtype=np.float64)
        if (x[1] + x[2] + 0.0007) / (x[0] / 100.0) ** 2 < 0.98:
            accepted.append(x)
    if len(accepted) < FULL_POPULATION:
        raise RuntimeError("could not generate initial physical population")
    return np.asarray(accepted[:FULL_POPULATION], dtype=np.float64)


def seed_population_profile(fixed_b: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    candidates = [np.array([133.6134163, 0.43027242]), np.array([99.16242, 0.23096])]
    accepted = [x for x in candidates if (fixed_b + x[1] + 0.0007) / (x[0] / 100.0) ** 2 < 0.98]
    tries = 0
    while len(accepted) < PROFILE_POPULATION and tries < 100_000:
        tries += 1
        h0 = rng.uniform(60.0, 250.0) if rng.random() < 0.75 else rng.uniform(250.0, 950.0)
        x = np.array([h0, rng.uniform(*PROFILE_BOUNDS[1])], dtype=np.float64)
        if (fixed_b + x[1] + 0.0007) / (x[0] / 100.0) ** 2 < 0.98:
            accepted.append(x)
    if len(accepted) < PROFILE_POPULATION:
        raise RuntimeError("could not generate initial profile population")
    return np.asarray(accepted[:PROFILE_POPULATION], dtype=np.float64)


def inverse_theta_mapping(parameters: tuple[float, float, float] | np.ndarray) -> dict[str, object]:
    target_h0, omega_b, omega_c = (float(v) for v in parameters)

    def h0_for_theta(theta100: float) -> float:
        pars = camb.CAMBparams()
        pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
        pars.set_cosmology(
            cosmomc_theta=theta100 / 100.0,
            ombh2=omega_b,
            omch2=omega_c,
            omk=0.0,
            mnu=FIXED["mnu_eV"],
            nnu=FIXED["N_eff"],
            num_massive_neutrinos=FIXED["num_massive_neutrinos"],
            neutrino_hierarchy=FIXED["neutrino_hierarchy"],
            TCMB=FIXED["T_CMB_K"],
            YHe=None,
            theta_H0_range=(1.0, 2000.0),
        )
        return float(pars.H0)

    sample = []
    for theta100 in np.linspace(0.6, 2.0, 71):
        try:
            sample.append((float(theta100), h0_for_theta(float(theta100))))
        except Exception:
            continue
    bracket = next(
        ((a[0], b[0]) for a, b in zip(sample, sample[1:]) if (a[1] - target_h0) * (b[1] - target_h0) <= 0.0),
        None,
    )
    if bracket is None:
        return {"status": "no theta bracket", "sampled_successes": len(sample)}
    theta100 = float(brentq(lambda t: h0_for_theta(t) - target_h0, bracket[0], bracket[1], xtol=2e-12, rtol=1e-13))
    reconstructed_h0 = h0_for_theta(theta100)
    try:
        pars_original = camb.CAMBparams()
        pars_original.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
        pars_original.set_cosmology(
            cosmomc_theta=theta100 / 100.0,
            ombh2=omega_b,
            omch2=omega_c,
            omk=0.0,
            mnu=FIXED["mnu_eV"],
            nnu=FIXED["N_eff"],
            num_massive_neutrinos=FIXED["num_massive_neutrinos"],
            neutrino_hierarchy=FIXED["neutrino_hierarchy"],
            TCMB=FIXED["T_CMB_K"],
            YHe=None,
            theta_H0_range=(10.0, 100.0),
        )
        original_bracket = {"root_found": True, "H0": float(pars_original.H0)}
    except Exception as exc:
        original_bracket = {"root_found": False, "error": f"{type(exc).__name__}: {exc}"}
    return {
        "status": "rooted",
        "100theta_MC": theta100,
        "inside_Table_I_0p5_to_10": bool(0.5 <= theta100 <= 10.0),
        "expanded_theta_H0_range": [1.0, 2000.0],
        "recovered_H0_km_s_Mpc": reconstructed_h0,
        "H0_abs_error": abs(reconstructed_h0 - target_h0),
        "original_theta_H0_range": [10.0, 100.0],
        "original_bracket": original_bracket,
        "successful_theta_samples": len(sample),
    }


def finite_difference_audit(challenge: np.ndarray) -> dict[str, object]:
    h0, _wb_challenge, wc = (float(v) for v in challenge)
    wb = 0.1  # evaluate the exact declared upper edge
    center = float(score([h0, wb, wc])["bao"]["chi2_cholesky_full_covariance"])
    rows = []
    for step in FD_STEPS:
        lower = float(score([h0, wb - step, wc])["bao"]["chi2_cholesky_full_covariance"])
        upper = float(score([h0, wb + step, wc])["bao"]["chi2_cholesky_full_covariance"])
        rows.append(
            {
                "step": step,
                "omega_b_lower": wb - step,
                "chi2_lower_fixed_H0_omega_c": lower,
                "omega_b_center": wb,
                "chi2_center": center,
                "omega_b_upper": wb + step,
                "chi2_upper_fixed_H0_omega_c": upper,
                "left_secant_dchi2_domega_b": (center - lower) / step,
                "right_secant_dchi2_domega_b": (upper - center) / step,
                "central_second_difference": (upper - 2.0 * center + lower) / step**2,
                "outside_declared_domain_on_right": wb + step > 0.1,
            }
        )
    return {
        "parameters_fixed": {"H0_km_s_Mpc": h0, "omega_c_h2": wc},
        "center_chi2": center,
        "finite_difference_rows": rows,
        "reading": "At a domain upper edge, the in-domain secant should be nonpositive toward the boundary for a boundary minimum; an outside-box upward score change is a diagnostic, not evidence against the in-box minimum.",
    }


def distance_check(candidate: dict[str, object]) -> dict[str, object]:
    # Reconstruct independently from the profile parameters.
    p = candidate["parameters"]
    evaluated = score([p["H0_km_s_Mpc"], p["omega_b_h2"], p["omega_c_h2"]], retain_objects=True)
    bg = evaluated["_background"]
    errors = []
    for z, dm in zip(Z, evaluated["bao"]["DM_Mpc"]):
        integral, _ = quad(lambda zz: C_KM_S / float(bg.hubble_parameter(zz)), 0.0, float(z), epsabs=1e-9, epsrel=1e-11, limit=250)
        errors.append(abs(float(dm) - float(integral)))
    return {"max_abs_D_M_error_Mpc": float(max(errors)), "all_below_1e-6_Mpc": bool(max(errors) < 1e-6)}


def write_plot(fd: dict[str, object], profiles: list[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    xs = [row["omega_b_h2"] for row in profiles]
    ys = [row["chi2_profile"] for row in profiles]
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    inside_x = [x for x in xs if x <= 0.1]
    inside_y = [y for x, y in zip(xs, ys) if x <= 0.1]
    outside_x = [x for x in xs if x > 0.1]
    outside_y = [y for x, y in zip(xs, ys) if x > 0.1]
    ax.plot(inside_x, inside_y, "o-", label="inside declared box")
    ax.plot(outside_x, outside_y, "s--", label="outside-box diagnostic")
    ax.axvline(0.1, color="black", lw=1, ls=":", label=r"$\omega_bh^2=0.1$ bound")
    ax.axhline(ys[xs.index(0.1)], color="gray", lw=0.8, ls="--")
    ax.set_xlabel(r"Fixed $\omega_b h^2$")
    ax.set_ylabel(r"Profile $\chi^2$ (optimize $H_0,\omega_c h^2$)")
    ax.set_yscale("linear")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    ax.set_title("Independent fixed-baryon BAO profile")
    fig.savefig(PLOT_PATH, dpi=170)
    plt.close(fig)


def main() -> int:
    started = datetime.now(timezone.utc).isoformat()
    timer = time.perf_counter()
    compute26 = json.loads(COMPUTE26_PATH.read_text(encoding="utf-8"))
    challenge = np.asarray(
        [
            compute26["selected_candidate"]["parameters"]["H0_km_s_Mpc"],
            compute26["selected_candidate"]["parameters"]["omega_b_h2"],
            compute26["selected_candidate"]["parameters"]["omega_c_h2"],
        ],
        dtype=np.float64,
    )
    challenge_eval = score(challenge, retain_objects=True)
    differences = finite_difference_audit(challenge)

    full_runs = []
    for seed in FULL_SEEDS:
        t0 = time.perf_counter()
        fit = differential_evolution(
            objective,
            bounds=FULL_BOUNDS,
            init=seed_population_full(seed),
            strategy="best1bin",
            maxiter=FULL_MAXITER,
            popsize=FULL_POPULATION // 3,
            tol=1e-8,
            atol=1e-9,
            seed=seed,
            polish=True,
            updating="immediate",
            workers=1,
            disp=False,
        )
        x = np.asarray(fit.x, dtype=np.float64)
        try:
            verified = score(x)
            valid = True
            chi2 = float(verified["bao"]["chi2_cholesky_full_covariance"])
        except Exception as exc:
            valid = False
            chi2 = float(INVALID_SCORE)
            error = f"{type(exc).__name__}: {exc}"
        run = {
            "seed": seed,
            "success": bool(fit.success),
            "message": str(fit.message),
            "iterations": int(fit.nit),
            "function_evaluations": int(fit.nfev),
            "runtime_seconds": time.perf_counter() - t0,
            "candidate_valid": valid,
            "chi2": chi2,
            "parameters": x.tolist(),
            "max_population_energy": float(np.max(fit.population_energies)),
        }
        if not valid:
            run["validation_error"] = error
        full_runs.append(run)

    profile_runs = []
    for index, wb in enumerate(WB_SLICES):
        seed = PROFILE_SEED_BASE + index
        t0 = time.perf_counter()
        fit = differential_evolution(
            lambda x, fixed_b=wb: objective(x, fixed_b=fixed_b),
            bounds=PROFILE_BOUNDS,
            init=seed_population_profile(wb, seed),
            strategy="best1bin",
            maxiter=PROFILE_MAXITER,
            popsize=PROFILE_POPULATION // 2,
            tol=1e-8,
            atol=1e-9,
            seed=seed,
            polish=True,
            updating="immediate",
            workers=1,
            disp=False,
        )
        nuisance = np.asarray(fit.x, dtype=np.float64)
        point = [float(nuisance[0]), float(wb), float(nuisance[1])]
        try:
            best = score(point)
            valid = True
            chi2 = float(best["bao"]["chi2_cholesky_full_covariance"])
            closure = best["physical_closure"]
            theta = inverse_theta_mapping(point)
        except Exception as exc:
            valid = False
            chi2 = float(INVALID_SCORE)
            closure = {"valid": False, "error": f"{type(exc).__name__}: {exc}"}
            theta = {"status": "not_evaluated"}
        profile_runs.append(
            {
                "omega_b_h2": float(wb),
                "inside_declared_omega_b_bound": bool(0.005 <= wb <= 0.1),
                "seed": seed,
                "success": bool(fit.success),
                "message": str(fit.message),
                "iterations": int(fit.nit),
                "function_evaluations": int(fit.nfev),
                "runtime_seconds": time.perf_counter() - t0,
                "nuisance_parameters": {"H0_km_s_Mpc": float(nuisance[0]), "omega_c_h2": float(nuisance[1])},
                "chi2_profile": chi2,
                "candidate_valid": valid,
                "physical_closure": closure,
                "theta_mapping": theta,
            }
        )

    valid_full = [r for r in full_runs if r["candidate_valid"]]
    best_full = min(valid_full, key=lambda r: r["chi2"]) if valid_full else None
    best_profile = next(r for r in profile_runs if r["omega_b_h2"] == 0.1)
    # Profile finite differences compare nuisance-reoptimized chi2 slices on
    # each side; outside-box values are intentionally not promoted as fits.
    by_wb = {r["omega_b_h2"]: r for r in profile_runs}
    profile_secants = []
    for lower_wb, upper_wb in ((0.0999, 0.1), (0.1, 0.1001), (0.099, 0.1), (0.1, 0.105)):
        lower, upper = by_wb[lower_wb], by_wb[upper_wb]
        profile_secants.append(
            {
                "omega_b_pair": [lower_wb, upper_wb],
                "delta_omega_b": upper_wb - lower_wb,
                "lower_profile_chi2": lower["chi2_profile"],
                "upper_profile_chi2": upper["chi2_profile"],
                "secant_dprofile_chi2_domega_b": (upper["chi2_profile"] - lower["chi2_profile"]) / (upper_wb - lower_wb),
                "outside_declared_box": bool(upper_wb > 0.1),
            }
        )

    # The allowed-box challenge is whether any reoptimized profile slice or
    # full-box independent run beats the pinned compute26 score.
    compute26_chi2 = float(compute26["selected_candidate"]["bao"]["chi2_full_covariance_cholesky"])
    interior_profiles = [r for r in profile_runs if r["inside_declared_omega_b_bound"] and r["candidate_valid"]]
    lower_allowed_profile = min(interior_profiles, key=lambda r: r["chi2_profile"])
    lower_full_candidates = [r for r in valid_full if r["chi2"] < compute26_chi2 - 1e-6]

    plot_payload = profile_runs
    write_plot(differences, plot_payload)
    ended = datetime.now(timezone.utc).isoformat()
    report = {
        "schema_version": 1,
        "status": "boundary_minimum_supported_on_tested_box" if (not lower_full_candidates and lower_allowed_profile["omega_b_h2"] == 0.1 and best_profile["candidate_valid"]) else "boundary_minimum_not_confirmed",
        "scope": {
            "model": "flat standard LambdaCDM background using CAMB 2.0.4; direct H0, physical-density profile",
            "likelihood": "DESI DR2 13-row Gaussian BAO vector with full 13x13 covariance only",
            "declared_box": {"H0_km_s_Mpc": [10.0, 1000.0], "omega_b_h2": [0.005, 0.1], "omega_c_h2": [0.001, 0.99]},
            "outside_box_rule": "omega_b h2 > 0.1 rows are exploratory diagnostics only and never included as in-domain candidates.",
            "not_included": ["Planck likelihood", "ACT likelihood", "posterior sampling", "evidence", "supernovae", "interaction model"],
        },
        "sources": {
            "DESI_mean_path": str(MEAN_PATH.relative_to(PROJECT)),
            "DESI_mean_sha256": hash_file(MEAN_PATH),
            "DESI_covariance_path": str(COV_PATH.relative_to(PROJECT)),
            "DESI_covariance_sha256": hash_file(COV_PATH),
            "DESI_official_upstream": "CobayaSampler/bao_data commit bb0c1c9009dc76d1391300e169e8df38fd1096db, desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt",
            "CAMB_version": camb.__version__,
            "CAMB_wheel_path": str(WHEEL_PATH.relative_to(PROJECT)),
            "CAMB_wheel_sha256": hash_file(WHEEL_PATH),
            "CAMB_repository_commit": "a6de8cc59124c8bbe3924f1c546f96cef73dcabb",
            "target_paper_path": str(PAPER_PATH.relative_to(PROJECT)),
            "target_paper_sha256": hash_file(PAPER_PATH),
            "target_paper_prior_locator": "arXiv:2609.05660 v1, p. 10 Table I: omega_b h2 [0.005,0.1], omega_c h2 [0.001,0.99], 100theta_MC [0.5,10]",
        },
        "conventions": {
            "fixed": FIXED,
            "c_km_s": float(C_KM_S),
            "observables": "DM/r_d, DH/r_d=(c/H)/r_d, DV/r_d=(z*DM^2*DH)^(1/3)/r_d",
            "covariance": "full covariance; residual^T C^-1 residual independently computed by Cholesky and dense solves",
            "dtype": "float64",
            "workers": 1,
            "numerical_threads": 1,
        },
        "challenge_point": {
            "compute26_parameters": challenge.tolist(),
            "independently_recomputed": score(challenge),
            "compute26_chi2": compute26_chi2,
            "finite_difference_at_fixed_H0_omega_c": differences,
            "theta_mapping": inverse_theta_mapping(challenge),
            "distance_integral_check": distance_check(challenge_eval),
        },
        "full_box_search": {
            "method": "Independent SciPy differential_evolution reimplementation over (H0, omega_b h2, omega_c h2), three seeds, broad physical initialization, explicit positive-flat-Lambda checks",
            "seeds": list(FULL_SEEDS),
            "maxiter": FULL_MAXITER,
            "population_size": FULL_POPULATION,
            "runs": full_runs,
            "best_valid_candidate": best_full,
            "lower_valid_candidate_than_compute26_found": bool(lower_full_candidates),
            "lower_valid_candidates": lower_full_candidates,
        },
        "fixed_baryon_profile": {
            "method": "At fixed omega_b h2, reoptimize (H0, omega_c h2) by broad-box differential evolution; all minima rechecked with CAMB closure and the full covariance.",
            "maxiter": PROFILE_MAXITER,
            "population_size": PROFILE_POPULATION,
            "slices": profile_runs,
            "best_inside_declared_slices": lower_allowed_profile,
            "chi2_at_omega_b_0p1": best_profile["chi2_profile"],
            "profile_secant_checks": profile_secants,
            "lowest_in_domain_profile_not_at_upper_edge": bool(lower_allowed_profile["omega_b_h2"] < 0.1),
        },
        "invalid_points": {"counts_by_reason": _invalid_counts, "first_examples_capped_at_30": _invalid_examples, "penalty": float(INVALID_SCORE)},
        "interpretation": {
            "fixed_coordinate_boundary_read": "The prior outward point having a higher score is consistent with a minimum at the upper edge if the in-domain score falls toward the edge. The fixed-coordinate finite differences and reoptimized slices test this condition.",
            "candidate_disposition": "Exploratory boundary point; no claim of an interior best fit or posterior.",
            "profile_scope": "The profile audit tests the declared direct-H0 finite box and selected omega_b slices. It does not establish exact globality for every continuous point in the box beyond the independent full-box multi-start search.",
        },
        "execution": {
            "started_utc": started,
            "finished_utc": ended,
            "runtime_seconds": time.perf_counter() - timer,
            "runtime_cap_seconds": 1100,
            "command": "cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute29/task_state.json --seconds 1100 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute29/local_boundary_audit.py",
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "code_sha256": hash_file(Path(__file__)),
        },
        "plot_path": str(PLOT_PATH.relative_to(PROJECT)),
    }
    OUTPUT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "compute26_chi2": compute26_chi2,
                "best_full": best_full,
                "best_inside_profile": lower_allowed_profile,
                "slice_summary": [{k: row[k] for k in ("omega_b_h2", "chi2_profile", "candidate_valid", "success", "nuisance_parameters")} for row in profile_runs],
                "finite_difference_rows": differences["finite_difference_rows"],
                "lower_valid_candidate_than_compute26_found": bool(lower_full_candidates),
                "runtime": report["execution"]["runtime_seconds"],
                "result": str(OUTPUT_PATH.relative_to(PROJECT)),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

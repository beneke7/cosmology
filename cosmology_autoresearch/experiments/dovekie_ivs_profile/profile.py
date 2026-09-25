#!/usr/bin/env python3
"""Bounded exploratory IVS profile on the exact DES-Dovekie SN likelihood.

This is a late-time background-only profile, not a joint BAO+SN fit, posterior,
or reproduction of the paper's modified-CAMB pipeline. The SN data parser and
GLS magnitude-offset algebra are reused from the audited Dovekie profile screen;
the interacting background/distance implementation below is separate from the
BAO optimizer.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from collections import Counter
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import platform
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RESULT_PATH = HERE / "result.json"
GRID_PATH = HERE / "profile_surface.npz"
PLOT_PATH = HERE / "profile_surface.png"
RECORD_PATH = HERE / "record.json"
sys.path.insert(0, str(ROOT / "work" / "inference"))
import fit_dovekie as sn  # noqa: E402

H0_GAUGE = 70.0
C_KM_S = 299792.458
U_MAX = float(np.log1p(1.14418))
INITIAL_OM_BOUNDS = (0.05, 0.60)
INITIAL_G_BOUNDS = (-3.0, 3.0)
MAX_OM = 0.95
MAX_ABS_G = 9.0
GRID_N_OM = 41
GRID_N_G = 41
N_STARTS = 16
SEED = 20260925
RTOL = 2.0e-10
ATOL = 2.0e-12
MAX_STEP_U = 0.025
INVALID = 1.0e80
_DATA: dict | None = None
_FAILURES: Counter[str] = Counter()


class DomainError(ValueError):
    """Raised when the declared late-time background is not physical/numeric."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integrate_background(omega_m: float, g: float, z_eval: np.ndarray, *,
                        rtol: float = RTOL, atol: float = ATOL,
                        max_step_u: float = MAX_STEP_U) -> dict:
    """Integrate m, x and dimensionless comoving distance in u=ln(1+z)."""
    om = float(omega_m)
    rate = float(g)
    z_eval = np.asarray(z_eval, dtype=np.float64)
    if not (np.isfinite(om) and np.isfinite(rate)):
        raise DomainError("nonfinite profile coordinate")
    if not 0.0 < om < 1.0:
        raise DomainError("present matter/vacuum densities must be positive")
    if np.any(~np.isfinite(z_eval)) or np.any(z_eval < 0.0) or np.any(z_eval > 1.14418 + 1e-12):
        raise DomainError("redshift outside the declared Dovekie interval")

    umax = float(np.log1p(max(float(np.max(z_eval)), 0.0)))

    def rhs(u: float, state: np.ndarray) -> tuple[float, float, float]:
        matter, vacuum, _distance = state
        e2 = matter + vacuum
        if not np.isfinite(e2) or e2 <= 0.0:
            raise DomainError("nonpositive/nonfinite E^2 in ODE RHS")
        e = math.sqrt(e2)
        source = rate * vacuum / e
        return (3.0 * matter + source, -source, math.exp(u) / e)

    om0 = om
    x0 = 1.0 - om
    sol = solve_ivp(
        rhs,
        (0.0, umax),
        (om0, x0, 0.0),
        method="DOP853",
        dense_output=True,
        rtol=rtol,
        atol=atol,
        max_step=max_step_u,
    )
    if not sol.success or sol.sol is None:
        raise DomainError(f"ODE solver failed: {sol.message}")

    # Dense trajectory feasibility gate over the full observed SN interval.
    # At least 257 points and delta-u <= 0.002; accepted endpoints receive a
    # second, finer pass in validate_point().
    n_check = max(257, int(math.ceil(umax / 0.002)) + 1)
    u_check = np.linspace(0.0, umax, n_check, dtype=np.float64)
    trajectory = np.asarray(sol.sol(u_check), dtype=np.float64)
    matter, vacuum, _distance = trajectory
    e2 = matter + vacuum
    if (
        not np.all(np.isfinite(trajectory))
        or np.any(matter <= 0.0)
        or np.any(vacuum <= 0.0)
        or np.any(e2 <= 0.0)
    ):
        raise DomainError("nonpositive/nonfinite m, x or E^2 on observed interval")

    # Exact baryon/CDM split ceiling. For r=m/(Omega_m0 exp(3u)),
    # r'=g*x/(E*Omega_m0 exp(3u)); with x,E>0 its minimum is at u=0 if
    # g>=0 and at u=U if g<0. This removes grid dependence from the split
    # feasibility check; the trajectory positivity check above remains sampled.
    if rate >= 0.0:
        fb_max = 1.0
    else:
        m_u = float(sol.sol(umax)[0])
        fb_max = m_u / (om0 * math.exp(3.0 * umax))
    if not np.isfinite(fb_max) or fb_max <= 0.0:
        raise DomainError("no positive constant baryon/CDM split exists")
    fb_witness = 0.5 * fb_max
    baryons = fb_witness * om0 * np.exp(3.0 * u_check)
    cdm = matter - baryons
    if np.any(baryons <= 0.0) or np.any(cdm <= 0.0):
        raise DomainError("interior baryon/CDM witness is not strictly positive")

    z_out = np.asarray(z_eval, dtype=np.float64)
    out = np.asarray(sol.sol(np.log1p(z_out)), dtype=np.float64)
    m_out, x_out, d_out = out
    e2_out = m_out + x_out
    if np.any(~np.isfinite(out)) or np.any(e2_out <= 0.0):
        raise DomainError("invalid requested-redshift solution")
    return {
        "matter": m_out,
        "vacuum": x_out,
        "distance": d_out,
        "e2": e2_out,
        "solution": sol.sol,
        "u_check": u_check,
        "trajectory": trajectory,
        "fb_max": float(fb_max),
        "fb_witness": float(fb_witness),
        "minimum_matter": float(np.min(matter)),
        "minimum_vacuum": float(np.min(vacuum)),
        "minimum_e2": float(np.min(e2)),
        "minimum_cdm_witness": float(np.min(cdm)),
    }


def prediction_and_profile(omega_m: float, g: float, *, diagnostics: bool = False,
                           ode_options: dict | None = None) -> dict:
    if _DATA is None:
        raise RuntimeError("Dovekie data not initialized")
    bg = integrate_background(omega_m, g, _DATA["z_hd"], **(ode_options or {}))
    dl = (1.0 + _DATA["z_hel"]) * (C_KM_S / H0_GAUGE) * bg["distance"]
    if np.any(~np.isfinite(dl)) or np.any(dl <= 0.0):
        raise DomainError("luminosity distance is not positive and finite")
    mu_base = 5.0 * np.log10(dl) + 25.0
    offset, chi2, residual = sn.profile_offset(
        mu_base,
        _DATA["mu_obs"],
        _DATA["precision"],
        _DATA["precision_ones"],
        _DATA["offset_information"],
    )
    if not np.isfinite(chi2) or chi2 < -1e-7:
        raise DomainError("profiled chi-square is nonfinite or negative")
    if not diagnostics:
        return {"chi2": float(chi2), "offset": float(offset)}
    diff = mu_base - _DATA["mu_obs"]
    alg = sn.independent_profile_check(
        diff,
        _DATA["precision_cholesky"],
        _DATA["offset_information"],
        offset,
        residual,
        chi2,
    )
    return {
        "chi2": float(chi2),
        "offset": float(offset),
        "normal_equation_abs_residual": float(abs(_DATA["precision_ones"] @ residual)),
        "independent_profile_algebra": alg,
        "background": bg,
        "mu_base": mu_base,
    }


def evaluate_grid_point(point: tuple[float, float]) -> tuple[float, str | None]:
    try:
        score = prediction_and_profile(float(point[0]), float(point[1]))["chi2"]
        return score, None
    except (DomainError, ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as exc:
        return INVALID, type(exc).__name__ + ": " + str(exc)


def make_grid(bounds: tuple[tuple[float, float], tuple[float, float]]) -> list[tuple[float, float]]:
    om_values = np.linspace(bounds[0][0], bounds[0][1], GRID_N_OM)
    g_values = np.linspace(bounds[1][0], bounds[1][1], GRID_N_G)
    return [(float(om), float(rate)) for rate in g_values for om in om_values]


def evaluate_grid(bounds: tuple[tuple[float, float], tuple[float, float]]) -> tuple[np.ndarray, dict]:
    points = make_grid(bounds)
    sample = points[::max(1, len(points) // 40)][:40]
    t0 = time.perf_counter()
    serial_probe = [evaluate_grid_point(p) for p in sample]
    serial_seconds = time.perf_counter() - t0
    workers = min(20, max(1, os.cpu_count() or 1), len(sample))
    parallel_seconds = None
    parallel_speedup = None
    use_parallel = False
    if workers > 1:
        t1 = time.perf_counter()
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("fork")) as pool:
            parallel_probe = list(pool.map(evaluate_grid_point, sample, chunksize=1))
        parallel_seconds = time.perf_counter() - t1
        if any(a[1] != b[1] or not np.isclose(a[0], b[0], rtol=2e-12, atol=2e-9)
               for a, b in zip(serial_probe, parallel_probe, strict=True)):
            raise ArithmeticError("serial/process grid probes disagree")
        parallel_speedup = serial_seconds / max(parallel_seconds, 1e-12)
        use_parallel = parallel_speedup >= 1.15 and len(points) >= 200

    t2 = time.perf_counter()
    if use_parallel:
        with ProcessPoolExecutor(max_workers=min(20, max(1, os.cpu_count() or 1)),
                                 mp_context=mp.get_context("fork")) as pool:
            records = list(pool.map(evaluate_grid_point, points, chunksize=4))
    else:
        records = [evaluate_grid_point(p) for p in points]
    surface_seconds = time.perf_counter() - t2
    values = np.asarray([row[0] for row in records], dtype=np.float64).reshape(GRID_N_G, GRID_N_OM)
    failures = Counter(reason for _, reason in records if reason is not None)
    return values, {
        "grid_point_count": len(points),
        "probe_point_count": len(sample),
        "probe_serial_seconds": float(serial_seconds),
        "probe_parallel_seconds": None if parallel_seconds is None else float(parallel_seconds),
        "probe_parallel_speedup": None if parallel_speedup is None else float(parallel_speedup),
        "probe_workers": int(workers),
        "grid_workers": int(min(20, max(1, os.cpu_count() or 1)) if use_parallel else 1),
        "parallel_grid_used": bool(use_parallel),
        "surface_seconds": float(surface_seconds),
        "physical_or_numeric_rejection_count": int(sum(failures.values())),
        "rejection_reasons": dict(failures),
    }


def build_starts(bounds: tuple[tuple[float, float], tuple[float, float]],
                 grid_minimum: tuple[float, float] | None = None) -> list[np.ndarray]:
    (om0, om1), (g0, g1) = bounds
    starts = [
        np.array([0.30, 0.0]),
        np.array([0.38697100769747017, -0.4666647299169934]),
        np.array([om0, 0.0]),
        np.array([om1, 0.0]),
        np.array([om0, g0]),
        np.array([om0, g1]),
        np.array([om1, g0]),
        np.array([om1, g1]),
    ]
    if grid_minimum is not None:
        starts.insert(0, np.asarray(grid_minimum, dtype=np.float64))
    rng = np.random.default_rng(SEED)
    for _ in range(max(0, N_STARTS - len(starts))):
        starts.append(np.array([rng.uniform(om0, om1), rng.uniform(g0, g1)], dtype=np.float64))
    return [np.clip(s, [om0, g0], [om1, g1]) for s in starts[:N_STARTS]]


def optimize_profile(bounds: tuple[tuple[float, float], tuple[float, float]],
                     grid_minimum: tuple[float, float] | None = None) -> list[dict]:
    records: list[dict] = []
    for index, start in enumerate(build_starts(bounds, grid_minimum)):
        failures = Counter()

        def objective(theta: np.ndarray) -> float:
            row, reason = evaluate_grid_point((float(theta[0]), float(theta[1])))
            if reason:
                failures[reason] += 1
            return row

        fit = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=list(bounds),
            options={"maxiter": 2000, "maxfun": 10000, "ftol": 1e-13, "gtol": 1e-8, "maxls": 40},
        )
        endpoint = evaluate_grid_point((float(fit.x[0]), float(fit.x[1])))
        recheck_difference = abs(float(fit.fun) - endpoint[0])
        if endpoint[1] is None and recheck_difference > 1e-6:
            raise ArithmeticError(
                f"optimizer endpoint disagrees with independent reevaluation by {recheck_difference}"
            )
        records.append({
            "start_index": index,
            "start": start.tolist(),
            "raw_success": bool(fit.success),
            "status": int(fit.status),
            "message": str(fit.message),
            "nfev": int(fit.nfev),
            "nit": int(fit.nit),
            "endpoint_parameters": fit.x.astype(float).tolist(),
            "reported_chi2": float(fit.fun),
            "independent_endpoint_recheck_chi2": float(endpoint[0]),
            "endpoint_recheck_difference": float(recheck_difference),
            "endpoint_valid_physical": bool(endpoint[1] is None and endpoint[0] < INVALID),
            "objective_rejection_counts": dict(failures),
        })
    return records


def boundary_hits(parameters: list[float], bounds: tuple[tuple[float, float], tuple[float, float]]) -> list[dict]:
    names = ["Omega_m0", "g"]
    hits = []
    for name, value, (lo, hi) in zip(names, parameters, bounds, strict=True):
        tol = 2e-5 * (hi - lo)
        if abs(value - lo) <= tol:
            hits.append({"parameter": name, "side": "lower", "bound": lo})
        if abs(value - hi) <= tol:
            hits.append({"parameter": name, "side": "upper", "bound": hi})
    return hits


def distance_quadrature_check(omega_m: float, g: float) -> dict:
    bg = integrate_background(omega_m, g, np.asarray([0.0, 1.14418]))
    solution = bg["solution"]
    points = [0.1, 0.5, 1.14418]
    checks = []
    for z in points:
        u = math.log1p(z)
        distance_ode = float(solution(u)[2])

        def inv_e(zz: float) -> float:
            m, x, _d = solution(math.log1p(zz))
            e2 = float(m + x)
            if e2 <= 0.0 or not np.isfinite(e2):
                raise DomainError("invalid E^2 in independent scalar quadrature")
            return 1.0 / math.sqrt(e2)

        distance_quad, quad_err = quad(inv_e, 0.0, z, epsabs=2e-12, epsrel=2e-12, limit=200)
        checks.append({
            "z": z,
            "ode_distance": distance_ode,
            "adaptive_scalar_distance": float(distance_quad),
            "absolute_difference": float(abs(distance_ode - distance_quad)),
            "quadrature_reported_error": float(quad_err),
        })
    return {"points": checks, "max_abs_difference": max(x["absolute_difference"] for x in checks)}


def validate_point(omega_m: float, g: float) -> dict:
    # Finer independent sampled trajectory around the selected/reference points.
    dense = integrate_background(omega_m, g, np.asarray([0.0, 1.14418]))
    u_fine = np.linspace(0.0, U_MAX, 4001)
    traj = np.asarray(dense["solution"](u_fine), dtype=np.float64)
    m, x, d = traj
    e2 = m + x
    fb_max = 1.0 if g >= 0.0 else float(m[-1] / (omega_m * math.exp(3.0 * U_MAX)))
    fb = 0.5 * fb_max
    cdm = m - fb * omega_m * np.exp(3.0 * u_fine)
    if np.any(m <= 0.0) or np.any(x <= 0.0) or np.any(e2 <= 0.0) or np.any(cdm <= 0.0):
        raise DomainError("refined accepted-point physicality validation failed")
    distance_check = distance_quadrature_check(omega_m, g)
    profile = prediction_and_profile(omega_m, g, diagnostics=True)
    return {
        "parameters": {"Omega_m0": float(omega_m), "g": float(g)},
        "profile": {"chi2": profile["chi2"], "magnitude_intercept": profile["offset"],
                    "normal_equation_abs_residual": profile["normal_equation_abs_residual"],
                    "independent_algebra": profile["independent_profile_algebra"]},
        "physicality": {
            "guard_z_interval": [0.0, 1.14418],
            "refined_u_grid_points": int(len(u_fine)),
            "fb_max_exact": float(fb_max),
            "fb_witness": float(fb),
            "minimum_matter": float(np.min(m)),
            "minimum_vacuum": float(np.min(x)),
            "minimum_e2": float(np.min(e2)),
            "minimum_cdm_witness": float(np.min(cdm)),
        },
        "distance_quadrature_check": distance_check,
    }


def plot_surface(om_bounds: tuple[float, float], g_bounds: tuple[float, float], surface: np.ndarray,
                 best: list[float], bao: list[float]) -> None:
    om_grid = np.linspace(*om_bounds, GRID_N_OM)
    g_grid = np.linspace(*g_bounds, GRID_N_G)
    delta = surface - float(np.nanmin(surface))
    fig, ax = plt.subplots(figsize=(9.5, 6.8), constrained_layout=True)
    masked = np.ma.masked_where(~np.isfinite(delta) | (delta >= INVALID / 2), delta)
    mesh = ax.pcolormesh(om_grid, g_grid, masked, shading="auto", cmap="viridis_r", rasterized=True)
    finite = delta[np.isfinite(delta) & (delta < INVALID / 2)]
    if finite.size and finite.min() <= 25:
        levels = [x for x in (1, 2.3, 4, 6.18, 9, 11.8, 20) if x <= max(float(np.max(finite)), 1.01)]
        if levels:
            ax.contour(om_grid, g_grid, masked, levels=levels, colors="white", linewidths=0.7, alpha=0.8)
    ax.scatter([best[0]], [best[1]], marker="*", s=140, color="#ffcc00", edgecolor="black", label="SN-only IVS profile minimum")
    ax.scatter([bao[0]], [bao[1]], marker="x", s=80, color="#ff4d6d", linewidth=2, label="BAO-selected IVS reference")
    ax.axhline(0.0, color="white", lw=0.8, alpha=0.8)
    ax.set_xlabel(r"$\Omega_{m0}$ (total baryons + CDM)")
    ax.set_ylabel(r"$g=Gamma/H_0$")
    ax.set_title("Dovekie-only IVS profile: Δχ² from the SN-only minimum\nExploratory search surface; not a posterior or BAO×SN combination")
    fig.colorbar(mesh, ax=ax, label=r"Profile $\Delta\chi^2_{\rm SN}$")
    ax.legend(frameon=True, loc="best")
    fig.savefig(PLOT_PATH, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    global _DATA
    t0 = time.perf_counter()
    _DATA = sn.load_inputs()
    data = _DATA
    if len(data["records"]) != 1820 or data["precision"].shape != (1820, 1820):
        raise ValueError("unexpected Dovekie row/covariance dimension")
    np.linalg.cholesky(data["precision"])

    # Baseline LCDM comes from the already independently audited Dovekie code;
    # the IVS model is separately integrated here and must recover its g=0 limit.
    lcdm = sn.fit_model(data, "LCDM", starts=8)
    test_om = 0.3
    ivs_g0 = prediction_and_profile(test_om, 0.0, diagnostics=True)
    lcdm_mu = sn.background_mu(data["z_hd"], data["z_hel"], np.asarray([test_om]), "LCDM")
    ivs_mu = ivs_g0["mu_base"]
    gamma0_distance_difference = float(np.max(np.abs(ivs_mu - lcdm_mu)))
    if gamma0_distance_difference > 2e-8:
        raise ArithmeticError(f"g=0 did not recover LCDM distances: {gamma0_distance_difference}")

    bounds = (INITIAL_OM_BOUNDS, INITIAL_G_BOUNDS)
    all_attempts: list[dict] = []
    bound_history: list[dict] = []
    surface = None
    scale_records = []
    best_fit = None
    for expansion in range(3):
        surface, scale = evaluate_grid(bounds)
        scale_records.append({"expansion_index": expansion, **scale, "bounds": [list(bounds[0]), list(bounds[1])]})
        finite_surface = np.isfinite(surface) & (surface < INVALID / 2)
        if not np.any(finite_surface):
            raise RuntimeError("coarse profile grid contains no valid physical points")
        grid_flat_index = int(np.argmin(np.where(finite_surface, surface, np.inf)))
        grid_g_index, grid_om_index = np.unravel_index(grid_flat_index, surface.shape)
        grid_minimum = (
            float(np.linspace(*bounds[0], GRID_N_OM)[grid_om_index]),
            float(np.linspace(*bounds[1], GRID_N_G)[grid_g_index]),
        )
        attempts = optimize_profile(bounds, grid_minimum)
        all_attempts.extend([{**row, "bounds": [list(bounds[0]), list(bounds[1])], "expansion_index": expansion}
                             for row in attempts])
        valid = [a for a in attempts if a["endpoint_valid_physical"]]
        if not valid:
            raise RuntimeError("no valid optimizer endpoints in declared box")
        best_fit = min(valid, key=lambda a: a["independent_endpoint_recheck_chi2"])
        point = best_fit["endpoint_parameters"]
        hits = boundary_hits(point, bounds)
        bound_history.append({"bounds": [list(bounds[0]), list(bounds[1])], "best_parameters": point, "boundary_hits": hits})
        if not hits or expansion == 2:
            break
        next_om, next_g = list(bounds[0]), list(bounds[1])
        for hit in hits:
            if hit["parameter"] == "Omega_m0":
                if hit["side"] == "lower":
                    next_om[0] = max(0.005, next_om[0] - 0.05)
                else:
                    next_om[1] = min(MAX_OM, next_om[1] + 0.10)
            elif hit["parameter"] == "g":
                if hit["side"] == "lower":
                    next_g[0] = max(-MAX_ABS_G, next_g[0] - 2.0)
                else:
                    next_g[1] = min(MAX_ABS_G, next_g[1] + 2.0)
        bounds = (tuple(next_om), tuple(next_g))

    assert surface is not None and best_fit is not None
    best_point = best_fit["endpoint_parameters"]
    best_om, best_g = map(float, best_point)
    ivs_min = validate_point(best_om, best_g)
    tight_options = {"rtol": 2.0e-12, "atol": 2.0e-14, "max_step_u": MAX_STEP_U / 2.0}
    tight_ivs = prediction_and_profile(best_om, best_g, diagnostics=True,
                                        ode_options=tight_options)
    tight_chi2_difference = abs(tight_ivs["chi2"] - ivs_min["profile"]["chi2"])
    # validate_point intentionally stores only scalar diagnostics; obtain the
    # reference-accuracy vector directly for a meaningful convergence check.
    default_ivs = prediction_and_profile(best_om, best_g, diagnostics=True)
    tight_mu_max_difference = float(np.max(np.abs(tight_ivs["mu_base"] - default_ivs["mu_base"])))
    if tight_chi2_difference > 1e-5 or tight_mu_max_difference > 1e-7:
        raise ArithmeticError("tightened ODE tolerance/step changed the selected profile materially")
    bao_point = (0.38697100769747017, -0.4666647299169934)
    bao_profile = validate_point(*bao_point)
    lcdm_point = (float(lcdm["parameters"]["Omega_m"]), 0.0)
    lcdm_profile = validate_point(*lcdm_point)
    # fit_model returns two row-level plotting arrays; they are not needed in
    # this compact profile artifact and NumPy arrays are not JSON-native.
    lcdm_record = {key: value for key, value in lcdm.items()
                   if not isinstance(value, np.ndarray)}
    bao_delta = float(bao_profile["profile"]["chi2"] - ivs_min["profile"]["chi2"])
    lcdm_delta = float(lcdm_profile["profile"]["chi2"] - ivs_min["profile"]["chi2"])
    final_hits = boundary_hits(best_point, bounds)
    flat_starts = []
    for start in np.linspace(bounds[0][0], bounds[0][1], 8):
        fit = minimize(lambda x: evaluate_grid_point((float(x[0]), 0.0))[0], np.asarray([start]),
                       method="L-BFGS-B", bounds=[bounds[0]],
                       options={"maxiter": 1000, "ftol": 1e-13, "gtol": 1e-8})
        flat_starts.append({"start_omega_m": float(start), "success": bool(fit.success),
                            "status": int(fit.status), "chi2": float(fit.fun), "omega_m": float(fit.x[0])})

    np.savez_compressed(
        GRID_PATH,
        omega_m=np.linspace(*bounds[0], GRID_N_OM),
        g=np.linspace(*bounds[1], GRID_N_G),
        chi2_surface=surface,
    )
    plot_surface(bounds[0], bounds[1], surface, best_point, list(bao_point))
    runtime = time.perf_counter() - t0
    result = {
        "status": "complete_exploratory_dovekie_only_ivs_profile",
        "objective": "Full STAT+SYS DES-Dovekie profile chi-square for flat late-time IVS background with an analytically profiled common magnitude intercept.",
        "scope_and_limits": [
            "SN-only reduced background; no BAO likelihood is multiplied in and no exact independence is asserted.",
            "No H0 or dimensional Gamma inference; H0=70 km/s/Mpc is only an absorbed distance-modulus gauge.",
            "Finite-search profile minima are not posterior/evidence, p-value, sigma, Wilks interval, or discovery evidence.",
            "Positive baryon/CDM split is existential and uncalibrated; no fb estimate. Physical checks cover only the observed Dovekie redshift range.",
            "Radiation, perturbations, CMB/BBN, early-time viability and the paper's modified-CAMB pipeline are not reproduced.",
        ],
        "likelihood_contract": {
            "row_count": 1820,
            "sample_order": "Pinned upstream release order, unchanged; no CID sorting or metadata join.",
            "covariance": "Full released STAT+SYS covariance represented by packed float32 upper-triangle precision P=C^-1; no separate MUERR addition.",
            "magnitude_intercept": "One unbounded common additive magnitude, profiled analytically by full-covariance GLS.",
            "distance": "DL=(1+zHEL)*(c/H0_gauge)*integral_0^zHD dz/E(z); zHD enters expansion and zHEL the prefactor.",
            "normalization": "Only profile chi-square differences are reported; fixed covariance logdet and Gaussian constant omitted from minima differences.",
            "H0_gauge_km_s_mpc": H0_GAUGE,
        },
        "background_model": {
            "source_convention": "Q=Gamma*rho_x>0 transfers CDM to vacuum; baryons separately conserved.",
            "u": "ln(1+z)=-ln(a)",
            "variables": "m=(rho_b+rho_c)/rho_crit0, x=rho_x/rho_crit0, E^2=m+x, g=Gamma/H0.",
            "ode": ["dx/du=-g*x/E", "dm/du=3*m+g*x/E", "dD/du=exp(u)/E"],
            "initial_conditions": ["m(0)=Omega_m0", "x(0)=1-Omega_m0", "D(0)=0"],
            "scope": "Flat late-time matter plus vacuum; radiation/neutrinos/perturbations omitted.",
        },
        "search": {
            "initial_bounds": {"Omega_m0": list(INITIAL_OM_BOUNDS), "g": list(INITIAL_G_BOUNDS)},
            "bounds_role": "Finite numerical search box; Omega_m bounds inherited from BAO code and are not a paper prior. g box matches paper Table I but no SN posterior measure is applied.",
            "final_bounds": {"Omega_m0": list(bounds[0]), "g": list(bounds[1])},
            "bound_history": bound_history,
            "grid_shape": [GRID_N_G, GRID_N_OM],
            "grid_scale_records": scale_records,
        "optimizer_method": "L-BFGS-B; deterministic reference/corner starts plus seeded uniform starts; every endpoint independently reevaluated.",
            "coarse_grid_minimum_included_as_optimizer_start": True,
            "optimizer_start_count_per_box": N_STARTS,
            "optimizer_attempts": all_attempts,
            "final_best_point_boundary_hits": final_hits,
            "profile_is_domain_limited": bool(final_hits),
            "flat_g0_multistarts": flat_starts,
        },
        "results": {
            "ivs_sn_only_minimum": ivs_min,
            "bao_selected_ivs_reference_on_sn_profile": bao_profile,
            "delta_chi2_bao_point_minus_sn_minimum": bao_delta,
            "dovekie_lcdm_reference": {"parameters": {"Omega_m0": lcdm_point[0], "g": 0.0},
                                         "chi2": lcdm_profile["profile"]["chi2"],
                                         "source_flat_profile_fit": lcdm_record},
            "delta_chi2_lcdm_minimum_minus_ivs_minimum": lcdm_delta,
            "interpretation_delta_sign": "Positive delta means the listed reference point has a higher SN-only profiled chi-square than the IVS SN-only minimum.",
        },
        "numerical_checks": {
            "precision_cholesky_pass": True,
            "precision_symmetric": bool(np.array_equal(data["precision"], data["precision"].T)),
            "precision_min_cholesky_diagonal": float(np.min(np.diag(data["precision_cholesky"]))),
            "g0_vs_flat_lcdm_max_distance_modulus_difference_mag": gamma0_distance_difference,
            "best_point_and_references_refined_physicality_pass": True,
            "best_point_tight_ode_recheck": {
                "default_rtol": RTOL,
                "default_atol": ATOL,
                "default_max_step_u": MAX_STEP_U,
                "tight_rtol": tight_options["rtol"],
                "tight_atol": tight_options["atol"],
                "tight_max_step_u": tight_options["max_step_u"],
                "absolute_profile_chi2_difference": tight_chi2_difference,
                "maximum_distance_modulus_difference_mag": tight_mu_max_difference,
                "pass_threshold_chi2": 1e-5,
                "pass_threshold_mu_mag": 1e-7,
            },
            "max_adaptive_distance_quadrature_difference": max(
                ivs_min["distance_quadrature_check"]["max_abs_difference"],
                bao_profile["distance_quadrature_check"]["max_abs_difference"],
                lcdm_profile["distance_quadrature_check"]["max_abs_difference"],
            ),
            "physicality_grid_step_u_max": float(U_MAX / (max(257, int(math.ceil(U_MAX / 0.002)) + 1) - 1)),
            "refined_physicality_grid_points": 4001,
        },
        "artifacts": {
            "result": str(RESULT_PATH.relative_to(ROOT)),
            "surface_npz": str(GRID_PATH.relative_to(ROOT)),
            "plot": str(PLOT_PATH.relative_to(ROOT)),
        },
        "inputs": {
            "hubble_diagram_sha256": data["hashes"]["des_dovekie_hd.csv"],
            "packed_precision_sha256": data["hashes"]["des_dovekie_stat_sys.npz"],
            "cid_order_sha256": hashlib.sha256(("\n".join(row["CID"] for row in data["records"]) + "\n").encode()).hexdigest(),
        },
        "runtime_seconds": float(runtime),
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cpu_count_visible": int(os.cpu_count() or 1),
            "blas_threads_requested": 1,
            "gpu_used": False,
        },
        "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --seconds 1200 -- .venv/bin/python experiments/dovekie_ivs_profile/profile.py",
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    record = {
        "id": "des_dovekie_ivs_profile_screen_20260925",
        "status": "exploratory_profile_screen_root_implementation_independent_review_pending",
        "objective": result["objective"],
        "contract": "work/theory7/dovekie_ivs_profile_contract.md/.json",
        "independent_physics_review": "work/theory8/dovekie_profile_physics_critic.md/.json",
        "command": result["command"],
        "runtime_seconds": float(runtime),
        "result": str(RESULT_PATH.relative_to(ROOT)),
        "source_hashes": {
            "profile_script": sha256(Path(__file__)),
            "dovekie_profile_helper": sha256(ROOT / "work/inference/fit_dovekie.py"),
            "bao_ivs_reference_adapter": sha256(ROOT / "experiments/interacting_vacuum_screen/interacting_vacuum_profile.py"),
            "physicality_contract": sha256(ROOT / "work/theory7/dovekie_ivs_profile_contract.json"),
            "independent_physics_review": sha256(ROOT / "work/theory8/dovekie_profile_physics_critic.json"),
        },
        "output_hashes": {
            "result": sha256(RESULT_PATH),
            "surface_npz": sha256(GRID_PATH),
            "plot": sha256(PLOT_PATH),
        },
        "limits": result["scope_and_limits"],
    }
    RECORD_PATH.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "runtime_seconds": runtime,
                      "best": best_point, "chi2_min": ivs_min["profile"]["chi2"],
                      "bao_point_delta": bao_delta, "lcdm_delta": lcdm_delta,
                      "parallel_grid_used": scale_records[-1]["parallel_grid_used"],
                      "profile_domain_limited": bool(final_hits)}, allow_nan=False))


if __name__ == "__main__":
    main()

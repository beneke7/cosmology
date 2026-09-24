#!/usr/bin/env python3
"""Late-time, flat interacting-vacuum profile screen against 13-row DESI DR2 BAO."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import solve_triangular
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import background_bao as bao  # noqa: E402

DT = np.float64
ALPHA_BOUNDS = (1e-6, 1e4)
OM_BOUNDS = (0.05, 0.60)
GAMMA_BOUNDS = (-3.0, 3.0)
ZMAX = 2.33
RTOL_MAIN = 2e-10
ATOL_MAIN = 2e-12
INVALID = 1e80


class DomainError(ValueError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def data_load():
    mean = ROOT / "context/data/desi_dr2_mean.txt"
    cov = ROOT / "context/data/desi_dr2_cov.txt"
    return bao.BAOData.from_files(mean, cov)


def integrate_profile(omega_m: float, g: float, z_eval: np.ndarray, *, rtol=RTOL_MAIN, atol=ATOL_MAIN):
    """Integrate total non-DE matter, vacuum, and dimensionless distance in u=ln(1+z).

    Returned density histories are normalized to today's critical density. Radiation
    and neutrinos are deliberately omitted per the explicitly limited low-z contract.
    """
    om, g = float(omega_m), float(g)
    z_eval = np.asarray(z_eval, dtype=DT)
    if not (np.isfinite(om) and np.isfinite(g)) or not OM_BOUNDS[0] <= om <= OM_BOUNDS[1]:
        raise DomainError("Omega_m outside finite profile domain")
    if not GAMMA_BOUNDS[0] <= g <= GAMMA_BOUNDS[1]:
        raise DomainError("Gamma/H0 outside finite profile domain")
    if np.any(~np.isfinite(z_eval)) or np.any(z_eval < 0.0) or np.any(z_eval > ZMAX):
        raise DomainError("redshift outside [0,2.33]")
    x0 = 1.0 - om
    if om <= 0.0 or x0 <= 0.0:
        raise DomainError("present matter and vacuum densities must be positive")
    umax = float(np.log1p(max(float(np.max(z_eval)), 0.0)))

    def rhs(u, state):
        matter, vacuum, distance = state
        e2 = matter + vacuum
        if not np.isfinite(e2) or e2 <= 0.0:
            raise DomainError("E^2 <= 0 during ODE integration")
        e = math.sqrt(e2)
        source = g * vacuum / e
        return (3.0 * matter + source, -source, math.exp(u) / e)

    sol = solve_ivp(
        rhs, (0.0, umax), (om, x0, 0.0), method="DOP853", dense_output=True,
        rtol=rtol, atol=atol, max_step=0.025,
    )
    if not sol.success or sol.sol is None:
        raise DomainError(f"ODE solver failed: {sol.message}")
    # Check the complete observed redshift interval, not only BAO row endpoints.
    ucheck = np.linspace(0.0, umax, max(257, int(math.ceil(umax / 0.002)) + 1), dtype=DT)
    hist = np.asarray(sol.sol(ucheck), dtype=DT)
    matter, vacuum = hist[0], hist[1]
    e2 = matter + vacuum
    if np.any(~np.isfinite(hist)) or np.any(e2 <= 0.0) or np.any(matter <= 0.0) or np.any(vacuum <= 0.0):
        raise DomainError("nonpositive/nonfinite matter, vacuum, or expansion on observed interval")
    # Since B(u)=f_b*Omega_m*exp(3u), this is the largest allowed baryon
    # fraction for which C(u)=M(u)-B(u) remains nonnegative everywhere.
    ratio = matter / (om * np.exp(3.0 * ucheck))
    fb_max = float(min(1.0, np.min(ratio)))
    if not np.isfinite(fb_max) or fb_max <= 0.0:
        raise DomainError("no positive baryon/CDM split keeps both components nonnegative")
    fb_witness = 0.5 * fb_max
    z_eval = np.asarray(z_eval, dtype=DT)
    out = np.asarray(sol.sol(np.log1p(z_eval)), dtype=DT)
    m_eval, x_eval, d_eval = out
    e2_eval = m_eval + x_eval
    if np.any(e2_eval <= 0.0) or np.any(~np.isfinite(out)):
        raise DomainError("invalid requested-redshift state")
    # Verify an interior, positive baryon/CDM decomposition as a physical witness.
    b_eval = fb_witness * om * (1.0 + z_eval) ** 3
    c_eval = m_eval - b_eval
    if np.any(b_eval <= 0.0) or np.any(c_eval <= 0.0):
        raise DomainError("interior baryon/CDM witness is not positive")
    return {
        "matter": m_eval, "vacuum": x_eval, "distance": d_eval,
        "e2": e2_eval, "e": np.sqrt(e2_eval), "fb_max": fb_max,
        "fb_witness": fb_witness, "minimum_e2": float(np.min(e2)),
        "minimum_matter": float(np.min(matter)), "minimum_vacuum": float(np.min(vacuum)),
        "solution": sol.sol,
    }


def unit_prediction(data, omega_m: float, g: float, **kwargs):
    bg = integrate_profile(omega_m, g, data.z, **kwargs)
    dm = bg["distance"]
    dh = 1.0 / bg["e"]
    dv = np.cbrt(data.z * dm * dm * dh)
    q = np.empty_like(data.value, dtype=DT)
    for label, values in (("DM_over_rs", dm), ("DH_over_rs", dh), ("DV_over_rs", dv)):
        q[data.observable == label] = values[data.observable == label]
    return q, bg


def profile_alpha(data, chol, y_white, omega_m: float, g: float, **kwargs):
    q, bg = unit_prediction(data, omega_m, g, **kwargs)
    qw = solve_triangular(chol, q, lower=True, check_finite=False)
    alpha = float(np.clip(np.dot(qw, y_white) / np.dot(qw, qw), *ALPHA_BOUNDS))
    residual = y_white - alpha * qw
    chi2 = float(np.dot(residual, residual))
    return chi2, alpha, alpha * q, q, bg


def fixed_split_check(omega_m: float, g: float, f_b: float) -> dict:
    """Check separately conserved baryons and the interacting CDM at a fixed split."""
    bg = integrate_profile(float(omega_m), float(g), np.array([0.0, ZMAX]))
    u = np.linspace(0.0, math.log1p(ZMAX), 1001)
    m, x, _ = np.asarray(bg["solution"](u), dtype=DT)
    b = float(f_b)*float(omega_m)*np.exp(3*u)
    c = m-b
    return {"f_b": float(f_b), "physical_on_0_to_zmax": bool(np.all(b > 0.0) and np.all(c > 0.0)),
            "minimum_baryon_density_fraction": float(np.min(b)),
            "minimum_CDM_density_fraction": float(np.min(c)),
            "f_b_max_for_some_positive_split": bg["fb_max"],
            "minimum_total_matter_fraction": float(np.min(m)),
            "minimum_vacuum_fraction": float(np.min(x))}


def optimizer_endpoint_audit(attempts: list[dict], bounds: list[tuple[float, float]], tolerance=0.005) -> dict:
    """Separate SciPy flags from physical objective endpoints and heuristically cluster them."""
    valid = [a for a in attempts if a.get("endpoint_valid_physical", False)
             and np.isfinite(a["chi2"]) and a["chi2"] < INVALID]
    invalid = [a for a in attempts if a not in valid]
    raw_success = [a for a in attempts if a["success"]]
    successful_valid = [a for a in valid if a["success"]]
    failed_valid = [a for a in valid if not a["success"]]
    raw_success_invalid = [a for a in invalid if a["success"]]

    # Connected components under a declared normalized-coordinate tolerance. This
    # describes endpoint grouping only; it is not proof of separate attraction basins.
    clusters: list[set[int]] = []
    lo = np.asarray([b[0] for b in bounds], dtype=DT)
    span = np.asarray([b[1]-b[0] for b in bounds], dtype=DT)
    points = np.asarray([a["x"] for a in valid], dtype=DT) if valid else np.empty((0, len(bounds)))
    scaled = (points-lo)/span if valid else points
    remaining = set(range(len(valid)))
    while remaining:
        seed = min(remaining)
        group = {seed}
        stack = [seed]
        remaining.remove(seed)
        while stack:
            current = stack.pop()
            neighbors = [j for j in list(remaining) if np.linalg.norm(scaled[current]-scaled[j]) <= tolerance]
            for j in neighbors:
                remaining.remove(j)
                group.add(j)
                stack.append(j)
        clusters.append(group)
    cluster_records = []
    for group in clusters:
        members = [valid[i] for i in group]
        representative = min(members, key=lambda a: a["chi2"])
        cluster_records.append({"member_count": len(members),
                                "raw_success_count": sum(a["success"] for a in members),
                                "best_chi2": float(representative["chi2"]),
                                "best_x": representative["x"]})
    cluster_records.sort(key=lambda a: a["best_chi2"])
    return {
        "start_count": len(attempts),
        "raw_scipy_success_count": len(raw_success),
        "raw_scipy_failure_count": len(attempts)-len(raw_success),
        "valid_finite_physical_endpoint_count": len(valid),
        "endpoint_rechecked_at_reported_coordinates": True,
        "invalid_penalty_endpoint_count": len(invalid),
        "raw_success_at_invalid_penalty_count": len(raw_success_invalid),
        "valid_endpoint_success_count": len(successful_valid),
        "valid_endpoint_failure_count": len(failed_valid),
        "valid_endpoint_cluster_count": len(cluster_records),
        "endpoint_cluster_method": "connected components using Euclidean distance <= 0.005 after each parameter is scaled to [0,1] in its search box; descriptive endpoint grouping only, not a basin/globality proof",
        "valid_endpoint_clusters_sorted_by_best_chi2": cluster_records,
    }


def run_fit(data, chol, y_white, *, gamma_fixed=None, gamma_bounds=GAMMA_BOUNDS,
            starts=32, rtol=RTOL_MAIN, atol=ATOL_MAIN):
    bounds = [OM_BOUNDS] if gamma_fixed is not None else [OM_BOUNDS, gamma_bounds]
    failures = Counter()

    def evaluate_endpoint(theta):
        om = float(theta[0])
        g = 0.0 if gamma_fixed is None and len(theta) == 1 else (float(gamma_fixed) if gamma_fixed is not None else float(theta[1]))
        try:
            score = profile_alpha(data, chol, y_white, om, g, rtol=rtol, atol=atol)[0]
            valid = bool(np.isfinite(score) and score < INVALID)
            return score, valid, None if valid else "nonfinite or penalty objective"
        except (DomainError, ValueError, FloatingPointError, OverflowError) as exc:
            return INVALID, False, str(exc)

    def objective(theta):
        score, valid, error = evaluate_endpoint(theta)
        if not valid:
            failures[str(error)] += 1
            return INVALID
        return score

    if gamma_fixed is None:
        om_starts = np.linspace(*OM_BOUNDS, 7)
        g_starts = np.linspace(*gamma_bounds, 5)
        starts_list = [np.array([om, g], dtype=DT) for om in om_starts for g in g_starts]
        # Deterministic extra diagonal/interior starts to supplement the 35-point mesh.
        diagonal_g = min(0.4, max(0.05, 0.25*(gamma_bounds[1]-gamma_bounds[0])))
        starts_list += [np.array([0.12, -diagonal_g]), np.array([0.22, diagonal_g]),
                        np.array([0.42, -diagonal_g]), np.array([0.52, diagonal_g])]
        if starts < len(starts_list):
            # Spread a smaller requested set over the full deterministic mesh.
            indices = np.linspace(0, len(starts_list)-1, starts).round().astype(int)
            starts_list = [starts_list[i] for i in indices]
    else:
        vals = np.linspace(*OM_BOUNDS, max(starts - 1, 2))
        starts_list = [np.array([float(v)], dtype=DT) for v in vals]
        starts_list.append(np.array([0.3], dtype=DT))
        starts_list = starts_list[:starts]
    attempts = []
    for i, x0 in enumerate(starts_list):
        before = failures.copy()
        opt = minimize(objective, x0, method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": 1200, "ftol": 1e-14, "gtol": 2e-8, "maxls": 40})
        endpoint_chi2, endpoint_valid, endpoint_error = evaluate_endpoint(opt.x)
        attempts.append({
            "start_index": i, "x0": [float(x) for x in x0], "x": [float(x) for x in opt.x],
            "chi2": float(endpoint_chi2), "scipy_fun": float(opt.fun),
            "endpoint_valid_physical": endpoint_valid, "endpoint_validation_error": endpoint_error,
            "success": bool(opt.success), "status": int(opt.status),
            "message": str(opt.message), "nit": int(opt.nit), "nfev": int(opt.nfev),
            "domain_failures_during_attempt": dict(failures - before),
        })
    finite = [a for a in attempts if a["endpoint_valid_physical"] and np.isfinite(a["chi2"]) and a["chi2"] < INVALID]
    if not finite:
        raise RuntimeError("no finite optimizer candidate")
    successful = [a for a in finite if a["success"]]
    best = min(successful or finite, key=lambda a: a["chi2"])
    failed_best = min((a for a in finite if not a["success"]), key=lambda a: a["chi2"], default=None)
    theta = np.asarray(best["x"], dtype=DT)
    best_g = 0.0 if gamma_fixed is not None else float(theta[1])
    score, alpha, pred, q, bg = profile_alpha(data, chol, y_white, float(theta[0]), best_g, rtol=rtol, atol=atol)
    return {
        "Omega_m": float(theta[0]), "Gamma_over_H0": best_g, "alpha": alpha,
        "alpha_bounds": ALPHA_BOUNDS, "shape_search_bounds": {"Omega_m": OM_BOUNDS, "Gamma_over_H0": (None if gamma_fixed is not None else gamma_bounds)},
        "alpha_bound_active": bool(abs(alpha-ALPHA_BOUNDS[0]) < 1e-12 or abs(alpha-ALPHA_BOUNDS[1]) < 1e-9),
        "chi2": score, "prediction": pred.tolist(), "unit_prediction": q.tolist(),
        "selected_start_index": best["start_index"], "selected_success": best["success"],
        "selected_status": best["status"], "starts_requested": len(starts_list),
        "optimizer_endpoint_audit": optimizer_endpoint_audit(attempts, bounds), "optimizer_attempts": attempts,
        "best_unsuccessful_candidate": (None if failed_best is None else {
            "start_index": failed_best["start_index"], "x": failed_best["x"], "chi2": failed_best["chi2"],
            "status": failed_best["status"], "message": failed_best["message"],
            "score_delta_vs_selected_success": float(failed_best["chi2"]-score)}),
        "domain_failure_evaluations_by_reason": dict(failures),
        "physical_decomposition_witness": {"f_b_open_interval": [0.0, bg["fb_max"]], "example_f_b": bg["fb_witness"]},
        "domain_minima": {k: bg[k] for k in ("minimum_e2", "minimum_matter", "minimum_vacuum")},
        "parameter_bound_hits": (
            (["Omega_m:lower"] if abs(theta[0] - OM_BOUNDS[0]) < 1e-8 else []) +
            (["Omega_m:upper"] if abs(theta[0] - OM_BOUNDS[1]) < 1e-8 else []) +
            ([] if gamma_fixed is not None else
             (["Gamma_over_H0:lower"] if abs(theta[1] - gamma_bounds[0]) < 1e-7 else []) +
             (["Gamma_over_H0:upper"] if abs(theta[1] - gamma_bounds[1]) < 1e-7 else []))
        ),
        "all_optimizer_statuses_retained": True,
    }


def adaptive_unit(data, omega_m, g, *, rtol=3e-12, atol=1e-14):
    """Independent scalar-distance check: dense ODE for E and adaptive QUAD for DM."""
    om, gg = float(omega_m), float(g)
    x0 = 1.0 - om
    umax = float(np.log1p(float(np.max(data.z))))

    def rhs(u, state):
        m, x = state
        e2 = m + x
        if e2 <= 0 or not np.isfinite(e2):
            raise DomainError("adaptive ODE E^2 failure")
        src = gg * x / math.sqrt(e2)
        return (3*m + src, -src)

    sol = solve_ivp(rhs, (0, umax), (om, x0), method="DOP853", dense_output=True,
                    rtol=rtol, atol=atol, max_step=0.0125)
    if not sol.success or sol.sol is None:
        raise DomainError("adaptive reference ODE failed")

    def inv_e(z):
        m, x = sol.sol(math.log1p(float(z)))
        e2 = float(m + x)
        if e2 <= 0 or not np.isfinite(e2):
            raise DomainError("adaptive reference E^2 failure")
        return 1.0 / math.sqrt(e2)

    dm = np.array([quad(inv_e, 0.0, float(z), epsabs=2e-12, epsrel=2e-12, limit=200)[0] for z in data.z])
    e = np.array([1.0 / inv_e(float(z)) for z in data.z])
    dh, dv = 1.0/e, np.cbrt(data.z*dm*dm/e)
    q = np.empty_like(data.value)
    for label, values in (("DM_over_rs", dm), ("DH_over_rs", dh), ("DV_over_rs", dv)):
        q[data.observable == label] = values[data.observable == label]
    return q


def profile_from_q(data, chol, y_white, q):
    qw = solve_triangular(chol, q, lower=True, check_finite=False)
    alpha = float(np.clip(np.dot(qw, y_white)/np.dot(qw, qw), *ALPHA_BOUNDS))
    pred = alpha*q
    r = y_white - alpha*qw
    return float(np.dot(r, r)), alpha, pred


def local_profile_curvature(data, chol, y_white, omega_m, g, h_om, h_g):
    """Finite-difference local curvature of profiled chi2; no interval interpretation."""
    f = lambda om, gg: profile_alpha(data, chol, y_white, float(om), float(gg))[0]
    f0 = f(omega_m, g)
    h11 = (f(omega_m+h_om,g)-2*f0+f(omega_m-h_om,g))/(h_om*h_om)
    h22 = (f(omega_m,g+h_g)-2*f0+f(omega_m,g-h_g))/(h_g*h_g)
    h12 = (f(omega_m+h_om,g+h_g)-f(omega_m+h_om,g-h_g)-
           f(omega_m-h_om,g+h_g)+f(omega_m-h_om,g-h_g))/(4*h_om*h_g)
    hess = np.array([[h11,h12],[h12,h22]], dtype=DT)
    eig = np.linalg.eigvalsh(hess)
    out = {"steps_Omega_m_Gamma_over_H0": [float(h_om),float(h_g)],
           "chi2_hessian": hess.tolist(), "eigenvalues": eig.tolist(),
           "positive_definite_local_curvature": bool(np.all(eig>0.0)),
           "ridge_slope_dOmega_m_dGamma_over_H0": float(-h12/h11) if h11 else None,
           "interpretation": "local curvature/correlation of the profiled chi2 surface only; not a posterior, confidence interval, or calibrated uncertainty"}
    if np.all(eig>0.0):
        inv = np.linalg.inv(hess)
        out["inverse_curvature_correlation"] = float(inv[0,1]/np.sqrt(inv[0,0]*inv[1,1]))
    else:
        out["inverse_curvature_correlation"] = None
    return out


def fixed_coupling_profile(data, chol, y_white, best_g, best_chi2):
    offsets = (-0.30,-0.20,-0.10,0.0,0.10,0.20,0.30)
    g_values = sorted(set([float(np.clip(best_g+d,*GAMMA_BOUNDS)) for d in offsets] + [0.0]))
    rows=[]
    for g in g_values:
        fit = run_fit(data, chol, y_white, gamma_fixed=g, starts=8)
        rows.append({"Gamma_over_H0_fixed":g,"Omega_m_profiled":fit["Omega_m"],"alpha_profiled":fit["alpha"],
                     "chi2_profile":fit["chi2"],"delta_chi2_vs_interacting_best":float(fit["chi2"]-best_chi2),
                     "delta_chi2_vs_gamma0":None,"optimizer_endpoint_audit":fit["optimizer_endpoint_audit"],
                     "optimizer_attempts":fit["optimizer_attempts"]})
    gamma0 = next((row["chi2_profile"] for row in rows if row["Gamma_over_H0_fixed"]==0.0), None)
    for row in rows:
        if gamma0 is not None: row["delta_chi2_vs_gamma0"] = float(row["chi2_profile"]-gamma0)
    return {"method":"fixed-gamma bounded Omega_m profile using 8 deterministic L-BFGS-B starts at each point",
            "rows":rows,"not_an_interval":"The profile curve is descriptive; no threshold crossing is interpreted as a posterior interval or calibrated confidence interval."}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--coarse-grid", type=int, default=41)
    p.add_argument("--starts", type=int, default=39)
    args = p.parse_args()
    started = time.time()
    data = data_load()
    chol = np.linalg.cholesky(data.covariance)
    y_white = solve_triangular(chol, data.value, lower=True, check_finite=False)

    flat = run_fit(data, chol, y_white, gamma_fixed=0.0, starts=16)
    fit = run_fit(data, chol, y_white, gamma_fixed=None, starts=args.starts)

    # Dense shape grid, followed by a bounded local refinement; grid is a search check,
    # not a proof of globality.
    grid_n = max(11, args.coarse_grid)
    oms = np.linspace(*OM_BOUNDS, grid_n)
    gs = np.linspace(*GAMMA_BOUNDS, grid_n)
    grid = np.full((grid_n, grid_n), np.nan, dtype=DT)
    grid_failures = Counter()
    for i, om in enumerate(oms):
        for j, g in enumerate(gs):
            try:
                grid[i, j] = profile_alpha(data, chol, y_white, float(om), float(g))[0]
            except (DomainError, ValueError, FloatingPointError, OverflowError) as exc:
                grid_failures[str(exc)] += 1
    grid_idx = np.unravel_index(np.nanargmin(grid), grid.shape)
    grid_x0 = np.array([oms[grid_idx[0]], gs[grid_idx[1]]], dtype=DT)
    grid_refine_failures = Counter()
    def grid_refine_objective(th):
        try:
            return profile_alpha(data, chol, y_white, float(th[0]), float(th[1]))[0]
        except (DomainError, ValueError, FloatingPointError, OverflowError) as exc:
            grid_refine_failures[str(exc)] += 1
            return INVALID
    multi_start_score = fit["chi2"]
    grid_opt = minimize(grid_refine_objective,
                        grid_x0, method="L-BFGS-B", bounds=[OM_BOUNDS, GAMMA_BOUNDS],
                        options={"maxiter": 1200, "ftol": 1e-14, "gtol": 2e-8})
    grid_refined = {
        "x0": grid_x0.tolist(), "x": [float(v) for v in grid_opt.x], "chi2": float(grid_opt.fun),
        "success": bool(grid_opt.success), "status": int(grid_opt.status), "message": str(grid_opt.message),
        "grid_n_per_axis": grid_n, "grid_valid_points": int(np.isfinite(grid).sum()),
        "grid_invalid_points": int((~np.isfinite(grid)).sum()), "grid_domain_failures": dict(grid_failures),
        "refinement_domain_failures": dict(grid_refine_failures),
        "grid_chi2_min": float(grid[grid_idx]), "globality_claim": False,
    }
    if grid_refined["chi2"] < fit["chi2"]:
        om_grid, g_grid = map(float, grid_opt.x)
        s_grid, a_grid, p_grid, q_grid, b_grid = profile_alpha(data, chol, y_white, om_grid, g_grid)
        fit.update({"Omega_m": om_grid, "Gamma_over_H0": g_grid, "alpha": a_grid,
                    "alpha_bounds": ALPHA_BOUNDS,
                    "alpha_bound_active": bool(abs(a_grid-ALPHA_BOUNDS[0]) < 1e-12 or abs(a_grid-ALPHA_BOUNDS[1]) < 1e-9),
                    "chi2": s_grid, "prediction": p_grid.tolist(), "unit_prediction": q_grid.tolist(),
                    "selected_start_index": None, "selected_success": bool(grid_opt.success),
                    "selected_status": int(grid_opt.status), "selected_source": "coarse-grid local refinement",
                    "physical_decomposition_witness": {"f_b_open_interval": [0.0, b_grid["fb_max"]],
                                                        "example_f_b": b_grid["fb_witness"]},
                    "domain_minima": {k: b_grid[k] for k in ("minimum_e2", "minimum_matter", "minimum_vacuum")},
                    "parameter_bound_hits": (["Omega_m:lower"] if abs(om_grid-OM_BOUNDS[0])<1e-8 else []) +
                                            (["Omega_m:upper"] if abs(om_grid-OM_BOUNDS[1])<1e-8 else []) +
                                            (["Gamma_over_H0:lower"] if abs(g_grid-GAMMA_BOUNDS[0])<1e-7 else []) +
                                            (["Gamma_over_H0:upper"] if abs(g_grid-GAMMA_BOUNDS[1])<1e-7 else [])})
    grid_refined["delta_chi2_vs_multistart"] = grid_refined["chi2"] - multi_start_score
    if fit.get("best_unsuccessful_candidate") is not None:
        fit["best_unsuccessful_candidate"]["score_delta_vs_selected_success"] = float(
            fit["best_unsuccessful_candidate"]["chi2"] - fit["chi2"])

    # Independent derivative-free optimizer, with adaptive scalar quadrature for D_M.
    adaptive_fail = Counter()
    def adaptive_obj(th):
        try:
            q = adaptive_unit(data, float(th[0]), float(th[1]))
            return profile_from_q(data, chol, y_white, q)[0]
        except (DomainError, ValueError, FloatingPointError, OverflowError) as exc:
            adaptive_fail[str(exc)] += 1
            return INVALID
    powell = minimize(adaptive_obj, [fit["Omega_m"], fit["Gamma_over_H0"]], method="Powell",
                      bounds=[OM_BOUNDS, GAMMA_BOUNDS],
                      options={"maxiter": 250, "xtol": 2e-6, "ftol": 1e-11})
    q_adapt = adaptive_unit(data, fit["Omega_m"], fit["Gamma_over_H0"])
    adaptive_score, adaptive_alpha, adaptive_pred = profile_from_q(data, chol, y_white, q_adapt)

    # Tight-tolerance convergence check at the selected primary optimum.
    score_tight, alpha_tight, pred_tight, _, _ = profile_alpha(
        data, chol, y_white, fit["Omega_m"], fit["Gamma_over_H0"], rtol=2e-12, atol=2e-14)

    # Gamma=0 must reproduce the background starter at identical (alpha, Omega_m).
    bkg = integrate_profile(flat["Omega_m"], 0.0, data.z)
    bkg_zero = integrate_profile(flat["Omega_m"], 0.0, np.array([0.0]))
    q0, _ = unit_prediction(data, flat["Omega_m"], 0.0)
    q_ref = bao.predict_bao(data.z, data.observable, 1.0, flat["Omega_m"])
    flat_limit = {"max_abs_unit_prediction_difference_vs_background_bao": float(np.max(np.abs(q0-q_ref))),
                  "E0": float(bkg_zero["e"][0]),
                  "flat_fit_score_difference_vs_seed": float(flat["chi2"] - 10.271041002565878)}

    # Zero-noise exact synthetic recovery and an explicit baryon-fraction degeneracy test.
    truth = {"Omega_m": 0.31, "Gamma_over_H0": 0.20, "alpha": 31.0, "f_b_generating_example": 0.16}
    q_truth, bg_truth = unit_prediction(data, truth["Omega_m"], truth["Gamma_over_H0"])
    synthetic_y = truth["alpha"] * q_truth
    synthetic_white = solve_triangular(chol, synthetic_y, lower=True, check_finite=False)
    synth = run_fit(data, chol, synthetic_white, gamma_fixed=None, starts=max(16, min(args.starts, 24)))
    # Distances depend on the sum B+C; several positive present-day splits share the same M.
    split_checks = []
    for fb in (0.01, 0.16, 0.50):
        cvals = bg_truth["matter"] - fb*truth["Omega_m"]*(1.0+data.z)**3
        split_checks.append({"f_b": fb, "minimum_CDM_at_data_rows": float(np.min(cvals)),
                             "same_total_background": True})
    split_invariance = float(np.max(np.abs(q_truth - unit_prediction(data, truth["Omega_m"], truth["Gamma_over_H0"])[0])))
    fixed_split_sensitivity = []
    for fb in (0.10, 0.15, 0.20):
        check = fixed_split_check(fit["Omega_m"], fit["Gamma_over_H0"], fb)
        check["likelihood_profile_score_if_physical"] = (fit["chi2"] if check["physical_on_0_to_zmax"] else None)
        check["not_used_as_prior_or_fit_parameter"] = True
        fixed_split_sensitivity.append(check)

    coupling_profile = fixed_coupling_profile(data, chol, y_white, fit["Gamma_over_H0"], fit["chi2"])
    curvature_scales = [local_profile_curvature(data, chol, y_white, fit["Omega_m"], fit["Gamma_over_H0"], 1e-4, 1e-3),
                        local_profile_curvature(data, chol, y_white, fit["Omega_m"], fit["Gamma_over_H0"], 2e-4, 2e-3)]
    nested_coupling_windows = []
    for g_bounds in ((-1.0, 1.0), (-0.5, 0.5)):
        window_fit = run_fit(data, chol, y_white, gamma_fixed=None, gamma_bounds=g_bounds, starts=25)
        window_fit["delta_chi2_vs_flat"] = float(flat["chi2"] - window_fit["chi2"])
        window_fit["delta_chi2_vs_full_box_best"] = float(window_fit["chi2"] - fit["chi2"])
        window_fit["coupling_boundary_hit"] = bool(any(
            name in window_fit["parameter_bound_hits"] for name in ("Gamma_over_H0:lower", "Gamma_over_H0:upper")))
        nested_coupling_windows.append({"Gamma_over_H0_bounds": g_bounds, "fit": window_fit,
                                        "minimum_acceptance": "selected finite objective (< INVALID) with a fully physical integrated state; raw SciPy statuses, including any invalid penalty endpoints, are retained"})

    # Sign convention and dimensional normalization checks.
    zsign = 1.0
    qplus, bgplus = unit_prediction(data, 0.3, 0.2)
    qzero, bgzero = unit_prediction(data, 0.3, 0.0)
    qminus, bgminus = unit_prediction(data, 0.3, -0.2)
    e_z1 = lambda bg, g: float(integrate_profile(0.3, g, np.array([zsign]))["e"][0])
    # alpha unit example: H0=70 km/s/Mpc, rd=147 Mpc.
    mpc_m = 648000.0/math.pi * 149597870700.0
    c_si = 299792458.0
    h0_si = 70.0*1000.0/mpc_m
    rd_m = 147.0*mpc_m
    alpha_unit = c_si/(h0_si*rd_m)
    alpha_rescaled = c_si/((2*h0_si)*(rd_m/2.0))

    # Boundary slices diagnose whether a nominal box face is hiding a better score.
    boundary_checks = []
    for gb in GAMMA_BOUNDS:
        vals = []
        for om in np.linspace(*OM_BOUNDS, 41):
            try: vals.append((profile_alpha(data, chol, y_white, float(om), float(gb))[0], float(om)))
            except DomainError: pass
        boundary_checks.append({"Gamma_over_H0": gb, "minimum_score": min(vals)[0] if vals else None,
                                "Omega_m_at_minimum": min(vals)[1] if vals else None})

    data_records = [{"z": float(z), "observable": str(o), "value": float(v)}
                    for z, o, v in zip(data.z, data.observable, data.value)]
    files = ["context/data/desi_dr2_mean.txt", "context/data/desi_dr2_cov.txt",
             "scripts/background_bao.py", "experiments/campaign_seed_bao/result.json",
             "context/papers/interacting_de_desi_dr2_2026.pdf", "BRIEF.md",
             "experiments/interacting_vacuum_screen/interacting_vacuum_profile.py"]
    hashes = {f: sha256(ROOT/f) for f in files}
    result = {
        "status": "complete_interacting_vacuum_late_time_profile_screen",
        "generated_unix_utc": time.time(),
        "scope": "13-row full-covariance DESI DR2 compressed BAO; local Q=Gamma*rho_x, w_x=-1; late-time matter+vacuum reduction; no early-time calibration/perturbations/posterior/evidence",
        "contract": {
            "paper_equations": "Eqs. (3)-(5),(8), pp.2-4; dot rho_x=Gamma rho_x and dot rho_c+3H rho_c=-Gamma rho_x for w_x=-1; positive Gamma is CDM -> DE",
            "late_time_reduction": "Radiation and neutrino terms in paper Eq.(3) explicitly omitted following BRIEF.md lane A and background_bao.py low-z screening contract; not a full paper-model reproduction",
            "state_equations_u_ln1pz": ["dX/du=-g*X/E", "dM/du=3*M+g*X/E", "E^2=M+X", "dD/du=exp(u)/E"],
            "initial_conditions": {"M0": "Omega_b0+Omega_c0=Omega_m0", "X0": "1-Omega_m0", "flat_closure": True},
            "bounds_are_search_domain_not_prior": {"Omega_m0": OM_BOUNDS, "Gamma_over_H0": GAMMA_BOUNDS, "alpha": ALPHA_BOUNDS, "source": "Omega_m and alpha from scripts/background_bao.py FIT_BOUNDS; Gamma/H0 from paper Table I IVS box"},
            "baryon_split": "No CMB/BBN baryon prior; f_b is profiled only by physical feasibility over 0<f_b<f_b_max, where f_b_max=min(1,min_z M(z)/(Omega_m*(1+z)^3)); no background BAO sensitivity to f_b because E depends only on B+C. f_b is not inferred.",
            "physical_domain": "Reject solver failures, nonfinite values, E^2<=0, M<=0, X<=0, or no nonempty positive baryon/CDM split on 0<=z<=2.33",
            "alpha_profile": "bounded GLS projection of unit-amplitude prediction using full-covariance Cholesky-whitened residual",
        },
        "inputs": {"mean_row_order": data_records, "covariance_shape": list(data.covariance.shape),
                   "covariance_cholesky_success": True, "covariance_min_eigenvalue": float(np.linalg.eigvalsh(data.covariance).min()),
                   "sha256": hashes},
        "fits": {"Gamma0_flat_LCDM": flat, "interacting_vacuum": fit},
        "flat_limit_check": flat_limit,
        "profile_comparison": {"delta_chi2_flat_minus_interacting": float(flat["chi2"]-fit["chi2"]),
                               "Gamma0_is_nested_limit": True, "interpretation": "finite profile improvement only; no statistical significance or model-selection claim"},
        "search_checks": {"coarse_grid_and_refinement": grid_refined, "boundary_slices_Gamma_pm3": boundary_checks,
                          "nested_coupling_windows": nested_coupling_windows,
                          "fixed_coupling_profile": coupling_profile,
                          "local_profile_curvature_correlation": {"scales": curvature_scales,
                              "not_a_posterior_interval": True},
                          "independent_adaptive_quad_Powell": {"x": [float(v) for v in powell.x], "chi2": float(powell.fun),
                              "success": bool(powell.success), "status": int(powell.status), "message": str(powell.message),
                              "nfev": int(powell.nfev), "domain_failures": dict(adaptive_fail),
                              "score_at_primary_optimum": adaptive_score, "alpha_at_primary_optimum": adaptive_alpha,
                              "max_abs_prediction_difference_vs_main": float(np.max(np.abs(adaptive_pred-fit["prediction"]))),
                              "delta_score_vs_main": float(adaptive_score-fit["chi2"])},
                          "tight_ODE_tolerance": {"rtol": 2e-12, "atol": 2e-14, "chi2": score_tight,
                              "alpha": alpha_tight, "delta_score_vs_main": float(score_tight-fit["chi2"]),
                              "max_abs_prediction_difference_vs_main": float(np.max(np.abs(pred_tight-fit["prediction"])))}},
        "synthetic_checks": {"truth": truth, "recovered": {k:synth[k] for k in ("Omega_m", "Gamma_over_H0", "alpha", "chi2", "selected_success")},
                             "optimizer_endpoint_audit": synth["optimizer_endpoint_audit"],
                             "baryon_split_examples": split_checks,
                             "max_unit_prediction_difference_repeated_split": split_invariance,
                             "fixed_split_sensitivity_at_observed_optimum": fixed_split_sensitivity,
                             "interpretation": "fixed f_b changes only component bookkeeping; when physically feasible, the exact BAO prediction and score are invariant"},
        "sign_and_units": {"source_at_today_positive_Gamma": {"Gamma_rho_x_sign": "positive", "CDM_source_minus_Gamma_rho_x_sign": "negative"},
                           "E_z1_gamma_minus_0p2": e_z1(bgminus,-0.2), "E_z1_gamma_0": e_z1(bgzero,0.0),
                           "E_z1_gamma_plus_0p2": e_z1(bgplus,0.2), "positive_gamma_sign_check": e_z1(bgplus,0.2)>e_z1(bgzero,0.0),
                           "negative_gamma_sign_check": e_z1(bgminus,-0.2)<e_z1(bgzero,0.0),
                           "alpha_dimensionless_example_H0_70_rd_147Mpc": alpha_unit,
                           "alpha_same_product_rescaling": alpha_rescaled, "alpha_rescaling_abs_difference": abs(alpha_unit-alpha_rescaled),
                           "alpha_unit_contract": "alpha=c/(H0*rd), H0 in s^-1 and rd in m; all returned BAO observables dimensionless"},
        "runtime": {"seconds": time.time()-started, "python": sys.version, "numpy": np.__version__,
                    "scipy": __import__("scipy").__version__, "platform": platform.platform(),
                    "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --seconds 1400 -- .venv/bin/python experiments/interacting_vacuum_screen/interacting_vacuum_profile.py --coarse-grid 41 --starts 39",
                    "workers": 1, "cpu_only": True, "gpu_used": False,
                    "blas_threads": {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS")}},
        "claim_limitations": ["compressed-BAO profile score only; no posterior/evidence/significance", "no rd or H0 inference because alpha is free",
                              "baryon/CDM present-day split exactly unidentifiable from this background distance likelihood", "no early-time calibration, radiation/neutrino transfer, or perturbations", "not a reproduction of published CMB+DESI+SNIa constraints"],
    }
    out = ROOT/"experiments/interacting_vacuum_screen/result.json"
    grid_out = ROOT/"experiments/interacting_vacuum_screen/profile_grid.npz"
    np.savez_compressed(grid_out, omega_m=oms, gamma_over_h0=gs, chi2=grid)
    # A compact profile surface is useful for auditing broad ridges and invalid regions.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
        zplot = np.where(np.isfinite(grid), grid - min(fit["chi2"], flat["chi2"]), np.nan)
        im = ax.pcolormesh(gs, oms, zplot, shading="auto", cmap="viridis", vmin=0, vmax=8)
        ax.contour(gs, oms, zplot, levels=[0.25,0.5,1,2,4], colors="white", linewidths=.7)
        ax.scatter([fit["Gamma_over_H0"]], [fit["Omega_m"]], marker="*", s=130, c="red", edgecolors="white", label="interacting best profile")
        ax.axvline(0,color="white",lw=1,ls="--",label="flat ΛCDM limit")
        ax.set(xlabel=r"$\Gamma/H_0$ (positive: CDM $\to$ DE)", ylabel=r"$\Omega_{m0}$", title=r"Profile $\Delta\chi^2$; full DESI DR2 BAO covariance")
        fig.colorbar(im, ax=ax, label=r"$\chi^2-\min(\chi^2_{flat},\chi^2_{int})$")
        ax.legend(loc="best", fontsize=8)
        for ext in ("png","svg"):
            fig.savefig(ROOT/f"experiments/interacting_vacuum_screen/profile_surface.{ext}", dpi=160)
        plt.close(fig)
    except ImportError:
        result["plot_note"] = "matplotlib unavailable; profile grid saved in NPZ"
    result["artifacts"] = {"profile_grid_npz_sha256": sha256(grid_out), "plots": {}}
    for ext in ("png", "svg"):
        plot = ROOT/f"experiments/interacting_vacuum_screen/profile_surface.{ext}"
        if plot.exists(): result["artifacts"]["plots"][plot.name] = sha256(plot)
    result["code_sha256"] = hashes["experiments/interacting_vacuum_screen/interacting_vacuum_profile.py"]
    out.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({"status":result["status"], "flat_chi2":flat["chi2"], "interacting_chi2":fit["chi2"],
                      "delta_chi2":flat["chi2"]-fit["chi2"], "gamma":fit["Gamma_over_H0"],
                      "omega_m":fit["Omega_m"], "elapsed":time.time()-started}, indent=2))


if __name__ == "__main__":
    main()

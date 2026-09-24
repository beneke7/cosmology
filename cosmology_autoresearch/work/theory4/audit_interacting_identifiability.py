#!/usr/bin/env python3
"""Independent scalar-distance audit of the interacting-vacuum BAO profile.

This is deliberately self-contained with respect to the experiment code: it reads
the DESI text files directly, integrates the background afresh, computes each
radial distance by adaptive scalar quadrature, and profiles the BAO amplitude
with the complete covariance.  All output stays in work/theory4/.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize, minimize_scalar


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "work/theory4"
ZMAX = 2.33
OM_BOUNDS = (0.05, 0.60)
G_BOUNDS = (-3.0, 3.0)
ALPHA_BOUNDS = (1e-6, 1e4)
FB_ILLUSTRATIVE = 0.16
RTOL = 2e-12
ATOL = 2e-14


class InvalidHistory(ValueError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_data():
    raw = np.loadtxt(ROOT / "context/data/desi_dr2_mean.txt", dtype=str, skiprows=1)
    z = raw[:, 0].astype(float)
    y = raw[:, 1].astype(float)
    labels = raw[:, 2].tolist()
    cov = np.loadtxt(ROOT / "context/data/desi_dr2_cov.txt", dtype=float)
    if z.size != 13 or cov.shape != (z.size, z.size):
        raise ValueError(f"expected 13 BAO rows and 13x13 covariance, got {z.size}, {cov.shape}")
    return z, y, labels, cov


def scalar_history(omega_m: float, g: float, zmax: float = ZMAX):
    """Return dense float64 (M,X)(u) solution, where u=ln(1+z)."""
    om, gg = float(omega_m), float(g)
    if not (OM_BOUNDS[0] <= om <= OM_BOUNDS[1] and G_BOUNDS[0] <= gg <= G_BOUNDS[1]):
        raise InvalidHistory("outside locked finite profile domain")
    if not (0.0 <= zmax <= ZMAX):
        raise InvalidHistory("redshift exceeds screen interval")
    x0 = 1.0 - om
    if om <= 0.0 or x0 <= 0.0:
        raise InvalidHistory("nonpositive present-day matter or vacuum")

    def rhs(u, state):
        matter, vacuum = state
        e2 = matter + vacuum
        if not np.isfinite(e2) or e2 <= 0.0:
            raise InvalidHistory("E^2 nonpositive during scalar ODE")
        e = math.sqrt(e2)
        transfer = gg * vacuum / e
        return (3.0 * matter + transfer, -transfer)

    umax = math.log1p(zmax)
    sol = solve_ivp(rhs, (0.0, umax), (om, x0), method="DOP853", dense_output=True,
                    rtol=RTOL, atol=ATOL, max_step=0.01)
    if not sol.success or sol.sol is None:
        raise InvalidHistory(f"DOP853 failed: {sol.message}")
    grid = np.linspace(0.0, umax, 257)
    hist = np.asarray(sol.sol(grid), dtype=np.float64)
    m, x = hist
    e2 = m + x
    if (not np.all(np.isfinite(hist)) or np.any(m <= 0.0) or np.any(x <= 0.0)
            or np.any(e2 <= 0.0)):
        raise InvalidHistory("nonpositive matter, vacuum, or expansion on sampled interval")
    return sol.sol, {
        "minimum_sampled_M": float(np.min(m)),
        "minimum_sampled_X": float(np.min(x)),
        "minimum_sampled_E2": float(np.min(e2)),
    }


def geometry_and_split(omega_m: float, g: float, z: np.ndarray, f_b: float = FB_ILLUSTRATIVE):
    """Compute E and DM by adaptive QUAD; independently test fixed baryon split."""
    sol, minima = scalar_history(omega_m, g, float(np.max(z)))

    def inverse_e(zp: float) -> float:
        u = math.log1p(float(zp))
        matter, vacuum = sol(u)
        e2 = float(matter + vacuum)
        if not math.isfinite(e2) or e2 <= 0.0:
            raise InvalidHistory("E^2 nonpositive in adaptive distance quadrature")
        return 1.0 / math.sqrt(e2)

    dm = np.asarray([
        quad(inverse_e, 0.0, float(zi), epsabs=5e-12, epsrel=5e-12, limit=200)[0]
        for zi in z
    ], dtype=np.float64)
    e = np.asarray([1.0 / inverse_e(float(zi)) for zi in z], dtype=np.float64)
    q = np.empty_like(dm)
    for i, label in enumerate(LABELS):
        if label == "DM_over_rs":
            q[i] = dm[i]
        elif label == "DH_over_rs":
            q[i] = 1.0 / e[i]
        elif label == "DV_over_rs":
            q[i] = np.cbrt(z[i] * dm[i] * dm[i] / e[i])
        else:
            raise ValueError(f"unknown observable {label!r}")

    # q_c=a^3 rho_c obeys dq_c/du=g*a^3*X/E.  Its monotonicity means
    # the minimum on [0,zmax] is at zmax for g<0, at z=0 for g>0.
    # Check the exact endpoint that can be limiting, plus a dense history check.
    zcheck = np.linspace(0.0, max(float(np.max(z)), 0.0), 1001)
    ucheck = np.log1p(zcheck)
    mcheck, xcheck = np.asarray(sol(ucheck), dtype=np.float64)
    bcheck = f_b * omega_m * np.exp(3.0 * ucheck)
    ccheck = mcheck - bcheck
    ratio = mcheck / (omega_m * np.exp(3.0 * ucheck))
    if g < 0.0:
        limiting_z = float(zcheck[-1])
    else:
        limiting_z = 0.0
    limiting_c = float(np.interp(limiting_z, zcheck, ccheck))
    split_ok = bool(np.all(bcheck > 0.0) and limiting_c > 0.0 and np.min(ccheck) > 0.0)
    if not split_ok:
        raise InvalidHistory("fixed illustrative baryon split has nonpositive CDM")
    return q, {
        **minima,
        "fb_max_sampled": float(min(1.0, np.min(ratio))),
        "fb_fixed": float(f_b),
        "fixed_split_positive": split_ok,
        "minimum_baryon_density_fraction": float(np.min(bcheck)),
        "minimum_CDM_density_fraction": float(np.min(ccheck)),
        "minimum_CDM_redshift": float(zcheck[int(np.argmin(ccheck))]),
        "limiting_comoving_CDM_density_at_zmax_or_today": limiting_c,
        "split_monotonic_endpoint": "zmax" if g < 0.0 else "today",
    }


def make_profile_function(z, y, chol, ywhite, cov):
    cache: dict[tuple[float, float], tuple[float, float, np.ndarray, dict]] = {}

    def evaluate(omega_m: float, g: float):
        key = (float(omega_m), float(g))
        if key in cache:
            return cache[key]
        try:
            q, phys = geometry_and_split(float(omega_m), float(g), z)
            qw = np.linalg.solve(chol, q)
            alpha = float(np.clip(np.dot(qw, ywhite) / np.dot(qw, qw), *ALPHA_BOUNDS))
            residual = ywhite - alpha * qw
            score = float(np.dot(residual, residual))
            cache[key] = (score, alpha, q, phys)
            return cache[key]
        except (InvalidHistory, ValueError, FloatingPointError, OverflowError):
            bad = (float("inf"), float("nan"), np.full_like(y, np.nan), {})
            cache[key] = bad
            return bad

    def omega_profile(g: float, n_omega: int):
        grid = np.linspace(*OM_BOUNDS, int(n_omega))
        vals = np.asarray([evaluate(float(om), float(g))[0] for om in grid])
        finite = np.flatnonzero(np.isfinite(vals))
        if not finite.size:
            return {"g": float(g), "valid": False, "chi2": None, "omega_m": None,
                    "alpha": None, "fb_max": None, "grid_valid_count": 0}
        ibest = int(finite[np.argmin(vals[finite])])
        best_om, best_chi = float(grid[ibest]), float(vals[ibest])
        # Refine only the grid cell around the best physically admissible point.
        left = float(grid[max(ibest - 1, 0)])
        right = float(grid[min(ibest + 1, grid.size - 1)])
        refine_bracket_ok = (
            right > left
            and np.isfinite(evaluate(left, float(g))[0])
            and np.isfinite(evaluate(right, float(g))[0])
        )
        if refine_bracket_ok:
            opt = minimize_scalar(lambda om: evaluate(float(om), float(g))[0],
                                  bounds=(left, right), method="bounded",
                                  options={"xatol": 2e-10, "maxiter": 160})
            if np.isfinite(opt.fun) and float(opt.fun) < best_chi:
                best_om, best_chi = float(opt.x), float(opt.fun)
        score, alpha, _q, phys = evaluate(best_om, float(g))
        return {"g": float(g), "valid": True, "chi2": float(score),
                "omega_m": float(best_om), "alpha": float(alpha),
                "fb_max": phys.get("fb_max_sampled"),
                "min_cdm": phys.get("minimum_CDM_density_fraction"),
                "grid_valid_count": int(finite.size)}

    return evaluate, omega_profile, cache


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gamma-grid", type=int, default=61)
    parser.add_argument("--omega-grid", type=int, default=41)
    parser.add_argument("--fine-gamma-grid", type=int, default=31)
    args = parser.parse_args()
    if args.gamma_grid < 21 or args.omega_grid < 21 or args.fine_gamma_grid < 11:
        parser.error("use gamma-grid>=21, omega-grid>=21, and fine-gamma-grid>=11")

    start = time.monotonic()
    global LABELS
    z, y, LABELS, cov = load_data()
    chol = np.linalg.cholesky(cov)
    ywhite = np.linalg.solve(chol, y)
    evaluate, omega_profile, cache = make_profile_function(z, y, chol, ywhite, cov)

    coarse_g = np.linspace(*G_BOUNDS, args.gamma_grid)
    coarse = [omega_profile(float(g), args.omega_grid) for g in coarse_g]
    valid_coarse = [r for r in coarse if r["valid"]]
    if not valid_coarse:
        raise RuntimeError("no fixed-split valid point on coarse gamma grid")
    coarse_best = min(valid_coarse, key=lambda r: r["chi2"])
    fine_g = np.linspace(max(G_BOUNDS[0], coarse_best["g"] - 0.30),
                         min(G_BOUNDS[1], coarse_best["g"] + 0.30), args.fine_gamma_grid)
    fine = [omega_profile(float(g), args.omega_grid) for g in fine_g]
    compact_g = np.linspace(-0.20, 0.20, 21)
    compact_profile = [omega_profile(float(g), args.omega_grid) for g in compact_g]

    # Independent two-shape local refinement, still using scalar quadrature.
    init = np.asarray([coarse_best["omega_m"], coarse_best["g"]], dtype=np.float64)
    local = minimize(lambda p: evaluate(float(p[0]), float(p[1]))[0], init,
                     method="Nelder-Mead", bounds=(OM_BOUNDS, G_BOUNDS),
                     options={"xatol": 2e-9, "fatol": 5e-10, "maxiter": 600})
    om_hat, g_hat = map(float, local.x)
    score_hat, alpha_hat, q_hat, phys_hat = evaluate(om_hat, g_hat)
    if not np.isfinite(score_hat):
        raise RuntimeError("independent local shape refinement ended outside the physical domain")
    pred_hat = alpha_hat * q_hat

    # Refine the exact nested null and compare the scalar implementation to both
    # archived profile fits.  At g=0 this should reproduce ordinary flat LCDM.
    grid0 = np.linspace(*OM_BOUNDS, args.omega_grid)
    scores0 = np.asarray([evaluate(float(om), 0.0)[0] for om in grid0])
    i0 = int(np.nanargmin(scores0))
    lo0, hi0 = grid0[max(i0-1, 0)], grid0[min(i0+1, len(grid0)-1)]
    opt0 = minimize_scalar(lambda om: evaluate(float(om), 0.0)[0], bounds=(lo0, hi0),
                           method="bounded", options={"xatol": 2e-10, "maxiter": 160})
    score0, alpha0, q0, phys0 = evaluate(float(opt0.x), 0.0)

    # Profile likelihood landmarks are diagnostics under this finite search box,
    # not posterior intervals. Linear interpolation is used only to locate crossings.
    curve_by_g = {round(r["g"], 12): r for r in coarse + fine + compact_profile}
    curve_by_g[round(g_hat, 12)] = {
        "g": g_hat, "valid": True, "chi2": float(score_hat), "omega_m": om_hat,
        "alpha": alpha_hat, "fb_max": phys_hat.get("fb_max_sampled"),
        "min_cdm": phys_hat.get("minimum_CDM_density_fraction"),
        "grid_valid_count": None, "point_type": "independent two-parameter local refinement",
    }
    curve = [curve_by_g[k] for k in sorted(curve_by_g)]
    mincurve = min((r for r in curve if r["valid"]), key=lambda r: r["chi2"])
    curve_min = float(score_hat)

    def crossing_intervals(delta: float):
        vals = [(r["g"], r["chi2"] - curve_min - delta) for r in curve if r["valid"]]
        g0 = float(mincurve["g"])
        left_cross = right_cross = None
        left = [(g, v) for g, v in vals if g < g0]
        right = [(g, v) for g, v in vals if g > g0]
        for seq, side in ((list(reversed(left)), "left"), (right, "right")):
            previous = (g0, -delta)
            for cur in seq:
                if cur[1] >= 0.0 and previous[1] <= 0.0:
                    frac = -previous[1] / (cur[1] - previous[1])
                    cross = previous[0] + frac * (cur[0] - previous[0])
                    if side == "left":
                        left_cross = float(cross)
                    else:
                        right_cross = float(cross)
                    break
                previous = cur
        return {"delta_chi2": delta, "left": left_cross, "right": right_cross,
                "classification": "profile-curve diagnostic only; not posterior or calibrated confidence limits"}

    # Local shape Jacobian. Whiten by L^{-1}; project out the free alpha direction.
    # Derivatives are checked at two centered step sizes for local numerical stability.
    def q_at(om, gg):
        return geometry_and_split(float(om), float(gg), z)[0]

    def jacobian_projected(step_scale: float):
        h_om = 2e-5 * step_scale
        h_g = 2e-5 * step_scale
        dq = np.column_stack((
            (q_at(om_hat+h_om, g_hat) - q_at(om_hat-h_om, g_hat))/(2*h_om),
            (q_at(om_hat, g_hat+h_g) - q_at(om_hat, g_hat-h_g))/(2*h_g),
        ))
        qwhite = np.linalg.solve(chol, q_hat)
        dqw = np.linalg.solve(chol, dq)
        unit = qwhite / np.linalg.norm(qwhite)
        projector = np.eye(len(y)) - np.outer(unit, unit)
        return alpha_hat * projector @ dqw

    jhalf = jacobian_projected(0.5)
    jnom = jacobian_projected(1.0)
    jdouble = jacobian_projected(2.0)
    u_svd, svals, vh = np.linalg.svd(jnom, full_matrices=False)
    rank = int(np.sum(svals > (svals[0] * 1e-10)))
    weak = vh[-1]
    fisher = jnom.T @ jnom
    fisher_cov = np.linalg.pinv(fisher, rcond=1e-12)
    corr = fisher_cov / np.sqrt(np.outer(np.diag(fisher_cov), np.diag(fisher_cov)))

    # This alternative is the exact finite-difference derivative of the profiled
    # whitened residual at the observed-data point; compare it with the standard
    # amplitude-projected mean Jacobian above.
    def profiled_residual(om, gg):
        score, aa, qq, _ = evaluate(float(om), float(gg))
        if not np.isfinite(score):
            raise InvalidHistory("profiled residual step outside physical domain")
        return np.linalg.solve(chol, y - aa*qq)

    h = 2e-5
    residual_j = np.column_stack((
        (profiled_residual(om_hat+h, g_hat)-profiled_residual(om_hat-h, g_hat))/(2*h),
        (profiled_residual(om_hat, g_hat+h)-profiled_residual(om_hat, g_hat-h))/(2*h),
    ))
    residual_svals = np.linalg.svd(residual_j, compute_uv=False)
    pred_diff = np.max(np.abs(pred_hat - np.asarray(
        json.loads((ROOT/"experiments/interacting_vacuum_screen/result.json").read_text())["fits"]["interacting_vacuum"]["prediction"])))
    archive = json.loads((ROOT / "experiments/interacting_vacuum_screen/result.json").read_text())
    reported = archive["fits"]["interacting_vacuum"]
    baseline_reported = archive["fits"]["Gamma0_flat_LCDM"]

    source_pdf = ROOT / "context/papers/interacting_de_desi_dr2_2026.pdf"
    current_code = ROOT / "experiments/interacting_vacuum_screen/interacting_vacuum_profile.py"
    current_result = ROOT / "experiments/interacting_vacuum_screen/result.json"
    baseline_code = ROOT / "scripts/background_bao.py"
    mean_path = ROOT / "context/data/desi_dr2_mean.txt"
    cov_path = ROOT / "context/data/desi_dr2_cov.txt"
    contract_path = ROOT / "work/theory3/interacting_vacuum_contract.md"
    report_path = ROOT / "work/compute3/interacting_vacuum_screen_report.md"

    result = {
        "classification": "independent numerical audit; local likelihood diagnostics, not posterior inference",
        "source_equation_audit": {
            "paper": "Do DESI-DR2 BAO data imply a coupling of dark matter and dark energy?",
            "local_pdf": str(source_pdf.relative_to(ROOT)),
            "pages_equations": "Eq. (3), p.2 for separate baryons/CDM/DE densities; Eqs. (4)-(5), p.2 for continuity; Eq. (8), p.3 for Q=Gamma*rho_x; text below Eq. (8) for sign",
            "equations_used": [
                "dot(rho_x)=Gamma*rho_x for w_x=-1",
                "dot(rho_c)+3H*rho_c=-Gamma*rho_x",
                "dot(rho_b)+3H*rho_b=0",
                "Q=Gamma*rho_x; Gamma>0 transfers CDM to DE",
                "u=ln(1+z), dX/du=-g X/E, dM/du=3M+gX/E, E^2=M+X, g=Gamma/H0",
            ],
            "consistency_notes": [
                "The sign of dM/du follows converting the paper's proper-time CDM continuity equation using du/dt=-H.",
                "M and X are fractions of today's critical density; E and u are dimensionless.",
                "Baryon and CDM densities sum to M. Separate baryon conservation is applied only in the fixed-split physicality check.",
                "Radiation/neutrinos and perturbations are omitted by the declared low-redshift BAO screen; no rd or H0 calibration is inferred.",
            ],
            "screen_contract_mismatch_to_note": "The older work/theory3 contract proposed a compact exploratory gamma window [-0.2,0.2], while the compute screen locked [-3,3] to the source paper IVS scan box (Table I); the independently audited fit at g=-0.4667 lies outside the compact initial sketch but inside the later recorded search domain.",
        },
        "inputs": {
            "data_rows": len(z), "covariance_shape": list(cov.shape),
            "covariance_min_eigenvalue": float(np.linalg.eigvalsh(cov)[0]),
            "covariance_cholesky_ok": True,
            "mean_order": [{"z": float(zi), "value": float(yi), "observable": str(lab)} for zi, yi, lab in zip(z,y,LABELS)],
            "sha256": {
                "source_paper_pdf": sha256(source_pdf),
                "mean": sha256(mean_path), "covariance": sha256(cov_path),
                "prior_experiment_script": sha256(current_code),
                "prior_experiment_result": sha256(current_result),
                "baseline_script": sha256(baseline_code),
                "theory_contract": sha256(contract_path),
                "compute_report": sha256(report_path),
            },
        },
        "profile": {
            "likelihood": "full 13x13 covariance; Cholesky whitening; bounded analytic GLS alpha profile; fixed illustrative f_b=0.16 positivity restriction",
            "bounds_search_domain_not_prior": {"omega_m": list(OM_BOUNDS), "g_gamma_over_H0": list(G_BOUNDS), "alpha": list(ALPHA_BOUNDS)},
            "fixed_baryon_fraction": FB_ILLUSTRATIVE,
            "independent_local_best": {"omega_m": om_hat, "g": g_hat, "alpha": alpha_hat,
                                        "chi2": float(score_hat), "delta_chi2_vs_flat": float(score0-score_hat),
                                        "optimizer_success": bool(local.success), "optimizer_message": str(local.message),
                                        "optimizer_nit": int(local.nit), "optimizer_nfev": int(local.nfev)},
            "flat_nested_null": {"omega_m": float(opt0.x), "g": 0.0, "alpha": float(alpha0), "chi2": float(score0)},
            "comparison_to_archived_screen": {
                "archived_omega_m": float(reported["Omega_m"]), "archived_g": float(reported["Gamma_over_H0"]),
                "archived_alpha": float(reported["alpha"]), "archived_chi2": float(reported["chi2"]),
                "abs_chi2_difference": float(abs(score_hat-reported["chi2"])),
                "max_prediction_difference": float(pred_diff),
                "archived_flat_chi2": float(baseline_reported["chi2"]),
            },
            "profile_curve": curve,
            "compact_initial_window_sensitivity": {
                "g_window": [-0.20, 0.20],
                "status": "initial theory contract's exploratory window, profiled separately on a 21-point g mesh with the same omega_m and alpha treatment",
                "best_grid_profile": min((r for r in compact_profile if r["valid"]), key=lambda r: r["chi2"]),
                "best_is_window_boundary": bool(abs(min((r for r in compact_profile if r["valid"]), key=lambda r: r["chi2"])["g"] + 0.20) < 1e-12),
                "profile_points": compact_profile,
            },
            "profile_delta_chi2_landmarks": [crossing_intervals(1.0), crossing_intervals(3.84)],
            "landmark_warning": "These finite-box profile crossings are descriptive only. Wilks thresholds are not asserted as calibrated confidence intervals; the screen excludes early physics and the candidate emerged within a broader adaptive campaign.",
            "fixed_split_physicality_at_best": phys_hat,
        },
        "local_identifiability_at_independent_best": {
            "definition": "J = alpha * (I - vv^T) * L^{-1} dq/d(omega_m,g), v is unit L^{-1}q; alpha amplitude direction projected out; columns are derivatives per unit dimensionless omega_m and g.",
            "finite_difference_steps": {"omega_m": [1e-5,2e-5,4e-5], "g": [1e-5,2e-5,4e-5]},
            "projected_shape_singular_values": svals.tolist(),
            "numerical_rank_relative_threshold_1e-10": rank,
            "singular_value_ratio_small_over_large": float(svals[-1]/svals[0]),
            "weak_shape_direction_omega_m_g": weak.tolist(),
            "weak_direction_delta_g_per_delta_omega_m": float(weak[1]/weak[0]) if abs(weak[0]) > 1e-14 else None,
            "local_linearized_shape_correlation": corr.tolist(),
            "finite_difference_sensitivity_max_abs_by_column": {
                "half_vs_nominal": np.max(np.abs(jhalf-jnom), axis=0).tolist(),
                "double_vs_nominal": np.max(np.abs(jdouble-jnom), axis=0).tolist(),
            },
            "exact_profiled_residual_jacobian_singular_values": residual_svals.tolist(),
            "interpretation": "A full local rank means the two shape parameters are locally distinguishable after allowing the common BAO amplitude to float; a small singular-value ratio indicates the weak compensating shape combination. Singular values and the local linearized covariance are diagnostics, not posterior constraints.",
        },
        "reproduction": {
            "command_from_project_root": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --seconds 1200 -- .venv/bin/python work/theory4/audit_interacting_identifiability.py --gamma-grid 61 --omega-grid 41 --fine-gamma-grid 31",
            "python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
            "platform": platform.platform(), "cpu_affinity_count": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
            "cpu_threads_requested": 1, "gpu": "not used; scalar ODE and adaptive quadrature reference is CPU-only",
            "wall_seconds": float(time.monotonic()-start),
            "script_sha256": sha256(Path(__file__)),
            "output_json": str((OUT/"interacting_vacuum_identifiability.json").relative_to(ROOT)),
            "limits": [
                "Profile is limited to the released 13-row compressed DESI DR2 BAO vector and its full covariance.",
                "No posterior sampling, Bayesian evidence, CMB/BBN calibration, early-time sound horizon, growth, perturbation closure, or covariant transfer-vector analysis.",
                "f_b=0.16 is an illustrative fixed split used only to test positivity, not a measured prior; BAO distances still have exact baryon/CDM split degeneracy.",
                "The selected paper IVS parameter box [-3,3] is a finite search region, not a probability prior.",
                "Profile crossings are not confidence intervals; null calibration remains separate.",
            ],
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    json_path = OUT / "interacting_vacuum_identifiability.json"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    # The report writer below consumes the saved result after this script completes.
    print(json.dumps({"output": str(json_path), "chi2": score_hat, "g": g_hat,
                      "omega_m": om_hat, "delta_chi2": score0-score_hat,
                      "singular_values": svals.tolist(), "profile_points": len(curve),
                      "wall_seconds": time.monotonic()-start}, indent=2))


if __name__ == "__main__":
    main()

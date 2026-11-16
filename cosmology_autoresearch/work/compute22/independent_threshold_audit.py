#!/usr/bin/env python3
"""Independent event-based audit of the compute19 benchmark threshold.

This implementation uses n=ln(a) and the state (x,q), where x=rho_x/rho_crit0
and q=a^3*rho_c/rho_crit0. It does not import or call the compute19 producer
or rerun driver. The saved compute19 JSON is read only as a comparison target.
"""
from __future__ import annotations

import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"
if hasattr(os, "sched_getaffinity"):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})

import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp


HERE = Path(__file__).resolve().parent
PROJECT = Path(__file__).resolve().parents[2]
PDF = PROJECT / "context/papers/interacting_de_desi_dr2_2026.pdf"
FIX_REPORT = PROJECT / "work/compute19/fix/physicality_rerun_report.md"
SAVED = PROJECT / "work/compute19/fix/physicality_result.json"
OUT = HERE / "independent_threshold_audit.json"

# Exact benchmark convention declared in compute19's fix report/JSON.
H0 = 67.4                         # km s^-1 Mpc^-1
HUBBLE_h = H0 / 100.0
TCMB = 2.7255                     # K
NEFF = 3.044                      # massless effective species
OMEGA_M = 0.315
OMEGA_B = 0.049
OMEGA_C = OMEGA_M - OMEGA_B
OMEGA_GAMMA_H2 = 2.469e-5         # fix report's stated convention
OMEGA_GAMMA = OMEGA_GAMMA_H2 / HUBBLE_h**2 * (TCMB / 2.7255)**4
NU_RATIO = (7.0 / 8.0) * (4.0 / 11.0)**(4.0 / 3.0) * NEFF
OMEGA_NU = OMEGA_GAMMA * NU_RATIO
OMEGA_RAD = OMEGA_GAMMA + OMEGA_NU
OMEGA_X = 1.0 - OMEGA_B - OMEGA_C - OMEGA_RAD
ZMAX = 1.0e10
N_MIN = -math.log1p(ZMAX)
U_MAX = -N_MIN
G_SELECTED = (-3.0, -1.25, -1.0, 0.0, 3.0)
N_SAMPLES = 4097
N_SAMPLE_GRID = np.linspace(0.0, N_MIN, N_SAMPLES)  # backward-time order
CFG = {"method": "RK45", "rtol": 3e-13, "atol": 2e-15, "max_step": 0.02}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def friedmann(n: float, x: float, q: float) -> float:
    """Flat Friedmann relation in independent comoving-CDM variables."""
    # rho_c = q*a^-3 = q*exp(-3n).
    return (OMEGA_RAD * math.exp(-4.0 * n)
            + OMEGA_B * math.exp(-3.0 * n)
            + q * math.exp(-3.0 * n) + x)


def crossing_event(n: float, y: np.ndarray) -> float:
    return float(y[1])

crossing_event.terminal = True
crossing_event.direction = 0


def integrate(g: float, samples: bool = True) -> dict:
    """Integrate x and q backward from a=1; stop at q=0 if reached."""
    def rhs(n: float, y: np.ndarray) -> tuple[float, float]:
        x, q = float(y[0]), float(y[1])
        e2 = friedmann(n, x, q)
        if not math.isfinite(e2) or e2 <= 0.0:
            raise FloatingPointError(f"E^2 invalid for g={g} at n={n}: {e2}")
        e = math.sqrt(e2)
        dx_dn = g * x / e
        dq_dn = -g * math.exp(3.0 * n) * x / e
        return dx_dn, dq_dn

    sol = solve_ivp(
        rhs, (0.0, N_MIN), (OMEGA_X, OMEGA_C),
        events=crossing_event,
        dense_output=True,
        t_eval=N_SAMPLE_GRID if samples else None,
        **CFG,
    )
    if not sol.success:
        raise RuntimeError(f"integration failed for g={g}: {sol.message}")
    hit = bool(sol.t_events[0].size)
    n_end = float(sol.t_events[0][0]) if hit else N_MIN
    x_end = float(sol.y_events[0][0, 0]) if hit else float(sol.y[0, -1])
    q_end = float(sol.y_events[0][0, 1]) if hit else float(sol.y[1, -1])
    e2_end = friedmann(n_end, x_end, q_end)
    z_end = math.expm1(-n_end)

    if samples:
        # Restrict diagnostics to the integrated interval, ending on the event.
        n_grid = np.linspace(0.0, n_end, N_SAMPLES)
        y_grid = sol.sol(n_grid)
        x_grid, q_grid = y_grid[0], y_grid[1]
        e2_grid = np.array([friedmann(float(n), float(x), float(q))
                            for n, x, q in zip(n_grid, x_grid, q_grid)])
        rc_grid = q_grid * np.exp(-3.0 * n_grid)
        x_min = float(np.min(x_grid))
        q_min = float(np.min(q_grid))
        rc_min = float(np.min(rc_grid))
        e2_min = float(np.min(e2_grid))
        positive = bool(np.all(x_grid > 0.0) and np.all(e2_grid > 0.0)
                        and np.all(q_grid >= -2e-13) and np.all(rc_grid >= -2e-13))
        if g < 0.0:
            q_u_monotonic = bool(np.all(np.diff(q_grid) <= 2e-13))
        elif g > 0.0:
            q_u_monotonic = bool(np.all(np.diff(q_grid) >= -2e-13))
        else:
            q_u_monotonic = bool(np.max(np.abs(q_grid - OMEGA_C)) < 1e-12)
    else:
        positive, q_u_monotonic = None, None
        x_min = q_min = rc_min = e2_min = None

    # Independent quadrature form of q(u)=Omega_c0+g*int_0^u x/E exp(-3v)dv.
    def q_integrand(n: float) -> float:
        x, q = (float(v) for v in sol.sol(n))
        e2 = friedmann(n, x, q)
        if e2 <= 0.0:
            raise FloatingPointError("nonpositive E2 in identity quadrature")
        return math.exp(3.0 * n) * x / math.sqrt(e2)

    integral, quad_error = quad(q_integrand, n_end, 0.0,
                                epsabs=2e-14, epsrel=2e-13, limit=500)
    q_from_integral = OMEGA_C + g * integral
    identity_abs = abs(q_end - q_from_integral)
    identity_rel = identity_abs / max(abs(q_end), OMEGA_C, abs(g * integral), 1e-300)

    return {
        "g": g,
        "crossed_zero_before_zmax": hit and (n_end > N_MIN + 1e-12),
        "event_or_endpoint": "q=0 event" if hit else "zmax endpoint",
        "n_end": n_end,
        "u_end": -n_end,
        "z_end": z_end,
        "x_end": x_end,
        "q_end": q_end,
        "rho_c_end": q_end * math.exp(-3.0 * n_end),
        "E2_end": e2_end,
        "x_min": x_min,
        "q_min": q_min,
        "omega_c_min": rc_min,
        "E2_min": e2_min,
        "q_from_independent_quad": q_from_integral,
        "q_identity_abs_residual": identity_abs,
        "q_identity_scaled_residual": identity_rel,
        "quadrature_reported_error": quad_error,
        "positivity_through_stop_sampled": positive,
        "q_monotonicity_in_u_sampled": q_u_monotonic,
        "sample_count": N_SAMPLES if samples else 0,
        "solver_nfev": int(sol.nfev),
    }


def find_threshold() -> dict:
    """Bisection on event/no-event classification, never continue past q=0."""
    lo, hi = -1.25, -1.0
    lo_result, hi_result = integrate(lo, samples=False), integrate(hi, samples=False)
    if not lo_result["crossed_zero_before_zmax"] or hi_result["crossed_zero_before_zmax"]:
        raise RuntimeError("initial threshold bracket does not straddle the first-zero event")
    history = []
    for i in range(35):
        mid = 0.5 * (lo + hi)
        current = integrate(mid, samples=False)
        history.append({"iteration": i + 1, "g": mid,
                        "crossed": current["crossed_zero_before_zmax"],
                        "n_end": current["n_end"], "q_end": current["q_end"]})
        if current["crossed_zero_before_zmax"]:
            lo, lo_result = mid, current
        else:
            hi, hi_result = mid, current
    # Retain sampled positivity and continuity checks for the final bracket sides.
    lo_result, hi_result = integrate(lo, samples=True), integrate(hi, samples=True)
    return {
        "crossing_side_g": lo,
        "survives_to_zmax_side_g": hi,
        "bracket_width": hi - lo,
        "crossing_side": lo_result,
        "survives_side": hi_result,
        "iterations": len(history),
        "bisection_history": history,
    }


def compare_to_report(independent: dict, saved: dict) -> dict:
    answer = {}
    pairs = (
        ("crossed_zero_before_zmax", "crossed_zero"),
        ("u_end", "u_end"), ("z_end", "z_end"),
        ("q_end", "q_end"), ("x_min", "x_min"),
        ("q_min", "q_min"), ("omega_c_min", "omega_c_min"),
        ("E2_min", "e2_min"),
    )
    for own_key, saved_key in pairs:
        left, right = independent[own_key], saved[saved_key]
        if isinstance(left, bool):
            answer[own_key] = {"independent": left, "saved": right, "match": left == right}
        else:
            scale = max(abs(left), abs(right))
            answer[own_key] = {"independent": left, "saved": right,
                           "absolute_difference": abs(left-right),
                           "relative_difference": abs(left-right)/scale if scale >= 1e-10 else None,
                           "near_zero_absolute_only": scale < 1e-10}
    return answer


def main():
    start = time.perf_counter()
    saved = json.loads(SAVED.read_text())
    saved_cases = saved["numerics"]["coarse_scan"]
    selected = []
    for g in G_SELECTED:
        independent = integrate(g)
        reference = next(c for c in saved_cases if c["g"] == g)
        selected.append({**independent, "saved_case_comparison": compare_to_report(independent, reference)})

    threshold = find_threshold()
    report_threshold = saved["numerics"]["threshold"]
    report_g_cross = report_threshold["g_crossed_side"]
    report_g_survive = report_threshold["g_positive_through_zmax_side"]
    threshold_compare = {
        "independent_crossing_side": threshold["crossing_side_g"],
        "saved_crossing_side": report_g_cross,
        "crossing_side_absolute_difference": abs(threshold["crossing_side_g"] - report_g_cross),
        "independent_survives_side": threshold["survives_to_zmax_side_g"],
        "saved_survives_side": report_g_survive,
        "survives_side_absolute_difference": abs(threshold["survives_to_zmax_side_g"] - report_g_survive),
        "saved_bracket_width": report_threshold["parameter_bracket_width"],
        "independent_bracket_width": threshold["bracket_width"],
        "saved_crossing_u": report_threshold["last_crossing_u"],
        "saved_crossing_z": report_threshold["last_crossing_z"],
        "saved_survives_q_at_zmax": report_threshold["physical_side_q_at_zmax"],
        "independent_bracket_is_inside_saved_bracket": (
            report_g_cross <= threshold["crossing_side_g"]
            <= threshold["survives_to_zmax_side_g"] <= report_g_survive
        ),
        "independent_threshold_event": threshold["crossing_side"],
        "independent_threshold_survives": threshold["survives_side"],
        "crossing_event_u_absolute_difference": abs(
            threshold["crossing_side"]["u_end"] - report_threshold["last_crossing_u"]
        ),
        "bisection_iterations": threshold["iterations"],
    }
    event_case = threshold["crossing_side"]
    dq_du_abs = abs(event_case["g"] * math.exp(3.0 * event_case["n_end"])
                    * event_case["x_end"] / math.sqrt(event_case["E2_end"]))
    threshold_compare["event_location_conditioning"] = {
        "abs_dq_du_at_event": dq_du_abs,
        "q_identity_abs_residual": event_case["q_identity_abs_residual"],
        "estimated_u_sensitivity_scale_residual_over_abs_dqdu": (
            event_case["q_identity_abs_residual"] / dq_du_abs
        ),
        "interpretation": "q is nearly tangent to zero here; event redshift is much less accurately determined than the g classification bracket",
    }

    result = {
        "assessment": "independent benchmark-specific event-threshold reproduction; not a paper-prior or posterior result",
        "source": {
            "primary_pdf": str(PDF.relative_to(PROJECT)),
            "primary_pdf_sha256": sha(PDF),
            "source_equations": "flat Friedmann Eq. (3), separate conservation Eqs. (4)-(5), Q=Gamma*rho_x Eq. (8); local arXiv v1 PDF",
            "paper_fix_report": str(FIX_REPORT.relative_to(PROJECT)),
        },
        "benchmark": {
            "H0_km_s_Mpc": H0,
            "Tcmb_K": TCMB,
            "Neff_massless": NEFF,
            "Omega_m0": OMEGA_M,
            "Omega_b0": OMEGA_B,
            "Omega_c0": OMEGA_C,
            "Omega_gamma_h2_convention": OMEGA_GAMMA_H2,
            "Omega_gamma0": OMEGA_GAMMA,
            "Omega_nu0": OMEGA_NU,
            "Omega_rad0": OMEGA_RAD,
            "Omega_x0_flat_closure": OMEGA_X,
            "neutrino_assumption": "massless, matching compute19 fixed benchmark; not exact paper's 0.06 eV neutrino background",
            "threshold_definition": "first backward q=a^3*rho_c/rho_crit0=0 event at or before zmax=1e10",
            "paper_prior_warning": "fixed Omega_m0 benchmark only; do not infer posterior, prior mass after other parameters, or likelihood support",
        },
        "independent_derivation": {
            "n": "ln(a), with n=0 today and n_min=-ln(1+zmax)",
            "x": "rho_x/rho_crit0",
            "q": "a^3*rho_c/rho_crit0",
            "E2": "Omega_rad*exp(-4n)+Omega_b*exp(-3n)+q*exp(-3n)+x",
            "dx_dn": "g*x/E",
            "dq_dn": "-g*exp(3n)*x/E",
            "comoving_identity": "q(n)=Omega_c0+g*integral_n^0[exp(3v)*x(v)/E(v)]dv",
            "paper_sign": "Gamma>0 transfers CDM to vacuum forward in cosmic time; g=Gamma/H0",
            "event_policy": "terminal q=0 event; no solution is continued to negative CDM",
            "redshift_mapping": "u=-n=ln(1+z); z=expm1(-n)",
        },
        "numerics": {
            "method": CFG["method"],
            "rtol": CFG["rtol"],
            "atol": CFG["atol"],
            "max_step_in_n": CFG["max_step"],
            "event_bisection_iterations": 35,
            "sample_points_for_positivity_and_monotonicity": N_SAMPLES,
            "quadrature": "adaptive scipy.integrate.quad over dense output for comoving-q identity",
        },
        "selected_couplings": selected,
        "threshold": threshold_compare,
        "comparison_input": {
            "saved_result_json": str(SAVED.relative_to(PROJECT)),
            "saved_result_json_sha256": sha(SAVED),
            "saved_fix_report_sha256": sha(FIX_REPORT),
            "producer_code_not_imported_or_called": True,
        },
        "execution": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "cpu_affinity_ids": sorted(os.sched_getaffinity(0)),
            "gpu_used": False,
            "runtime_seconds_inside_script": None,
            "code_path": str(Path(__file__).relative_to(PROJECT)),
            "code_sha256": sha(Path(__file__).resolve()),
            "bounded_command": "taskset -c 0 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --state work/compute22/no_campaign_deadline.json --seconds 600 -- .venv/bin/python work/compute22/independent_threshold_audit.py",
            "model_effort": "GPT-6 family; exact variant and reasoning-effort tier not exposed by worker runtime (inherited default)",
        },
    }
    result["execution"]["runtime_seconds_inside_script"] = time.perf_counter() - start
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({
        "output": str(OUT.relative_to(PROJECT)),
        "runtime": result["execution"]["runtime_seconds_inside_script"],
        "selected": [{"g": x["g"], "event": x["event_or_endpoint"], "z_end": x["z_end"],
                      "q_end": x["q_end"], "x_end": x["x_end"], "identity_rel": x["q_identity_scaled_residual"]}
                     for x in selected],
        "threshold": threshold_compare,
        "code_sha256": result["execution"]["code_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent fixed-metric implementation of the IVS paper's Eqs. 19-20.

This is a prescribed-metric two-component diagnostic, not a Boltzmann solver.
It independently rebuilds the Gate 1 massless-neutrino radiation background;
it does not import Gate 1 or production-profile code.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import cumulative_simpson, solve_ivp
from scipy.interpolate import CubicSpline


ROOT = Path(__file__).resolve().parents[2]
PDF_PATH = ROOT / "context/papers/interacting_de_desi_dr2_2026.pdf"
TEXT_PATH = ROOT / "context/papers/interacting_de_desi_dr2_2026.txt"
SOLVER_PATH = ROOT / "work/perturbation_audit/solver_path.md"
OUT_PATH = Path(__file__).resolve().with_name("perturbation_gate2.json")

# Same pinned homogeneous benchmark as Gate 1, recomputed here from constants.
H0_KM_S_MPC = 67.4
MPC_M = 3.0856775814913673e22
G_NEWTON_SI = 6.67430e-11
C_M_S = 299792458.0
A_RAD_J_M3_K4 = 7.5657e-16
T_CMB_K = 2.7255
N_EFF = 3.044
OMEGA_M0 = 0.315
OMEGA_B0 = 0.049
OMEGA_C0 = OMEGA_M0 - OMEGA_B0
Z_MAX = 1.0e10
N_MIN = -math.log1p(Z_MAX)
G_VALUES = (0.0, 0.46666, -0.46666)
N_EVAL = np.linspace(N_MIN, 0.0, 2049)

# Fixed synchronous-gauge trace source. h(n)=2*S0*(n-N_MIN), so
# h'=dh/dtau=2*S0*Hcal and h'/(2*Hcal)=S0 exactly.
METRIC_SOURCE = 1.0e-3
DELTA_INITIAL = 0.1
THETA_TODAY_OVER_H0 = 1.0e-6

ODE_CONFIG = {"rtol": 1.0e-12, "atol": 1.0e-14, "max_step": 0.02}
QUADRATURE_GRID_SIZES = (8193, 16385)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def pinned_present_densities() -> dict[str, float]:
    """Independent SI blackbody calculation for photons and massless neutrinos."""
    h0_s_inv = H0_KM_S_MPC * 1000.0 / MPC_M
    rho_critical = 3.0 * h0_s_inv**2 / (8.0 * math.pi * G_NEWTON_SI)
    rho_gamma = A_RAD_J_M3_K4 * T_CMB_K**4 / C_M_S**2
    omega_gamma = rho_gamma / rho_critical
    nu_to_gamma = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * N_EFF
    omega_nu = omega_gamma * nu_to_gamma
    omega_x = 1.0 - OMEGA_B0 - OMEGA_C0 - omega_gamma - omega_nu
    result = {
        "omega_b0": OMEGA_B0,
        "omega_c0": OMEGA_C0,
        "omega_gamma0": omega_gamma,
        "omega_nu0": omega_nu,
        "omega_x0": omega_x,
    }
    if not all(math.isfinite(x) and x > 0.0 for x in result.values()):
        raise RuntimeError(f"Pinned density is invalid: {result}")
    if not math.isclose(sum(result.values()), 1.0, rel_tol=0.0, abs_tol=3e-16):
        raise RuntimeError(f"Flat closure failed: {sum(result.values())}")
    return result


PRESENT = pinned_present_densities()


def background_terms(n: float | np.ndarray, rho_c: float | np.ndarray,
                     rho_x: float | np.ndarray) -> tuple[dict[str, np.ndarray], np.ndarray]:
    n_array = np.asarray(n)
    parts = {
        "rho_b": PRESENT["omega_b0"] * np.exp(-3.0 * n_array),
        "rho_gamma": PRESENT["omega_gamma0"] * np.exp(-4.0 * n_array),
        "rho_nu": PRESENT["omega_nu0"] * np.exp(-4.0 * n_array),
        "rho_c": np.asarray(rho_c),
        "rho_x": np.asarray(rho_x),
    }
    e2 = sum(parts.values())
    return parts, e2


def integrate_background(g: float):
    """Recompute rho_c,rho_x independently with n=ln(a) as the variable."""
    def rhs(n: float, state: np.ndarray) -> np.ndarray:
        rc, rx = state
        _, e2 = background_terms(n, rc, rx)
        if not math.isfinite(float(e2)) or e2 <= 0.0:
            raise FloatingPointError(f"Invalid E^2 at n={n}: {e2}")
        transfer = g * rx / math.sqrt(float(e2))
        return np.array([-3.0 * rc - transfer, transfer])

    solution = solve_ivp(
        rhs,
        (0.0, N_MIN),
        [PRESENT["omega_c0"], PRESENT["omega_x0"]],
        method="DOP853",
        t_eval=np.linspace(0.0, N_MIN, N_EVAL.size),
        dense_output=True,
        rtol=ODE_CONFIG["rtol"],
        atol=ODE_CONFIG["atol"],
        max_step=ODE_CONFIG["max_step"],
    )
    if not solution.success:
        raise RuntimeError(f"Background solve failed for g={g}: {solution.message}")
    rc_desc, rx_desc = solution.y
    _, e2_desc = background_terms(solution.t, rc_desc, rx_desc)
    if np.any(rc_desc <= 0.0) or np.any(rx_desc <= 0.0) or np.any(e2_desc <= 0.0):
        raise RuntimeError(f"Background positivity failure for g={g}")
    return solution


def eval_background(solution, n: float | np.ndarray):
    values = solution.sol(n)
    rc, rx = values[0], values[1]
    parts, e2 = background_terms(n, rc, rx)
    if np.any(rc <= 0.0) or np.any(rx <= 0.0) or np.any(e2 <= 0.0):
        raise FloatingPointError("Background density became non-positive")
    return rc, rx, np.sqrt(e2), parts, e2


def metric_source_over_Hcal(n: float | np.ndarray) -> np.ndarray:
    """Analytic h'/(2 Hcal) for h(n)=2*S*(n-N_MIN)."""
    return np.zeros_like(np.asarray(n), dtype=float) + METRIC_SOURCE


def theta_initial(gauge_test: str) -> float:
    if gauge_test == "theta_zero":
        return 0.0
    if gauge_test == "theta_nonzero_diagnostic":
        # Exact solution T(n)=T_today*exp(-n), with T=theta/H0.
        return THETA_TODAY_OVER_H0 * math.exp(-N_MIN)
    raise ValueError(gauge_test)


def perturbation_rhs(background, g: float):
    def rhs(n: float, y: np.ndarray) -> np.ndarray:
        delta, theta_over_h0 = y
        rc, rx, e, _, _ = eval_background(background, n)
        a = math.exp(n)
        hcal_over_h0 = a * e
        metric_ratio = float(metric_source_over_Hcal(n))
        # Paper Eq. (19), converted from conformal time to n=ln(a):
        # delta_n=-theta/Hcal-h'/(2Hcal)+(Gamma/H)*(rho_x/rho_c)*delta.
        theta_over_hcal = theta_over_h0 / hcal_over_h0
        gamma_over_h = g / e
        delta_n = (
            -theta_over_hcal
            - metric_ratio
            + gamma_over_h * (float(rx) / float(rc)) * delta
        )
        # Paper Eq. (20), in n: d(theta/H0)/dn=-(theta/H0).
        return np.array([delta_n, -theta_over_h0])
    return rhs


def solve_perturbations(background, g: float, gauge_test: str):
    theta_i = theta_initial(gauge_test)
    solution = solve_ivp(
        perturbation_rhs(background, g),
        (N_MIN, 0.0),
        [DELTA_INITIAL, theta_i],
        method="DOP853",
        t_eval=N_EVAL,
        dense_output=True,
        rtol=ODE_CONFIG["rtol"],
        atol=ODE_CONFIG["atol"],
        max_step=ODE_CONFIG["max_step"],
    )
    if not solution.success:
        raise RuntimeError(f"Perturbation solve failed for g={g}, {gauge_test}: {solution.message}")
    delta, theta_over_h0 = solution.y
    return solution, delta, theta_over_h0


def integrating_factor_quadrature(background, g: float, gauge_test: str,
                                  n_eval: np.ndarray, n_nodes: int) -> np.ndarray:
    """Composite-Simpson integrating-factor solution independent of delta ODE."""
    n_quad = np.linspace(N_MIN, 0.0, n_nodes)
    rc, rx, e, _, _ = eval_background(background, n_quad)
    a = np.exp(n_quad)
    t_i = theta_initial(gauge_test)
    theta_over_h0 = t_i * np.exp(-(n_quad - N_MIN))
    A = g / e * (rx / rc)
    forcing = -theta_over_h0 / (a * e) - metric_source_over_Hcal(n_quad)
    integrated_A = cumulative_simpson(A, x=n_quad, initial=0.0)
    integrating_factor = np.exp(-integrated_A)
    integrated_forcing = cumulative_simpson(
        integrating_factor * forcing, x=n_quad, initial=0.0
    )
    A_integral_spline = CubicSpline(n_quad, integrated_A)
    forcing_integral_spline = CubicSpline(n_quad, integrated_forcing)
    return np.exp(A_integral_spline(n_eval)) * (
        DELTA_INITIAL + forcing_integral_spline(n_eval)
    )


def solve_absolute_continuity(background, g: float, gauge_test: str):
    """Independently integrate U=rho_c*delta using source-free Eq. 19 form."""
    rc_i, _, _, _, _ = eval_background(background, N_MIN)
    u_i = float(rc_i) * DELTA_INITIAL
    t_i = theta_initial(gauge_test)

    def rhs(n: float, state: np.ndarray) -> np.ndarray:
        u = state[0]
        rc, _, e, _, _ = eval_background(background, n)
        a = math.exp(n)
        theta_over_h0 = t_i * math.exp(-(n - N_MIN))
        force = theta_over_h0 / (a * e) + float(metric_source_over_Hcal(n))
        # Equivalent absolute-CDM continuity: U_n+3U=-rho_c*(theta/Hcal+S).
        return np.array([-3.0 * u - float(rc) * force])

    solution = solve_ivp(
        rhs,
        (N_MIN, 0.0),
        [u_i],
        method="DOP853",
        t_eval=N_EVAL,
        dense_output=True,
        rtol=ODE_CONFIG["rtol"],
        atol=ODE_CONFIG["atol"],
        max_step=ODE_CONFIG["max_step"],
    )
    if not solution.success:
        raise RuntimeError(f"Absolute continuity solve failed: {solution.message}")
    return solution.y[0]


def relative_max(left: np.ndarray, right: np.ndarray) -> float:
    l = np.asarray(left, dtype=float)
    r = np.asarray(right, dtype=float)
    scale = np.maximum(np.maximum(np.abs(l), np.abs(r)), 1.0e-13)
    return float(np.max(np.abs(l - r) / scale))


def audit_variant(background, g: float, gauge_test: str) -> dict[str, object]:
    solution, delta, theta_over_h0 = solve_perturbations(background, g, gauge_test)
    theta_i = theta_initial(gauge_test)
    rc, rx, e, _, e2 = eval_background(background, N_EVAL)
    a = np.exp(N_EVAL)
    theta_over_hcal = theta_over_h0 / (a * e)
    s = metric_source_over_Hcal(N_EVAL)
    source_coefficient = (g / e) * (rx / rc)
    source_term = source_coefficient * delta

    delta_if_low = integrating_factor_quadrature(
        background, g, gauge_test, N_EVAL, QUADRATURE_GRID_SIZES[0]
    )
    delta_if_high = integrating_factor_quadrature(
        background, g, gauge_test, N_EVAL, QUADRATURE_GRID_SIZES[1]
    )

    # Dimensional check: Hcal=a*H0*E and Gamma_s=(Gamma/H0)*H0_s.
    h0_s_inv = H0_KM_S_MPC * 1000.0 / MPC_M
    gamma_s_inv = g * h0_s_inv
    hcal_s_inv = a * h0_s_inv * e
    conformal_source = a * gamma_s_inv / hcal_s_inv
    dimensionless_source = g / e
    unit_error = float(np.max(np.abs(conformal_source - dimensionless_source)))

    # From the contrast equation and rho_c'=-3Hcal*rho_c-a*Gamma*rho_x:
    # U'=d(rho_c delta_c)/d tau, and U'+3Hcal*U+
    # rho_c*(theta_c+h'/2)=0. Evaluate the residual in n units.
    rc_n = -3.0 * rc - (g / e) * rx
    delta_n = -theta_over_hcal - s + source_term
    u = rc * delta
    u_n_from_product = rc_n * delta + rc * delta_n
    absolute_residual_n = u_n_from_product + 3.0 * u + rc * (theta_over_hcal + s)
    residual_scale = np.maximum(
        np.abs(u_n_from_product) + 3.0 * np.abs(u) + np.abs(rc * (theta_over_hcal + s)),
        1.0e-300,
    )
    normalized_residual = np.abs(absolute_residual_n) / residual_scale

    absolute_solution = solve_absolute_continuity(background, g, gauge_test)
    absolute_continuity_delta = relative_max(u, absolute_solution)

    theta_a = a * theta_over_h0
    if gauge_test == "theta_zero":
        theta_check = bool(np.all(theta_over_h0 == 0.0))
        theta_invariant_error = 0.0
    else:
        invariant = theta_i * math.exp(N_MIN)
        theta_invariant_error = float(
            np.max(np.abs(theta_a - invariant)) / max(abs(invariant), 1.0e-300)
        )
        theta_check = theta_invariant_error < 2.0e-9

    if g == 0.0:
        if gauge_test == "theta_zero":
            delta_analytic = DELTA_INITIAL - METRIC_SOURCE * (N_EVAL - N_MIN)
            gamma_zero_error = relative_max(delta, delta_analytic)
        else:
            # The diagnostic has nonzero velocity forcing; it is checked by
            # the independent integrating-factor quadrature below.
            gamma_zero_error = None
        source_sign_check = bool(np.all(source_coefficient == 0.0))
    elif g > 0.0:
        gamma_zero_error = None
        source_sign_check = bool(np.all(source_coefficient > 0.0) and np.all(source_term > 0.0))
    else:
        gamma_zero_error = None
        source_sign_check = bool(np.all(source_coefficient < 0.0) and np.all(source_term < 0.0))

    result = {
        "gamma_over_H0": g,
        "Gamma_km_s_Mpc": g * H0_KM_S_MPC,
        "Gamma_s_inv": gamma_s_inv,
        "theta_test": gauge_test,
        "initial_conditions": {
            "n_initial": N_MIN,
            "delta_c_initial": DELTA_INITIAL,
            "theta_c_over_H0_initial": theta_i,
            "metric_h_of_n": "h(n)=2*S*(n-N_MIN)",
            "metric_hprime_over_2Hcal": METRIC_SOURCE,
        },
        "range_diagnostics": {
            "min_rho_c": float(np.min(rc)),
            "min_rho_x": float(np.min(rx)),
            "min_E2": float(np.min(e2)),
            "delta_c_min": float(np.min(delta)),
            "delta_c_max": float(np.max(delta)),
            "source_coefficient_min": float(np.min(source_coefficient)),
            "source_coefficient_max": float(np.max(source_coefficient)),
            "source_term_min": float(np.min(source_term)),
            "source_term_max": float(np.max(source_term)),
            "theta_a_invariant_max_relative_error": theta_invariant_error,
            "absolute_cdm_continuity_max_normalized_residual": float(
                np.max(normalized_residual)
            ),
            "absolute_continuity_independent_ode_vs_rho_c_delta_max_relative_error": absolute_continuity_delta,
        },
        "quadrature_comparison": {
            "simpson_grid_nodes_low_high": list(QUADRATURE_GRID_SIZES),
            "ode_delta_vs_high_grid_max_relative_error": relative_max(delta, delta_if_high),
            "high_vs_low_quadrature_delta_max_relative_difference": relative_max(
                delta_if_high, delta_if_low
            ),
        },
        "unit_conversion": {
            "H0_s_inv": h0_s_inv,
            "Hcal_equals_a_H0_E": True,
            "max_abs_difference_aGamma_over_Hcal_vs_g_over_E": unit_error,
        },
        "gamma_zero_analytic_delta_max_relative_error": gamma_zero_error,
        "checks": {
            "ode_success": bool(solution.success),
            "all_pinned_background_densities_positive": bool(
                np.all(rc > 0.0) and np.all(rx > 0.0) and np.all(e2 > 0.0)
            ),
            "source_sign_matches_Gamma": source_sign_check,
            "theta_zero_invariant_or_theta_a_constant": theta_check,
            "absolute_cdm_continuity_residual_small": bool(np.max(normalized_residual) < 2e-12),
            "independent_absolute_continuity_ode_agrees": absolute_continuity_delta < 2e-8,
            "integrating_factor_quadrature_agrees": relative_max(delta, delta_if_high) < 2e-7,
            "quadrature_grid_converges": relative_max(delta_if_high, delta_if_low) < 2e-7,
            "aGamma_over_Hcal_equals_Gamma_over_H": unit_error < 2e-15,
            "delta_stays_positive_for_declared_test": bool(np.all(delta > 0.0)),
        },
    }
    return result


def main() -> None:
    start = time.perf_counter()
    all_cases = []
    for g in G_VALUES:
        background = integrate_background(float(g))
        # Physical CDM-comoving solution, plus a distinct nonzero-theta
        # diagnostic of Eq. (20); the latter is explicitly not the chosen gauge.
        all_cases.append(audit_variant(background, float(g), "theta_zero"))
        all_cases.append(audit_variant(background, float(g), "theta_nonzero_diagnostic"))

    all_checks = all(
        all(bool(value) for value in case["checks"].values())
        for case in all_cases
    )
    zero_theta_zero = next(
        case for case in all_cases
        if case["gamma_over_H0"] == 0.0 and case["theta_test"] == "theta_zero"
    )
    all_checks = all_checks and (
        zero_theta_zero["gamma_zero_analytic_delta_max_relative_error"] < 2e-9
    )

    result = {
        "gate": "2: fixed-metric two-component perturbation audit",
        "status": "PASS" if all_checks else "FAIL",
        "primary_source": {
            "pdf": str(PDF_PATH.relative_to(ROOT)),
            "pdf_sha256": sha256(PDF_PATH),
            "text_copy": str(TEXT_PATH.relative_to(ROOT)),
            "text_copy_sha256": sha256(TEXT_PATH),
            "equations_checked": "Eqs. (19)-(20), PDF p. 4",
            "equation_contract": {
                "delta_c_prime_conformal": "-(theta_c+hprime/2)+a*Gamma*(rho_x/rho_c)*delta_c",
                "theta_c_prime_conformal": "-Hcal*theta_c",
                "paper_prime_variable": "conformal time tau",
            },
        },
        "solver_path": {
            "path": str(SOLVER_PATH.relative_to(ROOT)),
            "sha256": sha256(SOLVER_PATH),
            "declared_gate2_scope": "fixed-metric two-component perturbation test; no CAMB sampler or likelihood",
        },
        "implementation": {
            "script": str(Path(__file__).resolve().relative_to(ROOT)),
            "script_sha256": sha256(Path(__file__).resolve()),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cpu_count_visible": os.cpu_count(),
            "cpu_affinity_ids": sorted(os.sched_getaffinity(0)),
            "cpu_affinity_count": len(os.sched_getaffinity(0)),
            "workers": 1,
            "gpu_used": False,
            "elapsed_seconds_inside_script": time.perf_counter() - start,
            "bounded_command": "taskset -c 0 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --state work/compute18/no_campaign_deadline.json --seconds 600 -- .venv/bin/python work/compute18/perturbation_gate2.py",
            "state_path_was_missing": True,
            "explicit_runner_seconds": 600,
        },
        "pinned_background": {
            "independent_reimplementation": True,
            "H0_km_s_Mpc": H0_KM_S_MPC,
            "H0_s_inv": H0_KM_S_MPC * 1000.0 / MPC_M,
            "Omega_m0_baryons_plus_CDM": OMEGA_M0,
            "Omega_b0": OMEGA_B0,
            "Omega_c0": OMEGA_C0,
            "Omega_gamma0": PRESENT["omega_gamma0"],
            "Omega_nu0": PRESENT["omega_nu0"],
            "Omega_x0": PRESENT["omega_x0"],
            "T_CMB_K": T_CMB_K,
            "N_eff": N_EFF,
            "neutrino_approximation": "separately conserved massless effective neutrinos; the source paper's sum m_nu=0.06 eV is not reproduced",
            "flat_closure": "Omega_x0=1-Omega_b0-Omega_c0-Omega_gamma0-Omega_nu0",
            "densities_normalized_to": "today's critical density",
            "Q": "Gamma*rho_x",
            "w_x": -1,
            "sign_convention": "Gamma>0 transfers CDM energy to vacuum",
            "gamma_over_H0_values": list(G_VALUES),
            "redshift_range": [0.0, Z_MAX],
            "background_n_variable": "n=ln(a)",
            "background_equations": {
                "rho_b": "rho_b,n=-3*rho_b",
                "rho_gamma": "rho_gamma,n=-4*rho_gamma",
                "rho_nu": "rho_nu,n=-4*rho_nu",
                "rho_c": "rho_c,n=-3*rho_c-gamma*rho_x/E",
                "rho_x": "rho_x,n=gamma*rho_x/E",
                "E2": "sum of baryon, photon, neutrino, CDM, vacuum densities",
            },
        },
        "fixed_metric_and_frame": {
            "gauge": "synchronous; physical solution uses CDM-comoving theta_c=0",
            "closure": "geodesic CDM/no momentum transfer in CDM rest frame as in solver_path.md",
            "metric_source_choice": "S=hprime/(2*Hcal)=constant=1e-3; h(n)=2*S*(n-N_MIN), with hprime=dh/dtau",
            "theta_diagnostic": "A separate nonzero-theta run checks Eq. 20 and theta*a constancy; this is a diagnostic only and not the CDM-comoving physical run.",
            "initial_delta_c": DELTA_INITIAL,
            "nonzero_theta_today_over_H0": THETA_TODAY_OVER_H0,
            "nonzero_theta_initial_over_H0": theta_initial("theta_nonzero_diagnostic"),
            "conformal_conversion": "a*Gamma/Hcal=(Gamma/H0)/E=gamma/E because Hcal=a*H0*E",
            "delta_equation_in_ln_a": "delta_c,n=-theta_c/Hcal-S+(gamma/E)*(rho_x/rho_c)*delta_c",
            "theta_equation_in_ln_a": "(theta_c/H0),n=-(theta_c/H0)",
            "absolute_cdm_perturbation_identity": "U=rho_c*delta_c obeys U,n+3U+rho_c*(theta_c/Hcal+S)=0; interaction terms cancel",
            "independent_delta_method": "integrating factor with I(n)=integral A dn and delta=exp(I)*(delta_i+integral exp(-I)F dn), evaluated with composite Simpson quadrature on 8193 and 16385 nodes",
        },
        "integrator": {
            "method": "SciPy solve_ivp DOP853",
            "ode": ODE_CONFIG,
            "N_eval_points": int(N_EVAL.size),
            "quadrature_grid_sizes": list(QUADRATURE_GRID_SIZES),
        },
        "cases": all_cases,
        "scope_limit": "This fixed-metric two-component check omits Einstein metric feedback, photon/baryon/neutrino perturbation hierarchies, initial-condition generation, CMB spectra, and any stability classification or full-model reproduction.",
    }
    OUT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"status": result["status"], "output": str(OUT_PATH),
                      "case_count": len(all_cases)}, sort_keys=True))
    if not all_checks:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

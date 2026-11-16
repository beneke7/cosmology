#!/usr/bin/env python3
"""Independent audit of the fixed-metric interacting-vacuum perturbations.

Derivation follows the displayed equations on primary-paper PDF p. 4. This
script does not import the Gate 2 implementation. It uses an RK45 calculation
and only reads the Gate 2 JSON as a saved result for comparison.
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
from scipy.integrate import solve_ivp


HERE = Path(__file__).resolve().parent
PROJECT = Path(__file__).resolve().parents[2]
PDF = PROJECT / "context/papers/interacting_de_desi_dr2_2026.pdf"
GATE_CODE = PROJECT / "work/compute18/perturbation_gate2.py"
GATE_RESULT = PROJECT / "work/compute18/perturbation_gate2.json"
OUT = HERE / "independent_gate2_audit.json"

# Recompute the same declared massless-neutrino benchmark from constants.
H0_KM_S_MPC = 67.4
MPC_M = 3.0856775814913673e22
G_SI = 6.67430e-11
C_M_S = 299792458.0
A_RAD = 7.5657e-16
TCMB_K = 2.7255
NEFF = 3.044
OMEGA_M = 0.315
OMEGA_B = 0.049
OMEGA_C = OMEGA_M - OMEGA_B
ZMAX = 1e10
N0 = -math.log1p(ZMAX)
GAMMAS = (0.0, 0.46666, -0.46666)
N_GRID = np.linspace(N0, 0.0, 2049)
DELTA_I = 0.1
S_METRIC = 1e-3
T_TODAY = 1e-6
CFG = {"rtol": 1e-13, "atol": 1e-15, "max_step": 0.01}

H0_S = H0_KM_S_MPC * 1000.0 / MPC_M
RHOC0 = 3.0 * H0_S**2 / (8.0 * math.pi * G_SI)
RHO_GAMMA0 = A_RAD * TCMB_K**4 / C_M_S**2
OGAMMA = RHO_GAMMA0 / RHOC0
NU_RATIO = (7.0 / 8.0) * (4.0 / 11.0)**(4.0 / 3.0) * NEFF
ONU = OGAMMA * NU_RATIO
ORAD = OGAMMA + ONU
OX = 1.0 - OMEGA_B - OMEGA_C - ORAD


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def e2_at(n, rc, rx):
    return OMEGA_B * np.exp(-3.0 * np.asarray(n)) + ORAD * np.exp(-4.0 * np.asarray(n)) + rc + rx


def background(g: float):
    """Integrate paper's homogeneous continuity equations backward from today."""
    def f(n, y):
        rc, rx = y
        e2 = float(e2_at(n, rc, rx))
        if e2 <= 0.0 or not math.isfinite(e2):
            raise FloatingPointError(f"bad E2 at n={n}: {e2}")
        transfer = g * rx / math.sqrt(e2)
        return (-3.0 * rc - transfer, transfer)

    sol = solve_ivp(
        f, (0.0, N0), (OMEGA_C, OX), method="RK45",
        t_eval=np.linspace(0.0, N0, 2049), dense_output=True,
        **CFG,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def bg_at(sol, n):
    y = sol.sol(n)
    rc, rx = y[0], y[1]
    e2 = e2_at(n, rc, rx)
    if np.any(np.asarray(rc) <= 0) or np.any(np.asarray(rx) <= 0) or np.any(np.asarray(e2) <= 0):
        raise FloatingPointError("non-positive background component")
    return rc, rx, np.sqrt(e2), e2


def theta_init(kind: str) -> float:
    if kind == "theta_zero":
        return 0.0
    if kind == "theta_nonzero_diagnostic":
        # Choose T(0)=theta(today)/H0=1e-6; Eq. (20) implies T_n=-T.
        return T_TODAY * math.exp(-N0)
    raise ValueError(kind)


def solve_delta(bg, g: float, kind: str):
    t_i = theta_init(kind)
    def rhs(n, y):
        delta = y[0]
        # Integrate Eq. (20) analytically: T_n=-T implies T=T_i exp[-(n-n_i)].
        theta_h0 = t_i * math.exp(-(n - N0))
        rc, rx, e, _ = bg_at(bg, n)
        a = math.exp(n)
        hcal_over_h0 = a * float(e)
        # Eq. (19): delta' = -(theta+h'/2) + a Gamma (rho_x/rho_c) delta.
        # Since h'=2*S*Hcal and d/dn=(1/Hcal)d/dtau, its metric term is -S.
        delta_n = (
            -theta_h0 / hcal_over_h0 - S_METRIC
            + (g / float(e)) * (float(rx) / float(rc)) * delta
        )
        # Eq. (20) after division by Hcal: d(theta/H0)/dn=-(theta/H0).
        return (delta_n,)

    sol = solve_ivp(
        rhs, (N0, 0.0), (DELTA_I,), method="RK45",
        t_eval=N_GRID, dense_output=True, **CFG,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    theta_grid = t_i * np.exp(-(N_GRID - N0))
    return sol, theta_grid


def solve_absolute_u(bg, g: float, kind: str):
    """Integrate a^3 U using the separately derived absolute-CDM law."""
    t_i = theta_init(kind)
    rc_i = float(bg_at(bg, N0)[0])
    def rhs(n, y):
        rc, _, e, _ = bg_at(bg, n)
        a = math.exp(n)
        theta_h0 = t_i * math.exp(-(n - N0))
        theta_hcal = theta_h0 / (a * float(e))
        # For U_n+3U+rho_c*F=0, V=a^3*U obeys V_n=-a^3*rho_c*F.
        return (-a**3 * float(rc) * (theta_hcal + S_METRIC),)
    sol = solve_ivp(
        rhs, (N0, 0.0), (math.exp(3.0 * N0) * rc_i * DELTA_I,), method="RK45",
        t_eval=N_GRID, **CFG,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return np.exp(-3.0 * N_GRID) * sol.y[0]


def max_rel(a, b):
    aa, bb = np.asarray(a), np.asarray(b)
    return float(np.max(np.abs(aa - bb) / np.maximum(np.maximum(np.abs(aa), np.abs(bb)), 1e-300)))


def audit_case(bg, g: float, kind: str, reference: dict) -> dict:
    pert, theta_h0 = solve_delta(bg, g, kind)
    delta = pert.y[0]
    rc, rx, e, e2 = bg_at(bg, N_GRID)
    a = np.exp(N_GRID)
    hcal_h0 = a * e
    theta_hcal = theta_h0 / hcal_h0
    coefficient = (g / e) * (rx / rc)
    source_term = coefficient * delta
    # Build dU/dn independently from the product rule and check Eq. (19).
    rc_n = -3.0 * rc - (g / e) * rx
    delta_n = -theta_hcal - S_METRIC + source_term
    u = rc * delta
    u_n_product = rc_n * delta + rc * delta_n
    u_resid = u_n_product + 3.0 * u + rc * (theta_hcal + S_METRIC)
    denom = np.maximum(
        np.abs(u_n_product) + 3.0 * np.abs(u) + np.abs(rc * (theta_hcal + S_METRIC)),
        1e-300,
    )
    norm_resid = np.abs(u_resid) / denom
    u_separate = solve_absolute_u(bg, g, kind)
    u_rel = max_rel(u, u_separate)

    # Direct dimensional evaluation of a*Gamma/Hcal vs Gamma/H.
    gamma_km_s_mpc = g * H0_KM_S_MPC
    hcal_km_s_mpc = a * H0_KM_S_MPC * e
    lhs = a * gamma_km_s_mpc / hcal_km_s_mpc
    rhs = g / e
    unit_abs = float(np.max(np.abs(lhs - rhs)))

    theta_inv = a * theta_h0
    if kind == "theta_zero":
        theta_metric = float(np.max(np.abs(theta_h0)))
        theta_inv_error = 0.0
    else:
        theta_metric = float(np.max(np.abs(theta_h0)))
        theta_inv_error = max_rel(theta_inv, np.full_like(a, a[0] * theta_h0[0]))

    if g == 0.0 and kind == "theta_zero":
        analytic = DELTA_I - S_METRIC * (N_GRID - N0)
        analytic_error = max_rel(delta, analytic)
    else:
        analytic_error = None

    result = {
        "gamma_over_H0": g,
        "theta_test": kind,
        "initial_theta_over_H0": theta_init(kind),
        "range": {
            "delta_c_min": float(np.min(delta)),
            "delta_c_max": float(np.max(delta)),
            "rho_c_min": float(np.min(rc)),
            "rho_x_min": float(np.min(rx)),
            "E2_min": float(np.min(e2)),
            "theta_over_H0_max_abs": theta_metric,
            "source_coefficient_min": float(np.min(coefficient)),
            "source_coefficient_max": float(np.max(coefficient)),
            "source_term_min": float(np.min(source_term)),
            "source_term_max": float(np.max(source_term)),
            "theta_a_invariant_max_relative_error": theta_inv_error,
            "absolute_cdm_continuity_max_normalized_residual": float(np.max(norm_resid)),
            "absolute_continuity_independent_ode_vs_rho_c_delta_max_relative_error": u_rel,
        },
        "unit_conversion_max_abs_error": unit_abs,
        "gamma_zero_analytic_delta_max_relative_error": analytic_error,
        "compare_saved_gate2": {},
    }
    # Metrics saved in Gate 2 JSON are the comparison target; no Gate code is imported.
    for name in ("delta_c_min", "delta_c_max", "min_rho_c", "min_rho_x", "min_E2",
                 "theta_a_invariant_max_relative_error",
                 "absolute_cdm_continuity_max_normalized_residual",
                 "absolute_continuity_independent_ode_vs_rho_c_delta_max_relative_error"):
        ours_name = {"min_rho_c": "rho_c_min", "min_rho_x": "rho_x_min", "min_E2": "E2_min"}.get(name, name)
        ours = result["range"][ours_name]
        theirs = reference["range_diagnostics"][name]
        result["compare_saved_gate2"][name] = {
            "independent": ours,
            "saved": theirs,
            "absolute_difference": abs(ours - theirs),
            "relative_difference": abs(ours - theirs) / max(abs(ours), abs(theirs), 1e-300),
        }
    result["compare_saved_gate2"]["max_abs_unit_conversion_error_difference"] = {
        "independent": unit_abs,
        "saved": reference["unit_conversion"]["max_abs_difference_aGamma_over_Hcal_vs_g_over_E"],
        "absolute_difference": abs(unit_abs - reference["unit_conversion"]["max_abs_difference_aGamma_over_Hcal_vs_g_over_E"]),
    }
    return result


def main():
    started = time.perf_counter()
    saved = json.loads(GATE_RESULT.read_text())
    cases = []
    for g in GAMMAS:
        bg = background(g)
        bg_rc, bg_rx, bg_e, bg_e2 = bg_at(bg, N_GRID)
        if not (np.all(bg_rc > 0) and np.all(bg_rx > 0) and np.all(bg_e2 > 0)):
            raise RuntimeError(f"background positivity failure at gamma={g}")
        for kind in ("theta_zero", "theta_nonzero_diagnostic"):
            ref = next(c for c in saved["cases"]
                       if c["gamma_over_H0"] == g and c["theta_test"] == kind)
            cases.append(audit_case(bg, g, kind, ref))

    result = {
        "assessment": "independent derivation and numerical cross-check; not a Boltzmann or likelihood reproduction",
        "source": {
            "paper": str(PDF.relative_to(PROJECT)),
            "pdf_sha256": sha(PDF),
            "page": 4,
            "equations": "displayed IVS Eqs. (19)-(20); visually checked",
        },
        "derivation": {
            "conformal_equations": {
                "delta_c_prime": "-(theta_c+h_prime/2)+a*Gamma*(rho_x/rho_c)*delta_c",
                "theta_c_prime": "-Hcal*theta_c",
                "prime": "d/dtau; Hcal=a'/a=aH",
            },
            "conversion": "dn/dtau=Hcal; d/dn=(1/Hcal)d/dtau; a*Gamma/Hcal=Gamma/H=(Gamma/H0)/E=g/E",
            "ln_a_equations": {
                "metric_source": "S=h_prime/(2*Hcal)=h_n/2",
                "delta_c_n": "-theta_c/Hcal-S+(g/E)*(rho_x/rho_c)*delta_c",
                "theta_over_H0_n": "-(theta_c/H0)",
            },
            "theta_zero": "Eq. (20) is homogeneous; theta_c=0 is invariant. In synchronous gauge with geodesic CDM/no momentum transfer in the CDM rest frame this is the CDM-comoving choice.",
            "absolute_cdm_continuity": "rho_c,n=-3*rho_c-(g/E)*rho_x; U=rho_c*delta_c therefore U_n+3U+rho_c*(theta_c/Hcal+S)=0; the Gamma terms cancel.",
            "frame_scope": "paper states synchronous gauge, zero anisotropic stress, and no momentum transfer in the DM rest frame for its perturbation setup; Eq. (20) permits the theta_c=0 solution. This is a fixed-metric diagnostic, not a full gauge-invariant Einstein-fluid evolution.",
        },
        "benchmark": {
            "H0_km_s_Mpc": H0_KM_S_MPC,
            "Omega_m0_baryons_plus_CDM": OMEGA_M,
            "Omega_b0": OMEGA_B,
            "Omega_c0": OMEGA_C,
            "Omega_gamma0": OGAMMA,
            "Omega_nu0": ONU,
            "Omega_x0": OX,
            "T_CMB_K": TCMB_K,
            "Neff": NEFF,
            "neutrino_assumption": "massless effective neutrinos, as in saved Gate2 benchmark; paper footnote fixes sum m_nu=0.06 eV and is not reproduced",
            "gamma_values": list(GAMMAS),
            "redshift_range": [0.0, ZMAX],
            "fixed_metric": "S=1e-3; h(n)=2*S*(n-N0), so h_prime=2*S*Hcal",
            "delta_initial": DELTA_I,
            "nonzero_theta_today_over_H0": T_TODAY,
        },
        "integrator": {
            "method": "SciPy solve_ivp RK45 (different from saved Gate2 DOP853)",
            **CFG,
            "output_grid_points": int(N_GRID.size),
            "cpu_affinity_ids": sorted(os.sched_getaffinity(0)),
        },
        "comparison_reference": {
            "json_path": str(GATE_RESULT.relative_to(PROJECT)),
            "json_sha256": sha(GATE_RESULT),
            "code_path": str(GATE_CODE.relative_to(PROJECT)),
            "code_sha256": sha(GATE_CODE),
            "comparison_basis": "six saved case summaries from JSON; implementation inspected, not imported",
        },
        "execution": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "gpu_used": False,
            "elapsed_seconds": None,
            "code_path": str(Path(__file__).relative_to(PROJECT)),
            "code_sha256": sha(Path(__file__).resolve()),
            "bounded_command": "taskset -c 0 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --state work/compute20/no_campaign_deadline.json --seconds 600 -- .venv/bin/python work/compute20/independent_gate2_audit.py",
            "model_effort": "GPT-6; exact model variant and reasoning-effort tier are not exposed to this worker runtime (inherited default).",
        },
        "cases": cases,
    }
    result["execution"]["elapsed_seconds"] = time.perf_counter() - started
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({
        "output": str(OUT.relative_to(PROJECT)),
        "elapsed_seconds": result["execution"]["elapsed_seconds"],
        "case_count": len(cases),
        "cases": [{
            "g": c["gamma_over_H0"], "theta": c["theta_test"],
            "delta_min": c["range"]["delta_c_min"],
            "delta_max": c["range"]["delta_c_max"],
            "unit_error": c["unit_conversion_max_abs_error"],
            "U_identity": c["range"]["absolute_cdm_continuity_max_normalized_residual"],
            "saved_delta_min_absdiff": c["compare_saved_gate2"]["delta_c_min"]["absolute_difference"],
            "saved_delta_max_absdiff": c["compare_saved_gate2"]["delta_c_max"]["absolute_difference"],
        } for c in cases],
        "code_sha256": result["execution"]["code_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()

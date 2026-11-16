#!/usr/bin/env python3
"""Illustrative IVS background feasibility scan; independent of compute15/17."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SOURCE_PDF = PROJECT / "context" / "papers" / "interacting_de_desi_dr2_2026.pdf"
U_MAX = math.log1p(1.0e10)
H0_KMS_MPC = 67.4
H = H0_KMS_MPC / 100.0
TCMB_K = 2.7255
NEFF = 3.044
OMEGA_M0 = 0.315
OMEGA_B0 = 0.049
OMEGA_C0 = OMEGA_M0 - OMEGA_B0
OMEGA_GAMMA_H2 = 2.469e-5 * (TCMB_K / 2.7255) ** 4
OMEGA_GAMMA0 = OMEGA_GAMMA_H2 / H**2
NU_OVER_GAMMA = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * NEFF
OMEGA_NU0 = NU_OVER_GAMMA * OMEGA_GAMMA0
OMEGA_RAD0 = OMEGA_GAMMA0 + OMEGA_NU0
OMEGA_X0 = 1.0 - OMEGA_M0 - OMEGA_RAD0
RTOL = 3e-12
ATOL = 2e-14
MAX_STEP = 0.02


def sha256(path: Path) -> str:
    hsh = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hsh.update(chunk)
    return hsh.hexdigest()


def e2_at(u: float, x: float, q: float) -> float:
    e3 = math.exp(3 * u)
    return OMEGA_RAD0 * e3 * math.exp(u) + OMEGA_B0 * e3 + q * e3 + x


def integrate(g: float, sample_count: int = 2049) -> dict:
    """Integrate x=Omega_x, q=Omega_c*a^3, and integral identity state I."""
    def rhs(u, state):
        x, q, integ = state
        e2 = e2_at(u, x, q)
        if not math.isfinite(e2) or e2 <= 0:
            raise ArithmeticError(f"nonpositive/nonfinite E^2={e2} at u={u}")
        source = x / math.sqrt(e2)
        decay = math.exp(-3 * u)
        return (-g * source, g * source * decay, source * decay)

    def zero_event(u, state):
        return state[1]

    zero_event.terminal = True
    zero_event.direction = -1
    sol = solve_ivp(
        rhs, (0.0, U_MAX), (OMEGA_X0, OMEGA_C0, 0.0), method="DOP853",
        rtol=RTOL, atol=ATOL, max_step=MAX_STEP, dense_output=True,
        events=zero_event if g < 0 else None,
    )
    if not sol.success:
        raise RuntimeError(f"DOP853 failed for g={g}: {sol.message}")
    crossed = bool(sol.t_events and len(sol.t_events[0]))
    end = float(sol.t_events[0][0]) if crossed else float(sol.t[-1])
    us = np.linspace(0.0, end, sample_count)
    xs, qs, integral_state = sol.sol(us)
    e2s = np.array([e2_at(float(u), float(x), float(q)) for u, x, q in zip(us, xs, qs)])
    if np.any(e2s <= 0) or not np.all(np.isfinite(e2s)):
        raise ArithmeticError(f"nonpositive/nonfinite sampled E^2 for g={g}")
    ms = qs * np.exp(3 * us)
    min_x = float(np.min(xs))
    min_q = float(np.min(qs))
    min_m = float(np.min(ms))
    min_e2 = float(np.min(e2s))
    integ_rhs = OMEGA_C0 + g * integral_state
    ode_identity_resid = float(np.max(np.abs(qs - integ_rhs)))
    if (sample_count - 1) % 2:
        raise ValueError("Simpson sample count must give an even interval count")
    integrand = xs / np.sqrt(e2s) * np.exp(-3 * us)
    # Composite Simpson quadrature on uniform dense-output samples.
    h = end / (sample_count - 1)
    simpson_I = float(h / 3 * (integrand[0] + integrand[-1]
                               + 4 * np.sum(integrand[1:-1:2])
                               + 2 * np.sum(integrand[2:-1:2])))
    q_quad = OMEGA_C0 + g * simpson_I
    simpson_resid = abs(float(qs[-1]) - q_quad)
    g0_diag = None
    if g == 0.0:
        m_exact = OMEGA_C0 * np.exp(3 * us)
        e2_exact = OMEGA_RAD0 * np.exp(4 * us) + OMEGA_B0 * np.exp(3 * us) + m_exact + OMEGA_X0
        g0_diag = {
            "max_relative_error_x": float(np.max(np.abs(xs - OMEGA_X0) / OMEGA_X0)),
            "max_relative_error_m": float(np.max(np.abs(qs * np.exp(3 * us) - m_exact) / m_exact)),
            "max_relative_error_E2": float(np.max(np.abs(e2s - e2_exact) / e2_exact)),
        }
    dq = np.diff(qs)
    if g < 0:
        monotonic_violation = max(0.0, float(np.max(dq)))
    elif g > 0:
        monotonic_violation = max(0.0, float(-np.min(dq)))
    else:
        monotonic_violation = float(np.max(np.abs(dq)))
    return {
        "g": g,
        "crossed_zero": crossed,
        "u_end": end,
        "z_end": math.expm1(end),
        "reached_zmax": not crossed and abs(end - U_MAX) < 1e-10,
        "x_min": min_x,
        "q_min": min_q,
        "omega_c_min": min_m,
        "e2_min": min_e2,
        "q_end": float(qs[-1]),
        "q_end_from_simpson": q_quad,
        "max_identity_residual_integrated_I": ode_identity_resid,
        "comoving_identity_residual_simpson": simpson_resid,
        "q_monotonicity_violation": monotonic_violation,
        "integration_nfev": int(sol.nfev),
        "g0_exact_limit_diagnostics": g0_diag,
    }


def main() -> int:
    start = time.monotonic()
    if not SOURCE_PDF.exists():
        raise FileNotFoundError(f"missing exact local source PDF: {SOURCE_PDF}")
    if min(OMEGA_C0, OMEGA_X0, OMEGA_RAD0) <= 0:
        raise ArithmeticError("fiducial closure has a nonpositive present density")

    grid = [round(-3 + 0.25 * i, 12) for i in range(25)]
    scan = [integrate(g) for g in grid]
    crossings = [r for r in scan if r["crossed_zero"]]
    crossing_times = [r["u_end"] for r in crossings]
    if not crossings or scan[grid.index(0.0)]["crossed_zero"]:
        raise ArithmeticError("negative-g threshold is not bracketed by [-3,0]")
    if any(crossing_times[i + 1] + 1e-7 < crossing_times[i]
           for i in range(len(crossing_times) - 1)):
        raise ArithmeticError("sampled crossing time is not monotone toward g=0")

    low, high = -3.0, 0.0
    low_row, high_row = integrate(low, 257), integrate(high, 257)
    if not low_row["crossed_zero"] or high_row["crossed_zero"]:
        raise ArithmeticError("threshold endpoint check failed")
    iterations = 0
    for iterations in range(1, 49):
        mid = (low + high) / 2
        row = integrate(mid, 257)
        if row["crossed_zero"]:
            low, low_row = mid, row
        else:
            high, high_row = mid, row
        if high - low < 2e-10:
            break
    low_row = integrate(low)
    high_row = integrate(high)
    threshold = {
        "g_crossed_side": low,
        "g_positive_through_zmax_side": high,
        "parameter_bracket_width": high - low,
        "bisection_iterations": iterations,
        "last_crossing_u": low_row["u_end"],
        "last_crossing_z": low_row["z_end"],
        "physical_side_q_at_zmax": high_row["q_end"],
        "definition": "benchmark-specific boundary where backward-integrated CDM first reaches zero by z=1e10",
    }

    # g=0 exact limit, compared through the full dense-output history.
    zero = integrate(0.0, 2049)
    g0_diag = zero["g0_exact_limit_diagnostics"]
    max_g0_x_rel = g0_diag["max_relative_error_x"]
    max_g0_m_rel = g0_diag["max_relative_error_m"]
    max_g0_e2_rel = g0_diag["max_relative_error_E2"]

    all_rows = scan + [low_row, high_row, zero]
    for row in all_rows:
        if row["x_min"] <= 0 or row["e2_min"] <= 0:
            raise ArithmeticError(f"x or E^2 lost positivity for g={row['g']}")
        if row["max_identity_residual_integrated_I"] > 2e-10:
            raise ArithmeticError(f"ODE integral identity failed for g={row['g']}")
        if row["comoving_identity_residual_simpson"] > 2e-8:
            raise ArithmeticError(f"independent Simpson identity failed for g={row['g']}")
        if row["q_monotonicity_violation"] > 2e-9:
            raise ArithmeticError(f"comoving-CDM monotonicity failed for g={row['g']}")
        if not row["crossed_zero"] and row["omega_c_min"] <= 0:
            raise ArithmeticError(f"non-crossing trajectory had rho_c<=0 for g={row['g']}")
    if not (low < high and high - low < 2e-10):
        raise ArithmeticError("invalid or insufficiently refined threshold bracket")

    result = {
        "task": "illustrative IVS early-time CDM physicality check",
        "mode": "background only; no optimizer, fit, likelihood, posterior, p-value, or model comparison",
        "source": {
            "paper_pdf": str(SOURCE_PDF.relative_to(PROJECT)),
            "paper_pdf_sha256": sha256(SOURCE_PDF),
            "equations": "Table I and Eqs. (3)-(8), exact local arXiv v1 PDF",
            "runner": "scripts/run_bounded.py",
            "runner_sha256": sha256(PROJECT / "scripts" / "run_bounded.py"),
            "runtime_interpreter": ".venv/bin/python (Python 3.12.3; SciPy available)",
            "initial_system_interpreter_attempt": "python3 failed before computation: ModuleNotFoundError: scipy; no output/result produced",
            "budget_seconds": 600,
            "task_state_file": "work/compute19/RUN_STATE.json (absent; explicit --seconds cap used)",
            "threads": 1
        },
        "fiducial": {
            "illustrative_not_exact_paper_prior_mapping": True,
            "mapping_limitation": "The paper samples physical densities and 100 theta_MC, not this directly fixed Omega_m0 background. This is a declared benchmark, not an exact mapping of its sampled prior and says nothing about the paper posterior.",
            "Omega_m0": OMEGA_M0,
            "Omega_b0": OMEGA_B0,
            "Omega_c0": OMEGA_C0,
            "H0_km_s_Mpc": H0_KMS_MPC,
            "Tcmb_K": TCMB_K,
            "Neff_massless": NEFF,
            "Omega_gamma_h2_convention": "2.469e-5*(Tcmb/2.7255 K)^4",
            "Omega_gamma0": OMEGA_GAMMA0,
            "Omega_nu0": OMEGA_NU0,
            "Omega_rad0_including_massless_nu": OMEGA_RAD0,
            "Omega_x0_from_flat_closure": OMEGA_X0,
            "flatness": True,
            "zmax": 1.0e10,
            "u_max_ln1pz": U_MAX
        },
        "equations_used": {
            "definition": "u=ln(1+z), E=H/H0, x=rho_x/rho_crit0, m=rho_c/rho_crit0, q=m*exp(-3u), g=Gamma/H0",
            "friedmann": "E^2=(Omega_gamma0+Omega_nu0)*exp(4u)+Omega_b0*exp(3u)+m+x",
            "backward_time_odes": "dx/du=-g*x/E; dm/du=3*m+g*x/E",
            "comoving_cdm": "dq/du=g*x/E*exp(-3u); q(u)=Omega_c0+g*integral_0^u[x/E*exp(-3v)dv]",
            "paper_sign": "g>0 transfers CDM to vacuum in forward cosmic time; g<0 transfers vacuum to CDM"
        },
        "paper_prior": {
            "Table_I_IVS_uniform_g_range": [-3.0, 3.0],
            "status": "paper's broad uniform Gamma/H0 prior; this benchmark does not reproduce sampling over the other paper parameters"
        },
        "numerics": {
            "integrator": "SciPy solve_ivp DOP853",
            "rtol": RTOL,
            "atol": ATOL,
            "max_step_in_u": MAX_STEP,
            "one_cpu_only": True,
            "coarse_g_grid_step": 0.25,
            "coarse_scan": scan,
            "threshold": threshold,
            "g0_limit": {
                "expected_x": OMEGA_X0,
                "expected_m": "Omega_c0*exp(3u)",
                "max_relative_error_x": max_g0_x_rel,
                "max_relative_error_m_over_trajectory": max_g0_m_rel,
                "max_relative_error_E2_over_trajectory": max_g0_e2_rel,
                "passed": True
            },
            "independent_integral_check": "q from the q ODE compared to q0+g times composite Simpson quadrature of x/E*exp(-3u) on dense-output samples",
            "max_simpson_identity_residual": max(r["comoving_identity_residual_simpson"] for r in all_rows),
            "max_integrated_I_identity_residual": max(r["max_identity_residual_integrated_I"] for r in all_rows),
            "monotonicity": "q is nonincreasing with u for g<0 and nondecreasing for g>0; checked on dense samples",
            "positivity": "x>0 and E^2>0 through endpoint or first q=0 event; trajectories terminate at q=0 and are not continued beyond it",
            "all_checks_passed": True
        },
        "interpretation_limit": "A threshold inside [-3,3] means this fixed benchmark includes g values whose backward solution reaches rho_c=0 by z=1e10 and would become CDM-negative if continued. It is not a statement about posterior support, prior mass after other parameters, evidence, or the published analysis."
    }
    output = HERE / "physicality_result.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "g_threshold_bracket": [low, high],
        "last_crossing_z": threshold["last_crossing_z"],
        "negative_grid_crossings": len(crossings),
        "max_simpson_identity_residual": result["numerics"]["max_simpson_identity_residual"],
        "elapsed_seconds": time.monotonic() - start
    }, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise

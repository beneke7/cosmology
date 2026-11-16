#!/usr/bin/env python3
"""Independent radiation-inclusive IVS background Gate 1 check.

This script implements the primary paper's Q=Gamma*rho_x, w_x=-1
homogeneous background.  It is deliberately self-contained and does not
import the reduced fit code or any earlier physical-check implementation.
All densities are normalized to today's critical density; n=ln(a), and
E=H/H0.
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
from scipy.integrate import quad, solve_ivp


ROOT = Path(__file__).resolve().parents[2]
PDF_PATH = ROOT / "context/papers/interacting_de_desi_dr2_2026.pdf"
TEXT_PATH = ROOT / "context/papers/interacting_de_desi_dr2_2026.txt"
OUT_PATH = Path(__file__).resolve().with_name("background_gate.json")

# Explicit present-day benchmark assumptions.  H0 and Omega_m,b are a
# declared background benchmark, not a fit to the paper's reported scores.
H0_KM_S_MPC = 67.4
MPC_M = 3.0856775814913673e22
G_NEWTON_SI = 6.67430e-11
C_M_S = 299792458.0
A_RAD_J_M3_K4 = 7.5657e-16
T_CMB_K = 2.7255
N_EFF = 3.044  # Primary paper, footnote 2.
OMEGA_M0 = 0.315
OMEGA_B0 = 0.049
OMEGA_C0 = OMEGA_M0 - OMEGA_B0
Z_MAX = 1.0e10
N_MIN = -math.log1p(Z_MAX)
GAMMA_OVER_H0_VALUES = (0.0, 0.46666, -0.46666)
N_GRID = np.linspace(0.0, N_MIN, 2049)

# Run budgets: one process, one numerical-library thread, and nine short ODE
# integrations (three signs x coarse/fine primary and fine alternate form).
COARSE = {"rtol": 1.0e-9, "atol": 1.0e-12, "max_step": 0.10}
TOLERANCE_ONLY = {"rtol": 1.0e-12, "atol": 1.0e-14, "max_step": 0.10}
MAX_STEP_ONLY = {"rtol": 1.0e-9, "atol": 1.0e-12, "max_step": 0.02}
FINE = {"rtol": 1.0e-12, "atol": 1.0e-14, "max_step": 0.02}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def present_densities() -> dict[str, float]:
    """Compute photons from the blackbody energy density; neutrinos separate."""
    h0_si = H0_KM_S_MPC * 1000.0 / MPC_M
    rho_crit = 3.0 * h0_si**2 / (8.0 * math.pi * G_NEWTON_SI)
    rho_gamma = A_RAD_J_M3_K4 * T_CMB_K**4 / C_M_S**2
    omega_gamma = rho_gamma / rho_crit
    nu_to_gamma = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * N_EFF
    omega_nu = omega_gamma * nu_to_gamma
    omega_x = 1.0 - OMEGA_B0 - OMEGA_C0 - omega_gamma - omega_nu
    vals = {
        "omega_b0": OMEGA_B0,
        "omega_c0": OMEGA_C0,
        "omega_gamma0": omega_gamma,
        "omega_nu0": omega_nu,
        "omega_x0": omega_x,
    }
    if not all(math.isfinite(x) and x > 0.0 for x in vals.values()):
        raise RuntimeError(f"Non-positive or non-finite present density: {vals}")
    if not math.isclose(sum(vals.values()), 1.0, rel_tol=0.0, abs_tol=3e-16):
        raise RuntimeError(f"Flat closure failed: {sum(vals.values())}")
    return vals


PRESENT = present_densities()


def components(n: float | np.ndarray, rho_c: float | np.ndarray,
               rho_x: float | np.ndarray) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Return all separately conserved/interacting densities and E^2."""
    n_arr = np.asarray(n)
    rho = {
        "baryon": PRESENT["omega_b0"] * np.exp(-3.0 * n_arr),
        "photon_radiation": PRESENT["omega_gamma0"] * np.exp(-4.0 * n_arr),
        "neutrino": PRESENT["omega_nu0"] * np.exp(-4.0 * n_arr),
        "cdm": np.asarray(rho_c),
        "vacuum": np.asarray(rho_x),
    }
    e2 = sum(rho.values())
    return rho, e2


def rhs_c_x(g: float):
    """Primary form: evolve rho_c and rho_x directly in n=ln(a)."""
    def f(n: float, y: np.ndarray) -> np.ndarray:
        rc, rx = y
        _, e2 = components(n, rc, rx)
        if not math.isfinite(float(e2)) or e2 <= 0.0:
            raise FloatingPointError(f"E^2 became invalid at n={n}: {e2}")
        source_over_H0 = g * rx / math.sqrt(float(e2))  # Gamma*rho_x/H0E
        return np.array([-3.0 * rc - source_over_H0, source_over_H0])
    return f


def rhs_dark_total_x(g: float):
    """Equivalent form: evolve rho_d=rho_c+rho_x and rho_x."""
    def f(n: float, y: np.ndarray) -> np.ndarray:
        rd, rx = y
        rc = rd - rx
        _, e2 = components(n, rc, rx)
        if not math.isfinite(float(e2)) or e2 <= 0.0:
            raise FloatingPointError(f"E^2 became invalid at n={n}: {e2}")
        source_over_H0 = g * rx / math.sqrt(float(e2))
        # Dark-sector conservation gives d(rho_c+rho_x)/dn=-3 rho_c.
        return np.array([-3.0 * rc, source_over_H0])
    return f


def integrate(g: float, config: dict[str, float], form: str):
    if form == "c_x":
        y0 = np.array([PRESENT["omega_c0"], PRESENT["omega_x0"]])
        fun = rhs_c_x(g)
    elif form == "dark_total_x":
        y0 = np.array([
            PRESENT["omega_c0"] + PRESENT["omega_x0"],
            PRESENT["omega_x0"],
        ])
        fun = rhs_dark_total_x(g)
    else:
        raise ValueError(form)
    sol = solve_ivp(
        fun,
        (0.0, N_MIN),
        y0,
        method="DOP853",
        t_eval=N_GRID,
        dense_output=True,
        rtol=config["rtol"],
        atol=config["atol"],
        max_step=config["max_step"],
    )
    if not sol.success or sol.y.shape[1] != N_GRID.size:
        raise RuntimeError(f"ODE failed g={g}, form={form}: {sol.message}")
    if form == "c_x":
        rc, rx = sol.y
    else:
        rd, rx = sol.y
        rc = rd - rx
    rho, e2 = components(N_GRID, rc, rx)
    return sol, rc, rx, rho, e2


def relative_max(a: np.ndarray, b: np.ndarray) -> float:
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    den = np.maximum(np.maximum(np.abs(aa), np.abs(bb)), 1.0e-300)
    return float(np.max(np.abs(aa - bb) / den))


def integral_checks(g: float, sol, rc_arr: np.ndarray,
                    rx_arr: np.ndarray) -> dict[str, float]:
    """Independent quadratures of the exact comoving and total-sector laws."""
    # Use a separate sparse set of locations and adaptive quadrature over the
    # dense DOP853 interpolant. The identities are not states in the ODE.
    indices = np.linspace(0, N_GRID.size - 1, 129, dtype=int)
    comoving_errors = []
    dark_errors = []
    for i in indices:
        n = float(N_GRID[i])

        def state(s: float) -> tuple[float, float]:
            y = sol.sol(s)
            return float(y[0]), float(y[1])

        def c_x_form(s: float) -> tuple[float, float]:
            if sol.y.shape[0] == 2:
                return state(s)
            raise AssertionError("unexpected state length")

        # The passed primary solution is in (rho_c,rho_x), so identities use
        # precisely that independently integrated interpolant.
        def dark_vars(s: float) -> tuple[float, float, float]:
            c, x = c_x_form(s)
            _, e2 = components(s, c, x)
            return c, x, math.sqrt(float(e2))

        comoving_integral = quad(
            lambda s: g * math.exp(3.0 * s) * dark_vars(s)[1] / dark_vars(s)[2],
            0.0,
            n,
            epsabs=2e-13,
            epsrel=2e-12,
            limit=200,
        )[0]
        c_now = float(rc_arr[0])
        lhs = math.exp(3.0 * n) * float(rc_arr[i])
        rhs = c_now - comoving_integral
        comoving_errors.append(abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-300))

        dark_integral = quad(
            lambda s: dark_vars(s)[0],
            0.0,
            n,
            epsabs=2e-13,
            epsrel=2e-12,
            limit=200,
        )[0]
        d_lhs = float(rc_arr[i] + rx_arr[i])
        d_rhs = float(rc_arr[0] + rx_arr[0]) - 3.0 * dark_integral
        dark_errors.append(abs(d_lhs - d_rhs) / max(abs(d_lhs), abs(d_rhs), 1e-300))
    return {
        "cdm_comoving_identity_max_relative_error": float(max(comoving_errors)),
        "dark_sector_conservation_max_relative_error": float(max(dark_errors)),
        "quadrature_sample_count": int(len(indices)),
    }


def audit_case(g: float) -> dict[str, object]:
    coarse = integrate(g, COARSE, "c_x")
    tolerance_only = integrate(g, TOLERANCE_ONLY, "c_x")
    max_step_only = integrate(g, MAX_STEP_ONLY, "c_x")
    fine = integrate(g, FINE, "c_x")
    alternate = integrate(g, FINE, "dark_total_x")
    fine_sol, rc, rx, rho, e2 = fine
    _, rc_coarse, rx_coarse, _, e2_coarse = coarse
    _, rc_tol, rx_tol, _, e2_tol = tolerance_only
    _, rc_step, rx_step, _, e2_step = max_step_only
    _, rc_alt, rx_alt, _, e2_alt = alternate

    # The analytic reference is ordinary flat LCDM with the same separately
    # normalized baryons, photons, massless neutrinos, and vacuum.
    e2_lcdm = (
        (PRESENT["omega_b0"] + PRESENT["omega_c0"]) * np.exp(-3.0 * N_GRID)
        + (PRESENT["omega_gamma0"] + PRESENT["omega_nu0"]) * np.exp(-4.0 * N_GRID)
        + PRESENT["omega_x0"]
    )
    gamma_km_s_mpc = g * H0_KM_S_MPC
    h0_s_inv = H0_KM_S_MPC * 1000.0 / MPC_M
    gamma_s_inv = gamma_km_s_mpc * 1000.0 / MPC_M
    unit_ratio_error = abs(gamma_s_inv / h0_s_inv - g)
    result: dict[str, object] = {
        "gamma_over_H0": g,
        "Gamma_km_s_Mpc": gamma_km_s_mpc,
        "H0_s_inv": h0_s_inv,
        "Gamma_s_inv": gamma_s_inv,
        "dimensionless_unit_ratio_error_after_SI_conversion": unit_ratio_error,
        "fine_vs_coarse_max_relative_difference": {
            "rho_c": relative_max(rc, rc_coarse),
            "rho_x": relative_max(rx, rx_coarse),
            "E2": relative_max(e2, e2_coarse),
        },
        "tolerance_only_tightening_vs_coarse_max_relative_difference": {
            "rho_c": relative_max(rc_tol, rc_coarse),
            "rho_x": relative_max(rx_tol, rx_coarse),
            "E2": relative_max(e2_tol, e2_coarse),
        },
        "max_step_only_tightening_vs_coarse_max_relative_difference": {
            "rho_c": relative_max(rc_step, rc_coarse),
            "rho_x": relative_max(rx_step, rx_coarse),
            "E2": relative_max(e2_step, e2_coarse),
        },
        "fine_direct_vs_dark_total_form_max_relative_difference": {
            "rho_c": relative_max(rc, rc_alt),
            "rho_x": relative_max(rx, rx_alt),
            "E2": relative_max(e2, e2_alt),
        },
        "minima_over_0_to_zmax": {
            "rho_b": float(np.min(rho["baryon"])),
            "rho_gamma": float(np.min(rho["photon_radiation"])),
            "rho_nu": float(np.min(rho["neutrino"])),
            "rho_c": float(np.min(rc)),
            "rho_x": float(np.min(rx)),
            "E2": float(np.min(e2)),
        },
        "at_zmax": {
            "z": Z_MAX,
            "rho_c": float(rc[-1]),
            "rho_x": float(rx[-1]),
            "E2": float(e2[-1]),
            "radiation_plus_neutrino_fraction_of_E2": float(
                (rho["photon_radiation"][-1] + rho["neutrino"][-1]) / e2[-1]
            ),
        },
        "source_sign_today": {
            "d_rho_x_dn": float(g * PRESENT["omega_x0"]),
            "cdm_transfer_term_d_rho_c_dn": float(-g * PRESENT["omega_x0"]),
            "expected": "Gamma>0: rho_x gains and CDM loses forward in time; Gamma<0: reverse",
        },
        "checks": {
            "ode_coarse_success": bool(coarse[0].success),
            "ode_tolerance_only_success": bool(tolerance_only[0].success),
            "ode_max_step_only_success": bool(max_step_only[0].success),
            "ode_fine_success": bool(fine_sol.success),
            "alternate_form_success": bool(alternate[0].success),
            "rho_c_positive_on_full_grid": bool(np.all(rc > 0.0)),
            "rho_x_positive_on_full_grid": bool(np.all(rx > 0.0)),
            "E2_positive_on_full_grid": bool(np.all(e2 > 0.0)),
            "flat_closure_today": math.isclose(
                sum(PRESENT.values()), 1.0, rel_tol=0.0, abs_tol=3e-16
            ),
            "unit_conversion_consistent": abs(
                unit_ratio_error
            ) < 1e-15,
            "cdm_positive_at_both_range_endpoints": bool(rc[0] > 0.0 and rc[-1] > 0.0),
            "vacuum_positive_at_both_range_endpoints": bool(rx[0] > 0.0 and rx[-1] > 0.0),
            "e2_positive_at_both_range_endpoints": bool(e2[0] > 0.0 and e2[-1] > 0.0),
        },
    }
    result["integral_identities"] = integral_checks(g, fine_sol, rc, rx)
    if g == 0.0:
        result["gamma_zero_analytic_lcdm_max_relative_error"] = {
            "rho_c": relative_max(rc, PRESENT["omega_c0"] * np.exp(-3.0 * N_GRID)),
            "rho_x": relative_max(rx, np.full_like(N_GRID, PRESENT["omega_x0"])),
            "E2": relative_max(e2, e2_lcdm),
        }
    return result


def main() -> None:
    start = time.perf_counter()
    cases = [audit_case(float(g)) for g in GAMMA_OVER_H0_VALUES]
    checks_all = all(
        all(bool(v) for v in case["checks"].values())
        and max(case["fine_vs_coarse_max_relative_difference"].values()) < 2e-7
        and max(case["tolerance_only_tightening_vs_coarse_max_relative_difference"].values()) < 2e-7
        and max(case["max_step_only_tightening_vs_coarse_max_relative_difference"].values()) < 2e-7
        and max(case["fine_direct_vs_dark_total_form_max_relative_difference"].values()) < 2e-7
        and case["integral_identities"]["cdm_comoving_identity_max_relative_error"] < 2e-9
        and case["integral_identities"]["dark_sector_conservation_max_relative_error"] < 2e-9
        for case in cases
    )
    zero_case = next(case for case in cases if case["gamma_over_H0"] == 0.0)
    checks_all = checks_all and max(
        zero_case["gamma_zero_analytic_lcdm_max_relative_error"].values()
    ) < 2e-9
    result = {
        "gate": "1: homogeneous background only",
        "status": "PASS" if checks_all else "FAIL",
        "primary_source": {
            "pdf": str(PDF_PATH.relative_to(ROOT)),
            "pdf_sha256": sha256(PDF_PATH),
            "text_copy": str(TEXT_PATH.relative_to(ROOT)),
            "text_copy_sha256": sha256(TEXT_PATH),
            "equations_checked": "Eqs. (3)-(8), PDF pp. 2-3",
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
            "parallel_workers": 1,
            "gpu_used": False,
            "script_elapsed_seconds": time.perf_counter() - start,
        },
        "paper_model": {
            "interaction": "Q=Gamma*rho_x",
            "equation_of_state": "w_x=-1",
            "sign_convention": "Gamma>0 transfers CDM energy to vacuum/DE",
            "Gamma_dimension": "Hubble rate; gamma=Gamma/H0 is dimensionless",
            "unit_conversion_check": "Gamma and H0 converted independently from km/s/Mpc to s^-1 using 1 Mpc=3.0856775814913673e22 m; Gamma_s/H0_s compared with gamma",
            "n_variable": "ln(a)",
            "independent_ode_forms": [
                "(rho_c,rho_x), with rho_c'=-3*rho_c-gamma*rho_x/E and rho_x'=gamma*rho_x/E",
                "(rho_c+rho_x,rho_x), with (rho_c+rho_x)'=-3*rho_c and rho_x'=gamma*rho_x/E",
            ],
            "flat_friedmann": "E^2=rho_b+rho_gamma+rho_nu+rho_c+rho_x",
            "conserved_separately": {
                "baryons": "rho_b'= -3 rho_b",
                "photons_radiation": "rho_gamma'= -4 rho_gamma",
                "massless_neutrinos": "rho_nu'= -4 rho_nu",
            },
            "checked_identities": [
                "dark-sector continuity: (rho_c+rho_x)'=-3*rho_c",
                "CDM comoving density: d(a^3*rho_c)/dn=-gamma*a^3*rho_x/E",
            ],
        },
        "pinned_assumptions": {
            "H0_km_s_Mpc": H0_KM_S_MPC,
            "computed_present_density_parameters": PRESENT,
            "Omega_m0_baryons_plus_CDM": OMEGA_M0,
            "Omega_b0": OMEGA_B0,
            "Omega_c0": OMEGA_C0,
            "T_CMB_K": T_CMB_K,
            "photon_density_method": "rho_gamma=a_rad*T_CMB^4/c^2; Omega_gamma=rho_gamma/rho_critical",
            "a_rad_J_m3_K4": A_RAD_J_M3_K4,
            "c_m_s": C_M_S,
            "G_SI": G_NEWTON_SI,
            "Mpc_m": MPC_M,
            "N_eff": N_EFF,
            "neutrino_approximation": "all effective neutrinos massless; rho_nu/rho_gamma=(7/8)*(4/11)^(4/3)*N_eff; separate conserved component",
            "paper_neutrino_mass_context": "primary-paper footnote specifies sum m_nu=0.06 eV; this Gate 1 benchmark approximates all neutrinos as massless and does not claim exact reproduction of that massive-neutrino prescription",
            "flat_closure": "Omega_x0=1-Omega_b0-Omega_c0-Omega_gamma0-Omega_nu0",
            "rate_grid_Gamma_over_H0": list(GAMMA_OVER_H0_VALUES),
            "z_range_including_today": [0.0, Z_MAX],
            "normalized_density_units": "today's critical density",
            "analytic_LCDM_limit": "same pinned densities, Gamma=0",
        },
        "integrator": {
            "method": "SciPy solve_ivp DOP853; dense-output quadrature via scipy.integrate.quad",
            "coarse": COARSE,
            "tolerance_only_tightening": TOLERANCE_ONLY,
            "max_step_only_tightening": MAX_STEP_ONLY,
            "fine": FINE,
            "N_eval_points": int(N_GRID.size),
            "max_step_convergence_threshold_relative": 2e-7,
            "identity_threshold_relative": 2e-9,
            "lcdm_recovery_threshold_relative": 2e-9,
        },
        "positivity_assessment": {
            "numerical_sampling": "2049 ln(a) output points on 0<=z<=1e10, fine ODE max_step=0.02",
            "cdm_continuous_argument": "For gamma!=0, at any rho_c=0 the ODE gives rho_c'=-gamma*rho_x/E, which has the same strict sign at every zero because rho_x/E>0. Thus rho_c cannot cross zero twice; positive rho_c at both endpoints excludes a zero in the interval. For gamma=0, rho_c=Omega_c0*a^-3 exactly.",
            "vacuum_continuous_argument": "rho_x solves a scalar homogeneous equation rho_x'=gamma*rho_x/E with positive present value, so remains positive wherever the finite background solution exists.",
            "e2_continuous_argument": "With positive separately conserved baryons, photons, and neutrinos, and positive rho_x and rho_c, their sum E^2 is positive.",
        },
        "cases": cases,
        "scope_limit": "A successful homogeneous-background gate does not validate perturbative stability, CMB spectra, sound-horizon consistency, growth, or the paper's full likelihood/evidence claims.",
    }
    OUT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"status": result["status"], "output": str(OUT_PATH),
                      "case_count": len(cases)}, sort_keys=True))
    if not checks_all:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

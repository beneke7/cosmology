#!/usr/bin/env python3
"""Independent flat interacting-vacuum background audit.

Implements Yang et al. (2026), Eqs. (3)--(5), (8), for w_x=-1 in the
explicit all-massless-neutrino benchmark requested for Gate 1.  Densities are
in units of today's critical energy density.  The independent time coordinate
is s=ln(1+z)=-ln(a), so integration proceeds from today into the past.

This file has no dependency on work/compute15/background_gate.py.
"""
from __future__ import annotations

import os

# Bound this standalone audit to one CPU and one thread in numerical libraries.
for _thread_env in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_thread_env] = "1"
if hasattr(os, "sched_getaffinity"):
    _available = os.sched_getaffinity(0)
    os.sched_setaffinity(0, {min(_available)})

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
from scipy.constants import c as C_M_S, G as G_SI, Stefan_Boltzmann as SIGMA_SB


HERE = Path(__file__).resolve().parent
PROJECT = Path(__file__).resolve().parents[2]
PAPER = PROJECT / "context/papers/interacting_de_desi_dr2_2026.pdf"
GATE_JSON = PROJECT / "work/compute15/background_gate.json"
OUT = HERE / "independent_background.json"

# Pinned benchmark (Omega_m means baryons + CDM here; massive neutrino density
# is not folded into Omega_m because this benchmark treats neutrinos as massless).
H0_KM_S_MPC = 67.4
OMEGA_M0 = 0.315
OMEGA_B0 = 0.049
OMEGA_C0 = OMEGA_M0 - OMEGA_B0
TCMB_K = 2.7255
NEFF = 3.044
ZMAX = 1.0e10
GAMMAS = (0.0, 0.46666, -0.46666)
GRID_SIZE = 2049
S_END = math.log1p(ZMAX)
S_GRID = np.linspace(0.0, S_END, GRID_SIZE)

MPC_M = 3.0856775814913673e22
H0_S = H0_KM_S_MPC * 1000.0 / MPC_M
RHOCRIT_ENERGY_J_M3 = 3.0 * H0_S**2 * C_M_S**2 / (8.0 * math.pi * G_SI)
RADIATION_CONSTANT_J_M3_K4 = 4.0 * SIGMA_SB / C_M_S
PHOTON_ENERGY_J_M3 = RADIATION_CONSTANT_J_M3_K4 * TCMB_K**4
OMEGA_GAMMA0 = PHOTON_ENERGY_J_M3 / RHOCRIT_ENERGY_J_M3
NU_TO_GAMMA = (7.0 / 8.0) * (4.0 / 11.0)**(4.0 / 3.0) * NEFF
OMEGA_NU0 = NU_TO_GAMMA * OMEGA_GAMMA0
OMEGA_R0 = OMEGA_GAMMA0 + OMEGA_NU0
OMEGA_X0 = 1.0 - OMEGA_M0 - OMEGA_R0

# solve_ivp settings: independent coarse and fine resolutions, plus isolated
# maximum-step and error-tolerance refinements.
COARSE = {"rtol": 1e-9, "atol": 1e-12, "max_step": 0.10}
TOL_ONLY = {"rtol": 1e-12, "atol": 1e-14, "max_step": 0.10}
STEP_ONLY = {"rtol": 1e-9, "atol": 1e-12, "max_step": 0.02}
FINE = {"rtol": 1e-12, "atol": 1e-14, "max_step": 0.02}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def e2_at(s: float, rc: float, rx: float) -> float:
    """Flat Friedmann constraint, Eq. (3), written with a=exp(-s)."""
    return OMEGA_B0 * math.exp(3.0 * s) + OMEGA_R0 * math.exp(4.0 * s) + rc + rx


def rhs(gamma: float):
    """Continuity equations in s=-ln(a); E=sqrt(H^2/H0^2)."""
    def derivative(s: float, y: np.ndarray) -> np.ndarray:
        rc, rx = (float(y[0]), float(y[1]))
        e2 = e2_at(s, rc, rx)
        if not math.isfinite(e2) or e2 <= 0.0:
            raise FloatingPointError(f"invalid E^2={e2} at s={s}")
        transfer = gamma * rx / math.sqrt(e2)
        # Paper: rho_dot_x=Q, rho_dot_c+3H*rho_c=-Q, Q=Gamma*rho_x.
        # Since ds/dt=-H, the signs below reverse with the coordinate.
        return np.array([3.0 * rc + transfer, -transfer])
    return derivative


def integrate(gamma: float, cfg: dict[str, float]):
    solution = solve_ivp(
        rhs(gamma), (0.0, S_END), (OMEGA_C0, OMEGA_X0),
        method="DOP853", t_eval=S_GRID, dense_output=True,
        rtol=cfg["rtol"], atol=cfg["atol"], max_step=cfg["max_step"],
    )
    if not solution.success or solution.y.shape[1] != GRID_SIZE:
        raise RuntimeError(f"ODE failure for gamma={gamma}: {solution.message}")
    rc, rx = solution.y
    b = OMEGA_B0 * np.exp(3.0 * S_GRID)
    rad = OMEGA_R0 * np.exp(4.0 * S_GRID)
    e2 = b + rad + rc + rx
    return solution, rc, rx, e2, b, rad


def relative_max(a, b) -> float:
    aa, bb = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    denom = np.maximum(np.maximum(np.abs(aa), np.abs(bb)), 1e-300)
    return float(np.max(np.abs(aa - bb) / denom))


def integral_audits(sol, gamma: float) -> dict[str, float]:
    """Quadrature-check dark-sector continuity and CDM comoving transfer."""
    points = np.linspace(0.0, S_END, 65)
    continuity_rel, continuity_abs = [], []
    comoving_rel, comoving_abs = [], []

    def state(s: float):
        rc, rx = sol.sol(s)
        return float(rc), float(rx)

    def e2_sol(s: float, rc: float, rx: float) -> float:
        return e2_at(s, rc, rx)

    for s in points:
        rc_s, rx_s = state(float(s))
        integral_rc = quad(lambda u: state(u)[0], 0.0, float(s),
                           epsabs=1e-11, epsrel=2e-11, limit=300)[0]
        residual = rc_s + rx_s - (OMEGA_C0 + OMEGA_X0) - 3.0 * integral_rc
        scale = max(abs(rc_s + rx_s), abs(OMEGA_C0 + OMEGA_X0),
                    abs(3.0 * integral_rc), 1e-300)
        continuity_abs.append(abs(residual))
        continuity_rel.append(abs(residual) / scale)

        transfer_int = quad(
            lambda u: math.exp(-3.0 * u) * state(u)[1]
            / math.sqrt(e2_sol(u, *state(u))),
            0.0, float(s), epsabs=1e-13, epsrel=2e-11, limit=300,
        )[0]
        residual_f = math.exp(-3.0 * float(s)) * rc_s - OMEGA_C0 - gamma * transfer_int
        scale_f = max(abs(math.exp(-3.0 * float(s)) * rc_s), OMEGA_C0,
                       abs(gamma * transfer_int), 1e-300)
        comoving_abs.append(abs(residual_f))
        comoving_rel.append(abs(residual_f) / scale_f)

    return {
        "samples": int(len(points)),
        "dark_sector_continuity_max_abs": float(max(continuity_abs)),
        "dark_sector_continuity_max_rel": float(max(continuity_rel)),
        "cdm_comoving_transfer_max_abs": float(max(comoving_abs)),
        "cdm_comoving_transfer_max_rel": float(max(comoving_rel)),
    }


def compare_gate(case: dict, gamma: float) -> dict:
    """Read only the prior JSON measurements; never import its implementation."""
    if not GATE_JSON.exists():
        return {"available": False}
    gate = json.loads(GATE_JSON.read_text())
    gcase = next(c for c in gate["cases"] if c["gamma_over_H0"] == gamma)
    own = case["fine_values_at_reference_z"]
    other = gcase["at_zmax"]
    differences = {}
    for key in ("E2", "rho_c", "rho_x"):
        x, y = own[key], other[key]
        differences[key] = {
            "absolute": abs(x - y),
            "relative": abs(x - y) / max(abs(x), abs(y), 1e-300),
        }
    return {
        "available": True,
        "reference_script_sha256": gate["implementation"]["script_sha256"],
        "differences_at_zmax": differences,
    }


def run_case(gamma: float) -> dict:
    coarse = integrate(gamma, COARSE)
    tol_only = integrate(gamma, TOL_ONLY)
    step_only = integrate(gamma, STEP_ONLY)
    fine = integrate(gamma, FINE)
    _, rc, rx, e2, b, rad = fine

    analytic_c = OMEGA_C0 * np.exp(3.0 * S_GRID)
    analytic_x = np.full_like(S_GRID, OMEGA_X0)
    analytic_e2 = b + rad + analytic_c + analytic_x
    analytic_error = {
        "rho_c_max_relative": relative_max(rc, analytic_c),
        "rho_x_max_relative": relative_max(rx, analytic_x),
        "E2_max_relative": relative_max(e2, analytic_e2),
    } if gamma == 0.0 else None

    zmax_result = {
        "z": ZMAX,
        "E2": float(e2[-1]),
        "rho_c": float(rc[-1]),
        "rho_x": float(rx[-1]),
        "radiation_fraction_of_E2": float(rad[-1] / e2[-1]),
    }
    h0_rate = H0_KM_S_MPC
    out = {
        "gamma_over_H0": gamma,
        "Gamma_km_s_Mpc": gamma * h0_rate,
        "source_sign_today": {
            "d_rho_x_dln_a": gamma * OMEGA_X0,
            "interaction_part_d_rho_c_dln_a": -gamma * OMEGA_X0,
            "interpretation": "Gamma>0 transfers CDM to vacuum forward in time",
        },
        "minima_on_2049_point_grid": {
            "rho_c": float(np.min(rc)),
            "rho_x": float(np.min(rx)),
            "E2": float(np.min(e2)),
        },
        "positivity_claim": "sampled only on the 2049-point output grid",
        "fine_values_at_reference_z": zmax_result,
        "lcdm_analytic_recovery": analytic_error,
        "convergence_relative_to_fine": {
            "coarse": {"rho_c": relative_max(coarse[1], rc),
                       "rho_x": relative_max(coarse[2], rx),
                       "E2": relative_max(coarse[3], e2)},
            "tolerance_only": {"rho_c": relative_max(tol_only[1], rc),
                               "rho_x": relative_max(tol_only[2], rx),
                               "E2": relative_max(tol_only[3], e2)},
            "max_step_only": {"rho_c": relative_max(step_only[1], rc),
                              "rho_x": relative_max(step_only[2], rx),
                              "E2": relative_max(step_only[3], e2)},
        },
        "integral_identities": integral_audits(fine[0], gamma),
    }
    out["compute15_comparison"] = compare_gate(out, gamma)
    return out


def main() -> None:
    start = time.perf_counter()
    closure = OMEGA_B0 + OMEGA_C0 + OMEGA_GAMMA0 + OMEGA_NU0 + OMEGA_X0
    if abs(closure - 1.0) > 5e-16:
        raise RuntimeError(f"flat closure mismatch: {closure}")
    if not all(math.isfinite(q) and q > 0.0 for q in
               (OMEGA_B0, OMEGA_C0, OMEGA_GAMMA0, OMEGA_NU0, OMEGA_X0)):
        raise RuntimeError("non-positive or non-finite present component")

    affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    report = {
        "assessment": "independent background cross-check; not a likelihood reproduction",
        "source": {
            "paper": str(PAPER.relative_to(PROJECT)),
            "pdf_sha256": digest(PAPER),
            "equations": "flat Friedmann Eq. (3), continuity Eqs. (4)-(5), Q=Gamma*rho_x Eq. (8), vacuum w_x=-1",
            "page_2_3_neutrino_note": "paper footnote fixes sum m_nu=0.06 eV; current benchmark uses massless neutrinos, so this is not the exact paper neutrino background",
        },
        "parameters": {
            "H0_km_s_Mpc": H0_KM_S_MPC,
            "Omega_m0_baryon_plus_CDM": OMEGA_M0,
            "Omega_b0": OMEGA_B0,
            "Omega_c0": OMEGA_C0,
            "Tcmb_K": TCMB_K,
            "Neff": NEFF,
            "neutrino_approximation": "massless; rho_nu proportional to a^-4",
            "z_max": ZMAX,
            "gamma_over_H0": list(GAMMAS),
            "flat_closure_sum": closure,
        },
        "radiation_normalization_independent_SI": {
            "H0_per_s": H0_S,
            "critical_energy_density_J_m3": RHOCRIT_ENERGY_J_M3,
            "blackbody_radiation_constant_J_m3_K4": RADIATION_CONSTANT_J_M3_K4,
            "photon_energy_density_J_m3": PHOTON_ENERGY_J_M3,
            "omega_gamma0": OMEGA_GAMMA0,
            "nu_to_gamma_ratio": NU_TO_GAMMA,
            "omega_nu0": OMEGA_NU0,
            "omega_r0": OMEGA_R0,
            "omega_x0_from_flatness": OMEGA_X0,
            "units": "Gamma has H units; gamma=Gamma/H0 is dimensionless; all densities are energy densities divided by present critical energy density",
        },
        "integrator": {
            "method": "SciPy solve_ivp DOP853; dense output; independent s=ln(1+z) coordinate",
            "output_points": GRID_SIZE,
            "coarse": COARSE,
            "tolerance_only": TOL_ONLY,
            "max_step_only": STEP_ONLY,
            "fine": FINE,
            "quadrature": "adaptive scipy.integrate.quad on dense interpolant at 65 points",
        },
        "execution": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "cpu_affinity_ids": affinity,
            "gpu_used": False,
            "elapsed_seconds": None,
            "code_path": str(Path(__file__).relative_to(PROJECT)),
            "code_sha256": digest(Path(__file__).resolve()),
            "run_command": ".venv/bin/python work/compute17/independent_background.py",
        },
        "cases": [run_case(gamma) for gamma in GAMMAS],
    }
    report["execution"]["elapsed_seconds"] = time.perf_counter() - start
    OUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(OUT.relative_to(PROJECT)),
        "elapsed_seconds": report["execution"]["elapsed_seconds"],
        "radiation": report["radiation_normalization_independent_SI"],
        "cases": [{
            "gamma": c["gamma_over_H0"],
            "zmax": c["fine_values_at_reference_z"],
            "minima": c["minima_on_2049_point_grid"],
            "lcdm": c["lcdm_analytic_recovery"],
            "coarse_vs_fine": c["convergence_relative_to_fine"]["coarse"],
            "integrals": c["integral_identities"],
            "gate_comparison": c["compute15_comparison"],
        } for c in report["cases"]],
        "code_sha256": report["execution"]["code_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()

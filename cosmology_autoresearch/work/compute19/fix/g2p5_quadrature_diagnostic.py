#!/usr/bin/env python3
"""Isolate the g=2.5 comoving-CDM integral validation residual."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad, solve_ivp

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
PDF = PROJECT / "context/papers/interacting_de_desi_dr2_2026.pdf"
G = 2.5
H0 = 67.4
H = H0 / 100.0
TCMB = 2.7255
NEFF = 3.044
OM = 0.315
OB = 0.049
OC = OM - OB
OG = (2.469e-5 * (TCMB / 2.7255) ** 4) / H**2
ONU = OG * (7 / 8) * (4 / 11) ** (4 / 3) * NEFF
ORAD = OG + ONU
OX = 1 - OM - ORAD
UMAX = math.log1p(1e10)


def e2(u, x, q):
    e3 = math.exp(3 * u)
    return ORAD * e3 * math.exp(u) + OB * e3 + q * e3 + x


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    start = time.monotonic()

    def rhs(u, y):
        x, q, integ = y
        ee = e2(u, x, q)
        if ee <= 0 or not math.isfinite(ee):
            raise ArithmeticError(f"invalid E2={ee} at u={u}")
        f = x / math.sqrt(ee) * math.exp(-3 * u)
        return (-G * x / math.sqrt(ee), G * f, f)

    sol = solve_ivp(
        rhs, (0, UMAX), (OX, OC, 0.0), method="DOP853",
        rtol=3e-12, atol=2e-14, max_step=0.02, dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    x_end, q_end, i_ode = map(float, sol.y[:, -1])
    dense_end = tuple(float(v) for v in sol.sol(UMAX))
    ode_identity_end = abs(q_end - (OC + G * i_ode))
    ode_identity_dense_max = 0.0
    for u in np.linspace(0, UMAX, 4097):
        x, q, i_state = map(float, sol.sol(float(u)))
        ode_identity_dense_max = max(ode_identity_dense_max, abs(q - (OC + G * i_state)))

    simpson_rows = []
    for nodes in (2049, 4097, 8193, 16385):
        us = np.linspace(0, UMAX, nodes)
        state = sol.sol(us)
        xs, qs = state[0], state[1]
        e2s = np.array([e2(float(u), float(x), float(q)) for u, x, q in zip(us, xs, qs)])
        fs = xs / np.sqrt(e2s) * np.exp(-3 * us)
        n_intervals = nodes - 1
        du = UMAX / n_intervals
        integral = float(du / 3 * (fs[0] + fs[-1] + 4 * np.sum(fs[1:-1:2]) + 2 * np.sum(fs[2:-1:2])))
        q_from_int = OC + G * integral
        simpson_rows.append({
            "nodes": nodes,
            "intervals": n_intervals,
            "du": du,
            "integral": integral,
            "q_from_simpson": q_from_int,
            "absolute_q_residual": abs(q_end - q_from_int),
            "q_residual_vs_ode_I": abs(q_from_int - (OC + G * i_ode)),
        })

    def source(u):
        x, q, _ = sol.sol(u)
        ee = e2(u, float(x), float(q))
        if ee <= 0:
            raise ArithmeticError(f"invalid E2 in quad source at u={u}")
        return float(x) / math.sqrt(ee) * math.exp(-3 * u)

    quad_value, quad_error, quad_info = quad(
        source, 0.0, UMAX, epsabs=1e-13, epsrel=1e-13, limit=500, full_output=1
    )[:3]
    q_from_quad = OC + G * quad_value
    result = {
        "task": "diagnose g=2.5 quadrature residual only; no threshold/domain inference",
        "model": {
            "g": G,
            "Omega_m0": OM,
            "Omega_b0": OB,
            "Omega_c0": OC,
            "H0_km_s_Mpc": H0,
            "Tcmb_K": TCMB,
            "Neff_massless": NEFF,
            "Omega_gamma0": OG,
            "Omega_nu0": ONU,
            "Omega_rad0": ORAD,
            "Omega_x0": OX,
            "zmax": 1e10,
            "u_max": UMAX,
        },
        "ode": {
            "method": "SciPy solve_ivp DOP853",
            "rtol": 3e-12,
            "atol": 2e-14,
            "max_step_u": 0.02,
            "success": True,
            "nfev": int(sol.nfev),
            "dense_endpoint_state_x_q_I": dense_end,
            "ODE_identity_endpoint_abs_residual": ode_identity_end,
            "ODE_identity_dense_4097_sample_max_abs_residual": ode_identity_dense_max,
        },
        "simpson_convergence": simpson_rows,
        "adaptive_quad": {
            "integral": quad_value,
            "reported_error_estimate": quad_error,
            "function_evaluations": int(quad_info["neval"]),
            "q_from_integral": q_from_quad,
            "absolute_q_residual_vs_ode_endpoint": abs(q_end - q_from_quad),
            "absolute_integral_difference_vs_ODE_state_I": abs(quad_value - i_ode),
        },
        "equation_unit_identity": {
            "backward_u": "dx/du=-g*x/E; dm/du=3*m+g*x/E, with q=m*exp(-3u)",
            "comoving_identity": "dq/du=g*x/E*exp(-3u); q=q0+g*integral(x/E*exp(-3u)du)",
            "same_integrand_used_by_ode_and_adaptive_quad": True,
        },
        "diagnostic_runtime_seconds": time.monotonic() - start,
        "sha256": {
            "diagnostic_script": sha(Path(__file__)),
            "local_paper_pdf": sha(PDF),
        },
    }
    out = HERE / "g2p5_quadrature_diagnostic.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(out),
        "simpson_residuals": [r["absolute_q_residual"] for r in simpson_rows],
        "quad_residual": result["adaptive_quad"]["absolute_q_residual_vs_ode_endpoint"],
        "ode_state_residual": ode_identity_end,
        "runtime_seconds": result["diagnostic_runtime_seconds"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

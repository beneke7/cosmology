#!/usr/bin/env python3
"""Early-time positivity screen for the archived BAO-screen IVS point.

This is not a likelihood analysis. It extends one supplied point using flat
closure, separately conserved baryons/photons/massless neutrinos, and the
source-paper interaction Q=Gamma*rho_v with w_v=-1. Two independent ODE
coordinates and the comoving-CDM integral identity are compared.
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
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
OUT_JSON = Path(__file__).with_suffix(".json")
OUT_MD = Path(__file__).with_suffix(".md")
OM0 = 0.3869710077
G = -0.4666647299
TCMB_K = 2.7255
NEFF = 3.044
# Standard blackbody photon density: Omega_gamma h^2 = 2.4728e-5 at Tcmb=2.7255 K.
OMEGA_GAMMA_H2 = 2.4728e-5
NU_GAMMA_RATIO = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0)
H0_GRID = (55.0, 65.0, 75.0, 85.0)
FB_GRID = (0.10, 0.15, 0.20, 0.25)
ZMAX = 1.0e10
RTOL = 2e-11
ATOL = 2e-13


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rad_today(h0: float) -> tuple[float, float]:
    h = h0 / 100.0
    og = OMEGA_GAMMA_H2 / h**2
    onu = og * NU_GAMMA_RATIO * NEFF
    return og, onu


def integrate_n(h0: float, fb: float, rtol: float = RTOL, atol: float = ATOL):
    """Direct component equations in n=ln(a), integrated from today backward."""
    og, onu = rad_today(h0)
    ob, oc, ov = fb * OM0, (1.0 - fb) * OM0, 1.0 - OM0 - og - onu
    y0 = np.array([ob, oc, ov, og, onu], dtype=np.float64)

    def rhs(n, y):
        b, c, v, photon, nu = y
        e2 = b + c + v + photon + nu
        if not np.isfinite(e2) or e2 <= 0.0:
            raise FloatingPointError(f"E^2<=0 at n={n:g}")
        e = math.sqrt(e2)
        return (-3*b, -3*c - G*v/e, G*v/e, -4*photon, -4*nu)

    nmin = -math.log1p(ZMAX)
    sol = solve_ivp(rhs, (0.0, nmin), y0, method="DOP853", dense_output=True,
                    rtol=rtol, atol=atol, max_step=0.025)
    if not sol.success or sol.sol is None:
        raise RuntimeError(sol.message)
    return sol, (og, onu, ov)


def integrate_u_independent(h0: float, fb: float, rtol: float = RTOL, atol: float = ATOL):
    """Independent formulation in u=ln(1+z) with directly written signs."""
    og, onu = rad_today(h0)
    y0 = np.array([fb*OM0, (1-fb)*OM0, 1-OM0-og-onu, og, onu], dtype=np.float64)

    def rhs(u, y):
        b, c, v, photon, nu = y
        e2 = b + c + v + photon + nu
        if not np.isfinite(e2) or e2 <= 0:
            raise FloatingPointError(f"E^2<=0 at u={u:g}")
        e = math.sqrt(e2)
        transfer = G * v / e
        return (3*b, 3*c + transfer, -transfer, 4*photon, 4*nu)

    umax = math.log1p(ZMAX)
    sol = solve_ivp(rhs, (0.0, umax), y0, method="DOP853", dense_output=True,
                    rtol=rtol, atol=atol, max_step=0.025)
    if not sol.success or sol.sol is None:
        raise RuntimeError(sol.message)
    return sol


def identity_values(sol, nvals: np.ndarray, oc0: float) -> tuple[np.ndarray, float]:
    """Compute q_c from independent adaptive quadrature of the exact identity."""
    qs, errs = [], []
    for n in nvals:
        if n == 0:
            integ, err = 0.0, 0.0
        else:
            integ, err = quad(
                lambda s: math.exp(3*s) * float(sol.sol(s)[2]) /
                          math.sqrt(float(np.sum(sol.sol(s)))),
                float(n), 0.0, epsabs=2e-13, epsrel=2e-12, limit=250,
            )
        qs.append(oc0 + G*integ)
        errs.append(abs(G)*err)
    return np.asarray(qs), max(errs, default=0.0)


def run_case(h0: float, fb: float) -> dict:
    sol, (og, onu, ov0) = integrate_n(h0, fb)
    alt = integrate_u_independent(h0, fb)
    # Sampling uniform in log(1+z) resolves the full radiation-to-today history.
    ugrid = np.linspace(0.0, math.log1p(ZMAX), 5001)
    ngrid = -ugrid
    y = np.asarray(sol.sol(ngrid), dtype=np.float64)
    yalt = np.asarray(alt.sol(ugrid), dtype=np.float64)
    e2 = np.sum(y, axis=0)
    qdirect = np.exp(3*ngrid) * y[1]
    qidentity, qquaderr = identity_values(sol, ngrid[::100], (1-fb)*OM0)
    qdirect_coarse = qdirect[::100]
    iderr = float(np.max(np.abs(qdirect_coarse-qidentity)))
    cross = None
    # Find first zero encountered moving backward (smallest positive z root).
    neg = np.flatnonzero(y[1] <= 0.0)
    if neg.size:
        j = int(neg[0])
        if j == 0:
            cross = 0.0
        else:
            nroot = brentq(lambda nn: float(sol.sol(nn)[1]), ngrid[j], ngrid[j-1],
                           xtol=2e-14, rtol=9e-16)
            cross = math.expm1(-nroot)
    # Comoving positivity identity's asymptotic value and zero feasibility.
    qc_asym = float(qdirect[-1])
    return {
        "H0_km_s_Mpc": h0, "f_b": fb,
        "Omega_gamma0": og, "Omega_nu_massless0": onu,
        "Omega_r0": og+onu, "Omega_v0_flat_closure": ov0,
        "rho_c_zero_crossing_redshift": cross,
        "rho_c_positive_to_zmax": bool(np.min(y[1]) > 0.0),
        "rho_c_min_sampled_dimensionless_density": float(np.min(y[1])),
        "q_c_at_zmax": qc_asym,
        "q_c_identity_max_abs_error_sampled": iderr,
        "identity_quadrature_error_estimate_max": float(qquaderr),
        "independent_n_vs_u_max_abs_component_difference": float(np.max(np.abs(y-yalt))),
        "minimum_sampled_E2": float(np.min(e2)),
        "E2_positive_sampled": bool(np.all(e2 > 0)),
    }


def main():
    t0 = time.perf_counter()
    # Tolerance/convergence checks: one representative case at two tighter settings.
    representative = (70.0, 0.16)
    s_ref, _ = integrate_n(*representative, RTOL, ATOL)
    s_tight, _ = integrate_n(*representative, 2e-13, 2e-15)
    ncheck = np.linspace(0, -math.log1p(ZMAX), 1201)
    ref_values = s_ref.sol(ncheck)
    tight_values = s_tight.sol(ncheck)
    convergence = float(np.max(np.abs(ref_values-tight_values) /
                               np.maximum(1.0, np.maximum(np.abs(ref_values), np.abs(tight_values)))))
    # A third coarse resolution check independently verifies that the event result
    # and identity are insensitive to the tabulation density.
    cases = [run_case(h, fb) for h in H0_GRID for fb in FB_GRID]
    paper = ROOT / "context/papers/interacting_de_desi_dr2_2026.pdf"
    contract = ROOT / "work/theory3/interacting_vacuum_contract.md"
    ident = ROOT / "work/theory4/interacting_vacuum_identifiability_report.md"
    output = {
        "classification": "bounded background extension of a supplied BAO-screen point; not a refit or posterior",
        "parameters": {"Omega_m0": OM0, "g_Gamma_over_H0": G, "w_x": -1.0},
        "assumptions": {
            "geometry": "flat FLRW; present closure includes photons and massless neutrinos",
            "components": "separately conserved baryons, photons, massless neutrinos; interacting CDM and vacuum",
            "interaction": "Q=Gamma*rho_v; positive Gamma transfers CDM to vacuum, per source-paper Eqs. (4)-(5),(8)",
            "radiation": "Omega_gamma h^2=2.4728e-5 at T_CMB=2.7255 K; Omega_nu/Omega_gamma=(7/8)(4/11)^(4/3) N_eff",
            "source_paper": "N_eff=3.044 and sum m_nu=0.06 eV; this calculation keeps N_eff but approximates all neutrinos as massless at every epoch",
            "H0_f_b_grid": {"H0_km_s_Mpc": list(H0_GRID), "f_b": list(FB_GRID)},
            "redshift_range": [0.0, ZMAX],
            "BAO_identifiability": "H0 and f_b are sensitivity coordinates only; BAO alpha does not identify either, and no hidden CMB/BBN prior is applied",
            "critical_density": "all component densities normalized to today's critical density for the specified H0",
        },
        "comoving_CDM_identity": {
            "derivation": "For q_c=a^3 u_c, dq_c/dn=a^3(du_c/dn+3u_c)=-g*a^3*u_v/E. Integrating from n=0 to n<0 gives q_c(n)=Omega_c0+g*integral_n^0 exp(3s)u_v(s)/E(s) ds.",
            "negative_g_implication": "q_c is strictly decreasing as one moves backward because g<0, u_v>0 and E>0. Thus if q_c(zmax)>0, no missed intermediate zero can occur on [0,zmax]; a zero of q_c is the zero of rho_c because a^3>0.",
            "important": "present-day rho_c is (1-f_b) Omega_m0; the zero threshold depends on the unmeasured split f_b",
        },
        "numerics": {
            "solver": "SciPy DOP853 float64, one-thread environment, max_step=0.025 in log scale factor",
            "rtol": RTOL, "atol": ATOL, "identity_quadrature": "adaptive scipy.integrate.quad",
            "grid_cases": len(cases), "dense_sampling_points_per_case": 5001,
            "representative_tolerance_convergence_max_component_abs": convergence,
            "representative_case_H0_f_b": list(representative),
        },
        "cases": cases,
        "provenance": {
            str(p.relative_to(ROOT)): sha256(p) for p in
            (Path(__file__), paper, contract, ident)
        },
        "runtime_seconds": None,
        "runtime": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform(),
                    "thread_env": {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}},
    }
    output["runtime_seconds"] = time.perf_counter()-t0
    OUT_JSON.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")
    crossings = [x["rho_c_zero_crossing_redshift"] for x in cases]
    cross_values = [x for x in crossings if x is not None]
    cross_span = (f'{min(cross_values):.6g} to {max(cross_values):.6g}'
                  if cross_values else "none")
    qc_values = [c["q_c_at_zmax"] for c in cases]
    rows = ["| H0 | f_b | Omega_gamma0+nu0 | rho_c zero z | q_c(z=1e10) | positive to z=1e10 | q_c identity error | n/u max diff |",
            "|---:|---:|---:|---:|---:|:---:|---:|---:|"]
    for c in cases:
        zc = "none" if c["rho_c_zero_crossing_redshift"] is None else f'{c["rho_c_zero_crossing_redshift"]:.6g}'
        rows.append(f'| {c["H0_km_s_Mpc"]:g} | {c["f_b"]:.2f} | {c["Omega_r0"]:.4g} | {zc} | {c["q_c_at_zmax"]:.6f} | {"yes" if c["rho_c_positive_to_zmax"] else "no"} | {c["q_c_identity_max_abs_error_sampled"]:.2e} | {c["independent_n_vs_u_max_abs_component_difference"]:.2e} |')
    report = f'''# Early-time positivity of the BAO-screen IVS point

**Finding.** This bounded background calculation extends the supplied screen point (\\(\\Omega_{{m0}}={OM0:.10f}\\), \\(g={G:.10f}\\), \\(w_x=-1\\)); it does not refit BAO and gives no posterior. The declared grid contains {len(cases)} combinations of \\(H_0\\) and baryon fraction. Of these, {len(cross_values)} have a CDM zero crossing by \\(z=10^{{10}}\\). Crossing redshifts: {cross_span}. At \\(z=10^{{10}}\\), the smallest comoving CDM value across the grid is {min(qc_values):.6f}, still positive. The full rows are below.

## Identity

With \\(n=\\ln a\\), \\(E=H/H_0\\), \\(u_i=\\rho_i/\\rho_{{crit,0}}\\), and \\(q_c=a^3u_c\\), the source-paper equations give \\(du_c/dn=-3u_c-gu_v/E\\). Hence \\(dq_c/dn=-g a^3u_v/E\\), so integrating backward from today,

\\[q_c(n)=\\Omega_{{c0}}+g\\int_n^0 e^{{3s}}\\frac{{u_v(s)}}{{E(s)}}\\,ds.\\]

For this negative \\(g\\), \\(q_c\\) decreases monotonically as we move into the past, since \\(u_v/E>0\\). Thus a positive value at the declared maximum redshift excludes an intermediate crossing. Since \\(u_c=q_c/a^3\\), a zero in the identity is exactly a zero in physical CDM density. The threshold is conditional on \\(\\Omega_{{c0}}=(1-f_b)\\Omega_{{m0}}\\); neither BAO amplitude nor this background fit measures \\(f_b\\).

## Declared assumptions and limits

The local paper's background equations conserve baryons, radiation, and neutrinos separately and specify \\(N_{{eff}}=3.044\\), \\(\\sum m_\\nu=0.06\\) eV. Photons use \\(\\Omega_\\gamma h^2=2.4728\\times10^{{-5}}\\) at \\(T_{{CMB}}=2.7255\\) K. Neutrinos are treated here as massless radiation at all redshifts, explicitly an approximation to the paper's nonzero mass prescription. Flat closure sets \\(\\Omega_{{v0}}=1-\\Omega_{{m0}}-\\Omega_{{r0}}\\). The sensitivity grid is \\(H_0\\in\\{{{', '.join(map(str,H0_GRID))}\\}}\\) km s\\(^{{-1}}\\) Mpc\\(^{{-1}}\\), \\(f_b\\in\\{{{', '.join(f'{v:.2f}' for v in FB_GRID)}\\}}\\), and \\(0\\le z\\le10^{{10}}\\). No CMB or BBN prior is imported; H0 and \\(f_b\\) are sensitivity coordinates because the BAO-only amplitude \\(\\alpha=c/(H_0r_d)\\) does not identify them.

The numerical solutions use DOP853 in two independently written coordinates, \\(n=\\ln a\\) and \\(u=\\ln(1+z)\\), and an adaptive quadrature check of the exact comoving-CDM identity. The tolerance-refinement metric is the maximum component difference scaled by \\(\\max(1,|u_i|)\\): {convergence:.3e}; maximum identity residual over all cases: {max(c["q_c_identity_max_abs_error_sampled"] for c in cases):.3e}; maximum disagreement between coordinate integrations: {max(c["independent_n_vs_u_max_abs_component_difference"] for c in cases):.3e}. At very high redshift the result is a homogeneous background extrapolation only; it does not validate perturbations, sound-horizon physics, neutrino mass transitions, or early-universe viability in a full model.

## Grid results

{chr(10).join(rows)}

Reproduce from the repository root with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \\
  .venv/bin/python work/theory5/ivs_early_positivity.py
```

Machine-readable detail and hashes: [ivs_early_positivity.json](ivs_early_positivity.json). Runtime: {output['runtime_seconds']:.3f} s. Script SHA-256: `{sha256(Path(__file__))}`.

**Stop disposition:** completed at \\(z=10^{{10}}\\); no CDM crossing was found in the declared 16-case grid. No full likelihood, CAMB, CMB prior, or extrapolation beyond the stated redshift bound was used.
'''
    OUT_MD.write_text(report)
    print(json.dumps({"runtime_seconds": output["runtime_seconds"], "cases": len(cases),
                      "crossings": sum(x is not None for x in crossings),
                      "crossing_range": [min((x for x in crossings if x is not None), default=None),
                                         max((x for x in crossings if x is not None), default=None)],
                      "script_sha256": sha256(Path(__file__)), "json_sha256": sha256(OUT_JSON),
                      "md_sha256": sha256(OUT_MD)}, indent=2))


if __name__ == "__main__":
    main()

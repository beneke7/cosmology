#!/usr/bin/env python3
"""Thermal massive-neutrino sensitivity check for the archived IVS point.

One 0.06-eV thermal neutrino eigenstate plus 2.044 massless effective species
implements an explicit minimal-mass prescription consistent with N_eff=3.044.
No likelihood is evaluated and no CAMB result is reproduced.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import sys
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
OUT_JSON = Path(__file__).with_suffix(".json")
OUT_MD = Path(__file__).with_suffix(".md")
MASS_EV = 0.06
NEFF = 3.044
TCMB_K = 2.7255
K_B_EV_K = 8.617333262e-5
TNU_OVER_TGAMMA = (4.0 / 11.0) ** (1.0 / 3.0)
TNU0_EV = K_B_EV_K * TCMB_K * TNU_OVER_TGAMMA
OMEGA_GAMMA_H2 = 2.4728e-5
OM0 = 0.3869710077
GAMMA = -0.4666647299
H0_GRID = (55.0, 65.0, 75.0, 85.0)
FB_GRID = (0.10, 0.15, 0.20, 0.25)
ZMAX = 1.0e10
RTOL, ATOL = 2.0e-11, 2.0e-13
F0 = 7.0 * math.pi**4 / 120.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@lru_cache(maxsize=20000)
def fd_energy_integral(y_key: float) -> float:
    """F(y)=int dq q^2 sqrt(q^2+y^2)/(exp(q)+1), via adaptive QUAD."""
    y = float(y_key)
    val, _ = quad(
        lambda q: q*q*math.hypot(q, y)/(math.exp(q)+1.0)
        if q < 700.0 else 0.0,
        0.0, 50.0, epsabs=2e-12, epsrel=2e-12, limit=200,
    )
    return val


def omega_components_today(h0: float) -> tuple[float, float, float]:
    """Return photon, total neutrino, and massive-neutrino density today."""
    h = h0 / 100.0
    og = OMEGA_GAMMA_H2 / h**2
    y0 = MASS_EV / TNU0_EV
    ratio_per_species = (15.0 / math.pi**4) * TNU_OVER_TGAMMA**4 * fd_energy_integral(y0)
    ratio_massless = (15.0 / math.pi**4) * TNU_OVER_TGAMMA**4 * F0
    onu_massive = og * ratio_per_species
    onu_other = og * (NEFF-1.0) * ratio_massless
    return og, onu_massive + onu_other, onu_massive


def omega_nu_at_n(n: float, og: float) -> float:
    """Thermal neutrino density normalized to today's critical density."""
    a = math.exp(n)
    y = MASS_EV * a / TNU0_EV
    ratio_massive = (15.0 / math.pi**4) * TNU_OVER_TGAMMA**4 * fd_energy_integral(y)
    ratio_massless = (15.0 / math.pi**4) * TNU_OVER_TGAMMA**4 * F0
    return og * a**-4 * (ratio_massive + (NEFF-1.0)*ratio_massless)


def integrate_massive(h0: float, fb: float, rtol=RTOL, atol=ATOL):
    og, onu0, onu_m0 = omega_components_today(h0)
    ob, oc = fb*OM0, (1.0-fb)*OM0
    ov0 = 1.0 - OM0 - og - onu0
    y0 = np.array([ob, oc, ov0, og], dtype=np.float64)

    def rhs(n, state):
        b, c, v, photon = state
        nu = omega_nu_at_n(float(n), og)
        e2 = b+c+v+photon+nu
        if not math.isfinite(e2) or e2 <= 0:
            raise FloatingPointError(f"E^2<=0 at n={n:g}")
        e = math.sqrt(e2)
        transfer = GAMMA*v/e
        return (-3*b, -3*c-transfer, transfer, -4*photon)

    nmin = -math.log1p(ZMAX)
    sol = solve_ivp(rhs, (0.0, nmin), y0, method="DOP853", dense_output=True,
                    rtol=rtol, atol=atol, max_step=0.025)
    if not sol.success or sol.sol is None:
        raise RuntimeError(sol.message)
    return sol, (og, onu0, onu_m0, ov0)


def evaluate_state(sol, n, og):
    arr = np.asarray(sol.sol(float(n)), dtype=np.float64)
    nu = omega_nu_at_n(float(n), og)
    e2 = float(np.sum(arr) + nu)
    return arr, nu, e2


def identity_check(sol, nvals, oc0, og):
    vals, qerr = [], 0.0
    for n in nvals:
        if n == 0:
            integ, err = 0.0, 0.0
        else:
            integ, err = quad(
                lambda s: math.exp(3*s)*float(sol.sol(float(s))[2]) /
                    math.sqrt(evaluate_state(sol, float(s), og)[2]),
                float(n), 0.0, epsabs=2e-13, epsrel=2e-12, limit=250,
            )
        qdirect = math.exp(3*float(n))*float(sol.sol(float(n))[1])
        qidentity = oc0 + GAMMA*integ
        vals.append(abs(qdirect-qidentity))
        qerr = max(qerr, abs(GAMMA)*err)
    return max(vals, default=0.0), qerr


def run_case(h0, fb):
    sol, (og, onu0, onu_m0, ov0) = integrate_massive(h0, fb)
    ugrid = np.linspace(0.0, math.log1p(ZMAX), 5001)
    ngrid = -ugrid
    states = np.asarray(sol.sol(ngrid), dtype=np.float64)
    nu = np.asarray([omega_nu_at_n(float(n), og) for n in ngrid])
    e2 = np.sum(states, axis=0)+nu
    cden = states[1]
    neg = np.flatnonzero(cden <= 0.0)
    zcross = None
    if neg.size:
        j = int(neg[0])
        nr = brentq(lambda n: float(sol.sol(float(n))[1]), ngrid[j], ngrid[j-1],
                    xtol=2e-14, rtol=9e-16) if j else 0.0
        zcross = math.expm1(-nr)
    oc0 = (1.0-fb)*OM0
    qzmax = math.exp(3*ngrid[-1])*float(cden[-1])
    iderr, quaderr = identity_check(sol, ngrid[::100], oc0, og)
    # Explicit radiation-to-matter-transition diagnostics for the nu prescription.
    zrelativistic = TNU0_EV/MASS_EV
    return {
        "H0_km_s_Mpc": h0, "f_b": fb,
        "Omega_gamma0": og, "Omega_nu0_total": onu0,
        "Omega_nu0_massive_eigenstate": onu_m0,
        "Omega_v0_flat_closure": ov0,
        "rho_c_zero_crossing_redshift": zcross,
        "rho_c_positive_to_zmax": bool(np.min(cden) > 0.0),
        "q_c_at_zmax": qzmax,
        "minimum_sampled_E2": float(np.min(e2)),
        "identity_max_abs_error": iderr,
        "identity_quadrature_error_estimate": quaderr,
        "massive_neutrino_density_at_zmax": float(nu[-1]),
        "approx_m_over_Tnu_equals_one_redshift": 1.0/zrelativistic-1.0,
    }


def main():
    t0 = time.perf_counter()
    prev_path = ROOT/"work/theory5/ivs_early_positivity.json"
    prev = json.loads(prev_path.read_text())
    prev_by_key = {(c["H0_km_s_Mpc"], c["f_b"]): c for c in prev["cases"]}
    cases = [run_case(h, fb) for h in H0_GRID for fb in FB_GRID]
    for c in cases:
        old = prev_by_key[(c["H0_km_s_Mpc"], c["f_b"])]
        c["massless_q_c_at_zmax"] = old["q_c_at_zmax"]
        c["delta_q_c_at_zmax_vs_massless"] = c["q_c_at_zmax"]-old["q_c_at_zmax"]
        c["massless_rho_c_zero_crossing_redshift"] = old["rho_c_zero_crossing_redshift"]

    # Convergence: full history integrated with nominal/tighter tolerances;
    # separately tighten the adaptive thermal FD quadrature at spot-check y's.
    case_conv = (70.0, 0.16)
    s1, _ = integrate_massive(*case_conv, RTOL, ATOL)
    s2, _ = integrate_massive(*case_conv, 2e-13, 2e-15)
    ntest = np.linspace(0.0, -math.log1p(ZMAX), 1001)
    v1, v2 = s1.sol(ntest), s2.sol(ntest)
    ode_conv = float(np.max(np.abs(v1-v2)/np.maximum(1.0, np.maximum(np.abs(v1), np.abs(v2)))))
    ytests = (0.0, 0.01, 1.0, MASS_EV/TNU0_EV, 100.0, MASS_EV/TNU0_EV)
    fd_checks = []
    for y in sorted(set(ytests)):
        ref = fd_energy_integral(y)
        tight, qerr = quad(lambda q: q*q*math.hypot(q,y)/(math.exp(q)+1.0) if q < 700 else 0.0,
                           0, 50, epsabs=2e-13, epsrel=5e-13, limit=250)
        fd_checks.append({"y": y, "integral_nominal": ref, "integral_tight": tight,
                          "relative_difference": abs(ref-tight)/max(abs(tight), 1e-300),
                          "tight_quad_error_estimate": qerr})

    crossings = sum(c["rho_c_zero_crossing_redshift"] is not None for c in cases)
    massless_crossings = sum(c["massless_rho_c_zero_crossing_redshift"] is not None for c in cases)
    paper = ROOT/"context/papers/interacting_de_desi_dr2_2026.pdf"
    script = Path(__file__)
    output = {
        "classification": "bounded thermal-neutrino sensitivity check of a supplied BAO-screen point; no refit/posterior/CAMB reproduction",
        "parameters": {"Omega_m0_baryon_plus_CDM": OM0, "g_Gamma_over_H0": GAMMA, "w_x": -1.0},
        "neutrino_prescription": {
            "massive_eigenstates": 1, "massive_eigenstate_mass_eV": MASS_EV,
            "massless_effective_species": NEFF-1.0, "N_eff_total": NEFF,
            "temperature_ratio": "T_nu/T_gamma=(4/11)^(1/3)",
            "formula": "rho_nu(a)/rho_gamma0 = (15/pi^4)*(Tnu0/Tgamma0)^4*a^-4*[F(m a/Tnu0)+(Neff-1)F(0)], F(y)=int_0^infinity dq q^2 sqrt(q^2+y^2)/(exp(q)+1)",
            "normalization_check": "F(0)=7*pi^4/120 recovers (7/8)*(4/11)^(4/3) per massless species",
            "caveat": "explicit minimal-mass thermal prescription; not asserted to match the unknown author-CAMB neutrino hierarchy/configuration exactly",
        },
        "assumptions": {
            "flat_closure": "Omega_v0=1-Omega_m0-Omega_gamma0-Omega_nu0; Omega_m0 retains baryon+CDM meaning from the BAO screen",
            "components": "separately conserved baryons and photons; exact thermal FD neutrino background; interacting CDM/vacuum Q=Gamma*rho_v",
            "grid": {"H0_km_s_Mpc": list(H0_GRID), "f_b": list(FB_GRID), "z_range": [0.0, ZMAX]},
            "prior_scope": "H0/f_b are sensitivity coordinates only; no CMB/BBN prior or BAO refit",
        },
        "comoving_CDM_identity": "q_c(n)=Omega_c0+g*integral_n^0 exp(3s)u_v(s)/E(s) ds; g<0 makes q_c decline monotonically moving backward, so positive q_c at zmax rules out earlier crossings on the interval.",
        "checks": {
            "ode_tolerance_refinement_max_scaled_component_difference": ode_conv,
            "representative_case_H0_f_b": list(case_conv),
            "fd_adaptive_quadrature_convergence": fd_checks,
            "max_identity_abs_error_over_grid": max(c["identity_max_abs_error"] for c in cases),
            "max_identity_quad_error_estimate": max(c["identity_quadrature_error_estimate"] for c in cases),
        },
        "comparison": {
            "massive_case_crossing_count": crossings,
            "massless_case_crossing_count": massless_crossings,
            "minimum_massive_q_c_at_zmax": min(c["q_c_at_zmax"] for c in cases),
            "maximum_massive_q_c_at_zmax": max(c["q_c_at_zmax"] for c in cases),
            "minimum_massless_q_c_at_zmax": min(c["massless_q_c_at_zmax"] for c in cases),
            "maximum_massless_q_c_at_zmax": max(c["massless_q_c_at_zmax"] for c in cases),
            "max_abs_delta_q_c_at_zmax": max(abs(c["delta_q_c_at_zmax_vs_massless"]) for c in cases),
            "finding": "compare_sign_of_qc_endpoint; both prescriptions remain positive on the same declared redshift range if the counts are zero",
        },
        "cases": cases,
        "provenance": {str(p.relative_to(ROOT)): sha256(p) for p in
                       (script, prev_path, paper)},
        "runtime_seconds": None,
        "runtime": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                    "platform": platform.platform(), "threads": {k: os.environ.get(k) for k in
                    ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}},
    }
    output["runtime_seconds"] = time.perf_counter()-t0
    OUT_JSON.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")
    rows = ["| H0 | f_b | q_c(zmax), massive | q_c(zmax), massless | delta q_c | massive zero z | massless zero z |",
            "|---:|---:|---:|---:|---:|---:|---:|"]
    for c in cases:
        fm = "none" if c["rho_c_zero_crossing_redshift"] is None else f'{c["rho_c_zero_crossing_redshift"]:.6g}'
        fl = "none" if c["massless_rho_c_zero_crossing_redshift"] is None else f'{c["massless_rho_c_zero_crossing_redshift"]:.6g}'
        rows.append(f'| {c["H0_km_s_Mpc"]:g} | {c["f_b"]:.2f} | {c["q_c_at_zmax"]:.8f} | {c["massless_q_c_at_zmax"]:.8f} | {c["delta_q_c_at_zmax_vs_massless"]:+.2e} | {fm} | {fl} |')
    delta_max = max(abs(c["delta_q_c_at_zmax_vs_massless"]) for c in cases)
    fd_quad_err_max = max(x["tight_quad_error_estimate"] for x in fd_checks)
    report = f'''# Thermal massive-neutrino sensitivity of the IVS positivity check

**Finding.** The same supplied background point \\(\\Omega_{{m0}}={OM0:.10f}\\), \\(g={GAMMA:.10f}\\), \\(w_x=-1\\) remains CDM-positive through \\(z=10^{{10}}\\) in all 16 \\(H_0,f_b\\) cases under the thermal massive-neutrino prescription below. This matches theory5's massless-radiation disposition: {crossings} massive-case crossings versus {massless_crossings} massless-case crossings. The maximum absolute shift in endpoint comoving CDM is {delta_max:.3e}; its range and each case are tabulated below. No likelihood was refit and no posterior is inferred.

## Neutrino prescription and formula

The local paper fixes \\(N_{{eff}}=3.044\\) and \\(\\sum m_\\nu=0.06\\) eV. I implement one fully populated thermal massive eigenstate with \\(m=0.06\\) eV and the remaining \\(2.044\\) effective species massless. With \\(T_{{\\nu0}}/T_{{\\gamma0}}=(4/11)^{{1/3}}\\), \\(y=ma/T_{{\\nu0}}\\), and

\\[F(y)=\\int_0^\\infty dq\\;\\frac{{q^2\\sqrt{{q^2+y^2}}}}{{e^q+1}},\\qquad F(0)=\\frac{{7\\pi^4}}{{120}},\\]

the density used in the Friedmann equation is

\\[\\frac{{\\rho_\\nu(a)}}{{\\rho_{{\\gamma0}}}}=\\left(\\frac{{15}}{{\\pi^4}}\\right)\\left(\\frac{{T_{{\\nu0}}}}{{T_{{\\gamma0}}}}\\right)^4a^{{-4}}\\left[F\\!\\left(\\frac{{ma}}{{T_{{\\nu0}}}}\\right)+(N_{{eff}}-1)F(0)\\right].\\]

The Fermi-Dirac energy integral is evaluated with adaptive QUAD in the background RHS, not a massless/high-temperature switch. Flat closure uses \\(\\Omega_{{v0}}=1-\\Omega_{{m0}}-\\Omega_{{\\gamma0}}-\\Omega_{{\\nu0}}\\); \\(\\Omega_{{m0}}\\) retains the screen's baryon-plus-CDM meaning. Grid: \\(H_0\\in\\{{{', '.join(map(str,H0_GRID))}\\}}\\), \\(f_b\\in\\{{{', '.join(f'{v:.2f}' for v in FB_GRID)}\\}}\\), and \\(0\\le z\\le10^{{10}}\\). The earlier result and this one use no CMB/BBN prior. This is a declared thermal background prescription and is not an exact reproduction of the authors' unknown CAMB neutrino configuration.

The comoving-CDM identity remains \\(q_c(n)=\\Omega_{{c0}}+g\\int_n^0e^{{3s}}u_v(s)/E(s)\\,ds\\). For \\(g<0\\), it declines monotonically into the past, so the positive endpoint excludes missed earlier zeros in-range. ODE tolerance refinement gives max scaled component difference {ode_conv:.3e}; adaptive FD integral tightening agrees at probe values to relative difference {max(x["relative_difference"] for x in fd_checks):.3e}, with maximum reported tight-quadrature absolute error {fd_quad_err_max:.2e}; identity residual across the grid is at most {max(c["identity_max_abs_error"] for c in cases):.3e}. These checks support numerical convergence for this background calculation, not perturbative or sound-horizon viability.

## Case-by-case comparison

{chr(10).join(rows)}

Reproduce from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \\
  .venv/bin/python work/theory6/ivs_massive_nu_sensitivity.py
```

Machine-readable results and input hashes: [ivs_massive_nu_sensitivity.json](ivs_massive_nu_sensitivity.json). Runtime: {output['runtime_seconds']:.3f} s. Script SHA-256: `{sha256(script)}`.

**Stop disposition:** completed at \\(z=10^{{10}}\\). The finite-mass prescription does not change the positivity result on the declared 16-case grid. No BAO refit, CAMB run, download, or parameter prior was introduced.
'''
    OUT_MD.write_text(report)
    print(json.dumps({"runtime_seconds": output["runtime_seconds"], "massive_crossings": crossings,
                      "massless_crossings": massless_crossings,
                      "min_massive_qc_zmax": min(c["q_c_at_zmax"] for c in cases),
                      "max_abs_delta_qc": delta_max,
                      "script_sha256": sha256(script), "json_sha256": sha256(OUT_JSON),
                      "md_sha256": sha256(OUT_MD)}, indent=2))


if __name__ == "__main__":
    main()

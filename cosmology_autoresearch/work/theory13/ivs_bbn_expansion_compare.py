#!/usr/bin/env python3
"""Independent radiation-era H(z) check for the archived IVS screen point.

This computes no abundances. It compares the interacting-vacuum screen point
with flat Lambda-CDM at identical present-day component densities and closure,
using massless photons/neutrinos throughout. Two separate state choices
integrate the IVS background, and g=0 must recover the Lambda reference.
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
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
OUT_JSON = Path(__file__).with_suffix(".json")
OUT_MD = Path(__file__).with_suffix(".md")

# Declared BAO-screen point and sensitivity coordinates; not CMB/BBN priors.
OMEGA_M0 = 0.3869710077
G = -0.4666647299  # Gamma/H0; positive Gamma transfers CDM to vacuum.
H0_VALUES = (55.0, 65.0, 75.0, 85.0)
FB_VALUES = (0.10, 0.15, 0.20, 0.25)

# Radiation prescription from theory5: all neutrinos are approximated as
# massless at all epochs, with the same N_eff used for present closure.
TCMB_K = 2.7255
OMEGA_GAMMA_H2 = 2.4728e-5
N_EFF = 3.044
NU_OVER_GAMMA_PER_NEFF = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0)

Z_MAX = 1.0e10
Z_BBN = np.geomspace(1.0e8, 1.0e10, 401)
Z_D_PIVOT = 3.4e8 / TCMB_K - 1.0
Z_TARGETS = (1.0e8, 1.0e9, 3.0e9, 1.0e10, Z_D_PIVOT)
RTOL = 2.0e-12
ATOL = 2.0e-14
TIGHT_RTOL = 2.0e-13
TIGHT_ATOL = 2.0e-15
MAX_STEP = 0.04


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def radiation_today(h0_km_s_mpc: float) -> tuple[float, float, float]:
    h = h0_km_s_mpc / 100.0
    omega_gamma = OMEGA_GAMMA_H2 / h**2
    omega_nu = omega_gamma * NU_OVER_GAMMA_PER_NEFF * N_EFF
    return omega_gamma, omega_nu, omega_gamma + omega_nu


def component_parameters(h0_km_s_mpc: float, f_b: float) -> dict[str, float]:
    omega_gamma, omega_nu, omega_r = radiation_today(h0_km_s_mpc)
    omega_b = f_b * OMEGA_M0
    omega_c = (1.0 - f_b) * OMEGA_M0
    omega_v = 1.0 - OMEGA_M0 - omega_r
    if omega_v <= 0.0:
        raise ValueError("flat closure produced nonpositive present vacuum density")
    return {
        "Omega_b0": omega_b,
        "Omega_c0": omega_c,
        "Omega_gamma0": omega_gamma,
        "Omega_nu_massless0": omega_nu,
        "Omega_r0": omega_r,
        "Omega_v0": omega_v,
    }


def integrate_direct(
    h0_km_s_mpc: float, f_b: float, g: float = G,
    rtol: float = RTOL, atol: float = ATOL,
):
    """Integrate q_c=a^3 rho_c/rhocrit0 and v=rho_v/rhocrit0 in u=ln(1+z)."""
    p = component_parameters(h0_km_s_mpc, f_b)
    qc0, v0 = p["Omega_c0"], p["Omega_v0"]

    def rhs(u: float, y: np.ndarray) -> tuple[float, float]:
        qc, v = y
        e3, e4 = math.exp(3.0 * u), math.exp(4.0 * u)
        e2 = (p["Omega_b0"] + qc) * e3 + p["Omega_r0"] * e4 + v
        if not math.isfinite(e2) or e2 <= 0.0:
            raise FloatingPointError(f"E^2 <= 0 at u={u:g}")
        e = math.sqrt(e2)
        a3 = math.exp(-3.0 * u)
        # Source-paper convention: dot(rho_v)=Gamma*rho_v and
        # dot(rho_c)+3H*rho_c=-Gamma*rho_v; u increases to the past.
        return g * a3 * v / e, -g * v / e

    umax = math.log1p(Z_MAX)
    sol = solve_ivp(
        rhs, (0.0, umax), (qc0, v0), method="DOP853", dense_output=True,
        rtol=rtol, atol=atol, max_step=MAX_STEP,
    )
    if not sol.success or sol.sol is None:
        raise RuntimeError(sol.message)
    return sol, p


def integrate_dark_total(
    h0_km_s_mpc: float, f_b: float, g: float = G,
    rtol: float = RTOL, atol: float = ATOL,
):
    """Independent state choice: q_d=a^3(rho_c+rho_v) and v=rho_v.

    The summed continuity equation gives dq_d/du=-3*a^3*v, independent of
    Gamma. H^2 follows from rho_c+rho_v=q_d*a^-3, baryons, and radiation.
    """
    p = component_parameters(h0_km_s_mpc, f_b)
    qd0 = p["Omega_c0"] + p["Omega_v0"]
    v0 = p["Omega_v0"]

    def rhs(u: float, y: np.ndarray) -> tuple[float, float]:
        qd, v = y
        e3, e4 = math.exp(3.0 * u), math.exp(4.0 * u)
        e2 = (p["Omega_b0"] + qd) * e3 + p["Omega_r0"] * e4
        if not math.isfinite(e2) or e2 <= 0.0:
            raise FloatingPointError(f"E^2 <= 0 at u={u:g}")
        e = math.sqrt(e2)
        a3 = math.exp(-3.0 * u)
        return -3.0 * a3 * v, -g * v / e

    umax = math.log1p(Z_MAX)
    sol = solve_ivp(
        rhs, (0.0, umax), (qd0, v0), method="DOP853", dense_output=True,
        rtol=rtol, atol=atol, max_step=MAX_STEP,
    )
    if not sol.success or sol.sol is None:
        raise RuntimeError(sol.message)
    return sol, p


def evaluate_case(h0: float, f_b: float, rtol: float = RTOL, atol: float = ATOL) -> dict:
    direct, p = integrate_direct(h0, f_b, G, rtol, atol)
    total, p_total = integrate_dark_total(h0, f_b, G, rtol, atol)
    if p != p_total:
        raise AssertionError("present-day closure differs between formulations")

    z_eval = np.unique(np.concatenate((np.asarray(Z_TARGETS), Z_BBN)))
    u_eval = np.log1p(z_eval)
    qc, v = direct.sol(u_eval)
    qd, v_total = total.sol(u_eval)
    e3 = np.exp(3.0 * u_eval)
    e4 = np.exp(4.0 * u_eval)
    a3 = np.exp(-3.0 * u_eval)

    omega_c = qc * e3
    e2_direct = (p["Omega_b0"] + qc) * e3 + p["Omega_r0"] * e4 + v
    e2_total = (p["Omega_b0"] + qd) * e3 + p["Omega_r0"] * e4
    e2_lambda = OMEGA_M0 * e3 + p["Omega_r0"] * e4 + p["Omega_v0"]

    # Form the small differences analytically rather than subtracting two
    # radiation-dominated E^2 values of order 1e36 near z=1e10.
    delta_e2_direct = (qc - p["Omega_c0"]) * e3 + (v - p["Omega_v0"])
    delta_e2_total = (qd - p["Omega_c0"]) * e3 - p["Omega_v0"]
    rel_e2_direct = delta_e2_direct / e2_lambda
    rel_e2_total = delta_e2_total / e2_lambda
    if np.any(1.0 + rel_e2_direct <= 0) or np.any(1.0 + rel_e2_total <= 0):
        raise FloatingPointError("nonpositive IVS/Lambda H^2 ratio")
    delta_h_direct = np.expm1(0.5 * np.log1p(rel_e2_direct))
    delta_h_total = np.expm1(0.5 * np.log1p(rel_e2_total))
    delta_h_component = np.sqrt(e2_direct / e2_lambda) - 1.0

    target_rows = []
    for z in Z_TARGETS:
        ix = int(np.flatnonzero(z_eval == z)[0])
        target_rows.append({
            "z": float(z),
            "delta_H_over_H_direct": float(delta_h_direct[ix]),
            "delta_H_over_H_dark_total": float(delta_h_total[ix]),
            "delta_H_over_H_component_ratio_check": float(delta_h_component[ix]),
            "Omega_v_over_Omega_r_at_z": float(v[ix] / (p["Omega_r0"] * e4[ix])),
            "Omega_m_baryon_plus_cdm_over_Omega_r_at_z": float(
                (p["Omega_b0"] + omega_c[ix]) / (p["Omega_r0"] * math.exp(u_eval[ix]))
            ),
            "rho_c_comoving_q": float(qc[ix]),
            "rho_v_normalized_today_rhocrit": float(v[ix]),
        })

    bbn_rows = []
    bbn_idx = np.searchsorted(z_eval, Z_BBN)
    for j, z in enumerate(Z_BBN):
        bbn_rows.append({
            "z": float(z),
            "delta_H_over_H": float(delta_h_direct[bbn_idx[j]]),
            "delta_H_over_H_dark_total": float(delta_h_total[bbn_idx[j]]),
        })

    return {
        "H0_km_s_Mpc": h0,
        "f_b": f_b,
        **p,
        "closure_sum_today": float(sum(p.values())),
        "target_redshifts": target_rows,
        "bbn_grid": bbn_rows,
        "checks": {
            "max_abs_direct_vs_dark_total_delta_H_over_H": float(
                np.max(np.abs(delta_h_direct - delta_h_total))
            ),
            "max_abs_component_E2_vs_analytic_delta_H_over_H": float(
                np.max(np.abs(delta_h_direct - delta_h_component))
            ),
            "max_abs_dark_total_v_vs_direct_v": float(np.max(np.abs(v_total - v))),
            "minimum_q_c": float(np.min(qc)),
            "minimum_v": float(np.min(v)),
            "minimum_E2_direct": float(np.min(e2_direct)),
            "minimum_E2_dark_total": float(np.min(e2_total)),
            "minimum_E2_Lambda": float(np.min(e2_lambda)),
        },
    }


def evaluate_null_case(h0: float, f_b: float) -> dict:
    direct, p = integrate_direct(h0, f_b, 0.0)
    total, _ = integrate_dark_total(h0, f_b, 0.0)
    z = np.asarray(Z_TARGETS)
    u = np.log1p(z)
    qc, v = direct.sol(u)
    qd, vtot = total.sol(u)
    e3, e4 = np.exp(3 * u), np.exp(4 * u)
    e2_lambda = OMEGA_M0 * e3 + p["Omega_r0"] * e4 + p["Omega_v0"]
    e2_direct = (p["Omega_b0"] + qc) * e3 + p["Omega_r0"] * e4 + v
    e2_total = (p["Omega_b0"] + qd) * e3 + p["Omega_r0"] * e4
    return {
        "H0_km_s_Mpc": h0,
        "f_b": f_b,
        "max_abs_direct_delta_H_over_H": float(np.max(np.abs(np.sqrt(e2_direct / e2_lambda) - 1.0))),
        "max_abs_total_delta_H_over_H": float(np.max(np.abs(np.sqrt(e2_total / e2_lambda) - 1.0))),
        "max_abs_direct_vs_total_E2_relative": float(
            np.max(np.abs(e2_direct - e2_total) / e2_lambda)
        ),
        "max_abs_v_change": float(np.max(np.abs(v - vtot))),
    }


def main() -> None:
    start = time.perf_counter()
    cases = [evaluate_case(h0, fb) for h0 in H0_VALUES for fb in FB_VALUES]
    nulls = [evaluate_null_case(h0, fb) for h0 in H0_VALUES for fb in FB_VALUES]
    # Refine all 16 sensitivity cases; compare both independent formulations on
    # the four requested redshifts and the deuterium pivot.
    convergence_rows = []
    for h0 in H0_VALUES:
        for fb in FB_VALUES:
            default = evaluate_case(h0, fb, RTOL, ATOL)
            tight = evaluate_case(h0, fb, TIGHT_RTOL, TIGHT_ATOL)
            default_targets = default["target_redshifts"]
            tight_targets = tight["target_redshifts"]
            convergence_rows.append({
                "H0_km_s_Mpc": h0,
                "f_b": fb,
                "max_abs_delta_H_over_H_change_default_to_tight": max(
                    abs(a["delta_H_over_H_direct"] - b["delta_H_over_H_direct"])
                    for a, b in zip(default_targets, tight_targets)
                ),
                "max_abs_direct_vs_total_default": default["checks"][
                    "max_abs_direct_vs_dark_total_delta_H_over_H"
                ],
                "max_abs_direct_vs_total_tight": tight["checks"][
                    "max_abs_direct_vs_dark_total_delta_H_over_H"
                ],
            })

    bbn_values = [row["delta_H_over_H"] for case in cases for row in case["bbn_grid"]]
    target_values = [
        row["delta_H_over_H_direct"] for case in cases for row in case["target_redshifts"]
    ]
    max_direct_total = max(
        case["checks"]["max_abs_direct_vs_dark_total_delta_H_over_H"] for case in cases
    )
    max_component_check = max(
        case["checks"]["max_abs_component_E2_vs_analytic_delta_H_over_H"] for case in cases
    )
    max_null = max(
        max(row["max_abs_direct_delta_H_over_H"], row["max_abs_total_delta_H_over_H"])
        for row in nulls
    )
    max_convergence = max(
        row["max_abs_delta_H_over_H_change_default_to_tight"] for row in convergence_rows
    )

    provenance_paths = (
        Path(__file__),
        ROOT / "work/theory5/ivs_early_positivity.py",
        ROOT / "work/theory5/ivs_early_positivity.json",
        ROOT / "context/papers/interacting_de_desi_dr2_2026.pdf",
        ROOT / "work/theory11/ivs_idecamb_patch_contract.md",
    )
    output = {
        "classification": "independent bounded homogeneous-background comparison; no abundance computation",
        "source_point": {
            "Omega_m0": OMEGA_M0,
            "g_Gamma_over_H0": G,
            "w_v": -1.0,
            "Q_convention": "dot(rho_v)=Q=Gamma*rho_v; dot(rho_c)+3H*rho_c=-Q",
            "sign_interpretation": (
                "Gamma<0 makes rho_v decrease forward/increase backward; backward comoving CDM "
                "q_c=a^3*rho_c decreases from its present value"
            ),
        },
        "units_and_normalization": {
            "H0": "km s^-1 Mpc^-1; cancels in H_IVS/H_Lambda",
            "g": "dimensionless Gamma/H0",
            "density_variables": "dimensionless fractions normalized to the same present critical density within each case",
            "independent_variable": "u=ln(1+z), increasing toward the past",
            "omega_gamma_h2": OMEGA_GAMMA_H2,
            "T_CMB_K": TCMB_K,
            "N_eff": N_EFF,
        },
        "assumptions": {
            "geometry": "flat FLRW; same present-day Omega_b, Omega_c, Omega_v closure, photons and neutrinos in IVS and Lambda",
            "neutrinos": "all neutrinos treated as massless radiation at all redshifts",
            "radiation_density": "Omega_gamma h^2=2.4728e-5; Omega_nu/Omega_gamma=(7/8)(4/11)^(4/3) N_eff",
            "H0_values_km_s_Mpc": list(H0_VALUES),
            "f_b_values": list(FB_VALUES),
            "bbn_grid": "401 logarithmic points in 1e8 <= z <= 1e10, spanning the requested BBN comparison bracket",
            "target_redshifts": list(map(float, Z_TARGETS)),
            "deuterium_pivot_redshift": Z_D_PIVOT,
            "temperature_redshift_note": (
                "Redshift comparison is direct; converting T_CMB to photon temperature at high z "
                "would require e+e- entropy-transfer corrections, not modeled here."
            ),
        },
        "equations": {
            "direct_scaled_state": [
                "q_c=a^3*rho_c/rho_crit0; v=rho_v/rho_crit0; u=ln(1+z)",
                "dq_c/du=g*a^3*v/E",
                "dv/du=-g*v/E",
                "E^2=(Omega_b0+q_c)*(1+z)^3+Omega_r0*(1+z)^4+v",
            ],
            "alternate_dark_total_state": [
                "q_d=a^3*(rho_c+rho_v)/rho_crit0",
                "dq_d/du=-3*a^3*v",
                "dv/du=-g*v/E",
                "E^2=(Omega_b0+q_d)*(1+z)^3+Omega_r0*(1+z)^4",
                "q_c=q_d-a^3*v",
            ],
            "Lambda_reference": "E_Lambda^2=Omega_m0*(1+z)^3+Omega_r0*(1+z)^4+Omega_v0",
            "stable_difference": "delta(E^2)=[q_c-Omega_c0]*(1+z)^3+(v-Omega_v0)",
        },
        "numerics": {
            "solver": "SciPy DOP853, float64, one CPU thread, dense output",
            "default_rtol": RTOL,
            "default_atol": ATOL,
            "tight_rtol": TIGHT_RTOL,
            "tight_atol": TIGHT_ATOL,
            "max_step_in_log1pz": MAX_STEP,
            "n_cases": len(cases),
            "n_dense_bbn_points_per_case": int(len(Z_BBN)),
            "max_abs_direct_vs_dark_total_delta_H_over_H": max_direct_total,
            "max_abs_component_ratio_vs_stable_delta_H_formula": max_component_check,
            "g_zero_null_max_abs_delta_H_over_H": max_null,
            "max_abs_delta_H_over_H_change_default_to_tight": max_convergence,
        },
        "summary": {
            "delta_H_over_H_min_targets_all_cases": min(target_values),
            "delta_H_over_H_max_targets_all_cases": max(target_values),
            "delta_H_over_H_min_dense_BBN_all_cases": min(bbn_values),
            "delta_H_over_H_max_dense_BBN_all_cases": max(bbn_values),
            "max_abs_delta_H_over_H_dense_BBN": max(map(abs, bbn_values)),
            "max_abs_delta_H_ppm_dense_BBN": 1e6 * max(map(abs, bbn_values)),
            "max_abs_Omega_v_over_Omega_r_at_targets": max(
                row["Omega_v_over_Omega_r_at_z"]
                for case in cases for row in case["target_redshifts"]
            ),
            "min_Omega_v_over_Omega_r_at_targets": min(
                row["Omega_v_over_Omega_r_at_z"]
                for case in cases for row in case["target_redshifts"]
            ),
        },
        "cases": cases,
        "g_zero_null_cases": nulls,
        "tolerance_convergence_cases": convergence_rows,
        "provenance_sha256": {
            str(p.relative_to(ROOT)): sha256(p) for p in provenance_paths
        },
        "environment": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "thread_env": {
                k: os.environ.get(k)
                for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
            },
        },
        "runtime_seconds": None,
    }
    output["runtime_seconds"] = time.perf_counter() - start
    OUT_JSON.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    target_rows = []
    for z in Z_TARGETS:
        vals = [
            next(row["delta_H_over_H_direct"] for row in case["target_redshifts"] if row["z"] == z)
            for case in cases
        ]
        target_rows.append(
            f"| {z:.6g} | {min(vals):+.7e} | {max(vals):+.7e} | "
            f"{1e6 * max(map(abs, vals)):.5f} |"
        )
    null_ok = max_null < 5e-13
    convergence_ok = max_convergence < 5e-12
    report = f"""# IVS expansion-rate effect during standard BBN

**Finding.** For the archived IVS screen point, the largest absolute expansion-rate change relative to flat Lambda-CDM with identical present-day component densities is {1e6 * max(map(abs, bbn_values)):.4f} ppm over 1e8 <= z <= 1e10. The shift is negative for Gamma/H0={G}: H_IVS is lower than H_Lambda. It is far below a 1% BBN-era expansion change, so this background point does not invalidate use of a standard-SBBN baryon prior on expansion-rate grounds. This is a homogeneous expansion comparison only; no abundance response or likelihood was computed.

## Setup and sign

This independently reimplements the theory5 screen point Omega_m0={OMEGA_M0:.10f}, g=Gamma/H0={G:.10f}, w_v=-1 over all 16 cases H0={list(H0_VALUES)} km s^-1 Mpc^-1 and f_b={list(FB_VALUES)}. For each case, IVS and Lambda-CDM share Omega_b0, Omega_c0, Omega_r0 and flat-closure Omega_v0. Photons use Omega_gamma h^2={OMEGA_GAMMA_H2}; neutrinos use N_eff={N_EFF} and the massless approximation at all redshifts. This extends the source paper's nonzero-mass setup with the explicit approximation requested here.

The source sign convention is dot(rho_v)=Gamma*rho_v and dot(rho_c)+3H*rho_c=-Gamma*rho_v. With u=ln(1+z), q_c=a^3 rho_c/rho_crit, v=rho_v/rho_crit, and E=H/H0, the direct equations are

    dq_c/du = g*a^3*v/E
    dv/du   = -g*v/E
    E^2     = (Omega_b0+q_c)*(1+z)^3 + Omega_r0*(1+z)^4 + v.

Thus negative g makes v increase backward and q_c decrease backward from Omega_c0. The Lambda reference is E_Lambda^2=Omega_m0*(1+z)^3+Omega_r0*(1+z)^4+Omega_v0. The reported delta-H uses an analytic difference of E^2 terms to avoid subtracting two radiation-dominated totals directly.

An independent second integration uses q_d=a^3(rho_c+rho_v)/rho_crit and v, with dq_d/du=-3*a^3*v and the same dv/du. It reconstructs q_c=q_d-a^3*v and H from total dark-sector density. The g=0 solution is the exact Lambda null in both formulations.

## Results

Ranges below are across all 16 H0/f_b cases. Delta-H is (H_IVS/H_Lambda)-1.

| z | min delta-H/H | max delta-H/H | largest magnitude (ppm) |
|---:|---:|---:|---:|
{chr(10).join(target_rows)}

On the 401-point logarithmic BBN bracket, delta-H/H ranges from {min(bbn_values):+.7e} to {max(bbn_values):+.7e}; maximum |Omega_v/Omega_r| at requested target redshifts is {output["summary"]["max_abs_Omega_v_over_Omega_r_at_targets"]:.3e}. This shows the vacuum component itself is negligible relative to radiation at BBN. The rate difference is dominated by the altered early comoving-CDM normalization, still strongly radiation-suppressed.

Independent direct-vs-total agreement is {max_direct_total:.3e} in delta-H/H. The maximum component-ratio versus stable-difference check is {max_component_check:.3e}. Across all 16 cases, the g=0 null has max |delta-H/H|={max_null:.3e} ({'PASS' if null_ok else 'FAIL'}); default-to-tight tolerance refinement changes target delta-H/H by at most {max_convergence:.3e} ({'PASS' if convergence_ok else 'FAIL'}). The full 401-point-per-case results, all checks, runtime, environment, and input hashes are in the companion JSON.

## Scope and disposition

The comparison supports using the standard-SBBN prior for this one screened IVS point with respect to the homogeneous H(T) effect, because |delta-H/H| stays below {1e6*max(map(abs, bbn_values)):.4f} ppm over the declared BBN redshift bracket. It does not prove the prior valid for the full IVS posterior or different interaction points. It also does not account for electron-positron entropy transfer in the massless-radiation approximation, evolve neutrino masses, compute light-element abundances, or test perturbations. Reproducing PRIMAT with this modified expansion would be required to bound the actual abundance shift. No PRIMAT install, download, GPU, or shared-file edit was used.

Reproduce from the project root with:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \\
      .venv/bin/python scripts/run_bounded.py --seconds 300 -- \\
      .venv/bin/python work/theory13/ivs_bbn_expansion_compare.py

Runtime: {output["runtime_seconds"]:.3f} s. Script SHA-256: {sha256(Path(__file__))}.
"""
    OUT_MD.write_text(report, encoding="utf-8")
    print(json.dumps({
        "runtime_seconds": output["runtime_seconds"],
        "max_abs_delta_H_over_H_BBN": max(map(abs, bbn_values)),
        "delta_H_over_H_target_ranges": {
            str(z): [min(
                next(row["delta_H_over_H_direct"] for row in case["target_redshifts"] if row["z"] == z)
                for case in cases
            ), max(
                next(row["delta_H_over_H_direct"] for row in case["target_redshifts"] if row["z"] == z)
                for case in cases
            )]
            for z in Z_TARGETS
        },
        "max_abs_direct_vs_total": max_direct_total,
        "max_abs_component_check": max_component_check,
        "g0_null": max_null,
        "tolerance_refinement": max_convergence,
        "json": str(OUT_JSON.relative_to(ROOT)),
        "report": str(OUT_MD.relative_to(ROOT)),
        "script_sha256": sha256(Path(__file__)),
    }, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Small, self-contained asymptotic check for IVS radiation-era ICs.

This evaluates the leading radiation-era background series and the source-paper
CDM perturbation equation. It is deliberately not a Boltzmann solver or a full
Einstein-feedback calculation and imports no project model/production code.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import time


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PAPER_PDF = REPO / "context/papers/interacting_de_desi_dr2_2026.pdf"
PAPER_TEXT = REPO / "context/papers/interacting_de_desi_dr2_2026.txt"
SOLVER_PATH = REPO / "work/perturbation_audit/solver_path.md"
OUTPUT = HERE / "early_ic_asymptotic.json"

# Independently pinned benchmark constants (not claimed to be paper best-fit values).
H0_KM_S_MPC = 67.4
OMEGA_B0 = 0.049
OMEGA_C0 = 0.266
OMEGA_M0 = 0.315
TCMB_K = 2.7255
N_EFF = 3.044
SUM_MNU_EV_PAPER = 0.06
G_VALUES = (0.0, 0.46666, -0.46666)  # g = Gamma/H0

# SI constants for the independently evaluated radiation density.
MPC_M = 3.0856775814913673e22
C_M_S = 299792458.0
G_N_SI = 6.67430e-11
A_RAD_J_M3_K4 = 7.5657e-16
NEUTRINO_TO_PHOTON_FACTOR = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def log_slope(xs: list[float], ys: list[float]) -> float:
    """Least-squares slope of log(abs(y)) against log(a)."""
    xx = [math.log(x) for x in xs]
    yy = [math.log(abs(y)) for y in ys]
    xbar = sum(xx) / len(xx)
    ybar = sum(yy) / len(yy)
    return sum((x - xbar) * (y - ybar) for x, y in zip(xx, yy)) / sum(
        (x - xbar) ** 2 for x in xx
    )


def radiation_parameters() -> dict[str, float]:
    h0_s = H0_KM_S_MPC * 1000.0 / MPC_M
    rho_gamma_mass = A_RAD_J_M3_K4 * TCMB_K**4 / C_M_S**2
    rho_critical_mass = 3.0 * h0_s**2 / (8.0 * math.pi * G_N_SI)
    omega_gamma = rho_gamma_mass / rho_critical_mass
    omega_nu = omega_gamma * NEUTRINO_TO_PHOTON_FACTOR * N_EFF
    omega_r = omega_gamma + omega_nu
    omega_x0 = 1.0 - OMEGA_B0 - OMEGA_C0 - omega_r
    return {
        "H0_s_inverse": h0_s,
        "Omega_gamma0": omega_gamma,
        "Omega_nu0_relativistic_approx": omega_nu,
        "Omega_r0": omega_r,
        "Omega_x0_flat_closure": omega_x0,
        "R_nu": omega_nu / omega_r,
        "Hrad_over_H0": math.sqrt(omega_r),
    }


def series_audit(g: float, omega_r: float) -> dict[str, object]:
    """Check power laws using the leading radiation-dominated series.

    The test normalization X*=D=1 is arbitrary and used only for exponents;
    the powers and the C=0 limit do not depend on these positive amplitudes.
    """
    root_r = math.sqrt(omega_r)
    xstar = 1.0
    d_comoving = 1.0

    # In the ordinary D>0 branch, A=(Gamma/H)(rho_x/rho_c) ~ K a^5.
    a_grid = [1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2]
    u = g / (2.0 * root_r)
    ordinary_A: list[float] = []
    for a in a_grid:
        rho_x = xstar * (1.0 + u * a**2)
        rho_c = d_comoving * a**-3 - g * xstar / (5.0 * root_r) * a**2
        gamma_over_H = g / root_r * a**2
        ordinary_A.append(gamma_over_H * rho_x / rho_c)
    expected_coefficient = g * xstar / (root_r * d_comoving)
    ordinary_slope = None if g == 0.0 else log_slope(a_grid, ordinary_A)
    ordinary_sign_ok = all(
        (value == 0.0 if g == 0.0 else value * g > 0.0)
        for value in ordinary_A
    )

    result: dict[str, object] = {
        "gamma_over_H0": g,
        "ordinary_D_positive_branch": {
            "assumed_Xstar_over_rhocrit": xstar,
            "assumed_D_over_rhocrit": d_comoving,
            "rho_x_series": "X*[1 + g/(2 sqrt(Omega_r0))*a^2 + O(a^3)]",
            "rho_c_series": "D*a^-3 - g*X*/(5 sqrt(Omega_r0))*a^2 + ...",
            "A_leading": "g*X*/(sqrt(Omega_r0)*D)*a^5",
            "A_log_slope_fit": ordinary_slope,
            "A_sign_matches_Gamma": ordinary_sign_ok,
            "A_coefficient_for_unit_amplitudes": expected_coefficient,
        },
    }

    if g < 0.0:
        # D=0 is positive only for Gamma<0. Here A tends to -5.
        edge_grid = [1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2]
        rho_x_edge = [xstar * (1.0 + u * a**2) for a in edge_grid]
        rho_c_edge = [-g * xstar / (5.0 * root_r) * a**2 for a in edge_grid]
        edge_A = [
            (g / root_r * a**2) * rx / rc
            for a, rx, rc in zip(edge_grid, rho_x_edge, rho_c_edge)
        ]
        edge_delta = [-a**2 / 7.0 for a in edge_grid]  # h=a^2
        edge_residual = []
        adiabatic_residual = []
        for a, A, delta in zip(edge_grid, edge_A, edge_delta):
            h = a**2
            delta_n = 2.0 * delta
            edge_residual.append(delta_n + h - A * delta)
            delta_adiabatic = h / 3.0
            adiabatic_residual.append(2.0 * delta_adiabatic + h - A * delta_adiabatic)
        correction = [A + 5.0 for A in edge_A]
        result["D_zero_negative_Gamma_edge"] = {
            "rho_c_leading": "-g*X*/(5 sqrt(Omega_r0))*a^2 > 0",
            "beta_c": 2.0,
            "A_limit": -5.0,
            "A_plus_5_log_slope_fit": log_slope(edge_grid, correction),
            "regular_delta_for_h_a2": "-h/7",
            "Eq19_regular_mode_residual_over_h_at_A_limit": 0.0,
            "finite_series_regular_residual_over_h_log_slope": log_slope(
                edge_grid, [r / a**2 for r, a in zip(edge_residual, edge_grid)]
            ),
            "adiabatic_delta_from_equal_zeta": "+h/3",
            "equal_zeta_delta_Eq19_residual_over_h_at_A_limit": 10.0 / 3.0,
            "finite_series_equal_zeta_residual_over_h": max(
                abs(r / (a**2))
                for r, a in zip(adiabatic_residual, edge_grid)
            ),
            "zeta_c_minus_zeta_gamma_over_h": 5.0 / 21.0,
            "positive_rho_c_on_asymptotic_test_grid": all(
                rc > 0.0 for rc in rho_c_edge
            ),
        }
    elif g > 0.0:
        result["D_zero_positive_Gamma_edge"] = {
            "rho_c_leading": "-g*X*/(5 sqrt(Omega_r0))*a^2 < 0",
            "allowed_with_positive_rho_c": False,
        }

    return result


def perturbation_audit(omega_r: float) -> dict[str, object]:
    root_r = math.sqrt(omega_r)
    h_amplitude = 1.0
    a_grid = [1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3]
    positive_negative: dict[str, object] = {}

    for g in (0.46666, -0.46666):
        k = g / root_r  # unit X*=D=1
        residuals = []
        A_values = []
        for a in a_grid:
            h = h_amplitude * a**2
            A = k * a**5
            # The cancellation leaves A5^2*h/14; evaluate that analytic residual
            # directly so the small leading term is not lost to float subtraction.
            residuals.append(A * A * h / 14.0)
            A_values.append(A)
        positive_negative[str(g)] = {
            "A_sign": "positive" if g > 0 else "negative",
            "A_log_slope_fit": log_slope(a_grid, A_values),
            "delta_c_ideal_pure_A5_fixed_h_series": "-h/2 - A5*h/14; A5=K5*a^5",
            "full_background_correction_note": "A=A5+K6*a^6+... generally; the K6 response is O(a^8) for h~a^2 and is not evaluated by this pure-A5 test",
            "ideal_pure_A5_next_residual_order": "O(A5^2*h)=O(a^12)",
            "max_abs_Eq19_residual_over_h": max(
                abs(r / (a**2)) for r, a in zip(residuals, a_grid)
            ),
            "max_abs_residual_over_A5_squared_h": max(
                abs(r / ((k * a**5) ** 2 * a**2))
                for r, a in zip(residuals, a_grid)
            ),
        }

    # In the D=0, g<0 branch A=-5 at leading order. The regular solution is
    # delta=-h/7; equal-zeta instead requires +h/3 and cannot solve Eq. (19).
    h = 1.0
    delta_eq19 = -h / 7.0
    delta_ad = h / 3.0
    beta_c_edge = 2.0
    beta_gamma = -4.0
    delta_gamma = -2.0 * h / 3.0
    zeta_difference = -(
        delta_eq19 / beta_c_edge - delta_gamma / beta_gamma
    )
    return {
        "n_equation_theta_c_zero": "d(delta_c)/d ln(a) = -S + A*delta_c; S=h_,ln(a)/2",
        "A_definition": "(Gamma/H)*(rho_x/rho_c) = (g/E)*(rho_x/rho_c)",
        "gamma_zero_limit": {
            "A": 0.0,
            "regular_radiation_delta_c_over_h": -0.5,
            "Eq19_residual_over_h": 0.0,
            "interpretation": "recovers the standard CDM-comoving radiation-era limit",
        },
        "theta_equation_check": {
            "n_equation": "d(theta_c)/d ln(a)=-theta_c",
            "analytic_solution": "a*theta_c=constant",
            "theta_c_zero_invariant": True,
            "nonzero_solution_regular_as_a_to_zero": False,
        },
        "ordinary_branch_pure_A5_fixed_radiation_tests": positive_negative,
        "edge_branch_C_zero": {
            "S_over_h_for_h_proportional_a2": 1.0,
            "A_limit": -5.0,
            "delta_regular_over_h": delta_eq19 / h,
            "delta_equal_zeta_over_h": delta_ad / h,
            "equal_zeta_Eq19_residual_over_h": (
                2.0 * delta_ad + h - (-5.0) * delta_ad
            )
            / h,
            "zeta_c_minus_zeta_gamma_over_h": zeta_difference / h,
            "homogeneous_edge_solution": "proportional to a^-5; exclude as singular",
        },
    }


def main() -> None:
    started = time.perf_counter()
    radiation = radiation_parameters()
    omega_r = radiation["Omega_r0"]
    assert omega_r > 0.0 and radiation["Omega_x0_flat_closure"] > 0.0

    cases = [series_audit(g, omega_r) for g in G_VALUES]
    perturbations = perturbation_audit(omega_r)
    assert all(
        case["ordinary_D_positive_branch"]["A_sign_matches_Gamma"]
        for case in cases
    )
    assert all(
        case["ordinary_D_positive_branch"]["A_log_slope_fit"] is None
        or abs(case["ordinary_D_positive_branch"]["A_log_slope_fit"] - 5.0) < 1.0e-3
        for case in cases
    )
    edge = cases[2]["D_zero_negative_Gamma_edge"]
    assert abs(edge["A_limit"] + 5.0) == 0.0
    assert edge["positive_rho_c_on_asymptotic_test_grid"] is True
    assert abs(edge["A_plus_5_log_slope_fit"] - 2.0) < 1.0e-10
    assert abs(edge["equal_zeta_delta_Eq19_residual_over_h_at_A_limit"] - 10.0 / 3.0) < 1.0e-15
    assert abs(edge["zeta_c_minus_zeta_gamma_over_h"] - 5.0 / 21.0) < 1.0e-15
    assert abs(
        perturbations["edge_branch_C_zero"]["delta_regular_over_h"] + 1.0 / 7.0
    ) < 1.0e-15
    assert abs(
        perturbations["edge_branch_C_zero"]["delta_equal_zeta_over_h"] - 1.0 / 3.0
    ) < 1.0e-15
    for case in perturbations["ordinary_branch_pure_A5_fixed_radiation_tests"].values():
        assert abs(case["max_abs_residual_over_A5_squared_h"] - 1.0 / 14.0) < 1.0e-14

    report = {
        "status": "PASS_ASYMPTOTIC_ALGEBRA; C_ZERO_EDGE_IS_NOT_CONVENTIONALLY_ADIABATIC",
        "execution": {
            "command": "env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 taskset -c 0 python3 -B scripts/run_bounded.py --seconds 600 --state work/theory10/no_campaign_deadline.json -- /usr/bin/time -f 'inner_elapsed=%e sec' python3 -B work/theory10/early_ic_asymptotic_check.py",
            "task_cap_seconds": 600,
            "cpu_affinity": "CPU 0; child inherits single-core affinity",
            "state_path": "work/theory10/no_campaign_deadline.json (intentionally absent; --seconds supplies cap)",
            "runtime_seconds": None,
        },
        "scope": "Leading radiation-era series and fixed-metric evaluation of source-paper IVS Eqs. (19)-(20); not a Boltzmann solver or full Einstein-feedback test.",
        "source_equations_and_signs": {
            "background": "dot(rho_x)=Gamma*rho_x; dot(rho_c)+3H*rho_c=-Gamma*rho_x; Gamma>0 means CDM -> vacuum",
            "separate_conservation": "baryons, photons/radiation, and neutrinos are separately conserved",
            "perturbations_conformal_time": "delta_c'=-(theta_c+h'/2)+a*Gamma*(rho_x/rho_c)*delta_c; theta_c'=-Hcal*theta_c",
            "frame": "synchronous gauge; regular CDM-comoving solution theta_c=0",
            "dimensionless_conversion": "g=Gamma/H0; a*Gamma/Hcal=Gamma/H=g/E",
        },
        "pinned_benchmark_context": {
            "H0_km_s_Mpc": H0_KM_S_MPC,
            "Omega_b0": OMEGA_B0,
            "Omega_c0": OMEGA_C0,
            "Omega_m0_b_plus_c": OMEGA_M0,
            "Tcmb_K": TCMB_K,
            "N_eff": N_EFF,
            "paper_sum_mnu_eV": SUM_MNU_EV_PAPER,
            "early_neutrinos": "ultrarelativistic; rho_nu/rho_gamma=(7/8)*(4/11)^(4/3)*N_eff; massless effective approximation only for a->0 asymptote",
            "flat_closure": "Omega_x0=1-Omega_b0-Omega_c0-Omega_gamma0-Omega_nu0 (benchmark context only)",
            **radiation,
        },
        "background_asymptotics": {
            "radiation_domination": "E~sqrt(Omega_r0)*a^-2; Hrad=H0*sqrt(Omega_r0)",
            "vacuum": "rho_x=X*[1+g/(2 sqrt(Omega_r0))*a^2+...] with X*>0",
            "comoving_cdm": "D=lim(a^3*rho_c); a^3*rho_c=D-g*X*/(5 sqrt(Omega_r0))*a^5+...",
            "allowed_positive_branches": {
                "D_positive": "rho_c~D*a^-3 for either sign of g; A~g*X*/(sqrt(Omega_r0)*D)*a^5",
                "D_zero_g_negative": "rho_c~-g*X*/(5 sqrt(Omega_r0))*a^2>0; A->-5",
                "D_zero_g_positive": "rho_c is negative at sufficiently early times; exclude",
                "D_zero_g_zero": "rho_c=0; CDM contrast equation is undefined",
                "D_negative": "negative rho_c as a->0; exclude",
            },
            "radiation_constraint_scalings": {
                "baryon_over_radiation": "O(a)",
                "vacuum_over_radiation": "O(a^4)",
                "CDM_over_radiation_D_positive": "O(a)",
                "CDM_over_radiation_D_zero_g_negative": "O(a^6)",
            },
        },
        "standard_radiation_mode_reference": {
            "source": "Ma & Bertschinger (1995), arXiv:astro-ph/9506072, Eq. (96)",
            "mode": "physical fastest-growing superhorizon mode after removing synchronous gauge modes; a proportional to tau",
            "leading_relations": "h=C(k*tau)^2; delta_gamma=-2h/3; delta_nu=delta_gamma; delta_b=delta_c=3*delta_gamma/4=-h/2; theta_c=0",
            "neutrino_shear": "retain free-streaming neutrino hierarchy; R_nu changes eta/shear coefficients",
            "paper_IVS_total_radiation_stress_closure": "the paper does not specify the IVS total-radiation stress closure; its zero-anisotropic-stress sentence is attached to the wx != -1 perturbation system",
        },
        "case_checks": cases,
        "perturbation_checks": perturbations,
        "interpretation": {
            "D_positive": "The standard adiabatic CDM/baryon/photon/neutrino relations survive at leading order because A5=K5*a^5. In the idealized pure-A5, fixed-h~a^2 equation, delta_c=-h/2-A5*h/14 and the next interaction iteration is O(A5^2*h)=O(a^12). In the full background A generally has a K6*a^6 term from the O(a) matter correction to E/radiation; its response is generally O(a^8), so the pure-A5 residual is not the full-background remainder.",
            "D_zero_g_negative": "The radiation metric remains the leading Einstein solution because rho_c/rho_r=O(a^6), but the displayed Eq. (19) selects regular delta_c=-h/7. Imposing the usual equal-zeta CDM-radiation criterion instead requires delta_c=+h/3 (beta_c=+2, beta_gamma=-4); the residual relative entropy is (zeta_c-zeta_gamma)/h=5/21. This is an incompatibility between the displayed CDM perturbation equation and that imposed IC criterion on the tuned edge, not a demonstrated inconsistency of a fully specified covariant vacuum-transfer completion or a claim about the posterior.",
            "D_positive_subleading_caveat": "Because rho_c is not exactly a^-3, exact equal-zeta relations acquire higher-order terms. Eq. (19)'s fixed-metric correction is not by itself a proof of a full all-orders adiabatic mode; the source paper does not specify a full scalar IC prescription or a perturbed-transfer completion.",
            "vacuum_adiabaticity_caveat": "The paper gives no independent delta_x evolution and states that IVS DE has no perturbative contribution. If this is implemented as delta_x=0 in the CDM frame, then for nonzero Gamma beta_x=d ln(rho_x)/d ln(a)=Gamma/H=O(a^2), and strict equal-zeta with radiation would instead require delta_x=beta_x*h/6=O(a^4). Thus the standard relations established here are leading relations for radiation/baryons/CDM; a strict all-species adiabatic vacuum prescription is not supplied by the paper closure.",
            "theta_and_gauge": "Eq. (20) gives theta_c proportional to a^-1; the regular CDM-comoving radiation mode has theta_c=0. Synchronous residual gauge modes must be removed; do not call a singular homogeneous perturbation solution a physical adiabatic mode.",
            "neutrino_mass_subtlety": "The paper pins sum m_nu=0.06 eV and N_eff=3.044, but does not specify mass splitting/distribution details in the cited assumptions. As a->0 all are relativistic; a finite-start solver still needs its massive-neutrino phase-space hierarchy/thermal prescription and must retain free-streaming shear.",
            "scope_warning": "This homogeneous-background and fixed-metric asymptotic audit is NOT a perturbative stability proof, CMB validation, or full-model/Boltzmann result.",
        },
        "hashes": {
            "script_sha256": sha256(Path(__file__).resolve()),
            "paper_pdf_sha256": sha256(PAPER_PDF),
            "paper_text_sha256": sha256(PAPER_TEXT),
            "solver_path_sha256": sha256(SOLVER_PATH),
        },
    }
    report["execution"]["runtime_seconds"] = time.perf_counter() - started
    OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(OUTPUT), "status": report["status"], "Omega_r0": omega_r}, sort_keys=True))


if __name__ == "__main__":
    main()

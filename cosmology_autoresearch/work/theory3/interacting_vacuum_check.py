#!/usr/bin/env python3
"""Synthetic equation and positivity checks for Q=Gamma rho_x, w_x=-1.

This is not a DESI fit. It integrates the flat-FLRW component conservation
equations and compares them with their total-matter reduction.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp


ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "context/papers/interacting_de_desi_dr2_2026.pdf"
OUT = Path(__file__).resolve().with_suffix(".json")
Z_MAX = 2.33
N_MIN = -math.log1p(Z_MAX)
RTOL = 2.0e-11
ATOL = 1.0e-13
SAMPLE_N = np.linspace(0.0, N_MIN, 401, dtype=np.float64)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def integrate_components(omega_m: float, baryon_fraction: float, gamma: float):
    """Integrate (Omega_b(a), Omega_c(a), Omega_x(a)) vs N=ln(a)."""
    if not (0.0 < omega_m < 1.0 and 0.0 <= baryon_fraction <= 1.0):
        raise ValueError("unphysical present-day density split")
    y0 = np.asarray([
        baryon_fraction * omega_m,
        (1.0 - baryon_fraction) * omega_m,
        1.0 - omega_m,
    ], dtype=np.float64)

    def rhs(n: float, y: np.ndarray) -> np.ndarray:
        omega_b, omega_c, omega_x = map(float, y)
        e2 = omega_b + omega_c + omega_x
        if not math.isfinite(e2) or e2 <= 0.0:
            raise FloatingPointError(f"E^2<=0 at N={n:g}")
        e = math.sqrt(e2)
        return np.asarray([
            -3.0 * omega_b,
            -3.0 * omega_c - gamma * omega_x / e,
            gamma * omega_x / e,
        ], dtype=np.float64)

    result = solve_ivp(
        rhs, (0.0, N_MIN), y0, method="DOP853", dense_output=True,
        rtol=RTOL, atol=ATOL, max_step=0.01,
    )
    if not result.success or result.sol is None:
        raise RuntimeError(result.message)
    return result


def expansion_squared(solution, n: float) -> float:
    omega_b, omega_c, omega_x = map(float, solution.sol(float(n)))
    return omega_b + omega_c + omega_x


def main() -> None:
    # Gamma=0 must recover flat LCDM independently of the baryon/CDM split.
    z_test = np.asarray([0.0, 0.10, 0.50, 1.0, Z_MAX], dtype=np.float64)
    n_test = -np.log1p(z_test)
    omega_m0 = 0.30
    lcdm_solution = integrate_components(omega_m0, 0.37, 0.0)
    lcdm_computed = np.asarray([
        expansion_squared(lcdm_solution, float(n)) for n in n_test
    ])
    lcdm_exact = omega_m0 * (1.0 + z_test) ** 3 + (1.0 - omega_m0)
    lcdm_error = float(np.max(np.abs(lcdm_computed - lcdm_exact)))

    # The total matter+vacuum background is independent of baryon/CDM split.
    split_checks = []
    split_values = (0.10, 0.30, 0.60, 0.90)
    for gamma in (-0.05, 0.05):
        predictions = []
        split_minima = []
        for fb in split_values:
            solution = integrate_components(0.30, fb, gamma)
            values = np.asarray([
                expansion_squared(solution, float(n)) for n in SAMPLE_N
            ])
            predictions.append(values)
            states = solution.sol(SAMPLE_N)
            split_minima.append({
                "baryon_fraction": fb,
                "minimum_Omega_c": float(np.min(states[1])),
                "minimum_Omega_x": float(np.min(states[2])),
                "minimum_E2": float(np.min(values)),
            })
        reference = predictions[0]
        split_checks.append({
            "gamma": gamma,
            "split_values": list(split_values),
            "max_abs_E2_difference_across_splits": float(max(
                np.max(np.abs(values - reference)) for values in predictions
            )),
            "component_positivity_at_splits": split_minima,
        })

    # Check the exact comoving-CDM integral and the sign of energy flow.
    transfer_checks = []
    transfer_integral_max_error = 0.0
    for gamma in (-0.05, 0.05):
        solution = integrate_components(0.30, 0.30, gamma)
        rows = []
        for z in (0.0, 0.5, 1.0, Z_MAX):
            n = -math.log1p(z)
            omega_b, omega_c, omega_x = map(float, solution.sol(n))
            if n == 0.0:
                integral = 0.0
            else:
                integral, _ = quad(
                    lambda np_: math.exp(3.0 * np_)
                    * float(solution.sol(np_)[2])
                    / math.sqrt(expansion_squared(solution, np_)),
                    n, 0.0, epsabs=1.0e-13, epsrel=1.0e-13, limit=100,
                )
            q_c = math.exp(3.0 * n) * omega_c
            q_expected = (1.0 - 0.30) * 0.30 + gamma * integral
            error = abs(q_c - q_expected)
            transfer_integral_max_error = max(transfer_integral_max_error, error)
            rows.append({
                "z": z,
                "comoving_Omega_c": q_c,
                "integral_identity_rhs": q_expected,
                "absolute_identity_error": error,
                "Omega_x": omega_x,
            })
        transfer_checks.append({
            "gamma": gamma,
            "direction": "CDM_to_vacuum" if gamma > 0.0 else "vacuum_to_CDM",
            "rows": rows,
        })

    # Sample a compact numerical screen box. The negative-Gamma half is not
    # automatically physical for every baryon/CDM split: reject paths with
    # negative rho_c. Gamma>=0 has an analytic backward-time positivity proof.
    omega_m_grid = np.linspace(0.05, 0.60, 9)
    gamma_grid = np.linspace(-0.20, 0.20, 9)
    fb_grid = np.linspace(0.0, 1.0, 7)
    grid_min_e2 = math.inf
    grid_min_c = math.inf
    grid_min_x = math.inf
    invalid_cdm_count = 0
    invalid_examples = []
    failed_integrations = []
    for omega_m in omega_m_grid:
        for gamma in gamma_grid:
            for fb in fb_grid:
                try:
                    solution = integrate_components(
                        float(omega_m), float(fb), float(gamma))
                    states = solution.sol(SAMPLE_N)
                    e2_values = states[0] + states[1] + states[2]
                    min_b = float(np.min(states[0]))
                    min_c = float(np.min(states[1]))
                    min_x = float(np.min(states[2]))
                    min_e2 = float(np.min(e2_values))
                    grid_min_c = min(grid_min_c, min_c)
                    grid_min_x = min(grid_min_x, min_x)
                    grid_min_e2 = min(grid_min_e2, min_e2)
                    if min_b < -1.0e-10 or min_c < -1.0e-10 or min_x < -1.0e-10 or min_e2 <= 0.0:
                        invalid_cdm_count += 1
                        if len(invalid_examples) < 8:
                            invalid_examples.append({
                                "Omega_m0": float(omega_m),
                                "gamma": float(gamma),
                                "baryon_fraction": float(fb),
                                "min_Omega_b": min_b,
                                "min_Omega_c": min_c,
                                "min_Omega_x": min_x,
                                "min_E2": min_e2,
                            })
                except (FloatingPointError, RuntimeError, ValueError) as error:
                    failed_integrations.append({
                        "Omega_m0": float(omega_m),
                        "gamma": float(gamma),
                        "baryon_fraction": float(fb),
                        "error": str(error),
                    })

    script_hash = sha256(Path(__file__).resolve())
    output = {
        "status": "independent synthetic equation check; no DESI data fit",
        "model": {
            "interaction": "Q=Gamma*rho_x",
            "sign": "Gamma>0 transfers energy from CDM to vacuum/DE",
            "dimensionless_coupling": "gamma=Gamma/H0",
            "background_variables": "N=ln(a), u_i=rho_i/rho_crit, E=H/H0",
            "equations": [
                "du_b/dN=-3*u_b",
                "du_c/dN=-3*u_c-gamma*u_x/E",
                "du_x/dN=gamma*u_x/E",
                "E^2=u_b+u_c+u_x (flat; radiation omitted)",
            ],
            "initial_values": [
                "u_b(0)=f_b*Omega_m0",
                "u_c(0)=(1-f_b)*Omega_m0",
                "u_x(0)=1-Omega_m0",
            ],
        },
        "numerical_domain": {
            "redshift": [0.0, Z_MAX],
            "Omega_m0_grid": [float(omega_m_grid[0]), float(omega_m_grid[-1]), len(omega_m_grid)],
            "gamma_grid": [float(gamma_grid[0]), float(gamma_grid[-1]), len(gamma_grid)],
            "baryon_fraction_grid": [float(fb_grid[0]), float(fb_grid[-1]), len(fb_grid)],
            "grid_case_count": int(len(omega_m_grid) * len(gamma_grid) * len(fb_grid)),
            "integration": "SciPy solve_ivp DOP853, backward in N from today; dense sample of 401 N values",
            "rtol": RTOL,
            "atol": ATOL,
            "analytic_positivity": (
                "u_x>0; u_b>=0. For gamma>=0, "
                "a^3*u_c(N)=Omega_c0+gamma*int_N^0 exp(3n)u_x(n)/E(n)dn>=0."
            ),
            "negative_gamma_viability": (
                "Must impose u_c>=0 over the whole integration interval; "
                "the rectangular gamma/f_b box is not uniformly physical."
            ),
            "grid_minimum_E2": grid_min_e2,
            "grid_minimum_Omega_c": grid_min_c,
            "grid_minimum_Omega_x": grid_min_x,
            "negative_component_grid_cases": invalid_cdm_count,
            "integration_failures": failed_integrations,
            "invalid_examples": invalid_examples,
        },
        "synthetic_checks": {
            "Gamma_zero_LCDM_recovery": {
                "Omega_m0": omega_m0,
                "baryon_fraction": 0.37,
                "z": z_test.tolist(),
                "max_abs_E2_error_vs_Omega_m_a_minus3_plus_lambda": lcdm_error,
            },
            "background_invariance_to_baryon_CDM_split": split_checks,
            "comoving_CDM_integral_sign_and_identity": {
                "identity": "a^3*u_c=Omega_c0+gamma*int_N^0 exp(3n)u_x/E dn",
                "maximum_absolute_identity_error": transfer_integral_max_error,
                "checks": transfer_checks,
            },
        },
        "paper": {
            "title": "Do DESI-DR2 BAO data imply a coupling of dark matter and dark energy?",
            "doi": "10.1103/6kdf-vzq4",
            "online_primary": "https://journals.aps.org/prd/accepted/10.1103/6kdf-vzq4",
            "local_pdf": "context/papers/interacting_de_desi_dr2_2026.pdf",
            "local_pdf_sha256": sha256(PAPER),
        },
        "execution": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "cpu_only": True,
            "blas_thread_environment": {
                key: os.environ.get(key)
                for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
            },
            "script_sha256": script_hash,
        },
    }
    OUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUT),
        "lcdm_max_abs_E2_error": lcdm_error,
        "split_checks": split_checks,
        "transfer_integral_max_abs_error": transfer_integral_max_error,
        "positivity_grid": output["numerical_domain"],
        "script_sha256": script_hash,
    }, indent=2))


if __name__ == "__main__":
    main()

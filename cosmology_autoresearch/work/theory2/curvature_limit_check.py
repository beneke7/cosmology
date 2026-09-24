#!/usr/bin/env python3
"""Tiny float64 audit for the curved-LCDM domain and S_k flat limit."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "curvature_limit_check.json"


def s_kernel(chi: float, omega_k: float) -> float:
    """Dimensionless transverse-curvature kernel, stable near omega_k=0."""
    y = omega_k * chi * chi
    if abs(y) < 1.0e-4:
        return chi * (1.0 + y / 6.0 + y**2 / 120.0
                      + y**3 / 5040.0 + y**4 / 362880.0)
    if y > 0.0:
        root = np.sqrt(y)
        return float(chi * np.sinh(root) / root)
    root = np.sqrt(-y)
    return float(chi * np.sin(root) / root)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    om_grid = np.linspace(0.05, 0.60, 12, dtype=np.float64)
    ok_grid = np.linspace(-0.20, 0.20, 17, dtype=np.float64)
    z_grid = np.unique(np.concatenate((
        np.linspace(0.0, 2.33, 1001, dtype=np.float64),
        np.asarray([5.0 / 3.0], dtype=np.float64),
    )))
    x = 1.0 + z_grid
    e2 = (1.0 + om_grid[:, None, None] * (x[None, None, :]**3 - 1.0)
          + ok_grid[None, :, None] * (x[None, None, :]**2 - 1.0))
    min_index = np.unravel_index(int(np.argmin(e2)), e2.shape)
    min_e2 = float(e2[min_index])
    min_omega_lambda = float(np.min(
        1.0 - om_grid[:, None] - ok_grid[None, :]))

    chi = 1.2
    h = 1.0e-6
    fd_curvature_derivative = (
        s_kernel(chi, h) - s_kernel(chi, -h)) / (2.0 * h)
    analytic_curvature_derivative = chi**3 / 6.0
    small_values = []
    for omega_k in (-1.0e-12, -1.0e-8, 0.0, 1.0e-8, 1.0e-12):
        value = s_kernel(chi, omega_k)
        small_values.append({
            "Omega_k": omega_k,
            "S_k_at_chi_1p2": value,
            "flat_limit_chi": chi,
            "fractional_difference_from_flat": (value - chi) / chi,
        })

    output = {
        "status": "independent contract audit; no fit",
        "dtype": "float64",
        "cpu_only": True,
        "domain": {
            "Omega_m": [0.05, 0.60],
            "Omega_k": [-0.20, 0.20],
            "z": [0.0, 2.33],
            "grid_shape_Omega_m_Omega_k_z": list(e2.shape),
            "minimum_E2_on_grid": min_e2,
            "minimum_location": {
                "Omega_m": float(om_grid[min_index[0]]),
                "Omega_k": float(ok_grid[min_index[1]]),
                "z": float(z_grid[min_index[2]]),
            },
            "analytic_continuous_minimum_E2": 73.0 / 108.0,
            "analytic_minimum_x": 8.0 / 3.0,
            "minimum_Omega_Lambda_on_parameter_rectangle": min_omega_lambda,
            "closed_branch_argument_upper_bound": float(
                np.sqrt(0.20) * 2.33 / np.sqrt(73.0 / 108.0)),
            "closed_branch_argument_is_below_pi": bool(
                np.sqrt(0.20) * 2.33 / np.sqrt(73.0 / 108.0) < np.pi),
        },
        "flat_limit": {
            "S_k_chi_1p2_at_Omega_k_zero": s_kernel(chi, 0.0),
            "analytic_dS_dOmega_k_at_zero": analytic_curvature_derivative,
            "central_finite_difference_h": h,
            "central_finite_difference_derivative": fd_curvature_derivative,
            "absolute_derivative_error": abs(
                fd_curvature_derivative - analytic_curvature_derivative),
            "near_zero_examples": small_values,
            "series": "chi * (1 + y/6 + y^2/120 + y^3/5040 + y^4/362880), y=Omega_k*chi^2, used for |y|<1e-4",
        },
        "script_sha256": digest(Path(__file__).resolve()),
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "machine": platform.machine(),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

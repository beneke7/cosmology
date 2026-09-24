"""Independent fixed-source check of the printed constant-K transform.

Uses a collocation BVP for the transformed linear radial equation and an
initial-value solve of the original nonlinear radial equation.  This does not
import or call scripts/audit_tg.py.  Run from the project root with
    .venv/bin/python work/theory/independent_check.py
All quantities below are dimensionless synthetic values.
"""
import json

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp


def check():
    k = 0.4
    r0, r_outer = 1.0e-4, 5.0
    source = lambda r: np.exp(-r * r)

    # Solve u'' + 2u'/r + k S(r)u = 0 by collocation.  The first boundary
    # condition is the regular-center expansion u'/u = -k S(0) r/3+O(r^3);
    # u(r_outer)=1 fixes the potential gauge phi(r_outer)=0.
    def transformed_rhs(r, y):
        return np.vstack((y[1], -2.0 * y[1] / r - k * source(r) * y[0]))

    def transformed_bc(ya, yb):
        center_regular = ya[1] + k * source(r0) * r0 * ya[0] / 3.0
        outer_dirichlet = yb[0] - 1.0
        return np.array([center_regular, outer_dirichlet])

    mesh = np.linspace(r0, r_outer, 500)
    guess = np.vstack((np.ones_like(mesh), np.zeros_like(mesh)))
    linear = solve_bvp(
        transformed_rhs, transformed_bc, mesh, guess, tol=2e-10,
        max_nodes=20000,
    )
    if not linear.success:
        raise RuntimeError(linear.message)

    # Independently integrate g=-dphi/dr, which obeys the original stationary
    # Riccati equation g' + 2g/r + S + k g^2 = 0.  The center series follows by
    # substitution into that equation for S=exp(-r^2).
    a1 = -1.0 / 3.0
    a3 = 1.0 / 5.0 - k / 45.0
    a5 = -1.0 / 14.0 + 2.0 * k / 105.0 - 2.0 * k * k / 945.0
    g0 = a1 * r0 + a3 * r0**3 + a5 * r0**5

    def nonlinear_rhs(r, y):
        g = y[0]
        return [-2.0 * g / r - source(r) - k * g * g]

    nonlinear = solve_ivp(
        nonlinear_rhs, (r0, r_outer), [g0], rtol=2e-11, atol=2e-13,
        dense_output=True,
    )
    if not nonlinear.success:
        raise RuntimeError(nonlinear.message)

    r = np.linspace(0.02, r_outer, 1200)
    u, du = linear.sol(r)
    # From u=exp(-k phi), g=-phi'=u'/(k u).
    g_from_u = du / (k * u)
    g_direct = nonlinear.sol(r)[0]
    # Reconstruct the stationary transformed residual from the BVP ODE.
    second = -2.0 * du / r - k * source(r) * u
    transformed_residual = second + 2.0 * du / r + k * source(r) * u

    return {
        "scope": "synthetic prescribed source only; independent of audit_tg.py; no galaxy fit",
        "method": "solve_bvp for transformed linear equation vs solve_ivp for original Riccati equation",
        "source": "S(r)=exp(-r^2), dimensionless",
        "K": k,
        "domain": [r0, r_outer],
        "boundary_conditions": {
            "center": "regular-center series to O(r), applied at r=1e-4",
            "outer": "u(5)=1, equivalent to phi(5)=0",
        },
        "bvp_success": bool(linear.success),
        "bvp_max_rms_residual": float(np.max(linear.rms_residuals)),
        "minimum_sampled_u": float(np.min(u)),
        "outer_phi": float(-np.log(linear.sol(r_outer)[0]) / k),
        "max_abs_g_difference": float(np.max(np.abs(g_from_u - g_direct))),
        "max_relative_g_difference": float(np.max(np.abs((g_from_u - g_direct) / g_direct))),
        "max_abs_transformed_stationary_residual": float(np.max(np.abs(transformed_residual))),
    }


if __name__ == "__main__":
    print(json.dumps(check(), indent=2))

#!/usr/bin/env python3
"""Conditional isolated-boundary TG test using the paper's spherical surrogate.

This is deliberately not called a physical NGC 3198 prediction: the Rotmod
force curves do not specify a unique 3-D source, and Eq. (47)'s spherical
inversion is negative over a robust radial band. Negative inferred density is
retained; it is never clipped.

Run from the project root:
  .venv/bin/python scripts/run_bounded.py --state RUN_STATE.json --seconds 600 -- \
    env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    .venv/bin/python experiments/tg_ngc3198_isolated_surrogate/solve_isolated.py
"""

from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import time
import zipfile

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import Akima1DInterpolator, CubicSpline, PchipInterpolator
from scipy.optimize import brentq
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "context" / "data" / "Rotmod_LTG.zip"
MEMBER = "NGC3198_rotmod.dat"
G_KPC_KMS2_PER_MSUN = 4.30091e-6
K_PAPER = 3.40e-5
UPSILON_PAPER = 0.762
K_SIGMA_PAPER = 0.04e-5
UPSILON_SIGMA_PAPER = 0.007
EPSILON = 1.0e-7  # kpc; regular-center series seed, not a fitted inner boundary


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_data() -> dict[str, np.ndarray]:
    with zipfile.ZipFile(DATA) as archive:
        raw = archive.read(MEMBER).decode("ascii")
    values = np.loadtxt(io.StringIO(raw), comments="#")
    r, vobs, errv, vgas, vdisk, vbul = values[:, :6].T
    return {
        "r": r, "vobs": vobs, "errv": errv,
        "vgas": vgas, "vdisk": vdisk, "vbul": vbul,
    }


def build_density(data, method: str, upsilon: float, gas_signed: bool = True):
    r = data["r"]
    disk = data["vdisk"]
    gas = data["vgas"]
    if method == "pchip":
        disk_i = PchipInterpolator(r, disk)
        gas_i = PchipInterpolator(r, gas)
    elif method == "natural_cubic":
        disk_i = CubicSpline(r, disk, bc_type="natural")
        gas_i = CubicSpline(r, gas, bc_type="natural")
    elif method == "akima":
        disk_i = Akima1DInterpolator(r, disk)
        gas_i = Akima1DInterpolator(r, gas)
    else:
        raise ValueError(method)

    r0 = float(r[0])
    rend = float(r[-1])
    disk0 = float(disk[0])
    gas0 = float(gas[0])
    central_density = (
        3.0 * (upsilon * disk0**2 + gas0 * abs(gas0))
        / (4.0 * np.pi * G_KPC_KMS2_PER_MSUN * r0**2)
    )

    def measured_density(radius):
        radius = np.asarray(radius)
        dv = disk_i(radius)
        ddv = disk_i.derivative()(radius)
        gv = gas_i(radius)
        dgv = gas_i.derivative()(radius)
        stellar_term = upsilon * (dv**2 + 2.0 * radius * dv * ddv)
        if gas_signed:
            gas_term = gv * np.abs(gv) + 2.0 * radius * np.abs(gv) * dgv
        else:
            gas_term = gv**2 + 2.0 * radius * gv * dgv
        return (stellar_term + gas_term) / (
            4.0 * np.pi * G_KPC_KMS2_PER_MSUN * radius**2
        )

    rho_rend = float(measured_density(rend))

    def density(radius, taper_kpc):
        radius = float(radius)
        if radius < r0:
            # Declare solid-body component velocities inside the first datum:
            # V(r) proportional to r implies a constant sphericalized rho.
            return central_density
        if radius <= rend:
            return float(measured_density(radius))
        fraction = max(1.0 - (radius - rend) / taper_kpc, 0.0)
        # A linear taper to a vacuum exterior. Keep a negative endpoint negative.
        return rho_rend * fraction

    return density, central_density, rho_rend, r0, rend


def integrate_segments(rhs, breaks, y0, max_step, rtol, atol, dense_output=False):
    solutions = []
    y = np.asarray(y0, dtype=float)
    for left, right in zip(breaks[:-1], breaks[1:]):
        if right <= left:
            continue
        solution = solve_ivp(
            rhs, (left, right), y, method="DOP853", rtol=rtol, atol=atol,
            max_step=max_step, dense_output=dense_output,
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        solutions.append(solution)
        y = solution.y[:, -1]
    return solutions


def evaluate_segments(solutions, x):
    x = np.asarray(x)
    out = np.empty((2, x.size), dtype=float)
    for index, value in enumerate(x):
        chosen = None
        for solution in solutions:
            if solution.t[0] - 1e-12 <= value <= solution.t[-1] + 1e-12:
                chosen = solution
                break
        if chosen is None:
            raise ValueError(f"requested radius {value} outside integrated domain")
        out[:, index] = chosen.sol(float(value)) if chosen.sol else chosen.y[:, -1]
    return out


def solve_one(data, method, taper_kpc, k, upsilon, max_step=0.04, gas_signed=True):
    r = data["r"]
    r0, rend = float(r[0]), float(r[-1])
    rout = rend + float(taper_kpc)
    density, rho0, rhoend, _, _ = build_density(data, method, upsilon, gas_signed)

    def source(radius):
        return 4.0 * np.pi * G_KPC_KMS2_PER_MSUN * density(radius, taper_kpc)

    s0 = source(0.5 * r0)
    # Regular-center expansions for the transformed linear solution and the
    # original Riccati acceleration, both seeded at the same tiny epsilon.
    u0 = 1.0
    du0 = -k * s0 * EPSILON * u0 / 3.0

    def linear_rhs(radius, state):
        value, derivative = state
        return [derivative, -2.0 * derivative / radius - k * source(radius) * value]

    breaks = [EPSILON, r0, rend, rout]
    u_solutions = integrate_segments(
        linear_rhs, breaks, [u0, du0], max_step, rtol=2e-11, atol=2e-13,
        dense_output=True,
    )
    raw_edge = evaluate_segments(u_solutions, np.array([rout]))[:, 0]
    edge_u, edge_du = map(float, raw_edge)
    robin_denominator = edge_u + rout * edge_du
    denominator_scale = abs(edge_u) + rout * abs(edge_du)
    if abs(robin_denominator) <= 1e-12 * max(denominator_scale, 1.0):
        return {
            "valid": False, "reason": "isolated normalization denominator is singular",
            "robin_denominator": robin_denominator,
        }
    amplitude = 1.0 / robin_denominator
    dense_r = np.linspace(EPSILON, rout, 6001)
    raw = evaluate_segments(u_solutions, dense_r)
    u_dense = amplitude * raw[0]
    du_dense = amplitude * raw[1]
    u_at_data, du_at_data = evaluate_segments(u_solutions, r)
    u_at_data *= amplitude
    du_at_data *= amplitude
    positive = bool(np.all(u_dense > 0.0))
    if positive:
        v2 = -r * du_at_data / (k * u_at_data)
        velocity = np.sqrt(np.maximum(v2, 0.0))
        curve_real = bool(np.all(v2 >= 0.0))
        chi2 = float(np.sum(((data["vobs"] - velocity) / data["errv"])**2)) if curve_real else None
    else:
        v2 = np.full_like(r, np.nan)
        velocity = np.full_like(r, np.nan)
        curve_real = False
        chi2 = None

    # Independently integrate the original Riccati equation and compare g with
    # g=u'/(K u); this is not the same linear ODE implementation.
    a1 = -s0 / 3.0
    a3 = -k * s0**2 / 45.0
    g0 = a1 * EPSILON + a3 * EPSILON**3

    def riccati_rhs(radius, state):
        g = float(state[0])
        return [-2.0 * g / radius - source(radius) - k * g * g]

    g_solutions = integrate_segments(
        riccati_rhs, breaks, [g0], max_step, rtol=2e-11, atol=2e-13,
        dense_output=True,
    )
    g_at_data = np.array([
        evaluate_segments(g_solutions, np.array([radius]))[0, 0] for radius in r
    ])
    g_from_u = du_at_data / (k * u_at_data) if positive else np.full_like(r, np.nan)
    riccati_delta_v = (
        np.sqrt(np.maximum(-r * g_at_data, 0.0)) - velocity
        if curve_real else np.full_like(r, np.nan)
    )
    robin_residual = edge_du * amplitude + (edge_u * amplitude - 1.0) / rout
    exterior_b = rout * (amplitude * edge_u - 1.0)
    return {
        "valid": bool(positive and curve_real and np.isfinite(chi2)),
        "reason": None if positive and curve_real else "nonpositive u or non-real circular speed",
        "K_s2_per_km2": float(k),
        "Upsilon_disk": float(upsilon),
        "interpolator": method,
        "gas_convention": "signed Vgas*abs(Vgas)" if gas_signed else "literal positive Vgas^2",
        "outer_taper_kpc": float(taper_kpc),
        "numerical_max_step_kpc": float(max_step),
        "source_density_centre_Msun_kpc3": float(rho0),
        "source_density_at_last_datum_Msun_kpc3": float(rhoend),
        "robin_normalization_denominator": float(robin_denominator),
        "central_u_normalization": float(amplitude),
        "minimum_u_on_domain": float(np.min(u_dense)),
        "maximum_u_on_domain": float(np.max(u_dense)),
        "outer_u": float(amplitude * edge_u),
        "outer_exterior_B_kpc": float(exterior_b),
        "robin_boundary_residual_per_kpc": float(robin_residual),
        "min_circular_speed_squared_km2_s2": float(np.min(v2)) if curve_real else None,
        "max_abs_Riccati_vs_transformed_velocity_km_s": (
            float(np.max(np.abs(riccati_delta_v))) if curve_real else None
        ),
        "rms_Riccati_vs_transformed_velocity_km_s": (
            float(np.sqrt(np.mean(riccati_delta_v**2))) if curve_real else None
        ),
        "chi2_vs_all_43_observed_points": chi2,
        "dof_fixed_parameters": int(r.size),
        "max_abs_fractional_Riccati_g_difference": (
            float(np.max(np.abs((g_from_u - g_at_data) / np.maximum(np.abs(g_at_data), 1e-30))))
            if curve_real else None
        ),
        "predicted_vc_km_s": velocity.tolist() if curve_real else None,
    }


def isolated_normalization_denominator(
    data, method, taper_kpc, k, upsilon, max_step=0.08, gas_signed=True
):
    """Return F(K)=w(R)+R*w'(R) for the regular unit-center basis w."""
    r = data["r"]
    r0, rend = float(r[0]), float(r[-1])
    rout = rend + float(taper_kpc)
    density, _, _, _, _ = build_density(data, method, upsilon, gas_signed)

    def source(radius):
        return 4.0 * np.pi * G_KPC_KMS2_PER_MSUN * density(radius, taper_kpc)

    s0 = source(0.5 * r0)

    def rhs(radius, state):
        value, derivative = state
        return [derivative, -2.0 * derivative / radius - k * source(radius) * value]

    breaks = [EPSILON, r0, rend, rout]
    solutions = integrate_segments(
        rhs, breaks, [1.0, -k * s0 * EPSILON / 3.0], max_step,
        rtol=2e-10, atol=2e-12, dense_output=True,
    )
    w_r, dw_r = evaluate_segments(solutions, np.array([rout]))[:, 0]
    return float(w_r + rout * dw_r)


def first_normalization_threshold(
    data, method, taper_kpc, upsilon, k_max=1.0e-4, gas_signed=True
):
    """Find the first isolated-normalization singularity on K>=0 by bracketing."""
    k_grid = np.linspace(0.0, k_max, 101)
    f_grid = np.array([
        isolated_normalization_denominator(
            data, method, taper_kpc, float(k), upsilon,
            max_step=0.08, gas_signed=gas_signed,
        )
        for k in k_grid
    ])
    crossings = np.flatnonzero(f_grid[:-1] * f_grid[1:] <= 0.0)
    if crossings.size == 0:
        return {
            "found": False,
            "K_scan_max_s2_per_km2": k_max,
            "F_at_scan_max": float(f_grid[-1]),
            "scan_points": int(k_grid.size),
        }
    i = int(crossings[0])
    left, right = float(k_grid[i]), float(k_grid[i + 1])
    root = brentq(
        lambda k: isolated_normalization_denominator(
            data, method, taper_kpc, float(k), upsilon,
            max_step=0.04, gas_signed=gas_signed,
        ),
        left, right, xtol=1e-12, rtol=1e-12, maxiter=100,
    )
    fcrit = isolated_normalization_denominator(
        data, method, taper_kpc, root, upsilon,
        max_step=0.02, gas_signed=gas_signed,
    )
    below = isolated_normalization_denominator(
        data, method, taper_kpc, root * (1.0 - 1e-5), upsilon,
        max_step=0.02, gas_signed=gas_signed,
    )
    return {
        "found": True,
        "first_K_critical_s2_per_km2": float(root),
        "K_critical_over_K_paper": float(root / K_PAPER),
        "fractional_K_margin_from_paper": float(1.0 - K_PAPER / root),
        "F_at_K_critical": float(fcrit),
        "F_at_0.99999_Kcritical": float(below),
        "bracket_s2_per_km2": [left, right],
        "scan_points": int(k_grid.size),
        "K_scan_max_s2_per_km2": k_max,
    }


def leading_radial_eigenvalue(data, method, taper_kpc, k, upsilon, intervals):
    """P1 finite-element max eigenvalue of Delta_radial+K*S with Robin match.

    Under q=r*w perturbations, regularity gives w(0)=0 and the homogeneous
    static-exterior Robin condition gives w'(R)=0. The returned eigenvalue is
    per D (units kpc^-2); the actual growth rate is D times this value.
    """
    r0, rend = float(data["r"][0]), float(data["r"][-1])
    rout = rend + float(taper_kpc)
    density, _, _, _, _ = build_density(data, method, upsilon, gas_signed=True)
    h = rout / intervals
    nodes = np.linspace(0.0, rout, intervals + 1)
    midpoints = 0.5 * (nodes[:-1] + nodes[1:])
    v_mid = np.array([
        k * 4.0 * np.pi * G_KPC_KMS2_PER_MSUN * density(float(x), taper_kpc)
        for x in midpoints
    ])
    n = intervals  # nodes 1..N retained; node zero is Dirichlet
    stiffness_diag = np.full(n, -2.0 / h)
    stiffness_diag[-1] = -1.0 / h  # natural Neumann at the outer edge
    mass_diag = np.full(n, 2.0 * h / 3.0)
    mass_diag[-1] = h / 3.0
    potential_diag = np.zeros(n)
    potential_off = np.zeros(n - 1)
    mass_off = np.full(n - 1, h / 6.0)

    # Element e joins original nodes e,e+1; remove the Dirichlet node at e=0.
    for e, value in enumerate(v_mid):
        if e >= 1:
            potential_diag[e - 1] += value * h / 3.0
        potential_diag[e] += value * h / 3.0
        if e >= 1:
            potential_off[e - 1] += value * h / 6.0

    stiffness_off = np.full(n - 1, 1.0 / h)
    operator = diags(
        [stiffness_off + potential_off, stiffness_diag + potential_diag,
         stiffness_off + potential_off],
        offsets=[-1, 0, 1], format="csr",
    )
    mass = diags([mass_off, mass_diag, mass_off], offsets=[-1, 0, 1], format="csr")
    value = eigsh(
        operator, k=1, M=mass, which="LA", return_eigenvectors=False,
        tol=2e-10, maxiter=max(1000, 20 * intervals),
    )[0]
    return float(value)


def uniform_sphere_eigenvalue_check():
    """Independent FE check against the analytic qR=1.3 uniform-sphere mode."""
    rmax = 1.0
    x = 1.3
    exact = x**2 - (0.5 * np.pi)**2
    rows = []
    for intervals in (100, 200, 400):
        h = rmax / intervals
        n = intervals
        # V=q^2=x^2/R^2 is constant; same P1 matrices as in the nonuniform case.
        value = x**2 / rmax**2
        diag = np.full(n, -2.0 / h + value * 2.0 * h / 3.0)
        diag[-1] = -1.0 / h + value * h / 3.0
        off = np.full(n - 1, 1.0 / h + value * h / 6.0)
        mass_diag = np.full(n, 2.0 * h / 3.0)
        mass_diag[-1] = h / 3.0
        mass_off = np.full(n - 1, h / 6.0)
        operator = diags([off, diag, off], [-1, 0, 1], format="csr")
        mass = diags([mass_off, mass_diag, mass_off], [-1, 0, 1], format="csr")
        numeric = float(eigsh(
            operator, k=1, M=mass, which="LA", return_eigenvectors=False,
            tol=1e-11, maxiter=20 * intervals,
        )[0])
        rows.append({
            "intervals": intervals,
            "lambda_numeric_per_D": numeric,
            "lambda_exact_per_D": exact,
            "absolute_error_kpc2": abs(numeric - exact),
        })
    return rows


def main() -> None:
    started = time.monotonic()
    data = load_data()
    reference = solve_one(
        data, "pchip", 10.0, K_PAPER, UPSILON_PAPER, max_step=0.04,
        gas_signed=True,
    )

    # Numerical resolution convergence at fixed physical source and boundary.
    convergence = []
    for step in (0.08, 0.04, 0.02):
        convergence.append(solve_one(
            data, "pchip", 10.0, K_PAPER, UPSILON_PAPER, max_step=step,
            gas_signed=True,
        ))
    velocity_curves = [np.asarray(row["predicted_vc_km_s"]) for row in convergence]
    convergence_record = {
        "step_sizes_kpc": [0.08, 0.04, 0.02],
        "max_abs_delta_v_km_s_step_0.08_vs_0.04": float(
            np.max(np.abs(velocity_curves[0] - velocity_curves[1]))
        ),
        "max_abs_delta_v_km_s_step_0.04_vs_0.02": float(
            np.max(np.abs(velocity_curves[1] - velocity_curves[2]))
        ),
        "u_positive_classification": [bool(row["valid"]) for row in convergence],
        "pass_0.1_km_s_and_classification": bool(
            np.max(np.abs(velocity_curves[0] - velocity_curves[1])) < 0.1
            and np.max(np.abs(velocity_curves[1] - velocity_curves[2])) < 0.1
            and len({row["valid"] for row in convergence}) == 1
        ),
    }

    # Declare before inspection: interpolation and compact-support taper variations.
    source_sensitivity = []
    for method in ("pchip", "natural_cubic", "akima"):
        for taper in (5.0, 10.0, 20.0):
            source_sensitivity.append(solve_one(
                data, method, taper, K_PAPER, UPSILON_PAPER, max_step=0.04,
                gas_signed=True,
            ))

    # Isolate the sign ambiguity in Eq. (47); negative gas rows are small and
    # confined to the inner data, but the convention is still explicitly tested.
    gas_sign_sensitivity = solve_one(
        data, "pchip", 10.0, K_PAPER, UPSILON_PAPER, max_step=0.04,
        gas_signed=False,
    )

    # Bounded local grid around the paper's reported 1-sigma K and M/L ranges.
    local_grid = []
    for k in np.linspace(K_PAPER - K_SIGMA_PAPER, K_PAPER + K_SIGMA_PAPER, 9):
        for upsilon in (UPSILON_PAPER - UPSILON_SIGMA_PAPER, UPSILON_PAPER,
                        UPSILON_PAPER + UPSILON_SIGMA_PAPER):
            row = solve_one(
                data, "pchip", 10.0, float(k), float(upsilon), max_step=0.04,
                gas_signed=True,
            )
            local_grid.append(row)
    valid_grid = [row for row in local_grid if row["valid"]]
    best_local = min(valid_grid, key=lambda row: row["chi2_vs_all_43_observed_points"]) if valid_grid else None

    # Warden-requested distance to the first isolated normalization/stability
    # singularity. Scan the declared source family without using observed velocities.
    critical_cases = []
    cases = [
        ("pchip", 5.0), ("pchip", 10.0), ("pchip", 20.0),
        ("natural_cubic", 10.0), ("akima", 10.0),
    ]
    for method, taper in cases:
        critical_cases.append({
            "interpolator": method,
            "outer_taper_kpc": taper,
            **first_normalization_threshold(data, method, taper, UPSILON_PAPER),
        })

    eigenvalue_rows = []
    # P1 FE grid refinement for the paper point and just below the baseline first
    # singularity; D is factored out because its scale is not specified here.
    baseline_critical = critical_cases[1].get("first_K_critical_s2_per_km2")
    if baseline_critical is not None:
        for label, k_eval in (
            ("paper_K", K_PAPER),
            ("0.99_Kcritical", 0.99 * baseline_critical),
        ):
            eigenvalue_rows.append({
                "case": label,
                "K_s2_per_km2": float(k_eval),
                "grid_values_per_D_kpc2": [
                    leading_radial_eigenvalue(
                        data, "pchip", 10.0, float(k_eval), UPSILON_PAPER, intervals
                    )
                    for intervals in (400, 800, 1600)
                ],
                "intervals": [400, 800, 1600],
            })

    record = {
        "status": "conditional_isolated_spherical_surrogate_test_complete",
        "interpretation": (
            "No physical NGC 3198 prediction: the Eq. (47) pseudo-spherical source is negative "
            "over a robust band and SPARC Rotmod does not determine a unique 3-D source, center, "
            "or exterior taper. This tests only a predeclared mathematical surrogate."
        ),
        "inputs": {
            "paper": "user-provided Pszota & Van (2024), Eq. (47), Table 1",
            "paper_sha256": sha256(ROOT / "context/user_provided/pszota_van_2024.pdf"),
            "data_archive": "context/data/Rotmod_LTG.zip",
            "data_archive_sha256": sha256(DATA),
            "member": MEMBER,
            "rows": int(data["r"].size),
            "radial_range_kpc": [float(data["r"][0]), float(data["r"][-1])],
        },
        "equations_and_units": {
            "stationary_transformed": "u'' + 2u'/r + 4*pi*G*K*rho*u = 0; u=exp(-K*phi)>0",
            "isolated_robin": "u'(R)=-(u(R)-1)/R, with vacuum u=1+B/r outside R",
            "original_riccati": "g' + 2g/r + 4*pi*G*rho + K*g^2 = 0; g=-phi'",
            "velocity_map": "vc^2=-r*g=-r*u'/(K*u)",
            "G": f"{G_KPC_KMS2_PER_MSUN} kpc (km/s)^2 / Msun",
            "K": "s^2/km^2",
            "rho": "Msun/kpc^3",
            "dimensionless_operator_coefficient": "4*pi*G*K*rho has units kpc^-2",
        },
        "declared_source_and_boundary": {
            "source_recipe": "paper Eq. (47) pseudo-spherical inversion; interpolate Vdisk,Vgas then differentiate r*V_component^2",
            "gas_convention": "SPARC signed-force convention Vgas*abs(Vgas); literal positive square tested separately",
            "Upsilon_disk": UPSILON_PAPER,
            "interpolator_baseline": "PCHIP; natural cubic and Akima sensitivity runs",
            "center": "for r<0.32 kpc, linearly extrapolate each component velocity as V(r)=V(0.32)*r/0.32, giving constant central rho; no central data supplied",
            "exterior": "linearly taper inferred rho from its value at 44.08 kpc to zero across declared 5,10,20 kpc widths, then exact vacuum; inferred negative density is retained",
            "outer_boundary": "exact vacuum matching at R=44.08 kpc+taper; no observed outer velocity used as a boundary",
        },
        "reported_parameter_point": reference,
        "numerical_convergence": convergence_record,
        "uniform_sphere_eigenvalue_method_check": uniform_sphere_eigenvalue_check(),
        "first_isolated_normalization_thresholds": critical_cases,
        "leading_radial_growth_eigenvalue": {
            "method": "P1 finite elements for q=r*delta_u, Dirichlet q(0)=0, Neumann q'(R)=0",
            "finite_domain_quasistatic_exterior_only": True,
            "D_factor": "The physical growth rate is D times this per-D eigenvalue; D is not specified.",
            "grid_refinement": eigenvalue_rows,
            "interpretation": "negative is linearly decaying and positive growing only for the declared finite Robin realization; full-space diffusion has continuum spectrum reaching zero",
        },
        "resolution_records": convergence,
        "source_interpolation_taper_sensitivity": source_sensitivity,
        "literal_vs_signed_gas": {
            "literal_eq47_square": gas_sign_sensitivity,
            "signed_eq53_sparc_convention": reference,
            "max_delta_predicted_vc_km_s": (
                float(np.max(np.abs(
                    np.asarray(reference["predicted_vc_km_s"])
                    - np.asarray(gas_sign_sensitivity["predicted_vc_km_s"])
                )))
                if reference["valid"] and gas_sign_sensitivity["valid"] else None
            ),
        },
        "local_parameter_grid": {
            "K_range_s2_per_km2": [K_PAPER - K_SIGMA_PAPER, K_PAPER + K_SIGMA_PAPER],
            "Upsilon_values": [
                UPSILON_PAPER - UPSILON_SIGMA_PAPER,
                UPSILON_PAPER,
                UPSILON_PAPER + UPSILON_SIGMA_PAPER,
            ],
            "points": len(local_grid),
            "valid_points": len(valid_grid),
            "best_local_point": best_local,
            "all_points": local_grid,
        },
        "limitations": [
            "The spherical inversion gives negative total density over 14.5-16.0 kpc for multiple interpolators; no clipping is allowed.",
            "No local 3-D gas surface density, gas thickness, or source deprojection is staged.",
            "Center extrapolation and outer taper are declared assumptions, not observed data.",
            "Chi-square uses Rotmod errV only; no full covariance, inclination/distance systematics, or baryonic profile covariance is available.",
            "This NGC 3198 is a diagnostic, not a holdout; success would not establish thermodynamic gravity or a relativistic cosmology.",
        ],
        "runtime_seconds": time.monotonic() - started,
        "command": (
            ".venv/bin/python scripts/run_bounded.py --state RUN_STATE.json --seconds 600 -- "
            "env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python "
            "experiments/tg_ngc3198_isolated_surrogate/solve_isolated.py"
        ),
        "code_sha256": sha256(Path(__file__)),
    }
    (OUT / "record.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

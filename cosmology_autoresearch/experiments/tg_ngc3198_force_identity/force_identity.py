#!/usr/bin/env python3
"""Screen the stationary TG force identity against NGC 3198 Rotmod data.

Run from the project root with:
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
      .venv/bin/python experiments/tg_ngc3198_force_identity/force_identity.py

This is a spherical-surrogate, measurement-error-only consistency check, not a
three-dimensional baryon model or a physical detection.  It tests the necessary
integrated identity without using an observed outer velocity as a boundary value.
"""

from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import time
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2 as chi2_distribution


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "context" / "data" / "Rotmod_LTG.zip"
MEMBER = "NGC3198_rotmod.dat"
SEED = 20260925
N_DRAWS = 12000
K_PAPER = 3.40e-5  # s^2 / km^2
UPSILON_PAPER = 0.762
UPSILON_SIGMA = 0.007


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_linear_v4_integral(r_kpc: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    """Cumulative integral of v^4 under piecewise-linear interpolation of v."""
    left = velocity[:-1]
    right = velocity[1:]
    segment_mean_v4 = (
        left**4 + left**3 * right + left**2 * right**2
        + left * right**3 + right**4
    ) / 5.0
    increments = np.diff(r_kpc) * segment_mean_v4
    return np.concatenate(([0.0], np.cumsum(increments)))


def load_data() -> dict[str, np.ndarray]:
    with zipfile.ZipFile(DATA) as archive:
        raw = archive.read(MEMBER).decode("ascii")
    values = np.loadtxt(io.StringIO(raw), comments="#")
    if values.shape[1] != 8 or values.shape[0] < 10:
        raise RuntimeError(f"unexpected Rotmod layout: {values.shape}")
    radius, observed, error, gas, disk, bulge = values[:, :6].T
    if np.any(np.diff(radius) <= 0) or np.any(error <= 0):
        raise RuntimeError("non-monotone radii or non-positive velocity error")
    return {
        "radius_kpc": radius,
        "vobs_kms": observed,
        "errv_kms": error,
        "vgas_kms": gas,
        "vdisk_kms": disk,
        "vbul_kms": bulge,
    }


def cumulative_cross_covariance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Cov[a_i,b_j] from shared Monte Carlo draws, with variables on axis 1."""
    ac = a - a.mean(axis=0, keepdims=True)
    bc = b - b.mean(axis=0, keepdims=True)
    return ac.T @ bc / (a.shape[0] - 1)


def make_profile(y_observed, x_observed, y_draws, x_draws):
    """Build a generalized least-squares K profile with a free integration constant."""
    mean_y = y_draws.mean(axis=0)
    mean_x = x_draws.mean(axis=0)
    cov_yy = cumulative_cross_covariance(y_draws, y_draws)
    cov_xx = cumulative_cross_covariance(x_draws, x_draws)
    cov_yx = cumulative_cross_covariance(y_draws, x_draws)
    ones = np.ones(y_observed.size)

    def profile(k: float, return_details: bool = False):
        covariance = cov_yy - k * (cov_yx + cov_yx.T) + k * k * cov_xx
        covariance = 0.5 * (covariance + covariance.T)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        cutoff = max(eigenvalues[-1] * 1e-10, 0.0)
        keep = eigenvalues > cutoff
        if np.count_nonzero(keep) < 4:
            return (np.inf, {}) if return_details else np.inf
        precision = (eigenvectors[:, keep] / eigenvalues[keep]) @ eigenvectors[:, keep].T

        observed_residual = y_observed - k * x_observed
        # Correct small nonlinear transform bias estimated by the same draws
        # used for the shared-error covariance.
        expected_draw_residual = mean_y - k * mean_x
        adjusted = 2.0 * observed_residual - expected_draw_residual
        p_one = precision @ ones
        p_residual = precision @ adjusted
        constant = float((ones @ p_residual) / (ones @ p_one))
        residual = adjusted - constant
        value = float(residual @ precision @ residual)
        if not return_details:
            return value
        return value, {
            "intercept_kpc_km2_s2": constant,
            "effective_covariance_rank": int(np.count_nonzero(keep)),
            "covariance_min_kept_eigenvalue": float(eigenvalues[keep][0]),
            "covariance_max_eigenvalue": float(eigenvalues[-1]),
            "residual": residual,
            "covariance": covariance,
        }

    return profile


def main() -> None:
    started = time.monotonic()
    data = load_data()
    r = data["radius_kpc"]
    v = data["vobs_kms"]
    sigma_v = data["errv_kms"]
    gas = data["vgas_kms"]
    disk = data["vdisk_kms"]
    bulge = data["vbul_kms"]

    # SPARC stores gas force as a signed contribution to v^2.  In particular,
    # a negative Vgas is not to be squared into a positive source contribution.
    baryon_v2 = (
        UPSILON_PAPER * disk**2
        + gas * np.abs(gas)
        + UPSILON_PAPER * bulge**2
    )
    y_obs = r * (v**2 - baryon_v2)
    x_obs = exact_linear_v4_integral(r, v)

    rng = np.random.default_rng(SEED)
    v_draws = rng.normal(v, sigma_v, size=(N_DRAWS, v.size))
    upsilon_draws = rng.normal(UPSILON_PAPER, UPSILON_SIGMA, size=(N_DRAWS, 1))
    baryon_v2_draws = (
        upsilon_draws * (disk**2 + bulge**2)[None, :]
        + (gas * np.abs(gas))[None, :]
    )
    y_draws = r[None, :] * (v_draws**2 - baryon_v2_draws)
    x_draws = np.stack(
        [exact_linear_v4_integral(r, sample) for sample in v_draws], axis=0
    )

    # The inner table begins at 0.32 kpc, so marginalize a constant of
    # integration instead of inventing an unobserved central boundary value.
    profile = make_profile(y_obs, x_obs, y_draws, x_draws)

    result = minimize_scalar(
        profile, bounds=(0.0, 2.0e-4), method="bounded",
        options={"xatol": 1e-11, "maxiter": 1000},
    )
    if not result.success:
        raise RuntimeError(f"K profile minimization failed: {result.message}")
    k_best = float(result.x)
    chi2_best, details_best = profile(k_best, return_details=True)
    chi2_paper, details_paper = profile(K_PAPER, return_details=True)
    rank = details_best["effective_covariance_rank"]
    dof_best = max(rank - 2, 1)  # profiled K and integration constant
    dof_paper = max(details_paper["effective_covariance_rank"] - 1, 1)

    k_grid = np.linspace(0.0, 2.0e-4, 401)
    profile_grid = np.array([profile(float(k)) for k in k_grid])
    np.savetxt(
        OUT / "k_profile.csv",
        np.column_stack((k_grid, profile_grid)),
        delimiter=",", header="K_s2_per_km2,profile_chi2", comments="",
    )

    # Fixed, diagnostic radial-window checks: they test whether the shape result
    # is driven by inner points with large errV or by the last outer points. These
    # are not additional independent fits or a multiple-testing search.
    windows = {
        "r_ge_5_kpc": r >= 5.0,
        "5_le_r_le_35_kpc": (r >= 5.0) & (r <= 35.0),
        "r_le_35_kpc": r <= 35.0,
    }
    window_results = {}
    for label, selected in windows.items():
        window_profile = make_profile(
            y_obs[selected], x_obs[selected], y_draws[:, selected], x_draws[:, selected]
        )
        window_best = minimize_scalar(
            window_profile, bounds=(0.0, 2.0e-4), method="bounded",
            options={"xatol": 1e-11, "maxiter": 1000},
        )
        window_paper = window_profile(K_PAPER)
        window_chi2, window_details = window_profile(
            float(window_best.x), return_details=True
        )
        window_dof = max(window_details["effective_covariance_rank"] - 2, 1)
        window_results[label] = {
            "n_points": int(np.count_nonzero(selected)),
            "best_K_s2_per_km2": float(window_best.x),
            "profile_chi2_best": float(window_chi2),
            "approx_dof": window_dof,
            "approx_p_value_measurement_errors_only": float(
                chi2_distribution.sf(window_chi2, window_dof)
            ),
            "paper_K_profile_chi2": float(window_paper),
            "paper_minus_best_delta_profile_chi2": float(window_paper - window_chi2),
        }

    # A linear interpolation in v^4 instead of exact integration for linear v
    # supplies a second low-order quadrature check of the accumulated identity.
    x_trap = np.concatenate((
        [0.0], np.cumsum(np.diff(r) * 0.5 * (v[:-1]**4 + v[1:]**4))
    ))
    x_trap_draws = np.stack([
        np.concatenate((
            [0.0], np.cumsum(
                np.diff(r) * 0.5 * (sample[:-1]**4 + sample[1:]**4)
            )
        )) for sample in v_draws
    ])
    trapezoid_profile = make_profile(y_obs, x_trap, y_draws, x_trap_draws)
    trapezoid_best = minimize_scalar(
        trapezoid_profile, bounds=(0.0, 2.0e-4), method="bounded",
        options={"xatol": 1e-11, "maxiter": 1000},
    )
    trapezoid_chi2 = float(trapezoid_profile(float(trapezoid_best.x)))

    upsilon_sensitivity = {}
    for label, upsilon in (
        ("mean_minus_1sigma", UPSILON_PAPER - UPSILON_SIGMA),
        ("mean", UPSILON_PAPER),
        ("mean_plus_1sigma", UPSILON_PAPER + UPSILON_SIGMA),
    ):
        fixed_baryon_v2 = (
            upsilon * (disk**2 + bulge**2) + gas * np.abs(gas)
        )
        fixed_y_obs = r * (v**2 - fixed_baryon_v2)
        fixed_y_draws = r[None, :] * (
            v_draws**2 - fixed_baryon_v2[None, :]
        )
        fixed_profile = make_profile(fixed_y_obs, x_obs, fixed_y_draws, x_draws)
        fixed_best = minimize_scalar(
            fixed_profile, bounds=(0.0, 2.0e-4), method="bounded",
            options={"xatol": 1e-11, "maxiter": 1000},
        )
        fixed_chi2 = float(fixed_profile(float(fixed_best.x)))
        upsilon_sensitivity[label] = {
            "Upsilon_disk": upsilon,
            "best_K_s2_per_km2": float(fixed_best.x),
            "profile_chi2_best_velocity_errors_only": fixed_chi2,
            "paper_K_profile_chi2_velocity_errors_only": float(
                fixed_profile(K_PAPER)
            ),
        }

    # Approximate profile Δχ²=1 interval from the stored scan, if both crossings
    # lie inside the predeclared nonnegative K domain.
    delta_grid = profile_grid - chi2_best - 1.0
    crossings = []
    for left in range(len(k_grid) - 1):
        if delta_grid[left] == 0.0:
            crossings.append(float(k_grid[left]))
        elif delta_grid[left] * delta_grid[left + 1] < 0.0:
            fraction = -delta_grid[left] / (delta_grid[left + 1] - delta_grid[left])
            crossings.append(float(k_grid[left] + fraction * (k_grid[left + 1] - k_grid[left])))
    interval_68 = [crossings[0], crossings[-1]] if len(crossings) >= 2 else None

    fig, ax = plt.subplots(figsize=(7.4, 5.2), constrained_layout=True)
    ax.plot(x_obs, y_obs, "o", ms=4.3, color="#285f9e", label="SPARC central values")
    order = np.argsort(x_obs)
    fit_y = details_best["intercept_kpc_km2_s2"] + k_best * x_obs
    paper_y = details_paper["intercept_kpc_km2_s2"] + K_PAPER * x_obs
    ax.plot(x_obs[order], fit_y[order], color="#bd3b32", lw=2, label="profiled K")
    ax.plot(
        x_obs[order], paper_y[order], color="#ee9b24", lw=1.7, ls="--",
        label="paper K = 3.40 × 10⁻⁵",
    )
    ax.set_xlabel(r"$X(r)=\int_{r_{min}}^r v_{obs}^4\,dr$  [kpc (km/s)$^4$]")
    ax.set_ylabel(r"$Y(r)=r(v_{obs}^2-v_{bar}^2)$  [kpc (km/s)$^2$]")
    ax.set_title("NGC 3198: integrated stationary TG force identity")
    ax.legend(frameon=False)
    ax.text(
        0.02, 0.98,
        "full shared velocity-error covariance;\n"
        "measurement errors only; spherical surrogate",
        transform=ax.transAxes, va="top", fontsize=8.5,
    )
    fig.savefig(OUT / "force_identity.png", dpi=170)
    plt.close(fig)

    record = {
        "status": "measurement_error_screen_complete",
        "scope": (
            "Necessary integrated identity for the printed stationary TG equation; "
            "conditional pseudo-spherical SPARC source; not a physical 3D fit or detection."
        ),
        "source": {
            "archive": "context/data/Rotmod_LTG.zip",
            "member": MEMBER,
            "archive_sha256": sha256(DATA),
            "rows": int(r.size),
            "radial_range_kpc": [float(r[0]), float(r[-1])],
        },
        "model": {
            "identity": "d[r*(v_obs^2-v_bar^2)]/dr = K*v_obs^4",
            "integral": "Y(r)=C+K*X(r), X=integral(v_obs^4 dr)",
            "gas_convention": "SPARC signed v^2 contribution: Vgas*abs(Vgas)",
            "baryon_v2": "Upsilon_disk*(Vdisk^2+Vbul^2)+Vgas*abs(Vgas)",
            "upsilon_disk_mean": UPSILON_PAPER,
            "upsilon_disk_sigma": UPSILON_SIGMA,
            "integration": "exact segment integral under piecewise-linear interpolation of Vobs",
            "inner_boundary": (
                "profiled constant of integration because the first datum is at 0.32 kpc; "
                "no observed outer endpoint used as boundary"
            ),
            "K_units": "s^2/km^2",
        },
        "uncertainties": {
            "method": "parametric Monte Carlo covariance of shared Y and cumulative X errors",
            "seed": SEED,
            "draws": N_DRAWS,
            "observed_velocity_errors": "independent Gaussian Rotmod errV; no extra covariance supplied",
            "stellar_mass_to_light": "common Gaussian nuisance Upsilon=0.762 +/- 0.007",
            "omitted": [
                "gas/disk profile covariance",
                "distance and inclination systematics",
                "interpolation-model uncertainty",
                "non-spherical geometry uncertainty",
            ],
        },
        "fit": {
            "K_domain_s2_per_km2": [0.0, 2.0e-4],
            "best_K_s2_per_km2": k_best,
            "approx_profile_delta_chi2_1_interval": interval_68,
            "paper_K_s2_per_km2": K_PAPER,
            "profile_chi2_best": float(chi2_best),
            "profile_chi2_paper_K": float(chi2_paper),
            "delta_profile_chi2_paper_minus_best": float(chi2_paper - chi2_best),
            "effective_covariance_rank": rank,
            "dof_best_approx": dof_best,
            "p_value_best_approx": float(chi2_distribution.sf(chi2_best, dof_best)),
            "dof_paper_K_approx": dof_paper,
            "p_value_paper_K_approx": float(chi2_distribution.sf(chi2_paper, dof_paper)),
            "fixed_diagnostic_radial_windows": window_results,
            "Upsilon_fixed_value_sensitivity": upsilon_sensitivity,
            "quadrature_sensitivity": {
                "alternative": "trapezoidal integral of v^4 versus exact segment integral under linear v",
                "best_K_s2_per_km2_trapezoid": float(trapezoid_best.x),
                "profile_chi2_best_trapezoid": trapezoid_chi2,
                "delta_best_K_relative_to_exact": float(trapezoid_best.x - k_best),
            },
            "profiled_intercept_best_kpc_km2_s2": details_best[
                "intercept_kpc_km2_s2"
            ],
            "profiled_intercept_paper_kpc_km2_s2": details_paper[
                "intercept_kpc_km2_s2"
            ],
            "profile_csv": "experiments/tg_ngc3198_force_identity/k_profile.csv",
            "plot": "experiments/tg_ngc3198_force_identity/force_identity.png",
            "interpretation_caveat": (
                "Approximate chi-square tail uses a Monte Carlo-estimated, parameter-dependent "
                "covariance and omits substantial systematics; it is a screening statistic, "
                "not a calibrated rejection significance."
            ),
        },
        "checks": {
            "K_profile_optimizer_success": bool(result.success),
            "profile_chi2_grid_finite": bool(np.all(np.isfinite(profile_grid))),
            "min_K_nonnegative": bool(k_best >= 0.0),
            "all_Rotmod_errV_positive": bool(np.all(sigma_v > 0.0)),
        },
        "runtime_seconds": time.monotonic() - started,
        "command": (
            "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "
            ".venv/bin/python experiments/tg_ngc3198_force_identity/force_identity.py"
        ),
        "code_sha256": sha256(Path(__file__)),
    }
    (OUT / "record.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

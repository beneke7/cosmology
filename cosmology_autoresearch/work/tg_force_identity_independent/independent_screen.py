#!/usr/bin/env python3
"""Independent integrated stationary TG identity screen for NGC 3198.

This is a spherical-surrogate necessary-condition screen.  It does not solve a
disk geometry or make a physical rotation-curve prediction.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import zipfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ARCHIVE = REPO / "context/data/Rotmod_LTG.zip"
PDF = REPO / "context/user_provided/pszota_van_2024.pdf"
MEMBER = "NGC3198_rotmod.dat"
UPSILON = 0.762
UPSILON_SIGMA = 0.007
PUBLISHED_K = 3.40e-5  # Table 1: 3.40 +/- 0.04, in units 1e-5 s^2 km^-2.
PUBLISHED_K_SIGMA = 0.04e-5
SEED = 20260925
MONTE_CARLO_DRAWS = 100_000


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_data() -> tuple[dict[str, np.ndarray], bytes]:
    archive_bytes = ARCHIVE.read_bytes()
    with zipfile.ZipFile(io.BytesIO(archive_bytes)) as archive:
        raw = archive.read(MEMBER)
    values = np.genfromtxt(io.StringIO(raw.decode("ascii")), comments="#")
    names = ("r_kpc", "vobs_kms", "err_kms", "vgas_kms", "vdisk_kms", "vbul_kms", "sb_disk", "sb_bul")
    data = {name: values[:, index].astype(float) for index, name in enumerate(names)}
    return data, raw


def exact_linear_panel_integral_v4(v0: np.ndarray, v1: np.ndarray, dr: float) -> np.ndarray:
    """Exact integral of V^4 where V varies linearly between panel endpoints."""
    return dr * (
        v0**4 + v0**3 * v1 + v0**2 * v1**2 + v0 * v1**3 + v1**4
    ) / 5.0


def gauss_panel_integral_v4(v0: np.ndarray, v1: np.ndarray, dr: float, order: int = 8) -> np.ndarray:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = (nodes + 1.0) / 2.0
    w = weights / 2.0
    vals = v0[..., None] + (v1 - v0)[..., None] * t
    return dr * np.sum(w * vals**4, axis=-1)


def cumulative_integral(r: np.ndarray, v: np.ndarray, method: str = "exact") -> np.ndarray:
    panels = []
    for j, dr in enumerate(np.diff(r)):
        if method == "exact":
            q = exact_linear_panel_integral_v4(v[..., j], v[..., j + 1], float(dr))
        elif method == "gauss8":
            q = gauss_panel_integral_v4(v[..., j], v[..., j + 1], float(dr), order=8)
        else:
            raise ValueError(method)
        panels.append(q)
    panel_array = np.stack(panels, axis=-1)
    zero_shape = (*panel_array.shape[:-1], 1)
    zero = np.zeros(zero_shape, dtype=float)
    return np.concatenate((zero, np.cumsum(panel_array, axis=-1)), axis=-1)


def source_square(data: dict[str, np.ndarray], upsilon: np.ndarray | float) -> np.ndarray:
    # SPARC's gas term is signed: a negative Vgas means negative contribution
    # to total v^2.  NGC 3198's bulge column is zero throughout.
    gas = data["vgas_kms"] * np.abs(data["vgas_kms"])
    return upsilon * data["vdisk_kms"] ** 2 + gas


def line_fit(x: np.ndarray, y: np.ndarray) -> dict[str, float | np.ndarray]:
    xc = x - np.mean(x)
    yc = y - np.mean(y)
    slope = float(np.dot(xc, yc) / np.dot(xc, xc))
    k = -slope
    constant = float(np.mean(y) + k * np.mean(x))
    residual = y + k * x - constant
    scale = float(np.sqrt(np.mean((x * 0.0 + y) ** 2)))
    return {
        "K_s2_per_km2": k,
        "constant_kpc_km2_per_s2": constant,
        "residual": residual,
        "rms_residual_kpc_km2_per_s2": float(np.sqrt(np.mean(residual**2))),
        "rms_fraction_of_r_source_minus_observed": float(np.sqrt(np.mean(residual**2)) / scale),
        "max_abs_residual_kpc_km2_per_s2": float(np.max(np.abs(residual))),
        "r_squared": float(1.0 - np.sum(residual**2) / np.sum((y - np.mean(y)) ** 2)),
    }


def cumulative_jacobian_v4(r: np.ndarray, v: np.ndarray) -> np.ndarray:
    """d X_i / d V_j for X_i=int[r0,r_i] V^4 dr under linear panels."""
    n = len(r)
    out = np.zeros((n, n), dtype=float)
    for panel, dr in enumerate(np.diff(r)):
        a, b = v[panel], v[panel + 1]
        da = dr / 5.0 * sum((4 - m) * a ** (3 - m) * b**m for m in range(4))
        db = dr / 5.0 * sum(m * a ** (4 - m) * b ** (m - 1) for m in range(1, 5))
        out[panel + 1 :, panel] += da
        out[panel + 1 :, panel + 1] += db
    return out


def covariance_identity(
    data: dict[str, np.ndarray],
    k: float,
    jx: np.ndarray,
) -> np.ndarray:
    r = data["r_kpc"]
    v = data["vobs_kms"]
    err = data["err_kms"]
    jy = -2.0 * r * v
    jac = k * jx
    jac[np.diag_indices_from(jac)] += jy
    cov = (jac * err[None, :] ** 2) @ jac.T
    upsilon_jac = r * data["vdisk_kms"] ** 2
    cov += UPSILON_SIGMA**2 * np.outer(upsilon_jac, upsilon_jac)
    return cov


def feasible_gls(
    x: np.ndarray,
    y: np.ndarray,
    covariance: np.ndarray,
) -> dict[str, float]:
    """One-pass GLS line fit with covariance frozen at the central OLS K."""
    design = np.column_stack((np.ones_like(x), -x))
    # Cholesky is preferred; a small scale-relative diagonal regularizer is
    # used only if floating-point conditioning requires it.
    cov = covariance.copy()
    reg = 0.0
    for _ in range(8):
        try:
            chol = np.linalg.cholesky(cov + reg * np.eye(len(cov)))
            break
        except np.linalg.LinAlgError:
            reg = max(1e-14 * float(np.max(np.diag(cov))), 10.0 * reg if reg else 0.0)
    else:
        raise np.linalg.LinAlgError("identity covariance is not positive definite")
    inv_design = np.linalg.solve(chol.T, np.linalg.solve(chol, design))
    inv_y = np.linalg.solve(chol.T, np.linalg.solve(chol, y))
    beta_cov = np.linalg.inv(design.T @ inv_design)
    beta = beta_cov @ (design.T @ inv_y)
    residual = y - design @ beta
    whitened = np.linalg.solve(chol, residual)
    dof = len(y) - 2
    return {
        "constant_kpc_km2_per_s2": float(beta[0]),
        "K_s2_per_km2": float(beta[1]),
        "K_sigma_linearized_s2_per_km2": float(math.sqrt(beta_cov[1, 1])),
        "chi2": float(whitened @ whitened),
        "dof": int(dof),
        "reduced_chi2": float((whitened @ whitened) / dof),
        "covariance_regularizer": float(reg),
    }


def monte_carlo_fits(data: dict[str, np.ndarray], x_center: np.ndarray, y_center: np.ndarray, masks: dict[str, np.ndarray]) -> dict[str, dict[str, float]]:
    rng = np.random.default_rng(SEED)
    n_draw = MONTE_CARLO_DRAWS
    v_draws = rng.normal(data["vobs_kms"], data["err_kms"], size=(n_draw, len(data["r_kpc"])))
    upsilon_draws = rng.normal(UPSILON, UPSILON_SIGMA, size=n_draw)
    r = data["r_kpc"]
    b_draws = upsilon_draws[:, None] * data["vdisk_kms"][None, :] ** 2
    b_draws += data["vgas_kms"][None, :] * np.abs(data["vgas_kms"][None, :])
    y_draws = r[None, :] * (b_draws - v_draws**2)
    x_draws = cumulative_integral(r, v_draws, "exact")
    result: dict[str, dict[str, float]] = {}
    for name, mask in masks.items():
        xx = x_draws[:, mask]
        yy = y_draws[:, mask]
        xc = xx - np.mean(xx, axis=1)[:, None]
        yc = yy - np.mean(yy, axis=1)[:, None]
        k_draws = -(np.sum(xc * yc, axis=1) / np.sum(xc * xc, axis=1))
        c_draws = np.mean(yy, axis=1) + k_draws * np.mean(xx, axis=1)
        res = yy + k_draws[:, None] * xx - c_draws[:, None]
        rms = np.sqrt(np.mean(res**2, axis=1))
        result[name] = {
            "K_median_s2_per_km2": float(np.median(k_draws)),
            "K_mean_s2_per_km2": float(np.mean(k_draws)),
            "K_sd_s2_per_km2": float(np.std(k_draws, ddof=1)),
            "K_q16_s2_per_km2": float(np.quantile(k_draws, 0.16)),
            "K_q84_s2_per_km2": float(np.quantile(k_draws, 0.84)),
            "fraction_K_positive": float(np.mean(k_draws > 0.0)),
            "constant_median_kpc_km2_per_s2": float(np.median(c_draws)),
            "residual_rms_median_kpc_km2_per_s2": float(np.median(rms)),
            "residual_rms_q16_kpc_km2_per_s2": float(np.quantile(rms, 0.16)),
            "residual_rms_q84_kpc_km2_per_s2": float(np.quantile(rms, 0.84)),
            "draws": n_draw,
            "seed": SEED,
        }
    return result


def clean_fit(fit: dict[str, object]) -> dict[str, object]:
    return {k: v for k, v in fit.items() if k != "residual"}


def main() -> None:
    data, member_bytes = load_data()
    r = data["r_kpc"]
    v = data["vobs_kms"]
    b = source_square(data, UPSILON)
    x = cumulative_integral(r, v, "exact")
    x_alt = cumulative_integral(r, v, "gauss8")
    y = r * (b - v**2)
    mask_full = np.ones(len(r), dtype=bool)
    mask_outer = r >= 5.0
    masks = {"full_0p32_to_44p08_kpc": mask_full, "observed_points_r_ge_5_kpc": mask_outer}
    jx = cumulative_jacobian_v4(r, v)
    fits: dict[str, dict[str, object]] = {}
    for name, mask in masks.items():
        fit = line_fit(x[mask], y[mask])
        cov = covariance_identity(data, float(fit["K_s2_per_km2"]), jx)
        gls = feasible_gls(x[mask], y[mask], cov[np.ix_(mask, mask)])
        fit["linearized_error_covariance_gls"] = gls
        fit["N_points"] = int(mask.sum())
        fit["r_min_kpc"] = float(r[mask][0])
        fit["r_max_kpc"] = float(r[mask][-1])
        fit["published_K_difference_sigma"] = float((float(fit["K_s2_per_km2"]) - PUBLISHED_K) / PUBLISHED_K_SIGMA)
        fits[name] = fit
    mc = monte_carlo_fits(data, x, y, masks)
    max_abs_diff = float(np.max(np.abs(x - x_alt)))
    max_rel_diff = float(max_abs_diff / np.max(np.abs(x)))
    report = {
        "title": "Independent integrated stationary thermodynamic-gravity identity screen: NGC 3198",
        "classification": "spherical-surrogate necessary-condition screen; not a physical NGC 3198 prediction",
        "derivation": [
            "Paper Eq. (8): g' = -2g/r - 4 pi G rho - K g^2.",
            "Define b(r)=Upsilon Vdisk(r)^2 + Vgas(r)*abs(Vgas(r)); the paper's spherical surrogate gives 4 pi G rho = (r b)'/r^2.",
            "Multiply by r^2 and integrate: r^2 g + r b + K integral[rmin,r] s^2 g(s)^2 ds = C.",
            "Set g=-Vobs^2/r: r*(b - Vobs^2) + K integral[rmin,r] Vobs^4 ds = C.",
        ],
        "units": {
            "r": "kpc",
            "velocities": "km s^-1",
            "K": "s^2 km^-2",
            "each_identity_term": "kpc km^2 s^-2",
            "integral_intermediate": "kpc km^4 s^-4",
        },
        "input_hashes": {
            "Rotmod_LTG.zip_sha256": sha256_bytes(ARCHIVE.read_bytes()),
            "NGC3198_rotmod.dat_member_sha256": sha256_bytes(member_bytes),
            "pszota_van_2024.pdf_sha256": sha256_bytes(PDF.read_bytes()),
        },
        "input_data": {
            "source": "exact local Rotmod_LTG.zip member NGC3198_rotmod.dat",
            "rows": int(len(r)),
            "radial_range_kpc": [float(r[0]), float(r[-1])],
            "column_convention": "Rad, Vobs, errV, Vgas, Vdisk, Vbul, SBdisk, SBbul; units from archive header",
            "bulge_nonzero_rows": int(np.count_nonzero(data["vbul_kms"])),
        },
        "fixed_stellar_normalization": {
            "Upsilon": UPSILON,
            "sigma": UPSILON_SIGMA,
            "source": "Pszota & Van 2024 Table 1, NGC 3198 thermodynamic gravity fit",
            "treatment": "fixed at published central value; its quoted uncertainty is propagated as one shared nuisance in the linearized covariance and the Monte Carlo draws; not refit",
        },
        "gas_convention": "SPARC signed gas contribution is Vgas*abs(Vgas), so negative Vgas remains negative in b(r).",
        "integration": {
            "lower_limit_kpc": float(r[0]),
            "boundary_condition": "free integration constant C fitted on each selected range; no outer observed endpoint imposed as a boundary",
            "primary": "exact panel integral of Vobs^4 under piecewise-linear interpolation between adjacent observed points",
            "secondary": "independent 8-point Gauss-Legendre quadrature of the same piecewise-linear interpolant",
            "max_abs_cumulative_integral_difference_kpc_km4_per_s4": max_abs_diff,
            "max_relative_cumulative_integral_difference": max_rel_diff,
        },
        "residual_metric": {
            "primary_fit": "unweighted least-squares line Y=C-KX, with Y=r*(b-Vobs^2), X=int[rmin,r] Vobs^4 dr; report RMS(Y+KX-C) and its ratio to RMS(Y).",
            "uncertainty_fit": "single-pass feasible GLS with full first-order covariance from every independent observed-velocity error in both local Y and all shared cumulative X values, plus one shared Upsilon uncertainty; covariance frozen at central OLS K; chi2 uses two fitted coefficients.",
            "monte_carlo": "100,000 independent Gaussian draws using every row's errV jointly across all cumulative integrals and one shared Gaussian Upsilon draw; each draw refits C and K by unweighted OLS.",
            "limitations": "the propagated GLS covariance is linearized; the Monte Carlo distribution addresses nonlinear powers but assumes independent Gaussian input errors and fixed component curves.",
        },
        "published_comparison": {
            "paper_K_s2_per_km2": PUBLISHED_K,
            "paper_K_sigma_s2_per_km2": PUBLISHED_K_SIGMA,
            "reference": "Pszota & Van 2024 Table 1: K=(3.40 +/- 0.04) x 10^-5 s^2 km^-2; Upsilon=0.762 +/- 0.007.",
        },
        "ranges": {},
        "monte_carlo": mc,
        "assumptions_and_caveats": [
            "This only checks a necessary first-integral condition after replacing a disk galaxy with the paper's sphericalized source prescription.",
            "The exact interpolation assumption is piecewise linear Vobs(r) within each adjacent-data panel; it is numerical quadrature, not an asserted continuous observational model.",
            "No per-point errors are supplied for Vgas, Vdisk, or their geometry conversion, so these component curves are held fixed; only Upsilon's published uncertainty is propagated.",
            "Velocity errors are treated as independent Gaussian errors in Vobs; distance, inclination, and calibration uncertainties are not included.",
            "The inner point has errV greater than Vobs, so Gaussian propagation through Vobs^2 and Vobs^4 is strongly nonlinear there; the Monte Carlo summaries are more informative than the first-order GLS metric for the full range.",
            "The free integration constant is required because the first observed point is r_min=0.32 kpc; no condition is specified inside it. The last observed point is not used as a boundary condition.",
        ],
    }
    for name, fit in fits.items():
        report["ranges"][name] = clean_fit(fit)

    (HERE / "independent_result.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    lines = [
        "# Independent integrated stationary TG force-identity screen: NGC 3198",
        "",
        "Classification: spherical-surrogate necessary-condition screen, not a physical NGC 3198 prediction.",
        "",
        "## Derivation",
        "",
        "Pszota & Van (2024), Eq. (8), gives `g' = -2g/r - 4 pi G rho - K g^2`. Under the paper's sphericalized source construction, define `b(r)=Upsilon Vdisk^2 + Vgas*abs(Vgas)`, so `4 pi G rho = (r b)'/r^2`. Multiplication by `r^2` and integration from `r_min` gives `r^2 g + r b + K integral[r_min,r] s^2 g(s)^2 ds = C`. With `g=-Vobs^2/r`, the checked identity is:",
        "",
        "`r * (b - Vobs^2) + K * integral[r_min,r] Vobs^4 ds = C`.",
        "",
        "The constant is free because the first observed point is at 0.32 kpc. No outer observed endpoint is imposed as a boundary.",
        "",
        "## Frozen independent calculation",
        "",
        f"The exact local archive member has {len(r)} rows over {r[0]:.2f}–{r[-1]:.2f} kpc. I used the published fixed stellar normalization `Upsilon=0.762 +/- 0.007`, with the quoted uncertainty propagated as a shared nuisance, and the SPARC signed-gas term `Vgas*abs(Vgas)`. The bulge column is zero for all rows.",
        "",
        f"The cumulative `Vobs^4` integral uses the exact polynomial panel integral for linearly interpolated Vobs; an independent 8-point Gauss-Legendre calculation agrees to maximum relative difference {max_rel_diff:.3g}.",
        "",
        "Primary residual is the RMS of `Y+KX-C`, with `Y=r*(b-Vobs^2)`, `X=integral Vobs^4`, `C` free, and K fit as the slope; the reported fractional residual divides by RMS(Y). GLS uses the full covariance from shared observed-velocity errors in both Y and every cumulative X and includes the shared Upsilon error. A seeded 100,000-draw Monte Carlo refits C and K per draw to reflect nonlinear measurement propagation.",
        "",
        "| Range | N | OLS K (s²/km²) | C (kpc km²/s²) | RMS residual | RMS residual / RMS(Y) | GLS K | GLS reduced χ² |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, fit in fits.items():
        gls = fit["linearized_error_covariance_gls"]
        lines.append(
            f"| {name} | {fit['N_points']} | {fit['K_s2_per_km2']:.8g} | {fit['constant_kpc_km2_per_s2']:.8g} | {fit['rms_residual_kpc_km2_per_s2']:.8g} | {fit['rms_fraction_of_r_source_minus_observed']:.8g} | {gls['K_s2_per_km2']:.8g} | {gls['reduced_chi2']:.8g} |"
        )
    lines += ["", "Monte Carlo K and residual summaries:", ""]
    for name, stats in mc.items():
        lines.append(
            f"- {name}: K median {stats['K_median_s2_per_km2']:.8g}, mean {stats['K_mean_s2_per_km2']:.8g}, SD {stats['K_sd_s2_per_km2']:.8g}, 16–84% [{stats['K_q16_s2_per_km2']:.8g}, {stats['K_q84_s2_per_km2']:.8g}], P(K>0)={stats['fraction_K_positive']:.6g}; median residual RMS {stats['residual_rms_median_kpc_km2_per_s2']:.8g}.")
    lines += [
        "",
        f"Published Table 1 comparison: K={PUBLISHED_K:.8g} +/- {PUBLISHED_K_SIGMA:.2g} s²/km² and Upsilon={UPSILON} +/- {UPSILON_SIGMA}. The screen does not refit Upsilon.",
        "",
        "## Numerical and physical caveats",
        "",
        "- This checks a necessary integrated identity only after replacing the disk geometry by the paper's sphericalized source prescription; it is not a physical galaxy prediction.",
        "- Piecewise-linear interpolation of Vobs is used only to evaluate the cumulative integral; source values are evaluated at the data radii.",
        "- Per-row uncertainties for the gas/disk component curves are unavailable and those curves are held fixed. Distance, inclination, and calibration uncertainties are omitted.",
        "- Observed-velocity errors are independent Gaussian errors in Vobs. The innermost errV exceeds Vobs, so powers Vobs² and Vobs⁴ are nonlinear; Monte Carlo should be favored over the linearized GLS uncertainty for that range.",
        "- The integration constant is free from r_min=0.32 kpc. The outermost observed point is not a boundary condition.",
        "",
        "## Provenance and reproducibility",
        "",
        f"- ZIP SHA-256: `{report['input_hashes']['Rotmod_LTG.zip_sha256']}`",
        f"- NGC3198_rotmod.dat member SHA-256: `{report['input_hashes']['NGC3198_rotmod.dat_member_sha256']}`",
        f"- Supplied PDF SHA-256: `{report['input_hashes']['pszota_van_2024.pdf_sha256']}`",
        f"- Monte Carlo seed: `{SEED}`; draws: `{MONTE_CARLO_DRAWS}`.",
        "- Derivation reference: Pszota & Van (2024), Eqs. (8), (47), (53), Table 1; local PDF supplied in the workspace.",
    ]
    (HERE / "report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

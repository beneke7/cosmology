#!/usr/bin/env python3
"""Profile-only DES-Dovekie flat-background comparison.

This is an independent implementation of the pinned DES likelihood semantics.
The released file is SNANA whitespace text despite its .csv suffix, while the
released CosmoSIS module asks Astropy to parse comma-separated CSV.  This
script parses the official row format without reordering and uses the official
packed STAT+SYS inverse covariance directly.

The calculation is a late-time background screen: radiation and curvature are
omitted, the absolute magnitude is profiled analytically, and no H0 calibration
or posterior inference is attempted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "context" / "data"
HD_PATH = DATA_DIR / "des_dovekie_hd.csv"
PRECISION_PATH = DATA_DIR / "des_dovekie_stat_sys.npz"
OUTPUT_DIR = ROOT / "experiments" / "dovekie_screen"
PINNED_COMMIT = "c9a4fcafc4cbd19bd750dee47fc76194a45c181f"
EXPECTED_HASHES = {
    HD_PATH.name: "2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf",
    PRECISION_PATH.name: "ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b",
}
H0_GAUGE_KM_S_MPC = 70.0
C_KM_S = 299792.458
QUADRATURE_ORDER = 96
SEED = 20260924
BOUNDS: dict[str, list[tuple[float, float]]] = {
    "LCDM": [(0.05, 0.6)],
    "constant_w": [(0.05, 0.6), (-2.0, -0.3)],
    "CPL": [(0.05, 0.6), (-2.0, -0.3), (-3.0, 3.0)],
}
GL_QUADRATURE = {
    order: tuple(np.asarray(x, dtype=np.float64) for x in leggauss(order))
    for order in (QUADRATURE_ORDER, 160)
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_snana_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Parse VARNAMES/SN rows verbatim, preserving the covariance row order."""
    columns: list[str] | None = None
    records: list[dict[str, str]] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = raw.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "VARNAMES:":
            if columns is not None:
                raise ValueError(f"{path}:{line_number}: repeated VARNAMES row")
            columns = fields[1:]
        elif fields[0] == "SN:":
            if columns is None or len(fields[1:]) != len(columns):
                raise ValueError(f"{path}:{line_number}: malformed SN row")
            records.append(dict(zip(columns, fields[1:], strict=True)))
        else:
            raise ValueError(f"{path}:{line_number}: unexpected record {fields[0]!r}")
    if columns is None:
        raise ValueError(f"{path}: no VARNAMES header")
    return columns, records


def load_inputs() -> dict[str, Any]:
    hashes = {HD_PATH.name: sha256(HD_PATH), PRECISION_PATH.name: sha256(PRECISION_PATH)}
    if hashes != EXPECTED_HASHES:
        raise ValueError(f"input hash mismatch: {hashes!r}")

    columns, records = read_snana_table(HD_PATH)
    required = {"CID", "IDSURVEY", "zHD", "zHEL", "MU", "MUERR"}
    if not required.issubset(columns):
        raise ValueError(f"missing required HD columns: {sorted(required - set(columns))}")
    cids = [row["CID"] for row in records]
    if len(set(cids)) != len(cids):
        raise ValueError("duplicate CID in Hubble diagram")

    # This is the exact official likelihood cut. It retains every released row.
    selected = [row for row in records if float(row["zHD"]) > 0.0]
    z_hd = np.asarray([float(row["zHD"]) for row in selected], dtype=np.float64)
    z_hel = np.asarray([float(row["zHEL"]) for row in selected], dtype=np.float64)
    mu_obs = np.asarray([float(row["MU"]) for row in selected], dtype=np.float64)
    muerr = np.asarray([float(row["MUERR"]) for row in selected], dtype=np.float64)
    if not (
        np.all(np.isfinite(z_hd))
        and np.all(np.isfinite(z_hel))
        and np.all(np.isfinite(mu_obs))
        and np.all(np.isfinite(muerr))
    ):
        raise ValueError("non-finite redshift or MU")
    if np.any(z_hd <= 0) or np.any(z_hel <= -1) or np.any(muerr <= 0):
        raise ValueError("invalid released redshift or MUERR")
    if len(selected) != len(records):
        # The official covariance covers the full release, so a nontrivial cut
        # would require an explicitly justified marginal covariance operation.
        raise ValueError("official zHD>0 cut would change the 1820-row covariance")

    with np.load(PRECISION_PATH, allow_pickle=False) as archive:
        n = int(np.asarray(archive["nsn"]).reshape(-1)[0])
        packed = np.asarray(archive["cov"])
        archive_keys = list(archive.files)
    if n != len(records) or packed.size != n * (n + 1) // 2:
        raise ValueError(f"precision shape does not match HD rows: nsn={n}, rows={len(records)}")
    if packed.dtype != np.float32:
        raise ValueError(f"expected packed float32 precision, got {packed.dtype}")
    precision = np.zeros((n, n), dtype=np.float64)
    upper = np.triu_indices(n)
    precision[upper] = packed.astype(np.float64)
    lower = np.tril_indices(n, -1)
    precision[lower] = precision.T[lower]
    if not np.all(np.isfinite(precision)) or not np.array_equal(precision, precision.T):
        raise ValueError("precision is not finite and exactly symmetric")
    precision_cholesky = np.linalg.cholesky(precision)
    ones = np.ones(n, dtype=np.float64)
    precision_ones = precision @ ones
    offset_information = float(ones @ precision_ones)
    if not np.isfinite(offset_information) or offset_information <= 0:
        raise ValueError("constant-offset direction is not positively constrained")

    return {
        "columns": columns,
        "records": selected,
        "z_hd": z_hd,
        "z_hel": z_hel,
        "mu_obs": mu_obs,
        "muerr": muerr,
        "precision": precision,
        "precision_cholesky": precision_cholesky,
        "precision_ones": precision_ones,
        "offset_information": offset_information,
        "archive_keys": archive_keys,
        "hashes": hashes,
    }


def e2_flat(z: np.ndarray, omega_m: float, model: str, w0: float, wa: float) -> np.ndarray:
    one_plus_z = 1.0 + z
    matter = omega_m * one_plus_z**3
    omega_de = 1.0 - omega_m
    if model == "LCDM":
        de_factor = np.ones_like(z)
    elif model == "constant_w":
        de_factor = np.exp(3.0 * (1.0 + w0) * np.log1p(z))
    elif model == "CPL":
        log_de = 3.0 * (1.0 + w0 + wa) * np.log1p(z) - 3.0 * wa * z / one_plus_z
        de_factor = np.exp(log_de)
    else:
        raise ValueError(f"unknown model {model}")
    result = matter + omega_de * de_factor
    if not np.all(np.isfinite(result)) or np.any(result <= 0.0):
        raise ValueError("E(z)^2 must remain finite and positive")
    return result


def background_mu(
    z_hd: np.ndarray,
    z_hel: np.ndarray,
    parameters: np.ndarray,
    model: str,
    *,
    h0_gauge: float = H0_GAUGE_KM_S_MPC,
    quadrature_order: int = QUADRATURE_ORDER,
) -> np.ndarray:
    """Official SN distance convention with an arbitrary fixed H0 gauge.

    D_L = (1+zHEL) (c/H0) integral_0^zHD dz/E(z); no separate (1+zHD)
    remains after the official (1+zHD) D_A(zHD) construction.
    """
    omega_m = float(parameters[0])
    w0 = float(parameters[1]) if model != "LCDM" else -1.0
    wa = float(parameters[2]) if model == "CPL" else 0.0
    if quadrature_order in GL_QUADRATURE:
        nodes, weights = GL_QUADRATURE[quadrature_order]
    else:
        nodes, weights = leggauss(quadrature_order)
    z_eval = (z_hd[:, None] / 2.0) * (nodes[None, :] + 1.0)
    integral = (z_hd / 2.0) * np.sum(
        weights[None, :] / np.sqrt(e2_flat(z_eval, omega_m, model, w0, wa)), axis=1
    )
    dl_mpc = (1.0 + z_hel) * (C_KM_S / h0_gauge) * integral
    return 5.0 * np.log10(dl_mpc) + 25.0


def profile_offset(
    mu_base: np.ndarray,
    mu_obs: np.ndarray,
    precision: np.ndarray,
    precision_ones: np.ndarray,
    offset_information: float,
) -> tuple[float, float, np.ndarray]:
    """Return the profiled additive magnitude, chi-square, and residual."""
    difference = mu_base - mu_obs
    offset = -float(precision_ones @ difference) / offset_information
    residual = difference + offset
    weighted_residual = precision @ residual
    chi2 = float(residual @ weighted_residual)
    return offset, chi2, residual


def independent_profile_check(
    difference: np.ndarray,
    precision_cholesky: np.ndarray,
    offset_information: float,
    offset: float,
    residual: np.ndarray,
    direct_chi2: float,
) -> dict[str, float]:
    """Re-evaluate the profiled quadratic in Cholesky-whitened space."""
    whitened_difference = precision_cholesky.T @ difference
    whitened_ones = precision_cholesky.T @ np.ones_like(difference)
    q_profile = float(
        whitened_difference @ whitened_difference
        - (whitened_ones @ whitened_difference) ** 2 / offset_information
    )
    whitened_residual = precision_cholesky.T @ residual
    q_residual = float(whitened_residual @ whitened_residual)
    normal_equation_residual = abs(float(whitened_ones @ (whitened_difference + offset * whitened_ones)))
    return {
        "direct_residual_quadratic": direct_chi2,
        "schur_profile_quadratic": q_profile,
        "whitened_residual_quadratic": q_residual,
        "max_absolute_disagreement": max(abs(direct_chi2 - q_profile), abs(direct_chi2 - q_residual)),
        "offset_normal_equation_abs_residual": normal_equation_residual,
    }


def fit_model(data: dict[str, Any], model: str, starts: int) -> dict[str, Any]:
    bounds = BOUNDS[model]
    precision = data["precision"]
    p1 = data["precision_ones"]
    c0 = data["offset_information"]

    def objective(parameters: np.ndarray) -> float:
        try:
            base = background_mu(data["z_hd"], data["z_hel"], parameters, model)
            _, chi2, _ = profile_offset(base, data["mu_obs"], precision, p1, c0)
            if not np.isfinite(chi2):
                return 1e100
            return chi2
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return 1e100

    rng = np.random.default_rng(SEED + list(BOUNDS).index(model))
    midpoint = np.asarray([(lower + upper) / 2.0 for lower, upper in bounds])
    nominal = np.asarray([0.3, -1.0, 0.0][: len(bounds)], dtype=np.float64)
    initial = [midpoint, nominal]
    for _ in range(max(0, starts - len(initial))):
        initial.append(np.asarray([rng.uniform(lo, hi) for lo, hi in bounds]))
    initial = initial[:starts]

    results = []
    for start_index, start in enumerate(initial):
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000, "maxfun": 10000, "ftol": 1e-13, "gtol": 1e-8, "maxls": 40},
        )
        results.append(
            {
                "start_index": start_index,
                "success": bool(result.success),
                "message": str(result.message),
                "chi2": float(result.fun),
                "parameters": result.x.astype(float).tolist(),
                "nfev": int(result.nfev),
                "nit": int(result.nit),
            }
        )
    best = min(results, key=lambda item: item["chi2"])
    parameters = np.asarray(best["parameters"], dtype=np.float64)
    base = background_mu(data["z_hd"], data["z_hel"], parameters, model)
    offset, chi2, residual = profile_offset(
        base, data["mu_obs"], precision, p1, c0
    )
    difference = base - data["mu_obs"]
    independent = independent_profile_check(
        difference,
        data["precision_cholesky"],
        c0,
        offset,
        residual,
        chi2,
    )
    if independent["max_absolute_disagreement"] > 1e-6:
        raise ArithmeticError(f"independent profile algebra disagrees for {model}: {independent}")

    names = ["Omega_m"]
    if model != "LCDM":
        names.append("w0")
    if model == "CPL":
        names.append("wa")
    bound_hits = []
    for name, value, (lower, upper) in zip(names, parameters, bounds, strict=True):
        tolerance = 1e-6 * (1.0 + abs(float(value)))
        if abs(value - lower) <= tolerance:
            bound_hits.append({"parameter": name, "side": "lower", "bound": lower})
        if abs(value - upper) <= tolerance:
            bound_hits.append({"parameter": name, "side": "upper", "bound": upper})

    return {
        "model": model,
        "parameter_names": names,
        "parameters": {name: float(value) for name, value in zip(names, parameters, strict=True)},
        "parameter_bounds": {name: [float(lo), float(hi)] for name, (lo, hi) in zip(names, bounds, strict=True)},
        "profiled_additive_magnitude_offset_mag": offset,
        "chi2_profile": chi2,
        "n_data": int(len(data["mu_obs"])),
        "n_shape_parameters": len(parameters),
        "n_profiled_parameters": 1,
        "ndof_screening": int(len(data["mu_obs"]) - len(parameters) - 1),
        "chi2_per_ndof_screening": float(chi2 / (len(data["mu_obs"]) - len(parameters) - 1)),
        "bound_hits": bound_hits,
        "optimizer_best_start": best["start_index"],
        "optimizer_best_success": best["success"],
        "optimizer_best_message": best["message"],
        "optimizer_total_starts": starts,
        "optimizer_runs": results,
        "independent_algebra_check": independent,
        "residual_mag": residual,
        "mu_base_mag": base,
    }


def make_plot(data: dict[str, Any], fits: list[dict[str, Any]], output_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.2), sharex=True, constrained_layout=True)
    colors = {"LCDM": "#315a7d", "constant_w": "#c77c2c", "CPL": "#8d5a97"}
    for fit in fits:
        name = fit["model"]
        residual = fit["residual_mag"]
        axes[0].scatter(
            data["z_hd"], residual, s=7, alpha=0.20, color=colors[name],
            label=f"{name}: $\\chi^2_{{prof}}={fit['chi2_profile']:.1f}$",
            rasterized=True,
        )
        order = np.argsort(data["z_hd"])
        sorted_z = data["z_hd"][order]
        sorted_r = residual[order]
        bins = np.array_split(np.arange(len(sorted_z)), 32)
        binned_z = np.asarray([np.mean(sorted_z[chunk]) for chunk in bins])
        binned_r = np.asarray([np.median(sorted_r[chunk]) for chunk in bins])
        axes[0].plot(binned_z, binned_r, color=colors[name], linewidth=1.5)

    axes[0].axhline(0.0, color="black", linewidth=0.8, alpha=0.6)
    axes[0].set_ylabel("mu_base + Mhat - MU [mag]")
    axes[0].legend(loc="best", frameon=False, fontsize=8)
    axes[0].grid(alpha=0.15)

    lcdm = next(fit for fit in fits if fit["model"] == "LCDM")
    for fit in fits:
        if fit["model"] == "LCDM":
            continue
        delta = fit["mu_base_mag"] + fit["profiled_additive_magnitude_offset_mag"]
        delta -= lcdm["mu_base_mag"] + lcdm["profiled_additive_magnitude_offset_mag"]
        axes[1].scatter(data["z_hd"], delta, s=7, alpha=0.22, color=colors[fit["model"]], rasterized=True)
        order = np.argsort(data["z_hd"])
        bins = np.array_split(order, 32)
        binned_z = np.asarray([np.mean(data["z_hd"][chunk]) for chunk in bins])
        binned_delta = np.asarray([np.median(delta[chunk]) for chunk in bins])
        axes[1].plot(binned_z, binned_delta, color=colors[fit["model"]], linewidth=1.5, label=f"{fit['model']} minus LCDM")
    axes[1].axhline(0.0, color="black", linewidth=0.8, alpha=0.6)
    axes[1].set_xlabel("zHD")
    axes[1].set_ylabel("Profiled model difference from LCDM [mag]")
    axes[1].legend(loc="best", frameon=False, fontsize=8)
    axes[1].grid(alpha=0.15)
    fig.suptitle("DES-Dovekie flat-background profile screen\nvisual residual summaries; the fit uses the full correlated STAT+SYS precision")
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starts", type=int, default=12, help="deterministic multistarts per model")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    if args.starts < 2:
        raise ValueError("use at least two starts to check optimizer stability")

    command_started = time.perf_counter()
    data = load_inputs()
    fits = [fit_model(data, name, args.starts) for name in BOUNDS]
    fit_and_checks_elapsed = time.perf_counter() - command_started

    # At the minima, compare 96-point quadrature with an independent 160-point
    # Gauss-Legendre calculation and test arbitrary H0-gauge invariance.
    quadrature_checks = {}
    h0_gauge_check = {}
    for fit in fits:
        name = fit["model"]
        parameters = np.asarray([fit["parameters"][p] for p in fit["parameter_names"]])
        mu96 = background_mu(data["z_hd"], data["z_hel"], parameters, name, quadrature_order=96)
        mu160 = background_mu(data["z_hd"], data["z_hel"], parameters, name, quadrature_order=160)
        quadrature_checks[name] = float(np.max(np.abs(mu96 - mu160)))
        mu100 = background_mu(data["z_hd"], data["z_hel"], parameters, name, h0_gauge=100.0)
        off100, chi100, _ = profile_offset(
            mu100, data["mu_obs"], data["precision"], data["precision_ones"], data["offset_information"]
        )
        h0_gauge_check[name] = {
            "chi2_at_H0_70_gauge": fit["chi2_profile"],
            "chi2_at_H0_100_gauge": chi100,
            "chi2_abs_difference": abs(chi100 - fit["chi2_profile"]),
            "offset_shift_mag_H0_100_minus_70": off100 - fit["profiled_additive_magnitude_offset_mag"],
            "expected_offset_shift_mag": 5.0 * np.log10(100.0 / 70.0),
        }

    best_chi2 = min(fit["chi2_profile"] for fit in fits)
    for fit in fits:
        fit["delta_chi2_vs_LCDM"] = float(
            fit["chi2_profile"] - next(x["chi2_profile"] for x in fits if x["model"] == "LCDM")
        )
        fit["delta_chi2_vs_best_model"] = float(fit["chi2_profile"] - best_chi2)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = args.output_dir / "dovekie_hubble_residuals.png"
    make_plot(data, fits, plot_path)
    total_elapsed = time.perf_counter() - command_started
    code_hash = sha256(Path(__file__).resolve())
    count_by_survey: dict[str, int] = {}
    for record in data["records"]:
        survey = record["IDSURVEY"]
        count_by_survey[survey] = count_by_survey.get(survey, 0) + 1

    compact_fits = []
    for fit in fits:
        compact_fits.append({key: value for key, value in fit.items() if key not in {"residual_mag", "mu_base_mag"}})

    result = {
        "status": "exploratory_screening_calculation",
        "classification": "independent profile-likelihood implementation; not a published-posterior reproduction",
        "scientific_scope": {
            "dataset": "DES-Dovekie Hubble diagram only",
            "probe_combinations": [],
            "overlapping_supernova_compilations_used": [],
            "models": list(BOUNDS),
            "geometry": "flat, late-time matter plus dark energy",
            "radiation": "omitted for this z<=1.144 screening calculation",
            "curvature": "fixed to zero",
            "H0": "fixed to an arbitrary 70 km/s/Mpc distance-unit gauge and analytically absorbed by the free magnitude offset; not inferred or calibrated",
            "absolute_magnitude": "one unconstrained additive magnitude profiled analytically in every model",
            "likelihood_metric": "minimum profile chi-square using the complete released STAT+SYS inverse covariance",
            "muerr_handling": "MUERR and diagnostic MUERR_SYS are not added separately; released full STAT+SYS precision already supplies the Gaussian covariance",
        },
        "input_contract": {
            "repository": "https://github.com/des-science/DES-SN5YR",
            "repository_commit": PINNED_COMMIT,
            "hubble_diagram_url": f"https://raw.githubusercontent.com/des-science/DES-SN5YR/{PINNED_COMMIT}/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv",
            "precision_url": f"https://raw.githubusercontent.com/des-science/DES-SN5YR/{PINNED_COMMIT}/4_DISTANCES_COVMAT/STAT+SYS.npz",
            "hashes_sha256": data["hashes"],
            "hd_format": "SNANA VARNAMES/SN whitespace records despite .csv suffix",
            "official_reader_issue": "pinned Dovekie_cosmosis_likelihood.py calls Astropy Table.read(format='ascii.csv'), inconsistent with the pinned HD asset; this implementation uses the native official record format",
            "official_row_order": "unmodified file order, which the release README says is paired to the covariance; no sorting or metadata join",
            "official_filter": "retain zHD > 0.00, matching release code; all 1820 rows pass",
            "row_count": len(data["records"]),
            "unique_cid_count": len({row["CID"] for row in data["records"]}),
            "redshift": "zHD (CMB-frame with VPEC correction) enters E(z) integral; zHEL (heliocentric, no VPEC correction) enters the luminosity-distance prefactor",
            "calibration": "MU is used directly (released bias/contamination-corrected distance modulus assuming H0=70); no external absolute-magnitude or H0 calibration",
            "MUERR_range_mag_diagnostic_only": [float(data["muerr"].min()), float(data["muerr"].max())],
            "zHD_range": [float(data["z_hd"].min()), float(data["z_hd"].max())],
            "zHEL_range": [float(data["z_hel"].min()), float(data["z_hel"].max())],
            "IDSURVEY_counts": dict(sorted(count_by_survey.items())),
            "covariance_storage": "full STAT+SYS inverse covariance Covtot_inv, packed upper triangle, float32; unpacked to exactly symmetric float64 precision W",
            "archive_keys": data["archive_keys"],
            "precision_shape": list(data["precision"].shape),
            "packed_precision_entries": int(data["precision"].size // 2 + data["precision"].shape[0] // 2),
            "precision_cholesky": "pass",
            "offset_information_1TW1": data["offset_information"],
        },
        "likelihood_algebra": {
            "distance_modulus": "mu_base(zHD,zHEL)=5 log10[(1+zHEL)(c/H0_gauge) integral_0^zHD dz/E(z)] + 25",
            "chi2_before_profile": "(mu_base + M_offset - MU)^T W (mu_base + M_offset - MU)",
            "profile_solution": "Mhat = -1^T W d / (1^T W 1), d=mu_base-MU",
            "profile_chi2": "d^T W d - (1^T W d)^2/(1^T W 1), evaluated stably as r^T W r with r=d+Mhat*1",
            "released_likelihood_relation": "the release code adds log[(1^T W 1)/(2 pi)] for analytic integration over this same constant offset, plus a fixed Gaussian log determinant; both terms are model-independent for this fixed covariance, so profile Delta-chi-square is unchanged",
            "flat_expansion": {
                "LCDM": "E^2=Om(1+z)^3+(1-Om)",
                "constant_w": "E^2=Om(1+z)^3+(1-Om)(1+z)^(3(1+w))",
                "CPL": "E^2=Om(1+z)^3+(1-Om)(1+z)^(3(1+w0+wa))*exp[-3 wa z/(1+z)]",
            },
            "quadrature": f"{QUADRATURE_ORDER}-point Gauss-Legendre on each [0,zHD] interval; independent 160-point check at each reported optimum",
            "optimizer": "SciPy L-BFGS-B, 12 deterministic starts/model (midpoint, nominal nested point, seeded uniform starts), ftol=1e-13, gtol=1e-8, maxiter=2000, maxfun=10000",
            "parameter_domains": {model: {name: [lo, hi] for name, (lo, hi) in zip((["Omega_m"] + (["w0"] if model != "LCDM" else []) + (["wa"] if model == "CPL" else [])), box, strict=True)} for model, box in BOUNDS.items()},
            "domain_note": "bounds match the campaign's flat-background profile-screen domains; these are search domains, not posterior priors",
        },
        "checks": {
            "row_count_equals_precision_dimension": len(data["records"]) == data["precision"].shape[0],
            "all_rows_preserved_after_official_zHD_cut": True,
            "precision_cholesky_positive": True,
            "independent_cholesky_whitened_profile_algebra": {fit["model"]: fit["independent_algebra_check"] for fit in fits},
            "quadrature_max_abs_mu_difference_96_vs_160_mag": quadrature_checks,
            "arbitrary_H0_gauge_invariance": h0_gauge_check,
            "data_audit_condition_2_from_precision_eigenvalues": 59542388.186450124,
            "data_audit_precision_eigenvalue_range": [4.565498755832436e-06, 271.84069918452997],
        },
        "fits": compact_fits,
        "artifacts": {
            "code": str(Path(__file__).resolve().relative_to(ROOT)),
            "code_sha256": code_hash,
            "plot": str(plot_path.relative_to(ROOT)),
            "plot_sha256": sha256(plot_path),
            "result": str((args.output_dir / "result.json").relative_to(ROOT)),
        },
        "runtime": {
            "wall_seconds_fit_and_checks_excluding_plot_and_serialization": fit_and_checks_elapsed,
            "wall_seconds_through_plot_excluding_serialization": total_elapsed,
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "platform": platform.platform(),
            "logical_cpu_count_visible": os.cpu_count(),
            "device": "CPU double precision; no GPU calculation",
            "blas_thread_env": {name: os.environ.get(name) for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
            "random_seed": SEED,
            "command": f"OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python {Path(__file__).resolve().relative_to(ROOT)} --starts 12",
        },
        "run_history": [
            {
                "status": "failed_output_render",
                "command": f"OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python {Path(__file__).resolve().relative_to(ROOT)} --starts 12",
                "failure": "First invocation completed the three fits and numerical checks but Matplotlib rejected the plot's mathtext label (`\\mathcal`); it exited before writing result.json. The label was simplified and the full command rerun.",
            }
        ],
        "limits": [
            "This is an independent, reproducible profile comparison, not the official CosmoSIS pipeline run; the pinned official table reader has a format mismatch.",
            "No posterior sampling, Bayesian evidence, confidence intervals, or published-constraint reproduction is claimed.",
            "No calibrated H0 is measured; its distance normalization is exactly degenerate with the profiled magnitude offset in this radiation-free background screen.",
            "Radiation, curvature, perturbations, selection refits, and systematic nuisance parameters are not varied; the full released STAT+SYS covariance is retained.",
            "The three nested background families and bounded search were selected in advance; Delta-chi-square describes only these profile minima over recorded domains.",
        ],
    }

    output_path = args.output_dir / "result.json"
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "fits": compact_fits, "runtime": result["runtime"], "output": str(output_path)}, indent=2))


if __name__ == "__main__":
    main()

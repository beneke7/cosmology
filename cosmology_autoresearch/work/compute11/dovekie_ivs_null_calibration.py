#!/usr/bin/env python3
"""Finite-search Gaussian-null calibration for the Dovekie-only IVS profile.

The script independently reads the pinned DES-Dovekie Hubble diagram and
packed full STAT+SYS precision, generates mocks using P's Cholesky factor,
and repeats the audited finite LCDM/IVS searches. It uses one process pool
over realizations; each realization performs its grids and optimizations
serially. Run from the cosmology_autoresearch project root.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import time
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from scipy.stats import beta


ROOT = Path(__file__).resolve().parents[2]
HD_PATH = ROOT / "context/data/des_dovekie_hd.csv"
P_PATH = ROOT / "context/data/des_dovekie_stat_sys.npz"
PROFILE_RESULT_PATH = ROOT / "experiments/dovekie_ivs_profile/result.json"
PROFILE_SOURCE_PATH = ROOT / "experiments/dovekie_ivs_profile/profile.py"
INDEPENDENT_SOURCE_PATH = ROOT / "work/compute10/dovekie_ivs_profile_independent.py"
CONTRACT_PATH = ROOT / "work/theory7/dovekie_ivs_profile_contract.json"
OUT_DIR = ROOT / "work/compute11"
JSONL_PATH = OUT_DIR / "mock_results.jsonl"
RESULT_PATH = OUT_DIR / "result.json"
REPORT_PATH = OUT_DIR / "REPORT.md"

EXPECTED_HD_SHA = "2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf"
EXPECTED_P_SHA = "ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b"
EXPECTED_PROFILE_RESULT_SHA = "4068353752350739921f77c20754d54c85ae5024e0439280fa98e93656598817"
EXPECTED_PROFILE_SOURCE_SHA = "0dbc7acb0c4b5691e0d68144f86e82b2dcf0a1b17e8c07e8f49aa62a92a41f29"
EXPECTED_ROWS = 1820
EXPECTED_OMEGA_M = 0.3303167905483107
OBSERVED_T = 1.4213186471747576
MASTER_SEED_DEFAULT = 20260925
N_WORKERS = 20
PILOT_DEFAULT = 10
TARGET_DEFAULT = 200
CAP_DEFAULT_SECONDS = 1440.0
LCDM_STARTS = 8
IVS_STARTS = 16
GRID_N = 41
OM_BOUNDS = (0.05, 0.60)
G_BOUNDS = (-3.0, 3.0)
INVALID = 1.0e80
H0_GAUGE = 70.0
C_KM_S = 299792.458
U_MAX = math.log1p(1.14418)
RTOL = 2.0e-10
ATOL = 2.0e-12
MAX_STEP_U = 0.025
LCDM_GAUSS_N = 96
_GL_NODES, _GL_WEIGHTS = leggauss(LCDM_GAUSS_N)
_DATA: dict[str, Any] | None = None


class DomainError(ValueError):
    """Raised for nonphysical or numerically invalid IVS proposals."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_snana(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    columns: list[str] | None = None
    rows: list[dict[str, str]] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = raw.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "VARNAMES:":
            if columns is not None:
                raise ValueError(f"repeated VARNAMES at line {line_number}")
            columns = fields[1:]
        elif fields[0] == "SN:":
            if columns is None or len(fields[1:]) != len(columns):
                raise ValueError(f"malformed SN row at line {line_number}")
            rows.append(dict(zip(columns, fields[1:], strict=True)))
        else:
            raise ValueError(f"unexpected row at line {line_number}: {fields[0]!r}")
    if columns is None:
        raise ValueError("missing VARNAMES header")
    return columns, rows


def load_inputs() -> dict[str, Any]:
    hashes = {"hubble_diagram": sha256(HD_PATH), "packed_precision": sha256(P_PATH)}
    if hashes["hubble_diagram"] != EXPECTED_HD_SHA or hashes["packed_precision"] != EXPECTED_P_SHA:
        raise ValueError(f"pinned input hash mismatch: {hashes}")
    if sha256(PROFILE_RESULT_PATH) != EXPECTED_PROFILE_RESULT_SHA:
        raise ValueError("profile result identity changed from the independently reviewed version")
    if sha256(PROFILE_SOURCE_PATH) != EXPECTED_PROFILE_SOURCE_SHA:
        raise ValueError("profile source identity changed from the independently reviewed version")

    columns, rows = parse_snana(HD_PATH)
    if not {"CID", "zHD", "zHEL", "MU"}.issubset(columns):
        raise ValueError("pinned Hubble diagram lacks required columns")
    if len(rows) != EXPECTED_ROWS or len({row["CID"] for row in rows}) != EXPECTED_ROWS:
        raise ValueError("expected 1820 unique SN rows in release order")
    selected = [row for row in rows if float(row["zHD"]) > 0.0]
    if len(selected) != EXPECTED_ROWS:
        raise ValueError("official zHD > 0 cut changed the full covariance row count")
    z_hd = np.asarray([float(row["zHD"]) for row in selected], dtype=np.float64)
    z_hel = np.asarray([float(row["zHEL"]) for row in selected], dtype=np.float64)
    mu_obs = np.asarray([float(row["MU"]) for row in selected], dtype=np.float64)

    with np.load(P_PATH, allow_pickle=False) as archive:
        n = int(np.asarray(archive["nsn"]).reshape(-1)[0])
        packed = np.asarray(archive["cov"])
        archive_keys = list(archive.files)
    if n != EXPECTED_ROWS or packed.size != n * (n + 1) // 2 or packed.dtype != np.float32:
        raise ValueError("packed precision representation differs from the pinned contract")
    precision = np.zeros((n, n), dtype=np.float64)
    upper = np.triu_indices(n)
    precision[upper] = packed.astype(np.float64)
    precision[(upper[1], upper[0])] = precision[upper]
    if not np.array_equal(precision, precision.T) or not np.all(np.isfinite(precision)):
        raise ValueError("unpacked precision is not finite and exactly symmetric")
    precision_cholesky = np.linalg.cholesky(precision)
    one = np.ones(n, dtype=np.float64)
    precision_one = precision @ one
    q = float(one @ precision_one)
    if not np.isfinite(q) or q <= 0.0:
        raise ValueError("precision does not constrain the common intercept")

    profile_result = json.loads(PROFILE_RESULT_PATH.read_text(encoding="utf-8"))
    if profile_result.get("status") != "complete_exploratory_dovekie_only_ivs_profile":
        raise ValueError("source profile result is not complete")
    lcdm_ref = profile_result["results"]["dovekie_lcdm_reference"]
    omega_m = float(lcdm_ref["parameters"]["Omega_m0"])
    intercept = float(lcdm_ref["source_flat_profile_fit"]["profiled_additive_magnitude_offset_mag"])
    if omega_m != EXPECTED_OMEGA_M:
        raise ValueError(f"unexpected null Omega_m: {omega_m}")
    if abs(float(lcdm_ref["chi2"]) - 1631.4205355712409) > 1.0e-8:
        raise ValueError("source LCDM profile score changed")
    cid_hash = hashlib.sha256(("\n".join(row["CID"] for row in selected) + "\n").encode()).hexdigest()
    return {
        "z_hd": z_hd,
        "z_hel": z_hel,
        "mu_obs": mu_obs,
        "precision": precision,
        "precision_cholesky": precision_cholesky,
        "precision_one": precision_one,
        "offset_information": q,
        "omega_m_null": omega_m,
        "common_intercept": intercept,
        "source_lcdm_chi2": float(lcdm_ref["chi2"]),
        "cid_order_hash": cid_hash,
        "archive_keys": archive_keys,
        "hashes": hashes,
    }


def profile_score(mu_base: np.ndarray, mu_data: np.ndarray) -> tuple[float, float, float]:
    assert _DATA is not None
    difference = mu_base - mu_data
    precision_one = _DATA["precision_one"]
    q = float(_DATA["offset_information"])
    offset = -float(precision_one @ difference) / q
    residual = difference + offset
    chi2 = float(residual @ (_DATA["precision"] @ residual))
    normal = abs(float(precision_one @ residual))
    return offset, chi2, normal


def lcdm_mu(omega_m: float) -> np.ndarray:
    assert _DATA is not None
    z = _DATA["z_hd"]
    z_eval = (z[:, None] / 2.0) * (_GL_NODES[None, :] + 1.0)
    e2 = omega_m * (1.0 + z_eval) ** 3 + 1.0 - omega_m
    if np.any(~np.isfinite(e2)) or np.any(e2 <= 0.0):
        raise ValueError("LCDM expansion rate is invalid")
    integral = (z / 2.0) * np.sum(_GL_WEIGHTS[None, :] / np.sqrt(e2), axis=1)
    dl = (1.0 + _DATA["z_hel"]) * (C_KM_S / H0_GAUGE) * integral
    if np.any(~np.isfinite(dl)) or np.any(dl <= 0.0):
        raise ValueError("LCDM luminosity distance is invalid")
    return 5.0 * np.log10(dl) + 25.0


def integrate_ivs(omega_m: float, g: float) -> dict[str, Any]:
    assert _DATA is not None
    if not (math.isfinite(omega_m) and math.isfinite(g) and 0.0 < omega_m < 1.0):
        raise DomainError("nonfinite profile coordinate or invalid Omega_m")
    z_eval = _DATA["z_hd"]
    umax = math.log1p(float(np.max(z_eval)))

    def rhs(u: float, state: np.ndarray) -> tuple[float, float, float]:
        matter, vacuum, _distance = state
        e2 = matter + vacuum
        if not math.isfinite(e2) or e2 <= 0.0:
            raise DomainError("nonpositive/nonfinite E^2 in ODE RHS")
        e = math.sqrt(e2)
        source = g * vacuum / e
        return (3.0 * matter + source, -source, math.exp(u) / e)

    sol = solve_ivp(
        rhs,
        (0.0, umax),
        (omega_m, 1.0 - omega_m, 0.0),
        method="DOP853",
        dense_output=True,
        rtol=RTOL,
        atol=ATOL,
        max_step=MAX_STEP_U,
    )
    if not sol.success or sol.sol is None:
        raise DomainError(f"ODE solver failed: {sol.message}")

    n_check = max(257, int(math.ceil(umax / 0.002)) + 1)
    u_check = np.linspace(0.0, umax, n_check, dtype=np.float64)
    trajectory = np.asarray(sol.sol(u_check), dtype=np.float64)
    matter, vacuum, _distance = trajectory
    e2 = matter + vacuum
    if (
        np.any(~np.isfinite(trajectory))
        or np.any(matter <= 0.0)
        or np.any(vacuum <= 0.0)
        or np.any(e2 <= 0.0)
    ):
        raise DomainError("nonpositive/nonfinite m, x or E^2 on observed interval")
    if g >= 0.0:
        fb_max = 1.0
    else:
        fb_max = float(sol.sol(umax)[0]) / (omega_m * math.exp(3.0 * umax))
    if not math.isfinite(fb_max) or fb_max <= 0.0:
        raise DomainError("no positive constant baryon/CDM split exists")
    fb = 0.5 * fb_max
    cdm = matter - fb * omega_m * np.exp(3.0 * u_check)
    if np.any(cdm <= 0.0):
        raise DomainError("interior baryon/CDM witness is not positive")
    out = np.asarray(sol.sol(np.log1p(z_eval)), dtype=np.float64)
    if np.any(~np.isfinite(out)):
        raise DomainError("requested-redshift solution is invalid")
    distance = out[2]
    dl = (1.0 + _DATA["z_hel"]) * (C_KM_S / H0_GAUGE) * distance
    if np.any(~np.isfinite(dl)) or np.any(dl <= 0.0):
        raise DomainError("luminosity distance is invalid")
    mu_base = 5.0 * np.log10(dl) + 25.0
    return {
        "mu_base": mu_base,
        "min_matter": float(np.min(matter)),
        "min_vacuum": float(np.min(vacuum)),
        "min_e2": float(np.min(e2)),
        "fb_max": float(fb_max),
        "fb_witness": float(fb),
        "min_cdm_witness": float(np.min(cdm)),
    }


def score_ivs(point: tuple[float, float], mu_data: np.ndarray) -> tuple[float, str | None]:
    try:
        model = integrate_ivs(float(point[0]), float(point[1]))
        _offset, chi2, _normal = profile_score(model["mu_base"], mu_data)
        if not math.isfinite(chi2) or chi2 < -1.0e-7:
            raise DomainError("profiled chi-square is invalid")
        return float(chi2), None
    except (DomainError, ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as exc:
        return INVALID, type(exc).__name__ + ": " + str(exc)


def fit_lcdm(mu_data: np.ndarray) -> dict[str, Any]:
    bounds = [(OM_BOUNDS[0], OM_BOUNDS[1])]
    rng = np.random.default_rng(20260924)
    starts = [np.asarray([(OM_BOUNDS[0] + OM_BOUNDS[1]) / 2.0]), np.asarray([0.3])]
    for _ in range(LCDM_STARTS - len(starts)):
        starts.append(np.asarray([rng.uniform(*OM_BOUNDS)], dtype=np.float64))
    attempts: list[dict[str, Any]] = []

    def objective(theta: np.ndarray) -> float:
        try:
            mu = lcdm_mu(float(theta[0]))
            _offset, chi2, _normal = profile_score(mu, mu_data)
            return chi2 if math.isfinite(chi2) else 1.0e100
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return 1.0e100

    for index, start in enumerate(starts):
        fit = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000, "maxfun": 10000, "ftol": 1e-13, "gtol": 1e-8, "maxls": 40},
        )
        endpoint = objective(fit.x)
        attempts.append({
            "start_index": index,
            "start": start.astype(float).tolist(),
            "raw_success": bool(fit.success),
            "status": int(fit.status),
            "message": str(fit.message),
            "nfev": int(fit.nfev),
            "nit": int(fit.nit),
            "endpoint_omega_m": float(fit.x[0]),
            "reported_chi2": float(fit.fun),
            "independent_endpoint_recheck_chi2": float(endpoint),
            "endpoint_recheck_difference": float(abs(float(fit.fun) - endpoint)),
            "endpoint_valid": bool(math.isfinite(endpoint) and endpoint < 1.0e100),
        })
    valid = [row for row in attempts if row["endpoint_valid"]]
    if not valid:
        raise RuntimeError("all LCDM optimizer endpoints were invalid")
    best = min(valid, key=lambda row: row["independent_endpoint_recheck_chi2"])
    omega_m = float(best["endpoint_omega_m"])
    mu = lcdm_mu(omega_m)
    offset, chi2, normal = profile_score(mu, mu_data)
    return {
        "parameters": {"Omega_m0": omega_m},
        "chi2": chi2,
        "intercept": offset,
        "normal_equation_abs_residual": normal,
        "optimizer_best_start": best["start_index"],
        "attempts": attempts,
        "invalid_endpoint_count": len(attempts) - len(valid),
        "raw_unsuccessful_start_count": sum(not row["raw_success"] for row in attempts),
    }


def make_grid() -> list[tuple[float, float]]:
    oms = np.linspace(*OM_BOUNDS, GRID_N)
    gs = np.linspace(*G_BOUNDS, GRID_N)
    return [(float(om), float(g)) for g in gs for om in oms]


def build_ivs_starts(grid_minimum: tuple[float, float]) -> list[np.ndarray]:
    starts = [
        np.asarray(grid_minimum, dtype=np.float64),
        np.asarray([0.38697100769747017, -0.4666647299169934]),
        np.asarray([OM_BOUNDS[0], 0.0]),
        np.asarray([OM_BOUNDS[1], 0.0]),
        np.asarray([OM_BOUNDS[0], G_BOUNDS[0]]),
        np.asarray([OM_BOUNDS[0], G_BOUNDS[1]]),
        np.asarray([OM_BOUNDS[1], G_BOUNDS[0]]),
        np.asarray([OM_BOUNDS[1], G_BOUNDS[1]]),
    ]
    rng = np.random.default_rng(20260925)
    for _ in range(IVS_STARTS - len(starts)):
        starts.append(np.asarray([
            rng.uniform(*OM_BOUNDS), rng.uniform(*G_BOUNDS)
        ], dtype=np.float64))
    return [np.clip(row, [OM_BOUNDS[0], G_BOUNDS[0]], [OM_BOUNDS[1], G_BOUNDS[1]])
            for row in starts[:IVS_STARTS]]


def bound_hits(parameters: list[float]) -> list[dict[str, Any]]:
    hits = []
    for name, value, bounds in zip(("Omega_m0", "g"), parameters, (OM_BOUNDS, G_BOUNDS), strict=True):
        tol = 2e-5 * (bounds[1] - bounds[0])
        if abs(value - bounds[0]) <= tol:
            hits.append({"parameter": name, "side": "lower", "bound": bounds[0]})
        if abs(value - bounds[1]) <= tol:
            hits.append({"parameter": name, "side": "upper", "bound": bounds[1]})
    return hits


def fit_ivs(mu_data: np.ndarray) -> dict[str, Any]:
    grid = make_grid()
    grid_records = [score_ivs(point, mu_data) for point in grid]
    grid_invalid = Counter(reason for _score, reason in grid_records if reason is not None)
    grid_values = np.asarray([score for score, _reason in grid_records], dtype=np.float64).reshape(GRID_N, GRID_N)
    finite = np.isfinite(grid_values) & (grid_values < INVALID / 2.0)
    if not np.any(finite):
        raise RuntimeError("IVS coarse grid has no valid physical point")
    flat = int(np.argmin(np.where(finite, grid_values, np.inf)))
    g_index, om_index = np.unravel_index(flat, grid_values.shape)
    grid_minimum = (
        float(np.linspace(*OM_BOUNDS, GRID_N)[om_index]),
        float(np.linspace(*G_BOUNDS, GRID_N)[g_index]),
    )
    attempts: list[dict[str, Any]] = []
    for index, start in enumerate(build_ivs_starts(grid_minimum)):
        failures: Counter[str] = Counter()

        def objective(theta: np.ndarray) -> float:
            value, reason = score_ivs((float(theta[0]), float(theta[1])), mu_data)
            if reason is not None:
                failures[reason] += 1
            return value

        fit = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=[OM_BOUNDS, G_BOUNDS],
            options={"maxiter": 2000, "maxfun": 10000, "ftol": 1e-13, "gtol": 1e-8, "maxls": 40},
        )
        endpoint = score_ivs((float(fit.x[0]), float(fit.x[1])), mu_data)
        recheck_difference = abs(float(fit.fun) - endpoint[0])
        if endpoint[1] is None and recheck_difference > 1.0e-6:
            raise ArithmeticError(f"IVS optimizer endpoint recheck differs by {recheck_difference}")
        attempts.append({
            "start_index": index,
            "start": start.astype(float).tolist(),
            "raw_success": bool(fit.success),
            "status": int(fit.status),
            "message": str(fit.message),
            "nfev": int(fit.nfev),
            "nit": int(fit.nit),
            "endpoint_parameters": fit.x.astype(float).tolist(),
            "reported_chi2": float(fit.fun),
            "independent_endpoint_recheck_chi2": float(endpoint[0]),
            "endpoint_recheck_difference": float(recheck_difference),
            "endpoint_valid_physical": bool(endpoint[1] is None and endpoint[0] < INVALID),
            "objective_rejection_counts": dict(failures),
        })
    valid = [row for row in attempts if row["endpoint_valid_physical"]]
    if not valid:
        raise RuntimeError("all IVS optimizer endpoints were invalid")
    best = min(valid, key=lambda row: row["independent_endpoint_recheck_chi2"])
    om, g = map(float, best["endpoint_parameters"])
    point = integrate_ivs(om, g)
    offset, chi2, normal = profile_score(point["mu_base"], mu_data)
    return {
        "parameters": {"Omega_m0": om, "g": g},
        "chi2": chi2,
        "intercept": offset,
        "normal_equation_abs_residual": normal,
        "grid_minimum": list(grid_minimum),
        "grid_minimum_chi2": float(grid_values[g_index, om_index]),
        "grid_invalid_count": int(sum(grid_invalid.values())),
        "grid_invalid_reasons": dict(grid_invalid),
        "optimizer_attempts": attempts,
        "invalid_endpoint_count": len(attempts) - len(valid),
        "raw_unsuccessful_start_count": sum(not row["raw_success"] for row in attempts),
        "invalid_optimizer_objective_evaluations": int(sum(
            sum(row["objective_rejection_counts"].values()) for row in attempts
        )),
        "boundary_hits": bound_hits(best["endpoint_parameters"]),
        "physicality": {
            "minimum_matter": point["min_matter"],
            "minimum_vacuum": point["min_vacuum"],
            "minimum_e2": point["min_e2"],
            "fb_max_exact": point["fb_max"],
            "fb_witness": point["fb_witness"],
            "minimum_cdm_witness": point["min_cdm_witness"],
        },
    }


def run_mock(replicate_id: int, seed: int) -> dict[str, Any]:
    assert _DATA is not None
    start = time.perf_counter()
    try:
        rng = np.random.default_rng(seed)
        z = rng.standard_normal(EXPECTED_ROWS)
        # If P=L L^T, then epsilon=L^{-T} z has covariance
        # L^{-T} L^{-1}=P^{-1}=C.
        noise = solve_triangular(
            _DATA["precision_cholesky"].T,
            z,
            lower=False,
            check_finite=False,
            overwrite_b=False,
        )
        null_mean = _DATA["null_mu_base"] + _DATA["common_intercept"]
        mock_mu = null_mean + noise
        lcdm = fit_lcdm(mock_mu)
        ivs = fit_ivs(mock_mu)
        statistic = float(lcdm["chi2"] - ivs["chi2"])
        return {
            "replicate_id": replicate_id,
            "seed": seed,
            "status": "complete",
            "elapsed_seconds": float(time.perf_counter() - start),
            "sample_generation": {
                "noise_finite": bool(np.all(np.isfinite(noise))),
                "noise_norm2": float(noise @ noise),
                "whitened_recovery_max_abs": float(np.max(np.abs(
                    _DATA["precision_cholesky"].T @ noise - z
                ))),
            },
            "lcdm": lcdm,
            "ivs": ivs,
            "T_chi2_lcdm_minus_chi2_ivs": statistic,
            "exceeds_observed_T_inclusive": bool(statistic >= OBSERVED_T),
        }
    except Exception as exc:  # preserve failed realizations as explicit JSONL rows
        return {
            "replicate_id": replicate_id,
            "seed": seed,
            "status": "failed",
            "elapsed_seconds": float(time.perf_counter() - start),
            "error_type": type(exc).__name__,
            "error": str(exc),
            "T_chi2_lcdm_minus_chi2_ivs": None,
            "exceeds_observed_T_inclusive": None,
        }


def verify_sampling_algebra() -> dict[str, Any]:
    """Check the Cholesky sampling identity on a synthetic SPD matrix."""
    covariance = np.asarray([
        [1.3, 0.2, -0.1],
        [0.2, 0.9, 0.15],
        [-0.1, 0.15, 0.7],
    ], dtype=np.float64)
    precision = np.linalg.inv(covariance)
    lower = np.linalg.cholesky(precision)
    transform = solve_triangular(lower.T, np.eye(3), lower=False)
    generated_covariance = transform @ transform.T
    seed_z = np.asarray([0.4, -1.2, 2.1], dtype=np.float64)
    sample = solve_triangular(lower.T, seed_z, lower=False)
    covariance_error = float(np.max(np.abs(generated_covariance - covariance)))
    forward_error = float(np.max(np.abs(lower.T @ sample - seed_z)))
    if covariance_error > 2.0e-14 or forward_error > 2.0e-14:
        raise ArithmeticError("synthetic Cholesky sampling algebra check failed")
    return {
        "identity": "P=L L^T; epsilon=solve(L^T,z); Cov(epsilon)=L^{-T}L^{-1}=P^{-1}=C",
        "synthetic_spd_dimension": 3,
        "max_abs_generated_vs_target_covariance": covariance_error,
        "max_abs_whitened_recovery_error": forward_error,
        "passed": True,
    }


def append_row(stream: Any, row: dict[str, Any]) -> None:
    stream.write(json.dumps(row, allow_nan=False, separators=(",", ":")) + "\n")
    stream.flush()


def cp_interval(successes: int, total: int, confidence: float = 0.95) -> list[float] | None:
    if total <= 0:
        return None
    alpha = 1.0 - confidence
    lower = 0.0 if successes == 0 else float(beta.ppf(alpha / 2.0, successes, total - successes + 1))
    upper = 1.0 if successes == total else float(beta.ppf(1.0 - alpha / 2.0, successes + 1, total - successes))
    return [lower, upper]


def write_report(result: dict[str, Any]) -> None:
    interval = result["tail_probability_estimate"]["clopper_pearson_95_interval"]
    if interval is None:
        interval_text = "not defined (no valid mock fits)"
        fraction_text = "not defined"
    else:
        interval_text = f"[{interval[0]:.6g}, {interval[1]:.6g}]"
        fraction_text = (
            f"{result['tail_probability_estimate']['exceedance_count']}/"
            f"{result['tail_probability_estimate']['valid_mock_count']} = "
            f"{result['tail_probability_estimate']['empirical_fraction']:.6g}"
        )
    projection = result["pilot_projection"]
    lines = [
        "# Dovekie IVS finite-search null calibration",
        "",
        f"Status: `{result['status']}`. This is an empirical finite-search calibration under the fitted Dovekie LCDM shape; it is not a posterior p-value or discovery significance.",
        "",
        f"Generated {result['completed_realizations']} of {result['target_realizations']} target realizations with master seed `{result['master_seed']}`. Each mock used the 1,820-row released STAT+SYS precision and one common fitted magnitude intercept. The realized noise was generated as `solve(L.T, z)` for `P=L L.T`, which gives covariance `C=P^-1`.",
        "",
        f"Each valid mock refit LCDM with 8 starts and the same IVS finite search: Ωm ∈ [0.05, 0.60], g ∈ [-3, 3], a 41×41 grid, 16 L-BFGS-B starts including the grid minimum, the positive-history guard, invalid penalty, and endpoint rechecks. One 20-process pool ran realizations; searches within each realization ran serially, with BLAS threads set to one.",
        "",
        f"The count with `T = χ²_LCDM − χ²_IVS ≥ {OBSERVED_T:.16g}` was {fraction_text}. Exact two-sided 95% Clopper–Pearson interval: {interval_text}. The denominator includes completed mock fits; failed fits are counted and reported separately.",
        "",
        f"The 10-mock pilot took {projection['pilot_wall_seconds']:.3f} s wall time. Based on its measured per-mock runtime and 20 workers, the projected total with a 25% throughput allowance and 30 s finalization reserve was {projection['projected_total_seconds_with_safety']:.1f} s, against a {projection['command_cap_seconds']:.0f} s internal cap. Projection decision: `{projection['decision']}`.",
        "",
        f"Failed realization fits: {result['fit_diagnostics']['failed_realization_count']}; invalid IVS grid points summed over complete mocks: {result['fit_diagnostics']['invalid_grid_point_count_total']}; invalid IVS optimizer objective evaluations: {result['fit_diagnostics']['invalid_optimizer_objective_evaluation_count_total']}; IVS invalid endpoints: {result['fit_diagnostics']['invalid_ivs_endpoint_count_total']}; LCDM invalid endpoints: {result['fit_diagnostics']['invalid_lcdm_endpoint_count_total']}; raw unsuccessful starts (LCDM/IVS): {result['fit_diagnostics']['raw_unsuccessful_lcdm_start_count_total']}/{result['fit_diagnostics']['raw_unsuccessful_ivs_start_count_total']}.",
        "",
        f"Pinned input and code hashes, algebra check, full run command, per-realization statuses/statistics, and environment are recorded in [result.json](result.json) and [mock_results.jsonl](mock_results.jsonl). The generator, scoring and search are in [dovekie_ivs_null_calibration.py](dovekie_ivs_null_calibration.py).",
        "",
    ]
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    complete = [row for row in rows if row.get("status") == "complete"]
    failed = [row for row in rows if row.get("status") != "complete"]
    exceed = sum(bool(row["exceeds_observed_T_inclusive"]) for row in complete)
    lcdm_invalid_endpoints = sum(row["lcdm"]["invalid_endpoint_count"] for row in complete)
    ivs_invalid_endpoints = sum(row["ivs"]["invalid_endpoint_count"] for row in complete)
    lcdm_unsuccessful = sum(row["lcdm"]["raw_unsuccessful_start_count"] for row in complete)
    ivs_unsuccessful = sum(row["ivs"]["raw_unsuccessful_start_count"] for row in complete)
    return {
        "completed_mock_fit_count": len(complete),
        "failed_realization_count": len(failed),
        "exceedance_count": exceed,
        "empirical_fraction": (exceed / len(complete)) if complete else None,
        "clopper_pearson_95_interval": cp_interval(exceed, len(complete)),
        "invalid_grid_point_count_total": sum(row["ivs"]["grid_invalid_count"] for row in complete),
        "invalid_optimizer_objective_evaluation_count_total": sum(
            row["ivs"]["invalid_optimizer_objective_evaluations"] for row in complete
        ),
        "invalid_ivs_endpoint_count_total": ivs_invalid_endpoints,
        "invalid_lcdm_endpoint_count_total": lcdm_invalid_endpoints,
        "raw_unsuccessful_ivs_start_count_total": ivs_unsuccessful,
        "raw_unsuccessful_lcdm_start_count_total": lcdm_unsuccessful,
        "failed_realization_ids": [row["replicate_id"] for row in failed],
    }


def main() -> None:
    global _DATA
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=int, default=TARGET_DEFAULT)
    parser.add_argument("--pilot", type=int, default=PILOT_DEFAULT)
    parser.add_argument("--master-seed", type=int, default=MASTER_SEED_DEFAULT)
    parser.add_argument("--cap-seconds", type=float, default=CAP_DEFAULT_SECONDS)
    args = parser.parse_args()
    if args.target < 1 or args.pilot < 1 or args.pilot >= args.target:
        parser.error("require target > pilot >= 1")
    if args.cap_seconds <= 60.0 or not math.isfinite(args.cap_seconds):
        parser.error("cap-seconds must be finite and exceed 60 seconds")

    start_wall = time.perf_counter()
    deadline = start_wall + args.cap_seconds
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if JSONL_PATH.exists() or RESULT_PATH.exists() or REPORT_PATH.exists():
        raise FileExistsError("work/compute11 outputs already exist; preserving them and refusing overwrite")
    _DATA = load_inputs()
    _DATA["null_mu_base"] = lcdm_mu(_DATA["omega_m_null"])
    algebra_check = verify_sampling_algebra()
    affinity_count = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    if affinity_count < N_WORKERS:
        raise RuntimeError(f"20 workers requested but only {affinity_count} CPUs are available to this process")

    seeds = [
        int(child.generate_state(1, dtype=np.uint64)[0])
        for child in np.random.SeedSequence(args.master_seed).spawn(args.target)
    ]
    rows: list[dict[str, Any]] = []
    pilot_wall_seconds = None
    projected_total = None
    projection_decision = "not_calculated"
    command = (
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 "
        "PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --seconds 1470 -- "
        ".venv/bin/python work/compute11/dovekie_ivs_null_calibration.py "
        f"--target {args.target} --pilot {args.pilot} --master-seed {args.master_seed} "
        f"--cap-seconds {args.cap_seconds:g}"
    )

    with JSONL_PATH.open("x", encoding="utf-8") as out, ProcessPoolExecutor(
        max_workers=N_WORKERS,
        mp_context=mp.get_context("fork"),
    ) as pool:
        pilot_start = time.perf_counter()
        pilot_futures = {
            pool.submit(run_mock, replicate_id, seeds[replicate_id]): replicate_id
            for replicate_id in range(args.pilot)
        }
        for future in as_completed(pilot_futures):
            replicate_id = pilot_futures[future]
            try:
                row = future.result()
            except Exception as exc:
                row = {
                    "replicate_id": replicate_id,
                    "seed": seeds[replicate_id],
                    "status": "failed_worker_exception",
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "T_chi2_lcdm_minus_chi2_ivs": None,
                    "exceeds_observed_T_inclusive": None,
                }
            rows.append(row)
            append_row(out, row)
            print(f"completed realization {len(rows)}/{args.target}: id={replicate_id}, status={row['status']}", flush=True)
        pilot_wall_seconds = time.perf_counter() - pilot_start
        valid_pilot = [row for row in rows if row.get("status") == "complete"]
        if valid_pilot:
            mean_task = float(np.mean([row["elapsed_seconds"] for row in valid_pilot]))
            remaining = args.target - len(rows)
            projected_total = (
                time.perf_counter() - start_wall
                + remaining * mean_task / N_WORKERS * 1.25
                + 30.0
            )
            if projected_total <= args.cap_seconds:
                projection_decision = "continue_to_target"
            else:
                projection_decision = "pilot_only_projection_exceeds_cap"
        else:
            mean_task = None
            projection_decision = "pilot_failed_without_runtime_projection"

        print(json.dumps({
            "pilot_realizations": len(rows),
            "pilot_wall_seconds": pilot_wall_seconds,
            "mean_valid_mock_seconds": mean_task,
            "projected_total_seconds_with_safety": projected_total,
            "cap_seconds": args.cap_seconds,
            "decision": projection_decision,
        }, allow_nan=False), flush=True)

        next_id = len(rows)
        while (
            projection_decision == "continue_to_target"
            and next_id < args.target
        ):
            if time.perf_counter() >= deadline - 15.0:
                projection_decision = "stopped_at_cap_reserve"
                break
            batch_ids = list(range(next_id, min(next_id + N_WORKERS, args.target)))
            batch_futures = {pool.submit(run_mock, i, seeds[i]): i for i in batch_ids}
            for future in as_completed(batch_futures):
                replicate_id = batch_futures[future]
                try:
                    row = future.result()
                except Exception as exc:
                    row = {
                        "replicate_id": replicate_id,
                        "seed": seeds[replicate_id],
                        "status": "failed_worker_exception",
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                        "T_chi2_lcdm_minus_chi2_ivs": None,
                        "exceeds_observed_T_inclusive": None,
                    }
                rows.append(row)
                append_row(out, row)
                print(f"completed realization {len(rows)}/{args.target}: id={replicate_id}, status={row['status']}", flush=True)
            next_id += len(batch_ids)
            valid_rows = [row for row in rows if row.get("status") == "complete"]
            if valid_rows:
                mean_task = float(np.mean([row["elapsed_seconds"] for row in valid_rows]))
                remaining = args.target - len(rows)
                projected_total = (
                    time.perf_counter() - start_wall
                    + remaining * mean_task / N_WORKERS * 1.25
                    + 30.0
                )
                if projected_total > args.cap_seconds or projected_total > deadline - time.perf_counter():
                    projection_decision = "stopped_after_updated_projection_exceeds_cap"
                    break
            if time.perf_counter() >= deadline - 15.0 and next_id < args.target:
                projection_decision = "stopped_at_cap_reserve"
                break

    completed_count = len(rows)
    if completed_count == args.target and all(row.get("status") == "complete" for row in rows):
        status = "complete_200"
    elif completed_count == args.target:
        status = "target_realizations_with_fit_failures"
    else:
        status = "partial_within_command_cap"

    summarize_result = summarize(rows)
    source_hashes = {
        "hubble_diagram_sha256": _DATA["hashes"]["hubble_diagram"] if _DATA else EXPECTED_HD_SHA,
        "packed_stat_sys_precision_sha256": EXPECTED_P_SHA,
        "profile_result_sha256": EXPECTED_PROFILE_RESULT_SHA,
        "profile_source_sha256": EXPECTED_PROFILE_SOURCE_SHA,
        "independent_reproduction_script_sha256": sha256(INDEPENDENT_SOURCE_PATH),
        "profile_contract_sha256": sha256(CONTRACT_PATH),
        "calibration_script_sha256": sha256(Path(__file__).resolve()),
    }
    projection = {
        "pilot_realizations": args.pilot,
        "pilot_wall_seconds": pilot_wall_seconds,
        "mean_valid_mock_fit_seconds": mean_task,
        "projected_total_seconds_with_safety": projected_total,
        "command_cap_seconds": args.cap_seconds,
        "safety_factor_on_remaining_work": 1.25,
        "finalization_reserve_seconds": 30.0,
        "decision": projection_decision,
    }
    result = {
        "status": status,
        "objective": "Empirical finite-search ΛCDM-null calibration for the independently reproduced Dovekie-only IVS profile; not a posterior p-value or discovery significance.",
        "null_model": {
            "Omega_m0": EXPECTED_OMEGA_M,
            "g": 0.0,
            "common_intercept_mag": _DATA["common_intercept"] if _DATA else None,
            "common_intercept_role": "one shared additive magnitude included in every mock and profiled analytically in both fits",
            "true_model_chi2_profile": _DATA["source_lcdm_chi2"] if _DATA else None,
        },
        "statistic": {
            "definition": "T=chi2_LCDM-chi2_IVS",
            "inclusive_threshold": OBSERVED_T,
            "source_observed_delta": "work/compute10/dovekie_ivs_profile_independent_review.json",
        },
        "master_seed": args.master_seed,
        "target_realizations": args.target,
        "completed_realizations": completed_count,
        "seed_rule": "SeedSequence(master_seed).spawn(target); one uint64 generated per child and used with default_rng",
        "seed_values_by_replicate_id": seeds,
        "tail_probability_estimate": {
            "exceedance_count": summarize_result["exceedance_count"],
            "valid_mock_count": summarize_result["completed_mock_fit_count"],
            "empirical_fraction": summarize_result["empirical_fraction"],
            "confidence_level": 0.95,
            "interval_method": "two-sided exact Clopper-Pearson binomial interval",
            "clopper_pearson_95_interval": summarize_result["clopper_pearson_95_interval"],
        },
        "fit_diagnostics": {
            "failed_realization_count": summarize_result["failed_realization_count"],
            "failed_realization_ids": summarize_result["failed_realization_ids"],
            "invalid_grid_point_count_total": summarize_result["invalid_grid_point_count_total"],
            "invalid_optimizer_objective_evaluation_count_total": summarize_result["invalid_optimizer_objective_evaluation_count_total"],
            "invalid_ivs_endpoint_count_total": summarize_result["invalid_ivs_endpoint_count_total"],
            "invalid_lcdm_endpoint_count_total": summarize_result["invalid_lcdm_endpoint_count_total"],
            "raw_unsuccessful_ivs_start_count_total": summarize_result["raw_unsuccessful_ivs_start_count_total"],
            "raw_unsuccessful_lcdm_start_count_total": summarize_result["raw_unsuccessful_lcdm_start_count_total"],
        },
        "sampling_algebra_check": algebra_check,
        "search_contract": {
            "lcdm_shape_bounds": [OM_BOUNDS[0], OM_BOUNDS[1]],
            "lcdm_starts": LCDM_STARTS,
            "ivs_shape_bounds": {"Omega_m0": list(OM_BOUNDS), "g": list(G_BOUNDS)},
            "ivs_grid_shape": [GRID_N, GRID_N],
            "ivs_start_count": IVS_STARTS,
            "grid_minimum_included_as_start": True,
            "invalid_score_penalty": INVALID,
            "optimizer": "L-BFGS-B; maxiter=2000,maxfun=10000,ftol=1e-13,gtol=1e-8,maxls=40",
            "endpoint_independent_recheck": True,
            "positive_history_guard": "DOP853; rtol=2e-10,atol=2e-12,max_step_u=0.025; dense interval sample spacing <=0.002 in u plus positive interior baryon/CDM split witness",
        },
        "pilot_projection": projection,
        "inputs": {
            "rows": EXPECTED_ROWS,
            "precision": "full released STAT+SYS P=C^-1, unpacked float32 packed upper triangle into float64 in release order",
            "cid_order_sha256": "b7d5c6ad8dfdadf1e12443852b14006c0d9d0ceb4bcb46efb378f058ce5aa5a3",
            "zHD_range": [float(np.min(_DATA["z_hd"])) if _DATA else 0.02509, float(np.max(_DATA["z_hd"])) if _DATA else 1.14418],
            "zHEL_range": [float(np.min(_DATA["z_hel"])) if _DATA else 0.02385, float(np.max(_DATA["z_hel"])) if _DATA else 1.14501],
            "muerr_added_separately": False,
            "archive_keys": ["nsn", "cov", "allow_pickle"],
        },
        "hashes": source_hashes,
        "runtime_seconds": float(time.perf_counter() - start_wall),
        "environment": {
            "workers": N_WORKERS,
            "pool_scope": "across realizations only; no nested pools",
            "blas_threads_requested": 1,
            "cpu_count_visible": affinity_count,
            "gpu_used": False,
        },
        "command": command,
        "artifacts": {
            "jsonl": "work/compute11/mock_results.jsonl",
            "result": "work/compute11/result.json",
            "report": "work/compute11/REPORT.md",
        },
    }
    # Verify the pinned source artifacts remained unchanged during computation.
    if sha256(PROFILE_RESULT_PATH) != EXPECTED_PROFILE_RESULT_SHA or sha256(PROFILE_SOURCE_PATH) != EXPECTED_PROFILE_SOURCE_SHA:
        raise RuntimeError("source profile artifact changed during calibration")
    RESULT_PATH.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_report(result)
    print(json.dumps({
        "status": result["status"],
        "completed": result["completed_realizations"],
        "tail_count": result["tail_probability_estimate"]["exceedance_count"],
        "tail_fraction": result["tail_probability_estimate"]["empirical_fraction"],
        "clopper_pearson_95": result["tail_probability_estimate"]["clopper_pearson_95_interval"],
        "projection_decision": projection_decision,
        "runtime_seconds": result["runtime_seconds"],
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

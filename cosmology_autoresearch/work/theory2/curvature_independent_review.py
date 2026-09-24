#!/usr/bin/env python3
"""Independent scalar-quadrature curved-LCDM BAO profile cross-check.

No project prediction or likelihood implementation is imported.  Distances
are evaluated row by row with scipy.integrate.quad, then scored against the
full released covariance in the released mean-vector order.
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
from scipy.integrate import quad
from scipy.optimize import minimize, minimize_scalar


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().with_suffix(".json")
MEAN_PATH = ROOT / "context/data/desi_dr2_mean.txt"
COV_PATH = ROOT / "context/data/desi_dr2_cov.txt"
PAPER_PATH = ROOT / "context/papers/desi_dr2_bao.pdf"
STARTER_PATH = ROOT / "scripts/background_bao.py"
SCREEN_PATH = ROOT / "experiments/curvature_screen/result.json"

ALPHA_BOUNDS = (1.0e-6, 1.0e4)
OMEGA_M_BOUNDS = (0.05, 0.60)
OMEGA_K_LIMITS = (0.05, 0.10, 0.20)
Z_MAX = 2.33
FIT_EPSABS = 1.0e-12
FIT_EPSREL = 1.0e-12
TIGHT_EPSABS = 1.0e-13
TIGHT_EPSREL = 1.0e-13


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_mean(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    redshifts: list[float] = []
    values: list[float] = []
    observables: list[str] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"{path}:{line_number}: expected z value observable")
        redshifts.append(float(fields[0]))
        values.append(float(fields[1]))
        observables.append(fields[2])
    return (np.asarray(redshifts, dtype=np.float64),
            np.asarray(values, dtype=np.float64),
            np.asarray(observables, dtype=str))


def expansion_squared(z: float, omega_m: float, omega_k: float) -> float:
    """Curved LCDM E^2 with flat-closure Omega_Lambda."""
    x = 1.0 + z
    omega_lambda = 1.0 - omega_m - omega_k
    e2 = omega_m * x**3 + omega_k * x**2 + omega_lambda
    if not math.isfinite(e2) or e2 <= 0.0:
        raise ValueError(f"E^2 is not positive at z={z:g}: {e2!r}")
    return e2


def transverse_kernel(chi: float, omega_k: float) -> float:
    """Dimensionless S_k kernel with an analytic near-flat series."""
    y = omega_k * chi**2
    if abs(y) < 1.0e-4:
        return chi * (1.0 + y / 6.0 + y**2 / 120.0
                      + y**3 / 5040.0 + y**4 / 362880.0)
    if y > 0.0:
        root = math.sqrt(y)
        return math.sinh(root * 1.0) * chi / root
    root = math.sqrt(-y)
    return math.sin(root) * chi / root


def unit_amplitude_predictions(
    z: np.ndarray,
    observables: np.ndarray,
    omega_m: float,
    omega_k: float,
    *,
    epsabs: float = FIT_EPSABS,
    epsrel: float = FIT_EPSREL,
    return_quad_errors: bool = False,
) -> np.ndarray | tuple[np.ndarray, dict[str, float]]:
    """Return q with p=alpha*q, using scalar adaptive integration per row."""
    q = np.empty(len(z), dtype=np.float64)
    errors: dict[str, float] = {}
    omega_lambda = 1.0 - omega_m - omega_k
    if omega_lambda < 0.0:
        raise ValueError("profile domain requires nonnegative Omega_Lambda")
    for index, (redshift, observable) in enumerate(zip(z, observables, strict=True)):
        zi = float(redshift)
        label = str(observable)
        e2_z = expansion_squared(zi, omega_m, omega_k)
        e_z = math.sqrt(e2_z)
        if label == "DH_over_rs":
            q[index] = 1.0 / e_z
            continue

        def inverse_e(zp: float) -> float:
            return 1.0 / math.sqrt(expansion_squared(zp, omega_m, omega_k))

        chi, error = quad(inverse_e, 0.0, zi, epsabs=epsabs, epsrel=epsrel, limit=100)
        s_k = transverse_kernel(float(chi), omega_k)
        errors[f"row_{index}_chi_abs_error_estimate"] = float(error)
        if label == "DM_over_rs":
            q[index] = s_k
        elif label == "DV_over_rs":
            q[index] = (zi * s_k**2 / e_z) ** (1.0 / 3.0)
        else:
            raise ValueError(f"unknown observable {label!r}")
    if return_quad_errors:
        return q, errors
    return q


class FullCovarianceBAO:
    def __init__(self) -> None:
        self.z, self.y, self.observable = load_mean(MEAN_PATH)
        self.covariance = np.loadtxt(COV_PATH, dtype=np.float64)
        if self.covariance.shape != (len(self.y), len(self.y)):
            raise ValueError("covariance shape does not match the mean vector")
        self.cholesky = np.linalg.cholesky(self.covariance)
        self.cinv_y = np.linalg.solve(self.covariance, self.y)
        self.y_cinv_y = float(self.y @ self.cinv_y)

    def score(self, prediction: np.ndarray) -> float:
        residual = self.y - prediction
        return float(residual @ np.linalg.solve(self.covariance, residual))

    def profile_alpha(self, q: np.ndarray) -> dict[str, object]:
        """GLS profile for p=alpha*q; clip only to the existing alpha box."""
        cinv_q = np.linalg.solve(self.covariance, q)
        numerator = float(q @ self.cinv_y)
        denominator = float(q @ cinv_q)
        alpha_unbounded = numerator / denominator
        alpha = float(np.clip(alpha_unbounded, *ALPHA_BOUNDS))
        prediction = alpha * q
        return {
            "alpha_gls_unbounded": alpha_unbounded,
            "alpha_profiled": alpha,
            "alpha_is_box_clipped": alpha != alpha_unbounded,
            "chi2": self.score(prediction),
            "prediction": prediction,
            "quadratic_denominator": denominator,
            "gls_closed_form_chi2_unbounded": self.y_cinv_y - numerator**2 / denominator,
        }


DATA = FullCovarianceBAO()


def hit_bounds(parameters: dict[str, float], bounds: dict[str, tuple[float, float]]) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    for name, value in parameters.items():
        if name not in bounds:
            continue
        lo, hi = bounds[name]
        tolerance = max(1.0e-8, 1.0e-6 * (hi - lo))
        if abs(value - lo) <= tolerance:
            hits.append({"parameter": name, "side": "lower", "value": value})
        elif abs(value - hi) <= tolerance:
            hits.append({"parameter": name, "side": "upper", "value": value})
    return hits


def profile_at_shape(omega_m: float, omega_k: float, *, epsabs: float = FIT_EPSABS,
                     epsrel: float = FIT_EPSREL) -> dict[str, object]:
    q = unit_amplitude_predictions(DATA.z, DATA.observable, omega_m, omega_k,
                                   epsabs=epsabs, epsrel=epsrel)
    assert isinstance(q, np.ndarray)
    return DATA.profile_alpha(q)


def fit_flat() -> dict[str, object]:
    def f(omega_m: float) -> float:
        return float(profile_at_shape(float(omega_m), 0.0)["chi2"])

    result = minimize_scalar(f, bounds=OMEGA_M_BOUNDS, method="bounded",
                             options={"xatol": 1.0e-12, "maxiter": 500})
    candidates = [(f(OMEGA_M_BOUNDS[0]), OMEGA_M_BOUNDS[0]),
                  (f(OMEGA_M_BOUNDS[1]), OMEGA_M_BOUNDS[1]),
                  (float(result.fun), float(result.x))]
    chi2, omega_m = min(candidates)
    row = profile_at_shape(omega_m, 0.0)
    parameters = {"alpha": float(row["alpha_profiled"]),
                  "Omega_m": float(omega_m), "Omega_k": 0.0}
    bounds = {"alpha": ALPHA_BOUNDS, "Omega_m": OMEGA_M_BOUNDS}
    return {
        "model": "flat_LCDM",
        "profile_domain": {"Omega_m": list(OMEGA_M_BOUNDS), "Omega_k": [0.0]},
        "parameters": parameters,
        "Omega_Lambda": 1.0 - omega_m,
        "chi2": chi2,
        "n_data": len(DATA.y),
        "bound_hits": hit_bounds(parameters, bounds),
        "optimizer": {"method": "bounded scalar minimization on Omega_m",
                      "success": bool(result.success), "message": str(result.message),
                      "nfev": int(result.nfev)},
        "gls_alpha": {k: v for k, v in row.items() if k not in ("prediction",)},
    }


def fit_curved(limit: float) -> dict[str, object]:
    bounds = [OMEGA_M_BOUNDS, (-limit, limit)]

    def objective(shape: np.ndarray) -> float:
        try:
            return float(profile_at_shape(float(shape[0]), float(shape[1]))["chi2"])
        except (ValueError, FloatingPointError, OverflowError):
            return 1.0e100

    # Coarse deterministic starts, followed by bounded local refinement from
    # the six best grid cells. This is intentionally small and CPU-only.
    coarse_om = np.linspace(*OMEGA_M_BOUNDS, 7)
    coarse_ok = np.linspace(-limit, limit, 9)
    coarse: list[tuple[float, float, float]] = []
    for omega_m in coarse_om:
        for omega_k in coarse_ok:
            coarse.append((objective(np.asarray([omega_m, omega_k])),
                           float(omega_m), float(omega_k)))
    coarse.sort(key=lambda item: item[0])
    starts = [(om, ok) for _, om, ok in coarse[:6]]
    starts.append((float(flat_fit["parameters"]["Omega_m"]), 0.0))
    refined = []
    for start in starts:
        result = minimize(objective, np.asarray(start, dtype=np.float64),
                          method="L-BFGS-B", bounds=bounds,
                          options={"ftol": 1.0e-14, "gtol": 1.0e-8,
                                   "maxiter": 500, "maxls": 40})
        refined.append(result)
    best = min(refined, key=lambda item: float(item.fun))
    omega_m, omega_k = map(float, best.x)
    profile = profile_at_shape(omega_m, omega_k)
    parameters = {"alpha": float(profile["alpha_profiled"]),
                  "Omega_m": omega_m, "Omega_k": omega_k}
    fit_bounds = {"alpha": ALPHA_BOUNDS, "Omega_m": OMEGA_M_BOUNDS,
                  "Omega_k": (-limit, limit)}
    return {
        "model": "curved_LCDM",
        "profile_domain": {"Omega_m": list(OMEGA_M_BOUNDS),
                           "Omega_k": [-limit, limit]},
        "parameters": parameters,
        "Omega_Lambda": 1.0 - omega_m - omega_k,
        "chi2": float(profile["chi2"]),
        "n_data": len(DATA.y),
        "bound_hits": hit_bounds(parameters, fit_bounds),
        "optimizer": {
            "method": "6 best cells from 7x9 coarse shape grid, each refined with bounded L-BFGS-B",
            "coarse_shape_grid": [len(coarse_om), len(coarse_ok)],
            "local_refinement_starts": len(starts),
            "successful_refinements": sum(bool(item.success) for item in refined),
            "best_success": bool(best.success),
            "best_message": str(best.message),
            "best_nit": int(best.nit),
            "best_nfev": int(best.nfev),
        },
        "gls_alpha": {k: v for k, v in profile.items() if k != "prediction"},
    }


def direct_fit(model: dict[str, object], omega_k_limit: float | None) -> dict[str, object]:
    parameters = model["parameters"]
    assert isinstance(parameters, dict)
    if omega_k_limit is None:
        bounds = [ALPHA_BOUNDS, OMEGA_M_BOUNDS]

        def objective(theta: np.ndarray) -> float:
            alpha, omega_m = map(float, theta)
            q = unit_amplitude_predictions(DATA.z, DATA.observable, omega_m, 0.0)
            assert isinstance(q, np.ndarray)
            return DATA.score(alpha * q)

        optimum = [float(parameters["alpha"]), float(parameters["Omega_m"])]
        mid = [30.0, 0.30]
        names = ["alpha", "Omega_m"]
    else:
        bounds = [ALPHA_BOUNDS, OMEGA_M_BOUNDS, (-omega_k_limit, omega_k_limit)]

        def objective(theta: np.ndarray) -> float:
            alpha, omega_m, omega_k = map(float, theta)
            q = unit_amplitude_predictions(DATA.z, DATA.observable, omega_m, omega_k)
            assert isinstance(q, np.ndarray)
            return DATA.score(alpha * q)

        optimum = [float(parameters["alpha"]), float(parameters["Omega_m"]),
                   float(parameters["Omega_k"])]
        mid = [30.0, 0.30, 0.0]
        names = ["alpha", "Omega_m", "Omega_k"]
    seeds = [optimum, mid]
    candidates = []
    for seed in seeds:
        clipped = np.asarray([np.clip(value, *bound)
                              for value, bound in zip(seed, bounds, strict=True)],
                             dtype=np.float64)
        result = minimize(objective, clipped, method="L-BFGS-B", bounds=bounds,
                          options={"ftol": 1.0e-14, "gtol": 1.0e-8,
                                   "maxiter": 500, "maxls": 40})
        candidates.append(result)
    best = min(candidates, key=lambda item: float(item.fun))
    return {
        "method": "direct bounded 2-parameter fit" if omega_k_limit is None
                  else "direct bounded 3-parameter chi2 fit",
        "parameters": {name: float(value) for name, value in zip(names, best.x, strict=True)},
        "chi2": float(best.fun),
        "success": bool(best.success),
        "message": str(best.message),
        "iterations": int(best.nit),
        "function_evaluations": int(best.nfev),
    }


def diagnostic_predictions(fit: dict[str, object], *, epsabs: float,
                           epsrel: float) -> tuple[np.ndarray, dict[str, float]]:
    parameters = fit["parameters"]
    assert isinstance(parameters, dict)
    q, quad_errors = unit_amplitude_predictions(
        DATA.z, DATA.observable, float(parameters["Omega_m"]),
        float(parameters["Omega_k"]), epsabs=epsabs, epsrel=epsrel,
        return_quad_errors=True)
    assert isinstance(q, np.ndarray)
    return float(parameters["alpha"]) * q, quad_errors


def main() -> None:
    global flat_fit
    flat_fit = fit_flat()
    flat_fit["direct_optimizer_check"] = direct_fit(flat_fit, None)
    flat_prediction, flat_quad_errors = diagnostic_predictions(
        flat_fit, epsabs=FIT_EPSABS, epsrel=FIT_EPSREL)
    flat_fit["quadrature_diagnostics"] = {
        "max_reported_abs_error_estimate": max(flat_quad_errors.values(), default=0.0),
        "mean_vector_chi2_recomputed": DATA.score(flat_prediction),
    }

    curved_fits = []
    for limit in OMEGA_K_LIMITS:
        fit = fit_curved(limit)
        fit["direct_optimizer_check"] = direct_fit(fit, limit)
        prediction, quad_errors = diagnostic_predictions(
            fit, epsabs=FIT_EPSABS, epsrel=FIT_EPSREL)
        prediction_tight, tight_errors = diagnostic_predictions(
            fit, epsabs=TIGHT_EPSABS, epsrel=TIGHT_EPSREL)
        fit["quadrature_diagnostics"] = {
            "max_reported_abs_error_estimate": max(quad_errors.values(), default=0.0),
            "max_reported_abs_error_estimate_tight": max(tight_errors.values(), default=0.0),
            "max_abs_prediction_difference_standard_vs_tight": float(
                np.max(np.abs(prediction - prediction_tight))),
            "chi2_recomputed": DATA.score(prediction),
            "chi2_tight": DATA.score(prediction_tight),
        }
        curved_fits.append(fit)

    flat_shape = flat_fit["parameters"]
    assert isinstance(flat_shape, dict)
    q_flat = unit_amplitude_predictions(DATA.z, DATA.observable,
                                        float(flat_shape["Omega_m"]), 0.0)
    q_curved_zero = unit_amplitude_predictions(DATA.z, DATA.observable,
                                               float(flat_shape["Omega_m"]), 0.0)
    assert isinstance(q_flat, np.ndarray) and isinstance(q_curved_zero, np.ndarray)
    q_tiny_plus = unit_amplitude_predictions(
        DATA.z, DATA.observable, float(flat_shape["Omega_m"]), 1.0e-12)
    q_tiny_minus = unit_amplitude_predictions(
        DATA.z, DATA.observable, float(flat_shape["Omega_m"]), -1.0e-12)
    assert isinstance(q_tiny_plus, np.ndarray) and isinstance(q_tiny_minus, np.ndarray)

    # Independent audit of the bounded-background physical domain.
    om_scan = np.linspace(*OMEGA_M_BOUNDS, 12)
    ok_scan = np.linspace(-0.20, 0.20, 17)
    z_scan = np.unique(np.concatenate((np.linspace(0.0, Z_MAX, 1001),
                                       np.asarray([5.0 / 3.0]))))
    x_scan = 1.0 + z_scan
    e2_scan = (1.0 + om_scan[:, None, None] * (x_scan[None, None, :]**3 - 1.0)
               + ok_scan[None, :, None] * (x_scan[None, None, :]**2 - 1.0))
    min_idx = np.unravel_index(int(np.argmin(e2_scan)), e2_scan.shape)
    min_e2 = float(e2_scan[min_idx])
    closed_sine_argument_bound = math.sqrt(0.20) * Z_MAX / math.sqrt(73.0 / 108.0)

    source_hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (MEAN_PATH, COV_PATH, PAPER_PATH, STARTER_PATH)
    }
    script_hash = sha256(Path(__file__).resolve())
    all_fits = [flat_fit, *curved_fits]
    nested_monotonic = all(
        curved_fits[i + 1]["chi2"] <= curved_fits[i]["chi2"] + 1.0e-8
        for i in range(len(curved_fits) - 1)
    )
    screen_comparison: dict[str, object] | None = None
    if SCREEN_PATH.exists():
        screen_result = json.loads(SCREEN_PATH.read_text(encoding="utf-8"))
        observed = screen_result["observed_fits"]
        comparisons = []
        for fit in all_fits:
            if fit["model"] == "flat_LCDM":
                reference = observed["flat_LCDM"]
                domain_label = "flat"
            else:
                limit = float(fit["profile_domain"]["Omega_k"][1])
                reference = observed["curved_LCDM_by_abs_Omega_k_bound"][str(limit)]
                domain_label = f"Omega_k_abs_le_{limit:g}"
            parameters = fit["parameters"]
            comparisons.append({
                "domain": domain_label,
                "delta_independent_minus_screen": {
                    "alpha": float(parameters["alpha"] - reference["alpha"]),
                    "Omega_m": float(parameters["Omega_m"] - reference["Omega_m"]),
                    "Omega_k": float(parameters["Omega_k"] - reference["Omega_k"]),
                    "chi2": float(fit["chi2"] - reference["chi2"]),
                },
                "screen_bound_hits": reference["parameter_bound_hits"],
                "independent_bound_hits": fit["bound_hits"],
            })
        screen_comparison = {
            "reference": "experiments/curvature_screen/result.json observed_fits; read-only",
            "sha256": sha256(SCREEN_PATH),
            "comparisons": comparisons,
            "interpretation": "Agreement within small optimizer/float64 differences; this is an independent confirmation, not a refit of the screen artifacts.",
        }
    output = {
        "status": "independent_check; exploratory profile screen",
        "scope": (
            "Curved late-time LCDM background only; radiation omitted; free alpha=c/(H0 rd); "
            "released DESI DR2 13-row mean vector and full covariance; no posterior or evidence."
        ),
        "source": {
            "paper": "DESI DR2 Results II, Phys. Rev. D 112, 083515 (2025), DOI 10.1103/tr6y-kpc6",
            "equations_checked": "PDF page 4, Eqs. (3)-(6): curved DM S_k, flat limit, DH=c/H, and Friedmann curvature term +Omega_k(1+z)^2.",
            "curvature_sign": "Omega_k>0 open and sinh branch; Omega_k<0 closed and sin branch, matching the DESI Eq. (3),(6) convention.",
        },
        "conventions": {
            "E2": "Omega_m*(1+z)^3 + Omega_k*(1+z)^2 + Omega_Lambda",
            "flat_closure": "Omega_Lambda=1-Omega_m-Omega_k",
            "dimensionless_radial_coordinate": "chi=int_0^z dz'/E(z')",
            "transverse_kernel": "S_k=sinh(sqrt(Omega_k)*chi)/sqrt(Omega_k) for Omega_k>0; chi at 0; sin(sqrt(|Omega_k|)*chi)/sqrt(|Omega_k|) for Omega_k<0",
            "near_flat_series": "S_k=chi*(1+y/6+y^2/120+y^3/5040+y^4/362880), y=Omega_k*chi^2, used for |y|<1e-4",
            "BAO_predictions": "DM/rd=alpha*S_k; DH/rd=alpha/E; DV/rd=[z*(DM/rd)^2*(DH/rd)]^(1/3)",
            "amplitude": "alpha=c/(H0 rd), dimensionless and free; BAO-only screen does not identify H0 and rd separately.",
            "likelihood": "chi2=(y-p)^T C^-1 (y-p), using all 13 rows and the unchanged full covariance.",
            "analytic_profile_alpha": (
                "For p=alpha*q, alpha_GLS=(q^T C^-1 y)/(q^T C^-1 q), "
                "chi2_min=y^T C^-1 y-(q^T C^-1 y)^2/(q^T C^-1 q), "
                "then alpha is clipped only if outside the inherited [1e-6,1e4] box."
            ),
        },
        "domain": {
            "Omega_m": list(OMEGA_M_BOUNDS),
            "Omega_k_limits_profiled": list(OMEGA_K_LIMITS),
            "Omega_k_full_interval": [-0.20, 0.20],
            "alpha": list(ALPHA_BOUNDS),
            "redshift_max": Z_MAX,
            "bounds_are_profile_screen_not_priors": True,
            "positivity_audit": {
                "grid_shape_Omega_m_Omega_k_z": list(e2_scan.shape),
                "minimum_E2_grid": min_e2,
                "minimum_location": {
                    "Omega_m": float(om_scan[min_idx[0]]),
                    "Omega_k": float(ok_scan[min_idx[1]]),
                    "z": float(z_scan[min_idx[2]]),
                },
                "analytic_continuous_minimum_E2": 73.0 / 108.0,
                "analytic_minimum_x": 8.0 / 3.0,
                "minimum_Omega_Lambda": 0.20,
                "closed_branch_argument_upper_bound": closed_sine_argument_bound,
                "closed_branch_argument_below_pi": closed_sine_argument_bound < math.pi,
            },
        },
        "data": {
            "n_data": int(len(DATA.y)),
            "mean_row_order": [
                {"z": float(zi), "observable": str(obs)}
                for zi, obs in zip(DATA.z, DATA.observable, strict=True)
            ],
            "covariance_shape": list(DATA.covariance.shape),
            "covariance_cholesky_success": True,
            "sha256": source_hashes,
        },
        "fits": all_fits,
        "comparison_to_existing_screen": screen_comparison,
        "checks": {
            "profile_alpha_vs_direct_optimizer": [
                {
                    "model": fit["model"],
                    "profile_chi2": fit["chi2"],
                    "direct_chi2": fit["direct_optimizer_check"]["chi2"],
                    "delta_chi2_direct_minus_profile": (
                        fit["direct_optimizer_check"]["chi2"] - fit["chi2"]
                    ),
                    "alpha_profiled": fit["parameters"]["alpha"],
                    "alpha_direct": fit["direct_optimizer_check"]["parameters"]["alpha"],
                }
                for fit in all_fits
            ],
            "flat_limit_at_same_Omega_m": {
                "max_abs_q_difference_curvature_zero_vs_flat": float(
                    np.max(np.abs(q_curved_zero - q_flat))),
                "max_abs_q_delta_at_Omega_k_plus_1e-12": float(
                    np.max(np.abs(q_tiny_plus - q_flat))),
                "max_abs_q_delta_at_Omega_k_minus_1e-12": float(
                    np.max(np.abs(q_tiny_minus - q_flat))),
            },
            "nested_domain_chi2_nonincreasing": nested_monotonic,
            "quadrature": {
                "fit_epsabs": FIT_EPSABS,
                "fit_epsrel": FIT_EPSREL,
                "tight_epsabs": TIGHT_EPSABS,
                "tight_epsrel": TIGHT_EPSREL,
                "method": "scalar scipy.integrate.quad (QUADPACK), independent per non-radial BAO row",
            },
        },
        "limitations": [
            "The curvature screens are finite profile domains, not posterior priors or Bayesian evidence.",
            "The free alpha leaves H0 and rd separately unidentified.",
            "Radiation is omitted and rd is not predicted, so this is not a precision early-universe or sound-horizon calculation.",
            "This is a background-only geometry check with no perturbations, growth, or external probes.",
        ],
        "execution": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_only": True,
            "blas_thread_environment": {
                key: os.environ.get(key) for key in
                ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
            },
            "script_sha256": script_hash,
        },
    }
    OUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    compact = {
        "output": str(OUT),
        "fits": [
            {"model": fit["model"], "domain": fit["profile_domain"],
             "parameters": fit["parameters"], "chi2": fit["chi2"],
             "bound_hits": fit["bound_hits"],
             "direct_chi2": fit["direct_optimizer_check"]["chi2"]}
            for fit in all_fits
        ],
        "checks": output["checks"],
        "E2_min_grid": output["domain"]["positivity_audit"]["minimum_E2_grid"],
        "script_sha256": script_hash,
    }
    print(json.dumps(compact, indent=2))


flat_fit: dict[str, object]

if __name__ == "__main__":
    main()

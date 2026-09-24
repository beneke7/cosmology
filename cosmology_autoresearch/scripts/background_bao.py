#!/usr/bin/env python3
"""Late-time, flat CPL background and compressed BAO likelihood starter.

This intentionally stops at z <= 3 and omits radiation, curvature, perturbations,
and a prediction for the sound horizon.  It is a background-distance exercise,
not a complete cosmological inference pipeline.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys
import time
from typing import Iterable

import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.optimize import minimize


FLOAT = np.float64
Z_MAX = 3.0
OBSERVABLES = {"DM_over_rs", "DH_over_rs", "DV_over_rs"}
GL_NODES, GL_WEIGHTS = np.polynomial.legendre.leggauss(96)
GL_NODES = np.asarray(GL_NODES, dtype=FLOAT)
GL_WEIGHTS = np.asarray(GL_WEIGHTS, dtype=FLOAT)


def _redshifts(z: float | Iterable[float] | np.ndarray) -> np.ndarray:
    """Return a finite float64 redshift vector inside the stated low-z domain."""
    z_arr = np.asarray(z, dtype=FLOAT)
    if not np.all(np.isfinite(z_arr)) or np.any(z_arr < 0.0) or np.any(z_arr > Z_MAX):
        raise ValueError(f"redshifts must be finite and in [0, {Z_MAX:g}]")
    return z_arr


def e2_cpl(
    z: float | Iterable[float] | np.ndarray,
    omega_m: float,
    w0: float = -1.0,
    wa: float = 0.0,
) -> np.ndarray:
    """Flat, matter + CPL dark-energy E(z)^2; radiation is omitted by design."""
    z_arr = _redshifts(z)
    omega_m, w0, wa = FLOAT(omega_m), FLOAT(w0), FLOAT(wa)
    if not np.all(np.isfinite([omega_m, w0, wa])) or not 0.0 <= omega_m <= 1.0:
        raise ValueError("parameters must be finite and Omega_m must be in [0, 1]")
    one_plus_z = 1.0 + z_arr
    matter = omega_m * one_plus_z**3
    omega_de = 1.0 - omega_m
    if omega_de == 0.0:
        return np.asarray(matter, dtype=FLOAT)
    log_de = 3.0 * (1.0 + w0 + wa) * np.log1p(z_arr) - 3.0 * wa * z_arr / one_plus_z
    with np.errstate(over="ignore", invalid="ignore"):
        de = omega_de * np.exp(log_de)
    result = np.asarray(matter + de, dtype=FLOAT)
    if not np.all(np.isfinite(result)) or np.any(result <= 0.0):
        raise ValueError("CPL expansion rate is not finite and positive on this redshift grid")
    return result


def comoving_integral(
    z: float | Iterable[float] | np.ndarray,
    omega_m: float,
    w0: float = -1.0,
    wa: float = 0.0,
) -> np.ndarray:
    """Return integral_0^z dz'/E(z') using fixed 96-point Gauss-Legendre quadrature."""
    z_arr = _redshifts(z)
    nodes = (z_arr[..., None] / 2.0) * (GL_NODES + 1.0)
    e2 = e2_cpl(nodes, omega_m, w0, wa)
    integral = (z_arr / 2.0) * np.sum(GL_WEIGHTS / np.sqrt(e2), axis=-1, dtype=FLOAT)
    return np.asarray(integral, dtype=FLOAT)


def predict_bao(
    z: float | Iterable[float] | np.ndarray,
    observable: str | Iterable[str],
    alpha: float,
    omega_m: float,
    w0: float = -1.0,
    wa: float = 0.0,
) -> np.ndarray:
    """Predict DM/rd, DH/rd, or DV/rd using alpha = c/(H0 rd).

    DESI's released text labels use ``rs``.  The amplitude is free here, so this
    computes the same dimensionless distances without predicting rd itself.
    """
    z_arr = _redshifts(z)
    alpha = FLOAT(alpha)
    if not np.isfinite(alpha) or alpha <= 0.0:
        raise ValueError("alpha = c/(H0 rd) must be finite and positive")
    if isinstance(observable, str):
        obs = np.full(z_arr.shape, observable, dtype=object)
    else:
        obs = np.asarray(list(observable), dtype=object)
        if obs.shape != z_arr.shape:
            raise ValueError("observable and redshift arrays must have identical shapes")
    unknown = set(np.atleast_1d(obs).tolist()) - OBSERVABLES
    if unknown:
        raise ValueError(f"unknown BAO observable(s): {sorted(unknown)}")

    e_at_z = np.sqrt(e2_cpl(z_arr, omega_m, w0, wa))
    dm = alpha * comoving_integral(z_arr, omega_m, w0, wa)
    dh = alpha / e_at_z
    dv = np.cbrt(z_arr * dm**2 * dh)
    out = np.empty(z_arr.shape, dtype=FLOAT)
    for name, values in (("DM_over_rs", dm), ("DH_over_rs", dh), ("DV_over_rs", dv)):
        mask = obs == name
        out[mask] = values[mask]
    return out


@dataclass(frozen=True)
class BAOData:
    """Ordered BAO rows and their covariance in the same row order."""

    z: np.ndarray
    value: np.ndarray
    observable: np.ndarray
    covariance: np.ndarray

    def __post_init__(self) -> None:
        z = np.asarray(self.z, dtype=FLOAT)
        value = np.asarray(self.value, dtype=FLOAT)
        observable = np.asarray(self.observable, dtype=str)
        covariance = np.asarray(self.covariance, dtype=FLOAT)
        if z.ndim != 1 or value.shape != z.shape or observable.shape != z.shape:
            raise ValueError("mean data must contain aligned one-dimensional z, value, observable columns")
        if len(z) == 0:
            raise ValueError("the mean file contains no measurements")
        if not np.all(np.isfinite(z)) or not np.all(np.isfinite(value)):
            raise ValueError("redshifts and measurements must be finite")
        if np.any(z <= 0.0) or np.any(z > Z_MAX):
            raise ValueError(f"BAO measurement redshifts must be in (0, {Z_MAX:g}]")
        unknown = set(observable.tolist()) - OBSERVABLES
        if unknown:
            raise ValueError(f"unknown BAO observable(s): {sorted(unknown)}")
        n = len(z)
        if covariance.shape != (n, n) or not np.all(np.isfinite(covariance)):
            raise ValueError(f"covariance must be a finite {n} by {n} matrix in mean-row order")
        if not np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-12):
            raise ValueError("covariance matrix is not symmetric")
        covariance = (covariance + covariance.T) / 2.0
        try:
            np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as exc:
            raise ValueError("covariance matrix must be positive definite") from exc
        for array in (z, value, observable, covariance):
            array.setflags(write=False)
        object.__setattr__(self, "z", z)
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "observable", observable)
        object.__setattr__(self, "covariance", covariance)

    @classmethod
    def from_files(cls, mean_path: str | Path, cov_path: str | Path) -> "BAOData":
        """Read DESI-style three-column means plus a plain numeric covariance."""
        z_rows: list[float] = []
        value_rows: list[float] = []
        observable_rows: list[str] = []
        path = Path(mean_path)
        for line_no, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_no}: expected exactly: z value observable")
            try:
                z_value, measurement = FLOAT(fields[0]), FLOAT(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: redshift and value must be numeric") from exc
            label = fields[2]
            if label not in OBSERVABLES:
                raise ValueError(
                    f"{path}:{line_no}: unsupported observable {label!r}; expected one of "
                    f"{', '.join(sorted(OBSERVABLES))}"
                )
            z_rows.append(z_value)
            value_rows.append(measurement)
            observable_rows.append(label)
        covariance = np.loadtxt(cov_path, dtype=FLOAT, comments="#", ndmin=2)
        return cls(
            np.asarray(z_rows, dtype=FLOAT),
            np.asarray(value_rows, dtype=FLOAT),
            np.asarray(observable_rows, dtype=str),
            covariance,
        )


class GaussianBAOLikelihood:
    """Correlated Gaussian likelihood evaluated through a Cholesky solve."""

    def __init__(self, data: BAOData):
        self.data = data
        self._chol = np.linalg.cholesky(data.covariance)

    def chi2(self, alpha: float, omega_m: float, w0: float = -1.0, wa: float = 0.0) -> float:
        prediction = predict_bao(self.data.z, self.data.observable, alpha, omega_m, w0, wa)
        residual = self.data.value - prediction
        whitened = solve_triangular(self._chol, residual, lower=True, check_finite=False)
        return float(np.dot(whitened, whitened))


FIT_BOUNDS = {
    "lcdm": [(1e-6, 1e4), (0.05, 0.6)],
    "wcdm": [(1e-6, 1e4), (0.05, 0.6), (-2.0, -0.3)],
    "cpl": [(1e-6, 1e4), (0.05, 0.6), (-2.0, -0.3), (-3.0, 3.0)],
}


def _decode_parameters(model: str, vector: np.ndarray) -> tuple[float, float, float, float]:
    alpha, omega_m = map(FLOAT, vector[:2])
    if model == "lcdm":
        return alpha, omega_m, FLOAT(-1.0), FLOAT(0.0)
    if model == "wcdm":
        return alpha, omega_m, FLOAT(vector[2]), FLOAT(0.0)
    if model == "cpl":
        return alpha, omega_m, FLOAT(vector[2]), FLOAT(vector[3])
    raise ValueError(f"unknown model {model!r}")


def fit_model(data: BAOData, model: str = "cpl", starts: int = 16) -> dict[str, float | str]:
    """Bounded deterministic multi-start maximum-likelihood background fit."""
    if model not in FIT_BOUNDS:
        raise ValueError(f"model must be one of {', '.join(FIT_BOUNDS)}")
    if starts < 1:
        raise ValueError("starts must be at least one")
    bounds = FIT_BOUNDS[model]
    likelihood = GaussianBAOLikelihood(data)

    def objective(vector: np.ndarray) -> float:
        try:
            return likelihood.chi2(*_decode_parameters(model, vector))
        except (ValueError, FloatingPointError, OverflowError):
            return 1e100

    rng = np.random.default_rng(20260923)
    initial = [np.asarray([(lo + hi) / 2.0 for lo, hi in bounds], dtype=FLOAT)]
    # Physically ordinary amplitudes get an explicit seed; the rest are reproducible,
    # log-uniform in alpha and uniform in each dimensionless parameter range.
    nominal = [30.0, 0.3, -1.0, 0.0][: len(bounds)]
    initial.append(np.asarray(nominal, dtype=FLOAT))
    for _ in range(max(0, starts - len(initial))):
        row = []
        for i, (lo, hi) in enumerate(bounds):
            if i == 0:
                row.append(np.exp(rng.uniform(np.log(3.0), np.log(300.0))))
            else:
                row.append(rng.uniform(lo, hi))
        initial.append(np.asarray(row, dtype=FLOAT))
    initial = initial[:starts]

    best = None
    for x0 in initial:
        result = minimize(
            objective,
            x0,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8},
        )
        if best is None or result.fun < best.fun:
            best = result
    assert best is not None
    alpha, omega_m, w0, wa = _decode_parameters(model, best.x)
    bound_hits = []
    for name, value, (lower, upper) in zip(["alpha", "Omega_m", "w0", "wa"], best.x, bounds):
        tolerance = 1e-7 * (1 + abs(float(value)))
        if abs(value - lower) <= tolerance:
            bound_hits.append({"parameter": name, "side": "lower", "bound": lower})
        if abs(value - upper) <= tolerance:
            bound_hits.append({"parameter": name, "side": "upper", "bound": upper})
    return {
        "model": model,
        "alpha": float(alpha),
        "Omega_m": float(omega_m),
        "w0": float(w0),
        "wa": float(wa),
        "chi2": float(best.fun),
        "n_data": len(data.z),
        "success": bool(best.success),
        "message": str(best.message),
        "parameter_bound_hits": bound_hits,
    }


def _predict_batch_torch(z: np.ndarray, observable: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    """Optional CUDA batch evaluator, kept separate from the float64 CPU fit path."""
    import torch

    if not torch.cuda.is_available():
        raise RuntimeError("PyTorch is installed but no CUDA device is available")
    z_arr = _redshifts(z).reshape(-1)
    obs = np.asarray(observable, dtype=str).reshape(-1)
    pars = np.asarray(parameters, dtype=FLOAT)
    if pars.ndim == 1:
        pars = pars[None, :]
    if pars.ndim != 2 or pars.shape[1] != 4 or obs.shape != z_arr.shape:
        raise ValueError("batch parameters must have shape (batch, 4) and observables must match z")
    if not np.all(np.isfinite(pars)) or np.any(pars[:, 0] <= 0.0) or np.any((pars[:, 1] < 0.0) | (pars[:, 1] > 1.0)):
        raise ValueError("batch rows must be finite [alpha, Omega_m, w0, wa] with alpha>0 and Omega_m in [0,1]")
    unknown = set(obs.tolist()) - OBSERVABLES
    if unknown:
        raise ValueError(f"unknown BAO observable(s): {sorted(unknown)}")

    device = torch.device("cuda")
    t = lambda x: torch.as_tensor(x, dtype=torch.float64, device=device)
    z_t, p_t = t(z_arr), t(pars)
    x = (z_t[:, None] / 2.0) * (t(GL_NODES)[None, :] + 1.0)
    one_plus_x = 1.0 + x
    log_de = (3.0 * (1.0 + p_t[:, 2, None, None] + p_t[:, 3, None, None])
              * torch.log1p(x)[None, :, :]
              - 3.0 * p_t[:, 3, None, None] * x[None, :, :] / one_plus_x[None, :, :])
    with torch.no_grad():
        e2_x = (p_t[:, 1, None, None] * one_plus_x[None, :, :] ** 3
                + (1.0 - p_t[:, 1, None, None]) * torch.exp(log_de))
        dm = p_t[:, 0, None] * (z_t[None, :] / 2.0) * torch.sum(
            t(GL_WEIGHTS)[None, None, :] / torch.sqrt(e2_x), dim=-1
        )
        zp = 1.0 + z_t
        log_de_z = 3.0 * (1.0 + p_t[:, 2, None] + p_t[:, 3, None]) * torch.log(zp)[None, :]
        log_de_z = log_de_z - 3.0 * p_t[:, 3, None] * z_t[None, :] / zp[None, :]
        e2_z = p_t[:, 1, None] * zp[None, :] ** 3 + (1.0 - p_t[:, 1, None]) * torch.exp(log_de_z)
        dh = p_t[:, 0, None] / torch.sqrt(e2_z)
        dv = torch.pow(z_t[None, :] * dm**2 * dh, 1.0 / 3.0)
        pred = torch.zeros_like(dm)
        for name, values in (("DM_over_rs", dm), ("DH_over_rs", dh), ("DV_over_rs", dv)):
            mask = t(obs == name).to(dtype=torch.bool)
            pred = torch.where(mask[None, :], values, pred)
    return pred.cpu().numpy()


def _self_check() -> None:
    z = np.asarray([0.0, 0.05, 0.2, 0.7, 1.5, 2.33, 3.0], dtype=FLOAT)
    e2_lcdm = e2_cpl(z, 0.3, -1.0, 0.0)
    expected_lcdm = 0.3 * (1.0 + z) ** 3 + 0.7
    np.testing.assert_allclose(e2_lcdm, expected_lcdm, rtol=2e-15, atol=0.0)
    print("PASS flat LCDM limit of CPL E(z)^2")

    gauss = comoving_integral(z[1:], 0.3, -0.91, 0.22)
    quad_values = np.asarray(
        [quad(lambda zp: 1.0 / np.sqrt(float(e2_cpl(zp, 0.3, -0.91, 0.22))), 0.0, zi,
              epsabs=1e-13, epsrel=1e-13)[0] for zi in z[1:]], dtype=FLOAT
    )
    np.testing.assert_allclose(gauss, quad_values, rtol=3e-13, atol=3e-14)
    print("PASS vector Gauss-Legendre distances against independent scipy.quad")

    z_eds = z[1:]
    eds = comoving_integral(z_eds, 1.0, -1.0, 0.0)
    eds_exact = 2.0 * (1.0 - 1.0 / np.sqrt(1.0 + z_eds))
    np.testing.assert_allclose(eds, eds_exact, rtol=3e-13, atol=3e-14)
    print("PASS Einstein-de Sitter analytic comoving-distance integral")

    # Rows deliberately interleave observables; a plain covariance has no labels,
    # so its only defensible alignment is by the unchanged mean-file row index.
    z_order = np.asarray([0.3, 0.5, 0.5], dtype=FLOAT)
    obs_order = np.asarray(["DV_over_rs", "DH_over_rs", "DM_over_rs"])
    expected_order = np.asarray(
        [predict_bao([zv], [name], 30.0, 0.3)[0] for zv, name in zip(z_order, obs_order)], dtype=FLOAT
    )
    cov = np.asarray([[0.04, 0.002, 0.0], [0.002, 0.09, -0.01], [0.0, -0.01, 0.16]], dtype=FLOAT)
    data = BAOData(z_order, expected_order + [0.1, -0.2, 0.05], obs_order, cov)
    chol = np.linalg.cholesky(data.covariance)
    residual = data.value - expected_order
    whitened = solve_triangular(chol, residual, lower=True)
    expected_chi2 = float(whitened @ whitened)
    np.testing.assert_allclose(GaussianBAOLikelihood(data).chi2(30.0, 0.3), expected_chi2, rtol=1e-14)
    assert np.array_equal(data.observable, obs_order)
    assert data.z.dtype == data.value.dtype == data.covariance.dtype == np.dtype(FLOAT)
    print("PASS covariance symmetry, positive definiteness, mean-row order, and Cholesky likelihood")

    z_syn = np.asarray([0.2, 0.3, 0.5, 0.5, 0.7, 0.7, 0.9, 1.1, 1.3, 1.3,
                        1.5, 1.8, 2.0, 2.33, 2.33, 2.7], dtype=FLOAT)
    o_syn = np.asarray(["DV_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs",
                        "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs",
                        "DM_over_rs", "DH_over_rs", "DV_over_rs", "DM_over_rs",
                        "DH_over_rs", "DH_over_rs", "DM_over_rs", "DV_over_rs"])
    mock_values = predict_bao(z_syn, o_syn, 29.5, 0.31, -1.0, 0.0)
    mock_cov = np.diag(np.full(len(z_syn), 0.01**2, dtype=FLOAT))
    mock_data = BAOData(z_syn, mock_values, o_syn, mock_cov)
    recovered = fit_model(mock_data, "cpl", starts=8)
    assert recovered["chi2"] < 1e-6, recovered
    np.testing.assert_allclose([recovered["alpha"], recovered["Omega_m"], recovered["w0"], recovered["wa"]],
                               [29.5, 0.31, -1.0, 0.0], rtol=0.0, atol=2e-4)
    print("PASS deterministic CPL multi-start recovers a nested synthetic Lambda model")


def _benchmark() -> None:
    """Small synthetic CPU benchmark; torch/GPU is opportunistic, never required."""
    z = np.linspace(0.01, Z_MAX, 13, dtype=FLOAT)
    observable = np.asarray(["DM_over_rs", "DH_over_rs", "DV_over_rs"] * 4 + ["DM_over_rs"])
    start = time.perf_counter()
    checksum = 0.0
    for _ in range(3000):
        checksum += float(np.sum(predict_bao(z, observable, 30.0, 0.3, -0.95, 0.1)))
    elapsed = time.perf_counter() - start
    print(f"Synthetic NumPy CPU float64: 3000 x 13 predictions in {elapsed:.3f} s (checksum={checksum:.6g})")
    try:
        import torch  # optional; no dependency is required for the CPU path
    except ImportError:
        print("Optional torch GPU batch path: unavailable (torch is not installed)")
        return
    if not torch.cuda.is_available():
        print("Optional torch GPU batch path: unavailable (CUDA device not detected)")
        return
    parameters = np.tile(np.asarray([30.0, 0.3, -0.95, 0.1], dtype=FLOAT), (1024, 1))
    start = time.perf_counter()
    gpu_predictions = _predict_batch_torch(z, observable, parameters)
    gpu_elapsed = time.perf_counter() - start
    cpu_reference = predict_bao(z, observable, *parameters[0])
    np.testing.assert_allclose(gpu_predictions[0], cpu_reference, rtol=2e-12, atol=2e-12)
    print(f"Optional torch CUDA float64 batch: 1024 x 13 predictions in {gpu_elapsed:.3f} s; checked against CPU reference")


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mean", help="DESI-style mean file: z value observable")
    parser.add_argument("--cov", help="plain numeric covariance matrix in mean-file row order")
    parser.add_argument("--model", choices=tuple(FIT_BOUNDS), default="cpl")
    parser.add_argument("--starts", type=int, default=16, help="deterministic multi-start optimizations")
    parser.add_argument("--self-check", action="store_true", help="run numerical and synthetic recovery checks")
    parser.add_argument("--benchmark", action="store_true", help="run a synthetic float64 CPU benchmark")
    args = parser.parse_args(argv)
    if args.self_check:
        _self_check()
        return 0
    if args.benchmark:
        _benchmark()
        return 0
    if bool(args.mean) != bool(args.cov):
        parser.error("--mean and --cov must be supplied together")
    if not args.mean:
        parser.error("provide --self-check, --benchmark, or both --mean and --cov")
    data = BAOData.from_files(args.mean, args.cov)
    result = fit_model(data, args.model, args.starts)
    report = {
        "analysis": "exploratory_background_only_bao_fit",
        "scope": "flat late-time CPL background, radiation neglected, 0 < z <= 3",
        "sound_horizon": "not predicted; alpha = c/(H0*rd) is free",
        "early_universe_likelihood": False,
        "evidence_or_significance_claim": False,
        "mean_file": str(Path(args.mean)),
        "covariance_file": str(Path(args.cov)),
        "model": result["model"],
        "n_data": result["n_data"],
        "chi2": result["chi2"],
        "parameters": {
            "alpha": result["alpha"],
            "Omega_m": result["Omega_m"],
            "w0": result["w0"],
            "wa": result["wa"],
        },
        "optimizer": {"success": result["success"], "message": result["message"]},
        "parameter_bound_hits": result["parameter_bound_hits"],
        "starts": args.starts,
        "optimizer_seed": 20260923,
        "parameter_order": ["alpha", "Omega_m", "w0", "wa"][:len(FIT_BOUNDS[args.model])],
        "parameter_bounds": FIT_BOUNDS[args.model],
        "bounds_are_not_a_posterior_prior": True,
    }
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if result["success"] else 2


if __name__ == "__main__":
    sys.exit(_main())

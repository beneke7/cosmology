#!/usr/bin/env python3
"""Check the precision-Cholesky Gaussian sampler against released Dovekie P.

For P = L L.T with NumPy's lower Cholesky factor L, the required draw is
epsilon = L^{-T} z.  This script checks the whitening identity on the full
released matrix and compares a Monte Carlo set of scalar projections with
their exact covariance.  It writes no files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular


ROOT = Path(__file__).resolve().parents[2]
PRECISION_PATH = ROOT / "context" / "data" / "des_dovekie_stat_sys.npz"
N_DRAWS = 4096
SEED = 20260925


def main() -> None:
    with np.load(PRECISION_PATH, allow_pickle=False) as archive:
        n = int(np.asarray(archive["nsn"]).reshape(-1)[0])
        packed = np.asarray(archive["cov"])
    if packed.size != n * (n + 1) // 2 or packed.dtype != np.float32:
        raise ValueError("unexpected released packed precision shape or dtype")

    precision = np.zeros((n, n), dtype=np.float64)
    upper = np.triu_indices(n)
    precision[upper] = packed.astype(np.float64)
    lower = np.tril_indices(n, -1)
    precision[lower] = precision.T[lower]
    chol = np.linalg.cholesky(precision)

    rng = np.random.default_rng(SEED)
    z = rng.standard_normal((n, N_DRAWS), dtype=np.float64)
    # Solve L.T epsilon = z; each column is one draw in the exact release row order.
    epsilon = solve_triangular(chol.T, z, lower=False, check_finite=False)

    whitened = chol.T @ epsilon
    whitening_abs = float(np.max(np.abs(whitened - z)))
    whitening_rel = whitening_abs / max(float(np.max(np.abs(z))), 1.0)

    # The precision quadratic of each draw must equal the squared norm of its
    # generating standard-normal vector: epsilon.T P epsilon == z.T z.
    weighted = precision @ epsilon
    quadratic = np.sum(epsilon * weighted, axis=0)
    z_quadratic = np.sum(z * z, axis=0)
    quadratic_abs = float(np.max(np.abs(quadratic - z_quadratic)))
    quadratic_rel = quadratic_abs / max(float(np.max(z_quadratic)), 1.0)

    # Check sampled variances/covariances of fixed random contrasts against
    # a.T P^-1 a = ||L^-1 a||^2. Keep this low-dimensional so the check is
    # cheap while still using the actual 1820-dimensional precision matrix.
    contrasts = rng.standard_normal((n, 8), dtype=np.float64)
    contrasts /= np.linalg.norm(contrasts, axis=0, keepdims=True)
    samples = contrasts.T @ epsilon
    exact_whitened_contrasts = solve_triangular(
        chol, contrasts, lower=True, check_finite=False
    )
    exact_covariance = exact_whitened_contrasts.T @ exact_whitened_contrasts
    sample_covariance = np.cov(samples, ddof=1)
    variance_ratios = np.diag(sample_covariance) / np.diag(exact_covariance)

    print(f"n={n}; draws={N_DRAWS}; seed={SEED}")
    print(f"max |L.T @ epsilon - z| = {whitening_abs:.3e}")
    print(f"max whitened relative error = {whitening_rel:.3e}")
    print(f"max |epsilon.T @ P @ epsilon - z.T @ z| = {quadratic_abs:.3e}")
    print(f"max quadratic relative error = {quadratic_rel:.3e}")
    print("sample/exact contrast variance ratios =", np.array2string(variance_ratios, precision=3))
    print(
        "max absolute contrast covariance relative error = "
        f"{np.max(np.abs(sample_covariance - exact_covariance) / np.sqrt(np.outer(np.diag(exact_covariance), np.diag(exact_covariance)))):.3f}"
    )

    if whitening_rel > 2e-13 or quadratic_rel > 2e-12:
        raise AssertionError("precision-Cholesky sampling identity failed")
    if np.any((variance_ratios < 0.90) | (variance_ratios > 1.10)):
        raise AssertionError("Monte Carlo contrast variances outside broad 10% check")


if __name__ == "__main__":
    main()

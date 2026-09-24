#!/usr/bin/env python3
"""Synthetic reference check for Gaussian CV from a full precision matrix.

The data release supplies P=C^{-1}. For held-out H and training T:

    Cov(y_H | y_T)^{-1} = P_HH
    E[y_H | y_T] = mu_H - P_HH^{-1} P_HT (y_T-mu_T)
    Cov(y_T)^{-1} = P_TT - P_TH P_HH^{-1} P_HT.

This check compares those identities with a separately formed covariance C and
its direct marginal/conditional Gaussian calculations. It uses synthetic data
only; it does not score the DES-Dovekie observations.
"""

from __future__ import annotations

import json

import numpy as np


def main() -> None:
    rng = np.random.default_rng(20260924)
    n = 60
    design = rng.normal(size=(n, n))
    covariance = design @ design.T / n + 0.25 * np.eye(n)
    precision = np.linalg.inv(covariance)
    y = rng.normal(size=n)
    base = rng.normal(size=n)
    groups = np.array_split(np.arange(n), 4)

    max_offset_error = 0.0
    max_conditional_mean_error = 0.0
    max_conditional_covariance_error = 0.0
    max_conditional_score_error = 0.0
    max_intercept_marginal_covariance_error = 0.0
    max_intercept_marginal_score_error = 0.0
    max_training_precision_error = 0.0
    max_training_normal_equation_residual = 0.0
    min_precision_eigenvalue = float(np.linalg.eigvalsh(precision)[0])

    for held in groups:
        train = np.setdiff1d(np.arange(n), held)
        p_tt = precision[np.ix_(train, train)]
        p_th = precision[np.ix_(train, held)]
        p_ht = precision[np.ix_(held, train)]
        p_hh = precision[np.ix_(held, held)]

        # Schur complement gives the inverse marginal covariance on training.
        p_hh_inv_p_ht = np.linalg.solve(p_hh, p_ht)
        q_train = p_tt - p_th @ p_hh_inv_p_ht
        q_train = 0.5 * (q_train + q_train.T)
        np.linalg.cholesky(q_train)

        one = np.ones(train.size)
        residual0 = y[train] - base[train]
        offset = float(one @ (q_train @ residual0) / (one @ (q_train @ one)))

        # Reference training-only intercept from C_TT^{-1}.
        c_tt = covariance[np.ix_(train, train)]
        q_direct = np.linalg.inv(c_tt)
        offset_direct = float(one @ (q_direct @ residual0) / (one @ (q_direct @ one)))
        max_offset_error = max(max_offset_error, abs(offset - offset_direct))
        max_training_precision_error = max(
            max_training_precision_error, float(np.max(np.abs(q_train - q_direct)))
        )
        normal = float(one @ (q_train @ (residual0 - offset * one)))
        max_training_normal_equation_residual = max(
            max_training_normal_equation_residual, abs(normal)
        )

        mu_train = base[train] + offset
        mu_held = base[held] + offset
        mean_precision = mu_held - np.linalg.solve(p_hh, p_ht @ (y[train] - mu_train))

        c_ht = covariance[np.ix_(held, train)]
        c_th = covariance[np.ix_(train, held)]
        c_hh = covariance[np.ix_(held, held)]
        mean_direct = mu_held + c_ht @ np.linalg.solve(c_tt, y[train] - mu_train)
        s_cond = c_hh - c_ht @ np.linalg.solve(c_tt, c_th)
        max_conditional_mean_error = max(
            max_conditional_mean_error,
            float(np.max(np.abs(mean_precision - mean_direct))),
        )
        max_conditional_covariance_error = max(
            max_conditional_covariance_error,
            float(np.max(np.abs(s_cond - np.linalg.inv(p_hh)))),
        )

        residual_precision = y[held] - mean_precision
        residual_direct = y[held] - mean_direct
        score_precision = float(residual_precision @ (p_hh @ residual_precision))
        score_direct = float(residual_direct @ np.linalg.solve(s_cond, residual_direct))
        max_conditional_score_error = max(
            max_conditional_score_error, abs(score_precision - score_direct)
        )

        # Integrate the training-only flat-prior magnitude intercept. Its
        # posterior variance adds a rank-one term to the held-out covariance.
        sigma2 = 1.0 / float(one @ (q_train @ one))
        one_held = np.ones(held.size)
        v_precision = one_held + np.linalg.solve(p_hh, p_ht @ one)
        v_direct = one_held - c_ht @ np.linalg.solve(c_tt, one)
        predictive_precision_cov = np.linalg.inv(p_hh) + sigma2 * np.outer(v_precision, v_precision)
        sigma2_direct = 1.0 / float(one @ (q_direct @ one))
        predictive_direct_cov = s_cond + sigma2_direct * np.outer(v_direct, v_direct)
        max_intercept_marginal_covariance_error = max(
            max_intercept_marginal_covariance_error,
            float(np.max(np.abs(predictive_precision_cov - predictive_direct_cov))),
        )
        score_precision_marginal = float(
            residual_precision @ np.linalg.solve(predictive_precision_cov, residual_precision)
        )
        score_direct_marginal = float(
            residual_direct @ np.linalg.solve(predictive_direct_cov, residual_direct)
        )
        max_intercept_marginal_score_error = max(
            max_intercept_marginal_score_error,
            abs(score_precision_marginal - score_direct_marginal),
        )

    result = {
        "synthetic_only": True,
        "seed": 20260924,
        "dimension": n,
        "folds": len(groups),
        "all_training_schur_complements_cholesky_positive": True,
        "max_training_intercept_absolute_difference": max_offset_error,
        "max_training_precision_absolute_difference": max_training_precision_error,
        "max_training_normal_equation_residual": max_training_normal_equation_residual,
        "max_conditional_mean_absolute_difference": max_conditional_mean_error,
        "max_conditional_covariance_absolute_difference": max_conditional_covariance_error,
        "max_conditional_chi2_absolute_difference": max_conditional_score_error,
        "max_intercept_marginal_predictive_covariance_absolute_difference": max_intercept_marginal_covariance_error,
        "max_intercept_marginal_predictive_chi2_absolute_difference": max_intercept_marginal_score_error,
        "minimum_full_precision_eigenvalue": min_precision_eigenvalue,
        "disposition": "PASS" if max(
            max_offset_error,
            max_training_precision_error,
            max_training_normal_equation_residual,
            max_conditional_mean_error,
            max_conditional_covariance_error,
            max_conditional_score_error,
            max_intercept_marginal_covariance_error,
            max_intercept_marginal_score_error,
        ) < 1e-9 else "FAIL",
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["disposition"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

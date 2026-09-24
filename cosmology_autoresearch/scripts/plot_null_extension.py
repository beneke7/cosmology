#!/usr/bin/env python3
"""Recompute and plot the pooled DESI BAO profile-search null calibration."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta


def stats(values: np.ndarray, threshold: float) -> dict:
    values = np.asarray(values, dtype=np.float64)
    n = int(values.size)
    k = int(np.count_nonzero(values >= threshold))
    ci = [0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1)),
          1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))]
    return {"n": n, "exceedances": k, "fraction": k / n,
            "clopper_pearson_95pct": ci,
            "median": float(np.median(values)),
            "q90_linear": float(np.quantile(values, 0.90, method="linear")),
            "q95_linear": float(np.quantile(values, 0.95, method="linear")),
            "maximum": float(np.max(values))}


def collect(path: Path, statistic_key: str) -> np.ndarray:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload["adaptive_search_null_mock_calibration"]["realizations"]
    return np.asarray([row[statistic_key] for row in rows], dtype=np.float64)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", default="experiments/bao_robustness/result.json")
    parser.add_argument("--extension", default="experiments/bao_robustness/null_extension_1000_seed20260925.json")
    parser.add_argument("--output", default="experiments/bao_robustness/null_extension_pooled")
    args = parser.parse_args()
    original_path, extension_path = Path(args.original), Path(args.extension)
    original = json.loads(original_path.read_text(encoding="utf-8"))
    extension = json.loads(extension_path.read_text(encoding="utf-8"))
    labels = [("16-start search", "search_statistic_delta_chi2", "16_starts"),
              ("64-start convergence check", "expanded_search_statistic_delta_chi2", "64_starts")]
    pooled: dict[str, dict] = {}
    arrays: dict[str, np.ndarray] = {}
    thresholds: dict[str, float] = {}
    for _, key, budget in labels:
        threshold_a = original["adaptive_search_null_mock_calibration"]["observed_searches"][budget]["search_statistic_delta_chi2"]
        threshold_b = extension["adaptive_search_null_mock_calibration"]["observed_searches"][budget]["search_statistic_delta_chi2"]
        if threshold_a != threshold_b:
            raise ValueError(f"observed threshold differs across saved runs for {budget}")
        a = collect(original_path, key)
        b = collect(extension_path, key)
        values = np.concatenate((a, b))
        if not np.all(np.isfinite(values)):
            raise ValueError(f"non-finite mock statistic in {budget}")
        arrays[budget] = values
        thresholds[budget] = float(threshold_a)
        pooled[budget] = stats(values, threshold_a)

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.7), sharex=True, sharey=True,
                             constrained_layout=True)
    bins = np.linspace(0.0, 14.0, 36)
    for ax, (title, _, budget) in zip(axes, labels):
        values = arrays[budget]
        threshold = thresholds[budget]
        ax.hist(values[values < threshold], bins=bins, color="#3679a8", alpha=0.88,
                edgecolor="white", linewidth=0.45)
        ax.hist(values[values >= threshold], bins=bins, color="#d96b48", alpha=0.92,
                edgecolor="white", linewidth=0.45, label="at/above observed")
        ax.axvline(threshold, color="#8f2d1f", linestyle="--", linewidth=1.5,
                   label=fr"observed $\Delta\chi^2={threshold:.4f}$")
        summary = pooled[budget]
        low, high = summary["clopper_pearson_95pct"]
        ax.text(0.97, 0.96,
                f"N={summary['n']:,}; tail={summary['exceedances']}/{summary['n']}\n"
                f"fraction={summary['fraction']:.4f} [{low:.4f}, {high:.4f}] (exact 95%)\n"
                f"null q95={summary['q95_linear']:.3f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
                bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.92,
                       "edgecolor": "#bbbbbb"})
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.2)
        ax.set_xlabel(r"Selected $\Delta\chi^2$ under fitted flat-$\Lambda$CDM null")
        ax.legend(frameon=False, loc="upper left", fontsize=8)
    axes[0].set_ylabel("Mock realizations per bin")
    fig.suptitle("DESI DR2 BAO: finite-search null calibration (profile screen only)", fontsize=12)
    fig.text(0.5, -0.025,
             "Pooled disjoint PCG64 seeds 20260924 and 20260925; 13-row full-covariance Gaussian mocks. "
             "Not a posterior probability, evidence, or discovery significance.",
             ha="center", fontsize=8)
    base = Path(args.output)
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".png"), dpi=180,
                metadata={"Title": "DESI DR2 BAO finite-search null calibration",
                          "Description": "Pooled 1200-mock parametric bootstrap of a fixed profile search"})
    fig.savefig(base.with_suffix(".svg"),
                metadata={"Title": "DESI DR2 BAO finite-search null calibration"})
    print(json.dumps({"status": "recomputed", "pooled_by_start_budget": pooled,
                      "png": str(base.with_suffix(".png")),
                      "svg": str(base.with_suffix(".svg"))}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

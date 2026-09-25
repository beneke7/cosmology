#!/usr/bin/env python3
"""Plot the saved matched-LCDM Dovekie IVS null calibration."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "dovekie-ivs-null-calibration-v1"
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "work/compute11"
JSONL = RUN / "mock_results.jsonl"
SUMMARY = RUN / "result.json"
PNG = Path(__file__).resolve().parent / "null_calibration.png"
SVG = Path(__file__).resolve().parent / "null_calibration.svg"


def main() -> None:
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    rows = [json.loads(line) for line in JSONL.read_text(encoding="utf-8").splitlines() if line]
    complete = [row for row in rows if row.get("status") == "complete"]
    if len(complete) != summary["completed_realizations"] or len(complete) != 200:
        raise ValueError("saved JSONL does not contain the expected 200 completed mocks")
    values = np.asarray([row["T_chi2_lcdm_minus_chi2_ivs"] for row in complete], dtype=float)
    threshold = float(summary["statistic"]["inclusive_threshold"])
    exceedances = int(np.count_nonzero(values >= threshold))
    recorded = summary["tail_probability_estimate"]["exceedance_count"]
    if exceedances != recorded:
        raise ValueError(f"plot count {exceedances} differs from the saved result {recorded}")
    frac = exceedances / len(values)
    lo, hi = summary["tail_probability_estimate"]["clopper_pearson_95_interval"]

    fig, ax = plt.subplots(figsize=(9.2, 5.4), layout="constrained")
    weights = np.full(len(values), 1.0 / len(values))
    bin_edges = np.histogram_bin_edges(values, bins="fd")
    ax.hist(values, bins=bin_edges, weights=weights, color="#4C78A8", alpha=0.82,
            edgecolor="white", linewidth=0.8, label="Matched LCDM mock distribution")
    ax.axvline(threshold, color="#E45756", linewidth=2.2, linestyle="--",
               label=f"Observed statistic = {threshold:.3f}")
    ax.set_xlabel(r"$T = \chi^2_{\Lambda\mathrm{CDM}} - \chi^2_{\mathrm{IVS}}$")
    ax.set_ylabel("Fraction of mocks per bin")
    ax.set_title("Dovekie IVS profile: finite-search LCDM-null calibration")
    ax.text(
        0.98,
        0.96,
        f"Exceedances: {exceedances}/{len(values)} = {frac:.1%}\n"
        f"Exact 95% Clopper–Pearson interval: [{lo:.3f}, {hi:.3f}]\n"
        "200 full-covariance mocks; identical finite search",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "alpha": 0.9,
              "edgecolor": "#BBBBBB"},
    )
    ax.grid(axis="y", alpha=0.22)
    ax.legend(loc="upper left", frameon=False)
    fig.savefig(PNG, dpi=180)
    fig.savefig(SVG, metadata={"Date": None})
    # Matplotlib emits trailing blanks in multiline SVG path data. They are
    # semantically irrelevant but fail git's whitespace check; normalize them
    # so the checked-in vector artifact stays reproducible and reviewable.
    svg_text = SVG.read_text(encoding="utf-8")
    SVG.write_text("\n".join(line.rstrip(" \t") for line in svg_text.splitlines()) + "\n",
                   encoding="utf-8")
    plt.close(fig)
    print(json.dumps({
        "status": "complete",
        "mocks": len(values),
        "exceedances": exceedances,
        "empirical_fraction": frac,
        "clopper_pearson_95": [float(lo), float(hi)],
        "png": str(PNG.relative_to(ROOT)),
        "svg": str(SVG.relative_to(ROOT)),
    }, indent=2))


if __name__ == "__main__":
    main()

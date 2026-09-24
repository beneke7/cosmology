#!/usr/bin/env python3
"""Plot the predeclared frozen-shape Dovekie conditional-CV result."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RESULT = HERE / "result.json"


def main() -> None:
    data = json.loads(RESULT.read_text(encoding="utf-8"))
    if data.get("status") != "complete_exploratory_cross_probe_prediction_screen":
        raise ValueError("refusing to plot a result without a complete status")
    folds = data["fold_contract"]["folds"]
    variants = data["fold_contract"]["predeclared_edge_sensitivity"]["variants"]

    fig, (ax_folds, ax_edges) = plt.subplots(2, 1, figsize=(11.2, 8.8), constrained_layout=True)
    fold_y = np.arange(len(folds))
    fold_delta = np.asarray(
        [
            f["models"]["BAO_interacting_vacuum"]["integrated_intercept_conditional_chi2"]
            - f["models"]["BAO_flat_LCDM"]["integrated_intercept_conditional_chi2"]
            for f in folds
        ],
        dtype=np.float64,
    )
    fold_labels = [
        f"{f['zHD_min']:.3f}–{f['zHD_max']:.3f}  (n={f['heldout_count']})"
        for f in folds
    ]
    ax_folds.axvline(0.0, color="0.2", lw=1)
    ax_folds.barh(fold_y, fold_delta, color=np.where(fold_delta < 0.0, "#287c8e", "#d17a39"))
    ax_folds.set_yticks(fold_y, fold_labels)
    ax_folds.invert_yaxis()
    ax_folds.set_xlabel(r"$\Delta\chi^2_{\mathrm{cond}}$ (IVS $-$ $\Lambda$CDM)")
    ax_folds.set_title("Primary four-fold conditional scores")
    ax_folds.grid(axis="x", alpha=0.22)
    ax_folds.set_xlim(-6.7, 1.5)
    for y, value in zip(fold_y, fold_delta, strict=True):
        ax_folds.text(value + (0.10 if value >= 0 else -0.10), y, f"{value:+.2f}",
                      va="center", ha="left" if value >= 0 else "right", fontsize=9)

    edge_delta = np.asarray(
        [v["delta_IVS_minus_LCDM_integrated_conditional_neg2logpredictive"] for v in variants],
        dtype=np.float64,
    )
    edge_names = [v["variant"].replace("boundary_", "cut ").replace("_", " ") for v in variants]
    edge_y = np.arange(len(variants))
    primary = float(data["aggregate_predictive_chi2"][
        "delta_IVS_minus_LCDM_integrated_conditional_neg2logpredictive"
    ])
    ax_edges.axvline(0.0, color="0.2", lw=1)
    ax_edges.axvline(primary, color="#6a4c93", lw=1.6, ls="--", label=f"primary = {primary:+.2f}")
    ax_edges.scatter(edge_delta, edge_y, s=58, color="#287c8e", zorder=3)
    ax_edges.set_yticks(edge_y, edge_names)
    ax_edges.invert_yaxis()
    ax_edges.set_xlabel(r"Sum of conditional $\Delta(-2\log p)$ (IVS $-$ $\Lambda$CDM)")
    ax_edges.set_title("Predeclared one-edge ±1%-N perturbations")
    ax_edges.grid(axis="x", alpha=0.22)
    ax_edges.set_xlim(-4.65, 0.25)
    ax_edges.legend(frameon=False, loc="upper right")
    for y, value in zip(edge_y, edge_delta, strict=True):
        ax_edges.text(value + 0.05, y, f"{value:.2f}", va="center", ha="left", fontsize=8)

    fig.suptitle(
        "BAO-frozen shapes scored on DES-Dovekie SNe\n"
        "Full STAT+SYS covariance; training-only magnitude intercept integrated\n"
        "Negative favors IVS descriptively; leave-one-bin-out sums are not a joint likelihood or significance",
        fontsize=12,
    )
    fig.savefig(HERE / "conditional_cv.png", dpi=180, bbox_inches="tight")
    svg_path = HERE / "conditional_cv.svg"
    fig.savefig(svg_path, bbox_inches="tight")
    # Matplotlib wraps SVG path coordinates with trailing spaces; strip them so
    # the committed vector artifact passes whitespace checks without changing
    # its rendered geometry.
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text("\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n", encoding="utf-8")
    plt.close(fig)


if __name__ == "__main__":
    main()

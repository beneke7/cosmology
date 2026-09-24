#!/usr/bin/env python3
"""Render residual and parameter-boundary diagnostics from run_bao_robustness.py."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from background_bao import BAOData, FIT_BOUNDS, predict_bao


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", default="experiments/bao_robustness/result.json")
    parser.add_argument("--mean", default="context/data/desi_dr2_mean.txt")
    parser.add_argument("--cov", default="context/data/desi_dr2_cov.txt")
    parser.add_argument("--output", default="experiments/bao_robustness/bao_screen_diagnostics")
    args = parser.parse_args()

    result = json.loads(Path(args.result).read_text(encoding="utf-8"))
    data = BAOData.from_files(args.mean, args.cov)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), constrained_layout=True)
    ax = axes[0]
    x = np.arange(1, len(data.z) + 1)
    sigma = np.sqrt(np.diag(data.covariance))
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.axhspan(-1.0, 1.0, color="#eeeeee", zorder=0)
    lya_idx = np.flatnonzero(np.isclose(data.z, 2.33))
    if len(lya_idx):
        ax.axvspan(lya_idx[0] + 0.55, lya_idx[-1] + 1.45,
                   color="#e7f1fb", zorder=0, label="z = 2.33 block")
    colors = {"lcdm": "#2455a4", "wcdm": "#e07a16", "cpl": "#21835b"}
    offsets = {"lcdm": -0.12, "wcdm": 0.0, "cpl": 0.12}
    labels = {"lcdm": r"$\Lambda$CDM", "wcdm": r"constant $w$", "cpl": "CPL"}
    for model, fit in result["baseline_profile_fits"].items():
        params = [fit[name] for name in ("alpha", "Omega_m", "w0", "wa")]
        pred = predict_bao(data.z, data.observable, *params)
        pulls = (data.value - pred) / sigma
        ax.scatter(x + offsets[model], pulls, color=colors[model], s=30,
                   label=labels[model], zorder=3)
    ax.set_xlabel("DESI DR2 mean-vector row (released order)")
    ax.set_ylabel(r"Residual / marginal $1\sigma$")
    ax.set_title("BAO-only profile-fit residuals")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{i}\n{obs.replace('_over_rs', '')}" for i, obs in zip(x, data.observable)],
                       rotation=0, fontsize=7)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.text(0.02, 0.02,
            "Diagonal errors shown only for display; fits use the full covariance.",
            transform=ax.transAxes, fontsize=7, va="bottom")

    ax = axes[1]
    profiles = result["cpl_w0_upper_bound_sensitivity"]
    upper = np.asarray([item["w0_upper"] for item in profiles])
    chi2 = np.asarray([item["chi2"] for item in profiles])
    lcdm = result["baseline_profile_fits"]["lcdm"]["chi2"]
    ax.axhline(lcdm, color=colors["lcdm"], linestyle="--", linewidth=1.2,
               label=fr"$\Lambda$CDM minimum ($\chi^2={lcdm:.2f}$)")
    ax.plot(upper, chi2, "o-", color=colors["cpl"], linewidth=1.5,
            label="CPL profile minimum")
    for xval, yval, row in zip(upper, chi2, profiles):
        if row["parameter_bound_hits"]:
            ax.annotate("w₀ boundary", (xval, yval), xytext=(3, 7),
                        textcoords="offset points", fontsize=7)
    ax.set_xlabel(r"Allowed upper bound on $w_0$")
    ax.set_ylabel(r"Minimum profile $\chi^2$ (13 rows)")
    ax.set_title("CPL sensitivity to the chosen domain")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("DESI DR2 BAO screening diagnostics", fontsize=12)
    base = Path(args.output)
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".png"), dpi=180, metadata={"Title": "DESI DR2 BAO screening diagnostics"})
    fig.savefig(base.with_suffix(".svg"), metadata={"Title": "DESI DR2 BAO screening diagnostics"})
    print(json.dumps({"png": str(base.with_suffix('.png')),
                      "svg": str(base.with_suffix('.svg')),
                      "n_rows": len(data.z),
                      "fit_covariance": "full; residual display uses diagonal marginal standard deviations"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

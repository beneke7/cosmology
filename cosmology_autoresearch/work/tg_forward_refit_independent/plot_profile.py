#!/usr/bin/env python3
"""Render the verified NGC 3198 isolated-surrogate profile and rotation curve."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RESULT = HERE / "profile_checkpointed.json"
SURFACE = HERE / "profile_checkpointed_surface.csv"
OUT = HERE / "profile_summary.png"


def main():
    result = json.loads(RESULT.read_text())
    data = result["observed_rotmod_data"]
    r = np.asarray(data["r_kpc"])
    v = np.asarray(data["vobs_kms"])
    sigma = np.asarray(data["sigma_kms"])

    fig, (ax_curve, ax_profile) = plt.subplots(1, 2, figsize=(12.5, 5.1))
    ax_curve.errorbar(r, v, yerr=sigma, fmt="o", ms=3.2, color="0.25",
                      ecolor="0.65", capsize=1.5, label="SPARC Rotmod data")
    colors = {"1.0": "#0072B2", "1.5": "#D55E00", "2.0": "#009E73"}
    for p, entry in result["profiled_grid_minima"].items():
        model = entry.get("tight_model_v_kms")
        candidate = entry.get("selected_grid_point", {})
        if model is None:
            continue
        ax_curve.plot(r, model, "-", lw=1.6, color=colors.get(p),
                      label=fr"$p={p}$, $K={candidate['K_s2_km2']:.3g}$")
    ax_curve.set(xlabel="$r$ [kpc]", ylabel="$V_c$ [km s$^{-1}$]",
                 title="Isolated spherical-surrogate profile")
    ax_curve.grid(alpha=0.2)
    ax_curve.legend(fontsize=8, frameon=False)

    with SURFACE.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    for p, color in colors.items():
        selected = [row for row in rows
                    if float(row["p_inner"]) == float(p)
                    and abs(float(row["Upsilon"]) - 0.762) < 1e-10
                    and row["valid"] == "True"]
        if not selected:
            continue
        fraction = np.asarray([float(row["fraction"]) for row in selected])
        chi2 = np.asarray([float(row["chi2"]) for row in selected])
        ax_profile.plot(fraction, chi2, "o-", ms=3.5, lw=1.3, color=color,
                        label=fr"$p={p}$")
    ax_profile.axvline(1.0, ls="--", lw=1.0, color="0.25", label="first branch threshold")
    ax_profile.axvline(0.97, ls=":", lw=1.0, color="0.5", label="sampled minimum")
    ax_profile.set_yscale("log")
    ax_profile.set(xlim=(-0.025, 1.04), xlabel=r"$K/K_{\rm crit}$",
                   ylabel=r"direct velocity $\chi^2$",
                   title=r"$\Upsilon=0.762$ direct velocity score")
    ax_profile.grid(alpha=0.2, which="both")
    ax_profile.legend(fontsize=8, frameon=False)
    fig.suptitle("NGC 3198: conditional signed spherical pseudo-density; not a unique 3-D source")
    fig.tight_layout()
    fig.savefig(OUT, dpi=170, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

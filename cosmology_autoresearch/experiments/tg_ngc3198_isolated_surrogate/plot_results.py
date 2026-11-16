#!/usr/bin/env python3
"""Regenerate the NGC 3198 isolated-surrogate and density diagnostic plot."""
from __future__ import annotations

import io
import json
from pathlib import Path
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
with zipfile.ZipFile(ROOT / "context/data/Rotmod_LTG.zip") as archive:
    values = np.loadtxt(io.BytesIO(archive.read("NGC3198_rotmod.dat")), comments="#")
record = json.loads((HERE / "record.json").read_text())
r, vobs, errv, vgas, vdisk, _ = values[:, :6].T
reference = record["reported_parameter_point"]
best = record["local_parameter_grid"]["best_local_point"]
rho_diagnostic = record["source_interpolation_taper_sensitivity"][1]

# Recompute Eq. (47) PCHIP density for the lower panel; preserve signed Vgas.
disk_i = PchipInterpolator(r, vdisk)
gas_i = PchipInterpolator(r, vgas)
dense_r = np.linspace(r[0], r[-1], 3000)
dv, ddv = disk_i(dense_r), disk_i.derivative()(dense_r)
gv, dgv = gas_i(dense_r), gas_i.derivative()(dense_r)
upsilon = record["declared_source_and_boundary"]["Upsilon_disk"]
gkpc = 4.30091e-6
rho = (
    upsilon * (dv**2 + 2.0 * dense_r * dv * ddv)
    + gv * np.abs(gv) + 2.0 * dense_r * np.abs(gv) * dgv
) / (4.0 * np.pi * gkpc * dense_r**2)

fig, (ax, axrho) = plt.subplots(
    2, 1, figsize=(8.3, 8.1), sharex=True,
    gridspec_kw={"height_ratios": [1.55, 1.0]}, constrained_layout=True,
)
ax.errorbar(r, vobs, yerr=errv, fmt="o", ms=3.3, lw=0.8, capsize=1.5,
            color="#263f73", label="SPARC Rotmod observations")
ax.plot(r, reference["predicted_vc_km_s"], lw=2.1, color="#c23b32",
        label=(r"isolated, regular center: $K=3.40\times10^{-5}$, "
               r"$\Upsilon_*=0.762$"))
ax.plot(r, best["predicted_vc_km_s"], lw=1.8, ls="--", color="#dc9224",
        label=(r"best on local 1-$\sigma$ grid: $K=3.44\times10^{-5}$, "
               r"$\Upsilon_*=0.769$"))
ax.set_ylabel(r"$v_c$ [km s$^{-1}$]")
ax.legend(frameon=False, fontsize=8.5, loc="upper left")
ax.text(
    0.99, 0.06,
    "conditional spherical surrogate; not a 3D prediction",
    transform=ax.transAxes, ha="right", fontsize=8.5,
)

axrho.axhline(0.0, color="black", lw=0.9)
axrho.plot(dense_r, rho, color="#53458b", lw=1.7)
axrho.fill_between(dense_r, rho, 0.0, where=rho < 0.0,
                   color="#c23b32", alpha=0.22, interpolate=True,
                   label="negative inferred pseudo-density; retained")
axrho.set_ylabel(r"$\rho_{\rm pseudo}$ [$M_\odot$ kpc$^{-3}$]")
axrho.set_xlabel("radius [kpc]")
axrho.set_xlim(0, r[-1])
axrho.set_ylim(min(np.min(rho) * 1.3, -1.0), np.max(rho) * 1.05)
axrho.legend(frameon=False, fontsize=8.5, loc="lower left")
fig.suptitle("NGC 3198: isolated-boundary TG stress test and source audit")
output = HERE / "rotation_curve_and_source.png"
fig.savefig(output, dpi=180)
plt.close(fig)
print(output)

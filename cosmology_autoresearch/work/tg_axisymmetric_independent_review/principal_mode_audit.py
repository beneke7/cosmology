#!/usr/bin/env python3
"""Fine-grid generalized-eigenvalue and principal-mode audit.

Uses only the task-owned independent cylindrical FV assembly in
audit_axisymmetric.py; no production BVP or branch-edge code is imported.
"""
from __future__ import annotations

import hashlib
import io
import json
import time
import zipfile
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

from audit_axisymmetric import (DR, DZ, RMAX, ZMAX, G, UPSILON, PILOT,
                                source_on_grid, assemble)

HERE = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    t0w, t0c = time.monotonic(), time.process_time()
    source_path = PILOT / "source_profiles.csv"
    expected = "2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5"
    if sha(source_path) != expected:
        raise SystemExit("source profile hash mismatch; refusing calculation")
    dr, dz, rmax, zmax = 0.25, 0.125, 80.0, 60.0
    rc = (np.arange(round(rmax / dr)) + 0.5) * dr
    zc = (np.arange(round(zmax / dz)) + 0.5) * dz
    rho, rd = source_on_grid(rc, zc)
    B, D, volume = assemble(rho, rc, dr, dz, rmax, zmax)
    # Standard shift-invert generalized solve: nearest eigenvalue to sigma=0
    # is the smallest positive Kcrit since B and diag(D) are positive.
    vals, vecs = eigsh(B, k=1, M=diags(D, format="csr"), sigma=0.0,
                       which="LM", tol=2e-8, maxiter=4000)
    kcrit = float(vals[0])
    h = vecs[:, 0]
    if h.sum() < 0:
        h = -h
    h /= np.max(np.abs(h))
    residual = B @ h - kcrit * D * h
    residual_rel = float(np.linalg.norm(residual) / np.linalg.norm(B @ h))
    h2 = h.reshape(rho.shape)
    # Midplane extrapolation from even reflection about z=0.
    hmid = (9.0 * h2[:, 0] - h2[:, 1]) / 8.0
    rdat = np.loadtxt(io.StringIO(zipfile.ZipFile(
        Path(__file__).resolve().parents[2] / "context/data/Rotmod_LTG.zip"
    ).read("NGC3198_rotmod.dat").decode()), comments="#")
    radii, vobs, errv = rdat[:, 0], rdat[:, 1], rdat[:, 2]
    mode = PchipInterpolator(rc, hmid, extrapolate=True)
    h_eval = mode(radii)
    dh_eval = mode.derivative()(radii)
    v2 = -radii * dh_eval / (kcrit * h_eval)
    vmode = np.sqrt(np.maximum(v2, 0.0))
    chi2 = float(np.sum(((vobs - vmode) / errv) ** 2))
    # Source and operator checks: D is strictly positive for the frozen stellar
    # exponential tail; B is the symmetric positive exterior-regular operator.
    symmetry = float(np.max(np.abs((B - B.T).data))) if (B - B.T).nnz else 0.0
    out = {
        "status": "complete_principal_generalized_eigenmode",
        "definition": "smallest positive generalized eigenpair B h=Kcrit diag(D) h, where B is independently assembled integrated cylindrical -Laplacian with monopole Robin and D=4 pi G rho dV",
        "source_csv_sha256": sha(source_path),
        "independent_assembly_sha256": sha(HERE / "audit_axisymmetric.py"),
        "grid": {"dr_kpc": dr, "dz_kpc": dz, "rmax_kpc": rmax, "zmax_kpc": zmax,
                 "nr": len(rc), "nz": len(zc), "cells": int(rho.size)},
        "source": {"upsilon": UPSILON, "G": G, "rho_min_Msun_kpc3": float(rho.min()),
                   "rho_max_Msun_kpc3": float(rho.max()), "stellar_tail_Rd_kpc": rd,
                   "D_min": float(D.min()), "D_max": float(D.max())},
        "eigensolver": {"method": "scipy.sparse.linalg.eigsh", "k": 1, "M": "diag(D)",
                        "sigma": 0.0, "which": "LM", "tol": 2e-8, "maxiter": 4000,
                        "Kcrit_km_s_minus2": kcrit,
                        "reference_production_Kcrit": 5.3714678825e-5,
                        "difference_from_production": kcrit - 5.3714678825e-5,
                        "fractional_difference_from_production": kcrit / 5.3714678825e-5 - 1.0,
                        "generalized_residual_norm_over_Bh": residual_rel,
                        "eigenvector_min_after_positive_sign": float(h.min()),
                        "eigenvector_max": float(h.max()),
                        "B_symmetry_max_abs": symmetry},
        "limiting_principal_mode": {
            "meaning": "As K approaches Kcrit from below, w=(B-KD)^-1 D is dominated by this positive eigenmode. The reported mode curve is v^2=-R/Kcrit*dR(log h_mid). It is a branch-limit shape, not a finite-u solution or a verified physical solution.",
            "midplane_h_min": float(hmid.min()), "midplane_h_max": float(hmid.max()),
            "minimum_v2_km_s2": float(v2.min()),
            "diagonal_chi2_43": chi2,
            "radii_kpc": radii.tolist(), "observed_km_s": vobs.tolist(),
            "sigma_km_s": errv.tolist(), "mode_curve_km_s": vmode.tolist(),
            "mode_shape_hmid_normalized": hmid.tolist()},
        "cpu_seconds": time.process_time() - t0c,
        "wall_seconds": time.monotonic() - t0w,
    }
    path = HERE / "principal_mode_result.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"Kcrit": kcrit, "production_fractional_difference": out["eigensolver"]["fractional_difference_from_production"],
                      "residual": residual_rel, "h_min": float(h.min()),
                      "mode_speed_min_max": [float(vmode.min()), float(vmode.max())],
                      "chi2": chi2, "cpu_seconds": out["cpu_seconds"], "wall_seconds": out["wall_seconds"]}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent coarse FV check of the NGC 3198 transformed TG BVP.

This file deliberately builds its own cell-integrated operator and source;
it does not import the pilot solver or its build/solve routines.
"""
from __future__ import annotations

import csv
import argparse
import hashlib
import json
import math
import os
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import cg, eigsh

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent / "tg_axisymmetric_pilot"
G = 4.30091e-6
K = 3.9e-5
UPSILON = 0.762
DR, DZ, RMAX, ZMAX = 0.5, 0.25, 80.0, 60.0


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def source_on_grid(rc, zc):
    # Independent CSV parsing/interpolation and vertical source construction.
    with (PILOT / "source_profiles.csv").open(newline="") as f:
        p = list(csv.DictReader(f))
    rp = np.array([float(x["r_mid_kpc"]) for x in p])
    st = np.array([float(x["s4g_stars_faceon_Msun_pc2_Upsilon0p6"])
                   if x["s4g_stars_faceon_Msun_pc2_Upsilon0p6"] else np.nan
                   for x in p])
    ga = np.array([float(x["halogas_hihe_faceon_Msun_pc2"])
                   if x["halogas_hihe_faceon_Msun_pc2"] else np.nan for x in p])
    st[rp > 13.75] = np.nan
    ok = np.isfinite(st)
    rs, ss = rp[ok], st[ok]
    surf_s = np.empty_like(rc)
    measured = rc <= rs[-1]
    surf_s[measured] = PchipInterpolator(rs, ss)(np.maximum(rc[measured], rs[0]))
    tf = ok & (rp >= 5.0) & (rp <= 12.75) & (st > 0)
    slope, intercept = np.polyfit(rp[tf], np.log(st[tf]), 1)
    rd = -1.0 / slope
    surf_s[~measured] = ss[-1] * np.exp(-(rc[~measured] - rs[-1]) / rd)
    goodg = np.isfinite(ga)
    rg, sg = rp[goodg], ga[goodg]
    surf_g = np.zeros_like(rc)
    inside = rc <= rg[-1]
    surf_g[inside] = PchipInterpolator(rg, sg)(np.maximum(rc[inside], rg[0]))
    taper = (rc > rg[-1]) & (rc < rg[-1] + 5.0)
    surf_g[taper] = sg[-1] * (1.0 - (rc[taper] - rg[-1]) / 5.0)
    surf_s = np.maximum(surf_s, 0.0) * (UPSILON / 0.6)
    surf_g = np.maximum(surf_g, 0.0)
    zz = zc[None, :]
    # Surface density is Msun/pc^2; divide by 2h and convert pc^-2 to kpc^-2.
    rho_star = surf_s[:, None] * 1.0e6 / (2.0 * 0.3) * np.exp(-zz / 0.3)
    rho_gas = surf_g[:, None] * 1.0e6 * (
        0.85 / (2.0 * 0.2) * np.exp(-zz / 0.2)
        + 0.15 / (2.0 * 3.0) * np.exp(-zz / 3.0))
    return rho_star + rho_gas, float(rd)


def assemble(rho, rc, dr, dz, rmax, zmax):
    nr, nz = rho.shape
    n = nr * nz
    vol = (rc[:, None] * dr * dz) * np.ones((1, nz))
    diag = np.zeros((nr, nz))
    rows, cols, vals = [], [], []

    def links(a, b, c):
        np.add.at(diag.ravel(), a, c)
        np.add.at(diag.ravel(), b, c)
        rows.extend(np.r_[a, b].tolist())
        cols.extend(np.r_[b, a].tolist())
        vals.extend(np.r_[-c, -c].tolist())

    # Radial faces, including exact zero flux at R=0 by having no inner link.
    ii, jj = np.meshgrid(np.arange(nr - 1), np.arange(nz), indexing="ij")
    aa = (ii * nz + jj).ravel()
    bb = ((ii + 1) * nz + jj).ravel()
    c = (((ii + 1) * dr) * dz / dr).ravel()
    links(aa, bb, c)
    # Vertical internal faces; reflection at z=0 means no lower-face link.
    ii, jj = np.meshgrid(np.arange(nr), np.arange(nz - 1), indexing="ij")
    aa = (ii * nz + jj).ravel()
    bb = (ii * nz + jj + 1).ravel()
    c = (rc[:, None] * dr / dz * np.ones((1, nz - 1))).ravel()
    links(aa, bb, c)
    # Monopole exterior on w=(u-1)/K. Integrate Robin flux over the *boundary
    # face* with half-cell resistance: dn w=-kappa*w_b, w_b=w_c/(1+kappa*h/2).
    zc = (np.arange(nz) + 0.5) * dz
    rb = np.hypot(rmax, zc)
    kap = rmax / rb**2
    diag[-1, :] += kap / (1.0 + kap * dr / 2.0) * rmax * dz
    rb = np.hypot(rc, zmax)
    kap = zmax / rb**2
    diag[:, -1] += kap / (1.0 + kap * dz / 2.0) * rc * dr
    rows.extend(np.arange(n).tolist())
    cols.extend(np.arange(n).tolist())
    vals.extend(diag.ravel().tolist())
    B = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    D = (4.0 * math.pi * G * rho * vol).ravel()
    return B, D, vol


def solve(B, D, k):
    A = B - diags(k * D, format="csr")
    x, info = cg(A, D, rtol=1e-10, atol=0.0, maxiter=16000)
    res = A @ x - D
    return x, {"cg_info": int(info), "iterations_residual_l2_over_rhs":
               float(np.linalg.norm(res) / np.linalg.norm(D)),
               "residual_max_abs": float(np.max(np.abs(res)))}


def velocity(rc, w, radii, k):
    # Cell-centered z values; even symmetry gives a quadratic midplane estimate.
    wm = (9.0 * w[:, 0] - w[:, 1]) / 8.0
    f = PchipInterpolator(rc, wm, extrapolate=True)
    ww, dww = f(radii), f.derivative()(radii)
    vv2 = -radii * dww / (1.0 + k * ww)
    return np.sqrt(vv2), wm


def main():
    global DR, DZ, RMAX, ZMAX
    ap = argparse.ArgumentParser()
    ap.add_argument("--dr", type=float, default=DR)
    ap.add_argument("--dz", type=float, default=DZ)
    ap.add_argument("--rmax", type=float, default=RMAX)
    ap.add_argument("--zmax", type=float, default=ZMAX)
    ap.add_argument("--result", type=Path, default=HERE / "independent_axisymmetric_result.json")
    args = ap.parse_args()
    DR, DZ, RMAX, ZMAX = args.dr, args.dz, args.rmax, args.zmax
    t0, c0 = time.monotonic(), time.process_time()
    profile = PILOT / "source_profiles.csv"
    baseline = PILOT / "newtonian_k0_baseline.json"
    candidate_path = PILOT / "scan_u0762" / "axisymmetric_candidate_points.jsonl"
    summary_path = PILOT / "scan_u0762" / "axisymmetric_pilot_summary.json"
    expected = "2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5"
    if sha(profile) != expected:
        raise SystemExit("SOURCE HASH MISMATCH; audit refused")
    candidate = None
    with candidate_path.open() as f:
        for line in f:
            rec = json.loads(line)
            if abs(float(rec.get("K", -1)) - K) < 1e-14 and rec.get("status") == "ok":
                candidate = rec
                break
    if candidate is None:
        raise SystemExit("K=3.9e-5 production candidate missing")
    with summary_path.open() as f:
        production_summary = json.load(f)
    if production_summary["source_csv_sha256"] != expected:
        raise SystemExit("production summary source hash mismatch")
    import zipfile, io
    with zipfile.ZipFile(Path(__file__).resolve().parents[2] / "context/data/Rotmod_LTG.zip") as z:
        dat = np.loadtxt(io.StringIO(z.read("NGC3198_rotmod.dat").decode()), comments="#")
    radii = dat[:, 0]
    rc = (np.arange(round(RMAX / DR)) + 0.5) * DR
    zc = (np.arange(round(ZMAX / DZ)) + 0.5) * DZ
    rho, rd = source_on_grid(rc, zc)
    B, D, vol = assemble(rho, rc, DR, DZ, RMAX, ZMAX)
    # Principal generalized eigenvalue B h = Kcrit D h. Shift-invert at zero
    # identifies the smallest positive Kcrit rather than using CG breakdown.
    evals, evecs = eigsh(B, k=1, M=diags(D, format="csr"), sigma=0.0,
                         which="LM", tol=2e-8, maxiter=4000)
    kcrit = float(evals[0])
    h = evecs[:, 0]
    eig_res = float(np.linalg.norm(B @ h - kcrit * D * h) / np.linalg.norm(B @ h))
    # Mandatory nonlinear point and K=0 linear-reference solve.
    x, rn = solve(B, D, K)
    wk = x.reshape(rho.shape)
    u = 1.0 + K * wk
    vind, wm = velocity(rc, wk, radii, K)
    vprod = np.asarray(candidate["speeds_km_s"])
    x0, r0n = solve(B, D, 0.0)
    w0 = x0.reshape(rho.shape)
    v0, wm0 = velocity(rc, w0, radii, 0.0)
    bdat = np.genfromtxt(PILOT / "newtonian_k0_predictions.csv", delimiter=",", names=True)
    vstar06 = np.interp(radii, bdat["R_kpc"], bdat["Vstar_map_K0_km_s"])
    vgas = np.interp(radii, bdat["R_kpc"], bdat["Vgas_map_HALOGAS_thin_thick_K0_km_s"])
    # Compare at the requested Upsilon=0.762: baseline stars are normalized
    # to 0.6, while gas is unchanged; Newtonian component speeds add in V^2.
    vbase = np.sqrt((UPSILON / 0.6) * vstar06**2 + vgas**2)
    # Robin face mismatch using recovered boundary-face w values.
    wouter = wk[-1, :]
    zcent = zc
    kapR = RMAX / (RMAX**2 + zcent**2)
    wbR = wouter / (1.0 + kapR * DR / 2.0)
    dnR = -kapR * wbR
    gradR = (wbR - wouter) / (DR / 2.0)
    wtop = wk[:, -1]
    kapZ = ZMAX / (rc**2 + ZMAX**2)
    wbZ = wtop / (1.0 + kapZ * DZ / 2.0)
    dnZ = -kapZ * wbZ
    gradZ = (wbZ - wtop) / (DZ / 2.0)
    robin_rel = max(np.max(np.abs(gradR - dnR) / np.maximum(np.abs(dnR), 1e-30)),
                    np.max(np.abs(gradZ - dnZ) / np.maximum(np.abs(dnZ), 1e-30)))
    result = {
        "status": "complete", "source_csv_sha256": sha(profile),
        "baseline_json_sha256": sha(baseline), "candidate_jsonl_sha256": sha(candidate_path),
        "candidate_point": {"K": K, "upsilon": UPSILON, "production_grid": candidate["grid"],
                            "production_relative_residual": candidate["relative_residual"]},
        "independent_method": "cell-centered cylindrical FV; independently assembled conductance graph; even midplane reflection; regular zero-flux axis; face-integrated monopole Robin with half-cell boundary resistance; CG",
        "grid": {"dr_kpc": DR, "dz_kpc": DZ, "rmax_kpc": RMAX, "zmax_kpc": ZMAX,
                 "nr": len(rc), "nz": len(zc), "unknowns": int(rho.size)},
        "source": {"G_kpc_km2_s2_Msun": G, "Sigma_Msun_pc2_to_Msun_kpc2": 1e6,
                   "stellar_Upsilon": UPSILON, "stellar_h_kpc": 0.3,
                   "gas_HIHe": "HALOGAS; 85% h=.2 kpc + 15% h=3 kpc",
                   "stellar_tail_Rd_kpc": rd, "rho_min": float(rho.min()), "rho_max": float(rho.max()),
                   "coefficient_units": "4 pi G K rho = kpc^-2; with rho Msun/kpc^3 and G kpc (km/s)^2/Msun, K (km/s)^-2"},
        "K_point": {"solver": rn, "u_min_full_grid": float(u.min()), "u_max_full_grid": float(u.max()),
                    "speed_max_abs_difference_km_s": float(np.max(np.abs(vind-vprod))),
                    "speed_rms_difference_km_s": float(np.sqrt(np.mean((vind-vprod)**2))),
                    "speed_max_fractional_difference": float(np.max(np.abs(vind-vprod)/vprod)),
                    "production_speeds_min_max": [float(vprod.min()), float(vprod.max())],
                    "independent_speeds_min_max": [float(vind.min()), float(vind.max())],
                    "robin_relative_flux_residual": float(robin_rel)},
        "principal_eigenvalue": {"definition": "smallest positive generalized eigenvalue Kcrit of B h=Kcrit diag(D) h; B is independent integrated -Laplacian with monopole Robin, D=4 pi G rho times cylindrical cell volume",
                                 "Kcrit_km_s^-2": kcrit, "distance_K_over_Kcrit": K/kcrit,
                                 "generalized_eigen_residual": eig_res},
        "K0_limit": {"solver": r0n, "speed_max_abs_difference_from_Hankel_km_s": float(np.max(np.abs(v0-vbase))),
                     "speed_rms_difference_from_Hankel_km_s": float(np.sqrt(np.mean((v0-vbase)**2))),
                     "speed_max_fractional_difference_from_Hankel": float(np.max(np.abs(v0-vbase)/vbase)),
                     "comment": "Hankel stellar V^2 was rescaled from Upsilon=0.6 to 0.762 and gas V^2 held fixed. Difference includes grid/source quadrature, finite box and midplane/radial interpolation errors; it is not attributed solely to the TG discretization."},
        "radii_kpc": radii.tolist(), "production_speed_km_s": vprod.tolist(),
        "independent_speed_km_s": vind.tolist(), "K0_independent_speed_km_s": v0.tolist(),
        "K0_hankel_speed_km_s": vbase.tolist(),
        "cpu_seconds": time.process_time()-c0, "wall_seconds": time.monotonic()-t0,
        "thread_env": {x: os.environ.get(x) for x in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]}
    }
    out = args.result if args.result.is_absolute() else HERE / args.result
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ["grid", "K_point", "K0_limit", "cpu_seconds", "wall_seconds"]}, indent=2))


if __name__ == "__main__":
    main()

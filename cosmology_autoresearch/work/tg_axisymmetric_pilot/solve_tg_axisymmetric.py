#!/usr/bin/env python3
"""Finite-volume NGC 3198 pilot for the stationary thermodynamic-gravity field.

Solves (Delta + 4 pi G K rho) u = 0 with u -> 1 at infinity, using the exact
transformation u=exp(-K phi).  To avoid cancellation at small K, solve for
w=(u-1)/K instead: (B-KD)w=D, where B is the positive discrete -Laplacian
with monopole Robin exterior and D is the cell-integrated 4 pi G rho.

This is an axisymmetrized parametric disk-source pilot, not a direct 3D source
reconstruction and not a physical detection. See SOURCE_AUDIT.md for provenance
and declared source limitations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import cg

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PROFILE_PATH = OUT / "source_profiles.csv"
K0_PATH = OUT / "newtonian_k0_predictions.csv"
ROT_MOD = ROOT / "context/data/Rotmod_LTG.zip"
G = 4.30091e-6  # kpc (km/s)^2 / Msun
STAR_HZ = 0.3  # kpc
STAR_R_MAX = 13.75  # kpc; last fully covered measured annulus
GAS_HZ_THIN = 0.2  # kpc
GAS_HZ_THICK = 3.0  # kpc
GAS_THIN_FRACTION = 0.85
GAS_TAPER = 5.0  # kpc beyond final profile point


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_profile():
    with PROFILE_PATH.open(newline="") as stream:
        data = list(csv.DictReader(stream))
    r = np.array([float(x["r_mid_kpc"]) for x in data])
    def col(name):
        return np.array([float(x[name]) if x[name] else np.nan for x in data])
    stars = col("s4g_stars_faceon_Msun_pc2_Upsilon0p6")
    stars[r > STAR_R_MAX] = np.nan
    return r, stars, col("halogas_hihe_faceon_Msun_pc2"), col("things_hihe_faceon_Msun_pc2")


def source_profiles(rgrid, r, stars, gas, gas_alt, upsilon, *, star_tail_factor=1.0,
                    gas_name="halogas"):
    valid = np.isfinite(stars)
    rs, ss = r[valid], stars[valid]
    star_fn = PchipInterpolator(rs, ss, extrapolate=False)
    outstar = np.empty_like(rgrid)
    inside = rgrid <= rs[-1]
    outstar[inside] = star_fn(np.maximum(rgrid[inside], rs[0])) * (upsilon / 0.6)
    tail_fit = valid & (r >= 5.0) & (r <= 12.75) & (stars > 0)
    slope, intercept = np.polyfit(r[tail_fit], np.log(stars[tail_fit]), 1)
    if slope >= 0:
        raise ValueError("S4G measured stellar tail is not declining")
    rd = -1.0 / slope * star_tail_factor
    outstar[~inside] = ss[-1] * (upsilon / 0.6) * np.exp(-(rgrid[~inside] - rs[-1]) / rd)

    sigmagas = gas if gas_name == "halogas" else gas_alt
    validg = np.isfinite(sigmagas)
    rg, sg = r[validg], sigmagas[validg]
    gas_fn = PchipInterpolator(rg, sg, extrapolate=False)
    outgas = np.zeros_like(rgrid)
    ingas = rgrid <= rg[-1]
    outgas[ingas] = gas_fn(np.maximum(rgrid[ingas], rg[0]))
    taper = (rgrid > rg[-1]) & (rgrid < rg[-1] + GAS_TAPER)
    outgas[taper] = sg[-1] * (1.0 - (rgrid[taper] - rg[-1]) / GAS_TAPER)
    return np.maximum(outstar, 0), np.maximum(outgas, 0), {
        "star_tail_Rd_kpc": float(rd), "star_fit_RMS_log": float(np.sqrt(np.mean(
            (np.log(stars[tail_fit]) - (slope * r[tail_fit] + intercept))**2))),
        "star_fit_n": int(tail_fit.sum()), "gas_source": gas_name,
    }


def read_rotmod():
    import io
    import zipfile
    with zipfile.ZipFile(ROT_MOD) as zf:
        raw = zf.read("NGC3198_rotmod.dat")
    arr = np.loadtxt(io.BytesIO(raw), comments="#")
    # r, Vobs, sigmaV, Vgas, Vdisk, Vbul
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]


def build_operator(dr, dz, rmax, zmax, upsilon, *, gas_name="halogas", tail_factor=1.0):
    nr, nz = int(round(rmax / dr)), int(round(zmax / dz))
    dr_eff, dz_eff = rmax / nr, zmax / nz
    rc = (np.arange(nr) + 0.5) * dr_eff
    zc = (np.arange(nz) + 0.5) * dz_eff
    rr, zz = np.meshgrid(rc, zc, indexing="ij")
    rprof, stars, gas, gas_alt = load_profile()
    star_surf, gas_surf, src_cfg = source_profiles(rc, rprof, stars, gas, gas_alt, upsilon,
                                                   star_tail_factor=tail_factor,
                                                   gas_name=gas_name)
    rho = (star_surf[:, None] * 1.0e6 / (2.0 * STAR_HZ) * np.exp(-zz / STAR_HZ)
           + gas_surf[:, None] * 1.0e6 * (
               GAS_THIN_FRACTION / (2 * GAS_HZ_THIN) * np.exp(-zz / GAS_HZ_THIN)
               + (1 - GAS_THIN_FRACTION) / (2 * GAS_HZ_THICK) * np.exp(-zz / GAS_HZ_THICK)))

    # Cylindrical finite-volume weights omitting the common 2 pi factor:
    # volume=R_i dR dz; radial-face conductance=R_face dz/dR; vertical-face
    # conductance=R_i dR/dz. The inner axis has exactly zero radial flux.
    vol = rr * dr_eff * dz_eff
    n = nr * nz
    diagonal = np.zeros((nr, nz), dtype=np.float64)
    rows, cols, vals = [], [], []

    def add_link(i, j, ni, nj, conductance):
        a = i * nz + j
        b = ni * nz + nj
        diagonal[i, j] += conductance
        diagonal[ni, nj] += conductance
        rows.extend((a, b))
        cols.extend((b, a))
        vals.extend((-conductance, -conductance))

    for i in range(nr - 1):
        gface = (i + 1) * dr_eff * dz_eff / dr_eff
        for j in range(nz):
            add_link(i, j, i + 1, j, gface)
    for i in range(nr):
        gface = rc[i] * dr_eff / dz_eff
        for j in range(nz - 1):
            add_link(i, j, i, j + 1, gface)

    # At the symmetry plane z=0, Neumann flux is zero. On the two distant box
    # faces use the isolated monopole Robin operator: dn w + (n.r)/r^2 w=0.
    # It is an explicit exterior approximation, checked by enlarging the box.
    for j in range(nz):
        radius = math.hypot(rmax, zc[j])
        kappa = rmax / radius**2
        diagonal[-1, j] += kappa * rc[-1] * dz_eff
    for i in range(nr):
        radius = math.hypot(rc[i], zmax)
        kappa = zmax / radius**2
        diagonal[i, -1] += kappa * rc[i] * dr_eff

    rows.extend(range(n)); cols.extend(range(n)); vals.extend(diagonal.ravel())
    B = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))),
                   shape=(n, n)).tocsr()
    D = (4.0 * math.pi * G * rho * vol).ravel()
    return B, D, rho, vol, rr, zz, (dr_eff, dz_eff), src_cfg


def solve_point(B, D, shape, step, K, radii):
    # A = B - K diag(D/volume * volume) is the cell-integrated
    # ( -Delta - 4 pi G K rho ) operator. Right-hand side D has units
    # (km/s)^2 kpc; w=(u-1)/K has units (km/s)^2.
    A = B - diags(K * D, format="csr")
    rhs = D
    diagonal = A.diagonal()
    if np.any(diagonal <= 0):
        return {"status": "nonpositive_diagonal", "K": K}
    invdiag = 1.0 / diagonal
    M = diags(invdiag, format="csr")
    iteration = [0]
    def callback(_):
        iteration[0] += 1
    wvec, info = cg(A, rhs, M=M, rtol=2e-9, atol=0.0, maxiter=12000, callback=callback)
    residual = float(np.linalg.norm(A @ wvec - rhs) / np.linalg.norm(rhs))
    if info != 0 or not np.all(np.isfinite(wvec)):
        return {"status": "cg_failed", "K": K, "cg_info": int(info),
                "iterations": iteration[0], "relative_residual": residual}
    w = wvec.reshape(shape)
    dr, dz = step
    # Even reflection about the midplane; second-order quadratic extrapolation
    # from z=dz/2 and 3dz/2 removes the O(dz^2) value bias.
    wmid = (9.0 * w[:, 0] - w[:, 1]) / 8.0
    rw = (np.arange(shape[0]) + 0.5) * dr
    u_mid = 1.0 + K * wmid
    if np.any(u_mid <= 0):
        return {"status": "u_nonpositive", "K": K, "cg_info": int(info),
                "iterations": iteration[0], "relative_residual": residual,
                "minimum_u_mid": float(u_mid.min())}
    # The innermost SPARC point can lie inside the first cell center. The
    # regular axis condition is zero radial derivative, so a smooth one-sided
    # polynomial continuation to R=0 is the appropriate interpolation here.
    # Differentiate w directly: differentiating u=1+K w and dividing the
    # result by K is algebraically equivalent but loses precision as K->0.
    wfn = PchipInterpolator(rw, wmid, extrapolate=True)
    w_eval = wfn(radii)
    u_eval = 1.0 + K * w_eval
    v2 = -radii * wfn.derivative()(radii) / u_eval
    if np.any(~np.isfinite(v2)) or np.any(v2 <= 0):
        return {"status": "nonphysical_v2", "K": K, "cg_info": int(info),
                "iterations": iteration[0], "relative_residual": residual,
                "minimum_u_mid": float(u_mid.min()), "minimum_v2": float(np.nanmin(v2))}
    v = np.sqrt(v2)
    _, vobs, err, _, _ = read_rotmod()
    return {"status": "ok", "K": float(K), "iterations": int(iteration[0]),
            "relative_residual": residual, "minimum_u_mid": float(u_mid.min()),
            "maximum_u_mid": float(u_mid.max()),
            "chi2_diagonal_43": float(np.sum(((vobs - v) / err)**2)),
            "rms_resid_km_s": float(np.sqrt(np.mean((vobs - v)**2))),
            "speeds_km_s": v.tolist(), "u_mid": u_mid.tolist(),
            "r_grid_kpc": rw.tolist(), "w_mid": wmid.tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dr", type=float, default=0.5, help="radial grid spacing in kpc")
    parser.add_argument("--dz", type=float, default=0.25, help="vertical grid spacing in kpc")
    parser.add_argument("--rmax", type=float, default=80.0, help="outer cylindrical radius in kpc")
    parser.add_argument("--zmax", type=float, default=60.0, help="box half-height in kpc")
    parser.add_argument("--cpu-cap", type=float, default=1100.0, help="stop before 1200 aggregate CPU s")
    parser.add_argument("--upsilon", type=float, action="append",
                        help="stellar M/L to scan (repeatable; defaults to the declared grid)")
    parser.add_argument("--output-dir", type=Path, default=OUT,
                        help="isolated output directory for this bounded scan")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if os.environ.get("OPENBLAS_NUM_THREADS") != "1":
        print("NOTE: for reproducible sparse runs set OPENBLAS_NUM_THREADS=1")
    started_wall, started_cpu = time.perf_counter(), time.process_time()
    points_path = args.output_dir / "axisymmetric_candidate_points.jsonl"
    resume = []
    if args.resume and points_path.exists():
        with points_path.open() as stream:
            resume = [json.loads(line) for line in stream if line.strip()]
    done = {(round(x.get("K_requested", x["K"]), 14), round(x["upsilon"], 6), x["status"])
            for x in resume if "K" in x and "upsilon" in x}
    rprof, stars, gas, gas_alt = load_profile()
    radii, vobs, err, _, _ = read_rotmod()
    k0ref = np.genfromtxt(K0_PATH, delimiter=",", names=True)
    if radii.size != vobs.size or k0ref.size != radii.size:
        raise RuntimeError("SPARC and K=0 baseline row counts disagree")
    source_hash = sha256(PROFILE_PATH)

    # S4G stated stellar-mass uncertainty is about 0.1 dex; include the exact
    # paper's quoted NGC 3198 value (.762) as a labeled reference, not a prior.
    # This first scan stops at the published-parameter neighborhood. The
    # admissibility edge is computed separately from B h=Kcrit D h, not by
    # continuing CG past a singular matrix.
    K_grid = [0.0, 0.5e-5, 1.0e-5, 1.5e-5, 2.0e-5, 2.5e-5, 3.0e-5,
              3.3e-5, 3.5e-5, 3.7e-5, 3.9e-5]
    U_grid = args.upsilon or [10**(-0.1) * 0.6, 0.6, 10**0.1 * 0.6, 0.762]
    results = resume.copy()
    previous_ok = [x for x in results if x.get("status") == "ok" and
                    "chi2_diagonal_43" in x]
    best = min(previous_ok, key=lambda x: x["chi2_diagonal_43"]) if previous_ok else None
    solver_hash = sha256(Path(__file__))
    for upsilon in U_grid:
        if time.process_time() - started_cpu > args.cpu_cap:
            break
        B, D, rho, vol, rr, zz, step, src_cfg = build_operator(
            args.dr, args.dz, args.rmax, args.zmax, upsilon)
        shape = rho.shape
        for K in K_grid:
            if time.process_time() - started_cpu > args.cpu_cap:
                break
            if (round(K, 14), round(upsilon, 6), "ok") in done:
                continue
            if K == 0:
                # Newtonian limit is a separate Poisson solve, rather than a
                # cancellation-prone subtraction in u=1+K w.
                Keval = 1.0e-12
            else:
                Keval = K
            point = solve_point(B, D, shape, step, Keval, radii)
            point.update({"K_requested": float(K), "upsilon": float(upsilon),
                          "gas_source": "halogas", "stellar_tail": "exponential",
                          "grid": {"dr_kpc": step[0], "dz_kpc": step[1],
                                   "rmax_kpc": args.rmax, "zmax_kpc": args.zmax,
                                   "nr": shape[0], "nz": shape[1]},
                          "source_cfg": src_cfg, "source_csv_sha256": source_hash,
                          "solver_sha256": solver_hash,
                          "solver": "cylindrical_cell-centered FV; CG+Jacobi; outer monopole Robin",
                          "cg_rtol": 2e-9})
            results.append(point)
            with points_path.open("a") as stream:
                stream.write(json.dumps(point, sort_keys=True) + "\n")
            if point.get("status") == "ok" and (best is None or point["chi2_diagonal_43"] < best["chi2_diagonal_43"]):
                best = point
                # Compact reproducible field checkpoint for independent checking.
                record = dict(point)
                for key in ("w_mid",):
                    record.pop(key, None)
                (args.output_dir / "axisymmetric_pilot_best.json").write_text(json.dumps(record, indent=2) + "\n")
                pred = {"R_kpc": radii, "Vobs_km_s": vobs, "sigmaV_km_s": err,
                        "V_TG_km_s": point["speeds_km_s"]}
                np.savetxt(args.output_dir / "axisymmetric_pilot_best_curve.csv",
                           np.column_stack(list(pred.values())), delimiter=",",
                           header=",".join(pred), comments="")
                fig, ax = plt.subplots(figsize=(8, 5))
                ax.errorbar(radii, vobs, yerr=err, fmt="o", ms=3, color="black", label="SPARC NGC 3198")
                ax.plot(radii, point["speeds_km_s"], "-", lw=2, label="TG axisymmetric pilot")
                ax.plot(radii, np.sqrt(np.maximum(k0ref["Vbar_map_HALOGAS_K0_km_s"]**2, 0)),
                        "--", label="resolved-source K=0 baryons")
                ax.set(xlabel="R [kpc]", ylabel="circular speed [km/s]",
                       title=f"NGC 3198 pilot: K={K:g} (km/s)^-2, Υ*={upsilon:g}")
                ax.legend()
                ax.grid(alpha=.25)
                fig.tight_layout()
                fig.savefig(args.output_dir / "axisymmetric_pilot_best.png", dpi=170)
                plt.close(fig)
        print(f"completed Upsilon={upsilon:.6g}; CPU={time.process_time()-started_cpu:.1f}s; wall={time.perf_counter()-started_wall:.1f}s", flush=True)

    report = {
        "status": "bounded_scan_complete" if time.process_time() - started_cpu <= args.cpu_cap else "cpu_cap_checkpoint",
        "source_csv": str(PROFILE_PATH.relative_to(ROOT)), "source_csv_sha256": source_hash,
        "solver_sha256": solver_hash,
        "K0_reference_sha256": sha256(K0_PATH), "Rotmod_zip_sha256": sha256(ROT_MOD),
        "python": os.sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
        "grid": {"dr_kpc": args.dr, "dz_kpc": args.dz, "rmax_kpc": args.rmax,
                 "zmax_kpc": args.zmax},
        "K_grid_requested": K_grid, "upsilon_grid_requested": U_grid,
        "boundary": "monopole Robin dn(w)+(n.r)/r^2 w=0 on R=Rmax and z=Zmax; Neumann z=0; regular zero-flux axis",
        "source": "S4G stars inside complete rings then exponential tail; HALOGAS HI+He, .85 h=.2 kpc + .15 h=3 kpc; stellar h=.3 kpc",
        "chi2_scope": "43 SPARC NGC3198 pointwise errors, diagonal descriptive score only; no covariance",
        "positive_u": "checked at interpolated midplane only; candidate independent report must check full field",
        "wall_seconds": time.perf_counter() - started_wall, "cpu_seconds": time.process_time() - started_cpu,
        "completed_points": len(results),
        "ok_points": sum(x.get("status") == "ok" for x in results),
        "best_ok": {k: v for k, v in (best or {}).items() if k not in ("w_mid", "u_mid", "r_grid_kpc")},
    }
    (args.output_dir / "axisymmetric_pilot_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

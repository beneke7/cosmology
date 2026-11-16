#!/usr/bin/env python3
"""Fixed-physical-K NGC 3198 finite-volume refinement check.

Reconstruct the axisymmetrized S4G + HALOGAS source at D=16.6 Mpc, i=72 deg,
Upsilon*=0.762, h*=0.3 kpc and an all-thin h_gas=0.2 kpc profile.  At every
mesh K is held at the supplied physical value; it is never rescaled to Kcrit.
The source and outer monopole Robin boundary are the same family documented
in work/tg_axisymmetric_nuisance_envelope/nuisance_envelope.py and
work/tg_axisymmetric_pilot/solve_tg_axisymmetric.py.

Run from cosmology_autoresearch with one-thread BLAS, a 180 CPU-second and
420 wall-second cap (the exact command is in REPORT.md).
"""
from __future__ import annotations

import csv
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, splu
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PILOT = ROOT / "work/tg_axisymmetric_pilot"
NE_DIR = ROOT / "work/tg_axisymmetric_nuisance_envelope"
SOLVER_PATH = PILOT / "solve_tg_axisymmetric.py"
SOURCE_BUILDER_PATH = NE_DIR / "nuisance_envelope.py"
K_FIXED = 4.398380672762078e-5  # (km/s)^-2, fixed in every solve
D_MPC, INC_DEG, UPSILON, HSTAR = 16.6, 72.0, 0.762, 0.3
D_REF_MPC, I_REF_DEG = 13.8, 73.0
G = 4.30091e-6  # kpc (km/s)^2 / Msun
HGAS = 0.2
CPU_CAP_S = 180.0
WALL_CAP_S = 420.0
CASE_SPECS = [
    {"name": "reference", "dr": 0.25, "dz": 0.125, "rmax": 80.0, "zmax": 60.0},
    {"name": "finer_mesh", "dr": 0.20, "dz": 0.10, "rmax": 80.0, "zmax": 60.0},
    {"name": "half_step_mesh", "dr": 0.125, "dz": 0.0625, "rmax": 80.0, "zmax": 60.0},
    {"name": "larger_domain", "dr": 0.25, "dz": 0.125, "rmax": 100.0, "zmax": 80.0},
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def build_geometry_source(builder, maps, rcent: np.ndarray, zcent: np.ndarray,
                          rmax: float) -> tuple[np.ndarray, dict]:
    profile = builder.transformed_sources(maps, D_MPC, INC_DEG)
    # Evaluate the exact stated P5 exponential stellar continuation and
    # 5-kpc HALOGAS taper across the whole requested radial box. This makes
    # enlarged-domain source support explicit instead of clamping at 50 kpc.
    source_r = np.arange(0.0, rmax + 0.0200001, 0.02)
    stars, star_cfg = builder.surface_density_on_grid(source_r, profile["r"],
                                                       profile["stars"], "stars")
    gas, gas_cfg = builder.surface_density_on_grid(source_r, profile["r"],
                                                    profile["gas"], "gas")
    star_sigma = np.interp(rcent, source_r, stars)
    gas_sigma = np.interp(rcent, source_r, gas)
    rho_stars = (star_sigma[:, None] * 1.0e6 / (2.0 * HSTAR)
                 * np.exp(-zcent[None, :] / HSTAR))
    rho_gas = (gas_sigma[:, None] * 1.0e6 / (2.0 * HGAS)
               * np.exp(-zcent[None, :] / HGAS))
    rho = UPSILON / 0.6 * rho_stars + rho_gas
    return rho, {
        "distance_Mpc": D_MPC, "inclination_deg": INC_DEG,
        "upsilon": UPSILON, "stellar_h_kpc": HSTAR,
        "gas_model": "all thin HI+He exponential, h=0.2 kpc",
        "profile_bin_kpc": 0.5,
        "map_profile_coverage_minimum": 0.8,
        "stellar_profile_support": "P5 ellipse; exponential tail fitted to supported 5-12.75 kpc bins",
        "gas_profile_support": "HALOGAS HI+He; linear 5 kpc outer taper",
        "star_profile": star_cfg, "gas_profile": gas_cfg,
        "source_radial_evaluation_max_kpc": float(source_r[-1]),
    }


def run_case(spec, solver, builder, maps, data):
    wall0, cpu0 = time.perf_counter(), time.process_time()
    B, _ignored_D, _ignored_rho, vol, rr, zz, step, _ignored_cfg = solver.build_operator(
        spec["dr"], spec["dz"], spec["rmax"], spec["zmax"], UPSILON)
    nr, nz = rr.shape
    dr, dz = step
    rcent = (np.arange(nr) + 0.5) * dr
    zcent = (np.arange(nz) + 0.5) * dz
    rho, source_cfg = build_geometry_source(builder, maps, rcent, zcent, spec["rmax"])
    dvec = (4.0 * math.pi * G * rho * vol).ravel()

    # Generalized principal eigenvalue B h = Kcrit diag(D) h, using the same
    # shift-invert method as the parent branch-edge audit. Residual is evaluated
    # against the assembled finite-volume matrices.
    eig_values, eig_vectors = eigsh(B, k=1, M=diags(dvec, format="csr"),
                                    sigma=0.0, which="LM", tol=2e-8,
                                    maxiter=3000)
    kcrit = float(eig_values[0])
    h = eig_vectors[:, 0]
    if np.sum(h) < 0.0:
        h *= -1.0
    eig_resid = float(np.linalg.norm(B @ h - kcrit * dvec * h) /
                      np.linalg.norm(B @ h))
    del eig_vectors
    row = {
        "case": spec["name"], **spec, "nr": nr, "nz": nz,
        "effective_dr_kpc": dr, "effective_dz_kpc": dz,
        "K_fixed_km_s_minus2": K_FIXED, "Kcrit_km_s_minus2": kcrit,
        "K_less_than_Kcrit": bool(K_FIXED < kcrit),
        "K_over_Kcrit": K_FIXED / kcrit,
        "eigen_relative_residual": eig_resid,
        "source_min_rho_Msun_kpc3": float(rho.min()),
        "source_max_rho_Msun_kpc3": float(rho.max()),
        "grid_cells": int(nr * nz), "source": source_cfg,
    }
    if K_FIXED >= kcrit:
        row.update({"status": "supercritical_fixed_K", "minimum_u_full": None,
                    "minimum_u_midplane": None, "linear_relative_residual": None,
                    "chi2_diagonal_43": None, "rms_residual_km_s": None,
                    "predicted_speeds_km_s": None})
    else:
        A = B - diags(K_FIXED * dvec, format="csr")
        lu = splu(A.tocsc())
        wvec = lu.solve(dvec)
        residual = float(np.linalg.norm(A @ wvec - dvec) / np.linalg.norm(dvec))
        w = wvec.reshape((nr, nz))
        u = 1.0 + K_FIXED * w
        wmid = (9.0 * w[:, 0] - w[:, 1]) / 8.0
        umid = (9.0 * u[:, 0] - u[:, 1]) / 8.0
        wfn = PchipInterpolator(rcent, wmid, extrapolate=True)
        radii = data["r"] * (D_MPC / D_REF_MPC)
        v2 = -radii * wfn.derivative()(radii) / (1.0 + K_FIXED * wfn(radii))
        speed = np.sqrt(v2) if np.all(np.isfinite(v2)) and np.all(v2 > 0.0) else None
        if speed is None:
            row.update({"status": "subcritical_but_nonphysical_curve",
                        "minimum_u_full": float(np.min(u)),
                        "maximum_u_full": float(np.max(u)),
                        "minimum_u_midplane": float(np.min(umid)),
                        "linear_relative_residual": residual,
                        "chi2_diagonal_43": None, "rms_residual_km_s": None,
                        "predicted_speeds_km_s": None})
        else:
            obs_factor = math.sin(math.radians(I_REF_DEG)) / math.sin(math.radians(INC_DEG))
            obs, err = data["vobs"] * obs_factor, data["err"] * obs_factor
            resid = obs - speed
            row.update({"status": "ok", "minimum_u_full": float(np.min(u)),
                        "maximum_u_full": float(np.max(u)),
                        "minimum_u_midplane": float(np.min(umid)),
                        "linear_relative_residual": residual,
                        "chi2_diagonal_43": float(np.sum((resid / err) ** 2)),
                        "rms_residual_km_s": float(np.sqrt(np.mean(resid**2))),
                        "predicted_speeds_km_s": speed.tolist(),
                        "observed_velocity_factor_from_i73": obs_factor})
    row["cpu_seconds"] = time.process_time() - cpu0
    row["wall_seconds"] = time.perf_counter() - wall0
    return row


def main():
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=[x["name"] for x in CASE_SPECS],
                        help="run only this independent mesh/domain case")
    parser.add_argument("--output-dir", type=Path, default=HERE,
                        help="write this run's artifacts here")
    args = parser.parse_args()
    out_dir = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    case_specs = [x for x in CASE_SPECS if args.case is None or x["name"] == args.case]
    solver = load_module("tg_solver_fixedk", SOLVER_PATH)
    builder = load_module("tg_geometry_builder_fixedk", SOURCE_BUILDER_PATH)
    maps = builder.load_maps()
    r, vobs, err, *_ = builder.read_rotmod()
    data = {"r": r, "vobs": vobs, "err": err}
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    rows = []
    progress_path = out_dir / "progress.jsonl"
    progress_path.write_text("")
    for spec in case_specs:
        if time.process_time() - start_cpu >= CPU_CAP_S or time.perf_counter() - start_wall >= WALL_CAP_S:
            break
        row = run_case(spec, solver, builder, maps, data)
        rows.append(row)
        with progress_path.open("a") as stream:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
        if rows and rows[0].get("status") == "ok" and row.get("status") == "ok":
            ref = np.asarray(rows[0]["predicted_speeds_km_s"])
            current = np.asarray(row["predicted_speeds_km_s"])
            delta = current - ref
            row["delta_chi2_vs_reference"] = row["chi2_diagonal_43"] - rows[0]["chi2_diagonal_43"]
            row["delta_speed_rms_vs_reference_km_s"] = float(np.sqrt(np.mean(delta**2)))
            row["delta_speed_max_abs_vs_reference_km_s"] = float(np.max(np.abs(delta)))
            # Refresh progress with comparison fields after the reference exists.
            lines = progress_path.read_text().splitlines()
            lines[-1] = json.dumps(row, sort_keys=True)
            progress_path.write_text("\n".join(lines) + "\n")
        write_outputs(rows, data, start_cpu, start_wall, out_dir)
        print(f"{spec['name']}: Kcrit={row['Kcrit_km_s_minus2']:.12g}, "
              f"status={row['status']}, CPU={row['cpu_seconds']:.1f}s, "
              f"wall={row['wall_seconds']:.1f}s", flush=True)
    write_outputs(rows, data, start_cpu, start_wall, out_dir)


def write_outputs(rows, data, start_cpu, start_wall, out_dir):
    summaries = []
    for r in rows:
        summaries.append({k: v for k, v in r.items() if k != "predicted_speeds_km_s"})
    result = {
        "classification": "fixed physical-K NGC3198 robustness check; not a detection",
        "status": "complete" if len(rows) == len(CASE_SPECS) else "bounded_subset_or_partial",
        "fixed_parameters": {"K_km_s_minus2": K_FIXED, "distance_Mpc": D_MPC,
                             "inclination_deg": INC_DEG, "upsilon": UPSILON,
                             "stellar_h_kpc": HSTAR, "gas": "all thin, h=0.2 kpc"},
        "source_scope": "NGC3198 S4G P5 stellar image and HALOGAS HI+He map, azimuthal means; no NGC4559 inputs",
        "reference_case": CASE_SPECS[0], "cases": summaries,
        "elapsed_cpu_seconds": time.process_time() - start_cpu,
        "elapsed_wall_seconds": time.perf_counter() - start_wall,
        "caps_seconds": {"cpu": CPU_CAP_S, "wall": WALL_CAP_S},
    }
    result_path = out_dir / "results.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n")
    if rows and rows[0].get("predicted_speeds_km_s") is not None:
        reference = np.asarray(rows[0]["predicted_speeds_km_s"])
        csv_path = out_dir / "curves.csv"
        with csv_path.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["R_catalog_kpc", "R_physical_kpc", "Vobs_i72_km_s",
                             "sigmaV_i72_km_s"] + [f"V_{r['case']}_km_s" for r in rows])
            fac = math.sin(math.radians(I_REF_DEG)) / math.sin(math.radians(INC_DEG))
            for j in range(data["r"].size):
                values = [data["r"][j], data["r"][j] * D_MPC / D_REF_MPC,
                          data["vobs"][j] * fac, data["err"][j] * fac]
                for row in rows:
                    predictions = row.get("predicted_speeds_km_s")
                    values.append(predictions[j] if predictions is not None else "")
                writer.writerow(values)
    # Plot only the reconstructed NGC3198 curves and its NGC3198 Rotmod data.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8.5, 5.2))
        fac = math.sin(math.radians(I_REF_DEG)) / math.sin(math.radians(INC_DEG))
        radii = data["r"] * D_MPC / D_REF_MPC
        ax.errorbar(radii, data["vobs"] * fac, yerr=data["err"] * fac,
                    fmt="o", ms=3, color="black", label="NGC 3198 Rotmod, i=72°")
        for row in rows:
            if row.get("predicted_speeds_km_s") is not None:
                ax.plot(radii, row["predicted_speeds_km_s"], lw=1.7,
                        label=f"{row['case']}: K fixed")
        ax.set(xlabel="R [kpc at D=16.6 Mpc]", ylabel="V [km/s]",
               title="NGC 3198 fixed-K FV refinement")
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / "rotation_curve.png", dpi=170)
        plt.close(fig)
    except Exception as exc:
        result["plot_error"] = repr(exc)
        result_path.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Find the first TG positive-u boundary and its limiting circular curve.

For the assembled positive finite-volume Laplacian B and source diagonal D,
the transformed problem is (B-KD)u=b. Its first pole is the smallest
generalized eigenvalue B h=Kcrit D h. Below it, a positive branch exists;
at/above it no strictly positive solution can satisfy the same positive
exterior boundary source b. Although u's normalization diverges as K->Kcrit-,
the logarithmic radial force has the finite limit
V^2=-R/Kcrit d_R log(h_mid).

This script reports that limiting curve separately from finite subcritical
solutions. It reuses the pilot's frozen source/operator builder; independent
discretization and residual review is in the neighboring independent-review
directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh, splu

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOLVER_PATH = HERE / "solve_tg_axisymmetric.py"
SPEC = importlib.util.spec_from_file_location("tg_solver", SOLVER_PATH)
tg = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tg)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def critical_mode(dr, dz, rmax, zmax, upsilon, *, gas_name="halogas", tail_factor=1.0):
    B, D, rho, vol, rr, zz, step, source_cfg = tg.build_operator(
        dr, dz, rmax, zmax, upsilon, gas_name=gas_name, tail_factor=tail_factor)
    # Shift-invert at zero solves the smallest positive generalized eigenvalue
    # without attempting the poorly conditioned M-inverse formulation.
    values, vectors = eigsh(B, k=1, M=diags(D, format="csr"), sigma=0.0,
                            which="LM", tol=2e-8, maxiter=5000)
    kcrit = float(values[0])
    h = vectors[:, 0].reshape(rho.shape)
    if float(h.sum()) < 0:
        h *= -1.0
    h /= float(np.max(h))
    residual = float(np.linalg.norm(B @ vectors[:, 0] - kcrit * D * vectors[:, 0]) /
                      np.linalg.norm(B @ vectors[:, 0]))
    dr_eff, dz_eff = step
    rcent = (np.arange(rho.shape[0]) + 0.5) * dr_eff
    hmid = (9.0 * h[:, 0] - h[:, 1]) / 8.0
    return {
        "B": B, "D": D, "rho": rho, "shape": rho.shape, "step": step,
        "rcent": rcent, "h": h, "hmid": hmid, "kcrit": kcrit,
        "eig_residual": residual, "min_h_full": float(h.min()),
        "max_h_full": float(h.max()), "source_cfg": source_cfg,
    }


def edge_curve(mode, radii):
    hm = mode["hmid"]
    if np.any(hm <= 0):
        return None, {"status": "nonpositive_midplane_eigenmode"}
    rcent = mode["rcent"]
    hfn = PchipInterpolator(rcent, hm, extrapolate=True)
    h_eval = hfn(radii)
    v2 = -radii / mode["kcrit"] * hfn.derivative()(radii) / h_eval
    if np.any(h_eval <= 0) or np.any(v2 <= 0) or np.any(~np.isfinite(v2)):
        return None, {"status": "nonphysical_edge_curve",
                      "minimum_h_at_data": float(np.min(h_eval)),
                      "minimum_v2": float(np.nanmin(v2))}
    return np.sqrt(v2), {"status": "ok"}


def finite_subcritical(mode, K, radii, vobs, errors):
    A = mode["B"] - diags(K * mode["D"], format="csr")
    lu = splu(A.tocsc())
    wvec = lu.solve(mode["D"])
    linear_resid = float(np.linalg.norm(A @ wvec - mode["D"]) /
                         np.linalg.norm(mode["D"]))
    w = wvec.reshape(mode["shape"])
    _, dz = mode["step"]
    wmid = (9.0 * w[:, 0] - w[:, 1]) / 8.0
    u = 1.0 + K * w
    u_mid = (9.0 * u[:, 0] - u[:, 1]) / 8.0
    if np.any(u <= 0) or np.any(u_mid <= 0):
        return {"status": "u_nonpositive", "minimum_u": float(u.min()),
                "linear_residual": linear_resid}
    wfn = PchipInterpolator(mode["rcent"], wmid, extrapolate=True)
    v2 = -radii * wfn.derivative()(radii) / (1.0 + K * wfn(radii))
    if np.any(v2 <= 0) or np.any(~np.isfinite(v2)):
        return {"status": "nonphysical_speed", "linear_residual": linear_resid,
                "minimum_u": float(u.min()), "minimum_v2": float(np.nanmin(v2))}
    speed = np.sqrt(v2)
    return {
        "status": "ok", "K": K, "fractional_distance_to_critical": 1-K/mode["kcrit"],
        "minimum_u_full": float(u.min()), "maximum_u_full": float(u.max()),
        "linear_residual": linear_resid,
        "chi2_diagonal_43": float(np.sum(((vobs-speed)/errors)**2)),
        "rms_resid_km_s": float(np.sqrt(np.mean((vobs-speed)**2))),
        "speeds_km_s": speed.tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu-cap", type=float, default=240.0)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    start_cpu = time.process_time()
    start_wall = time.perf_counter()
    radii, vobs, errors, _, _ = tg.read_rotmod()
    k0 = np.genfromtxt(HERE / "newtonian_k0_predictions.csv", delimiter=",", names=True)
    # Include the S4G +/-0.1 dex envelope, its central value, and the paper's
    # quoted NGC3198 M/L as a reference rather than as a prior.
    upsilons = [0.6*10**-0.1, 0.6, 0.6*10**0.1, 0.762]
    summaries, curves, best_mode = [], {}, None
    for upsilon in upsilons:
        if time.process_time() - start_cpu > args.cpu_cap:
            break
        mode = critical_mode(.25, .125, 80., 60., upsilon)
        edge, edge_status = edge_curve(mode, radii)
        record = {
            "upsilon": upsilon, "Kcrit_km_s_minus2": mode["kcrit"],
            "eigen_residual": mode["eig_residual"],
            "minimum_principal_mode_full": mode["min_h_full"],
            "maximum_principal_mode_full": mode["max_h_full"],
            "edge_curve_status": edge_status["status"],
            "grid": {"dr_kpc": .25, "dz_kpc": .125, "Rmax_kpc": 80.,
                     "Zmax_kpc": 60., "nr": mode["shape"][0], "nz": mode["shape"][1]},
            "source_cfg": mode["source_cfg"],
        }
        if edge is not None:
            record.update({
                "chi2_diagonal_43_limiting_curve": float(np.sum(((vobs-edge)/errors)**2)),
                "rms_resid_km_s_limiting_curve": float(np.sqrt(np.mean((vobs-edge)**2))),
                "min_limiting_speed_km_s": float(edge.min()),
                "max_limiting_speed_km_s": float(edge.max()),
            })
            curves[f"edge_Upsilon_{upsilon:.6f}"] = edge.tolist()
        summaries.append(record)
        if upsilon == .762:
            best_mode = mode
            np.savez_compressed(args.output_dir / "axisymmetric_principal_mode_U0p762.npz",
                                h=mode["h"], r_kpc=mode["rcent"],
                                kcrit=mode["kcrit"], upsilon=upsilon)

    # Limited but physically motivated source sensitivities at the high-M/L
    # reference: independent gas-map conversion and +/-25% outer stellar tail
    # scale. These are alternatives, not a profile fit to the observed curve.
    source_sensitivities = []
    baseline_edge = curves.get("edge_Upsilon_0.762000")
    variant_specs = [("THINGS_gas_map", "things", 1.0),
                     ("shorter_stellar_tail_Rd_x0.75", "halogas", .75),
                     ("longer_stellar_tail_Rd_x1.25", "halogas", 1.25)]
    for label, gas_name, tail_factor in variant_specs:
        if time.process_time() - start_cpu > args.cpu_cap:
            break
        variant = critical_mode(.25, .125, 80., 60., .762,
                                gas_name=gas_name, tail_factor=tail_factor)
        variant_speed, variant_status = edge_curve(variant, radii)
        item = {
            "label": label, "gas_name": gas_name, "stellar_tail_factor": tail_factor,
            "upsilon": .762, "Kcrit_km_s_minus2": variant["kcrit"],
            "eigen_residual": variant["eig_residual"],
            "principal_mode_min": variant["min_h_full"],
            "edge_curve_status": variant_status["status"],
            "source_cfg": variant["source_cfg"],
        }
        if variant_speed is not None:
            item["chi2_diagonal_43_limiting_curve"] = float(np.sum(((vobs-variant_speed)/errors)**2))
            item["rms_resid_km_s_limiting_curve"] = float(np.sqrt(np.mean((vobs-variant_speed)**2)))
            item["max_speed_difference_vs_baseline_km_s"] = float(
                np.max(np.abs(variant_speed-baseline_edge)))
            item["rms_speed_difference_vs_baseline_km_s"] = float(
                np.sqrt(np.mean((variant_speed-baseline_edge)**2)))
            curves["source_"+label] = variant_speed.tolist()
        source_sensitivities.append(item)

    finite = []
    if best_mode is not None:
        for fraction in (0.95, 0.99, 0.998):
            K = best_mode["kcrit"] * fraction
            item = finite_subcritical(best_mode, K, radii, vobs, errors)
            item["upsilon"] = .762
            finite.append(item)
            if item.get("status") == "ok":
                curves[f"finite_{fraction:.3f}_Kcrit_Upsilon_0.762"] = item["speeds_km_s"]

    # Separate resolution from box-boundary sensitivity: same coarse spacing
    # on the base and enlarged boxes, then compare both to the production grid.
    coarse_mode = critical_mode(.5, .25, 80., 60., .762)
    coarse_edge, coarse_status = edge_curve(coarse_mode, radii)
    domain_mode = critical_mode(.5, .25, 100., 80., .762)
    domain_edge, domain_status = edge_curve(domain_mode, radii)
    domain_record = {
        "upsilon": .762, "Kcrit_km_s_minus2": domain_mode["kcrit"],
        "eigen_residual": domain_mode["eig_residual"],
        "edge_curve_status": domain_status["status"],
        "grid": {"dr_kpc": .5, "dz_kpc": .25, "Rmax_kpc": 100.,
                 "Zmax_kpc": 80., "nr": domain_mode["shape"][0],
                 "nz": domain_mode["shape"][1]},
    }
    coarse_record = {
        "upsilon": .762, "Kcrit_km_s_minus2": coarse_mode["kcrit"],
        "eigen_residual": coarse_mode["eig_residual"],
        "edge_curve_status": coarse_status["status"],
        "grid": {"dr_kpc": .5, "dz_kpc": .25, "Rmax_kpc": 80.,
                 "Zmax_kpc": 60., "nr": coarse_mode["shape"][0],
                 "nz": coarse_mode["shape"][1]},
    }
    fine_u762 = next(x for x in summaries if abs(x["upsilon"]-.762) < 1e-12)
    fine_edge = curves.get("edge_Upsilon_0.762000")
    for label, candidate, candidate_edge in (("coarse_same_domain", coarse_record, coarse_edge),
                                               ("coarse_larger_domain", domain_record, domain_edge)):
        candidate["Kcrit_fractional_difference_vs_fine"] = (
            candidate["Kcrit_km_s_minus2"] / fine_u762["Kcrit_km_s_minus2"] - 1.0)
        if candidate_edge is not None and fine_edge is not None:
            delta = np.asarray(candidate_edge) - np.asarray(fine_edge)
            candidate["edge_speed_max_abs_difference_vs_fine_km_s"] = float(np.max(np.abs(delta)))
            candidate["edge_speed_RMS_difference_vs_fine_km_s"] = float(np.sqrt(np.mean(delta**2)))
            candidate["edge_chi2_diagonal_43"] = float(np.sum(((vobs-candidate_edge)/errors)**2))
            curves[label + "_edge_Upsilon_0.762"] = candidate_edge.tolist()
    summaries.append(coarse_record)
    summaries.append(domain_record)

    output_rows = []
    for idx, radius in enumerate(radii):
        row = {"R_kpc": radius, "Vobs_km_s": vobs[idx], "sigmaV_km_s": errors[idx],
               "K0_resolved_baryons_U0p6_km_s": k0["Vbar_map_HALOGAS_K0_km_s"][idx],
               "K0_resolved_baryons_U0p762_km_s": float(np.sqrt(max(
                   k0["Vstar_map_K0_km_s"][idx]**2 * (.762/.6) +
                   (k0["Vbar_map_HALOGAS_K0_km_s"][idx]**2 -
                    k0["Vstar_map_K0_km_s"][idx]**2), 0.0))) }
        for key, speeds in curves.items():
            row[key + "_km_s"] = speeds[idx]
        output_rows.append(row)
    csv_path = args.output_dir / "axisymmetric_branch_edge_curves.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    fig, ax = plt.subplots(figsize=(9, 5.6))
    ax.errorbar(radii, vobs, yerr=errors, fmt="o", ms=3, color="black", label="SPARC NGC 3198")
    for key, speeds in curves.items():
        if key.startswith("edge_"):
            ax.plot(radii, speeds, lw=2, label=key.replace("edge_", "K→Kcrit: "))
    ax.plot(radii, k0["Vbar_map_HALOGAS_K0_km_s"], "--", color="0.45",
            label="resolved baryons, K=0 (Υ*=0.6)")
    ax.set(xlabel="R [kpc]", ylabel="circular speed [km/s]",
           title="NGC 3198: TG first positive-u boundary and limiting curve")
    ax.grid(alpha=.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    plot_path = args.output_dir / "axisymmetric_branch_edge.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    summary = {
        "classification": "principal-mode limiting curve is not a finite-normalization positive-u solution",
        "source_csv_sha256": file_hash(HERE / "source_profiles.csv"),
        "solver_sha256": file_hash(SOLVER_PATH),
        "analysis_script_sha256": file_hash(Path(__file__)),
        "Rotmod_sha256": file_hash(ROOT / "context/data/Rotmod_LTG.zip"),
        "K0_csv_sha256": file_hash(HERE / "newtonian_k0_predictions.csv"),
        "eigen_method": "scipy eigsh(B,k=1,M=diag(D),sigma=0,which=LM); D_i=4piG rho_i R_i dR dz",
        "field_equation": "(B-KD)u=b; first positive-solution edge is smallest generalized eigenvalue B h=Kcrit D h",
        "velocity_limit": "v_c^2=-R/Kcrit * d_R log(h_midplane); compare separately from finite-K solutions",
        "chi2_scope": "43 pointwise SPARC errV values, diagonal descriptive score only, no full covariance",
        "fine_grid_results": [x for x in summaries if x["grid"]["dr_kpc"] == .25],
        "resolution_and_domain_sensitivity": {
            "same_domain_coarse": coarse_record,
            "coarse_extended_domain": domain_record,
        },
        "source_sensitivities_at_U0p762": source_sensitivities,
        "finite_subcritical_at_U0p762": finite,
        "cpu_seconds": time.process_time()-start_cpu, "wall_seconds": time.perf_counter()-start_wall,
        "cpu_cap_seconds": args.cpu_cap,
    }
    (args.output_dir / "axisymmetric_branch_edge_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

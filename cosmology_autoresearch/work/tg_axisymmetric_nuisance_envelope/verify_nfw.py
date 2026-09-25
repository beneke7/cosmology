#!/usr/bin/env python3
"""Independent all-node optimizer audit for the saved NFW nuisance grid.

Reads the existing envelope progress rows but does not alter them. Rebuilds the
same map-derived source curves, reruns each of the three deterministic NFW
starts at all 450 geometry/vertical/M-L nodes, and records every fit status and
objective. A single process with BLAS threads fixed to one is used; a 120 CPU
second cap is enforced internally.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = HERE/"nuisance_envelope.py"
PROGRESS = HERE/"envelope_progress.jsonl"
OUT = HERE/"nfw_verification.jsonl"
SUMMARY = HERE/"nfw_verification_summary.json"


def file_hash(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module():
    spec = importlib.util.spec_from_file_location("tg_nuisance_envelope", SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fit_all_starts(envelope, r, obs, err, baryon_v2):
    lo = np.array([9.0, 0.0])
    hi = np.array([14.0, math.log10(50.0)])

    def predict(theta):
        return np.sqrt(np.maximum(baryon_v2+envelope.nfw_v2(
            r, float(theta[0]), float(theta[1])), 0.0))

    def residual(theta):
        return (predict(theta)-obs)/err

    starts = [(11.0, 0.9), (12.0, 0.7), (11.5, 1.1)]
    results = []
    for initial in starts:
        fit = least_squares(residual, np.asarray(initial), bounds=(lo, hi),
                            x_scale=np.array([1.0, 0.5]), max_nfev=1500,
                            ftol=1e-10, xtol=1e-10, gtol=1e-10)
        pred = predict(fit.x)
        chi2 = float(np.dot(fit.fun, fit.fun))
        results.append({"initial_log10_M200_c200": list(initial),
                        "success": bool(fit.success), "status": int(fit.status),
                        "message": str(fit.message), "nfev": int(fit.nfev),
                        "njev": int(fit.njev) if fit.njev is not None else None,
                        "optimality": float(fit.optimality),
                        "active_mask": fit.active_mask.astype(int).tolist(),
                        "objective_chi2_data": chi2,
                        "log10_M200_Msun": float(fit.x[0]),
                        "M200_Msun": float(10.0**fit.x[0]),
                        "log10_c200": float(fit.x[1]),
                        "c200": float(10.0**fit.x[1]),
                        "prediction_min_km_s": float(np.min(pred)),
                        "prediction_max_km_s": float(np.max(pred))})
    selected_idx = min(range(len(results)), key=lambda idx: results[idx]["objective_chi2_data"])
    return results, selected_idx


def main(cpu_cap: float = 120.0):
    started_cpu, started_wall = time.process_time(), time.perf_counter()
    env = load_module()
    progress_by_key = {}
    with PROGRESS.open() as stream:
        for line in stream:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("model") == "NFW":
                key = (row["vertical_case"], round(row["distance_Mpc"], 8),
                       round(row["inclination_deg"], 8), round(row["upsilon"], 10))
                progress_by_key[key] = row
    expected_keys = len(env.D_SIGMA_GRID)*len(env.I_SIGMA_GRID)*len(env.VERTICAL_CASES)*len(env.Q_GRID)
    if len(progress_by_key) != expected_keys:
        raise RuntimeError(f"expected {expected_keys} stored NFW rows; found {len(progress_by_key)}")

    maps = env.load_maps()
    r_catalog, v_catalog, e_catalog, _, _ = env.read_rotmod()
    source_r = np.arange(0.0, env.R_SOURCE_MAX+0.5*env.R_SOURCE_STEP, env.R_SOURCE_STEP)
    d_grid = [env.D_REF+x*env.DISTANCE_SIGMA_MPC for x in env.D_SIGMA_GRID]
    i_grid = [env.I_REF+x*env.INCLINATION_SIGMA_DEG for x in env.I_SIGMA_GRID]
    u_grid = [env.UPSILON_REF*10.0**q for q in env.Q_GRID]
    out_rows = []
    best_global = None
    max_chi2_delta = 0.0
    max_objective_delta = 0.0
    max_logm_delta = 0.0
    max_logc_delta = 0.0
    selected_success_count = 0
    start_success_count = 0
    complete = True
    OUT.write_text("")
    for d in d_grid:
        for inc in i_grid:
            if time.process_time()-started_cpu > cpu_cap:
                complete = False
                break
            prof = env.transformed_sources(maps, d, inc)
            reval = r_catalog*(d/env.D_REF)
            obs, err, factor = env.observed_for_inclination(v_catalog, e_catalog, inc)
            star_sigma, star_cfg = env.surface_density_on_grid(source_r, prof["r"], prof["stars"], "stars")
            gas_sigma, gas_cfg = env.surface_density_on_grid(source_r, prof["r"], prof["gas"], "gas")
            ks, ts = env.hankel_transform(source_r, star_sigma)
            kg, tgtr = env.hankel_transform(source_r, gas_sigma)
            star_v2 = {h: env.hankel_v2_from_transform(reval, ks, ts, h) for h in (0.3, 0.6)}
            gas_thin = env.hankel_v2_from_transform(reval, kg, tgtr, 0.2)
            gas_thick = env.hankel_v2_from_transform(reval, kg, tgtr, 3.0)
            for vertical, vc in env.VERTICAL_CASES.items():
                gas_v2 = vc["f_thin"]*gas_thin+(1.0-vc["f_thin"])*gas_thick
                for q, upsilon in zip(env.Q_GRID, u_grid):
                    if time.process_time()-started_cpu > cpu_cap:
                        complete = False
                        break
                    key = (vertical, round(d, 8), round(inc, 8), round(upsilon, 10))
                    old = progress_by_key[key]
                    baryon_v2 = (upsilon/env.UPSILON_REF)*star_v2[vc["hstar"]]+gas_v2
                    starts, selected_idx = fit_all_starts(env, reval, obs, err, baryon_v2)
                    selected = starts[selected_idx]
                    prior_d, prior_i = env.geometry_prior(d, inc)
                    prior_ml = (q/env.MSTAR_SIGMA_DEX)**2
                    prior_total = prior_d+prior_i+prior_ml
                    objective = selected["objective_chi2_data"]+prior_total
                    d_chi2 = selected["objective_chi2_data"]-old["chi2_data"]
                    d_obj = objective-old["objective_chi2_plus_priors"]
                    d_logm = selected["log10_M200_Msun"]-old["NFW_log10_M200"]
                    d_logc = selected["log10_c200"]-old["NFW_log10_c200"]
                    max_chi2_delta = max(max_chi2_delta, abs(d_chi2))
                    max_objective_delta = max(max_objective_delta, abs(d_obj))
                    max_logm_delta = max(max_logm_delta, abs(d_logm))
                    max_logc_delta = max(max_logc_delta, abs(d_logc))
                    selected_success_count += int(selected["success"])
                    start_success_count += sum(int(x["success"]) for x in starts)
                    row = {"vertical_case": vertical, "distance_Mpc": d,
                           "inclination_deg": inc, "D_sigma": (d-env.D_REF)/env.DISTANCE_SIGMA_MPC,
                           "i_sigma": (inc-env.I_REF)/env.INCLINATION_SIGMA_DEG,
                           "upsilon": upsilon, "q_log10_ML": q,
                           "prior_D": prior_d, "prior_i": prior_i, "prior_ML": prior_ml,
                           "stored_chi2_data": old["chi2_data"],
                           "verified_chi2_data": selected["objective_chi2_data"],
                           "delta_chi2_data": d_chi2,
                           "stored_objective_chi2_plus_priors": old["objective_chi2_plus_priors"],
                           "verified_objective_chi2_plus_priors": objective,
                           "delta_objective": d_obj,
                           "stored_log10_M200": old["NFW_log10_M200"],
                           "verified_log10_M200": selected["log10_M200_Msun"],
                           "delta_log10_M200": d_logm,
                           "stored_log10_c200": old["NFW_log10_c200"],
                           "verified_log10_c200": selected["log10_c200"],
                           "delta_log10_c200": d_logc,
                           "selected_start_index": selected_idx,
                           "selected_fit_success": selected["success"],
                           "selected_fit_message": selected["message"],
                           "selected_fit_nfev": selected["nfev"],
                           "all_starts": starts,
                           "source_fit": {"stars_tail": star_cfg, "gas_tail": gas_cfg},
                           "obs_velocity_factor_from_i73": factor}
                    out_rows.append(row)
                    with OUT.open("a") as stream:
                        stream.write(json.dumps(row, sort_keys=True)+"\n")
                    if (best_global is None or objective < best_global["verified_objective_chi2_plus_priors"]):
                        best_global = {k: row[k] for k in row if k != "all_starts"}
                if not complete:
                    break
            if not complete:
                break
        if not complete:
            break

    summary = {"status": "complete" if complete and len(out_rows)==expected_keys else "partial_cpu_cap",
               "node_count_expected": expected_keys, "node_count_verified": len(out_rows),
               "starts_per_node": 3, "total_start_fits": len(out_rows)*3,
               "selected_fit_success_count": selected_success_count,
               "total_successful_start_fits": start_success_count,
               "all_selected_fits_successful": selected_success_count==len(out_rows),
               "all_450_node_scores_match_existing": max_chi2_delta < 1e-8 and len(out_rows)==expected_keys,
               "max_abs_delta_chi2_data": max_chi2_delta,
               "max_abs_delta_total_objective": max_objective_delta,
               "max_abs_delta_log10_M200": max_logm_delta,
               "max_abs_delta_log10_c200": max_logc_delta,
               "selected_start_selection": "minimum objective across the same three deterministic starts; bounds [9,14] log10(M200/Msun), [0,log10(50)] log10(c200)",
               "optimizer": {"method": "scipy.optimize.least_squares", "x_scale": [1.0,0.5],
                             "ftol": 1e-10, "xtol": 1e-10, "gtol": 1e-10, "max_nfev": 1500,
                             "processes": 1, "BLAS_threads": 1},
               "best_verified_global_NFW_node": best_global,
               "runtime": {"cpu_seconds": time.process_time()-started_cpu,
                           "wall_seconds": time.perf_counter()-started_wall,
                           "cpu_cap_seconds": cpu_cap},
               "inputs_sha256": {str(path.relative_to(ROOT)): file_hash(path) for path in
                                 [SOURCE, PROGRESS, ROOT/"work/tg_axisymmetric_pilot/newtonian_k0_predictions.csv",
                                  ROOT/"context/data/Rotmod_LTG.zip"]},
               "unchanged_original_results": "results.json and envelope_progress.jsonl are read-only inputs; this verifier writes separate nfw_verification* artifacts"}
    SUMMARY.write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu-cap", type=float, default=120.0)
    args = parser.parse_args()
    main(args.cpu_cap)

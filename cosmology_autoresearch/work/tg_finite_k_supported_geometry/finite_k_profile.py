#!/usr/bin/env python3
"""Bounded finite-u TG profile on H I-supported NGC 3198 geometry nodes.

All scores are exploratory, diagonal, in-sample summaries. K/Kcrit never
includes the singular Kcrit point; K is an added profiled parameter without a
prior. The output directory is task-owned and separate from the old envelope.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import platform
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator
from scipy.sparse import diags
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
ENV_PATH = ROOT / "work/tg_axisymmetric_nuisance_envelope/nuisance_envelope.py"
SOLVER_PATH = ROOT / "work/tg_axisymmetric_pilot/solve_tg_axisymmetric.py"
ENVELOPE_CSV = ROOT / "work/tg_axisymmetric_nuisance_envelope/envelope_grid.csv"
MODE_JSON = ROOT / "work/tg_axisymmetric_independent_review/principal_mode_result.json"
FRACTIONS = (0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.998)
DISTANCES = (13.8, 16.6)  # SPARC prior center and +2 sigma (13.8 +/- 1.4 Mpc)
INCLINATION = 72.0  # THINGS source geometry; HALOGAS adopts about 70 degrees.
UPSILONS = (0.4766, 0.6, 0.762)  # existing grid: -1 sigma, center, published reference
VERTICALS = {
    "gas_mix_hstar0p3": {"hstar": 0.3, "f_thin": 0.85},
    "gas_thin_hstar0p3": {"hstar": 0.3, "f_thin": 1.0},
    "gas_mix_hstar0p6": {"hstar": 0.6, "f_thin": 0.85},
}
MAX_CPU_DEFAULT = 175.0
MAX_WALL_DEFAULT = 410.0


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def write_json(path: Path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def source_case(env, tg, maps, B, volume, rcent, zcent, step, rotmod,
                distance, inc, upsilon, vertical, cpu_started, cpu_cap, wall_started, wall_cap):
    """Compute one source-specific Kcrit, eight finite-u curves, and matched controls."""
    prof = env.transformed_sources(maps, distance, inc)
    source_r = np.arange(0.0, env.R_SOURCE_MAX + 0.5 * env.R_SOURCE_STEP, env.R_SOURCE_STEP)
    star_sigma, star_cfg = env.surface_density_on_grid(source_r, prof["r"], prof["stars"], "stars")
    gas_sigma, gas_cfg = env.surface_density_on_grid(source_r, prof["r"], prof["gas"], "gas")
    star_on_grid = np.interp(rcent, source_r, star_sigma)
    gas_on_grid = np.interp(rcent, source_r, gas_sigma)
    vc = VERTICALS[vertical]
    star_rho, gas_rho = env.tg_vertical_density(star_on_grid, gas_on_grid, rcent, zcent,
                                                 vc["hstar"], vc["f_thin"])
    rho = (upsilon / env.UPSILON_REF) * star_rho + gas_rho
    dvec = (4.0 * math.pi * env.G * rho * volume).ravel()

    # Recompute the source-specific principal mode using the checked FV operator.
    edge = env.solve_tg_edge(B, volume, rcent, zcent, star_rho, gas_rho, upsilon,
                             rotmod[:, 0] * distance / env.D_REF,
                             *env.observed_for_inclination(rotmod[:, 1], rotmod[:, 2], inc)[:2])
    kcrit = float(edge["kcrit"])
    reval = rotmod[:, 0] * distance / env.D_REF
    obs, err, obs_factor = env.observed_for_inclination(rotmod[:, 1], rotmod[:, 2], inc)
    prior_d, prior_i = env.geometry_prior(distance, inc)
    q = math.log10(upsilon / env.UPSILON_REF)
    prior_ml = (q / env.MSTAR_SIGMA_DEX) ** 2
    prior_sum = prior_d + prior_i + prior_ml

    # Same-source, same-nuisance controls at this exact D/i/M/L/vertical node.
    ks, ts = env.hankel_transform(source_r, star_sigma)
    kg, gt = env.hankel_transform(source_r, gas_sigma)
    star_v2 = env.hankel_v2_from_transform(reval, ks, ts, vc["hstar"])
    gas_v2 = (vc["f_thin"] * env.hankel_v2_from_transform(reval, kg, gt, 0.2)
              + (1.0-vc["f_thin"]) * env.hankel_v2_from_transform(reval, kg, gt, 3.0))
    baryon_v2 = (upsilon / env.UPSILON_REF) * star_v2 + gas_v2
    nfw = env.fit_nfw_for_geometry(reval, obs, err, baryon_v2, upsilon)
    gnewton = baryon_v2 / reval
    gmond = 0.5 * (gnewton + np.sqrt(gnewton*gnewton + 4.0*env.A0*gnewton))
    mond_speed = np.sqrt(reval * gmond)
    mond_chi2 = float(np.sum(((obs-mond_speed)/err)**2))
    controls = {
        "NFW": {"chi2_diag": nfw["chi2_data"], "objective_chi2_plus_priors": nfw["chi2_data"]+prior_sum,
                "optimizer_success": nfw["optimizer_success"], "optimizer_message": nfw["optimizer_message"],
                "log10_M200_Msun": nfw["log10_M200_Msun"], "c200": nfw["concentration_c200"]},
        "MOND": {"chi2_diag": mond_chi2, "objective_chi2_plus_priors": mond_chi2+prior_sum,
                 "a0_m_s2": 1.2e-10, "mu": "x/(1+x)", "EFE": False},
    }
    common = {
        "D_Mpc": distance, "D_sigma": (distance-env.D_REF)/env.DISTANCE_SIGMA_MPC,
        "i_deg": inc, "i_sigma": (inc-env.I_REF)/env.INCLINATION_SIGMA_DEG,
        "upsilon": upsilon, "q_log10_ML_dex": q, "vertical_case": vertical,
        "Kcrit_km_s_minus2": kcrit, "principal_mode_relative_residual": edge["mode_residual"],
        "edge_limit_chi2_diag_only_for_context_not_finite": edge["chi2_data"],
        "edge_limit_objective_only_for_context_not_finite": edge["chi2_data"]+prior_sum,
        "prior_D": prior_d, "prior_i": prior_i, "prior_ML": prior_ml,
        "prior_K": None, "prior_total_existing_terms": prior_sum,
        "obs_velocity_factor_from_i73": obs_factor,
        "star_tail": star_cfg, "gas_tail": gas_cfg, "matched_controls_same_source_node": controls,
    }
    case_id = f"D{distance:.1f}_i{inc:.0f}_U{upsilon:.4f}_{vertical}"
    records, curves = [], []
    for fraction in FRACTIONS:
        if time.process_time()-cpu_started >= cpu_cap-1.0 or time.perf_counter()-wall_started >= wall_cap-3.0:
            records.append({**common, "case_id": case_id, "K_over_Kcrit": fraction,
                            "K_km_s_minus2": fraction*kcrit, "status": "not_run_cpu_or_wall_guard"})
            continue
        kval = fraction * kcrit
        record = {**common, "case_id": case_id, "K_over_Kcrit": fraction,
                  "K_km_s_minus2": kval, "status": "factorization_failed"}
        try:
            A = B - diags(kval*dvec, format="csr")
            lu = splu(A.tocsc())
            wvec = lu.solve(dvec)
            residual = float(np.linalg.norm(A @ wvec-dvec) / np.linalg.norm(dvec))
            w = wvec.reshape(rho.shape)
            u = 1.0 + kval*w
            wmid = (9.0*w[:, 0]-w[:, 1])/8.0
            umid = (9.0*u[:, 0]-u[:, 1])/8.0
            wf = PchipInterpolator(rcent, wmid, extrapolate=True)
            u_eval = 1.0 + kval*wf(reval)
            v2 = -reval*wf.derivative()(reval)/u_eval
            record.update({
                "relative_linear_solve_residual": residual,
                "min_u_full_half_domain": float(np.min(u)), "max_u_full_half_domain": float(np.max(u)),
                "min_u_midplane_extrapolated": float(np.min(umid)),
                "max_u_midplane_extrapolated": float(np.max(umid)),
                "min_u_at_rotation_radii": float(np.min(u_eval)), "max_u_at_rotation_radii": float(np.max(u_eval)),
                "full_domain_u_positive": bool(np.all(np.isfinite(u)) and np.min(u)>0.0),
                "midplane_u_positive": bool(np.all(np.isfinite(umid)) and np.min(umid)>0.0),
                "speed_squared_positive_all_points": bool(np.all(np.isfinite(v2)) and np.all(v2>0.0)),
                "min_v2_km2_s2": float(np.nanmin(v2)),
            })
            if residual > 1.0e-8:
                record["status"] = "linear_residual_failure"
            elif not record["full_domain_u_positive"] or not record["midplane_u_positive"]:
                record["status"] = "u_nonpositive"
            elif not record["speed_squared_positive_all_points"]:
                record["status"] = "nonphysical_v2"
            else:
                speed = np.sqrt(v2)
                resid = obs-speed
                chi2 = float(np.sum((resid/err)**2))
                record.update({"status": "ok", "chi2_diag": chi2,
                               "objective_chi2_plus_existing_priors": chi2+prior_sum,
                               "rms_residual_km_s": float(np.sqrt(np.mean(resid**2))),
                               "residuals_km_s": resid.tolist(), "speeds_km_s": speed.tolist()})
                for ix in range(len(reval)):
                    curves.append({"case_id": case_id, "K_over_Kcrit": fraction,
                                   "index": ix, "R_kpc_at_D": float(reval[ix]),
                                   "Vobs_corrected_km_s": float(obs[ix]), "sigmaV_corrected_km_s": float(err[ix]),
                                   "V_TG_finite_km_s": float(speed[ix]), "residual_km_s": float(resid[ix])})
        except Exception as exc:
            record.update({"status": "factorization_or_evaluation_failure", "error": repr(exc)})
        records.append(record)
        with (OUT/"progress.jsonl").open("a") as stream:
            stream.write(json.dumps(record, allow_nan=False, sort_keys=True)+"\n")
            stream.flush()
    return common, records, curves


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu-cap", type=float, default=MAX_CPU_DEFAULT)
    ap.add_argument("--wall-cap", type=float, default=MAX_WALL_DEFAULT)
    args = ap.parse_args()
    started_wall, started_cpu = time.perf_counter(), time.process_time()
    env = load_module("tg_env_finite", ENV_PATH)
    tg = load_module("tg_solver_finite", SOLVER_PATH)
    maps = env.load_maps()
    rotmod = np.asarray(env.read_rotmod()).T
    # The exact Laplacian is source-independent. Its grid matches the checked envelope.
    B, _, _, volume, rr, zz, step, _ = tg.build_operator(
        env.TG_DR, env.TG_DZ, env.TG_RMAX, env.TG_ZMAX, env.UPSILON_REF)
    rcent = (np.arange(rr.shape[0])+0.5)*step[0]
    zcent = (np.arange(rr.shape[1])+0.5)*step[1]
    (OUT/"progress.jsonl").write_text("")
    cases = [(d, INCLINATION, u, v) for d in DISTANCES for v in VERTICALS for u in UPSILONS]
    all_records, all_curves, source_cases = [], [], []
    hard_stop = False
    for d, inc, upsilon, vertical in cases:
        if time.process_time()-started_cpu >= args.cpu_cap-2.0 or time.perf_counter()-started_wall >= args.wall_cap-5.0:
            hard_stop = True
            break
        common, records, curves = source_case(env, tg, maps, B, volume, rcent, zcent, step,
                                              rotmod, d, inc, upsilon, vertical,
                                              started_cpu, args.cpu_cap, started_wall, args.wall_cap)
        source_cases.append(common)
        all_records.extend(records)
        all_curves.extend(curves)

    good = [r for r in all_records if r.get("status") == "ok"]
    best = min(good, key=lambda r:r["objective_chi2_plus_existing_priors"]) if good else None
    by_vertical = {}
    for vertical in VERTICALS:
        group = [r for r in good if r["vertical_case"] == vertical]
        by_vertical[vertical] = min(group, key=lambda r:r["objective_chi2_plus_existing_priors"]) if group else None
    counts = {}
    for r in all_records:
        counts[r.get("status", "missing")] = counts.get(r.get("status", "missing"), 0)+1
    results = {
        "status": "complete" if len(all_records)==len(cases)*len(FRACTIONS) and counts.get("not_run_cpu_or_wall_guard",0)==0 else "partial_or_failures",
        "classification": "exploratory_in_sample_diagonal_screen_no_detection_or_model_selection_significance",
        "predeclared_grid": {"D_Mpc": list(DISTANCES), "D_prior": "13.8 +/- 1.4 Mpc; selected 0 and +2 sigma nodes",
                             "inclination_deg": [INCLINATION], "source_support": "THINGS adopts 72 deg; HALOGAS describes about 70 deg",
                             "upsilon": list(UPSILONS), "vertical_cases": VERTICALS,
                             "K_over_Kcrit": list(FRACTIONS), "singular_Kcrit_sampled_as_finite": False},
        "scope": "NGC 3198 source profile only. NGC 4559 holdout not read or changed. K is an added in-sample profiled galaxy parameter and has no prior.",
        "priors": {"P_D": "((D-13.8)/1.4)^2", "P_i": "((i-73)/3)^2",
                   "P_ML": "(log10(Upsilon/0.6)/0.1)^2", "P_K": None},
        "likelihood": "43 pointwise SPARC Vobs/errorV diagonal chi-square; no radial covariance; descriptive only",
        "matched_controls": "NFW and fixed-a0 simple MOND computed directly at every exact TG D/i/M/L/vertical node from the same transformed source and deprojected observations. The saved envelope has no exact i=72 rows.",
        "n_source_cases_expected": len(cases), "n_source_cases_completed": len(source_cases),
        "n_finite_rows_expected": len(cases)*len(FRACTIONS), "n_finite_rows": len(all_records),
        "status_counts": counts, "best_finite": best, "best_finite_by_vertical_case": by_vertical,
        "source_case_summaries": source_cases, "rows": all_records,
        "runtime": {"cpu_seconds": time.process_time()-started_cpu,
                    "wall_seconds": time.perf_counter()-started_wall,
                    "cpu_cap_seconds": args.cpu_cap, "wall_cap_seconds": args.wall_cap,
                    "threads": {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS")},
                    "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_inputs_sha256": {},
        "references": {"existing_envelope_csv_sha256": sha(ENVELOPE_CSV),
                       "checked_principal_mode_json_sha256": sha(MODE_JSON),
                       "source_builder_sha256": sha(ENV_PATH), "solver_sha256": sha(SOLVER_PATH)},
    }
    input_paths = {
        "S4G_stellar_map": ROOT/"work/tg_axisymmetric_pilot/data/NGC3198.stellar.fits",
        "S4G_mask": ROOT/"work/tg_axisymmetric_pilot/data/NGC3198.ICAmask.fits",
        "HALOGAS_HI_map": ROOT/"work/tg_resolved_source_audit/NGC3198-HR_coldens.fits",
        "THINGS_moment0_map": ROOT/"work/tg_resolved_source_audit/NGC_3198_NA_MOM0_THINGS.FITS",
        "SPARC_Rotmod_zip": ROOT/"context/data/Rotmod_LTG.zip",
        "old_envelope_results_json": ROOT/"work/tg_axisymmetric_nuisance_envelope/results.json",
    }
    results["source_inputs_sha256"] = {k:sha(v) for k,v in input_paths.items()}
    results["runner_sha256"] = sha(Path(__file__).resolve())

    write_json(OUT/"results.json", results)
    scalar_fields = ["case_id","status","D_Mpc","D_sigma","i_deg","i_sigma","upsilon","q_log10_ML_dex","vertical_case",
                     "Kcrit_km_s_minus2","K_over_Kcrit","K_km_s_minus2","chi2_diag","objective_chi2_plus_existing_priors",
                     "prior_D","prior_i","prior_ML","prior_K","prior_total_existing_terms",
                     "relative_linear_solve_residual","min_u_full_half_domain","max_u_full_half_domain",
                     "min_u_midplane_extrapolated","max_u_midplane_extrapolated","full_domain_u_positive",
                     "midplane_u_positive","speed_squared_positive_all_points","rms_residual_km_s"]
    with (OUT/"finite_profile.csv").open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=scalar_fields, extrasaction="ignore"); w.writeheader(); w.writerows(all_records)
    curve_fields = list(all_curves[0]) if all_curves else ["case_id","K_over_Kcrit","index","R_kpc_at_D","Vobs_corrected_km_s","sigmaV_corrected_km_s","V_TG_finite_km_s","residual_km_s"]
    with (OUT/"finite_curves.csv").open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=curve_fields); w.writeheader(); w.writerows(all_curves)

    # Compact plot: best finite score over D and M/L at each finite K fraction,
    # alongside the best directly matched same-source controls for each vertical case.
    fig, axes=plt.subplots(1,3,figsize=(14,4.6),sharey=True,constrained_layout=True)
    cols={"gas_mix_hstar0p3":"tab:blue","gas_thin_hstar0p3":"tab:orange","gas_mix_hstar0p6":"tab:green"}
    for ax,v in zip(axes,VERTICALS):
        group=[r for r in good if r["vertical_case"]==v]
        prof=[min((r for r in group if r["K_over_Kcrit"]==fr),key=lambda r:r["objective_chi2_plus_existing_priors"],default=None) for fr in FRACTIONS]
        ax.plot([r["K_over_Kcrit"] for r in prof if r],[r["objective_chi2_plus_existing_priors"] for r in prof if r],"o-",color=cols[v],label="TG finite-u profile")
        s=[x for x in source_cases if x["vertical_case"]==v]
        if s:
            ax.axhline(min(x["matched_controls_same_source_node"]["NFW"]["objective_chi2_plus_priors"] for x in s),ls="--",color="tab:red",label="NFW, same-node minimum")
            ax.axhline(min(x["matched_controls_same_source_node"]["MOND"]["objective_chi2_plus_priors"] for x in s),ls=":",color="tab:purple",label="MOND, same-node minimum")
        ax.set_title(v); ax.set_xlabel("K/Kcrit (finite only)"); ax.grid(alpha=.25)
    axes[0].set_ylabel("diagonal chi-square + D/i/M/L priors")
    axes[-1].legend(fontsize=7,loc="best")
    fig.savefig(OUT/"finite_k_profile.png",dpi=170); plt.close(fig)

    def best_desc(r):
        if r is None: return "No successful finite-u point."
        c=next(x for x in source_cases if x["D_Mpc"]==r["D_Mpc"] and x["i_deg"]==r["i_deg"] and x["upsilon"]==r["upsilon"] and x["vertical_case"]==r["vertical_case"])
        return (f"D={r['D_Mpc']:.1f} Mpc, i={r['i_deg']:.0f}°, Υ*={r['upsilon']:.4f}, "
                f"{r['vertical_case']}, K/Kcrit={r['K_over_Kcrit']:.3f}; objective={r['objective_chi2_plus_existing_priors']:.3f} "
                f"(χ²diag={r['chi2_diag']:.3f}, priors={r['prior_total_existing_terms']:.3f}). Same-node controls: "
                f"NFW={c['matched_controls_same_source_node']['NFW']['objective_chi2_plus_priors']:.3f} "
                f"(optimizer_success={c['matched_controls_same_source_node']['NFW']['optimizer_success']}), "
                f"MOND={c['matched_controls_same_source_node']['MOND']['objective_chi2_plus_priors']:.3f}.")
    grow=[]
    for fr in FRACTIONS:
        x=[r["max_u_full_half_domain"] for r in good if r["K_over_Kcrit"]==fr]
        if x: grow.append((fr,float(np.median(x)),float(np.min(x)),float(np.max(x))))
    reshash=sha(OUT/"results.json")
    report=("# Finite-u TG profile at H I-supported NGC 3198 geometry\n\n"
            "**Finding — exploratory finite-u in-sample screen.** The run profiles actual positive-u solutions over the declared finite K/Kcrit grid at i=72° and D={13.8,16.6} Mpc. "
            "The singular Kcrit limit is used only to set the scale and is never scored as a finite point. The result is descriptive and establishes neither a detection nor model-selection significance.\n\n"
            "## Grid and source support\n\n"
            "The source builder is the geometry-aware HALOGAS H I plus S4G old-stellar map pipeline used by the existing nuisance envelope. THINGS adopts i=72° for NGC 3198; HALOGAS describes i≈70°. This reduced run uses the THINGS-supported 72° node only, so the HALOGAS 70° node remains a geometry limitation. D uses only the existing 13.8±1.4 Mpc SPARC prior nodes at 0 and +2σ (13.8 and 16.6 Mpc). Υ* uses existing prior/grid nodes {0.4766, 0.6, 0.762}; all three established vertical source cases are included. There are 18 source cases and eight finite K/Kcrit fractions per case.\n\n"
            "Every source case recomputes Kcrit from the principal generalized mode of the checked finite-volume operator. At each K fraction the actual B−K·diag(D) system is factored, u is checked on the full half-domain and extrapolated midplane, and the curve is evaluated at all 43 SPARC radii. The diagonal score is χ² plus the same D/i/log-M/L prior terms as the envelope. K has no prior. NFW and fixed-a0 simple MOND controls are recomputed directly from the same transformed source, distance-scaled radii, and inclination-corrected velocities at each exact nuisance node; NFW optimizer success is recorded.\n\n"
            "## Results\n\n"
            f"Completed {len(source_cases)}/18 source cases and {len(all_records)}/144 finite rows; status counts: `{json.dumps(counts,sort_keys=True)}`. "
            f"Best finite-u point: {best_desc(best)}\n\n"
            "Median/max u growth by finite K fraction (median, min, max of full-domain maximum u across successful rows):\n\n"
            + "\n".join(f"- {fr:.3f}: {med:.4g}, {mn:.4g}, {mx:.4g}" for fr,med,mn,mx in grow)+"\n\n"
            "The plot profiles over the selected D/M-L nodes separately for each vertical case; its NFW/MOND lines are the minima of same-node direct controls over those same nodes. The best point is grid-bounded if it selects K/Kcrit=0.998, D=16.6 Mpc, or an endpoint of the selected M/L subset; no interpolation beyond the declared grid is implied.\n\n"
            "## Interpretation and limits\n\n"
            "K is an added profiled galaxy parameter. The comparisons are in-sample diagonal-score screens with no radial covariance, no K prior, and a finite nuisance subset; they do not establish model preference, significance, or detection. The saved nuisance envelope has no exact i=72° rows, so the NFW/MOND controls were recomputed directly using the same transformed maps and geometry. These are not holdout predictions. NGC 4559 holdout data were not read or modified.\n\n"
            "## Runtime, hashes, and outputs\n\n"
            f"Measured CPU {results['runtime']['cpu_seconds']:.2f} s; wall {results['runtime']['wall_seconds']:.2f} s; single-thread BLAS environment: `{results['runtime']['threads']}`. Code SHA-256: `{results['runner_sha256']}`. Results SHA-256: `{reshash}`. Input hashes and all scalar/curve results are in `results.json`; `SHA256SUMS.txt` hashes task outputs.\n\n"
            "Outputs: `results.json`, `finite_profile.csv`, `finite_curves.csv`, and `finite_k_profile.png`.\n")
    (OUT/"REPORT.md").write_text(report)
    files=["finite_k_profile.py","results.json","finite_profile.csv","finite_curves.csv","finite_k_profile.png","REPORT.md"]
    (OUT/"SHA256SUMS.txt").write_text("".join(f"{sha(OUT/name)}  {name}\n" for name in files))
    print(json.dumps({"status":results["status"],"n_rows":len(all_records),"status_counts":counts,
                      "best_finite":best,"runtime":results["runtime"],"results_sha256":reshash},indent=2))


if __name__ == "__main__":
    main()

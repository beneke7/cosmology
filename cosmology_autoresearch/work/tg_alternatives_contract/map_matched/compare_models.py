#!/usr/bin/env python3
"""Same-resolved-source NGC 3198 Newtonian, NFW, MOND comparison."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import time
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares, minimize_scalar

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
ROT_ZIP = ROOT / "context/data/Rotmod_LTG.zip"
K0_CSV = ROOT / "work/tg_axisymmetric_pilot/newtonian_k0_predictions.csv"
EDGE_CSV = ROOT / "work/tg_axisymmetric_pilot/axisymmetric_branch_edge_curves.csv"
EDGE_JSON = ROOT / "work/tg_axisymmetric_pilot/axisymmetric_branch_edge_summary.json"
U_LO, U_HI = 0.6 * 10**(-0.1), 0.762
U_FIXED = (0.6, 0.762)
G_KPC = 4.30091e-6  # kpc (km/s)^2 / Msun
H0_KPC = 73.0 / 1000.0  # km/s/kpc
RHO_CRIT = 3 * H0_KPC**2 / (8 * math.pi * G_KPC)
A0_SI = 1.2e-10
KPC_M = 3.085677581491367e19
A0 = A0_SI * KPC_M / 1.0e6  # (km/s)^2/kpc


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rotmod_rows():
    with zipfile.ZipFile(ROT_ZIP) as zf:
        raw = zf.read("NGC3198_rotmod.dat").decode("utf-8")
    rows = []
    for line in raw.splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            x = [float(t) for t in line.split()]
            rows.append(x)
    return raw, np.asarray(rows, dtype=float)


def read_csv(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def chi2(obs, err, pred):
    return float(np.sum(((obs - pred) / err) ** 2))


def rms(obs, pred):
    return float(np.sqrt(np.mean((obs - pred) ** 2)))


def nfw_v2(r, logm, logc):
    m200, c = 10.0**logm, 10.0**logc
    r200 = (3 * m200 / (4 * math.pi * 200 * RHO_CRIT)) ** (1 / 3)
    rs = r200 / c
    fc = math.log1p(c) - c / (1 + c)
    x = np.maximum(r / rs, 1e-14)
    fx = np.log1p(x) - x / (1 + x)
    return G_KPC * m200 * fx / (fc * r)


def fit_1d(name, obs, err, r, gstar06, ggas, kind):
    def pred(u):
        gb = (u / 0.6) * gstar06 + ggas
        if kind == "newtonian":
            v2 = r * gb
        else:
            if np.any(gb <= 0):
                raise ValueError("MOND requires positive inward Newtonian acceleration")
            gm = 0.5 * (gb + np.sqrt(gb * gb + 4 * A0 * gb))
            v2 = r * gm
        return np.sqrt(v2)

    def objective(u):
        return chi2(obs, err, pred(float(u)))

    opt = minimize_scalar(objective, bounds=(U_LO, U_HI), method="bounded",
                          options={"xatol": 1e-12, "maxiter": 1000})
    candidates = [(float(opt.x), float(opt.fun))]
    candidates.extend((u, objective(u)) for u in (*U_FIXED, U_LO))
    u_best, score = min(candidates, key=lambda x: x[1])
    return {"name": name, "n_free": 1, "upsilon": u_best,
            "chi2_diag": score, "rms_km_s": rms(obs, pred(u_best)),
            "predictions": pred(u_best), "fixed": {
                str(u): {"upsilon": u, "chi2_diag": objective(u),
                         "rms_km_s": rms(obs, pred(u)), "predictions": pred(u)}
                for u in U_FIXED}}


def fit_nfw(obs, err, r, gstar06, ggas):
    # Optimization domains are explicit boxes in log10(M200/Msun), log10(c),
    # and linear Upsilon. No concentration-mass relation or prior is imposed.
    lo = np.array([9.0, 0.0, U_LO])
    hi = np.array([14.0, math.log10(50.0), U_HI])

    def predict(theta):
        logm, logc, u = theta
        gb = (u / 0.6) * gstar06 + ggas
        v2 = r * gb + nfw_v2(r, logm, logc)
        return np.sqrt(np.maximum(v2, 0.0))

    def resid(theta):
        return (predict(theta) - obs) / err

    # Deterministic multistart least-squares protects against a local basin.
    starts = [(10.0, 0.7, 0.6), (11.0, 0.9, 0.6), (12.0, 1.0, 0.6),
              (13.0, 0.8, 0.7), (11.5, 1.4, 0.55), (13.5, 1.4, 0.75)]
    fits = [least_squares(resid, np.asarray(s), bounds=(lo, hi),
                          x_scale=np.array([1.0, 0.5, 0.1]),
                          ftol=1e-12, xtol=1e-12, gtol=1e-12,
                          max_nfev=4000) for s in starts]
    best = min(fits, key=lambda x: float(np.dot(x.fun, x.fun)))
    p = predict(best.x)
    logm, logc, u = map(float, best.x)
    return {"name": "Newtonian baryons + NFW halo", "n_free": 3,
            "upsilon": u, "log10_M200_Msun": logm,
            "M200_Msun": 10**logm, "concentration_c200": 10**logc,
            "log10_c200": logc, "chi2_diag": chi2(obs, err, p),
            "rms_km_s": rms(obs, p), "predictions": p,
            "optimizer_success": bool(best.success),
            "optimizer_message": best.message,
            "optimizer_starts": len(starts)}


def pack_model(model):
    return {k: v for k, v in model.items() if k != "predictions" and k != "fixed"}


def main():
    wall0, cpu0 = time.perf_counter(), time.process_time()
    _, rot = rotmod_rows()
    k0 = read_csv(K0_CSV)
    edge = read_csv(EDGE_CSV)
    if len(rot) != 43 or len(k0) != 43 or len(edge) != 43:
        raise RuntimeError(f"expected 43 matched rows, got {len(rot)}, {len(k0)}, {len(edge)}")
    r, obs, err = rot[:, 0], rot[:, 1], rot[:, 2]
    kr = np.asarray([float(x["R_kpc"]) for x in k0])
    er = np.asarray([float(x["R_kpc"]) for x in edge])
    if not (np.allclose(r, kr, rtol=0, atol=1e-9) and np.allclose(r, er, rtol=0, atol=1e-9)):
        raise RuntimeError("SPARC and map/branch radii are not identical")
    for col, arr in [("Vobs_km_s", obs), ("sigmaV_km_s", err)]:
        got = np.asarray([float(x[col]) for x in k0])
        if not np.allclose(got, arr, rtol=0, atol=1e-9):
            raise RuntimeError(f"map CSV {col} fails exact SPARC row check")

    vstar06 = np.asarray([float(x["Vstar_map_K0_km_s"]) for x in k0])
    vtot06 = np.asarray([float(x["Vbar_map_HALOGAS_K0_km_s"]) for x in k0])
    vtot762_map = np.asarray([float(x["K0_resolved_baryons_U0p762_km_s"]) for x in edge])
    # Frozen file's unit-M/L axisymmetric star speed is at U=0.6. Since the
    # component total is positive, recover gas acceleration as a signed v^2
    # difference. The negative values are retained, never abs'ed away.
    gas_v2 = vtot06**2 - vstar06**2
    gstar06 = vstar06**2 / r
    ggas = gas_v2 / r
    gbar06 = gstar06 + ggas
    # Independent second-M/L total checks linear stellar acceleration scaling.
    vtot_hi_rebuilt = np.sqrt((U_HI / 0.6) * vstar06**2 + gas_v2)
    map_scaling_max_abs = float(np.max(np.abs(vtot_hi_rebuilt - vtot762_map)))
    if np.any(vtot06 <= 0) or np.any(vtot762_map <= 0):
        raise RuntimeError("total baryonic speed is not positive; signed recovery ambiguous")

    models = []
    models.append(fit_1d("Resolved-map Newtonian baryons only", obs, err, r,
                         gstar06, ggas, "newtonian"))
    models.append(fit_nfw(obs, err, r, gstar06, ggas))
    models.append(fit_1d("Isolated simple-MOND, fixed a0", obs, err, r,
                         gstar06, ggas, "mond"))

    # Edge is a non-finite-normalization TG limiting curve; summarize all
    # reported fine-grid edge curves and select lowest diagonal chi-square.
    edge_cols = [c for c in edge[0] if c.startswith("edge_Upsilon_")]
    edge_models = []
    for col in edge_cols:
        curve = np.asarray([float(x[col]) for x in edge])
        u = float(col.removeprefix("edge_Upsilon_").removesuffix("_km_s"))
        edge_models.append({"name": "TG branch edge limit", "upsilon": u,
                            "curve_column": col,
                            "chi2_diag": chi2(obs, err, curve),
                            "rms_km_s": rms(obs, curve), "predictions": curve})
    best_edge = min(edge_models, key=lambda m: m["chi2_diag"])

    # Include fixed-M/L NFW sensitivity at both requested values.
    nfw_fixed = []
    for uf in U_FIXED:
        # Profile only the halo's two parameters with stellar M/L held fixed.
        def vpred(theta):
            gb = (uf / 0.6) * gstar06 + ggas
            return np.sqrt(np.maximum(r * gb + nfw_v2(r, *theta), 0.0))
        def rr(theta): return (vpred(theta) - obs) / err
        sols = [least_squares(rr, np.asarray(s), bounds=([9.0, 0.0], [14.0, math.log10(50)]),
                              x_scale=[1.0, 0.5], max_nfev=4000,
                              ftol=1e-12, xtol=1e-12, gtol=1e-12)
                for s in [(10, .7), (11, 1), (12, 1.4), (13, .8)]]
        sol = min(sols, key=lambda x: float(np.dot(x.fun, x.fun)))
        pp = vpred(sol.x)
        nfw_fixed.append({"upsilon_fixed": uf, "n_free": 2,
                          "log10_M200_Msun": float(sol.x[0]),
                          "M200_Msun": float(10**sol.x[0]),
                          "log10_c200": float(sol.x[1]),
                          "concentration_c200": float(10**sol.x[1]),
                          "chi2_diag": chi2(obs, err, pp), "rms_km_s": rms(obs, pp)})

    # Ensure algebraic MOND has the required inward Newtonian field across all
    # allowed stellar M/L values, not merely at the fitted point.
    gb_lo = (U_LO / .6) * gstar06 + ggas
    gb_hi = (U_HI / .6) * gstar06 + ggas
    mon_diagnostics = {"minimum_net_inward_gN_over_full_U_domain_km2_s2_kpc": float(min(gb_lo.min(), gb_hi.min())),
                       "minimum_gN_at_U_lower_km2_s2_kpc": float(gb_lo.min()),
                       "minimum_gN_at_U_upper_km2_s2_kpc": float(gb_hi.min()),
                       "all_positive_for_entire_domain": bool(np.all(gb_lo > 0) and np.all(gb_hi > 0)),
                       "signed_gas_acceleration_negative_point_count": int(np.count_nonzero(ggas < 0)),
                       "minimum_signed_gas_g_km2_s2_kpc": float(ggas.min()),
                       "maximum_signed_gas_g_km2_s2_kpc": float(ggas.max())}
    if not mon_diagnostics["all_positive_for_entire_domain"]:
        raise RuntimeError("MOND source has non-positive inward Newtonian gN")

    # Machine-readable per-model predictions.
    pred_models = [(m["name"] + f" U={m['upsilon']:.6f}", m["predictions"]) for m in models]
    for m in models:
        if "fixed" in m:
            pred_models += [(m["name"] + f" U={float(u):.3f} fixed", d["predictions"])
                            for u, d in m["fixed"].items()]
    pred_models.append((f"TG branch edge limit U={best_edge['upsilon']:.6f}", best_edge["predictions"]))
    with (OUT / "predictions.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["R_kpc", "Vobs_km_s", "sigmaV_km_s", "signed_gas_v2_km2_s2",
                    "gN_best_baryon_km2_s2_kpc", *[n for n, _ in pred_models]])
        bestb = models[0]["upsilon"]
        bestg = (bestb/.6)*gstar06 + ggas
        for i in range(43):
            w.writerow([r[i], obs[i], err[i], gas_v2[i], bestg[i],
                        *[p[i] for _, p in pred_models]])

    results = {
        "status": "completed_same_resolved_source_diagonal_error_descriptive_comparison_no_preference_claim",
        "n_observations": 43,
        "inputs": {"rotmod_zip": str(ROT_ZIP.relative_to(ROOT)), "rotmod_sha256": sha256(ROT_ZIP),
                   "k0_predictions_csv": str(K0_CSV.relative_to(ROOT)), "k0_predictions_sha256": sha256(K0_CSV),
                   "branch_edge_csv": str(EDGE_CSV.relative_to(ROOT)), "branch_edge_csv_sha256": sha256(EDGE_CSV),
                   "branch_edge_summary_sha256": sha256(EDGE_JSON)},
        "source_recovery": {"stellar_reference_upsilon": 0.6,
                            "gas_v2_equation": "Vgas_signed^2 = Vbar_map_HALOGAS_K0^2 - Vstar_map_K0^2",
                            "gas_v2_negative_count": int(np.count_nonzero(gas_v2 < 0)),
                            "gas_v2_min_km2_s2": float(gas_v2.min()),
                            "gas_v2_max_km2_s2": float(gas_v2.max()),
                            "second_ML_scaling_check_max_abs_km_s": map_scaling_max_abs,
                            "reference_profile": "same resolved HALOGAS map-derived K=0 stellar+gas field for all; stellar term scales linearly in Upsilon"},
        "likelihood": {"form": "chi2=sum(((Vobs-Vmodel)/errV)^2)", "error_model": "diagonal SPARC errV only; no full covariance",
                       "dof_not_used_for_claim": True, "radii_kpc_min": float(r.min()), "radii_kpc_max": float(r.max())},
        "units_and_constants": {"r": "kpc", "velocity": "km/s", "g": "(km/s)^2/kpc", "G": G_KPC,
                                "H0_km_s_Mpc": 73.0, "rho_critical_Msun_kpc3": RHO_CRIT,
                                "MOND_a0_m_s2": A0_SI, "MOND_a0_km2_s2_kpc": A0},
        "domains": {"upsilon_Msun_Lsun_3p6": [U_LO, U_HI], "fixed_upsilon_values": list(U_FIXED),
                    "NFW_log10_M200_Msun": [9.0, 14.0], "NFW_c200": [1.0, 50.0],
                    "NFW_nuisance_counts": {"profile_upsilon": 3, "fixed_upsilon": 2},
                    "Newtonian_baryons_nuisance_count": 1,
                    "isolated_simple_MOND_nuisance_count": 1},
        "equations": {"Newtonian_baryons": "V^2 = r*(U/0.6*gstar_0.6 + ggas_signed)",
                      "NFW": "V^2 = Vbar^2 + G*M200*f(r/rs)/(r*f(c)); r200=(3M200/(800*pi*rho_crit))^(1/3), rs=r200/c; f(x)=ln(1+x)-x/(1+x)",
                      "MOND": "gN=(U/0.6*gstar_0.6+ggas_signed)>0; mu=x/(1+x), x=g/a0; g=(gN+sqrt(gN^2+4*a0*gN))/2; V^2=r*g"},
        "models": [pack_model(m) | {"fixed_ML_sensitivity": {k: {kk: vv for kk, vv in v.items() if kk != "predictions"}
                                                        for k, v in m.get("fixed", {}).items()}}
                   for m in models],
        "nfw_fixed_ML_sensitivity": nfw_fixed,
        "tg_branch_edge_limit_only": {k: v for k, v in best_edge.items() if k != "predictions"} |
            {"upsilon_selected_from_grid_points": len(edge_models),
             "K_treatment": "fixed to Kcrit; eigenmode limiting curve, not fitted"},
        "all_fine_grid_branch_edges": [{k: v for k, v in m.items() if k != "predictions"} for m in edge_models],
        "mond_net_inward_check": mon_diagnostics,
        "runtime": {"cpu_seconds": time.process_time() - cpu0,
                    "wall_seconds": time.perf_counter() - wall0,
                    "single_thread_env": {k: os.environ.get(k) for k in
                        ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}},
        "limitations": ["diagonal velocity errors only; no full SPARC covariance",
                        "the map-derived source is an axisymmetric pilot source, not a unique 3D baryon density",
                        "MOND is the isolated one-dimensional algebraic simple-mu relation; no EFE",
                        "TG edge is a non-finite-normalization branch limit, shown as a limiting curve only",
                        "scores are descriptive; no statistical model preference or detection claim"]}
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.errorbar(r, obs, yerr=err, fmt="o", ms=3, color="black", alpha=.75, label="SPARC NGC 3198")
    for m, style in zip(models, ["-", "-", "-"]):
        ax.plot(r, m["predictions"], style, lw=2, label=f"{m['name']} (U={m['upsilon']:.3f})")
    ax.plot(r, best_edge["predictions"], "--", lw=2,
            label=f"TG branch edge limit only (U={best_edge['upsilon']:.3f})")
    ax.set(xlabel="Radius (kpc)", ylabel="Circular speed (km s$^{-1}$)",
           title="NGC 3198: same resolved baryon source, diagonal-error fits")
    ax.grid(alpha=.25); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(OUT / "map_matched_comparison.png", dpi=180)
    results["runtime"] = {"cpu_seconds": time.process_time() - cpu0,
                          "wall_seconds": time.perf_counter() - wall0,
                          "single_thread_env": {k: os.environ.get(k) for k in
                              ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}}
    (OUT / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    print(json.dumps(results["runtime"], sort_keys=True))
    print(json.dumps({"models": results["models"], "nfw_fixed": nfw_fixed,
                      "edge": results["tg_branch_edge_limit_only"],
                      "mond_check": mon_diagnostics, "map_scaling_max": map_scaling_max_abs}, indent=2))


if __name__ == "__main__":
    main()

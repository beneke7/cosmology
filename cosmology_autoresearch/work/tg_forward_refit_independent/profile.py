#!/usr/bin/env python3
"""Independent Pszota–Ván NGC 3198 spherical-surrogate profile.

Reads only the local public Rotmod archive and the paper. It does not import
or inspect the campaign's other isolated-solver implementation.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import platform
import resource
import time
import zipfile
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
ARCHIVE = PROJECT / "context" / "data" / "Rotmod_LTG.zip"
PAPER = PROJECT / "context" / "user_provided" / "pszota_van_2024.pdf"
MEMBER = "NGC3198_rotmod.dat"
G = 4.30091e-6  # kpc (km/s)^2 / Msun
UPS_RANGE = (0.2, 1.0)
UPS_GRID = sorted(set(np.round(np.linspace(*UPS_RANGE, 9), 12).tolist() + [0.762]))
K_DOMAIN = (0.0, 1.0e-2)  # s^2/km^2; 294 x the paper's reported central K
INNER_POWERS = (1.0, 1.5, 2.0)  # regular v ~ r^p continuations below first datum
TAIL_LENGTH_KPC = 5.0
TAIL_EXTENTS = 12.0  # impose exterior Robin match after 12 taper lengths
FRACTIONS = sorted(set([0.0, 0.9999] + np.geomspace(1e-5, 0.99, 42).tolist()
                       + np.linspace(0.02, 0.98, 19).tolist()))
CPU_LIMIT_SECONDS = 1100.0
CPU_START = time.process_time()
WALL_START = time.perf_counter()
N_SOLVES = 0
DATA = None
PAPER_SHA = None
ARCHIVE_SHA = None


def check_budget() -> None:
    if time.process_time() - CPU_START > CPU_LIMIT_SECONDS:
        raise TimeoutError("numerical CPU budget reached")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_local_rotmod() -> dict[str, np.ndarray]:
    """Separate line-oriented parser; assert the expected eight-column schema."""
    with zipfile.ZipFile(ARCHIVE) as zf:
        raw = zf.read(MEMBER)
    rows = []
    for line_no, line in enumerate(raw.decode("ascii").splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 8:
            raise ValueError(f"{MEMBER}:{line_no}: expected 8 fields, found {len(fields)}")
        rows.append([float(field) for field in fields])
    a = np.asarray(rows, dtype=float)
    if a.shape != (43, 8):
        raise ValueError(f"Expected all 43 SPARC points, found shape {a.shape}")
    if np.any(np.diff(a[:, 0]) <= 0) or np.any(a[:, 2] <= 0):
        raise ValueError("radii must increase and every velocity uncertainty must be positive")
    return {
        "raw": raw,
        "r": a[:, 0], "vobs": a[:, 1], "sigma": a[:, 2],
        "vgas": a[:, 3], "vdisk": a[:, 4],
        "sbdisk": a[:, 6], "sbbul": a[:, 7],
    }


class Source:
    """Signed spherical inversion of the disk and gas force components."""

    def __init__(self, p_inner: float, tail_length: float = TAIL_LENGTH_KPC):
        self.p = float(p_inner)
        self.L = float(tail_length)
        self.r = DATA["r"]
        self.r0, self.rd = float(self.r[0]), float(self.r[-1])
        self.disk = PchipInterpolator(self.r, DATA["vdisk"], extrapolate=False)
        self.gas = PchipInterpolator(self.r, DATA["vgas"], extrapolate=False)
        self.vd0 = float(DATA["vdisk"][0])
        self.vg0 = float(DATA["vgas"][0])
        self.m_center = 2.0 * self.p - 2.0
        self.R = self.rd + TAIL_EXTENTS * self.L
        # Coefficient at r0 in rho_center(r) = rho0 * (r/r0)^m.
        mult = 1.0 + 2.0 * self.p
        self.rho_d0 = mult * self.vd0**2 / (4.0 * np.pi * G * self.r0**2)
        self.rho_g0 = mult * self.vg0 * abs(self.vg0) / (4.0 * np.pi * G * self.r0**2)

    def rho_parts(self, radius: float) -> tuple[float, float]:
        r = float(radius)
        if r < self.r0:
            scale = (max(r, 0.0) / self.r0) ** self.m_center
            return self.rho_d0 * scale, self.rho_g0 * scale
        if r <= self.rd:
            vd, vdp = float(self.disk(r)), float(self.disk.derivative()(r))
            vg, vgp = float(self.gas(r)), float(self.gas.derivative()(r))
            # Signed gas mass proxy q(v)=v|v|; the disk curve is positive.
            rdens = (vd * vd + 2.0 * r * vd * vdp) / (4.0 * np.pi * G * r * r)
            rgden = (vg * abs(vg) + 2.0 * r * abs(vg) * vgp) / (4.0 * np.pi * G * r * r)
            return rdens, rgden
        # Explicit signed exponential continuation of the endpoint pseudo-density.
        dd, gg = self.rho_parts(self.rd)
        factor = np.exp(-(r - self.rd) / self.L)
        return dd * factor, gg * factor

    def rho(self, radius: float, ups: float) -> float:
        d, g = self.rho_parts(radius)
        return ups * d + g

    def center_coefficient(self, ups: float) -> float:
        return ups * self.rho_d0 + self.rho_g0

    def integration_boundaries(self) -> list[float]:
        eps = min(1.0e-7, self.r0 * 1.0e-6)
        return [eps, self.r0, self.rd, self.R]


def eval_segments(segments, radius: float) -> np.ndarray:
    r = float(radius)
    for a, b, sol in segments:
        if a - 1e-13 <= r <= b + 1e-13:
            return np.asarray(sol.sol(min(max(r, a), b)), dtype=float)
    raise ValueError(f"radius {r} is not covered by the integrated source domain")


def solve_x(source: Source, kappa: float, ups: float, *, rtol=2e-9,
            atol=1e-9, max_step=0.25):
    """Integrate x=(1-w)/K, whose K=0 equation is exactly Poisson."""
    global N_SOLVES
    N_SOLVES += 1
    if not (K_DOMAIN[0] <= kappa <= K_DOMAIN[1]):
        return {"valid": False, "reason": "outside_declared_K_domain"}
    bounds = source.integration_boundaries()
    eps, r0, rd, R = bounds
    m = source.m_center
    q0 = 4.0 * np.pi * G * source.center_coefficient(ups)
    denom = (m + 3.0) * (m + 2.0)
    x0 = q0 * eps ** (m + 2.0) / (r0**m * denom)
    xp0 = q0 * eps ** (m + 1.0) / (r0**m * (m + 3.0))
    state = np.array([x0, xp0], dtype=float)
    segments = []

    def rhs(rad, y):
        rho = source.rho(rad, ups)
        return np.array([y[1], 4.0 * np.pi * G * rho * (1.0 - kappa * y[0]) - 2.0 * y[1] / rad])

    for left, right in zip(bounds[:-1], bounds[1:]):
        sol = solve_ivp(rhs, (left, right), state, method="DOP853",
                        rtol=rtol, atol=atol, max_step=max_step, dense_output=True)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return {"valid": False, "reason": "integration_failure", "message": sol.message}
        segments.append((left, right, sol))
        state = sol.y[:, -1]
    # For K>0, w=1-Kx and the exterior u->1 normalization is fixed at R.
    xR, xpR = state
    norm = 1.0 - kappa * (xR + R * xpR)
    if not np.isfinite(norm):
        return {"valid": False, "reason": "nonfinite_exterior_normalization"}
    dense_r = np.unique(np.concatenate([
        np.linspace(eps, r0, 120), np.linspace(r0, rd, 1200),
        np.linspace(rd, R, 240),
    ]))
    dense_state = np.array([eval_segments(segments, rr) for rr in dense_r])
    den = 1.0 - kappa * dense_state[:, 0]
    if kappa == 0.0:
        min_u_numerator = 1.0
        min_u = 1.0
    else:
        min_u_numerator = float(np.min(den))
        min_u = min_u_numerator / norm if norm != 0.0 else float("nan")
        if norm <= 0.0:
            return {"valid": False, "reason": "normalization_threshold", "normalization": float(norm)}
        if min_u_numerator <= 0.0:
            return {"valid": False, "reason": "u_positivity_threshold",
                    "normalization": float(norm), "min_w": min_u_numerator}
    states = np.array([eval_segments(segments, rr) for rr in DATA["r"]])
    xd = states[:, 1]
    denominator = 1.0 - kappa * states[:, 0]
    model_v2 = DATA["r"] * xd / denominator if kappa != 0.0 else DATA["r"] * xd
    if np.any(~np.isfinite(model_v2)):
        return {"valid": False, "reason": "nonfinite_model_velocity"}
    if np.any(model_v2 <= 0.0):
        return {"valid": False, "reason": "nonpositive_model_v2",
                "minimum_model_v2": float(np.min(model_v2))}
    model_v = np.sqrt(model_v2)
    chi2 = float(np.sum(((DATA["vobs"] - model_v) / DATA["sigma"]) ** 2))
    return {
        "valid": True, "chi2": chi2, "model_v_kms": model_v,
        "model_v2_kms2": model_v2, "segments": segments,
        "normalization": float(norm), "min_w": float(min_u_numerator),
        "min_normalized_u": float(min_u), "R_kpc": R,
        "nfev": int(sum(sol.nfev for _, _, sol in segments)),
    }


def find_first_threshold(source: Source, ups: float) -> dict:
    """Find the first connected-branch positivity/normalization loss in K."""
    scan = np.geomspace(1e-11, K_DOMAIN[1], 66)
    last_valid = 0.0
    last_good = solve_x(source, 0.0, ups)
    if not last_good["valid"]:
        return {"threshold_found": True, "kind": last_good["reason"], "K_lower": 0.0, "K_upper": 0.0}
    for k in scan:
        check_budget()
        result = solve_x(source, float(k), ups)
        if not result["valid"]:
            reason = result["reason"]
            if reason == "integration_failure":
                return {"threshold_found": False, "status": "numerical_failure_before_threshold",
                        "last_valid_K": last_valid, "first_failed_K": float(k), "reason": reason}
            lo, hi = last_valid, float(k)
            # Geometric bisection once a positive lower bracket exists.
            for _ in range(28):
                check_budget()
                mid = np.sqrt(lo * hi) if lo > 0.0 else hi * 0.5
                trial = solve_x(source, float(mid), ups)
                if trial["valid"]:
                    lo = mid
                elif trial["reason"] in ("normalization_threshold", "u_positivity_threshold", "nonpositive_model_v2"):
                    hi = mid
                else:
                    return {"threshold_found": False, "status": "threshold_refinement_numerical_failure",
                            "last_valid_K": lo, "first_failed_K": hi, "reason": trial["reason"]}
            near = solve_x(source, float(lo), ups, rtol=1e-11, atol=1e-11, max_step=0.08)
            return {
                "threshold_found": True, "kind": reason, "K_just_below": float(lo),
                "K_just_above": float(hi), "relative_bracket_width": float((hi - lo) / hi),
                "near_threshold_refined_valid": bool(near["valid"]),
                "near_threshold_normalization": near.get("normalization"),
                "near_threshold_min_w": near.get("min_w"),
            }
        last_valid = float(k)
    return {"threshold_found": False, "status": "no_threshold_in_declared_domain",
            "K_domain_max": K_DOMAIN[1]}


def score(source: Source, kappa: float, ups: float, *, rtol=2e-9,
          atol=1e-9, max_step=0.25) -> dict:
    out = solve_x(source, kappa, ups, rtol=rtol, atol=atol, max_step=max_step)
    if out["valid"]:
        return {key: out[key] for key in ("valid", "chi2", "normalization", "min_w",
                                          "min_normalized_u", "R_kpc", "nfev")}
    return {key: out.get(key) for key in ("valid", "reason", "normalization",
                                          "min_w", "minimum_model_v2")}


def k0_recovery(source: Source, ups: float) -> dict:
    result = solve_x(source, 0.0, ups)
    if not result["valid"]:
        return result
    target_v2 = ups * DATA["vdisk"]**2 + DATA["vgas"] * np.abs(DATA["vgas"])
    err = result["model_v2_kms2"] - target_v2
    return {
        "valid": True,
        "max_abs_v2_recovery_error_km2_s2": float(np.max(np.abs(err))),
        "rms_v2_recovery_error_km2_s2": float(np.sqrt(np.mean(err**2))),
        "chi2": result["chi2"],
    }


def riccati_crosscheck(source: Source, kappa: float, ups: float, x_result: dict) -> dict:
    """Independent direct integration of the original acceleration Riccati ODE."""
    eps, r0, rd, R = source.integration_boundaries()
    m = source.m_center
    q0 = 4.0 * np.pi * G * source.center_coefficient(ups)
    g0 = -q0 * eps ** (m + 1.0) / (r0**m * (m + 3.0))
    state = np.array([g0])
    segments = []

    def rhs(rad, y):
        g = y[0]
        return np.array([-2.0 * g / rad - 4.0 * np.pi * G * source.rho(rad, ups) - kappa * g * g])

    for left, right in zip((eps, r0, rd), (r0, rd, R)):
        sol = solve_ivp(rhs, (left, right), state, method="DOP853",
                        rtol=2e-11, atol=1e-11, max_step=0.08, dense_output=True)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return {"status": "riccati_integration_failure", "message": sol.message}
        segments.append((left, right, sol))
        state = sol.y[:, -1]
    radii = DATA["r"]
    g_ric = np.array([eval_segments(segments, rr)[0] for rr in radii])
    g_x = []
    for rr in radii:
        x, xp = eval_segments(x_result["segments"], rr)
        g_x.append(-xp / (1.0 - kappa * x))
    g_x = np.asarray(g_x)
    absdiff = np.abs(g_ric - g_x)
    rel = absdiff / np.maximum(np.abs(g_x), 1e-14)
    gR = float(state[0])
    C = float(-1.0 / (gR * R * R) - kappa / R) if gR != 0 else float("inf")
    return {
        "status": "independent Riccati integration completed",
        "max_abs_g_difference_km_s2": float(np.max(absdiff)),
        "max_relative_g_difference": float(np.max(rel)),
        "g_at_R_km_s2": gR,
        "vacuum_C_km2_s2_per_kpc2": C,
        "C_positive_attractive_vacuum_branch": bool(C > 0.0),
    }


def main() -> None:
    global DATA, PAPER_SHA, ARCHIVE_SHA
    t0 = time.perf_counter()
    DATA = read_local_rotmod()
    PAPER_SHA = sha(PAPER.read_bytes())
    ARCHIVE_SHA = sha(ARCHIVE.read_bytes())
    source_hash = sha(DATA["raw"])
    paper_k = 3.40e-5
    paper_ups = 0.762
    k0_tests = {}

    sources = {p: Source(p, TAIL_LENGTH_KPC) for p in INNER_POWERS}
    density_diag = {}
    for p, source in sources.items():
        rr = np.linspace(source.r0, source.rd, 20_001)
        rho = np.array([source.rho(x, paper_ups) for x in rr])
        neg = np.flatnonzero(rho < 0.0)
        intervals = []
        if len(neg):
            start = prev = int(neg[0])
            for idx in neg[1:]:
                idx = int(idx)
                if idx > prev + 1:
                    intervals.append([float(rr[start]), float(rr[prev])])
                    start = idx
                prev = idx
            intervals.append([float(rr[start]), float(rr[prev])])
        density_diag[str(p)] = {
            "negative_intervals_kpc": intervals,
            "rho_min_Msun_kpc3": float(np.min(rho)),
            "r_at_min_kpc": float(rr[int(np.argmin(rho))]),
            "negative_fraction_on_grid": float(np.mean(rho < 0.0)),
        }
        k0_tests[str(p)] = {
            "K0_exact_recovery": k0_recovery(source, paper_ups),
            "small_K": {
                str(k): score(source, k, paper_ups)
                for k in (1e-12, 1e-10, 1e-8, paper_k)
            },
        }

    # Locate first normalization / positivity threshold on the connected K>=0 branch.
    thresholds = {}
    for p, source in sources.items():
        thresholds[str(p)] = {}
        for ups in UPS_GRID:
            check_budget()
            thresholds[str(p)][str(ups)] = find_first_threshold(source, ups)

    # Broad deterministic two-parameter profile over the entire positive branch.
    profile_rows = []
    best_by_core = {}
    for p, source in sources.items():
        local_best = None
        for ups in UPS_GRID:
            check_budget()
            t = thresholds[str(p)][str(ups)]
            if t.get("threshold_found"):
                kmax = float(t["K_just_below"])
                branch_status = "first_positive_branch"
            elif t.get("status") in ("numerical_failure_before_threshold", "threshold_refinement_numerical_failure"):
                kmax = float(t.get("last_valid_K", t.get("K_lower", 0.0)))
                branch_status = "numerically_verified_branch_only"
            else:
                kmax = K_DOMAIN[1]
                branch_status = t.get("status", "domain_cap")
            for frac in FRACTIONS:
                check_budget()
                k = float(frac * kmax)
                res = score(source, k, ups)
                row = {
                    "p_inner": p, "Upsilon": ups, "K_s2_per_km2": k,
                    "fraction_of_first_branch": float(frac), "branch_status": branch_status,
                    **res,
                }
                profile_rows.append(row)
                if row["valid"] and (local_best is None or row["chi2"] < local_best["chi2"]):
                    local_best = dict(row)
        best_by_core[str(p)] = local_best

    # Refine around each coarse profile minimum in branch-fraction / Upsilon coordinates.
    refined = {}
    for p, source in sources.items():
        check_budget()
        best = best_by_core[str(p)]
        crit_u = np.array(UPS_GRID, dtype=float)
        crit_k = np.array([
            float(thresholds[str(p)][str(u)].get("K_just_below", K_DOMAIN[1]))
            if thresholds[str(p)][str(u)].get("threshold_found") else K_DOMAIN[1]
            for u in UPS_GRID
        ])
        if best is None:
            refined[str(p)] = {"status": "no_valid_profile_point"}
            continue

        def kcrit_interp(u):
            return float(np.interp(u, crit_u, crit_k))

        def objective(theta):
            fraction, ups = float(theta[0]), float(theta[1])
            k = max(0.0, fraction) * kcrit_interp(ups)
            res = solve_x(source, k, ups)
            if not res["valid"]:
                return 1e12
            return res["chi2"]

        frac0 = float(best["fraction_of_first_branch"])
        result = minimize(objective, np.array([frac0, best["Upsilon"]]),
                          method="Powell", bounds=[(0.0, 0.9999), UPS_RANGE],
                          options={"xtol": 2e-5, "ftol": 1e-9, "maxiter": 90})
        k_ref = float(result.x[0] * kcrit_interp(result.x[1]))
        fit_ref = solve_x(source, k_ref, float(result.x[1]))
        if not fit_ref["valid"]:
            ref_candidate = best
            refine_status = "local_refinement_invalid; retained best grid point"
            k_ref = float(best["K_s2_per_km2"])
            ups_ref = float(best["Upsilon"])
        else:
            ref_candidate = {
                "p_inner": p, "Upsilon": float(result.x[1]), "K_s2_per_km2": k_ref,
                "chi2": fit_ref["chi2"], "fraction_of_first_branch": float(result.x[0]),
            }
            refine_status = "Powell local refinement completed"
            ups_ref = float(result.x[1])

        # Independently refine the exact Newtonian K=0 slice in Upsilon.
        zero_opt = minimize_scalar(lambda u: score(source, 0.0, float(u))["chi2"],
                                   bounds=UPS_RANGE, method="bounded",
                                   options={"xatol": 1e-9, "maxiter": 100})
        zero_point = {"K_s2_per_km2": 0.0, "Upsilon": float(zero_opt.x),
                      "chi2": float(zero_opt.fun), "success": bool(zero_opt.success)}
        candidates = [best, ref_candidate, zero_point]
        candidate = min(candidates, key=lambda item: item["chi2"])
        kc, uc = float(candidate["K_s2_per_km2"]), float(candidate["Upsilon"])
        tight = solve_x(source, kc, uc, rtol=2e-11, atol=1e-11, max_step=0.06)
        if not tight["valid"]:
            tight_record = {"valid": False, "reason": tight["reason"]}
            ric = {"status": "skipped_invalid_candidate"}
        else:
            tight_record = {
                "valid": True, "chi2_tight": float(tight["chi2"]),
                "normalization_tight": float(tight["normalization"]),
                "min_w_tight": float(tight["min_w"]),
                "max_nfev": int(tight["nfev"]),
                "delta_chi2_vs_standard_tolerance": float(tight["chi2"] - candidate["chi2"]),
            }
            ric = riccati_crosscheck(source, kc, uc, tight) if kc > 0.0 else {"status": "K=0 Newtonian limit; Riccati check not needed"}
        # Explicit fixed-parameter taper sensitivity from 5 to 10 kpc.
        source_10 = Source(p, 10.0)
        tail = solve_x(source_10, kc, uc, rtol=2e-11, atol=1e-11, max_step=0.08)
        tail_sensitivity = {
            "taper_length_kpc": 10.0,
            "valid": bool(tail["valid"]),
            "chi2_same_K_Upsilon": float(tail["chi2"]) if tail["valid"] else None,
            "chi2_change_from_5kpc": float(tail["chi2"] - tight["chi2"]) if tail["valid"] and tight["valid"] else None,
            "reason_if_invalid": tail.get("reason"),
        }
        refined[str(p)] = {
            "coarse_grid_best": best,
            "local_refinement_status": refine_status,
            "refined_candidate": ref_candidate,
            "K0_slice_minimum": zero_point,
            "selected_candidate": candidate,
            "tight_tolerance_verification": tight_record,
            "independent_Riccati_crosscheck": ric,
            "fixed_parameter_outer_taper_sensitivity": tail_sensitivity,
            "positivity_threshold_at_paper_Upsilon": thresholds[str(p)][str(paper_ups)],
        }

    # Save profile surface for full inspection, including invalid points.
    csv_path = HERE / "profile_surface.csv"
    fields = ["p_inner", "Upsilon", "K_s2_per_km2", "fraction_of_first_branch",
              "branch_status", "valid", "chi2", "reason", "normalization", "min_w",
              "min_normalized_u", "R_kpc", "nfev"]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(profile_rows)

    cpu_s = time.process_time() - CPU_START
    wall_s = time.perf_counter() - WALL_START
    result = {
        "status": "completed independent fixed-source spherical-surrogate profile",
        "scope": "all 43 NGC3198 Rotmod points; signed Vgas force; no density clipping; direct velocity chi-square; no observed-endpoint boundary condition",
        "source": {
            "paper_pdf_sha256": PAPER_SHA,
            "rotmod_archive_sha256": ARCHIVE_SHA,
            "rotmod_member": MEMBER,
            "rotmod_member_sha256": source_hash,
            "N_data": int(len(DATA["r"])),
            "r_min_kpc": float(DATA["r"][0]),
            "r_max_kpc": float(DATA["r"][-1]),
            "outer_observed_velocity_km_s": float(DATA["vobs"][-1]),
            "outer_velocity_error_km_s": float(DATA["sigma"][-1]),
            "outer_vgas_km_s": float(DATA["vgas"][-1]),
            "outer_vdisk_km_s": float(DATA["vdisk"][-1]),
            "signed_gas_rows": int(np.count_nonzero(DATA["vgas"] < 0)),
            "source_inversion": "rho=[Upsilon d(r Vdisk^2)/dr + d(r Vgas*abs(Vgas))/dr]/(4 pi G r^2); negative intervals retained",
            "interpolation": "PCHIP of tabulated component velocities; differentiate mass proxy analytically",
            "inner_continuations": {str(p): "V_component(r)=V_component(r0)*(r/r0)^p, p=" + str(p) for p in INNER_POWERS},
            "outer_taper": "each signed pseudo-density component continues as rho(rd)*exp(-(r-rd)/L); L=5 kpc baseline; integrate to R=rd+12L then harmonic Robin match to u(infinity)=1",
            "outer_boundary": "u'(R)=-(u(R)-1)/R, no observed outer V used",
            "physical_source_limitation": "pseudo-spherical rho has negative intervals; local Rotmod does not uniquely determine a 3D source, inner core or exterior taper",
        },
        "equation": {
            "stationary_phi": "phi''+2phi'/r=4 pi G rho + K(phi')^2",
            "transform": "u=exp(-K phi)>0; u''+2u'/r+4 pi G K rho u=0",
            "stable_variable": "x=(1-w)/K; x''+2x'/r=4 pi G rho*(1-K*x), with exact K=0 Poisson limit x''+2x'/r=4 pi G rho",
            "model_velocity": "v_model^2=r*x'/(1-K*x) for K>0; v_model^2=r*x' at K=0; nonpositive values reject the model (no clipping)",
            "K_units": "s^2/km^2",
            "G_kpc_km2_s2_Msun": G,
        },
        "predeclared_domain": {
            "K_s2_per_km2": list(K_DOMAIN),
            "K_upper_over_paper_K": K_DOMAIN[1] / paper_k,
            "Upsilon_Msun_per_Lsun": list(UPS_RANGE),
            "Upsilon_grid": UPS_GRID,
            "central_velocity_powers": list(INNER_POWERS),
            "outer_taper_length_kpc_baseline": TAIL_LENGTH_KPC,
            "taper_tail_extent": TAIL_EXTENTS,
            "fit_statistic": "sum_i [(Vobs_i - Vmodel_i)/errV_i]^2 over all 43 points",
            "random_seed": None,
            "optimizer": "deterministic threshold scan, dense grid, bounded Powell refinement, bounded Newtonian-slice minimization",
            "stopping_rule": "stop K scan at first normalization or u-positivity failure connected to K=0; profile is finite and no larger than 1e-2",
        },
        "pseudodensity_diagnostics_at_paper_Upsilon": density_diag,
        "K0_limit_stability_at_paper_Upsilon": k0_tests,
        "first_branch_thresholds_by_Upsilon": thresholds,
        "profiled_minima": refined,
        "profile_surface_csv": csv_path.name,
        "profile_grid_point_count": len(profile_rows),
        "run": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "logical_cpu_count": os.cpu_count(),
            "numerical_cpu_seconds": cpu_s,
            "wall_seconds": wall_s,
            "ode_solve_count": N_SOLVES,
            "cpu_cap_seconds": CPU_LIMIT_SECONDS,
            "solver": "SciPy solve_ivp DOP853, rtol 2e-9 / atol 1e-9 / max_step 0.25 kpc; candidate check 2e-11 / 1e-11 / 0.06 kpc",
        },
        "code_sha256": sha(Path(__file__).read_bytes()),
    }
    out = HERE / "profile.json"
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"profile_json": str(out), "profile_csv": str(csv_path),
                      "cpu_seconds": cpu_s, "wall_seconds": wall_s,
                      "solve_count": N_SOLVES,
                      "selected": {p: v.get("selected_candidate") for p, v in refined.items()}}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except TimeoutError as exc:
        print(str(exc))
        raise

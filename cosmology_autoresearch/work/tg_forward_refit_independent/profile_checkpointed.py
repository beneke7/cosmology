#!/usr/bin/env python3
"""Checkpointed parallel profile with an independent strict parser and spherical ODE."""
from __future__ import annotations

import concurrent.futures as cf
import csv
import hashlib
import json
import multiprocessing as mp
import os
import platform
import resource
import time
import zipfile
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
INITIAL_SCRIPT_PATH = HERE / "profile.py"  # preserved first independent implementation
ARCHIVE = PROJECT / "context" / "data" / "Rotmod_LTG.zip"
PAPER = PROJECT / "context" / "user_provided" / "pszota_van_2024.pdf"
MEMBER = "NGC3198_rotmod.dat"
GROUP_DIR = HERE / "groups"
OUT_JSON = HERE / "profile_checkpointed.json"
OUT_CSV = HERE / "profile_checkpointed_surface.csv"
CHECKPOINT_LOG = HERE / "checkpoint.jsonl"

K_DOMAIN = (0.0, 1.0e-2)
UPS_RANGE = (0.2, 1.0)
UPS_GRID = (0.2, 0.5, 0.7, 0.762, 1.0)
INNER_POWERS = (1.0, 1.5, 2.0)
PAPER_K = 3.40e-5
PAPER_UPS = 0.762
BASE_TAPER_KPC = 5.0
TAIL_EXTENTS = 12.0
# A deterministic stratified profile over the connected positive branch.
K_FRACTIONS = (0.0, 0.001, 0.01, 0.03, 0.06, 0.10, 0.18, 0.28,
               0.40, 0.55, 0.70, 0.82, 0.92, 0.97, 0.99, 0.995, 0.999)
K_SCAN = np.geomspace(1.0e-11, K_DOMAIN[1], 12)
BISECTION_STEPS = 12
GROUP_CPU_CAP_S = 8.5
SOFT_COMPLETED_CPU_S = 135.0
INITIAL_SERIAL_CPU_S = 979.0  # process snapshot at the explicit stop, rounded to seconds
TOTAL_NUMERIC_CPU_CAP_S = 1200.0
PRIOR_BENCHMARK_CPU_S = 2.657  # measured five-solve speed tests, step-size check, and one tight-solve timing
SOLVER_VERSION = "v2-regular-source-isolated-robin-u-and-speed-events"
PARENT_CPU_START = time.process_time()
G = 4.30091e-6  # kpc (km/s)^2 / Msun


def parse_rotmod() -> dict[str, np.ndarray]:
    """Independent strict parser for the exact public SPARC Rotmod member."""
    with zipfile.ZipFile(ARCHIVE) as zf:
        raw = zf.read(MEMBER)
    rows = []
    for line_no, line in enumerate(raw.decode("ascii").splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 8:
            raise ValueError(f"{MEMBER}:{line_no}: expected 8 columns, found {len(fields)}")
        rows.append([float(item) for item in fields])
    a = np.asarray(rows, dtype=float)
    if a.shape != (43, 8) or np.any(np.diff(a[:, 0]) <= 0) or np.any(a[:, 2] <= 0):
        raise ValueError(f"invalid NGC3198 Rotmod table shape/order/errors: {a.shape}")
    return {"raw": raw, "r": a[:, 0], "vobs": a[:, 1], "sigma": a[:, 2],
            "vgas": a[:, 3], "vdisk": a[:, 4], "sbdisk": a[:, 6], "sbbul": a[:, 7]}


class FastSource:
    """Signed spherical inversion with cached PCHIP derivatives and declared tails."""

    def __init__(self, p_inner: float, tail_length: float = BASE_TAPER_KPC):
        self.p = float(p_inner)
        self.L = float(tail_length)
        self.r = DATA["r"]
        self.r0, self.rd = float(self.r[0]), float(self.r[-1])
        self.disk = PchipInterpolator(self.r, DATA["vdisk"], extrapolate=False)
        self.gas = PchipInterpolator(self.r, DATA["vgas"], extrapolate=False)
        self.disk_derivative = self.disk.derivative()
        self.gas_derivative = self.gas.derivative()
        self.vd0, self.vg0 = float(DATA["vdisk"][0]), float(DATA["vgas"][0])
        self.m_center = 2.0 * self.p - 2.0
        mult = 1.0 + 2.0 * self.p
        self.rho_d0 = mult * self.vd0**2 / (4.0 * np.pi * G * self.r0**2)
        self.rho_g0 = mult * self.vg0 * abs(self.vg0) / (4.0 * np.pi * G * self.r0**2)
        self.rho_d_edge, self.rho_g_edge = self._density_parts_at_edge()
        self.R = self.rd + TAIL_EXTENTS * self.L

    def _density_parts_at_edge(self):
        r = self.rd
        vd, vdp = float(self.disk(r)), float(self.disk_derivative(r))
        vg, vgp = float(self.gas(r)), float(self.gas_derivative(r))
        den = 4.0 * np.pi * G * r * r
        return ((vd * vd + 2.0 * r * vd * vdp) / den,
                (vg * abs(vg) + 2.0 * r * abs(vg) * vgp) / den)

    def rho_parts(self, radius: float):
        r = float(radius)
        if r < self.r0:
            scale = (max(r, 0.0) / self.r0) ** self.m_center
            return self.rho_d0 * scale, self.rho_g0 * scale
        if r <= self.rd:
            vd, vdp = float(self.disk(r)), float(self.disk_derivative(r))
            vg, vgp = float(self.gas(r)), float(self.gas_derivative(r))
            den = 4.0 * np.pi * G * r * r
            return ((vd * vd + 2.0 * r * vd * vdp) / den,
                    (vg * abs(vg) + 2.0 * r * abs(vg) * vgp) / den)
        factor = np.exp(-(r - self.rd) / self.L)
        return self.rho_d_edge * factor, self.rho_g_edge * factor

    def rho(self, radius: float, ups: float) -> float:
        d, g = self.rho_parts(radius)
        return float(ups * d + g)

    def center_coefficient(self, ups: float) -> float:
        return float(ups * self.rho_d0 + self.rho_g0)

    def integration_boundaries(self):
        eps = min(1.0e-7, self.r0 * 1.0e-6)
        return [eps, self.r0, self.rd, self.R]


DATA = parse_rotmod()
GROUP_DIR.mkdir(exist_ok=True)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def atomic_json(path: Path, obj: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def eval_segments(segments, radius: float) -> np.ndarray:
    rr = float(radius)
    for left, right, sol in segments:
        if left - 1e-12 <= rr <= right + 1e-12:
            return np.asarray(sol.sol(min(max(rr, left), right)), dtype=float)
    raise ValueError(f"radius {rr} outside integrated source domain")


def solve_fast(source, kappa: float, ups: float, *, rtol=2e-8,
               atol=1e-8, max_step=0.35) -> dict:
    """Independent fast x-ODE solve, with first-zero event and isolated Robin match."""
    if not (K_DOMAIN[0] <= kappa <= K_DOMAIN[1] and UPS_RANGE[0] <= ups <= UPS_RANGE[1]):
        return {"valid": False, "reason": "outside_declared_domain"}
    eps, r0, rd, R = source.integration_boundaries()
    m = source.m_center
    q0 = 4.0 * np.pi * G * source.center_coefficient(ups)
    x0 = q0 * eps**(m + 2.0) / (r0**m * (m + 3.0) * (m + 2.0))
    xp0 = q0 * eps**(m + 1.0) / (r0**m * (m + 3.0))
    y0 = np.array([x0, xp0], dtype=float)
    cuts = (eps, r0, rd, R)
    segments = []
    rhs_calls = 0

    def rhs(rr, y):
        nonlocal rhs_calls
        rhs_calls += 1
        rho = source.rho(rr, ups)
        return np.array([y[1], 4.0 * np.pi * G * rho * (1.0 - kappa * y[0]) - 2.0 * y[1] / rr])

    def speed_event(rr, y):
        return y[1]
    speed_event.terminal = True
    speed_event.direction = -1.0

    if kappa == 0.0:
        events = [speed_event]
    else:
        def u_event(rr, y):
            return 1.0 - kappa * y[0]
        u_event.terminal = True
        u_event.direction = -1.0
        events = [u_event, speed_event]

    state = y0
    for left, right in zip(cuts[:-1], cuts[1:]):
        sol = solve_ivp(rhs, (left, right), state, method="DOP853",
                        rtol=rtol, atol=atol, max_step=max_step,
                        dense_output=True, events=events)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return {"valid": False, "reason": "integration_failure", "rhs_calls": rhs_calls}
        if kappa > 0.0 and len(sol.t_events[0]):
            return {"valid": False, "reason": "u_positivity_threshold",
                    "first_zero_r_kpc": float(sol.t_events[0][0]), "rhs_calls": rhs_calls}
        if len(sol.t_events[-1]):
            return {"valid": False, "reason": "nonpositive_model_v2",
                    "first_zero_r_kpc": float(sol.t_events[-1][0]), "rhs_calls": rhs_calls}
        segments.append((left, right, sol))
        state = sol.y[:, -1]
    xR, xpR = state
    norm = 1.0 - kappa * (xR + R * xpR)
    if not np.isfinite(norm):
        return {"valid": False, "reason": "nonfinite_normalization", "rhs_calls": rhs_calls}
    if kappa > 0.0 and norm <= 0.0:
        return {"valid": False, "reason": "normalization_threshold",
                "normalization": float(norm), "rhs_calls": rhs_calls}
    # Accepted integrator states independently protect against an event missed between steps.
    if kappa > 0.0:
        wmin = min(float(np.min(1.0 - kappa * sol.y[0])) for _, _, sol in segments)
        if wmin <= 0.0:
            return {"valid": False, "reason": "u_positivity_threshold",
                    "normalization": float(norm), "min_w": wmin, "rhs_calls": rhs_calls}
    xpmin = min(float(np.min(sol.y[1])) for _, _, sol in segments)
    if xpmin <= 0.0:
        return {"valid": False, "reason": "nonpositive_model_v2",
                "min_model_v2": xpmin, "rhs_calls": rhs_calls}
    states = np.array([eval_segments(segments, r) for r in DATA["r"]])
    denom = 1.0 - kappa * states[:, 0]
    v2 = DATA["r"] * states[:, 1] / denom if kappa > 0.0 else DATA["r"] * states[:, 1]
    if np.any(~np.isfinite(v2)):
        return {"valid": False, "reason": "nonfinite_model_v2", "rhs_calls": rhs_calls}
    if np.any(v2 <= 0.0):
        return {"valid": False, "reason": "nonpositive_model_v2",
                "min_model_v2": float(np.min(v2)), "rhs_calls": rhs_calls}
    vm = np.sqrt(v2)
    chi2 = float(np.sum(((DATA["vobs"] - vm) / DATA["sigma"])**2))
    uR = float((1.0 - kappa * xR) / norm) if kappa > 0.0 else 1.0
    return {
        "valid": True, "chi2": chi2, "model_v_kms": vm, "model_v2_km2_s2": v2,
        "segments": segments, "normalization": float(norm),
        "min_w": float(wmin) if kappa > 0.0 else 1.0, "u_at_R": uR,
        "R_kpc": float(R), "rhs_calls": rhs_calls,
    }


def score_only(source, kappa, ups, *, tight=False) -> dict:
    if tight:
        res = solve_fast(source, kappa, ups, rtol=1e-11, atol=1e-11, max_step=0.06)
    else:
        res = solve_fast(source, kappa, ups)
    if res["valid"]:
        return {key: res[key] for key in ("valid", "chi2", "normalization", "min_w",
                                          "u_at_R", "R_kpc", "rhs_calls")}
    return {key: res.get(key) for key in ("valid", "reason", "normalization",
                                          "min_w", "first_zero_r_kpc", "min_model_v2",
                                          "rhs_calls")}


def group_file(p, ups) -> Path:
    return GROUP_DIR / f"p{str(p).replace('.', 'p')}_u{str(ups).replace('.', 'p')}.json"


def run_group(task: tuple[float, float]) -> dict:
    """One resumable central-continuation/M-L group; unique file per group."""
    p, ups = task
    start_cpu = time.process_time()
    src = FastSource(p, BASE_TAPER_KPC)
    path = group_file(p, ups)
    record = {
        "status": "running",
        "p_inner": p, "Upsilon": ups,
        "K_domain_s2_km2": list(K_DOMAIN),
        "Upsilon_domain": list(UPS_RANGE),
        "threshold_scan_K": K_SCAN.tolist(),
        "profile_fractions": list(K_FRACTIONS),
        "profile": [],
        "threshold": None,
        "source_taper_kpc": BASE_TAPER_KPC,
        "score": "direct sum of squared velocity residuals over all 43 Rotmod points",
        "boundary": "regular center; source tapered; Robin match at R=rd+12L; never uses observed endpoint velocity",
    }

    # Reuse a durable completed group, or resume a previously thresholded profile.
    existing = None
    if path.exists():
        try:
            existing = json.loads(path.read_text())
        except Exception:
            existing = None
    if existing and existing.get("status") == "complete":
        existing["worker_cpu_seconds_prior"] = float(existing.get("worker_cpu_seconds_total",
                                                    existing.get("worker_cpu_seconds", 0.0)))
        existing["worker_cpu_seconds_current"] = 0.0
        atomic_json(path, existing)
        return existing
    previous_cpu = float(existing.get("worker_cpu_seconds_total",
                                      existing.get("worker_cpu_seconds", 0.0))) if existing else 0.0
    if existing:
        record["profile"] = list(existing.get("profile", []))
        prior_threshold = existing.get("threshold")
        if prior_threshold and prior_threshold.get("status") in (
                "first_connected_branch_threshold_bracketed", "no_threshold_within_domain_cap"):
            record["threshold"] = prior_threshold

    def checkpoint(status: str):
        record["status"] = status
        elapsed = time.process_time() - start_cpu
        record["worker_cpu_seconds"] = elapsed
        record["worker_cpu_seconds_current"] = elapsed
        record["worker_cpu_seconds_prior"] = previous_cpu
        record["worker_cpu_seconds_total"] = previous_cpu + elapsed
        atomic_json(path, record)

    def group_budget():
        return time.process_time() - start_cpu >= GROUP_CPU_CAP_S

    # K=0 is the stable Newtonian slice and start of the connected physical branch.
    if record["threshold"] is None:
        zero = solve_fast(src, 0.0, ups)
        if not zero["valid"]:
            record["threshold"] = {"status": "Newtonian_slice_invalid", "details": score_only(src, 0.0, ups)}
            checkpoint("complete_invalid_newtonian")
            return record
        last = 0.0
        bracket = None
        threshold_kind = None
        for kval in K_SCAN:
            if group_budget():
                record["threshold"] = {"status": "partial_before_threshold", "last_valid_K": last}
                checkpoint("partial_group_cpu_cap")
                return record
            trial = solve_fast(src, float(kval), ups)
            if not trial["valid"]:
                if trial["reason"] == "integration_failure":
                    record["threshold"] = {"status": "numerical_failure_before_threshold",
                                           "last_valid_K": last, "first_failed_K": float(kval)}
                    checkpoint("partial_numerical_failure")
                    return record
                bracket = [last, float(kval)]
                threshold_kind = trial["reason"]
                break
            last = float(kval)
        if bracket is None:
            kcrit = K_DOMAIN[1]
            tstatus = "no_threshold_within_domain_cap"
            hi = None
        else:
            lo, hi = bracket
            for _ in range(BISECTION_STEPS):
                if group_budget():
                    record["threshold"] = {"status": "partial_threshold_bisection",
                                           "bracket_K": [lo, hi]}
                    checkpoint("partial_group_cpu_cap")
                    return record
                mid = np.sqrt(lo * hi) if lo > 0.0 else 0.5 * hi
                trial = solve_fast(src, float(mid), ups)
                if trial["valid"]:
                    lo = float(mid)
                elif trial["reason"] in ("normalization_threshold", "u_positivity_threshold", "nonpositive_model_v2"):
                    hi = float(mid)
                else:
                    record["threshold"] = {"status": "threshold_refinement_numerical_failure",
                                           "bracket_K": [lo, hi], "reason": trial["reason"]}
                    checkpoint("partial_numerical_failure")
                    return record
            kcrit = lo
            tstatus = "first_connected_branch_threshold_bracketed"
        record["threshold"] = {
            "status": tstatus, "K_just_below_s2_km2": float(kcrit),
            "first_invalid_K_s2_km2": float(hi) if bracket else None,
            "threshold_kind": threshold_kind,
        }
        checkpoint("threshold_complete")
    else:
        kcrit = float(record["threshold"]["K_just_below_s2_km2"])

    already_profiled = {float(row["fraction"]) for row in record["profile"]}
    for frac in K_FRACTIONS:
        if float(frac) in already_profiled:
            continue
        if group_budget():
            checkpoint("partial_group_cpu_cap")
            return record
        kval = float(frac * kcrit)
        res = score_only(src, kval, ups)
        row = {"fraction": float(frac), "K_s2_km2": kval, **res}
        record["profile"].append(row)
        checkpoint("profiling")
    valid_rows = [x for x in record["profile"] if x.get("valid")]
    record["best_group_point"] = min(valid_rows, key=lambda x: x["chi2"]) if valid_rows else None
    record["status"] = "complete"
    elapsed = time.process_time() - start_cpu
    record["worker_cpu_seconds"] = elapsed
    record["worker_cpu_seconds_current"] = elapsed
    record["worker_cpu_seconds_prior"] = previous_cpu
    record["worker_cpu_seconds_total"] = previous_cpu + elapsed
    atomic_json(path, record)
    return record


def riccati_check(source, kappa, ups, x_res) -> dict:
    eps, r0, rd, R = source.integration_boundaries()
    m = source.m_center
    q0 = 4.0 * np.pi * G * source.center_coefficient(ups)
    g0 = -q0 * eps**(m + 1.0) / (r0**m * (m + 3.0))
    y = np.array([g0])
    segments = []

    def rhs(rr, state):
        g = state[0]
        return np.array([-2.0 * g / rr - 4.0 * np.pi * G * source.rho(rr, ups) - kappa * g * g])

    for left, right in zip((eps, r0, rd), (r0, rd, R)):
        sol = solve_ivp(rhs, (left, right), y, method="DOP853",
                        rtol=1e-10, atol=1e-10, max_step=0.08, dense_output=True)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return {"status": "riccati_integration_failure", "message": sol.message}
        segments.append((left, right, sol))
        y = sol.y[:, -1]
    g_direct = np.array([eval_segments(segments, r)[0] for r in DATA["r"]])
    g_transformed = []
    for r in DATA["r"]:
        x, xp = eval_segments(x_res["segments"], r)
        g_transformed.append(-xp / (1.0 - kappa * x))
    g_transformed = np.asarray(g_transformed)
    diff = np.abs(g_direct - g_transformed)
    rel = diff / np.maximum(np.abs(g_transformed), 1e-14)
    gr = float(y[0])
    C = float(-1.0 / (gr * R * R) - kappa / R) if gr != 0.0 else float("inf")
    return {
        "status": "completed",
        "max_abs_acceleration_difference_km_s2": float(np.max(diff)),
        "max_relative_acceleration_difference": float(np.max(rel)),
        "g_at_R_km_s2": gr,
        "vacuum_C": C,
        "attractive_C_positive": bool(C > 0.0),
    }


def verify_profile_group(task: tuple[dict, str]) -> dict:
    """Recheck all profile points and the first physical threshold with speed/u events."""
    record, _unused = task
    started = time.process_time()
    p, ups = float(record["p_inner"]), float(record["Upsilon"])
    path = group_file(p, ups)
    if (record.get("verification", {}).get("status") == "complete" and
            record.get("verification", {}).get("solver_version") == SOLVER_VERSION and
            record.get("verification", {}).get("profile_rows") == len(K_FRACTIONS)):
        record["verification"]["worker_cpu_seconds_current"] = 0.0
        atomic_json(path, record)
        return record

    src = FastSource(p, BASE_TAPER_KPC)
    threshold = record.get("threshold") or {}
    # Start from K=0 and independently scan the full declared interval with both
    # positivity and circular-speed events enabled.
    kzero = solve_fast(src, 0.0, ups)
    if not kzero["valid"]:
        record["verification"] = {
            "status": "newtonian_slice_invalid", "solver_version": SOLVER_VERSION,
            "k0": score_only(src, 0.0, ups),
        }
        atomic_json(path, record)
        return record
    last, bracket, initial_reason = 0.0, None, None
    for kval in K_SCAN:
        trial = solve_fast(src, float(kval), ups)
        if trial["valid"]:
            last = float(kval)
            continue
        if trial["reason"] == "integration_failure":
            record["verification"] = {
                "status": "threshold_scan_numerical_failure", "solver_version": SOLVER_VERSION,
                "last_valid_K": last, "first_failed_K": float(kval), "reason": trial["reason"],
            }
            atomic_json(path, record)
            return record
        bracket, initial_reason = [last, float(kval)], trial["reason"]
        break

    if bracket is None:
        kcrit, final_hi, threshold_status, tight_reason = K_DOMAIN[1], None, "no_threshold_within_domain_cap", None
        tight_below = solve_fast(src, kcrit, ups, rtol=1e-11, atol=1e-11, max_step=0.06)
        tight_above = None
    else:
        lo, hi = bracket
        for _ in range(BISECTION_STEPS):
            mid = np.sqrt(lo * hi) if lo > 0.0 else 0.5 * hi
            trial = solve_fast(src, float(mid), ups)
            if trial["valid"]:
                lo = float(mid)
            elif trial["reason"] in ("normalization_threshold", "u_positivity_threshold", "nonpositive_model_v2"):
                hi = float(mid)
            else:
                record["verification"] = {
                    "status": "normal_threshold_refinement_failure", "solver_version": SOLVER_VERSION,
                    "bracket_K": [lo, hi], "reason": trial["reason"],
                }
                atomic_json(path, record)
                return record
        # Tighter bisection resolves whether the last normal-tolerance valid solution
        # was displaced by numerical error. Six logarithmic steps give ~1e-6 relative
        # bracketing after the normal-tolerance stage.
        for _ in range(6):
            mid = np.sqrt(lo * hi) if lo > 0.0 else 0.5 * hi
            trial = solve_fast(src, float(mid), ups, rtol=1e-11, atol=1e-11, max_step=0.06)
            if trial["valid"]:
                lo = float(mid)
            elif trial["reason"] in ("normalization_threshold", "u_positivity_threshold", "nonpositive_model_v2"):
                hi = float(mid)
            else:
                record["verification"] = {
                    "status": "tight_threshold_refinement_failure", "solver_version": SOLVER_VERSION,
                    "bracket_K": [lo, hi], "reason": trial["reason"],
                }
                atomic_json(path, record)
                return record
        kcrit, final_hi = lo, hi
        tight_below = solve_fast(src, kcrit, ups, rtol=1e-11, atol=1e-11, max_step=0.06)
        tight_above = solve_fast(src, final_hi, ups, rtol=1e-11, atol=1e-11, max_step=0.06)
        threshold_status = "first_connected_branch_threshold_bracketed"
        tight_reason = tight_above.get("reason")

    lower_ok = bool(tight_below and tight_below.get("valid"))
    upper_ok = (tight_above is None or not tight_above.get("valid"))
    record["threshold"] = {
        "status": threshold_status,
        "K_just_below_s2_km2": float(kcrit),
        "first_invalid_K_s2_km2": float(final_hi) if final_hi is not None else None,
        "tight_bracket_K_s2_km2": [float(kcrit), float(final_hi)] if final_hi is not None else None,
        "initial_12_node_scan_reason": initial_reason,
        "tight_upper_failure_reason": tight_reason,
        "tight_bisection_steps": 6,
        "tight_lower_valid": lower_ok,
        "tight_upper_invalid": upper_ok,
        "tight_bracket_confirmed": bool(lower_ok and upper_ok),
        "profile_threshold_before_refinement_K": threshold.get("K_just_below_s2_km2"),
        "threshold_kind": tight_reason or initial_reason,
    }

    new_rows = []
    max_delta = 0.0
    for fraction in K_FRACTIONS:
        kval = float(fraction * kcrit)
        checked = score_only(src, kval, ups)
        old = next((row for row in record.get("profile", [])
                    if abs(float(row["fraction"]) - float(fraction)) < 1e-14), {})
        if old.get("valid") and checked.get("valid"):
            max_delta = max(max_delta, abs(float(checked["chi2"]) - float(old["chi2"])))
        new_rows.append({"fraction": float(fraction), "K_s2_km2": kval, **checked})
    record["profile"] = new_rows
    valid_rows = [row for row in new_rows if row.get("valid")]
    record["best_group_point"] = min(valid_rows, key=lambda row: row["chi2"]) if valid_rows else None
    record["verification"] = {
        "status": "complete", "solver_version": SOLVER_VERSION,
        "profile_rows": len(new_rows), "profile_re_evaluation_max_abs_delta_chi2": float(max_delta),
        "threshold_tight_bracket_confirmed": bool(lower_ok and upper_ok),
        "threshold_tight_lower_valid": lower_ok, "threshold_tight_upper_invalid": upper_ok,
        "worker_cpu_seconds_current": time.process_time() - started,
    }
    record["threshold"]["scan_and_profile_solver_version"] = SOLVER_VERSION
    record["profile_solver_version"] = SOLVER_VERSION
    atomic_json(path, record)
    return record


def append_checkpoint(record):
    with CHECKPOINT_LOG.open("a") as f:
        f.write(json.dumps({
            "group": [record["p_inner"], record["Upsilon"]],
            "status": record["status"],
            "worker_cpu_seconds": record.get("worker_cpu_seconds"),
            "profile_rows": len(record.get("profile", [])),
            "group_file": group_file(record["p_inner"], record["Upsilon"]).name,
        }, allow_nan=False) + "\n")


def main():
    wall0 = time.perf_counter()
    # Priority ordering yields all three central continuations at the paper value first.
    tasks = [(p, PAPER_UPS) for p in INNER_POWERS]
    tasks += [(p, u) for u in (0.5, 0.7, 0.2, 1.0) for p in INNER_POWERS]
    records = []
    stop_submission = False
    completed_group_cpu = 0.0
    affinity = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    max_workers = min(len(tasks), max(1, affinity))
    ctx = mp.get_context("fork")
    child_cpu_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    child_cpu_before_s = child_cpu_before.ru_utime + child_cpu_before.ru_stime
    with cf.ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as pool:
        futures = {pool.submit(run_group, task): task for task in tasks}
        for fut in cf.as_completed(futures):
            try:
                rec = fut.result()
            except Exception as exc:
                task = futures[fut]
                rec = {"status": "worker_exception", "p_inner": task[0], "Upsilon": task[1],
                       "error": repr(exc), "worker_cpu_seconds": 0.0, "profile": []}
                atomic_json(group_file(*task), rec)
            records.append(rec)
            completed_group_cpu += float(rec.get("worker_cpu_seconds", 0.0))
            append_checkpoint(rec)
            if completed_group_cpu >= SOFT_COMPLETED_CPU_S and not stop_submission:
                stop_submission = True
                for pending in futures:
                    if not pending.done() and not pending.running():
                        pending.cancel()
                # Each submitted group has its own 8.5 CPU-s checkpoint guard.
    # Robustly read every durable group artifact, including any finished just before cancellation.
    all_records = []
    for path in sorted(GROUP_DIR.glob("*.json")):
        try:
            all_records.append(json.loads(path.read_text()))
        except Exception:
            continue

    # Independently re-evaluate the surface with a circular-speed-zero event as
    # well as u-positivity, then tighten each connected-branch threshold.
    verification_tasks = [(rec, "") for rec in all_records if rec.get("status") == "complete"]
    verification_records = []
    with cf.ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as pool:
        futures = {pool.submit(verify_profile_group, task): task for task in verification_tasks}
        for future in cf.as_completed(futures):
            try:
                verification_records.append(future.result())
            except Exception as exc:
                failed_record = dict(futures[future][0])
                failed_record["verification"] = {"status": "verification_worker_exception", "error": repr(exc)}
                verification_records.append(failed_record)
    verified_keys = {(r.get("p_inner"), r.get("Upsilon")) for r in verification_records}
    untouched = [r for r in all_records if (r.get("p_inner"), r.get("Upsilon")) not in verified_keys]
    all_records = verification_records + untouched

    grid_rows = []
    for rec in all_records:
        for row in rec.get("profile", []):
            grid_rows.append({
                "p_inner": rec.get("p_inner"), "Upsilon": rec.get("Upsilon"),
                "K_s2_km2": row.get("K_s2_km2"), "fraction": row.get("fraction"),
                "valid": row.get("valid"), "chi2": row.get("chi2"), "reason": row.get("reason"),
                "normalization": row.get("normalization"), "min_w": row.get("min_w"),
                "group_status": rec.get("status"),
            })
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(grid_rows[0].keys()) if grid_rows else
                                ["p_inner", "Upsilon", "K_s2_km2", "fraction", "valid", "chi2", "reason", "normalization", "min_w", "group_status"])
        writer.writeheader()
        writer.writerows(grid_rows)

    # Summarize minima separately for each continuation; verify only best candidate with tight solves.
    minima = {}
    tight_cpu_start = time.process_time()
    for p in INNER_POWERS:
        candidates = [r for r in grid_rows if r["p_inner"] == p and r.get("valid") is True]
        if not candidates:
            minima[str(p)] = {"status": "no_valid_profile_point"}
            continue
        candidate = min(candidates, key=lambda r: r["chi2"])
        source = FastSource(p, BASE_TAPER_KPC)
        tight = solve_fast(source, candidate["K_s2_km2"], candidate["Upsilon"],
                           rtol=1e-11, atol=1e-11, max_step=0.06)
        verification = {
            "valid": bool(tight["valid"]),
            "chi2_tight": float(tight["chi2"]) if tight["valid"] else None,
            "delta_chi2_tight_minus_profile": float(tight["chi2"] - candidate["chi2"]) if tight["valid"] else None,
            "normalization_tight": tight.get("normalization"),
            "min_w_tight": tight.get("min_w"),
        }
        if tight["valid"] and candidate["K_s2_km2"] > 0.0:
            ric = riccati_check(source, candidate["K_s2_km2"], candidate["Upsilon"], tight)
        else:
            ric = {"status": "K=0 branch or invalid tight solution; direct Riccati comparison omitted"}
        source10 = FastSource(p, 10.0)
        tail = solve_fast(source10, candidate["K_s2_km2"], candidate["Upsilon"],
                          rtol=1e-11, atol=1e-11, max_step=0.08)
        tail_sensitivity = {
            "L_kpc": 10.0, "valid": bool(tail["valid"]),
            "chi2_same_K_Upsilon": float(tail["chi2"]) if tail["valid"] else None,
            "delta_chi2_from_L5": float(tail["chi2"] - tight["chi2"]) if tail["valid"] and tight["valid"] else None,
            "reason": tail.get("reason"),
        }
        minima[str(p)] = {
            "status": "grid minimum over completed groups, not claimed continuous global optimum",
            "selected_grid_point": candidate,
            "tight_tolerance_verification": verification,
            "tight_model_v_kms": [float(v) for v in tight["model_v_kms"]] if tight["valid"] else None,
            "independent_original_Riccati_check": ric,
            "fixed_parameter_outer_taper_sensitivity": tail_sensitivity,
        }

    # Exact K=0 recovery and the small-K limit: direct Newtonian baryon curve must be recovered.
    k0_checks = {}
    for p in INNER_POWERS:
        source = FastSource(p, BASE_TAPER_KPC)
        values = {}
        for k in (0.0, 1e-12, 1e-11, 1e-10, 1e-8, PAPER_K):
            res = solve_fast(source, k, PAPER_UPS, rtol=1e-11, atol=1e-11, max_step=0.06)
            if not res["valid"]:
                values[str(k)] = {"valid": False, "reason": res.get("reason")}
                continue
            if k == 0.0:
                target = PAPER_UPS * DATA["vdisk"]**2 + DATA["vgas"] * np.abs(DATA["vgas"])
                err = res["model_v2_km2_s2"] - target
                values[str(k)] = {
                    "valid": True,
                    "max_abs_baryon_v2_recovery_error_km2_s2": float(np.max(np.abs(err))),
                    "max_fractional_baryon_v2_recovery_error": float(np.max(np.abs(err) / np.maximum(np.abs(target), 1.0))),
                    "chi2": float(res["chi2"]),
                }
            else:
                values[str(k)] = {"valid": True, "chi2": float(res["chi2"]),
                                  "delta_chi2_from_K0": float(res["chi2"] - values["0.0"]["chi2"])}
        k0_checks[str(p)] = values

    threshold_summary = {}
    for p in INNER_POWERS:
        vals = [x for x in all_records if x.get("p_inner") == p and x.get("threshold")]
        threshold_summary[str(p)] = {
            str(r["Upsilon"]): r["threshold"] for r in vals
        }
    worker_cpu = sum(float(r.get("worker_cpu_seconds_total", r.get("worker_cpu_seconds", 0.0))) for r in all_records)
    prior_group_cpu = sum(float(r.get("worker_cpu_seconds_prior", 0.0)) for r in all_records)
    verification_worker_cpu = sum(float(r.get("verification", {}).get("worker_cpu_seconds_current", 0.0))
                                  for r in all_records)
    # Dense sign audit of the raw signed pseudo-spherical inversion, without clipping.
    sign_r = np.linspace(float(DATA["r"][0]), float(DATA["r"][-1]), 5001)
    sign_source = FastSource(1.0)
    sign_vd, sign_vdp = sign_source.disk(sign_r), sign_source.disk_derivative(sign_r)
    sign_vg, sign_vgp = sign_source.gas(sign_r), sign_source.gas_derivative(sign_r)
    sign_den = 4.0 * np.pi * G * sign_r * sign_r
    rho_d = (sign_vd**2 + 2.0 * sign_r * sign_vd * sign_vdp) / sign_den
    rho_g = (sign_vg * np.abs(sign_vg) + 2.0 * sign_r * np.abs(sign_vg) * sign_vgp) / sign_den
    density_sign = {}
    for ups in UPS_GRID:
        rho = ups * rho_d + rho_g
        negative = rho < 0.0
        edges = np.flatnonzero(np.diff(np.r_[False, negative, False].astype(np.int8)))
        intervals = []
        for start, stop in zip(edges[::2], edges[1::2]):
            left = float(sign_r[start])
            right = float(sign_r[stop - 1])
            if start > 0:
                a, b = rho[start - 1], rho[start]
                left = float(sign_r[start - 1] + (sign_r[start] - sign_r[start - 1]) * (-a) / (b - a))
            if stop < len(sign_r):
                a, b = rho[stop - 1], rho[stop]
                right = float(sign_r[stop - 1] + (sign_r[stop] - sign_r[stop - 1]) * (-a) / (b - a))
            intervals.append([left, right])
        density_sign[str(ups)] = {
            "negative_intervals_kpc_linear_crossing_estimate": intervals,
            "minimum_rho_Msun_kpc3": float(np.min(rho)),
            "rho_at_first_and_last_data_Msun_kpc3": [float(rho[0]), float(rho[-1])],
            "grid_points": int(len(sign_r)),
        }

    parent_cpu = time.process_time() - PARENT_CPU_START
    child_cpu_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    child_pool_cpu = child_cpu_after.ru_utime + child_cpu_after.ru_stime - child_cpu_before_s
    cpu_accounting_guard = 2.0
    total_campaign_cpu = (INITIAL_SERIAL_CPU_S + PRIOR_BENCHMARK_CPU_S + prior_group_cpu +
                          child_pool_cpu + parent_cpu + cpu_accounting_guard)

    aggregate = {
        "status": "checkpointed independent profile; report scope is exactly the completed grid rows",
        "source_contract": {
            "parser": "independent strict line parser in checkpointed v2; verifies 43 rows/eight columns",
            "member_sha256": sha(DATA["raw"]),
            "paper_pdf_sha256": sha(PAPER.read_bytes()),
            "archive_sha256": sha(ARCHIVE.read_bytes()),
            "initial_v1_script_sha256": sha(INITIAL_SCRIPT_PATH.read_bytes()),
            "all_43_points_used": True,
            "signed_gas": True,
            "negative_density_clipped": False,
            "pseudo_spherical_source_limitation": True,
            "density_sign_audit": "signed analytic PCHIP derivative, 5001 radii over observed Rotmod support; negative values retained",
        },
        "equation": {
            "x_transform": "x=(1-w)/K gives x''+2x'/r=4piG rho(1-Kx); exact K=0 Poisson limit",
            "velocity": "v^2=r x'/(1-Kx), K>0; v^2=r x' at K=0; any nonpositive v^2 invalid",
            "outer_match": "u'(R)=-(u(R)-1)/R at R=rd+12L; no observed V(rd) used",
            "physical_branch": "event-checks u=1-Kx>0 and x'>0 throughout integration; positive exterior normalization; v^2 positive at all observed radii",
            "inner": "three regular continuations V_i(r)=V_i(r0)*(r/r0)^p for r<r0; p=1,1.5,2",
            "outer_taper": "signed pseudo-density terms exponentially tapered with L=5 kpc, fixed-parameter sensitivity L=10 kpc",
        },
        "predeclared_domain": {
            "K_s2_per_km2": list(K_DOMAIN), "K_upper_over_paper_value": K_DOMAIN[1] / PAPER_K,
            "Upsilon": list(UPS_RANGE), "Upsilon_grid": list(UPS_GRID),
            "inner_powers": list(INNER_POWERS), "K_fractions_per_branch": list(K_FRACTIONS),
            "K_scan_nodes_per_group": int(len(K_SCAN)), "bisection_steps": BISECTION_STEPS,
            "profile_method": "grid profile of direct all-point velocity chi-square, one threshold/profile group per (inner power,Upsilon)",
            "random_seed": None,
            "threshold_rule": "first loss of regular positive u, exterior normalization, or positive circular-speed solution connected to K=0",
        },
        "thresholds_by_completed_group": threshold_summary,
        "pseudo_density_sign_by_Upsilon": density_sign,
        "profiled_grid_minima": minima,
        "observed_rotmod_data": {
            "r_kpc": [float(v) for v in DATA["r"]],
            "vobs_kms": [float(v) for v in DATA["vobs"]],
            "sigma_kms": [float(v) for v in DATA["sigma"]],
            "vgas_kms_signed": [float(v) for v in DATA["vgas"]],
            "vdisk_kms": [float(v) for v in DATA["vdisk"]],
        },
        "K0_limit_checks": k0_checks,
        "completed_group_count": len(all_records),
        "completed_groups": [{"p_inner": x.get("p_inner"), "Upsilon": x.get("Upsilon"),
                              "status": x.get("status"), "profile_rows": len(x.get("profile", [])),
                              "worker_cpu_seconds": x.get("worker_cpu_seconds")} for x in all_records],
        "profile_surface_csv": OUT_CSV.name,
        "checkpoint_log": CHECKPOINT_LOG.name,
        "run_accounting": {
            "initial_serial_cpu_seconds_snapshot": INITIAL_SERIAL_CPU_S,
            "prior_solver_benchmark_validation_cpu_seconds": PRIOR_BENCHMARK_CPU_S,
            "checkpointed_child_group_cpu_seconds": worker_cpu,
            "checkpointed_prior_attempt_group_cpu_seconds": prior_group_cpu,
            "current_profile_verification_worker_cpu_seconds": verification_worker_cpu,
            "current_process_pool_user_system_cpu_seconds": child_pool_cpu,
            "orchestrator_cpu_seconds": parent_cpu,
            "conservative_cpu_accounting_guard_seconds": cpu_accounting_guard,
            "aggregate_numerical_cpu_seconds": total_campaign_cpu,
            "cap_seconds": TOTAL_NUMERIC_CPU_CAP_S,
            "cap_satisfied": bool(total_campaign_cpu <= TOTAL_NUMERIC_CPU_CAP_S),
            "wall_seconds_checkpointed_attempt": time.perf_counter() - wall0,
            "workers": max_workers,
            "verification_workers": max_workers,
            "cpu_affinity_count": affinity,
            "blas_thread_caps": {key: os.environ.get(key) for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
            "per_group_cpu_guard_seconds": GROUP_CPU_CAP_S,
            "soft_group_cpu_stop_seconds": SOFT_COMPLETED_CPU_S,
            "initial_serial_run_stopped_after_no-artifacts": True,
            "initial_serial_run_cpu_seconds_snapshot_source": "parent process monitor; 16:19 CPU time, rounded to 979 s",
        },
        "runtime": {
            "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "logical_cpu_count": os.cpu_count(), "solver": "DOP853",
            "normal_rtol": 2e-8, "normal_atol": 1e-8, "normal_max_step_kpc": 0.35,
            "candidate_rtol": 1e-11, "candidate_atol": 1e-11, "candidate_max_step_kpc": 0.06,
            "profile_validation": "recomputed every K point with an x'=0 event and recomputed thresholds with a tight bracket",
            "benchmark_initial_v1_five_solves_cpu_seconds": 0.917659,
            "benchmark_cached_v2_five_solves_cpu_seconds": 0.503968,
            "benchmark_interpretation": "same 1e-5, Upsilon=0.762 solve; cached derivative implementation is ~1.82x faster in this microbenchmark",
        },
        "code_sha256": sha(Path(__file__).read_bytes()),
    }
    atomic_json(OUT_JSON, aggregate)
    print(json.dumps({"json": str(OUT_JSON), "csv": str(OUT_CSV),
                      "completed_groups": len(all_records),
                      "profile_rows": len(grid_rows), "worker_cpu_seconds": worker_cpu,
                      "aggregate_numerical_cpu_seconds": total_campaign_cpu,
                      "minima": minima}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent score reproduction for the Dovekie IVS profile screen.

This checker intentionally does not import the production profile, its SN
scorer, or its ODE/distance helpers. It parses the pinned SNANA table itself,
unpacks the released full STAT+SYS precision itself, integrates the contract's
source-sign ODE with SciPy RK45, and evaluates the one-offset GLS profile.
Run from the cosmology_autoresearch project root after the profile result is
finalized.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad, solve_ivp


ROOT = Path(__file__).resolve().parents[2]
HD_PATH = ROOT / "context/data/des_dovekie_hd.csv"
P_PATH = ROOT / "context/data/des_dovekie_stat_sys.npz"
CONTRACT_JSON = ROOT / "work/theory7/dovekie_ivs_profile_contract.json"
CONTRACT_MD = ROOT / "work/theory7/dovekie_ivs_profile_contract.md"
RESULT_PATH = ROOT / "experiments/dovekie_ivs_profile/result.json"
PROFILE_PATH = ROOT / "experiments/dovekie_ivs_profile/profile.py"
PROFILE_HELPER = ROOT / "work/inference/fit_dovekie.py"
BAO_RESULT = ROOT / "experiments/interacting_vacuum_screen/result.json"
BAO_CODE = ROOT / "experiments/interacting_vacuum_screen/interacting_vacuum_profile.py"
OUT_JSON = ROOT / "work/compute10/dovekie_ivs_profile_independent_review.json"
OUT_MD = ROOT / "work/compute10/dovekie_ivs_profile_independent_review.md"

EXPECTED_HD_SHA = "2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf"
EXPECTED_P_SHA = "ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b"
EXPECTED_ROWS = 1820
H0_GAUGE = 70.0
C_KM_S = 299792.458
SCORE_TOL = 1.0e-4
DISTANCE_MODULUS_TOL = 2.0e-8


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_snana(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read VARNAMES/SN whitespace rows verbatim and preserve release order."""
    columns: list[str] | None = None
    records: list[dict[str, str]] = []
    for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "VARNAMES:":
            if columns is not None:
                raise ValueError(f"repeated VARNAMES at line {line_no}")
            columns = fields[1:]
        elif fields[0] == "SN:":
            if columns is None or len(fields[1:]) != len(columns):
                raise ValueError(f"bad SN row at line {line_no}")
            records.append(dict(zip(columns, fields[1:], strict=True)))
        else:
            raise ValueError(f"unexpected row at line {line_no}: {fields[0]}")
    if columns is None:
        raise ValueError("missing VARNAMES header")
    return columns, records


def load_release() -> dict:
    hd_sha, p_sha = file_sha(HD_PATH), file_sha(P_PATH)
    if hd_sha != EXPECTED_HD_SHA or p_sha != EXPECTED_P_SHA:
        raise ValueError(f"input hash mismatch: HD={hd_sha}, precision={p_sha}")
    columns, rows = read_snana(HD_PATH)
    needed = {"CID", "zHD", "zHEL", "MU"}
    if not needed.issubset(columns):
        raise ValueError(f"missing columns: {sorted(needed - set(columns))}")
    if len(rows) != EXPECTED_ROWS:
        raise ValueError(f"expected {EXPECTED_ROWS} rows, found {len(rows)}")
    if len({r["CID"] for r in rows}) != EXPECTED_ROWS:
        raise ValueError("CID is not unique")

    selected = [r for r in rows if float(r["zHD"]) > 0.0]
    if len(selected) != EXPECTED_ROWS:
        raise ValueError("official zHD>0 selection did not retain the whole release")
    z_hd = np.array([float(r["zHD"]) for r in selected], dtype=np.float64)
    z_hel = np.array([float(r["zHEL"]) for r in selected], dtype=np.float64)
    mu = np.array([float(r["MU"]) for r in selected], dtype=np.float64)

    with np.load(P_PATH, allow_pickle=False) as arc:
        n = int(np.asarray(arc["nsn"]).reshape(-1)[0])
        packed = np.asarray(arc["cov"])
        archive_keys = list(arc.files)
    if n != EXPECTED_ROWS or packed.size != n * (n + 1) // 2 or packed.dtype != np.float32:
        raise ValueError("packed precision format/dimensions differ from the pinned contract")
    precision = np.zeros((n, n), dtype=np.float64)
    upper = np.triu_indices(n)
    precision[upper] = packed.astype(np.float64)
    precision[(upper[1], upper[0])] = precision[upper]
    if not np.array_equal(precision, precision.T) or not np.all(np.isfinite(precision)):
        raise ValueError("unpacked precision is nonfinite or asymmetric")
    chol = np.linalg.cholesky(precision)
    one = np.ones(n, dtype=np.float64)
    p_one = precision @ one
    q = float(one @ p_one)
    if not (np.isfinite(q) and q > 0.0):
        raise ValueError("constant intercept direction is not positively constrained")
    cid_hash = hashlib.sha256(("\n".join(r["CID"] for r in rows) + "\n").encode()).hexdigest()
    return {
        "rows": rows,
        "z_hd": z_hd,
        "z_hel": z_hel,
        "mu": mu,
        "precision": precision,
        "chol": chol,
        "p_one": p_one,
        "q": q,
        "cid_hash": cid_hash,
        "archive_keys": archive_keys,
        "precision_packed_dtype": str(packed.dtype),
    }


def rhs(u: float, state: np.ndarray, omega_m: float, g: float) -> np.ndarray:
    """Contract equations in u=ln(1+z), with g=Gamma/H0."""
    matter, vacuum, _distance = state
    e2 = matter + vacuum
    if not (np.isfinite(e2) and e2 > 0.0):
        raise FloatingPointError("E^2 became nonpositive or nonfinite")
    e = math.sqrt(e2)
    transfer_over_h = g / e  # (Gamma/H0)/(H/H0) = Gamma/H
    source = transfer_over_h * vacuum
    return np.array((3.0 * matter + source, -source, math.exp(u) / e), dtype=np.float64)


def integrate(omega_m: float, g: float, z_eval: np.ndarray, *, tight: bool = False) -> dict:
    z_eval = np.asarray(z_eval, dtype=np.float64)
    if not (0.0 < omega_m < 1.0 and np.all(np.isfinite(z_eval)) and np.all(z_eval >= 0.0)):
        raise ValueError("invalid profile coordinate or redshift")
    umax = math.log1p(float(np.max(z_eval)))
    options = {"rtol": 2e-13, "atol": 2e-15, "max_step": 0.001} if tight else {
        "rtol": 2e-11, "atol": 2e-13, "max_step": 0.003
    }
    sol = solve_ivp(
        lambda u, y: rhs(u, y, omega_m, g),
        (0.0, umax),
        (omega_m, 1.0 - omega_m, 0.0),
        method="RK45",
        dense_output=True,
        **options,
    )
    if not sol.success or sol.sol is None:
        raise RuntimeError(f"independent RK45 integration failed: {sol.message}")
    u_grid = np.linspace(0.0, umax, 4001, dtype=np.float64)
    trajectory = np.asarray(sol.sol(u_grid), dtype=np.float64)
    m, x, d = trajectory
    e2 = m + x
    if np.any(~np.isfinite(trajectory)) or np.any(m <= 0) or np.any(x <= 0) or np.any(e2 <= 0):
        raise ValueError("m, x or E^2 failed the 4001-point positivity guard")
    fb_max_exact = 1.0 if g >= 0.0 else float(m[-1] / (omega_m * math.exp(3.0 * umax)))
    fb_witness = 0.5 * fb_max_exact
    baryons = fb_witness * omega_m * np.exp(3.0 * u_grid)
    cdm = m - baryons
    if not (0.0 < fb_witness < fb_max_exact) or np.any(baryons <= 0) or np.any(cdm <= 0):
        raise ValueError("no strictly positive baryon/CDM split witness on the guard grid")
    at_z = np.asarray(sol.sol(np.log1p(z_eval)), dtype=np.float64)
    return {
        "sol": sol.sol,
        "distance": at_z[2],
        "trajectory": trajectory,
        "u_grid": u_grid,
        "minimum_matter": float(np.min(m)),
        "minimum_vacuum": float(np.min(x)),
        "minimum_e2": float(np.min(e2)),
        "minimum_cdm_witness": float(np.min(cdm)),
        "fb_max_exact": fb_max_exact,
        "fb_witness": fb_witness,
        "u_max": umax,
    }


def flat_lcdm_distance(omega_m: float, z: float) -> float:
    value, _error = quad(
        lambda zz: 1.0 / math.sqrt(omega_m * (1.0 + zz) ** 3 + 1.0 - omega_m),
        0.0,
        float(z),
        epsabs=2e-13,
        epsrel=2e-13,
        limit=200,
    )
    return float(value)


def profile_at(data: dict, omega_m: float, g: float, *, tight: bool = False) -> dict:
    bg = integrate(omega_m, g, data["z_hd"], tight=tight)
    dl_mpc = (1.0 + data["z_hel"]) * (C_KM_S / H0_GAUGE) * bg["distance"]
    if np.any(~np.isfinite(dl_mpc)) or np.any(dl_mpc <= 0.0):
        raise ValueError("luminosity distance is not finite and positive")
    mu0 = 5.0 * np.log10(dl_mpc) + 25.0
    delta = mu0 - data["mu"]
    p_delta = data["precision"] @ delta
    offset = -float(data["p_one"] @ delta) / data["q"]
    residual = delta + offset
    p_residual = data["precision"] @ residual
    chi_direct = float(residual @ p_residual)
    whitened = data["chol"].T @ residual
    chi_whitened = float(whitened @ whitened)
    chi_schur = float(delta @ p_delta - (data["p_one"] @ delta) ** 2 / data["q"])
    normal = abs(float(data["p_one"] @ residual))
    return {
        "omega_m": float(omega_m),
        "g": float(g),
        "chi2": chi_direct,
        "offset": offset,
        "normal_equation_abs_residual": normal,
        "whitened_chi2": chi_whitened,
        "schur_chi2": chi_schur,
        "max_profile_algebra_disagreement": max(abs(chi_direct - chi_whitened), abs(chi_direct - chi_schur)),
        "background": bg,
    }


def compare_score(name: str, own: dict, reference: dict, result_delta: float | None = None) -> dict:
    root_chi = float(reference["chi2"])
    item = {
        "name": name,
        "parameters": {"Omega_m0": own["omega_m"], "g": own["g"]},
        "independent_chi2": own["chi2"],
        "reported_chi2": root_chi,
        "absolute_score_difference": abs(own["chi2"] - root_chi),
        "independent_profile_offset_mag": own["offset"],
        "normal_equation_abs_residual": own["normal_equation_abs_residual"],
        "max_gls_algebra_disagreement": own["max_profile_algebra_disagreement"],
    }
    if result_delta is not None:
        item["reported_delta_from_ivs_minimum"] = float(result_delta)
    return item


def main() -> None:
    start = time.perf_counter()
    if not RESULT_PATH.exists():
        raise FileNotFoundError(f"finalized profile result not found: {RESULT_PATH}")
    root_result_sha_before = file_sha(RESULT_PATH)
    result = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    if result.get("status") != "complete_exploratory_dovekie_only_ivs_profile":
        raise ValueError(f"profile result is not complete: {result.get('status')!r}")
    contract = json.loads(CONTRACT_JSON.read_text(encoding="utf-8"))
    data = load_release()

    reported_best = result["results"]["ivs_sn_only_minimum"]
    best_parameters = reported_best["parameters"]
    points = {
        "ivs_minimum": (float(best_parameters["Omega_m0"]), float(best_parameters["g"])),
    }
    # Avoid dict ordering assumptions when taking the fixed contract coordinate.
    points["bao_ivs_reference"] = (
        float(contract["bao_point_comparison"]["bao_ivs_shape"]["Omega_m0"]),
        float(contract["bao_point_comparison"]["bao_ivs_shape"]["g"]),
    )
    result_bao_parameters = result["results"]["bao_selected_ivs_reference_on_sn_profile"]["parameters"]
    if not np.allclose(points["bao_ivs_reference"],
                       (float(result_bao_parameters["Omega_m0"]), float(result_bao_parameters["g"])),
                       rtol=0.0, atol=1e-14):
        raise ValueError("profile's BAO reference coordinate differs from the frozen contract")

    lcdm_profile_ref = result["results"]["dovekie_lcdm_reference"]
    lcdm_profile_parameters = lcdm_profile_ref["parameters"]
    points["sn_lcdm_reference"] = (
        float(lcdm_profile_parameters["Omega_m0"]), float(lcdm_profile_parameters["g"])
    )
    points["bao_lcdm_reference"] = (
        float(contract["bao_point_comparison"]["bao_lcdm_shape"]["Omega_m0"]),
        float(contract["bao_point_comparison"]["bao_lcdm_shape"]["g"]),
    )

    root_refs = {
        "ivs_minimum": {"chi2": float(reported_best["profile"]["chi2"])},
        "bao_ivs_reference": {
            "chi2": float(result["results"]["bao_selected_ivs_reference_on_sn_profile"]["profile"]["chi2"])
        },
        "sn_lcdm_reference": {"chi2": float(lcdm_profile_ref["chi2"])},
    }
    own_scores = {name: profile_at(data, *coordinate) for name, coordinate in points.items()}
    comparisons = [compare_score("ivs_minimum", own_scores["ivs_minimum"], root_refs["ivs_minimum"])]
    comparisons.append(compare_score(
        "bao_ivs_reference", own_scores["bao_ivs_reference"], root_refs["bao_ivs_reference"],
        result["results"]["delta_chi2_bao_point_minus_sn_minimum"],
    ))
    comparisons.append(compare_score(
        "sn_lcdm_reference", own_scores["sn_lcdm_reference"], root_refs["sn_lcdm_reference"],
        result["results"]["delta_chi2_lcdm_minimum_minus_ivs_minimum"],
    ))
    comparisons.append(compare_score(
        "bao_lcdm_reference", own_scores["bao_lcdm_reference"],
        {"chi2": float("nan")},
    ))
    best_chi = own_scores["ivs_minimum"]["chi2"]
    for row in comparisons:
        row["independent_delta_from_ivs_minimum"] = row["independent_chi2"] - best_chi
    comparisons[-1]["reported_chi2"] = None
    comparisons[-1]["absolute_score_difference"] = None
    comparisons[-1].pop("reported_delta_from_ivs_minimum", None)

    # The explicit LCDM limit uses an analytic E(z) and adaptive scalar
    # quadrature, independent of the coupled IVS ODE solver.
    test_om = 0.3
    lcdm_limit_bg = integrate(test_om, 0.0, data["z_hd"])
    lcdm_limit_distance = np.array([flat_lcdm_distance(test_om, float(z)) for z in data["z_hd"]])
    lcdm_limit_mu = 5.0 * np.log10(
        (1.0 + data["z_hel"]) * (C_KM_S / H0_GAUGE) * lcdm_limit_distance
    ) + 25.0
    ivs_g0_mu = 5.0 * np.log10(
        (1.0 + data["z_hel"]) * (C_KM_S / H0_GAUGE) * lcdm_limit_bg["distance"]
    ) + 25.0
    g0_max_mu_diff = float(np.max(np.abs(ivs_g0_mu - lcdm_limit_mu)))

    # Cross-tolerance reproduction of every scored point using tighter RK45.
    convergence = []
    for name in points:
        standard = own_scores[name]
        tight = profile_at(data, *points[name], tight=True)
        convergence.append({
            "name": name,
            "standard_rk45_chi2": standard["chi2"],
            "tight_rk45_chi2": tight["chi2"],
            "absolute_chi2_difference": abs(standard["chi2"] - tight["chi2"]),
            "max_distance_modulus_difference_mag": float(np.max(np.abs(
                (5.0 * np.log10((1.0 + data["z_hel"]) * (C_KM_S / H0_GAUGE) * standard["background"]["distance"]) + 25.0)
                - (5.0 * np.log10((1.0 + data["z_hel"]) * (C_KM_S / H0_GAUGE) * tight["background"]["distance"]) + 25.0)
            ))),
        })

    for item in comparisons[:3]:
        if item["absolute_score_difference"] > SCORE_TOL:
            raise ArithmeticError(f"score mismatch for {item['name']}: {item['absolute_score_difference']}")
    if g0_max_mu_diff > DISTANCE_MODULUS_TOL:
        raise ArithmeticError(f"g=0 failed LCDM distance limit by {g0_max_mu_diff} mag")
    if max(r["absolute_chi2_difference"] for r in convergence) > SCORE_TOL:
        raise ArithmeticError("independent RK45 tolerance-refinement score check failed")
    if max(r["max_distance_modulus_difference_mag"] for r in convergence) > DISTANCE_MODULUS_TOL:
        raise ArithmeticError("independent RK45 tolerance-refinement distance check failed")

    # Check the implemented source direction numerically at u=0 for one point
    # with nonzero g; the source signs are fixed by Q=Gamma*rho_x>0 (CDM -> DE).
    sign_omega, sign_g = points["bao_ivs_reference"]
    sign_derivative = rhs(0.0, np.array([sign_omega, 1.0 - sign_omega, 0.0]), sign_omega, sign_g)
    sign_check = {
        "reference_g": sign_g,
        "dimensionless_rate_definition": "g=Gamma/H0; Gamma/H=g/E",
        "dx_du_at_today": float(sign_derivative[1]),
        "dm_du_minus_3m_at_today": float(sign_derivative[0] - 3.0 * sign_omega),
        "expected_x_source_sign": "opposite sign of g for x>0",
        "expected_m_source_sign": "same sign as g for x>0",
        "passed": bool(sign_derivative[1] * sign_g < 0 and
                        (sign_derivative[0] - 3.0 * sign_omega) * sign_g > 0),
    }
    if not sign_check["passed"]:
        raise ArithmeticError("source-sign check failed")

    point_physicality = {}
    for name, score in own_scores.items():
        bg = score["background"]
        point_physicality[name] = {
            "guard_z_range": [0.0, float(np.max(data["z_hd"]))],
            "guard_grid_points": int(len(bg["u_grid"])),
            "minimum_matter": bg["minimum_matter"],
            "minimum_vacuum": bg["minimum_vacuum"],
            "minimum_e2": bg["minimum_e2"],
            "fb_max_exact_from_monotone_ratio": bg["fb_max_exact"],
            "fb_witness": bg["fb_witness"],
            "minimum_cdm_witness": bg["minimum_cdm_witness"],
            "all_positive_on_sampled_guard": bool(bg["minimum_matter"] > 0 and
                                                     bg["minimum_vacuum"] > 0 and
                                                     bg["minimum_e2"] > 0 and
                                                     bg["minimum_cdm_witness"] > 0),
        }

    # Confirm the profile result did not change while this checker was running.
    result_sha_after = file_sha(RESULT_PATH)
    if result_sha_after != root_result_sha_before:
        raise RuntimeError("production result changed during independent reproduction; rerun against final hash")

    checker_hash = file_sha(Path(__file__).resolve())
    profile_hash = file_sha(PROFILE_PATH)
    hashes = {
        "hubble_diagram_sha256": file_sha(HD_PATH),
        "packed_stat_sys_precision_sha256": file_sha(P_PATH),
        "production_profile_sha256": profile_hash,
        "production_result_sha256": root_result_sha_before,
        "profile_contract_json_sha256": file_sha(CONTRACT_JSON),
        "profile_contract_markdown_sha256": file_sha(CONTRACT_MD),
        "frozen_bao_result_sha256": file_sha(BAO_RESULT),
        "bao_profile_source_sha256": file_sha(BAO_CODE),
        "production_sn_helper_sha256": file_sha(PROFILE_HELPER),
        "independent_checker_sha256": checker_hash,
    }
    command = (
        "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
        "NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 "
        ".venv/bin/python work/compute10/dovekie_ivs_profile_independent.py"
    )
    report = {
        "status": "independent_reproduction_passed",
        "method": "Standalone SNANA parser, independent full-precision unpack, contract-sign IVS ODE integrated with RK45, and independent full-covariance one-offset GLS profile.",
        "nonreuse_statement": "This checker imports none of profile.py, fit_dovekie.py, or their scoring/ODE/distance functions. Its parser, matrix unpacking, integration, distance conversion, intercept profile, and quadratic checks are implemented locally.",
        "source_contract": {
            "paper_convention": "Q=Gamma*rho_x>0 transfers CDM to vacuum; baryons are separately conserved.",
            "u": "ln(1+z)=-ln(a)",
            "ode": ["dx/du=-g*x/E", "dm/du=3*m+g*x/E", "E^2=m+x", "dD/du=exp(u)/E"],
            "rate_units": "g=Gamma/H0 is dimensionless and Gamma/H=g/E; c/H0_gauge is in Mpc.",
            "distance": "DL=(1+zHEL)*(c/H0_gauge)*D(zHD), H0_gauge=70 km/s/Mpc.",
        },
        "inputs": {
            "rows": EXPECTED_ROWS,
            "precision_representation": "Full released STAT+SYS inverse covariance P, unpacked from packed upper float32 values into float64 in original release order.",
            "precision_archive_keys": data["archive_keys"],
            "precision_cholesky_pass": True,
            "precision_q_one_p_one": data["q"],
            "zHD_range": [float(np.min(data["z_hd"])), float(np.max(data["z_hd"]))],
            "zHEL_range": [float(np.min(data["z_hel"])), float(np.max(data["z_hel"]))],
            "official_zHD_positive_cut_keeps_all_rows": True,
            "cid_order_sha256_lf_joined_trailing_lf": data["cid_hash"],
            "muerr_or_muerr_sys_added": False,
        },
        "hashes": hashes,
        "exact_command_from_project_root": command,
        "score_tolerance_absolute": SCORE_TOL,
        "g0_lcdm_limit": {
            "test_Omega_m0": test_om,
            "max_abs_distance_modulus_difference_mag_vs_analytic_LCDM_quadrature": g0_max_mu_diff,
            "tolerance_mag": DISTANCE_MODULUS_TOL,
            "passed": g0_max_mu_diff <= DISTANCE_MODULUS_TOL,
        },
        "source_sign_check": sign_check,
        "scores": comparisons,
        "physicality": point_physicality,
        "rk45_tolerance_refinement": {
            "standard": "rtol=2e-11, atol=2e-13, max_step=0.003",
            "tight": "rtol=2e-13, atol=2e-15, max_step=0.001",
            "points": convergence,
            "passed": True,
        },
        "interpretation_caveats": [
            "The checker reproduces finite-search profile scores at the reported minimum and fixed reference coordinates; it does not establish global optimality.",
            "The physicality check samples 4,001 points across the observed interval and checks the exact monotonic-ratio endpoint for f_b,max; it is not a formal continuous-ODE positivity proof.",
            "The BAO IVS and BAO LCDM coordinates are frozen external reference points. These are SN-only scores/differences, not a BAO×SN combined likelihood, calibrated interval, or model evidence.",
            "The H0 value is a common distance gauge absorbed by the profiled magnitude intercept; Dovekie alone does not identify dimensional Gamma.",
        ],
        "runtime_seconds": time.perf_counter() - start,
        "exact_root_result_hash_after_check": result_sha_after,
    }
    report["scores"][0]["independent_delta_from_ivs_minimum"] = 0.0
    report["production_delta_reproduction"] = {
        "bao_reference_delta_difference": abs(
            comparisons[1]["independent_delta_from_ivs_minimum"]
            - float(result["results"]["delta_chi2_bao_point_minus_sn_minimum"])
        ),
        "sn_lcdm_reference_delta_difference": abs(
            comparisons[2]["independent_delta_from_ivs_minimum"]
            - float(result["results"]["delta_chi2_lcdm_minimum_minus_ivs_minimum"])
        ),
    }
    if max(report["production_delta_reproduction"].values()) > SCORE_TOL:
        raise ArithmeticError("profile delta reproduction exceeded score tolerance")
    OUT_JSON.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    score_lines = []
    for item in comparisons:
        delta_text = item["independent_delta_from_ivs_minimum"]
        delta_text = "—" if delta_text is None else f"{delta_text:.12f}"
        difference = item["absolute_score_difference"]
        difference_text = "—" if difference is None else f"{difference:.3g}"
        score_lines.append(
            f"| {item['name']} | {item['parameters']['Omega_m0']:.12g} | {item['parameters']['g']:.12g} "
            f"| {item['independent_chi2']:.12f} | "
            f"{delta_text} | {difference_text} |"
        )
    md = "\n".join([
        "# Independent Dovekie IVS profile reproduction",
        "",
        "Status: **passed**. This standalone checker parsed the SNANA release order and unpacked the full STAT+SYS precision itself. It imports none of the production profile, SN helper, ODE, distance, or scoring functions.",
        "",
        "The local ODE uses the source convention `Q=Gamma*rho_x>0` (CDM to vacuum), with `g=Gamma/H0`, `Gamma/H=g/E`, `dx/du=-g*x/E`, `dm/du=3*m+g*x/E`, and `dD/du=exp(u)/E`. Distances use `zHD` in the integral and `zHEL` in the luminosity prefactor at a fixed 70 km/s/Mpc gauge. One unbounded intercept is profiled by GLS over the full 1,820-row precision; no MUERR term is added.",
        "",
        f"Input CID-order fingerprint: `{data['cid_hash']}`. Both pinned data files and all code/result identities are in the JSON report. Full precision Cholesky passed; the independent intercept normal equation and direct/whitened/Schur quadratic evaluations are recorded per point.",
        "",
        "| Point | Ωm0 | g | Independent χ²profile | Δχ² from IVS minimum | Abs. χ² difference from production |",
        "|---|---:|---:|---:|---:|---:|",
        *score_lines,
        "",
        f"The analytic flat-LCDM distance limit at g=0 agrees with the independent IVS ODE to {g0_max_mu_diff:.3g} mag (tolerance {DISTANCE_MODULUS_TOL:g}); the source-term direction and standard-vs-tight RK45 refinements pass. Root-reported score and delta discrepancies are below {SCORE_TOL:g}.",
        "",
        "Physicality was checked over the observed 0 ≤ zHD ≤ 1.14418 range using 4,001 trajectory samples, with positive m, x, E² and an interior constant baryon/CDM split witness. This is a dense numerical check, not a formal continuous positivity proof. The reported finite-search minimum is reproduced at its coordinate; this check does not establish global optimality or turn the fixed BAO references into a joint likelihood.",
        "",
        f"Exact command from `{ROOT}`: `{command}`",
        "",
    ])
    OUT_MD.write_text(md, encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "result_sha256": root_result_sha_before,
        "ivs_minimum": [points["ivs_minimum"][0], points["ivs_minimum"][1], own_scores["ivs_minimum"]["chi2"]],
        "bao_delta": comparisons[1]["independent_delta_from_ivs_minimum"],
        "sn_lcdm_delta": comparisons[2]["independent_delta_from_ivs_minimum"],
        "bao_lcdm_reference_chi2": own_scores["bao_lcdm_reference"]["chi2"],
        "g0_max_mu_diff": g0_max_mu_diff,
        "outputs": [str(OUT_JSON.relative_to(ROOT)), str(OUT_MD.relative_to(ROOT))],
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

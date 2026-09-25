#!/usr/bin/env python3
"""No-fit local GLS identifiability audit at the pinned Dovekie IVS point.

This is a standalone response-derivative calculation. It does not import the
production profile/scorer/ODE, run an optimizer, or generate mocks. A forward
sensitivity system for (m, x, D, d/dOmega_m, d/dg) is integrated alongside the
background; central finite differences of a separate base ODE are used only as
a convergence check. All artifacts stay in work/compute13/.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parents[2]
HD_PATH = ROOT / "context/data/des_dovekie_hd.csv"
P_PATH = ROOT / "context/data/des_dovekie_stat_sys.npz"
PROFILE_RESULT_PATH = ROOT / "experiments/dovekie_ivs_profile/result.json"
OUT_JSON = ROOT / "work/compute13/dovekie_profile_identifiability.json"
OUT_MD = ROOT / "work/compute13/dovekie_profile_identifiability.md"

OM0 = 0.39301742026988545
G0 = -0.6293878145169555
H0_GAUGE = 70.0  # km s^-1 Mpc^-1
C_KM_S = 299792.458
EXPECTED_HD_SHA = "2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf"
EXPECTED_P_SHA = "ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b"
N_EXPECTED = 1820
FD_STEPS = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_hd(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    columns: list[str] | None = None
    records: list[dict[str, str]] = []
    for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        f = line.split()
        if not f or f[0].startswith("#"):
            continue
        if f[0] == "VARNAMES:":
            if columns is not None:
                raise ValueError(f"duplicate VARNAMES header at line {line_no}")
            columns = f[1:]
        elif f[0] == "SN:":
            if columns is None or len(f[1:]) != len(columns):
                raise ValueError(f"bad SN row at line {line_no}")
            records.append(dict(zip(columns, f[1:], strict=True)))
        else:
            raise ValueError(f"unexpected row type {f[0]!r} at line {line_no}")
    if columns is None:
        raise ValueError("missing VARNAMES header")
    return columns, records


def load_inputs() -> dict:
    hd_sha, p_sha = sha256(HD_PATH), sha256(P_PATH)
    if (hd_sha, p_sha) != (EXPECTED_HD_SHA, EXPECTED_P_SHA):
        raise ValueError(f"pinned input hash mismatch: HD={hd_sha}, precision={p_sha}")
    columns, rows = read_hd(HD_PATH)
    required = {"CID", "zHD", "zHEL", "MU"}
    if not required.issubset(columns):
        raise ValueError(f"missing columns {sorted(required - set(columns))}")
    if len(rows) != N_EXPECTED or len({r["CID"] for r in rows}) != N_EXPECTED:
        raise ValueError("Hubble diagram row count or CID uniqueness mismatch")
    selected = [r for r in rows if float(r["zHD"]) > 0.0]
    if len(selected) != N_EXPECTED:
        raise ValueError("official zHD>0 selection does not retain all rows")
    zhd = np.array([float(r["zHD"]) for r in selected], dtype=np.float64)
    zhel = np.array([float(r["zHEL"]) for r in selected], dtype=np.float64)
    mu = np.array([float(r["MU"]) for r in selected], dtype=np.float64)
    cid_hash = hashlib.sha256(("\n".join(r["CID"] for r in rows) + "\n").encode()).hexdigest()

    with np.load(P_PATH, allow_pickle=False) as archive:
        n = int(np.asarray(archive["nsn"]).reshape(-1)[0])
        packed = np.asarray(archive["cov"])
        keys = list(archive.files)
    if n != N_EXPECTED or packed.size != n * (n + 1) // 2 or packed.dtype != np.float32:
        raise ValueError("released precision archive differs from expected packed format")
    precision = np.zeros((n, n), dtype=np.float64)
    upper = np.triu_indices(n)
    precision[upper] = packed.astype(np.float64)
    precision[(upper[1], upper[0])] = precision[upper]
    if not np.array_equal(precision, precision.T) or not np.all(np.isfinite(precision)):
        raise ValueError("unpacked precision is nonfinite or asymmetric")
    chol = np.linalg.cholesky(precision)
    one = np.ones(n, dtype=np.float64)
    p_one = precision @ one
    q1 = float(one @ p_one)
    if not np.isfinite(q1) or q1 <= 0:
        raise ValueError("intercept direction is not positively constrained")
    return {"zHD": zhd, "zHEL": zhel, "MU": mu, "CID_sha256": cid_hash,
            "P": precision, "L": chol, "q1": q1, "archive_keys": keys,
            "HD_sha256": hd_sha, "P_sha256": p_sha}


def rhs_sens(u: float, y: np.ndarray, om: float, g: float) -> np.ndarray:
    """ODE for base background and its analytic forward parameter sensitivities."""
    m, x, _d = y[:3]
    e2 = m + x
    if not (math.isfinite(e2) and e2 > 0.0):
        raise FloatingPointError("E^2 <= 0 in sensitivity integration")
    e = math.sqrt(e2)
    source = g * x / e
    out = np.empty_like(y)
    out[0] = 3.0 * m + source
    out[1] = -source
    out[2] = math.exp(u) / e
    # State blocks are (m_p, x_p, D_p), first p=Omega_m0, then p=g.
    for j, is_g in ((3, False), (6, True)):
        mp, xp, _dp = y[j:j + 3]
        source_p = (x / e if is_g else 0.0) + g * (
            xp / e - x * (mp + xp) / (2.0 * e**3)
        )
        out[j] = 3.0 * mp + source_p
        out[j + 1] = -source_p
        out[j + 2] = -math.exp(u) * (mp + xp) / (2.0 * e**3)
    return out


def solve_sens(om: float, g: float, zhd: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    umax = math.log1p(float(np.max(zhd)))
    # m(0)=om, x(0)=1-om, D(0)=0 and their exact initial derivatives.
    y0 = np.array([om, 1.0 - om, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    sol = solve_ivp(lambda u, y: rhs_sens(u, y, om, g), (0.0, umax), y0,
                    method="RK45", dense_output=True, rtol=3e-13, atol=3e-15,
                    max_step=0.003)
    if not sol.success or sol.sol is None:
        raise RuntimeError(f"sensitivity ODE failed: {sol.message}")
    state = np.asarray(sol.sol(np.log1p(zhd)), dtype=np.float64)
    d = state[2]
    dp = np.vstack((state[5] / d, state[8] / d)).T * (5.0 / math.log(10.0))
    # Return dmu/d(om,g), after the log derivative of D; zHEL and H0 are common factors.
    return d, dp


def solve_base_mu(om: float, g: float, zhd: np.ndarray, zhel: np.ndarray) -> np.ndarray:
    """Independent base-only integration used for the central-FD convergence sweep."""
    umax = math.log1p(float(np.max(zhd)))

    def rhs(u: float, state: np.ndarray) -> np.ndarray:
        m, x, _d = state
        e2 = m + x
        if not (math.isfinite(e2) and e2 > 0.0):
            raise FloatingPointError("E^2 <= 0 in base integration")
        e = math.sqrt(e2)
        src = g * x / e
        return np.array([3.0 * m + src, -src, math.exp(u) / e], dtype=np.float64)

    sol = solve_ivp(rhs, (0.0, umax), np.array([om, 1.0 - om, 0.0]),
                    method="RK45", dense_output=True, rtol=3e-13, atol=3e-15,
                    max_step=0.003)
    if not sol.success or sol.sol is None:
        raise RuntimeError(f"base ODE failed: {sol.message}")
    dist = np.asarray(sol.sol(np.log1p(zhd))[2], dtype=np.float64)
    dl_mpc = (1.0 + zhel) * (C_KM_S / H0_GAUGE) * dist
    return 5.0 * np.log10(dl_mpc) + 25.0


def rank_aware_basis(a: np.ndarray) -> tuple[np.ndarray, int, list[float], float]:
    """Return left singular basis, numerical rank, singular values and cutoff."""
    u, s, _vh = np.linalg.svd(a, full_matrices=False)
    s0 = float(s[0]) if s.size else 0.0
    tol = float(np.finfo(np.float64).eps * max(a.shape, default=0) * s0)
    rank = int(np.sum(s > tol))
    return u[:, :rank], rank, s.tolist(), tol


def project_out(v: np.ndarray, a: np.ndarray) -> tuple[np.ndarray, int, list[float], float]:
    basis, rank, singular_values, tol = rank_aware_basis(a)
    return v - basis @ (basis.T @ v), rank, singular_values, tol


def main() -> None:
    dat = load_inputs()
    zhd, zhel, y = dat["zHD"], dat["zHEL"], dat["MU"]
    p, l = dat["P"], dat["L"]
    n = len(y)

    d, dmu = solve_sens(OM0, G0, zhd)
    if np.any(d <= 0) or np.any(~np.isfinite(dmu)):
        raise FloatingPointError("invalid distance or sensitivity")
    mu0 = 5.0 * np.log10((1.0 + zhel) * (C_KM_S / H0_GAUGE) * d) + 25.0

    # Profile one common additive magnitude intercept by its exact GLS normal equation.
    delta = mu0 - y
    mhat = -float(np.ones(n) @ (p @ delta)) / dat["q1"]
    residual = delta + mhat * np.ones(n)
    normal_residual = float(np.ones(n) @ (p @ residual))
    chi2 = float(residual @ (p @ residual))

    # If P=L L^T, the Euclidean-whitened derivative is L^T dmu, not L^-1 dmu.
    w_intercept = l.T @ np.ones(n)
    w_shape = l.T @ dmu
    raw_cosine = float((w_shape[:, 0] @ w_shape[:, 1]) /
                       (np.linalg.norm(w_shape[:, 0]) * np.linalg.norm(w_shape[:, 1])))
    projected_shape, intercept_rank, intercept_sv, intercept_tol = project_out(
        w_shape, w_intercept[:, None])
    post_cosine = float((projected_shape[:, 0] @ projected_shape[:, 1]) /
                        (np.linalg.norm(projected_shape[:, 0]) * np.linalg.norm(projected_shape[:, 1])))
    sv_shape = np.linalg.svd(projected_shape, compute_uv=False)
    cond_shape = float(sv_shape[0] / sv_shape[-1])
    gram_after_intercept = projected_shape.T @ projected_shape

    # Efficient local information in g with Omega_m0 and the magnitude intercept
    # both treated as free nuisance coordinates.
    nuisance = np.column_stack((w_intercept, w_shape[:, 0]))
    wg_cond, nuisance_rank, nuisance_sv, nuisance_tol = project_out(
        w_shape[:, 1], nuisance)
    info_g_cond = float(wg_cond @ wg_cond)
    schur_info = float(gram_after_intercept[1, 1] -
                       gram_after_intercept[0, 1]**2 / gram_after_intercept[0, 0])
    info_gap = abs(info_g_cond - schur_info)
    g_info_before_om = float(gram_after_intercept[1, 1])

    # Central finite differences of the independent base-only RK45 ODE.
    fd_rows = []
    for h in FD_STEPS:
        derivs = []
        for col, center in ((0, OM0), (1, G0)):
            plus = [OM0, G0]
            minus = [OM0, G0]
            plus[col] = center + h
            minus[col] = center - h
            mup = solve_base_mu(plus[0], plus[1], zhd, zhel)
            mum = solve_base_mu(minus[0], minus[1], zhd, zhel)
            derivs.append((mup - mum) / (2.0 * h))
        fd = np.column_stack(derivs)
        wfd = l.T @ fd
        errs = []
        raw_errs = []
        cosines = []
        for j in range(2):
            diff = wfd[:, j] - w_shape[:, j]
            errs.append(float(np.linalg.norm(diff) / np.linalg.norm(w_shape[:, j])))
            raw_errs.append(float(np.linalg.norm(fd[:, j] - dmu[:, j]) / np.linalg.norm(dmu[:, j])))
            cosines.append(float((wfd[:, j] @ w_shape[:, j]) /
                                 (np.linalg.norm(wfd[:, j]) * np.linalg.norm(w_shape[:, j]))))
        fd_projected, fd_rank, _fd_isv, _fd_itol = project_out(wfd, w_intercept[:, None])
        fd_shape_sv = np.linalg.svd(fd_projected, compute_uv=False)
        fd_rows.append({
            "step_both_parameters": h,
            "whitened_relative_derivative_error": errs,
            "raw_relative_derivative_error": raw_errs,
            "whitened_derivative_cosine_with_sensitivity": cosines,
            "intercept_rank": fd_rank,
            "projected_shape_singular_values": fd_shape_sv.tolist(),
            "projected_shape_condition": float(fd_shape_sv[0] / fd_shape_sv[-1]),
        })

    # Directly verify the GLS-whitening and projection algebra numerically.
    test_v = np.column_stack((dmu[:, 0], dmu[:, 1], residual))
    whitened_test = l.T @ test_v
    direct_gram = test_v.T @ p @ test_v
    whitened_gram = whitened_test.T @ whitened_test
    gls_identity_abs = float(np.max(np.abs(direct_gram - whitened_gram)))
    residual_white = l.T @ residual
    residual_white_norm2 = float(residual_white @ residual_white)
    schur_projector_matrix = p - np.outer(p @ np.ones(n), p @ np.ones(n)) / dat["q1"]
    direct_profile_gram = dmu.T @ schur_projector_matrix @ dmu
    profile_gram_abs = float(np.max(np.abs(direct_profile_gram - gram_after_intercept)))
    normal_residual_white = float(w_intercept @ residual_white)

    # Descriptive residual alignment with the nuisance-orthogonal g direction.
    residual_g_dot = float(residual_white @ wg_cond)
    residual_g_cos = float(residual_g_dot /
                           (np.linalg.norm(residual_white) * np.linalg.norm(wg_cond)))
    residual_g_projection = float(residual_g_dot / math.sqrt(info_g_cond))

    source_score = None
    if PROFILE_RESULT_PATH.exists():
        source = json.loads(PROFILE_RESULT_PATH.read_text())
        results_section = source.get("results", {})
        nested_minimum = results_section.get("ivs_sn_only_minimum", {}) if isinstance(results_section, dict) else {}
        nested_profile = nested_minimum.get("profile") if isinstance(nested_minimum, dict) else None
        for value in (nested_profile, source.get("ivs_sn_only_minimum"), source.get("fit"),
                      source.get("ivs_fit"), source.get("profile_fit")):
            if isinstance(value, dict) and "chi2" in value:
                source_score = float(value["chi2"])
                break
        if source_score is None:
            # Locate the profile point by its top-level naming without traversing arrays.
            for value in source.values():
                if isinstance(value, dict) and "chi2" in value and (
                        abs(float(value.get("Omega_m0", value.get("Omega_m", math.inf))) - OM0) < 1e-14):
                    source_score = float(value["chi2"])
                    break

    result = {
        "status": "completed_no_fit_local_identifiability_audit",
        "scope": "Local response geometry at the pinned exploratory Dovekie-only IVS profile coordinate; not a posterior, uncertainty, significance, or evidence calculation.",
        "pinned_point": {"Omega_m0": OM0, "g": G0, "H0_gauge_km_s_Mpc": H0_GAUGE},
        "inputs": {
            "hubble_diagram": str(HD_PATH.relative_to(ROOT)), "hubble_diagram_sha256": dat["HD_sha256"],
            "precision_archive": str(P_PATH.relative_to(ROOT)), "precision_archive_sha256": dat["P_sha256"],
            "rows": n, "cid_order_sha256": dat["CID_sha256"], "precision_archive_keys": dat["archive_keys"],
            "zHD_range": [float(np.min(zhd)), float(np.max(zhd))],
            "zHEL_range": [float(np.min(zhel)), float(np.max(zhel))],
            "MUERR_or_MUERR_SYS_added": False,
        },
        "independent_model": {
            "integrator": "SciPy solve_ivp RK45; standalone implementation, no production profile/scorer imports",
            "background": ["dm/du=3m+g*x/E", "dx/du=-g*x/E", "dD/du=exp(u)/E", "E^2=m+x"],
            "sensitivity": "Analytic forward sensitivities integrated in the same ODE, with dE/du dependence differentiated explicitly.",
            "mu_definition": "5 log10((1+zHEL)*(c/H0_gauge)*D(zHD)/Mpc)+25",
            "derivative_units": "mag per unit dimensionless Omega_m0 or g; H0 and zHEL prefactors are common and parameter-independent.",
        },
        "profiled_intercept_check": {
            "Mhat_mag": mhat, "chi2_at_pinned_point_independently_recomputed": chi2,
            "reference_profile_chi2_if_read": source_score,
            "abs_difference_from_reference_profile_chi2": None if source_score is None else abs(chi2 - source_score),
            "normal_equation_abs_1TP_r": abs(normal_residual),
            "normal_equation_abs_after_whitening": abs(normal_residual_white),
            "residual_chi2_whitened": residual_white_norm2,
        },
        "gls_algebra": {
            "precision_factorization": "P=L L^T, with L lower Cholesky",
            "whitening_map": "v -> L^T v, since ||L^T v||^2=v^T P v",
            "max_abs_direct_vs_whitened_3x3_gram": gls_identity_abs,
            "intercept_projector_rank": intercept_rank,
            "intercept_rank_tolerance": intercept_tol,
            "intercept_whitened_singular_values": intercept_sv,
            "max_abs_direct_schur_vs_whitened_projected_shape_gram": profile_gram_abs,
        },
        "derivative_alignment": {
            "raw_whitened_norms_before_intercept_projection": [float(np.linalg.norm(w_shape[:, 0])), float(np.linalg.norm(w_shape[:, 1]))],
            "cosine_before_intercept_projection": raw_cosine,
            "whitened_norms_after_intercept_projection": [float(np.linalg.norm(projected_shape[:, 0])), float(np.linalg.norm(projected_shape[:, 1]))],
            "cosine_after_intercept_projection": post_cosine,
            "absolute_cosine_after_intercept_projection": abs(post_cosine),
            "interpretation": "Signed cosine near +1/-1 means local response vectors are nearly parallel after removing the common magnitude direction; it is descriptive local shape geometry.",
        },
        "local_shape_information": {
            "parameters": ["Omega_m0", "g"],
            "projected_response_singular_values": sv_shape.tolist(),
            "projected_response_condition_number": cond_shape,
            "precision_weighted_information_matrix_after_intercept": gram_after_intercept.tolist(),
            "information_units": "inverse squared parameter units; both shape coordinates are dimensionless.",
        },
        "g_conditional_information": {
            "nuisance_columns_projected_out": ["common magnitude intercept", "Omega_m0 response"],
            "nuisance_design_rank": nuisance_rank,
            "nuisance_design_singular_values": nuisance_sv,
            "nuisance_rank_tolerance": nuisance_tol,
            "efficient_local_information_for_g": info_g_cond,
            "schur_complement_information_for_g": schur_info,
            "abs_projector_vs_schur_difference": info_gap,
            "raw_g_information_after_intercept_before_Omega_projection": g_info_before_om,
            "fraction_of_after_intercept_g_information_retained_after_Omega_projection": info_g_cond / g_info_before_om,
        },
        "finite_difference_convergence": {
            "method": "Central differences of independently integrated base-only distance-modulus vectors; same tight RK45 tolerance, max step 0.003 in u.",
            "sensitivity_reference": "Analytic forward-sensitivity ODE at the pinned point.",
            "parameter_step_values": list(FD_STEPS),
            "rows": fd_rows,
            "best_whitened_relative_error_by_parameter": [
                min(row["whitened_relative_derivative_error"][j] for row in fd_rows) for j in range(2)
            ],
        },
        "optional_residual_alignment": {
            "scope": "Descriptive diagnostic only; it is not a significance or evidence calculation.",
            "residual_cosine_with_g_after_intercept_and_Omega_projection": residual_g_cos,
            "signed_whitened_residual_projection_on_unit_conditional_g_direction": residual_g_projection,
        },
        "interpretation_limits": [
            "Local derivative geometry at one fixed fitted coordinate; no optimization, mocks, posterior, calibrated interval, significance, or evidence.",
            "The condition number depends on the stated dimensionless parameter coordinates and is not by itself an uncertainty.",
            "This SN-only late-time background diagnostic does not establish early-time or perturbation viability and does not combine BAO with SN.",
        ],
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    report = f"""# Dovekie IVS local identifiability audit

**Status:** completed standalone no-fit local response audit at the pinned point
\\((\\Omega_{{m0}},g)=({OM0:.14f},{G0:.14f})\\). It is local linear GLS geometry,
not a posterior uncertainty, significance, or evidence calculation.

The SN model response derivatives were obtained by an independent RK45 forward-
sensitivity ODE. Central finite differences of a separately written base-only
ODE were checked over steps from `1e-2` to `1e-6`. The full 1,820 by 1,820
released STAT+SYS precision was unpacked in release order and factored as
`P = L L^T`; whitening used `L.T @ v`. One common magnitude intercept was
removed by an SVD rank-aware projection.

After removing the intercept, the signed cosine between the whitened
\\(\\Omega_{{m0}}\\) and \\(g\\) responses is **{post_cosine:.8f}**
(raw, before intercept removal: {raw_cosine:.8f}). Singular values of the
two-column projected response are **{sv_shape[0]:.6g}, {sv_shape[1]:.6g}**,
with condition number **{cond_shape:.6g}** in the stated dimensionless
parameter coordinates.

After projecting out both the intercept and the \\(\\Omega_{{m0}}\\) response,
the conditional local information for \\(g\\) is **{info_g_cond:.8g}**. The
Schur-complement value is {schur_info:.8g}; their absolute difference is
{info_gap:.2e}. The finite-difference sweep's best whitened relative derivative
errors are {min(row['whitened_relative_derivative_error'][0] for row in fd_rows):.2e}
for \\(\\Omega_{{m0}}\\) and {min(row['whitened_relative_derivative_error'][1] for row in fd_rows):.2e}
for \\(g\\).

The independent GLS intercept is {mhat:.8g} mag and the normal-equation
residual is {abs(normal_residual):.2e}; this re-evaluation only supports the
descriptive residual alignment and is not a new fit. The optional whitened
residual cosine with the nuisance-orthogonal \\(g\\) direction is
{residual_g_cos:.6g}, reported descriptively only.
The independently recomputed fixed-point \\(\\chi^2_{{\\rm prof}}\\) is
{chi2:.12f}, matching the pinned profile record to
{abs(chi2 - source_score):.2e} if that record is available.

The independent central-difference convergence sweep is:

| Step | Whitened relative error, \\(\\Omega_{{m0}}\\) | Whitened relative error, \\(g\\) | Projected response condition |
|---:|---:|---:|---:|
{chr(10).join(f"| {row['step_both_parameters']:.0e} | {row['whitened_relative_derivative_error'][0]:.2e} | {row['whitened_relative_derivative_error'][1]:.2e} | {row['projected_shape_condition']:.6g} |" for row in fd_rows)}

The alignment, singular values, condition number, and conditional information
describe the chosen local response coordinates. They do not establish a
posterior width, calibrated constraint, evidence, or physical viability beyond
the scoped low-redshift model. Full numerical values and every finite-difference
step are in `dovekie_profile_identifiability.json`.
"""
    OUT_MD.write_text(report, encoding="utf-8")
    print(json.dumps({
        "output_json": str(OUT_JSON.relative_to(ROOT)),
        "output_md": str(OUT_MD.relative_to(ROOT)),
        "cosine_after_intercept": post_cosine,
        "singular_values": sv_shape.tolist(),
        "condition_number": cond_shape,
        "g_conditional_information": info_g_cond,
        "finite_difference_best_relative_errors": result["finite_difference_convergence"]["best_whitened_relative_error_by_parameter"],
        "max_gls_gram_identity_abs": gls_identity_abs,
        "max_projected_information_identity_abs": profile_gram_abs,
        "chi2_at_pinned_point": chi2,
    }, indent=2))


if __name__ == "__main__":
    main()

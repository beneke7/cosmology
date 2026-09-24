#!/usr/bin/env python3
"""Independent fixed-point and bounded-profile audit of DES-Dovekie.

This implementation intentionally does not import any work/inference module.
It parses the release SNANA records in source order, unpacks the release's
float32 upper-triangle precision, and evaluates distances and profile
quadratics directly.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.optimize import differential_evolution


ROOT = Path(__file__).resolve().parents[2]
HD_FILE = ROOT / "context/data/des_dovekie_hd.csv"
PRECISION_FILE = ROOT / "context/data/des_dovekie_stat_sys.npz"
OUT_FILE = ROOT / "work/critic/independent_result.json"
EXPECTED = {
    "des_dovekie_hd.csv": "2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf",
    "des_dovekie_stat_sys.npz": "ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b",
}
C_KM_S = 299792.458
H0_GAUGE = 70.0
GL_X, GL_W = leggauss(128)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_release(path: Path):
    names = None
    rows = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = raw.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "VARNAMES:":
            if names is not None:
                raise ValueError(f"repeated VARNAMES at line {number}")
            names = fields[1:]
        elif fields[0] == "SN:":
            if names is None or len(fields[1:]) != len(names):
                raise ValueError(f"bad SN record at line {number}")
            rows.append(dict(zip(names, fields[1:], strict=True)))
        else:
            raise ValueError(f"unknown record at line {number}: {fields[0]}")
    if names is None:
        raise ValueError("missing VARNAMES")
    return names, rows


def unpack_precision(path: Path, nrows: int):
    with np.load(path, allow_pickle=False) as z:
        if set(z.files) != {"nsn", "cov", "allow_pickle"}:
            raise ValueError(f"unexpected NPZ members: {z.files}")
        n = int(z["nsn"][0])
        packed = np.asarray(z["cov"])
    if n != nrows or packed.dtype != np.float32 or packed.size != n * (n + 1) // 2:
        raise ValueError(f"packed precision mismatch: n={n}, rows={nrows}, dtype={packed.dtype}, size={packed.size}")
    p = np.zeros((n, n), dtype=np.float64)
    ij = np.triu_indices(n)
    p[ij] = packed.astype(np.float64)
    p[(ij[1], ij[0])] = packed.astype(np.float64)
    if not np.array_equal(p, p.T) or not np.all(np.isfinite(p)):
        raise ValueError("unpacked precision is not finite symmetric")
    # A factorization is a direct positive-definiteness check of the released P.
    np.linalg.cholesky(p)
    return p, packed


def e_squared(z: np.ndarray, om: float, kind: str, w0: float, wa: float) -> np.ndarray:
    zp1 = 1.0 + z
    matter = om * zp1**3
    if kind == "LCDM":
        log_rho_de = np.zeros_like(z)
    elif kind == "constant_w":
        log_rho_de = 3.0 * (1.0 + w0) * np.log1p(z)
    elif kind == "CPL":
        log_rho_de = 3.0 * (1.0 + w0 + wa) * np.log1p(z) - 3.0 * wa * z / zp1
    else:
        raise ValueError(kind)
    answer = matter + (1.0 - om) * np.exp(log_rho_de)
    if np.any(answer <= 0.0) or not np.all(np.isfinite(answer)):
        raise FloatingPointError("nonpositive or nonfinite E(z)^2")
    return answer


def mu_theory(zhd: np.ndarray, zhel: np.ndarray, params, kind: str, prefactor: str = "zHEL"):
    om = float(params[0])
    w0 = -1.0 if kind == "LCDM" else float(params[1])
    wa = 0.0 if kind != "CPL" else float(params[2])
    zz = zhd[:, None] * (GL_X[None, :] + 1.0) * 0.5
    integral = zhd * 0.5 * np.sum(GL_W[None, :] / np.sqrt(e_squared(zz, om, kind, w0, wa)), axis=1)
    zouter = zhel if prefactor == "zHEL" else zhd
    dl = (1.0 + zouter) * (C_KM_S / H0_GAUGE) * integral
    return 5.0 * np.log10(dl) + 25.0


def profile(model: np.ndarray, observed: np.ndarray, precision: np.ndarray, p_one: np.ndarray, c_one: float):
    delta = model - observed
    # Center the vector before applying the ill-conditioned full precision. A
    # constant is absorbed exactly by the analytically profiled magnitude.
    anchor = float(delta[0])
    centered = delta - anchor
    shift = -float(p_one @ centered) / c_one
    resid = centered + shift
    qresid = precision @ resid
    chi2_resid = float(resid @ qresid)
    offset = shift - anchor
    return offset, chi2_resid, centered, resid


def main():
    started = time.perf_counter()
    hashes = {HD_FILE.name: sha256(HD_FILE), PRECISION_FILE.name: sha256(PRECISION_FILE)}
    if hashes != EXPECTED:
        raise ValueError(f"pinned input hash mismatch: {hashes}")
    columns, rows = parse_release(HD_FILE)
    required = {"CID", "zHD", "zHEL", "MU", "MUERR"}
    if not required.issubset(columns):
        raise ValueError(f"missing columns: {sorted(required - set(columns))}")
    cids = [x["CID"] for x in rows]
    if len(set(cids)) != len(cids):
        raise ValueError("duplicate CID")
    zhd = np.asarray([float(x["zHD"]) for x in rows], dtype=np.float64)
    zhel = np.asarray([float(x["zHEL"]) for x in rows], dtype=np.float64)
    muobs = np.asarray([float(x["MU"]) for x in rows], dtype=np.float64)
    muerr = np.asarray([float(x["MUERR"]) for x in rows], dtype=np.float64)
    if not (np.all(np.isfinite(zhd)) and np.all(np.isfinite(zhel)) and np.all(np.isfinite(muobs)) and np.all(muerr > 0.0)):
        raise ValueError("invalid data values")
    keep = zhd > 0.0
    if not np.all(keep):
        raise ValueError("official zHD>0 selection would invalidate positional covariance pairing")
    precision, packed = unpack_precision(PRECISION_FILE, len(rows))
    precision_cholesky = np.linalg.cholesky(precision)
    one = np.ones(len(rows), dtype=np.float64)
    p_one = precision @ one
    c_one = float(one @ p_one)

    points = {
        "LCDM": ("LCDM", np.array([0.3303167905483107]), [(0.05, 0.6)]),
        "constant_w": ("constant_w", np.array([0.2610419, -0.83268]), [(0.05, 0.6), (-2.0, -0.3)]),
        "CPL_wa_minus_3": ("CPL", np.array([0.409838, -0.788386, -3.0]), [(0.05, 0.6), (-2.0, -0.3), (-3.0, 3.0)]),
        "CPL_widened": ("CPL", np.array([0.465666, -0.538859, -6.7818]), [(0.05, 0.6), (-2.0, -0.3), (-10.0, 3.0)]),
    }
    fixed = {}
    for label, (kind, pars, _) in points.items():
        model = mu_theory(zhd, zhel, pars, kind)
        offset, chi2, centered, resid = profile(model, muobs, precision, p_one, c_one)
        p_centered = precision @ centered
        schur = float(centered @ p_centered - (p_one @ centered) ** 2 / c_one)
        normal = float(p_one @ resid)
        whiten = np.linalg.norm(precision_cholesky.T @ resid) ** 2
        wrong_zhd_model = mu_theory(zhd, zhel, pars, kind, prefactor="zHD")
        wrong_offset, wrong_chi2, _, _ = profile(wrong_zhd_model, muobs, precision, p_one, c_one)
        fixed[label] = {
            "model": kind,
            "parameters": [float(v) for v in pars],
            "profile_offset_mag": offset,
            "profile_chi2": chi2,
            "schur_chi2": schur,
            "cholesky_whitened_chi2": float(whiten),
            "offset_normal_equation": normal,
            "max_algebra_difference": max(abs(chi2 - schur), abs(chi2 - whiten)),
            "if_zHD_used_in_prefactor_chi2": wrong_chi2,
            "if_zHD_used_in_prefactor_offset_mag": wrong_offset,
        }

    # Independent bounded global search: seeded differential evolution followed
    # by its bounded polishing stage. Search boxes match the reported screens.
    fit_domains = {
        "LCDM": ("LCDM", [(0.05, 0.6)]),
        "constant_w": ("constant_w", [(0.05, 0.6), (-2.0, -0.3)]),
        "CPL_wa_ge_minus3": ("CPL", [(0.05, 0.6), (-2.0, -0.3), (-3.0, 3.0)]),
        "CPL_wa_ge_minus10": ("CPL", [(0.05, 0.6), (-2.0, -0.3), (-10.0, 3.0)]),
    }
    fits = {}
    for ix, (label, (kind, bounds)) in enumerate(fit_domains.items()):
        def objective(x):
            try:
                theory = mu_theory(zhd, zhel, x, kind)
                return profile(theory, muobs, precision, p_one, c_one)[1]
            except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
                return 1.0e100

        fit = differential_evolution(
            objective,
            bounds,
            strategy="best1bin",
            maxiter=65,
            popsize=8,
            tol=2.0e-10,
            atol=1.0e-9,
            mutation=(0.5, 1.0),
            recombination=0.7,
            seed=90210 + ix,
            polish=True,
            updating="immediate",
            workers=1,
        )
        fits[label] = {
            "parameters": [float(v) for v in fit.x],
            "profile_chi2": float(fit.fun),
            "success": bool(fit.success),
            "message": str(fit.message),
            "function_evaluations": int(fit.nfev),
            "iterations": int(fit.nit),
            "bounds": [[float(a), float(b)] for a, b in bounds],
        }

    # Sensitivity probe: preserve the distance integral in zHD but substitute
    # zHD for zHEL in the luminosity-distance prefactor, a common convention error.
    zhd_probe = {}
    for label, (kind, pars, bounds) in points.items():
        def objective_wrong(x):
            theory = mu_theory(zhd, zhel, x, kind, prefactor="zHD")
            return profile(theory, muobs, precision, p_one, c_one)[1]

        wrong_fit = differential_evolution(
            objective_wrong,
            bounds,
            maxiter=40,
            popsize=7,
            tol=1.0e-8,
            seed=8711 + len(zhd_probe),
            polish=True,
            workers=1,
        )
        zhd_probe[label] = {
            "parameters": [float(v) for v in wrong_fit.x],
            "profile_chi2": float(wrong_fit.fun),
            "success": bool(wrong_fit.success),
        }

    packed_order_digest = hashlib.sha256("\n".join(cids).encode("utf-8")).hexdigest()
    result = {
        "status": "independent_profile_validation",
        "inputs": {
            "sha256": hashes,
            "pinned_release_commit": "c9a4fcafc4cbd19bd750dee47fc76194a45c181f",
            "rows": int(len(rows)),
            "selected_zHD_positive_rows": int(np.count_nonzero(keep)),
            "all_rows_positive_zHD": bool(np.all(zhd > 0.0)),
            "unique_CIDs": int(len(set(cids))),
            "CID_first": cids[0],
            "CID_last": cids[-1],
            "CID_order_sha256": packed_order_digest,
            "zHD_range": [float(zhd.min()), float(zhd.max())],
            "zHEL_range": [float(zhel.min()), float(zhel.max())],
            "muerr_range": [float(muerr.min()), float(muerr.max())],
            "precision_source": "NPZ cov array, upper triangle float32, symmetrized; used as P=C_STAT+SYS^-1",
            "precision_shape": list(precision.shape),
            "precision_dtype_on_disk": str(packed.dtype),
            "precision_cholesky": "pass",
            "precision_offset_information_1TP1": c_one,
            "muerr_added_separately": False,
        },
        "likelihood": {
            "distance_redshift": "zHD",
            "luminosity_prefactor_redshift": "zHEL",
            "distance_modulus": "5 log10[(1+zHEL) (c/H0_gauge) integral_0^zHD dz/E(z)] +25",
            "offset": "Mhat=-(1^T P d)/(1^T P 1), d=mu_theory-mu_obs",
            "chi2": "(d+Mhat*1)^T P (d+Mhat*1); full dense P, no diagonal approximation",
            "H0_gauge_km_s_Mpc": H0_GAUGE,
            "quadrature": "128-point Gauss-Legendre on [0,zHD] per row",
            "rationale_for_recentering": "subtract d[0] before precision products; the profiled constant absorbs this exactly",
        },
        "fixed_points": fixed,
        "bounded_differential_evolution": fits,
        "zHD_prefactor_error_probe": zhd_probe,
        "runtime_seconds": time.perf_counter() - started,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cpu_affinity_count": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
            "thread_env": {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
            "optimizer": "scipy.optimize.differential_evolution; seeded independently by fit; one worker",
        },
    }
    OUT_FILE.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Test whether the theta-solver H0 ceiling truncates the CAMB BAO profile.

This is a conditional support diagnostic at the valid best point from the
first broad CAMB screen, not a full profile or posterior. It reparameterizes
the same flat-LCDM physical densities directly by H0 and recomputes DESI DR2
BAO with an independent parser/scorer.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import platform
import time

import camb
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "camb_h0_support_test.json"
C_KM_S = 299792.458
OMEGA_B_H2 = 0.060691736587168255
OMEGA_C_H2 = 0.23095938241922825
H0_GRID = [80.0, 90.0, 95.0, 99.0, 99.158774458752, 100.0, 105.0, 110.0, 125.0, 150.0, 200.0]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_data() -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    rows = []
    for line in (ROOT / "context/data/desi_dr2_mean.txt").read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            z_s, mean_s, kind = line.split()
            rows.append((float(z_s), float(mean_s), kind))
    z = np.array([r[0] for r in rows])
    y = np.array([r[1] for r in rows])
    cov = np.loadtxt(ROOT / "context/data/desi_dr2_cov.txt", dtype=np.float64)
    if cov.shape != (len(rows), len(rows)) or not np.array_equal(cov, cov.T):
        raise ValueError("DESI input/covariance shape or symmetry check failed")
    L = np.linalg.cholesky(cov)
    return z, y, cov, [r[2] for r in rows]


def evaluate(h0: float, z: np.ndarray, y: np.ndarray, cov: np.ndarray, kinds: list[str]) -> dict:
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
    pars.set_cosmology(
        H0=float(h0), ombh2=OMEGA_B_H2, omch2=OMEGA_C_H2, omk=0.0,
        mnu=0.06, nnu=3.044, num_massive_neutrinos=1,
        neutrino_hierarchy="degenerate", TCMB=2.7255, YHe=None,
    )
    bg = camb.get_background(pars)
    derived = bg.get_derived_params()
    dm = np.asarray(bg.comoving_radial_distance(z), dtype=np.float64)
    hz = np.asarray(bg.hubble_parameter(z), dtype=np.float64)
    dh = C_KM_S / hz
    pred = np.empty(len(z), dtype=np.float64)
    for i, kind in enumerate(kinds):
        if kind == "DM_over_rs":
            pred[i] = dm[i] / derived["rdrag"]
        elif kind == "DH_over_rs":
            pred[i] = dh[i] / derived["rdrag"]
        elif kind == "DV_over_rs":
            pred[i] = (z[i] * dm[i] ** 2 * dh[i]) ** (1.0 / 3.0) / derived["rdrag"]
        else:
            raise ValueError(kind)
    residual = pred - y
    white = np.linalg.solve(np.linalg.cholesky(cov), residual)
    chi2_chol = float(white @ white)
    chi2_dense = float(residual @ np.linalg.solve(cov, residual))
    if abs(chi2_chol - chi2_dense) > 1e-10:
        raise ArithmeticError("dense and Cholesky covariance scores disagree")
    dist_err = []
    for zi, dmi in zip(z, dm):
        integ, _ = quad(lambda zz: C_KM_S / float(bg.hubble_parameter(zz)), 0.0, float(zi), epsabs=1e-9, epsrel=1e-11, limit=250)
        dist_err.append(abs(float(dmi) - integ))
    return {
        "H0_km_s_Mpc": h0,
        "derived_Omega_m0_including_neutrino": float(pars.omegam),
        "derived_Omega_b0": float(pars.omegab),
        "derived_Omega_c0": float(pars.omegac),
        "r_drag_Mpc": float(derived["rdrag"]),
        "100theta_star_CAMB": float(derived["thetastar"]),
        "chi2_cholesky": chi2_chol,
        "chi2_dense": chi2_dense,
        "max_abs_distance_integral_error_Mpc": float(max(dist_err)),
        "max_abs_bao_residual": float(np.max(np.abs(residual))),
    }


def main() -> None:
    t0 = time.perf_counter()
    z, y, cov, kinds = load_data()
    rows = [evaluate(h0, z, y, cov, kinds) for h0 in H0_GRID]
    conditional = minimize_scalar(
        lambda h0: evaluate(float(h0), z, y, cov, kinds)["chi2_cholesky"],
        bounds=(95.0, 105.0), method="bounded",
        options={"xatol": 1e-7, "maxiter": 80},
    )
    conditional_row = evaluate(float(conditional.x), z, y, cov, kinds)
    result = {
        "status": "complete_conditional_solver_support_diagnostic",
        "interpretation": "Fixed omega_b h2 and omega_c h2 at the previous screen's selected candidate; varies H0 directly, avoiding CAMB's theta-to-H0 root upper bracket of 100. This tests whether the reported point is close to an artificial solver support boundary. Not a re-fit over all parameters, not a posterior/prior result.",
        "parameterization": {"omega_b_h2_fixed": OMEGA_B_H2, "omega_c_h2_fixed": OMEGA_C_H2, "H0_km_s_Mpc_grid": H0_GRID, "flat": True, "w": -1.0, "mnu_eV": 0.06, "N_eff": 3.044, "Tcmb_K": 2.7255, "YHe": "CAMB BBN consistency"},
        "input_hashes": {
            "mean_sha256": sha256(ROOT / "context/data/desi_dr2_mean.txt"),
            "covariance_sha256": sha256(ROOT / "context/data/desi_dr2_cov.txt"),
            "CAMB_wheel_sha256": sha256(ROOT / "context/code/camb-2.0.4-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl"),
        },
        "script_sha256": sha256(Path(__file__)),
        "CAMB_version": camb.__version__,
        "NumPy_version": np.__version__,
        "SciPy_version": scipy.__version__,
        "python": os.sys.version,
        "platform": platform.platform(),
        "grid_results": rows,
        "conditional_H0_minimum": {
            "bounds_km_s_Mpc": [95.0, 105.0],
            "success": bool(conditional.success),
            "message": str(conditional.message),
            "nfev": int(conditional.nfev),
            **conditional_row,
        },
        "duration_seconds": time.perf_counter() - t0,
    }
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "results": rows, "duration_seconds": result["duration_seconds"], "output": str(OUT.relative_to(ROOT))}, sort_keys=True))


if __name__ == "__main__":
    main()

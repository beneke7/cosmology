#!/usr/bin/env python3
"""Independent one-point check of the compute26 CAMB BAO candidate."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PINNED_SITE = ROOT / "work/compute24/site"
if str(PINNED_SITE) not in sys.path:
    sys.path.insert(0, str(PINNED_SITE))
import camb

MEAN_FILE = ROOT / "context/data/desi_dr2_mean.txt"
COV_FILE = ROOT / "context/data/desi_dr2_cov.txt"
REFERENCE_FILE = ROOT / "work/compute26/independent_lcdm_bao_profile.json"
WHEEL_FILE = ROOT / "context/code/camb-2.0.4-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl"
OUT = HERE / "independent_candidate_check.json"
C_KM_S = 299792.458
TOL_CHI2 = 2.0e-8


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_mean(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    parsed = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        z, y, name = raw.split()
        parsed.append((float(z), float(y), name))
    return (
        np.asarray([r[0] for r in parsed], dtype=np.float64),
        np.asarray([r[1] for r in parsed], dtype=np.float64),
        [r[2] for r in parsed],
    )


def main() -> None:
    started = time.perf_counter()
    reference_bytes = REFERENCE_FILE.read_bytes()
    reference = json.loads(reference_bytes)
    candidate = reference["selected_candidate"]
    pars_in = candidate["parameters"]
    h0 = float(pars_in["H0_km_s_Mpc"])
    omega_b = float(pars_in["omega_b_h2"])
    omega_c = float(pars_in["omega_c_h2"])

    z, observed, kinds = read_mean(MEAN_FILE)
    covariance = np.asarray(np.loadtxt(COV_FILE, dtype=np.float64), dtype=np.float64)
    if len(z) != 13 or covariance.shape != (13, 13):
        raise ValueError("DESI DR2 vector/covariance does not have the locked 13-row shape")
    if not np.all(np.isfinite(covariance)) or not np.allclose(covariance, covariance.T, atol=1e-12, rtol=0):
        raise ValueError("DESI covariance symmetry/finiteness failed")

    pars = camb.CAMBparams()
    pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
    pars.set_cosmology(
        H0=h0, ombh2=omega_b, omch2=omega_c, omk=0.0,
        mnu=0.06, nnu=3.044, num_massive_neutrinos=1,
        neutrino_hierarchy="degenerate", TCMB=2.7255, YHe=None,
        tau=0.054,
    )
    background = camb.get_background(pars)
    derived = background.get_derived_params()
    rdrag = float(derived["rdrag"])
    dm = np.asarray(background.comoving_radial_distance(z), dtype=np.float64)
    hz = np.asarray(background.hubble_parameter(z), dtype=np.float64)
    dh = C_KM_S / hz
    pred = np.asarray(
        [
            dm[i] / rdrag if kind == "DM_over_rs" else
            dh[i] / rdrag if kind == "DH_over_rs" else
            (z[i] * dm[i] ** 2 * dh[i]) ** (1.0 / 3.0) / rdrag if kind == "DV_over_rs" else
            np.nan
            for i, kind in enumerate(kinds)
        ],
        dtype=np.float64,
    )
    if not np.all(np.isfinite(pred)):
        raise ValueError("Unknown BAO observable or non-finite prediction")
    residual = pred - observed
    precision_residual = np.linalg.solve(covariance, residual)
    chi2_dense = float(residual @ precision_residual)
    lower = np.linalg.cholesky(covariance)
    whitened = np.linalg.solve(lower, residual)
    chi2_cholesky = float(whitened @ whitened)

    distance_errors = []
    for zi, dmi in zip(z, dm):
        integrated, _ = quad(
            lambda zz: C_KM_S / float(background.hubble_parameter(zz)),
            0.0, float(zi), epsabs=1e-9, epsrel=1e-11, limit=250,
        )
        distance_errors.append(abs(float(dmi) - integrated))

    claimed = float(candidate["bao"]["chi2_full_covariance_cholesky"])
    delta = abs(chi2_cholesky - claimed)
    result = {
        "status": "complete_independent_one_point_score_and_distance_audit",
        "implementation_note": "Standalone parser, CAMB setup, observable construction, dense precision solve, Cholesky whitening, and c/H quadrature; imports no compute26 scoring code.",
        "reference_candidate_parameters": {"H0_km_s_Mpc": h0, "omega_b_h2": omega_b, "omega_c_h2": omega_c},
        "recomputed_derived": {
            "Omega_m0_including_massive_nu": float(pars.omegam),
            "Omega_b0": float(pars.omegab),
            "Omega_c0": float(pars.omegac),
            "r_drag_Mpc": rdrag,
            "100theta_star_CAMB": float(derived["thetastar"]),
        },
        "reference_chi2": claimed,
        "chi2_dense_precision": chi2_dense,
        "chi2_cholesky_whitening": chi2_cholesky,
        "abs_delta_to_reference": delta,
        "chi2_match_pass": bool(delta < TOL_CHI2),
        "dense_vs_cholesky_abs_delta": abs(chi2_dense - chi2_cholesky),
        "max_abs_CAMB_distance_minus_c_over_H_integral_Mpc": float(max(distance_errors)),
        "all_DM_integrals_below_1e-6_Mpc": bool(max(distance_errors) < 1e-6),
        "predictions_in_data_order": pred.tolist(),
        "residuals_in_data_order": residual.tolist(),
        "whitened_residuals": whitened.tolist(),
        "source_hashes": {
            "mean": sha256(MEAN_FILE),
            "covariance": sha256(COV_FILE),
            "compute26_reference": sha256(REFERENCE_FILE),
            "CAMB_wheel": sha256(WHEEL_FILE),
        },
        "versions": {
            "CAMB": camb.__version__,
            "NumPy": np.__version__,
            "SciPy": scipy.__version__,
            "Python": sys.version,
            "platform": platform.platform(),
        },
        "runtime_seconds": time.perf_counter() - started,
    }
    result["script_sha256"] = sha256(Path(__file__))
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("status", "reference_candidate_parameters", "reference_chi2", "chi2_dense_precision", "chi2_cholesky_whitening", "abs_delta_to_reference", "chi2_match_pass", "max_abs_CAMB_distance_minus_c_over_H_integral_Mpc", "runtime_seconds")}, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent full-covariance point-score audit of work/compute31.

This deliberately does not import compute31's scorer or any earlier BAO
likelihood module. It reads the released DESI vector/covariance and reruns
CAMB for each saved point, then forms the distance observables and solves the
covariance system directly.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "work/compute24/site"))
import camb  # noqa: E402

MEAN = ROOT / "context/data/desi_dr2_mean.txt"
COV = ROOT / "context/data/desi_dr2_cov.txt"
SOURCE = ROOT / "work/compute31/bbn_bao_profile.py"
RESULT = ROOT / "work/compute31/bbn_bao_profile.json"
OUTPUT = Path(__file__).with_name("independent_bbn_bao_score_audit.json")
EXPECTED_NAMES = [
    "DV_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs",
    "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs",
    "DH_over_rs", "DH_over_rs", "DM_over_rs",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_inputs() -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray]:
    rows = []
    for line in MEAN.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if text and not text.startswith("#"):
            z_text, value_text, name = text.split()
            rows.append((float(z_text), float(value_text), name))
    redshifts = np.array([row[0] for row in rows], dtype=np.float64)
    observed = np.array([row[1] for row in rows], dtype=np.float64)
    names = [row[2] for row in rows]
    covariance = np.loadtxt(COV, dtype=np.float64)
    if redshifts.shape != (13,) or observed.shape != (13,) or names != EXPECTED_NAMES:
        raise AssertionError("released DR2 BAO vector shape or observable order changed")
    if covariance.shape != (13, 13) or not np.allclose(covariance, covariance.T, atol=1e-12, rtol=0):
        raise AssertionError("released DR2 covariance shape/symmetry check failed")
    np.linalg.cholesky(covariance)
    return redshifts, observed, names, covariance


def predict(point: dict, redshifts: np.ndarray, names: list[str]) -> dict:
    pars = camb.CAMBparams()
    pars.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model="fluid")
    pars.set_cosmology(
        H0=float(point["H0_km_s_Mpc"]),
        ombh2=float(point["omega_b_h2"]),
        omch2=float(point["omega_c_h2"]),
        omk=0.0,
        mnu=0.06,
        nnu=3.044,
        num_massive_neutrinos=1,
        neutrino_hierarchy="degenerate",
        TCMB=2.7255,
        YHe=None,
    )
    omega_m0 = float(pars.omegam)
    if not math.isfinite(omega_m0) or not 0.0 < omega_m0 < 1.0:
        raise AssertionError(f"non-positive flat Lambda closure Omega_m0={omega_m0}")
    background = camb.get_background(pars)
    omega_de0 = float(background.get_Omega("de", 0.0))
    if not math.isfinite(omega_de0) or omega_de0 <= 0.0:
        raise AssertionError(f"non-positive Lambda closure Omega_de0={omega_de0}")
    rd = float(background.get_derived_params()["rdrag"])
    dm = np.asarray(background.comoving_radial_distance(redshifts), dtype=np.float64)
    dh = 299792.458 / np.asarray(background.hubble_parameter(redshifts), dtype=np.float64)
    predicted = []
    for z, d_m, d_h, name in zip(redshifts, dm, dh, names):
        if name == "DV_over_rs":
            value = (z * d_m * d_m * d_h) ** (1.0 / 3.0) / rd
        elif name == "DM_over_rs":
            value = d_m / rd
        elif name == "DH_over_rs":
            value = d_h / rd
        else:
            raise AssertionError(f"unrecognized observable {name}")
        predicted.append(value)
    return {
        "prediction": np.asarray(predicted, dtype=np.float64),
        "r_drag_Mpc": rd,
        "Omega_m0": omega_m0,
        "Omega_de0": omega_de0,
        "theta_star": float(background.get_derived_params()["thetastar"]),
    }


def main() -> None:
    started = time.perf_counter()
    z, observed, names, covariance = load_inputs()
    cholesky = np.linalg.cholesky(covariance)
    record = json.loads(RESULT.read_text(encoding="utf-8"))
    points = [("pure_BAO_best_valid_fit", record["pure_BAO_best_valid_fit"]["profile_record"]["best_valid_trial"],
               record["pure_BAO_best_valid_fit"]["profile_record"]["selected_score_checks"])]
    for item in record["BBN_prior_sensitivity_fits"]:
        points.append((item["prior"]["id"], item["best_valid_trial"], item["selected_score_checks"]))

    audits = []
    for label, point, saved in points:
        calc = predict(point, z, names)
        residual = calc["prediction"] - observed
        dense_score = float(residual @ np.linalg.solve(covariance, residual))
        white = np.linalg.solve(cholesky, residual)
        whitened_score = float(white @ white)
        delta_prediction = float(np.max(np.abs(calc["prediction"] - np.asarray(saved["predictions_in_released_order"]))))
        saved_score = float(saved["chi2_dense_solve"])
        score_delta = abs(dense_score - saved_score)
        if delta_prediction > 2e-12 or score_delta > 2e-10 or abs(dense_score - whitened_score) > 2e-10:
            raise AssertionError(f"independent point audit failed for {label}")
        audits.append({
            "label": label,
            "point": {key: point[key] for key in ("H0_km_s_Mpc", "omega_b_h2", "omega_c_h2")},
            "independent_chi2_dense": dense_score,
            "independent_chi2_cholesky": whitened_score,
            "saved_chi2": saved_score,
            "absolute_chi2_delta": score_delta,
            "maximum_prediction_delta": delta_prediction,
            "r_drag_Mpc": calc["r_drag_Mpc"],
            "Omega_m0": calc["Omega_m0"],
            "Omega_de0": calc["Omega_de0"],
            "100theta_star_CAMB": calc["theta_star"],
        })

    output = {
        "status": "independent_one_point_full_covariance_score_audit_pass",
        "scope": "Four saved CAMB points re-evaluated from DESI DR2 mean/covariance; not an independent optimizer/globality audit.",
        "implementation": "Standalone parser, CAMB prediction assembly, dense covariance solve and Cholesky whitening; imports no compute31 or other BAO scoring module.",
        "versions": {"Python": platform.python_version(), "NumPy": np.__version__, "CAMB": getattr(camb, "__version__", "unknown")},
        "elapsed_seconds": time.perf_counter() - started,
        "tolerances": {"max_prediction_abs": 2e-12, "max_chi2_abs": 2e-10},
        "audits": audits,
        "hashes": {
            "mean_vector_sha256": sha256(MEAN),
            "covariance_sha256": sha256(COV),
            "compute31_script_sha256": sha256(SOURCE),
            "compute31_result_sha256": sha256(RESULT),
            "this_script_sha256": sha256(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

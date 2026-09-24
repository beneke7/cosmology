#!/usr/bin/env python3
"""Independent scalar reference for the radiation sensitivity lane.

No campaign prediction, radiation-worker, or likelihood code is imported.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import scipy
from scipy import constants
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[2]
MEAN = ROOT / "context/data/desi_dr2_mean.txt"
COV = ROOT / "context/data/desi_dr2_cov.txt"
SEED = ROOT / "experiments/campaign_seed_bao/result.json"
T_CMB = 2.7255
N_EFF = 3.046


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physical_radiation_density() -> dict:
    """Derive photon mass density directly from the blackbody energy density."""
    a_rad = (8 * math.pi**5 * constants.k**4
             / (15 * constants.h**3 * constants.c**3))
    rho_gamma = a_rad * T_CMB**4 / constants.c**2
    # Use the stated IAU angular convention explicitly and independently verify
    # the parsec exported by scipy.constants, avoiding an AU/tan(arcsec) mix-up.
    iau_parsec = constants.au * 648_000 / math.pi
    hubble100 = 100_000 / (1_000_000 * iau_parsec)
    rho_critical100 = 3 * hubble100**2 / (8 * math.pi * constants.G)
    omega_gamma = rho_gamma / rho_critical100
    neutrino_factor = (7 / 8) * (4 / 11)**(4 / 3) * N_EFF
    return {
        "T_CMB_K": T_CMB,
        "N_eff_massless_reference": N_EFF,
        "constants_from_scipy_SI": {
            "k_B": constants.k,
            "h_P": constants.h,
            "c": constants.c,
            "G": constants.G,
            "AU": constants.au,
            "parsec_scipy_reference_only": constants.parsec,
        },
        "parsec_used_exact_IAU_m": iau_parsec,
        "relative_difference_scipy_parsec_from_IAU": constants.parsec / iau_parsec - 1,
        "radiation_energy_density_constant_J_m3_K4": a_rad,
        "photon_mass_density_kg_m3": rho_gamma,
        "H100_per_s": hubble100,
        "rho_critical_h2_kg_m3": rho_critical100,
        "omega_gamma": omega_gamma,
        "omega_r": omega_gamma * (1 + neutrino_factor),
    }


def prediction(zs, labels, parameters, omega_r):
    alpha, om, w0, wa = parameters
    ode = 1 - om - omega_r
    assert 0 <= omega_r < 1 - om

    def e2(z):
        zp = 1 + z
        f = zp**(3 * (1 + w0 + wa)) * math.exp(-3 * wa * z / zp)
        return om * zp**3 + omega_r * zp**4 + ode * f

    assert abs(e2(0) - 1) < 4e-16
    output = []
    integration_errors = []
    for z, label in zip(zs, labels):
        integral, error = quad(lambda x: 1 / math.sqrt(e2(x)), 0, float(z),
                               epsabs=1e-12, epsrel=1e-12)
        dm = alpha * integral
        dh = alpha / math.sqrt(e2(z))
        values = {"DM_over_rs": dm, "DH_over_rs": dh,
                  "DV_over_rs": (z * dm * dm * dh)**(1 / 3)}
        output.append(values[label])
        integration_errors.append(error)
    return np.array(output), max(integration_errors)


def score(data, covariance, model):
    residual = data - model
    return float(residual @ np.linalg.solve(covariance, residual))


def displacement(delta, covariance):
    norm2 = float(delta @ np.linalg.solve(covariance, delta))
    return {
        "delta_prediction": delta.tolist(),
        "max_abs_shift_over_marginal_sigma": float(np.max(np.abs(delta) / np.sqrt(np.diag(covariance)))),
        "full_covariance_displacement_squared": norm2,
        "full_covariance_displacement_norm": math.sqrt(max(norm2, 0)),
    }


def main():
    started = time.perf_counter()
    numeric = np.loadtxt(MEAN, usecols=(0, 1))
    labels = np.loadtxt(MEAN, usecols=(2,), dtype=str)
    zs, observed = numeric.T
    covariance = np.loadtxt(COV)
    assert len(observed) == 13 and covariance.shape == (13, 13)
    np.linalg.cholesky(covariance)
    assert labels[-2:].tolist() == ["DH_over_rs", "DM_over_rs"]
    radiation = physical_radiation_density()
    baseline = json.loads(SEED.read_text())
    records = []
    maximum_baseline_prediction_error = 0.0
    maximum_baseline_score_error = 0.0
    max_profile_error = 0.0
    for fit in baseline["models"]:
        parameters = [fit[key] for key in ["alpha", "Omega_m", "w0", "wa"]]
        no_rad, quadrature_error = prediction(zs, labels, parameters, 0)
        chi0 = score(observed, covariance, no_rad)
        maximum_baseline_prediction_error = max(maximum_baseline_prediction_error,
            float(np.max(np.abs(no_rad - fit["prediction"]))))
        maximum_baseline_score_error = max(maximum_baseline_score_error,
            abs(chi0 - fit["chi2"]))
        unit_template = no_rad / parameters[0]
        best_alpha = float(unit_template @ np.linalg.solve(covariance, observed)
                           / (unit_template @ np.linalg.solve(covariance, unit_template)))
        best_alpha = float(np.clip(best_alpha, 1e-6, 1e4))
        max_profile_error = max(max_profile_error, abs(best_alpha - parameters[0]))
        fixed = []
        for h0 in [50, 55, 60, 65, 67.4, 70, 75, 80, 85, 90]:
            omega_r = radiation["omega_r"] / (h0 / 100)**2
            model, error = prediction(zs, labels, parameters, omega_r)
            delta = model - no_rad
            chisq = score(observed, covariance, model)
            linear = float(-2 * (observed - no_rad) @ np.linalg.solve(covariance, delta))
            record = {"H0_km_s_Mpc": h0, "Omega_r": omega_r,
                      "chi2_fixed": chisq, "delta_chi2_fixed": chisq - chi0,
                      "first_order_residual_term": linear,
                      "max_abs_fractional_prediction_shift": float(np.max(np.abs(delta / no_rad))),
                      "quad_reported_max_abs_error": error,
                      **displacement(delta, covariance)}
            assert abs(record["delta_chi2_fixed"] - linear
                       - record["full_covariance_displacement_squared"]) < 1e-10
            fixed.append(record)
        records.append({"model": fit["model"], "parameters": parameters,
                        "no_radiation_chi2": chi0, "prediction": no_rad.tolist(),
                        "profile_alpha_at_seed_shape": best_alpha,
                        "quad_reported_max_abs_error": quadrature_error,
                        "fixed_baseline_radiation_scan": fixed})
    assert maximum_baseline_prediction_error < 1e-10
    assert maximum_baseline_score_error < 1e-9
    result = {
        "status": "independently_checked_fixed_baseline",
        "scope": "Separate scalar reference; sensitivity domain, no H0 inference, no posterior/evidence/significance.",
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "command": "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ./.venv/bin/python scripts/run_bounded.py --seconds 120 -- ./.venv/bin/python work/orchestrator2_radiation/scalar_reference.py",
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__},
        "sha256": {str(path.relative_to(ROOT)): digest(path) for path in [MEAN, COV, SEED, Path(__file__)]},
        "constants": radiation,
        "checks": {"max_seed_prediction_error": maximum_baseline_prediction_error,
                   "max_seed_chi2_error": maximum_baseline_score_error,
                   "max_seed_alpha_difference_from_exact_profile": max_profile_error,
                   "quadratic_score_identity_abs_tolerance": 1e-10,
                   "full_covariance_in_release_order": True},
        "records": records,
        "runtime_seconds": time.perf_counter() - started,
    }
    output = Path(__file__).with_name("scalar_reference.json")
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(output.relative_to(ROOT)), "checks": result["checks"],
                      "omega_gamma": radiation["omega_gamma"], "omega_r": radiation["omega_r"],
                      "runtime_seconds": result["runtime_seconds"],
                      "H0_50_fixed_corrections": [{"model": item["model"], **item["fixed_baseline_radiation_scan"][0]}
                                                for item in records]}, indent=2))


if __name__ == "__main__":
    main()

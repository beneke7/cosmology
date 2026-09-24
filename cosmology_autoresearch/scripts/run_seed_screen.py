#!/usr/bin/env python3
"""Reproducible BAO-only seed run; no posterior, evidence or discovery inference."""
import argparse
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad

from background_bao import BAOData, FIT_BOUNDS, fit_model, predict_bao


def independent_prediction(data, pars, omega_r=0.0):
    """Scalar adaptive integration; independent of the vectorized kernel."""
    alpha, om, w0, wa = pars
    ode = 1 - om - omega_r
    def expansion(z):
        a_inv = 1 + z
        de = math.exp(3 * ((1 + w0 + wa) * math.log(a_inv) - wa * z / a_inv))
        return math.sqrt(om * a_inv**3 + omega_r * a_inv**4 + ode * de)
    pred = []
    for z, name in zip(data.z, data.observable):
        dm = alpha * quad(lambda x: 1 / expansion(x), 0, float(z), epsabs=1e-12, epsrel=1e-12)[0]
        dh = alpha / expansion(float(z))
        pred.append({"DM_over_rs": dm, "DH_over_rs": dh,
                     "DV_over_rs": (float(z) * dm * dm * dh)**(1 / 3)}[name])
    return np.asarray(pred)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mean", default="context/data/desi_dr2_mean.txt")
    parser.add_argument("--cov", default="context/data/desi_dr2_cov.txt")
    parser.add_argument("--output", default="experiments/seed_bao/result.json")
    args = parser.parse_args()
    data = BAOData.from_files(args.mean, args.cov)
    started = time.perf_counter()
    result = {"status": "exploratory_screen", "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "scope": "13-row compressed DESI DR2 BAO only; free c/(H0 rd), flat late-time background, no radiation in fit; no posterior/evidence/significance claim",
              "data_sha256": {args.mean: digest(args.mean), args.cov: digest(args.cov)},
              "code_sha256": {str(p): digest(p) for p in [Path(__file__), Path(__file__).with_name("background_bao.py")]},
              "mean_row_order": [{"z": float(z), "observable": str(o)} for z, o in zip(data.z, data.observable)],
              "models": []}
    for model in ("lcdm", "wcdm", "cpl"):
        fit = fit_model(data, model, starts=16)
        pars = [fit[key] for key in ["alpha", "Omega_m", "w0", "wa"]]
        vector = predict_bao(data.z, data.observable, *pars)
        independent = independent_prediction(data, pars)
        np.testing.assert_allclose(vector, independent, rtol=1e-11, atol=1e-11)
        residual = data.value - independent
        chi2_independent = float(residual @ np.linalg.solve(data.covariance, residual))
        assert abs(chi2_independent - fit["chi2"]) < 1e-7
        # Sensitivity exercise at fixed best-fit parameters, NOT a radiation refit.
        omega_r = 9e-5
        radiation = independent_prediction(data, pars, omega_r=omega_r)
        rr = data.value - radiation
        radiation_chi2 = float(rr @ np.linalg.solve(data.covariance, rr))
        fit.update({"parameter_bounds": FIT_BOUNDS[model], "starts": 16, "optimizer_seed": 20260923,
                    "independent_max_prediction_difference": float(np.max(np.abs(vector - independent))),
                    "independent_chi2": chi2_independent,
                    "prediction": vector.tolist(),
                    "radiation_sensitivity": {"Omega_r": omega_r, "refitted": False,
                      "assumption": "constant present-day radiation parameter, flatness restored by reducing Omega_de; illustrative only",
                      "delta_chi2_at_fixed_parameters": radiation_chi2 - fit["chi2"],
                      "max_shift_over_marginal_measurement_sigma": float(np.max(np.abs(radiation-vector) / np.sqrt(np.diag(data.covariance))))}})
        result["models"].append(fit)
    if not all(m["success"] for m in result["models"]):
        result["status"] = "optimizer_failure_needs_review"
    result["runtime_seconds"] = time.perf_counter() - started
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(output), "status": result["status"],
                      "models": [{key: model[key] for key in ["model", "chi2", "success"]} for model in result["models"]]}, indent=2))
    return 0 if result["status"] == "exploratory_screen" else 2


if __name__ == "__main__":
    raise SystemExit(main())

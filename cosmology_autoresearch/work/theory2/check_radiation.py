#!/usr/bin/env python3
"""Independent float64 check of the late-time radiation omission screen.

This script does not import the project BAO background implementation.  It
evaluates the flat LCDM/wCDM/CPL distance equations independently at the
campaign-seed fit points, treating alpha = c/(H0 rd) as an independent fixed
amplitude while adding a small massless-radiation reference term.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "radiation_check.json"
MEAN = ROOT / "context/data/desi_dr2_mean.txt"
COV = ROOT / "context/data/desi_dr2_cov.txt"
SEED = ROOT / "experiments/campaign_seed_bao/result.json"
STARTER = ROOT / "scripts/background_bao.py"

# CODATA 2022 SI constants. c, h_P and k_B are exact in the current SI;
# G carries the quoted 2022 CODATA relative standard uncertainty 2.2e-5.
C = np.float64(299_792_458.0)  # m s^-1
H_PLANCK = np.float64(6.626_070_15e-34)  # J s
K_B = np.float64(1.380_649e-23)  # J K^-1
G_NEWTON = np.float64(6.674_30e-11)  # m^3 kg^-1 s^-2
AU_M = np.float64(149_597_870_700.0)  # exact IAU conventional length
# IAU definition: 1 pc = (648000/pi) au exactly.
MPC_M = np.float64(1e6) * AU_M * (648_000.0 / np.pi)
T_CMB_K = np.float64(2.7255)
N_EFF = np.float64(3.046)
H100_SI = np.float64(100_000.0) / MPC_M


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_mean(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z, y, obs = [], [], []
    for line in path.read_text(encoding="utf-8").splitlines():
        row = line.strip()
        if not row or row.startswith("#"):
            continue
        a, b, c = row.split()
        z.append(float(a))
        y.append(float(b))
        obs.append(c)
    return (np.asarray(z, dtype=np.float64), np.asarray(y, dtype=np.float64),
            np.asarray(obs, dtype=str))


def de_factor(z: np.ndarray, w0: float, wa: float) -> np.ndarray:
    """rho_DE(z)/rho_DE(0) for CPL; wa=0 gives constant-w."""
    x = 1.0 + z
    log_f = 3.0 * (1.0 + w0 + wa) * np.log1p(z) - 3.0 * wa * z / x
    return np.exp(log_f)


def predictions_and_derivative(z: np.ndarray, obs: np.ndarray, alpha: float,
                               om: float, w0: float, wa: float,
                               omega_r: float, n_quad: int) -> tuple[np.ndarray, np.ndarray]:
    """Return predictions and analytic d(predictions)/d(Omega_r) at fixed shape."""
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    qz = 0.5 * z[:, None] * (nodes[None, :] + 1.0)
    qx = 1.0 + qz
    f_q = de_factor(qz, w0, wa)
    e2_q = om * qx**3 + omega_r * qx**4 + (1.0 - om - omega_r) * f_q
    if np.any(e2_q <= 0.0) or not np.all(np.isfinite(e2_q)):
        raise ValueError("non-positive or non-finite E(z)^2")
    q_term = qx**4 - f_q  # derivative of E^2 after enforcing flat closure
    integ = 0.5 * z * np.sum(weights[None, :] / np.sqrt(e2_q), axis=1,
                              dtype=np.float64)
    d_integ = 0.5 * z * np.sum(-0.5 * weights[None, :] * q_term / e2_q**1.5,
                               axis=1, dtype=np.float64)

    x = 1.0 + z
    f_z = de_factor(z, w0, wa)
    e2_z = om * x**3 + omega_r * x**4 + (1.0 - om - omega_r) * f_z
    q_z = x**4 - f_z
    dm = alpha * integ
    dh = alpha / np.sqrt(e2_z)
    dv = np.cbrt(z * dm**2 * dh)
    d_dm = alpha * d_integ
    d_dh = -0.5 * alpha * q_z / e2_z**1.5
    d_dv = dv * ((2.0 / 3.0) * d_integ / integ - q_z / (6.0 * e2_z))

    prediction = np.empty_like(z, dtype=np.float64)
    derivative = np.empty_like(z, dtype=np.float64)
    selectors = {
        "DM_over_rs": (dm, d_dm),
        "DH_over_rs": (dh, d_dh),
        "DV_over_rs": (dv, d_dv),
    }
    for label, (value, deriv) in selectors.items():
        mask = obs == label
        prediction[mask] = value[mask]
        derivative[mask] = deriv[mask]
    return prediction, derivative


def chi2(residual: np.ndarray, covariance: np.ndarray) -> float:
    """Full-covariance quadratic form, via a CPU float64 solve."""
    return float(residual @ np.linalg.solve(covariance, residual))


def main() -> None:
    # Blackbody energy density from integrating Planck's spectrum:
    # u_gamma = a_rad T^4, a_rad = 8 pi^5 k_B^4/(15 h_P^3 c^3) = 4 sigma_SB/c.
    a_rad = 8.0 * np.pi**5 * K_B**4 / (15.0 * H_PLANCK**3 * C**3)
    sigma_sb = 2.0 * np.pi**5 * K_B**4 / (15.0 * H_PLANCK**3 * C**2)
    u_gamma = a_rad * T_CMB_K**4
    rho_gamma_mass = u_gamma / C**2
    rho_crit_h1 = 3.0 * H100_SI**2 / (8.0 * np.pi * G_NEWTON)
    omega_gamma_h2 = rho_gamma_mass / rho_crit_h1
    nu_factor_per_neff = (7.0 / 8.0) * (4.0 / 11.0)**(4.0 / 3.0)
    omega_r_h2 = omega_gamma_h2 * (1.0 + nu_factor_per_neff * N_EFF)

    z, y, obs = read_mean(MEAN)
    cov = np.loadtxt(COV, dtype=np.float64)
    seed = json.loads(SEED.read_text(encoding="utf-8"))
    sigmas = np.sqrt(np.diag(cov))
    seed_models = {row["model"]: row for row in seed["models"]}
    models = []
    for name in ("lcdm", "wcdm", "cpl"):
        row = seed_models[name]
        alpha, om = float(row["alpha"]), float(row["Omega_m"])
        w0, wa = float(row["w0"]), float(row["wa"])
        base, derivative = predictions_and_derivative(z, obs, alpha, om, w0, wa, 0.0, 512)
        # Independent reproduction of the zero-radiation campaign seed predictions.
        source_predictions = np.asarray(row["prediction"], dtype=np.float64)
        model_record = {
            "model": name,
            "reference_parameters": {"alpha": alpha, "Omega_m": om, "w0": w0, "wa": wa},
            "independent_zero_radiation_max_abs_prediction_difference": float(np.max(np.abs(base-source_predictions))),
            "results_by_H0_screen_endpoint": [],
        }
        for H0 in (50.0, 90.0):
            h_dimless = H0 / 100.0
            omega_r = float(omega_r_h2 / h_dimless**2)
            first_order_delta = derivative * omega_r
            rad_exact, _ = predictions_and_derivative(z, obs, alpha, om, w0, wa, omega_r, 512)
            exact_delta = rad_exact - base
            linear_error = exact_delta - first_order_delta
            baseline_chi2 = chi2(y - base, cov)
            rad_chi2 = chi2(y - rad_exact, cov)
            shift_precision2 = chi2(exact_delta, cov)
            model_record["results_by_H0_screen_endpoint"].append({
                "H0_km_s_Mpc": H0,
                "h_dimensionless_screen_value": h_dimless,
                "Omega_r_today": omega_r,
                "first_order_prediction_shift": first_order_delta.tolist(),
                "exact_prediction_shift": exact_delta.tolist(),
                "exact_shift_over_marginal_row_sigma": (exact_delta / sigmas).tolist(),
                "max_abs_shift_over_marginal_row_sigma": float(np.max(np.abs(exact_delta) / sigmas)),
                "full_covariance_shift_mahalanobis2": shift_precision2,
                "full_covariance_shift_mahalanobis": float(np.sqrt(shift_precision2)),
                "baseline_chi2_recomputed": baseline_chi2,
                "chi2_with_radiation_at_same_fit_parameters": rad_chi2,
                "delta_chi2_at_same_fit_parameters": rad_chi2 - baseline_chi2,
                "max_abs_exact_minus_first_order_over_marginal_sigma": float(np.max(np.abs(linear_error) / sigmas)),
                "first_order_relative_l2_error_in_prediction_shift": float(
                    np.linalg.norm(linear_error) / np.linalg.norm(exact_delta)
                ),
            })
        grid_summary = []
        for H0 in np.linspace(50.0, 90.0, 41):
            omega_r = float(omega_r_h2 / (H0 / 100.0)**2)
            rad_exact, _ = predictions_and_derivative(z, obs, alpha, om, w0, wa, omega_r, 512)
            delta = rad_exact - base
            grid_summary.append({
                "H0_km_s_Mpc": float(H0),
                "max_abs_shift_over_marginal_row_sigma": float(np.max(np.abs(delta) / sigmas)),
                "full_covariance_shift_mahalanobis": float(np.sqrt(chi2(delta, cov))),
                "abs_delta_chi2_at_same_fit_parameters": abs(chi2(y-rad_exact, cov) - chi2(y-base, cov)),
            })
        model_record["H0_screen_41_point_grid_summary"] = {
            "grid_endpoints_km_s_Mpc": [50.0, 90.0],
            "point_count": len(grid_summary),
            "max_marginal_sigma_shift": max(grid_summary, key=lambda x: x["max_abs_shift_over_marginal_row_sigma"]),
            "max_full_covariance_displacement": max(grid_summary, key=lambda x: x["full_covariance_shift_mahalanobis"]),
            "max_abs_delta_chi2_at_same_fit_parameters": max(grid_summary, key=lambda x: x["abs_delta_chi2_at_same_fit_parameters"]),
        }
        models.append(model_record)

    # Independent quadrature convergence check for the largest radiation endpoint.
    convergence = {}
    for name in ("lcdm", "wcdm", "cpl"):
        row = seed_models[name]
        H0 = 50.0
        omega_r = float(omega_r_h2 / (H0 / 100.0)**2)
        predictions = {}
        for n_quad in (64, 128, 256, 512):
            predictions[n_quad], _ = predictions_and_derivative(
                z, obs, float(row["alpha"]), float(row["Omega_m"]),
                float(row["w0"]), float(row["wa"]), omega_r, n_quad)
        convergence[name] = {
            "max_abs_prediction_difference_64_vs_512": float(np.max(np.abs(predictions[64]-predictions[512]))),
            "max_abs_prediction_difference_128_vs_512": float(np.max(np.abs(predictions[128]-predictions[512]))),
            "max_abs_prediction_difference_256_vs_512": float(np.max(np.abs(predictions[256]-predictions[512]))),
        }

    source_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in (MEAN, COV, STARTER, SEED)}
    output = {
        "status": "independent_check; exploratory fixed-shape/fixed-amplitude screen",
        "purpose": "Quantify the omitted radiation term against the 13-row supplied DESI DR2 BAO covariance for flat LCDM, wCDM and CPL at campaign-seed fit points.",
        "scientific_scope": {
            "radiation_reference": "photons plus effectively massless radiation normalized with N_eff=3.046 relative to the nominal (4/11)^(1/3) neutrino-to-photon temperature ratio; N_eff is the effective-density parameterization, not a species-by-species thermal-history calculation",
            "neutrino_caveat": "This massless-neutrino reference extension does not reproduce the DESI massive-neutrino density transition; Omega_m and fitted shape parameters are held fixed while flat closure is restored through Omega_de.",
            "amplitude": "alpha=c/(H0 rd) remains a free BAO amplitude and is held fixed at each seed fit point; no H0 or rd inference is made.",
            "H0_domain": "50-90 km s^-1 Mpc^-1 is a sensitivity screen only, not a prior or inference.",
            "fit_points": "Supplied experiments/campaign_seed_bao/result.json exploratory optimizer fit points; models are not re-fit after radiation is added.",
            "precision_metrics": "Per-row shifts use sqrt(diag(C)); full-covariance displacement is delta_prediction^T C^-1 delta_prediction. Delta-chi2 compares the data residual at the same fit parameters, before and after radiation.",
        },
        "constants_and_physical_density": {
            "T_CMB_K_input": float(T_CMB_K),
            "h_P_J_s": float(H_PLANCK),
            "k_B_J_K": float(K_B),
            "c_m_s": float(C),
            "G_m3_kg_s2": float(G_NEWTON),
            "G_relative_standard_uncertainty_CODATA_2022": 2.2e-5,
            "AU_m_exact_IAU": float(AU_M),
            "parsec_m_from_pc_definition": float(MPC_M/1e6),
            "Mpc_m": float(MPC_M),
            "H100_s_inv": float(H100_SI),
            "a_rad_J_m3_K4_from_h_k_c": float(a_rad),
            "sigma_SB_4sigma_over_c_J_m2_s_K4": float(sigma_sb),
            "relative_a_rad_sigma_identity_residual": float(abs(a_rad-4.0*sigma_sb/C)/a_rad),
            "photon_energy_density_J_m3_at_input_temperature": float(u_gamma),
            "photon_mass_density_kg_m3_at_input_temperature": float(rho_gamma_mass),
            "critical_mass_density_kg_m3_at_h_equals_1": float(rho_crit_h1),
            "Omega_gamma_h2": float(omega_gamma_h2),
            "N_eff": float(N_EFF),
            "neutrino_to_photon_density_factor_per_Neff": float(nu_factor_per_neff),
            "Omega_r_h2_massless_reference": float(omega_r_h2),
            "Omega_r_screen_endpoints": {
                "H0_50": float(omega_r_h2/(0.5**2)),
                "H0_90": float(omega_r_h2/(0.9**2)),
            },
            "formulae": {
                "u_gamma": "(8*pi^5/15) (k_B^4/(h_P^3*c^3)) T_CMB^4 = a_rad*T_CMB^4",
                "Omega_gamma_h2": "(8*pi*G/(3*c^2*H100^2))*u_gamma, with H100=100 km s^-1 Mpc^-1",
                "Omega_r_h2": "Omega_gamma_h2 * [1+(7/8)*(4/11)^(4/3)*N_eff]",
                "Omega_r_H0": "Omega_r_h2 / h^2, h=H0/(100 km s^-1 Mpc^-1)",
            },
        },
        "input_sha256": source_hashes,
        "script_sha256": sha256(Path(__file__).resolve()),
        "execution": {
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_only": True,
            "dtype": "float64",
            "blas_thread_environment": {key: os.environ.get(key) for key in (
                "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")},
            "quadrature": "Independent Gauss-Legendre integration, 512-point production values; convergence compared at 64/128/256 points.",
        },
        "quadrature_convergence_at_H0_50": convergence,
        "models": models,
    }
    OUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUT),
        "Omega_gamma_h2": float(omega_gamma_h2),
        "Omega_r_h2": float(omega_r_h2),
        "Omega_r_H0_50": output["constants_and_physical_density"]["Omega_r_screen_endpoints"]["H0_50"],
        "Omega_r_H0_90": output["constants_and_physical_density"]["Omega_r_screen_endpoints"]["H0_90"],
        "models": [{
            "model": m["model"],
            "max_sigma_at_H0_50": m["results_by_H0_screen_endpoint"][0]["max_abs_shift_over_marginal_row_sigma"],
            "mahalanobis_at_H0_50": m["results_by_H0_screen_endpoint"][0]["full_covariance_shift_mahalanobis"],
            "delta_chi2_at_H0_50": m["results_by_H0_screen_endpoint"][0]["delta_chi2_at_same_fit_parameters"],
            "max_sigma_at_H0_90": m["results_by_H0_screen_endpoint"][1]["max_abs_shift_over_marginal_row_sigma"],
            "mahalanobis_at_H0_90": m["results_by_H0_screen_endpoint"][1]["full_covariance_shift_mahalanobis"],
            "delta_chi2_at_H0_90": m["results_by_H0_screen_endpoint"][1]["delta_chi2_at_same_fit_parameters"],
        } for m in models],
    }, indent=2))


if __name__ == "__main__":
    main()

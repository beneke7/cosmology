#!/usr/bin/env python3
"""K=0 Newtonian circular-speed baseline from the resolved NGC 3198 profiles.

The axisymmetric Hankel integral is independently checked against the analytic
Freeman exponential disk. This is a baryon-source and forward-force baseline,
not a TG fit or a halo/MOND comparison.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import platform
import time
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator
from scipy.special import iv, j0, j1, kv


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PROFILE_PATH = OUT / "source_profiles.csv"
ROTMOD_ZIP = ROOT / "context/data/Rotmod_LTG.zip"
ROTMOD_MEMBER = "NGC3198_rotmod.dat"
G = 4.30091e-6  # kpc (km/s)^2 / Msun
R_SRC_MAX = 50.0  # kpc; all maps end by 35 kpc, stellar tail is explicit to 50
R_SRC_STEP = 0.01
K_MAX = 80.0  # kpc^-1; safely below the radial source-grid Nyquist frequency
N_K = 8001
R_GAS_TAPER = 5.0
STAR_HZ_KPC = 0.3
STAR_MEASURED_RMAX_KPC = 13.75  # outer full S4G annulus; 14.25 kpc is only 92% covered
GAS_THIN_HZ_KPC = 0.2
GAS_THICK_HZ_KPC = 3.0
GAS_THIN_FRACTION = 0.85


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def load_profiles() -> dict[str, np.ndarray]:
    with PROFILE_PATH.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    radii = np.array([float(row["r_mid_kpc"]) for row in rows])

    def field(name: str) -> np.ndarray:
        return np.array([float(row[name]) if row[name] else np.nan for row in rows])

    # The P5 analysis ellipse clips the 14.0--14.5 kpc ring. Preserve the
    # 14.25-kpc value in the source CSV as a coverage sensitivity, but do not
    # interpolate through it as if it were a complete observed azimuthal ring.
    stars = field("s4g_stars_faceon_Msun_pc2_Upsilon0p6")
    stars[radii > STAR_MEASURED_RMAX_KPC] = np.nan

    return {
        "r": radii,
        "things": field("things_hihe_faceon_Msun_pc2"),
        "halogas": field("halogas_hihe_faceon_Msun_pc2"),
        "stars": stars,
    }


def load_rotmod() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with zipfile.ZipFile(ROTMOD_ZIP) as archive:
        raw = archive.read(ROTMOD_MEMBER)
    rows = np.loadtxt(io.BytesIO(raw), comments="#")
    return rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3], rows[:, 4]


def gas_surface_density(r_eval: np.ndarray, radii: np.ndarray, sigma: np.ndarray,
                        taper_kpc: float = R_GAS_TAPER) -> np.ndarray:
    valid = np.isfinite(sigma)
    x, y = radii[valid], sigma[valid]
    fn = PchipInterpolator(x, y, extrapolate=False)
    out = np.zeros_like(r_eval)
    inside = r_eval <= x[-1]
    out[inside] = fn(np.maximum(r_eval[inside], x[0]))
    outer = (r_eval > x[-1]) & (r_eval < x[-1] + taper_kpc)
    out[outer] = y[-1] * (1.0 - (r_eval[outer] - x[-1]) / taper_kpc)
    return np.maximum(out, 0.0)


def fit_stellar_outer_scale(radii: np.ndarray, sigma: np.ndarray) -> tuple[float, float, float]:
    valid = np.isfinite(sigma) & (radii >= 5.0) & (radii <= 12.75) & (sigma > 0)
    slope, intercept = np.polyfit(radii[valid], np.log(sigma[valid]), 1)
    if slope >= 0:
        raise ValueError(f"stellar outer profile did not decline: slope={slope}")
    residual = np.log(sigma[valid]) - (slope * radii[valid] + intercept)
    return -1.0 / slope, float(np.sqrt(np.mean(residual**2))), int(np.count_nonzero(valid))


def stellar_surface_density(r_eval: np.ndarray, radii: np.ndarray, sigma: np.ndarray,
                            *, tail_scale_factor: float = 1.0,
                            tail: str = "exponential") -> np.ndarray:
    valid = np.isfinite(sigma)
    x, y = radii[valid], sigma[valid]
    fn = PchipInterpolator(x, y, extrapolate=False)
    out = np.zeros_like(r_eval)
    inside = r_eval <= x[-1]
    out[inside] = fn(np.maximum(r_eval[inside], x[0]))
    if tail == "exponential":
        rd, _, _ = fit_stellar_outer_scale(radii, sigma)
        rd *= tail_scale_factor
        outer = r_eval > x[-1]
        out[outer] = y[-1] * np.exp(-(r_eval[outer] - x[-1]) / rd)
    elif tail == "truncate_2kpc":
        outer = (r_eval > x[-1]) & (r_eval < x[-1] + 2.0)
        out[outer] = y[-1] * (1.0 - (r_eval[outer] - x[-1]) / 2.0)
    elif tail != "hard_truncate":
        raise ValueError(tail)
    return np.maximum(out, 0.0)


def hankel_v2(r_eval: np.ndarray, r_source: np.ndarray, sigma_msun_pc2: np.ndarray,
              h_z_kpc: float, *, k_max: float = K_MAX, n_k: int = N_K) -> np.ndarray:
    """Return circular-speed squared for exp(-|z|/h)/(2h) vertical structure.

    V^2(R) = 2 pi G R int dk k J1(kR) SigmaTilde(k)/(1+k h),
    SigmaTilde(k) = int dR' R' Sigma(R') J0(k R').
    Sigma is converted from Msun/pc^2 to Msun/kpc^2 before integration.
    """
    sigma_kpc2 = sigma_msun_pc2 * 1.0e6
    k = np.linspace(0.0, k_max, n_k)
    transform = np.empty_like(k)
    # Chunk k to bound memory while doing the two independent Bessel integrals.
    weighted_source = r_source * sigma_kpc2
    chunk = 256
    for start in range(0, n_k, chunk):
        stop = min(start + chunk, n_k)
        kernel = j0(np.outer(k[start:stop], r_source))
        transform[start:stop] = np.trapezoid(kernel * weighted_source[None, :], r_source, axis=1)
    transform /= 1.0 + k * h_z_kpc
    out = np.zeros_like(r_eval)
    for idx, radius in enumerate(r_eval):
        if radius <= 0:
            continue
        integrand = k * j1(k * radius) * transform
        out[idx] = 2.0 * math.pi * G * radius * np.trapezoid(integrand, k)
    return out


def analytic_freeman_v2(radius: np.ndarray, sigma0_msun_pc2: float,
                        rd_kpc: float) -> np.ndarray:
    y = radius / (2.0 * rd_kpc)
    return 4.0 * math.pi * G * sigma0_msun_pc2 * 1.0e6 * rd_kpc * y**2 * (
        iv(0, y) * kv(0, y) - iv(1, y) * kv(1, y)
    )


def main() -> None:
    t0_wall = time.perf_counter()
    t0_cpu = time.process_time()
    profile = load_profiles()
    rdata = profile["r"]
    rsrc = np.arange(0.0, R_SRC_MAX + 0.5 * R_SRC_STEP, R_SRC_STEP)
    reval, vobs, err, vgas_sparc, vdisk_sparc = load_rotmod()
    sigma_gas_halo = gas_surface_density(rsrc, rdata, profile["halogas"])
    sigma_gas_things = gas_surface_density(rsrc, rdata, profile["things"])
    sigma_stars = stellar_surface_density(rsrc, rdata, profile["stars"], tail="exponential")
    sigma_stars_short = stellar_surface_density(rsrc, rdata, profile["stars"], tail="exponential",
                                                tail_scale_factor=0.75)
    sigma_stars_long = stellar_surface_density(rsrc, rdata, profile["stars"], tail="exponential",
                                                tail_scale_factor=1.25)
    sigma_stars_truncated = stellar_surface_density(rsrc, rdata, profile["stars"], tail="truncate_2kpc")
    rd_star, fit_rms_log, n_fit = fit_stellar_outer_scale(rdata, profile["stars"])

    v2_star = hankel_v2(reval, rsrc, sigma_stars, STAR_HZ_KPC)
    v2_star_short = hankel_v2(reval, rsrc, sigma_stars_short, STAR_HZ_KPC)
    v2_star_long = hankel_v2(reval, rsrc, sigma_stars_long, STAR_HZ_KPC)
    v2_star_truncated = hankel_v2(reval, rsrc, sigma_stars_truncated, STAR_HZ_KPC)
    v2_gas_thin = hankel_v2(reval, rsrc, sigma_gas_halo, GAS_THIN_HZ_KPC)
    v2_gas_thick = hankel_v2(reval, rsrc, sigma_gas_halo, GAS_THICK_HZ_KPC)
    v2_gas_halo = GAS_THIN_FRACTION * v2_gas_thin + (1.0 - GAS_THIN_FRACTION) * v2_gas_thick
    v2_gas_things = hankel_v2(reval, rsrc, sigma_gas_things, GAS_THIN_HZ_KPC)

    speed_star = np.sqrt(np.maximum(v2_star, 0.0))
    speed_gas_halo = np.sqrt(np.maximum(v2_gas_halo, 0.0))
    signed_speed_gas_halo = np.sign(v2_gas_halo) * np.sqrt(np.abs(v2_gas_halo))
    speed_total = np.sqrt(np.maximum(v2_star + v2_gas_halo, 0.0))
    speed_total_thin_gas = np.sqrt(np.maximum(v2_star + v2_gas_thin, 0.0))
    speed_total_things = np.sqrt(np.maximum(v2_star + v2_gas_things, 0.0))
    speed_total_star_short = np.sqrt(np.maximum(v2_star_short + v2_gas_halo, 0.0))
    speed_total_star_long = np.sqrt(np.maximum(v2_star_long + v2_gas_halo, 0.0))
    speed_total_star_truncated = np.sqrt(np.maximum(v2_star_truncated + v2_gas_halo, 0.0))
    speed_sparc_baryon = np.sqrt(np.maximum(vgas_sparc * np.abs(vgas_sparc) +
                                            0.6 * vdisk_sparc * np.abs(vdisk_sparc), 0.0))

    # Analytic thin exponential disk validation and spectral/grid convergence.
    analytic_r = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 15.0])
    analytic_sigma0, analytic_rd = 800.0, 3.0
    analytic_source = analytic_sigma0 * np.exp(-rsrc / analytic_rd)
    v2_analytic_num = hankel_v2(analytic_r, rsrc, analytic_source, 0.0,
                                k_max=K_MAX, n_k=N_K)
    v2_analytic_exact = analytic_freeman_v2(analytic_r, analytic_sigma0, analytic_rd)
    analytic_rel = np.abs(v2_analytic_num - v2_analytic_exact) / v2_analytic_exact

    conv_v2_star = hankel_v2(reval, rsrc, sigma_stars, STAR_HZ_KPC, k_max=60.0, n_k=6001)
    conv_v2_gas = 0.85 * hankel_v2(reval, rsrc, sigma_gas_halo, GAS_THIN_HZ_KPC,
                                   k_max=60.0, n_k=6001) + 0.15 * hankel_v2(
                                       reval, rsrc, sigma_gas_halo, GAS_THICK_HZ_KPC,
                                       k_max=60.0, n_k=6001)
    conv_speed = np.sqrt(np.maximum(conv_v2_star + conv_v2_gas, 0.0))
    convergence_fraction = np.abs(conv_speed - speed_total) / np.maximum(speed_total, 1.0)

    # The pointwise SPARC errors are not a full covariance; this score is a
    # descriptive diagonal comparison only.
    chi2_k0 = float(np.sum(((vobs - speed_total) / err) ** 2))
    chi2_k0_things = float(np.sum(((vobs - speed_total_things) / err) ** 2))
    chi2_rotmod_baryon = float(np.sum(((vobs - speed_sparc_baryon) / err) ** 2))
    valid_stars = reval <= STAR_MEASURED_RMAX_KPC  # only complete P5-supported rings
    component_compare = {
        "stars_R_le_13p75_kpc_full_rings": {
            "n": int(valid_stars.sum()),
            "rms_km_s_map_minus_sparc_Upsilon0p6": float(np.sqrt(np.mean((speed_star[valid_stars] -
                                                                            math.sqrt(0.6) * vdisk_sparc[valid_stars]) ** 2))),
        },
        "gas_R_le_44p08_kpc": {
            "n": int(reval.size),
            "rms_km_s_signed_map_minus_sparc_gas": float(np.sqrt(np.mean((signed_speed_gas_halo - vgas_sparc) ** 2))),
        },
        "baryon_total_R_le_44p08_kpc": {
            "rms_km_s_map_minus_sparc_Upsilon0p6": float(np.sqrt(np.mean((speed_total - speed_sparc_baryon) ** 2))),
            "max_abs_km_s_map_minus_sparc_Upsilon0p6": float(np.max(np.abs(speed_total - speed_sparc_baryon))),
        },
    }

    rows = []
    for i, radius in enumerate(reval):
        rows.append({
            "R_kpc": radius,
            "Vobs_km_s": vobs[i],
            "sigmaV_km_s": err[i],
            "Vstar_map_K0_km_s": speed_star[i],
            "Vgas_map_HALOGAS_thin_thick_K0_km_s": speed_gas_halo[i],
            "Vbar_map_HALOGAS_K0_km_s": speed_total[i],
            "Vbar_map_THINGS_thin_K0_km_s": speed_total_things[i],
            "Vbar_map_HALOGAS_thin_only_K0_km_s": speed_total_thin_gas[i],
            "Vbar_SPARC_Upsilon0p6_km_s": speed_sparc_baryon[i],
            "Vbar_residual_km_s": vobs[i] - speed_total[i],
        })
    with (OUT / "newtonian_k0_predictions.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    data = np.column_stack([reval, vobs, err, speed_star, speed_gas_halo, speed_total,
                            speed_total_things, speed_total_thin_gas, speed_sparc_baryon])
    figure, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
    ax.errorbar(reval, vobs, yerr=err, fmt="o", ms=3, color="black", alpha=0.75, label="SPARC observed")
    ax.plot(reval, speed_total, color="tab:blue", label="Resolved maps, Newtonian K=0")
    ax.plot(reval, speed_total_things, "--", color="tab:green", label="THINGS gas, thin")
    ax.plot(reval, speed_total_thin_gas, ":", color="tab:cyan", label="HALOGAS thin-only")
    ax.plot(reval, speed_sparc_baryon, "-.", color="tab:orange", label="SPARC baryons, Υ=0.6")
    ax.axvline(14.4, color="gray", alpha=0.45, lw=1, label="S4G measured aperture")
    ax.set(xlabel="R [kpc]", ylabel="circular speed [km s$^{-1}$]", xlim=(0, 46), ylim=(0, 180))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    figure.savefig(OUT / "newtonian_k0_baseline.png", dpi=180)
    plt.close(figure)

    record = {
        "status": "complete_source_derived_Newtonian_K0_baseline; no TG fit and no halo_or_MOND_fit",
        "interpretation": "Exploratory conditional forward-force comparison. The SPARC error column is treated diagonally because the release supplies no per-galaxy covariance here.",
        "inputs": {
            "source_profiles_csv_sha256": sha256_file(PROFILE_PATH),
            "rotmod_archive_sha256": sha256_file(ROTMOD_ZIP),
            "rotmod_member": ROTMOD_MEMBER,
            "rotmod_member_sha256": sha256_bytes(zipfile.ZipFile(ROTMOD_ZIP).read(ROTMOD_MEMBER)),
            "n_source_profile_bins": int(np.count_nonzero(np.isfinite(profile["halogas"]))),
            "n_sparc_points": int(reval.size),
        },
        "model": {
            "surface_density_units": "Msun/pc^2 face-on",
            "axisymmetric_assumption": "azimuthally averaged resolved 2D maps",
            "stars": "S4G measured to its 14.4-kpc ICA ellipse, M/L=0.6; use only complete rings through 13.75 kpc (the 14.25-kpc ring is partial); exponential tail beyond last full ring, fitted Rd from 5-12.75 kpc",
            "stellar_full_ring_support_kpc": STAR_MEASURED_RMAX_KPC,
            "star_scale_length_kpc": rd_star,
            "star_outer_fit_rms_ln_sigma": fit_rms_log,
            "star_outer_fit_bins": n_fit,
            "gas": "HALOGAS primary, 1.36 helium factor already included; 5-kpc linear taper beyond R=34.75 kpc",
            "gas_vertical": "85% exp(|z|/h=0.2 kpc) + 15% exp(|z|/h=3 kpc); thin-only is separately tested",
            "stellar_vertical_scale_kpc": STAR_HZ_KPC,
            "gravity": "Newtonian Poisson equation, K=0; no halo and no non-baryonic component",
            "velocity_kernel": "V^2=2*pi*G*R*integral dk*k*J1(kR)*SigmaTilde(k)/(1+k*h), SigmaTilde=integral dRprime*Rprime*Sigma(Rprime)*J0(kRprime)",
        },
        "numerics": {
            "G_kpc_km2_s2_Msun": G,
            "R_source_step_kpc": R_SRC_STEP,
            "R_source_max_kpc": R_SRC_MAX,
            "k_max_kpc_inverse": K_MAX,
            "n_k": N_K,
            "freeman_exponential_disk_test_max_relative_v2_error": float(np.max(analytic_rel)),
            "freeman_test_radii_kpc": analytic_r.tolist(),
            "freeman_relative_errors": analytic_rel.tolist(),
            "k_spectral_convergence_max_fractional_speed_change": float(np.max(convergence_fraction)),
            "no_random_seed": True,
        },
        "results": {
            "chi2_resolved_baryons_vs_SPARC_diagonal": chi2_k0,
            "chi2_THINGS_gas_thin_vs_SPARC_diagonal": chi2_k0_things,
            "chi2_SPARC_baryons_Upsilon0p6_vs_SPARC_diagonal": chi2_rotmod_baryon,
        "component_comparisons": component_compare,
            "mean_abs_residual_km_s": float(np.mean(np.abs(vobs - speed_total))),
            "max_abs_residual_km_s": float(np.max(np.abs(vobs - speed_total))),
            "source_tail_sensitivity_max_speed_delta_km_s": float(max(
                np.max(np.abs(speed_total - speed_total_star_short)),
                np.max(np.abs(speed_total - speed_total_star_long)),
                np.max(np.abs(speed_total - speed_total_star_truncated)),
            )),
            "source_tail_sensitivity": {
                "star_Rd_x0p75_max_delta_km_s": float(np.max(np.abs(speed_total - speed_total_star_short))),
                "star_Rd_x1p25_max_delta_km_s": float(np.max(np.abs(speed_total - speed_total_star_long))),
                "star_2kpc_taper_max_delta_km_s": float(np.max(np.abs(speed_total - speed_total_star_truncated))),
                "gas_thin_only_max_delta_km_s": float(np.max(np.abs(speed_total - speed_total_thin_gas))),
                "THINGS_gas_thin_max_delta_km_s": float(np.max(np.abs(speed_total - speed_total_things))),
            },
        },
        "runtime": {
            "wall_seconds": time.perf_counter() - t0_wall,
            "process_cpu_seconds": time.process_time() - t0_cpu,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "files": {
            "code_sha256": sha256_file(Path(__file__)),
            "result_csv": "work/tg_axisymmetric_pilot/newtonian_k0_predictions.csv",
            "plot": "work/tg_axisymmetric_pilot/newtonian_k0_baseline.png",
        },
    }
    (OUT / "newtonian_k0_baseline.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["results"], indent=2))


if __name__ == "__main__":
    main()

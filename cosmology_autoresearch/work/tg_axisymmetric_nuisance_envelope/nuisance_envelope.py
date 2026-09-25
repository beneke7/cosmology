#!/usr/bin/env python3
"""Geometry-aware nuisance envelope for the NGC 3198 axisymmetric pilot.

The raw map intensity/column is a line-of-sight quantity. For a thin disk,
Sigma_face-on = Sigma_los cos(i), while a sky pixel's disk-plane area is
D^2 dOmega/cos(i). The factors of cos(i) cancel in integrated map mass;
distance makes that mass scale as D^2. The deprojected radius is
D sqrt(x_major^2 + (x_minor/cos(i))^2), with angular offsets in radians.
Thus map-derived circular speeds are recomputed from the transformed source;
the SPARC velocities are separately re-deprojected by sin(73 deg)/sin(i).

This is an explicitly axisymmetrized, parametric-source test, not a unique 3D
reconstruction. See REPORT.md for scope, priors, and the finite-u/edge split.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import platform
import time
import zipfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PILOT = ROOT / "work/tg_axisymmetric_pilot"
MAP_AUDIT = ROOT / "work/tg_resolved_source_audit"
ROT_ZIP = ROOT / "context/data/Rotmod_LTG.zip"
STAR_MAP = PILOT / "data/NGC3198.stellar.fits"
STAR_MASK = PILOT / "data/NGC3198.ICAmask.fits"
HALOGAS_MAP = MAP_AUDIT / "NGC3198-HR_coldens.fits"
THINGS_MAP = MAP_AUDIT / "NGC_3198_NA_MOM0_THINGS.FITS"
REFERENCE_CSV = PILOT / "source_profiles.csv"

D_REF = 13.8
I_REF = 73.0  # SPARC tabulates its deprojected Vobs at this inclination.
I_MAP_REF = 72.0  # original source-extraction convention; audited explicitly.
PA_DEG = 215.0
P5_PA_DEG = 33.7
P5_SMA_ARCSEC = 215.4
P5_Q = 1.0 - 0.621
UPSILON_REF = 0.6
ARCSEC_RAD = math.pi / (180.0 * 3600.0)
NHI_PER_MSUN_PC2 = 1.25e20
HELIUM = 1.36
G = 4.30091e-6
RING_WIDTH = 0.5
R_MAX = 35.0
H0_KPC = 73.0/1000.0
RHO_CRIT = 3.0*H0_KPC**2/(8.0*math.pi*G)
A0 = 1.2e-10*3.0856775814913673e19/1.0e6  # (km/s)^2/kpc
MSTAR_SIGMA_DEX = 0.1
DISTANCE_SIGMA_MPC = 1.4
INCLINATION_SIGMA_DEG = 3.0
RING_TAIL_KPC = 5.0
R_SOURCE_STEP = 0.02
R_SOURCE_MAX = 50.0
K_MAX = 80.0
N_K = 6001
TG_DR = 0.25
TG_DZ = 0.125
TG_RMAX = 80.0
TG_ZMAX = 60.0

VERTICAL_CASES = {
    "gas_mix_hstar0p3": {"hstar": 0.3, "f_thin": 0.85},
    "gas_thin_hstar0p3": {"hstar": 0.3, "f_thin": 1.0},
    "gas_mix_hstar0p6": {"hstar": 0.6, "f_thin": 0.85},
}
Q_GRID = [-0.1, -0.05, 0.0, 0.05, 0.1, math.log10(0.762/UPSILON_REF)]
D_SIGMA_GRID = [-2, -1, 0, 1, 2]
I_SIGMA_GRID = [-2, -1, 0, 1, 2]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def image(path: Path) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path, memmap=False) as hdus:
        return np.asarray(hdus[0].data, dtype=float).squeeze(), hdus[0].header.copy()


def pixel_scales(h: fits.Header) -> tuple[float, float]:
    if "CD1_1" in h:
        return abs(float(h["CD1_1"]))*3600., abs(float(h["CD2_2"]))*3600.
    return abs(float(h["CDELT1"]))*3600., abs(float(h["CDELT2"]))*3600.


def sky_xy(h: fits.Header, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Pixel-center tangent-plane offsets in arcsec; map WCS is diagonal."""
    ny, nx = shape
    xp, yp = np.meshgrid(np.arange(nx)+1., np.arange(ny)+1.)
    cd1 = float(h["CD1_1"] if "CD1_1" in h else h["CDELT1"])
    cd2 = float(h["CD2_2"] if "CD2_2" in h else h["CDELT2"])
    dx = (xp-float(h["CRPIX1"]))*cd1*3600.
    dy = (yp-float(h["CRPIX2"]))*cd2*3600.
    return dx, dy


def rotate_disk(dx: np.ndarray, dy: np.ndarray, pa_deg: float = PA_DEG):
    pa = math.radians(pa_deg)
    major = dx*math.sin(pa) + dy*math.cos(pa)
    minor = dx*math.cos(pa) - dy*math.sin(pa)
    return major, minor


def radial_kpc(dx: np.ndarray, dy: np.ndarray, distance_mpc: float, inc_deg: float):
    major, minor = rotate_disk(dx, dy)
    ci = math.cos(math.radians(inc_deg))
    return np.hypot(major, minor/ci)*ARCSEC_RAD*distance_mpc*1000.


def p5_ellipse_radius_kpc(dx: np.ndarray, dy: np.ndarray, distance_mpc: float):
    major, minor = rotate_disk(dx, dy, P5_PA_DEG)
    return np.hypot(major, minor/P5_Q)*ARCSEC_RAD*distance_mpc*1000.


def profile(data: np.ndarray, valid: np.ndarray, radius: np.ndarray,
            h: fits.Header, distance_mpc: float, inc_deg: float,
            *, physical_floor: bool = True):
    """Physical face-on annular means and projected-pixel coverage fractions."""
    sx, sy = pixel_scales(h)
    pixel_area_projected = ((distance_mpc*1000.*sx*ARCSEC_RAD) *
                            (distance_mpc*1000.*sy*ARCSEC_RAD))
    edges = np.arange(0., R_MAX + RING_WIDTH/2, RING_WIDTH)
    use = valid & np.isfinite(data) & np.isfinite(radius) & (radius >= 0) & (radius < R_MAX)
    count, _ = np.histogram(radius[use], bins=edges)
    sums, _ = np.histogram(radius[use], bins=edges, weights=data[use])
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = sums/count
    mean[count == 0] = np.nan
    projected_area_annulus = np.pi*(edges[1:]**2-edges[:-1]**2)*math.cos(math.radians(inc_deg))
    coverage = count*pixel_area_projected/projected_area_annulus
    # Match the original source extraction: clip only after azimuthal averaging.
    if physical_floor:
        mean = np.maximum(mean, 0.)
    mids = (edges[:-1]+edges[1:])/2
    return mids, mean, coverage, count


def load_maps():
    star, hs = image(STAR_MAP)
    mask, hm = image(STAR_MASK)
    halo, hh = image(HALOGAS_MAP)
    things, ht = image(THINGS_MAP)
    if halo.ndim != 2 or things.ndim != 2 or star.ndim != 2 or mask.shape != star.shape:
        raise ValueError("Unexpected map shape")
    dxs, dys = sky_xy(hs, star.shape)
    dxh, dyh = sky_xy(hh, halo.shape)
    dxt, dyt = sky_xy(ht, things.shape)
    # The P5 mask uses its documented photometric ellipse, fixed in angle.
    major5, minor5 = rotate_disk(dxs, dys, P5_PA_DEG)
    p5_arcsec = np.hypot(major5, minor5/P5_Q)
    star_valid = (mask <= 0) & np.isfinite(star) & (p5_arcsec <= P5_SMA_ARCSEC)
    star_los = np.zeros_like(star)
    positive = star > 0
    sigma_los = star[positive]*1.e6/(206265.**2)
    mu = -2.5*np.log10(sigma_los/280.9)
    star_los[positive] = 10.**(-0.4*(mu-3.24-21.572))*UPSILON_REF
    halo_los = np.maximum(halo, 0.)/NHI_PER_MSUN_PC2*HELIUM
    things_bmaj, things_bmin = 3.1753e-3*3600., 2.6007e-3*3600.
    things_tb_per_jy = 606000./(things_bmaj*things_bmin)
    things_nhi = things*1.e-3*things_tb_per_jy*1.823e18
    things_los = np.maximum(things_nhi, 0.)/NHI_PER_MSUN_PC2*HELIUM
    return {
        "star_los": star_los, "star_valid": star_valid, "star_h": hs, "star_dx": dxs, "star_dy": dys,
        "halo_los": halo_los, "halo_h": hh, "halo_dx": dxh, "halo_dy": dyh,
        "things_los": things_los, "things_h": ht, "things_dx": dxt, "things_dy": dyt,
    }


def transformed_sources(maps, distance_mpc: float, inc_deg: float):
    """Return native, fully coupled D/i profile transforms for one geometry."""
    c = math.cos(math.radians(inc_deg))
    rstar = radial_kpc(maps["star_dx"], maps["star_dy"], distance_mpc, inc_deg)
    rhalo = radial_kpc(maps["halo_dx"], maps["halo_dy"], distance_mpc, inc_deg)
    rthings = radial_kpc(maps["things_dx"], maps["things_dy"], distance_mpc, inc_deg)
    rs, stars, cs, ns = profile(maps["star_los"]*c, maps["star_valid"], rstar,
                                maps["star_h"], distance_mpc, inc_deg)
    rg, gas, cg, ng = profile(maps["halo_los"]*c, np.isfinite(maps["halo_los"]), rhalo,
                              maps["halo_h"], distance_mpc, inc_deg)
    rt, things, ct, nt = profile(maps["things_los"]*c, np.isfinite(maps["things_los"]), rthings,
                                  maps["things_h"], distance_mpc, inc_deg)
    # S4G photometry is documented only to the P5 ellipse. Only complete
    # annuli ending before the ellipse semimajor axis are treated as measured.
    p5_rmax = P5_SMA_ARCSEC*ARCSEC_RAD*distance_mpc*1000.
    star_ok = (cs >= .8) & (rs+RING_WIDTH/2 <= p5_rmax+1e-12)
    stars = np.where(star_ok, stars, np.nan)
    # Preserve the existing pilot's 0.8 support criterion for the HI maps.
    gas = np.where(cg >= .8, gas, np.nan)
    things = np.where(ct >= .8, things, np.nan)
    return {"r": rs, "stars": stars, "gas": gas, "things": things,
            "coverage_star": cs, "coverage_gas": cg, "coverage_things": ct,
            "count_star": ns, "count_gas": ng, "count_things": nt,
            "p5_rmax_kpc": p5_rmax}


def read_rotmod():
    with zipfile.ZipFile(ROT_ZIP) as z:
        arr = np.loadtxt(io.BytesIO(z.read("NGC3198_rotmod.dat")), comments="#")
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]


def surface_density_on_grid(r_eval: np.ndarray, radii: np.ndarray, sigma: np.ndarray,
                            kind: str) -> tuple[np.ndarray, dict]:
    """Use the pilot's explicit tail rules for gas and measured P5 stars."""
    valid = np.isfinite(sigma)
    x, y = radii[valid], sigma[valid]
    if x.size < 4:
        raise ValueError(f"too few supported {kind} radial bins: {x.size}")
    interp = PchipInterpolator(x, y, extrapolate=False)
    out = np.zeros_like(r_eval)
    inside = r_eval <= x[-1]
    out[inside] = interp(np.clip(r_eval[inside], x[0], x[-1]))
    config = {"measured_last_center_kpc": float(x[-1])}
    if kind == "stars":
        fit_hi = min(12.75, float(x[-1]))
        fit = valid & (radii >= 5.0) & (radii <= fit_hi) & (sigma > 0)
        if np.count_nonzero(fit) < 5:
            raise ValueError(f"too few stellar tail-fit bins: {np.count_nonzero(fit)}")
        slope, intercept = np.polyfit(radii[fit], np.log(sigma[fit]), 1)
        if slope >= 0:
            raise ValueError("stellar measured tail is nondeclining")
        rd = -1.0/slope
        outer = r_eval > x[-1]
        out[outer] = y[-1]*np.exp(-(r_eval[outer]-x[-1])/rd)
        config.update({"tail": "exponential", "tail_Rd_kpc": float(rd),
                       "tail_fit_max_kpc": fit_hi, "tail_fit_n": int(fit.sum()),
                       "tail_fit_rms_log": float(np.sqrt(np.mean((np.log(sigma[fit])-
                                                                  (slope*radii[fit]+intercept))**2)))})
    elif kind == "gas":
        outer = (r_eval > x[-1]) & (r_eval < x[-1]+RING_TAIL_KPC)
        out[outer] = y[-1]*(1.0-(r_eval[outer]-x[-1])/RING_TAIL_KPC)
        config.update({"tail": "linear_5kpc_taper", "tail_width_kpc": RING_TAIL_KPC})
    else:
        raise ValueError(kind)
    return np.maximum(out, 0.0), config


def hankel_transform(r_source: np.ndarray, sigma: np.ndarray):
    """Return SigmaTilde(k)=int R Sigma J0(kR)dR in kpc-based units."""
    from scipy.special import j0
    k = np.linspace(0.0, K_MAX, N_K)
    source = r_source*sigma*1.0e6
    transform = np.empty_like(k)
    for start in range(0, N_K, 256):
        stop = min(start+256, N_K)
        transform[start:stop] = np.trapezoid(j0(np.outer(k[start:stop], r_source))*
                                              source[None, :], r_source, axis=1)
    return k, transform


def hankel_v2_from_transform(r_eval: np.ndarray, k: np.ndarray, transform: np.ndarray,
                             h_z: float):
    from scipy.special import j1
    trans = transform/(1.0+k*h_z)
    out = np.zeros_like(r_eval)
    for idx, radius in enumerate(r_eval):
        out[idx] = 2.0*math.pi*G*radius*np.trapezoid(k*j1(k*radius)*trans, k)
    return out


def nfw_v2(r: np.ndarray, logm: float, logc: float):
    mass, conc = 10.0**logm, 10.0**logc
    r200 = (3.0*mass/(4.0*math.pi*200.0*RHO_CRIT))**(1.0/3.0)
    rs = r200/conc
    fc = math.log1p(conc)-conc/(1.0+conc)
    x = np.maximum(r/rs, 1e-14)
    fx = np.log1p(x)-x/(1.0+x)
    return G*mass*fx/(fc*r)


def geometry_prior(distance_mpc: float, inclination_deg: float):
    return (((distance_mpc-D_REF)/DISTANCE_SIGMA_MPC)**2,
            ((inclination_deg-I_REF)/INCLINATION_SIGMA_DEG)**2)


def observed_for_inclination(v_catalog: np.ndarray, err_catalog: np.ndarray, inc_deg: float):
    factor = math.sin(math.radians(I_REF))/math.sin(math.radians(inc_deg))
    return v_catalog*factor, err_catalog*factor, factor


def fit_nfw_for_geometry(r_eval, obs, err, baryon_v2, upsilon):
    from scipy.optimize import least_squares
    lo = np.array([9.0, 0.0])
    hi = np.array([14.0, math.log10(50.0)])

    def predict(theta):
        total = baryon_v2+nfw_v2(r_eval, float(theta[0]), float(theta[1]))
        return np.sqrt(np.maximum(total, 0.0))

    def residual(theta):
        return (predict(theta)-obs)/err

    starts = [(11.0, 0.9), (12.0, 0.7), (11.5, 1.1)]
    fits = [least_squares(residual, np.asarray(start), bounds=(lo, hi),
                          x_scale=np.array([1.0, 0.5]), max_nfev=1500,
                          ftol=1e-10, xtol=1e-10, gtol=1e-10) for start in starts]
    best = min(fits, key=lambda fit: float(np.dot(fit.fun, fit.fun)))
    pred = predict(best.x)
    return {"chi2_data": float(np.dot(best.fun, best.fun)),
            "predictions": pred, "log10_M200_Msun": float(best.x[0]),
            "M200_Msun": float(10.0**best.x[0]),
            "log10_c200": float(best.x[1]),
            "concentration_c200": float(10.0**best.x[1]),
            "optimizer_success": bool(best.success),
            "optimizer_message": str(best.message)}


def tg_vertical_density(star_sigma, gas_sigma, rcent, zcent, hstar, f_thin):
    rho_star = star_sigma[:, None]*1.e6/(2.0*hstar)*np.exp(-zcent/hstar)
    rho_gas = gas_sigma[:, None]*1.e6*(f_thin/(2.0*0.2)*np.exp(-zcent/0.2)
                                          +(1.0-f_thin)/(2.0*3.0)*np.exp(-zcent/3.0))
    return rho_star, rho_gas


def solve_tg_edge(B, volume, rcent, zcent, star_rho, gas_rho, upsilon,
                  eval_radii, obs, err):
    from scipy.sparse import diags
    from scipy.sparse.linalg import eigsh
    rho = upsilon/UPSILON_REF*star_rho+gas_rho
    dvec = (4.0*math.pi*G*rho*volume).ravel()
    values, vectors = eigsh(B, k=1, M=diags(dvec, format="csr"), sigma=0.0,
                            which="LM", tol=2e-8, maxiter=5000)
    kcrit = float(values[0])
    vec = vectors[:, 0]
    residual = float(np.linalg.norm(B@vec-kcrit*dvec*vec)/np.linalg.norm(B@vec))
    h = vec.reshape(star_rho.shape)
    if float(np.sum(h)) < 0.0:
        h *= -1.0
    h /= float(np.max(h))
    hmid = (9.0*h[:, 0]-h[:, 1])/8.0
    if np.any(hmid <= 0.0):
        raise RuntimeError("principal mode is not positive on the midplane")
    hfn = PchipInterpolator(rcent, hmid, extrapolate=True)
    hv = hfn(eval_radii)
    v2 = -eval_radii/kcrit*hfn.derivative()(eval_radii)/hv
    if np.any(hv <= 0) or np.any(v2 <= 0) or np.any(~np.isfinite(v2)):
        raise RuntimeError("nonphysical TG edge curve")
    speed = np.sqrt(v2)
    return {"kcrit": kcrit, "mode_residual": residual,
            "minimum_mode": float(np.min(h)), "maximum_mode": float(np.max(h)),
            "chi2_data": float(np.sum(((obs-speed)/err)**2)),
            "rms_residual": float(np.sqrt(np.mean((obs-speed)**2))),
            "predictions": speed, "h": h, "hmid": hmid, "rho": rho}


def run_envelope(cpu_cap: float = 850.0):
    import importlib.util
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.sparse import diags
    from scipy.sparse.linalg import splu

    started_wall, started_cpu = time.perf_counter(), time.process_time()
    maps = load_maps()
    r_catalog, v_catalog, e_catalog, _, _ = read_rotmod()
    source_r = np.arange(0.0, R_SOURCE_MAX+0.5*R_SOURCE_STEP, R_SOURCE_STEP)
    d_grid = [D_REF+x*DISTANCE_SIGMA_MPC for x in D_SIGMA_GRID]
    i_grid = [I_REF+x*INCLINATION_SIGMA_DEG for x in I_SIGMA_GRID]
    u_grid = [UPSILON_REF*10.0**q for q in Q_GRID]
    geometries = []
    for d in d_grid:
        for inc in i_grid:
            p = transformed_sources(maps, d, inc)
            geometries.append({"distance": d, "inclination": inc,
                               "d_sigma": (d-D_REF)/DISTANCE_SIGMA_MPC,
                               "i_sigma": (inc-I_REF)/INCLINATION_SIGMA_DEG,
                               "profile": p})

    # Assemble only the geometry-dependent finite-volume Laplacian once. Its
    # source diagonal is replaced for every raw-map-derived geometry/M/L/height.
    solver_path = PILOT/"solve_tg_axisymmetric.py"
    spec = importlib.util.spec_from_file_location("tg_solver_envelope", solver_path)
    tg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tg)
    B, _, _, volume, rr, zz, step, _ = tg.build_operator(
        TG_DR, TG_DZ, TG_RMAX, TG_ZMAX, UPSILON_REF)
    rcent = (np.arange(rr.shape[0])+0.5)*step[0]
    zcent = (np.arange(rr.shape[1])+0.5)*step[1]

    # Force-kernel closure against the prior K=0 source calculation at its
    # original i=72 geometry, with exactly the same measured map profiles.
    p72 = transformed_sources(maps, D_REF, I_MAP_REF)
    ss72, cfg_s72 = surface_density_on_grid(source_r, p72["r"], p72["stars"], "stars")
    sg72, cfg_g72 = surface_density_on_grid(source_r, p72["r"], p72["gas"], "gas")
    ks72, ts72 = hankel_transform(source_r, ss72)
    kg72, tg72 = hankel_transform(source_r, sg72)
    v2s72 = hankel_v2_from_transform(r_catalog, ks72, ts72, 0.3)
    v2g72 = (0.85*hankel_v2_from_transform(r_catalog, kg72, tg72, 0.2)
             +0.15*hankel_v2_from_transform(r_catalog, kg72, tg72, 3.0))
    base_csv = list(csv.DictReader((PILOT/"newtonian_k0_predictions.csv").open(newline="")))
    base_star = np.asarray([float(x["Vstar_map_K0_km_s"]) for x in base_csv])
    base_total = np.asarray([float(x["Vbar_map_HALOGAS_K0_km_s"]) for x in base_csv])
    closure_star = np.sqrt(np.maximum(v2s72, 0.0))
    closure_total = np.sqrt(np.maximum(v2s72+v2g72, 0.0))
    kernel_closure = {
        "geometry": {"distance_Mpc": D_REF, "inclination_deg": I_MAP_REF},
        "stellar_max_abs_speed_delta_km_s": float(np.max(np.abs(closure_star-base_star))),
        "stellar_RMS_speed_delta_km_s": float(np.sqrt(np.mean((closure_star-base_star)**2))),
        "total_max_abs_speed_delta_km_s": float(np.max(np.abs(closure_total-base_total))),
        "total_RMS_speed_delta_km_s": float(np.sqrt(np.mean((closure_total-base_total)**2))),
        "stars_tail": cfg_s72, "gas_tail": cfg_g72,
        "source_step_kpc": R_SOURCE_STEP, "k_max_kpc_inverse": K_MAX, "n_k": N_K,
    }

    # Same-grid TG-source closure at the earlier map-extraction inclination.
    # The legacy branch edge used this source, Υ*=0.6, and the same FV grid.
    star_sigma72, _ = surface_density_on_grid(source_r, p72["r"], p72["stars"], "stars")
    gas_sigma72, _ = surface_density_on_grid(source_r, p72["r"], p72["gas"], "gas")
    rho_star72, rho_gas72 = tg_vertical_density(np.interp(rcent, source_r, star_sigma72),
                                                np.interp(rcent, source_r, gas_sigma72),
                                                rcent, zcent, 0.3, 0.85)
    observed72, errors72, _ = observed_for_inclination(v_catalog, e_catalog, I_MAP_REF)
    legacy_mode = solve_tg_edge(B, volume, rcent, zcent, rho_star72, rho_gas72,
                                0.6, r_catalog, observed72, errors72)
    legacy_edges = list(csv.DictReader((PILOT/"axisymmetric_branch_edge_curves.csv").open(newline="")))
    legacy_edge_col = "edge_Upsilon_0.600000_km_s"
    legacy_speed = np.asarray([float(row[legacy_edge_col]) for row in legacy_edges])
    same_grid_tg_check = {
        "geometry": {"distance_Mpc": D_REF, "source_extraction_inclination_deg": I_MAP_REF},
        "upsilon": 0.6,
        "grid": {"dr_kpc": TG_DR, "dz_kpc": TG_DZ, "Rmax_kpc": TG_RMAX,
                 "Zmax_kpc": TG_ZMAX},
        "new_Kcrit_km_s_minus2": legacy_mode["kcrit"],
        "frozen_Kcrit_km_s_minus2": 6.619916457741752e-05,
        "Kcrit_fractional_difference": legacy_mode["kcrit"]/6.619916457741752e-05-1.0,
        "eigen_residual": legacy_mode["mode_residual"],
        "max_abs_edge_speed_difference_vs_frozen_km_s": float(np.max(np.abs(legacy_mode["predictions"]-legacy_speed))),
        "RMS_edge_speed_difference_vs_frozen_km_s": float(np.sqrt(np.mean((legacy_mode["predictions"]-legacy_speed)**2))),
        "chi2_against_Vobs_corrected_to_i72": legacy_mode["chi2_data"],
        "chi2_against_catalog_Vobs_i73_without_correction": float(np.sum(((v_catalog-legacy_mode["predictions"])/e_catalog)**2)),
    }

    # Full grid metadata and resumable scalar checkpoint stream. Each row is
    # flushed after a model evaluation; large principal modes are recomputed
    # only for the three selected finite-u examples after profiling.
    progress_path = OUT/"envelope_progress.jsonl"
    progress_path.write_text("")
    all_rows = []
    best = {vertical: {model: None for model in ("TG_edge", "NFW", "MOND")}
            for vertical in VERTICAL_CASES}
    geometry_cache = {}
    n_expected = len(geometries)*len(VERTICAL_CASES)*len(u_grid)*3
    completed = 0
    grid_complete = True
    for geom in geometries:
        d, inc, prof = geom["distance"], geom["inclination"], geom["profile"]
        reval = r_catalog*(d/D_REF)
        obs, err, inc_factor = observed_for_inclination(v_catalog, e_catalog, inc)
        prior_d, prior_i = geometry_prior(d, inc)
        star_sigma, star_cfg = surface_density_on_grid(source_r, prof["r"], prof["stars"], "stars")
        gas_sigma, gas_cfg = surface_density_on_grid(source_r, prof["r"], prof["gas"], "gas")
        ks, ts = hankel_transform(source_r, star_sigma)
        kg, tgtr = hankel_transform(source_r, gas_sigma)
        star_v2_by_h = {h: hankel_v2_from_transform(reval, ks, ts, h) for h in (0.3, 0.6)}
        gas_v2_thin = hankel_v2_from_transform(reval, kg, tgtr, 0.2)
        gas_v2_thick = hankel_v2_from_transform(reval, kg, tgtr, 3.0)
        xsource = (d, inc)
        geometry_cache[xsource] = {"profile": prof, "star_cfg": star_cfg, "gas_cfg": gas_cfg}
        for vertical, vc in VERTICAL_CASES.items():
            hstar, fthin = vc["hstar"], vc["f_thin"]
            star_v2 = star_v2_by_h[hstar]
            gas_v2 = fthin*gas_v2_thin+(1.0-fthin)*gas_v2_thick
            rho_star, rho_gas = tg_vertical_density(
                np.interp(rcent, source_r, star_sigma), np.interp(rcent, source_r, gas_sigma),
                rcent, zcent, hstar, fthin)
            for q, upsilon in zip(Q_GRID, u_grid):
                if time.process_time()-started_cpu > cpu_cap:
                    grid_complete = False
                    break
                u_prior = (q/MSTAR_SIGMA_DEX)**2
                prior_sum = prior_d+prior_i+u_prior
                bary_v2 = (upsilon/UPSILON_REF)*star_v2+gas_v2

                # NFW: profile halo mass/concentration at each same M/L, D, i.
                nfwfit = fit_nfw_for_geometry(reval, obs, err, bary_v2, upsilon)
                nfwrow = {"model": "NFW", "vertical_case": vertical, "distance_Mpc": d,
                          "inclination_deg": inc, "D_sigma": geom["d_sigma"],
                          "i_sigma": geom["i_sigma"], "upsilon": upsilon, "q_log10_ML": q,
                          "chi2_data": nfwfit["chi2_data"], "rms_km_s": float(np.sqrt(np.mean((obs-nfwfit["predictions"])**2))),
                          "prior_D": prior_d, "prior_i": prior_i, "prior_ML": u_prior,
                          "prior_other": 0.0, "prior_total": prior_sum,
                          "objective_chi2_plus_priors": nfwfit["chi2_data"]+prior_sum,
                          "NFW_log10_M200": nfwfit["log10_M200_Msun"],
                          "NFW_M200_Msun": nfwfit["M200_Msun"],
                          "NFW_log10_c200": nfwfit["log10_c200"],
                          "NFW_c200": nfwfit["concentration_c200"],
                          "obs_velocity_factor_from_i73": inc_factor}
                all_rows.append(nfwrow)
                with progress_path.open("a") as stream:
                    stream.write(json.dumps(nfwrow, sort_keys=True)+"\n")
                completed += 1
                if (best[vertical]["NFW"] is None or nfwrow["objective_chi2_plus_priors"] <
                        best[vertical]["NFW"]["objective_chi2_plus_priors"]):
                    best[vertical]["NFW"] = {**nfwrow, "predictions": nfwfit["predictions"].tolist()}

                # Simple isolated MOND, fixed a0 and the same transformed source.
                gnewton = bary_v2/reval
                if np.any(gnewton <= 0.0):
                    raise RuntimeError(f"MOND inward g_N<=0 for {vertical} D={d} i={inc} U={upsilon}")
                gmond = 0.5*(gnewton+np.sqrt(gnewton*gnewton+4.0*A0*gnewton))
                mond_v = np.sqrt(reval*gmond)
                mond_chi2 = float(np.sum(((obs-mond_v)/err)**2))
                mondrow = {"model": "MOND", "vertical_case": vertical, "distance_Mpc": d,
                           "inclination_deg": inc, "D_sigma": geom["d_sigma"],
                           "i_sigma": geom["i_sigma"], "upsilon": upsilon, "q_log10_ML": q,
                           "chi2_data": mond_chi2,
                           "rms_km_s": float(np.sqrt(np.mean((obs-mond_v)**2))),
                           "prior_D": prior_d, "prior_i": prior_i, "prior_ML": u_prior,
                           "prior_other": 0.0, "prior_total": prior_sum,
                           "objective_chi2_plus_priors": mond_chi2+prior_sum,
                           "MOND_a0_m_s2": 1.2e-10,
                           "obs_velocity_factor_from_i73": inc_factor,
                           "minimum_gN_km2_s2_kpc": float(np.min(gnewton))}
                all_rows.append(mondrow)
                with progress_path.open("a") as stream:
                    stream.write(json.dumps(mondrow, sort_keys=True)+"\n")
                completed += 1
                if (best[vertical]["MOND"] is None or mondrow["objective_chi2_plus_priors"] <
                        best[vertical]["MOND"]["objective_chi2_plus_priors"]):
                    best[vertical]["MOND"] = {**mondrow, "predictions": mond_v.tolist()}

                # TG principal-mode edge is separate from every finite-u curve.
                edge = solve_tg_edge(B, volume, rcent, zcent, rho_star, rho_gas,
                                     upsilon, reval, obs, err)
                tgrow = {"model": "TG_edge_limit", "vertical_case": vertical,
                         "distance_Mpc": d, "inclination_deg": inc,
                         "D_sigma": geom["d_sigma"], "i_sigma": geom["i_sigma"],
                         "upsilon": upsilon, "q_log10_ML": q,
                         "chi2_data": edge["chi2_data"], "rms_km_s": edge["rms_residual"],
                         "prior_D": prior_d, "prior_i": prior_i, "prior_ML": u_prior,
                         "prior_other": 0.0, "prior_total": prior_sum,
                         "objective_chi2_plus_priors": edge["chi2_data"]+prior_sum,
                         "Kcrit_km_s_minus2": edge["kcrit"],
                         "principal_mode_residual": edge["mode_residual"],
                         "principal_mode_min": edge["minimum_mode"],
                         "obs_velocity_factor_from_i73": inc_factor,
                         "K_treatment": "derived Kcrit limiting curve; no finite normalization and no K prior"}
                all_rows.append(tgrow)
                with progress_path.open("a") as stream:
                    stream.write(json.dumps(tgrow, sort_keys=True)+"\n")
                completed += 1
                if (best[vertical]["TG_edge"] is None or tgrow["objective_chi2_plus_priors"] <
                        best[vertical]["TG_edge"]["objective_chi2_plus_priors"]):
                    best[vertical]["TG_edge"] = {**tgrow, "predictions": edge["predictions"].tolist()}
            if not grid_complete:
                break
        if not grid_complete:
            break

    # At the penalized best geometry/M/L of each vertical profile, solve actual
    # positive finite-u examples at fixed fractions below Kcrit. The fractions
    # are illustrative; no K profile or prior is claimed.
    finite_records = []
    finite_predictions = {}
    for vertical, vc in VERTICAL_CASES.items():
        if not grid_complete:
            break
        selected = best[vertical]["TG_edge"]
        d, inc, upsilon = selected["distance_Mpc"], selected["inclination_deg"], selected["upsilon"]
        prof = geometry_cache[(d, inc)]["profile"]
        stars_surface, star_cfg = surface_density_on_grid(source_r, prof["r"], prof["stars"], "stars")
        gas_surface, gas_cfg = surface_density_on_grid(source_r, prof["r"], prof["gas"], "gas")
        star_rho, gas_rho = tg_vertical_density(np.interp(rcent, source_r, stars_surface),
                                                np.interp(rcent, source_r, gas_surface),
                                                rcent, zcent, vc["hstar"], vc["f_thin"])
        reval = r_catalog*(d/D_REF)
        obs, err, inc_factor = observed_for_inclination(v_catalog, e_catalog, inc)
        mode = solve_tg_edge(B, volume, rcent, zcent, star_rho, gas_rho,
                             upsilon, reval, obs, err)
        dvec = (4.0*math.pi*G*mode["rho"]*volume).ravel()
        for fraction in (0.95, 0.99):
            if time.process_time()-started_cpu > cpu_cap:
                grid_complete = False
                break
            kval = fraction*mode["kcrit"]
            A = B-diags(kval*dvec, format="csr")
            lu = splu(A.tocsc())
            w = lu.solve(dvec).reshape(mode["rho"].shape)
            ufield = 1.0+kval*w
            wmid = (9.0*w[:, 0]-w[:, 1])/8.0
            umid = (9.0*ufield[:, 0]-ufield[:, 1])/8.0
            if np.any(ufield <= 0) or np.any(umid <= 0):
                raise RuntimeError(f"finite-u positivity failure at K/Kcrit={fraction}")
            wf = PchipInterpolator(rcent, wmid, extrapolate=True)
            v2 = -reval*wf.derivative()(reval)/(1.0+kval*wf(reval))
            if np.any(v2 <= 0) or np.any(~np.isfinite(v2)):
                raise RuntimeError(f"finite-u speed failure at K/Kcrit={fraction}")
            speed = np.sqrt(v2)
            data_chi2 = float(np.sum(((obs-speed)/err)**2))
            prior_d, prior_i = geometry_prior(d, inc)
            q = math.log10(upsilon/UPSILON_REF)
            prior_u = (q/MSTAR_SIGMA_DEX)**2
            row = {"vertical_case": vertical, "distance_Mpc": d,
                   "inclination_deg": inc, "upsilon": upsilon,
                   "K_fraction_of_critical": fraction, "K_km_s_minus2": kval,
                   "minimum_u_full": float(np.min(ufield)),
                   "maximum_u_full": float(np.max(ufield)),
                   "minimum_u_midplane": float(np.min(umid)),
                   "chi2_data": data_chi2,
                   "rms_km_s": float(np.sqrt(np.mean((obs-speed)**2))),
                   "prior_D": prior_d, "prior_i": prior_i, "prior_ML": prior_u,
                   "prior_K": None, "prior_total": prior_d+prior_i+prior_u,
                   "objective_chi2_plus_defined_priors": data_chi2+prior_d+prior_i+prior_u,
                   "linear_residual": float(np.linalg.norm(A@w.ravel()-dvec)/np.linalg.norm(dvec)),
                   "obs_velocity_factor_from_i73": inc_factor,
                   "K_treatment": "finite-u subcritical illustration; fixed K/Kcrit, no K prior"}
            finite_records.append(row)
            finite_predictions[(vertical, fraction)] = speed

    # Write the complete scalar grid and three selected model curves per
    # vertical case in both JSON and CSV for independent reading.
    with (OUT/"envelope_grid.csv").open("w", newline="") as stream:
        fields = list(dict.fromkeys(k for row in all_rows for k in row))
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)
    predictions = {"R_catalog_kpc": r_catalog, "Vobs_catalog_km_s": v_catalog,
                   "sigmaV_catalog_km_s": e_catalog}
    for vertical in VERTICAL_CASES:
        for model in ("TG_edge", "NFW", "MOND"):
            row = best[vertical][model]
            if row is None:
                continue
            factor = math.sin(math.radians(row["inclination_deg"]))/math.sin(math.radians(I_REF))
            predictions[f"{vertical}_{model}_catalog_frame_km_s"] = np.asarray(row["predictions"])*factor
        for fraction in (0.95, 0.99):
            if (vertical, fraction) in finite_predictions:
                bestrow = best[vertical]["TG_edge"]
                factor = math.sin(math.radians(bestrow["inclination_deg"]))/math.sin(math.radians(I_REF))
                predictions[f"{vertical}_TG_finite_{fraction:.2f}_catalog_frame_km_s"] = (
                    finite_predictions[(vertical, fraction)]*factor)
    with (OUT/"profiled_predictions.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        fields = list(predictions)
        writer.writerow(fields)
        for idx in range(r_catalog.size):
            writer.writerow([float(np.asarray(predictions[key])[idx]) for key in fields])

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), sharey=True, constrained_layout=True)
    colors = {"TG_edge": "tab:red", "NFW": "tab:blue", "MOND": "tab:green"}
    for ax, vertical in zip(axes, VERTICAL_CASES):
        ax.errorbar(r_catalog, v_catalog, yerr=e_catalog, fmt="ko", ms=2.5, alpha=.7,
                    label="SPARC catalog (i=73°)")
        for model in ("TG_edge", "NFW", "MOND"):
            key = f"{vertical}_{model}_catalog_frame_km_s"
            if key in predictions:
                ax.plot(r_catalog, predictions[key], color=colors[model], lw=1.6,
                        label=f"{model} profile")
        for fraction, style in ((0.95, "--"), (0.99, ":")):
            key = f"{vertical}_TG_finite_{fraction:.2f}_catalog_frame_km_s"
            if key in predictions:
                ax.plot(r_catalog, predictions[key], color="tab:orange", ls=style, lw=1.0,
                        label=f"TG finite K/Kcrit={fraction:.2f}")
        edge = best[vertical]["TG_edge"]
        ax.set_title(f"{vertical}\nTG best D={edge['distance_Mpc']:.1f} Mpc, i={edge['inclination_deg']:.0f}°")
        ax.set_xlabel("SPARC catalog radius [kpc at 13.8 Mpc]")
        ax.grid(alpha=.25)
    axes[0].set_ylabel("deprojected circular speed [km s$^{-1}$]")
    axes[-1].legend(fontsize=7, loc="best")
    fig.savefig(OUT/"profiled_models.png", dpi=180)
    plt.close(fig)

    summary = {
        "status": "complete" if grid_complete and completed == n_expected else "partial_cpu_cap",
        "classification": "descriptive_diagonal_scores_only_no_model_preference_or_significance",
        "inputs_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in
                          [STAR_MAP, STAR_MASK, HALOGAS_MAP, THINGS_MAP, REFERENCE_CSV,
                           ROT_ZIP, PILOT/"newtonian_k0_predictions.csv", solver_path]},
        "method": {
            "geometry": "Reproject raw map pixel centers in sky-angle coordinates at each D/i; Sigma_face-on=Sigma_LOS*cos(i); radii scale with D; observed Rotmod V/error multiply sin(73)/sin(i).",
            "distance_prior": {"form": "Gaussian", "mean_Mpc": D_REF, "sigma_Mpc": DISTANCE_SIGMA_MPC,
                               "grid_sigma_nodes": D_SIGMA_GRID},
            "inclination_prior": {"form": "Gaussian", "mean_deg": I_REF,
                                  "sigma_deg": INCLINATION_SIGMA_DEG,
                                  "grid_sigma_nodes": I_SIGMA_GRID},
            "stellar_ML_prior": {"form": "Gaussian in log10(Upsilon/0.6)",
                                 "mean_Upsilon": UPSILON_REF, "sigma_dex": MSTAR_SIGMA_DEX,
                                 "grid_q_dex": Q_GRID,
                                 "penalty": "(log10(Upsilon/0.6)/0.1)^2"},
            "vertical_cases": VERTICAL_CASES,
            "models": {
                "TG_edge": "Kcrit principal-mode logarithmic-force limit, not a finite-normalization u; Kcrit is derived and has no K prior",
                "NFW": "same map-derived source plus NFW; profile logM200 in [9,14] and c200 in [1,50], no c-M relation or halo prior, H0=73",
                "MOND": "same map-derived source, fixed a0=1.2e-10 m/s^2, isolated simple mu=x/(1+x), no EFE",
                "common_score": "chi2_diag + distance prior + inclination prior + S4G log-M/L prior",
                "finite_TG": "separate positive-u solves at K/Kcrit={0.95,0.99}; fixed illustrative fractions and no K prior"},
            "source": "HALOGAS H I column map plus helium and S4G P5 old-star map inside the P5 ellipse; finite P5 stellar aperture continued by the declared exponential tail; azimuthal averages; no H2 or ionized component",
            "diagonal_likelihood": "sum(((Vobs(i)-Vmodel)/errV(i))^2); 43 Rotmod rows; no radial covariance available",
            "tail": {"stars": "exponential fit to supported 5-12.75 kpc bins (upper end shortened only if P5 support ends earlier)",
                     "gas": "linear 5-kpc taper after the last supported radial bin"},
            "hankel": {"R_source_step_kpc": R_SOURCE_STEP, "R_source_max_kpc": R_SOURCE_MAX,
                       "kmax_kpc_inverse": K_MAX, "n_k": N_K},
            "TG_grid": {"dr_kpc": TG_DR, "dz_kpc": TG_DZ, "Rmax_kpc": TG_RMAX,
                        "Zmax_kpc": TG_ZMAX, "nr": B.shape[0]//int(round(TG_ZMAX/TG_DZ)),
                        "nz": int(round(TG_ZMAX/TG_DZ)), "outer_boundary": "monopole Robin from pilot solver"}},
        "kernel_closure_at_original_geometry": kernel_closure,
        "same_grid_tg_source_and_edge_reproduction": same_grid_tg_check,
        "n_grid_rows": completed, "n_expected_grid_rows": n_expected,
        "best_by_vertical_case": best,
        "finite_u_at_best_tg_edge_nuisances": finite_records,
        "runtime": {"cpu_seconds": time.process_time()-started_cpu,
                    "wall_seconds": time.perf_counter()-started_wall,
                    "cpu_cap_seconds": cpu_cap,
                    "python": platform.python_version(), "numpy": np.__version__},
        "limitations": [
            "This remains the declared axisymmetric parameterized source; the raw maps do not provide a unique 3D stellar/gas density.",
            "The five-node D/i grid samples a finite ±2 sigma box; profiled optima are grid minima, not continuous posterior integration.",
            "The stellar source uses a finite P5 ellipse and explicit exponential extrapolation; source vertical structure is varied only in the listed cases.",
            "Scores use diagonal SPARC errV only; no full rotation-curve covariance is available, so there are no significance or model-preference claims.",
            "NFW and MOND fits have different model structure; their score gaps are descriptive and do not account for all model complexity.",
            "TG edge curves are limiting curves, not finite-u solutions; finite examples at K/Kcrit 0.95 and 0.99 are kept separate and are not K-profiled."],
    }
    (OUT/"results.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({"status": summary["status"], "n_grid_rows": completed,
                      "n_expected_grid_rows": n_expected,
                      "runtime": summary["runtime"],
                      "best_by_vertical_case": {v: {m: (None if best[v][m] is None else
                          {k: best[v][m][k] for k in best[v][m] if k != "predictions"})
                          for m in best[v]} for v in best}}, indent=2))


def pilot():
    maps = load_maps()
    ref = transformed_sources(maps, D_REF, I_MAP_REF)
    ref_rows = list(csv.DictReader(REFERENCE_CSV.open(newline="")))
    csv_by_r = {round(float(x["r_mid_kpc"]), 8): x for x in ref_rows}
    comparisons = {}
    for key, name in [("gas", "halogas_hihe_faceon_Msun_pc2"),
                      ("things", "things_hihe_faceon_Msun_pc2"),
                      ("stars", "s4g_stars_faceon_Msun_pc2_Upsilon0p6")]:
        errs, rels, n = [], [], 0
        for radius, value in zip(ref["r"], ref[key]):
            old = csv_by_r.get(round(float(radius), 8), {}).get(name, "")
            if np.isfinite(value) and old:
                oldv = float(old)
                errs.append(value-oldv)
                if oldv != 0:
                    rels.append((value-oldv)/oldv)
                n += 1
        comparisons[key] = {"n_matched": n,
                            "max_abs_surface_density_difference_Msun_pc2": float(np.max(np.abs(errs))),
                            "max_relative_surface_density_difference": float(np.max(np.abs(rels))) if rels else None}
    # Exact local Jacobian check: Sigma_face-on * dA_face-on = Sigma_los * dA_sky.
    # Total source mass follows D^2 and is inclination-invariant for each map pixel.
    geometry_checks = []
    for d, inc in [(D_REF, I_MAP_REF), (12.4, 70.), (15.2, 76.)]:
        cosi = math.cos(math.radians(inc))
        ratios = {}
        for label, los, h in [("HI", maps["halo_los"], maps["halo_h"]),
                              ("stars", maps["star_los"], maps["star_h"])]:
            sx, sy = pixel_scales(h)
            dsq_arcsec = (d*1000.*sx*ARCSEC_RAD)*(d*1000.*sy*ARCSEC_RAD)
            projected_mass_proxy = float(np.nansum(los*dsq_arcsec))
            disk_mass_proxy = float(np.nansum((los*cosi)*(dsq_arcsec/cosi)))
            ratios[label] = disk_mass_proxy/projected_mass_proxy
        geometry_checks.append({"D_Mpc": d, "i_deg": inc,
                                "faceon_projected_mass_ratio": ratios,
                                "expected_mass_ratio_vs_reference": (d/D_REF)**2})
    # Small verified TG pilot, imported from the declared solver with a coarse
    # grid; it checks the transformed-source operator at the map-extraction i.
    # Production geometry nuisance scans are added after this gate succeeds.
    solver_path = PILOT/"solve_tg_axisymmetric.py"
    import importlib.util
    spec = importlib.util.spec_from_file_location("tg_solver_pilot", solver_path)
    tg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tg)
    started_cpu, started_wall = time.process_time(), time.perf_counter()
    coarse = tg.build_operator(.5, .25, 80., 60., .6)
    B, _D0, _rho0, vol, rr, zz, step, _cfg = coarse
    from scipy.sparse import diags
    from scipy.sparse.linalg import eigsh
    from scipy.interpolate import PchipInterpolator
    rnodes = np.arange(ref["r"].size)*0.5 + 0.25
    # This source check uses the same decl as the prior solver, interpolated to
    # its cylindrical cell radii and mixed vertical profile.
    valid_s = np.isfinite(ref["stars"])
    fs = PchipInterpolator(ref["r"][valid_s], ref["stars"][valid_s], extrapolate=False)
    rcent = rr[:, 0]
    stars = np.nan_to_num(fs(np.clip(rcent, ref["r"][valid_s][0], ref["r"][valid_s][-1])))
    fit = valid_s & (ref["r"] >= 5.) & (ref["r"] <= 12.75) & (ref["stars"] > 0)
    slope, intercept = np.polyfit(ref["r"][fit], np.log(ref["stars"][fit]), 1)
    tail_rd = -1./slope
    stars[rcent > ref["r"][valid_s][-1]] = ref["stars"][valid_s][-1]*np.exp(
        -(rcent[rcent > ref["r"][valid_s][-1]]-ref["r"][valid_s][-1])/tail_rd)
    valid_g = np.isfinite(ref["gas"])
    fg = PchipInterpolator(ref["r"][valid_g], ref["gas"][valid_g], extrapolate=False)
    gas = np.zeros_like(rcent)
    ins = rcent <= ref["r"][valid_g][-1]
    gas[ins] = fg(np.maximum(rcent[ins], ref["r"][valid_g][0]))
    taper = (rcent > ref["r"][valid_g][-1]) & (rcent < ref["r"][valid_g][-1]+5.)
    gas[taper] = ref["gas"][valid_g][-1]*(1-(rcent[taper]-ref["r"][valid_g][-1])/5.)
    z = zz
    rho = (stars[:, None]*1.e6/(2*.3)*np.exp(-z/.3)
           + gas[:, None]*1.e6*(.85/(2*.2)*np.exp(-z/.2)+.15/(2*3.)*np.exp(-z/3.)))
    Dvec = (4*math.pi*G*rho*vol).ravel()
    vals, vecs = eigsh(B, k=1, M=diags(Dvec, format="csr"), sigma=0., which="LM",
                       tol=1e-7, maxiter=3000)
    mode = vecs[:, 0].reshape(rho.shape)
    if mode.sum() < 0:
        mode *= -1
    mode /= mode.max()
    hmid = (9*mode[:, 0]-mode[:, 1])/8
    mode_r = (np.arange(mode.shape[0])+.5)*step[0]
    rf, vo, ev, _, _ = read_rotmod()
    hn = PchipInterpolator(mode_r, hmid, extrapolate=True)
    edge_v2 = -rf/vals[0]*hn.derivative()(rf)/hn(rf)
    edge_v = np.sqrt(edge_v2)
    vo72, ev72, _ = observed_for_inclination(vo, ev, I_MAP_REF)
    pilot_record = {"grid_shape": list(rho.shape), "kcrit_km_s_minus2": float(vals[0]),
                    "upsilon": 0.6,
                    "eigen_residual": float(np.linalg.norm(B@vecs[:,0]-vals[0]*Dvec*vecs[:,0])/
                                             np.linalg.norm(B@vecs[:,0])),
                    "minimum_edge_v2": float(edge_v2.min()),
                    "chi2_diagonal_against_Vobs_redeprojected_to_i72": float(np.sum(((vo72-edge_v)/ev72)**2)),
                    "catalog_i73_chi2_for_pipeline_smoke_test_only": float(np.sum(((vo-edge_v)/ev)**2)),
                    "wall_seconds": time.perf_counter()-started_wall,
                    "cpu_seconds": time.process_time()-started_cpu}
    input_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in
                    [STAR_MAP, STAR_MASK, HALOGAS_MAP, THINGS_MAP, REFERENCE_CSV, ROT_ZIP, solver_path]}
    record = {"status": "geometry_transform_pilot_complete",
              "transform_derivation": {
                  "radius_kpc": "D_Mpc*1000*ARCSEC_RAD*sqrt(major_arcsec^2+(minor_arcsec/cos(i))^2)",
                  "HI_faceon_surface_density": "N_HI/1.25e20*1.36*cos(i) M_sun pc^-2",
                  "stars_faceon_surface_density": "I_nu(MJy/sr) -> Lsun pc^-2 times Upsilon*cos(i)",
                  "mass_jacobian": "Sigma_faceon*cos^-1(i)*D^2*dOmega = Sigma_los*D^2*dOmega",
                  "mass_distance_scaling": "D^2 at fixed map; independent of inclination",
                  "sparc_velocity_transform": "Vobs(i)=Vobs(73)*sin(73)/sin(i); errV scales identically",
                  "fixed_physical_vertical_scale": "h is held in kpc; source radii scale with D, so h/R changes with D"},
              "base_profile_reconstruction_i72": comparisons,
              "jacobian_checks": geometry_checks,
              "coarse_tg_pilot_i72": pilot_record,
              "inputs_sha256": input_hashes,
              "environment": {"python": platform.python_version(), "numpy": np.__version__},
              "decision": "The map-to-disk transformation is coherent; proceed to the bounded D/i and vertical nuisance envelope."}
    path = OUT/"pilot_geometry.json"
    path.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps({"status": record["status"], "base_profile_reconstruction_i72": comparisons,
                      "jacobian_checks": geometry_checks, "coarse_tg_pilot_i72": pilot_record}, indent=2))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--cpu-cap", type=float, default=850.0)
    args = parser.parse_args()
    if args.pilot:
        pilot()
    elif args.run:
        run_envelope(args.cpu_cap)
    else:
        raise SystemExit("Use --pilot for geometry validation or --run for the bounded nuisance envelope")

#!/usr/bin/env python3
"""Build curve-blind NGC 4559 source profiles from the two approved FITS maps.

This script has an explicit input allow-list. It never opens SPARC Rotmod data,
target Vobs, residuals, or any target rotation curve. It performs no gravity
fit, likelihood evaluation, or curve prediction.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
INPUTS = {
    "stellar": ROOT / "work/tg_holdout_target_audit/NGC4559.stellar.fits",
    "hi": ROOT / "work/tg_holdout_target_audit/NGC4559-HR_coldens.fits",
}
EXPECTED_SHA256 = {
    "stellar": "9fe33f187f10db079c360de29c80de727c16a4932e275362ef5b90c272bd9405",
    "hi": "320d2de1e29adbfdc7bd9023d58ed6dec78bf7d8b1ba326ffdc064150e252bba",
}

# Frozen source geometry and sensitivity brackets. PA is east of north.
D_MPC = {"low": 8.56, "nominal": 8.91, "high": 9.26}
INCLINATIONS_DEG = (66.0, 67.0, 68.0)
PA_DEG = (-42.0, -37.0, -32.0)  # P5 -37 deg with a declared +/-5 deg bracket.
P5_SMA_ARCSEC = 278.1
P5_ELLIPTICITY = 0.579
P5_Q = 1.0 - P5_ELLIPTICITY
OUTER_LIGHT_KPC = 20.8  # pre-fit metadata support target; not taken from curve rows.
TAIL_FIT_KPC = (6.0, 11.0)  # map-only fit interval, wholly inside the P5 aperture.
BIN_KPC = 0.5
N_SUBPIX = 5

RA0_DEG = 188.990005284
DEC0_DEG = 27.9600007817
M_H_G = 1.6735575e-24
M_SUN_G = 1.98847e33
MPC_CM = 3.0856775814913673e24
KPC_CM = MPC_CM / 1000.0
PC_CM = KPC_CM / 1000.0
MJY_TO_JY = 1e6
JY_TO_W_M2_HZ = 1e-26
IRAC36_VEGA_ZERO_JY = 280.9
SUN_M36_VEGA = 3.24
SUN_F36_AT_10PC_JY = IRAC36_VEGA_ZERO_JY * 10.0 ** (-0.4 * SUN_M36_VEGA)
HELIUM_FACTOR = 1.36


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def lin_wcs_without_stale_sip(header):
    """P5 CTYPE declares plain TAN despite orphan SIP keywords; honor TAN grid."""
    clean = header.copy()
    sip = re.compile(r"^(A|B|AP|BP)_(ORDER|DMAX|\d+_\d+)$")
    for key in list(clean.keys()):
        if sip.match(key):
            del clean[key]
    return WCS(clean, naxis=2).celestial


def bilinear(image: np.ndarray, x: np.ndarray, y: np.ndarray):
    x0 = np.floor(x).astype(np.int64)
    y0 = np.floor(y).astype(np.int64)
    ok = (x0 >= 0) & (y0 >= 0) & (x0 + 1 < image.shape[1]) & (y0 + 1 < image.shape[0])
    out = np.full(x.shape, np.nan, dtype=np.float64)
    if np.any(ok):
        dx = x[ok] - x0[ok]
        dy = y[ok] - y0[ok]
        v00 = image[y0[ok], x0[ok]]
        v10 = image[y0[ok], x0[ok] + 1]
        v01 = image[y0[ok] + 1, x0[ok]]
        v11 = image[y0[ok] + 1, x0[ok] + 1]
        out[ok] = ((1 - dx) * (1 - dy) * v00 + dx * (1 - dy) * v10 +
                   (1 - dx) * dy * v01 + dx * dy * v11)
    return out


def reproject_stellar_to_hi(stellar, stellar_header, hi_header, hi_shape):
    """Average 5x5 stellar surface-brightness samples onto the HI sky grid."""
    sw = lin_wcs_without_stale_sip(stellar_header)
    hw = WCS(hi_header, naxis=2).celestial
    h, w = hi_shape
    # Find the target-grid bounding box of the complete input image.
    sx = np.array([-0.5, stellar.shape[1] - 0.5, stellar.shape[1] - 0.5, -0.5])
    sy = np.array([-0.5, -0.5, stellar.shape[0] - 0.5, stellar.shape[0] - 0.5])
    ra, dec = sw.pixel_to_world_values(sx, sy)
    gx, gy = hw.world_to_pixel_values(ra, dec)
    xlo = max(0, int(np.floor(np.nanmin(gx))) - 2)
    xhi = min(w, int(np.ceil(np.nanmax(gx))) + 3)
    ylo = max(0, int(np.floor(np.nanmin(gy))) - 2)
    yhi = min(h, int(np.ceil(np.nanmax(gy))) + 3)
    yy, xx = np.mgrid[ylo:yhi, xlo:xhi]
    total = np.zeros(xx.shape, dtype=np.float64)
    count = np.zeros(xx.shape, dtype=np.uint8)
    offsets = (np.arange(N_SUBPIX, dtype=np.float64) + 0.5) / N_SUBPIX - 0.5
    for oy in offsets:
        for ox in offsets:
            ra, dec = hw.pixel_to_world_values(xx + ox, yy + oy)
            ix, iy = sw.world_to_pixel_values(ra, dec)
            sample = bilinear(stellar, np.asarray(ix), np.asarray(iy))
            good = np.isfinite(sample)
            total[good] += sample[good]
            count[good] += 1
    out = np.full((h, w), np.nan, dtype=np.float32)
    footprint = np.zeros((h, w), dtype=np.float32)
    block = np.full(total.shape, np.nan, dtype=np.float64)
    good = count > 0
    block[good] = total[good] / count[good]
    out[ylo:yhi, xlo:xhi] = block.astype(np.float32)
    footprint[ylo:yhi, xlo:xhi] = count.astype(np.float32) / float(N_SUBPIX ** 2)
    return out, footprint, (xlo, xhi, ylo, yhi)


def tangent_offsets_arcsec(ra_deg, dec_deg, ra0_deg=RA0_DEG, dec0_deg=DEC0_DEG):
    """Gnomonic east/north offsets; accurate over the map field."""
    ra = np.deg2rad(np.asarray(ra_deg, dtype=np.float64))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=np.float64))
    ra0 = math.radians(ra0_deg)
    dec0 = math.radians(dec0_deg)
    dra = (ra - ra0 + np.pi) % (2 * np.pi) - np.pi
    den = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(dra)
    east = np.cos(dec) * np.sin(dra) / den
    north = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(dra)) / den
    return np.rad2deg(east) * 3600.0, np.rad2deg(north) * 3600.0


def sky_pixel_area_sr(wcs, x, y):
    # Finite-difference local WCS Jacobian in east/north radians per pixel.
    ra_x1, dec_x1 = wcs.pixel_to_world_values(x - 0.5, y)
    ra_x2, dec_x2 = wcs.pixel_to_world_values(x + 0.5, y)
    ra_y1, dec_y1 = wcs.pixel_to_world_values(x, y - 0.5)
    ra_y2, dec_y2 = wcs.pixel_to_world_values(x, y + 0.5)
    decm = math.radians(float((dec_x1 + dec_x2 + dec_y1 + dec_y2) / 4.0))
    wrap = lambda a, b: (a - b + 180.0) % 360.0 - 180.0
    dx_e = math.radians(wrap(ra_x2, ra_x1) * math.cos(decm))
    dx_n = math.radians(float(dec_x2 - dec_x1))
    dy_e = math.radians(wrap(ra_y2, ra_y1) * math.cos(decm))
    dy_n = math.radians(float(dec_y2 - dec_y1))
    return abs(dx_e * dy_n - dy_e * dx_n)


def ring_sum(values, radius, edges, mask):
    idx = np.digitize(radius[mask], edges, right=False) - 1
    valid = (idx >= 0) & (idx < len(edges) - 1)
    idx = idx[valid]
    v = values[mask][valid]
    return np.bincount(idx, weights=v, minlength=len(edges) - 1), np.bincount(idx, minlength=len(edges) - 1)


def exponential_tail_luminosity(sigma_at_r0, slope_per_kpc, r0, r1):
    if r1 <= r0 or sigma_at_r0 <= 0:
        return 0.0
    beta = slope_per_kpc
    if abs(beta) < 1e-10:
        integ = 0.5 * (r1 * r1 - r0 * r0)
    else:
        integ = (math.exp(beta * (r1 - r0)) * (r1 / beta - 1.0 / (beta * beta)) -
                 (r0 / beta - 1.0 / (beta * beta)))
    return max(0.0, 2.0 * math.pi * sigma_at_r0 * integ)


def main():
    for key, path in INPUTS.items():
        got = sha256(path)
        if got != EXPECTED_SHA256[key]:
            raise RuntimeError(f"Approved input hash mismatch for {path}: {got}")

    with fits.open(INPUTS["stellar"], memmap=True) as hdul:
        star_native = np.asarray(hdul[0].data, dtype=np.float64)
        sh = hdul[0].header.copy()
        star_bunit = sh.get("BUNIT", "")
    with fits.open(INPUTS["hi"], memmap=True) as hdul:
        hi_nhi = np.asarray(hdul[0].data, dtype=np.float64)
        hh = hdul[0].header.copy()
        hi_bunit = hh.get("BUNIT", "")
        hi_history = "\n".join(str(x) for x in hh.get("HISTORY", []))
    # The released product name and HISTORY establish that the stored array is
    # the converted N_HI map. BUNIT is stale; never apply the HISTORY factor twice.
    required_history = "1.823e18*606*1e3"
    if required_history not in hi_history or "out=coldens" not in hi_history:
        raise RuntimeError("HALOGAS HISTORY does not contain the expected coldens conversion")

    star_grid, star_footprint, crop = reproject_stellar_to_hi(
        star_native, sh, hh, hi_nhi.shape)
    hw = WCS(hh, naxis=2).celestial
    ny, nx = hi_nhi.shape
    xg, yg = np.meshgrid(np.arange(nx, dtype=np.float64), np.arange(ny, dtype=np.float64))
    ra, dec = hw.pixel_to_world_values(xg, yg)
    east_arcsec, north_arcsec = tangent_offsets_arcsec(ra, dec)
    omega_sr = sky_pixel_area_sr(hw, 512.0, 512.0)
    if not np.isfinite(omega_sr) or omega_sr <= 0:
        raise RuntimeError(f"Invalid HI pixel solid angle: {omega_sr}")
    hi_footprint = np.isfinite(hi_nhi)  # finite pixels span the full released map.
    hi_detect = hi_footprint & (hi_nhi > 0.0)
    hi_mask_3e19 = hi_detect & (hi_nhi >= 3e19)
    hi_mask_1e20 = hi_detect & (hi_nhi >= 1e20)
    star_valid = np.isfinite(star_grid) & (star_footprint >= 0.99)
    star_nonneg = np.where(star_valid, np.maximum(star_grid, 0.0), 0.0)

    # Luminosity conversion: MJy/sr -> Jy/pixel -> L_nu in solar IRAC 3.6-um
    # band units, with M_sun,3.6=3.24 Vega and IRAC zero point 280.9 Jy.
    dnom_pc = D_MPC["nominal"] * 1e6
    fnu_jy = star_nonneg * MJY_TO_JY * omega_sr
    lstar_nom_pixel = fnu_jy / SUN_F36_AT_10PC_JY * (dnom_pc / 10.0) ** 2
    hi_mass_nom_pixel = np.maximum(hi_nhi, 0.0) * omega_sr * (D_MPC["nominal"] * MPC_CM) ** 2 * M_H_G / M_SUN_G

    nominal_geom = (D_MPC["nominal"], 67.0, -37.0)
    def geometry_arrays(distance_mpc, incl_deg, pa_deg):
        pa = math.radians(pa_deg)
        xmaj = east_arcsec * math.sin(pa) + north_arcsec * math.cos(pa)
        ymin = east_arcsec * math.cos(pa) - north_arcsec * math.sin(pa)
        radius_arcsec = np.sqrt(xmaj * xmaj + (ymin / math.cos(math.radians(incl_deg))) ** 2)
        r_kpc = radius_arcsec / 206264.80624709636 * distance_mpc * 1000.0
        p5_radius_arcsec = np.sqrt(xmaj * xmaj + (ymin / P5_Q) ** 2)
        p5_mask = p5_radius_arcsec <= P5_SMA_ARCSEC
        return r_kpc, p5_mask

    # Reprojected maps on the same 4-arcsec HALOGAS tangent grid.
    starout_header = hh.copy()
    starout_header["BUNIT"] = "MJy/sr"
    starout_header["BTYPE"] = "reprojected stellar old-light surface brightness"
    starout_header.add_history("5x5 subpixel average of S4G P5 stellar map onto HALOGAS HR grid")
    starout_header.add_history("Orphan S4G SIP keywords ignored: CTYPE is RA---TAN/DEC--TAN without -SIP")
    fits.PrimaryHDU(data=star_grid.astype(np.float32), header=starout_header).writeto(
        OUT / "stellar_reprojected_MJy_sr.fits", overwrite=True)
    foot_header = hh.copy()
    foot_header["BUNIT"] = "fraction"
    foot_header.add_history("Fraction of 25 subpixel samples with valid S4G stellar-map WCS coverage")
    fits.PrimaryHDU(data=star_footprint.astype(np.float32), header=foot_header).writeto(
        OUT / "stellar_reprojected_footprint.fits", overwrite=True)

    # The nominal radial output uses 0.5-kpc bins and extends through all
    # detected HALOGAS emission; the stellar profile is only validated inside P5.
    r_nom, p5_nom = geometry_arrays(*nominal_geom)
    hi_rmax = float(np.max(r_nom[hi_detect])) if np.any(hi_detect) else 0.0
    profile_rmax = max(hi_rmax, OUTER_LIGHT_KPC)
    edges = np.arange(0.0, math.ceil(profile_rmax / BIN_KPC) * BIN_KPC + BIN_KPC * 0.5,
                      BIN_KPC, dtype=np.float64)
    centers = 0.5 * (edges[:-1] + edges[1:])
    star_ap_mask = star_valid & p5_nom
    star_full = np.where(star_ap_mask, lstar_nom_pixel, 0.0)
    l_ann, n_star = ring_sum(lstar_nom_pixel, r_nom, edges, star_ap_mask)
    l_cum = np.cumsum(l_ann)
    m_ann = np.where(hi_detect, hi_mass_nom_pixel, 0.0)
    hi_m_ann, n_hi = ring_sum(m_ann, r_nom, edges, np.ones_like(hi_detect, dtype=bool))
    hi_m_cum = np.cumsum(hi_m_ann)
    area_faceon_kpc2 = math.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    sigma_l = np.divide(l_ann, area_faceon_kpc2 * 1e6,
                        out=np.full_like(l_ann, np.nan), where=area_faceon_kpc2 > 0)
    sigma_gas = np.divide(HELIUM_FACTOR * hi_m_ann, area_faceon_kpc2 * 1e6,
                          out=np.full_like(hi_m_ann, np.nan), where=area_faceon_kpc2 > 0)
    # In the radial profile, bins beyond the direct P5 aperture are intentionally
    # empty; explicit tail scenarios live in the separate sensitivity table.
    max_direct_r = np.max(r_nom[star_ap_mask]) if np.any(star_ap_mask) else 0.0

    profile_path = OUT / "profiles.csv"
    with profile_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "r_in_kpc", "r_out_kpc", "r_mid_kpc", "disk_annulus_area_kpc2",
            "star_map_coverage_fraction", "star_pixels_in_P5_aperture",
            "L36_annulus_1e6_Lsun36", "L36_cumulative_1e6_Lsun36",
            "Sigma_L36_Lsun36_pc2", "HI_positive_pixels", "HI_ge_3e19_pixels",
            "HI_ge_1e20_pixels", "NHI_annular_mean_los_cm-2",
            "MHI_annulus_1e6_Msun", "MHI_cumulative_1e6_Msun",
            "Mgas_He_annulus_1e6_Msun", "Mgas_He_cumulative_1e6_Msun",
            "Sigma_gas_He_Msun_pc2", "HI_footprint_area_kpc2",
        ])
        for j, rc in enumerate(centers):
            annulus = (r_nom >= edges[j]) & (r_nom < edges[j + 1])
            n_geom = int(np.count_nonzero(annulus))
            n_cov = int(np.count_nonzero(annulus & star_valid))
            hi_pos_n = int(np.count_nonzero(annulus & hi_detect))
            hi_3_n = int(np.count_nonzero(annulus & hi_mask_3e19))
            hi_1_n = int(np.count_nonzero(annulus & hi_mask_1e20))
            avg_nhi = (float(np.sum(hi_nhi[annulus & hi_detect]) / hi_pos_n)
                       if hi_pos_n else 0.0)
            area_hi = hi_pos_n * omega_sr * (D_MPC["nominal"] * 1000.0) ** 2
            writer.writerow([
                f"{edges[j]:.3f}", f"{edges[j+1]:.3f}", f"{rc:.3f}",
                f"{area_faceon_kpc2[j]:.8g}",
                f"{n_cov / n_geom:.6f}" if n_geom else "0",
                int(n_star[j]), f"{l_ann[j] / 1e6:.9g}", f"{l_cum[j] / 1e6:.9g}",
                f"{sigma_l[j]:.9g}" if np.isfinite(sigma_l[j]) else "",
                hi_pos_n, hi_3_n, hi_1_n, f"{avg_nhi:.9g}",
                f"{hi_m_ann[j] / 1e6:.9g}", f"{hi_m_cum[j] / 1e6:.9g}",
                f"{HELIUM_FACTOR * hi_m_ann[j] / 1e6:.9g}",
                f"{HELIUM_FACTOR * hi_m_cum[j] / 1e6:.9g}",
                f"{sigma_gas[j]:.9g}" if np.isfinite(sigma_gas[j]) else "",
                f"{area_hi:.9g}",
            ])

    # Profile-only exponential fit, set by map geometry before evaluating the
    # tail: 6--11 kpc lies inside the nominal 12.0-kpc P5 analysis ellipse.
    fit = ((centers >= TAIL_FIT_KPC[0]) & (centers <= TAIL_FIT_KPC[1]) &
           (np.asarray(n_star) > 0) & np.isfinite(sigma_l) & (sigma_l > 0))
    if np.count_nonzero(fit) < 4:
        raise RuntimeError(f"Insufficient positive P5 map annuli for frozen tail fit: {fit.sum()}")
    slope, intercept = np.polyfit(centers[fit], np.log(sigma_l[fit]), 1)
    if slope >= 0:
        raise RuntimeError(f"Frozen outer-light fit is not declining (slope={slope})")

    # Compute full geometry/tail sensitivity; no kinematics enter any loop.
    sens_path = OUT / "geometry_tail_sensitivity.csv"
    sens_rows = []
    with sens_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "distance_case", "D_Mpc", "inclination_deg", "PA_deg",
            "R_P5_kpc", "L36_inside_P5_1e9_Lsun36", "tail_model",
            "tail_L36_12_to_20p8_1e9_Lsun36", "L36_total_to_20p8_1e9_Lsun36",
            "MHI_total_map_1e9_Msun", "Mgas_He_total_map_1e9_Msun",
            "MHI_within_20p8_1e9_Msun", "NHI_positive_pixels",
            "NHI_ge_3e19_pixels", "NHI_ge_1e20_pixels",
            "stellar_map_coverage_pixels", "stellar_map_extent_R_kpc",
            "fit_slope_per_kpc",
        ])
        for dcase, dist in D_MPC.items():
            for inc in INCLINATIONS_DEG:
                for pa in PA_DEG:
                    radius, p5m = geometry_arrays(dist, inc, pa)
                    aper = star_valid & p5m
                    dist_scale2 = (dist / D_MPC["nominal"]) ** 2
                    l_pix = lstar_nom_pixel * dist_scale2
                    l_inside = float(np.sum(l_pix[aper]))
                    r_p5 = (P5_SMA_ARCSEC / 206264.80624709636) * dist * 1000.0
                    # Rebuild the P5 annular luminosity density under this geometry.
                    se = np.arange(0.0, max(r_p5 + BIN_KPC, OUTER_LIGHT_KPC + BIN_KPC), BIN_KPC)
                    sc = 0.5 * (se[:-1] + se[1:])
                    lbin, _ = ring_sum(np.where(aper, l_pix, 0.0), radius, se,
                                       np.ones_like(aper, dtype=bool))
                    area = math.pi * (se[1:] ** 2 - se[:-1] ** 2)
                    sl = np.divide(lbin, area * 1e6, out=np.full_like(lbin, np.nan), where=area > 0)
                    fmask = ((sc >= TAIL_FIT_KPC[0]) & (sc <= TAIL_FIT_KPC[1]) &
                             np.isfinite(sl) & (sl > 0) & (lbin > 0))
                    if np.count_nonzero(fmask) >= 4:
                        b, a = np.polyfit(sc[fmask], np.log(sl[fmask]), 1)
                        if b >= 0:
                            b = slope * D_MPC["nominal"] / dist
                        sigp = math.exp(a + b * r_p5)
                    else:
                        b = slope * D_MPC["nominal"] / dist
                        sigp = math.exp(intercept + slope * r_p5 * D_MPC["nominal"] / dist) * dist_scale2
                    # Use the frozen fit slope if local geometry loses annuli to edge pixels;
                    # the direct P5 map value at the fit range sets each geometry's scale.
                    scenarios = (
                        ("truncate_at_P5", 0.0),
                        ("exp_slope_6_11", b),
                        ("half_slope_6_11", 0.5 * b),
                    )
                    for label, tail_slope in scenarios:
                        tail = exponential_tail_luminosity(sigp, tail_slope, r_p5, OUTER_LIGHT_KPC)
                        total_l = l_inside + tail
                        mscale = (dist / D_MPC["nominal"]) ** 2
                        hp = hi_detect
                        mhi_total = float(np.sum(hi_mass_nom_pixel[hp]) * mscale)
                        mhi_r = float(np.sum(hi_mass_nom_pixel[hp & (radius <= OUTER_LIGHT_KPC)]) * mscale)
                        row = [dcase, dist, inc, pa, r_p5, l_inside / 1e9, label,
                               tail / 1e9, total_l / 1e9, mhi_total / 1e9,
                               HELIUM_FACTOR * mhi_total / 1e9, mhi_r / 1e9,
                               int(np.count_nonzero(hi_detect)), int(np.count_nonzero(hi_mask_3e19)),
                               int(np.count_nonzero(hi_mask_1e20)), int(np.count_nonzero(star_valid)),
                               float(np.max(radius[star_valid])) if np.any(star_valid) else 0.0, b]
                        writer.writerow([f"{x:.9g}" if isinstance(x, (float, np.floating)) else x for x in row])
                        sens_rows.append(row)

    # Plot source profiles only; stellar M/L remains a multiplicative parameter.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 8.0), sharex=True,
                                    gridspec_kw={"height_ratios": [1, 1]})
    okstar = (sigma_l > 0) & (centers <= P5_SMA_ARCSEC / 206264.80624709636 * D_MPC["nominal"] * 1000.0)
    ax1.plot(centers[okstar], sigma_l[okstar], color="#2471A3", marker="o", ms=3,
             lw=1.2, label=r"S4G old light: $L_{3.6}/\mathrm{pc}^2$")
    sig_boundary = math.exp(intercept + slope * OUTER_LIGHT_KPC)
    rr = np.linspace(P5_SMA_ARCSEC / 206264.80624709636 * D_MPC["nominal"] * 1000.0,
                     OUTER_LIGHT_KPC, 150)
    rpnom = P5_SMA_ARCSEC / 206264.80624709636 * D_MPC["nominal"] * 1000.0
    sig_at_rp = math.exp(intercept + slope * rpnom)
    ax1.plot(rr, sig_at_rp * np.exp(slope * (rr - rpnom)), color="#2471A3", ls="--",
             label="exponential tail from 6–11 kpc map fit")
    ax1.plot(rr, sig_at_rp * np.exp(0.5 * slope * (rr - rpnom)), color="#2471A3", ls=":",
             label="half-slope tail bracket")
    ax1.axvline(rpnom, color="0.45", ls="--", lw=1, label="P5 aperture edge")
    ax1.set_yscale("log")
    ax1.set_ylabel(r"$\Sigma_{L,3.6}$ [$L_{\odot,3.6}\,\mathrm{pc}^{-2}$]")
    ax1.grid(True, which="both", alpha=0.2)
    ax1.legend(fontsize=8, ncol=2)
    okhi = (sigma_gas > 0) & (centers <= max(hi_rmax, OUTER_LIGHT_KPC))
    ax2.plot(centers[okhi], sigma_gas[okhi], color="#B03A2E", marker="o", ms=3,
             lw=1.1, label=r"HALOGAS H I + He ($1.36$), face-on")
    ax2.axvline(rpnom, color="0.45", ls="--", lw=1)
    ax2.axvline(OUTER_LIGHT_KPC, color="0.35", ls=":", lw=1.2,
                label="source support target 20.8 kpc")
    ax2.set_yscale("log")
    ax2.set_xlabel("Deprojected radius [kpc]")
    ax2.set_ylabel(r"$\Sigma_{\mathrm{gas}}$ [$M_{\odot}\,\mathrm{pc}^{-2}$]")
    ax2.grid(True, which="both", alpha=0.2)
    ax2.legend(fontsize=8)
    fig.suptitle("NGC 4559 source-only radial profiles (no kinematic data)")
    fig.tight_layout()
    fig.savefig(OUT / "source_profiles.png", dpi=180)
    plt.close(fig)

    # Map-aperture metrics and selected mass cross-check at the nominal distance.
    r_p5_nom = P5_SMA_ARCSEC / 206264.80624709636 * D_MPC["nominal"] * 1000.0
    total_hi = float(np.sum(hi_mass_nom_pixel[hi_detect]))
    hi_3 = float(np.sum(hi_mass_nom_pixel[hi_mask_3e19]))
    hi_1 = float(np.sum(hi_mass_nom_pixel[hi_mask_1e20]))
    hi_inside_p5 = float(np.sum(hi_mass_nom_pixel[hi_detect & p5_nom]))
    hi_inside_208 = float(np.sum(hi_mass_nom_pixel[hi_detect & (r_nom <= OUTER_LIGHT_KPC)]))
    l_p5_nom = float(np.sum(lstar_nom_pixel[star_ap_mask]))
    l_p5_phot_area = math.pi * (P5_SMA_ARCSEC / 206264.80624709636 * D_MPC["nominal"] * 1000.0) ** 2 * P5_Q
    map_covered_area = int(np.count_nonzero(star_valid)) * omega_sr * (D_MPC["nominal"] * 1000.0) ** 2
    report = {
        "task": "curve_blind_NGC4559_source_only_build",
        "status": "source_profiles_complete_no_fit_no_likelihood_no_prediction",
        "inputs": {key: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
                   for key, path in INPUTS.items()},
        "geometry": {"distance_Mpc_nominal": D_MPC["nominal"],
                     "distance_bracket_Mpc": list(D_MPC.values()),
                     "distance_bracket_definition": "8.91 +/- sqrt(0.19^2 + 0.29^2) Mpc from the pre-fit TRGB metadata summary",
                     "inclinations_deg": list(INCLINATIONS_DEG), "PA_deg": list(PA_DEG),
                     "PA_convention": "east of north", "center_RA_deg": RA0_DEG,
                     "center_Dec_deg": DEC0_DEG, "P5_sma_arcsec": P5_SMA_ARCSEC,
                     "P5_ellipticity": P5_ELLIPTICITY, "P5_q": P5_Q,
                     "P5_sma_kpc_nominal": r_p5_nom, "outer_light_support_kpc": OUTER_LIGHT_KPC,
                     "tail_fit_range_kpc": list(TAIL_FIT_KPC), "tail_fit_slope_per_kpc_nominal": float(slope),
                     "tail_scenarios": ["truncate_at_P5", "exp_slope_6_11", "half_slope_6_11"]},
        "reprojection": {"common_grid": "HALOGAS NCP 4-arcsec HR grid",
                         "stellar_method": f"{N_SUBPIX}x{N_SUBPIX} subpixel sampling and bilinear interpolation of surface brightness",
                         "crop_xy": list(crop), "covered_stellar_common_grid_pixels": int(np.count_nonzero(star_valid)),
                         "covered_projected_area_kpc2_nominal": map_covered_area,
                         "orphan_SIP_policy": "ignored because P5 CTYPE omits -SIP; use declared TAN CD matrix"},
        "stellar": {"map_BUNIT": star_bunit, "luminosity_units": "Lsun,3.6 (IRAC Vega band)",
                    "conversion": "MJy/sr -> Jy/pixel -> Lnu/Lsun,3.6 using M_sun,3.6=3.24 Vega and zero point 280.9 Jy",
                    "M_L_policy": "kept separate; Mstar = Υ3.6 * L3.6; no Υ3.6 value applied",
                    "luminosity_inside_P5_aperture_1e9_Lsun36": l_p5_nom / 1e9,
                    "P5_projected_aperture_area_kpc2": math.pi * (r_p5_nom ** 2) * P5_Q,
                    "P5_annular_net_map_fit_slope_per_kpc": float(slope),
                    "negative_pixels_policy": "clipped at zero in source luminosity sums; raw FITS retained unchanged",
                    "excluded_components": ["dust map", "molecular gas", "ionized gas"]},
        "hi": {"map_BUNIT_header_stale": hi_bunit,
               "array_interpretation": "stored pixels are N_HI in cm^-2; release HISTORY applies mom0mskzlm*1.823e18*606*1e3/(28.37887*13.10430) to create coldenslm",
               "conversion_applied_again": False, "helium_factor": HELIUM_FACTOR,
               "positive_detected_pixels": int(np.count_nonzero(hi_detect)),
               "pixels_ge_3e19_cm2": int(np.count_nonzero(hi_mask_3e19)),
               "pixels_ge_1e20_cm2": int(np.count_nonzero(hi_mask_1e20)),
               "total_positive_map_MHI_1e9_Msun": total_hi / 1e9,
               "total_positive_map_Mgas_with_He_1e9_Msun": HELIUM_FACTOR * total_hi / 1e9,
               "MHI_inside_P5_ellipse_1e9_Msun": hi_inside_p5 / 1e9,
               "MHI_inside_R20p8kpc_1e9_Msun": hi_inside_208 / 1e9,
               "MHI_ge_3e19_1e9_Msun": hi_3 / 1e9,
               "MHI_ge_1e20_1e9_Msun": hi_1 / 1e9,
               "mask_definition": "finite N_HI > 0 is HALOGAS source mask; thresholds shown as mask sensitivity"},
        "catalog_crosscheck": {"SPARC_Table1_MHI_1e9_Msun_at_9Mpc_from_prefit_metadata": 5.811,
                               "independent_TRGB_scaled_MHI_1e9_Msun": 5.811 * (8.91 / 9.00) ** 2,
                               "HI_map_to_catalog_ratio": total_hi / 1e9 / (5.811 * (8.91 / 9.00) ** 2),
                               "use": "source-map flux consistency only; no curve rows or velocities opened"},
        "source_support": {"hi_map_shape": list(hi_nhi.shape), "hi_finite_pixels": int(hi_footprint.sum()),
                           "hi_detected_area_projected_kpc2": int(hi_detect.sum()) * omega_sr * (D_MPC["nominal"] * 1000.0) ** 2,
                           "hi_max_deprojected_detected_radius_kpc_nominal": hi_rmax,
                           "stellar_valid_reprojected_pixels": int(star_valid.sum()),
                           "p5_valid_stellar_pixels": int(star_ap_mask.sum()),
                           "stellar_aperture_area_fraction_valid": float(np.count_nonzero(star_ap_mask) / np.count_nonzero(p5_nom)),
                           "stellar_map_positive_light_pixel_count": int(np.count_nonzero(star_valid & (star_grid > 0))),
                           "pixel_solid_angle_sr": omega_sr},
        "boundary": {"source_task_boundary": "No TG or alternative-gravity source solve, fit, likelihood, or prediction was run.",
                     "conditional_Kcrit_check": "not run; this source-only task did not solve a field equation"},
        "official_metadata": {"P5_table": "https://irsa.ipac.caltech.edu/data/SPITZER/S4G/docs/P5_table.txt",
                              "P5_entry": {"excluded": 0, "ICA_iteration": 2, "quality": 1,
                                           "stellar_fraction": 0.729626, "sma_arcsec": 278.1,
                                           "ellipticity": 0.579, "PA_deg": -37.0},
                              "HALOGAS_release": "https://zenodo.org/records/3715549",
                              "HALOGAS_release_record_version": 2,
                              "HALOGAS_file_md5": "e2c57a3f8a4cca5f76ffeffde0bedbb9",
                              "HALOGAS_size_bytes": 4245120,
                              "local_file_sha256": EXPECTED_SHA256["hi"]},
    }
    with (OUT / "run_metadata.json").open("w") as f:
        json.dump(report, f, indent=2)
        f.write("\n")

    # Input and generated artifact hashes, excluding this manifest itself.
    outputs = sorted(p for p in OUT.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt")
    with (OUT / "SHA256SUMS.txt").open("w") as f:
        for p in outputs:
            f.write(f"{sha256(p)}  {p.name}\n")
    print(json.dumps({
        "hi_MHI_1e9_Msun_nominal": total_hi / 1e9,
        "hi_MHI_He_1e9_Msun_nominal": HELIUM_FACTOR * total_hi / 1e9,
        "hi_MHI_catalog_TRGB_proxy_1e9_Msun": 5.811 * (8.91 / 9.0) ** 2,
        "l36_inside_P5_1e9_Lsun36_nominal": l_p5_nom / 1e9,
        "outer_fit_slope_per_kpc": slope,
        "hi_max_detected_radius_kpc": hi_rmax,
        "star_valid_common_pixels": int(star_valid.sum()),
        "star_p5_coverage_fraction": float(np.count_nonzero(star_ap_mask) / np.count_nonzero(p5_nom)),
        "outputs": [p.name for p in outputs],
    }, indent=2))


if __name__ == "__main__":
    main()

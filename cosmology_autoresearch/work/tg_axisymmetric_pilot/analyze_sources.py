#!/usr/bin/env python3
"""Build unit-checked radial source profiles for the NGC 3198 TG pilot.

This is a source-construction/closure stage, not a TG fit.  Inputs are the
target-specific THINGS and HALOGAS H I maps plus the S4G P5 old-stellar-light
map.  The script preserves signed measurements for closure, then clips only
the final azimuthal mean at zero when defining a physical source profile.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = OUT / "data"
THINGS = ROOT / "work/tg_resolved_source_audit/NGC_3198_NA_MOM0_THINGS.FITS"
HALOGAS = ROOT / "work/tg_resolved_source_audit/NGC3198-HR_coldens.fits"
STELLAR = DATA / "NGC3198.stellar.fits"
ICA_MASK = DATA / "NGC3198.ICAmask.fits"

# Walter et al. (2008), Table 2; also inspect the map's converted AIPS history.
D_MPC = 13.8
INCLINATION_DEG = 72.0
PA_DEG = 215.0
# The selected FITS product records the restoring CLEAN beam in AIPS history.
# Its Table-2 nominal natural-weighted beam is retained as a sensitivity case.
THINGS_CLEAN_BEAM_MAJOR_DEG = 3.1753e-3
THINGS_CLEAN_BEAM_MINOR_DEG = 2.6007e-3
THINGS_TABLE_BEAM_MAJOR_ARCSEC = 13.01
THINGS_TABLE_BEAM_MINOR_ARCSEC = 11.56
H0_CGS = 1.6735575e-24
M_SUN_CGS = 1.98847e33
KPC_CM = 3.0856775814913673e21
ARCSEC_RAD = math.pi / (180.0 * 3600.0)
NHI_PER_MSUN_PC2 = 1.25e20
HELIUM_FACTOR = 1.36
RING_WIDTH_KPC = 0.5
R_MAX_KPC = 35.0

# S4G/IRAC conventions.  Querejeta et al. use M_sun,3.6(Vega)=3.24 and
# adopt Υ_3.6=0.6 M_sun/L_sun for the old stellar population.
IRAC_CH1_ZEROPOINT_JY = 280.9
M_SUN_36_VEGA = 3.24
UPSILON_36_REFERENCE = 0.6
S4G_PIXEL_ARCSEC = 0.75
S4G_COMMON_FWHM_ARCSEC = 19.0
S4G_P5_SMA_ARCSEC = 215.4
S4G_P5_ELLIPTICITY = 0.621
S4G_P5_PA_DEG = 33.7


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def image_and_header(path: Path) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path, memmap=False) as hdus:
        header = hdus[0].header.copy()
        data = np.asarray(hdus[0].data, dtype=float).squeeze()
    if data.ndim != 2:
        raise ValueError(f"expected a 2D map after squeeze, got {data.shape}: {path}")
    return data, header


def pixel_scales_arcsec(header: fits.Header) -> tuple[float, float]:
    if "CD1_1" in header:
        sx = abs(float(header["CD1_1"])) * 3600.0
        sy = abs(float(header["CD2_2"])) * 3600.0
    else:
        sx = abs(float(header["CDELT1"])) * 3600.0
        sy = abs(float(header["CDELT2"])) * 3600.0
    if not (sx > 0 and sy > 0):
        raise ValueError(f"non-positive pixel scale: {sx}, {sy}")
    return sx, sy


def disk_coordinates(header: fits.Header, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return projected offsets (arcsec) and deprojected radius (kpc)."""
    sx, sy = pixel_scales_arcsec(header)
    ny, nx = shape
    xpix, ypix = np.meshgrid(np.arange(nx) + 1.0, np.arange(ny) + 1.0)
    # These released maps have diagonal tangent-plane CD/CDELT matrices.
    cd1 = float(header["CD1_1"] if "CD1_1" in header else header["CDELT1"])
    cd2 = float(header["CD2_2"] if "CD2_2" in header else header["CDELT2"])
    dx = (xpix - float(header["CRPIX1"])) * cd1 * 3600.0
    dy = (ypix - float(header["CRPIX2"])) * cd2 * 3600.0
    pa = math.radians(PA_DEG)
    major = dx * math.sin(pa) + dy * math.cos(pa)
    minor = dx * math.cos(pa) - dy * math.sin(pa)
    cosi = math.cos(math.radians(INCLINATION_DEG))
    r_arcsec = np.sqrt(major**2 + (minor / cosi) ** 2)
    r_kpc = r_arcsec * ARCSEC_RAD * D_MPC * 1000.0
    return dx, dy, r_kpc


def elliptical_radius_kpc(dx_arcsec: np.ndarray, dy_arcsec: np.ndarray, pa_deg: float,
                          axis_ratio: float) -> np.ndarray:
    """Projected elliptical radius for an image-space photometric aperture."""
    pa = math.radians(pa_deg)
    major = dx_arcsec * math.sin(pa) + dy_arcsec * math.cos(pa)
    minor = dx_arcsec * math.cos(pa) - dy_arcsec * math.sin(pa)
    return np.sqrt(major**2 + (minor / axis_ratio) ** 2) * ARCSEC_RAD * D_MPC * 1000.0


def ring_profile(data: np.ndarray, radius_kpc: np.ndarray, header: fits.Header, *,
                 physical_floor: bool = False) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    edges = np.arange(0.0, R_MAX_KPC + RING_WIDTH_KPC, RING_WIDTH_KPC)
    sx, sy = pixel_scales_arcsec(header)
    pixel_area_kpc2 = (D_MPC * 1000.0 * sx * ARCSEC_RAD) * (D_MPC * 1000.0 * sy * ARCSEC_RAD)
    projected_area_factor = math.cos(math.radians(INCLINATION_DEG))
    for lo, hi in zip(edges[:-1], edges[1:]):
        use = (radius_kpc >= lo) & (radius_kpc < hi) & np.isfinite(data)
        if not np.any(use):
            continue
        vals = data[use]
        mean = float(np.mean(vals))
        expected_projected_area = math.pi * (hi**2 - lo**2) * projected_area_factor
        coverage = float(np.count_nonzero(use) * pixel_area_kpc2 / expected_projected_area)
        rows.append({
            "r_mid_kpc": 0.5 * (lo + hi),
            "mean": mean,
            "median": float(np.median(vals)),
            "scatter": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
            "n_pixels": int(vals.size),
            "coverage_fraction": coverage,
            "physical_mean": max(0.0, mean) if physical_floor else mean,
        })
    return rows


def sum_in_aperture(data: np.ndarray, radius_kpc: np.ndarray, rmax_kpc: float) -> float:
    use = np.isfinite(data) & (radius_kpc < rmax_kpc)
    return float(np.sum(data[use]))


def hi_mass_from_column(column_cm2: np.ndarray, header: fits.Header, radius_kpc: np.ndarray,
                        rmax_kpc: float) -> float:
    sx, sy = pixel_scales_arcsec(header)
    omega = sx * sy * ARCSEC_RAD**2
    d_cm = D_MPC * 1.0e6 * 3.0856775814913673e18
    area_cm2 = d_cm**2 * omega
    use = np.isfinite(column_cm2) & (radius_kpc < rmax_kpc)
    return float(H0_CGS * area_cm2 * np.sum(column_cm2[use]) / M_SUN_CGS)


def stellar_flux_and_mass(intensity_mjy_sr: np.ndarray, valid: np.ndarray, header: fits.Header,
                          radius_kpc: np.ndarray, rmax_kpc: float) -> dict[str, float]:
    sx, sy = pixel_scales_arcsec(header)
    omega = sx * sy * ARCSEC_RAD**2
    d_pc = D_MPC * 1.0e6
    fnu_jy = intensity_mjy_sr * 1.0e6 * omega
    aperture = valid & np.isfinite(intensity_mjy_sr) & (radius_kpc < rmax_kpc)
    signed_flux = float(np.sum(fnu_jy[aperture]))
    positive_flux = float(np.sum(np.maximum(fnu_jy[aperture], 0.0)))
    solar_flux_10pc = IRAC_CH1_ZEROPOINT_JY * 10.0 ** (-0.4 * M_SUN_36_VEGA)
    lum_signed = signed_flux / solar_flux_10pc * (d_pc / 10.0) ** 2
    lum_positive = positive_flux / solar_flux_10pc * (d_pc / 10.0) ** 2
    return {
        "flux_density_signed_Jy": signed_flux,
        "flux_density_positive_only_Jy": positive_flux,
        "luminosity_signed_Lsun36": lum_signed,
        "luminosity_positive_only_Lsun36": lum_positive,
        "mass_signed_Msun_at_Upsilon_0p6": UPSILON_36_REFERENCE * lum_signed,
        "mass_positive_only_Msun_at_Upsilon_0p6": UPSILON_36_REFERENCE * lum_positive,
    }


def main() -> None:
    things_map, things_h = image_and_header(THINGS)
    halo_map, halo_h = image_and_header(HALOGAS)
    star_map, star_h = image_and_header(STELLAR)
    ica_mask, mask_h = image_and_header(ICA_MASK)
    if ica_mask.shape != star_map.shape:
        raise ValueError(f"S4G mask shape {ica_mask.shape} differs from stellar map {star_map.shape}")
    # P5's published row defines its ICA solution area as a finite P3 ellipse.
    # Outside it the map does not guarantee ICA-cleaned light, so do not silently
    # treat the full FITS cutout as resolved stellar photometry.
    dx_star, dy_star, r_star = disk_coordinates(star_h, star_map.shape)
    p5_q = 1.0 - S4G_P5_ELLIPTICITY
    p5_ellipse_r = elliptical_radius_kpc(dx_star, dy_star, S4G_P5_PA_DEG, p5_q)
    p5_sma_kpc = S4G_P5_SMA_ARCSEC * ARCSEC_RAD * D_MPC * 1000.0

    # P5's positive masks are inherited from P4 and identify excluded
    # foreground/background/artifact pixels. Negative labels flag regions used
    # by the recursive second ICA pass; Querejeta et al. recommend interpolation
    # across these small areas in the released old-star map. Keep those filled
    # values in the physical source, and retain mask==0 as a strict sensitivity.
    star_valid = (ica_mask <= 0) & np.isfinite(star_map) & (p5_ellipse_r <= p5_sma_kpc)

    # THINGS M0 is Jy beam^-1 m s^-1.  Convert m/s -> km/s, then use the
    # standard 21-cm beam-temperature relation and optically thin N_HI law.
    things_beam_major = THINGS_CLEAN_BEAM_MAJOR_DEG * 3600.0
    things_beam_minor = THINGS_CLEAN_BEAM_MINOR_DEG * 3600.0
    beam_area = math.pi / (4.0 * math.log(2.0)) * things_beam_major * things_beam_minor
    t_b_per_jy = 606000.0 / (things_beam_major * things_beam_minor)
    things_nhi = things_map * 1.0e-3 * t_b_per_jy * 1.823e18
    # Independent integrated-flux route.  This also exposes a metadata trap:
    # using only the survey table beam for this FITS map underestimates flux.
    sx_things, sy_things = pixel_scales_arcsec(things_h)
    pixel_area_arcsec2 = sx_things * sy_things
    things_flux_jykms = float(np.sum(things_map) * pixel_area_arcsec2 / beam_area / 1000.0)
    things_flux_tablebeam_jykms = float(
        np.sum(things_map) * pixel_area_arcsec2 /
        (math.pi / (4.0 * math.log(2.0)) * THINGS_TABLE_BEAM_MAJOR_ARCSEC * THINGS_TABLE_BEAM_MINOR_ARCSEC) / 1000.0
    )

    # HALOGAS HR coldens is already N_HI (cm^-2); the BUNIT is stale, but
    # the FITS HISTORY documents the conversion from Jy/beam km/s to N_HI.
    halo_nhi = halo_map.copy()

    _, _, r_things = disk_coordinates(things_h, things_nhi.shape)
    _, _, r_halo = disk_coordinates(halo_h, halo_nhi.shape)
    star_sigma_los = np.maximum(star_map, 0.0) * 10.0**6 / (206265.0**2)
    star_mu_vega = np.full_like(star_sigma_los, np.inf)
    good = star_sigma_los > 0
    star_mu_vega[good] = -2.5 * np.log10(star_sigma_los[good] / IRAC_CH1_ZEROPOINT_JY)
    star_lsun_pc2_los = 10.0 ** (-0.4 * (star_mu_vega - M_SUN_36_VEGA - 21.572))
    star_sigma_msun_pc2 = star_lsun_pc2_los * UPSILON_36_REFERENCE * math.cos(math.radians(INCLINATION_DEG))
    # Smooth the source image to approximately the coarser HALOGAS beam before
    # radial binning.  This changes neither aperture-integrated flux nor mass.
    sigma_pix = (S4G_COMMON_FWHM_ARCSEC / 2.354820045) / S4G_PIXEL_ARCSEC
    star_sigma_num = gaussian_filter(np.where(star_valid, star_sigma_msun_pc2, 0.0),
                                     sigma=sigma_pix, mode="constant", cval=0.0)
    star_sigma_den = gaussian_filter(star_valid.astype(float), sigma=sigma_pix, mode="constant", cval=0.0)
    star_sigma_smooth = np.divide(star_sigma_num, star_sigma_den,
                                  out=np.zeros_like(star_sigma_num), where=star_sigma_den > 0.2)

    things_sigma_gas = np.maximum(things_nhi, 0.0) / NHI_PER_MSUN_PC2 * HELIUM_FACTOR * math.cos(math.radians(INCLINATION_DEG))
    halo_sigma_gas = np.maximum(halo_nhi, 0.0) / NHI_PER_MSUN_PC2 * HELIUM_FACTOR * math.cos(math.radians(INCLINATION_DEG))

    profiles = {
        "THINGS_HIplusHe_faceon_Msun_pc2": ring_profile(things_sigma_gas, r_things, things_h, physical_floor=True),
        "HALOGAS_HIplusHe_faceon_Msun_pc2": ring_profile(halo_sigma_gas, r_halo, halo_h, physical_floor=True),
        "S4G_stars_native_faceon_Msun_pc2_Upsilon0p6": ring_profile(
            np.where(star_valid, star_sigma_msun_pc2, np.nan), r_star, star_h, physical_floor=True),
        "S4G_stars_smoothed_faceon_Msun_pc2_Upsilon0p6": ring_profile(
            np.where(star_valid, star_sigma_smooth, np.nan), r_star, star_h, physical_floor=True),
    }

    apertures = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
    masses: dict[str, dict[str, float]] = {}
    for radius in apertures:
        masses[f"R_lt_{radius:g}_kpc"] = {
            "THINGS_MHI_Msun": hi_mass_from_column(things_nhi, things_h, r_things, radius),
            "HALOGAS_MHI_Msun": hi_mass_from_column(halo_nhi, halo_h, r_halo, radius),
            **{f"S4G_{k}": v for k, v in stellar_flux_and_mass(star_map, star_valid, star_h, r_star, radius).items()},
        }

    meta = {
        "status": "source-construction and closure only; no TG solution or parameter fit",
        "data": {
            str(p.relative_to(ROOT)): {"bytes": p.stat().st_size, "sha256": sha256(p)}
            for p in (THINGS, HALOGAS, STELLAR, ICA_MASK)
        },
        "geometry": {"distance_Mpc": D_MPC, "inclination_deg": INCLINATION_DEG,
                     "position_angle_deg": PA_DEG, "R_max_kpc": R_MAX_KPC,
                     "ring_width_kpc": RING_WIDTH_KPC},
        "THINGS_conversion": {
            "BUNIT": things_h.get("BUNIT"),
            "moment0_to_Tb": "M0[Jy/beam m/s]/1000 * 606000/(11.43108*9.36252) K per Jy/beam",
            "column": "N_HI[cm^-2]=1.823e18 * integral(Tb dv)[K km/s]",
            "selected_FITS_AIPS_CLEAN_beam_arcsec": [things_beam_major, things_beam_minor],
            "Walter_2008_Table_2_NA_beam_arcsec": [THINGS_TABLE_BEAM_MAJOR_ARCSEC, THINGS_TABLE_BEAM_MINOR_ARCSEC],
            "beam_area_arcsec2": beam_area,
            "integrated_flux_Jy_km_s_using_FITS_CLEAN_beam": things_flux_jykms,
            "integrated_flux_Jy_km_s_using_Table_2_beam_sensitivity": things_flux_tablebeam_jykms,
            "published_Walter_2008_flux_Jy_km_s": 227.0,
            "metadata_note": "The selected MOM0 FITS HISTORY says AIPS CLEAN BMAJ=3.1753e-3 deg, BMIN=2.6007e-3 deg (11.43x9.36 arcsec), whereas Walter et al. Table 2 lists the nominal NA cube beam as 13.01x11.56 arcsec. The product-history beam makes integrated flux/mass close to the published total; table-beam-only conversion is retained as sensitivity, not silently substituted.",
        },
        "HALOGAS_conversion": {
            "header_BUNIT": halo_h.get("BUNIT"),
            "used_data_as": "N_HI [cm^-2], as recorded in FITS HISTORY conversion",
            "caveat": "BUNIT is stale Jy/beam km/s; HISTORY documents the column-density conversion",
        },
        "HI_mass_method": "m_H * sum(N_HI * D^2 * pixel_solid_angle) / M_sun; no cos(i) in total mass because projected area and LOS column cancel for a thin disk",
        "gas_source": "face-on Sigma=(N_HI/1.25e20)*1.36*cos(i); thin exponential vertical profile h=0.2 kpc is a later explicit pilot assumption; HALOGAS thick component sensitivity is separate",
        "stellar_conversion": {
            "BUNIT": star_h.get("BUNIT"),
            "ICA_mask_BUNIT": mask_h.get("BUNIT"),
            "ICA_mask_rule": "Use mask<=0 for baseline: positive labels are inherited P4 excluded contaminants; negative labels mark the recursive second-ICA treatment and are retained because P5 recommends interpolation across those regions. Strict mask==0 and unmasked-map alternatives are reported as sensitivities.",
            "P5_quality": {"excluded_from_pipeline": 0, "ICA_iteration": 2, "quality_flag": 2,
                           "s1_over_total_flux_fraction": 0.727999,
                           "analysis_semimajor_arcsec": 215.4,
                           "analysis_ellipticity": 0.621,
                           "analysis_position_angle_deg": 33.7,
                           "analysis_semimajor_kpc_at_adopted_distance": p5_sma_kpc},
            "ICA_mask_pixel_counts": {
                "negative_recursive_ICA": int(np.count_nonzero(ica_mask < 0)),
                "zero_unmasked": int(np.count_nonzero(ica_mask == 0)),
                "positive_inherited_P4": int(np.count_nonzero(ica_mask > 0)),
            },
            "flux": "F_nu[pixel]=I[MJy/sr]*1e6*pixel_solid_angle[sr] Jy",
            "luminosity": "L/Lsun=(Fnu/Fnu_sun_at_10pc)*(D/10pc)^2, Fnu_sun=280.9 Jy*10^(-0.4*3.24)",
            "mass_to_light": UPSILON_36_REFERENCE,
            "face_on_correction": "Sigma_faceon=Sigma_los*cos(i)",
            "radial_profile_smoothing_FWHM_arcsec": S4G_COMMON_FWHM_ARCSEC,
            "aperture_rule": "Limit baseline stellar source to the published P5 analysis ellipse; pixels outside that P3-defined ellipse are not claimed as ICA-cleaned photometry. Any stellar tail beyond it must be an explicit model/extrapolation.",
            "caveat": "S4G P5 is a quality-2 cleaned old-stellar-light map; total mass still depends on M/L, the finite analysis aperture and outer-disk extrapolation",
        },
        "aperture_masses_and_fluxes": masses,
        "S4G_mask_sensitivity_R_lt_35_kpc": {
            "mask_le_0_inside_P5_analysis_ellipse_baseline": stellar_flux_and_mass(star_map, star_valid, star_h, r_star, 35.0),
            "mask_eq_0_inside_P5_analysis_ellipse_strict": stellar_flux_and_mass(
                star_map, (ica_mask == 0) & np.isfinite(star_map) & (p5_ellipse_r <= p5_sma_kpc),
                star_h, r_star, 35.0),
            "all_finite_inside_P5_analysis_ellipse_mask_sensitivity_only": stellar_flux_and_mass(
                star_map, np.isfinite(star_map) & (p5_ellipse_r <= p5_sma_kpc), star_h, r_star, 35.0),
        },
        "full_map_HI_mass_Msun": {
            "THINGS_from_NHI_map": hi_mass_from_column(things_nhi, things_h, r_things, float("inf")),
            "THINGS_from_integrated_flux": 2.356e5 * D_MPC**2 * things_flux_jykms,
            "THINGS_using_Walter_Table_2_beam": 2.356e5 * D_MPC**2 * things_flux_tablebeam_jykms,
            "HALOGAS_from_NHI_map": hi_mass_from_column(halo_nhi, halo_h, r_halo, float("inf")),
        },
        "profile_sample_counts": {k: sum(int(row["n_pixels"]) for row in v) for k, v in profiles.items()},
        "profile_complete_to_kpc_coverage_ge_0p8": {
            k: max((row["r_mid_kpc"] for row in v if row["coverage_fraction"] >= 0.8), default=None)
            for k, v in profiles.items()
        },
    }
    (OUT / "source_summary.json").write_text(json.dumps(meta, indent=2) + "\n")

    # Intersect radial centers and export all profiles to a single CSV.
    by_profile = {
        name: {row["r_mid_kpc"]: row for row in rows if row["coverage_fraction"] >= 0.8}
        for name, rows in profiles.items()
    }
    radial_centers = sorted(set.union(*(set(x) for x in by_profile.values())))
    with (OUT / "source_profiles.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["r_mid_kpc", "things_hihe_faceon_Msun_pc2", "halogas_hihe_faceon_Msun_pc2",
                         "s4g_stars_faceon_Msun_pc2_Upsilon0p6",
                         "s4g_stars_smoothed_faceon_Msun_pc2_Upsilon0p6",
                         "things_npix", "halogas_npix", "s4g_npix", "s4g_smoothed_npix",
                         "things_coverage", "halogas_coverage", "s4g_coverage"])
        for radius in radial_centers:
            things_row = by_profile["THINGS_HIplusHe_faceon_Msun_pc2"].get(radius)
            halo_row = by_profile["HALOGAS_HIplusHe_faceon_Msun_pc2"].get(radius)
            native_star = by_profile["S4G_stars_native_faceon_Msun_pc2_Upsilon0p6"].get(radius)
            smooth_star = by_profile["S4G_stars_smoothed_faceon_Msun_pc2_Upsilon0p6"].get(radius)
            writer.writerow([radius,
                             things_row["physical_mean"] if things_row else "",
                             halo_row["physical_mean"] if halo_row else "",
                             native_star["physical_mean"] if native_star else "",
                             smooth_star["physical_mean"] if smooth_star else "",
                             things_row["n_pixels"] if things_row else "",
                             halo_row["n_pixels"] if halo_row else "",
                             native_star["n_pixels"] if native_star else "",
                             smooth_star["n_pixels"] if smooth_star else "",
                             things_row["coverage_fraction"] if things_row else "",
                             halo_row["coverage_fraction"] if halo_row else "",
                             native_star["coverage_fraction"] if native_star else ""])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
    for ax, image, radius, title, cmap in [
        (axes[0], np.where(star_valid, star_map, np.nan), r_star, "S4G P5 old-star light (MJy/sr)", "magma"),
        (axes[1], np.maximum(halo_nhi, 1.0), r_halo, "HALOGAS N_HI (cm$^{-2}$)", "viridis"),
    ]:
        if "S4G" in title:
            vmin, vmax = 0.005, 10.0
            plotted = np.maximum(image, vmin)
            norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            vmin, vmax = 1.0e19, 1.0e21
            plotted = np.maximum(image, vmin)
            norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
        ax.imshow(plotted, origin="lower", cmap=cmap, norm=norm)
        ax.set_title(title)
        ax.set_xlabel("pixel x")
        ax.set_ylabel("pixel y")
    fig.savefig(OUT / "resolved_source_maps.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for name, label, style in [
        ("THINGS_HIplusHe_faceon_Msun_pc2", "THINGS H I + He", "-"),
        ("HALOGAS_HIplusHe_faceon_Msun_pc2", "HALOGAS H I + He", "--"),
        ("S4G_stars_native_faceon_Msun_pc2_Upsilon0p6", "S4G stars native (Υ=0.6)", "-"),
        ("S4G_stars_smoothed_faceon_Msun_pc2_Upsilon0p6", "S4G stars smoothed to gas beam", ":"),
    ]:
        rows = [row for row in profiles[name] if row["coverage_fraction"] >= 0.8]
        x = np.array([row["r_mid_kpc"] for row in rows])
        y = np.array([row["physical_mean"] for row in rows])
        ax.plot(x, np.maximum(y, 1e-3), style, label=label)
    ax.set_yscale("log")
    ax.set_xlim(0, R_MAX_KPC)
    ax.set_ylim(1e-2, 1e3)
    ax.set_xlabel("deprojected R [kpc]")
    ax.set_ylabel("face-on surface density [M$_\\odot$ pc$^{-2}$]")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(OUT / "radial_baryon_profiles.png", dpi=180)
    plt.close(fig)

    print(json.dumps({"status": meta["status"], "data": meta["data"],
                      "aperture_masses_and_fluxes": masses}, indent=2))


if __name__ == "__main__":
    main()

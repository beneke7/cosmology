#!/usr/bin/env python3
"""Reproduce a local Eq. (47) source-contract audit; no fit or TG solve."""
from __future__ import annotations

import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
from scipy.interpolate import Akima1DInterpolator, CubicSpline, PchipInterpolator


PROJECT = Path(__file__).resolve().parents[2]
ZIP_PATH = PROJECT / "context/data/Rotmod_LTG.zip"
PDF_PATH = PROJECT / "context/user_provided/pszota_van_2024.pdf"
MEMBER = "NGC3198_rotmod.dat"
G_KPC_KMS2_PER_MSUN = 4.30091e-6
UPSILON = 0.762  # Pszota & Ván (2024), Table 1, dimensionless M/L in solar units
NGRID = 20_001


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def negative_intervals(x: np.ndarray, y: np.ndarray) -> list[dict[str, float]]:
    neg = np.flatnonzero(y < 0)
    if not len(neg):
        return []
    ranges = []
    start = previous = int(neg[0])
    for index in neg[1:]:
        index = int(index)
        if index > previous + 1:
            ranges.append({"r_start_kpc": float(x[start]), "r_end_kpc": float(x[previous])})
            start = index
        previous = index
    ranges.append({"r_start_kpc": float(x[start]), "r_end_kpc": float(x[previous])})
    return ranges


def main() -> None:
    zip_bytes = ZIP_PATH.read_bytes()
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        raw = zf.read(MEMBER)
    d = np.loadtxt(io.BytesIO(raw), comments="#")
    r, vobs, err, vgas, vdisk = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4]
    dense_r = np.linspace(r[0], r[-1], NGRID)
    result = {
        "status": "reproduced local source diagnostic; no fit and no TG boundary-value solve",
        "inputs": {
            "paper": "user-provided Pszota & Van, Physics of the Dark Universe 46 (2024) 101660, printed Eq. (47), p. 5; Table 1 M/L",
            "paper_sha256": sha256(PDF_PATH.read_bytes()),
            "rotmod_archive": "SPARC Rotmod_LTG.zip",
            "rotmod_archive_sha256": sha256(zip_bytes),
            "member": MEMBER,
            "member_sha256": sha256(raw),
            "member_bytes": len(raw),
            "rotmod_rows": int(len(r)),
        },
        "equation_contract": {
            "literal_paper_eq47": "rho = [Upsilon d(r Vdisk^2)/dr + d(r Vgas^2)/dr] / (4 pi G r^2)",
            "signed_gas_comparison": "replace Vgas^2 by Vgas*abs(Vgas), preserving the SPARC signed gas-force convention also explicit in paper Eq. (53)",
            "G_kpc_km2_s2_per_Msun": G_KPC_KMS2_PER_MSUN,
            "upsilon_Msun_per_Lsun": UPSILON,
            "radial_data_range_kpc": [float(r[0]), float(r[-1])],
            "uniform_grid_points_for_negative_fraction": NGRID,
            "density_units": "Msun/kpc^3",
            "interpolated_quantity": "velocity components Vdisk(r), Vgas(r), followed by analytic derivative of r*V(r)^2; signed alternative differentiates r*Vgas*abs(Vgas)",
        },
        "endpoint_rows": {
            "first": {"r_kpc": float(r[0]), "Vobs_km_s": float(vobs[0]), "errV_km_s": float(err[0]), "Vgas_km_s": float(vgas[0]), "Vdisk_km_s": float(vdisk[0])},
            "last": {"r_kpc": float(r[-1]), "Vobs_km_s": float(vobs[-1]), "errV_km_s": float(err[-1]), "Vgas_km_s": float(vgas[-1]), "Vdisk_km_s": float(vdisk[-1])},
        },
        "negative_gas_rows": [{"r_kpc": float(r[i]), "Vgas_km_s": float(vgas[i])} for i in np.flatnonzero(vgas < 0)],
        "interpolation_results": {},
        "interpretation": {
            "nonnegative": False,
            "robust_negative_band_kpc": [14.507, 15.993],
            "physical_source_status": "The literal pseudo-spherical inversion yields negative total density robustly near 14.5-16.0 kpc; do not clip. Eq. (47)'s written Vgas^2 also loses SPARC's sign for negative gas velocities.",
            "isolated_exterior_status": "A conditional spherical numerical problem can be posed by explicit center extrapolation, source taper beyond the measured interval, and vacuum Robin matching. These choices are not determined by the local files. Thus a physical isolated-source prediction is not data-closed; a declared surrogate stress test is possible but retains source/model limitation.",
        },
    }

    interpolators = {
        "PCHIP": PchipInterpolator,
        "natural_cubic_spline": lambda x, y: CubicSpline(x, y, bc_type="natural"),
        "Akima": Akima1DInterpolator,
    }
    for name, build in interpolators.items():
        disk_i = build(r, vdisk)
        gas_i = build(r, vgas)
        disk_v = disk_i(dense_r)
        disk_dv = disk_i.derivative()(dense_r)
        gas_v = gas_i(dense_r)
        gas_dv = gas_i.derivative()(dense_r)
        disk_term = UPSILON * (disk_v**2 + 2 * dense_r * disk_v * disk_dv)
        rows = {}
        for convention in ("literal_vgas_squared", "signed_vgas_abs_vgas"):
            if convention == "literal_vgas_squared":
                gas_term = gas_v**2 + 2 * dense_r * gas_v * gas_dv
                mass_proxy = UPSILON * dense_r * disk_v**2 + dense_r * gas_v**2
            else:
                gas_term = gas_v * np.abs(gas_v) + 2 * dense_r * np.abs(gas_v) * gas_dv
                mass_proxy = UPSILON * dense_r * disk_v**2 + dense_r * gas_v * np.abs(gas_v)
            rho = (disk_term + gas_term) / (4 * np.pi * G_KPC_KMS2_PER_MSUN * dense_r**2)
            imin = int(np.argmin(rho))
            # Independent derivative route: finite-difference the enclosed-mass proxy
            # rather than reusing the analytic product-rule derivative above.
            rho_fd = np.gradient(mass_proxy, dense_r, edge_order=2) / (
                4 * np.pi * G_KPC_KMS2_PER_MSUN * dense_r**2
            )
            fd_interior = rho_fd[1:-1]
            x_interior = dense_r[1:-1]
            fd_min_index = int(np.argmin(fd_interior))
            rows[convention] = {
                "rho_min_Msun_kpc3": float(rho[imin]),
                "r_at_min_kpc": float(dense_r[imin]),
                "rho_max_Msun_kpc3": float(np.max(rho)),
                "negative_grid_count": int(np.count_nonzero(rho < 0)),
                "negative_grid_fraction": float(np.mean(rho < 0)),
                "negative_intervals_kpc": negative_intervals(dense_r, rho),
                "independent_finite_difference": {
                    "edge_cells_excluded": 1,
                    "rho_min_Msun_kpc3": float(fd_interior[fd_min_index]),
                    "r_at_min_kpc": float(x_interior[fd_min_index]),
                    "negative_intervals_kpc": negative_intervals(x_interior, fd_interior),
                },
            }
        result["interpolation_results"][name] = rows

    out = Path(__file__).with_name("audit.json")
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(out)


if __name__ == "__main__":
    main()

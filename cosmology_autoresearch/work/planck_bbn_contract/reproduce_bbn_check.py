#!/usr/bin/env python3
"""Reproduce the local CAMB PArthENoPE-grid audit at Planck's reference point."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path

import camb
import numpy as np
from scipy.interpolate import RegularGridInterpolator


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
PAPER = HERE / "1807.06209v4.pdf"
TABLE_NAMES = (
    "PArthENoPE_880.2_standard.dat",
    "PArthENoPE_880.2_marcucci.dat",
    "PRIMAT_Yp_DH_ErrorMC_2021.dat",
)
OMEGA_B = 0.02236
N_EFF_PLANCK = 3.046
T_CMB_PLANCK = 2.7255


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path) -> tuple[list[str], np.ndarray]:
    header = None
    with path.open(encoding="ascii") as stream:
        for line in stream:
            if line.startswith("#") and line[1:].strip().startswith("ombh2"):
                header = line[1:].split()
                break
    if header is None:
        raise ValueError(f"No column-header row in {path}")
    data = np.loadtxt(path)
    if data.shape[1] != len(header):
        raise ValueError(f"Header/data column mismatch for {path}")
    return header, data


def table_result(name: str) -> dict[str, object]:
    path = Path(camb.__file__).resolve().parent / name
    columns, data = read_table(path)
    i_omega, i_delta = columns.index("ombh2"), columns.index("DeltaN")
    x = np.unique(data[:, i_omega])
    y = np.unique(data[:, i_delta])
    outputs = ("Yp", "Yp^BBN", "D/H")
    splined: dict[str, float] = {}
    linear: dict[str, float] = {}
    for output in outputs:
        i = columns.index(output)
        grid = np.array(
            [
                [data[(data[:, i_omega] == xx) & (data[:, i_delta] == yy), i][0] for yy in y]
                for xx in x
            ]
        )
        splined[output] = float(camb.bbn.get_predictor(str(path)).get(output, OMEGA_B, 0.0))
        linear[output] = float(
            RegularGridInterpolator((x, y), grid, method="linear")([[OMEGA_B, 0.0]])[0]
        )
    pars = camb.set_params(
        H0=67.27,
        ombh2=OMEGA_B,
        omch2=0.1202,
        mnu=0.06,
        nnu=N_EFF_PLANCK,
        standard_neutrino_neff=N_EFF_PLANCK,
        TCMB=T_CMB_PLANCK,
        bbn_predictor=name,
    )
    pars_default_reference = camb.set_params(
        H0=67.27,
        ombh2=OMEGA_B,
        omch2=0.1202,
        mnu=0.06,
        nnu=N_EFF_PLANCK,
        standard_neutrino_neff=3.044,
        TCMB=T_CMB_PLANCK,
        bbn_predictor=name,
    )
    return {
        "path": str(path),
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
        "columns": columns,
        "grid": {
            "omega_b_h2_min": float(x[0]),
            "omega_b_h2_max": float(x[-1]),
            "omega_b_h2_points": int(x.size),
            "delta_n_eff_min": float(y[0]),
            "delta_n_eff_max": float(y[-1]),
            "delta_n_eff_points": int(y.size),
        },
        "at_omega_b_h2_0_02236_delta_n_eff_0": {
            "camb_bicubic_interpolation": splined,
            "independent_scipy_bilinear_interpolation": linear,
            "baryon_number_fraction_to_mass_fraction_via_CAMB": float(
                camb.bbn.ypBBN_to_yhe(splined["Yp^BBN"])
            ),
        },
        "direct_camb_set_params": {
            "standard_neutrino_neff_3_046": {
                "YHe_mass_fraction": float(pars.YHe),
                "omega_nu_h2": float(pars.omnuh2),
            },
            "standard_neutrino_neff_3_044_CAMB_default": {
                "YHe_mass_fraction": float(pars_default_reference.YHe),
                "omega_nu_h2": float(pars_default_reference.omnuh2),
            },
        },
    }


def main() -> None:
    package = Path(camb.__file__).resolve().parent
    source_files = [
        package / "bbn.py",
        package / "model.py",
        package / "constants.py",
        PROJECT / ".venv/lib/python3.12/site-packages/cobaya/theories/camb/camb.py",
        PAPER,
    ]
    artifacts = {
        str(path): {"sha256": sha256(path), "size_bytes": path.stat().st_size}
        for path in source_files
        if path.is_file()
    }
    for name in TABLE_NAMES:
        path = package / name
        if path.is_file():
            artifacts[str(path)] = {"sha256": sha256(path), "size_bytes": path.stat().st_size}
    result = {
        "status": "reproduced_local_interpolation_and_CAMB_parameter_mapping_only",
        "camb_version": camb.__version__,
        "cobaya_version": importlib.metadata.version("cobaya"),
        "inputs": {
            "omega_b_h2": OMEGA_B,
            "N_eff": N_EFF_PLANCK,
            "standard_neutrino_neff_for_planck_matching": N_EFF_PLANCK,
            "TCMB_K": T_CMB_PLANCK,
            "delta_neff_table_coordinate": 0.0,
            "H0_initialization_km_s_Mpc": 67.27,
            "omega_c_h2_initialization": 0.1202,
            "mnu_eV_initialization": 0.06,
        },
        "tables": {name: table_result(name) for name in TABLE_NAMES},
        "artifacts": artifacts,
        "limitations": [
            "No PArthENoPE executable/source was run; calculation interpolates the shipped tables.",
            "The CAMB table provenance states PArthENoPE 2017 / May 2017 and tau_n=880.2 s, but does not explicitly identify the executable as PArthENoPE 1.1; exact code-version equivalence is not established.",
            "This is not a Planck likelihood evaluation, parameter fit, or posterior reproduction.",
        ],
    }
    marcucci = result["tables"][TABLE_NAMES[1]]["at_omega_b_h2_0_02236_delta_n_eff_0"]
    if abs(marcucci["camb_bicubic_interpolation"]["Yp^BBN"] - 0.24672) > 1e-5:
        raise AssertionError("Marcucci interpolation does not reproduce Planck's rounded Yp^BBN")
    if abs(marcucci["baryon_number_fraction_to_mass_fraction_via_CAMB"] - 0.2454) > 1e-5:
        raise AssertionError("Marcucci conversion does not reproduce Planck's rounded Y_P")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

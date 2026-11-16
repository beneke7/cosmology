#!/usr/bin/env python3
"""Build the pinned Planck PR3 no-reconstruction-lensing Cobaya input."""

from __future__ import annotations

import hashlib
from importlib.metadata import version as package_version
import json
from pathlib import Path

import cobaya
import camb
import yaml
from cobaya.cosmo_input.create_input import create_input


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "work" / "planck_data" / "plc_3.0"

EXPECTED_COBAYA = "3.6.2"
EXPECTED_CAMB = "1.6.6"
EXPECTED_CLIPY = "0.15"
PLANCK_BBN_TABLE = "PArthENoPE_880.2_marcucci.dat"
PLANCK_NEFF = 3.046
PLANCK_TAU_N_SECONDS = 880.2
PLANCK_YHE_AT_TABLE2_CENTRE = 0.2454
PLANCK_1E5_DH_AT_TABLE2_CENTRE = 2.455
EXPECTED_COMPONENTS = {
    "planck_2018_highl_plik.TTTEEE": (
        DATA / "hi_l" / "plik" / "plik_rd12_HM_v22b_TTTEEE.clik"
    ),
    "planck_2018_lowl.TT_clik": (
        DATA / "low_l" / "commander" / "commander_dx12_v3_2_29.clik"
    ),
    "planck_2018_lowl.EE_clik": (
        DATA / "low_l" / "simall" / "simall_100x143_offlike5_EE_Aplanck_B.clik"
    ),
}

PUBLISHED_CENTERS = {
    "ombh2": 0.02236,
    "omch2": 0.1202,
    "theta_MC_100": 1.04090,
    "tau": 0.0544,
    "logA": 3.045,
    "ns": 0.9649,
}

PUBLISHED_ERRORS = {
    "ombh2": 0.00015,
    "omch2": 0.0014,
    "theta_MC_100": 0.00031,
    "tau": 0.0075,
    "logA": 0.016,
    "ns": 0.0044,
}

PLANCK_PRIORS = {
    "ombh2": {"min": 0.005, "max": 0.1},
    "omch2": {"min": 0.001, "max": 0.99},
    "theta_MC_100": {"min": 0.5, "max": 10.0},
    "tau": {"min": 0.01, "max": 0.8},
    "logA": {"min": 1.61, "max": 3.91},
    "ns": {"min": 0.9, "max": 1.1},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if cobaya.__version__ != EXPECTED_COBAYA:
        raise RuntimeError(
            f"expected Cobaya {EXPECTED_COBAYA}, got {cobaya.__version__}"
        )
    if camb.__version__ != EXPECTED_CAMB:
        raise RuntimeError(f"expected CAMB {EXPECTED_CAMB}, got {camb.__version__}")
    clipy_version = package_version("clipy-like")
    if clipy_version != EXPECTED_CLIPY:
        raise RuntimeError(f"expected clipy-like {EXPECTED_CLIPY}, got {clipy_version}")
    for key, product in EXPECTED_COMPONENTS.items():
        if not product.is_dir():
            raise FileNotFoundError(f"missing selected Planck product for {key}: {product}")

    # Use Cobaya's version-pinned Planck 2018 CAMB preset for the standard cosmology,
    # nuisance defaults and shared likelihood priors. Replace only its likelihood
    # combination with the exact no-lensing-reconstruction PR3 products.
    info = create_input(preset="planck_2018_camb")
    info["likelihood"] = {
        key: {"clik_file": str(product)}
        for key, product in EXPECTED_COMPONENTS.items()
    }
    info["theory"]["camb"]["version"] = "1.6.6"
    info["theory"]["camb"].setdefault("extra_args", {}).update(
        {
            "lens_potential_accuracy": 1,
            "num_massive_neutrinos": 1,
            # Planck VI's baseline BBN branch (its case b) used PArthENoPE 1.1 with
            # the Marcucci et al. (2016) d(p,gamma)3He rate and tau_n=880.2 s. Pin
            # that CAMB table, not CAMB 1.6.6's PRIMAT-2021 default or PArthENoPE's
            # alternate experimental-rate table. CAMB defines Delta N_eff relative
            # to standard_neutrino_neff, so both must match Planck's 3.046 to make
            # the BBN-table interpolation use DeltaN=0.
            "bbn_predictor": PLANCK_BBN_TABLE,
            "nnu": PLANCK_NEFF,
            "standard_neutrino_neff": PLANCK_NEFF,
            "theta_H0_range": [20, 100],
        }
    )

    bbn_predictor = camb.bbn.get_predictor(PLANCK_BBN_TABLE)
    yhe_at_table2_centre = bbn_predictor.Y_He(
        PUBLISHED_CENTERS["ombh2"], delta_neff=0.0
    )
    onee5_dh_at_table2_centre = 1e5 * bbn_predictor.DH(
        PUBLISHED_CENTERS["ombh2"], delta_neff=0.0
    )
    if abs(yhe_at_table2_centre - PLANCK_YHE_AT_TABLE2_CENTRE) > 5e-5:
        raise AssertionError(
            "Pinned Planck BBN table does not reproduce the published reference "
            f"Y_P≈0.2454: got {yhe_at_table2_centre:.9f}"
        )
    if abs(onee5_dh_at_table2_centre - PLANCK_1E5_DH_AT_TABLE2_CENTRE) > 0.005:
        raise AssertionError(
            "Pinned Planck BBN table does not reproduce the case-(b) deuterium "
            f"reference y_DP≈2.455: got {onee5_dh_at_table2_centre:.9f}"
        )

    for name, bounds in PLANCK_PRIORS.items():
        info["params"][name]["prior"] = bounds
        info["params"][name]["ref"] = {
            "dist": "norm",
            "loc": PUBLISHED_CENTERS[name],
            "scale": PUBLISHED_ERRORS[name],
        }
    for name, proposal in {
        "ombh2": 0.00010,
        "omch2": 0.00080,
        "theta_MC_100": 0.00020,
        "tau": 0.0040,
        "logA": 0.010,
        "ns": 0.0025,
    }.items():
        info["params"][name]["proposal"] = proposal

    info["sampler"]["mcmc"].update(
        {
            "Rminus1_stop": 0.01,
            "Rminus1_cl_stop": 0.2,
        }
    )
    info["output"] = str(OUT / "chains" / "planck2018_nolens")

    # Guard against accidental inclusion of the separate lensing reconstruction or a
    # compressed/lite substitute, and verify the exact six sampled cosmological priors.
    if set(info["likelihood"]) != set(EXPECTED_COMPONENTS):
        raise AssertionError(f"unexpected likelihood combination: {info['likelihood'].keys()}")
    if "planck_2018_lensing.clik" in info["likelihood"]:
        raise AssertionError("separate Planck lensing likelihood must not be included")
    for name, bounds in PLANCK_PRIORS.items():
        if info["params"][name]["prior"] != bounds:
            raise AssertionError(f"prior mismatch for {name}")

    yaml_path = OUT / "planck2018_nolens_camb1166.yaml"
    yaml_path.write_text(
        "# Generated by make_config.py; standard Planck PR3 TT,TE,EE+lowE only.\n"
        "# Physical CMB lensing remains on in CAMB; no separate lensing reconstruction.\n"
        + yaml.safe_dump(info, sort_keys=False, width=100),
        encoding="utf-8",
    )
    config_hash = sha256(yaml_path)
    smoke_path = OUT / "smoke_result.json"
    smoke_record = None
    if smoke_path.is_file():
        candidate = json.loads(smoke_path.read_text(encoding="utf-8"))
        if (
            candidate.get("status") == "full_likelihood_one_point_smoke_pass_not_fit"
            and candidate.get("model", {}).get("config_sha256") == config_hash
        ):
            smoke_record = candidate
    manifest = {
        "status": (
            "full_likelihood_one_point_smoke_pass_not_fit"
            if smoke_record is not None
            else "configuration_generated_not_yet_evaluated"
        ),
        "cobaya_version": cobaya.__version__,
        "camb_version": camb.__version__,
        "clipy_like_version": clipy_version,
        "source_preset": "Cobaya 3.6.2 planck_2018_camb; edited only for exact TT,TE,EE+lowE data combination and Planck 2018 N_eff/prior contract",
        "paper_reference": "Planck Collaboration 2018 VI, Table 2, TT,TE,EE+lowE column",
        "separate_lensing_reconstruction": False,
        "physical_cmb_lensing": True,
        "components": {
            key: {"path": str(product), "product_present": product.is_dir()}
            for key, product in EXPECTED_COMPONENTS.items()
        },
        "base_parameters": {
            name: {
                "reference_center": PUBLISHED_CENTERS[name],
                "reference_scale": PUBLISHED_ERRORS[name],
                "prior": PLANCK_PRIORS[name],
            }
            for name in PLANCK_PRIORS
        },
        "fixed_model_settings": {
            "Omega_k": 0.0,
            "w": -1.0,
            "mnu_eV": 0.06,
            "massive_neutrino_species": 1,
            "N_eff": 3.046,
            "helium": {
                "prescription": "PArthENoPE 1.1 BBN consistency interpolation",
                "camb_table": PLANCK_BBN_TABLE,
                "neutron_lifetime_s": PLANCK_TAU_N_SECONDS,
                "standard_neutrino_neff": PLANCK_NEFF,
                "delta_neff_at_reference": 0.0,
                "YHe_at_ombh2_0.02236": yhe_at_table2_centre,
                "1e5_DH_at_ombh2_0.02236": onee5_dh_at_table2_centre,
                "table_sha256": sha256(Path(camb.__file__).parent / PLANCK_BBN_TABLE),
                "dpgamma3He_rate": "Marcucci et al. 2016 theoretical rate, Planck VI case b",
            },
        },
        "config_path": str(yaml_path),
        "config_sha256": config_hash,
        "smoke_evaluation": (
            {
                "record_path": str(smoke_path),
                "log_likelihood_components": smoke_record["log_likelihood_components"],
                "log_likelihood_sum": smoke_record["log_likelihood_sum"],
                "seed_for_nuisance_reference_draw": smoke_record[
                    "seed_for_nuisance_reference_draw"
                ],
                "runtime_seconds": smoke_record["runtime_seconds"],
                "interpretation": smoke_record["interpretation"],
            }
            if smoke_record is not None
            else None
        ),
        "nuisance_policy": "Use full Cobaya 3.6.2 component nuisance defaults and SZ combination prior; no nuisance is frozen by this script.",
    }
    manifest_path = OUT / "config_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {yaml_path}")
    print(f"sha256={manifest['config_sha256']}")
    print(f"wrote {manifest_path}")


if __name__ == "__main__":
    main()

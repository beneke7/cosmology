#!/usr/bin/env python3
"""Initialize and evaluate the full Planck PR3 likelihood at one diagnostic point."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time

import camb
import cobaya
import numpy as np
from cobaya.model import get_model
import yaml


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONFIG = HERE / "planck2018_nolens_camb1166.yaml"
RESULT = HERE / "smoke_result.json"
SEED = 20260925
CENTRE = {
    "ombh2": 0.02236,
    "omch2": 0.1202,
    "theta_MC_100": 1.04090,
    "tau": 0.0544,
    "logA": 3.045,
    "ns": 0.9649,
}
PRODUCTS = {
    "planck_2018_highl_plik.TTTEEE": ROOT
    / "work/planck_data/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik",
    "planck_2018_lowl.TT_clik": ROOT
    / "work/planck_data/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik",
    "planck_2018_lowl.EE_clik": ROOT
    / "work/planck_data/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_sha256(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rel = path.relative_to(root).as_posix().encode()
        digest.update(len(rel).to_bytes(8, "big"))
        digest.update(rel)
        digest.update(path.stat().st_size.to_bytes(8, "big"))
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def main() -> None:
    os.environ["CLIPY_NOJAX"] = "1"
    t0 = time.monotonic()
    info = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    model = get_model(info)

    names = list(model.parameterization.sampled_params())
    rng = np.random.default_rng(SEED)
    reference = model.prior.reference(random_state=rng)
    point = dict(zip(names, map(float, reference), strict=True))
    point.update(CENTRE)
    prior_parts = model.logpriors(point, as_dict=True)
    logprior = float(sum(prior_parts.values()))
    if not np.isfinite(logprior):
        raise RuntimeError(f"diagnostic point has non-finite prior: {prior_parts}")

    loglikes = model.loglikes(point, as_dict=True, return_derived=False)
    if set(loglikes) != set(PRODUCTS):
        raise RuntimeError(f"likelihood-component mismatch: {sorted(loglikes)}")
    if not all(np.isfinite(value) for value in loglikes.values()):
        raise RuntimeError(f"non-finite Planck component score(s): {loglikes}")

    result = {
        "status": "full_likelihood_one_point_smoke_pass_not_fit",
        "interpretation": (
            "One exact likelihood evaluation at the published Table 2 marginalized "
            "cosmological centers with one deterministic draw from the configured "
            "nuisance reference/prior measure. This is not the published best fit "
            "or a reproduction of marginalized constraints."
        ),
        "seed_for_nuisance_reference_draw": SEED,
        "model": {
            "cobaya": cobaya.__version__,
            "camb": camb.__version__,
            "clipy_like": __import__("importlib.metadata", fromlist=["version"])
            .version("clipy-like"),
            "config_path": str(CONFIG),
            "config_sha256": file_sha256(CONFIG),
            "sampled_parameter_count": len(names),
            "sampled_parameter_names": names,
        },
        "inputs": {
            "published_cosmological_centres": CENTRE,
            "complete_sampled_point": point,
            "likelihood_product_tree_sha256": {
                name: tree_sha256(path) for name, path in PRODUCTS.items()
            },
        },
        "log_prior_components": prior_parts,
        "log_prior_sum": logprior,
        "log_likelihood_components": loglikes,
        "log_likelihood_sum": float(sum(loglikes.values())),
        "runtime_seconds": time.monotonic() - t0,
        "thread_environment": {
            key: os.environ.get(key)
            for key in (
                "CLIPY_NOJAX",
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
    }
    RESULT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

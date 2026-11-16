#!/usr/bin/env python3
"""One-point, no-sampler Planck PR3 likelihood smoke test."""

from __future__ import annotations

import hashlib
import importlib.metadata as metadata
import json
import os
import platform
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import yaml


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONFIG = ROOT / "work" / "planck_baseline" / "planck2018_nolens_camb1166.yaml"
MANIFEST = ROOT / "work" / "planck_baseline" / "config_manifest.json"
DATA_ROOT = ROOT / "work" / "planck_data" / "plc_3.0"
RESULT = HERE / "result.json"
SEED = 20260925

COMPONENT_PATHS = {
    "planck_2018_highl_plik.TTTEEE": DATA_ROOT / "hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik",
    "planck_2018_lowl.TT_clik": DATA_ROOT / "low_l/commander/commander_dx12_v3_2_29.clik",
    "planck_2018_lowl.EE_clik": DATA_ROOT / "low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik",
}
TABLE2_CENTERS = {
    "ombh2": 0.02236,
    "omch2": 0.1202,
    "theta_MC_100": 1.04090,
    "tau": 0.0544,
    "logA": 3.045,
    "ns": 0.9649,
}
EXPECTED_VERSIONS = {
    "cobaya": "3.6.2",
    "CAMB": "1.6.6",
    "clipy-like": "0.15",
    "numpy": "2.5.3",
    "scipy": "1.18.1",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha256_tree(path: Path) -> tuple[str, int, int]:
    """Hash sorted relative file names and file hashes, independent of mtimes."""
    files = sorted(p for p in path.rglob("*") if p.is_file())
    h = hashlib.sha256()
    total = 0
    for p in files:
        total += p.stat().st_size
        h.update(p.relative_to(path).as_posix().encode())
        h.update(b"\0")
        h.update(sha256_file(p).encode())
        h.update(b"\n")
    return h.hexdigest(), len(files), total


def dependency_inventory() -> dict[str, dict[str, str]]:
    output = {}
    for name, expected in EXPECTED_VERSIONS.items():
        dist = metadata.distribution(name)
        actual = dist.version
        if actual != expected:
            raise RuntimeError(f"{name} version {actual} != pinned {expected}")
        metadata_path = Path(dist._path) / "METADATA"
        record_path = Path(dist._path) / "RECORD"
        output[name] = {
            "version": actual,
            "metadata_sha256": sha256_file(metadata_path),
            "record_sha256": sha256_file(record_path),
        }
    return output


def write_result(payload: dict) -> None:
    RESULT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    started = time.perf_counter()
    payload: dict = {
        "status": "started",
        "experiment": "one-point full Planck PR3 TTTEEE+lowE likelihood smoke",
        "interpretation": (
            "The Table 2 cosmological marginal centers are evaluated at one seeded "
            "nuisance reference/prior draw. This is not a Planck published best-fit, "
            "maximum-likelihood point, posterior sample, or reproduction of parameter intervals."
        ),
        "seed": SEED,
        "table2_cosmological_centers": TABLE2_CENTERS,
        "python": sys.version,
        "platform": platform.platform(),
        "thread_environment": {
            k: os.environ.get(k)
            for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
        },
        "cpu_affinity_count": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
    }
    try:
        versions = dependency_inventory()
        payload["dependencies"] = versions
        config_sha = sha256_file(CONFIG)
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        if manifest["config_sha256"] != config_sha:
            raise RuntimeError("config SHA256 disagrees with config_manifest.json")
        raw_config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
        config_likelihoods = raw_config.get("likelihood", {})
        if set(config_likelihoods) != set(COMPONENT_PATHS):
            raise RuntimeError(f"likelihood set mismatch: {sorted(config_likelihoods)}")
        for name, expected_path in COMPONENT_PATHS.items():
            configured = Path(config_likelihoods[name]["clik_file"]).resolve()
            if configured != expected_path.resolve() or not configured.is_dir():
                raise RuntimeError(f"likelihood path mismatch or missing product: {name}: {configured}")
        products = {}
        for name, path in COMPONENT_PATHS.items():
            tree_hash, nfiles, nbytes = sha256_tree(path)
            products[name] = {
                "path": str(path),
                "tree_sha256": tree_hash,
                "file_count": nfiles,
                "bytes": nbytes,
            }
        payload["config"] = {"path": str(CONFIG), "sha256": config_sha}
        payload["products"] = products
        payload["likelihood_components_expected"] = list(COMPONENT_PATHS)
        payload["bbn_contract"] = manifest["fixed_model_settings"]["helium"]

        from cobaya.model import get_model

        init_started = time.perf_counter()
        model = get_model(CONFIG, stop_at_error=True)
        init_seconds = time.perf_counter() - init_started
        likelihood_names = list(model.likelihood)
        if set(likelihood_names) != set(COMPONENT_PATHS):
            raise RuntimeError(f"initialized components differ: {likelihood_names}")

        sampled_names = list(model.parameterization.sampled_params())
        fixed_parameters = model.parameterization.constant_params()
        fixed_parameters = {
            k: (v.item() if isinstance(v, np.generic) else v)
            for k, v in fixed_parameters.items()
        }
        nuisance_ref = model.prior.reference(random_state=np.random.default_rng(SEED))
        point = dict(zip(sampled_names, map(float, nuisance_ref), strict=True))
        missing_cosmo = set(TABLE2_CENTERS).difference(point)
        if missing_cosmo:
            raise RuntimeError(f"Table 2 parameter names absent among sampled parameters: {sorted(missing_cosmo)}")
        point.update(TABLE2_CENTERS)
        # Cobaya's ``Prior`` consumes the ordered numeric vector; the parameterization
        # helper intentionally returns a dict when given one.
        point_vector = np.asarray([point[name] for name in model.prior.params], dtype=float)
        logprior_parts = np.asarray(model.prior.logps(point_vector), dtype=float)
        if not np.all(np.isfinite(logprior_parts)):
            raise RuntimeError(f"seeded central-cosmology point has non-finite prior parts: {logprior_parts}")

        eval_started = time.perf_counter()
        loglikes = model.loglikes(point, as_dict=True, return_derived=False)
        eval_seconds = time.perf_counter() - eval_started
        if not isinstance(loglikes, dict) or set(loglikes) != set(COMPONENT_PATHS):
            raise RuntimeError(f"unexpected returned likelihood map: {loglikes}")
        total_loglike = float(sum(loglikes.values()))
        if not np.isfinite(total_loglike):
            raise RuntimeError(f"non-finite total log-likelihood: {total_loglike}")

        nuisance_values = {k: v for k, v in point.items() if k not in TABLE2_CENTERS}
        payload.update(
            {
                "status": "evaluated",
                "initialized_components": likelihood_names,
                "fixed_parameter_values": fixed_parameters,
                "sampled_parameter_count": len(sampled_names),
                "sampled_parameter_names": sampled_names,
                "reference_nuisance_seed": SEED,
                "nuisance_point_and_cosmology": point,
                "nuisance_parameters_only": nuisance_values,
                "logprior_parts": logprior_parts.tolist(),
                "summed_logprior": float(np.sum(logprior_parts)),
                "individual_loglikes": {k: float(v) for k, v in loglikes.items()},
                "summed_loglike": total_loglike,
                "initialization_seconds": init_seconds,
                "evaluation_seconds": eval_seconds,
            }
        )
    except Exception as exc:  # retain failure evidence for the parent
        payload["status"] = "error"
        payload["error_type"] = type(exc).__name__
        payload["error"] = str(exc)
        payload["traceback"] = traceback.format_exc()
    finally:
        payload["total_wall_seconds"] = time.perf_counter() - started
        payload["source_sha256"] = sha256_file(Path(__file__).resolve())
        write_result(payload)
    print(json.dumps({k: payload.get(k) for k in (
        "status", "error_type", "error", "initialized_components", "individual_loglikes",
        "summed_loglike", "summed_logprior", "initialization_seconds", "evaluation_seconds",
        "total_wall_seconds", "source_sha256"
    )}, indent=2, sort_keys=True))
    return 0 if payload["status"] == "evaluated" else 1


if __name__ == "__main__":
    raise SystemExit(main())

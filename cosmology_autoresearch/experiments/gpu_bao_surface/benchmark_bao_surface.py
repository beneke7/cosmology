#!/usr/bin/env python3
"""Benchmark a float64 DESI DR2 CPL profile-chi2 surface on CPU and CUDA.

The three Sobol coordinates are (Omega_m, w0, wa).  At each point alpha =
c/(H0 rd) is profiled analytically over the declared fit bounds, using all
13 rows and the full released covariance.  This is a background-only,
radiation-free screening calculation, not a posterior or evidence estimate.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from background_bao import BAOData, FIT_BOUNDS, GL_NODES, GL_WEIGHTS, predict_bao  # noqa: E402

FLOAT = np.float64
PARAMETER_ORDER = ["Omega_m", "w0", "wa"]
DOMAINS = {"Omega_m": [0.05, 0.6], "w0": [-2.0, -0.3], "wa": [-3.0, 3.0]}
ALPHA_BOUNDS = list(FIT_BOUNDS["cpl"][0])
SEED = 20260924
QUAD_ORDER = len(GL_NODES)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def system_info() -> dict:
    try:
        nvidia = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,driver_version,memory.total,memory.used,memory.free",
             "--format=csv,noheader,nounits"],
            check=True, capture_output=True, text=True, timeout=10,
        ).stdout.strip()
    except Exception as exc:  # diagnostic only; benchmark can still run without the query.
        nvidia = f"unavailable: {type(exc).__name__}: {exc}"
    cpu_model = "unknown"
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu_model = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return {
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "platform": platform.platform(),
        "processor": cpu_model,
        "affinity_cpus": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "thread_environment": {k: os.environ.get(k) for k in
                               ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
        "gpu_nvidia_smi_before_cuda_context": nvidia,
    }


def make_context(mean_path: Path, cov_path: Path) -> dict:
    data = BAOData.from_files(mean_path, cov_path)
    chol = np.linalg.cholesky(data.covariance)
    whiten = solve_triangular(chol, np.eye(len(data.z), dtype=FLOAT), lower=True,
                              check_finite=False)
    y_white = whiten @ data.value
    z = np.asarray(data.z, dtype=FLOAT)
    zq = z[:, None] * (GL_NODES[None, :] + 1.0) / 2.0
    one_plus_zq = 1.0 + zq
    return {
        "data": data,
        "z": z,
        "observable": np.asarray(data.observable),
        "obs_dm": np.asarray(data.observable == "DM_over_rs"),
        "obs_dh": np.asarray(data.observable == "DH_over_rs"),
        "obs_dv": np.asarray(data.observable == "DV_over_rs"),
        "covariance": np.asarray(data.covariance, dtype=FLOAT),
        "chol": chol,
        "whiten": whiten,
        "y_white": y_white,
        "y_norm2": float(y_white @ y_white),
        "zq": zq,
        "one_plus_zq": one_plus_zq,
        "log_one_plus_zq": np.log1p(zq),
        "z_over_one_plus_zq": zq / one_plus_zq,
        "quadrature_scale": z / 2.0,
        "weights": np.asarray(GL_WEIGHTS, dtype=FLOAT),
        "one_plus_z": 1.0 + z,
        "log_one_plus_z": np.log1p(z),
        "z_over_one_plus_z": z / (1.0 + z),
    }


def sobol_parameters(max_power: int, seed: int) -> np.ndarray:
    unit = qmc.Sobol(d=3, scramble=True, seed=seed).random_base2(max_power)
    lower = np.asarray([DOMAINS[name][0] for name in PARAMETER_ORDER], dtype=FLOAT)
    upper = np.asarray([DOMAINS[name][1] for name in PARAMETER_ORDER], dtype=FLOAT)
    return np.asarray(qmc.scale(unit, lower, upper), dtype=FLOAT)


def unit_predictions(xp, parameters, ctx):
    """Return dimensionless predictions at alpha=1 for [B,3] parameters."""
    om = parameters[:, 0, None, None]
    w0 = parameters[:, 1, None, None]
    wa = parameters[:, 2, None, None]

    log_de = 3.0 * (1.0 + w0 + wa) * ctx["log_one_plus_zq"][None, :, :]
    log_de = log_de - 3.0 * wa * ctx["z_over_one_plus_zq"][None, :, :]
    de = (1.0 - om) * xp.exp(log_de)
    e2_q = om * (ctx["one_plus_zq"][None, :, :] ** 3) + de
    inv_e_q = 1.0 / xp.sqrt(e2_q)
    integrals = ctx["quadrature_scale"][None, :] * xp.sum(
        inv_e_q * ctx["weights"][None, None, :], axis=2, dtype=xp.float64
    )

    om_at = parameters[:, 0, None]
    w0_at = parameters[:, 1, None]
    wa_at = parameters[:, 2, None]
    log_de_at = 3.0 * (1.0 + w0_at + wa_at) * ctx["log_one_plus_z"][None, :]
    log_de_at = log_de_at - 3.0 * wa_at * ctx["z_over_one_plus_z"][None, :]
    e2_at = om_at * (ctx["one_plus_z"][None, :] ** 3)
    e2_at = e2_at + (1.0 - om_at) * xp.exp(log_de_at)
    dh = 1.0 / xp.sqrt(e2_at)
    dv = xp.cbrt(ctx["z"][None, :] * integrals**2 * dh)

    return xp.where(ctx["obs_dm"][None, :], integrals,
                    xp.where(ctx["obs_dh"][None, :], dh, dv))


def profile_chunk(xp, parameters, ctx):
    shape = unit_predictions(xp, parameters, ctx)
    shape_white = shape @ ctx["whiten"].T
    numerator = shape_white @ ctx["y_white"]
    denominator = xp.sum(shape_white * shape_white, axis=1, dtype=xp.float64)
    alpha = xp.clip(numerator / denominator, ALPHA_BOUNDS[0], ALPHA_BOUNDS[1])
    residual_white = ctx["y_white"][None, :] - alpha[:, None] * shape_white
    chi2 = xp.sum(residual_white * residual_white, axis=1, dtype=xp.float64)
    return alpha, chi2, shape


def cpu_eval(parameters: np.ndarray, ctx: dict, chunk_size: int) -> tuple[np.ndarray, np.ndarray]:
    alphas = np.empty(len(parameters), dtype=FLOAT)
    chi2 = np.empty(len(parameters), dtype=FLOAT)
    for start in range(0, len(parameters), chunk_size):
        end = min(start + chunk_size, len(parameters))
        alpha, score, _ = profile_chunk(np, parameters[start:end], ctx)
        alphas[start:end] = alpha
        chi2[start:end] = score
    return alphas, chi2


def gpu_eval(parameters: np.ndarray, ctx: dict, chunk_size: int, pool_limit_gib: float) -> dict:
    # Import, CUDA context/device initialization, host-to-device staging, calculation,
    # output transfer, and an explicit final synchronization are all inside wall time.
    start_time = time.perf_counter()
    import cupy as cp

    cp.cuda.Device(0).use()
    pool = cp.get_default_memory_pool()
    pool.set_limit(size=int(pool_limit_gib * 1024**3))
    free_after_context, total_bytes = cp.cuda.runtime.memGetInfo()
    device_context = {key: cp.asarray(value) for key, value in ctx.items()
                      if key in ("whiten", "y_white", "zq", "one_plus_zq", "log_one_plus_zq",
                                 "z_over_one_plus_zq", "quadrature_scale", "weights",
                                 "one_plus_z", "log_one_plus_z", "z_over_one_plus_z",
                                 "z", "obs_dm", "obs_dh", "obs_dv")}
    device_context["y_norm2"] = cp.float64(ctx["y_norm2"])
    n = len(parameters)
    alphas = np.empty(n, dtype=FLOAT)
    chi2 = np.empty(n, dtype=FLOAT)
    peak_pool_reserved = pool.total_bytes()
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        pars_gpu = cp.asarray(parameters[start:end], dtype=cp.float64)
        alpha_gpu, chi2_gpu, _ = profile_chunk(cp, pars_gpu, device_context)
        alphas[start:end] = cp.asnumpy(alpha_gpu)
        chi2[start:end] = cp.asnumpy(chi2_gpu)
        cp.cuda.Stream.null.synchronize()
        peak_pool_reserved = max(peak_pool_reserved, pool.total_bytes())
        del pars_gpu, alpha_gpu, chi2_gpu
    cp.cuda.Stream.null.synchronize()
    elapsed = time.perf_counter() - start_time
    result = {
        "alpha": alphas,
        "chi2": chi2,
        "wall_seconds": elapsed,
        "cupy_version": cp.__version__,
        "cuda_runtime_version": cp.cuda.runtime.runtimeGetVersion(),
        "cuda_driver_version": cp.cuda.runtime.driverGetVersion(),
        "device_name": cp.cuda.runtime.getDeviceProperties(0)["name"].decode(),
        "free_bytes_after_cuda_context": int(free_after_context),
        "total_vram_bytes": int(total_bytes),
        "peak_cupy_pool_reserved_bytes": int(peak_pool_reserved),
        "pool_limit_bytes": int(pool_limit_gib * 1024**3),
    }
    del device_context
    pool.free_all_blocks()
    return result


def scalar_adaptive_shape(parameters: np.ndarray, ctx: dict) -> np.ndarray:
    om, w0, wa = map(float, parameters)
    out: list[float] = []

    def expansion(z: float) -> float:
        one = 1.0 + z
        log_de = 3.0 * ((1.0 + w0 + wa) * math.log1p(z) - wa * z / one)
        return math.sqrt(om * one**3 + (1.0 - om) * math.exp(log_de))

    for z, name in zip(ctx["data"].z, ctx["data"].observable, strict=True):
        zf = float(z)
        dm_unit = quad(lambda zp: 1.0 / expansion(zp), 0.0, zf,
                       epsabs=1e-13, epsrel=1e-13, limit=100)[0]
        dh_unit = 1.0 / expansion(zf)
        values = {"DM_over_rs": dm_unit, "DH_over_rs": dh_unit,
                  "DV_over_rs": (zf * dm_unit * dm_unit * dh_unit) ** (1.0 / 3.0)}
        out.append(values[str(name)])
    return np.asarray(out, dtype=FLOAT)


def adaptive_subset_check(parameters: np.ndarray, cpu_alpha: np.ndarray,
                          cpu_chi2: np.ndarray, ctx: dict, requested: int) -> dict:
    count = min(requested, len(parameters))
    indices = np.unique(np.linspace(0, len(parameters) - 1, count, dtype=np.int64))
    prediction_max_abs = 0.0
    chi2_max_abs = 0.0
    direct_chi2_max_abs = 0.0
    adaptive_runtime = 0.0
    for index in indices:
        t0 = time.perf_counter()
        shape = scalar_adaptive_shape(parameters[index], ctx)
        adaptive_runtime += time.perf_counter() - t0
        alpha = cpu_alpha[index]
        prediction = alpha * shape
        _, _, gl_shape_batch = profile_chunk(np, parameters[index:index + 1], ctx)
        gl_prediction = cpu_alpha[index] * gl_shape_batch[0]
        prediction_max_abs = max(prediction_max_abs, float(np.max(np.abs(prediction - gl_prediction))))

        residual = ctx["data"].value - prediction
        direct_chi2 = float(residual @ np.linalg.solve(ctx["covariance"], residual))
        chi2_max_abs = max(chi2_max_abs, abs(direct_chi2 - float(cpu_chi2[index])))
        direct_chi2_max_abs = max(direct_chi2_max_abs,
                                  abs(direct_chi2 - float(cpu_chi2[index])))
    return {
        "count": int(len(indices)),
        "independent_method": "scalar Python CPL E(z) plus scipy.integrate.quad per observed redshift; direct residual^T C^-1 residual via np.linalg.solve",
        "quad_epsabs": 1e-13,
        "quad_epsrel": 1e-13,
        "max_abs_prediction_difference_vs_96_point_GL": prediction_max_abs,
        "max_abs_chi2_difference_vs_float64_CPU_whitened": chi2_max_abs,
        "max_abs_chi2_difference_direct_full_covariance": direct_chi2_max_abs,
        "total_adaptive_integral_seconds": adaptive_runtime,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mean", type=Path, default=ROOT / "context/data/desi_dr2_mean.txt")
    parser.add_argument("--cov", type=Path, default=ROOT / "context/data/desi_dr2_cov.txt")
    parser.add_argument("--output", type=Path, default=ROOT / "experiments/gpu_bao_surface/result.json")
    parser.add_argument("--surface", type=Path, default=ROOT / "experiments/gpu_bao_surface/surface.npz")
    parser.add_argument("--powers", default="10,12,14,18", help="comma-separated Sobol 2^power batch sizes")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--adaptive-check-points", type=int, default=128)
    parser.add_argument("--pool-limit-gib", type=float, default=2.0)
    args = parser.parse_args()

    powers = [int(value) for value in args.powers.split(",")]
    if not powers or powers != sorted(set(powers)) or min(powers) < 1:
        parser.error("--powers must be unique increasing positive integers")
    if args.chunk_size < 1 or args.pool_limit_gib <= 0:
        parser.error("chunk size and GPU pool limit must be positive")
    max_power = max(powers)

    setup_start = time.perf_counter()
    info = system_info()
    ctx = make_context(args.mean, args.cov)
    parameters = sobol_parameters(max_power, args.seed)
    common_setup_seconds = time.perf_counter() - setup_start

    # One vectorized implementation check against the repository's single-point API.
    check_index = min(17, len(parameters) - 1)
    _, _, shape = profile_chunk(np, parameters[check_index:check_index + 1], ctx)
    prof_alpha, _, _ = profile_chunk(np, parameters[check_index:check_index + 1], ctx)
    actual_alpha = float(prof_alpha[0])
    scalar_prediction = predict_bao(ctx["z"], ctx["observable"], actual_alpha,
                                    *map(float, parameters[check_index]))
    kernel_prediction = actual_alpha * shape[0]
    seed_api_max_abs = float(np.max(np.abs(scalar_prediction - kernel_prediction)))
    if not np.allclose(scalar_prediction, kernel_prediction, rtol=2e-13, atol=2e-13):
        raise AssertionError(f"batch CPU kernel disagrees with repository prediction API: {seed_api_max_abs}")

    result = {
        "status": "reproduced benchmark; exploratory profile surface",
        "compute_agent_model_effort": "gpt-6-luna/max (assigned and active for this worker)",
        "scope": "DESI DR2 13-row compressed BAO likelihood; flat late-time CPL, radiation omitted; free alpha=c/(H0 rd) analytically profiled; no posterior/evidence/significance claim",
        "parameter_order": PARAMETER_ORDER,
        "sobol_domains": DOMAINS,
        "alpha_profile_bounds": ALPHA_BOUNDS,
        "seed": args.seed,
        "sobol": "SciPy scrambled Sobol; deterministic prefix of one 3D sequence; affine uniform mapping on listed domains",
        "quadrature": {"method": "96-point Gauss-Legendre on each [0,z] interval", "order": QUAD_ORDER,
                       "matches_repository_kernel": True},
        "likelihood": "all 13 ordered means and full 13x13 covariance; Cholesky-whitened float64 residual norm; alpha profiled by whitened least squares, then clipped to original CPL alpha bounds",
        "mean_row_order": [{"z": float(z), "observable": str(o)}
                           for z, o in zip(ctx["data"].z, ctx["data"].observable, strict=True)],
        "data_sha256": {str(args.mean): sha256(args.mean), str(args.cov): sha256(args.cov)},
        "code_sha256": {str(Path(__file__)): sha256(Path(__file__)),
                        str(ROOT / "scripts/background_bao.py"): sha256(ROOT / "scripts/background_bao.py")},
        "git_head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                    capture_output=True, text=True, check=True).stdout.strip(),
        "execution": {
            "script_command": shlex.join([sys.executable, *sys.argv]),
            "bounded_wrapper_seconds": 900,
            "thread_limits": {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                              "MKL_NUM_THREADS": "1"},
        },
        "system": info,
        "common_setup_seconds": common_setup_seconds,
        "shared_covariance_details": {
            "shape": list(ctx["data"].covariance.shape),
            "condition_number_2": float(np.linalg.cond(ctx["data"].covariance)),
            "cholesky_min_diagonal": float(np.min(np.diag(ctx["chol"]))),
        },
        "api_prediction_check": {"method": "repository predict_bao() vs independent batched NumPy kernel",
                                 "max_abs_prediction_difference": seed_api_max_abs},
        "cpu_chunk_size": args.chunk_size,
        "gpu_chunk_size": args.chunk_size,
        "gpu_pool_limit_gib": args.pool_limit_gib,
        "adaptive_check_requested": args.adaptive_check_points,
        "batches": [],
    }

    last_cpu_alpha = last_cpu_chi2 = None
    last_gpu_alpha = last_gpu_chi2 = None
    for power in powers:
        n = 1 << power
        subset = parameters[:n]
        cpu_process_start = time.process_time()
        cpu_start = time.perf_counter()
        cpu_alpha, cpu_chi2 = cpu_eval(subset, ctx, args.chunk_size)
        cpu_seconds = time.perf_counter() - cpu_start
        cpu_process_seconds = time.process_time() - cpu_process_start

        gpu_process_start = time.process_time()
        gpu = gpu_eval(subset, ctx, args.chunk_size, args.pool_limit_gib)
        gpu_process_seconds = time.process_time() - gpu_process_start
        gpu_alpha, gpu_chi2 = gpu["alpha"], gpu["chi2"]

        abs_alpha = np.abs(cpu_alpha - gpu_alpha)
        abs_chi2 = np.abs(cpu_chi2 - gpu_chi2)
        cpu_min_index = int(np.argmin(cpu_chi2))
        batch = {
            "power": power,
            "samples": n,
            "cpu_float64_seconds": cpu_seconds,
            "cpu_process_seconds": cpu_process_seconds,
            "cpu_core_equivalent_utilization": cpu_process_seconds / cpu_seconds,
            "cpu_throughput_samples_per_second": n / cpu_seconds,
            "cpu_wall_including_common_setup_seconds": common_setup_seconds + cpu_seconds,
            "gpu_end_to_end_float64_seconds": gpu["wall_seconds"],
            "gpu_process_cpu_seconds": gpu_process_seconds,
            "gpu_process_core_equivalent_utilization": gpu_process_seconds / gpu["wall_seconds"],
            "gpu_throughput_samples_per_second": n / gpu["wall_seconds"],
            "gpu_wall_including_common_setup_seconds": common_setup_seconds + gpu["wall_seconds"],
            "gpu_speedup_over_cpu_compute_with_gpu_setup_and_transfers": cpu_seconds / gpu["wall_seconds"],
            "gpu_speedup_over_full_pipeline_including_common_setup":
                (common_setup_seconds + cpu_seconds) / (common_setup_seconds + gpu["wall_seconds"]),
            "gpu_includes": ["CuPy import on first measured batch", "CUDA context/device setup",
                             "host-to-device data and parameter transfers", "96-point distance integration",
                             "full-covariance whitening and alpha profile", "device-to-host outputs",
                             "explicit CUDA synchronization"],
            "max_abs_alpha_difference_gpu_vs_cpu": float(np.max(abs_alpha)),
            "max_abs_chi2_difference_gpu_vs_cpu": float(np.max(abs_chi2)),
            "max_rel_chi2_difference_gpu_vs_cpu": float(np.max(abs_chi2 / np.maximum(np.abs(cpu_chi2), 1.0))),
            "max_cpu_chi2": float(np.max(cpu_chi2)),
            "min_sampled_profile_chi2": float(cpu_chi2[cpu_min_index]),
            "min_sampled_profile_parameters": {
                "Omega_m": float(subset[cpu_min_index, 0]),
                "w0": float(subset[cpu_min_index, 1]),
                "wa": float(subset[cpu_min_index, 2]),
                "alpha": float(cpu_alpha[cpu_min_index]),
            },
            "fraction_alpha_at_profile_bound": float(np.mean(
                np.isclose(cpu_alpha, ALPHA_BOUNDS[0], rtol=0.0, atol=1e-12) |
                np.isclose(cpu_alpha, ALPHA_BOUNDS[1], rtol=0.0, atol=1e-12))),
            "gpu": {key: value for key, value in gpu.items()
                    if key not in ("alpha", "chi2", "wall_seconds")},
        }
        if not np.all(np.isfinite(cpu_chi2)) or not np.all(np.isfinite(gpu_chi2)):
            raise AssertionError(f"non-finite likelihood output at 2^{power}")
        if batch["max_abs_chi2_difference_gpu_vs_cpu"] > 2e-8:
            raise AssertionError(f"GPU/CPU float64 chi2 mismatch at 2^{power}: {batch['max_abs_chi2_difference_gpu_vs_cpu']}")
        result["batches"].append(batch)
        last_cpu_alpha, last_cpu_chi2 = cpu_alpha, cpu_chi2
        last_gpu_alpha, last_gpu_chi2 = gpu_alpha, gpu_chi2
        print(json.dumps({key: batch[key] for key in
                          ("power", "samples", "cpu_float64_seconds", "cpu_process_seconds",
                           "gpu_end_to_end_float64_seconds", "gpu_process_cpu_seconds",
                           "gpu_speedup_over_full_pipeline_including_common_setup",
                           "max_abs_chi2_difference_gpu_vs_cpu",
                           "min_sampled_profile_chi2")}, sort_keys=True), flush=True)

    assert last_cpu_alpha is not None and last_gpu_alpha is not None
    result["adaptive_scalar_cpu_check"] = adaptive_subset_check(
        parameters[:len(last_cpu_chi2)], last_cpu_alpha, last_cpu_chi2,
        ctx, args.adaptive_check_points,
    )
    result["adaptive_scalar_cpu_check"]["gpu_subset_max_abs_chi2_difference"] = float(
        np.max(np.abs(last_gpu_chi2 - last_cpu_chi2)))
    result["surface_file"] = str(args.surface)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.surface.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.surface,
                        parameters=parameters[:len(last_cpu_chi2)],
                        alpha_profile=last_cpu_alpha,
                        chi2_cpu_float64=last_cpu_chi2,
                        chi2_gpu_float64=last_gpu_chi2)
    result["surface_sha256"] = sha256(args.surface)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "surface": str(args.surface),
                      "surface_sha256": result["surface_sha256"],
                      "adaptive_scalar_cpu_check": result["adaptive_scalar_cpu_check"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

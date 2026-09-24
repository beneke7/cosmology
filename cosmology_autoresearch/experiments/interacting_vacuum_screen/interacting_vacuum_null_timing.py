#!/usr/bin/env python3
"""40-mock timing tranche for the matched interacting-vacuum BAO null search.

This deliberately evaluates only indices 0--19 from each of two pre-existing
PCG64 streams. It is not the full calibration and is not intended to estimate
the tail probability precisely.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
from scipy.optimize import minimize
from scipy.linalg import solve_triangular

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import interacting_vacuum_profile as iv  # noqa: E402

SEEDS = {20260924: 200, 20260925: 1000}
N_PER_SEED = 20
N_STARTS_INTERACTION = 39
GRID_N = 41
WORKER_DATA = None
WORKER_CHOL = None

JSONL_PATH = HERE / "null_timing_40.jsonl"
RESULT_PATH = HERE / "null_timing_40.json"
CHECKPOINT_PATH = HERE / "null_timing_40_checkpoint.json"
PROVENANCE_PATH = ROOT / "work/data_audit2/provenance_check.json"
SEED_MEAN_PATH = ROOT / "experiments/campaign_seed_bao/result.json"
SEED_RUNS = {
    20260924: ROOT / "experiments/bao_robustness/result.json",
    20260925: ROOT / "experiments/bao_robustness/null_extension_1000_seed20260925.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def vector_sha256(y: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(y, dtype="<f8").tobytes()).hexdigest()


def worker_init(z, observable, covariance):
    global WORKER_DATA, WORKER_CHOL
    WORKER_DATA = iv.bao.BAOData(np.asarray(z, dtype=np.float64), np.zeros(len(z), dtype=np.float64),
                                 np.asarray(observable, dtype=str), np.asarray(covariance, dtype=np.float64))
    WORKER_CHOL = np.linalg.cholesky(WORKER_DATA.covariance)


def grid_and_refinement(data, chol, y_white):
    """Repeat the observed-fit 41x41 scan and one bounded L-BFGS-B refinement."""
    oms = np.linspace(*iv.OM_BOUNDS, GRID_N, dtype=np.float64)
    gs = np.linspace(*iv.GAMMA_BOUNDS, GRID_N, dtype=np.float64)
    scores = np.full((GRID_N, GRID_N), np.nan, dtype=np.float64)
    grid_failures = {}
    for i, om in enumerate(oms):
        for j, g in enumerate(gs):
            try:
                scores[i, j] = iv.profile_alpha(data, chol, y_white, float(om), float(g))[0]
            except (iv.DomainError, ValueError, FloatingPointError, OverflowError) as exc:
                name = str(exc)
                grid_failures[name] = grid_failures.get(name, 0) + 1
    idx = np.unravel_index(np.nanargmin(scores), scores.shape)
    x0 = np.asarray([oms[idx[0]], gs[idx[1]]], dtype=np.float64)
    refine_failures = {}

    def objective(theta):
        try:
            return iv.profile_alpha(data, chol, y_white, float(theta[0]), float(theta[1]))[0]
        except (iv.DomainError, ValueError, FloatingPointError, OverflowError) as exc:
            name = str(exc)
            refine_failures[name] = refine_failures.get(name, 0) + 1
            return iv.INVALID

    opt = minimize(objective, x0, method="L-BFGS-B", bounds=[iv.OM_BOUNDS, iv.GAMMA_BOUNDS],
                   options={"maxiter": 1200, "ftol": 1e-14, "gtol": 2e-8})
    try:
        score, alpha, pred, unit, bg = iv.profile_alpha(
            data, chol, y_white, float(opt.x[0]), float(opt.x[1]))
        physical = bool(np.isfinite(score) and score < iv.INVALID)
        validation_error = None if physical else "nonfinite or penalty objective"
    except (iv.DomainError, ValueError, FloatingPointError, OverflowError) as exc:
        score, alpha, pred, unit, bg = iv.INVALID, None, None, None, None
        physical, validation_error = False, str(exc)
    return {
        "x0": x0.tolist(), "grid_best_x": [float(oms[idx[0]]), float(gs[idx[1]])],
        "grid_min_chi2": float(scores[idx]), "grid_valid_points": int(np.isfinite(scores).sum()),
        "grid_invalid_points": int(scores.size - np.isfinite(scores).sum()),
        "grid_domain_failures": grid_failures,
        "refinement": {"x": [float(v) for v in opt.x], "scipy_fun": float(opt.fun),
                       "endpoint_chi2": float(score), "endpoint_valid_physical": physical,
                       "endpoint_validation_error": validation_error,
                       "success": bool(opt.success), "status": int(opt.status),
                       "message": str(opt.message), "nit": int(opt.nit), "nfev": int(opt.nfev),
                       "domain_failures": refine_failures},
        "candidate": (None if not physical else {"Omega_m": float(opt.x[0]), "Gamma_over_H0": float(opt.x[1]),
                     "alpha": float(alpha), "chi2": float(score), "prediction": pred.tolist(),
                     "unit_prediction": unit.tolist(), "status": int(opt.status),
                     "success": bool(opt.success), "source": "41x41 grid plus L-BFGS-B refinement",
                     "parameter_bound_hits": (
                         (["Omega_m:lower"] if abs(float(opt.x[0])-iv.OM_BOUNDS[0]) < 1e-8 else []) +
                         (["Omega_m:upper"] if abs(float(opt.x[0])-iv.OM_BOUNDS[1]) < 1e-8 else []) +
                         (["Gamma_over_H0:lower"] if abs(float(opt.x[1])-iv.GAMMA_BOUNDS[0]) < 1e-7 else []) +
                         (["Gamma_over_H0:upper"] if abs(float(opt.x[1])-iv.GAMMA_BOUNDS[1]) < 1e-7 else [])),
                     "domain_minima": {k: bg[k] for k in ("minimum_e2", "minimum_matter", "minimum_vacuum")}}),
    }


def fit_mock(task):
    """Worker function. Every endpoint/status is retained in the returned row."""
    seed, index, y_list, mock_hash = task
    cpu_start = time.process_time()
    wall_start = time.monotonic()
    y = np.asarray(y_list, dtype=np.float64)
    data = iv.bao.BAOData(WORKER_DATA.z, y, WORKER_DATA.observable, WORKER_DATA.covariance)
    y_white = solve_triangular(WORKER_CHOL, y, lower=True, check_finite=False)
    row = {"seed": int(seed), "index": int(index), "mock_vector_sha256_little_endian_float64": mock_hash,
           "mock_vector": y.tolist(), "status": "complete"}
    try:
        flat = iv.run_fit(data, WORKER_CHOL, y_white, gamma_fixed=0.0, starts=16)
        interaction_multistart = iv.run_fit(data, WORKER_CHOL, y_white, gamma_fixed=None,
                                            gamma_bounds=iv.GAMMA_BOUNDS, starts=N_STARTS_INTERACTION)
        grid = grid_and_refinement(data, WORKER_CHOL, y_white)
        candidates = [{"Omega_m": interaction_multistart["Omega_m"],
                       "Gamma_over_H0": interaction_multistart["Gamma_over_H0"],
                       "alpha": interaction_multistart["alpha"], "chi2": interaction_multistart["chi2"],
                       "prediction": interaction_multistart["prediction"],
                       "unit_prediction": interaction_multistart["unit_prediction"],
                       "status": interaction_multistart["selected_status"],
                       "success": interaction_multistart["selected_success"],
                       "source": "39-start bounded multistart",
                       "parameter_bound_hits": interaction_multistart["parameter_bound_hits"],
                       "physical_decomposition_witness": interaction_multistart["physical_decomposition_witness"],
                       "domain_minima": interaction_multistart["domain_minima"]}]
        if grid["candidate"] is not None:
            candidates.append(grid["candidate"])
        raw_best = min(candidates, key=lambda item: item["chi2"])
        fallback_used = raw_best["chi2"] > flat["chi2"]
        selected = ({"Omega_m": flat["Omega_m"], "Gamma_over_H0": 0.0, "alpha": flat["alpha"],
                     "chi2": flat["chi2"], "prediction": flat["prediction"],
                     "unit_prediction": flat["unit_prediction"], "source": "exact refitted nested flat-LCDM fallback",
                     "status": flat["selected_status"], "success": flat["selected_success"],
                     "parameter_bound_hits": flat["parameter_bound_hits"]} if fallback_used else raw_best)
        row.update({
            "flat_fit": {k: flat[k] for k in ("Omega_m", "Gamma_over_H0", "alpha", "chi2", "prediction",
                                                "selected_start_index", "selected_success", "selected_status",
                                                "parameter_bound_hits", "optimizer_endpoint_audit",
                                                "optimizer_attempts", "best_unsuccessful_candidate",
                                                "domain_failure_evaluations_by_reason")},
            "interaction_raw_search_optimum": raw_best,
            "interaction_multistart": {"Omega_m": interaction_multistart["Omega_m"],
                "Gamma_over_H0": interaction_multistart["Gamma_over_H0"], "alpha": interaction_multistart["alpha"],
                "chi2": interaction_multistart["chi2"], "selected_start_index": interaction_multistart["selected_start_index"],
                "selected_success": interaction_multistart["selected_success"], "selected_status": interaction_multistart["selected_status"],
                "parameter_bound_hits": interaction_multistart["parameter_bound_hits"],
                "optimizer_endpoint_audit": interaction_multistart["optimizer_endpoint_audit"],
                "optimizer_attempts": interaction_multistart["optimizer_attempts"],
                "best_unsuccessful_candidate": interaction_multistart["best_unsuccessful_candidate"],
                "domain_failure_evaluations_by_reason": interaction_multistart["domain_failure_evaluations_by_reason"]},
            "interaction_grid_refinement": grid,
            "nested_flat_fallback_used": bool(fallback_used),
            "selected_interaction_candidate": selected,
            "raw_delta_chi2_flat_minus_interaction": float(flat["chi2"] - raw_best["chi2"]),
            "selection_delta_chi2_flat_minus_interaction": float(flat["chi2"] - selected["chi2"]),
            "fit_runtime_wall_seconds": time.monotonic() - wall_start,
        })
    except Exception as exc:
        row.update({"status": "fit_exception", "error_type": type(exc).__name__, "error": str(exc)})
    row["worker_cpu_seconds"] = time.process_time() - cpu_start
    return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--checkpoint-every", type=int, default=25)
    args = parser.parse_args()
    if args.workers != 20:
        raise ValueError("the locked timing tranche requires exactly 20 worker processes")
    for output in (JSONL_PATH, RESULT_PATH, CHECKPOINT_PATH):
        if output.exists():
            raise FileExistsError(f"refusing to overwrite existing tranche artifact: {output}")
    affinity = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    if affinity < args.workers:
        raise RuntimeError(f"20 workers requested but only {affinity} CPUs are in process affinity")

    mean_path = ROOT / "context/data/desi_dr2_mean.txt"
    cov_path = ROOT / "context/data/desi_dr2_cov.txt"
    provenance = json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))
    if provenance.get("status") != "pass":
        raise RuntimeError("data_audit2 provenance record is not in pass state")
    provenance_streams = provenance["mock_generator_contract"]["streams"]
    provenance_sources = provenance["source_hash_audit"]
    for seed, path in SEED_RUNS.items():
        audited_source = provenance_sources[str(seed)]
        if not audited_source["data_hashes_match_current"] or not audited_source["run_bao_robustness_hash_matches_current"]:
            raise RuntimeError(f"source/data hashes for prior seed {seed} fail work/data_audit2 verification")
        if sha256(path) != audited_source["result_sha256"]:
            raise RuntimeError(f"prior result hash changed after work/data_audit2 for seed {seed}")
    seed_doc = json.loads(SEED_MEAN_PATH.read_text(encoding="utf-8"))
    flat_seed = next(item for item in seed_doc["models"] if item["model"] == "lcdm")
    null_mean = np.asarray(flat_seed["prediction"], dtype=np.float64)
    base = iv.data_load()
    covariance = np.asarray(base.covariance, dtype=np.float64)
    chol = np.linalg.cholesky(covariance)
    null_mean_recheck = flat_seed["alpha"] * iv.unit_prediction(base, flat_seed["Omega_m"], 0.0)[0]
    null_mean_max_abs_difference = float(np.max(np.abs(null_mean - null_mean_recheck)))
    if null_mean_max_abs_difference > 1e-10:
        raise RuntimeError(f"fitted flat null mean mismatch: {null_mean_max_abs_difference:.4g}")

    tasks = []
    stream_records = {}
    for seed, count in SEEDS.items():
        rng = np.random.default_rng(seed)
        if type(rng.bit_generator).__name__ != "PCG64":
            raise RuntimeError("NumPy default_rng did not resolve to PCG64")
        normals = rng.standard_normal((count, len(base.z)))
        digest = hashlib.sha256(np.ascontiguousarray(normals, dtype="<f8").tobytes()).hexdigest()
        audited = provenance_streams[str(seed)]
        if digest != audited["normal_stream_sha256_little_endian_float64"]:
            raise RuntimeError(f"PCG64 stream digest mismatch for seed {seed}")
        mocks = null_mean[None, :] + normals @ chol.T
        indices = list(range(N_PER_SEED))
        stream_records[str(seed)] = {"full_stream_n": count, "selected_zero_based_indices": indices,
                                     "bit_generator": "PCG64", "normals_shape": list(normals.shape),
                                     "normal_stream_sha256_little_endian_float64": digest,
                                     "normal_digest_matches_data_audit2": True,
                                     "selected_vector_sha256_little_endian_float64":
                                         [vector_sha256(mocks[i]) for i in indices]}
        for i in indices:
            tasks.append((seed, i, mocks[i].tolist(), vector_sha256(mocks[i])))

    fit_source = HERE / "interacting_vacuum_profile.py"
    input_paths = [mean_path, cov_path, SEED_MEAN_PATH, *SEED_RUNS.values(), PROVENANCE_PATH, fit_source,
                   ROOT / "scripts/background_bao.py"]
    input_hashes = {str(path.relative_to(ROOT)): sha256(path) for path in input_paths}
    observed = json.loads((HERE / "result.json").read_text(encoding="utf-8"))
    observed_delta = float(observed["profile_comparison"]["delta_chi2_flat_minus_interacting"])
    run_started = time.monotonic()
    rows = []
    result_idx = 0
    with JSONL_PATH.open("x", encoding="utf-8") as output_stream:
        with ProcessPoolExecutor(max_workers=args.workers, initializer=worker_init,
                                 initargs=(base.z.tolist(), base.observable.tolist(), covariance.tolist())) as pool:
            futures = [pool.submit(fit_mock, task) for task in tasks]
            for future in as_completed(futures):
                row = future.result()
                output_stream.write(json.dumps(row, separators=(",", ":"), allow_nan=False) + "\n")
                output_stream.flush()
                os.fsync(output_stream.fileno())
                rows.append(row)
                result_idx += 1
                if result_idx % args.checkpoint_every == 0 or result_idx == len(tasks):
                    checkpoint = {"status": "running" if result_idx < len(tasks) else "complete",
                                  "completed_rows": result_idx, "requested_rows": len(tasks),
                                  "completed_seed_indices": [[r["seed"], r["index"]] for r in rows],
                                  "jsonl_sha256": sha256(JSONL_PATH),
                                  "last_checkpoint_wall_seconds": time.monotonic() - run_started,
                                  "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --seconds 1800 -- .venv/bin/python experiments/interacting_vacuum_screen/interacting_vacuum_null_timing.py --workers 20 --checkpoint-every 25"}
                    CHECKPOINT_PATH.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    elapsed = time.monotonic() - run_started
    rows.sort(key=lambda r: (r["seed"], r["index"]))
    completed = [row for row in rows if row["status"] == "complete"]
    deltas = np.asarray([row["selection_delta_chi2_flat_minus_interaction"] for row in completed], dtype=np.float64)
    raw_deltas = np.asarray([row["raw_delta_chi2_flat_minus_interaction"] for row in completed], dtype=np.float64)
    mock_failures = [row for row in rows if row["status"] != "complete"]
    flat_raw_failures = sum(row["flat_fit"]["optimizer_endpoint_audit"]["raw_scipy_failure_count"] for row in completed)
    interaction_raw_failures = sum(row["interaction_multistart"]["optimizer_endpoint_audit"]["raw_scipy_failure_count"] for row in completed)
    flat_invalid = sum(row["flat_fit"]["optimizer_endpoint_audit"]["invalid_penalty_endpoint_count"] for row in completed)
    interaction_invalid = sum(row["interaction_multistart"]["optimizer_endpoint_audit"]["invalid_penalty_endpoint_count"] for row in completed)
    result = {
        "status": "complete_40_mock_timing_tranche" if len(completed) == len(tasks) else "partial_40_mock_timing_tranche",
        "scope": "Timing tranche only: the first 20 PCG64 mock indices from each previously audited fitted-flat-LCDM seed stream; not the remaining 1160 mocks and not a stable tail calibration.",
        "null_definition": {"mean": {"source": "campaign_seed_bao/result.json fitted flat LCDM prediction",
                        "Omega_m": flat_seed["Omega_m"], "alpha": flat_seed["alpha"], "prediction_sha256_little_endian_float64": vector_sha256(null_mean)},
                        "draw": "np.random.default_rng(seed).standard_normal((N,13)); y_i=mean+normal_i@chol(full covariance).T",
                        "seed_streams": stream_records,
                        "mean_prediction_max_abs_recheck_difference": null_mean_max_abs_difference,
                        "stream_separation_note": "The two seeds are independently initialized streams; selected exact mock vectors are kept separate by seed/index and were regenerated from complete N-row draws."},
        "contract": {"Omega_m": iv.OM_BOUNDS, "Gamma_over_H0": iv.GAMMA_BOUNDS, "alpha": iv.ALPHA_BOUNDS,
                     "interaction_fit": "39 deterministic L-BFGS-B starts plus the exact observed 41x41 grid and bounded L-BFGS-B refinement; each returned optimizer endpoint is re-evaluated with the physical ODE/positive-component feasibility mask",
                     "flat_fit": "16 deterministic L-BFGS-B starts with the exact nested flat LCDM limit",
                     "accepted_endpoint_policy": "For multistart, select the lowest-scoring physically valid finite endpoint among raw-success endpoints when any exist, otherwise among valid endpoints; grid refinement is a separate valid candidate whether or not raw SciPy success is true, matching observed search procedure. Raw SciPy status is never conflated with physical endpoint validity.",
                     "nested_flat_fallback": "If the best accepted raw interaction-search candidate is worse than the independently refitted flat score, include that exact flat fit at Gamma/H0=0 as the selected nested fallback while retaining the raw interaction result separately.",
                     "statistic": "raw delta=chi2_flat-refit-chi2_interaction-search; selected delta=chi2_flat-refit-chi2_selected-interaction-including-nested-flat-fallback"},
        "checks": {"data_audit2_status": provenance["status"], "data_audit2_normal_stream_digests_match": True,
                   "exact_mock_count": len(tasks), "completed_mock_count": len(completed),
                   "failed_mock_records": mock_failures,
                   "tail_count_at_or_above_observed_delta_chi2": int(np.sum(deltas >= observed_delta)),
                   "tail_fraction_descriptive_only": float(np.mean(deltas >= observed_delta)) if len(deltas) else None,
                   "observed_delta_chi2_reference": observed_delta,
                   "not_a_stable_tail_estimate": True,
                   "flat_optimizer_raw_scipy_failures_total": int(flat_raw_failures),
                   "interaction_multistart_raw_scipy_failures_total": int(interaction_raw_failures),
                   "flat_invalid_penalty_endpoints_total": int(flat_invalid),
                   "interaction_invalid_penalty_endpoints_total": int(interaction_invalid),
                   "grid_refinement_raw_scipy_failures_total": int(sum(not r["interaction_grid_refinement"]["refinement"]["success"] for r in completed)),
                   "grid_invalid_evaluations_total": int(sum(r["interaction_grid_refinement"]["grid_invalid_points"] for r in completed)),
                   "nested_fallback_count": int(sum(r["nested_flat_fallback_used"] for r in completed)),
                   "selection_delta_summary_descriptive": ({"minimum": float(np.min(deltas)), "median": float(np.median(deltas)),
                       "maximum": float(np.max(deltas)), "values": deltas.tolist()} if len(deltas) else None),
                   "raw_delta_summary_descriptive": ({"minimum": float(np.min(raw_deltas)), "median": float(np.median(raw_deltas)),
                       "maximum": float(np.max(raw_deltas)), "values": raw_deltas.tolist()} if len(raw_deltas) else None)},
        "runtime": {"wall_seconds": elapsed, "workers_requested": args.workers, "workers_used": args.workers,
                    "cpu_affinity_count": affinity, "worker_cpu_seconds_sum": float(sum(r.get("worker_cpu_seconds", 0.0) for r in rows)),
                    "cpu_only": True, "gpu_used": False, "blas_threads": {k: os.environ.get(k) for k in
                         ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
                    "python": sys.version, "numpy": np.__version__, "scipy": __import__("scipy").__version__,
                    "platform": platform.platform(),
                    "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --seconds 1800 -- .venv/bin/python experiments/interacting_vacuum_screen/interacting_vacuum_null_timing.py --workers 20 --checkpoint-every 25"},
        "sha256": {**input_hashes, "experiments/interacting_vacuum_screen/interacting_vacuum_null_timing.py": sha256(Path(__file__)),
                   "experiments/interacting_vacuum_screen/null_timing_40.jsonl": sha256(JSONL_PATH)},
        "row_count": len(rows),
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    CHECKPOINT_PATH.write_text(json.dumps({"status": result["status"], "completed_rows": len(rows),
                              "requested_rows": len(tasks), "jsonl_sha256": sha256(JSONL_PATH),
                              "result_sha256": sha256(RESULT_PATH)}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "row_count": len(rows), "elapsed_wall_seconds": elapsed,
                      "worker_cpu_seconds_sum": result["runtime"]["worker_cpu_seconds_sum"],
                      "descriptive_tail_count": result["checks"]["tail_count_at_or_above_observed_delta_chi2"]}, indent=2))
    return 0 if len(completed) == len(tasks) else 2


if __name__ == "__main__":
    raise SystemExit(main())

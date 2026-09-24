#!/usr/bin/env python3
"""Resumable matched interacting-vacuum null extension (indices 20 onward).

This driver is prepared for review but must only be run after explicit approval.
It never writes to or modifies the separate 40-mock timing-tranche files.
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

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import interacting_vacuum_profile as iv  # noqa: E402
import interacting_vacuum_null_timing as common  # noqa: E402

SEEDS = {20260924: 200, 20260925: 1000}
FIRST_REMAINING_INDEX = 20
WORKERS_REQUIRED = 20
CHECKPOINT_DEFAULT = 25

JSONL_PATH = HERE / "null_extension_1160.jsonl"
SUMMARY_PATH = HERE / "null_extension_1160_summary.json"
CHECKPOINT_PATH = HERE / "null_extension_1160_checkpoint.json"
PROVENANCE_PATH = ROOT / "work/data_audit2/provenance_check.json"
SEED_MEAN_PATH = ROOT / "experiments/campaign_seed_bao/result.json"
SEED_RUNS = {
    20260924: ROOT / "experiments/bao_robustness/result.json",
    20260925: ROOT / "experiments/bao_robustness/null_extension_1000_seed20260925.json",
}

# Hard guards protect the already-reviewed 40-mock pilot. The extension is
# deliberately output-only and must not overwrite any of these files.
PILOT_EXPECTED_SHA256 = {
    "experiments/interacting_vacuum_screen/interacting_vacuum_null_timing.py": "35f8edf052a7d76eaa63ea1e2390d79123c1111e7c605f8b2b49973f5b79d5d4",
    "experiments/interacting_vacuum_screen/null_timing_40.jsonl": "a5138b807ff21baf776e389848177db617e08b6f59d9318b3ff1e421f9dab09e",
    "experiments/interacting_vacuum_screen/null_timing_40.json": "402ed8fe5c9864be6629aa0022f7e666d3a919674b8f75b3333951363cb27b97",
    "experiments/interacting_vacuum_screen/null_timing_40_checkpoint.json": "3f283ee1083f8353ebb8c00b8ae63695702174704afef8be97946d6bee5961ca",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, document: dict) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temp.replace(path)


def verify_pilot_unchanged() -> dict:
    observed = {}
    for relative, expected in PILOT_EXPECTED_SHA256.items():
        path = ROOT / relative
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"reviewed pilot artifact changed: {relative}")
        observed[relative] = {"sha256": actual, "matches_reviewed_hash": True}
    return observed


def load_verified_vectors():
    provenance = json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))
    if provenance.get("status") != "pass":
        raise RuntimeError("data_audit2 provenance status is not pass")
    sources = provenance["source_hash_audit"]
    streams = provenance["mock_generator_contract"]["streams"]
    for seed, path in SEED_RUNS.items():
        item = sources[str(seed)]
        if not item["data_hashes_match_current"] or not item["run_bao_robustness_hash_matches_current"]:
            raise RuntimeError(f"data/source hash audit failed for seed {seed}")
        if sha256(path) != item["result_sha256"]:
            raise RuntimeError(f"prior seed result changed since audit for seed {seed}")

    seed_doc = json.loads(SEED_MEAN_PATH.read_text(encoding="utf-8"))
    flat = next(model for model in seed_doc["models"] if model["model"] == "lcdm")
    mean = np.asarray(flat["prediction"], dtype=np.float64)
    data = iv.data_load()
    covariance = np.asarray(data.covariance, dtype=np.float64)
    chol = np.linalg.cholesky(covariance)
    reconstructed_mean = float(flat["alpha"]) * iv.unit_prediction(data, float(flat["Omega_m"]), 0.0)[0]
    mean_difference = float(np.max(np.abs(mean-reconstructed_mean)))
    if mean_difference > 1e-10:
        raise RuntimeError(f"fitted null mean recheck failed: {mean_difference:.4g}")

    tasks_by_key = {}
    stream_info = {}
    for seed, size in SEEDS.items():
        rng = np.random.default_rng(seed)
        if type(rng.bit_generator).__name__ != "PCG64":
            raise RuntimeError("default_rng is not PCG64 in this environment")
        normals = rng.standard_normal((size, len(data.z)))
        digest = hashlib.sha256(np.ascontiguousarray(normals, dtype="<f8").tobytes()).hexdigest()
        if digest != streams[str(seed)]["normal_stream_sha256_little_endian_float64"]:
            raise RuntimeError(f"full PCG64 normal digest mismatch for seed {seed}")
        mock_vectors = mean[None, :] + normals @ chol.T
        indices = range(FIRST_REMAINING_INDEX, size)
        for index in indices:
            vector = mock_vectors[index]
            tasks_by_key[(seed, index)] = (seed, index, vector.tolist(), common.vector_sha256(vector))
        stream_info[str(seed)] = {"full_stream_n": size, "selected_indices": [FIRST_REMAINING_INDEX, size-1],
                                  "selected_count": size-FIRST_REMAINING_INDEX,
                                  "normal_stream_sha256_little_endian_float64": digest,
                                  "normal_digest_matches_data_audit2": True}
    return data, covariance, tasks_by_key, stream_info, mean_difference, provenance


def read_existing_rows(tasks_by_key: dict) -> tuple[list[dict], dict[tuple[int, int], list[dict]]]:
    if not JSONL_PATH.exists():
        return [], {}
    raw = JSONL_PATH.read_bytes()
    if raw and not raw.endswith(b"\n"):
        last_newline = raw.rfind(b"\n")
        with JSONL_PATH.open("r+b") as stream:
            stream.truncate(last_newline + 1 if last_newline >= 0 else 0)
        raw = JSONL_PATH.read_bytes()
    rows = []
    by_key: dict[tuple[int, int], list[dict]] = {}
    for line_number, line in enumerate(raw.splitlines(), 1):
        try:
            row = json.loads(line)
        except json.JSONDecodeError as exc:
            raise RuntimeError(f"invalid JSONL row {line_number}: {exc}") from exc
        key = (int(row["seed"]), int(row["index"]))
        if key not in tasks_by_key:
            raise RuntimeError(f"extension JSONL contains out-of-domain key {key}")
        expected_vector = np.asarray(tasks_by_key[key][2], dtype=np.float64)
        stored_vector = np.asarray(row.get("mock_vector"), dtype=np.float64)
        if stored_vector.shape != (13,) or expected_vector.shape != (13,):
            raise RuntimeError(f"mock vector must have exact shape (13,) for {key}")
        if not np.all(np.isfinite(stored_vector)):
            raise RuntimeError(f"stored mock vector contains nonfinite values for {key}")
        recomputed_hash = common.vector_sha256(stored_vector)
        if recomputed_hash != row.get("mock_vector_sha256_little_endian_float64"):
            raise RuntimeError(f"stored mock vector does not match its recorded digest for {key}")
        if recomputed_hash != tasks_by_key[key][3] or not np.array_equal(stored_vector, expected_vector):
            raise RuntimeError(f"stored mock vector differs from regenerated seed/index vector for {key}")
        rows.append(row)
        by_key.setdefault(key, []).append(row)
    for key, attempts in by_key.items():
        if sum(row.get("status") == "complete" for row in attempts) > 1:
            raise RuntimeError(f"duplicate completed fit records for key {key}")
    return rows, by_key


def summarize_rows(rows: list[dict], expected_count: int, *, status: str, run_seconds: float,
                   prior_checkpoint: dict | None, provenance: dict, stream_info: dict,
                   mean_difference: float, input_hashes: dict, pilot_hashes: dict,
                   resume_audit: dict) -> dict:
    by_key: dict[tuple[int, int], list[dict]] = {}
    for row in rows:
        by_key.setdefault((int(row["seed"]), int(row["index"])), []).append(row)
    latest = {key: attempts[-1] for key, attempts in by_key.items()}
    complete = [row for row in latest.values() if row.get("status") == "complete"]
    failed_latest = [row for row in latest.values() if row.get("status") != "complete"]
    deltas = np.asarray([row["selection_delta_chi2_flat_minus_interaction"] for row in complete], dtype=np.float64)
    raw_deltas = np.asarray([row["raw_delta_chi2_flat_minus_interaction"] for row in complete], dtype=np.float64)
    prior_cum = float((prior_checkpoint or {}).get("cumulative_wall_seconds_checkpointed",
                       (prior_checkpoint or {}).get("cumulative_wall_seconds_checkpointed_lower_bound", 0.0)))
    prior_invocations = int((prior_checkpoint or {}).get("checkpointed_invocation_count", 0))
    observed = json.loads((HERE / "result.json").read_text(encoding="utf-8"))
    observed_delta = float(observed["profile_comparison"]["delta_chi2_flat_minus_interacting"])
    result = {
        "status": status,
        "scope": "Extension only: seed 20260924 indices 20-199 and seed 20260925 indices 20-999; preserves the separately reviewed 40-mock timing tranche.",
        "expected_unique_mock_count": expected_count,
        "unique_keys_seen": len(latest), "complete_unique_mock_count": len(complete),
        "failed_latest_records": failed_latest,
        "attempt_record_count": len(rows),
        "null_definition": {"mean": "saved campaign_seed_bao fitted flat-LCDM prediction",
                            "random_generator": "np.random.default_rng(seed), NumPy PCG64, full standard_normal((N,13)) generated separately per seed",
                            "draw": "mock_i = fitted_flat_prediction + normals[i] @ cholesky(full covariance).T",
                            "stream_records": stream_info,
                            "mean_prediction_max_abs_difference_from_independent_profile": mean_difference},
        "contract": {"Omega_m": iv.OM_BOUNDS, "Gamma_over_H0": iv.GAMMA_BOUNDS, "alpha": iv.ALPHA_BOUNDS,
                     "flat": "16 deterministic L-BFGS-B starts per mock",
                     "interaction": "39 deterministic L-BFGS-B starts plus 41x41 grid and bounded L-BFGS-B refinement per mock",
                     "physical_domain": "same profile ODE and baryon/CDM positive-density feasibility test through z=2.33; endpoints re-evaluated at the reported coordinates",
                     "endpoint_policy": "preserve SciPy success/status/fun separately from endpoint physical validity; select lowest valid raw-success multistart endpoint where available (otherwise lowest valid endpoint); allow valid grid refinement candidate as observed; include independently refit exact flat fit as selected nested fallback only if search result is worse, while keeping raw interaction optimum distinct"},
        "pilot_protection": {"unchanged_hash_checks": pilot_hashes, "pilot_output_paths_modified": False},
        "resume_integrity": resume_audit,
        "provenance": {"data_audit2_status": provenance.get("status"),
                       "all_normal_stream_digests_match_audit": all(v["normal_digest_matches_data_audit2"] for v in stream_info.values()),
                       "seed_sources_and_data_hashes_match_audit": True},
        "checks": {
            "tail_count_at_or_above_observed_delta_chi2": int(np.sum(deltas >= observed_delta)),
            "tail_fraction_descriptive_only": float(np.mean(deltas >= observed_delta)) if len(deltas) else None,
            "observed_delta_chi2": observed_delta,
            "flat_raw_scipy_failures": int(sum(r["flat_fit"]["optimizer_endpoint_audit"]["raw_scipy_failure_count"] for r in complete)),
            "flat_invalid_penalty_endpoints": int(sum(r["flat_fit"]["optimizer_endpoint_audit"]["invalid_penalty_endpoint_count"] for r in complete)),
            "interaction_raw_scipy_failures": int(sum(r["interaction_multistart"]["optimizer_endpoint_audit"]["raw_scipy_failure_count"] for r in complete)),
            "interaction_valid_endpoints": int(sum(r["interaction_multistart"]["optimizer_endpoint_audit"]["valid_finite_physical_endpoint_count"] for r in complete)),
            "interaction_invalid_penalty_endpoints": int(sum(r["interaction_multistart"]["optimizer_endpoint_audit"]["invalid_penalty_endpoint_count"] for r in complete)),
            "grid_refinement_raw_scipy_failures": int(sum(not r["interaction_grid_refinement"]["refinement"]["success"] for r in complete)),
            "grid_invalid_evaluations": int(sum(r["interaction_grid_refinement"]["grid_invalid_points"] for r in complete)),
            "nested_flat_fallback_count": int(sum(r["nested_flat_fallback_used"] for r in complete)),
            "selected_delta_summary_descriptive": ({"minimum": float(np.min(deltas)), "median": float(np.median(deltas)),
                "maximum": float(np.max(deltas)), "values": deltas.tolist()} if len(deltas) else None),
            "raw_delta_summary_descriptive": ({"minimum": float(np.min(raw_deltas)), "median": float(np.median(raw_deltas)),
                "maximum": float(np.max(raw_deltas)), "values": raw_deltas.tolist()} if len(raw_deltas) else None),
            "not_posterior_p_value_evidence_or_discovery_significance": True},
        "runtime": {"wall_seconds_this_invocation": run_seconds,
                    "cumulative_wall_seconds_checkpointed_lower_bound": prior_cum + run_seconds,
                    "checkpointed_invocation_count": prior_invocations + 1,
                    "workers": WORKERS_REQUIRED, "cpu_only": True, "gpu_used": False,
                    "blas_thread_environment": {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
                    "python": sys.version, "numpy": np.__version__, "scipy": __import__("scipy").__version__,
                    "platform": platform.platform(),
                    "worker_cpu_seconds_sum_all_attempt_rows": float(sum(r.get("worker_cpu_seconds", 0.0) for r in rows)),
                    "command": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --seconds 3600 -- .venv/bin/python experiments/interacting_vacuum_screen/interacting_vacuum_null_extension.py --workers 20 --checkpoint-every 25"},
        "sha256": input_hashes,
        "jsonl_sha256": sha256(JSONL_PATH) if JSONL_PATH.exists() else None,
        "summary_code_sha256": sha256(Path(__file__)),
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=WORKERS_REQUIRED)
    parser.add_argument("--checkpoint-every", type=int, default=CHECKPOINT_DEFAULT)
    parser.add_argument("--preflight", action="store_true", help="validate hashes/streams and exit without any mock fitting")
    args = parser.parse_args()
    if args.workers != WORKERS_REQUIRED:
        raise ValueError("locked extension contract requires exactly 20 worker processes")
    if args.checkpoint_every != CHECKPOINT_DEFAULT:
        raise ValueError("locked extension contract checkpoints every 25 completed rows")
    affinity = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    if affinity < WORKERS_REQUIRED:
        raise RuntimeError(f"20 workers requested but only {affinity} CPUs available in affinity")
    pilot_hashes = verify_pilot_unchanged()
    data, covariance, tasks_by_key, stream_info, mean_difference, provenance = load_verified_vectors()
    input_paths = [ROOT / "context/data/desi_dr2_mean.txt", ROOT / "context/data/desi_dr2_cov.txt",
                   SEED_MEAN_PATH, *SEED_RUNS.values(), PROVENANCE_PATH,
                   ROOT / "scripts/background_bao.py", HERE / "interacting_vacuum_profile.py",
                   HERE / "interacting_vacuum_null_timing.py", HERE / "result.json", Path(__file__).resolve()]
    input_hashes = {str(path.relative_to(ROOT)): sha256(path) for path in input_paths}
    if input_hashes["experiments/interacting_vacuum_screen/interacting_vacuum_null_timing.py"] != PILOT_EXPECTED_SHA256["experiments/interacting_vacuum_screen/interacting_vacuum_null_timing.py"]:
        raise RuntimeError("shared mock-fit helper changed since pilot review")
    if input_hashes["experiments/interacting_vacuum_screen/interacting_vacuum_profile.py"] != "ae99ae789aa02885667e27aa2cbbaa0b035323e69d1a4bf7bbd3848949fc5282":
        raise RuntimeError("interacting profile fitter changed since pilot review")
    if input_hashes["experiments/campaign_seed_bao/result.json"] != "92662f3ee42500d2be34d9f87b485fe56c82e9550fa674e42f4dec2bf2cab1f0":
        raise RuntimeError("fitted null mean source changed since pilot review")
    if input_hashes["experiments/interacting_vacuum_screen/result.json"] != "e347b8934b2afdc23d23ef7d39c70f963b73888d17137f55caf7e51af87d43ae":
        raise RuntimeError("observed fit result changed since pilot review")

    if args.preflight:
        print(json.dumps({"status": "provenance_preflight_passed", "fit_started": False,
                          "remaining_mock_count": len(tasks_by_key), "stream_records": stream_info,
                          "fitted_mean_max_abs_recheck_difference": mean_difference,
                          "data_audit2_status": provenance.get("status"),
                          "pilot_artifacts_unchanged": pilot_hashes,
                          "input_sha256": input_hashes}, indent=2))
        return 0

    old_rows, old_by_key = read_existing_rows(tasks_by_key)
    existing_rows_verified_count = len(old_rows)
    completed_keys = {key for key, attempts in old_by_key.items() if any(row.get("status") == "complete" for row in attempts)}
    tasks = [tasks_by_key[key] for key in tasks_by_key if key not in completed_keys]
    previous_checkpoint = None
    jsonl_ahead_of_prior_checkpoint = False
    if CHECKPOINT_PATH.exists():
        previous_checkpoint = json.loads(CHECKPOINT_PATH.read_text(encoding="utf-8"))
        if previous_checkpoint.get("jsonl_sha256") and previous_checkpoint["jsonl_sha256"] != sha256(JSONL_PATH):
            # Each result row is fsync'd before the checkpoint cadence. After an
            # interruption the JSONL may legitimately be ahead of its checkpoint;
            # reconcile rows by seed/index/vector hash and resume missing keys.
            previous_checkpoint["jsonl_ahead_of_checkpoint_on_resume"] = True
            jsonl_ahead_of_prior_checkpoint = True
    resume_audit = {"prior_checkpoint_present": previous_checkpoint is not None,
                    "jsonl_ahead_of_prior_checkpoint_on_resume": jsonl_ahead_of_prior_checkpoint,
                    "existing_jsonl_attempt_rows_vector_reverified": existing_rows_verified_count,
                    "vector_shape_digest_and_exact_seed_index_value_checked": True,
                    "prior_checkpoint_jsonl_sha256": (None if previous_checkpoint is None else previous_checkpoint.get("jsonl_sha256"))}

    session_started = time.monotonic()
    new_rows = []
    processed = 0
    attempts_per_key = {key: len(value) for key, value in old_by_key.items()}
    mode = "a" if JSONL_PATH.exists() else "x"
    with JSONL_PATH.open(mode, encoding="utf-8") as output_stream:
        if tasks:
            with ProcessPoolExecutor(max_workers=args.workers, initializer=common.worker_init,
                                     initargs=(data.z.tolist(), data.observable.tolist(), covariance.tolist())) as pool:
                futures = [pool.submit(common.fit_mock, task) for task in tasks]
                for future in as_completed(futures):
                    row = future.result()
                    key = (int(row["seed"]), int(row["index"]))
                    row["extension_attempt_number"] = attempts_per_key.get(key, 0) + 1
                    row["extension_only"] = True
                    output_stream.write(json.dumps(row, separators=(",", ":"), allow_nan=False) + "\n")
                    output_stream.flush()
                    os.fsync(output_stream.fileno())
                    new_rows.append(row)
                    old_rows.append(row)
                    attempts_per_key[key] = attempts_per_key.get(key, 0) + 1
                    processed += 1
                    if processed % args.checkpoint_every == 0:
                        elapsed = time.monotonic() - session_started
                        current = summarize_rows(old_rows, len(tasks_by_key), status="running", run_seconds=elapsed,
                            prior_checkpoint=previous_checkpoint, provenance=provenance, stream_info=stream_info,
                            mean_difference=mean_difference, input_hashes=input_hashes, pilot_hashes=pilot_hashes,
                            resume_audit=resume_audit)
                        atomic_json(SUMMARY_PATH, current)
                        atomic_json(CHECKPOINT_PATH, {"status": "running", "processed_this_invocation": processed,
                            "requested_this_invocation": len(tasks), "complete_unique_count": current["complete_unique_mock_count"],
                            "required_unique_count": len(tasks_by_key), "checkpointed_invocation_count": current["runtime"]["checkpointed_invocation_count"],
                            "cumulative_wall_seconds_checkpointed_lower_bound": current["runtime"]["cumulative_wall_seconds_checkpointed_lower_bound"],
                            "jsonl_ahead_of_prior_checkpoint_on_resume": jsonl_ahead_of_prior_checkpoint,
                            "existing_jsonl_attempt_rows_vector_reverified": existing_rows_verified_count,
                            "jsonl_sha256": sha256(JSONL_PATH), "summary_sha256": sha256(SUMMARY_PATH),
                            "pilot_artifacts_unchanged": pilot_hashes})
    elapsed = time.monotonic() - session_started
    # Verify the prior pilot still matches the exact review hashes before publishing summary.
    final_pilot_hashes = verify_pilot_unchanged()
    final = summarize_rows(old_rows, len(tasks_by_key), status="complete" if len({(r["seed"],r["index"]) for r in old_rows if r.get("status")=="complete"}) == len(tasks_by_key) else "partial",
                           run_seconds=elapsed, prior_checkpoint=previous_checkpoint, provenance=provenance,
                           stream_info=stream_info, mean_difference=mean_difference,
                           input_hashes=input_hashes, pilot_hashes=final_pilot_hashes,
                           resume_audit=resume_audit)
    final["runtime"]["cpu_affinity_count"] = affinity
    atomic_json(SUMMARY_PATH, final)
    atomic_json(CHECKPOINT_PATH, {"status": final["status"], "processed_this_invocation": processed,
        "requested_this_invocation": len(tasks), "complete_unique_count": final["complete_unique_mock_count"],
        "required_unique_count": len(tasks_by_key), "checkpointed_invocation_count": final["runtime"]["checkpointed_invocation_count"],
        "cumulative_wall_seconds_checkpointed_lower_bound": final["runtime"]["cumulative_wall_seconds_checkpointed_lower_bound"],
        "jsonl_ahead_of_prior_checkpoint_on_resume": jsonl_ahead_of_prior_checkpoint,
        "existing_jsonl_attempt_rows_vector_reverified": existing_rows_verified_count,
        "jsonl_sha256": sha256(JSONL_PATH), "summary_sha256": sha256(SUMMARY_PATH),
        "pilot_artifacts_unchanged": final_pilot_hashes})
    print(json.dumps({"status": final["status"], "new_rows": processed,
                      "complete_unique_mock_count": final["complete_unique_mock_count"],
                      "expected_unique_mock_count": len(tasks_by_key), "elapsed_wall_seconds": elapsed}, indent=2))
    return 0 if final["status"] == "complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())

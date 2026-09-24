#!/usr/bin/env python3
"""Independent audit of the completed 1,200-row interacting-vacuum null screen.

This reads and verifies saved records; it does not optimize/refit any realization.
It independently reevaluates eight deterministic rows at their archived fit
parameters using the adaptive scalar-distance implementation in work/theory4.

Run only through the project's time-bounded wrapper, with one numerical thread:
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --seconds 900 -- \
  .venv/bin/python work/data_audit3/audit_interaction_null.py
"""
from __future__ import annotations

from collections import Counter
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.linalg import solve_triangular
from scipy.stats import beta as beta_distribution


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PILOT_ROWS = ROOT / "experiments/interacting_vacuum_screen/null_timing_40.jsonl"
PILOT_RESULT = ROOT / "experiments/interacting_vacuum_screen/null_timing_40.json"
PILOT_RECORD = ROOT / "experiments/interacting_vacuum_screen/null_timing_40_record.json"
PILOT_CHECKPOINT = ROOT / "experiments/interacting_vacuum_screen/null_timing_40_checkpoint.json"
EXTENSION_ROWS = ROOT / "experiments/interacting_vacuum_screen/null_extension_1160.jsonl"
EXTENSION_SUMMARY = ROOT / "experiments/interacting_vacuum_screen/null_extension_1160_summary.json"
EXTENSION_CHECKPOINT = ROOT / "experiments/interacting_vacuum_screen/null_extension_1160_checkpoint.json"
MEAN_PATH = ROOT / "context/data/desi_dr2_mean.txt"
COV_PATH = ROOT / "context/data/desi_dr2_cov.txt"
SEED_MEAN_PATH = ROOT / "experiments/campaign_seed_bao/result.json"
DATA_AUDIT2_PATH = ROOT / "work/data_audit2/provenance_check.json"
OBSERVED_PATH = ROOT / "experiments/interacting_vacuum_screen/result.json"
EXTENSION_SCRIPT = ROOT / "experiments/interacting_vacuum_screen/interacting_vacuum_null_extension.py"
TIMING_SCRIPT = ROOT / "experiments/interacting_vacuum_screen/interacting_vacuum_null_timing.py"
THEORY4_SCRIPT = ROOT / "work/theory4/audit_interacting_identifiability.py"
OUTPUT_PATH = HERE / "interaction_null_audit.json"

SEED_SIZES = {20260924: 200, 20260925: 1000}
PILOT_LIMIT = 20
OBSERVED_DELTA_REFERENCE = 1.5171630696615157
INVALID_PENALTY = 1e80
FIT_ABS_TOL = 2e-7
PREDICTION_ABS_TOL = 2e-5


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def vector_digest(vector: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(vector, dtype="<f8").tobytes()).hexdigest()


def read_mean() -> tuple[np.ndarray, np.ndarray, list[str]]:
    rows = []
    for line_no, line in enumerate(MEAN_PATH.read_text(encoding="utf-8").splitlines(), 1):
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) != 3:
            raise ValueError(f"mean file line {line_no}: expected three columns")
        rows.append((float(fields[0]), float(fields[1]), fields[2]))
    return (np.asarray([r[0] for r in rows], dtype=np.float64),
            np.asarray([r[1] for r in rows], dtype=np.float64),
            [r[2] for r in rows])


def load_theory4():
    spec = importlib.util.spec_from_file_location("theory4_scalar_audit", THEORY4_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import independent theory4 implementation at {THEORY4_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_jsonl(path: Path):
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid JSON: {exc}") from exc


def close_enough(a: float, b: float, atol: float = 2e-9) -> bool:
    return math.isfinite(float(a)) and math.isfinite(float(b)) and abs(float(a) - float(b)) <= atol


def add_attempt_audit(fit: dict, prefix: str, counts: Counter, row_mismatch: Counter) -> None:
    attempts = fit.get("optimizer_attempts")
    audit = fit.get("optimizer_endpoint_audit")
    if not isinstance(attempts, list) or not isinstance(audit, dict):
        row_mismatch[f"{prefix}_missing_attempt_audit"] += 1
        return
    total = len(attempts)
    raw_success = sum(bool(a.get("success")) for a in attempts)
    valid = [bool(a.get("endpoint_valid_physical")) and math.isfinite(float(a.get("chi2", math.nan)))
             and float(a.get("chi2", math.inf)) < INVALID_PENALTY for a in attempts]
    valid_count = sum(valid)
    counts[f"{prefix}_attempts"] += total
    counts[f"{prefix}_raw_success"] += raw_success
    counts[f"{prefix}_raw_failure"] += total - raw_success
    counts[f"{prefix}_valid_physical"] += valid_count
    counts[f"{prefix}_invalid_or_penalty"] += total - valid_count
    counts[f"{prefix}_success_and_valid"] += sum(bool(a.get("success")) and is_valid for a, is_valid in zip(attempts, valid))
    counts[f"{prefix}_failure_but_valid"] += sum(not bool(a.get("success")) and is_valid for a, is_valid in zip(attempts, valid))
    counts[f"{prefix}_success_but_invalid"] += sum(bool(a.get("success")) and not is_valid for a, is_valid in zip(attempts, valid))
    counts[f"{prefix}_failure_and_invalid"] += sum(not bool(a.get("success")) and not is_valid for a, is_valid in zip(attempts, valid))
    counts[f"{prefix}_invalid_penalty_endpoints"] += sum(
        (not is_valid) and math.isfinite(float(a.get("chi2", math.nan)))
        and float(a.get("chi2", 0.0)) >= INVALID_PENALTY for a, is_valid in zip(attempts, valid)
    )
    expected = {
        "start_count": total,
        "raw_scipy_success_count": raw_success,
        "raw_scipy_failure_count": total - raw_success,
        "valid_finite_physical_endpoint_count": valid_count,
        "invalid_penalty_endpoint_count": total - valid_count,
        "raw_success_at_invalid_penalty_count": sum(bool(a.get("success")) and not is_valid for a, is_valid in zip(attempts, valid)),
        "valid_endpoint_success_count": sum(bool(a.get("success")) and is_valid for a, is_valid in zip(attempts, valid)),
        "valid_endpoint_failure_count": sum(not bool(a.get("success")) and is_valid for a, is_valid in zip(attempts, valid)),
    }
    for name, expected_value in expected.items():
        if audit.get(name) != expected_value:
            row_mismatch[f"{prefix}_audit_{name}"] += 1
    for attempt in attempts:
        for field in ("chi2", "scipy_fun"):
            if not math.isfinite(float(attempt.get(field, math.nan))):
                row_mismatch[f"{prefix}_nonfinite_{field}"] += 1
        if not isinstance(attempt.get("success"), bool) or not isinstance(attempt.get("endpoint_valid_physical"), bool):
            row_mismatch[f"{prefix}_missing_separate_status_flag"] += 1
        if not isinstance(attempt.get("status"), int):
            row_mismatch[f"{prefix}_missing_integer_scipy_status"] += 1


def scalar_unit_prediction(theory4, z: np.ndarray, labels: list[str], omega_m: float, g: float,
                            preferred_fb: float = 0.16) -> tuple[np.ndarray, dict]:
    """Use theory4's independent ODE plus adaptive scalar-quadrature geometry."""
    theory4.LABELS = labels
    try:
        return theory4.geometry_and_split(omega_m, g, z, f_b=preferred_fb)
    except theory4.InvalidHistory as first_error:
        # Background distances do not depend on the baryon/CDM split. If the
        # illustrative f_b=0.16 check alone fails, retry at an explicit positive
        # split below the independently sampled viability limit; retain both facts.
        sol, minima = theory4.scalar_history(omega_m, g, float(np.max(z)))
        zcheck = np.linspace(0.0, float(np.max(z)), 1001, dtype=np.float64)
        ucheck = np.log1p(zcheck)
        matter_check, _vacuum_check = np.asarray(sol(ucheck), dtype=np.float64)
        ratio = matter_check / (omega_m * np.exp(3.0 * ucheck))
        fb_max = float(min(1.0, np.min(ratio)))
        if not (math.isfinite(fb_max) and fb_max > 0.0):
            raise RuntimeError("theory4 scalar history has no positive baryon split") from first_error
        witness_fb = min(preferred_fb, 0.5 * fb_max)
        if witness_fb <= 0.0:
            raise RuntimeError("nonpositive theory4 baryon witness") from first_error
        q, phys = theory4.geometry_and_split(omega_m, g, z, f_b=witness_fb)
        phys["preferred_fb_016_physical"] = False
        phys["adaptive_distance_witness_fb"] = witness_fb
        phys["preferred_fb_016_failure"] = str(first_error)
        return q, phys


def evaluate_saved_fit(theory4, z, labels, chol, mock, fit: dict, role: str) -> dict:
    om = float(fit["Omega_m"])
    g = float(fit.get("Gamma_over_H0", 0.0))
    alpha = float(fit["alpha"])
    q, physics = scalar_unit_prediction(theory4, z, labels, om, g)
    pred = alpha * q
    residual_white = solve_triangular(chol, mock - pred, lower=True, check_finite=False)
    chi2 = float(np.dot(residual_white, residual_white))
    saved_prediction = np.asarray(fit["prediction"], dtype=np.float64)
    saved_score = float(fit["chi2"])
    return {
        "role": role,
        "parameters": {"Omega_m": om, "Gamma_over_H0": g, "alpha": alpha},
        "independent_chi2": chi2,
        "stored_chi2": saved_score,
        "chi2_difference_independent_minus_stored": chi2 - saved_score,
        "independent_prediction": pred.tolist(),
        "max_abs_prediction_difference": float(np.max(np.abs(pred - saved_prediction))),
        "unit_prediction_max_abs_difference": float(np.max(np.abs(q - np.asarray(fit.get("unit_prediction", q), dtype=np.float64)))),
        "theory4_physicality": physics,
    }


def main() -> int:
    started = time.monotonic()
    mean_doc = json.loads(SEED_MEAN_PATH.read_text(encoding="utf-8"))
    flat_seed = next(m for m in mean_doc["models"] if m["model"] == "lcdm")
    null_mean = np.asarray(flat_seed["prediction"], dtype=np.float64)
    z, published_mean, labels = read_mean()
    covariance = np.loadtxt(COV_PATH, dtype=np.float64, ndmin=2)
    if len(z) != 13 or covariance.shape != (13, 13) or null_mean.shape != (13,):
        raise ValueError("expected a 13-row mean, 13x13 covariance, and 13-value null mean")
    if not np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-12):
        raise ValueError("full covariance is not symmetric")
    chol = np.linalg.cholesky(covariance)
    eig = np.linalg.eigvalsh(covariance)

    a2 = json.loads(DATA_AUDIT2_PATH.read_text(encoding="utf-8"))
    if a2.get("status") != "pass":
        raise RuntimeError("prior data_audit2 provenance contract does not pass")
    stream_info = a2["mock_generator_contract"]["streams"]
    rng_data = a2["mock_generator_contract"]
    observed = json.loads(OBSERVED_PATH.read_text(encoding="utf-8"))
    observed_delta = float(observed["profile_comparison"]["delta_chi2_flat_minus_interacting"])
    if not close_enough(observed_delta, OBSERVED_DELTA_REFERENCE, atol=1e-13):
        raise RuntimeError("observed-delta reference does not match the user-supplied locked value")

    pilot_summary = json.loads(PILOT_RESULT.read_text(encoding="utf-8"))
    pilot_record = json.loads(PILOT_RECORD.read_text(encoding="utf-8"))
    pilot_checkpoint = json.loads(PILOT_CHECKPOINT.read_text(encoding="utf-8"))
    extension_summary = json.loads(EXTENSION_SUMMARY.read_text(encoding="utf-8"))
    extension_checkpoint = json.loads(EXTENSION_CHECKPOINT.read_text(encoding="utf-8"))

    # Independently regenerate both complete PCG64 streams and all Cholesky mocks.
    rng_rebuild = {}
    regenerated_vectors = {}
    for seed, count in SEED_SIZES.items():
        rng = np.random.default_rng(seed)
        if type(rng.bit_generator).__name__ != "PCG64":
            raise RuntimeError("NumPy default_rng is not using PCG64")
        normals = rng.standard_normal((count, 13))
        normal_hash = hashlib.sha256(np.ascontiguousarray(normals, dtype="<f8").tobytes()).hexdigest()
        saved = stream_info[str(seed)]
        if normal_hash != saved["normal_stream_sha256_little_endian_float64"]:
            raise RuntimeError(f"full PCG64 stream digest differs from data_audit2 for seed {seed}")
        regenerated_vectors[seed] = null_mean[None, :] + normals @ chol.T
        rng_rebuild[str(seed)] = {
            "count": count, "bit_generator": "PCG64", "shape": list(normals.shape),
            "normal_stream_sha256_little_endian_float64": normal_hash,
            "matches_passing_data_audit2_contract": True,
        }

    expected_pilot = {(seed, i) for seed in SEED_SIZES for i in range(PILOT_LIMIT)}
    expected_extension = {(20260924, i) for i in range(20, 200)} | {(20260925, i) for i in range(20, 1000)}
    expected_all = {(seed, i) for seed, n in SEED_SIZES.items() for i in range(n)}
    seen_pilot: set[tuple[int, int]] = set()
    seen_extension: set[tuple[int, int]] = set()
    sample_rows: dict[str, list[tuple[tuple[int, int], dict]]] = {"exceedance": [], "non_exceedance": []}
    all_deltas: list[float] = []
    all_raw_deltas: list[float] = []
    selected_g: list[float] = []
    raw_g: list[float] = []
    flat_omega: list[float] = []
    counters: Counter = Counter()
    mismatches: Counter = Counter()
    max_vector_diff = 0.0
    max_delta_formula_error = 0.0
    max_raw_delta_formula_error = 0.0
    record_statuses = Counter()
    candidate_sources = Counter()
    selected_bound_hits = Counter()
    selected_param_bound_counts = Counter()
    interaction_raw_param_bound_counts = Counter()
    flat_param_bound_counts = Counter()
    selected_success = Counter()
    selected_gamma_min = math.inf
    selected_gamma_max = -math.inf

    def examine(row: dict, source: str) -> None:
        nonlocal max_vector_diff, max_delta_formula_error, max_raw_delta_formula_error
        nonlocal selected_gamma_min, selected_gamma_max
        key = (int(row["seed"]), int(row["index"]))
        target_keys = seen_pilot if source == "pilot" else seen_extension
        if key in target_keys:
            mismatches[f"duplicate_key_within_{source}"] += 1
        target_keys.add(key)
        if key not in (expected_pilot if source == "pilot" else expected_extension):
            mismatches[f"out_of_domain_key_in_{source}"] += 1

        stored_vector = np.asarray(row.get("mock_vector"), dtype=np.float64)
        if stored_vector.shape != (13,):
            mismatches["mock_vector_bad_shape"] += 1
            return
        expected_vector = regenerated_vectors[key[0]][key[1]]
        diff = float(np.max(np.abs(stored_vector - expected_vector)))
        max_vector_diff = max(max_vector_diff, diff)
        actual_vhash = vector_digest(stored_vector)
        if (not np.array_equal(stored_vector, expected_vector)
                or actual_vhash != row.get("mock_vector_sha256_little_endian_float64")):
            mismatches["mock_vector_value_or_hash"] += 1
        if not np.all(np.isfinite(stored_vector)):
            mismatches["nonfinite_mock_vector"] += 1
        record_statuses[str(row.get("status"))] += 1
        if row.get("status") != "complete":
            mismatches["noncomplete_record"] += 1
        for field in ("flat_fit", "interaction_multistart", "interaction_grid_refinement",
                      "interaction_raw_search_optimum", "selected_interaction_candidate"):
            if field not in row:
                mismatches[f"missing_{field}"] += 1

        flat = row["flat_fit"]
        multi = row["interaction_multistart"]
        grid = row["interaction_grid_refinement"]
        raw = row["interaction_raw_search_optimum"]
        selected = row["selected_interaction_candidate"]
        for name, value in (("flat", flat), ("multistart", multi), ("raw", raw), ("selected", selected)):
            for scalar in ("chi2", "alpha", "Omega_m"):
                if scalar not in value or not math.isfinite(float(value[scalar])):
                    mismatches[f"{name}_nonfinite_{scalar}"] += 1
            if float(value.get("chi2", math.nan)) < 0.0:
                mismatches[f"{name}_negative_chi2"] += 1
        for name, value in (("flat", flat), ("raw", raw), ("selected", selected)):
            pred = np.asarray(value.get("prediction", []), dtype=np.float64)
            if pred.shape != (13,) or not np.all(np.isfinite(pred)):
                mismatches[f"{name}_prediction_shape_or_nonfinite"] += 1
        for value in (multi, raw, selected):
            if not math.isfinite(float(value.get("Gamma_over_H0", math.nan))):
                mismatches["interaction_gamma_nonfinite"] += 1

        add_attempt_audit(flat, "flat", counters, mismatches)
        add_attempt_audit(multi, "interaction", counters, mismatches)
        for hit in flat.get("parameter_bound_hits", []):
            flat_param_bound_counts[str(hit.get("parameter")) + ":" + str(hit.get("side"))] += 1
        for hit in multi.get("parameter_bound_hits", []):
            interaction_raw_param_bound_counts[str(hit.get("parameter")) + ":" + str(hit.get("side"))] += 1

        ref = grid.get("refinement", {})
        grid_valid = bool(ref.get("endpoint_valid_physical"))
        grid_success = bool(ref.get("success"))
        counters["grid_refinement_total"] += 1
        counters["grid_refinement_physically_valid"] += int(grid_valid)
        counters["grid_refinement_raw_success"] += int(grid_success)
        counters["grid_refinement_raw_failure"] += int(not grid_success)
        counters["grid_refinement_success_but_invalid"] += int(grid_success and not grid_valid)
        counters["grid_refinement_failure_but_valid"] += int((not grid_success) and grid_valid)
        if (not math.isfinite(float(ref.get("endpoint_chi2", math.nan)))
                or not math.isfinite(float(ref.get("scipy_fun", math.nan)))):
            mismatches["grid_refinement_nonfinite_score"] += 1
        candidate = grid.get("candidate")
        if grid_valid != (candidate is not None):
            mismatches["grid_validity_candidate_presence"] += 1
        if candidate is not None:
            candidate_sources[str(candidate.get("source"))] += 1
            if not close_enough(float(candidate["chi2"]), float(ref["endpoint_chi2"]), atol=2e-8):
                mismatches["grid_candidate_endpoint_score"] += 1

        for metric_name in ("grid_invalid_points",):
            value = grid.get(metric_name)
            if not isinstance(value, int) or value < 0:
                mismatches[f"invalid_{metric_name}"] += 1
            else:
                counters["grid_invalid_evaluations"] += value

        raw_candidates = [{"chi2": float(multi["chi2"]), "source": "39-start bounded multistart"}]
        if candidate is not None:
            raw_candidates.append({"chi2": float(candidate["chi2"]), "source": candidate["source"]})
        expected_raw = min(raw_candidates, key=lambda c: c["chi2"])
        if (not close_enough(float(raw["chi2"]), expected_raw["chi2"], atol=2e-9)
                or raw.get("source") != expected_raw["source"]):
            mismatches["raw_candidate_minimum_or_source"] += 1

        flat_score = float(flat["chi2"])
        raw_score = float(raw["chi2"])
        selected_score = float(selected["chi2"])
        raw_delta = float(row["raw_delta_chi2_flat_minus_interaction"])
        delta = float(row["selection_delta_chi2_flat_minus_interaction"])
        expected_raw_delta = flat_score - raw_score
        expected_delta = flat_score - selected_score
        max_raw_delta_formula_error = max(max_raw_delta_formula_error, abs(raw_delta - expected_raw_delta))
        max_delta_formula_error = max(max_delta_formula_error, abs(delta - expected_delta))
        if not close_enough(raw_delta, expected_raw_delta, atol=2e-8):
            mismatches["raw_delta_formula"] += 1
        if not close_enough(delta, expected_delta, atol=2e-8):
            mismatches["selection_delta_formula"] += 1
        if delta < -2e-8:
            mismatches["selected_delta_negative_after_fallback"] += 1
        fallback = bool(row["nested_flat_fallback_used"])
        expected_fallback = raw_score > flat_score
        counters["nested_fallback_total"] += int(fallback)
        if fallback != expected_fallback:
            mismatches["nested_fallback_decision"] += 1
        if fallback:
            if (not close_enough(selected_score, flat_score, atol=2e-9)
                    or not close_enough(float(selected.get("Gamma_over_H0", math.nan)), 0.0, atol=1e-14)
                    or not str(selected.get("source", "")).startswith("exact refitted nested flat")):
                mismatches["nested_fallback_payload"] += 1
        elif (not close_enough(selected_score, raw_score, atol=2e-9)
              or selected.get("source") != raw.get("source")):
            mismatches["selected_not_raw_best_without_fallback"] += 1

        all_deltas.append(delta)
        all_raw_deltas.append(raw_delta)
        flat_omega.append(float(flat["Omega_m"]))
        selected_gamma = float(selected.get("Gamma_over_H0", 0.0))
        raw_gamma = float(raw.get("Gamma_over_H0", math.nan))
        selected_g.append(selected_gamma)
        raw_g.append(raw_gamma)
        selected_gamma_min = min(selected_gamma_min, selected_gamma)
        selected_gamma_max = max(selected_gamma_max, selected_gamma)
        counters["selected_gamma_negative"] += int(selected_gamma < -1e-12)
        counters["selected_gamma_zero"] += int(abs(selected_gamma) <= 1e-12)
        counters["selected_gamma_positive"] += int(selected_gamma > 1e-12)
        counters["raw_gamma_negative"] += int(raw_gamma < -1e-12)
        counters["raw_gamma_zero"] += int(abs(raw_gamma) <= 1e-12)
        counters["raw_gamma_positive"] += int(raw_gamma > 1e-12)
        selected_success[str(bool(selected.get("success"))) + ":" + str(selected.get("source"))] += 1
        for hit in selected.get("parameter_bound_hits", []):
            label = str(hit.get("parameter")) + ":" + str(hit.get("side"))
            selected_param_bound_counts[label] += 1
            selected_bound_hits[label] += 1
        all_values = sample_rows["exceedance" if delta >= observed_delta else "non_exceedance"]
        all_values.append((key, row))
        all_values.sort(key=lambda item: item[0])
        del all_values[4:]

    pilot_count = 0
    for row in read_jsonl(PILOT_ROWS):
        pilot_count += 1
        examine(row, "pilot")
    extension_count = 0
    for row in read_jsonl(EXTENSION_ROWS):
        extension_count += 1
        examine(row, "extension")

    if seen_pilot != expected_pilot:
        mismatches["pilot_key_set"] += len(expected_pilot.symmetric_difference(seen_pilot)) or 1
    if seen_extension != expected_extension:
        mismatches["extension_key_set"] += len(expected_extension.symmetric_difference(seen_extension)) or 1
    if seen_pilot & seen_extension:
        mismatches["pilot_extension_overlap"] += len(seen_pilot & seen_extension)
    if seen_pilot | seen_extension != expected_all:
        mismatches["pooled_key_set"] += len(expected_all.symmetric_difference(seen_pilot | seen_extension)) or 1

    # Verify the result files and the hard pilot-protection/input hashes.
    hash_audit = {"pilot_record_sha256": sha256_path(PILOT_RECORD), "extension_checkpoint_sha256": sha256_path(EXTENSION_CHECKPOINT)}
    pilot_expected_hashes = pilot_summary.get("sha256", {})
    pilot_hash_checks = {}
    for relative, expected_hash in pilot_expected_hashes.items():
        actual_path = ROOT / relative
        actual_hash = sha256_path(actual_path) if actual_path.is_file() else None
        pilot_hash_checks[relative] = {"expected": expected_hash, "actual": actual_hash, "matches": actual_hash == expected_hash}
    ext_expected_hashes = extension_summary.get("sha256", {})
    extension_hash_checks = {}
    for relative, expected_hash in ext_expected_hashes.items():
        actual_path = ROOT / relative
        actual_hash = sha256_path(actual_path) if actual_path.is_file() else None
        extension_hash_checks[relative] = {"expected": expected_hash, "actual": actual_hash, "matches": actual_hash == expected_hash}
    pilot_record_ok = (
        pilot_record.get("status") == "complete_timing_pilot_only"
        and pilot_record.get("rows_sha256") == sha256_path(PILOT_ROWS)
        and pilot_record.get("result_sha256") == sha256_path(PILOT_RESULT)
        and pilot_record.get("rows") == str(PILOT_ROWS.relative_to(ROOT))
    )
    pilot_report_ok = (
        pilot_summary.get("status") == "complete_40_mock_timing_tranche"
        and int(pilot_summary.get("row_count", -1)) == 40
        and int(pilot_summary.get("checks", {}).get("completed_mock_count", -1)) == 40
    )
    pilot_hashes_ok = all(item["matches"] for item in pilot_hash_checks.values())
    extension_hashes_ok = all(item["matches"] for item in extension_hash_checks.values())
    pprotect = extension_summary.get("pilot_protection", {})
    pprotect_checks = pprotect.get("unchanged_hash_checks", {})
    pprotect_ok = bool(pprotect) and pprotect.get("pilot_output_paths_modified") is False and all(
        bool(item.get("matches_reviewed_hash")) and sha256_path(ROOT / path) == item.get("sha256")
        for path, item in pprotect_checks.items()
    )
    data_audit2_hash = sha256_path(DATA_AUDIT2_PATH)
    provenance_hash_ok = (
        extension_summary.get("provenance", {}).get("data_audit2_status") == "pass"
        and extension_summary.get("sha256", {}).get("work/data_audit2/provenance_check.json") == data_audit2_hash
        and a2.get("status") == "pass"
    )
    extension_rows_hash_ok = (
        extension_summary.get("jsonl_sha256") == sha256_path(EXTENSION_ROWS)
        and extension_checkpoint.get("jsonl_sha256") == sha256_path(EXTENSION_ROWS)
        and extension_checkpoint.get("status") == "complete"
        and int(extension_checkpoint.get("complete_unique_count", -1)) == 1160
        and extension_checkpoint.get("summary_sha256") == sha256_path(EXTENSION_SUMMARY)
    )

    deltas = np.asarray(all_deltas, dtype=np.float64)
    raw_deltas = np.asarray(all_raw_deltas, dtype=np.float64)
    if len(deltas) != 1200 or not np.all(np.isfinite(deltas)) or not np.all(np.isfinite(raw_deltas)):
        mismatches["pooled_delta_count_or_finiteness"] += 1
    k = int(np.count_nonzero(deltas >= observed_delta))
    n = int(len(deltas))
    cp = [float(beta_distribution.ppf(0.025, k, n-k+1)) if k else 0.0,
          float(beta_distribution.ppf(0.975, k+1, n-k)) if k < n else 1.0]
    q50, q95, q99 = map(float, np.quantile(deltas, [0.50, 0.95, 0.99], method="linear"))
    q_raw = list(map(float, np.quantile(raw_deltas, [0.50, 0.95, 0.99], method="linear")))

    endpoint_summary_values = {
        "flat_raw_scipy_failures": int(pilot_summary["checks"]["flat_optimizer_raw_scipy_failures_total"])
            + int(extension_summary["checks"]["flat_raw_scipy_failures"]),
        "flat_invalid_penalty_endpoints": int(pilot_summary["checks"]["flat_invalid_penalty_endpoints_total"])
            + int(extension_summary["checks"]["flat_invalid_penalty_endpoints"]),
        "interaction_raw_scipy_failures": int(pilot_summary["checks"]["interaction_multistart_raw_scipy_failures_total"])
            + int(extension_summary["checks"]["interaction_raw_scipy_failures"]),
        "interaction_invalid_penalty_endpoints": int(pilot_summary["checks"]["interaction_invalid_penalty_endpoints_total"])
            + int(extension_summary["checks"]["interaction_invalid_penalty_endpoints"]),
        "grid_refinement_raw_scipy_failures": int(pilot_summary["checks"]["grid_refinement_raw_scipy_failures_total"])
            + int(extension_summary["checks"]["grid_refinement_raw_scipy_failures"]),
        "grid_invalid_evaluations": int(pilot_summary["checks"]["grid_invalid_evaluations_total"])
            + int(extension_summary["checks"]["grid_invalid_evaluations"]),
        "nested_flat_fallback_count": int(pilot_summary["checks"]["nested_fallback_count"])
            + int(extension_summary["checks"]["nested_flat_fallback_count"]),
    }
    endpoint_summary_reconciliation = {
        name: {"record_recount": int(counters[counter_name]), "stored_tranche_summaries": expected,
               "matches": int(counters[counter_name]) == expected}
        for name, counter_name, expected in (
            ("flat_raw_scipy_failures", "flat_raw_failure", endpoint_summary_values["flat_raw_scipy_failures"]),
            ("flat_invalid_penalty_endpoints", "flat_invalid_penalty_endpoints", endpoint_summary_values["flat_invalid_penalty_endpoints"]),
            ("interaction_raw_scipy_failures", "interaction_raw_failure", endpoint_summary_values["interaction_raw_scipy_failures"]),
            ("interaction_invalid_penalty_endpoints", "interaction_invalid_penalty_endpoints", endpoint_summary_values["interaction_invalid_penalty_endpoints"]),
            ("grid_refinement_raw_scipy_failures", "grid_refinement_raw_failure", endpoint_summary_values["grid_refinement_raw_scipy_failures"]),
            ("grid_invalid_evaluations", "grid_invalid_evaluations", endpoint_summary_values["grid_invalid_evaluations"]),
            ("nested_flat_fallback_count", "nested_fallback_total", endpoint_summary_values["nested_flat_fallback_count"]),
        )
    }

    theory4 = load_theory4()
    theory4.LABELS = labels
    sample_audit = []
    scalar_failures = []
    for class_name in ("exceedance", "non_exceedance"):
        chosen = sample_rows[class_name]
        if len(chosen) != 4:
            mismatches[f"insufficient_{class_name}_sample"] += 1
            continue
        for key, row in chosen:
            mock = np.asarray(row["mock_vector"], dtype=np.float64)
            fits = [
                ("flat_LCDM", row["flat_fit"]),
                ("raw_interacting_optimum", row["interaction_raw_search_optimum"]),
                ("selected_interaction_candidate", row["selected_interaction_candidate"]),
            ]
            fit_results = []
            for role, fit in fits:
                try:
                    checked = evaluate_saved_fit(theory4, z, labels, chol, mock, fit, role)
                    fit_results.append(checked)
                except Exception as exc:
                    failure = {"seed": key[0], "index": key[1], "role": role,
                               "error_type": type(exc).__name__, "error": str(exc)}
                    scalar_failures.append(failure)
                    fit_results.append(failure)
            sample_audit.append({
                "seed": key[0], "index": key[1], "class": class_name,
                "selected_delta_chi2": float(row["selection_delta_chi2_flat_minus_interaction"]),
                "nested_flat_fallback_used": bool(row["nested_flat_fallback_used"]),
                "selected_source": row["selected_interaction_candidate"].get("source"),
                "fit_reevaluations": fit_results,
            })

    fit_checks = []
    for sample in sample_audit:
        for reevaluation in sample["fit_reevaluations"]:
            if "independent_chi2" not in reevaluation:
                continue
            fit_checks.append(reevaluation)
    max_obj_difference = max((abs(x["chi2_difference_independent_minus_stored"]) for x in fit_checks), default=math.inf)
    max_pred_difference = max((x["max_abs_prediction_difference"] for x in fit_checks), default=math.inf)

    current_input_paths = [MEAN_PATH, COV_PATH, SEED_MEAN_PATH, DATA_AUDIT2_PATH, OBSERVED_PATH,
                           PILOT_ROWS, PILOT_RESULT, PILOT_RECORD, PILOT_CHECKPOINT,
                           EXTENSION_ROWS, EXTENSION_SUMMARY, EXTENSION_CHECKPOINT,
                           EXTENSION_SCRIPT, TIMING_SCRIPT, ROOT / "experiments/interacting_vacuum_screen/interacting_vacuum_profile.py",
                           THEORY4_SCRIPT]
    input_hashes = {str(path.relative_to(ROOT)): sha256_path(path) for path in current_input_paths}
    input_hashes["work/theory4/interacting_vacuum_identifiability.json"] = sha256_path(ROOT / "work/theory4/interacting_vacuum_identifiability.json")

    checks = {
        "pilot_exact_40_rows_and_expected_indices": pilot_count == 40 and seen_pilot == expected_pilot,
        "extension_exact_1160_rows_and_expected_indices": extension_count == 1160 and seen_extension == expected_extension,
        "pooled_exact_1200_unique_expected_keys_no_overlap": len(seen_pilot | seen_extension) == 1200 and not (seen_pilot & seen_extension) and seen_pilot | seen_extension == expected_all,
        "all_1200_recorded_vectors_exactly_regenerated_and_hash_matched": mismatches["mock_vector_value_or_hash"] == 0 and mismatches["mock_vector_bad_shape"] == 0,
        "all_source_data_and_full_stream_hashes_match_prior_passing_provenance": provenance_hash_ok and all(x["matches_passing_data_audit2_contract"] for x in rng_rebuild.values()),
        "pilot_output_hashes_unchanged_and_protected": pilot_record_ok and pilot_report_ok and pilot_hashes_ok and pprotect_ok,
        "extension_input_hashes_and_checkpoint_integrity_match": extension_hashes_ok and extension_rows_hash_ok,
        "all_records_complete_scores_and_predictions_finite": mismatches["noncomplete_record"] == 0 and not any("nonfinite" in key for key in mismatches),
        "saved_delta_formulas_and_nested_fallbacks_consistent": not any(key in mismatches for key in ("raw_delta_formula", "selection_delta_formula", "nested_fallback_decision", "nested_fallback_payload", "selected_not_raw_best_without_fallback", "raw_candidate_minimum_or_source", "selected_delta_negative_after_fallback")),
        "optimizer_endpoint_status_counts_reconcile": not any("audit_" in key or "missing_attempt_audit" in key or "missing_separate_status" in key for key in mismatches)
            and all(value["matches"] for value in endpoint_summary_reconciliation.values()),
        "eight_deterministic_scalar_quad_rows_recomputed": len(sample_audit) == 8 and not scalar_failures and max_obj_difference < FIT_ABS_TOL and max_pred_difference < PREDICTION_ABS_TOL,
    }
    failures_by_code = {key: int(value) for key, value in mismatches.items() if value}
    status = "pass" if all(checks.values()) else "fail"
    summary = {
        "status": status,
        "scope": "Fitted-flat, finite-search, 13-row DESI DR2 BAO background null calibration; provenance/status audit and eight saved-parameter scalar-distance reevaluations; no full search/refits",
        "environment": {"python": sys.version.split()[0], "numpy": np.__version__,
                        "platform": platform.platform(), "cpu_only": True, "gpu_used": False,
                        "numerical_threads": {name: os.environ.get(name) for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
                        "audit_runtime_seconds": float(time.monotonic() - started)},
        "commands": {
            "audit": "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --seconds 900 -- .venv/bin/python work/data_audit3/audit_interaction_null.py",
            "no_refit_statement": "All 1,200 archived optimizer records were inspected; only the eight selected rows were reevaluated at saved parameter values. No optimizer was called in this audit.",
        },
        "provenance": {
            "input_sha256": input_hashes,
            "prior_data_audit2_status": a2["status"],
            "prior_stream_reconstruction": rng_rebuild,
            "fitted_flat_null_mean_source": "experiments/campaign_seed_bao/result.json, model=lcdm saved prediction",
            "fitted_flat_null_mean_vector_sha256_little_endian_float64": vector_digest(null_mean),
            "covariance_sha256": input_hashes["context/data/desi_dr2_cov.txt"],
            "covariance_shape": list(covariance.shape),
            "covariance_cholesky_positive": True,
            "covariance_min_eigenvalue": float(eig[0]),
            "covariance_condition_number_2": float(eig[-1] / eig[0]),
            "pilot_result_hash_audit": {"record": pilot_record_ok, "result": pilot_report_ok,
                                         "all_result_input_hashes_match_current": pilot_hashes_ok,
                                         "file_hash_checks": pilot_hash_checks,
                                         "recorded_file_sha256": {"pilot_rows": pilot_record.get("rows_sha256"), "pilot_result": pilot_record.get("result_sha256")}},
            "extension_provenance_audit": {"all_summary_input_hashes_match_current": extension_hashes_ok,
                                           "pilot_protection_hashes_match": pprotect_ok,
                                           "pilot_protection": pprotect,
                                           "checkpoint_and_jsonl_hashes_match": extension_rows_hash_ok,
                                           "summary_code_sha256": extension_summary.get("summary_code_sha256"),
                                           "extension_code_sha256_current": sha256_path(EXTENSION_SCRIPT),
                                           "extension_attempt_record_count": extension_summary.get("attempt_record_count"),
                                           "extension_unique_complete_count": extension_summary.get("complete_unique_mock_count")},
            "record_counts": {"pilot_file_rows": pilot_count, "extension_file_rows": extension_count,
                              "pilot_expected_unique_keys": len(expected_pilot), "extension_expected_unique_keys": len(expected_extension),
                              "pooled_unique_keys": len(seen_pilot | seen_extension),
                              "pilot_extension_key_overlap": len(seen_pilot & seen_extension),
                              "complete_status_counts": dict(record_statuses),
                              "missing_or_duplicated_keys": {"pilot": len(expected_pilot.symmetric_difference(seen_pilot)),
                                                             "extension": len(expected_extension.symmetric_difference(seen_extension)),
                                                             "pooled": len(expected_all.symmetric_difference(seen_pilot | seen_extension))}},
            "mock_vector_reconstruction": {"generator": "np.random.default_rng(seed) (PCG64); standard_normal((N,13)); mean + normals[i] @ L.T with L=lower Cholesky(full covariance)",
                                           "stored_vector_field": "mock_vector, float64 values checked array-exactly; little-endian float64 SHA-256 also checked",
                                           "max_abs_value_difference": max_vector_diff,
                                           "vectors_checked": int(pilot_count + extension_count),
                                           "sampled_vector_digests_are_not_substitutes_for_full_stream_digest": True},
        },
        "pooled_null_statistics": {
            "null": "fixed fitted flat-LCDM mean; full correlated Gaussian 13-vector mocks; search rerun under a finite multi-start plus 41x41/refinement interacting-vacuum BAO background profile",
            "statistic": "selection_delta_chi2_flat_minus_selected_interaction_including_exact_nested_flat_fallback",
            "observed_delta_chi2_reference": observed_delta,
            "n": n,
            "exceedance_count_delta_ge_observed": k,
            "exceedance_fraction_descriptive_only": k / n,
            "exact_two_sided_clopper_pearson_95_interval": cp,
            "quantiles_linear_interpolation": {"q50": q50, "q95": q95, "q99": q99},
            "raw_delta_quantiles_linear_interpolation": {"q50": q_raw[0], "q95": q_raw[1], "q99": q_raw[2]},
            "selected_delta_min_median_max": [float(np.min(deltas)), float(np.median(deltas)), float(np.max(deltas))],
            "raw_delta_min_median_max": [float(np.min(raw_deltas)), float(np.median(raw_deltas)), float(np.max(raw_deltas))],
            "classification": "fitted-flat finite-search BAO-background null frequency, not a p-value, posterior probability, evidence, discovery significance, or validation of the full physical interacting cosmology",
        },
        "search_integrity_and_negative_checks": {
            "mismatches": failures_by_code,
            "selected_delta_negative_count": int(np.count_nonzero(deltas < -2e-8)),
            "raw_delta_negative_count": int(np.count_nonzero(raw_deltas < -2e-8)),
            "nested_flat_fallback_count_independent": int(counters["nested_fallback_total"]),
            "stored_summary_nested_fallback_count": int(extension_summary.get("checks", {}).get("nested_flat_fallback_count", -1) + pilot_summary.get("checks", {}).get("nested_fallback_count", -1)),
            "selected_delta_formula_max_abs_error": max_delta_formula_error,
            "raw_delta_formula_max_abs_error": max_raw_delta_formula_error,
            "endpoint_optimizer_audit": {key: int(value) for key, value in counters.items()},
            "endpoint_status_count_reconciliation_vs_tranche_summaries": endpoint_summary_reconciliation,
            "endpoint_audit_mismatch_counts": {key: int(value) for key, value in mismatches.items() if "audit_" in key or "attempt_audit" in key},
            "flat_selected_bound_hits": dict(flat_param_bound_counts),
            "interaction_multistart_selected_bound_hits": dict(interaction_raw_param_bound_counts),
            "final_selected_bound_hits": dict(selected_param_bound_counts),
            "selected_candidate_success_and_source_counts": dict(selected_success),
            "raw_and_selected_gamma_sign_counts": {
                "selected": {"negative": int(counters["selected_gamma_negative"]), "zero": int(counters["selected_gamma_zero"]), "positive": int(counters["selected_gamma_positive"]),
                             "minimum": selected_gamma_min, "maximum": selected_gamma_max},
                "raw_best": {"negative": int(counters["raw_gamma_negative"]), "zero": int(counters["raw_gamma_zero"]), "positive": int(counters["raw_gamma_positive"]),
                             "minimum": float(np.min(raw_g)), "maximum": float(np.max(raw_g))},
            },
            "grid_valid_candidate_source_counts": dict(candidate_sources),
            "known_numerical_behavior": "Raw SciPy success is recorded separately from endpoint physical validity. Invalid penalty endpoints that carry success=true remain invalid and are excluded from profile candidate selection; valid endpoints with success=false remain separately visible. The grid-refinement candidate is retained when its endpoint is physically revalidated, matching the observed-data search contract.",
        },
        "eight_deterministic_theory4_checks": {
            "selection_rule": "sort all complete rows by (seed,index); choose the first four with selected delta >= observed delta and first four below it",
            "independent_code": str(THEORY4_SCRIPT.relative_to(ROOT)),
            "method": "theory4 scalar_history DOP853 plus adaptive scalar quadrature for D_M and full-covariance Cholesky likelihood at saved alpha, Omega_m and Gamma/H0; flat, raw-interaction and selected-interaction candidates evaluated per row",
            "fixed_baryon_split": "try theory4 f_b=0.16; if that illustrative split is not viable, repeat geometry at an explicit positive split below the independently computed feasibility ceiling and retain the f_b=0.16 failure flag",
            "rows": sample_audit,
            "scalar_failures": scalar_failures,
            "maximum_absolute_score_difference": max_obj_difference,
            "maximum_absolute_prediction_difference": max_pred_difference,
            "acceptance_tolerances": {"score": FIT_ABS_TOL, "prediction_component": PREDICTION_ABS_TOL},
        },
        "checks": checks,
    }
    OUTPUT_PATH.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "audit_json": str(OUTPUT_PATH.relative_to(ROOT)),
                      "checks": checks, "n": n, "exceedances": k, "fraction": k/n,
                      "cp95": cp, "q50_q95_q99": [q50, q95, q99],
                      "sample_keys": {name: [[key[0], key[1]] for key, _ in values] for name, values in sample_rows.items()},
                      "scalar_max_objective_diff": max_obj_difference,
                      "scalar_max_prediction_diff": max_pred_difference,
                      "mismatches": failures_by_code}, indent=2))
    return 0 if status == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

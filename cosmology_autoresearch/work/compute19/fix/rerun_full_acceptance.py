#!/usr/bin/env python3
"""Run the preserved producer with a quadrature-resolution override."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
FAILED_PRODUCER = HERE.parent / "ivs_physicality.py"
SIMPSON_NODES = 8193

spec = importlib.util.spec_from_file_location("ivs_failed_producer", FAILED_PRODUCER)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot load preserved producer {FAILED_PRODUCER}")
producer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(producer)
original_integrate = producer.integrate


def integrate_at_resolved_default(g, sample_count=SIMPSON_NODES):
    return original_integrate(g, sample_count)


# Run all existing domain, threshold, g=0, identity, monotonicity and positivity
# checks unchanged, changing only the default dense-output Simpson node count.
producer.HERE = HERE
producer.integrate = integrate_at_resolved_default
exit_status = producer.main()
if exit_status != 0:
    raise SystemExit(exit_status)

result_path = HERE / "physicality_result.json"
result = json.loads(result_path.read_text(encoding="utf-8"))
result["numerics"]["comoving_identity_simpson_sample_count_default"] = SIMPSON_NODES
result["numerics"]["quadrature_diagnosis"] = (
    "The preserved 2049-node attempt exceeded the residual threshold at g=2.5; "
    "convergence at 4097/8193/16385 nodes and adaptive quad showed resolution-only error. "
    "This full rerun uses 8193 Simpson nodes by default; threshold and all other acceptance checks are unchanged."
)
result["numerics"]["preserved_failed_attempt"] = "../ivs_physicality.py and ../physicality_result.json"
result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps({
    "full_acceptance_rerun": "passed all producer checks",
    "simpson_default_nodes": SIMPSON_NODES,
    "output": str(result_path),
}, indent=2))

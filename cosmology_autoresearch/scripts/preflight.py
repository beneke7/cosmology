#!/usr/bin/env python3
"""Record available local resources; never inspect credentials or change settings."""
import argparse
import datetime as dt
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys


def probe(command):
    if not shutil.which(command[0]):
        return {"available": False}
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=10)
        return {"available": True, "returncode": result.returncode,
                "stdout": result.stdout[-10000:], "stderr": result.stderr[-2000:]}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"available": True, "error": str(exc)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="runs/preflight.json")
    parser.add_argument("--start-run-hours", type=float)
    parser.add_argument("--state", default="RUN_STATE.json")
    args = parser.parse_args()
    if args.start_run_hours is not None and args.start_run_hours <= 0:
        parser.error("--start-run-hours must be positive")
    now = dt.datetime.now(dt.timezone.utc)
    versions = {}
    for name in ["numpy", "scipy", "sympy", "torch", "jax", "jaxlib", "cobaya", "camb", "classy", "pymupdf"]:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    memory = None
    try:
        memory = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        pass
    limits = {}
    for rel in ["cpu.max", "memory.max"]:
        path = Path("/sys/fs/cgroup") / rel
        if path.is_file():
            limits[rel] = path.read_text().strip()
    disk = shutil.disk_usage(Path.cwd())
    result = {
        "utc": now.isoformat(), "platform": platform.platform(),
        "python": sys.version, "python_executable": sys.executable,
        "logical_cpu_count": os.cpu_count(),
        "affinity_cpu_count": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        "physical_memory_bytes": memory, "cgroup_limits": limits,
        "free_disk_bytes": disk.free, "package_versions": versions,
        "codex": probe(["codex", "--version"]),
        "gpu": probe(["nvidia-smi", "--query-gpu=name,driver_version,memory.total,memory.free,utilization.gpu", "--format=csv,noheader"]),
        "tools": {name: shutil.which(name) for name in ["git", "pdftotext", "tmux"]},
        "note": "GPU/CPU visibility does not establish an exclusive allocation. No model availability or CUDA kernel execution was tested."
    }
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2) + "\n")
    if args.start_run_hours is not None:
        state = {"status": "started", "start_utc": now.isoformat(),
                 "deadline_utc": (now + dt.timedelta(hours=args.start_run_hours)).isoformat(),
                 "preflight": str(target), "actual_agent_models": {}, "tasks": []}
        state_path = Path(args.state)
        state_path.parent.mkdir(parents=True, exist_ok=True)
        # Exclusive creation protects an existing run's deadline during resume.
        with state_path.open("x") as stream:
            json.dump(state, stream, indent=2)
            stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

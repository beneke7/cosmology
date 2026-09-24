#!/usr/bin/env python3
"""Check whether the DES-Dovekie CPL profile minimum is set by wa's bound."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

from work.inference import fit_dovekie as screen


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "experiments" / "dovekie_screen" / "boundary_sensitivity.json"
WA_INTERVALS = [(-5.0, 3.0), (-10.0, 3.0), (-20.0, 3.0)]
STARTS = 12


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    started = time.perf_counter()
    data = screen.load_inputs()
    base_path = ROOT / "experiments" / "dovekie_screen" / "result.json"
    base = json.loads(base_path.read_text(encoding="utf-8"))
    lcdm_chi2 = next(fit["chi2_profile"] for fit in base["fits"] if fit["model"] == "LCDM")

    fits = []
    for wa_low, wa_high in WA_INTERVALS:
        screen.BOUNDS["CPL"] = [(0.05, 0.6), (-2.0, -0.3), (wa_low, wa_high)]
        fit = screen.fit_model(data, "CPL", STARTS)
        fit.pop("residual_mag")
        fit.pop("mu_base_mag")
        fit["delta_chi2_vs_LCDM"] = float(fit["chi2_profile"] - lcdm_chi2)
        fits.append(fit)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "exploratory_boundary_sensitivity",
        "purpose": "CPL profile minimum in the base comparison is at wa=-3; extend only that lower bound while holding other domains fixed.",
        "baseline_domains": {
            "Omega_m": [0.05, 0.6],
            "w0": [-2.0, -0.3],
            "wa": [-3.0, 3.0],
        },
        "extended_wa_domains": [[-5.0, 3.0], [-10.0, 3.0], [-20.0, 3.0]],
        "parameters_held_fixed": {"Omega_m": [0.05, 0.6], "w0": [-2.0, -0.3]},
        "data_sha256": data["hashes"],
        "main_result_sha256": file_sha256(base_path),
        "main_code_sha256": file_sha256(Path(screen.__file__).resolve()),
        "boundary_check_code_sha256": file_sha256(Path(__file__).resolve()),
        "seed": screen.SEED + list(screen.BOUNDS).index("CPL"),
        "starts_per_model": STARTS,
        "reference_LCDM_chi2_profile": lcdm_chi2,
        "fits": fits,
        "runtime": {
            "wall_seconds": time.perf_counter() - started,
            "command": "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python work/inference/check_cpl_boundary.py",
            "optimizer": "same L-BFGS-B settings and deterministic start recipe as the main profile calculation",
        },
        "interpretation_limit": "An extended-domain minimum at the new lower edge remains boundary-limited; it does not establish an unconstrained CPL best fit or evidence for evolving dark energy.",
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "fits": fits, "runtime": result["runtime"], "output": str(OUTPUT)}, indent=2))


if __name__ == "__main__":
    main()

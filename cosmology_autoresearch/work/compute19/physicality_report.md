# IVS early-time physicality check — stopped on failed validation

Status: **no validated numerical result**. The run was stopped as instructed when the independent comoving-CDM Simpson check failed. No threshold, exclusion, or posterior statement is made.

The source was the exact local arXiv v1 PDF at `context/papers/interacting_de_desi_dr2_2026.pdf` (SHA-256 below). Direct text extraction confirmed Table I's IVS uniform `g=Gamma/H0` interval `[-3,3]` and Eqs. (3)–(8): flat Friedmann density sum; `dot(rho_x)+3H(1+w_x)rho_x=Q`; `dot(rho_c)+3H rho_c=-Q`; `Q=Gamma rho_x`; hence positive Gamma is forward-time CDM→vacuum transfer. The standalone script derives the backward-u system and comoving identity in its header/result schema. This is not a mapping of the paper's sampled physical-density and `100 theta_MC` priors to fixed `Omega_m0`.

Fiducial requested: `Omega_m0=0.315`, `Omega_b0=0.049`, `H0=67.4`, `Tcmb=2.7255 K`, `Neff=3.044`, massless neutrinos, flat closure, integrate toward `z=1e10`. Implementation is an independent one-CPU DOP853 background script; it does not import compute15/17 code and does not evaluate a likelihood.

The first bounded command used system `python3` and failed immediately because that interpreter lacked SciPy (`ModuleNotFoundError`); no integration was started. After confirming `.venv/bin/python` was available, the bounded DOP853 run reached validation and stopped with `ArithmeticError: independent Simpson identity failed for g=2.5` at the configured `2e-8` tolerance. The script writes its result JSON only after every numerical check passes, so no result JSON was produced. I did not loosen the tolerance or rerun after this validation failure.

Both commands used `scripts/run_bounded.py --state work/compute19/RUN_STATE.json --seconds 600`, where the task-local state file was absent, and environment limits `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`. The second command was:

```text
python3 scripts/run_bounded.py --state work/compute19/RUN_STATE.json --seconds 600 -- env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python work/compute19/ivs_physicality.py
```

The scan/threshold diagnostics computed before the failed check are deliberately not reported as findings. Required next step is to resolve or independently strengthen the quadrature validation before any physical threshold can be promoted. Any future threshold would be benchmark-specific only and would not imply anything about the published analysis's posterior.

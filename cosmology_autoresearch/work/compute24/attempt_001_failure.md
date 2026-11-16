# Attempt 001 — CAMB physical-parameter profile

Status: **stopped before any fit result**. The first differential-evolution population included a trial for which CAMB's `set_H0_for_theta` could not bracket an H0 solution in the explicitly recorded [10,100] km/s/Mpc interval. The resulting `CAMBParamRangeError` escaped the objective and stopped the search. No score or parameter estimate from this attempt is retained as a finding.

- Command: `env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=work/compute24/site OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --state work/compute24/no_campaign_deadline.json --seconds 600 -- .venv/bin/python work/compute24/camb_bao_lcdm.py --workers 20 --maxiter 60 --popsize 12 --seeds 20260925 20260926`
- Source-script SHA-256 before correction: `5216d15c014078c9cf00eede5c3cf5e6e48187652cfa4bce5829b704c290be11`.
- Failure class: `camb.baseconfig.CAMBParamRangeError` (subclass of `CAMBError`), raised for at least one point in the broad Table-I `100theta_MC`, `omega_b h^2`, `omega_c h^2` box.
- Correction: map only documented CAMB/parameter-domain errors to the predeclared invalid-point penalty. Keep the Table-I box and the H0 bracketing interval unchanged; preserve invalid support as a reported limitation rather than altering bounds to help the fit.
- No output JSON or plot was produced by attempt 001.

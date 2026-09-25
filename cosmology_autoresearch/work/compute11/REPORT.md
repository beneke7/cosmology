# Dovekie IVS finite-search null calibration

Status: `complete_200`. This is an empirical finite-search calibration under the fitted Dovekie LCDM shape; it is not a posterior p-value or discovery significance.

Generated 200 of 200 target realizations with master seed `20260925`. Each mock used the 1,820-row released STAT+SYS precision and one common fitted magnitude intercept. The realized noise was generated as `solve(L.T, z)` for `P=L L.T`, which gives covariance `C=P^-1`.

Each valid mock refit LCDM with 8 starts and the same IVS finite search: Ωm ∈ [0.05, 0.60], g ∈ [-3, 3], a 41×41 grid, 16 L-BFGS-B starts including the grid minimum, the positive-history guard, invalid penalty, and endpoint rechecks. One 20-process pool ran realizations; searches within each realization ran serially, with BLAS threads set to one.

The count with `T = χ²_LCDM − χ²_IVS ≥ 1.421318647174758` was 50/200 = 0.25. Exact two-sided 95% Clopper–Pearson interval: [0.191607, 0.315963]. The denominator includes completed mock fits; failed fits are counted and reported separately.

The 10-mock pilot took 10.898 s wall time. Based on its measured per-mock runtime and 20 workers, the projected total with a 25% throughput allowance and 30 s finalization reserve was 158.0 s, against a 1440 s internal cap. Projection decision: `continue_to_target`.

Failed realization fits: 0; invalid IVS grid points summed over complete mocks: 91600; invalid IVS optimizer objective evaluations: 7293; IVS invalid endpoints: 1400; LCDM invalid endpoints: 0; raw unsuccessful starts (LCDM/IVS): 9/2.

Pinned input and code hashes, algebra check, full run command, per-realization statuses/statistics, and environment are recorded in [result.json](result.json) and [mock_results.jsonl](mock_results.jsonl). SHA-256 checksums for the saved artifacts are in [SHA256SUMS](SHA256SUMS). The generator, scoring and search are in [dovekie_ivs_null_calibration.py](dovekie_ivs_null_calibration.py).

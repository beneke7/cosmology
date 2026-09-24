# Independent audit: interacting-vacuum fitted-flat null

**Disposition: PASS.** I independently verified the pooled 1,200 archived mock records, their complete PCG64 vector streams, per-record selection arithmetic and saved status tallies. I also reevaluated eight predeclared rows at their saved parameters with the separate adaptive scalar-distance implementation. No mock was refit; recorded audit runtime was 0.48 s, CPU-only, with one numerical thread.

## Pooled result

The fixed null is the saved fitted flat-LCDM mean, with full-covariance Gaussian 13-vectors. Each mock reruns the finite interacting-vacuum BAO-background search (bounded multistart plus 41×41 grid/refinement); the reported statistic is

`Δχ² = χ²_flat − χ²_selected`,

where selection includes the exact nested-flat fallback. Against the observed `Δχ² = 1.5171630696615157`, **256/1,200** mocks meet or exceed the threshold: a descriptive frequency of **0.21333**, with exact two-sided Clopper–Pearson 95% interval **[0.19046, 0.23762]**. Linear-interpolation quantiles are q50 = **0.47047**, q95 = **3.78610**, q99 = **6.15132**. The selected statistic ranges from `2.3343e-7` to `7.48847`; no selected negative deltas were found.

This is a finite-search, fitted-flat, DESI DR2 BAO-background null frequency—not a p-value, posterior probability, evidence, discovery significance, or validation of the full physical interacting cosmology. A cross-model selection claim still requires the complete selection rule and its calibration.

## Provenance and vector reconstruction

The exact row keys are complete and disjoint: seed **20260924**, indices 0–199 (pilot 0–19 plus extension 20–199); seed **20260925**, indices 0–999 (pilot 0–19 plus extension 20–999). Both separately initialized `np.random.default_rng(seed)` streams used **PCG64** and `standard_normal((N,13))`; each full little-endian float64 normal stream matches the prior passing `data_audit2` digest. With the saved fitted-flat prediction `μ`, and `L = cholesky(C)` for the full 13×13 covariance, zero-based draw `i` is reproduced as `μ + normals[i] @ L.T`. All **1,200/1,200** saved vectors are shape-checked, array-exact against regenerated values (maximum absolute difference **0**), and their saved vector hashes match.

| Input | SHA-256 |
|---|---|
| DESI mean source `context/data/desi_dr2_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| DESI full covariance `context/data/desi_dr2_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| Saved fitted-flat seed result | `92662f3ee42500d2be34d9f87b485fe56c82e9550fa674e42f4dec2bf2cab1f0` |
| Prior passing `data_audit2/provenance_check.json` | `8c4fed3713b9670c37133bfef13ca1a315acbc7d7880e23547f35673fdd7f581` |
| Observed interacting-vacuum result | `e347b8934b2afdc23d23ef7d39c70f963b73888d17137f55caf7e51af87d43ae` |
| Pilot 40 JSONL | `a5138b807ff21baf776e389848177db617e08b6f59d9318b3ff1e421f9dab09e` |
| Pilot 40 result JSON | `402ed8fe5c9864be6629aa0022f7e666d3a919674b8f75b3333951363cb27b97` |
| Extension 1,160 JSONL | `caf01a78e711ede965f98b5e389e28affe74e1e292e888c1483fc8fcc85feed1` |
| Extension summary JSON | `8bcc93486846237f020b3178403518a8d7f668da5761c0bf7d4c4e83ec8c287a` |
| Extension checkpoint JSON | `55b616637a08c0cf47e2559969dc15beedb34359e4f81f50ce8b44dd3e9d2a1d` |
| Independent scalar implementation `work/theory4/audit_interacting_identifiability.py` | `6cb55d20d50baa4520970296fdb74bf94f66273c1db3cf0bf2438018194a3052` |

The saved fitted-flat mean vector's little-endian float64 digest is `707f6894421b38233ed79459c814cdf23a76df29e28973ac5dbebdf8851096b0`. Full normal-stream digests are `a5265eba65345219facd2ce66e1edf653f7a2a872f51a79905db228fe1f3de21` (seed 20260924, N=200) and `b9cf266dee55dd1f67366c030273a27ae1fb107691baf8bfc3eab07fa407c985` (seed 20260925, N=1,000). The seeds were independently initialized, and the audited key partitions do not overlap; within each stream pilot and extension indices partition the complete zero-based range.

The pilot record, checkpoint, result inputs and extension's pilot-protection hashes all match the current files; the extension reports no pilot output paths modified. The covariance Cholesky succeeds; minimum eigenvalue is **0.00578999** and 2-norm condition number **116.824**. The JSON artifact contains the full file-hash manifest, both stream digests and all row-level audit data.

## Search/status checks and negative results

Every record is marked complete, has finite saved scores and predictions, and passes both delta identities (`flat − raw` and `flat − selected`); the largest formula discrepancy is **0**. There were **zero** nested-flat fallbacks, zero negative raw/selected deltas, and no final selected parameter-bound hits. The selected `Γ/H₀` signs are nearly balanced: 597 negative / 603 positive (range **−1.09104 to +1.09452**), with no selected boundary hit in the declared primary `[-3,3]` domain.

Status and physical validity were counted separately and reconcile against the archived tranche summaries. Across 1,200 records: flat search had 19,200 starts, 742 raw SciPy failures but 19,200 physically valid finite endpoints; interaction multistart had 46,800 starts, 88 raw failures, 32,400 valid endpoints, and 14,400 invalid penalty endpoints. Notably, all 14,400 invalid interaction endpoints carry raw SciPy `success=true`, so success alone is not a validity test. Grid refinement had 80 raw failures but all 1,200 endpoints were physically revalidated; it encountered 584,400 invalid grid evaluations. Those are search diagnostics, not mismatches in the saved selected candidates.

The eight deterministic rows were the first four exceedances and first four non-exceedances in sorted `(seed,index)` order: exceedances `(20260924,3)`, `(20260924,16)`, `(20260924,29)`, `(20260924,31)`; non-exceedances `(20260924,0)`, `(20260924,1)`, `(20260924,2)`, `(20260924,4)`. For each row I reevaluated flat, raw-interaction and selected-interaction fits (24 fixed-parameter evaluations total) using theory4's scalar ODE plus adaptive quadrature and the full-covariance Cholesky likelihood. Maximum absolute objective discrepancy was **5.24e-13** and maximum componentwise prediction discrepancy **4.27e-14**. No optimizer was called; all eight passed the `f_b=0.16` split check.

## Reproduction

Run from the repository root, CPU-only and with one numerical thread:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --seconds 900 -- .venv/bin/python work/data_audit3/audit_interaction_null.py
```

It completed in 0.48 s of recorded audit time with all 11 checks passing. The checker SHA-256 is `3a897a9c039198d0ea15db7512a8df4e151dbebb5833dec02c7b302f18b3859c`; the machine-readable result is `interaction_null_audit.json` (SHA-256 `418bdad5321db200ee54bde121f7f8ef42c117b4a5dc9caf60ac502c1534723d`).

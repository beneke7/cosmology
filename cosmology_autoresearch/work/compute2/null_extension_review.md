# Audit of the DESI BAO null-search extension

Status: independently audited from the saved per-realization search statistics. I did not rerun the mock-fitting ensemble. The files contain 200 realizations from seed 20260924 and 1,000 from seed 20260925.

## Reproduction and input identity

The original and extension records use the same search contract: ΛCDM, constant-w, and CPL with CPL `w0` upper bounds `{-0.3, -0.1, 0, 0.1}`; 16-start primary and 64-start expanded searches; the other bounds and per-fit optimizer seed `20260923` are identical. Both select `cpl_w0_upper_0.1` for the observed extension. Recomputing `χ²_LCDM − min(χ²_extension)` from the stored observed candidate fits reproduces the recorded statistic exactly: 4.652440971036443 for 16 starts and 4.652440971083040 for 64 starts. The two output files have identical source data and code hashes, and each recorded hash matches the live file.

The generator in `scripts/run_bao_robustness.py` draws `standard_normal((n_mocks, 13))` from `np.random.default_rng(seed)` and transforms with the released covariance Cholesky factor. The two recorded seeds instantiate PCG64 streams with distinct reproducible draws. I regenerated the normal arrays only: identical seeds reproduce their arrays bit for bit, while the two seed streams differ. The mock data vectors themselves are not stored, so this check validates the deterministic stream definitions but does not independently refit each realization from its mock vector.

## Tail counts and pooled calibration

Counts below come directly from each saved realization's `search_statistic_delta_chi2` (16 starts) or `expanded_search_statistic_delta_chi2` (64 starts), using `statistic >= observed`. Intervals are two-sided 95% Clopper–Pearson intervals. The pooled sample concatenates the two files' disjoint seed runs; `q95` uses NumPy's linear quantile convention.

| Search | Observed Δχ² | Extension, seed 20260925 | Pooled seeds 20260924+20260925 | Pooled q95 |
|---|---:|---:|---:|---:|
| 16 starts | 4.652440971036443 | 66/1000 = 0.066; 95% CI [0.051409, 0.083206] | 83/1200 = 0.069167; 95% CI [0.055464, 0.085024] | 5.169980759925 |
| 64 starts | 4.652440971083040 | 66/1000 = 0.066; 95% CI [0.051409, 0.083206] | 83/1200 = 0.069167; 95% CI [0.055464, 0.085024] | 5.169980759959 |

For context, the original 200-mock run has 17 exceedances for each search budget (17/200 = 0.085; 95% CI [0.050296, 0.132605]). The pooled q95 exceeds the observed statistic for both budgets. These are Monte Carlo tail frequencies and finite-count intervals under the fitted ΛCDM null and this finite profile-search procedure. They are not posterior p-values, model evidence, or discovery significance.

## Abnormal optimizer candidates

The extension marks 282 of 12,000 candidate fits as nonconverged across the 16- and 64-start searches. All 282 are present in the diagnostic recheck file and marked abnormal in the original records. Across those records, directly recomputed original objective values differ from reported values by at most `2.3448×10⁻¹³`. The 846 alternative optimizer runs all succeeded. Best checked χ² minus original direct χ² ranges from `−4.1795×10⁻¹⁰` to zero, and none improves by more than `1×10⁻⁷`. Rechecking therefore does not change the saved exceedance counts at the stated precision.

## Files and hashes

The machine-readable audit is [null_extension_review.json](null_extension_review.json). Input and code hashes checked before interpreting results:

- `null_extension_1000_seed20260925.json`: `5063ae50f267477b41febd30eabe00c3167b0ccb2d826cd94f3fc966dc025693`
- `result.json`: `0f477285029024fa359740e726b46f2abf8c87ea8704f488d8ad6658c3fc2a3e`
- `null_extension_1000_optimizer_recheck.json`: `69829e76be7139edceb7ee7e9699bcc058fa2f8b8ed75bf056c7bc15fd305c5b`
- Mean data: `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`
- Covariance: `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`
- `run_bao_robustness.py`: `54cb6d78b45f0b9a495bdf8733459ccef6b6a9b2cfc6cc5e8c835a602a09ebcb`
- `background_bao.py`: `34bdf1d81397b39b740faa3e2f4a79748d5fa9d3bfdd8303619d631dbf8aa2e2`

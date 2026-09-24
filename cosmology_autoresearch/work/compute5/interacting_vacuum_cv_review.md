# Independent review: seven-fold IVS vs flat-LCDM block CV

**Disposition: PASS.** Dense-solve reconstruction, independent scalar-distance predictions, row/covariance order, and all published input-hash claims check out. I ran no optimizer and did not refit any parameter. The audit took 0.14 s, CPU-only, with one numerical thread.

| Held-out z | Rows (0-based) | Flat conditional χ² | IVS conditional χ² | Δ(IVS − flat) |
|---:|---:|---:|---:|---:|
| 0.295 | 0 | 0.503274 | 0.021391 | −0.481882 |
| 0.510 | 1–2 | 5.025057 | 4.190968 | −0.834089 |
| 0.706 | 3–4 | 5.137111 | 4.225458 | −0.911653 |
| 0.934 | 5–6 | 0.974987 | 3.652959 | +2.677971 |
| 1.321 | 7–8 | 0.474356 | 0.622664 | +0.148308 |
| 1.484 | 9–10 | 0.698714 | 0.780546 | +0.081832 |
| 2.330 | 11–12 | 0.153520 | 22.223332 | +22.069812 |
| Sum of seven folds | 13 held-out appearances | 12.967018 | 35.717318 | **+22.750300** |

My dense solve gives Δ = **+22.75029962434854**, within `8.1e-13` of the result artifact. The record file rounds this to `22.750300` (difference from full precision `3.76e-7`, exactly consistent with six-decimal rounding). At z=2.33 the independently reconstructed contribution is **+22.069812381804915**, within `9.3e-13` of the saved value. This single fold contributes about 97% of the aggregate difference; the aggregate is therefore a localized high-redshift stress result, not a broad uniform separation.

## Independent numerical checks

I read the original mean file in non-comment row order and recomputed all seven partitions. Held-out rows are exactly `[0]`, `[1,2]`, `[3,4]`, `[5,6]`, `[7,8]`, `[9,10]`, `[11,12]`; every training set is the exact complement. Row metadata (z, observable, value) also matches. The z=2.33 pair is retained in source order as `DH/rs`, then `DM/rs`.

For each fold I formed `C_TT`, `C_HT`, `C_TH`, and `C_HH` from the released full covariance and used general dense `numpy.linalg.solve` calls to calculate `μ_H + C_HT solve(C_TT, y_T−μ_T)` and `S = C_HH − C_HT solve(C_TT,C_TH)`. I then scored with another dense solve through `S`. The maximum difference versus saved conditional means was **2.13e-14**, versus saved Schur covariances **0**, and versus saved conditional χ² **9.42e-13**. All seven Schur matrices and the full 13×13 covariance are positive definite; the minimum eigenvalue is **0.00578998687** (full covariance condition number **116.824**).

The original result stores training predictions and conditional held-out means/covariances, but not the raw held-out model-only predictions. I independently generated all **14 model predictions** (two models × seven saved fold fits) at the saved `alpha`, `Omega_m`, and `Gamma/H0`: the separate `work/theory4/audit_interacting_identifiability.py` scalar ODE for `(M,X)` plus adaptive QUAD for distances, with flat LCDM represented by `g=0`. No refitting occurred. Maximum component difference from stored training predictions is **2.84e-14**; using regenerated held-out predictions with the saved training predictions reproduces the stored conditional means above. The IVS histories have positive `M`, vacuum and `E²`, with a positive baryon/CDM split witness throughout the 0–2.33 interval.

The aggregate totals independently reconstruct as flat **12.96701832203225**, IVS **35.71731794638079**, difference **22.75029962434854**. Both reported fold-level delta arithmetic and the record/result aggregate agree at the reported precision.

## CV semantics and interpretation

Code inspection confirms that each fold sends only the complement rows to `fit_lcdm(train)` and `fit_ivs(train)`. Amplitude profiling, optimizer/grid/refinement candidate selection, and the training prediction all use training data/covariance only. Held-out values enter only after parameters are fixed, in the conditional residual. The search's z=2.33 physical-domain ceiling is a model-domain check, not use of held-out measurements.

No arithmetic or leakage defect was found. A few interpretation limits matter:

- These are plug-in conditional chi-squares, not posterior-predictive CV; parameter uncertainty is not integrated.
- Omitting the conditional log determinant is harmless for the *within-fold model difference*: `S` depends on the shared data covariance and fold, not on which model is scored. It would not be an absolute predictive log score. Summing seven fold differences is descriptive only because the folds have overlapping training sets and share one survey/covariance; it is not a joint likelihood, evidence, or independent seven-fold sample.
- The aggregate is dominated by the two-row z=2.33 holdout. The result supports “this saved IVS fit predicts that held-out high-z pair poorly under this contract,” not universal IVS failure or a physical exclusion.
- This is compressed late-time BAO background only; no sound-horizon calibration, CMB, perturbation closure/stability, or growth test is performed.

## Provenance and reproduction

All report, record, and result SHA claims were checked against current files. Current SHA-256 values:

| File | SHA-256 |
|---|---|
| `experiments/interacting_vacuum_block_cv/record.json` | `8ad8df6c4f2187fe11d386f4d88198d42442e9f8cf5bc3fe5bf657aff8478907` |
| `experiments/interacting_vacuum_block_cv/result.json` | `9d1bf8ac6536328d36d5e26d3a1dd2e2a8ddf6c2761ad6472a1ac9adc9b69a55` |
| `experiments/interacting_vacuum_block_cv/leave_bin_out.py` | `81c5ce44ab6ec27539290a9d3754e85d5dd809559dc04d8558a4a6d31cb3367c` |
| `context/data/desi_dr2_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `context/data/desi_dr2_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| `scripts/background_bao.py` | `34bdf1d81397b39b740faa3e2f4a79748d5fa9d3bfdd8303619d631dbf8aa2e2` |
| `experiments/interacting_vacuum_screen/interacting_vacuum_profile.py` | `ae99ae789aa02885667e27aa2cbbaa0b035323e69d1a4bf7bbd3848949fc5282` |
| Published report `work/compute4/interacting_vacuum_block_cv.md` | `2bea08327f663287b918e459a9e228463103d4e0ecfc4a80ac13d182cc2c94df` |
| Independent theory4 scalar code | `6cb55d20d50baa4520970296fdb74bf94f66273c1db3cf0bf2438018194a3052` |

Run from repository root with one numerical thread; `run_bounded.py` enforces the six-minute ceiling:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --seconds 360 -- .venv/bin/python work/compute5/interacting_vacuum_cv_review.py
```

All seven audit checks pass. Machine-readable fold reconstructions and hash checks are in `interacting_vacuum_cv_review.json`; checker SHA-256 is `f118995d103cd1a27763224c0e75928208e1c1f19867b862f5e0f1ebb108ea88`, result JSON SHA-256 is `e448a6e53c4c9a34a8b4088b30c4ef22e44a4bf1a3147d21e1a6f23c954bbd86`.

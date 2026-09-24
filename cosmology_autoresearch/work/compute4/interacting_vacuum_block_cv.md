# Leave-one-redshift-bin-out BAO predictive check

## Result

For the seven leave-one-redshift-group-out folds, the sum of conditional predictive chi-square values was 12.9670 for flat LCDM and 35.7173 for the late-time IVS screen. With Δ defined as IVS minus LCDM, the aggregate is +22.7503, so this plug-in predictive score favors LCDM overall. IVS scores better in the three lowest-redshift folds, but its held-out score is worse by 22.0698 at z=2.33. Fold scores are:

| Held-out z | Rows | LCDM conditional χ² | IVS conditional χ² | Δ(IVS − LCDM) |
|---:|---:|---:|---:|---:|
| 0.295 | 1 | 0.503274 | 0.021391 | -0.481882 |
| 0.510 | 2 | 5.025057 | 4.190968 | -0.834089 |
| 0.706 | 2 | 5.137111 | 4.225458 | -0.911653 |
| 0.934 | 2 | 0.974987 | 3.652959 | +2.677971 |
| 1.321 | 2 | 0.474356 | 0.622664 | +0.148308 |
| 1.484 | 2 | 0.698714 | 0.780546 | +0.081832 |
| 2.330 | 2 | 0.153520 | 22.223332 | +22.069812 |
| Sum across folds | 13 held-out appearances | 12.967018 | 35.717318 | +22.750300 |

The sum is a convenient summary of seven overlapping conditional assessments; it is not the chi-square of one joint prediction. These are plug-in fits, not posterior predictive cross-validation or an evidence calculation. There are only seven redshift groups, and the folds come from one DESI DR2 release with its shared covariance. This is not an independent survey, a physical consistency test, or a CMB test.

## Rows and covariance

The input order is the 13 non-comment rows of `context/data/desi_dr2_mean.txt`; row indices below are zero-based. Each distinct z is held out as one group, with all remaining rows used for training:

| z | Input rows | Measurements |
|---:|---:|---|
| 0.295 | 0 | DV/rs |
| 0.510 | 1–2 | DM/rs, DH/rs |
| 0.706 | 3–4 | DM/rs, DH/rs |
| 0.934 | 5–6 | DM/rs, DH/rs |
| 1.321 | 7–8 | DM/rs, DH/rs |
| 1.484 | 9–10 | DM/rs, DH/rs |
| 2.330 | 11–12 | DH/rs, DM/rs |

For training indices T and held-out indices H, the calculation used the exact released covariance blocks in this same row order. With model mean μ fitted only on T,

\[
\mu_{H|T}=\mu_H+C_{HT}C_{TT}^{-1}(y_T-\mu_T),\qquad
S_{H|T}=C_{HH}-C_{HT}C_{TT}^{-1}C_{TH},
\]

and scored \((y_H-\mu_{H|T})^T S_{H|T}^{-1}(y_H-\mu_{H|T})\). The inverses are implemented as Cholesky solves, not diagonal or block-diagonal approximations. Each model fits alpha and Omega_m to the training rows; IVS additionally fits g=Gamma/H0. Alpha is analytically profiled using the training covariance. Search boxes are alpha [1e-6, 1e4], Omega_m [0.05, 0.60], and g [-3, 3]; they are numerical search domains, not priors.

The IVS BAO likelihood uses the existing late-time flat background implementation. Every optimizer and grid candidate was checked for positive E^2, matter, vacuum, and a feasible positive baryon/CDM split densely across 0 <= z <= 2.33, even on folds omitting z=2.33.

## Numerical and optimizer checks

The conditional means and Schur covariances match independently recomputed dense linear solves exactly at stored float precision. The largest difference in held-out chi-square versus a dense solve was 3.6e-15. All seven conditional covariance matrices were symmetric at stored precision and positive definite; the smallest reported eigenvalue was 0.00579.

The IVS used eight deterministic L-BFGS-B starts and a 7x7 Omega_m/g grid per fold, with Powell refinement from each fold's three best finite grid points. All 56 multistart status records and all 21 refinement records are retained in `experiments/interacting_vacuum_block_cv/result.json`. The multistart returned 56 SciPy successes (status 0); 35 endpoints were physically finite, while invalid endpoints and their domain failures are retained. Grid refinement returned 21 successes (status 0). Across multistart endpoints, parameter-bound hits were g lower: 21, Omega_m lower: 20, and g upper: 14; there were no alpha endpoint hits. The selected best fit in each fold was finite, status 0, and had no active parameter bound. The grid/refinement selected training scores agreed with or improved upon the multistart results in every fold.

For LCDM, 88 of 91 starts converged with status 0; three returned status 2 (`ABNORMAL`). The failed candidates' endpoint scores and all messages/statuses are retained; each selected fold result came from a successful, finite endpoint. No selected LCDM endpoint hit the alpha or Omega_m search boundary. The complete per-start optimizer records, boundary hits, selected fold parameters, and physical-domain minima are in the machine-readable result.

## Reproduction and provenance

Run from `cosmology_autoresearch/` with a single numerical thread:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 timeout 90s .venv/bin/python experiments/interacting_vacuum_block_cv/leave_bin_out.py
```

The final run completed all seven folds in 34.6 seconds (wall time approximately 35 seconds). The run used no GPU and downloaded no data. A first, incomplete timing attempt with a larger grid/start count was interrupted after 4.5 minutes; it produced no result artifact. The final run completed after reducing those counts while retaining deterministic multistart and grid/local cross-checks. Stop disposition: requested folds and checks complete; compute stopped normally within the 15-minute ceiling.

SHA-256:

| File | SHA-256 |
|---|---|
| `experiments/interacting_vacuum_block_cv/leave_bin_out.py` | `81c5ce44ab6ec27539290a9d3754e85d5dd809559dc04d8558a4a6d31cb3367c` |
| `experiments/interacting_vacuum_block_cv/result.json` | `9d1bf8ac6536328d36d5e26d3a1dd2e2a8ddf6c2761ad6472a1ac9adc9b69a55` |
| `context/data/desi_dr2_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `context/data/desi_dr2_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| `scripts/background_bao.py` | `34bdf1d81397b39b740faa3e2f4a79748d5fa9d3bfdd8303619d631dbf8aa2e2` |
| `experiments/interacting_vacuum_screen/interacting_vacuum_profile.py` | `ae99ae789aa02885667e27aa2cbbaa0b035323e69d1a4bf7bbd3848949fc5282` |

The machine-readable result records software/environment information and input hashes. No shared report, ledger, or `RUN_STATE` file was edited.

# Independent Dovekie frozen-shape prediction review

**Independent result:** using dense covariance-domain conditioning on the pinned SN rows, the IVS-minus-flat-ΛCDM descriptive sum of the four primary normalized conditional scores is **−4.0365877404**. The training-offset plug-in sensitivity is **−5.1975327992**. All six predeclared edge variants retain a negative aggregate IVS-minus-ΛCDM difference. These overlapping conditional sums are descriptive, not one additive likelihood or significance measure.

## Independent reconstruction

I parsed the SNANA `VARNAMES`/`SN` records directly in release order, unpacked the archived upper-triangle float32 precision, Cholesky-solved the dense full covariance `C = P^-1`, and performed conditioning from covariance submatrices. I did not import `conditional_cv.py` or any of its helper modules. My independent prediction integrated the continuity equations in `u = ln(1+z)`: `dM/du = 3M + gX/E`, `dX/du = -gX/E`, and `E^2 = M+X`. I then used adaptive scalar quadrature for `DC(zHD) = integral from 0 to zHD of dz/E(z)` and `DL = (1+zHEL)(c/H0)DC`. The exact BAO-result points were ΛCDM `(Omega_m,g) = (0.2974618149542753, 0)` and IVS `(0.38697100769747017, -0.4666647299169934)`. No shape fit or refit was done.

For each fold, the covariance-domain calculation used `Mhat_T = (1_T^T C_TT^-1 (y_T-g_T))/(1_T^T C_TT^-1 1_T)`, `S = C_HH-C_HT C_TT^-1 C_TH`, and `a = 1_H-C_HT C_TT^-1 1_T`. The conditional mean was `g_H+Mhat_T 1_H+C_HT C_TT^-1 (y_T-g_T-Mhat_T 1_T)`. Primary covariance was `S + a a^T/(1_T^T C_TT^-1 1_T)`; the plug-in sensitivity used `S` alone.

## Primary four-fold result

The locked edges are `zHD = (0.32433, 0.479345, 0.626935)`; each fold holds out 455 SNe and trains on 1,365. Per-fold differences are IVS minus flat ΛCDM, so negative values favor IVS on that conditional score.

| Held-out fold | zHD range | Integrated Δ(−2 log predictive) | Plug-in Δ(−2 log predictive) |
|---:|---:|---:|---:|
| 0 | 0.02509–0.32394 | −5.5418304001 | −6.7825833325 |
| 1 | 0.32472–0.47905 | +0.7875784746 | +0.7668307624 |
| 2 | 0.47964–0.62684 | +0.3389940653 | +0.2926124516 |
| 3 | 0.62703–1.14418 | +0.3786701199 | +0.5256073193 |
| Descriptive sum | — | **−4.0365877404** | **−5.1975327992** |

Summed integrated scores were −513.2567841977 (ΛCDM) and −517.2933719381 (IVS); summed plug-in scores were −512.6967574233 and −517.8942902225. Fold 0 supplies the main IVS improvement; the other three folds favor ΛCDM modestly. The four conditional training sets overlap, so neither sum is an additive joint predictive density.

## Predeclared edge variants

Each internal quartile target was moved independently by ±1% of `N`, with the other two targets fixed; equal-z groups remained intact. I recomputed all six maps and their fold evaluations. Their edges and held-out row maps match the recorded run exactly.

| Variant | Edges zHD | Aggregate integrated Δ | Aggregate plug-in Δ |
|---|---|---:|---:|
| boundary 1 left | 0.31484, 0.479345, 0.626935 | −3.0485688103 | −3.5226631279 |
| boundary 1 right | 0.33044, 0.479345, 0.626935 | −3.9865401339 | −5.1962456871 |
| boundary 2 left | 0.32433, 0.473285, 0.626935 | −3.9300529621 | −5.0385096222 |
| boundary 2 right | 0.32433, 0.48994, 0.626935 | −4.0890347025 | −5.2883397599 |
| boundary 3 left | 0.32433, 0.479345, 0.621895 | −3.8525614719 | −4.9392660704 |
| boundary 3 right | 0.32433, 0.479345, 0.63399 | −4.1179544579 | −5.3665242491 |

Across the variants, integrated Δ ranges from **−4.1179544579 to −3.0485688103** and plug-in Δ from **−5.3665242491 to −3.5226631280**. The aggregate sign is stable over this perturbation set. This is a descriptive sensitivity of the fold scheme, not a model-selection correction.

## Cross-checks against the recorded run

Input SHA-256 values were `2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf` (Hubble diagram) and `ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b` (precision archive). Original-order CID SHA-256 was `b7d5c6ad8dfdadf1e12443852b14006c0d9d0ceb4bcb46efb378f058ce5aa5a3`. The full covariance was symmetric to `9.71e-17` before roundoff symmetrization; Cholesky succeeded with minimum diagonal `0.06750`; max absolute `P*C-I` was `7.52e-15`. The independent Gamma=0 distance check agreed with flat ΛCDM at the same Omega_m to `9.64e-16` mag; adaptive-quadrature estimated relative error was at most `1.12e-14`.

I compared against `experiments/dovekie_frozen_shape_cv/result.json` only after completing the independent calculation. The baseline and all six variant fold maps/edges match exactly. Maximum absolute scalar discrepancies across model-by-fold offset, offset variance, integrated/plugin quadratic, log determinant and normalized score were:

- offset: `4.72e-16` mag;
- offset variance: `3.39e-21` mag²;
- integrated quadratic: `7.96e-13`; integrated log determinant: `6.82e-13`; integrated normalized score: `6.82e-13`;
- plug-in quadratic: `8.53e-13`; plug-in log determinant: `4.55e-13`; plug-in normalized score: `9.09e-13`.

Maximum aggregate primary/plugin Δ discrepancies were `1.25e-12` and `1.82e-12`. Across the six variants, the largest Δ discrepancy was `2.84e-12` integrated and `3.07e-12` plug-in. The recorded result JSON SHA-256 was `c28468a4ab9b20a71822eb2c7d971bb78a3e780e587992b2dc6c76ea0c8cd64f`. This is an independent numerical reproduction to floating-point agreement, not a reproduction of the paper's sampler or full likelihood pipeline.

## Reproducibility and limits

The calculation completed in 5.785 seconds with Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, and one numerical thread. Command prefix: `timeout 590s env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /home/v/proj/bene/cosmology/cosmology_autoresearch/.venv/bin/python -`. The independent implementation was passed inline and was not saved as a separate script because this assignment owns only these report files. No other files were written by this worker.

The distances are late-time background predictions with radiation and perturbations omitted. `H0=70` is only a distance-modulus gauge absorbed by the training intercept; no H0 or absolute-calibration inference is made. BAO and SN remain separate; no cross-probe independence claim or joint score is made. The BAO points are fixed points with no propagated shape uncertainty. Dovekie overlaps other SN compilations; none was added. These scores are not a posterior, evidence, calibrated p-value, physical-stability test, or full paper reproduction.

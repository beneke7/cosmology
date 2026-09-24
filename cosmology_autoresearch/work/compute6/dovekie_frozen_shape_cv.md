# Frozen BAO-shape prediction on DES-Dovekie

**Status:** exploratory result computed and independently reproduced PASS by a dense covariance-domain implementation (no production-score import or refit).

## Question and result

With the BAO-only shape points held fixed, does the interacting-vacuum background predict the separate DES-Dovekie supernova distance data better or worse than flat ΛCDM? On the predeclared four contiguous zHD holdouts, each with 455 held-out SNe and the released full STAT+SYS covariance, the sum of flat-prior-intercept-integrated conditional quadratic scores is 1640.3725 for ΛCDM and 1636.3359 for IVS. Thus

\[
\Delta_{\rm IVS-LCDM}=\sum_k[-2\log p_{\rm IVS}(y_{H_k}\mid y_{T_k})+2\log p_{\rm LCDM}(y_{H_k}\mid y_{T_k})]=-4.0366,
\]

where negative favors IVS in this **descriptive, overlapping leave-one-bin-out screen**. Conditional log determinants and held-out normalization are the same for the two frozen shapes within a fold, so the full normalized-score difference equals the reported Δχ² difference. The four conditional distributions are not independent and this sum is not a joint likelihood or significance.

The result is localized: the lowest-redshift holdout (zHD 0.02509–0.32394) contributes −5.5418; the other three folds contribute +0.7876, +0.3390, and +0.3787. Fixed-offset plug-in scores give a larger IVS advantage (Δ=−5.1975); integrating the training-offset uncertainty is the primary contract and softens that comparison. Six predeclared one-edge ±1%-of-N perturbations retain an IVS-favoring descriptive sum from −4.1180 to −3.0486. The split is stable to these small edge moves, but the gain is driven by the first fold rather than spread over all redshift blocks.

## Contract and numerical checks

- Source data and likelihood contract: the exact local Hubble diagram and packed precision match the official DES-SN5YR repository tree at commit `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`; see the [pinned input/likelihood audit](../data_audit4/dovekie_frozen_shape_contract.md). Its 1,820 SNANA rows are kept in release order, and the packed float32 product is the upper triangle of P=C⁻¹; no `MUERR` is added separately. The official reader’s CSV default is incompatible with the paired whitespace-format Dovekie file, so an order-preserving adapter is used.
- BAO shape points are loaded at full precision from `experiments/interacting_vacuum_screen/result.json` and checked against predeclared rounded coordinates: ΛCDM Ωm=0.2974618149542753, g=0; IVS Ωm=0.38697100769747017, g=−0.4666647299169934. Both recorded selected optimizer statuses are 0. There is no SN cosmological refit.
- For every fold the training GLS offset uses Q_T=P_TT−P_TH P_HH⁻¹P_HT, not P_TT. The primary predictive covariance includes the training-only flat-prior intercept term σ²_M vvᵀ; the fixed-offset result using S_H alone is secondary. Bins are contiguous in zHD and equal-redshift groups are never split.
- The actual 1,820×1,820 covariance was reconstructed by Cholesky solves. It is positive definite; max asymmetry before symmetrization is 9.7×10⁻¹⁷ and max |PC−I| is 7.52×10⁻¹⁵. On reference fold 1, the pre-score covariance-domain reference agrees at ≤1.14×10⁻¹³ in conditional χ². A separate implementation independently parsed the SNANA rows, reconstructed C=P⁻¹, integrated distances, and reproduced all fold maps and edge variants: max fold-score discrepancy 9.1×10⁻¹³; max variant-delta discrepancy 3.1×10⁻¹². See [independent result review](../compute8/dovekie_result_independent_review.md).
- Independently coded adaptive scalar distances agree with the main LCDM and IVS distances within 7.2×10⁻¹⁵ mag; the Γ=0 limit agrees with LCDM to 1.5×10⁻¹⁴ mag. An illustrative f_b=0.16 split remains positive through z=2.33. These checks concern the stated late-time background only—not perturbations, CMB, growth, or the paper’s modified-CAMB solver.

## Interpretation and limitations

This is a useful but weak cross-probe prediction result: the BAO-frozen IVS shape is not uniformly worse on the SN data and modestly improves the conditional screen, robustly under the specified edge perturbations. It does not establish a detection, posterior preference, calibrated p-value, physical model validation, or exact independence of BAO and SNe. No combined BAO×SN likelihood is formed because the public products provide no cross-probe covariance; the probes are reported separately. The Dovekie release substantially overlaps original DES-SN5YR, and no second SN compilation is included.

This also does not erase the prior BAO-block stress test: that separate seven-fold check found IVS worse by +22.7503 in its summed plug-in conditional χ², of which +22.0698 came from the z=2.33 pair. Evidence is therefore mixed across predictive diagnostics, not a broad victory for IVS.

The first run reached result assembly but failed only at JSON serialization (`ndarray` type); no comparative result was persisted from that attempt. The failure is preserved as `experiments/dovekie_frozen_shape_cv/attempt_001_failure.json`. Attempt 2 changed only JSON conversion and completed successfully in 1.554 s on one numerical thread.

## Reproduction

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
  .venv/bin/python scripts/run_bounded.py --seconds 1200 -- \
  .venv/bin/python experiments/dovekie_frozen_shape_cv/conditional_cv.py

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
  .venv/bin/python scripts/run_bounded.py --seconds 30 -- \
  .venv/bin/python experiments/dovekie_frozen_shape_cv/plot_result.py
```

The primary figure is [conditional_cv.png](../../experiments/dovekie_frozen_shape_cv/conditional_cv.png) (vector version: `conditional_cv.svg`). Machine-readable fold maps, all scores, covariance checks, input hashes and source hashes are in [result.json](../../experiments/dovekie_frozen_shape_cv/result.json); the independent disposition is in [independent_review_record.json](../../experiments/dovekie_frozen_shape_cv/independent_review_record.json). The root test suite passes 17/17 tests; synthetic precision-vs-covariance algebra also passes separately.

## Next discriminating test

After review of the no-score physical/domain contract, profile the same two-parameter late-time IVS background on Dovekie alone, retaining full STAT+SYS covariance and analytically profiling only the common magnitude offset. Compare its likelihood contour/profile with the BAO-selected point; do not combine probe likelihoods or interpret a profile minimum as evidence. This will show whether the modest frozen-shape prediction is compatible with an SN-preferred background shape or is mostly a single low-z block effect. Keep the z=2.33 BAO conditional failure as an explicit counter-signal.

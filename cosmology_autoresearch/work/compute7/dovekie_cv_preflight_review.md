# Dovekie frozen-shape CV preflight re-review

**Verdict: PASS for the pre-score implementation audit.** The revision addresses the previous holds: the actual-matrix covariance reference runs before the fold/scenario score loop, its tolerance failures produce a diagnostic-only record and exit, and the tolerance keys cover every returned reference check. No scoring script, fit, refit, or score-producing calculation was run. Thus the code gates are reviewed, while their numerical outcomes and the edge-perturbation stability outcome remain unobserved.

## Re-review of the prior HOLD items

| Item | Result | Evidence |
|---|---|---|
| Frozen BAO point provenance | PASS | `conditional_cv.py:92-134` loads exact coordinates from the BAO result, checks selected optimizer success/status and locked rounded targets, and records exact points, selected chi-square and status. The BAO result is hashed in the source list at lines 377-388. |
| Predeclared bin-edge perturbations | PASS (implemented; outcomes unobserved) | `EDGE_SHIFT_FRACTION=0.01` at line 46; lines 463-536 independently shift each of the three quantile targets by ±1% of N, preserve equal-z groups through `z_quantile_bins`, and record all six maps, fold results and comparative deltas. Lines 624-643 record the variants and delta range. I did not run these calculations and make no claim about their observed stability. |
| Normalized conditional scores/log determinants | PASS | Lines 271-288 calculate per-fold primary `chi2 + logdet(V) + |H|log(2π)` and plug-in `chi2 + logdet(S) + |H|log(2π)`. Fold sums are labelled descriptive, not a joint score, at lines 619-623. |
| Actual-matrix covariance-domain reconstruction and stop rule | PASS (static gate) | Lines 389-395 form `C=P^-1`, check symmetry and `PC≈I`, and Cholesky-factor C. Lines 297-355 reconstruct the reference fold in covariance space, factor `C_TT`, `S_H`, and the predictive covariance, and return offset, mean, rank-one, quadratic, logdet and pre-symmetrization asymmetry checks. The seven tolerance keys at lines 397-405 match the seven returned keys exactly. If any check fails, lines 414-439 write a diagnostic record with `comparative_scoring_performed=false` and exit before the primary/scenario fold loop at line 441; that record has no comparative score fields. |
| Numerical status, hashes, command | PASS (static reporting) | The failure branch records the failed status, tolerance/check values, input/source hashes and disposition at lines 419-439. The success path is reachable only after `reference_pass`; it records the bounded one-thread command, Python/NumPy/SciPy versions, code/data/contract hashes, and `pass_after_pre-score_covariance_reference` status at lines 545-665. |

The other audited items remain PASS by static inspection: precision-block signs and Schur-complement training precision; training-only GLS intercept and positive rank-one predictive variance; `zHD`/`zHEL` luminosity-distance convention and H0 cancellation; within-fold log-determinant comparability; and separate-probe/sample-overlap caveats.

This is a code-path review only. Before interpreting the experiment, its bounded run must still demonstrate that the covariance reference passes and report the actual edge-perturbation delta range; no result or comparative direction is asserted here.

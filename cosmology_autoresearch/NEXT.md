# Next experiment: profile IVS on DES-Dovekie alone

## Current evidence

The frozen BAO-shape prediction has passed independent covariance-domain reproduction. Four overlapping 455-SN conditional folds give integrated-intercept Δ(−2 log predictive) = −4.0366 (IVS−ΛCDM), with six locked edge variants all negative (−4.1180 to −3.0486). The lowest-redshift fold contributes −5.5418; the other three collectively favor ΛCDM by +1.5052. The result is a modest descriptive prediction diagnostic, not a joint BAO×SN likelihood, significance, posterior, or evidence. See [analysis](work/compute6/dovekie_frozen_shape_cv.md), [independent review](work/compute8/dovekie_result_independent_review.md), and [machine-readable result](experiments/dovekie_frozen_shape_cv/result.json).

Prior counter-signals remain material: the IVS BAO profile gain is exceeded by 256/1200 matched flat-null searches, and the BAO redshift-block conditional score penalizes IVS by +22.7503, dominated by the z=2.33 pair. A Dovekie-only profile is the next useful bounded check: determine whether the SN likelihood itself prefers a compatible background shape or whether the fixed-point result is mostly local to one low-z block.

## Gate before scoring

A no-score model/likelihood contract is complete and root-reviewed in `work/theory7/dovekie_ivs_profile_contract.md` and `.json`. The gate passes for a bounded exploratory run: the paper's g=Γ/H0 range is distinguished from the inherited code Ωm search interval; neither is silently interpreted as an SN posterior prior. The late-time background variables, analytic intercept profile, and positive-density plus existential baryon/CDM split checks over the observed SN redshift interval are specified. No optimizer or observed profile score has yet been run.

The next run can proceed under this reviewed contract:

1. Use only the pinned 1,820-row DES-Dovekie branch and exact full STAT+SYS covariance; no second SN compilation and no extra `MUERR` term.
2. Fit Ωm and g in the explicitly justified domain, profiling the one common additive magnitude intercept analytically by GLS. Keep the BAO profile separate; do not multiply likelihoods or claim exact probe independence.
3. Verify Γ=0 recovers flat ΛCDM, units/sign convention, covariance positivity, optimizer convergence from independent starts, physical feasibility, and numerical tolerance. Compare a direct covariance-domain score with an independent implementation.
4. Record the profile surface/contours, best fit, Δχ² relative to the Dovekie ΛCDM fit, and sensitivity to search bounds. Label the result exploratory; no evidence or calibrated significance without a separately designed null procedure.
5. Start one numerical thread and measure runtime first. Scale CPU parallelism only if runtime and independent-work structure justify it; do not use the GPU for this small fit absent an end-to-end benchmark.

## Standing limitations and source requests

- The IVS computation remains a background-only independent implementation; exact author-modified CAMB source/base commit, full Cobaya settings, and explicit perturbation initial conditions were not located. The paper PDF is already local (SHA-256 `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`). If the user has the missing author artifacts, request them and hash/audit before source-faithful claims; otherwise retain the [solver-path gates](work/perturbation_audit/solver_path.md).
- Do not append DESI Results IV Lyα values to the current BAO vector without an official replacement likelihood and cross-covariance; see [release audit](work/literature/desi_lya_release_probe.md).
- The campaign's root-owned durable state is `RUN_STATE.json`, `experiments/ledger.jsonl`, and `MORNING_REPORT.md`. After a reviewed checkpoint is committed and pushed, invoke the configured one-shot Astra Max Warden and include its brief advice in the running state.

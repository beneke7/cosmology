# Next experiment: test whether the BAO-selected IVS shape predicts an independent probe

## Question

Does the simple interacting-vacuum background selected by DESI DR2 BAO predict the separate DES-Dovekie Type Ia supernova distance-redshift data at least as well as the BAO-selected flat-ΛCDM shape?

This follows the evidence already collected: the IVS profile gain is common in 1,200 matched BAO null searches, its seven-fold BAO conditional score is worse overall and dominated by the z=2.33 pair, but it does not suffer a CDM-positivity crossing in the declared early-time sensitivity domain. The unanswered issue is whether the shape improvement travels to a genuinely different probe.

## Locked first-pass contract

1. Audit the existing DES-Dovekie files and likelihood definition, data/sample provenance, redshift ordering, covariance symmetry/positive-definiteness and hashes. Treat Dovekie as the one SN branch; do not add Pantheon+ or another overlapping compilation.
2. Freeze the BAO-only flat-ΛCDM shape Ωm=0.2974618150 and IVS shape (Ωm=0.3869710077, g=Γ/H0=−0.4666647299). Do not refit shape parameters on the SN data. Use the source-paper sign convention and the same explicitly background-only late-time model.
3. For each model and each predeclared contiguous redshift holdout, fit only the shared SN magnitude intercept on the training rows. Evaluate held-out distances and scores with the exact full covariance using the Gaussian conditional mean and Schur complement; never diagonalize or drop cross-block terms. Choose the number/bin edges from the observed redshift design before scoring. Store the fold map and all optimizer/linear-solve status.
4. Use luminosity distance with a freely fitted/marginalized magnitude intercept. This cancels the absolute H0 scale; report no H0 or absolute-calibration inference. Compare the held-out conditional scores descriptively and keep the SN result separate from BAO—no joint likelihood or combined evidence.
5. Validate the SN intercept solution against its analytic generalized-least-squares expression; independently reconstruct at least one fold by dense solves; check Γ=0 and distance units; verify numerical stability under equivalent covariance factorizations and bin-edge perturbations fixed before looking at comparative scores.

## Stop rules and interpretation

- Stop before any score if the released Dovekie likelihood/covariance, sample definition, or row mapping cannot be verified from local source records.
- Stop if full-covariance conditional calculations fail their independent reconstruction or positive-definiteness checks.
- A better/worse descriptive score is a predictive diagnostic, not a posterior, evidence, calibrated p-value, physical stability test, or full paper reproduction. Do not claim the two probes are statistically independent without documenting the measurement/systematic contract.
- First implementation/audit budget: 30 minutes CPU wall time with one numerical thread; explicitly cap every command. Use more CPU parallelism only if measured runtime warrants it. No GPU is expected to help this matrix/GLS calculation; benchmark before assigning one.
- Do not proceed to a Dovekie-wide IVS parameter refit or a mock-calibrated claim unless the frozen-shape prediction itself is stable and informative.

## Source-artifact gate in parallel

The paper identifies modified CAMB and a linear perturbation closure, but the exact patch/version, complete Cobaya configuration, and initial-condition prescription remain unverified in the public source trail. If the user supplies those exact artifacts, hash and audit them before any source-reproduction claim. Otherwise, any solver implementation is an independent model implementation and must pass the background, perturbation-source, and Γ=0/spectrum gates in `work/perturbation_audit/solver_path.md`.

Do not append DESI Results IV Lyα values to the existing BAO vector unless the official replacement/joint likelihood and cross-covariance become available. The current blocked-release evidence is `work/literature/desi_lya_release_probe.md`.

## References

- Current status and caveats: `MORNING_REPORT.md`.
- BAO-selected profile and full-covariance null: `work/compute3/interacting_vacuum_screen_report.md`, `work/compute3/interacting_vacuum_null_extension_report.md`.
- BAO block holdout and independent review: `work/compute4/interacting_vacuum_block_cv.md`, `work/compute5/interacting_vacuum_cv_review.md`.
- Existing DES-Dovekie profile/source analysis: `work/inference/REPORT.md`, `work/critic/REPORT.md`, and `experiments/dovekie_screen/`.

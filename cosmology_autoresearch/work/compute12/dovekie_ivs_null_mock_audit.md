# Independent audit: Dovekie IVS finite-search null mocks

Audit scope: read-only review of `work/compute11/dovekie_ivs_null_calibration.py` and its completed 200-realization output. This audit's only artifacts are in `work/compute12/`.

## Sampler derivation and numerical check

The release stores precision `P = C^-1`. For its lower Cholesky factor `P = L L^T`, a standard-normal vector `z` must be transformed as `epsilon = L^-T z`, i.e. solve `L^T epsilon = z`. Then

`Cov(epsilon) = L^-T I L^-1 = P^-1 = C`.

Solving `L epsilon = z` would generally be the wrong covariance. The driver uses `solve_triangular(L.T, z, lower=False)`, so its orientation is correct. [The independent full-matrix check](check_dovekie_precision_sampler.py) loaded the exact packed Dovekie precision, factored the 1,820-dimensional matrix, and generated 4,096 seeded draws. It found maximum relative whitening error `2.49e-15`, maximum relative precision-quadratic error `5.75e-15`, contrast sample/exact variance ratios from `0.979` to `1.017`, and maximum standardized contrast covariance error `0.027`.

## Driver review

The null mean uses the fitted LCDM shape `(Omega_m0, g) = (0.3303167905483107, 0)` and the fitted shared offset `0.007449125568204219 mag`. Adding this offset to each mock mean is valid; both model fits profile the same unbounded GLS intercept again. The score uses the full released precision, without an extra `MUERR` diagonal, and the statistic is `chi2_LCDM - chi2_IVS` against the requested observed value `1.4213186471747576`.

I found the searches aligned with the source profile: LCDM's `[0.05, 0.60]` bound and eight seeded starts, plus IVS's `[0.05, 0.60] x [-3, 3]` box, `41 x 41` grid, 16 starts (including that grid's minimum), physical-history check, `1e80` invalid score, and endpoint re-evaluation. Invalid IVS grid points and invalid optimizer endpoints are recorded; valid endpoints are selected even when an optimizer status is unsuccessful, matching the source's finite-search procedure. The per-replicate seed comes from a master `SeedSequence` child, so worker completion order does not change a mock's data or fit. JSONL rows are written in nondeterministic completion order, with replicate IDs and seeds preserved.

The analytic Clopper-Pearson interval implementation is the standard two-sided beta-quantile interval, and exceedances use the inclusive `T >= T_obs` rule. The finished run had 200 complete mocks and no failed realizations: 50 exceeded the observed statistic, for an empirical finite-search exceedance fraction of `0.25`, with interval `[0.191607, 0.315963]`. Invalid endpoints do not enter the selected statistic; the diagnostics report 1,400 invalid IVS endpoints and two raw unsuccessful IVS starts, as well as nine raw unsuccessful LCDM starts. This remains a finite-search null calibration only.

## Physicality spot check and wording

The mock driver stores the same sampled physicality check used by the search, while the source contract asks for a finer check at accepted fit points. I independently ran the source profile's 4,001-point `validate_point` check on all 200 accepted IVS mock fits: every point passed, with minimum `m=0.1872`, `x=0.3962`, `E^2=1`, and positive interior CDM witness (minimum `0.0936`). Maximum independent distance-integral discrepancy was `6.7e-16`. The numerical results therefore pass the finer check, though the driver/result does not currently record that per-mock validation.

For exact alignment with the requested interpretation, describe `50/200` as an “empirical null exceedance fraction” or “finite-search calibration rate” throughout; avoid naming the reported quantity `tail_probability_estimate` or implying a p-value. The current report says “not a posterior p-value,” which leaves unnecessary ambiguity about frequentist interpretation.

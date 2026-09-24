# Interacting-vacuum background screen

## Question and scope

Test whether a constant local interaction `Q = Γ ρ_x` in the vacuum case `w_x = −1` improves the 13-row DESI DR2 compressed-BAO profile fit relative to its `Γ = 0` limit. Use the supplied mean vector and full covariance unchanged. This is a low-redshift background-distance screen, not a reproduction of the paper's full CAMB/Cobaya analysis, perturbations, posterior constraints, or evidence.

## Source equations and conventions

The local source `context/papers/interacting_de_desi_dr2_2026.pdf` states in Eqs. (3)–(5), (8), pp. 2–4, that the flat expansion includes separate baryon, CDM, DE, radiation and neutrino terms; baryons, radiation and neutrinos are conserved separately; and

```text
dot(rho_x) + 3 H (1+w_x) rho_x = Q
dot(rho_c) + 3 H rho_c = -Q
Q = Gamma rho_x
```

Thus for `w_x = −1`, `dot(rho_x) = Gamma rho_x`; positive `Gamma` means CDM → DE, exactly as the paper specifies below Eq. (8). With `u = ln(1+z) = −ln(a)`, `g = Gamma/H0`, `M = Omega_b + Omega_c`, and `X = Omega_x`, the paper's low-redshift matter-plus-vacuum subsystem is

```text
dX/du = -g X/E
dM/du = 3 M + g X/E
E^2   = M + X
M(0) = Omega_m0 = Omega_b0 + Omega_c0
X(0) = 1 - Omega_m0
```

The separate baryon equation is `dOmega_b/du = 3 Omega_b`, and `Omega_c(u) = M(u) - Omega_b(u)`. This form follows by summing the separately conserved baryon and interacting-CDM equations; it preserves the paper's positive-`Gamma` convention.

### Explicit component approximation

Although Eq. (3) includes radiation and neutrinos, this screen omits them, following the repository's documented low-z BAO screening contract (`BRIEF.md`, lane A; `scripts/background_bao.py`). At `z ≤ 2.33`, no early-time or `r_d` prediction is attempted. Inserting fixed fractional radiation or massive-neutrino density would require an explicit `H0` / physical-density calibration, which the free BAO amplitude `alpha = c/(H0 r_d)` does not provide. The omission is explicit, not a hidden CMB prior; results are a late-time truncation of the paper's background equations, not a reproduction of its full model constraints.

## Locked profile domain and nuisance treatment

- Search `Omega_m0 ∈ [0.05, 0.60]`, exactly the flat-LCDM shape interval in `scripts/background_bao.py`, for a like-for-like BAO baseline. Flat closure gives `Omega_x0 = 1 - Omega_m0 > 0`.
- Search `g = Gamma/H0 ∈ [−3, 3]`, matching the explicit IVS range in the paper's Table I; here this is a finite profile-search box, not a prior weighting or posterior range.
- Profile `alpha = c/(H0 r_d)` analytically with exact baseline bounds `[1e−6, 1e4]`. `alpha` is dimensionless; `H0` and `r_d` remain separately unidentified.
- There is no CMB/BBN baryon-density prior. Present baryon fraction `f_b = Omega_b0/Omega_m0` is considered over all strictly positive values below 1. The observable background depends on `M = Omega_b + Omega_c`, not on the split. For each candidate `M(u), X(u)`, require a nonempty physical split interval `0 < f_b < f_b,max`, where

  ```text
  f_b,max = min(1, min_{0 <= u <= ln(3.33)} M(u)/(Omega_m0 exp(3u)))
  ```

  This is a mathematical existence/positivity check, not a measured or realistic baryon abundance. It profiles over the split by feasibility only; it does not infer `f_b`. An illustrative fixed-split check at `f_b ∈ [0.10, 0.20]` will be reported separately. This range is not an empirical prior, an inferred abundance range, or used in the profile. If it is not physically feasible at the selected fit, no realistic-split robustness will be claimed.
- Reject any ODE/integration failure, nonfinite state, `E^2 <= 0`, `M <= 0`, `X <= 0`, or candidate with no physically allowed baryon/CDM split. Keep such objective failures visible in the optimizer records.

## BAO likelihood and checks planned before inspecting the interacting fit

Use the exact ordered 13-row mean and full `13×13` covariance. For each shape pair `(Omega_m0, g)`, calculate unit-amplitude predictions `q` for

```text
DM/r_d = alpha * integral_0^z dz'/E(z')
DH/r_d = alpha/E(z)
DV/r_d = [z (DM/r_d)^2 (DH/r_d)]^(1/3)
```

Profile `alpha` as the bounded GLS solution

```text
alpha_hat = clip[(q^T C^-1 y)/(q^T C^-1 q), 1e-6, 1e4]
```

scoring the full-covariance whitened residual. Fit `Gamma = 0` independently on the same domain and compare its predictions/score to the existing flat-LCDM baseline. Before interpreting an observed extension fit, run zero-coupling recovery, synthetic nonzero-coupling recovery, baryon-split invariance and CDM-positivity checks, interaction-sign and distance-unit checks, ODE tolerance/convergence comparison, and an independent scalar-distance/derivative-free check. This contract precedes the interacting observed-data profile.

## Reproduction status

The contract was locked before the observed-data optimization. The component-split domain is a mathematical feasibility envelope, not a measured or independently validated baryon abundance; the BAO-only setup is conditional on omitting CMB/BBN calibration.

## Profile result

| Profile | Omega_m0 | Gamma/H0 | alpha | chi-squared | Notes |
|---|---:|---:|---:|---:|---|
| Gamma = 0, flat LCDM | 0.2974618150 | 0 | 29.52463331 | 10.271041002565 | same Omega_m interval and alpha profile as baseline |
| Interacting vacuum | 0.3869710077 | -0.4666647299 | 30.21858932 | 8.753877932904 | interior point; alpha and shape bounds inactive |

The profiled score difference is `chi2_flat - chi2_interacting = 1.51716306966`. The fitted coupling is negative, corresponding to DE → CDM with the paper's sign convention. The improvement is a finite profile difference for this one pre-specified background extension; no null mocks were run in this first task, and no significance/model-selection claim follows.

The main `g ∈ [−3,3]` box is the paper's Table-I search box, not a prior. Two nested search-window checks, `g ∈ [−1,1]` and `g ∈ [−0.5,0.5]`, independently return `chi2=8.753877932901` and `8.753877932908`, with `g=−0.4666653` and `−0.4666644`; neither optimum is on a coupling boundary. Relative to the full-box result, the score differences are below `4.1e−12`. The wider nested fit has 25 raw SciPy success flags, 20 valid finite physical endpoints, and 5 success-flagged invalid-penalty endpoints; the narrower has 25 raw successes, 23 valid endpoints, and 2 invalid-penalty endpoints. Both selected valid minima have status 0. These are finite search-window robustness checks, not prior restrictions or interval estimates.

At the best fit, `f_b,max = 0.78449`; the allowed positive-split interval is `0 < f_b < 0.78449`. For the separate illustrative diagnostic `f_b = 0.10, 0.15, 0.20`, both baryon and CDM densities stay positive through `z=2.33`, and the profile score remains exactly `8.753877932904` in each case. This only shows that these selected decompositions are physically feasible in the reduced background. BAO distances do not distinguish them, so this is not an abundance measurement, calibration, or evidence that the model's baryon sector is realistic. There is an exact `f_b` degeneracy in the background likelihood; the fitted profile is conditional on leaving the split uncalibrated and omitting BBN/CMB information.

## Numerical checks and optimizer record

- The flat `Gamma=0` predictions agree with `scripts/background_bao.py` to maximum absolute unit-prediction difference `1.11e-15`, with `E(0)=1`; its fitted score differs from the saved seed baseline by `-7.46e-13`.
- A 39-start L-BFGS-B interacting profile and a 41×41 coarse shape grid plus refinement agree within `4.95e-12` in chi-squared. The grid's best point refined to the same optimum. This does not prove globality. A separate Powell optimization using an independently adaptive scalar `quad` distance calculation agrees with the main score at the optimum within `2.1e-13` and the prediction vector within `1.1e-14`; its optimizer success/status is saved.
- Tightening the DOP853 ODE tolerances from `rtol=2e-10, atol=2e-12` to `2e-12, 2e-14` changes chi-squared by `2.72e-13` and predictions by at most `1.43e-14`.
- A zero-noise synthetic vector at `(Omega_m0, Gamma/H0, alpha) = (0.31, 0.20, 31)` recovers `(0.30999970, 0.20000153, 30.99999771)` with chi-squared `1.08e-11`. This is a self-consistency test, not an independent data validation.
- The sign checks give `E(z=1)=1.81627` at `Gamma/H0=+0.2`, `1.76068` at zero and `1.69887` at `-0.2`. Positive coupling sources DE and removes CDM in forward time as required by Eqs. (4), (5), and (8).
- With `H0=70 km/s/Mpc`, `r_d=147 Mpc`, and the explicit SI/Mpc conversion, `alpha=c/(H0 r_d)=29.13435` and is unchanged under `H0 -> 2 H0, r_d -> r_d/2`; all BAO predictions are dimensionless.
- The full covariance is Cholesky positive definite; its minimum eigenvalue is `0.00578998687`. The exact 13-row mean/covariance order and hashes are in the JSON.
- The optimizer statuses for every start are retained in `result.json`, and raw SciPy `success` flags are explicitly separated from accepted endpoints. The flat fit had 15 raw successes and 16 valid physical endpoints; its single status-2 abnormal endpoint was valid and lower than the best successful score by only `4.97e-14`, so the successful candidate was selected. For the main 39-start interacting fit, SciPy marked all 39 runs successful, but only 27 ended at a finite valid physical objective; the remaining 12 returned `success=True` at the `1e80` invalid-domain penalty. Thus “39/39 succeeded” is only a raw API status count, not a scientific convergence count. All 27 valid endpoints were raw successes. Under the declared normalized-coordinate connected-component tolerance of 0.005, those endpoints form 24 descriptive endpoint clusters, with 4 members in the best cluster; these are not certified attraction basins or proof of globality. The interacting optimizer recorded 105 invalid evaluations where `E^2 <= 0` and 3 with nonpositive states; all invalid endpoints/evaluations were rejected. The coarse grid marked 487 of 1,681 evaluations invalid (449 `E^2` failures, 38 nonpositive-state failures). No selected fit bound was active. The `Gamma/H0=-3` face had no valid point in its 41-point Omega_m boundary scan; the `+3` face minimum was `chi2=485.31`.
- A local finite-difference Hessian diagnostic at two step sizes gives inverse-curvature correlation `−0.99430` and ridge slope `dOmega_m/dg ≈ −0.1816`; both Hessians are positive definite and agree closely. This is only local geometry of the profiled score, not a posterior correlation, uncertainty interval, or calibrated uncertainty.

## Reproduction and artifacts

Command (bounded at 1,400 seconds; it completed well within the 1,500-second task cap):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python scripts/run_bounded.py --seconds 1400 -- \
  .venv/bin/python experiments/interacting_vacuum_screen/interacting_vacuum_profile.py \
  --coarse-grid 41 --starts 39
```

Recorded runtime was 24.17 s wall, one CPU process, one BLAS/OpenMP thread, no GPU. Environment: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, Linux x86_64.

- `experiments/interacting_vacuum_screen/interacting_vacuum_profile.py` — implementation; SHA-256 `ae99ae789aa02885667e27aa2cbbaa0b035323e69d1a4bf7bbd3848949fc5282`.
- `experiments/interacting_vacuum_screen/result.json` — full fit vectors, predictions, profile/domain details, every optimizer status, checks, input/artifact hashes; SHA-256 `e347b8934b2afdc23d23ef7d39c70f963b73888d17137f55caf7e51af87d43ae`.
- `experiments/interacting_vacuum_screen/profile_grid.npz` — 41×41 surface including invalid cells as NaN; SHA-256 `28171fc8f495da602686d58d47ecf84cc516413f3cfa2573d13dec2ccf09bca7`.
- `experiments/interacting_vacuum_screen/profile_surface.png` — SHA-256 `47670329c96fdb2c849cdc5f7ace20d6fd2dcf8fab59bd8e66ed26051e39a31a`; SVG — SHA-256 `4f814d09321b29a06e44bd3913be95a5ee10a8925d9384aa7ed69e9dba97855e`.

Key input hashes: paper PDF `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`; mean `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`; covariance `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`; baseline code `34bdf1d81397b39b740faa3e2f4a79748d5fa9d3bfdd8303619d631dbf8aa2e2`; saved baseline fit `92662f3ee42500d2be34d9f87b485fe56c82e9550fa674e42f4dec2bf2cab1f0`.

# Dovekie IVS local identifiability audit

**Status:** completed standalone no-fit local response audit at the pinned point
\((\Omega_{m0},g)=(0.39301742026989,-0.62938781451696)\). It is local linear GLS geometry,
not a posterior uncertainty, significance, or evidence calculation.

The SN model response derivatives were obtained by an independent RK45 forward-
sensitivity ODE. Central finite differences of a separately written base-only
ODE were checked over steps from `1e-2` to `1e-6`. The full 1,820 by 1,820
released STAT+SYS precision was unpacked in release order and factored as
`P = L L^T`; whitening used `L.T @ v`. One common magnitude intercept was
removed by an SVD rank-aware projection.

After removing the intercept, the signed cosine between the whitened
\(\Omega_{m0}\) and \(g\) responses is **0.96938919**
(raw, before intercept removal: 0.98032051). Singular values of the
two-column projected response are **74.1801, 1.8118**,
with condition number **40.9428** in the stated dimensionless
parameter coordinates.

After projecting out both the intercept and the \(\Omega_{m0}\) response,
the conditional local information for \(g\) is **3.3137374**. The
Schur-complement value is 3.3137374; their absolute difference is
2.80e-14. The finite-difference sweep's best whitened relative derivative
errors are 3.50e-10
for \(\Omega_{m0}\) and 6.86e-10
for \(g\).

The independent GLS intercept is 0.017690169 mag and the normal-equation
residual is 1.02e-12; this re-evaluation only supports the
descriptive residual alignment and is not a new fit. The optional whitened
residual cosine with the nuisance-orthogonal \(g\) direction is
1.58181e-07, reported descriptively only.
The independently recomputed fixed-point \(\chi^2_{\rm prof}\) is
1629.999216924067, matching the pinned profile record to
6.82e-13 if that record is available.

The independent central-difference convergence sweep is:

| Step | Whitened relative error, \(\Omega_{m0}\) | Whitened relative error, \(g\) | Projected response condition |
|---:|---:|---:|---:|
| 1e-02 | 8.54e-05 | 5.56e-06 | 40.9581 |
| 3e-03 | 7.69e-06 | 5.00e-07 | 40.9442 |
| 1e-03 | 8.54e-07 | 5.56e-08 | 40.9429 |
| 3e-04 | 7.69e-08 | 5.01e-09 | 40.9428 |
| 1e-04 | 8.55e-09 | 6.86e-10 | 40.9428 |
| 3e-05 | 7.82e-10 | 1.41e-09 | 40.9428 |
| 1e-05 | 3.50e-10 | 4.19e-09 | 40.9428 |
| 3e-06 | 1.08e-09 | 1.36e-08 | 40.9428 |
| 1e-06 | 3.14e-09 | 4.30e-08 | 40.9428 |

The alignment, singular values, condition number, and conditional information
describe the chosen local response coordinates. They do not establish a
posterior width, calibrated constraint, evidence, or physical viability beyond
the scoped low-redshift model. Full numerical values and every finite-difference
step are in `dovekie_profile_identifiability.json`.

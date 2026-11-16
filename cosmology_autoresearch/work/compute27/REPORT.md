# CAMB BAO conditional H0-support diagnostic

## Question and disposition

Does the nonconverged broad CAMB profile's candidate owe its high (H_0\simeq99.16\) to the theta-to-H0 root solver's upper bracket at 100 km s\(^{-1}\) Mpc\(^{-1}\)? **Not at fixed candidate physical densities.** Directly specifying (H_0) in flat CAMB and profiling it conditionally gives a nearby interior minimum with effectively unchanged full-covariance BAO score. This removes one suspected solver-bracket artifact, but it does **not** establish convergence or global optimality of the three-parameter profile.

## Method and checks

Held ω\(_b h^2\)=0.060691736587168255 and ω\(_c h^2\)=0.23095938241922825 at the candidate from `work/compute24/camb_bao_lcdm_result.json`; fixed flat geometry, (w=-1), (m_\nu=0.06\) eV, (N_{\rm eff}=3.044\), (T_{\rm CMB}=2.7255\) K, one degenerate massive neutrino and CAMB BBN-consistency helium. Evaluated (H_0=80,90,95,99,99.1588,100,105,110,125,150,200) directly, then used bounded scalar minimization on [95,105]. Read the released 13-row DESI DR2 mean and full 13x13 covariance independently of the production scorer; checked dense covariance solve against Cholesky and CAMB (D_M(z)) against adaptive integration of (c/H(z)).

The conditional minimum is (H_0=99.16241760\) km s\(^{-1}\) Mpc\(^{-1}\), ω\(_b h^2\), ω\(_c h^2\) fixed as above, with χ\(^2=10.27630436\). The previous theta-root point scores 10.27636166; at H\(_0=100\) the score is 13.2775, and it rises rapidly over the rest of the sampled range. This is a conditional, non-posterior diagnostic with broad uncalibrated BAO-only physical densities, not an inference of (H_0\).

Numerical checks: maximum radial-distance integration discrepancy (1.37\times10^{-12}\) Mpc across the 13 rows; dense-vs-Cholesky score differs by (1.8\times10^{-15}\). Runtime 0.077 s. CAMB 2.0.4 (wheel SHA-256 `9b7e6558b2049757e96496b8a33461d8ff11acd79017cfbd00a573dbc8d1e366`), NumPy 2.5.3, SciPy 1.18.1. DESI mean SHA-256 `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`; covariance SHA-256 `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`.

## Reproduction

From the repository root:

```sh
PYTHONPATH=work/compute24/site OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_bounded.py --seconds 600 --state RUN_STATE.json -- \
  .venv/bin/python work/compute27/camb_h0_support_test.py
```

The output is `camb_h0_support_test.json`. Script SHA-256: `3158e2b8ec5e9186d567ac0d657bdc3c073757f9338197a93ca23284156388a6`; result SHA-256 is recorded in the experiment ledger after independent profile review.

## Limitation / next discriminating check

The full ω\(_b h^2\), ω\(_c h^2\), and (H_0\) profile still needs independent converged multi-start reproduction with support/boundary checks. The fixed-density scalar minimum cannot substitute for that. The broad Table-I boxes are search limits here, not justified posterior priors; no CMB or BBN observational likelihood is included.

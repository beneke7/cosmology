# Independent compute26 boundary audit

Status: independently checked. The compute26 point is supported as the constrained profile minimum at the upper \(\omega_bh^2=0.1\) edge of its declared direct-\(H_0\) box. It is not an interior fit and remains unaccepted under the campaign's interior-minimum rule. Reoptimized points above \(0.1\) have lower scores, but lie outside the declared density bound and are exploratory only.

## What resolved the apparent inconsistency

The compute26 report's fixed-\((H_0,\omega_ch^2)\) outward stress changes \(\omega_bh^2\) while holding nuisance coordinates fixed. At the exact upper edge and the compute26 \(H_0,\omega_ch^2\), the independently rebuilt score is 10.273886272854. Two-sided finite differences show a minimum along that fixed-coordinate slice:

| Step in \(\omega_bh^2\) | \(\chi^2(0.1-h)\), in box | \(\chi^2(0.1)\) | \(\chi^2(0.1+h)\), outside box | Left/right secant slopes |
|---:|---:|---:|---:|---:|
| 0.001 | 11.08151608 | 10.27388627 | 11.07531828 | −807.630 / +801.432 |
| 0.0001 | 10.28194016 | 10.27388627 | 10.28192530 | −80.539 / +80.390 |
| 0.00001 | 10.27396709 | 10.27388627 | 10.27396641 | −8.082 / +8.014 |

After reoptimizing \(H_0\) and \(\omega_ch^2\) at each fixed baryon density, the profile moves along a shallow degeneracy and continues downward through the bound:

| Fixed \(\omega_bh^2\) | Status | Profile \(\chi^2\) | Reoptimized \(H_0\) | Reoptimized \(\omega_ch^2\) |
|---:|---|---:|---:|---:|
| 0.0800 | In box | 10.27485540 | 115.43574 | 0.31560768 |
| 0.0950 | In box | 10.27409767 | 128.91942 | 0.39861540 |
| 0.0990 | In box | 10.27392714 | 132.66600 | 0.42376827 |
| 0.0999 | In box | 10.27389034 | 133.51864 | 0.42962049 |
| 0.1000 | Upper edge | 10.27388627 | 133.61342 | 0.43027242 |
| 0.1001 | Outside box | 10.27388222 | 133.70843 | 0.43092815 |
| 0.1050 | Outside box | 10.27369183 | 138.41814 | 0.46415161 |
| 0.1200 | Outside box | 10.27319394 | 153.56600 | 0.58071696 |

Thus the fixed-coordinate rise beyond the boundary is not contradictory. Along the nuisance-reoptimized profile, the score falls on both sides because \(H_0\) and \(\omega_ch^2\) adjust together. The in-box profile approaches its minimum at the allowed upper edge; extending the baryon range would shift the minimum outside the declared box.

## Search and checks

I reimplemented CAMB parameter construction and the DESI Gaussian score in `local_boundary_audit.py`; the script does not import compute26's likelihood code. It uses the local pinned CAMB 2.0.4 site, float64, the full 13×13 covariance, flat positive-\(\Lambda\) closure, and the same fixed neutrino/temperature/BBN conventions. Three independent full-box differential-evolution starts (seeds 20260929–20260931) converged to the same in-domain point, \(\chi^2=10.273886272853535\), \((H_0,\omega_bh^2,\omega_ch^2)=(133.61341630,0.1,0.43027242194)\). None found a valid in-box candidate with lower \(\chi^2\). The fixed-\(\omega_bh^2\) nuisance profiles at all tested in-box slices also decrease toward 0.1; no lower in-box slice was found.

Every reported profile minimum passed CAMB closure: \(\Omega_{m0}\simeq0.2974\), \(\Omega_{de0}\simeq0.7026>0\). Inverse CAMB theta mapping at all eight slices recovers \(100\theta_{MC}\) between 1.0944 and 1.1041, inside the target paper's broad [0.5,10] interval; expanded \(H_0\in[1,2000]\) roots reproduce the fitted \(H_0\) to at most about \(3\times10^{-10}\) km/s/Mpc. The original [10,100] root bracket fails for every profiled point because all fitted \(H_0>100\).

The independent full-covariance dense and Cholesky scores agree. Independently integrating \(c/H(z)\) at the challenge point reproduces CAMB \(D_M(z)\) to \(9.1\times10^{-13}\) Mpc. Invalid search points were retained: 95 positive-flat-\(\Lambda\) closure failures and 15 CAMB recombination failures. They were never eligible as candidates.

The outside-box rows are diagnostics, not an extension of the declared fit range. The optimizer agreement supports, but does not mathematically prove, globality over every continuous point in the box. Per the standing criterion, the boundary candidate is not an accepted interior minimum. This remains a background-only DESI DR2 BAO profile; it is not a Planck/ACT likelihood, posterior, or evidence calculation.

## Reproduction

Run from the workspace root; the wrapper capped runtime at 1100 seconds. Actual runtime was 86.2 seconds using one CPU worker and one numerical thread:

```sh
cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute29/task_state.json --seconds 1100 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute29/local_boundary_audit.py
```

Versions: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, CAMB 2.0.4; Linux x86_64. No install or network request was used. Complete seeds, optimizer traces, profile values, theta roots, invalid-point samples, checks, and hashes are in `local_boundary_audit.json`; `local_boundary_audit.png` plots the fixed-density profile and marks out-of-box diagnostics.

SHA-256 inputs:

- DESI mean: `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`
- DESI covariance: `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`
- CAMB 2.0.4 wheel: `9b7e6558b2049757e96496b8a33461d8ff11acd79017cfbd00a573dbc8d1e366`
- Target paper PDF: `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`


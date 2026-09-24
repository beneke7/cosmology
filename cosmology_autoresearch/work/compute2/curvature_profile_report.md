# Curved-LCDM compressed-BAO profile screen

## Bottom line

On the unchanged 13-row DESI DR2 compressed-BAO mean vector and full covariance, profiling a standard FRW curvature parameter lowers the best-fit score by only

\[
\Delta\chi^2 \equiv \chi^2_{\rm flat}-\chi^2_{\rm curved}=0.3173631.
\]

Under 1,200 full-covariance parametric mocks generated from the saved fitted flat-\(\Lambda\)CDM null, this much or more finite-search improvement occurred in 692/1,200 runs (0.5767; exact two-sided 95% Clopper–Pearson interval 0.5481–0.6048). The pooled mock \(q_{95}\) is 3.822. Thus curvature improves this particular profile minimum, but the amount is common under this fitted-null calibration. This is a compressed-BAO profile/bootstrap result, not a posterior p-value, model evidence, or discovery significance.

## Fit and domains

Used the full SPD 13×13 covariance without changing the data or row order. The minimum covariance eigenvalue is 0.00578998687. The model is

\[
E^2(z)=\Omega_m(1+z)^3+\Omega_k(1+z)^2+\Omega_\Lambda,\qquad
\Omega_\Lambda=1-\Omega_m-\Omega_k,
\]

with \(\Omega_m\in[0.05,0.60]\), primary \(\Omega_k\in[-0.20,0.20]\), and nested \(\Omega_k\) domains ±0.10 and ±0.05. Positive curvature parameter means open geometry: \(S_k(\chi)=\sinh(\sqrt{\Omega_k}\chi)/\sqrt{\Omega_k}\); negative means closed and uses sine. Flat closure and radial/transverse/volume BAO distances follow DESI DR2 Results II, Eqs. (3)–(6), and `work/theory2/curvature_contract.md`. The independent BAO amplitude \(\alpha=c/(H_0r_d)\) was profiled analytically and clipped to the exact [1e−6, 1e4] bounds. All calculations use float64.

| Fit | \(\Omega_m\) | \(\Omega_k\) | \(\Omega_\Lambda\) | \(\alpha\) | \(\chi^2\) | \(\Delta\chi^2\) vs flat |
|---|---:|---:|---:|---:|---:|---:|
| Flat ΛCDM | 0.2974618173 | 0 | 0.7025381827 | 29.52463336 | 10.27104100256 | — |
| Curved ΛCDM, ±0.20 | 0.2931444835 | +0.0225169635 | 0.6843385531 | 29.60748965 | 9.95367789913 | 0.31736310343 |
| Curved ΛCDM, ±0.10 | 0.2931444825 | +0.0225169684 | 0.6843385490 | 29.60748967 | 9.95367789913 | 0.31736310343 |
| Curved ΛCDM, ±0.05 | 0.2931444810 | +0.0225169760 | 0.6843385430 | 29.60748969 | 9.95367789913 | 0.31736310343 |

No fitted parameters hit bounds. The flat result reproduces the saved baseline score to \(-9.54\times10^{-13}\), with exact prediction agreement when evaluated at the saved seed parameters. The 16-start full-domain result was stable against 64 starts (score difference \(1.39\times10^{-13}\)); a 61×61 coarse grid located a basin whose local refinement returned the same minimum. These checks support numerical stability but do not establish globality.

## Null calibration

Each mock is a fitted flat-ΛCDM prediction plus a standard-normal vector transformed by the Cholesky factor of the complete observed covariance. The 200- and 1,000-mock standard-normal arrays were regenerated using NumPy `default_rng` / PCG64 at seeds 20260924 and 20260925. They are bitwise reproducible per seed and have distinct stream digests. Both flat and curved fits were rerun for every mock using the same fixed search contract: 16 deterministic L-BFGS-B starts per model, \(\Omega_k\in[-0.20,0.20]\); the curved starts include the corresponding flat optimum at \(\Omega_k=0\). The criterion is \(\Delta\chi^2\ge0.31736310343\).

Specifically, the observed flat profile and all 1,200 flat-null and curved mock fits were independently recomputed for this curvature statistic. Prior mock vectors were recreated from the recorded seeds and covariance rather than taking prior fits' chi-square values as given. The earlier unrelated flat/CPL robustness ensembles and their model comparisons were not rerun.

| Seed | Mocks | Exceedances | Fraction | Exact CP 95% interval | Mock \(q_{95}(\Delta\chi^2)\) |
|---:|---:|---:|---:|---:|---:|
| 20260924 | 200 | 113 | 0.5650 | [0.4933, 0.6348] | 3.8752 |
| 20260925 | 1,000 | 579 | 0.5790 | [0.5477, 0.6098] | 3.7519 |
| Pooled disjoint streams | 1,200 | 692 | 0.5767 | [0.5481, 0.6048] | 3.8220 |

The pooled median, 90th, 95th, and 99th percentiles are respectively 0.4277, 2.7040, 3.8220, and 6.7919. The 200- and 1,000-mock summaries and exact pooled values are in `result.json`; each generated mock's statistic and optimizer records are in `mock_realizations.jsonl`.

There were 52 unsuccessful starts among 19,200 flat starts and 8 among 19,200 curved starts. The unsuccessful-start statuses remain recorded. The selected best candidate was unsuccessful in 4 flat and 2 curved mock fits; all six received successful bounded Powell rechecks. The most favorable recheck changed the score by only about \(2\times10^{-13}\), with no retry improving by more than \(10^{-7}\). No claim that every non-selected failed start is an optimum is made.

## Numerical and domain checks

- At \(\Omega_k=0\), predictions match `scripts/background_bao.py` to max absolute error \(3.55\times10^{-15}\); \(E(0)=1\). For \(\chi=1.2\), open/closed curvature kernels lie respectively above/below the flat value, consistent with the sign convention.
- Analytic profiling of \(\alpha\) agrees with bounded scalar minimization to below \(10^{-13}\). At both observed optima, 96- versus 192-node Gauss–Legendre predictions differ by less than \(10^{-13}\); an independent scalar adaptive-quadrature Powell reoptimization agrees in score at this scale.
- Under the specified domains and \(z\le2.33\), the contract's continuous lower bound is \(E^2_{\min}=73/108\simeq0.675926>0\); checked numerically at the analytic minimum, and sampled over the domain. Closure has \(\Omega_\Lambda\ge0.2\).
- The 13 rows include the Lyα BAO block at \(z=2.33\). DESI DR2 Lyα full-shape work reports BAO information from the same dataset; it must not be appended as an independent likelihood without the relevant joint covariance.

This screen concerns background distances inferred from the compressed BAO vector only. It does not implement early-time calibration or perturbations, and it does not provide posterior sampling or Bayesian evidence.

## Reproduction and artifacts

Command (the bounded runner's 1,450 s ceiling was not approached):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python scripts/run_bounded.py --seconds 1450 -- \
  .venv/bin/python experiments/curvature_screen/curvature_profile.py --workers 20
```

Recorded runtime was 6.364 s wall, with 20 worker processes, 20 CPUs in affinity, one BLAS/OpenMP thread, CPU only. Runtime environment: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, Linux x86_64.

Review artifacts:

- `experiments/curvature_screen/result.json` — complete fit results, bounds, per-start observed search records, checks, source hashes, pooled and per-seed summaries. SHA-256: `5457e02073b51a1c3905007f8263698a403c910e983e54ce632437a03d99350c`.
- `experiments/curvature_screen/mock_realizations.jsonl` — all 1,200 mock statistics and flat/curved per-start/recheck statuses. SHA-256: `812167e6b5f9dbc8d164e20180f421a980e1c793ef5afa35662da89b15c3752f`.
- `experiments/curvature_screen/curvature_profile.py` — standalone implementation. SHA-256: `191d933c6c263a767e9d87c2aee294066e3e3c023c64e6de54886acc59491108`.
- `experiments/curvature_screen/curvature_null_calibration.png` SHA-256: `fe0711d42bee977fe4f285f24a6a5c31e125c5df59896d8a5c46cc3cf02c6afb`; SVG SHA-256: `681f09c126e8918b54e5aed39bd42e58fc7e051074ab181877251f2da4345f63`.

Mean/covariance hashes are `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` and `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`. The result JSON records hashes for the seed fit, both input bootstrap results, code, DESI source papers, and curvature contract; it confirms the input data/code hashes stored in the seed bootstrap files match the current files.

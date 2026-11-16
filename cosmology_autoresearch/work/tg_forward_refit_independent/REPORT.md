# NGC 3198 isolated spherical-surrogate profile

**Disposition: exploratory, independently implemented conditional test; not a physical detection or a complete falsification of thermodynamic gravity.** The numerical branch and score are reproducible. Their interpretation is limited by the source model: SPARC rotation-curve components do not define a unique three-dimensional baryon density, and the spherical inversion used here contains robust negative-density intervals.

## Question and implementation

I tested whether the printed constant- K equation admits a regular, isolated exterior solution for an explicitly declared spherical surrogate built from all 43 local NGC 3198 SPARC Rotmod points. The fit statistic is the direct diagonal velocity score, `sum(((Vobs - Vmodel)/sigma)**2)`. No cumulative force-identity statistic, endpoint velocity boundary condition, density clipping, or covariance matrix was used. The data product supplies pointwise errors, not a covariance; any unknown correlated/systematic errors remain outside this screening likelihood.

The input Pszota–Ván PDF is the user-provided local copy (SHA-256 below). Its Eq. (47), visually checked in the companion source audit, has the plus-sign spherical inversion

`rho = [Upsilon d(r Vdisk^2)/dr + d(r Vgas |Vgas|)/dr] / (4 pi G r^2)`.

Here `Upsilon` multiplies the stellar-disk term only; `Vgas |Vgas|` preserves SPARC's signed gas convention, including its negative gas rows. With `r` in kpc, velocities in km/s, and `G = 4.30091e-6 kpc (km/s)^2/Msun`, the density unit is Msun/kpc^3 and K is s^2/km^2. The paper prints an unsigned square in Eq. (47), while its later discussion (Eq. 53) preserves negative gas contributions; this implementation follows the latter convention. See [the source-contract audit](../tg_sparc_source_audit/REPORT.md) for equation-page, interpolation, and sign details.

For `x=(1-u)/K`, I independently integrated

`x'' + 2 x'/r = 4 pi G rho (1 - K x)`, with `v^2 = r x'/(1-Kx)` for K>0 and the exact Poisson/Newtonian limit `v^2 = r x'` at K=0. The center is regularized with three declared continuations `V_i(r)=V_i(r0)(r/r0)^p`, p=1, 1.5, 2, below the first datum at 0.32 kpc. Beyond the last datum at 44.08 kpc, each signed pseudo-density component is exponentially tapered with L=5 kpc. The finite exterior match is at R=104.08 kpc (12 taper lengths beyond the last point), using `u'(R)=-(u(R)-1)/R`, equivalent to `u -> 1` at infinity in a vacuum exterior. No observed outer velocity enters this condition. L=10 kpc was checked at each selected candidate without refitting.

The predeclared domain was K in [0, 0.01] s^2/km^2 (including K=0; upper edge 294 times the paper's quoted K), Upsilon in [0.2,1.0] with a five-value grid {0.2,0.5,0.7,0.762,1.0}, and the three central continuations. Each of the 15 (p,Upsilon) groups profiled 17 deterministic fractions of its connected physical K branch. First-branch thresholds were found by a 12-node logarithmic scan, 12 normal-tolerance bisections, then six tight-tolerance bisections. Every point was re-evaluated with an event for `x'=0`; a candidate with nonpositive speed anywhere in the integrated source region is rejected. Positivity of `u`, the exterior normalization, and positive observed `v^2` were also enforced. All 15 tight lower/upper threshold checks passed.

## Results

The minimum of the completed coarse profile surface is at the reference mass-to-light value Upsilon=0.762 and the 0.97 threshold fraction for all three center continuations:

| Inner continuation p | Grid-minimum K (s^2/km^2) | First invalid K (s^2/km^2) | K/Kcrit | Direct chi-square | Tight chi-square |
|---:|---:|---:|---:|---:|---:|
| 1.0 | 3.71051e-5 | 3.82529e-5 | 0.970 | 92.1742 | 92.1741 |
| 1.5 | 3.71112e-5 | 3.82592e-5 | 0.970 | 92.3791 | 92.3769 |
| 2.0 | 3.71147e-5 | 3.82628e-5 | 0.970 | 92.4924 | 92.4982 |

The lower side of each tight threshold bracket is valid; the upper side first fails the exterior normalization. The initial coarse scan's first rejected node instead reports `u_positivity_threshold`, because it lies well above the refined branch edge. At Upsilon=0.762, the threshold is about 1.125 times the paper's quoted K and the best sampled K is about 1.091 times that value. The score minimum is only about 3% below the isolated-branch edge. Across the five-point Upsilon grid, the best score occurs at 0.762; this is not a continuous Upsilon profile and the reference value was included deliberately.

The K=0 baryon-only score is 19002.23 for the same data and pointwise errors; the paper K=3.40e-5 gives scores around 885–890 for the three center continuations at Upsilon=0.762. Thus this source prescription produces a very large descriptive improvement over baryons alone, but its best sampled total chi-square is still about 92 on 43 points (roughly twice the point count). This is not a satisfactory precision fit by a naive independent-error reading, and neither a formal significance nor evidence for the theory follows: the source is conditional, Upsilon is only coarsely profiled, there is no reported covariance, and Newtonian baryons alone are not a competitive halo/MOND baseline.

The fine profile is sharply sensitive near the branch edge: at p=1.5 and Upsilon=0.762, chi-square is 396.1 at K/Kcrit=0.92, 92.38 at 0.97, 260.1 at 0.99, and 419.3 at 0.999. The 0.97 grid point is therefore a sampled local minimum, not proof of the continuous optimum. Tight candidate solves shift chi-square by at most 0.0059. Re-anchoring the whole grid to the refined thresholds changed some high-Upsilon near-edge scores by as much as 41.6; this is recorded, not hidden, and demonstrates the conditioning problem close to the singular normalization edge.

The L=5 to fixed-parameter L=10 kpc tail change shifts the candidate chi-square by at most 4.3e-6. Changing the inner continuation shifts the best grid score by only 0.32, but neither test removes the outer normalization limitation. A separate Riccati integration of the original acceleration equation agrees with the transformed ODE to at most 3.6e-8 relative on the candidate data radii; the matched exterior constants are positive (attractive branch).

The K=0 Poisson solution recovers the signed baryon `v^2` curve to at worst 1.32e-7 fractional error across the three center continuations. Tight scores approach their exact K=0 values continuously: at K=1e-12 s^2/km^2, the largest absolute score shift is 0.0013. This checks the singular-looking K transform without dividing by K.

## Source physicality and negative result

The signed PCHIP spherical inversion at the reference Upsilon=0.762 has negative pseudo-density over approximately 14.498–15.997 kpc and again over 43.057–44.08 kpc; the minimum sampled value is -1.56e5 Msun/kpc^3. The inner negative interval agrees with natural-cubic, Akima, and finite-difference checks in the companion audit. The outer interval is more sensitive to endpoint interpolation. No negative density was clipped. This is a serious limitation of interpreting this curve inversion as a physical spherical matter source; it is not a unique 3-D reconstruction and does not rule out another geometry/source prescription.

**Negative result / boundary:** the public local files do not supply the resolved 3-D stellar-plus-gas source, its inner continuation, or its outer taper. The calculation therefore establishes only that the declared signed spherical surrogate admits a narrow, low-score, positive isolated branch. Its near-threshold sensitivity and negative source intervals prevent promoting the result to a physical detection. This test cannot validate a cosmological or relativistic completion.

## Provenance, checksums, compute, and reproduction

- Exact user-provided Pszota–Ván PDF SHA-256: `55cbfbbb0a916ea2ce0b70523334ee301897695860fb4bdb5bfd2278ab47b62e`.
- Local `Rotmod_LTG.zip` SHA-256: `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`; exact `NGC3198_rotmod.dat` member SHA-256: `17b774ad11e7dd745a067073b76f1909d4de32fa87616df12b64da3f225cf953` (43 rows, 2,075 bytes).
- Frozen first-version independent implementation [profile.py](profile.py) SHA-256: `b7e097e6bfe020b967074d66fb8d4835bd6a57c8f6bb6e190bf6bbec7b5b462a`; unchanged after the initial serial run. Final parser/solver [profile_checkpointed.py](profile_checkpointed.py) SHA-256: `7508b62025efd9c3ba5b00c8a64c36a440add66042d38461679141fc947d0df8`. Its hash is also embedded in the JSON record.
- No random numbers were used (seed: none). Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, DOP853; normal tolerances rtol=2e-8, atol=1e-8, max step 0.35 kpc; tight rtol=atol=1e-11, max step 0.06 kpc.
- The failed-to-checkpoint first v1 attempt was stopped only at its task-owned PID after 979 s process CPU (16:19 snapshot). The replacement used 15 independent worker processes over 20 allowed CPUs and capped BLAS/OpenMP at one thread. Its derivative-cache microbenchmark reduced five equivalent solves from 0.918 to 0.504 CPU s (~1.82x); the small sequential adaptive ODE has no demonstrated GPU benefit, so GPU was not used.
- Conservative cumulative numerical CPU accounting is 1181.08/1200 s, including the first stopped attempt, profile and verification workers, parent CPU, measured solver benchmarks, and a 2 s guard. The checkpointed attempt's wall time was 17.34 s. All 15 groups contain 17 rows; [checkpoint.jsonl](checkpoint.jsonl) and `groups/*.json` preserve resumable results.

From `cosmology_autoresearch/`, the reproducible command is:

```sh
.venv/bin/python scripts/run_bounded.py --state /tmp/tg-no-campaign-deadline.json --seconds 120 -- env OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python work/tg_forward_refit_independent/profile_checkpointed.py
MPLBACKEND=Agg .venv/bin/python work/tg_forward_refit_independent/plot_profile.py
```

The machine-readable results are [profile_checkpointed.json](profile_checkpointed.json) and [profile_checkpointed_surface.csv](profile_checkpointed_surface.csv); the plot is [profile_summary.png](profile_summary.png), generated by [plot_profile.py](plot_profile.py) (SHA-256 `92c239971e9300f85c2844646514487b1065cbf460b48bb9d078aa04a69eb25c`).

## Ranked next experiment

1. **Most informative:** replace the spherical inversion with an explicit, resolved disk-plus-HI source and a defensible geometry, then redo a dense K/Upsilon profile with an independently justified isolated boundary and whole-galaxy holdout. First ask whether that source is available; the local Rotmod archive is not enough. This determines whether the near-normalization-edge score is a physical prediction or a surrogate artifact.
2. Refine the current conditional profile around K/Kcrit=0.93–0.995 and continuously profile Upsilon near 0.762, including numerical-mesh/exterior-radius convergence. This resolves the sampled minimum, but cannot repair the missing physical source.
3. Compare only after the source is physical: matched Newtonian+halo and MOND fits, and then cross-galaxy prediction. A rotation-curve result would remain separate from a relativistic/cosmological viability claim.

# Finite-u TG profile at H I-supported NGC 3198 geometry

**Finding — exploratory finite-u in-sample screen.** The run profiles actual positive-u solutions over the declared finite K/Kcrit grid at i=72° and D={13.8,16.6} Mpc. The singular Kcrit limit is used only to set the scale and is never scored as a finite point. The result is descriptive and establishes neither a detection nor model-selection significance.

## Grid and source support

The source builder is the geometry-aware HALOGAS H I plus S4G old-stellar map pipeline used by the existing nuisance envelope. THINGS adopts i=72° for NGC 3198; HALOGAS describes i≈70°. This reduced run uses the THINGS-supported 72° node only, so the HALOGAS 70° node remains a geometry limitation. D uses only the existing 13.8±1.4 Mpc SPARC prior nodes at 0 and +2σ (13.8 and 16.6 Mpc). Υ* uses existing prior/grid nodes {0.4766, 0.6, 0.762}; all three established vertical source cases are included. There are 18 source cases and eight finite K/Kcrit fractions per case.

Every source case recomputes Kcrit from the principal generalized mode of the checked finite-volume operator. At each K fraction the actual B−K·diag(D) system is factored, u is checked on the full half-domain and extrapolated midplane, and the curve is evaluated at all 43 SPARC radii. The diagonal score is χ² plus the same D/i/log-M/L prior terms as the envelope. K has no prior. NFW and fixed-a0 simple MOND controls are recomputed directly from the same transformed source, distance-scaled radii, and inclination-corrected velocities at each exact nuisance node; NFW optimizer success is recorded.

## Results

Completed 18/18 source cases and 144/144 finite rows; status counts: `{"ok": 144}`. Best finite-u point: D=16.6 Mpc, i=72°, Υ*=0.7620, gas_thin_hstar0p3, K/Kcrit=0.995; objective=864.164 (χ²diag=858.975, priors=5.189). Same-node controls: NFW=83.395 (optimizer_success=True), MOND=5101.728.

At this source node, the Kcrit edge-limit objective is 868.919; it is reported only as a separate limit for context. The finite score falls from 13,734.8 at K/Kcrit=0.5 to 864.164 at 0.995, then rises slightly to 864.887 at 0.998. At the finite best, u spans 15.03–394.48 over the half-domain; at 0.998 it spans 36.13–988.62. Thus the best sampled finite point is numerically positive but remains near the critical branch, with a large field amplitude. This screen does not find a distinct moderate-u optimum.

As a source/eigensolver check, the D=13.8 Mpc, i=72°, Υ*=0.762, gas-mix/h*=0.3 case gives Kcrit=5.3713855702×10⁻⁵ (km/s)⁻². The independently assembled checked-mode result is 5.3713874645×10⁻⁵ (km/s)⁻², a fractional difference of −3.53×10⁻⁷. Across all 18 source cases, the largest principal-mode residual is 1.29×10⁻¹¹; all 18 NFW optimizers report success. The maximum finite linear residual is 2.75×10⁻⁹.

Median/max u growth by finite K fraction (median, min, max of full-domain maximum u across successful rows):

- 0.500: 2.735, 2.697, 2.758
- 0.800: 8.472, 8.234, 8.592
- 0.900: 18.3, 17.68, 18.59
- 0.950: 38.08, 36.66, 38.73
- 0.980: 97.54, 93.69, 99.27
- 0.990: 196.7, 188.8, 200.2
- 0.995: 395, 378.9, 402.1
- 0.998: 989.8, 949.4, 1008

The plot profiles over the selected D/M-L nodes separately for each vertical case; its NFW/MOND lines are the minima of same-node direct controls over those same nodes. The overall finite best is at the sampled +2σ distance boundary and the upper endpoint of the selected M/L subset; its K fraction is 0.995 rather than the largest sampled 0.998. Inclination is fixed at 72°, so this run does not profile it.

## Interpretation and limits

K is an added profiled galaxy parameter. The comparisons are in-sample diagonal-score screens with no radial covariance, no K prior, and a finite nuisance subset; they do not establish model preference, significance, or detection. The saved nuisance envelope has no exact i=72° rows, so the NFW/MOND controls were recomputed directly using the same transformed maps and geometry. These are not holdout predictions. NGC 4559 holdout data were not read or modified.

## Runtime, hashes, and outputs

Measured CPU 86.48 s; wall 86.51 s; single-thread BLAS environment: `{'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1'}`. Code SHA-256: `434ae2a5c4d2e6a2094ba7f343d339f46eaf8fcde2209233ca8e0ab6891414e0`. Results SHA-256: `a7287d97d6f03351a3902ac1180beeab32cbaa79d60511346e63881bb7362076`. Input hashes and all scalar/curve results are in `results.json`; `SHA256SUMS.txt` hashes task outputs.

Outputs: `results.json`, `finite_profile.csv`, `finite_curves.csv`, and `finite_k_profile.png`.

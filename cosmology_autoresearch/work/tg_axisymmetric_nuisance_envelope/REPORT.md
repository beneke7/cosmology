# NGC 3198 TG source and geometry nuisance envelope

**Finding — independently checked source geometry; exploratory model envelope.** Reprojecting the raw gas and S4G maps with a common distance and inclination materially improves the conditional TG limiting curve relative to the earlier fixed-source score. At the best sampled nuisance point, the TG edge is close to simple MOND on the diagonal score, but remains well above the same-source NFW score. The TG edge minimum lands at the edge of the sampled distance/inclination grid, so it is not a continuous optimum under the Gaussian priors. No model-preference or significance claim is justified.

## Reproduce

Run from the `cosmology_autoresearch/` project directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python work/tg_axisymmetric_nuisance_envelope/nuisance_envelope.py --pilot

.venv/bin/python scripts/run_bounded.py --state RUN_STATE.json --seconds 1700 -- \
  env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python work/tg_axisymmetric_nuisance_envelope/nuisance_envelope.py --run --cpu-cap 850

.venv/bin/python scripts/run_bounded.py --state RUN_STATE.json --seconds 300 -- \
  env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python work/tg_axisymmetric_nuisance_envelope/verify_nfw.py --cpu-cap 120
```

The production run completed all 1,350 model-grid rows in 387.66 CPU seconds and 387.83 wall seconds, below the 850 CPU / 1,700 wall-second bounds. It used one CPU thread, no GPU, and NumPy 2.5.3 / SciPy 1.18.1 / Python 3.12.3. The separate NFW verification completed in 13.61 CPU seconds and 13.63 wall seconds under its 120 CPU / 300 wall-second cap. The row-by-row checkpoint is [`envelope_progress.jsonl`](envelope_progress.jsonl); scalar rows are also in [`envelope_grid.csv`](envelope_grid.csv), and the summary is [`results.json`](results.json).

## Geometry and shared source

The map-to-disk transformation is applied to the raw maps at every grid point. For sky offsets rotated to major/minor axes, the deprojected radius is

`R = D × 1000 × ARCSEC_RAD × sqrt(x_major² + (x_minor/cos i)²)`.

The H I column and stellar surface brightness are line-of-sight quantities. Their face-on surface densities are multiplied by `cos(i)`, while disk-plane pixel area is the projected physical pixel area divided by `cos(i)`. The factors cancel in integrated map mass; at fixed angular map, mass scales as `D²` and is independent of inclination. The radial distribution still changes with both `D` and `i`, so the disk force was recomputed from the transformed sources. The source vertical heights remain in kpc; changing distance therefore changes the ratio of height to radial scale.

SPARC Rotmod velocities and errors are tabulated deprojected at 73°. At each inclination node, both are multiplied by `sin(73°)/sin(i)`, and the Rotmod radii are multiplied by `D/13.8` before model evaluation. Thus the run does not apply a velocity-only scale to the map model.

The source is the same axisymmetric parametric source family used in the prior pilot: HALOGAS H I column multiplied by 1.36 for helium, plus the S4G P5 old-stellar map within its published ellipse. The measured P5 profile is kept only for complete supported annuli; its outer stellar source is continued by the declared exponential fit. Gas uses the existing 5-kpc taper beyond its last supported annulus. There is no H₂ or ionized component. The maps are azimuthally averaged, not a unique 3D baryon reconstruction.

The HI transformation holds a flat disk and fixed position angle while changing `i`. The NGC 3198 gas map and the prior source extraction are centered near `i≈72°`; excursions beyond about 2° from that extraction are not independently validated against a tilted-ring or warp model. The ±2σ nuisance grid extends as far as 67°–79°, and its TG minimum uses 67° (5° below the map-extraction value). Therefore that high-leverage point is conditional on a flat-disk deprojection that the present source audit does not establish at that inclination. The optical P5 ellipse and HI morphology need not share one exact inclination at all radii.

The original extraction geometry was reproduced before varying nuisances. At `D=13.8 Mpc, i=72°`, the new profile builder matches the frozen HALOGAS and THINGS radial profiles at all 70 gas radii and the 28 complete supported stellar rings, with maximum relative surface-density differences below `1.9×10⁻¹³`. The independent K=0 Hankel implementation agrees with the frozen Newtonian curve to `0.0103 km/s` maximum speed difference. On the same 320×480 TG grid, `D=13.8 Mpc, i=72°, Υ*=0.6`, the new `Kcrit` differs from the frozen value by `−1.47×10⁻⁵` fractionally; the edge curve differs by at most `0.00245 km/s`. Eigen residual is `5.1×10⁻¹²`. These checks validate the source transformation and the reused solver conventions.

At that i=72° geometry, correcting the observed velocities from their catalog 73° deprojection changes the frozen `Υ*=0.6` edge score from `9,820.1` (uncorrected catalog comparison) to `10,175.9` (i=72°-corrected). This is a velocity-convention change, separate from the source and eigensolver agreement.

## Priors, model definitions, and scores

The common nuisance treatment uses Gaussian penalties on the `−2 log prior` scale:

- `P_D = ((D−13.8)/1.4)²` with D in Mpc;
- `P_i = ((i−73)/3)²` with i in degrees;
- `P_Υ = (log10(Υ*/0.6)/0.1)²`.

Distance and inclination are sampled at `−2, −1, 0, +1, +2σ`. The stellar M/L grid spans `−1, −0.5, 0, +0.5, +1σ` in log space and also includes the quoted NGC 3198 reference `Υ*=0.762` (1.038σ above 0.6). Every model uses these same nuisance nodes and penalties. Each reported objective is `χ²diag + P_D + P_i + P_Υ`; the NFW halo has no concentration–mass relation or halo prior and has only its stated bounds (`10⁹≤M200/M☉≤10¹⁴`, `1≤c200≤50`). Simple MOND uses fixed `a0=1.2×10⁻¹⁰ m/s²`, `μ(x)=x/(1+x)`, and no external-field effect. H0 for NFW is fixed at 73 km/s/Mpc.

TG's main row is the principal-mode `Kcrit` edge. It is a limiting curve without a finite-normalization positive-u solution; `Kcrit` is derived from the source and no K prior is applied. At each vertical case's best edge nuisance point, positive finite-u curves were also solved separately at `K/Kcrit=0.95` and `0.99`. These fixed fractions are illustrative and have no K prior or K profile.

| Vertical source case | Model | Best sampled D, i, Υ* | χ²diag | P_D / P_i / P_Υ | χ²diag + priors |
|---|---|---:|---:|---:|---:|
| Gas 85/15, h* = 0.3 kpc | TG edge limit | 16.6 Mpc, 67°, 0.762 | 300.10 | 4 / 4 / 1.08 | 309.18 |
|  | NFW | 13.8 Mpc, 73°, 0.762 | 59.31 | 0 / 0 / 1.08 | 60.38 |
|  | Simple MOND | 12.4 Mpc, 67°, 0.762 | 295.72 | 1 / 4 / 1.08 | 301.80 |
| Gas thin-only, h* = 0.3 kpc | TG edge limit | 16.6 Mpc, 67°, 0.762 | 315.30 | 4 / 4 / 1.08 | 324.38 |
|  | NFW | 13.8 Mpc, 73°, 0.762 | 59.94 | 0 / 0 / 1.08 | 61.02 |
|  | Simple MOND | 12.4 Mpc, 67°, 0.762 | 307.03 | 1 / 4 / 1.08 | 313.11 |
| Gas 85/15, h* = 0.6 kpc | TG edge limit | 16.6 Mpc, 67°, 0.762 | 648.33 | 4 / 4 / 1.08 | 657.41 |
|  | NFW | 13.8 Mpc, 73°, 0.762 | 62.70 | 0 / 0 / 1.08 | 63.78 |
|  | Simple MOND | 12.4 Mpc, 67°, 0.762 | 515.57 | 1 / 4 / 1.08 | 521.65 |

For the baseline gas mixture and `h*=0.3 kpc`, the NFW profile has `M200=3.71×10¹¹ M☉`, `c200=10.96`. Its best nuisance point is the prior center. The TG edge and MOND choose the high-M/L node, but different distances: the TG sampled minimum selects `D=+2σ, i=−2σ`, while MOND selects `D=−1σ, i=−2σ`. The TG edge minimum hits both distance/inclination grid boundaries; expanding the grid or doing continuous refinement could change it. Treat its objective as the best point in this declared grid, not a fully profiled Gaussian-prior optimum.

The NFW optimizer was independently rerun for all 450 geometry/vertical/M-L nodes, with all three deterministic least-squares starts retained. All 1,350 starts reported convergence, every selected fit succeeded, and all saved NFW `χ²diag`, prior-penalized objectives, `log10(M200)`, and `log10(c200)` match the independent rerun exactly (maximum absolute difference 0). The per-start status, message, objective, and parameters are in [`nfw_verification.jsonl`](nfw_verification.jsonl); totals are in [`nfw_verification_summary.json`](nfw_verification_summary.json). The original [`results.json`](results.json) and [`envelope_progress.jsonl`](envelope_progress.jsonl) were read-only inputs to this audit.

### Map-supported inclination screen

The full ±2σ SPARC grid contains inclinations 67°, 70°, 73°, 76°, and 79°. The independent H I source audit recommends a nominal geometry near 71° and at least a ±2° check; more extreme flat-disk re-deprojections are not validated by the tilted-ring/warp source models. Filtering the *existing coarse grid* to 69°–73° therefore retains only its 70° and 73° nodes; this is a post-hoc robustness screen, not a new continuous profile or posterior. For the baseline 85/15 gas mixture and `h*=0.3 kpc`, the best retained nodes are:

| Model | Best retained D, i, Υ* | χ²diag | χ²diag + priors |
|---|---:|---:|---:|
| TG edge limit | 16.6 Mpc, 70°, 0.762 | 511.08 | 517.15 |
| NFW | 13.8 Mpc, 73°, 0.762 | 59.31 | 60.38 |
| Simple MOND | 12.4 Mpc, 70°, 0.762 | 490.45 | 493.52 |

The TG edge and simple MOND are near one another on this diagonal screen; both select the high-M/L node and low-inclination side. NFW remains descriptively much lower. None of these score gaps is a calibrated model-selection result. The `i=67°` TG edge minimum and all six finite-u examples in the next table lie outside this map-supported screen. A separate finite-`K` profile at the source-supported inclinations is being run before interpreting whether any regular TG solution remains competitive.

At the TG edge-selected nuisance point, the finite-u results are:

| Vertical source case | K/Kcrit | χ²diag | Prior terms D / i / Υ* | χ²diag + priors | Minimum u |
|---|---:|---:|---:|---:|---:|
| Gas 85/15, h* = 0.3 kpc | 0.95 | 514.36 | 4 / 4 / 1.08 | 523.44 | 2.24 |
|  | 0.99 | 207.61 | 4 / 4 / 1.08 | 216.69 | 7.35 |
| Gas thin-only, h* = 0.3 kpc | 0.95 | 488.95 | 4 / 4 / 1.08 | 498.03 | 2.24 |
|  | 0.99 | 212.74 | 4 / 4 / 1.08 | 221.82 | 7.34 |
| Gas 85/15, h* = 0.6 kpc | 0.95 | 1,176.09 | 4 / 4 / 1.08 | 1,185.17 | 2.33 |
|  | 0.99 | 650.45 | 4 / 4 / 1.08 | 659.53 | 7.84 |

All six finite examples have positive u throughout the finite-volume domain and linear-solve residuals below `5.3×10⁻¹⁰`. Their fixed K fractions are not a physical posterior profile. The limiting edge and finite-u values remain separate in [`results.json`](results.json), [`envelope_grid.csv`](envelope_grid.csv), and [`profiled_predictions.csv`](profiled_predictions.csv). The curves are plotted in [`profiled_models.png`](profiled_models.png).

## Interpretation and limits

The fixed-source TG edge score near 3,600–4,500 at high M/L is not stable to the coupled source geometry and velocity corrections. In the sampled ±2σ D/i box, the best `h*=0.3` edge score falls to 300.1, nearly the simple-MOND score of 295.7 after the respective prior terms (309.2 versus 301.8). Thus the earlier large conditional gap from MOND does not survive this sampled geometry envelope. The same comparison still has an NFW score near 60.4, substantially lower as a descriptive diagonal score; this experiment does not turn that gap into a statistical preference. Increasing the stellar scale height to 0.6 kpc worsens both TG and MOND scores substantially, while the NFW score changes only modestly. Gas thin-only versus the 85/15 mixture has small effects here.

The strongest result is a nuisance sensitivity: treating map inclination and the SPARC velocity deprojection consistently can move the TG edge by a large amount. The TG score minimum at the ±2σ corner means the sampled envelope is not closed under further nuisance profiling. No radial covariance is available, so all `χ²diag` values are descriptive pointwise-error summaries only. The NFW and MOND families also have different fitted flexibility, and the TG edge itself is not a finite solution.

The model remains an axisymmetric radial source surrogate built from imperfect source maps. Stellar light is the S4G P5 old-star component with a finite ellipse and an exponential outer continuation; gas excludes unmodeled molecular and ionized components. Non-axisymmetric structure, vertical density beyond the tested exponential forms, and the map-derived source systematic against SPARC's total-light disk are not resolved here. This work tests the conditional NGC 3198 setup only; it establishes neither a universal TG parameter nor a galaxy-wide prediction.

No conventional test suite was run. The checks here are numerical scientific checks: map-profile reconstruction, surface-density/area Jacobians, K=0 kernel closure, same-grid TG edge/eigenvalue reproduction, positive-u finite solves with residuals, and the all-node NFW optimizer audit.

## Evidence and hashes

Inputs and code versions are recorded in [`results.json`](results.json) and [`nfw_verification_summary.json`](nfw_verification_summary.json). The task files and plots are listed in [`SHA256SUMS.txt`](SHA256SUMS.txt); source-map and SPARC input hashes are embedded in the result records.

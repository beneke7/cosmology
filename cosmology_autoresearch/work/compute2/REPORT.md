# Radiation sensitivity profile scan

Status: exploratory reference-model calculation; independently checked against a separate scalar implementation. This is a background-only profile-χ² scan, not posterior inference, evidence, a significance test, or a reproduction of the DESI collaboration's physical model.

## Setup and reproducibility

The fit uses all 13 rows of the DESI DR2 mean vector and the supplied full covariance in unchanged row order. The expansion history is

`E²(z) = Ωm(1+z)³ + Ωr(1+z)⁴ + Ωde fDE(z)`, with `Ωde=1−Ωm−Ωr`.

Photons use `a_rad=8π⁵k_B⁴/(15h_P³c³)` and `ρcrit,100=3H100²/(8πG)`. With the supplied SI constants, `T_CMB=2.7255 K` and massless-neutrino reference `N_eff=3.046`, the calculation gives `ωγ=2.472975328714088×10⁻⁵`, `ωr=4.183702725849734×10⁻⁵`, and `Ωr(H0=70)=8.538168828264765×10⁻⁵`. The sampled H0 envelope gives `Ωr=5.1650650936416465×10⁻⁵` to `1.6734810903398935×10⁻⁴`.

The neutrino approximation is a reference background, not an exact DESI physical model with massive neutrinos. H0 is only used to map `ωr` to `Ωr`; it is a sensitivity grid, not a prior or an H0 measurement. The sound horizon is not predicted. The free amplitude `α=c/(H0 rd)` is analytically profiled within `[1e-6,1e4]`; shape bounds are imported unchanged from `background_bao.FIT_BOUNDS`.

Command:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python scripts/run_bounded.py --seconds 1200 -- \
  .venv/bin/python experiments/radiation_sensitivity/radiation_scan.py
```

Runtime was 6.835 s on Python 3.12.3, NumPy 2.5.3 and SciPy 1.18.1, with 20 CPUs in affinity and one numerical-library thread. No GPU was used. Arithmetic and optimization used float64.

Input SHA-256 values:

- `context/data/desi_dr2_mean.txt`: `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`
- `context/data/desi_dr2_cov.txt`: `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`
- `experiments/campaign_seed_bao/result.json`: `92662f3ee42500d2be34d9f87b485fe56c82e9550fa674e42f4dec2bf2cab1f0`
- `scripts/background_bao.py`: `34bdf1d81397b39b740faa3e2f4a79748d5fa9d3bfdd8303619d631dbf8aa2e2`
- `experiments/radiation_sensitivity/radiation_scan.py`: `b825ee41df8296c16bb332535acb2ae0b33eadf15b51227328d067292f109e21`

The machine-readable record at [result.json](../../experiments/radiation_sensitivity/result.json) includes the exact H0 grid, rowwise predictions and shifts, parameters, hashes, full optimizer starts/statuses, convergence checks, and run metadata. All three models are refitted at ten H0 values: 50,55,60,65,67.4,70,75,80,85,90 km/s/Mpc. The independent theory worker's separate 41-value scan holds the seed parameters fixed; it is not 41 radiation refits. An independent review of all saved optima is recorded in [saved_fit_review.json](../orchestrator2_radiation/saved_fit_review.json); [scalar_reference.json](../orchestrator2_radiation/scalar_reference.json) separately checks fixed seed parameters.

## Results

The independently optimized no-radiation fits reproduce the seed χ² values to better than `1e-11`; the Ωr→0 limit gives `E(0)=1`, exact agreement in `E²` with the starter background, and predictions within `3.6×10⁻¹⁵` of its implementation.

| Model | No-radiation χ² | Radiation-refit χ² over H0 grid | Refitted Δχ² vs no radiation |
|---|---:|---:|---:|
| flat ΛCDM | 10.27104100 | 10.27829435–10.29479686 | +0.00725335 to +0.02375586 |
| wCDM | 9.04104661 | 9.02534381–9.03619183 | −0.01570280 to −0.00485477 |
| CPL | 5.70061594 | 5.69818103–5.69985957 | −0.00243490 to −0.00075637 |

The ordering remains CPL, wCDM, ΛCDM at every sampled H0. These are best-fit comparisons with different model dimensions and selected parameter domains, not significance or evidence statements. The extrema above cover only the stated finite H0 grid and the corresponding baseline/refitted optima; they are not a bound over the full shape-parameter domain.

At H0=50, the fixed-baseline prediction shifts (holding the no-radiation best-fit shape and α fixed) are:

| Model | Fixed-point Δχ² | Largest `|Δprediction|/sqrt(Cii)` | `Δpredictionᵀ C⁻¹ Δprediction` |
|---|---:|---:|---:|
| flat ΛCDM | +0.03962229 | 0.07452707 | 0.01628546 |
| wCDM | −0.00066883 | 0.07307493 | 0.01502564 |
| CPL | +0.00842568 | 0.06236726 | 0.01087807 |

The largest fixed-point fractional row shift over the models is below `8.8×10⁻⁴`. After refitting, the largest marginal-sigma row shift is at most 0.0112; the largest full-covariance displacement quadratic is `3.62×10⁻⁴`. Radiation therefore causes small prediction changes for these fitted points on this data vector. This statement does not assert small changes over every allowed shape.

## Numerical and optimizer checks

- At the fitted points, 96-node vs adaptive quadrature changes predictions by at most `4.3×10⁻¹⁴`; 96-node vs 192-node changes them by at most `1.21×10⁻¹³`. The largest adaptive-quadrature χ² difference is `2.03×10⁻¹³`.
- The analytic profiled α agrees with a bounded scalar minimization to `2.3×10⁻¹³` at worst. Its bound is inactive at all reported fits.
- The bounded multi-start search retained 16 L-BFGS-B starts for each of 33 fits. Nineteen of 528 starts returned unsuccessful optimizer statuses; all are kept in the JSON. All 33 independent differential-evolution checks converged, with maximum χ² difference `3.10×10⁻⁹`. Eleven active-face Powell checks converged. The CPL optimum has `w0=-0.3` on the stated upper bound; feasible inward steps are recorded for the boundary checks.
- The independent scalar implementation checked all 33 saved optima; its largest prediction and χ² discrepancies were below `5.0×10⁻¹⁴` and `2.1×10⁻¹³`, respectively.

The fixed-baseline shifts are evaluated at the independently optimized no-radiation best-fit parameters in this scan. Their differences from calculations anchored on the full-precision saved seed parameter tuple arise from finite optimizer tolerances; the no-radiation χ² reproduction remains much tighter than the stated `1e-6` tolerance.

# Fixed-Υ thermodynamic-gravity axisymmetric pilot: Υ*=0.762

## Result

Ran the pre-existing `solve_tg_axisymmetric.py` without code changes on the declared 320×480 cylindrical grid (`dr=0.25 kpc`, `dz=0.125 kpc`, `Rmax=80 kpc`, `Zmax=60 kpc`), with only the paper-reference stellar mass-to-light ratio Υ*=0.762. All 11 requested K points completed with solver status `ok`. The lowest sampled diagonal score is at the upper edge of this scan, K=3.9×10⁻⁵ (km/s)⁻²: χ²diag=13,693.18 for 43 SPARC rotation-curve points and RMS speed residual 35.79 km/s. This is an exploratory descriptive score, not a likelihood or model preference; it is still a poor curve match, uses no covariance, and reaches the scan boundary. No fit or detection claim follows.

## Per-K outcomes

| K [(km/s)⁻²] | Status | χ²diag (43 points) | CG relative residual | min midplane u |
|---:|:---|---:|---:|---:|
| 0 | ok (Newtonian-limit evaluation) | 30,762.18 | 1.9592e-9 | 1.00000000 |
| 5.0e-6 | ok | 29,141.92 | 1.9956e-9 | 1.01051 |
| 1.0e-5 | ok | 27,403.45 | 1.9704e-9 | 1.02293 |
| 1.5e-5 | ok | 25,528.01 | 1.9620e-9 | 1.03797 |
| 2.0e-5 | ok | 23,492.33 | 1.9877e-9 | 1.05676 |
| 2.5e-5 | ok | 21,267.38 | 1.9877e-9 | 1.08120 |
| 3.0e-5 | ok | 18,816.92 | 1.9910e-9 | 1.11478 |
| 3.3e-5 | ok | 17,220.31 | 1.9516e-9 | 1.14203 |
| 3.5e-5 | ok | 16,096.73 | 1.9663e-9 | 1.16472 |
| 3.7e-5 | ok | 14,922.04 | 1.9572e-9 | 1.19254 |
| 3.9e-5 | ok; best sampled | 13,693.18 | 2.0089e-9 | 1.22757 |

The score improves monotonically over this sampled K range; therefore the best sampled point does not locate an interior optimum or a robust preferred K.

## Numerical and source checks

- Main scan measured 23.849 s process CPU and 23.872 s wall. A separate one-point full-field positivity/residual check at the best K took 2.61 s user CPU / 2.68 s wall; total measured compute was about 26.5 CPU seconds, well below the 300 CPU-second ceiling and 1800 s wall bound. Numerical-library threads were each set to one. No GPU was used.
- At the best K, a repeated solve of the same assembled finite-volume operator gave CG info=0, relative residual 2.0089e-9, full-grid min(u)=1.180391 and max(u)=5.847963 on all 153,600 cells; zero cells had u≤0. Midplane range in the production output is 1.227574–5.858022. This verifies positivity only for this discretized candidate/box, not all admissible solutions or continuums.
- Source profile SHA-256: `2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5`.
- Other fixed source inputs: K=0 reference CSV SHA-256 `5090562be95b6b757cbb9cb7e623d1e42510f7319a9b35e1ddc13e15fa849eff`; SPARC Rotmod archive SHA-256 `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`.
- Source construction is S4G stars through the last fully covered 13.75-kpc annulus, with an exponential tail (fitted scale 4.0553 kpc; log-profile RMS 0.0552), plus HALOGAS HI+He with the declared thin/thick vertical mixture; stellar scale height is 0.3 kpc. Outer boundaries use the solver's monopole Robin approximation, with midplane reflection and zero radial flux at the axis.

## Interpretation, caveats, and issues

Υ*=0.762 is included because it is quoted in the paper, as a reference value only—not an independent prior or a preferred mass calibration. At fixed source assumptions, the tested positive K values improve the diagonal curve score, but even the best is extremely large and has a 35.8 km/s RMS residual. This is not a significant positive physical result. It also is not a global negative result: this scan does not establish the first singular/admissibility edge, vary source geometry or boundary approximation, include full observational/systematic covariance, or compare matched Newtonian+halo and MOND alternatives. In particular, the optimum-at-scan-edge behavior is a warning against interpreting the score trend as parameter evidence.

No solver failures or other execution issues occurred. The score uses the 43 tabulated point errors diagonally; the reported χ² is not covariance-aware. Stability under resolution/domain changes has not been tested in this fixed-Υ task.

## Reproduction

From the repository root, the bounded invocation was:

```sh
.venv/bin/python scripts/run_bounded.py --state work/tg_axisymmetric_pilot/scan_u0762/run_state.json --seconds 1800 -- prlimit --cpu=300 -- env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python work/tg_axisymmetric_pilot/solve_tg_axisymmetric.py --dr 0.25 --dz 0.125 --rmax 80 --zmax 60 --cpu-cap 280 --upsilon 0.762 --output-dir work/tg_axisymmetric_pilot/scan_u0762
```

The production summary, per-K records, best-candidate curve, and plot are saved alongside this report. The full-field positivity check repeated the best-K operator assembly and conjugate-gradient solve using the same implementation; it was a field-coverage check, not an independent solver.

# NGC 3198 same resolved-source alternatives

Completed 2026-09-25. This is a bounded descriptive fit to the 43 SPARC points, using the frozen HALOGAS axisymmetric map-derived K=0 baryon components for all alternatives. Errors are diagonal `errV`; there is no full covariance. Scores are not a statistical preference or detection claim.

## Method and provenance

Reproduction entry point: [`compare_models.py`](compare_models.py). Machine-readable fit results and 43-row predictions are in [`results.json`](results.json) and [`predictions.csv`](predictions.csv); the comparison plot is [`map_matched_comparison.png`](map_matched_comparison.png). It reads exactly 43 rows from `context/data/Rotmod_LTG.zip::NGC3198_rotmod.dat` and matches every radius, observed speed, and error against the frozen K=0 and branch-edge files before fitting.

Input SHA-256 values:

- SPARC Rotmod archive: `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`
- [`newtonian_k0_predictions.csv`](../../../tg_axisymmetric_pilot/newtonian_k0_predictions.csv): `5090562be95b6b757cbb9cb7e623d1e42510f7319a9b35e1ddc13e15fa849eff`
- [`axisymmetric_branch_edge_curves.csv`](../../../tg_axisymmetric_pilot/axisymmetric_branch_edge_curves.csv), corrected signed-gas output: `3a7ced618b539180d0f50b4ab4fd6f26f30e234bdb247d97fb6de43b31c568dd`
- [`axisymmetric_branch_edge_summary.json`](../../../tg_axisymmetric_pilot/axisymmetric_branch_edge_summary.json): `566f95a36fe993f9c02ebba6a6c8faf4bac69f3cc1c2ed94126395349ce596b2`
- This reproduction script: `00d7c63e1b69641986bf8da1f841805e410be8e01288011341e383c56657b688`

At reference Υ*=0.6, the frozen map table gives stellar-only speed and total baryonic speed. I recovered the *signed* gas term as `Vgas_signed² = Vbar_map² − Vstar_map²`, then used `g_N(r;Υ*) = [(Υ*/0.6)Vstar_map² + Vgas_signed²]/r`. The recovered gas acceleration term is negative at 16/43 radii (minimum −25.56 `(km/s)²/kpc`); none were replaced with absolute values. The independent axisymmetric K=0 total at Υ*=0.762 agrees with rescaling the frozen Hankel K=0 components to floating-point precision (maximum difference 0.000 km/s) after correcting the branch-edge CSV to preserve signed gas V². Its earlier gas-speed squaring had produced a superseded 0.583 km/s mismatch. The fits use the signed decomposition above.

The common stellar policy is a continuous bounded profile over Υ* ∈ [0.6×10⁻⁰·¹, 0.762] = [0.47659694, 0.762] M☉/L☉ at 3.6 μm, with no extra M/L likelihood penalty. Also reported are fixed Υ*=0.6 and 0.762 results. Newtonian baryons have one free parameter (Υ*); isolated simple-MOND has one (Υ*, with fixed `a0=1.2×10⁻¹⁰ m/s²`); Newtonian+NFW has three (Υ*, `M200`, `c200`). NFW uses log-box bounds `9≤log10(M200/M☉)≤14` and `1≤c200≤50`, `H0=73 km/s/Mpc`, and no concentration–mass relation. Fixed-Υ* NFW sensitivity has two halo parameters.

Equations used (r in kpc, speed in km/s, acceleration in `(km/s)²/kpc`):

- Newtonian baryons: `V² = r g_N`.
- NFW: `V² = r g_N + G M200 f(r/rs)/(r f(c200))`, with `f(x)=ln(1+x)−x/(1+x)`, `r200=[3M200/(800πρcrit)]^(1/3)`, `rs=r200/c200`, and `ρcrit=3H0²/(8πG)`.
- Simple isolated MOND: `μ(x)=x/(1+x)`, `x=g/a0`; hence `g=[g_N+sqrt(g_N²+4a0g_N)]/2` and `V²=rg`. The minimum net inward `g_N` over every radius and the whole M/L interval is positive: `71.50 (km/s)²/kpc` at the lower M/L bound. The algebraic MOND relation is therefore defined for all fitted points and allowed M/L values. This prescription has no external-field effect.

The best branch-edge curve is the `Υ*=0.762` TG principal-mode edge (lowest diagonal score among four fine-grid M/L nodes). Per the branch summary it is a limiting curve, not a finite-normalization positive-u solution. `K` is fixed to `Kcrit` rather than optimized; the best M/L is selected from the declared grid. This is a source-constrained envelope, not a conventional fit or an equivalent parameter-count comparison.

## Descriptive results

`χ²diag = Σ[(Vobs−Vmodel)/errV]²` over the same 43 radii; RMS is unweighted velocity RMS.

| Curve | Parameter treatment | Υ* | Other fitted values | χ²diag | RMS (km/s) |
|---|---:|---:|---|---:|---:|
| Resolved-map Newtonian baryons | Υ* profiled (1) | 0.762 (upper bound) | — | 30397.96 | 53.342 |
| Newtonian + NFW | Υ*, M200, c200 profiled (3) | 0.762 (upper bound) | M200=3.812×10¹¹ M☉, c200=10.493 | 58.49 | 7.683 |
| Isolated simple-MOND | Υ* profiled (1), fixed a0 | 0.656626 | no EFE | 1032.23 | 10.589 |
| TG branch edge limit | Υ* selected from 4-node grid; K=Kcrit limit | 0.762 | not a finite-u solution | 3585.85 | 17.203 |

At fixed Υ*=0.6: baryons-only χ²=39286.30, RMS=60.152 km/s; MOND χ²=1131.90, RMS=10.433 km/s; NFW with two halo parameters has M200=3.638×10¹¹ M☉, c200=12.599, χ²=74.57, RMS=8.440 km/s. At fixed Υ*=0.762, NFW is χ²=58.49, RMS=7.683 km/s (halo fit above), and MOND is χ²=1334.74, RMS=12.665 km/s.

## Runtime and limits

Run command used single-thread settings for OpenBLAS, OMP, MKL and NumExpr; no GPU. The rerun against the corrected branch-edge input took 0.46 CPU seconds and 0.50 wall seconds including Python startup and plot writing (`time.log`). The calculation's internal timed section took 0.156 CPU seconds and 0.156 wall seconds (`results.json`), below the 200 CPU-second and 1200 wall-second caps.

The result is conditional on this pilot's axisymmetric source reconstruction and map force calculation; SPARC force components do not uniquely determine a 3D density. It uses no covariance beyond the pointwise diagonal errors. The NFW and MOND scores quantify only these specified families and boxes. No inference of a detection, physical model preference, or equivalence between the TG limiting curve and a finite TG solution is made.

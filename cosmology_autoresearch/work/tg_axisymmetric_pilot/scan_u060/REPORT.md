# Bounded NGC 3198 axisymmetric pilot: Υ*=0.6

**Classification: exploratory numerical branch, not a detection or model comparison.**

## Run and provenance

Ran the existing `solve_tg_axisymmetric.py` unchanged at `dr=0.25 kpc`, `dz=0.125 kpc`, `Rmax=80 kpc`, `Zmax=60 kpc`; `Υ*=0.6`; declared K grid `0` through `3.9e-5 (km/s)^-2` (11 points). Threads were limited to one; no GPU was used. `run_bounded.py` enforced an 1800 s wall ceiling; solver CPU cap was 240 s, below the 300 s task ceiling. `/usr/bin/time -v` measured **39.90 s wall, 39.90 s CPU** (solver-reported 39.44 s each), peak RSS 217 MB, exit 0.

- Source profile SHA-256: `2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5`
- Solver SHA-256: `db4caf73ae80ecb680559607fbeb7d41c7fcb040426009b3b77fa3b8de9161fe`
- SPARC Rotmod archive SHA-256: `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`
- All candidates, including residuals and scores: `axisymmetric_candidate_points.jsonl`
- Run summary: `axisymmetric_pilot_summary.json`

The model source is an axisymmetrized, parametric disk: audited S4G stellar profile inside complete rings then an exponential tail, plus HALOGAS H I+He split into thin/thick vertical components. It is not a direct 3D mass reconstruction. The boundary uses the declared distant-box monopole Robin approximation.

## Numerical outcomes

All 11 grid entries returned `status=ok`, with CG relative residuals from `1.9666e-9` to `1.9977e-9` (requested `2e-9`). Interpolated midplane `u` was positive at every candidate and SPARC radius; its minimum ranged from `1.0000000017` (near-zero proxy) to `1.13750` (K=`3.9e-5`). The script does **not** retain the full 2D field, so this run does not independently establish positivity everywhere in the computational domain.

The minimum sampled descriptive diagonal score was **χ²=26479.80 for 43 SPARC points at the top edge K=`3.9e-5`**, with velocity RMS residual 49.34 km/s. Scores decrease monotonically over the tested K grid. Thus this scan has not located an interior best fit or an admissibility edge; the edge minimum is only a reason for a separately bounded follow-up if the root reviewer agrees. The score uses pointwise errors and no covariance, and must not be interpreted as evidence for TG or compared as if nuisance/prior choices were matched.

| K [ (km/s)^-2 ] | status | diagonal χ² (43 points) | relative residual | minimum midplane u |
|---:|:---|---:|---:|---:|
| 0 | ok* | 39674.70 | 1.9977e-9 | 1.0000000017 |
| 0.5e-5 | ok | 38346.38 | 1.9746e-9 | 1.00897 |
| 1.0e-5 | ok | 36937.71 | 1.9840e-9 | 1.01921 |
| 1.5e-5 | ok | 35438.19 | 1.9666e-9 | 1.03107 |
| 2.0e-5 | ok | 33835.16 | 1.9715e-9 | 1.04506 |
| 2.5e-5 | ok | 32113.21 | 1.9886e-9 | 1.06192 |
| 3.0e-5 | ok | 30253.44 | 1.9744e-9 | 1.08281 |
| 3.3e-5 | ok | 29061.87 | 1.9726e-9 | 1.09802 |
| 3.5e-5 | ok | 28232.34 | 1.9968e-9 | 1.10962 |
| 3.7e-5 | ok | 27372.38 | 1.9934e-9 | 1.12268 |
| 3.9e-5 | ok | 26479.80 | 1.9923e-9 | 1.13750 |

`*` **K=0 caveat:** the script evaluates K=0 with `K=1e-12` and then differentiates `u=1+K w`; this is cancellation-prone. Its χ²=39674.70 differs from the separately computed Hankel K=0 baseline (χ²=39286.3, same nominal resolved baryonic source) by 388.4. Treat this scan's K=0 value as a numerical proxy, not an independent Poisson reproduction. The positive-K results do not use this limit substitution.

## Interpretation / next check

No positive physical detection is established. Even the best sampled branch leaves a very poor descriptive residual, and the finite K grid ends while scores are still falling. Before interpreting that trend physically, prioritize a task-owned independent full-field `u>0` check and a tighter upper-K/admissibility scan with source/domain/grid convergence; independently verify the midplane force extraction. Keep K=0 controlled by the existing Hankel baseline, not the cancellation-prone proxy. Any later halo/MOND comparison must use matched source, radii, errors/covariance, and nuisance treatment.

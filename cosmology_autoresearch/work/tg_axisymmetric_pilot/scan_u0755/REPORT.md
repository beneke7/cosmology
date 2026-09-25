# Fixed stellar M/L axisymmetric TG scan: Υ*=0.7546788

## Result

The assigned bounded scan completed normally: 11/11 requested K nodes returned
`ok`, using 35.018 CPU seconds and 35.025 wall seconds, with BLAS/OpenMP
thread counts set to one. This is well below the 300 CPU-second task ceiling
and 1,800-second wall ceiling. No GPU was used; the sparse CPU run completed
quickly.

At this fixed stellar mass-to-light ratio (0.6×10^0.1 = 0.7546788), the best
sampled descriptive diagonal score is χ²=14,213.41 at the upper edge of the
tested interval, K=3.9×10⁻⁵ (km/s)⁻². The score declines monotonically over the
sampled K grid, so this run finds no interior best fit and does not locate the
admissibility/singularity boundary. It is not evidence of model preference or
a detection. The score uses 43 pointwise SPARC errors without a covariance
matrix and is only a descriptive statistic.

| K [(km/s)⁻²] | Status | Diagonal χ² (43 points) | RMS residual [km/s] | CG iterations | Relative linear residual | Minimum sampled midplane u |
|---:|:---|---:|---:|---:|---:|---:|
| 0 | ok | 31,116.78 | 53.941 | 2708 | 1.96e-9 | 1.000000002 |
| 5.0e-6 | ok | 29,508.56 | 52.555 | 2711 | 1.99e-9 | 1.01044 |
| 1.0e-5 | ok | 27,783.82 | 51.022 | 2714 | 1.96e-9 | 1.02276 |
| 1.5e-5 | ok | 25,924.20 | 49.309 | 2718 | 1.95e-9 | 1.03764 |
| 2.0e-5 | ok | 23,906.96 | 47.371 | 2722 | 1.98e-9 | 1.05618 |
| 2.5e-5 | ok | 21,703.69 | 45.147 | 2728 | 1.99e-9 | 1.08020 |
| 3.0e-5 | ok | 19,278.95 | 42.548 | 2749 | 1.99e-9 | 1.11302 |
| 3.3e-5 | ok | 17,700.01 | 40.755 | 2755 | 1.92e-9 | 1.13950 |
| 3.5e-5 | ok | 16,589.21 | 39.439 | 2758 | 1.93e-9 | 1.16144 |
| 3.7e-5 | ok | 15,428.07 | 38.009 | 2761 | 1.98e-9 | 1.18819 |
| 3.9e-5 | ok | 14,213.41 | 36.447 | 2770 | 1.99e-9 | 1.22165 |

The K=0 row in this solver is evaluated through its 10⁻¹² small-K limit. Its
predicted curve differs from the separately implemented Hankel K=0 reference,
after rescaling the reference stellar component to this Υ*, by 0.489 km/s RMS
(maximum absolute difference 1.359 km/s). The corresponding diagonal scores
are 31,116.8 and 30,741.8. The speed-level agreement is fairly close, but the
score difference is amplified by the data uncertainties; the two baselines
should not be silently treated as identical.

## Numerical and physical scope

- Finite-volume cylindrical quadrant grid: ΔR=0.25 kpc, Δz=0.125 kpc,
  Rmax=80 kpc, Zmax=60 kpc (320×480 cells). Axis regularity and midplane
  reflection are imposed; the outer box uses a monopole Robin approximation.
- The reported CG residual is the relative residual of the assembled linear
  system, below 2.0e-9 for every K. All 11 points returned finite positive
  speeds and positive sampled midplane u. The candidate summary only checks
  u positivity on the interpolated midplane, not throughout the full 2D field;
  full-field positivity is therefore not established by this run.
- Source is the audited axisymmetrized parametric profile: S4G stellar surface
  density through the last complete ring, an exponential stellar tail
  (Rd=4.05533 kpc), and HALOGAS H I plus helium with the declared thin/thick
  vertical mixture. This is not a direct 3D source-map reconstruction.
- A fixed Υ* branch cannot assess the stated stellar M/L uncertainty. The
  single highest-K score is a grid-edge result, not a fitted optimum. No
  source, nuisance, covariance, halo, or MOND alternatives were fit here.

## Reproducibility and provenance

Run from the project root:

```sh
.venv/bin/python scripts/run_bounded.py \
  --state work/tg_axisymmetric_pilot/scan_u0755/run_state.json \
  --seconds 1800 -- env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python work/tg_axisymmetric_pilot/solve_tg_axisymmetric.py \
  --dr 0.25 --dz 0.125 --rmax 80 --zmax 60 --cpu-cap 300 \
  --upsilon 0.7546788 \
  --output-dir work/tg_axisymmetric_pilot/scan_u0755
```

Input `source_profiles.csv` SHA-256:
`2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5`.
K=0 comparison CSV SHA-256:
`5090562be95b6b757cbb9cb7e623d1e42510f7319a9b35e1ddc13e15fa849eff`.
SPARC Rotmod archive SHA-256:
`0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`.

Machine-readable per-K results are in `axisymmetric_candidate_points.jsonl`;
the solver summary, best-curve CSV, plot, and selected candidate JSON are in
this same task-owned directory. No shared code, ledger, report, or campaign
state was changed.

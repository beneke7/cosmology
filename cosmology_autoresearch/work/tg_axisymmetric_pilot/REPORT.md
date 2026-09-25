# Thermodynamic gravity: NGC 3198 physical-evidence checkpoint

**Checkpoint:** 2026-09-25. **Finding:** reproduced calculation; independently
checked numerics; conditional negative for one declared source/boundary family.
No detection, formal model selection, or TG-wide falsification is established.

## The two papers are different theories

Both full PDFs are local and hash-pinned:

| Branch | Paper | Local full PDF | SHA-256 |
|---|---|---|---|
| Galaxy, 2024 | Pszota & Ván, [“Field equation of thermodynamic gravity and galactic rotational curves”](https://arxiv.org/abs/2306.01825), arXiv:2306.01825v3; *Physics of the Dark Universe* 46 (2024) 101660 | [`context/user_provided/pszota_van_2024.pdf`](../../context/user_provided/pszota_van_2024.pdf) | `55cbfbbb0a916ea2ce0b70523334ee301897695860fb4bdb5bfd2278ab47b62e` |
| Short-range, 2026 | [“Ghost-free higher-gradient Newtonian gravity from the Second Law of Thermodynamics”](https://arxiv.org/abs/2609.00317), arXiv:2609.00317v1 | [`context/papers/tg_higher_gradient_2026.pdf`](../../context/papers/tg_higher_gradient_2026.pdf) | `c0992fc15699bcd2c7abd5a4b044ace25203c2793c0cca42cb775888945fe165` |

The 2024 paper's stationary galaxy equation is
`Δφ − K|∇φ|² = 4πGρ`; for fixed density, `u=exp(−Kφ)` gives the exact linear
equation `(Δ + 4πGKρ)u=0`, with the isolated boundary `u→1`. This is a
fixed-source field calculation, not coupled matter evolution or a relativistic
cosmology. The 2026 paper is a separate linear higher-gradient/relaxation
model with two length scales and a correlated short-range potential. Its lab
parameters and experimental implications are not tests of the 2024 constant-K
galaxy equation. The local theory-critical report checks dimensions, operator
algebra and the coupled Jeans condition, and records caveats about neutral
zero modes and boundary-dependent entropy claims:
[`work/tg_thermo_theory_critical/REPORT.md`](../tg_thermo_theory_critical/REPORT.md).

The 2026 paper itself states `ℓ1 ≲ 4×10⁻⁵ m` from the UW single-Yukawa
`|α|=1` contour and its regularity sum rule. Our independent algebra audit
finds that the quoted `λ<38.6 μm` contour directly yields this length bound
only in the single-range limit `ℓ2=0`. In the real two-root case,
`ℓ1²=λa²+λb²` with linked amplitudes; the contour is not by itself a joint
two-range limit. This is an existing experimental constraint in a restricted
limit, not a detection, and we have not independently reproduced the full
torque likelihood.
The revised short-range review ([arXiv:2605.18212v2](https://arxiv.org/abs/2605.18212v2))
is also local at [`work/tg_short_range_literature/REPORT.md`](../tg_short_range_literature/REPORT.md).
Its displayed Yukawa limits are restricted to attractive `α>0`, so those plots
cannot be transferred directly to the 2026 TG model's `α<0`/correlated
two-range sector. The sign-specific 2020 UW crossing is a contour read-off,
not a newly calibrated limit; the model-specific torsion-torque re-fit remains
unfinished.

## NGC 3198 result

The source is built from audited THINGS and HALOGAS H I maps plus S4G P5 old
stars. Main-disk THINGS/HALOGAS H I masses agree at the 1.7% level. S4G gives
`1.6318e10 M☉` inside its official analysis ellipse for `Υ*=0.6`; only fully
covered stellar annuli through 13.75 kpc are measured and the outer tail is
explicitly extrapolated. The map-derived stellar force still differs from the
SPARC Rotmod stellar component by 18.39 km/s RMS over the fully covered
support, so the source is a parametric axisymmetrization, not a unique 3-D mass
map. The complete source conversion and hashes are in
[`SOURCE_AUDIT.md`](SOURCE_AUDIT.md).

The analytic Freeman exponential-disk check and k-spectrum refinement pass at
`2.12e-4` relative error in `V²` and `2.11e-4` fractional speed, respectively.
The finite-volume TG solution was run at fixed `Υ*=0.4766, 0.6, 0.7547, 0.762`
on a `320×480` cylindrical grid (`dR=0.25`, `dz=0.125`, `Rmax=80`,
`Zmax=60 kpc`), with an observed-curve-independent monopole Robin exterior.
The fixed-source transformed solve avoids cancellation in the `K→0` force.

At `Υ*=0.762`, the first positive-solution edge is
`Kcrit=5.3714678825e-5 (km/s)^−2`. A separately assembled conductance-graph
solver gives `5.3713874645e-5` (fractional difference `1.50e-5`) with eigenpair
residual `4.53e-12`; the production residual is `7.2e-12` and the principal
mode is positive on the full grid. The finite-`u` curve approaches a limiting
eigenmode curve, but that limit itself is **not** a finite-normalization
solution. Its best tested descriptive score is `χ²diag=3585.85`, `17.20 km/s`
RMS on 43 SPARC points. At `0.998 Kcrit`, the field remains positive but
`u_max≈987`; the finite curve has `χ²diag=3641.76` and `17.36 km/s` RMS.
Thus the improving score toward a pole is not evidence of a regular physical
fit.

The tested coarse/fine shift in `Kcrit` is 1.85%; extending the box at the
coarse spacing changes it negligibly. Replacing HALOGAS by THINGS shifts the
edge curve by 0.12 km/s RMS, and ±25% stellar-tail scale changes it by
0.28–0.31 km/s RMS. Stellar `M/L` is the much larger tested source nuisance.
These checks do not yet cover distance/inclination or vertical-profile
uncertainty.

## Same-source alternatives

Each model uses the same 43 radii, observed velocities and pointwise diagonal
`errV` values, with signed gas acceleration preserved. These are descriptive
scores only: no full radial covariance is released, so neither differences nor
their scale are model-selection significance.

| Model | Parameter treatment | `χ²diag` | Velocity RMS |
|---|---|---:|---:|
| Newtonian baryons + NFW | `Υ*, M200, c200` profiled (3); no c–M relation | 58.49 | 7.68 km/s |
| Isolated algebraic simple MOND | `Υ*` profiled (1), fixed `a0`; no EFE; approximate disk control | 1032.23 | 10.59 km/s |
| 2024 TG branch-edge limit | `Υ*` selected from four declared grid nodes; `K=Kcrit` | 3585.85 | 17.20 km/s |

The TG edge is a singular limit, unlike either finite alternative. NFW also has
more fitted parameters. This is not a formal evidence comparison, but the gap
is large enough to say that the tested TG branch is not competitive as an
NGC 3198 curve under this source family and screen. The lower-M/L nodes are
worse. The algebraic MOND control is not an exact AQUAL/QUMOND solve.

**Interpretation:** a meaningful conditional negative for the isolated
positive-`K`, fixed-source NGC 3198 construction, not a TG-wide falsification.
There is no physical detection. The limiting branch, missing covariance,
map-vs-Rotmod stellar discrepancy, untested distance/inclination/vertical
source variants, and absence of a second-galaxy prediction all constrain the
claim.

## Reproduction and provenance

Run from the repository root in the project environment; one numerical thread
per process was used. The exact scan command and runtime are preserved in
[`scan_u0762/REPORT.md`](scan_u0762/REPORT.md). The principal-mode envelope
entry point is `work/tg_axisymmetric_pilot/analyze_branch_edge.py`; it writes
to this task directory by default, so select a new output directory to retain
the pinned outputs:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python work/tg_axisymmetric_pilot/analyze_branch_edge.py \
  --cpu-cap 240 --output-dir /tmp/tg-edge-reproduction
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python work/tg_alternatives_contract/map_matched/compare_models.py
```

The independent principal-mode implementation and numerical checks are in
[`work/tg_axisymmetric_independent_review/PRINCIPAL_MODE_REPORT.md`](../tg_axisymmetric_independent_review/PRINCIPAL_MODE_REPORT.md).
Machine-readable branch results and curve CSV, plot, scans, source maps,
independent result, and same-source comparison inputs/outputs remain beside
these reports. Key frozen hashes:

- Source profiles: `2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5`.
- Production TG solver: `ee8f24d57954b304bdd524f3f62fd27f12a6f95e1e7ddc4b1e293911fa25f7bc`.
- Branch-edge script / result: `37118cc8efbe692bb09af9c11ddc3fef0bfaabc208aa103a11244ed5b7893ee8` /
  `566f95a36fe993f9c02ebba6a6c8faf4bac69f3cc1c2ed94126395349ce596b2`.
- Independent eigenmode code / result:
  `18eb3766ce69b2b32b3bd6cfc4de9f85d0510544a433982910a7bda19a126371` /
  `b013f8887c7ca12bb440f40ee6cb9e527014d103978ce417047a79cfc337b158`.
- Same-source comparator / result:
  `00d7c63e1b69641986bf8da1f841805e410be8e01288011341e383c56657b688` /
  `fa4b2cd8f1e628413a86e06bb0899c8db7c4d4c8068bbbb46d22ab3555491b21`.
- SPARC Rotmod archive: `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`.

The complete thermodynamics tranche used an estimated 280 CPU seconds across
the bounded scans, eigenmode and independent audits, sensitivity calculations
and matched-model screen, below its 1,200 CPU-second task ceiling. The GPU was
not used: an end-to-end sparse-solve benchmark was slower on GPU (3.7 s) than
CPU (2.1 s). No performance claim depends on using the accelerator.

## Ranked next tests

1. **Finish the physical nuisance envelope:** propagate distance/inclination,
   vertical thickness, stellar `M/L` and measured/extrapolated light separately
   through the same TG, NFW and MOND source pipelines. Freeze the prior/range
   before looking at each score. This determines whether the present conditional
   negative survives realistic source geometry; do not use the rotation curve
   to retune the photometric source.
2. **Make an out-of-sample galaxy prediction:** predeclare a shared-`K` or
   mass-scaling law from theory/units, then hold out a whole galaxy and compare
   its full curve against the same source-aware alternatives. A refit of NGC
   3198 cannot be a detection.
3. **Close the 2026 laboratory branch:** use the published torsion apparatus
   geometry and torque data to fit the model's correlated two-range potential,
   profiling measured nuisance parameters and the actual covariance. Do not
   translate an attractive-only Yukawa exclusion into the repulsive/correlated
   TG sector. This is a direct null-test opportunity, not a claimed signal.

The Warden's advisory may reorder these tests, but it does not replace the
underlying physical or numerical checks.

# Independent audit of Gate 2

**Disposition: independently checked; no unresolved sign, frame, unit, or implementation defect found in the saved six-case result.** The agreement applies to the declared fixed-metric, two-component diagnostic, not a complete cosmological perturbation calculation.

## Source equation and transformation

I read and visually checked the primary PDF, page 4, displayed Eqs. (19)–(20). PDF SHA-256: `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`.

The paper's IVS equations, with prime denoting conformal time, are:

```text
delta_c' = -(theta_c + h'/2) + a*Gamma*(rho_x/rho_c)*delta_c
theta_c' = -Hcal*theta_c
```

Here `Hcal = a'/a = a H`. Set `n=ln(a)`, `E=H/H0`, and `g=Gamma/H0`. Since `dn/dtau=Hcal`, division by `Hcal` gives:

```text
S = h'/(2*Hcal) = (1/2)*dh/dn
d(delta_c)/dn = -theta_c/Hcal - S + (g/E)*(rho_x/rho_c)*delta_c
d(theta_c/H0)/dn = -(theta_c/H0)
```

The rate conversion is exact dimensionally: `a*Gamma/Hcal = a*Gamma/(a*H0*E) = g/E`. Numerical evaluation across the 2,049-point grid differed by at most `1.11e-16` absolute for nonzero rates and exactly zero for `g=0`.

The paper's perturbation setup specifies synchronous gauge, zero anisotropic stress, and no momentum transfer in the DM rest frame. In this geodesic-CDM frame, setting `theta_c=0` is consistent: Eq. (20) is homogeneous, so zero remains zero. The separate nonzero-theta run is correctly labeled a diagnostic rather than a CDM-comoving solution. The prescribed metric `h(n)=2*S*(n-n_initial)` has `h'=2*S*Hcal`, so the metric source is exactly the stated constant `S=1e-3`.

## Absolute-CDM continuity check

The homogeneous CDM equation in `n` is `rho_c,n = -3*rho_c - (g/E)*rho_x`. With `U=rho_c*delta_c`, the product rule and the density-contrast equation give:

```text
U_n = [-3*rho_c - (g/E)*rho_x]*delta_c
      + rho_c*[-theta_c/Hcal - S + (g/E)*(rho_x/rho_c)*delta_c]
    = -3*U - rho_c*(theta_c/Hcal + S)
```

The interaction terms cancel, yielding `U_n + 3*U + rho_c*(theta_c/Hcal + S)=0`. I checked this pointwise from the independently computed RHS, then separately integrated the comoving quantity `a^3 U`, which obeys `(a^3 U)_n = -a^3*rho_c*(theta_c/Hcal+S)`.

## Independent numerical comparison

The implementation [independent_gate2_audit.py](independent_gate2_audit.py) recomputes the pinned flat massless-neutrino background and all six `(g, theta test)` cases. It uses SciPy RK45 at `rtol=1e-13`, `atol=1e-15`, `max_step=0.01`, different from Gate 2's DOP853 solver. Eq. (20) is evaluated from its analytic solution for the independent perturbation run. The script reads the saved Gate 2 JSON only to compare outputs; it does not import the Gate 2 code.

Across all six cases:

- Independent `delta_c` minima differ from the saved values by at most `1.72e-14` absolute; all maxima are `0.1` and agree exactly at saved precision.
- Sampled minimum `rho_c`, `rho_x`, and `E^2` agree within `1.14e-14` absolute; all components remain positive on the sampled grid.
- `a*theta_c` is constant to at most `2.12e-16` relative in the nonzero-theta diagnostic; `theta_c` stays identically zero in the physical CDM-comoving cases.
- The pointwise absolute-CDM identity has worst normalized residual `2.22e-16` (Gate 2: `2.09e-16`). The independent RK45 integration of `a^3 U` agrees with `rho_c*delta_c` to at most `3.20e-12` relative. This is below the audit's `1e-9` numerical comparison tolerance; Gate 2's DOP853 version reports a smaller maximum, `1.01e-13`.
- At `g=0`, `theta_c=0`, the analytic solution `delta_c=delta_i-S*(n-n_i)` is recovered to `9.70e-14` maximum relative error (saved: `9.68e-14`). For `g=+0.46666`, the source coefficient and contribution are positive; for `g=-0.46666`, both are negative, as required by the paper's sign convention.

The saved summary JSON retains extrema and diagnostic metrics, not full trajectories. Therefore, the quantitative code comparison is over saved case summaries; the independent trajectory was checked internally against the analytic `g=0` solution and the absolute-CDM identity. The source paper fixes `sum m_nu=0.06 eV`, while Gate 2 declares massless effective neutrinos; the audit reproduces that benchmark approximation, not the paper's exact late-time neutrino background.

The preserved first Gate 2 attempt failed before producing a physical result because its forward `solve_ivp` call received a descending `t_eval` grid. The corrected current script uses the ascending grid and the saved successful JSON. This is a recorded, corrected execution defect, not an equation or convention mismatch.

## Reproduction record

Successful command, from the project root (600-second cap; CPU 0 only):

```sh
taskset -c 0 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --state work/compute20/no_campaign_deadline.json --seconds 600 -- .venv/bin/python work/compute20/independent_gate2_audit.py
```

It completed in 3.32 seconds; no GPU was used. Environment: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1. Model/effort metadata available to this worker: GPT-6 family; the exact variant and reasoning-effort tier are not exposed by the runtime, so the inherited default was used.

- Independent implementation SHA-256: `4ac12933b867e31dc19d4969095df1712955d5dedb8089983bdc4566fd892341`
- Independent result JSON SHA-256: `fe7d7642a1b4b11ea68287f9afb916f935b233c913566eee01db7ba026f1bccd`
- Compared Gate 2 script SHA-256: `55090793914ccecefe9bdcbfcb2bcb9dd8314d3493294975a094316095728a89`
- Compared Gate 2 result JSON SHA-256: `c4d1e39f3ee368aed296515ee11d06bb4dad1a1d99e5818fa3d5553154beb58d`
- Primary paper PDF SHA-256: `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`

The machine-readable per-case comparison, code hash, parameters, command, and run metadata are in [independent_gate2_audit.json](independent_gate2_audit.json). This remains a fixed-metric dark-sector check: it omits Einstein feedback, radiation/baryon perturbation hierarchies, initial-condition generation, CMB spectra, and a stability classification.

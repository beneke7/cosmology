# Independent Gate 1 audit

**Disposition: independently checked on the declared massless-neutrino background benchmark.** The paper's source sign and dimensional convention match Gate 1. This does not reproduce the paper's massive-neutrino background or its likelihood results.

## Source equations and convention

I read and visually checked the primary paper PDF, pages 2–3: Yang et al., *Do DESI-DR2 BAO data imply a coupling of dark matter and dark energy?*, Eqs. (3)–(5), (8). The source PDF SHA-256 is `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`.

For the interacting-vacuum case, (w_x=-1), the paper gives

\[
\dot\rho_x=Q,\qquad \dot\rho_c+3H\rho_c=-Q,\qquad Q=\Gamma\rho_x.
\]

With (n=\ln a), (E=H/H_0), and (g=\Gamma/H_0), the implemented form is

\[
\frac{d\rho_x}{dn}=\frac{g\rho_x}{E},\qquad
\frac{d\rho_c}{dn}=-3\rho_c-\frac{g\rho_x}{E},\qquad
E^2=\Omega_ba^{-3}+\Omega_ra^{-4}+\rho_c+\rho_x.
\]

Positive (\Gamma) therefore adds to vacuum and removes from CDM forward in time, agreeing with the paper's explicit text. The code integrates in the independent coordinate (s=\ln(1+z)=-n), reversing the two interaction-term signs. (\Gamma) has Hubble-rate units; (g) is dimensionless. Densities are energy densities divided by today's critical energy density. The SI normalization uses (\rho_{\rm crit,0}=3H_0^2c^2/(8\pi G)); the paper writes the same relation in (c=1) units.

One declared model gap remains: the paper's page 2 footnote fixes \(\sum m_\nu=0.06\,\mathrm{eV}\), whereas the pinned benchmark treats neutrinos as massless, (\rho_\nu\propto a^{-4}\). I retained that benchmark for like-to-like checking. The result must not be described as the exact paper neutrino background.

## Independent radiation normalization

Using (a_{\rm rad}=4\sigma_{\rm SB}/c), (u_\gamma=a_{\rm rad}T_{\rm CMB}^4), and the massless-neutrino ratio (u_\nu/u_\gamma=(7/8)(4/11)^{4/3}N_{\rm eff}), I obtain:

| Quantity | Independent value |
| --- | ---: |
| (H_0) | (2.1842852410855023\times10^{-18}\,\mathrm{s}^{-1}) |
| (\rho_{\rm crit,0}) (energy) | (7.668947767821907\times10^{-10}\,\mathrm{J\,m}^{-3}) |
| (a_{\rm rad}) | (7.565733250280007\times10^{-16}\,\mathrm{J\,m}^{-3}\,\mathrm{K}^{-4}) |
| (u_\gamma) | (4.1748009199298555\times10^{-14}\,\mathrm{J\,m}^{-3}) |
| (\Omega_{\gamma0}) | (5.443772791681902\times10^{-5}) |
| (\Omega_{\nu0}) | (3.763360018025512\times10^{-5}) |
| (\Omega_{r0}=\Omega_{\gamma0}+\Omega_{\nu0}) | (9.207132809707413\times10^{-5}) |
| (\Omega_{x0}) from flat closure | (0.684907928671903) |

The independent photon density parameter differs from compute15 by ΔΩγ0 = 2.3924577253e-10 absolute, or 4.395e-6 relative. This is a rounded-constant sensitivity, not an implementation failure: compute15 uses a_rad = 7.5657e-16, while this calculation uses 4σSB/c = 7.56573325028e-16, a relative difference of about 4.4e-6. The same fractional offset appears in Ων0. At z=1e10, radiation supplies more than 0.9999995 of E², explaining the measured 4.39485e-6 relative E² offset.

## Numerical checks

The implementation is [independent_background.py](independent_background.py); it does not import Gate 1 code. It integrates 2,049 points from (z=0) to (10^{10}), with DOP853, a coarse run ((rtol,atol,h_{\max})=(10^{-9},10^{-12},0.10)), and a fine run ((10^{-12},10^{-14},0.02)). Isolated tolerance-only and maximum-step-only refinements are also recorded in the JSON output.

- At (g=0), the analytic flat LCDM solution is recovered to maximum relative errors (1.4423\times10^{-14}) in (\rho_c), (0) in constant (\rho_x), and (1.9164\times10^{-15}) in (E^2).
- Across the three rates, coarse-to-fine maximum relative changes over the output grid are (1.063\times10^{-10}) in (E^2), (3.239\times10^{-10}) in (\rho_c), and (4.081\times10^{-12}) in (\rho_x). Tightening only tolerance changes them by at most (3.168\times10^{-12}), (9.781\times10^{-12}), and (1.238\times10^{-13}), respectively. Tightening only maximum step changes them by at most (1.956\times10^{-15}), (5.594\times10^{-15}), and (3.007\times10^{-15}).
- (E^2,\rho_c,\rho_x) are positive at every one of the 2,049 output points for every rate. The minimum sampled values across all runs are (1), (0.266), and (0.44938827967897693), respectively. This is a grid-sampled check, not a proof of positivity at every continuum point.
- Two separately quadrature-evaluated identities hold over 65 locations on the dense solution: (d(\rho_c+\rho_x)/ds=3\rho_c), with worst scaled residual (7.17\times10^{-15}), and the integrated CDM comoving-density relation
  \[
  e^{-3s}\rho_c(s)-\Omega_{c0}-g\int_0^s e^{-3u}\frac{\rho_x(u)}{E(u)}du=0,
  \]
  with worst absolute residual (4.32\times10^{-15}) and worst scaled residual (1.44\times10^{-14}).

## Comparison with compute15

The following are absolute and relative endpoint differences at (z=10^{10}), using compute15's final JSON (its script SHA-256: `c07ede03539d08702f24db69140badd8361c859e76be608fb246ad1867c91d9c`). Relative acceptance limits were (5\times10^{-6}) for (E^2), (10^{-9}) for (\rho_c), and (2\times10^{-9}) for (\rho_x). The larger (E^2) allowance accounts for the independently identified radiation-constant rounding above.

| (g) | (E^2: |\Delta|\) (relative) | (\rho_c: |\Delta|\) (relative) | (\rho_x: |\Delta|\) (relative) |
| ---: | ---: | ---: | ---: |
| 0 | (4.0464\times10^{30}\) ((4.3949\times10^{-6})) | (0) ((0)) | (4.0464\times10^{-10}) ((5.9079\times10^{-10})) |
| +0.46666 | (4.0464\times10^{30}\) ((4.3949\times10^{-6})) | (6.2296\times10^{19}) ((1.8186\times10^{-10})) | (8.6641\times10^{-12}) ((1.9280\times10^{-11})) |
| −0.46666 | (4.0464\times10^{30}\) ((4.3949\times10^{-6})) | (1.0276\times10^{20}) ((6.1064\times10^{-10})) | (2.0601\times10^{-9}) ((1.8598\times10^{-9})) |

All endpoint comparisons fall within the stated relative tolerances. The endpoint (E^2) offset is the expected radiation normalization difference; CDM and vacuum agree more closely. This comparison is endpoint-only because compute15's JSON does not retain its full trajectory.

## Reproduction

From the project root, run:

```sh
.venv/bin/python work/compute17/independent_background.py
```

Successful run: 0.90 s, CPU affinity `[0]`, one numerical-library thread, no GPU. Environment was Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1. Code SHA-256: `51741fe15312ede9960cdee3902313387468816944fdb824e8c3e9c50d2cbff0`. Machine-readable output, including the exact code/source hashes, tolerances, measured values, and command, is [independent_background.json](independent_background.json).

Limitations: this is a homogeneous background check only. It does not include massive-neutrino evolution, perturbations, sound-horizon/likelihood consistency, or a likelihood reproduction. Positivity is sampled on the fixed output grid. The compute15 comparison uses endpoint values saved in its JSON, not a full-trajectory pointwise comparison.

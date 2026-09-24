# Early-time positivity of the BAO-screen IVS point

**Finding.** This bounded background calculation extends the supplied screen point (\(\Omega_{m0}=0.3869710077\), \(g=-0.4666647299\), \(w_x=-1\)); it does not refit BAO and gives no posterior. The declared grid contains 16 combinations of \(H_0\) and baryon fraction. Of these, 0 have a CDM zero crossing by \(z=10^{10}\). Crossing redshifts: none. At \(z=10^{10}\), the smallest comoving CDM value across the grid is 0.206114, still positive. The full rows are below.

## Identity

With \(n=\ln a\), \(E=H/H_0\), \(u_i=\rho_i/\rho_{crit,0}\), and \(q_c=a^3u_c\), the source-paper equations give \(du_c/dn=-3u_c-gu_v/E\). Hence \(dq_c/dn=-g a^3u_v/E\), so integrating backward from today,

\[q_c(n)=\Omega_{c0}+g\int_n^0 e^{3s}\frac{u_v(s)}{E(s)}\,ds.\]

For this negative \(g\), \(q_c\) decreases monotonically as we move into the past, since \(u_v/E>0\). Thus a positive value at the declared maximum redshift excludes an intermediate crossing. Since \(u_c=q_c/a^3\), a zero in the identity is exactly a zero in physical CDM density. The threshold is conditional on \(\Omega_{c0}=(1-f_b)\Omega_{m0}\); neither BAO amplitude nor this background fit measures \(f_b\).

## Declared assumptions and limits

The local paper's background equations conserve baryons, radiation, and neutrinos separately and specify \(N_{eff}=3.044\), \(\sum m_\nu=0.06\) eV. Photons use \(\Omega_\gamma h^2=2.4728\times10^{-5}\) at \(T_{CMB}=2.7255\) K. Neutrinos are treated here as massless radiation at all redshifts, explicitly an approximation to the paper's nonzero mass prescription. Flat closure sets \(\Omega_{v0}=1-\Omega_{m0}-\Omega_{r0}\). The sensitivity grid is \(H_0\in\{55.0, 65.0, 75.0, 85.0\}\) km s\(^{-1}\) Mpc\(^{-1}\), \(f_b\in\{0.10, 0.15, 0.20, 0.25\}\), and \(0\le z\le10^{10}\). No CMB or BBN prior is imported; H0 and \(f_b\) are sensitivity coordinates because the BAO-only amplitude \(\alpha=c/(H_0r_d)\) does not identify them.

The numerical solutions use DOP853 in two independently written coordinates, \(n=\ln a\) and \(u=\ln(1+z)\), and an adaptive quadrature check of the exact comoving-CDM identity. The tolerance-refinement metric is the maximum component difference scaled by \(\max(1,|u_i|)\): 1.928e-14; maximum identity residual over all cases: 3.747e-15; maximum disagreement between coordinate integrations: 0.000e+00. At very high redshift the result is a homogeneous background extrapolation only; it does not validate perturbations, sound-horizon physics, neutrino mass transitions, or early-universe viability in a full model.

## Grid results

| H0 | f_b | Omega_gamma0+nu0 | rho_c zero z | q_c(z=1e10) | positive to z=1e10 | q_c identity error | n/u max diff |
|---:|---:|---:|---:|---:|:---:|---:|---:|
| 55 | 0.10 | 0.0001383 | none | 0.264177 | yes | 2.16e-15 | 0.00e+00 |
| 55 | 0.15 | 0.0001383 | none | 0.244829 | yes | 2.03e-15 | 0.00e+00 |
| 55 | 0.20 | 0.0001383 | none | 0.225480 | yes | 3.08e-15 | 0.00e+00 |
| 55 | 0.25 | 0.0001383 | none | 0.206131 | yes | 1.78e-15 | 0.00e+00 |
| 65 | 0.10 | 9.899e-05 | none | 0.264169 | yes | 2.78e-15 | 0.00e+00 |
| 65 | 0.15 | 9.899e-05 | none | 0.244820 | yes | 2.05e-15 | 0.00e+00 |
| 65 | 0.20 | 9.899e-05 | none | 0.225472 | yes | 3.75e-15 | 0.00e+00 |
| 65 | 0.25 | 9.899e-05 | none | 0.206123 | yes | 1.39e-15 | 0.00e+00 |
| 75 | 0.10 | 7.435e-05 | none | 0.264163 | yes | 2.28e-15 | 0.00e+00 |
| 75 | 0.15 | 7.435e-05 | none | 0.244815 | yes | 3.19e-15 | 0.00e+00 |
| 75 | 0.20 | 7.435e-05 | none | 0.225466 | yes | 2.03e-15 | 0.00e+00 |
| 75 | 0.25 | 7.435e-05 | none | 0.206118 | yes | 2.00e-15 | 0.00e+00 |
| 85 | 0.10 | 5.789e-05 | none | 0.264160 | yes | 3.00e-15 | 0.00e+00 |
| 85 | 0.15 | 5.789e-05 | none | 0.244811 | yes | 2.14e-15 | 0.00e+00 |
| 85 | 0.20 | 5.789e-05 | none | 0.225463 | yes | 1.72e-15 | 0.00e+00 |
| 85 | 0.25 | 5.789e-05 | none | 0.206114 | yes | 1.78e-15 | 0.00e+00 |

Reproduce from the repository root with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python work/theory5/ivs_early_positivity.py
```

Machine-readable detail and hashes: [ivs_early_positivity.json](ivs_early_positivity.json). Runtime: 3.128 s. Script SHA-256: `004431c201daae6f638663502986068578e3289882545e0f8f2c53681392bba7`.

**Stop disposition:** completed at \(z=10^{10}\); no CDM crossing was found in the declared 16-case grid. No full likelihood, CAMB, CMB prior, or extrapolation beyond the stated redshift bound was used.

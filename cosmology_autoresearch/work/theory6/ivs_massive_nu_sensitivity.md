# Thermal massive-neutrino sensitivity of the IVS positivity check

**Finding.** The same supplied background point \(\Omega_{m0}=0.3869710077\), \(g=-0.4666647299\), \(w_x=-1\) remains CDM-positive through \(z=10^{10}\) in all 16 \(H_0,f_b\) cases under the thermal massive-neutrino prescription below. This matches theory5's massless-radiation disposition: 0 massive-case crossings versus 0 massless-case crossings. The maximum absolute shift in endpoint comoving CDM is 3.729e-04; its range and each case are tabulated below. No likelihood was refit and no posterior is inferred.

## Neutrino prescription and formula

The local paper fixes \(N_{eff}=3.044\) and \(\sum m_\nu=0.06\) eV. I implement one fully populated thermal massive eigenstate with \(m=0.06\) eV and the remaining \(2.044\) effective species massless. With \(T_{\nu0}/T_{\gamma0}=(4/11)^{1/3}\), \(y=ma/T_{\nu0}\), and

\[F(y)=\int_0^\infty dq\;\frac{q^2\sqrt{q^2+y^2}}{e^q+1},\qquad F(0)=\frac{7\pi^4}{120},\]

the density used in the Friedmann equation is

\[\frac{\rho_\nu(a)}{\rho_{\gamma0}}=\left(\frac{15}{\pi^4}\right)\left(\frac{T_{\nu0}}{T_{\gamma0}}\right)^4a^{-4}\left[F\!\left(\frac{ma}{T_{\nu0}}\right)+(N_{eff}-1)F(0)\right].\]

The Fermi-Dirac energy integral is evaluated with adaptive QUAD in the background RHS, not a massless/high-temperature switch. Flat closure uses \(\Omega_{v0}=1-\Omega_{m0}-\Omega_{\gamma0}-\Omega_{\nu0}\); \(\Omega_{m0}\) retains the screen's baryon-plus-CDM meaning. Grid: \(H_0\in\{55.0, 65.0, 75.0, 85.0\}\), \(f_b\in\{0.10, 0.15, 0.20, 0.25\}\), and \(0\le z\le10^{10}\). The earlier result and this one use no CMB/BBN prior. This is a declared thermal background prescription and is not an exact reproduction of the authors' unknown CAMB neutrino configuration.

The comoving-CDM identity remains \(q_c(n)=\Omega_{c0}+g\int_n^0e^{3s}u_v(s)/E(s)\,ds\). For \(g<0\), it declines monotonically into the past, so the positive endpoint excludes missed earlier zeros in-range. ODE tolerance refinement gives max scaled component difference 1.900e-14; adaptive FD integral tightening agrees at probe values to relative difference 0.000e+00, with maximum reported tight-quadrature absolute error 1.18e-11; identity residual across the grid is at most 2.998e-15. These checks support numerical convergence for this background calculation, not perturbative or sound-horizon viability.

## Case-by-case comparison

| H0 | f_b | q_c(zmax), massive | q_c(zmax), massless | delta q_c | massive zero z | massless zero z |
|---:|---:|---:|---:|---:|---:|---:|
| 55 | 0.10 | 0.26454991 | 0.26417706 | +3.73e-04 | none | none |
| 55 | 0.15 | 0.24520136 | 0.24482851 | +3.73e-04 | none | none |
| 55 | 0.20 | 0.22585281 | 0.22547996 | +3.73e-04 | none | none |
| 55 | 0.25 | 0.20650426 | 0.20613141 | +3.73e-04 | none | none |
| 65 | 0.10 | 0.26443572 | 0.26416863 | +2.67e-04 | none | none |
| 65 | 0.15 | 0.24508717 | 0.24482008 | +2.67e-04 | none | none |
| 65 | 0.20 | 0.22573862 | 0.22547153 | +2.67e-04 | none | none |
| 65 | 0.25 | 0.20639007 | 0.20612298 | +2.67e-04 | none | none |
| 75 | 0.10 | 0.26436402 | 0.26416335 | +2.01e-04 | none | none |
| 75 | 0.15 | 0.24501547 | 0.24481480 | +2.01e-04 | none | none |
| 75 | 0.20 | 0.22566692 | 0.22546625 | +2.01e-04 | none | none |
| 75 | 0.25 | 0.20631837 | 0.20611770 | +2.01e-04 | none | none |
| 85 | 0.10 | 0.26431608 | 0.26415982 | +1.56e-04 | none | none |
| 85 | 0.15 | 0.24496753 | 0.24481126 | +1.56e-04 | none | none |
| 85 | 0.20 | 0.22561898 | 0.22546271 | +1.56e-04 | none | none |
| 85 | 0.25 | 0.20627043 | 0.20611416 | +1.56e-04 | none | none |

Reproduce from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python work/theory6/ivs_massive_nu_sensitivity.py
```

Machine-readable results and input hashes: [ivs_massive_nu_sensitivity.json](ivs_massive_nu_sensitivity.json). Runtime: 9.516 s. Script SHA-256: `f0dae2ba45f4c8785822dd5db3cb76741462e204b4d658b47277d3961ae14347`.

**Stop disposition:** completed at \(z=10^{10}\). The finite-mass prescription does not change the positivity result on the declared 16-case grid. No BAO refit, CAMB run, download, or parameter prior was introduced.

# Gate 2: fixed-metric IVS perturbation equations

**Result: PASS for the declared fixed-metric two-component tests.** This checks the primary paper's Eqs. (19)–(20), PDF p. 4, not a coupled Einstein–Boltzmann solution. The exact PDF SHA-256 is `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`; the local text copy SHA-256 is `a99a4afdd858cdce3cbc2c41f0a13098c30c361dc67891331ee7d9d5f824868f`.

The standalone [perturbation_gate2.py](perturbation_gate2.py) independently rebuilds the pinned Gate 1 background without importing compute15: H0 = 67.4 km/s/Mpc, Omega_m0 = 0.315, Omega_b0 = 0.049, T_CMB = 2.7255 K, N_eff = 3.044, separately conserved massless neutrinos, and flat closure. It uses Q = Gamma*rho_x, w_x = -1, and the paper convention Gamma > 0 transfers energy from CDM to vacuum. Rates are Gamma/H0 = 0, +0.46666, and -0.46666; the checked interval is z = 10^10 to today.

With n = ln(a), E = H/H0, and Hcal = a*H0*E, the paper equations become:

```text
S = h'/(2*Hcal)
d(delta_c)/dn = -S - theta_c/Hcal + (g/E)*(rho_x/rho_c)*delta_c
d(theta_c/H0)/dn = -(theta_c/H0)
```

The fixed metric is analytic: S = 1e-3 and h(n) = 2*S*(n - n_min). The physical gauge run sets theta_c = 0. A separate nonzero-theta diagnostic checks Eq. (20) and a*theta_c = constant; it is not treated as a CDM-comoving solution.

Across all six rate/gauge cases, the perturbation ODE agrees with an independent integrating-factor solution using composite Simpson quadrature on 16,385 nodes to at most 6.52e-12 relative error. Changing the quadrature grid from 8,193 to 16,385 nodes changes delta_c by at most 9.72e-11. At Gamma/H0 = 0 and theta_c = 0, the analytic solution delta_c = delta_i - S*(n - n_i) is recovered to 9.69e-14. Theta_c = 0 remains exactly zero; in the nonzero diagnostic, a*theta_c is constant to 3.60e-15 relative error.

The unit conversion a*Gamma/Hcal = (Gamma/H0)/E agrees to 1.12e-16 absolute error. For positive Gamma, the source coefficient and its delta_c contribution are positive; for negative Gamma, both are negative. From the density-contrast equation and background continuity, U = rho_c*delta_c obeys `dU/dn + 3*U + rho_c*(theta_c/Hcal + S) = 0`: the interaction terms cancel. Its maximum normalized residual is 2.09e-16, and an independent integration of this absolute-CDM equation agrees with rho_c*delta_c to 1.01e-13 relative error.

The successful second invocation used one CPU (affinity pinned to CPU 0), no GPU, and the required 600-second cap with a missing state path. Script runtime was 4.1571 seconds; bounded-command wall time was about 4.406 seconds. Exact command:

```text
taskset -c 0 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --state work/compute18/no_campaign_deadline.json --seconds 600 -- .venv/bin/python work/compute18/perturbation_gate2.py
```

The first invocation failed before producing a physical result because the forward solver received a descending evaluation grid; the exact command, exception, and correction are retained in [attempt_001_failure.json](attempt_001_failure.json). Final measurements and assumptions are in [perturbation_gate2.json](perturbation_gate2.json).

Artifact SHA-256: script `55090793914ccecefe9bdcbfcb2bcb9dd8314d3493294975a094316095728a89`; result JSON `c4d1e39f3ee368aed296515ee11d06bb4dad1a1d99e5818fa3d5553154beb58d`; preserved first-attempt record `f466161e40795abf4062529e6b260b973c3cb0f75abeeb77eb5308ec3921b3a2`.

This is **not** a perturbative stability classification, CMB validation, or full-model reproduction. The metric is prescribed, so Einstein feedback is absent; photon, baryon, and neutrino perturbation hierarchies, initial-condition generation, and CMB spectra are not included.

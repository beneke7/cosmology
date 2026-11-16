# Compute30 boundary-gradient audit

Status: **independently checked; no inward \(\omega_b h^2\) downhill direction found in the tested neighborhood**. A same-edge \(H_0\)/\(\omega_c h^2\) polish gives \(\Delta\chi^2=-5.05e-11\) relative to the supplied point, far below compute26's declared \(2\times10^{-5}\) optimizer agreement tolerance. This is a BAO-only flat-LCDM finite-box/profile result, not physical cosmological inference.

The standalone implementation reproduces the candidate score as \(10.273886272857\), matching compute26 \(10.273886272857\). At fixed \(H_0\) and \(\omega_ch^2\), moving either inward or outward from \(\omega_b h^2=0.1\) raises the score at all tested steps \(10^{-5}\), \(10^{-4}\), and \(5\times10^{-4}\). For example, the inward \(10^{-4}\) point is higher by \(0.00805388275\); the outward point is outside the declared box and is diagnostic only. The compute26 outward stress therefore did not imply an inward feasible descent.

Re-optimizing \(H_0\) and \(\omega_ch^2\) from three starts at 9 fixed \(\omega_b h^2\) slices puts the best sampled score at the upper edge: \(\chi^2=10.273886272807\). Every tested inward slice is higher; the nearest, \(\omega_b h^2=0.0999999\), is higher by \(4.05e-09\), and the lowest inward-slice increase is \(4.05e-09\). The polished edge point has \(H_0=133.613428198\), \(\omega_c h^2=0.430272600279\), and \(100\theta_{MC}=1.10012396\), all inside the declared domain. The tested profile supports a local minimum at the \(\omega_b h^2\) boundary. The tiny literal improvement over the supplied point reflects final same-edge polishing; discrete slices and local Powell solves do not establish global optimality.

The score uses the full released 13-row order and full 13x13 covariance; both dense precision solve and Cholesky whitening are checked at every selected point. CAMB 2.0.4 is configured for flat fluid \(\Lambda\)CDM with \(\Omega_k=0\), \(m_\nu=0.06\) eV, \(N_{eff}=3.044\), \(T_{CMB}=2.7255\) K, one degenerate massive species and BBN-consistent helium. Every retained point has \(0<\Omega_m<1\), \(\Omega_\Lambda>0\), and finite positive background distances and Hubble rates. The CAMB \(\theta_{MC}\) setter, with expanded internal H0 root bracket \([1,2000]\), independently maps each selected \(H_0\) back to the target theta interval.

## Reproduction

```sh
cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute30/no_campaign_deadline.json --seconds 560 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute30/boundary_gradient_audit.py
```

The run is bounded at 560 seconds and fixes numerical-library threads to one. Detailed finite differences, all optimizer starts, score vectors, closure checks, versions and source hashes are in [boundary_gradient_audit.json](boundary_gradient_audit.json). Script SHA-256: `6085f7b281fae413a4c9a00bac5ece6444e250ba0a6a667e4ac56445bff18fc9`. Result SHA-256: `1550c66aaccf637e7428268168cf40181b82351620f69b9f9a1815e789897cca`.

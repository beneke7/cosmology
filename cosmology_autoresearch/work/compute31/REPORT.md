# Compute31: standard-BBN-region DESI DR2 BAO profile

Status: **exploratory, independently implemented CAMB × DESI DR2 BAO profile**. The data-only minimum is pinned to the upper scan edge \(\omega_bh^2=0.0240\), so the scan does not establish an interior baryon-density minimum. Three separate BBN prior sensitivities are reported. The user-directed conservative prior remains provisional; the source-matched 2026 input is a preprint translated to a local piecewise Gaussianized penalty rather than importing its likelihood. This is not a posterior, evidence calculation, or DESI collaboration parameter reproduction.

The released 13-row BAO mean vector and full 13×13 covariance were parsed and scored in their pinned order. The prediction code independently evaluates (D_V/r_d), (D_M/r_d), and (D_H/r_d=(c/H)/r_d) with local CAMB 2.0.4. The fixed setup is flat fluid ΛCDM, (m_\nu=0.06\,\mathrm{eV}), (N_{eff}=3.044), one degenerate massive species, (T_{CMB}=2.7255\,\mathrm{K}), and CAMB's BBN-consistent helium setting (`YHe=None`).

## Fit and checks

For each fixed \(\omega_bh^2\) in [0.0205, 0.0240], the code reoptimized \(H_0\) and \(\omega_ch^2\) in \([10,1000]\,\mathrm{km\,s^{-1}\,Mpc^{-1}}\) and [0.001, 0.99] using four or five deterministic Powell starts. The coarse profile has 18 points; the local refinement has 9 points. At the lowest sampled data-only score, \(\chi^2=10.2814618861\), \(\omega_bh^2=0.0240000\), \(H_0=69.993182\), and \(\omega_ch^2=0.1209550\). The refined-versus-coarse profile change is \(\Delta\chi^2=6.75e-14\) at the respective minima; it is a resolution check, not proof of globality.

The coarse and refined sampled minima agree at the upper boundary with \(\Delta\chi^2=6.75e-14\); this is a local resolution check. The multistart audit found 10 valid but much worse bounded solutions, generally ending at the \(\omega_ch^2=0.99\) upper bound from broad high-\(H_0\) starts. These alternatives remain in the JSON; only 103 of 125 fixed-slice starts landed within \(10^{-7}\) of that slice's best start. The lowest-scoring basin was independently recovered from standard-scale starts, but globality over the full continuous box is not proved. Powell success flags and every start are retained, including non-success outcomes. Selected fixed-slice dense-precision and Cholesky-whitening scores agree to at most 7.11e-15. Every selected profile point was re-mapped through CAMB's `cosmomc_theta` setter and checked against broad \(100\theta_{MC}\in[0.5,10]\) support; all selected minima have positive flat-Λ closure.

Three distinct prior sensitivities were evaluated as alternatives and were never multiplied together. Poulin et al., [arXiv:2607.20635v1](https://arxiv.org/html/2607.20635), Appendix A/Table 2, report a BBN-only LambdaCDM result based on updated D/H = 2.508 +/- 0.030, 100 omega_b = 2.2017 (+0.0221/-0.0217) at 68% (preprint). I use a piecewise Gaussianized penalty with omega_b h2 = 0.022017 (+0.000221/-0.000217); this is not the source paper's exact profile. The separate PDG 2025 SBBN sensitivity supplied by the orchestrator is 0.02205 +/- 0.00043. The provisional conservative sensitivity is 0.0222 +/- 0.0005; its source provenance remains unreviewed. The reoptimized omega_b h2 values are 0.0220170 (preprint), 0.0220500 (PDG 2025), and 0.0222000 (provisional), respectively. The PDG entry's page/URL was not present in this worker's local inputs and remains to be attached by the orchestrator.

## Reproduction

The isolated script ran in 125.80 seconds with one numerical worker, numerical-library threads set to one, and no GPU. No installs or downloads were used. The 850-second command cap is below the task's 900-second limit:

```sh
cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute31/no_campaign_deadline.json --seconds 850 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute31/bbn_bao_profile.py
```

The machine-readable record contains exact bounds and prior form, all profile slices and optimizer starts, invalid probes and non-success outcomes, software versions, run time, input and code hashes, and the score/closure/theta checks. Result SHA-256: `790250c7948afec861c8d9aef1c2c5ef1028e11b8c3144b246fc4485fc6b3326`. Software: CAMB 2.0.4, NumPy 2.5.3, SciPy 1.18.1, Python 3.12.3.

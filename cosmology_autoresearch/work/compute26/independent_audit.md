# Independent CAMB ΛCDM × DESI DR2 BAO profile

Status: reproduced as a bounded, single-probe profile calculation; incomplete as an accepted minimum because the converged profile is pinned to the upper prior boundary in \(\omega_bh^2\). The optimizer candidate is exploratory and is not reported as an accepted interior minimum.

## Audit of compute24

The physical mapping in `work/compute24/camb_bao_lcdm.py` is internally consistent: it passes `100theta_MC / 100` to CAMB's `cosmomc_theta`, sets flat geometry, fixes standard Λ (`w=-1`) before CAMB solves for \(H_0\), and uses `omk=0`. Its background observables use \(D_M/r_d\), \(D_H/r_d=(c/H)/r_d\), and \(D_V/r_d=[zD_M^2D_H]^{1/3}/r_d\). The massive-neutrino, temperature, and BBN-consistency choices are explicit. `As`, `ns`, and `tau` do not affect this background-only BAO calculation.

The table limits \(100\theta_{MC}\in[0.5,10]\), \(\omega_bh^2\in[0.005,0.1]\), and \(\omega_ch^2\in[0.001,0.99]\) are the target paper's broad Table-I parameter priors (local PDF p. 10, Table I), originally used in joint CMB+BAO+SN analyses. They can be declared as a finite BAO-only search box, but do not reproduce those joint constraints or a posterior. The hard `theta_H0_range=(10,100)` is much narrower in practice: for \(\omega_bh^2=0.0224,\omega_ch^2=0.12\), `100theta_MC=1.2` has a valid CAMB root at \(H_0\simeq129.48\) only when the range is expanded. My independently selected profile point has inferred `100theta_MC=1.10012392067`, \(H_0=133.6134163\), lies inside the paper's theta box, and fails CAMB's original \([10,100]\) root bracket. Expanding the root bracket to \([1,2000]\) recovers the direct-\(H_0\) value with absolute error \(2.3\times10^{-13}\) km/s/Mpc. Thus compute24's root bracket does omit valid parts of the stated box.

Invalid-point handling in compute24 uses `1e100`. That value is finite, and its seed 20260925 stopped after one generation with every population energy equal to `1e100`, while SciPy reported “Optimization terminated successfully.” The second seed reached its 60-iteration cap at \(\chi^2=10.27636166\), so its `success` flag was false. The seed marked successful was not a valid candidate. The compute24 script does check the full covariance with both Cholesky whitening and a dense solve; those agree at the selected valid point to \(3.6\times10^{-15}\). Its 13-row ordering, including high-redshift `DH` then `DM`, agrees with the pinned input.

## Independent profile

`independent_lcdm_bao_profile.py` uses the pinned CAMB 2.0.4 wheel and a direct-\(H_0\) parameterization to remove the nested theta root from optimizer calls. It uses IEEE float64, all 13 rows and their full 13×13 covariance, flat ΛCDM, \(m_\nu=0.06\) eV, \(N_{eff}=3.044\), \(T_{CMB}=2.7255\) K, one degenerate massive species, and CAMB's BBN-consistency default. Bounds are \(H_0\in[10,1000]\) km/s/Mpc and the two Table-I physical-density bounds. The \(H_0\) interval is an explicitly broad profile coordinate, not an Hubble prior. This direct-\(H_0\) finite box is not asserted to be exactly the same domain as a uniform finite box in \(100\theta_{MC}\); the selected point itself is mapped back and verified inside the paper's theta interval.

Four independent differential-evolution seeds (20260925–20260928), each with up to 600 iterations, 48 starting points and one CPU worker, converged to the same point. The pairwise scaled \(L_\infty\) parameter distances were below \(2.7\times10^{-7}\), and the largest pairwise \(\chi^2\) difference was below \(3.4\times10^{-10}\), within declared tolerances of \(2\times10^{-3}\) and \(2\times10^{-5}\). However, every seed converged to \(\omega_bh^2=0.1\), the upper Table-I boundary. This fails the required interior-minimum gate, so the machine record labels it `incomplete_boundary_pinned_profile` and `candidate_acceptance.accepted=false`.

At that boundary point, \(H_0=133.6134163\) km/s/Mpc, \(\omega_ch^2=0.43027242\), \(\Omega_{m0}=0.29739027\), \(r_d=75.99673\) Mpc, and the full-covariance score is \(\chi^2=10.27388627\). This is a reproducible BAO-only boundary profile, not a claimed best-fit cosmology. A fixed-\(H_0\), fixed-\(\omega_ch^2\) edge check has \(\chi^2=10.27388627\) at \(\omega_bh^2=0.1\) and \(10.28192530\) at 0.1001 (outside the declared box); that is a one-dimensional diagnostic, not an extension of the accepted domain.

The independent distance check integrated \(c/H(z)\) by adaptive quadrature at all 13 redshifts. Its maximum absolute difference from CAMB's \(D_M(z)\) is \(9.1\times10^{-13}\) Mpc. At the profile point, dense and Cholesky full-covariance scores are both 10.273886272857048 (difference zero). Invalid evaluations were retained as categories and examples: 158 points failed the flat-positive-Λ closure check and 27 raised CAMB errors; they received a separate \(10^{80}\) sentinel and were never eligible for candidate acceptance.

This is only the standard flat-ΛCDM background against released DESI DR2 compressed BAO. It includes no Planck or ACT likelihood, supernova likelihood, posterior sampling, evidence calculation, or interacting-dark-energy model.

## Reproduction record

Run from the workspace root; the wrapper bounded the process at 540 seconds. Actual runtime was about 92 seconds with one worker and one numerical-library thread:

```sh
cd cosmology_autoresearch && .venv/bin/python scripts/run_bounded.py --state work/compute26/task_state.json --seconds 540 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute26/independent_lcdm_bao_profile.py
```

Versions: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, CAMB 2.0.4 (local unpacked site from the pinned wheel); Linux x86_64/glibc 2.39. Exact command, UTC timestamps, optimizer seed records, criteria, hashes, distances and covariance checks are in `independent_lcdm_bao_profile.json`.

Input hashes:

- DESI DR2 BAO mean: `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`
- DESI DR2 BAO covariance: `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`
- CAMB 2.0.4 wheel: `9b7e6558b2049757e96496b8a33461d8ff11acd79017cfbd00a573dbc8d1e366`
- Target-paper PDF: `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`

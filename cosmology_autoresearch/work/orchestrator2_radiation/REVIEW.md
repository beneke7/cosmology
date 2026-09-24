# Independent radiation-lane review

Status: independently checked; exploratory reference-model sensitivity screen complete, 2026-09-24. Lane owner: `/root/orchestrator2_radiation`, gpt-6-astra/max as assigned by the campaign root. Specialist workers explicitly launched as gpt-6-luna/max: `theory2` and `compute2`.

## Scope and input contract

This lane checks the radiation-neglect approximation in the existing 13-row DESI DR2 BAO background screen. It does not replace the DESI collaboration inference or create posterior, evidence, significance, or H0 constraints. The released covariance, mean ordering, model bounds, and independently free distance amplitude must be preserved.

Baseline files reviewed before work:

| File | SHA-256 |
|---|---|
| `scripts/background_bao.py` | `34bdf1d81397b39b740faa3e2f4a79748d5fa9d3bfdd8303619d631dbf8aa2e2` |
| `context/data/desi_dr2_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `context/data/desi_dr2_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| `experiments/campaign_seed_bao/result.json` | `92662f3ee42500d2be34d9f87b485fe56c82e9550fa674e42f4dec2bf2cab1f0` |

The contract is one Gaussian statistic `(data-prediction)^T C^-1 (data-prediction)` using the complete 13 by 13 covariance in release row order, including the final `DH`, then `DM` pair at z=2.33. The domains are alpha in [1e-6,1e4], Omega_m in [0.05,0.6], w0 in [-2,-0.3], and wa in [-3,3], with unused parameters fixed to LCDM/wCDM limits. Analytic profiling of alpha is equivalent to the original bounded optimization because every distance ratio is linear in positive alpha.

For this sensitivity screen H0 is held at declared values from 50 to 90 km/s/Mpc. These values specify the radiation fraction, not a sampled posterior prior; alpha=c/(H0 r_d) remains independently free. Radiation is parameterized by a fixed physical density omega_r=Omega_r h^2 derived from T_CMB=2.7255 K and reference N_eff=3.046, with flat closure Omega_de=1-Omega_m-Omega_r. The original fixed Omega_r=9e-5 seed calculation was explicitly illustrative and did not refit; it is not used as a precision validation.

## Physical-model limitation

The all-massless-neutrino prescription is a declared reference extension. The DESI DR2 paper, local `context/papers/desi_dr2_bao.pdf`, printed p.5, Eqs. (6)-(7), separates photons, baryons/CDM, and neutrino density with its relativistic-to-nonrelativistic transition. Its Omega_m definition includes neutrinos when nonrelativistic. Accordingly, the reference calculation must not be described as reproducing the DESI massive-neutrino background or proving a bound for arbitrary physical neutrino assumptions. CPL density evolution follows Eq. (10), and the distance definitions follow Eqs. (4)-(5).

The claimed measurement-precision comparison must identify where it was evaluated: baseline best fits and radiation-refitted optima over the stated H0 scan. A claim over the entire wide shape domain would require a separate calculation.

## Work ownership and review plan

- `work/theory2/`: independent constants, dimensional derivation, first-order radiation response, and scalar numerical checks.
- `work/compute2/`, `experiments/radiation_sensitivity/`: exact covariance profile refits and computation record.
- `work/orchestrator2_radiation/`: lane review, independent reproduction of saved optima, and handoff.

Both workers use CPU FP64 and one numerical-library thread; neither uses the GPU. No baseline/shared campaign state, source maps, manifest, ledger, morning report, or continuation report is owned by this lane. No commits, pushes, or warden calls are made here.

## Independent fixed-baseline result

`scalar_reference.py` uses scalar `scipy.integrate.quad` integrals and dense solves with the complete raw covariance; it imports no campaign or worker prediction/likelihood code. Across all three supplied seed fits, predictions reproduce within 4.27e-14 and chi-squared within 7.29e-14. The independently solved amplitude differs from the saved optimizer amplitude by at most 6.02e-8.

From the SI blackbody expression

`u_gamma = 8*pi^5*k_B^4*T^4/(15*h_P^3*c^3)`,

`omega_gamma = (u_gamma/c^2)/(3*H100^2/(8*pi*G))`,

and `omega_r = omega_gamma*[1+(7/8)*(4/11)^(4/3)*N_eff]`,

the independent reference yields omega_gamma = 2.472975328714088e-5 and omega_r = 4.183702725849734e-5. Here lowercase omega denotes Omega*h^2, so no further h^2 factor belongs on the left-hand side. The unit convention is exactly 1 pc = 648000/pi AU with AU = 149597870700 m; `scipy.constants.parsec` agrees. Omega_r therefore ranges from 1.6734810903398935e-4 at H0=50 to 5.165065093641646e-5 at H0=90.

Holding each saved baseline's alpha and shape fixed, the H0=50 endpoint gives:

| Model | Maximum absolute fractional distance shift | Maximum shift / marginal measurement error | Full-covariance displacement norm | Delta chi-squared at fixed parameters |
|---|---:|---:|---:|---:|
| LCDM | 0.0008720953 | 0.07452707 | 0.12761449 | +0.03962239 |
| wCDM | 0.0008497292 | 0.07307494 | 0.12257911 | -0.00066914 |
| CPL | 0.0007291893 | 0.06236726 | 0.10429800 | +0.00842548 |

The norm is `sqrt(delta_prediction^T C^-1 delta_prediction)`, a measurement-space displacement diagnostic, not a discovery significance. The independent H0 grid is 50,55,60,65,67.4,70,75,80,85,90. The refitted minima below independently confirm that the profile model comparison changes little under this reference correction.

Reproduction, from the project root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ./.venv/bin/python scripts/run_bounded.py --seconds 120 -- ./.venv/bin/python work/orchestrator2_radiation/scalar_reference.py
```

The numerical work takes about 0.005 seconds (about 0.09 seconds for the wrapped process in the review environment). The script emits `scalar_reference.json`, including exact source/code hashes, constants, all predictions/shifts, dependency versions, command, and timing. It verifies the full-covariance identity `Delta_chi2 = -2*r^T C^-1 delta_prediction + delta_prediction^T C^-1 delta_prediction` within 1e-10.

## Profile refits and independent review

`review_saved_fits.py` independently verified all 33 saved optima: three no-radiation fits plus three models at each of ten H0 values. It enforces the agreed domains and flat closure, re-evaluates scalar adaptive distances with raw full-covariance solves, verifies the analytic amplitude against an independent dense solve, checks code/data/result hashes, and checks the fixed-parameter shifts separately from refitted predictions.

| Model | No-radiation chi-squared | Radiation chi-squared at H0=50 | Radiation chi-squared at H0=90 | Refitted delta chi-squared range over grid |
|---|---:|---:|---:|---:|
| LCDM | 10.271041003 | 10.294796863 | 10.278294351 | +0.00725335 to +0.02375586 |
| wCDM | 9.041046606 | 9.025343808 | 9.036191833 | -0.01570280 to -0.00485477 |
| CPL | 5.700615936 | 5.698181033 | 5.699859565 | -0.00243490 to -0.00075637 |

The lowest-to-highest chi-squared order is CPL, wCDM, LCDM at every grid value. CPL retains w0=-0.3 at the upper box boundary throughout; radiation does not remove the pre-existing domain sensitivity. The largest absolute parameter changes at a refitted optimum are |delta alpha|=0.004926, |delta Omega_m|=0.00061293, |delta w0|=0.00125213, and |delta wa|=0.00364184 (maxima can belong to different models). Alpha stays free; these are not H0 measurements.

After refitting, the maximum change of any predicted datum is 0.0111742 marginal measurement errors, and the largest full-covariance displacement norm is 0.0190223. The maximum absolute profile chi-squared change is 0.0237559. Before refitting the largest norm was 0.127615, already below unit measurement precision. Thus radiation omission is small for this background profile screen at the tested optima and reference assumptions. This does not upgrade it to the official physical-model inference, a posterior comparison, or a precision massive-neutrino treatment.

Independent checks agree with every saved prediction within 4.98e-14, chi-squared within 2.10e-13, and analytic amplitude within 1.78e-14. Fixed-parameter shift vectors agree within 1.43e-14. The worker's GL96/GL192/adaptive-quadrature checks differ by at most 1.21e-13 in predictions and 2.03e-13 in chi-squared. This numerical error is far smaller than the radiation correction.

The worker used 16-start bounded L-BFGS-B with a separate differential-evolution search for every fit, plus Powell checks on active faces. Nineteen of 528 L-BFGS-B start attempts have failed/abnormal status and are retained in the result; they are not silently counted as successful fits. All 33 differential-evolution checks and all 11 active-face checks succeeded. The largest discrepancy between the selected radiation optimum and its independent differential-evolution check is 3.10e-9 in chi-squared. Feasible inward steps from the CPL boundary increase the objective. These are reproducible optimizer checks, not a global mathematical certificate.

Main scan (6.835 seconds measured scientific runtime, CPU FP64, one numerical-library thread):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --seconds 1200 -- .venv/bin/python experiments/radiation_sensitivity/radiation_scan.py
```

Independent saved-fit review, including all 41 fixed-point H0 values from the theory check (0.027 seconds measured scientific runtime, approximately 0.14 seconds wrapped process):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ./.venv/bin/python scripts/run_bounded.py --seconds 120 -- ./.venv/bin/python work/orchestrator2_radiation/review_saved_fits.py
```

`saved_fit_review.json` records the exact checked result/code hashes. Main code SHA-256: `b825ee41df8296c16bb332535acb2ae0b33eadf15b51227328d067292f109e21`. Main result SHA-256: `07a99215b4998e73a39c4d99c22bbe04e61cc9fa87ae86b161bb78ce8b842e2b`. `scripts/background_bao.py` retains its pre-lane hash.

## Source and derivation review

The theory worker's [report](../theory2/report.md) supplies the first-order response and exact source locators. Differentiating the flat-closure model gives `dE^2/dOmega_r=(1+z)^4-f_DE(z)`; the subtracted dark-energy term is essential. The integral derivative is `dI/dOmega_r=-0.5*integral[Q/E^3]dz`, yielding the reported DM, DH and DV responses. These equations agree with the independent finite-radiation calculations. At H0=50, exact-minus-linear prediction shifts are below 9.8e-5 marginal measurement errors. The theory result has detailed endpoint vectors at H0=50 and90 plus extrema over 41 fixed-point H0 values (50 through90 in steps of1). The main refits use the ten-point grid stated above. The independent saved-fit reviewer also recomputes all 41 theory values and reproduces their reported extrema within 3.20e-14, including the wCDM fixed-point absolute chi-squared change peaking at H0=69. This distinction avoids treating the fixed-point grid as 41 independent refits.

The lane reviewer also re-opened and checked these primary-source details on 2026-09-24:

- [CODATA 2022 recommended constants, Table XXXII, printed p.44](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf): c, Planck's constant, and G; the numerical assumptions match the calculation. No uncertainty sampling of constants is performed.
- [Fixsen 2009, arXiv v2 abstract](https://arxiv.org/abs/0911.1955v2): the reported temperature combination is 2.72548 +/-0.00057 K, consistent with using the rounded reference input 2.7255 K.
- [IAU 2015 Resolution B2, note4, PDF p.5](https://iauarchive.eso.org/static/resolutions/IAU2015_English.pdf): exact parsec convention `648000/pi AU`. The original `www.iau.org/static/resolutions/IAU2015_English.pdf` link failed in the review browser; the official archive link works.
- [PDG 2024 Neutrinos in Cosmology, Eq.(26.1), printed p.3](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-neutrinos-in-cosmology.pdf): the effective relativistic density normalization and its nominal temperature-ratio convention. N_eff=3.046 is the declared conventional reference used here; newer decoupling values are not silently substituted. This equation does not treat all physical neutrino mass states as radiation at low redshift.
- Local DESI DR2 v3 paper, printed pp.4-5, Eqs.(4)-(10): distances, CPL evolution, and separate physical neutrino term as already described above. The archived data audit pins the compressed BAO files to the exact release commit.

Source access here is targeted equation/constant review, not a claim to have read every cited paper in full. The source-derived definitions and independent numerical reproduction support the stated conditional screen result. Root owns any promotion to the shared source map, ledger, or morning report.

# Current interaction literature check

**Status:** Literature-only. I checked each public arXiv abstract/version page, downloaded and read both full v1 PDFs, and recorded local PDF SHA-256 values below and in crosscheck.json. No fit or result from either paper was independently reproduced.

## Wang, Yu & Wu, arXiv:2609.03062v1

[Reassessing Evidence for Dark-Sector Interactions with Dynamical Dark Energy and DESI DR2](https://arxiv.org/abs/2609.03062), submitted 2 Sep 2026. Local copy: arXiv_2609.03062v1.pdf, SHA-256 0e2b97a2fd988a6657e0c563bbbde7843529dea9c9ab8ac63dff88fefde356e5. Public full-text v1; 17 pages.

The paper uses a flat model with
\[
Q=\beta H\rho_{\rm de},\qquad
\dot\rho_c+3H\rho_c=-Q,\quad
\dot\rho_{\rm de}+3H(1+w)\rho_{\rm de}=Q,
\]
and \(Q^\nu=Qu_c^\nu\), so there is no momentum transfer in the CDM rest frame. Thus \(\beta>0\) transfers energy from CDM to dark energy. They evolve perturbations with interacting generalized PPF, which permits crossing \(w=-1\). Model and perturbation equations are on printed pp. 3–4 (Eqs. 1–10).

Their primary fit combines Planck low-\(\ell\) temperature/polarization and restricted high-\(\ell\) Plik-lite spectra, compressed ACT DR6 bandpowers, ACT+Planck lensing, DESI DR2 BAO, and DES-Dovekie SNe (printed pp. 4–5). They use a modified CAMB with Cobaya MCMC and GetDist. Comparison is by MAP \(\Delta\chi^2\) and AIC, not Bayes evidence (Eqs. 11–12, printed p. 4). The Dovekie description identifies the recalibrated DES five-year sample and its photometric/light-curve updates, but does not spell out the SN intercept or nuisance/covariance treatment.

For the CPL fit, \(Q=\beta H\rho_{\rm de}\) gives \(\beta=-0.35^{+0.77}_{-0.62}\), consistent with zero; \(w_0=-0.93^{+0.24}_{-0.19}\), \(w_a=-0.60\pm0.21\). The noninteracting CPL and interacting CPL MAP improvements over \(\Lambda\)CDM are \(\Delta\chi^2=-10.72\) and \(-10.75\): adding the coupling gains only 0.03 in MAP \(\chi^2\). Their \(\Delta{\rm AIC}\) values versus \(\Lambda\)CDM are \(-6.72\) and \(-4.75\), respectively (Table 1, printed p. 5). The interaction therefore adds no useful fit improvement once CPL freedom is present, while the dynamical-DE preference remains.

This is parameterization-dependent: relative to each corresponding noninteracting trajectory, the interaction changes MAP \(\chi^2\) by \(-6.07\) for thawing (\(\beta=+0.53^{+0.24}_{-0.13}\)), \(-0.25\) for mirage (\(\beta=-0.034\pm0.074\)), and \(-8.86\) for GEDE (\(\beta=-0.89^{+0.11}_{-0.19}\)); Table 5, printed p. 13. These are AIC/MAP comparisons, not evidence ratios. The paper cautions that its plotted growth spectra use representative parameter sets, not posterior propagation, and calls for direct low-redshift growth tests (printed pp. 6–15).

## Yang et al., arXiv:2609.05660v1

[Do DESI-DR2 BAO data imply a coupling of dark matter and dark energy?](https://arxiv.org/abs/2609.05660), submitted 4 Sep 2026; accepted for *Physical Review D*. Local copy: arXiv_2609.05660v1.pdf, SHA-256 0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c. Public full-text v1; 25 pages.

The interaction is
\[
Q=\Gamma\rho_x,\qquad
\dot\rho_x+3H(1+w_x)\rho_x=Q,\quad
\dot\rho_c+3H\rho_c=-Q.
\]
Here \(\Gamma>0\) transfers CDM energy to DE; their sampled parameter is \(\Gamma/H_0\). They assume flat FLRW, separately conserved baryons/radiation/neutrinos, constant \(w_x\), and fix \(\sum m_\nu=0.06\) eV and \(N_{\rm eff}=3.044\) (printed pp. 2–3, Eqs. 2–13). Perturbations use synchronous gauge, zero anisotropic stress, and no momentum transfer in the DM frame. For \(w_x\ne-1\) they evolve DE perturbations and fix \(c_s^2=1\) in the inference; IVS (\(w_x=-1\)) has no DE perturbation contribution (Eqs. 14–21, printed pp. 4–5). Stable phantom/quintessence branches use \(\Gamma/H_0<0\) and \(>0\), respectively.

The sampled likelihoods are Planck 2018 plikTTTEEE+lowl+lowE, DESI DR2 BAO (BGS, ELG, LRG, QSO and Ly-\(\alpha\) auto/cross measurements), and one of PantheonPlus (1701 light curves/1550 SNe), Union3 (2087 SNe), or DES-Dovekie (printed pp. 8–10). They use modified CAMB, Cobaya MCMC, Gelman–Rubin convergence checks, GetDist, and MCEvidence for \(\ln B_{i\Lambda}\), with uniform priors listed in Table I: \(\Gamma/H_0\in[-3,3]\) for IVS, \([-3,0]\) for phantom, and \([0,3]\) for quintessence; the corresponding constant-\(w_x\) ranges are \([-3,-1]\) and \([-1,0]\) (printed pp. 9–11). The article names the Dovekie sample but does not detail its SN intercept, covariance, or calibration-nuisance treatment. Its listed likelihoods do not include a low-redshift RSD or weak-lensing likelihood.

For the directly relevant CMB+DESI+DES-Dovekie IVS case, Table II (printed p. 11) gives \(\Gamma/H_0\simeq-0.053\): zero is outside the 68% interval but inside the 95% interval. The reported \(\ln B_{i\Lambda}=-4.4\) favors \(\Lambda\)CDM (moderate evidence on the paper's adopted Jeffreys scale). For the phantom branch, Table III gives roughly \(\Gamma/H_0=-0.075\) for the Dovekie combination and \(\ln B_{i\Lambda}=-7.5\); for quintessence, the coupling remains an upper limit and \(\ln B_{i\Lambda}=-8.1\) (Tables III–IV, printed p. 12). Their overall conclusion is that coupling indications depend on the assumed \(w_x\), interaction, and dataset, while Bayesian evidence favors \(\Lambda\)CDM for all tested combinations (printed pp. 16–17). This is not a CPL dynamical-\(w\) test.

## Relation to the current Dovekie-only profile

The local screen is a reduced, flat late-time vacuum model with \(Q=\Gamma\rho_{\rm vac}\), fit to the 1,820-row full STAT+SYS Dovekie covariance while profiling one shared additive magnitude intercept. Its exploratory minimum is \((\Omega_{m0},g)=(0.3930,-0.6294)\), \(g=\Gamma/H_0\), \(\chi^2_{\rm prof}=1629.9992\); the best profiled \(g=0\) reference is \(\chi^2=1631.4205\), so the descriptive minimum-to-null difference is 1.4213. It is a finite-domain profile, not a posterior, evidence, significance, or a joint fit; it omits radiation, perturbations, CMB, and BAO. The profile and limits are in experiments/dovekie_ivs_profile/result.json and work/theory7/dovekie_ivs_profile_contract.md.

For \(w_x=-1\), the local screen and Yang et al. use the same \(Q=\Gamma\rho_x\) law and sign convention, but their likelihoods and physical model scope differ substantially: the paper adds Planck and DESI and fits full cosmology, while the screen is Dovekie-only and late-time reduced. Its negative-\(g\) best-fit therefore neither conflicts with the paper's negative-\(g\), 68%-only posterior indication nor establishes the interaction; the paper's \(\ln B\) favors \(\Lambda\)CDM. Wang et al. use the different \(Q=\beta H\rho_{\rm de}\) law and add CPL freedom; their result directly supports the claim that interaction preference can disappear under dynamical-DE freedom.

**Current verdict:** Both full texts reinforce model- and dataset-dependence. The local result is a shallow exploratory Dovekie-only shape preference, not literature-confirmed evidence for an interaction. Neither paper is a like-for-like validation of that profile.

**Implication for the next test:** Reproduce the closest published comparison first: a full Planck 2018 + DESI DR2 + DES-Dovekie IVS inference with \(Q=\Gamma\rho_x\), the paper's \(\Gamma/H_0\in[-3,3]\) prior, and the same Dovekie likelihood. Report posterior constraints and \(\ln B\). This will show how the reduced SN-only optimum moves when the external likelihoods are added.

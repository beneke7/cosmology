# vEDE-at-deuterium-burning physicalness review

**Finding.** This is a coherent phenomenological demonstration that a transient increase in H during late deuterium burning can raise surviving D/H with little effect on Yp. It does not establish a BBN-only preference for new physics or a complete scalar-field cosmology. The headline preference is conditional on the PRIMAT rate prescription and a compressed, high-omega_b EDE CMB target. Treat it as a useful conditional mechanism and hypothesis-generating calculation, not evidence that a second scalar field is present.

**Status:** equations and likelihood algebra independently checked; source/model limitations literature-only (the exact preprint, 2026 data-driven-rate paper, and PDG 2025); exact paper-level numerical reproduction not performed.

## Equation and likelihood check

For Eq. (7), w_n=(n-1)/(n+1), so n=3 gives w=1/2 and alpha=3(1+w)=9/2. The interpolant

    rho_v(T) = 2 rho_v(Tc) / [1 + (Tc/T)^alpha]

is continuous at T=Tc: it tends to 2 rho_v(Tc) on the frozen high-temperature side and falls as T^alpha proportional to a^(-alpha) on the post-transition side. Eq. (8) is the usual fractional density at Tc, and Eq. (9) is the flat-GR Friedmann relation for the stipulated standard plus vEDE densities. These are dimensionally consistent.

Eqs. (10)–(12) also map correctly. With H_std^2=8 pi G rho_std/3,

    rho_v(TD) = [(1+delta_D)^2 - 1] rho_std(TD)
    rho_v(Tc) = 0.5 rho_v(TD) [1 + (Tc/TD)^alpha].

Thus delta_D is dimensionless and rho_v is nonnegative on the stated prior 0 <= delta_D <= 0.20. Eq. (13)'s conversion is sound: kT=0.0293 MeV at 3.4e8 K, and TD/T0-1 is about 1.25e8. Caveat: Tc=T0(1+zc) assumes post-e+-annihilation photon-temperature scaling, while the flat prior 7<log10(zc)<10 reaches Tc about 2.35 MeV, before that scaling applies. The reported combined posterior's log10(zc)<8.3 (95%) is after annihilation, but prior support at high z should not be read as a physically calibrated transition redshift without entropy-history mapping.

Eq. (14)'s independent-Gaussian abundance likelihood has the correct -2 log-likelihood:

    chi2_BBN = sum_X [(X-Xobs)^2/(sigma_obs,X^2 + sigma_th,X(theta)^2)
                      + ln(2 pi [sigma_obs,X^2 + sigma_th,X(theta)^2])],

for X=Yp and 10^5 D/H. Since the theory variance depends on parameters, retaining the Gaussian normalization is necessary; the text explicitly says it is retained. I found no algebraic normalization error. Limits: PRIMAT Monte Carlo errors are compressed to abundance variances rather than sampled rate nuisances, and the two abundance likelihoods omit any theory covariance between Yp and D/H. The source says it uses independent observational Gaussians and tabulated PRIMAT rate/lifetime errors (pp. 5–6, Eq. 14).

The sampled priors are uniform 0.0205<omega_b<0.0240, 0<delta_D<0.20, and 7<log10(zc)<10. Sampling delta_D is sensible for the BBN response; the mapping to f_vEDE(zc) is correctly derived but nonlinear. Appendix C finds f_vEDE essentially unconstrained, log10(zc)<8.3 one-sided at 95%, and posterior pile-up near f~1 from broad prior volume. The abundances therefore constrain local H(TD), not a unique peak scalar fraction or transition history.

## Inference scope and fragility

- The BBN-only profile's vEDE best fit is at the no-boost boundary, delta_D=0, as for Lambda-CDM; profiling over delta_D and zc makes its high-omega_b profile nearly flat (Fig. 2 and Table I). This shows compatibility can be restored, not that BBN alone detects a departure.
- The quoted delta_D=0.087 (+0.036/-0.037) and Q_DMAP=0.7 sigma are from the **BBN plus compressed EDE-omega_b Gaussian** fit, omega_b=0.02272 +/- 0.00014, not the BBN-only profile. Q_DMAP is a dataset-compatibility statistic, not model selection. Separately, delta-AIC about -5.3 counts two added parameters for that same combined compressed likelihood; it is not a detection significance or Bayesian evidence. The Gaussian is reused for all three BBN models as a proxy; the authors state a full joint EDE+vEDE fit to CMB data is beyond scope (pp. 5–7).
- The D/H source matters. With updated 10^5 D/H=2.508 +/- 0.030, Appendix A reports vEDE delta_D=0.073 (+0.035/-0.036), Q_DMAP about 0.7 sigma, and delta-AIC about -3.0, versus -5.3 using Cooke 2018's 2.527 +/- 0.030. The lower value reduces the Lambda-CDM discrepancy to about 2.7 sigma. PDG 2025 independently recommends 10^6 D/H=25.08 +/- 0.29 (Table 24.1 mean of 12 systems, with scatter inflation), close to the Appendix A value. The updated-source branch deserves primary weight over the older baseline.
- The EDE target matters too. Appendix B's lower Lambda-CDM-SPA target omega_b=0.022398 +/- 0.000095 makes vEDE nearly neutral in AIC (delta-AIC=+0.02); SPA+DESI gives only mild preference (-1.5). The stronger -5.3 result is specific to the higher EDE omega_b target, not generic support for BBN-era vEDE.
- Nuclear source dependence is load-bearing. The v1 paper itself says PArthENoPE and Ref. [3] predict D/H closer to observation than PRIMAT at comparable baryon density, attributing the offset to D(D,n)3He and D(D,p)3H rate conversions (pp. 2–3). The independent 2026 GP rate analysis, arXiv:2604.16600v1, finds 10^5 D/H=2.442 +/- 0.040 with Planck omega_b and 2.437 +/- 0.034 with SPA (Table I); it identifies those two D+D S-factors as needing measurement in the 0.1–0.6 MeV range. It does not reproduce the EDE-target fit, but confirms nuclear-rate methodology remains a leading theory systematic. The vEDE likelihood propagates PRIMAT's own tabulated errors, not between-code/rate-method discrepancy.

The vEDE density is an axion-motivated fluid interpolation, not a solved scalar field. No field equation, initial condition, mass-to-zc mapping, or perturbation evolution is fitted. The background-only BBN calculation can demonstrate the chemistry response under prescribed H(T), but it does not establish a field realization. The authors acknowledge possible small-scale transfer-function effects and omit them; the full CMB+BBN analysis is unperformed (pp. 5–6, 9–10). This is a completion gate, not a minor refinement.

**SBBN prior caveat for the separate IVS lane:** the PDG BBN baryon prior presumes standard expansion during nucleosynthesis. The current IVS equations conserve baryons, photons, and neutrinos separately and exchange energy only between CDM and vacuum. Under the source convention dot(rho_v)=Gamma*rho_v, the reported negative coupling means vacuum density increases backward; the early-time radiation-era contribution is still asymptotically subdominant. Whether the resulting H_IVS/H_std shift is negligible over BBN was not quantified in the earlier screen, so the SBBN prior is conditional pending that check; otherwise use an IVS-conditioned abundance likelihood. The vEDE departure studied here is explicitly not in that negligible-departure case.

## Reproducibility and next gate

The exact local input is arXiv:2607.20635v1 (22 July 2026). Its final paragraph says code, data, and notebooks needed to reproduce figures **will be public upon acceptance** at Zenodo record 21494819. The arXiv record still lists v1 only (checked 25 September 2026); I found no exact paper analysis bundle in the local manifest/cache, and the linked Zenodo record was not accessible through read-only lookup. The paper links the separately public PRIMAT and primat_tools projects, but those do not establish availability of the paper-specific vEDE modification, chains, or figure data. Exact figure/likelihood reproduction is therefore **not verified/available from checked materials**; do not call it reproduced.

**Recommended gate:** before implementing the scalar/CMB model, evaluate the same updated-D/H likelihood with PRIMAT and independent GP/current PArthENoPE rate alternatives, and with both EDE and lower Lambda-CDM-SPA omega_b targets, retaining identifiable rate covariance. If a material high-omega_b residual survives, the discriminating experimental priority is improved D(D,n)3He and D(D,p)3H S-factors over 0.1–0.6 MeV, as the 2026 rate paper recommends. Then replace the fluid with an explicit scalar solution and test perturbations, transfer functions, and a full joint CMB+BBN likelihood.

**Source locators:** Poulin et al., arXiv:2607.20635v1, Eqs. (7)–(14), (16)–(20), Figs. 1–2, Tables I–III, Appendices A–C; https://arxiv.org/abs/2607.20635. Launders et al., arXiv:2604.16600v1, Eqs. (4)–(12), Table I; https://arxiv.org/abs/2604.16600. PDG 2025, *Big-Bang Nucleosynthesis*, Secs. 24.2–24.4, especially Table 24.1 / Eq. (24.2); local exact PDF/text at context/papers/pdg_2025_big_bang_nucleosynthesis.*.

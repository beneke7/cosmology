# Interacting-vacuum perturbation/physicality audit

**Status:** advisory source/equation audit, 2026-09-24; no solver or fit run. This does not revise the campaign’s BAO-only result.

## Bottom line

The BAO ansatz is not, by itself, a covariant perturbation theory. However, the cited current paper does specify an operational linear-fluid closure for its interacting-vacuum case: synchronous gauge, no momentum transfer in the CDM rest frame, and the resulting CDM perturbation equations. This is the geodesic-CDM closure. Under that choice the vacuum is spatially homogeneous on CDM-comoving slices, not unperturbed in every gauge. The paper does not give a dark-sector action/microphysical origin, and its published text does not spell out the primordial initial-condition prescription or a versioned CAMB patch. Thus a CMB/growth calculation is implementable as a *phenomenological model*, but reproducing one requires explicit code/configuration choices and early-time physicality checks.

The exact current source is Yang et al., [arXiv:2609.05660v1](https://arxiv.org/abs/2609.05660v1), accepted for PRD. Its local PDF is `context/papers/interacting_de_desi_dr2_2026.pdf`, SHA-256 `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`. I read the PDF/text and visually checked the printed perturbation equations on pp. 4–5. Equation locators below refer to the PDF numbering.

## What is fixed, and what a Boltzmann implementation must state

| Ingredient | Source-paper status | Reproducible implementation contract |
|---|---|---|
| Background split and signs | Eqs. (3)–(5), (8), pp. 2–3: baryons/radiation/neutrinos evolve separately; \(\dot V=Q\), \(\dot\rho_c+3H\rho_c=-Q\), \(Q=\Gamma V\). Positive \(\Gamma\) is CDM→vacuum. | Keep this sign in every background and perturbation source. \(\Gamma\) has proper-time units; use \(g=\Gamma/H_0\) only as a dimensionless parameter. Specify whether the covariant scalar is the local \(V(x)\), not merely the FLRW background value. |
| Transfer four-vector | The paper states “no momentum transfer in the DM rest frame” before its fluid equations and supplies Eqs. (19)–(20), pp. 4–5, for IVS via Ref. [141]. It does not print a standalone covariant \(Q^\mu\) equation for this exact ansatz or derive it from an action. | The standard completion consistent with its signs is \(\nabla_\mu T_c^{\mu\nu}=-Q u_c^\nu\), \(\nabla_\mu T_V^{\mu\nu}=+Q u_c^\nu\), with \(T_V^{\mu\nu}=-Vg^{\mu\nu}\), \(Q=\Gamma V\). This is pure energy transfer in the CDM frame: CDM is geodesic. Because vacuum stress implies \(-\nabla^\nu V=Q u_c^\nu\), its spatial gradient vanishes in that frame. Martinelli et al. note this geodesic/potential-flow approximation is suited to linear scalar perturbations; it is not a general nonlinear halo closure. Treat this as the closure implied by the paper’s equations, not a unique consequence of the background ansatz. |
| Vacuum perturbation | The paper says IVS has no DE contribution at perturbative level and omits a vacuum density-contrast equation. | Set \(\delta V_{(c)}=0\) on CDM-comoving slices; if \(Q=\Gamma V\) locally, then \(\delta Q_{(c)}=0\) there. Do not state that \(\delta V=0\) in every gauge. Martinelli et al. explicitly note it is homogeneous only in the geodesic-CDM frame [arXiv:1902.10694v2, §2.3](https://arxiv.org/abs/1902.10694). |
| Gauge/variables | Synchronous-gauge metric (Eq. 14), Fourier-space \(\delta_c,\theta_c,h\), conformal-time primes and \(\mathcal H=a'/a\); CDM Euler equation is \(\theta_c'=-\mathcal H\theta_c\), Eqs. (19)–(20). | Fix the residual synchronous gauge by choosing CDM-comoving \(\theta_c=0\) (the geodesic closure preserves it). Then evolve \(\delta_c'=-h'/2+a\Gamma(V/\rho_c)\delta_c\), together with Einstein constraints and the standard species. Check gauge-invariant observables independently; do not compare raw gauge-dependent density contrasts across codes. |
| Radiation/neutrinos and calibration | Eq. (3), p. 2, includes photons/radiation, baryons, CDM, vacuum and neutrinos. The paper fixes \(\sum m_\nu=0.06\,\mathrm{eV}\), \(N_{\rm eff}=3.044\); it analyzes Planck 2018 TT/TE/EE+low-\(\ell\)/lowE and DESI DR2. §III says a modified CAMB is used. | For an independent CMB run, record \(T_{\rm CMB}\), neutrino mass/species split and hierarchy, \(N_{\rm eff}\), \(\omega_b,\omega_c,H_0\), helium/BBN and recombination settings, primordial \(A_s,n_s\), reionization, and the exact likelihood/config/code version. Retain photon/baryon evolution, neutrino hierarchy/free-streaming shear and metric feedback. A low-z BAO amplitude \(c/(H_0r_d)\) cannot supply these calibrations. |
| Initial conditions | The paper’s prose gives the modified-CAMB/MCMC method and standard parameter set (including \(A_s,n_s\)) but no explicit IVS superhorizon initial-condition formula, start epoch, isocurvature choice, or vacuum-mode initialization in the accessible article text. | State and test a regular growing adiabatic initial mode for photons, baryons, CDM and neutrinos, and explicitly match it to \(\delta V_{(c)}=0\); specify whether any dark-sector isocurvature mode is excluded. Verify convergence with starting redshift and \(k\tau\), and check that the chosen mode satisfies the modified Einstein/continuity constraints. Do not silently assume a vacuum sound-speed prescription. |

The source’s \(w_x\ne-1\) fluid branches are distinct: it quotes large-scale stable regions \(w_x<-1,\Gamma/H_0<0\) (phantom) and \(w_x>-1,\Gamma/H_0>0\) (quintessence), and commonly fixes their rest-frame \(c_s^2=1\). Those sound-speed/sign rules are not a stability theorem for IVS. At exactly \(w_x=-1\), the source’s dark-energy fluid perturbations are removed under its chosen vacuum closure; there is no independent vacuum-fluid \(c_s^2\) to set. The paper’s uniform IVS scan range is \(\Gamma/H_0\in[-3,3]\) (Table I), without a stated IVS sign cut. This does not establish ghost/gradient stability of a microphysical theory: no such action is supplied.

## Exact sign and positivity checks

For the source convention and \(V_0=\rho_{x0}>0\),

\[
\dot V=\Gamma V,\qquad V(t)=V_0 e^{\Gamma(t-t_0)}>0,
\]

and

\[
\frac{d}{dt}(a^3\rho_c)=-\Gamma a^3V,\qquad
a^3(t)\rho_c(t)=\rho_{c0}-\Gamma\!\int_{t_0}^{t}a^3(t')V(t')\,dt'.
\]

These are exact consequences of Eqs. (4), (5), (8), not numerical fits. For a past-time interval \(t<t_0\), positive \(\Gamma\) makes the comoving CDM density larger than today and hence positive if \(\rho_{c0}>0\). Negative \(\Gamma\) has a finite past-positivity condition:

\[
\rho_{c0}>|\Gamma|\int_{t}^{t_0}a^3V\,dt'.
\]

For a future interval, the corresponding possible failure is for positive \(\Gamma\): require \(\rho_{c0}>\Gamma\int_{t_0}^{t}a^3V\,dt'\). Vacuum positivity alone therefore does not guarantee physical CDM. In the linear source equation, \(\rho_c\) is also in the denominator, so approaching zero is both a physical boundary and a perturbation-equation singularity. Check \(\rho_c>0\), \(V>0\), expanding \(H^2>0\), and regular finite perturbations over the *entire* integration range, not just observed BAO redshifts.

The signs are convention-specific. For example, the older [Yang et al. 2020 paper, arXiv:2001.04307](https://arxiv.org/abs/2001.04307), uses the same-looking \(Q=\Gamma\rho_x\) rate but describes positive \(\Gamma\) as vacuum→CDM in its Model II (PDF pp. 3–4, Eqs. 13–14); its CDM perturbation source consequently has the opposite sign. It also states energy flow parallel to CDM four-velocity and vacuum perturbations vanish in that frame. Convert the continuity-equation convention before reusing any sign or source term.

## Current-literature boundary

The exact source above (submitted 2026-09-04, currently arXiv v1 and accepted for PRD) is the directly matched recent analysis. The covariant explanation in [Martinelli et al., arXiv:1902.10694v2](https://arxiv.org/abs/1902.10694) is a useful primary derivation of geodesic-CDM vacuum perturbations, but its fitted coupling is \(Q\propto H V\), **not** the present constant proper-time \(Q=\Gamma V\). Likewise, the recent nonlinear study [Zhai et al., arXiv:2606.11368](https://arxiv.org/abs/2606.11368) uses \(Q=\xi\mathcal H\rho_x\); its abstract reports model-dependent nonlinear departures and warns against assuming \(\Lambda\)CDM-calibrated nonlinear prescriptions, but it does not validate this exact \(\Gamma V\) model. Its full text was not read in this bounded lane.

## Campaign implication / next physicality gate

The campaign profile remains a low-redshift background screen (preferred \(g\simeq-0.467\), checked only through \(z=2.33\)); its existing report explicitly omits radiation, neutrinos, perturbations and sound-horizon calibration. No sign, fit score or campaign claim is changed here. Before giving that point physical/CMB meaning, the next gate should be a small, versioned linear-Boltzmann reproduction of the exact geodesic closure above: first test \(\Gamma=0\) recovery and the sign in a two-component perturbation check; then integrate early enough to verify CDM positivity, regular adiabatic initial conditions and tolerance/start-redshift convergence with the stated radiation/neutrino sector. Only if it passes should a full CMB+growth likelihood or null calibration be considered.

See [equation/source checklist](CHECKLIST.json) for machine-readable locators, provenance and gates.

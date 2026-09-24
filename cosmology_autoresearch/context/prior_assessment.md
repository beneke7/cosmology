**Minimal cosmology, thermodynamic gravity, and an agentic research programme**

Research brief for Benedek. Literature checked on 23 September 2026. This is a targeted assessment of observational tensions, minimal extensions, and computational opportunities, with an independent equation-level audit of the supplied Pszota–Ván paper. It is not an exhaustive review, a new observational fit, or a demonstrated discovery of new physics.

**The most promising starting point is a predictive test of the supplied gravity model.** A separate, more ambitious cosmology project could ask which minimal changes to expansion and gravity remain identifiable after accounting for observational systematics. Both are feasible on a strong workstation. The scientific bottleneck is distinguishing explanations that fit the same data, and demonstrating a gain on data excluded from model development.

The premise should remain open: ΛCDM faces serious tensions, but has not been decisively displaced by an alternative that explains the complete evidence better. Evolving dark energy, changes to dark matter, and modified gravity are distinct hypotheses. Evidence against a constant dark-energy equation of state is not evidence that dark matter is unnecessary.

**The observational situation has changed since the original DESI headlines.** Numbers below belong to specific combinations and statistical definitions. They should not be combined as independent significance measurements.

| Evidence | Result checked | Research implication |
|---|---|---|
| Primary CMB | ACT DR6 temperature and polarization remain well described by ΛCDM. Its joint analysis with Planck, lensing and DESI DR2 gives H0 = 68.43 ± 0.27 km/s/Mpc. [1] | The early-universe prediction machinery remains a demanding baseline. A successful alternative must preserve much more than a distance–redshift curve. |
| Local expansion | The H0 Distance Network reports 73.50 ± 0.81 km/s/Mpc, with covariance between distance indicators included. The collaboration reports 7.1σ relative to its chosen combined-CMB ΛCDM comparison, and 5.0σ relative to BBN + DESI BAO. [2] | A serious discrepancy. Its interpretation still depends on calibration, likelihood and cosmological assumptions; a discrepancy is not an identification of its cause. |
| Original DESI DR2 | The 2025 BAO paper reports a 2.8–4.2σ preference for evolving dark energy when CMB and different supernova compilations are included. [3] | A reason to test simple extensions, with explicit sample and prior dependence. |
| Recalibrated supernovae | DES-Dovekie reduces the earlier 4.2σ preference to 3.2σ, and reports approximately 5:1 Bayesian model preference under its assumptions. [4] | Calibration can materially move the apparent evidence. Bayesian model preference and a frequentist significance answer different questions. |
| July 2026 DESI update | Lyα full-shape Alcock–Paczyński information shifts the high-redshift anchor toward ΛCDM. The joint analyses report 2.7σ for DESI + CMB and 3.2σ with supernovae. [5] | The data have not moved monotonically toward rejection of ΛCDM. This is a valuable new discriminator for proposed late-time deformations. |
| September 2026 supernova update | The Unite preprint combines 2,884 likely SNe Ia. It reports 3.3σ using a MAP comparison, 3.1σ using maximum likelihood, and weak Bayesian preference for time variation. Its likelihood/data are promised upon acceptance, not currently assumed available. [6] | Include the result in the literature assessment, but use released likelihoods for an immediately executable project. |
| Structure growth | KiDS-Legacy gives S8 = 0.815(+0.016/−0.021), consistent with Planck at 0.73σ. DES Y6 3×2pt gives S8 = 0.789 ± 0.012; its difference from combined CMB is 2.6σ in S8 and 1.8σ in the full parameter space. [7,8] | There is no uniform, survey-independent verdict from “the S8 tension.” Photometric redshifts, intrinsic alignments, baryonic feedback and parameter projection matter. |

These examples also explain why extra compute is not a shortcut to certainty. A larger search can find increasingly convincing fits to a calibration error or a statistical fluctuation unless the entire search procedure is tested on mock observations.

**Occam’s razor is useful, provided complexity is counted at the level of predictions.** The standard six-parameter model contains strong structural assumptions: FLRW geometry, general relativity, a cosmological constant, cold collisionless dark matter, a simple primordial spectrum, and specified neutrino/recombination physics. A short equation that introduces an arbitrary function, a boundary prescription for each galaxy, or independent nuisance parameters for each object need not be simpler in the statistical sense.

For a family M, the evidence is

\[
p(D\mid M)=\int p(D\mid\theta,M)\,p(\theta\mid M)\,d\theta.
\]

Report best fit, predictive performance, evidence with prior sensitivity, and individual-dataset residuals. AIC/BIC are useful diagnostics, but cannot substitute for proper evidence and mock-calibrated model selection, especially near boundaries or after searching many functional forms. Broader error bars can reduce a quoted tension without increasing explanatory power.

Small changes are already a major research programme. “Small” should mean a controlled, physically specified family, not merely a small number of symbols. A correction may need to become order unity in the low-acceleration regime, even if it is negligible in the Solar System.

| Family | What is changed | Existing capability and limitation |
|---|---|---|
| wCDM and CPL | Constant w, or w(a) = w0 + wa(1−a) | Essential baselines. DESI has already compared parametrizations, bins and Gaussian-process reconstructions; its extended analysis finds two parameters sufficient to capture the principal trend. Repeating this with another unconstrained curve is weak novelty. [9] |
| Scalar-field dark energy | A specified field and potential | Canonical minimally coupled quintessence with positive energy density has w ≥ −1. A reconstructed crossing requires more structure or a different effective interpretation. Background fits must be accompanied by consistent perturbations. |
| Stable scalar–tensor extensions | Time-dependent gravitational coupling and scalar interactions | hi_class/EFTCAMB implement broad families. mochi_class uses a basis designed to enforce ghost/gradient stability. In June 2026, Cataneo and Koyama already scanned stable kinetic-gravity-braiding models with phantom crossing against several probes. [10,11,12] |
| Early dark energy | A transient pre-recombination component | Potentially changes the sound horizon and inferred H0. ACT DR6 does not independently prefer EDE; a 2025 combined analysis nevertheless finds it remains viable, with important prior-volume effects. This is active and contested, not a settled solution. [13] |
| Dark-sector interactions and neutrino changes | Couplings, radiation, masses, or decay | Good small-family comparisons when perturbations, BBN, CMB and growth are treated consistently. One attractive late-time fit is insufficient. |
| Self-interacting dark matter | Scattering and gravothermal halo evolution | A physically motivated galaxy-scale programme. A 2026 SPARC analysis already improves a semi-analytic SIDM Jeans model, so “fit SIDM to SPARC” alone is not an opening. Baryons and core-collapse modelling create substantive uncertainty. [14,15] |
| MOND/AQUAL/QUMOND and relativistic variants | The acceleration law and, in full theories, additional fields | Strong galaxy phenomenology deserves fair comparison. Modern relativistic MOND is more developed than a simple claim that “MOND cannot fit the CMB”: Skordis–Złośnik demonstrated agreement with linear CMB and matter spectra in their theory. This does not validate every nonrelativistic force law. [16] |
| Supplied thermodynamic gravity | A nonlinear correction to the potential equation | A concrete equation with a small fitting space, but the supplied work is a one-galaxy demonstration, not a complete replacement cosmology. [P] |

**The attachment provides a particularly tractable mathematical target.** Pszota and Ván study NGC 3198 with a spherical surrogate for the baryonic mass distribution. They fit K and a mass-to-light conversion, and set the two radial derivative boundaries using observed rotation velocities. They explicitly count four fitted/conditioned quantities in their reduced chi-squared. Their Table 1 reports reduced chi-squared 1.29 for TG, 1.26 for DC14 and 1.62 for MOND EFE, under different fitting treatments. Those numbers establish that the curve can be represented well; they do not establish superior out-of-sample predictions. [P, pp. 4–8]

The spherical conversion is also substantive. Disk geometry permits signed contributions from a gas distribution, whereas a positive spherical enclosed-mass model obeys different constraints. SPARC’s tabulated Newtonian baryonic velocities are not automatically the source density for a different field equation. An axisymmetric calculation needs a consistently reconstructed stellar and gas distribution, with vertical structure and deprojection uncertainty. Verify resolved gas-profile availability galaxy by galaxy; the tabulated gas circular-speed contribution alone is insufficient to assume a unique three-dimensional density. Restrict the pilot to objects with suitable source-survey data if necessary.

The paper allows K to differ among galaxies because it describes a constitutive coupling. That is a legitimate modelling choice. Its scientific value becomes much greater if K can be predicted from independently measured galaxy properties, with uncertainty, instead of fitted separately to every rotation curve.

**An exact change of variables substantially simplifies the PDE.** Starting with the printed equation, write D = l²/τ and S = 4πGρ:

\[
\partial_t\phi=D\left(\nabla^2\phi-K|\nabla\phi|^2-S\right).
\]

Assume K is constant in space and time, and the source density is prescribed. Set

\[
u=\exp(-K\phi)>0.
\]

Then the chain rule gives

\[
\nabla^2u=-Ku\nabla^2\phi+K^2u|\nabla\phi|^2,
\qquad
\partial_tu=-Ku\partial_t\phi,
\]

so that

\[
\boxed{\partial_tu=D(\nabla^2u+KSu)},
\qquad
\boxed{(\nabla^2+KS)u=0\quad\text{at stationarity}}.
\]

This is a linear equation for u, including nonspherical prescribed densities. It is a standard exponential/Cole–Hopf-type manipulation, independently derived here; no novelty claim is attached to the transformation. It does not linearize a fully coupled matter–gravity system, and a variable K introduces extra derivative terms.

The practical consequence is large: use a conventional sparse linear boundary-value solver as a reference before considering a neural PDE solver. Recover acceleration as g = ∇u/(Ku). Enforce u > 0, transform the physical boundary conditions, and test domain convergence. A Neumann condition on φ becomes a Robin condition on u. The operator is not automatically positive definite; near a zero eigenvalue it can become ill-conditioned. A formal change of variables does not establish existence, uniqueness or dynamical stability for every source and boundary choice. Treat K = 0 with the Poisson formulation or a controlled limiting expansion to avoid division and cancellation errors. Thermodynamic consistency also does not establish that a chosen constitutive closure is empirically correct.

For a perturbative starting point, φ = φN + Kφ1 + O(K²) gives

\[
\nabla^2\phi_N=4\pi G\rho,
\qquad
\nabla^2\phi_1=|\nabla\phi_N|^2.
\]

Thus even the first correction is a pair of Poisson solves. Its regime must be checked using dimensionless combinations such as KGM/R. Large galaxy mass discrepancies need not lie in that perturbative regime.

**There is a concrete sign-convention issue to reconcile.** The paper’s Eq. (8), with g = −φ′, is

\[
g'=-2g/r-S-Kg^2.
\]

For constant S > 0 and K > 0, its regular central expansion is

\[
g=-Sr/3-KS^2r^3/45+O(r^5).
\]

The printed Eq. (23) instead has a positive cubic coefficient and the hyperbolic-cotangent branch

\[
g_p=\frac1{Kr}-\sqrt{S/K}\,\coth(\sqrt{KS}\,r).
\]

Direct differentiation gives

\[
g_p'+2g_p/r+S+Kg_p^2=2Kg_p^2,
\]

rather than zero. The regular positive-K branch of Eq. (8) is

\[
g_c=-\frac1{Kr}+\sqrt{S/K}\,\cot(\sqrt{KS}\,r),
\]

on intervals before its poles. The paper itself notes an opposite-sign K comparison between numerical and analytic solutions in Section 2.2.3. This is therefore a documented convention/consistency problem to resolve explicitly, not an accusation of an undisclosed error. It does not, by itself, establish that the numerical galaxy fit is wrong. No author correspondence or implementation inspection was performed.

The supplied Python audit independently checks finite-difference residuals and integrates the nonlinear g equation and the transformed u equation. These are numerical checks, not a Lean proof and not astronomical evidence.

| Check actually executed | Result |
|---|---:|
| Printed coth branch residual, K = S = 1, 0.2 ≤ r ≤ 1.2 | Maximum absolute residual 0.2682 |
| Correct cot branch, same test | Maximum absolute residual 4.25 × 10⁻¹² |
| Numerical check that the printed-branch residual equals 2Kg² | Maximum difference 4.28 × 10⁻¹² |
| Independent nonlinear/transformed integrations, S(r) = exp(−r²), K = 0.4 | Maximum absolute acceleration difference 4.80 × 10⁻¹¹ |

**The transformed problem exposes a useful physical sensitivity test.** Consider an isolated uniform-density sphere of radius R, with a regular centre and φ → 0 at infinity. Let q² = KS and x = qR. On the positive branch connected to Newtonian gravity,

\[
u_\mathrm{in}(r)=\frac{\sin(qr)}{qr\cos x},
\qquad
u_\mathrm{out}(r)=1+\frac{B}{r},
\qquad
B=R\left(\frac{\tan x}{x}-1\right).
\]

Matching the far-field acceleration to −GMapp/r² gives

\[
\frac{M_\mathrm{app}}{M_b}
=\frac{3}{x^2}\left(\frac{\tan x}{x}-1\right),
\qquad 0<x<\pi/2.
\]

| x | Apparent-to-baryonic mass ratio |
|---:|---:|
| 0.2 | 1.016 |
| 0.5 | 1.111 |
| 1.0 | 1.672 |
| 1.3 | 3.144 |
| 1.5 | 11.201 |

This is an analytic toy example, not a realistic galaxy model. It shows why a viable benchmark should measure sensitivity to K, density and boundaries: substantial amplification can approach a critical branch. Whether actual galaxy solutions occupy such a sensitive regime is an empirical question. A spherical example does not prove fine-tuning for all disk galaxies.

In vacuum, the paper’s solution gives

\[
g=-\frac1{Kr+Cr^2},\qquad v_c^2=\frac1{K+Cr}.
\]

For K,C > 0, a flat portion can occur where Cr ≪ K; sufficiently far away the curve becomes Keplerian. It is not generically flat to infinite radius. When this flat-regime mechanism is invoked, vflat⁴ ≈ 1/K². Reproducing the baryonic Tully–Fisher scaling vflat⁴ ∝ Mb therefore demands an appropriate K–mass relation, roughly K ∝ Mb⁻¹/² in that regime. Fitting a different K to each galaxy does not by itself explain the relation.

**The first proposed experiment is a galaxy-level predictive benchmark.** The question is: can a single low-complexity constitutive prescription for K, combined with independently estimated baryons, predict rotation curves of galaxies excluded from model construction?

1. Freeze the equation, units, signs, source-density mapping, boundary prescription, quality cuts and scoring rule. Reproduce the paper as a conditional fit, clearly recording its observed endpoint conditions. Separately implement a genuinely predictive boundary-value problem with regularity and specified exterior/environment conditions.
2. Use the public SPARC mass models and source references. Begin with a small, prespecified subset spanning gas-rich dwarfs and high-surface-brightness disks. Validate the Newtonian limit and a spherical manufactured solution, then implement axisymmetry. Do not silently turn signed disk force contributions into negative spherical mass density.
3. Compare K = 0, individually fitted K, universal K, a simple mass-scaled law, and a hierarchical relation K = f(baryonic properties). Include intrinsic scatter. Infer the shared law through the forward model, propagating uncertainty in mass, distance, inclination and mass-to-light ratios; regressing noisy fitted K values on noisy masses can create spurious relations.
4. Use MOND with an explicitly stated interpolation and external-field treatment, and CDM halo models with documented population priors. Equal treatment of distance, inclination and stellar-population uncertainties matters more than forcing equal parameter counts. Include both flexible and population-constrained halo comparisons, because the former tests curve representation while the latter tests stronger predictions.
5. Hold out entire galaxies. Randomly splitting radial points leaks object-level information. Since NGC 3198 was already used to motivate the theory, treat it as development data. Keep a final sample unavailable to the proposing agent. Public-data reuse limits how strong a discovery claim can be, even with this split.
6. Evaluate log predictive density and calibrated intervals, along with residual dependence on mass, surface brightness, gas fraction, radius and environment. Record negative results, including the possibility that K cannot be predicted without near-free per-galaxy flexibility.
7. On mock observations, inject both CDM-like galaxy diversity and alternative-gravity signals through the same observational assumptions. Rerun the whole fitting and model-selection process. A signal should survive nuisance changes and improve prediction of new galaxies.
8. If this passes, test a new observable: outer rotation curves or vertical dynamics can be informative. Lensing requires a specified relativistic relation between metric potentials. The scalar nonrelativistic equation alone does not predict light deflection, CMB anisotropies or cosmic expansion.

The key novelty would be a fair geometry-aware, uncertainty-aware predictive comparison and a constrained K closure. “More SPARC fits” is not sufficient. I did not verify the absence of an existing full-sample TG analysis, and a new Pszota–Ván higher-gradient preprint, arXiv:2609.00317, appeared in September 2026. Its metadata was identified but its full text could not be retrieved in this session; it should be checked before fixing publication novelty.

**The cosmology experiment should search for identifiable modifications.** A useful core question is: which one- or two-mode changes to background expansion or gravitational response cannot be absorbed by allowed nuisance shifts, and predict an independent probe better?

For the background, start with w(z) = −1 + δw(z), and use the exact continuity relation

\[
\rho_\mathrm{de}(z)=\rho_{\mathrm{de},0}
\exp\left[3\int_0^z\frac{\delta w(z')}{1+z'}dz'\right].
\]

Fit actual BAO distance ratios and supernova distance data with their full released covariance. Keep ΛCDM, wCDM and CPL as baselines. A spline, Gaussian process or symbolic expression should be compared with these at comparable effective complexity. Do not fit an independently reconstructed w(z) curve as if its points had independent errors.

For gravity, use a stable, finite-dimensional subset of an established EFT solver. Include both expansion and perturbation predictions. General changes to H(z), μ(k,z) or gravitational slip are useful phenomenological diagnostics but do not automatically correspond to a stable covariant theory. Canonical-quintessence constraints do not apply identically to an effective w reconstructed in modified gravity.

A cheap analytic stage can guide the search. Locally write

\[
d\simeq f_\Lambda(\theta_0)+J\,\delta\theta+B\alpha+N\eta+\epsilon,
\]

where J contains baseline-parameter responses, B contains candidate-physics responses and N contains nuisance responses. Whiten with the data covariance. Project candidate responses away from the directions spanned by the baseline and allowed nuisance shifts, or use the corresponding marginalized Fisher matrix with informative priors. Singular values then identify directions that the observations can distinguish. This is a local sensitivity calculation, not a universal proof of identifiability or a novel method by itself.

After this screening, run exact nonlinear inference for a small number of candidates. For example, fit geometry data and ask whether a candidate predicts held-out growth/lensing measurements, with their own nuisance parameters properly marginalized. Evaluate the inverse direction as a robustness check only if it was prespecified. Correlated probes need a joint covariance or a conditional predictive test, not naive multiplication of likelihoods.

A concrete first campaign could compare three physical families and three calibration/systematic alternatives, with fixed priors and an initial budget of 10,000 exact theory evaluations. Promote at most three candidates to full inference. Profile likelihood can diagnose prior-volume effects; evidence and predictive scores quantify other aspects. Every promoted point must be re-evaluated with the exact solver at higher numerical accuracy.

The potential paper would be an identifiability and false-discovery study: where minimal changes become distinguishable, which data separate them, and how often an adaptive search invents a departure in a ΛCDM mock universe. This is a plausible contribution, not a confirmed gap in the literature. The 2026 KGB work is a particularly close comparator. [12]

**Existing ML results constrain the novelty claim.** Desmond, Bartlett and Ferreira already applied exhaustive symbolic regression to the SPARC radial-acceleration relation. Their mock tests show that the data’s uncertainties and dynamic range limit recovery even of the generating law. A shorter formula alone would therefore be weak evidence for a new force law. [17]

Symbolic regression has also been applied directly to DESI DR2 and supernova dark-energy reconstructions. Its outputs depend on the search grammar, data treatment and external assumptions; the term “model independent” does not remove those choices. [18]

CosmoPower and CosmoPower-JAX already accelerate power-spectrum prediction and differentiable inference. Use these as baselines for inference acceleration, and never use a ΛCDM-trained surrogate to validate an untrained gravity extension. Any extended emulator needs exact training data in its new domain and measured error in likelihood space. [19,20]

SimBIG demonstrates simulation-based inference from nonlinear galaxy clustering, including information beyond standard two-point summaries. This is a mature direction with substantial forward-model validation. CAMELS and Quijote offer existing simulations for further work; they do not automatically provide training examples for arbitrary new gravity. [21,22,23]

**An agentic loop should produce auditable comparisons.** The following division of responsibilities is a proposed future workflow; no autonomous research campaign or remote hardware job was launched during this review.

| Component | Bounded responsibility | Required output |
|---|---|---|
| Literature and proposal agent | Propose a physical family, identify nearest prior work, specify a discriminating prediction | Versioned hypothesis with assumptions, units and prior |
| Mathematics checker | Derive limits, transformations, regularity and stability conditions | Reproducible algebra or proof plus counterexamples attempted |
| Experiment builder | Implement the candidate against a frozen likelihood contract | Exact executable, input hashes, baseline-reproduction result |
| Scheduler | Allocate calls to informative and sufficiently validated candidates | Compute ledger and a stop/promote decision |
| Inference engine | Fit fixed models with standard samplers or validated surrogates | Chains, convergence diagnostics, predictive intervals and evidence uncertainty |
| Critic/reproduction step | Recompute claims independently and inspect nuisance dependence | Confirmed result, narrowed claim or rejection |
| Locked evaluator | Run prespecified held-out and mock tests | Selection-calibrated performance unavailable to the proposing stage |

Start with a deterministic scheduler and ordinary Bayesian optimization or active learning. A learned scheduler is justified only if it beats those under the same total compute and correctness requirements. An LLM’s confidence is not an acquisition function or a scientific score.

Preserve every explored branch. Adaptive model search creates a multiple-testing problem even if each individual candidate has an attractive chi-squared. A Bayes factor for the final selected formula does not automatically account for all the formulae tried unless the whole search distribution is represented in the prior/model space.

Run the complete adaptive pipeline on synthetic null data, including selection and promotion. An initial 200 null realizations can estimate the frequency of fairly common false positives and catch obvious failure. It cannot establish a 5σ tail. For a rare-event claim, design a much larger or otherwise statistically justified calibration. A permanent holdout reduces selection bias but does not replace systematics modelling.

**Lean is useful for a narrow mathematical contract.** If “Lean” means the theorem prover, use it for reusable identities such as the transformed field equation, branch/sign statements with explicit assumptions, dimensional encodings, and selected conservation arguments. Lean’s kernel checks formal proofs; it does not check that the physical axioms are true or that a numerical likelihood is an adequate description of data. [24]

For this project, formalization should follow a successful analytic and numerical prototype. A short exact result is a good target. Formalizing a complete Boltzmann hierarchy, inference code or coupled hydrodynamic simulator is a separate large project. No Lean proof was executed here. If “lean” instead means a compact workflow, the same priorities apply: few families, exact baselines, frozen scores, and aggressive stopping of uninformative branches.

**The hardware is sufficient for the first two tracks.** These are planning judgments, not measurements on the user’s workstation.

| Stage | Sensible use of hardware | First budget |
|---|---|---|
| Equation audit and spherical TG | Laptop/CPU; avoid GPU setup overhead | Hours to a few focused days including interpretation |
| Axisymmetric galaxy benchmark | CPU sparse solvers, parallel independent galaxies/chains; GPU only after profiling | A two-week pilot; an article-quality study likely takes additional weeks |
| Background BAO/SN screening | CPU or batched JAX calculations | A few days to reproduce released baselines and nuisance handling |
| Full linear cosmology | CPU CLASS/CAMB/EFT solver evaluations; 5090 for validated emulator training and batches | Pilot 10⁴ exact calls, then size the campaign from measured cost |
| Simulation-based robustness | Existing CAMELS/Quijote summaries before snapshots; 5090 for networks | Select one observable and one known simulation domain |
| New cosmological hydrodynamic suite | Much larger engineering and validation effort | Defer until a cheap discriminating test identifies the need |

For planning, T_wall ≈ N_calls × t_call / N_effective_workers, plus inference, I/O and orchestration overhead. If a measured exact call costs 5 seconds, 10⁴ calls require about 14 core-hours; this is an illustration, not a benchmark. Some models, likelihoods and accuracy settings cost much more. Chain mixing and stability-boundary rejection can dominate.

The RTX 5090 is ample for small emulators and batched inference. Use CPU reference calculations at appropriate precision. If the 1,300-core allocation is also available to this project, independent exact solves and mock campaigns are a natural use; do not assume a single solve scales to all cores. Access to the additional Blackwell machines would mainly help larger surrogate/SBI campaigns and repeated training. High utilization is useful only after the scientific workload is well specified.

**The public-data starting set is already strong.** Use released summary statistics and likelihoods initially. Their availability differs from that of raw catalogs and from that of newly announced analyses.

| Resource | Immediate role | Access and cautions |
|---|---|---|
| [SPARC](https://astroweb.case.edu/SPARC/) | Galaxy rotation curves, baryonic mass models, photometry and sample metadata | Public. Source references and geometry conventions need to travel with the data. [25] |
| [DESI DR2 products](https://www.desi.lbl.gov/2025/10/06/desi-dr2-cosmology-chains-and-data-products-released/) and [Cobaya BAO vectors](https://github.com/CobayaSampler/bao_data) | Reproduce BAO inference and published chains | Use official covariance. Pin the release/commit; account for overlap if adding full-shape or updated Lyα constraints. [26,27] |
| [DES-Dovekie](https://github.com/des-science/DES-SN5YR) | Recalibrated SN likelihood, statistical/systematic covariance and mocks | Current repository supersedes the original DES-SN5YR analysis; an older tag remains available. Do not combine these as independent samples. [28] |
| [Pantheon+ and SH0ES](https://github.com/PantheonPlusSH0ES/DataRelease) | A separate SN robustness branch or correctly calibrated ladder analysis | Different releases share objects/calibrators. Account for overlap; uncalibrated SN distances alone do not measure H0. [29] |
| [NASA LAMBDA ACT products](https://lambda.gsfc.nasa.gov/product/act/actadv_prod_table.html) | CMB spectra, lensing and released likelihood products | Respect supplied combinations, foreground treatment and overlap. [30] |
| [KiDS science products](https://kids.strw.leidenuniv.nl/sciencedata.php) and [DES Y6](https://www.darkenergysurvey.org/des-y6-cosmology-results-papers/) | Growth and lensing tests | Use survey likelihoods, scale cuts and nuisance parameters, rather than treating an S8 summary as a generic likelihood. [7,31] |
| [CAMELS](https://camels.readthedocs.io/en/latest/) | Baryonic-model robustness and simulation-based inference | Public release products and documented access procedures; start with summaries. Large raw archives are unnecessary for the pilot. [22] |
| [Quijote](https://quijote-simulations.readthedocs.io/en/latest/) | Cosmology response, covariance and learning benchmarks | Public simulations and summaries within their simulated model families. [23] |
| Unite | Latest SN result for literature context | The checked September preprint says data/likelihood release is deferred until acceptance. [6] |

The cosmology stack can be CLASS/CAMB plus Cobaya, adding mochi_class or EFTCAMB for a specific supported extension. NumPy/SciPy is sufficient for the first TG reference solver; JAX is useful after a differentiable implementation has an independent numerical reference. [10,11,32]

**My ranking is based on tractability and potential information gain.** It is not a probability of discovering new physics.

| Rank | Project | Why pursue it | Principal risk |
|---:|---|---|---|
| 1 | Predictive, geometry-aware TG/SPARC benchmark with a shared K law | Directly uses the supplied physical idea; exact simplification; inexpensive checks; interpretable outcomes | The theory may require unconstrained per-galaxy freedom or fail after fair boundary/geometry treatment |
| 2 | Minimal cosmological modifications versus nuisance directions | Directly addresses ΛCDM tensions; public likelihoods; strong role for analytic response and calibrated active learning | Crowded prior art; apparent gains may reflect priors or calibration |
| 3 | Stable scalar–tensor active-learning benchmark | Physically constrained search, with clear existing software and June 2026 comparators | Requires greater cosmology expertise; “automating the scan” is not sufficient novelty |
| 4 | SIDM surrogate with baryonic and core-collapse uncertainty | Connects computational physics, simulation and galaxy data | Serious astrophysical modelling, with recent semi-analytic competitors |
| 5 | Unrestricted symbolic search for a replacement cosmology | Broad exploratory scope | Enormous search bias, poorly defined complexity, and insufficient independent validation |

For a first two-week pilot, I would choose rank 1. Days 1–3: resolve equation conventions and physical boundaries, reproduce the published conditional setup, freeze a quality-selected sample protocol. Days 4–7: construct an axisymmetric reference and matched nuisance treatment on a small prespecified set. Days 8–10: compare constant, independently fitted and simple shared K models. Days 11–14: run whole-galaxy predictive tests and mock checks. If a stage fails, the outcome is a documented limitation and a decision about whether the failure is numerical, observational or theoretical.

A useful project can succeed by finding a simple law that predicts new data, proving that a claimed distinction is not identifiable from current data, or demonstrating that an apparent fit advantage comes from boundaries or nuisance flexibility. Those outcomes make the agentic loop scientifically valuable even if ΛCDM remains preferred.

**Sources and reading order.** The supplied paper was read in full as local extracted text, with its equation pages also visually inspected. Other entries are primary papers or official project/software sources; many were assessed through abstracts or release documentation rather than full-paper reproduction. Source links are enough to reconstruct the literature path, but do not imply every cited analysis was rerun.

[P] Pszota & Ván, *Field equation of thermodynamic gravity and galactic rotational curves*, Physics of the Dark Universe 46 (2024), 101660. Supplied PDF. [DOI](https://doi.org/10.1016/j.dark.2024.101660); [preprint](https://arxiv.org/abs/2306.01825).

[1] Louis et al., [ACT DR6 power spectra, likelihoods and ΛCDM parameters](https://arxiv.org/abs/2503.14452).

[2] H0DN Collaboration, [The Local Distance Network](https://arxiv.org/abs/2510.23823), A&A 708, A166 (2026).

[3] DESI Collaboration, [DR2 results II: BAO and cosmological constraints](https://arxiv.org/abs/2503.14738).

[4] Popovic et al., [DES-Dovekie reanalysis](https://arxiv.org/abs/2511.07517), revised March 2026.

[5] DESI Collaboration, [DR2 results IV: Lyα Alcock–Paczyński measurements](https://arxiv.org/abs/2607.27410), July/August 2026. [Official explanation](https://www.desi.lbl.gov/2026/07/30/new-desi-dr2-lyman-alpha-results-shed-light-on-dark-energy/).

[6] Camilleri et al., [Supernovae Unite](https://arxiv.org/abs/2609.05053), September 2026 preprint.

[7] KiDS, [official science-data and paper index](https://kids.strw.leidenuniv.nl/sciencedata.php).

[8] DES Collaboration, [Y6 galaxy clustering and weak-lensing constraints](https://arxiv.org/abs/2601.14559).

[9] DESI Collaboration, [Extended dark-energy analysis using DR2 BAO](https://arxiv.org/abs/2503.14743).

[10] Cataneo & Bellini, [mochi_class](https://arxiv.org/abs/2407.11968); [source code](https://github.com/mcataneo/mochi_class_public).

[11] [hi_class](https://hiclass-code.net/) and [EFTCAMB](https://github.com/EFTCAMB/EFTCAMB), official projects.

[12] Cataneo & Koyama, [Nonparametric exploration of minimally coupled gravity with phantom crossing](https://doi.org/10.1103/ckns-5ytv), Phys. Rev. D 113, 122006 (June 2026).

[13] Poulin et al., [Impact of ACT DR6 and DESI DR2 for Early Dark Energy and the Hubble tension](https://arxiv.org/abs/2505.08051).

[14] [Astrophysical tests of dark matter self-interactions](https://doi.org/10.1103/m2vm-59y3), Rev. Mod. Phys. 97, 045004 (2025).

[15] [An enhanced isothermal Jeans approach to constraining dark matter self-interactions from galactic kinematics](https://academic.oup.com/mnras/article/549/3/stag969/8691061), MNRAS (2026).

[16] Skordis & Złośnik, [A new relativistic theory for MOND](https://arxiv.org/abs/2007.00082), Phys. Rev. Lett. 127, 161302 (2021).

[17] Desmond, Bartlett & Ferreira, [On the functional form of the radial acceleration relation](https://arxiv.org/abs/2301.04368); [ESR source code](https://github.com/DeaglanBartlett/ESR).

[18] Sousa-Neto et al., [DESI DR2/SN symbolic regression](https://arxiv.org/abs/2502.10506), revised June 2025; [journal publication](https://doi.org/10.1016/j.dark.2025.102108).

[19] Piras & Spurio Mancini, [CosmoPower-JAX](https://arxiv.org/abs/2305.06347); [source code](https://github.com/dpiras/cosmopower-jax).

[20] [CosmoPower](https://github.com/alessiospuriomancini/cosmopower), official implementation.

[21] [SimBIG nonlinear/non-Gaussian cosmological inference](https://www.nature.com/articles/s41550-024-02344-2); [project documentation](https://changhoonhahn.github.io/simbig/).

[22] [CAMELS documentation](https://camels.readthedocs.io/en/latest/) and [data-access instructions](https://camels.readthedocs.io/en/latest/data_access.html).

[23] [Quijote documentation](https://quijote-simulations.readthedocs.io/en/latest/).

[24] [Lean reference](https://lean-lang.org/doc/reference/latest) and [propositions and proofs](https://docs.lean-lang.org/theorem_proving_in_lean4/Propositions-and-Proofs/).

[25] [SPARC official data page](https://astroweb.case.edu/SPARC/).

[26] [DESI DR2 chains and products release](https://www.desi.lbl.gov/2025/10/06/desi-dr2-cosmology-chains-and-data-products-released/).

[27] [Cobaya BAO data](https://github.com/CobayaSampler/bao_data).

[28] [DES-Dovekie release repository](https://github.com/des-science/DES-SN5YR).

[29] [Pantheon+/SH0ES release repository](https://github.com/PantheonPlusSH0ES/DataRelease).

[30] [ACT data on NASA LAMBDA](https://lambda.gsfc.nasa.gov/product/act/actadv_prod_table.html).

[31] [DES Y6 cosmology release](https://www.darkenergysurvey.org/des-y6-cosmology-results-papers/).

[32] [Cobaya theory and likelihood documentation](https://cobaya.readthedocs.io/en/latest/cosmo_theories_likes.html).

**Executable appendix.** Save the code below as `audit.py`, install NumPy and SciPy if necessary, and run `python audit.py`. All tests use synthetic dimensionless inputs. The printed JSON is the result of checks, not fitted observational parameters. The comparison uses the literal sign conventions of the supplied PDF.

```python
"""Equation-level audit of Pszota & Van (2024), not an observational fit.
Dependencies: numpy, scipy. Run: python audit.py
Conventions: S(r) = 4*pi*G*rho(r), g = -d(phi)/dr, constant K > 0.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp


def g_coth(r, K=1.0, S=1.0):
    q = np.sqrt(K*S)
    return 1/(K*r) - q/(K*np.tanh(q*r))


def g_cot(r, K=1.0, S=1.0):
    q = np.sqrt(K*S)
    return -1/(K*r) + q/(K*np.tan(q*r))


def derivative_5point(f, r, h=2e-4):
    return (f(r-2*h)-8*f(r-h)+8*f(r+h)-f(r+2*h))/(12*h)


def equation_residual(f, r, K=1.0, S=1.0):
    g = f(r)
    return derivative_5point(f,r) + 2*g/r + S + K*g*g


def run_audit():
    grid = np.linspace(0.2,1.2,201)
    paper_residual = equation_residual(g_coth,grid)
    cot_residual = equation_residual(g_cot,grid)
    algebra_residual = paper_residual - 2*g_coth(grid)**2

    # Independent integrations of nonlinear g and transformed linear u.
    K = 0.4
    r0, rmax = 1e-4, 5.0
    source = lambda r: np.exp(-r*r)
    g0 = -r0/3 + (1/5-K/45)*r0**3
    u0 = 1-K*r0*r0/6+(K/20+K*K/120)*r0**4
    p0 = -K*r0/3+(K/5+K*K/30)*r0**3
    nonlinear = solve_ivp(
        lambda r,y: [-2*y[0]/r-source(r)-K*y[0]*y[0]],
        (r0,rmax),[g0],rtol=2e-11,atol=2e-13,dense_output=True)
    linear = solve_ivp(
        lambda r,y: [y[1],-2*y[1]/r-K*source(r)*y[0]],
        (r0,rmax),[u0,p0],rtol=2e-11,atol=2e-13,dense_output=True)
    if not (nonlinear.success and linear.success):
        raise RuntimeError('Integration failed')
    radii = np.linspace(0.01,rmax,1000)
    u,p = linear.sol(radii)
    direct = nonlinear.sol(radii)[0]
    transformed = p/(K*u)

    # Isolated constant-density sphere, qR=x, u -> 1 at infinity.
    # Interior u = sinc(qr)/cos(x); requires 0<x<pi/2 for
    # the everywhere-positive, regular branch connected to Newtonian gravity.
    # Exterior mass amplification is analytic, not a galaxy fit.
    xs = np.array([0.2,0.5,1.0,1.3,1.5])
    amp = 3*(np.tan(xs)/xs-1)/xs**2
    result = {
        'scope':'Equation checks and synthetic density only; no SPARC/cosmology fit; no Lean proof',
        'uniform_density_test': {
            'K':1.0,'S':1.0,'radius_range':[0.2,1.2],
            'paper_Eq23_max_absolute_residual':float(np.max(np.abs(paper_residual))),
            'corrected_cot_max_absolute_residual':float(np.max(np.abs(cot_residual))),
            'max_error_in_identity_R_paper_equals_2K_g_squared':float(np.max(np.abs(algebra_residual)))
        },
        'gaussian_density_crosscheck': {
            'K':K,'S':'exp(-r^2)','radius_range':[r0,rmax],
            'minimum_sampled_u':float(np.min(u)),
            'max_absolute_g_difference':float(np.max(np.abs(direct-transformed))),
            'max_relative_g_difference':float(np.max(np.abs((direct-transformed)/direct)))
        },
        'isolated_uniform_sphere': {
            'assumptions':'constant K>0, S>0; density truncates at R; regular centre; phi(infinity)=0; Newtonian-connected positive branch',
            'x_critical':float(np.pi/2),
            'x_qR':xs.tolist(),'M_app_over_M_b':amp.tolist()
        },
        'paper_K_flat_regime_speed_km_s':float((3.4e-5)**(-0.5))
    }
    assert np.max(np.abs(cot_residual)) < 1e-7
    assert np.max(np.abs(algebra_residual)) < 1e-7
    assert np.max(np.abs(direct-transformed)) < 1e-8
    return result


if __name__ == '__main__':
    print(json.dumps(run_audit(),indent=2))
```

**Execution record.** This audit ran locally during this review.

```json
{
  "scope": "Equation checks and synthetic density only; no SPARC/cosmology fit; no Lean proof",
  "uniform_density_test": {
    "K": 1.0,
    "S": 1.0,
    "radius_range": [
      0.2,
      1.2
    ],
    "paper_Eq23_max_absolute_residual": 0.2682110481022828,
    "corrected_cot_max_absolute_residual": 4.247784415878364e-12,
    "max_error_in_identity_R_paper_equals_2K_g_squared": 4.279392812334137e-12
  },
  "gaussian_density_crosscheck": {
    "K": 0.4,
    "S": "exp(-r^2)",
    "radius_range": [
      0.0001,
      5.0
    ],
    "minimum_sampled_u": 0.8418442408438386,
    "max_absolute_g_difference": 4.800362884971321e-11,
    "max_relative_g_difference": 5.010181446444847e-10
  },
  "isolated_uniform_sphere": {
    "assumptions": "constant K>0, S>0; density truncates at R; regular centre; phi(infinity)=0; Newtonian-connected positive branch",
    "x_critical": 1.5707963267948966,
    "x_qR": [
      0.2,
      0.5,
      1.0,
      1.3,
      1.5
    ],
    "M_app_over_M_b": [
      1.0162633157521803,
      1.1112597562509716,
      1.6722231739647069,
      3.143517225263512,
      11.20126217526375
    ]
  },
  "paper_K_flat_regime_speed_km_s": 171.49858514250883
}
```

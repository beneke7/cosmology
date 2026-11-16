# Next experiment: physical thermodynamic-gravity pilot

## September 25 user steering: prioritize physical tests

The user's current priority is thermodynamic gravity, physical detections, and
matched alternative theories. The branch contains two distinct models and
must not blend their parameters or evidence: (1) the 2024 nonlinear constant-K
galaxy model, tested conditionally on NGC 3198; and (2) the August 2026
linear higher-gradient Newtonian model with short-range Yukawa scales.
The PDFs are local at [`2024 galaxy paper`](context/user_provided/pszota_van_2024.pdf)
and [`2026 higher-gradient paper`](context/papers/tg_higher_gradient_2026.pdf).

Current physical evidence includes the audited **K=0** source/force baseline,
the nonzero-K axisymmetric stationary solve, its independent PDE/eigenmode
audit, and same-resolved-source NFW/MOND controls. The positive branch closes
at a fine-grid `Kcrit=5.3714e-5 (km/s)^-2` for `Υ*=0.762`; the principal-mode
limit has `χ²diag=3585.9`, 17.2 km/s RMS on the 43 SPARC points. This limit is
not a finite-normalization solution. The nearest tested finite solution at
`0.998 Kcrit` has `u_max≈987` and `χ²diag=3641.8`. These remain descriptive
diagonal scores, not covariant model selection.

The best matched-source controls are substantially closer to the curve:
Newtonian+NFW profiles `Υ*, M200, c200` and gets `χ²diag=58.5` (7.68 km/s RMS);
isolated simple-MOND with fixed `a0` and profiled `Υ*` gets `1032.2` (10.59
km/s RMS). TG's branch-edge envelope at `Υ*=0.762` gets `3585.9`, but is a
singular limit, not a realization. Signed gas acceleration is preserved; all
43 radii/velocity errors were matched exactly. There is no released radial
covariance, and the resolved profile is a parametric axisymmetrization rather
than a unique 3D source. Distance/inclination variation and vertical-profile
alternatives remain open, so this is a conditional negative for this NGC 3198
source/boundary family—not a TG-wide falsification or a physical detection.

Source facts remain important: THINGS and HALOGAS main-disk H I masses agree to
about 1.7%; S4G P5 gives `1.6318e10 Msun` inside its official analysis ellipse
at M/L=0.6. Only complete stellar rings through 13.75 kpc are measured; the
outer source is explicitly extrapolated. The independent Freeman-disk and
spectral checks pass at approximately `2e-4`; the independent fine-grid solver
agrees on the selected curve to `0.00067 km/s` RMS and on Kcrit to `1.5e-5`
fractionally. THINGS-vs-HALOGAS gas and ±25% stellar-tail changes shift the
branch-edge curve by only about 0.1–0.3 km/s RMS, while the stated stellar M/L
range changes it materially. The map-derived vs SPARC Rotmod stellar-force RMS
discrepancy of 18.39 km/s remains a source-definition caveat. See
[`source audit`](work/tg_axisymmetric_pilot/SOURCE_AUDIT.md),
[`branch-edge summary`](work/tg_axisymmetric_pilot/axisymmetric_branch_edge_summary.json),
[`independent principal-mode audit`](work/tg_axisymmetric_independent_review/PRINCIPAL_MODE_REPORT.md),
[`matched alternatives`](work/tg_alternatives_contract/map_matched/REPORT.md),
[`K=0 result`](work/tg_axisymmetric_pilot/newtonian_k0_baseline.json), and
[`independent S4G review`](work/tg_s4g_independent/REPORT.md).

**Next: propagate shared distance/inclination and vertical-source uncertainty,
then prioritize a whole-galaxy holdout rather than refitting NGC 3198.** A
predictive TG relation needs a predeclared shared-K/scaling law and a clean
external galaxy test. Separately, the 2026 higher-gradient paper's physical
short-range prediction still needs a model-specific two-range fit to published
torsion data; no direct detection is present in the reviewed literature.
NGC 3198 alone is not a new-field detection or a relativistic/cosmological
completion.

The 2026 short-range paper is separately under audit. The 2026 review
`arXiv:2605.18212v2` has been acquired locally; its plotted Yukawa limits
restrict to attractive alpha>0 and do not directly test the TG model's
alpha=-1 one-range or correlated two-range sector. The UW/HUST primary curves
remain the current relevant published results located, with no new
model-specific two-range contour. Preserve the sign-specific UW 25–28 μm
published-contour crossing as a read-off, not a new calibrated limit. See
[`short-range literature update`](work/tg_short_range_literature/REPORT.md)
and [`lab audit`](work/tg_lab_bound_independent/REPORT.md).

The CMB/IVS lane remains useful but is lower priority for this user-directed
branch; prior results and safeguards are retained below for continuity.

## September 25 BBN-calibration update

The BBN source audit, profile and forward feasibility pass are now complete. PDG 2025's broad SBBN mapping (ω_bh²=0.02205±0.00043) is suitable as an independent *probe-level* calibration beside DESI only if the model leaves the BBN thermal history standard. With CAMB 2.0.4 and the full released 13×13 DESI DR2 covariance, a 25-node ω_b profile changes the data-only χ² by <9.1×10⁻⁴ across [0.0205,0.024]; 600 deterministic starts, five differential-evolution searches and a standalone point-score audit find no lower tested basin. Thus this combination removes the pathological high-ω_b BAO-only edge only conditionally; it is not a posterior or H₀ determination. See [compute31](work/compute31/REPORT.md), [compute33](work/compute33/REPORT.md), and [point audit](work/compute32/independent_bbn_bao_score_audit.json).

The generic public PRIMAT v0.3.2 standard-background grid and five-point repeat are bit-identical, and the PDG-center D/H prediction is descriptively within 0.16 observation-only σ. No abundance/rate likelihood or posterior was run, and the abundance data behind PDG's ω_b mapping cannot be counted again independently. For the screened late-time IVS point, a separate radiation-era comparison gives at most 7.27 ppm in H over 10⁸≤z≤10¹⁰, supporting the standard-BBN prior only at that point on expansion grounds—not across the IVS posterior or at abundance level. Details and exact hashes are in [the source/prior audit](work/bbn_prior_audit/REPORT.md), [PRIMAT report](work/data9_primat/REPORT.md), and [IVS–BBN background comparison](work/theory13/ivs_bbn_expansion_compare.md).

## Current evidence

The frozen-BAO-shape prediction has passed independent covariance-domain reproduction. Four overlapping 455-SN conditional folds give integrated-intercept Δ(−2 log predictive) = −4.0366 (IVS−ΛCDM), with six locked edge variants all negative (−4.1180 to −3.0486). The lowest-redshift fold contributes −5.5418; the other three collectively favor ΛCDM by +1.5052. This is descriptive same-sample prediction, not a joint BAO×SN likelihood or evidence. See [result](experiments/dovekie_frozen_shape_cv/result.json) and [independent review](work/compute8/dovekie_result_independent_review.md).

The Dovekie-only full-covariance IVS profile is now complete and independently reproduced. Its finite-box minimum is χ²=1629.999216924 at (Ωm0,g)=(0.39301742,−0.62938781), only Δχ²=1.421319 below the Dovekie ΛCDM profile despite one additional shape parameter (g) beyond the shared (Omega_{m0}). The frozen BAO-selected IVS point is +0.598982 above the SN-only minimum. No boundary was touched. Independent RK45/full-covariance GLS scores agree within 4.55×10⁻¹², and ODE tolerance refinements, the g=0 limit, adaptive-distance quadrature, and the sampled low-z positive-density/split checks pass. The result remains an exploratory finite-box background profile, not significance or physical-model validation. See [the experiment](experiments/dovekie_ivs_profile/result.json), [plot](experiments/dovekie_ivs_profile/profile_surface.png), and [independent scorer/review](work/compute10/dovekie_ivs_profile_independent_review.md).

The declared wider-domain check is complete: a 61×61 grid on Ωm0=[0.005,0.95], g=[−5,5] changes the minimum by only 5.68×10⁻¹¹ and finds no new boundary contact. This reuses the production scorer and is a finite-domain sensitivity check, not an independent proof of the global minimum. See [record](experiments/dovekie_ivs_profile/domain_sensitivity_record.json).

The exact finite-search null calibration is complete and independently audited. Of 200 full-covariance Gaussian SN mocks generated under the fitted Dovekie ΛCDM shape and refit with the same search, 50 reach or exceed the observed Δχ²=1.421318647: empirical finite-search exceedance fraction 25%, exact Clopper–Pearson 95% interval [19.16%,31.60%]. There were no failed realizations. The independent sampler and 4,001-point accepted-endpoint physicality audit pass. This is a coarse fitted-null calibration, not a universal p-value, model evidence, or full physical-model result. See [experiment record](experiments/dovekie_ivs_profile/null_calibration_record.json), [mock report](work/compute11/REPORT.md), [independent audit](work/compute12/dovekie_ivs_null_mock_audit.md), and [plot](experiments/dovekie_ivs_profile/null_calibration.png).

The standalone local identifiability audit is now complete. At the pinned profile point, the intercept-projected whitened Ωm and g response cosine is +0.969389, and the two-column response condition number is 40.94; only 6.0% of the post-intercept g information remains after projecting Ωm. The residual cosine with the remaining g direction is only 1.58×10⁻⁷, which is expected from first-order stationarity at the joint profile minimum and is a fit-consistency diagnostic, not independent evidence against IVS. The forward-sensitivity ODE agrees with a separate central-difference ODE sweep; direct GLS and projected/Schur algebra agree to better than 1e-11. These are local response-geometry diagnostics, not posterior errors or significance. In combination with the matched null, they make flexible SN-only background extensions low priority. See [independent report](work/compute13/dovekie_profile_identifiability.md), [code/results](work/compute13/dovekie_profile_identifiability.py), and [durable record](experiments/dovekie_ivs_profile/identifiability_record.json).

Counter-signals remain material: the IVS BAO gain is exceeded by 256/1200 matched flat-null searches, and its BAO redshift-block conditional score is worse by +22.7503, dominated by the z=2.33 pair. Frozen-shape SN gains concentrate at low z. Recent primary-source work likewise warns that the coupling preference depends on the chosen DE dynamics and SN sample; it is not safe to infer an interaction from a restricted background fit alone.

## September 25 physical-solver and baseline update

The radiation-inclusive homogeneous IVS benchmark (compute15/17) and the fixed-metric two-component perturbation equations (compute18/20) both pass their stated numerical checks, but neither is a coupled Einstein–Boltzmann solve. A benchmark-specific early-CDM event audit brackets (g\simeq-1.022557) for a zero by (z=10^{10}); it uses fixed late-time density parameters and massless neutrinos, and does not imply posterior support/exclusion under the paper's physical-density/100\(\theta_{MC}\) prior.

A source-level IDECAMB/CAMB contract is now complete: `H0rde` matches the target background with \(\beta=\Gamma/H_0\), but IDECAMB's generic interacting PPF perturbations are not the target's \(\delta\rho_v=0\), Eqs. 19–20. The exact author patch, full Cobaya configuration, and initial conditions were not located in bounded public-source checks. Any code we build is an independent implementation, not a paper reproduction. See [adaptation contract](work/theory11/ivs_idecamb_patch_contract.md), [solver provenance](work/literature7/solver_provenance.md), and [CAMB source gate](work/compute23/camb_2_0_4_ivs_feasibility.md).

The independent flat-ΛCDM × DESI DR2 BAO profile finds a constrained boundary minimum at \((H_0,\omega_bh^2,\omega_ch^2)\simeq(133.61,0.1,0.43027)\) with \(\chi^2=10.2738863\). Three direct-H0 seeds agree, but the solution pins the broad \(\omega_bh^2\) upper bound; nuisance-reoptimized profiles continue lower outside it. Compute30 independently reproduced the score and found the best of nine reoptimized in-box slices at the boundary; the nearest inward slice is higher by \(4.05\times10^{-9}\). A same-edge polish changes the score by only \(5.05\times10^{-11}\), well below the declared optimizer tolerance. Thus the tested profile supports a local constrained boundary minimum, not global optimality. This is an identifiability/domain diagnostic, not a physical baseline or H0 constraint: it lacks CMB/BBN likelihood information. See [compute29 profile](work/compute29/boundary_audit_report.md) and [compute30 independent review](work/compute30/REPORT.md).

Planck 2018 full likelihood is publicly mapped to Plik HM TTTEEE + Commander low-TT + SimAll lowE, with no lensing component. The exact archive size is undisclosed; an official range probe transferred no bytes. No data or packages have been staged. The BRIEF.md permits expanding the selected research cache to 20 GiB after recording the rationale; this exact public archive is necessary for the next calibrated baseline and to resolve the BAO-only physical-density boundary. I am proceeding with a hard total-cache cap of 20 GiB, staging only the official baseline archive and later selecting only the named likelihood products; preserve CC BY-NC attribution and full hashes. Standard CAMB + Planck ΛCDM validation remains the highest-information next baseline. The data/config contract is [here](work/data7_planck_gate/planck_2018_full_likelihood_gate.md).

## Recommended next test

Do not add more freedom to the SN-only IVS background fit: its local interaction response is mostly an \(\Omega_m\) displacement and matched null mocks routinely achieve the observed gain. Do not spend compute relaxing the DESI-only CAMB density bounds without an early-universe prior: the profile demonstrably moves outside the declared box. The most informative next baseline is standard CAMB + the exact public Planck PR3 Plik TTTEEE, Commander TT, and SimAll EE likelihoods with declared cosmological and nuisance settings. If the pending user decision allows staging the full archive, build a hash-locked project-local environment and validate a ΛCDM run first. Do not add lensing, substitute Plik-lite/chains/distance priors, or exceed 20 GiB total campaign-cache expansion. If the archive is not authorized or cannot be staged reproducibly, record the limitation and continue with source-only work.

For IVS, exact author-run reproduction is blocked unless the user supplies the modified-CAMB patch/base commit, full Cobaya YAMLs/nuisance setup, and initial-condition prescription. Do not repeat the same public search. Once a standard CMB baseline exists, a separate independent IVS pilot needs a recursive pinned CAMB checkout with the `forutils` gitlink resolved, project-local Fortran compiler, joint coupled background/perturbation implementation, and regression checks including \(\Gamma=0\), both signs, radiation shear, conservation, early-time IC and convergence. Never call it source-faithful. Keep BAO and SN separate without justified cross-covariance. The RTX 5090 is not indicated for present scalar background gates; benchmark any full likelihood before scheduling GPU work.

## Physical and data limits

- The BAO/SN IVS fit is still a flat, late-time background analysis with uncalibrated physical-density split and no full perturbation/CMB likelihood. Its observed-range positive split is existential; it is not a full early-time or microphysical viability certificate.
- Separate radiation-background and prescribed-metric perturbation gates now pass, but they are not Einstein–Boltzmann feedback, a perturbative stability test, CMB spectra, or reproduction of the target's massive-neutrino setup. The exact author patch/config/IC remain unavailable in checked sources; if the user has them, request original artifacts, otherwise label future solver work as independent.
- The fixed benchmark's negative-coupling CDM-zero threshold is not a posterior statement. The BAO-only high-density boundary solution is not a physical ΛCDM best fit or H0 measurement.
- Keep DESI BAO and Dovekie as separate probes; no BAO×SN product or exact independence claim is justified by the missing cross-covariance. Do not add DESI Results IV Lyα to the current BAO vector without an official replacement likelihood and cross-covariance.

The campaign state and durable record are `RUN_STATE.json`, `experiments/ledger.jsonl`, `MORNING_REPORT.md`, and `experiments/dovekie_ivs_profile/`. After the next reviewed checkpoint is committed and pushed, invoke the configured one-shot Astra Max Warden after each event.

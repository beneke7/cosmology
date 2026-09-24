# IVS solver/reproducibility path (advisory audit)

Scope: implementation path for the exact (Q=\Gamma\rho_x,\;w_x=-1) interacting-vacuum scenario in arXiv:2609.05660v1. This was a source/local-paper audit only; no software or data were downloaded, installed, cloned, or run. It is not a revised fit or a claim that the source paper’s numerical results have been reproduced.

## Public artifact status

The article reports a modified CAMB + Cobaya analysis and names Planck 2018 `plikTTTEEE+lowl+lowE`. It fixes \(\sum m_\nu=0.06\,\mathrm{eV}\), \(N_\mathrm{eff}=3.044\), and gives model priors, but the visible article/arXiv metadata do not identify a CAMB release/commit, IVS patch, complete Cobaya YAML/INI, likelihood build/version, initial-condition setting, or chain/Cls/transfer products. Bibliography [157] points to the public [`wgcosmo` repository](https://github.com/williamgiare/wgcosmo), which the paper associates with its Cobaya/MCEvidence interface. Its browsable root/README describe a general cosmology toolbox with a `yamls/` collection; no exact IVS CAMB patch or run config was verifiable here. Conclusion: *public general toolbox confirmed; exact model implementation/config not confirmed*. A separate public [IDECAMB](https://github.com/liaocrane/IDECAMB) project patches CAMB/CosmoMC for interacting dark-energy families, but its README describes coupled-quintessence/coupled-fluid models (including PPF), not a verified implementation of this paper’s `Q=Gamma*rho_x`, geodesic-vacuum closure. It is not provenance for this paper’s code.

The named Planck likelihood data are publicly obtainable: the official [Planck Legacy Archive](https://pla.esac.esa.int/) serves public mission products; the [ESA Planck release page](https://www.cosmos.esa.int/web/planck) records the 2018 likelihood-code release. The primary [Cobaya Planck-likelihood documentation](https://github.com/CobayaSampler/cobaya/blob/master/docs/likelihood_planck.rst) identifies `COM_Likelihood_Data-baseline_R3.00.tar.gz`, the official PLA route, and example Commander and Plik members. Thus the likelihood family is available, but not the authors’ end-to-end config. Software/data versions must still be pinned, acquired, hashed, and smoke-tested. Do not substitute distance priors or a different PR4/NPIPE/CamSpec likelihood when claiming replication.

At audit date (2026-09-24), PyPI lists CAMB **2.0.4** (25 Aug 2026); its verified publishing attestation points to upstream commit `a6de8cc59124c8bbe3924f1c546f96cef73dcabb` (sdist SHA-256 `2f1f5f3b3964a3746693dc88a54218dca4ea049a04726f8177664c11d2d2616e`). The upstream repository separately labels 1.6.6 as its last “Planck-era precision” release and retains `CAMB_v0` as the Fortran-oriented source used in the Planck 2018 analysis. This identifies reproducible upstream choices, **not** the version modified by the paper. First ask for the authors’ exact patch/base version; if unavailable, pin a fresh implementation and do not call it a source-paper reproduction.

## Closure and equation checks

The paper states separate conservation for baryons/radiation/neutrinos; `Q=Gamma*rho_x`, with positive `Gamma` transferring energy out of CDM into vacuum; synchronous-gauge scalar perturbations; no momentum transfer in the DM rest frame; and, for `w_x=-1`, no DE perturbative contribution, giving Eqs. (19)-(20). Under the usual covariant split this is the geodesic-CDM choice, i.e. the transfer four-vector is parallel to CDM four-velocity (the sign depends on which sector’s `Q^mu` is named). In CDM-comoving synchronous gauge, `theta_c=0` for the regular CDM-frame solution; vacuum is homogeneous on those slices. This is a closure assumption, not a universal property of every `w=-1` interaction.

With the paper’s signs, the background equations are

\[
\dot\rho_x=\Gamma\rho_x,\qquad \dot\rho_c+3H\rho_c=-\Gamma\rho_x.
\]

For positive today’s vacuum density and expanding `H>0`, exact integral identities are

\[
\rho_x(a)=\rho_{x0}\exp\!\left[-\Gamma\int_a^1\frac{d\ln a'}{H(a')}\right],\quad
\rho_c(a)=a^{-3}\left[\rho_{c0}+\int_a^1\frac{\Gamma}{H(a')}\rho_x(a')a'^3\,d\ln a'\right].
\]

They provide independent unit/sign oracles; `Gamma` has inverse-time units and the perturbation source is the conformal-time rate `a*Gamma`. For past `a<=1`, positive `Gamma` with positive `rho_x,H,rho_c0` makes the CDM bracket larger than `rho_c0`; negative `Gamma` can drive that bracket through zero. Require the bracket stay strictly positive across the whole integration range, and `H^2>0`. Given `rho_x0>0`, the exponential solution keeps `rho_x>0` at finite times. These are background admissibility tests, not perturbative-stability proofs. In the given density-contrast equation, `+a*Gamma*(rho_x/rho_c)*delta_c` enhances the derivative of positive `delta_c` for positive `Gamma`, and damps it for negative `Gamma` at fixed positive densities; metric feedback means that term’s sign alone is not a stability classification.

Paper Eqs. (19)-(20) are

\[
\delta_c'=-(\theta_c+h'/2)+a\Gamma(\rho_x/\rho_c)\delta_c,\qquad
\theta_c'=-\mathcal H\theta_c.
\]

Choices that still need to be explicit in a reproducer or obtained from the authors: (1) write the covariant transfer sign convention and `Q^mu || u_c^mu`; (2) clarify whether the paper’s phrase “zero anisotropic stress” applies only to dark-sector stresses—standard photon/neutrino shear must be retained or this is not a full CMB Boltzmann model; (3) state scalar initial conditions (adiabatic or specified correlated isocurvature; not named in the visible paper/config); (4) set early radiation, massive-neutrino hierarchy, helium/BBN, recombination and primordial-spectrum settings; (5) define growth/RSD density and velocity convention. CAMB standard transfer outputs explicitly weight baryon and CDM perturbations; do not silently treat a CDM-only perturbation as total matter. (6) verify conversion of sampled `Gamma/H0` to CAMB internal time and wavenumber units. Do not introduce an early-time interaction switch except as a clearly separate robustness model.

## CAMB adaptation map

The verified CAMB 2.0.4 [`Fortran readme`](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/readme.html) maps background and perturbation evolution, `dtauda(a)` and `init_background` to `equations.f90`. Its pinned [`equations.f90`](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/equations.f90) shows `dtauda` obtains dark energy from the separate DE background interface; scalar `derivs` uses fixed `grhoc_t=State%grhoc/a` and current CDM derivative `ayprime(ix_clxc)=-k*z`. This means a complete extension cannot be a late-time `H(z)` wrapper alone.

1. Replace the stock independent CDM/constant-vacuum background by coupled \(\rho_c(a),\rho_x(a),H(a)\) at every scale factor used by thermodynamics, distances, sound horizon and perturbations. Update background initialization and density outputs as well.
2. In scalar synchronous/CDM-frame `derivs`, implement the paper’s CDM source with the physical-unit conversion checked (`Gamma/H0` may be the sampled dimensionless parameter); keep CDM geodesic and do not evolve an independent vacuum density-contrast hierarchy for this closure. Ensure metric constraints and Einstein sources use the coupled densities consistently.
3. Expose `Gamma` in `CAMBparams` and the Python wrapper. CAMB 2.0.4 [parameter docs](https://camb.readthedocs.io/en/latest/model.html) explicitly require a Fortran `CAMBparams` field plus the corresponding Python `model.py` `_fields_` entry. Then test Cobaya parameter transfer end-to-end.
4. Audit derived/background quantities and outputs (`Omega_c`, equality, `zstar`/`zdrag`, `r_drag`, transfers, lensing, total-matter and baryon+CDM growth); current result methods use present CDM normalization and standard `a`-scaling at several points. Add exact `Gamma=0` recovery and an independent background ODE/integral comparison. Check CAMB’s CDM-frame `z`/`etak` variables against the paper’s `h'/2` before mapping the perturbation source.

For model closure, retain photons/neutrinos and their standard perturbation hierarchies, thermal history, and radiation-era initial conditions. Start linear: evaluate TT/TE/EE (and lensing if included in the chosen likelihood), transfers/growth at named redshifts, and recompute `r_drag`. Verify `Gamma=0` recovers stock background/spectra within declared tolerances; independently compare the background and source term; vary integration/accuracy controls and demonstrate convergence. Any non-linear `P(k)`/HALOFIT conclusion needs separate model-validity checks.

## Minimal staged gates and resource limits

| Gate | Pass evidence | Bound / do not do |
|---|---|---|
| 0 — source contract | Request modified-CAMB diff/commit + base version, complete Cobaya config, likelihood version, scalar IC and fiducial output. If unavailable, record limitation. | No download, code edit, or fits. |
| 1 — background | Independent ODE and quadrature agree with integral identities; `Gamma=0` recovers standard densities/distances; both signs scanned from radiation era to today with `rho_c>0`, `rho_x>0`, `H^2>0`. | ≤1 CPU core, ≤10 min, no GPU/data. Stop on first sign/unit/conservation/positivity failure. |
| 2 — two-component perturbation | Eqs. 19-20 in CDM-comoving synchronous gauge pass `Gamma=0` and ±Gamma fixed-metric tests, conformal `a*Gamma` units, dark-sector conservation/gauge diagnostics and independent source implementation. | ≤1 CPU core, ≤10 min; no CAMB sampler or likelihood. Stop on unresolved frame/sign/units. |
| 3 — full linear CAMB CMB + growth | Pinned and tested solver; coupled densities feed thermo/distances/transfers; standard early species/IC/BBN/recombination explicit; `Gamma=0` spectra recovery; convergence for spectra, sound horizon and growth; exact Planck likelihood smoke test before bounded posterior work. | First solver smoke one process, checkpoint within 30 min; benchmark before scaling. No GPU absent end-to-end speed evidence. Begin inference with ≤4 chains, one BLAS/OpenMP thread each; scale only on measured throughput and RAM. |

For Gate 3 acquire only the needed Planck baseline likelihood code/data, not sky maps, all chains, or replacement products. Inspect the official PLA package size first. Keep within the campaign’s 200 MiB bootstrap cap unless root records a specific reason and invokes the documented 20 GiB research-cache expansion; hash the archive and list exact internal members. The precise MCMC/CPU budget should follow a one-point CAMB timing benchmark. If source artifacts remain unavailable, explicitly name results an independent implementation.

## Provenance and primary links

- Local source PDF: `context/papers/interacting_de_desi_dr2_2026.pdf`, SHA-256 `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`; acquisition record is in `context/manifest.lock.json` (`verified_existing`, HTTP 200). Metadata: [arXiv 2609.05660](https://arxiv.org/abs/2609.05660), DOI 10.1103/6kdf-vzq4.
- Exact model predecessor cited as [141]: Yang et al., [arXiv:1808.01669](https://arxiv.org/abs/1808.01669). Covariant/geodesic-CDM closure precedent, **not the same `Q(a)` law**: Martinelli et al., [arXiv:1902.10694](https://arxiv.org/abs/1902.10694).
- CAMB version/provenance: [PyPI 2.0.4](https://pypi.org/project/camb/) and [upstream repository/version notes](https://github.com/cmbant/CAMB); release-pinned [equations.f90](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/equations.f90), [Fortran readme](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/readme.html), [results.f90](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/results.f90); [parameter docs](https://camb.readthedocs.io/en/latest/model.html).
- Planck: [official PLA](https://pla.esac.esa.int/); [ESA Planck release page](https://www.cosmos.esa.int/web/planck); [Cobaya Planck likelihood docs](https://github.com/CobayaSampler/cobaya/blob/master/docs/likelihood_planck.rst); Planck Collaboration, [2018 results V, arXiv:1907.12875](https://arxiv.org/abs/1907.12875).
- Non-equivalent public IDE patch: [IDECAMB README](https://github.com/liaocrane/IDECAMB).

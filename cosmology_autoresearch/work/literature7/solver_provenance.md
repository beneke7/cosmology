# Solver and reproduction provenance: interacting-DE papers

**Status:** Literature-only, bounded primary-source audit, as of 2026-09-25. No solver was built or run; no parameter fit was reproduced.

## Finding

Both papers explicitly identify CAMB as the Boltzmann solver. Neither publishes an identifiable author-modified CAMB patch, its base CAMB commit/version, a complete run-specific Cobaya configuration, or initial-condition details. The Yang et al. IVS paper does publish the vacuum perturbation equations and states that vacuum DE has no perturbative contribution; this is enough to implement that prescription independently, but does not recover the authors' exact solver. Wang et al. name generalized interacting PPF and give the fluid equations, but do not give a complete PPF implementation/closure or initial conditions.

The only repository explicitly linked by the target paper is William Giarè's `wgcosmo`, cited for the Cobaya interface used with MCEvidence. At the inspected public tree (`80fa9742f10616bcaf55d5888a15e8f07e74cd25`, 2026-08-05), it contains reusable experiment/model YAMLs and evidence tooling, but no target-specific IVS/PPF configuration or CAMB patch; the repository has no releases or tags. Its YAML README describes files collected for the maintainer's personal organization. This is not the target authors' modified CAMB release.

| Item | Yang et al., arXiv:2609.05660v1 | Wang et al., arXiv:2609.03062v1 |
|---|---|---|
| Solver named by the paper | Modified CAMB, printed p. 9, Sec. III; refs. [153] (CAMB) and [154] (Cobaya). | Modified CAMB, printed p. 4, Sec. 3; CAMB citation follows the claim. |
| Exact modified patch/base commit | Not identified in the paper, its arXiv record, or the linked repository trail. | Not identified in the paper or its arXiv record; no paper-linked code repository found. |
| Complete Cobaya configuration | Not found. Paper gives component priors (Table I, p. 10), datasets and broad sampler description (pp. 8–10), but not the complete analysis YAML/settings. The `wgcosmo` link at p. 10 is specifically in the MCEvidence discussion. | Not found. The paper says it uses Cobaya MCMC, but supplies no run YAML or code/chain link. |
| Perturbation prescription | Partly explicit: synchronous gauge, zero anisotropic stress, no momentum transfer in the DM rest frame; for IVS, no DE perturbative contribution and CDM equations (19)–(20), printed pp. 4–5. | Partly explicit: fluid equations (8)–(10), printed pp. 3–4, and a generalized interacting PPF prescription described in prose on p. 4. The actual PPF closure/implementation choices are not specified to source-level completeness. |
| Initial conditions | No explicit initial-condition prescription found in the full text. | No explicit initial-condition prescription found in the full text. |
| Reproduction label | A new implementation using the stated equations and public CAMB would be independent, not source-faithful. | A new implementation of the stated model and PPF prescription would be independent unless the authors' code/config is obtained. |

The statement that the solver is CAMB is direct textual evidence, not a guess based on typical practice: Yang et al. say “modified version of the cosmology code CAMB” (p. 9); Wang et al. say “modified version of CAMB” (p. 4). No switch to CLASS is reported.

## Yang et al. IVS equations and missing run details

The exact public v1 PDF is available at [arXiv:2609.05660v1](https://arxiv.org/abs/2609.05660v1) and [its PDF endpoint](https://arxiv.org/pdf/2609.05660v1). The local copy at `work/literature5/arXiv_2609.05660v1.pdf` and the context copy `context/papers/interacting_de_desi_dr2_2026.pdf` are byte-identical, SHA-256 `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c` (4,524,530 bytes, 25 pages). Streaming the versioned arXiv PDF endpoint through `sha256sum` returned the same digest.

For IVS, the paper states (w_x=-1), “the DE has no contribution at the perturbative level,” and gives the synchronous-gauge CDM equations

```text
δ′c = −(θc + h′/2) + Γ a (ρx/ρc) δc
θ′c = −ℋ θc
```

as Eqs. (19)–(20), printed p. 4. The definitions and the sound-speed equation follow in Eq. (21), p. 5; the text says IVS is independent of (c^2_{s,x}). For the non-vacuum branches, Eqs. (15)–(18) give the fluid perturbations; (c^2_{s,x}=1) is fixed for the later statistical analyses (printed pp. 5–6). These statements specify a useful equation-level target, but the full text does not give the starting time, regular-mode/adiabatic initial values, or other code initial-condition settings.

The public prior table and data descriptions do not replace a complete Cobaya config. Table I (p. 10) gives the IVS prior Γ/H0 ∈ [−3, 3] and other cosmological parameter ranges. The paper names Planck 2018, DESI DR2, and PantheonPlus/Union3/DES-Dovekie combinations (pp. 8–10), CAMB, Cobaya, Gelman–Rubin, and GetDist (p. 9), and MCEvidence for evidences. It does not provide all likelihood YAMLs, exact likelihood versions and nuisance settings, sampler proposal/run settings, chains, or a solver patch.

The arXiv record links its TeX source, but that is not a solver release; an attempted source-archive fetch returned HTTP 406. The official APS accepted-paper page currently shows no supplementary-material download and says author-provided supplements will appear with the version of record, if provided. Neither observation establishes that no later release could appear.

## Wang et al. close comparison

The exact public v1 PDF is at [arXiv:2609.03062v1](https://arxiv.org/abs/2609.03062v1) and [its PDF endpoint](https://arxiv.org/pdf/2609.03062v1). Local copy: `work/literature5/arXiv_2609.03062v1.pdf`, SHA-256 `0e2b97a2fd988a6657e0c563bbbde7843529dea9c9ab8ac63dff88fefde356e5` (1,690,567 bytes, 17 pages); the versioned arXiv PDF endpoint returned the same digest.

This paper uses a different interaction, (Q=\beta H\rho_{de}), with transfer vector (Q^\nu=Q u_c^\nu), and evolving (w(a)). Eqs. (8)–(10), printed pp. 3–4, give the CDM/fluid perturbation equations and expose the conventional fluid singularity at (w=-1). Its p. 4 text says it uses generalized interacting PPF to cross the phantom divide and describes a large-scale relation smoothly connecting to small-scale behavior. It does not spell out the paper-specific PPF source modifications, closure settings, or initial conditions, and it names no base CAMB commit or complete Cobaya configuration. Its equations and PPF description are enough to define an independent implementation target only after implementation choices absent from the paper are supplied.

## Comparison with pinned public IDECAMB patch

The local archive `context/code/idecamb_4f1093d9efe46f28cf7e2acb4d07ae116ad5e075.tar.gz` is SHA-256 `c167aa0e0eae369a5b94816343ac6bdc8e3bbee32c5cf66fd116ad71450a6e44`. It is public [IDECAMB](https://github.com/liaocrane/IDECAMB), commit [`4f1093d9efe46f28cf7e2acb4d07ae116ad5e075`](https://github.com/liaocrane/IDECAMB/commit/4f1093d9efe46f28cf7e2acb4d07ae116ad5e075), dated 2023-06-20. Its README calls it a CAMB/CosmoMC patch, instructs copying into the Planck-2018 CosmoMC package, and says coupled-fluid perturbations use extended PPF; this is a separate, older public implementation, not either target paper's patch.

The background interaction and transfer direction can be mapped: `camb/equations_ppfi.f90` defines `H0rde` as (Q=\beta H_0\rho_{de}) (archive member lines 122–140), and `uc` as (Q_\mu=Q u_{\mu,c}) (lines 161–175). With β = Γ/H0 these match Yang et al.'s (Q=\Gamma\rho_x) and geodesic-CDM transfer choice. However, it is not equivalent unchanged at the perturbation level. IDECAMB's `LogicalPertur_CF` enables DE perturbations for the coupled (w=-1, \beta\ne0) fluid (lines 177–193); its default PPF branch at (w_{eff}\simeq-1) includes

```text
dgrhoe = −2 k² Kf(1) Gamma_PPF + gQ vT/k
```

in the DE density perturbation and adds this to the metric density source (`equations_ppfi.f90`, lines 2901–2917; the source calls the auxiliary PPF variable `Gamma`). Yang et al. instead state no IVS DE perturbative contribution and give Eqs. (19)–(20). IDECAMB initializes its PPF variable to zero and the interacting CDM velocity to zero in the documented CAMB initial routine (`modification_manual.f90`, lines 177–197), but that PPF initialization does not remove the interaction-sourced `gQ*vT/k` term. Its sample `test_ide.ini` selects `QForm_CF=4` (βH0ρc) and (w_0=-1.05), so the sample itself is not the target IVS setup.

**Disposition:** IDECAMB cannot be used unchanged as a source-faithful IVS reproduction. Its background law and transfer-vector option match after parameter mapping, but it requires an IVS-specific perturbation adaptation to implement Yang et al.'s Eqs. (19)–(20) and target closure. As-is, it is a different perturbation realization of the same background coupling. A hand-adapted IDECAMB or current public CAMB build would be an independent implementation, not the authors' reproduction.

## Bounded source trail and negative checks

- Read both arXiv version pages, PDFs, and local searchable texts. The 2609.05660 PDF itself contains the only article-linked repository URL: [wgcosmo](https://github.com/williamgiare/wgcosmo), ref. [157], printed p. 24; the body points to that repository for the Cobaya/MCEvidence interface, printed p. 10. The target paper itself describes the repository in that specific context.
- Inspected the `wgcosmo` public repository at main commit `80fa9742f10616bcaf55d5888a15e8f07e74cd25` (2026-08-05), the [`yamls/README.md`](https://github.com/williamgiare/wgcosmo/blob/80fa9742f10616bcaf55d5888a15e8f07e74cd25/yamls/README.md), its recursive tree, and release/tag endpoints. The tree has component YAMLs for DESI DR2, Planck/ACT, generic SN, and standard models, plus MCEvidence/Cobaya-wrapper files; it has no file path for IVS/interacting DE/PPF or either arXiv ID. Release and tag counts were both zero. No whole repository was downloaded.
- Searched the local full-text extracts and official GitHub web search for the two arXiv IDs/titles together with CAMB, code, GitHub, and interaction terms. No exact patch/config/chain release surfaced for either paper. The 2609.03062 PDF contains no GitHub/code URL; its CAMB modification is stated without provenance. These are bounded negative searches, not proof that an unindexed/private release does not exist.
- Opened the APS accepted-paper page for 2609.05660; it exposes no current supplement link. No public run chains were identified in the sources above.

**Next step if exact reproduction is needed:** request the authors' CAMB patch or a pinned commit, complete Cobaya YAMLs (including custom likelihood and nuisance configuration), chain settings/output, and initial-condition prescription. Until then, mark any new CAMB calculation as independent and compare its IVS equations explicitly to Eqs. (19)–(20).

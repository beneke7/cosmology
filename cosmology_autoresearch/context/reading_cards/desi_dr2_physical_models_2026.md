# Physical interpretations of DESI DR2: 2026 model-paper audit

- Audit date: 2026-09-24. Four current papers were read in full from local PDF/text; exact versions, hashes, equation/page locators, dataset combinations and access notes are consolidated in [`work/literature/physical_model_update.md`](../../work/literature/physical_model_update.md).
- Shared-data caution: every paper uses the same DESI DR2 BAO release; model/sample branches are not independent survey confirmations. Their CMB and SN combinations also overlap or are presented as alternatives.

## What is physically specified

| Model class | Explicit physical content | What a BAO-only profile can say | Main caveat |
|---|---|---|---|
| CPL/BA/EXP/LOG/JBP evolving (w(z)), [arXiv:2609.10567](https://arxiv.org/abs/2609.10567) | Kinematic flat-FLRW distance histories; the paper does not specify a dark-energy perturbation closure or phantom-crossing stability treatment in the reviewed text. | Compare background distances under chosen priors. | Do not promote preference for a crossing history into a stable scalar-field realization. |
| (Q=\Gamma\rho_x) interacting DM-DE, [arXiv:2609.05660](https://arxiv.org/abs/2609.05660) | Continuity equations and linear synchronous-gauge prescription are explicit; (c_{s,x}^2>0), often set to 1, with branch-dependent early stability conditions. | Screen the background response of the stated interaction. | In this paper (+\Gamma) transfers energy **from CDM to DE**. BAO does not test the perturbation closure or identify a microphysical interaction. |
| Metastable/decaying DE, [arXiv:2608.01844](https://arxiv.org/abs/2608.01844) | Three distinct channels: effective decaying-DE history; transfer to DM; transfer to dark radiation with a relativistic hierarchy. | Test the effective background history. | For transfer models (+\Gamma) means DE→daughter sector, opposite to the previous paper. Full-shape input is DESI DR1, not DR2. |
| Mass-varying neutrinos, [arXiv:2609.25090](https://arxiv.org/abs/2609.25090) | Conformal scalar-neutrino coupling motivates three background cases; one assumes adiabatic minimum tracking. | Constrain only the chosen background distance curve. | Neutrino/scalar perturbations, clustering and possible small-scale instabilities are omitted; adiabatic tracking needs posterior-wide verification. Figs. 1–2 label off-diagonal correlations synthetic, not sampler-derived. |

The interacting model is the most explicit one-rate perturbation closure among these papers, but implementing its BAO background alone would still not test that closure. The sign convention for “decay” is not portable between authors: preserve each continuity equation, not just the symbol Γ. The MaVaN paper's synthetic plot correlations do not automatically invalidate separately tabulated 1D constraints/evidence, but its joint-posterior visual degeneracies are not empirical evidence.

## Current Lyα boundary

DESI DR2 Results IV [arXiv:2607.27410](https://arxiv.org/abs/2607.27410) reports (D_M/D_H=4.572\pm0.046) from full-shape AP and (D_M/r_d=39.32\pm0.33, D_H/r_d=8.600\pm0.066), correlation +0.225, at (z_{\rm eff}=2.33). It analyzes the same DR2 Lyα sample family as Results I [arXiv:2503.14739](https://arxiv.org/abs/2503.14739), whose older compressed BAO row is already in our 13-vector. Results IV must replace/condition that row, not be appended independently absent the joint likelihood/covariance.

The Results IV validation companion [arXiv:2607.27411](https://arxiv.org/abs/2607.27411) validates BAO/AP behavior on mocks/splits, but reports significant mock bias in (f\sigma_8), which is excluded from its final analysis. Thus current DR2 Lyα provides strong geometry, not a validated growth/stability test. The official release page did not yield a directly usable replacement joint likelihood in this pass.

## Campaign decision

Keep the current 13-row full-covariance BAO analysis labeled as a background-distance screen. Use one dataset contract per fit. Before a physical-growth claim, acquire and audit a non-overlapping growth/full-shape likelihood and implement the exact model perturbations. See the detailed source/version ledger in the linked audit before reusing any paper's combined-data result.

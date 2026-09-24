# Next experiment

## 1. Physically explicit interacting-vacuum background (highest executable value)

Derive and implement the paper's (Q=\Gamma\rho_x), (w_x=-1) continuity system on the same 13-row DESI DR2 BAO mean/full covariance. The derivation must preserve the paper's sign convention (positive Γ transfers CDM energy to DE), keep baryons separately conserved, and state whether the baryon/CDM split is fitted as a nuisance or fixed from an external prior. No hidden CMB calibration: keep α=c/(H0 rd) free. Scope is background geometry only.

Pre-fit checks: Γ=0 recovers flat ΛCDM, units and signs match the paper's continuity equations, densities and H² stay positive throughout 0≤z≤2.33, and two independent integration/quadrature implementations agree. Compare finite bounds and one sensitivity to the baryon split. Fit with all covariance terms, retain optimizer statuses, and if the profile adds useful leverage run the same selected search on the existing 1,200 disjoint-seed fitted-flat-null vectors. Stop if the physical domain is unstable or the rate is unidentifiable without an unjustified prior.

Existing source: [Yang et al., arXiv:2609.05660](https://arxiv.org/abs/2609.05660); exact local PDF/text and equation-level notes in `context/papers/interacting_de_desi_dr2_2026.pdf` and `work/literature/physical_model_update.md`.

## 2. Official Lyα replacement likelihood (blocked pending a public artifact)

The completed 2026-09-24 official-page/directory, Results IV paper, Zenodo and linked-code audit found no machine-readable replacement/joint likelihood or cross-covariance. Eq. (26) has rounded summary values for an explicitly approximate 2D Gaussian; do not mistake it for an exact likelihood. The August 2026 unified-tracer supplement concerns galaxy/quasar BAO and is not a substitute. Full evidence is in `work/literature/desi_lya_release_probe.md`. Recheck official release links if they change, then pin/hash the artifact and validate overlap/order before use.

## Stop/quality criteria

- No BAO-only H0 claim: α leaves H0 and rd degenerate.
- A likelihood gain is exploratory until search-selection effects are calibrated on the same flat null.
- No perturbation/growth/stability claim without implementing and validating the paper's full perturbation closure.
- Every command stays bounded by an explicit task-level runtime; GPU only if an end-to-end benchmark shows a real advantage.

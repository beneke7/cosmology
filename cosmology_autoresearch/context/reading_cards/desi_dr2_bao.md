# DESI DR2 BAO — Results II

- Source/version: arXiv:2503.14738 v3; local full PDF/text read at `context/papers/desi_dr2_bao.pdf` (SHA-256 `1e82f26e4cc3901b16168cd147f252bfa804f9c3caad3f4f7e3532640d237841`).
- Question: what do the DR2 galaxy, quasar and Ly-alpha BAO measurements imply for late-time expansion and extensions to LCDM?
- Governing relations: Eqs. (1)–(5), pp. 4–5, use DM/rd and DH/rd; Eqs. (8)–(10) define evolving-DE/CPL density. The released compressed measurement is a 13-entry vector with full covariance.
- Data/likelihood: the local mean and covariance are pinned to CobayaSampler/bao_data commit `bb0c1c9009dc76d1391300e169e8df38fd1096db`; exact row order is required. The vector includes a Ly-alpha block at z_eff=2.33.
- Reported result: Table VI, p. 24, gives conditional profile Delta-chi-square comparisons for specified CMB/SN combinations. These are not Bayesian odds or a stand-alone detection.
- Assumptions/limits: the default CMB bundle includes Planck spectra and ACT/Planck lensing. SN branches overlap and are not independent; a free BAO scale constrains H0*rd, not H0 alone.
- Experiment implication: compare LCDM, constant-w and CPL on the same vector/full covariance; retain a free scale and call background-only minima profile screens, not posterior reproductions.
- Related files: `experiments/campaign_seed_bao/result.json`, `experiments/bao_robustness/result.json`, `work/literature/reading_map.md`.

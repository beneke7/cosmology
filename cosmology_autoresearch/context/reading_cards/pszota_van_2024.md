# Pszota–Van thermodynamic-gravity galaxy application

- Source/version: user-provided Physics of the Dark Universe 46 (2024) 101660, DOI `10.1016/j.dark.2024.101660`; local full PDF `context/user_provided/pszota_van_2024.pdf` (SHA-256 `55cbfbbb0a916ea2ce0b70523334ee301897695860fb4bdb5fd2278ab47b62e`). Eq. (1) was visually checked. Matching arXiv:2306.01825 v3 is also locally read but versions are not asserted byte-identical.
- Equation: Eq. (1), p. 2, is partial_t phi = D(laplacian(phi) - K|grad(phi)|^2 - 4 pi G rho), D=ell^2/tau. Spherical/stationary equations are Eqs. (2), (5)–(8), pp. 3–4; vacuum branch Eqs. (9)–(17), pp. 4–5.
- Exact transform: for prescribed rho, constant K and positive u=exp(-K phi), chain rule yields partial_t u = D[laplacian(u)+4 pi G K rho u]; the stationary fixed-source equation is linear. This is not a novelty claim, stability proof or coupled-matter solution.
- Galaxy fit/limits: the NGC 3198 numerical scheme fits K and stellar mass-to-light normalization in a sphericalized source; Eqs. (41)–(42), p. 8, impose observed endpoint velocities as derivative boundary conditions. It is conditional, not an isolated-boundary/holdout prediction. SPARC force components do not uniquely specify a 3D density.
- Experiment implication: derive transformed boundary data, test positivity and K-to-zero limits, and specify a physical outer boundary independent of held-out observed velocities before predictive fitting.
- Local checks: `runs/tg_audit.json`, `work/theory/derivation.md`, `work/theory/independent_check.py`.

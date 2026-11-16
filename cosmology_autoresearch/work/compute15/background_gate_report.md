# Gate 1: radiation-inclusive homogeneous background

**Result: PASS for the declared background benchmark.** The standalone implementation checks the equations in Yang et al., *Do DESI-DR2 BAO data imply a coupling of dark matter and dark energy?*, Eqs. (3)–(8), PDF pp. 2–3. It does not compare or reproduce the paper's likelihood scores.

The primary PDF is `context/papers/interacting_de_desi_dr2_2026.pdf` (SHA-256 `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`); the local text copy hash is `a99a4afdd858cdce3cbc2c41f0a13098c30c361dc67891331ee7d9d5f824868f`. The implementation is independent and self-contained in [background_gate.py](background_gate.py); full machine-readable measurements, assumptions, and hashes are in [background_gate.json](background_gate.json).

With (n=\ln a), (g=\Gamma/H_0), and (E=H/H_0), I integrated

\[
\rho_x'=g\rho_x/E,\quad
\rho_c'=-3\rho_c-g\rho_x/E,\quad
\rho_b'=-3\rho_b,\quad
\rho_\gamma'=-4\rho_\gamma,\quad
\rho_\nu'=-4\rho_\nu,
\]

using (E^2=\rho_b+\rho_\gamma+\rho_\nu+\rho_c+\rho_x). The transfer sign follows the paper: Γ>0 moves energy from CDM to vacuum. As an independent equivalent ODE, the script evolves (\rho_d=\rho_c+\rho_x), using (\rho_d'=-3\rho_c). Separate adaptive quadratures check both this dark-sector continuity identity and

\[
\frac{d(a^3\rho_c)}{d\ln a}=-g\,a^3\rho_x/E.
\]

Assumptions are pinned in the JSON: (H_0=67.4\) km s⁻¹ Mpc⁻¹, (\Omega_{m0}=0.315), (\Omega_{b0}=0.049), (\Omega_{c0}=0.266), (T_{\rm CMB}=2.7255\) K, (N_{\rm eff}=3.044), photon density from blackbody radiation, flat closure, and densities normalized to today's critical density. Neutrinos are treated as a separately conserved massless component, (\rho_\nu\propto a^{-4}). The paper's footnote specifies \(\sum m_\nu=0.06\) eV; this benchmark therefore does not claim its exact massive-neutrino background. The tested rates are (g=0,+0.46666,-0.46666), corresponding to (\Gamma=0,+31.452884,-31.452884\) km s⁻¹ Mpc⁻¹.

Results over (0\le z\le10^{10}) (2,049 output points, DOP853, maximum step 0.02 in the fine run):

- (\rho_c,\rho_x,E^2) stayed positive on all 2,049 output points for all three rates; minimum sampled (\rho_c=0.266), minimum (\rho_x=0.4493883), and minimum (E^2=1). This also implies continuous positivity in the integrated interval: at any zero of (\rho_c), (\rho_c'=-g\rho_x/E) has the same strict sign, so a second crossing is impossible; positive endpoint values exclude a single crossing. For (g=0), (\rho_c=\Omega_{c0}a^{-3}) exactly. The scalar vacuum equation preserves positive (\rho_x), and positive component densities give positive (E^2). At (z=10^{10}), photons plus massless neutrinos supply more than 0.9999995 of (E^2).
- At (g=0), the analytic flat LCDM solution is recovered with maximum relative errors (1.45\times10^{-14}) in (\rho_c), (2.04\times10^{-15}) in (E^2), and zero in constant (\rho_x).
- Across all rates, the independent ρd formulation agrees with direct ((\rho_c,\rho_x)) integration to at most (5.04\times10^{-15}). Dark-sector continuity and the exact CDM comoving-density identity have maximum relative quadrature residuals (7.34\times10^{-15}) and (1.44\times10^{-14}), respectively.
- Tightening tolerance alone, maximum step alone, and both together changed (\rho_c,\rho_x,E^2) by at most (3.24\times10^{-10}) relative to the coarse run. (\Gamma/H_0) is dimensionless; (\Gamma) and (H_0) were also independently converted from km s⁻¹ Mpc⁻¹ to s⁻¹ and their ratio checked against (g).

The run used one CPU (affinity pinned to CPU 0), one worker, no GPU, and completed in about 2.3 seconds. Script SHA-256: `c07ede03539d08702f24db69140badd8361c859e76be608fb246ad1867c91d9c`.

This passing homogeneous-background check is **not** perturbative stability or CMB validation. It does not test perturbation evolution, CMB spectra, massive-neutrino evolution, sound-horizon consistency, structure growth, or the paper's full likelihood/evidence analysis.

# IVS early-time perturbation IC audit

## Result

For the usual radiation-era branch with positive early comoving CDM density, the standard synchronous-gauge adiabatic mode remains valid at leading order. The interaction term in the paper's CDM contrast equation is suppressed as \(A_5=K_5a^5\) in \(d/d\ln a\); on a fixed leading radiation metric its first correction to the regular CDM solution is \(-A_5 h/14=O(a^7)\). In the full background, matter corrections to \(E\) generally generate \(A_6=K_6a^6\), whose response is \(O(a^8)\); this is not included in the pure-\(A_5\) truncation. This is an asymptotic statement, not a full Einstein-Boltzmann solution.

There is also a tuned *early-time-positive* edge branch with zero early comoving CDM density, possible only for \(\Gamma<0\). There the interaction coefficient tends to \(-5\), not zero. The regular solution of the displayed Eq. (19) is \(\delta_c=-h/7\), whereas imposing the usual CDM-radiation equal-\(\zeta\) criterion requires \(\delta_c=+h/3\). Thus Eq. (19) and that imposed IC criterion are incompatible on this edge branch. The radiation-era Einstein constraints themselves remain radiation-dominated: \(\rho_c/\rho_r=O(a^6)\) there, so either CDM contrast contributes only at still higher order. Whether this asymptotic branch stays positive through the whole history requires a separate global background check. This result is limited to the displayed perturbation equation and imposed IC criterion; it is not a demonstrated inconsistency of a fully specified covariant vacuum-transfer completion and is not a claim about the authors' posterior.

## Paper-specified facts and conventions

The source paper takes \(w_x=-1\), \(Q=\Gamma\rho_x\), and separately conserves baryons, radiation, and neutrinos. Its Eqs. (4)-(5) give

\[
\dot\rho_x=\Gamma\rho_x,\qquad \dot\rho_c+3H\rho_c=-\Gamma\rho_x,
\]

so \(\Gamma>0\) means CDM transfers energy to vacuum. Its synchronous-gauge Eqs. (19)-(20), with primes denoting conformal time, are

\[
\delta_c'=-(\theta_c+h'/2)+a\Gamma\frac{\rho_x}{\rho_c}\delta_c,
\qquad \theta_c'=-\mathcal H\theta_c.
\]

The regular CDM-comoving solution has \(\theta_c=0\). With \(n=\ln a\), \(E=H/H_0\), and \(g=\Gamma/H_0\), the exact conversion is \(a\Gamma/\mathcal H=\Gamma/H=g/E\), and the density equation becomes

\[
\frac{d\delta_c}{dn}=-S+A\delta_c,\quad
S=\frac{1}{2}\frac{dh}{dn},\quad
A=\frac{g}{E}\frac{\rho_x}{\rho_c}.
\]

These are facts from the paper, distinct from the asymptotic derivation below. Local PDF SHA-256: 0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c.

## Derived background and perturbation scalings

Assume a nonzero conserved relativistic density, expanding flat background, and positive vacuum density. Then \(E\sim\sqrt{\Omega_{r0}}a^{-2}\), so \(\Gamma/H=O(a^2)\) and the past integral of \(\Gamma/H\) converges. For \(X_*\equiv\lim_{a\to0}\rho_x>0\),

\[
\rho_x=X_*\left[1+\frac{g}{2\sqrt{\Omega_{r0}}}a^2+O(a^3)\right].
\]

Writing \(D=\lim_{a\to0}a^3\rho_c\), the CDM background equation gives

\[
a^3\rho_c=D-\frac{gX_*}{5\sqrt{\Omega_{r0}}}a^5+\cdots,
\qquad
\rho_c=D a^{-3}-\frac{gX_*}{5\sqrt{\Omega_{r0}}}a^2+\cdots.
\]

Positivity excludes \(D<0\). If \(D>0\), either sign of \(g\) is regular and

\[
A\sim\frac{gX_*}{\sqrt{\Omega_{r0}}D}a^5.
\]

Thus the standard growing radiation mode \(h=\mathcal C(k\tau)^2\propto a^2\) has \(S=h\), and its leading CDM value \(-h/2\) is unchanged. Solving Eq. (19) one order farther on this fixed leading metric gives

\[
\delta_c=-\frac{h}{2}-\frac{A_5 h}{14}+\cdots,\qquad A_5=K_5a^5.
\]

The displayed \(O(a^7)\) correction is the leading interaction response. If the full-background coefficient is written \(A=A_5+A_6+\cdots\), with \(A_6=K_6a^6\) generally induced by the \(O(a)\) matter correction to \(E/(\sqrt{\Omega_{r0}}a^{-2})\), then its response is generally \(O(a^8)\). Only in the idealized pure-\(A_5\), fixed-\(h\propto a^2\) truncation is the next iteration residual \(O(A_5^2h)=O(a^{12})\). That \(O(A_5^2h)\) remainder applies only to this truncation and is not a bound on the full-background remainder. This is not an assertion that the fixed-metric series alone supplies all higher-order adiabatic ICs for every species.

If \(D=0\), then \(g>0\) makes the leading CDM density negative and is excluded; \(g=0\) leaves no CDM, so Eq. (19)'s density-contrast form is undefined. For \(D=0,g<0\),

\[
\rho_c\sim-\frac{gX_*}{5\sqrt{\Omega_{r0}}}a^2>0,
\qquad A\longrightarrow -5.
\]

With \(S=h\propto a^2\), Eq. (19) then has a singular homogeneous solution \(\propto a^{-5}\) and a unique regular particular solution \(\delta_c=-h/7\). The standard massless-radiation synchronous adiabatic relations give \(\delta_\gamma=-2h/3\) and \(\beta_\gamma=d\ln\rho_\gamma/d\ln a=-4\). Here \(\beta_c=+2\). Using \(\zeta_i=-\psi_{\rm slice}-\delta_i/\beta_i\) (the common metric term cancels in differences), equal \(\zeta_i\) demands

\[
\delta_c=\frac{\beta_c}{\beta_\gamma}\delta_\gamma=+\frac{h}{3},
\]

not \(-h/7\). For the regular Eq. (19) solution,

\[
\frac{\zeta_c-\zeta_\gamma}{h}
=-\frac{1}{h}\left(\frac{\delta_c}{\beta_c}-\frac{\delta_\gamma}{\beta_\gamma}\right)=\frac{5}{21}.
\]

This mismatch is at leading order on the \(D=0,g<0\) edge; it cannot be repaired by the subleading CDM contribution to the radiation metric. The CDM-comoving frame should be understood as a limit from \(a>0\), since \(\rho_c\to0\) at the boundary.

There is a separate vacuum subtlety even on the ordinary \(D>0\) branch. The paper supplies no independent \(\delta_x\) equation and says the IVS vacuum has no perturbative contribution. If that closure is implemented as \(\delta_x=0\), then for nonzero \(\Gamma\), \(\beta_x=d\ln\rho_x/d\ln a=\Gamma/H=O(a^2)\). A strict all-species equal-\(\zeta\) time shift would require \(\delta_x=\beta_x h/6=O(a^4)\), not zero. Accordingly, the standard IC result established above is the leading adiabatic relation for radiation, baryons, and the ordinary CDM branch; the source paper does not specify a strict all-species vacuum IC prescription. At \(\Gamma=0\), as for a true cosmological constant, \(\zeta_x\) is undefined and is not an independent adiabatic constraint.

## Standard radiation mode, gauge, and neutrinos

The primary reference is Ma & Bertschinger, [*Cosmological Perturbation Theory in the Synchronous and Conformal Newtonian Gauges*](https://arxiv.org/abs/astro-ph/9506072), especially Eq. (96). Their fastest-growing physical super-horizon synchronous solution has \(h=\mathcal C(k\tau)^2\), \(\delta_\gamma=\delta_\nu=-2h/3\), \(\delta_b=\delta_c=3\delta_\gamma/4=-h/2\), and \(\theta_c=0\). They also give

\[
\sigma_\nu=\frac{4\mathcal C}{3(15+4R_\nu)}(k\tau)^2,\quad
\theta_\nu=\frac{23+4R_\nu}{15+4R_\nu}\theta_\gamma,\quad
\eta=2\mathcal C-\frac{5+4R_\nu}{6(15+4R_\nu)}\mathcal C(k\tau)^2.
\]

They identify the other first two synchronous solutions as gauge modes and the second physical radiation-era mode as subdominant/decaying later; an IC implementation must remove residual gauge modes and not mistake an added CDM entropy mode for the physical adiabatic mode.

Retain the free-streaming neutrino hierarchy and its shear in a CMB solver unless the authors specify otherwise. The paper's “zero anisotropic stress” phrase is attached to the \(w_x\ne-1\) perturbation system; the paper does not specify the IVS *total-radiation stress closure*. With the benchmark \(N_{\rm eff}=3.044\), the early relativistic fraction is \(R_\nu\simeq0.40874\), which affects the \(\eta\), shear, and neutrino-velocity coefficients even though the leading \(h,\delta_\gamma,\delta_\nu\) powers above are unchanged. The source paper pins \(\sum m_\nu=0.06\,\mathrm{eV}\) and \(N_{\rm eff}=3.044\); as \(a\to0\) the masses are negligible. A finite-start Boltzmann implementation still needs a pinned massive-neutrino phase-space/thermal prescription and mass splitting. The paper does not state that prescription or a scalar IC config in the visible material.

The radiation normalization used only to anchor the asymptotics is \(H_0=67.4\,\mathrm{km\,s^{-1}\,Mpc^{-1}}\), \(\Omega_{b0}=0.049\), \(\Omega_{c0}=0.266\), \(T_{\rm CMB}=2.7255\,\mathrm K\), flat closure, and \(\rho_\nu/\rho_\gamma=(7/8)(4/11)^{4/3}N_{\rm eff}\) in the asymptotic relativistic limit. These are benchmark conventions, not a fit to the paper's headline scores; the exponents do not depend on these late-time values.

## Executable check and scope

The script early_ic_asymptotic_check.py independently evaluates the pinned radiation density using SI constants, the \(D>0\) and \(D=0\) leading series, log-slopes \(A\sim a^5\) and \(A+5\sim a^2\), signs, and the regular/same-slice equal-\(\zeta\) residuals. It uses no production-profile, compute15, or Boltzmann code. The machine-readable results are in early_ic_asymptotic.json.

The bounded invocation was:

    env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 taskset -c 0 python3 -B scripts/run_bounded.py --seconds 600 --state work/theory10/no_campaign_deadline.json -- /usr/bin/time -f 'inner_elapsed=%e sec' python3 -B work/theory10/early_ic_asymptotic_check.py

The process finished successfully in 0.02 s of inner wall time (the script's own timer was about 0.0025 s). The state path was intentionally absent, so the explicit 600 s task cap governed the run.

SHA-256 provenance:

- Script: cad766208f030d67fa15fefef92b306efc1d111a2dd3a293ccd1bc9eb9634516
- Machine-readable JSON: fff0703c5fa8036054b0226f2e9ea78a30c35031d7cea862812ce0c286bee063
- Source PDF: 0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c
- Local source text: a99a4afdd858cdce3cbc2c41f0a13098c30c361dc67891331ee7d9d5f824868f
- Solver-path context: 810ebc06aa9cbf22aeaf9fdb974fce6b9b7c35c41a6b70336a0392ff824e15b3

This is an asymptotic algebra and fixed-metric check only. It is **not** a perturbative stability proof, CMB validation, or full-model/Boltzmann result. No CAMB/Cobaya build, posterior, likelihood, or GPU was run.

# Independent critical audit: 2026 higher-gradient thermodynamic gravity

**Scope and source identity.** Audited the complete four-page local v1 PDF, `context/papers/tg_higher_gradient_2026.pdf` (arXiv:2609.00317v1; SHA-256 `c0992fc15699bcd2c7abd5a4b044ace25203c2793c0cca42cb775888945fe165`). Page and equation locators below refer to that version. I also checked the relevant cited primary sources: Lazar's static gradient-gravity derivation [arXiv:2009.09846](https://arxiv.org/abs/2009.09846), Szücs & Ván's weakly nonlocal thermodynamic framework [arXiv:2504.07296v2](https://arxiv.org/abs/2504.07296v2), Lee et al.'s torsion-balance result [arXiv:2002.11761](https://arxiv.org/abs/2002.11761), and the full supplied 2024 galaxy paper [arXiv:2306.01825v3](https://arxiv.org/abs/2306.01825v3). Ván & Abe's earlier paper [arXiv:1912.00252](https://arxiv.org/abs/1912.00252) was checked at abstract/source-index level, not independently re-derived. No fit or external write occurred.

## Independently checked / confirmed

- **Static operator, signs, and dimensions (p. 1, Eqs. 2–3; p. 2, Eq. 7).** Varying
  \[
  \varepsilon=\rho\phi+\frac{(\nabla\phi)^2+\ell_1^2(\Delta\phi)^2+\ell_2^4\lVert\nabla^3\phi\rVert^2}{8\pi G_N}
  \]
  gives \(\mathcal G=-\delta_\phi\varepsilon=[\Delta\phi-\ell_1^2\Delta^2\phi+\ell_2^4\Delta^3\phi]/(4\pi G_N)-\rho\); the signs in the sixth-order stationary equation are consistent with Lazar's Eq. (9). With \([\phi]=L^2T^{-2}\), each term in \(\varepsilon\) has energy-density units; \(\ell_1,\ell_2\) are lengths, while the separate transport scale \(\ell^2/\tau\) has units \(L^2/T\). In Eq. (7), \([\mathcal G]=M/L^3\), \([\Pi]=M/(LT^2)\), so \([K]=T^2/L^2\), matching the 2024 paper.

- **Entropy-production closure and sign of the relaxation (p. 2, Eqs. 5–10).** Let \(\theta=\nabla\cdot v\). Eq. (6) yields \(\theta=(-\Pi-l_{21}\mathcal G)/l_2\), hence
  \[
  \dot\phi=\frac{\det L}{l_2}\,\mathcal G-\frac{l_{12}}{l_2}\Pi.
  \]
  The stated definitions in Eq. (7) reproduce both coefficients and signs. The entropy bound is positivity of the symmetric part of \(L\), \(l_1l_2\ge(l_{12}+l_{21})^2/4\). For \(\Pi=0\), \(D\equiv\ell^2/\tau=\det L/(4\pi G_N l_2)>0\), and \(\mathcal G=-\delta_\phi\varepsilon\), giving the energy gradient flow \(\dot\phi=-4\pi G_ND\,\delta_\phi\varepsilon\). This is first order in time; the Ostrogradsky theorem for nondegenerate higher-time-derivative Lagrangians does not apply to this chosen dynamics. That is a valid conditional “no ghost in this relaxation model” result.

- **Pressure identity (p. 2, Eqs. 3–4, 8).** Direct 1-D differentiation of the printed \(P^g\), using \(M=A-B'+C''\), \(N=B-C'\), cancels the \(\phi''\), \(\phi'''\), and \(\phi''''\) terms and leaves \((P^g)'=\rho(\partial_\rho\varepsilon)'+\mathcal G\phi'\), as Eq. (8) states. The same cancellation is the standard higher-gradient stress identity. On shell for Eq. (2), \(\partial_\rho\varepsilon=\phi\), so \(\nabla\cdot P^g=\rho\nabla\phi\). The paper correctly adds (p. 2, after Eq. 8) that test-particle Euler motion also requires mechanical perfection and an integrable thermal-pressure/entropy force.

- **Field-mode spectrum (p. 2, Eq. 10).** Linearization at fixed \(\rho\) gives
  \[
  \lambda(k)=-D Q(k),\qquad Q(k)=k^2+\ell_1^2k^4+\ell_2^4k^6.
  \]
  For positive mobility and nonnegative gradient-energy coefficients this is damped for every **nonzero** real wavenumber, including the complex-pole regime: the Fourier symbol \(Q(k)\) stays positive. I independently checked the claimed full-fluid Jeans result for a perfect barotropic background with \(c_s^2>0\), \(\Pi=0\), and the printed pressure identity. Its characteristic polynomial is
  \[
  (\lambda+DQ)(\lambda^2+c_s^2k^2)-4\pi G_N\rho_0Dk^2=0.
  \]
  **Correction/retraction:** the earlier draft's “above \(k_J^2=4\pi G_N\rho_0/c_s^2\)” stability threshold was overbroad; it omitted the gradient terms. Expanding the cubic gives
  \[
  \lambda^3+DQ\lambda^2+c_s^2k^2\lambda+Dk^2(c_s^2Q-4\pi G_N\rho_0)=0.
  \]
  For \(k>0\), \(D>0\), \(c_s^2>0\), and \(\rho_0>0\), set \(A=DQ\), \(B=c_s^2k^2\), and \(C=Dk^2(c_s^2Q-4\pi G_N\rho_0)\). Then \(A,B>0\) and direct Routh–Hurwitz gives stability iff \(C>0\): \(AB-C=4\pi G_N\rho_0Dk^2>0\), so the exact threshold is \(c_s^2Q(k)=4\pi G_N\rho_0\). If \(C<0\), the coefficient signs are \(+,+,+,-\), hence Descartes' rule gives exactly one positive real (growing) root. The classical \(k_J^2=4\pi G_N\rho_0/c_s^2\) boundary follows only in the zero-gradient limit \(\ell_1=\ell_2=0\); for nonzero gradient lengths, \(Q=k^2+\ell_1^2k^4+\ell_2^4k^6\) shifts the boundary. At equality there is a neutral root and the remaining quadratic has roots with negative real parts. This verifies the result for that perfect barotropic closure, not every possible dissipative/thermal completion.

- **Regularity sum rule (p. 3, Eqs. 11–13).** For the decaying point-source Green function of the linear stationary operator, \(\phi=-GM/r[1+\sum_i\alpha_i e^{-r/\lambda_i}]\). Finiteness at \(r=0\) requires \(1+\sum_i\alpha_i=0\). In the real-root case \(a+b=\ell_1^2\), \(ab=\ell_2^4\), \(\lambda_{a,b}=\sqrt{a,b}\), with \(\alpha_a=-a/(a-b)\le-1\), \(\alpha_b=b/(a-b)\ge0\); these follow by partial fractions. The amplitudes are correlated, not independent. The “Lee–Wick-like” label is supported only as an analogy for complex-conjugate poles and oscillatory static corrections; it does not establish a Lee–Wick quantum theory or its unitarity.

## Qualifications / likely overstatements

- **Strict negativity and concavity (p. 2, Eqs. 9–10).** Eq. (10) prints \(\lambda(k)<0\) “for all \(k\),” but \(\lambda(0)=0\): the constant-potential mode is neutral. Strict damping holds for \(k>0\) if \(D>0\). The calculation is a fixed-density field-sector test; it does not prove concavity or nonlinear stability of the entire fluid state space. Also, strict \(\ell_2^4>0\) is not necessary in the single-length limit discussed on p. 3: \(\ell_2=0\) still gives \(Q(k)>0\) for nonzero \(k\).

- **Lyapunov wording (p. 2, discussion after Eq. 9).** Local nonnegative production \(\sigma_s\ge0\) is confirmed from Eq. (5), but a local entropy balance has an entropy-flux divergence; this alone does not make the entropy *density* pointwise monotone or eliminate boundary conditions. For fixed source and appropriate periodic/decay or variational boundary conditions, the global field functional does obey
  \[
  \frac{d}{dt}\int\varepsilon\,d^3x=-4\pi G_ND\int(\delta_\phi\varepsilon)^2d^3x\le0.
  \]
  At fixed \(\rho,e\) and positive \(T\), the corresponding global entropy functional increases, subject to its boundary flux. This supports a conditional global Lyapunov statement, not the paper's stronger “entropy density”/“no boundary condition” formulation.

- **Onsager parameters and setting \(K=0\) (p. 2, Eqs. 6–7).** Under the paper's highlighted zero-direct-field-dissipation branch \(l_1=0\), positivity forces \(l_{21}=-l_{12}\), so \(\det L=l_{12}^2\) and the printed \(K=l_{12}/(6\det L)=1/(6l_{12})\) for finite nonzero coupling. Thus \(K\) cannot independently be set to zero while retaining that branch's finite positive mobility. A separate branch \(l_{12}=0,l_1>0\) permits \(K=0\) and relaxation, but then field relaxation is directly dissipative. Since the subsequent field-only stability discussion also sets \(\Pi=0\), the \(K\Pi\) term vanishes there either way; the paper should specify which constitutive branch it means. This is an ambiguity in the claimed origin/independence of the cross-effect, not a sign error in Eq. (7).

- **Laboratory number (p. 3, Eq. 14).** Lee et al. report that a *single* gravitational-strength Yukawa term (including both signs in their analysis) with \(|\alpha|=1\) must have \(\lambda<38.6\,\mu\mathrm m\) at 95% confidence; they measured gaps from 52 μm to 3.0 mm. This directly gives \(\ell_1<38.6\,\mu\mathrm m\) in the single-scale \(\ell_2=0\) limit, where \(\alpha=-1,\lambda=\ell_1\). It is not by itself a published joint limit on the paper's two-range or complex-root force template. In the real-root case, \(\ell_1^2=\lambda_a^2+\lambda_b^2\), so bounding the largest range \(\lambda_a\) does not set \(\ell_1=\lambda_a\) (at most \(\ell_1\le\sqrt2\lambda_a\)); correlated two-range torques and complex-root oscillations need a direct fit to the apparatus data. Thus “about 40 μm” is a sound single-scale/decoupled-short-range estimate, not an established bound over the full parameter space.

- **Novelty and physical scope.** The static sixth-order operator, regular point-source Green function, and correlated bi-Yukawa amplitudes are already in Lazar 2020 (see its point-mass Eqs. 38–39); \(\sum_i\alpha_i=-1\) is the algebraic cancellation condition behind that result, not a new static prediction. The plausible new contribution is the thermodynamic derivation/relaxational dynamics with third spatial derivatives of \(\phi\): the cited Szücs–Ván v2 framework uses second-order weakly nonlocal state variables, while this Letter includes third-order dependence. Priority beyond this cited lineage is not established by the targeted check. “No variational principle can produce” is too broad without restricting it to an ordinary conservative single-field action: nonconservative actions with doubled fields are established for dissipative systems [Galley, Tsang & Stein](https://arxiv.org/abs/1412.3082). The first-order model has no Ostrogradsky mode, but this establishes neither a relativistic completion nor universal physical stability. The assumptions are Newtonian instantaneous gravity, the selected positive quadratic gradient energy, a local entropy/equation of state with positive temperature, Onsager closure, and a specified matter/thermal boundary-value problem.

## Keep the 2026 and 2024 theories distinct

The 2024 galaxy model (Pszota–Ván, arXiv:2306.01825v3, Eqs. 1–2) is a **first-gradient, nonlinear** relaxation equation, \(\dot\phi=(\ell^2/\tau)[\Delta\phi-4\pi G_N\rho-K(\nabla\phi)^2]\), with fitted material cross-coupling \(K\); it was a one-galaxy NGC 3198 proof of concept. The 2026 paper instead introduces a **linear sixth-order spatial operator** using \(\ell_1,\ell_2\), derives a first-order field relaxation from weak nonlocal thermodynamics, and explicitly sets the bulk-pressure coupling term aside for its main analysis. It is not a validation or extension of the 2024 galaxy fit, and its approximately 40 μm laboratory discussion tests the short-range higher-gradient sector, not the galaxy-scale K mechanism.

## Single most discriminating physical experiment

Run an improved Eöt–Wash-style patterned-mass torsion balance with a minimum calibrated gap in the 5–10 μm range and fit the measured separation-dependent torque directly to the **full correlated two-range/complex-root potential** (not a one-Yukawa limit curve), profiling both \(\ell_1\) and \(\ell_2\) and nuisance geometry. The single-scale model predicts a fixed-sign Yukawa correction with \(\alpha=-1\), which weakens Newtonian attraction below \(\lambda\); the two-scale model predicts the linked sum above. A direct differential torque scan across 5–100 μm most cleanly tests those shape predictions against Newtonian gravity.

### Primary sources checked

- Pszota & Ván 2026, exact v1: [arXiv:2609.00317v1](https://arxiv.org/abs/2609.00317v1).
- Lazar 2020, sixth-order static gradient gravity and Green functions: [arXiv:2009.09846](https://arxiv.org/abs/2009.09846).
- Szücs & Ván 2026 v2, thermodynamic weakly nonlocal fluid framework: [arXiv:2504.07296v2](https://arxiv.org/abs/2504.07296v2).
- Lee et al. 2020, 52 μm torsion-balance test and 38.6 μm single-Yukawa bound: [arXiv:2002.11761](https://arxiv.org/abs/2002.11761), [Phys. Rev. Lett. 124, 101101](https://doi.org/10.1103/PhysRevLett.124.101101).
- Pszota & Ván 2024, distinct galaxy (K) model: [arXiv:2306.01825](https://arxiv.org/abs/2306.01825), [journal article](https://doi.org/10.1016/j.dark.2024.101660).
- Ván & Abe 2022, earlier thermodynamic gravity: [arXiv:1912.00252](https://arxiv.org/abs/1912.00252).
- 2026 source cited for Newtonian energy density distinction: Trasarti-Battistoni, Pszota & Ván, [arXiv:2605.02976](https://arxiv.org/abs/2605.02976) (bibliographic/abstract check only).
- Galley, Tsang & Stein 2015, nonconservative variational principle for classical fields: [arXiv:1412.3082](https://arxiv.org/abs/1412.3082).

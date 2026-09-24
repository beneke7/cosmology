# Constant-\(K\) transform and NGC 3198 readiness

## Status and scope

- **Reproduced:** `scripts/audit_tg.py` runs in the project `.venv`. Its printed-convention uniform-density residual for the paper's Eq. 23 coth profile is `0.268211`, versus `4.25e-12` for the cot profile. Its Gaussian-source nonlinear/transformed comparison differs by at most `4.80e-11` in $g$.
- **Independently checked:** `independent_check.py` solves the transformed radial boundary-value problem by collocation and the original nonlinear equation by a separate initial-value solve. Their $g$ fields agree to $2.84\times10^{-11}$ absolute and $3.09\times10^{-10}$ relative on the sampled domain; $u$ stays positive.
- **Blocked for prediction:** the supplied NGC 3198 material supports at most a conditional spherical-surrogate reproduction. Its three-dimensional baryon source and a data-independent outer boundary are not specified. No galaxy fit was attempted.

Sources inspected: Pszota & Ván, *Physics of the Dark Universe* 46 (2024) 101660, DOI [10.1016/j.dark.2024.101660](https://doi.org/10.1016/j.dark.2024.101660), supplied as `context/user_provided/pszota_van_2024.pdf` (SHA-256 `55cbfbbb0a916ea2ce0b70523334ee301897695860fb4bdb5bfd2278ab47b62e`); and the provided `scripts/audit_tg.py` (SHA-256 `a980a8bbc544f2463b9df75d081867e9b4f2b4ac5b2159730843a77f3fb1868b`). The paper's equation and page references below are to its printed numbering/pages.

## Transformation, signs, and units

Write the printed equation (paper Eq. 1, p. 2) as

\[
\partial_t\phi=D\bigl(\Delta\phi-K|\nabla\phi|^2-S\bigr),
\qquad S=4\pi G\rho,
\qquad u=e^{-K\phi},
\]

where $D=\ell^2/\tau$ and $K$ is constant in space and time. Direct differentiation gives

\[
\partial_tu=-Ku\partial_t\phi,
\qquad \nabla u=-Ku\nabla\phi,
\qquad \Delta u=-Ku\Delta\phi+K^2u|\nabla\phi|^2.
\]

Therefore

\[
\boxed{\partial_tu=D(\Delta u+KSu)}.
\]

The plus sign in $+KSu$ follows from the minus source in the printed $\phi$ equation and the minus sign in the exponent. At stationarity with prescribed $S(\mathbf{x})$, the equation is the linear homogeneous boundary-value problem

\[
(\Delta+KS)u=0.
\]

It is homogeneous rather than a Poisson equation with an additive source: density acts as a spatially varying coefficient multiplying $u$. A boundary condition must fix the remaining normalization, which corresponds to the additive zero of $\phi$.

The dimensions are $[\phi]=\mathrm{m^2\,s^{-2}}$, $[S]=\mathrm{s^{-2}}$, $[K]=\mathrm{s^2\,m^{-2}}=[\phi]^{-1}$, and $[D]=\mathrm{m^2\,s^{-1}}$. Thus $K\phi$ is dimensionless, $KS$ has units $\mathrm{m^{-2}}$, and both terms in the transformed bracket have units $\mathrm{m^{-2}}$. In the paper's NGC 3198 table, $K=3.40\times10^{-5}\,\mathrm{s^2\,km^{-2}}$, so its nominal flat-regime speed $K^{-1/2}$ is $171.5\,\mathrm{km\,s^{-1}}$.

The $K\to0$ limit is regular in the original variable, not in a formula that divides by $K$. Expanding $u=1-K\phi+O(K^2)$ in the stationary equation gives $\Delta\phi=S+O(K)$, hence $\Delta\phi=S$ at $K=0$. The time-dependent original equation becomes $\partial_t\phi=D(\Delta\phi-S)$. Numerically, near zero $K$, solve in $\phi$ or use $(u-1)/K\approx-\phi$ with a stable exponential evaluation; do not divide noisy $u-1$ by $K$ without care.

## Boundaries, positivity, and stationarity

Boundary data must be transformed with the field. Dirichlet data $\phi|_{\partial\Omega}=h$ become $u|_{\partial\Omega}=e^{-Kh}>0$. If the normal derivative is given, $\partial_n\phi=h$, then $\partial_nu=-Khu$; it is not generally a fixed Neumann value for $u$. A shift $\phi\mapsto\phi+c$ rescales $u$ by $e^{-Kc}$, so an isolated gauge $\phi(\infty)=0$ means $u(\infty)=1$.

For a compact spherical source, the exterior stationary solution obeys $\Delta u=0$, hence $u=1+B/r$ for the isolated gauge. Matching at a finite radius $R$ gives $u'(R)=-(u(R)-1)/R$, an exterior-matching Robin condition. For attractive gravity on the regular branch $B>0$,

\[
g=-\phi'=\frac{u'}{Ku}=-\frac{B}{Kr(r+B)},
\qquad v_c^2=r\phi'=-rg=\frac{B}{K(r+B)}.
\]

This gives a flat interval when (r\ll B), with (v_c^2\simeq1/K), and a Kepler-like decline for (r\gg B). The observed velocity at a finite endpoint does not supply the isolated boundary condition; replacing the matching condition with that velocity makes the result conditional on the data being modeled.

The exponential guarantees (u>0) when it is constructed from a finite real \(\phi\). An unconstrained linear solve does not guarantee a positive solution. For (K>0,S>0), the stationary operator is Helmholtz-like; it can become poorly conditioned near a zero eigenvalue, and solutions can cross zero. For a uniform-density sphere truncated at (R), regular-center and isolated exterior matching give \(u_{\rm in}=\sec(x)\,\sin(qr)/(qr)\), with \(q=\sqrt{KS_0}\), \(x=qR\), and \(B/R=\tan(x)/x-1\). The positive Newtonian-connected branch requires \(0<x<\pi/2\); the edge (x\to\pi/2) is near-singular and should not be treated as a robust fit regime. More generally, stationarity is not a stability proof: at fixed source the time equation is linear, but its modes have growth rates set by (D(\Delta+KS)), which need not all be negative under the chosen boundaries.

These statements assume fixed constant (K). If (S) is prescribed and static, the transformed stationary problem is linear. In the full material problem, \(\rho\) evolves through matter/hydrodynamic equations and may depend on \(\phi\); the joint system is then coupled and is not made linear by transforming only \(\phi\). The paper invokes a short gravitational relaxation time to justify its quasistatic treatment, but that does not establish stability of every stationary branch or solve the coupled evolution.

## Paper's constant-density sign check

The paper's spherical field equation is Eq. 8, p. 2:

\[
g'+\frac{2g}{r}+S+Kg^2=0,\qquad g=-\phi'.
\]

For (K>0), constant (S>0), the transform yields \(u''+2u'/r+KSu=0\). The regular solution is \(u\propto\sin(qr)/(qr)\), so

\[
g(r)=-\frac{1}{Kr}+\frac{q}{K}\cot(qr),\qquad q=\sqrt{KS}.
\]

Near the center this is \(g=-Sr/3+O(r^3)\). The positive-(K) coth profile printed as Eq. 23 instead has the audit's nonzero residual; the cot profile satisfies Eq. 8 to numerical precision. This sign issue is already acknowledged in the supplied paper: p. 4 says its nonlinear numerical and theoretical solutions can only be compared using the opposite sign of (K). It should be reported as a printed-sign convention issue, not as an undisclosed discovery or a refutation of the broader theory.

## Independent synthetic calculation

`independent_check.py` uses the fixed synthetic source \(S(r)=e^{-r^2}\), dimensionless (K=0.4), and domain \([10^{-4},5]\). It obtains (u) from `solve_bvp` for the transformed linear equation with regular-center series data and (u(5)=1), then obtains (g) independently from `solve_ivp` applied to the original nonlinear Riccati equation. It does not import the project's audit. The collocation solve succeeds; its maximum interval RMS residual is \(1.98\times10^{-10}\), sampled (u\ge1), \(\phi(5)=0\), the maximum transformed stationary residual is \(3.47\times10^{-18}\), and the two (g) solutions agree to the errors reported above. This independently checks the printed-sign transform and positivity for this finite-domain example; it is not a proof for arbitrary sources or boundaries.

## NGC 3198: conditional reproduction versus prediction

The supplied SPARC archive contains `NGC3198_rotmod.dat` with radii, observed velocities/errors, and stellar-disk and HI rotation contributions; the photometry archive includes `NGC3198.sfb`. That is enough to attempt the paper's *conditional spherical-surrogate* recipe. The paper (Sec. 3.1, Eqs. 46–49, pp. 5–6) converts disk-plane component curves to an equivalent pseudo-spherical density through a radial derivative and a fitted stellar mass-to-light ratio. It explicitly acknowledges the spherical simplification. The SPARC curves and surface photometry do not uniquely determine the stellar and gas three-dimensional density, vertical structure, or a deprojection prescription; the gas force contribution can also encode nonspherical geometry. So they do not define a unique physical (S(\mathbf{x})).

More decisively for prediction, the paper's numerical domain ends at the outermost observed point (r_d), and its boundary derivatives at (r_{\min}) and (r_d) are set from the observed endpoint velocities (Eqs. 41–42, pp. 3–4). It fits (K) and \(\Upsilon_1\) using the observed curve; Table 1 reports \(K=(3.40\pm0.04)\times10^{-5}\,\mathrm{s^2\,km^{-2}}\), \(\Upsilon_1=0.762\pm0.007\), and \(\chi_\nu^2=1.29\). Those choices are adequate to describe a conditional fit, but they use held-in observed velocities to close the finite-domain problem and provide no isolated exterior matching. A prediction for an excluded galaxy therefore cannot use this boundary prescription. The exact interpolation/smoothing used before differentiating the tabulated curves is also not fully specified in the paper, so a bitwise reproduction is not assured from the article alone.

**Decision:** stop before fitting NGC 3198. A physical prediction needs either a resolved/deprojected three-dimensional source with stated geometry and finite-domain exterior matching, or an explicit limitation to a spherical surrogate together with a boundary rule fixed independently of the held-out velocities. The current files are insufficient for the former; the paper's endpoint rule only supports a conditional reproduction of its in-sample method.

## (K\) scaling and a falsifiable holdout

In the exterior spherical solution written as \(g=-1/(Kr+Cr^2)\), the corresponding circular speed is \(v_c^2=1/(K+Cr)\). If a measured radial interval is truly source-free and (Kr\gg Cr^2) in the denominator, then \(v_{\rm flat}^2\simeq1/K\). Combining this with a baryonic Tully–Fisher law \(v_{\rm flat}^4=A M_b\) gives

\[
K\simeq (A M_b)^{-1/2}\quad\Longrightarrow\quad K\propto M_b^{-1/2}.
\]

That exponent relies on a flat interval governed by the (1/K) term, the same BTFR normalization across galaxies, a consistent velocity and baryonic-mass calibration, and no unmodeled environmental or geometry dependence. It also assumes that a per-galaxy constant (K) is meaningful. If a different radial branch, mass-to-light freedom, or endpoint condition sets (v_{\rm flat}), the scaling does not follow. Fitting a distinct (K) for every curve and then noticing the inverse-square-root trend is not a predictive test.

A falsifiable next test is a whole-galaxy holdout. First define the quality cuts, fixed baryon/geometry prescription, and uncertainty model. On training galaxies only, calibrate the normalization in (K=A_K(M_b/M_0)^{-1/2}) (with exponent fixed at \(-1/2\)); propagate distance, inclination, stellar mass-to-light, and gas-mass uncertainty under predeclared priors. For each held-out galaxy, infer (M_b) from photometry and gas data and predict (K) without using its rotation curve. Solve the stationary problem with regular central conditions and an exterior-matched isolated boundary, then predict the full held-out curve through (v_c^2=r\partial_r\phi). Do not use held-out endpoint velocities as boundary data. Score the frozen predictions against held-out curves and compare matched nuisance freedom with Newtonian+halo and MOND models. Systematic mass-dependent residuals, excess between-galaxy scatter, or failure of predictive scores under reasonable source/geometry uncertainty would falsify the short (K(M_b)) law. Universal (K) is a separate, stricter hypothesis and should be tested separately.

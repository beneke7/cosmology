# Independent audit: isolated spherical TG boundary problem

## Result

For the *printed* Pszota–Ván equation, the Warden's stationary transformed equation and isolated-exterior Robin condition are correct, provided the source is compactly supported at radius \(R\), \(K\) is constant, and the gauge is \(\phi(\infty)=0\). For a uniform positive-density sphere and \(K>0\), the positive, finite, Newtonian-connected solution exists only for

\[
x\equiv qR=R\sqrt{4\pi G K\rho}<\frac{\pi}{2}.
\]

The equality is both the first singularity of the *normalized stationary boundary-value problem* and the zero-growth threshold for radial perturbations. It is not a blow-up of the local Riccati acceleration: that field remains finite at \(R\) as the threshold is approached. Proceed with the Warden's NGC 3198 Test 1 as a bounded **spherical-surrogate viability test**, not as a physical galaxy prediction or TG detection. A nonuniform source requires computing its principal radial growth eigenvalue and checking \(u>0\); the uniform-density \(\pi/2\) formula is a check, not an average-density shortcut.

## Equation, units, and sign

The supplied paper's Eq. (1), p. 2, with \(D=\ell^2/\tau>0\) and \(S=4\pi G\rho\), is

\[
\phi_t=D(\Delta\phi-K|\nabla\phi|^2-S).
\]

For constant \(K\), setting \(u=e^{-K\phi}\) gives \(u_t=D(\Delta u+KSu)\). Thus stationarity is \((\Delta+KS)u=0\). With \(g=-\phi'\), the original spherical stationary equation (paper Eq. 8, p. 4) is

\[
g'+\frac{2g}{r}+S+Kg^2=0,\qquad \frac{u'}u=Kg.
\]

Signs check: the plus \(KS\) in the transformed equation follows from the paper's minus source in the \(\phi\) equation; a positive central density gives \(g=-Sr/3+O(r^3)\), an attractive acceleration. This matches the Warden's displayed BVP. It does **not** make the paper's positive-\(K\) Eq. (23) coth profile a solution of the printed Eq. (8): the independent project audit already found a nonzero residual for that branch, while the cot branch satisfies the printed convention. The paper itself notes (p. 9) that its relaxation and analytic curves require opposite signs of \(K\) for comparison. Keep that caveat explicit; do not splice the coth branch into this BVP.

Dimensions in SI are \([\phi]=\mathrm{m^2\,s^{-2}}\), \([D]=\mathrm{m^2\,s^{-1}}\), \([K]=\mathrm{s^2\,m^{-2}}\), \([S]=\mathrm{s^{-2}}\), and \([q]=\mathrm{m^{-1}}\), where \(q^2=KS\). Hence \(x=qR\) and \(K\phi\) are dimensionless; \(K g^2\), \(S\), and \(g'\) all have units \(\mathrm{s^{-2}}\). The paper's fitted \(K=3.40\times10^{-5}\,\mathrm{s^2\,km^{-2}}\) is \(3.40\times10^{-11}\,\mathrm{s^2\,m^{-2}}\); do not combine it with SI density/radius without conversion.

## Uniform sphere: matching, positivity, and Riccati cross-check

For \(S=S_0>0\) on \(0\le r<R\), \(S=0\) outside, and \(K>0\), let \(q=\sqrt{KS_0}\). Regularity at the centre gives \(u_{\rm in}=A\sin(qr)/(qr)\). The isolated exterior solution with \(u(\infty)=1\) is \(u_{\rm out}=1+B/r\). Continuity of \(u,u'\) at \(R\) is equivalently

\[
R u'(R)+u(R)=1\quad\Longleftrightarrow\quad u'(R)=-\frac{u(R)-1}{R}.
\]

Since \(f(x)+xf'(x)=\cos x\) for \(f(x)=\sin x/x\), matching fixes

\[
A=\sec x,\qquad \frac BR=\frac{\tan x}{x}-1,
\quad u_{\rm in}(r)=\sec x\,\frac{\sin(qr)}{qr}.
\]

For \(0<x<\pi/2\), \(A>0\), the interior has no node, \(B>0\), and \(u>0\) everywhere. At \(x\uparrow\pi/2\), \(A\) and \(B\) diverge: the affine isolated-gauge boundary condition is incompatible with any finite amplitude at equality. Immediately above it, the unique matched interior has \(A<0\) and \(u(0)<0\), so it is not an admissible exponential transform. Later oscillatory branches have interior nodes and are also excluded. Equivalently,

\[
\rho_{\rm crit}=\frac{\pi}{16GKR^2}\quad(K>0),
\]

with units of mass density. The weak-coupling limit is regular: \(B\simeq R x^2/3=KG M_b\), where \(M_b=4\pi\rho R^3/3\), and \(g\to-GM_b/r^2\) outside.

The independent Riccati solution is

\[
g(r)=\frac{q\cot(qr)-1/r}{K},
\]

which obeys \(g'+2g/r+S_0+Kg^2=0\) and \(g=-S_0r/3+O(r^3)\). It stays finite at \(x=\pi/2\), with \(g(R)=-1/(KR)\). Thus the threshold is a global normalization/isolated-matching singularity (and stability boundary), not a local Riccati pole. A solver monitoring only \(g\) or \(v_c\) can miss the diverging \(u\), potential offset, and apparent exterior mass \(B/(KG)\).

For \(K<0\), putting \(y=R\sqrt{|K|S_0}\), the regular matched branch is hyperbolic:

\[
u_{\rm in}=\operatorname{sech}(y)\frac{\sinh(yr/R)}{yr/R},\qquad
\frac BR=\frac{\tanh y}{y}-1.
\]

It remains positive for every finite \(y\); there is no oscillatory \(\pi/2\) threshold. At exactly \(K=0\), the transform degenerates and one should use the original Poisson equation, not divide an expression by \(K\).

## Stationarity versus relaxation

Under the **same instantaneous exterior-matched finite-domain evolution**, perturbations \(w=\delta u\) obey \(w_t=D(\Delta+q^2)w\) inside and \(w'(R)=-w(R)/R\). A regular radial eigenfunction \(w\propto\sin(kr)/(kr)\) satisfies \(\cos(kR)=0\); hence the leading eigenvalue is

\[
\lambda_0=\frac{D}{R^2}\left[x^2-\left(\frac\pi2\right)^2\right].
\]

It is negative on the positive stationary branch, zero at its singular endpoint, and positive beyond. This checks the Warden's proposed radial stability threshold. The rate has units \(\mathrm{s^{-1}}\), and \(D\) changes the rate, not the threshold.

**Boundary caveat:** this discrete \(\lambda_0\) assumes the exterior is forced to adjust quasi-statically, so the Robin condition is imposed at every time. The paper's field equation also evolves the vacuum exterior. For the actual all-space diffusion problem, perturbations outside obey \(w_t=D\Delta w\); the spectrum has a continuum up to zero, not this isolated finite-domain decay rate. Matching a zero-growth exterior \(1/r\) mode still gives \(\cos x=0\), so the first radial instability threshold remains \(x=\pi/2\), but the subcritical long-time relaxation is not a single exponential mode. Report which evolution/boundary realization is used; stationarity alone is never a relaxation proof.

## Audit of the Warden claims and NGC 3198 recommendation

- **Correct as stated, with assumptions:** the transformed stationary BVP and Robin condition are exactly the compact-source, \(\phi(\infty)=0\) matching problem. Requiring both \(u>0\) and convergence under source taper/domain changes is appropriate.
- **Correct for a uniform sphere:** the claimed first threshold \(\sqrt{4\pi G K\rho}\,R=\pi/2\). For real SPARC-derived nonuniform \(\rho(r)\), use \(\lambda_{\max}[D(\Delta+KS)]\) with homogeneous perturbation matching plus direct positivity; do not substitute a mean density into this formula.
- **Important limitation:** Robin matching is exact for a stationary compact source, but applying it to every relaxation timestep is a quasi-static-exterior approximation. If testing the paper's dissipative evolution itself, evolve the vacuum exterior or clearly label the Robin eigenvalue as a bounded-domain proxy.
- **Source physicalness:** SPARC component rotation curves do not uniquely specify a 3-D spherical density. The paper's Eq. (47) uses a pseudo-spherical inversion and explicitly acknowledges a spherical simplification; retain signed gas-force contributions and never clip a negative inferred density. A failure there indicts that source prescription, not by itself TG.
- **Independent mass-scale check:** paper Eq. (17) gives crossover radius \(R_K=KGM_{aa}\). The text's Eq. (66) instead inserts \(M_{aa,K=0}\). With its Table 2 masses and fitted \(K\), these give respectively 135.7 kpc (\(M_{aa}=927.8\times10^9M_\odot\)) and 33.3 kpc (\(M_{aa,K=0}=227.5\times10^9M_\odot\)). The Warden correctly flags this as an unresolved mass-definition substitution; do not treat 33.3 kpc as the fitted branch's exact crossover scale.
- **Recommendation:** proceed with Warden Test 1, limited to an explicitly labelled conditional spherical-surrogate test. Use the published \((K,\Upsilon_1)\) as a point and bounded local profile; declare interpolation/taper first; solve the isolated problem without observed endpoint velocities; report the minimum \(u\), principal growth eigenvalue (and which exterior dynamics it represents), convergence, and signed-density failures. This is worth doing because a nonpositive or unstable isolated branch would directly expose the limitation of the paper's endpoint-closed fit. Do not promote it to a physical NGC 3198 verdict absent a defensible 3-D source. This audit does not run that fit.

## Reproducibility

Machine-readable assumptions, formulas, sample values, and checks are in [`checks.json`](checks.json). The numerical spot checks used elementary double-precision evaluation: at \(x=1.3\), direct Riccati residuals at \(r/R=0.23,0.61,0.91\) were at most \(2.1\times10^{-15}\) in normalized units; matching residuals were below \(5\times10^{-16}\). Source: local supplied paper `context/user_provided/pszota_van_2024.pdf`, SHA-256 `55cbfbbb0a916ea2ce0b70523334ee301897695860fb4bdb5bfd2278ab47b62e`; equations/pages are paper Eq. (1) p. 2, Eq. (8) p. 4, Eqs. (21)–(24) pp. 4–5, Eqs. (41)–(42) p. 8, Eqs. (61)–(66) pp. 13–17.

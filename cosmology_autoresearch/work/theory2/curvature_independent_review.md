# Independent curved-LCDM BAO cross-check

**Status: independent check; exploratory finite-domain profile screen.** This is
not a posterior, Bayesian evidence, or an independent cosmological constraint.
It checks standard curved-FRW distance equations against the released DESI DR2
13-row BAO mean and its full covariance. No mocks were generated.

## Equations and units

I checked the exact published DESI DR2 Results II PDF, pp. 4–5, Eqs. (3)–(6):
[official APS PDF](https://journals.aps.org/prd/pdf/10.1103/tr6y-kpc6),
[article DOI](https://doi.org/10.1103/tr6y-kpc6). Equation (3) maps the radial
integral to the curved transverse comoving distance with the `sinh` branch for
positive curvature density; Eq. (4) is its flat limit; Eq. (5) defines
$D_H=c/H(z)$; Eq. (6) has the curvature contribution
$+\Omega_K(1+z)^2$ in $E^2$. I use the conventional
$\Omega_k=-Kc^2/(a_0^2H_0^2)$ sign: positive $\Omega_k$ is open ($K<0$) and
uses `sinh`, while negative $\Omega_k$ is closed ($K>0$) and uses `sin`.

For late-time curved LCDM, set $x=1+z$ and impose flat closure including
curvature:

$$
\Omega_\Lambda=1-\Omega_m-\Omega_k,\qquad
E^2(z)=\Omega_m x^3+\Omega_k x^2+\Omega_\Lambda
=1+\Omega_m(x^3-1)+\Omega_k(x^2-1).
$$

With $\chi(z)=\int_0^z dz'/E(z')$, define the dimensionless transverse
kernel

$$
S_k(\chi)=\begin{cases}
\sinh(\sqrt{\Omega_k}\,\chi)/\sqrt{\Omega_k},&\Omega_k>0,\\
\chi,&\Omega_k=0,\\
\sin(\sqrt{|\Omega_k|}\,\chi)/\sqrt{|\Omega_k|},&\Omega_k<0.
\end{cases}
$$

Then $D_M=(c/H_0)S_k$, $D_H=(c/H_0)/E$, and
$D_V=[zD_M^2D_H]^{1/3}$. Writing
$\alpha=c/(H_0r_d)$ gives dimensionless BAO observables

$$
\frac{D_M}{r_d}=\alpha S_k,\qquad
\frac{D_H}{r_d}=\frac{\alpha}{E},\qquad
\frac{D_V}{r_d}=\left[z\left(\frac{D_M}{r_d}\right)^2
\left(\frac{D_H}{r_d}\right)\right]^{1/3}.
$$

Here $E,\chi,S_k,\alpha$ are dimensionless; $D_M,D_H,D_V,r_d$ have length
units. Since the screen uses free $\alpha$, it cannot separately identify
$H_0$ and $r_d$; it makes no separate inference about either.

For stable near-flat evaluation, with $y=\Omega_k\chi^2$, both branches use
the common analytic continuation

$$
S_k=\chi\left(1+\frac{y}{6}+\frac{y^2}{120}
+\frac{y^3}{5040}+\frac{y^4}{362880}+O(y^5)\right),
$$

so $S_k\to\chi$ as $\Omega_k\to0$ and
$\partial S_k/\partial\Omega_k|_0=\chi^3/6$. The standalone implementation
uses this polynomial for $|y|<10^{-4}$ and otherwise the appropriate
`sinh`/`sin` expression.

For $p=\alpha q$, the full-covariance objective is

$$
\chi^2(\alpha)=y^TC^{-1}y-2\alpha q^TC^{-1}y
+\alpha^2q^TC^{-1}q.
$$

Setting its derivative to zero gives the GLS profile
$\alpha_*=(q^TC^{-1}y)/(q^TC^{-1}q)$ and
$\chi^2_* = y^TC^{-1}y-(q^TC^{-1}y)^2/(q^TC^{-1}q)$. The script clips this
only if needed to the inherited $[10^{-6},10^4]$ amplitude interval; it was
interior in all reported fits. I also attempted direct joint optimization of
$\chi^2$ over $(\alpha,\Omega_m,\Omega_k)$ as an independent check of the
profiling step; convergence status is given below.

## Independent fits

Each distance integral is evaluated row-by-row with scalar adaptive
`scipy.integrate.quad` and float64 arithmetic. The full 13-by-13 covariance is
used without diagonalization or row reordering. Mean/covariance row order is
$D_V/r_d$ at $z=0.295$; then $(D_M/r_d,D_H/r_d)$ at $z=0.510,0.706,0.934,
1.321,1.484$; and $(D_H/r_d,D_M/r_d)$ at $z=2.330$.

| Model / $\Omega_k$ profile box | $\alpha$ | $\Omega_m$ | $\Omega_k$ | $\chi^2_{\min}$ | Active bounds |
|---|---:|---:|---:|---:|---|
| Flat LCDM ($\Omega_k=0$) | 29.52463342 | 0.29746182 | 0 | 10.271041003 | none |
| Curved LCDM, $[-0.05,0.05]$ | 29.60748971 | 0.29314448 | 0.02251697 | 9.953677899 | none |
| Curved LCDM, $[-0.10,0.10]$ | 29.60748963 | 0.29314448 | 0.02251698 | 9.953677899 | none |
| Curved LCDM, $[-0.20,0.20]$ | 29.60748962 | 0.29314448 | 0.02251697 | 9.953677899 | none |

The three nested curved domains give non-increasing profile minima (equal to
displayed precision); all remain away from the shape and amplitude boundaries.
I also checked the profiled result against direct fits in $(\alpha,\Omega_m,
\Omega_k)$ (with the flat fit omitting $\Omega_k$). Three of the four
direct-optimizer calls reported convergence; their $\chi^2$ values agree
with the profile within $2.4\times10^{-14}$. The $[-0.10,0.10]$ direct call
returned the same $\chi^2$ as GLS at its initial point, but reported
success=false, ABNORMAL at zero iterations. Treat that case as an
objective-value consistency check only, not a successful independent
optimization. The independent GLS profile/grid fit for this domain remains
the reported fit. At the flat solution, the GLS and direct objective values
are $10.271041002564933$ and $10.271041002564917$. The three curved-domain
profile/direct $\chi^2$ differences are respectively $-1.60\times10^{-14}$,
$0$, and $-2.31\times10^{-14}$.

Flat-limit recovery is exact in the computed prediction vector at
$\Omega_k=0$ (maximum $|\Delta q|=0$). Perturbing to either
$\Omega_k=\pm10^{-12}$ changes it by only $2.97\times10^{-13}$ at the flat
best-fit $\Omega_m$. At the curved best fits, tighter quadrature tolerances
($10^{-13}$ versus $10^{-12}$ absolute and relative) produced no prediction
change in float64; QUADPACK's maximum reported absolute integral-error
estimate was $2.84\times10^{-14}$.

The existing `experiments/curvature_screen/result.json` was compared
read-only. Across flat and all three curved domains, independent-minus-screen
differences are at most $3.1\times10^{-13}$ in $\chi^2$,
$5.7\times10^{-8}$ in $\alpha$, $4.3\times10^{-9}$ in $\Omega_m$, and
$7.2\times10^{-9}$ in $\Omega_k$; neither set reports a bound hit. This is
consistent with small optimizer/float64 differences. Exact deltas and the
reference-file hash are recorded in the JSON result.

## Positivity and closed-branch audit

Across $\Omega_m\in[0.05,0.60]$, $\Omega_k\in[-0.20,0.20]$, and
$0\le z\le2.33$, both derivatives of $E^2$ with respect to the shape
parameters are nonnegative:
$\partial E^2/\partial\Omega_m=x^3-1\ge0$ and
$\partial E^2/\partial\Omega_k=x^2-1\ge0$. The minimum is therefore at
$(\Omega_m,\Omega_k)=(0.05,-0.20)$. There,
$E^2=1+0.05(x^3-1)-0.20(x^2-1)$ has its continuous minimum at $x=8/3$
($z=5/3$), with $E^2_{\min}=73/108\simeq0.675926>0$.
Also $\Omega_\Lambda\ge0.20$. Thus $\chi\le z/\sqrt{73/108}$ and, on the
closed edge, $\sqrt{|\Omega_k|}\chi\le1.26743<\pi$; the sine argument
remains on its positive first lobe throughout the domain.

The independent 12-by-17-by-1002 parameter/redshift grid attains minimum
$E^2=0.6759259259259258$ at $(0.05,-0.20,5/3)$, agreeing with the analytic
continuous result to rounding. This independently audits the positivity and
closed-distance-branch claims in [`curvature_contract.md`](curvature_contract.md).

## Scope and reproduction

This is a late-time curved-LCDM geometry screen with radiation omitted; matter
is represented by $\Omega_m$, dark energy is a cosmological constant, and
$r_d$ is not computed. It does not model early-time physics, massive-neutrino
transitions, perturbations, growth, external probes, or posterior volume.
Finite profile bounds are a numerical search domain, not priors or evidence.
The reported best-fit curvature and profile-score improvement must not be
interpreted as a calibrated significance.

Run from the repository root with one BLAS thread:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python work/theory2/curvature_independent_review.py
```

This writes `curvature_independent_review.json`. It records source-data hashes,
software/runtime details, optimizer diagnostics, quadrature checks, and the
read-only comparison to the existing screen. The run was CPU-only; no GPU was
used. The script SHA-256 for the final run is
`897a7ce63d7ff1deb96b9635f7c595cfec1767906954234b8fc6a4102efc6fe0`.

# Curved-LCDM extension contract for the DESI DR2 BAO screen

**Status: derived model contract; no fit performed.** This note specifies the next background-only extension of the supplied 13-row flat LCDM profile screen. It preserves the existing mean vector, full covariance, row ordering, and free BAO amplitude.

## Conventions and predictions

Use the standard curvature density parameter

$$
\Omega_k \equiv -\frac{Kc^2}{a_0^2H_0^2}=1-\Omega_{\rm tot}.
$$

Thus $\Omega_k>0$ denotes an open universe ($K<0$); $\Omega_k<0$ denotes a closed universe ($K>0$). This matches DESI DR2 Results II Eq. (3), which uses $\sinh$ for $\Omega_k>0$, and its Friedmann Eq. (6), which adds $+\Omega_k(1+z)^2$ to $E^2$. [DESI DR2 Results II, Eqs. (3)–(6)](https://doi.org/10.1103/tr6y-kpc6).

For late-time curved LCDM, with radiation omitted as in the existing screen, define $x=1+z$ and

$$
\Omega_\Lambda=1-\Omega_m-\Omega_k,\qquad
E^2(z)=\Omega_m x^3+\Omega_k x^2+\Omega_\Lambda
=1+\Omega_m(x^3-1)+\Omega_k(x^2-1).
$$

The radial dimensionless comoving coordinate and curvature map are

$$
\chi(z)=\int_0^z\frac{dz'}{E(z')},\qquad
S_k(\chi)=
\begin{cases}
\sinh(\sqrt{\Omega_k}\,\chi)/\sqrt{\Omega_k},&\Omega_k>0,\\
\chi,&\Omega_k=0,\\
\sin(\sqrt{|\Omega_k|}\,\chi)/\sqrt{|\Omega_k|},&\Omega_k<0.
\end{cases}
$$

Then $D_M=(c/H_0)S_k(\chi)$, $D_H=c/H(z)=(c/H_0)/E(z)$, and

$$
\frac{D_M}{r_d}=\alpha S_k(\chi),\qquad
\frac{D_H}{r_d}=\frac{\alpha}{E(z)},\qquad
\frac{D_V}{r_d}=\left[z\left(\frac{D_M}{r_d}\right)^2
\left(\frac{D_H}{r_d}\right)\right]^{1/3},
\quad
\alpha\equiv\frac{c}{H_0r_d}.
$$

In the volume-averaged distance, the two bracketed factors multiply:
$D_V/r_d=[z(D_M/r_d)^2(D_H/r_d)]^{1/3}$, equivalent to
$D_V=[zD_M^2D_H]^{1/3}$ divided by $r_d$. These definitions agree with
DESI's curved transverse distance and line-of-sight distance definitions;
the volume-averaged distance is the standard isotropic BAO combination.
[DESI DR2 Results II, Eqs. (3)–(6)](https://doi.org/10.1103/tr6y-kpc6).

All three BAO predictions are dimensionless. $E$, $\chi$, and $S_k$ are
dimensionless, while $D_M,D_H,D_V,r_d$ have length units. The separate
values of $H_0$ and $r_d$ are not identifiable in this screen: predictions
depend on them only through $\alpha$. Under $H_0\mapsto\lambda H_0$,
$r_d\mapsto r_d/\lambda$, $\alpha$ is unchanged. The screen therefore
profiles $\alpha$ freely and makes no separate $H_0$ or $r_d$ inference.

## Flat limit and stable evaluation

For $y=\Omega_k\chi^2$, both curvature branches have the common analytic
expansion

$$
S_k(\chi)=\chi\left(1+\frac{y}{6}+\frac{y^2}{120}
+\frac{y^3}{5040}+\frac{y^4}{362880}+O(y^5)\right).
$$

Taking $\Omega_k\to0$ gives $S_k\to\chi$, hence $D_M/r_d\to\alpha\chi$,
exactly the flat-LCDM screen. The derivative at flatness is

$$
\left.\frac{\partial S_k}{\partial\Omega_k}\right|_0=\frac{\chi^3}{6}.
$$

For implementation, evaluate the polynomial when $|y|<10^{-4}$ and the
appropriate $\sinh(\sqrt{y})/\sqrt{y}$ or
$\sin(\sqrt{-y})/\sqrt{-y}$ branch otherwise. This avoids the removable
$0/0$ form and retains the correct sign of the first curvature correction
on both sides of zero.

## Proposed bounded profile domain and positivity

Use $\Omega_m\in[0.05,0.60]$, inherited from the starter LCDM bounds, and
profile $\Omega_k\in[-0.20,0.20]$ with the existing positive-amplitude
domain $\alpha\in[10^{-6},10^4]$. The curvature interval is a finite
numerical screen, not a prior or a claim about observational constraints.
Keep the redshift scope $0<z\le2.33$ (the starter background routine
permits $z\le3$).

For $x\ge1$,

$$
\frac{\partial E^2}{\partial\Omega_m}=x^3-1\ge0,\qquad
\frac{\partial E^2}{\partial\Omega_k}=x^2-1\ge0.
$$

Therefore the minimum across the rectangular $(\Omega_m,\Omega_k)$ domain
occurs at $\Omega_m=0.05,\Omega_k=-0.20$. There,

$$
E^2_{\min}(x)=1+0.05(x^3-1)-0.20(x^2-1).
$$

On $1\le x\le3.33$, the derivative is $x(0.15x-0.40)$; the minimum occurs
at $x=8/3$ and equals $73/108\simeq0.675926>0$. Thus $E^2$ is positive
over the full continuous domain, not only at sampled optimizer points.
Also $\Omega_\Lambda=1-\Omega_m-\Omega_k\ge0.20$. A conservative bound
$\chi\le2.33/\sqrt{73/108}$ gives $\sqrt{0.20}\chi<1.27<\pi$, so the
closed-branch $S_k$ remains positive throughout the screen.

The accompanying [CPU checker](curvature_limit_check.py) samples the domain
and evaluates the near-zero series and flat-limit derivative; its complete
output is [curvature_limit_check.json](curvature_limit_check.json). The grid
audits the analytic bounds and does not fit the data. It used Python 3.12.3,
NumPy 2.5.3, float64, and no GPU. The 12-by-17-by-1002 grid gives
$\min E^2=0.6759259259259258$ at
$(\Omega_m,\Omega_k,z)=(0.05,-0.20,5/3)$, matching the analytic
$73/108$ bound. At $\chi=1.2$, the central finite-difference check gives
$\partial S_k/\partial\Omega_k|_0=0.2880000000704541$ against the analytic
value $0.288$.

Reproduce it from the project root with:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python work/theory2/curvature_limit_check.py

## Likelihood and row-order contract

Use the unchanged values in context/data/desi_dr2_mean.txt and the
unchanged full 13-by-13 covariance in context/data/desi_dr2_cov.txt. Pair
row $i$ only with covariance row/column $i$. The released order is:

1. $z=0.295$: $D_V/r_d$.
2. $z=0.510$: $D_M/r_d$, then $D_H/r_d$.
3. $z=0.706$: $D_M/r_d$, then $D_H/r_d$.
4. $z=0.934$: $D_M/r_d$, then $D_H/r_d$.
5. $z=1.321$: $D_M/r_d$, then $D_H/r_d$.
6. $z=1.484$: $D_M/r_d$, then $D_H/r_d$.
7. $z=2.330$: $D_H/r_d$, then $D_M/r_d$.

For prediction vector $\mathbf p(\alpha,\Omega_m,\Omega_k)$ and mean
$\mathbf y$, retain

$$
\chi^2=(\mathbf y-\mathbf p)^T C^{-1}(\mathbf y-\mathbf p)
$$

with the full released covariance. The curved-LCDM profile parameter
vector is $(\alpha,\Omega_m,\Omega_k)$. Record optimizer starts,
convergence, active boundaries and profile values. A finite bounded
profile search does not define a posterior or Bayesian evidence; those
require an explicit measure/prior and integration over parameter volume.
This contract covers late-time BAO distances only, with radiation omitted
and no sound-horizon calibration or perturbation calculation.

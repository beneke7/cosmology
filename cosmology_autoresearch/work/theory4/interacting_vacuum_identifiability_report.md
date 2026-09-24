# Independent identifiability audit: interacting vacuum BAO screen

## Finding

The independent implementation reproduces the archived interacting-vacuum
profile, but the 13-row BAO vector only gives a modest best-fit reduction. The
shape parameters are locally full-rank after profiling the free BAO amplitude,
yet their whitened response is strongly degenerate. For an illustrative fixed
baryon split (f_b=0.16), the selected history has positive baryon and CDM
density through (z=2.33). This is a useful background-screen result; it is
not evidence for a physical interaction or a calibrated constraint on
(\Gamma).

## Reproduction and identifiability

The independently profiled best fit is

| Model | \(\Omega_{m0}\) | \(g=\Gamma/H_0\) | \(\alpha=c/(H_0r_d)\) | \(\chi^2\) |
|---|---:|---:|---:|---:|
| Flat \(\Lambda\)CDM, \(g=0\) | 0.297462 | 0 | 29.524633 | 10.271041 |
| Interacting vacuum | 0.386971 | −0.466665 | 30.218590 | 8.753878 |

Thus \(\Delta\chi^2=1.517163\), with the same full 13×13 covariance and
bounded GLS amplitude profile. The score differs from the previously recorded
screen by only \(3.4\times10^{-12}\); the largest component difference in the
predicted BAO vector is \(2.2\times10^{-7}\). The exact nested \(g=0\) limit
reproduces the independently fitted flat baseline.

After whitening with the Cholesky factor of the complete covariance and
projecting out the free amplitude direction, the local shape Jacobian for
\((\Omega_{m0},g)\) has singular values 136.40 and 2.501, numerical rank 2
(relative threshold \(10^{-10}\)), and a small/large ratio of 0.0183. Its
linearized parameter correlation is −0.9946; the weak direction has
\(\delta g/\delta\Omega_{m0}\simeq-5.51\). So BAO can distinguish the two
shapes locally in exact arithmetic, but the response is dominated by a long
compensating direction. Step sizes of \(10^{-5},2\times10^{-5},4\times10^{-5}\)
give stable derivatives at the displayed precision.

Profile landmarks, obtained by linear interpolation of the independently
profiled curve, are \(g\in[-0.849,-0.088]\) at \(\Delta\chi^2=1\) and
\(g\in[-1.229,0.271]\) at \(\Delta\chi^2=3.84\). These are descriptive
finite-search profile crossings, **not** posterior or calibrated confidence
intervals; no Wilks calibration is asserted.

The initial theory contract sketched a narrower exploratory window
\([-0.2,0.2]\), while the completed compute contract used the source paper's
IVS scan box \([-3,3]\). A 21-point reprofile restricted to the initial
window puts its best value on the lower boundary (g=-0.2), with
\(\chi^2=9.245618\) (\(\Delta\chi^2=1.025423\) versus flat). Thus the narrower
window would conceal the later broad-box minimum and supplies no interior
measurement of the rate.

## Physical check and scope

The source PDF's Eqs. (3)–(5), pp. 2–3, gives separately conserved baryons
and the interacting continuity equations
\[
\dot\rho_x+3H(1+w_x)\rho_x=Q,\qquad
\dot\rho_c+3H\rho_c=-Q.
\]
Its Eq. (8), p. 3, specifies \(Q=\Gamma\rho_x\), with the text below it
stating that positive \(\Gamma\) transfers energy from CDM to DE. Therefore,
for \(w_x=-1\), \(\dot\rho_x=\Gamma\rho_x\). Changing variables with
\(u=\ln(1+z)\), \(du/dt=-H\), and defining
\(M=\Omega_b+\Omega_c\), \(X=\Omega_x\), and \(g=\Gamma/H_0\), gives
\[
\frac{dX}{du}=-\frac{gX}{E},\qquad
\frac{dM}{du}=3M+\frac{gX}{E},\qquad E^2=M+X.
\]
This checks the transfer sign used by the screen: its preferred negative
\(g\) is DE-to-CDM transfer under the paper's convention. All integrated
densities are normalized to today's critical density; \(u,E,g\), the
dimensionless distances, and \(\alpha\) are dimensionless. BAO does not
separate \(H_0\) and \(r_d\), so it cannot calibrate dimensional \(\Gamma\).

At the selected fit and fixed illustrative \(f_b=0.16\), both components
remain positive over the full interval: the sampled upper split allowing
positive CDM is \(f_{b,\max}=0.7845\), and the minimum CDM density fraction
is 0.3251 today. This verifies viability only for that declared split within
the reduced background. The baryon fraction is exactly absent from the
homogeneous distance likelihood at fixed total matter; \(f_b=0.16\) is not
measured or used as an empirical prior.

The profile is late-time background-only. It omits radiation, neutrinos,
early-time sound-horizon physics, perturbations, growth, and the covariant
transfer four-vector. It therefore tests neither the paper's full
CMB+BAO+supernova analysis nor perturbative stability. The source paper is
locally available as an exact PDF, but the present audit checks its
background equations only.

## Recommendation on null calibration

Do not spend the next large matched-null ensemble on this screen. The gain is
only \(\Delta\chi^2=1.52\), and the candidate sits on a strong
\(\Omega_m\)-rate degeneracy while lacking early-time calibration and
perturbation physics. A null ensemble would be appropriate if this model
survives a physically complete follow-up and the project then needs a
search-calibrated model-selection claim. At that stage, each null realization
must repeat the same model-selection procedure used on the observations.

## Reproduction record

The script is independent of the screen implementation: it reads the mean
and covariance text directly, integrates \((M,X)\) with float64 DOP853, uses
adaptive scalar quadrature for every comoving distance, profiles \(\alpha\)
with the full covariance, and profiles \(\Omega_m\) on 102 distinct \(g\)
points spanning the broad scan, local refinement and narrow-window sensitivity
meshes. CPU only; one CPU thread; no GPU.

From the project root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python scripts/run_bounded.py --seconds 1200 -- \
  .venv/bin/python work/theory4/audit_interacting_identifiability.py \
  --gamma-grid 61 --omega-grid 41 --fine-gamma-grid 31
```

The latest run completed in 62.05 s under Python 3.12.3, NumPy 2.5.3 and
SciPy 1.18.1 on Linux x86_64. The finite search domains are
\(\Omega_{m0}\in[0.05,0.60]\), \(g\in[-3,3]\), and
\(\alpha\in[10^{-6},10^4]\); they are search bounds, not priors. The
artifact records every curve point, input hashes, machine/runtime metadata,
and Jacobian calculation.

Key SHA-256 provenance:

- Paper PDF: `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`
- BAO mean: `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`
- BAO covariance: `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`
- Audit script: `6cb55d20d50baa4520970296fdb74bf94f66273c1db3cf0bf2438018194a3052`

The machine-readable result is [interacting_vacuum_identifiability.json](interacting_vacuum_identifiability.json); the reproducible code is [audit_interacting_identifiability.py](audit_interacting_identifiability.py).

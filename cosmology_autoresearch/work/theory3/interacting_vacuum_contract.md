# Interacting-vacuum background contract for a BAO-only profile screen

**Status: independently derived theory contract plus synthetic equation check;
no DESI data were fit.** This defines a limited late-time background
extension, not a full interacting-cosmology inference.

## Source equations and conventions

I read the local paper *Do DESI-DR2 BAO data imply a coupling of dark matter
and dark energy?*, context/papers/interacting_de_desi_dr2_2026.pdf
(SHA-256
0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c).
The local copy is the September 2026 preprint, identified with the accepted
APS paper, [DOI 10.1103/6kdf-vzq4](https://doi.org/10.1103/6kdf-vzq4)
([official APS accepted-paper page](https://journals.aps.org/prd/accepted/10.1103/6kdf-vzq4);
[arXiv record](https://arxiv.org/abs/2609.05660)). Equation locators below
refer to the local PDF.

The paper's flat Friedmann equation (Eq. 3, p. 2) includes radiation,
baryons, interacting CDM, DE and neutrinos. Its dark-sector continuity
equations (Eqs. 4–5, p. 2) are

$$
\dot\rho_x+3H(1+w_x)\rho_x=Q,\qquad
\dot\rho_c+3H\rho_c=-Q.
$$

Thus \(Q>0\) transfers energy from CDM to DE. Eqs. (6)–(7), p. 3, rewrite
these as separate effective-fluid continuity equations with
\(w_{x,\mathrm{eff}}=w_x-Q/(3H\rho_x)\) and
\(w_{c,\mathrm{eff}}=Q/(3H\rho_c)\). The interaction studied is
\(Q=\Gamma\rho_x\) (Eq. 8); \(\Gamma\) has units of inverse time, the same
as \(H\), and \(\Gamma>0\) is CDM \(\to\) DE. Eq. (9) then gives
\(w_{x,\mathrm{eff}}=w_x-\Gamma/(3H)\) and
\(w_{c,\mathrm{eff}}=\Gamma\rho_x/(3H\rho_c)\).

For \(\rho_t=\rho_c+\rho_x\), Eqs. (10)–(11), p. 4, express the components as

$$
\rho_x=-\frac{1}{w_x}\left(\frac{d\rho_t}{dN}+\rho_t\right),\qquad
\rho_c=\rho_t+\frac{1}{w_x}\left(\frac{d\rho_t}{dN}+\rho_t\right),
\quad N=3\ln a.
$$

These follow from total dark-sector conservation and are consistent with
Eqs. (4)–(5). For \(w_x=-1\), the more direct first-order equations below
are simpler than the paper's second-order total-density form.

## Derived \(w_x=-1\) flat-FLRW background

Use dimensionless densities \(u_i=\rho_i/\rho_{\mathrm{crit},0}\),
expansion \(E=H/H_0\), time coordinate \(n=\ln a\), and coupling
\(\gamma=\Gamma/H_0\). Then \(d/dt=H\,d/dn\), and the vacuum and CDM equations
become

$$
\frac{du_x}{dn}=\frac{\gamma}{E}u_x,\qquad
\frac{du_c}{dn}=-3u_c-\frac{\gamma}{E}u_x.
$$

Baryons are separately conserved:

$$
\frac{du_b}{dn}=-3u_b,\qquad u_b(n)=\Omega_{b0}e^{-3n}.
$$

For this deliberately late-time screen I omit radiation explicitly, giving
flat closure and the Friedmann constraint

$$
E^2=u_b+u_c+u_x,\qquad
\Omega_{b0}+\Omega_{c0}+\Omega_{x0}=1.
$$

If radiation is retained, it must be added as a declared component (and
neutrino evolution specified); it is not silently fixed from CMB assumptions.

Let \(\Omega_{m0}=\Omega_{b0}+\Omega_{c0}\) and
\(\Omega_{x0}=1-\Omega_{m0}\). The total dust density \(u_m=u_b+u_c\)
satisfies

$$
\frac{du_m}{dn}=-3u_m-\frac{\gamma}{E}u_x,\qquad
E^2=u_m+u_x,\qquad
u_m(0)=\Omega_{m0},\quad u_x(0)=1-\Omega_{m0}.
$$

This closed two-component system is all a homogeneous BAO distance prediction
needs. The vacuum also has the exact proper-time form
\(\rho_x(t)=\rho_{x0}\exp[\Gamma(t-t_0)]\); in \(n=\ln a\) its evolution is
implicit through \(E(a)\), so I recommend integrating the coupled first-order
system rather than treating this as a closed-form \(E(z)\).

The \(\Gamma=0\) limit is exactly flat LCDM:

$$
u_m=\Omega_{m0}a^{-3},\qquad u_x=1-\Omega_{m0},\qquad
E^2=\Omega_{m0}(1+z)^3+1-\Omega_{m0}.
$$

With \(\chi(z)=\int_0^z dz'/E(z')\) and
\(\alpha=c/(H_0r_d)\), a BAO-only prediction is

$$
\frac{D_M}{r_d}=\alpha\chi,\qquad
\frac{D_H}{r_d}=\frac{\alpha}{E(z)},\qquad
\frac{D_V}{r_d}=\left[z\left(\frac{D_M}{r_d}\right)^2
\left(\frac{D_H}{r_d}\right)\right]^{1/3}.
$$

All three ratios and \(\alpha\) are dimensionless. With \(c\) in km/s,
\(H_0\) in km s\(^{-1}\) Mpc\(^{-1}\), and \(r_d\) in Mpc, \(\alpha\) has no
units. BAO alone with free \(\alpha\) cannot distinguish \(H_0\) from \(r_d\),
nor recover a dimensional \(\Gamma\) from \(\gamma=\Gamma/H_0\).

## BAO-only parameterization and identifiability

For a future screen, a transparent finite parameter box is

$$
(\alpha,\Omega_{m0},\gamma),\quad
\alpha\in[10^{-6},10^4],\quad
\Omega_{m0}\in[0.05,0.60],\quad
\gamma\in[-0.20,0.20].
$$

The \(\Omega_{m0}\) and \(\alpha\) intervals inherit the existing background
screen; the compact \(\gamma\) interval is an exploratory search window (the
paper illustrates IVS histories over roughly this range), not a prior or
constraint. Nested \(\gamma\) windows and wider-window checks should be
reported. Restrict distance calculations to the released BAO range
\(z\le2.33\), retain the 13-row mean order and full released covariance, and
profile \(\alpha\) freely. No Planck \(H_0\), \(\Omega_bh^2\), \(\Omega_ch^2\),
neutrino prior or sound-horizon calibration is imported. In particular,
\(r_d\) is not recomputed from this late-time model.

The baryon/CDM split is **exactly unidentifiable from this homogeneous
background**, not merely weakly constrained. At fixed
\(\Omega_{m0}=\Omega_{b0}+\Omega_{c0}\), the equation for \(u_m\) above is
independent of \(f_b=\Omega_{b0}/\Omega_{m0}\). Consequently \(E(z)\) and all
BAO distance ratios are invariant under changing \(f_b\); the synthetic
calculation below verifies this numerically. Do not fit \(f_b\) as though BAO
measured it or report a CMB/BBN-calibrated interacting cosmology.

The split still matters for physicality when \(\gamma<0\). Defining the
comoving CDM density \(q_c=a^3u_c\), integration gives

$$
q_c(n)=\Omega_{c0}
+\gamma\int_n^0 e^{3s}\frac{u_x(s)}{E(s)}\,ds.
$$

For \(\gamma\ge0\), the integral is nonnegative on the past-directed interval
\(n\le0\), so any nonnegative present CDM density stays nonnegative backward
in time; \(u_x>0\) as long as \(u_{x0}>0\) and \(E>0\). For \(\gamma<0\), the
integral reduces \(q_c\) into the past. It must be checked over the whole
interval; \(E^2>0\) alone does not suffice.

Equivalently, \(q_m=a^3u_m=\Omega_{m0}+\gamma I(n)\), where
\(I(n)=\int_n^0e^{3s}u_x(s)/E(s)\,ds\). For a specified split,
\(q_c=q_m-f_b\Omega_{m0}\), so physical histories require
\[
0<f_b<\min\left(1,\frac{\min_{n\in[n_{\min},0]}q_m(n)}{\Omega_{m0}}\right),
\quad q_m>0,
\]
with equality the zero-CDM boundary. If the split is left fully free, the
existence of some positive split is only a feasibility condition, not a
measurement: \(\min_n q_m(n)>0\) admits sufficiently small positive \(f_b\).
Any negative-\(\gamma\) BAO profile is therefore conditional on this explicit
viability domain. A physical baryon fraction requires independent information.

## Perturbation and stability scope

The source paper's perturbation calculation is not encoded by the background
equations alone. For its non-vacuum cases, it assumes synchronous gauge, zero
anisotropic stress and no momentum transfer in the DM rest frame; it defines
a DE rest-frame sound speed \(c_{s,x}^2=(\delta p_x/\delta\rho_x)_{\rm ref}\).
The authors require \(c_{s,x}^2>0\) to avoid pressure/large-scale instabilities
and commonly set \(c_{s,x}^2=1\). They discuss early-time instabilities and
report the large-scale-stable sectors \(w_x<-1,\ \gamma<0\) for phantom DE
and \(w_x>-1,\ \gamma>0\) for quintessence. Those \(w_x\ne-1\) conditions must
not be transferred uncritically to the vacuum case.

For \(w_x=-1\), the paper treats IVS separately: its Eqs. (19)–(20), p. 4,
have no DE perturbation contribution under the adopted transfer prescription;
the CDM density/velocity perturbations remain and depend on the specified
interaction. In this special case the background model has no \(c_{s,x}^2\)
parameter. A background-only BAO profile does not test these perturbations,
growth, early-time stability, or the covariant form of the transfer
four-vector. It is not a complete physical model comparison. The paper itself
fixes \(\sum m_\nu=0.06\) eV and \(N_{\mathrm{eff}}=3.044\) in its stated
analysis setup and combines CMB, BAO and supernovae; those assumptions and
results are not carried into this BAO-only contract. Its IVS treatment
considers both signs of \(\gamma\); a background positivity cut is distinct
from the paper's perturbation prescription.

## Independent synthetic check and recommendation

The accompanying CPU-only checker integrates the three separate component
equations backward from \(a=1\) using float64 DOP853, then tests the
LCDM limit, split invariance, the integrated CDM transfer identity, and a
sampled physicality box. It ran with one BLAS thread. Results in
[interacting_vacuum_check.json](interacting_vacuum_check.json):

- At \(\gamma=0\), maximum \(|\Delta E^2|\) against analytic flat LCDM at
  \(z=(0,0.1,0.5,1,2.33)\) is \(8.9\times10^{-15}\).
- For \(\Omega_{m0}=0.30\), changing \(f_b\) over \(0.10,0.30,0.60,0.90\)
  changes sampled \(E^2(z)\) by at most \(1.5\times10^{-14}\), for either
  \(\gamma=\pm0.05\).
- The independent quadrature check of the \(a^3u_c\) integral identity has
  maximum absolute error \(2.5\times10^{-16}\); the comoving CDM density
  moves up into the past for positive \(\gamma\), and down for negative
  \(\gamma\), consistent with the sign convention.
- On a 567-case sampled box
  \(\Omega_{m0}\in[0.05,0.60]\), \(\gamma\in[-0.20,0.20]\),
  \(f_b\in[0,1]\), \(z\in[0,2.33]\), no ODE solve failed and sampled
  \(E^2\) stayed positive, but 59 split/parameter cases had negative CDM.
  Thus the full rectangle is not a physical domain for both signs. The
  analytic \(\gamma\ge0\) result and explicit negative-\(\gamma\) viability
  mask are preferable to relying on sampled \(E^2\) alone.

**Recommendation:** a DESI-BAO profile test is defensible as a clearly
labelled, late-time geometry screen if it uses the free dimensionless
\(\alpha,\gamma\), the full 13-row covariance, the exact \(\Gamma=0\) null,
and a component-positivity filter for signed coupling. Do not interpret it as
a measurement of \(f_b\), a calibrated \(H_0\) or \(r_d\), evidence for
interaction, or a stability test. For a physical interacting-vacuum
constraint, add explicit early-time sound-horizon physics and a declared
baryon/neutrino calibration, and evaluate the corresponding perturbation
prescription; those are beyond the present scope.

Reproduce the synthetic check from the repository root with:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      .venv/bin/python work/theory3/interacting_vacuum_check.py

The script SHA-256 for the recorded run is
5713cef34f98b6486cd78478862e743734f126bffef222caa40c9b5d51814ff4.

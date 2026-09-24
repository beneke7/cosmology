# DES-Dovekie-only interacting-vacuum profile contract

**Status:** pre-fit contract complete. This specifies a two-shape-parameter, full-covariance fit to DES-Dovekie alone, with one shared magnitude intercept profiled analytically. No optimizer or new comparative score was run for this contract.

## Scope and source contract

Use the verified Dovekie release pairing and order in [the data audit](../data_audit4/dovekie_frozen_shape_contract.md): the 1,820-row SNANA whitespace Hubble diagram, its 1,820 by 1,820 full STAT+SYS precision matrix, and the official \(z_{\rm HD}>0\) selection, which keeps every row. Preserve the Hubble-diagram/covariance row map; do not sort or join on CID. The packed upper triangle in the NPZ is \(P=C^{-1}\). Unpack it in release order, reflect it, and use the full matrix. Do not add MUERR or MUERR_SYS to this covariance.

The observable is \(y_i={\tt MU}_i\). At shape \((\Omega_{m0},g)\), form

\[
D(z)=\int_0^z \frac{dz'}{E(z')},\qquad
D_L(z_{\rm HD},z_{\rm HEL})=(1+z_{\rm HEL})\frac{c}{H_{0,\rm gauge}}D(z_{\rm HD}),
\]

\[
\mu_i^0=5\log_{10}\!\left(D_{L,i}/{\rm Mpc}\right)+25.
\]

Use \(z_{\rm HD}\) in the expansion-history integral and \(z_{\rm HEL}\) only in its luminosity-distance prefactor. The shared fixed gauge is \(H_{0,\rm gauge}=70\ {\rm km\,s^{-1}\,Mpc^{-1}}\), matching the release convention and current adapter.

## Late-time background and sign convention

The local interacting-DE paper is *Do DESI-DR2 BAO data imply a coupling of dark matter and dark energy?*, arXiv:2609.05660v1, local copy at context/papers/interacting_de_desi_dr2_2026.pdf. Its Eqs. (4), (5), and (8) specify

\[
\dot\rho_x=\Gamma\rho_x,\qquad
\dot\rho_c+3H\rho_c=-\Gamma\rho_x,\qquad
Q=\Gamma\rho_x,\quad w_x=-1.
\]

The paper defines \(Q>0\) as energy transfer from CDM to vacuum. Baryons are separately conserved. Let \(u=\ln(1+z)=-\ln a\), \(g=\Gamma/H_0\), \(m=(\rho_b+\rho_c)/\rho_{\rm cr,0}\), \(x=\rho_x/\rho_{\rm cr,0}\), and \(E=H/H_0\). The paper's sign convention gives \(\Gamma/H=g/E\) and the implemented low-redshift equations are

\[
\frac{dx}{du}=-\frac{g x}{E},\qquad
\frac{dm}{du}=3m+\frac{g x}{E},\qquad
E^2=m+x,\qquad
\frac{dD}{du}=\frac{e^u}{E}.
\]

Initial conditions are \(m(0)=\Omega_{m0}\), \(x(0)=1-\Omega_{m0}\), \(D(0)=0\), so the background is flat and \(E(0)=1\). In lookback-redshift \(u\), the source terms have opposite signs: positive \(g\) reduces vacuum density and increases total matter relative to the uncoupled past history. At \(g=0\), \(x=1-\Omega_{m0}\) and \(m=\Omega_{m0}(1+z)^3\), exactly recovering flat \(\Lambda\)CDM at the same \(\Omega_{m0}\).

This is the current late-time reduction: radiation, neutrinos, perturbations, and early-time calibration are omitted. It is not a reproduction of the paper's CMB+BAO+SN analysis.

## Analytic intercept profile and normalization

For each \((\Omega_{m0},g)\), define \(\delta=\mu^0-y\), \(\mathbf 1\) as the 1,820-vector of ones, and \(q=\mathbf 1^TP\mathbf 1>0\). The single common additive SN magnitude intercept \(\mathcal M\) is profiled over all real values:

\[
\widehat{\mathcal M}=-\frac{\mathbf 1^TP\delta}{q}
=\frac{\mathbf 1^TP(y-\mu^0)}{q},
\qquad
r=\delta+\widehat{\mathcal M}\mathbf 1,
\]

\[
\chi^2_{\rm prof}=r^TPr
=\delta^TP\delta-\frac{(\mathbf 1^TP\delta)^2}{q}.
\]

The normalized fixed-covariance profile likelihood is

\[
-2\log L_{\rm prof}=\chi^2_{\rm prof}+\log|C|+N\log(2\pi),\qquad N=1820.
\]

Because \(C\) is fixed, its log determinant and the Gaussian normalization do not change the shape profile. Compute \(\log|C|=-\log|P|\) from a Cholesky factorization if an absolute normalized score is recorded. Keep the likelihood semantics as profiling. The release helper instead integrates over a flat intercept and adds \(\log(q/(2\pi))\); with this fixed covariance that is a shape-independent constant, but it is not part of the requested profiled likelihood.

The normal-equation diagnostic is \(|\mathbf 1^TPr|\), which should be near floating-point tolerance. A synthetic 3-vector SPD-covariance algebra check gave zero reported discrepancy between the direct residual quadratic and its closed-form profile expression, and zero reported normal-equation residual. This check involved no cosmological model or fit.

## Identifiability and baryon/CDM physicality

Changing the physical \(H_0\) multiplies every luminosity distance by the same inverse scale and shifts every \(\mu_i^0\) by one additive magnitude. That shift is exactly absorbed by \(\mathcal M\). Therefore DES-Dovekie alone cannot identify \(H_0\) separately from the SN absolute-magnitude/calibration intercept. The fit may constrain the dimensionless expansion-shape parameter \(g\), but cannot identify the dimensional rate \(\Gamma\) without \(H_0\) calibration.

The background ODE evolves total non-DE matter \(m=b+c\), while the physical model requires baryons to be separately conserved. For a constant present-day baryon fraction \(f_b=\Omega_{b0}/\Omega_{m0}\),

\[
b(u)=f_b\Omega_{m0}e^{3u},\qquad c(u)=m(u)-b(u).
\]

On the Dovekie interval, a positive split exists only if

\[
0<f_b<f_{b,\max},\qquad
f_{b,\max}=\min\!\left(1,\min_{u\in[0,\ln(1+z_{\rm HD,max})]}
\frac{m(u)}{\Omega_{m0}e^{3u}}\right).
\]

Set the fit's physical-feasibility guard to \(0\le z\le z_{\rm HD,max}=1.14418\), the full observed Dovekie interval. At every trial shape require finite positive \(E^2\), \(m\), and \(x\) on a dense integration grid across that interval, plus a nonempty open interval \(0<f_b<f_{b,\max}\); reject trials that fail. The current integrator samples the full interval (at least 257 points and no more than 0.002 spacing in \(u\)) and checks an interior witness \(f_b=f_{b,\max}/2\). Refine this check at accepted fit/profile points so a narrow invalid region is not missed.

This is only an existence test for a physical split. There is no CMB/BBN baryon calibration in a Dovekie-only likelihood, so \(f_b\) must not be reported as fitted or constrained. In particular, \(f_b=0.16\) is an illustrative fixed-split diagnostic in prior work, not a measurement or a prior for this fit. The code's maximum supported requested redshift remains 2.33; the explicit positivity guard for this SN-only fit is the observed Dovekie range above. The previously selected BAO point was checked by its own BAO fit through \(z=2.33\).

## Search box and interpretation

Use \(\Omega_{m0}\in[0.05,0.60]\) and \(g\in[-3,3]\) as the initial finite numerical search box, matching the current BAO-screen implementation. The local paper's Table I gives a uniform IVS prior \(g\in[-3,3]\), which agrees with the \(g\) box. It does not give a direct \(\Omega_{m0}\) prior: it samples \(\omega_b=\Omega_bh^2\in[0.005,0.1]\), \(\omega_c=\Omega_ch^2\in[0.001,0.99]\), and other parameters. Thus \([0.05,0.60]\) is an inherited search domain from scripts/background_bao.py, not a paper prior and not a probability measure for this SN-only profile. The magnitude intercept is unbounded. If the selected maximum or relevant profile surface touches either shape-box face, flag it as domain-limited and extend/recheck the domain before interpreting the profile as closed.

For the future fit, record the SN-only minimum and a descriptive \(\Delta\chi^2_{\rm SN,prof}(\Omega_{m0},g)\) surface. The already selected BAO-only IVS coordinate is \((0.38697100769747017,-0.4666647299169934)\), from experiments/interacting_vacuum_screen/result.json; the nested BAO \(\Lambda\)CDM coordinate is \((0.2974618149542753,0)\). Overlay these as fixed reference coordinates and evaluate the SN-only profile at the IVS point:

\[
\Delta\chi^2_{\rm SN,BAO\ point}
=\chi^2_{\rm SN,prof}(\Omega_{m,\rm BAO},g_{\rm BAO})
-\min_{\Omega_{m0},g}\chi^2_{\rm SN,prof}.
\]

This is a displacement within the SN-only profile, not a joint BAO×SN likelihood, evidence ratio, p-value, or significance. Do not add the BAO and SN likelihoods or claim exact probe independence; the public products contain no BAO–SN cross-covariance. Do not apply Wilks calibration to profile differences or contour thresholds. The existing [frozen-prediction summary](../compute6/dovekie_frozen_shape_cv.md) remains a separate predictive diagnostic and is not a score component of this fit.

## Required fit diagnostics

The eventual fit record should include:

- The pinned Hubble-diagram, precision archive, profile source, data-audit contract, paper, and BAO-point source identities/hashes; row count, row-order fingerprint, and \(z_{\rm HD}\)/\(z_{\rm HEL}\) ranges.
- Successful positive-definite Cholesky checks for the unpacked full precision, stable evaluation of the profiled quadratic, and the fixed-covariance normalization convention. Record that no separate MUERR term was added.
- For the selected solution and reported surface points: \((\Omega_{m0},g)\), \(\widehat{\mathcal M}\), \(\chi^2_{\rm prof}\), normal-equation residual, physicality minima, \(f_{b,\max}\), and the witness fraction. Record which parameters touch search-box edges.
- All optimizer statuses and endpoint objective rechecks, multiple starting points, a coarse full-box surface/refinement check, and a clear statement that these checks do not prove globality. A one-dimensional \(g=0\) slice should recover the independently evaluated flat-\(\Lambda\)CDM distance curve.
- An independent distance calculation at the best point and BAO reference point, with a declared tolerance; solver failure and every physical-domain rejection should be counted and reported.
- The SN profile at the BAO-only IVS coordinate and its displacement from the SN-only minimum, reported as descriptive profile geometry only. No joint score or calibrated interval is authorized by this contract.

## Review status, choices, and blockers

The independent frozen-shape review is not a prerequisite for running this separate full-Dovekie profile: it audits the conditional redshift-fold calculation, which this likelihood does not use. The parent reports that the independent empirical reproduction in [the compute8 review](../compute8/dovekie_result_independent_review.md) has now passed. That closes the outstanding review of the earlier frozen-shape result; it does not substitute for the new profile's own covariance, distance, and physicality diagnostics.

There are no blockers to implementing the specified likelihood. Remaining conditions are explicit: use the observed Dovekie interval for the fit's positivity/split guard; treat the two parameter ranges as finite search domains; expand the box if the best point or relevant surface reaches its boundary; and interpret \(f_b\) only as an uncalibrated feasibility witness. A Dovekie-only background profile does not test perturbations, CMB physics, early-time behavior, or an externally calibrated baryon fraction.

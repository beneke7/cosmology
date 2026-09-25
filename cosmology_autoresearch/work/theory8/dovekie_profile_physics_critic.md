# Independent physicality review: Dovekie-only IVS profile

**Recommendation: PASS for the specified SN-only late-time profile.** No sign, normalization, or baryon/CDM-split error found. The dense-grid split calculation can be replaced by an exact endpoint expression derived below. Keep the trajectory checks for finite, real (E) and positive (m,x,E^2) separate; this review did not run the profile or inspect SN scores.

## Exchange equations and limiting case

The local paper states ρ̇x+3H(1+w_x)ρx=Q and ρ̇c+3Hρc=−Q, with positive (Q) transferring energy from CDM to DE (Eqs. (4)–(5), printed p. 2); it sets (Q=Γρ_x) (Eq. (8), printed p. 3). Thus for (w_x=-1), ρ̇x=Γρx. With (u=-\ln a=\ln(1+z)), (du/dt=-H), (g=Γ/H_0), and (E=H/H_0), division by (-H) gives

\[
\frac{dx}{du}=-\frac{g}{E}x,\qquad
\frac{dc}{du}=3c+\frac{g}{E}x.
\]

Baryons satisfy (db/du=3b), so (m=b+c) obeys (dm/du=3m+gx/E). Under the stated flat, radiation-free late-time reduction, (E^2=m+x), and (dD/du=e^u/E) follows from (D(z)=\int_0^z dz'/E(z')). For (g=0), (x=1-\Omega_{m0}) and (m=\Omega_{m0}e^{3u}), exactly flat ΛCDM. This agrees with the equations in [the pre-fit contract](../theory7/dovekie_ivs_profile_contract.md) and the RHS in [the current profile adapter](../../experiments/interacting_vacuum_screen/interacting_vacuum_profile.py:68).

## Exact condition for a positive baryon/CDM split

Let (U=\ln(1+1.14418)), (b=f_b\Omega_{m0}e^{3u}), and (c=m-b), with Ωm0>0. Define

\[
r(u)=\frac{m(u)}{\Omega_{m0}e^{3u}},\qquad
r'(u)=\frac{m'-3m}{\Omega_{m0}e^{3u}}
=\frac{g x}{E\Omega_{m0}e^{3u}}.
\]

Conditional on a regular solution with (E>0) and (x>0), (r) is monotone: increasing for (g>0), constant for (g=0), and decreasing for (g<0). Since (r(0)=1), the exact value is

\[
f_{b,\max}=\min_{u\in[0,U]}r(u)=
\begin{cases}
1,&g\ge 0,\\
m(U)/(\Omega_{m0}e^{3U}),&g<0.
\end{cases}
\]

There exists a constant (0<f_b<1) giving (b>0) and (c>0) throughout the closed interval if and only if (f_{b,\max}>0); every such split is (0<f_b<f_{b,\max}). The upper endpoint must be excluded because it makes CDM zero where the minimum is attained. Thus the contract's dense-grid minimum is a sound approximation to this quantity, but the ODE makes the minimum an endpoint calculation. Encoding the piecewise expression would remove grid-resolution dependence for the split witness.

This proof concerns the baryon/CDM decomposition only and is conditional on a regular positive-expansion trajectory. It does not replace the separate numerical checks that the ODE remains finite and that (m,x,E^2) stay positive along the whole interval. The existing adapter checks the RHS (E^2) and samples dense output at spacing no larger than 0.002 in (u) (lines 68–89); retain those trajectory checks and refine accepted points as the contract requires.

## Scope, parameters, and remaining limitations

The contract's (\Omega_{m0}) is today's *total baryon plus CDM* fraction, not a CDM-only fraction. (g=\Gamma/H_0) is dimensionless and Γ/H=g/E; SN distance shapes can at most identify this ratio under the model assumptions. Changing physical (H_0) rescales every (D_L\propto H_0^{-1}) by a common factor, which shifts all distance moduli by the same additive constant. The freely profiled intercept absorbs that shift exactly. The fixed 70 km s⁻¹ Mpc⁻¹ gauge therefore sets the intercept convention; Dovekie alone cannot determine (H_0) or dimensional Γ. The local data audit likewise identifies (H_0) as fully degenerate with (M) ([data contract](../data_audit4/dovekie_frozen_shape_contract.md:49)).

The (\Omega_{m0}\in[0.05,0.60]) interval is a finite search domain inherited from the BAO screen, not a probability prior. The (g\in[-3,3]) interval also matches the paper's IVS prior range, but an SN-only profile over that box is not a posterior unless a measure/prior is separately specified. Boundary contact must be treated as domain-limited, as the contract says.

Checking (0\le z\le1.14418) is sufficient for feasibility of the background and positive split on the redshift interval used by this SN-only distance calculation, under the stated reduced model. It is not a certificate of an externally calibrated baryon fraction: any positive witness is existential, and arbitrarily small (f_b>0) may pass. Nor does it establish positivity beyond this interval, radiation/neutrino-era consistency, perturbation stability, CMB/BBN viability, or equivalence to the paper's full CMB+BAO+SN model. The paper explicitly evolves baryons separately and includes radiation/neutrinos (printed p. 2); the profile contract transparently omits these for its late-time scope.

Prior positivity sensitivities are consistent but narrower: theory5 and theory6 tested only the supplied BAO IVS point, over 16 (H_0,f_b) combinations, through (z=10^{10}), with no CDM crossing; theory6 reports a maximum endpoint shift (3.729\times10^{-4}) relative to its massless-neutrino version ([theory5 report](../theory5/ivs_early_positivity.md:3), [theory6 report](../theory6/ivs_massive_nu_sensitivity.md:3)). These do not certify every point in the new SN profile and do not constrain (f_b).

**Must-fix items:** none before running the explicitly scoped SN-only profile. Prefer the exact (f_{b,\max}) endpoint formula for the split check. Any eventual physicality claim must remain conditional on the SN redshift interval and an existential, uncalibrated (f_b); it must not be presented as full-model viability.

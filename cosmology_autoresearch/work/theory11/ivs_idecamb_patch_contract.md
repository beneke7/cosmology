# IVS-specific IDECAMB adaptation contract

Status: source-level contract and independent algebra check complete. This is not a code patch, build, run, or reproduction of the target paper's unpublished solver. All source line numbers below refer to the exact pinned archive members named in the JSON companion.

## Finding

The target paper's IVS background is exactly the `H0 × rho_DE` interaction available in the archived IDECAMB coupled-fluid code after identifying

\[
\beta=\Gamma/H_0,\qquad Q=\Gamma\rho_v=\beta H_0\rho_v.
\]

The target text gives Eqs. (4), (5), and (8) as \(\dot\rho_v=Q\), \(\dot\rho_c+3H\rho_c=-Q\), and \(Q=\Gamma\rho_x\), with positive \(Q\) transferring energy from CDM to DE. It explicitly distinguishes its local law from an alternative proportional to \(H\rho_x\). Thus IDECAMB's `H0rde` branch matches; its `Hrde` branch does not.

The vector choice \(Q_v^\mu=Q u_c^\mu\) is compatible with the target IVS equations and is the geodesic-CDM closure that yields homogeneous vacuum density on CDM-comoving slices. The visible IVS text, however, cites Eq. (19)-(20) from Ref. [141] without spelling out this covariant vector. The nearby explicit “no momentum transfer in the DM rest frame” assumption is introduced for the \(w_x\ne-1\) fluid equations. Therefore record \(Q_v^\mu\parallel u_c^\mu\) as the closure required by this adaptation and consistent with the IVS equations, not as a verified quotation of the target's exact source implementation.

The archived IDECAMB PPF path is not the target IVS perturbation closure. When its coupled-fluid model has nonzero coupling and \(w=-1\), it still marks DE perturbations active, evolves a PPF state, adds effective-DE density/momentum to the Einstein sources, and feeds a DE perturbation term back into CDM continuity. An IVS-specific branch must explicitly suppress that PPF stress while retaining the coupled background and the standard radiation hierarchies.

## Sign, vector, and unit cross-check

Use signature \((-+++ )\), future-directed \(u_c^\mu\), and the target's split

\[
\nabla_\nu T_v^{\nu\mu}=Q_v^\mu=-\nabla_\nu T_c^{\nu\mu},
\qquad Q_v^\mu=Q u_c^\mu,
\qquad T_v^{\mu\nu}=-\rho_v g^{\mu\nu}.
\]

Then \(\nabla_\nu T_v^{\nu\mu}=-\nabla^\mu\rho_v\). Its time projection gives \(\dot\rho_v=Q\); the CDM time projection gives \(\dot\rho_c+3H\rho_c=-Q\). Projecting orthogonally to \(u_c^\mu\) gives zero CDM force, hence geodesic CDM and \(\theta_c'=-\mathcal H\theta_c\). It also gives zero spatial gradient of \(\rho_v\) on CDM-orthogonal slices, so \(\delta\rho_v=0\) in CDM-comoving synchronous gauge. In this frame \(\delta Q=\Gamma\delta\rho_v=0\).

Dividing the perturbed CDM energy equation by its interacting background explains the positive contrast source:

\[
\delta_c'=-(\theta_c+h'/2)+\frac{aQ}{\rho_c}\delta_c
              -\frac{a\,\delta Q}{\rho_c}
            =-(\theta_c+h'/2)+a\Gamma\frac{\rho_v}{\rho_c}\delta_c.
\]

Thus positive \(\Gamma\) means CDM loses background energy, while the normalized contrast equation has a positive \(+a\Gamma\rho_v/\rho_c\) term. With \(\mathcal H=aH\), the term has units of inverse conformal time; in \(n=\ln a\) it becomes \((\Gamma/H)(\rho_v/\rho_c)\delta_c=(\beta/E)(\rho_v/\rho_c)\delta_c\), where \(E=H/H_0\). \(\Gamma\) and \(H_0\) must be expressed in the same time units.

IDECAMB documents `gQ = 8*pi*G*a^3*Q` and `grhoc_t = 8*pi*G*a^2*rho_c` in `equations_ppfi.f90:122-125`; therefore

\[
\frac{gQ}{\texttt{grhoc\_t}}=\frac{aQ}{\rho_c}
 =a\Gamma\frac{\rho_v}{\rho_c}.
\]

Its `H0rde` implementation is `gQ=beta*sqrt(grhom/3)*a*grhov_t` (`:132-134`), the archive's explicit \(\beta H_0\rho_{DE}\) branch. On port to CAMB 2.0.4, rederive this coefficient in that release's units: its `State%grhocrit` is \(3H_0^2/c^2\) in Mpc\(^{-2}\) (`results.f90:462-465`), while the IDECAMB archive uses a legacy `grhom` global not defined by the pinned 2.0.4 Fortran source. Do not copy that legacy symbol or pass dimensionless \(\beta\) as a dimensional rate without conversion.

## Narrow implementation contract

1. **Background law.** In the coupled-fluid background, use physical \(w_v=-1\) and `QForm=H0rde` semantics. Integrate both \(\rho_v\) and \(\rho_c\) from today with the same sign as Eqs. (4)-(5), feed both into \(H(a)\), and use those densities at every background/thermodynamics/transfer query. IDECAMB's `Eqs_CF` at `equations_ppfi.f90:243-255` evolves the scaled densities `y(i_de4)=8*pi*G*a^4*rho_v` and `y(i_c4)=8*pi*G*a^4*rho_c`; for \(w=-1\), its derivatives reduce to \(d y_v/da=4a\,\texttt{grhov\_t}+gQ/H\) and \(d y_c/da=a\,\texttt{grhoc\_t}-gQ/H\). This matches the target background when `gQ` is the H0rde form. In current CAMB the background entry points include `equations.f90:4-24` (`dtauda`, including `grho_no_de`) and `results.f90:1393-1410`; current density initialization/closure and equality are at `results.f90:462-498`, and tabulated species outputs at `results.f90:975-1005`. A late-time `H(z)` replacement alone is incomplete.

2. **Dedicated perturbation semantics.** Add an IVS/geodesic-vacuum branch separate from the generic PPF switch and from the stock constant-\(\Lambda\) flag. The vacuum background evolves, so marking the whole model as a true constant-\(\Lambda\) background would be wrong. In the IVS branch, allocate/evolve no DE `w_ix`/PPF perturbation state; set the DE perturbative density, momentum, pressure perturbation, and anisotropic stress to zero at every call; and add no DE perturbative terms to `dgrho`, `dgq`, or `dgpi`. Retain the DE background density and pressure in the Friedmann and acceleration equations.

3. **CDM frame and source.** In synchronous CDM-comoving gauge enforce \(\theta_c=0\) for the regular solution. CAMB's stock scalar line is `clxcdot=-k*z` at `equations.f90:2441-2443`; it represents the \(-h'/2\) source under its synchronous conventions. The equivalent IDECAMB generic line is `equations_ppfi.f90:3024-3027`:

   `clxcdot = -k*z - k*vc + (gQ*clxc-gC1*clxq-gC2*clxc-gC3*vq/k)/grhoc_t`.

   For the target branch, `vc=0`, `clxq=0`, and no DE velocity exists; the result must reduce to

   `clxcdot = -k*z + (gQ/grhoc_t)*clxc`,

   i.e. Eq. (19). If a non-comoving CDM velocity variable is retained for diagnostics, preserve the geodesic equation and its `-k*vc` term, then initialize/set it to zero for the target CDM-comoving solution. Verify CAMB's `z`/`etak` convention against `h'/2` before committing a port.

4. **Keep ordinary species standard.** Continue the stock photon and massless-neutrino hierarchies and their shear: CAMB obtains neutrino quadrupole `pir` at `equations.f90:2383-2387`, photon quadrupole `pig` at `:2389-2405`, and derives tight-coupling photon shear in `:2485-2526`. The paper's “zero anisotropic stress” sentence precedes its \(w_x\ne-1\) dark-fluid system; it does not specify deleting photon/free-streaming-neutrino shear from a CMB calculation. Dark-sector zero stress in the IVS must not zero the radiation hierarchy.

5. **Outputs and gauge.** Preserve coupled \(\rho_c(a)\) in equality, background-density exports, transfers, matter weighting, and any growth observable. Report the gauge and whether growth uses CDM-only or baryon+CDM total matter. CAMB's `results.f90:998-1004` currently exports standard CDM \(a^{-3}\) scaling, so it also needs an IVS-aware value.

## Why generic IDECAMB PPF is a different model

In the exact IDECAMB archive:

- `equations_ppfi.f90:67-84` defines selectors (`H0rde=3`, `uc=1`); `:132-134` implements the matching H0rde background; `:151-159` maps this interaction's `gC=(gQ,0,0)`; `:167-175` gives the `uc` covariant choice `gD=(0,gQ)` (the `ude` choice instead gives `(gQ,0)`). For `uc`, `LogicalPertur_CF` leaves `evolveVc` false because `gD1=0`, at `:179-193`.
- Nevertheless, for coupled CPL fluid with \(w_0=-1,w_1=0,\beta\ne0\), `LogicalPertur_CF` keeps `perturDE=true`; it disables that only when coupling is also zero (`:179-193`). PPF is the default (`equations_ppfi.f90:582-588`; readme says CF perturbations use extended PPF, README `:1-5`). The initial-state allocation and initialization are at `:1198-1218` and `:177-197`.
- The PPF stress is added to total density and momentum in `derivs` at `:2897-2971`, stored in `EV%dgrho_e_ppf`/`EV%dgq_e_ppf`, and re-added to output stress at `:1971-1982`. Even its exact \(w=-1\) special branch computes `dgrhoe=-2*k2*Kf(1)*Gamma + gQ*vT/k` (`:2901-2917`); with a nonzero interaction the `gQ*vT/k` term can create nonzero effective-DE density perturbation even if the PPF state starts at zero. For \(w\ne-1\), the extended PPF source includes interaction terms at `:2934-2968`.
- The generic CDM source at `:3024-3027` includes `-gC1*clxq`; for H0rde with `uc`, `gC1=gQ`, so this term is generally present in PPF mode. It disappears only under the IVS-specific imposed `clxq=0` semantics. Merely choosing `uc`, `w=-1`, or a zero initial PPF variable is not sufficient.

These are IDECAMB's implementation terms, not evidence that the target paper used this patch, PPF, or its initial-condition/configuration choices. The pinned IDECAMB archive is an independent implementation reference; it is not the target authors' source diff.

## Initial conditions and unresolved target choices

The target paper specifies \(\sum m_\nu=0.06\,\mathrm{eV}\) and \(N_{\rm eff}=3.044\), but its visible source does not name scalar IC mode/configuration, the exact IVS covariant vector, or the code-level vacuum perturbation closure. The independent implementation should state its choices.

For the ordinary radiation-era branch with positive early comoving CDM abundance \(D=\lim_{a\to0}a^3\rho_c>0\), the reviewed theory10 asymptotic result supports the usual leading regular synchronous adiabatic seed: \(h\propto a^2\), \(\delta_\gamma=\delta_\nu=-2h/3\), \(\delta_b=\delta_c=-h/2\), and \(\theta_c=0\). The IVS source coefficient \(A=(\Gamma/H)(\rho_v/\rho_c)\sim a^5\), so it does not alter this leading seed; on a fixed leading radiation metric its first correction is \(-A_5h/14=O(a^7)\). This check does not supply every higher-order species IC coefficient for a full Einstein-Boltzmann solve.

There are two explicit ambiguities. First, for nonzero \(\Gamma\), setting \(\delta\rho_v=0\) is not the strict all-species equal-curvature time-shift condition: \(d\ln\rho_v/d\ln a=\Gamma/H=O(a^2)\), for which that condition would give \(\delta_v=(\Gamma/H)h/6=O(a^4)\). Treat the paper's vacuum-homogeneity closure as an explicit exception and do not add an independent vacuum adiabatic IC. Second, a tuned \(D=0,\Gamma<0\) edge has \(A\to-5\); Eq. (19)'s regular solution is \(\delta_c=-h/7\), whereas equal \(\zeta_c=\zeta_\gamma\) would require \(+h/3\). This is an edge-branch incompatibility in the displayed equation/IC criterion, not proof the authors sampled that branch or used an inconsistent full covariant solver. The background must be checked globally for \(\rho_c>0\); at \(D=0,\Gamma=0\), the CDM contrast is undefined.

The recommended starting mode is CAMB's ordinary `initial_adiabatic` for the standard species, with CDM comoving, plus the explicit no-vacuum-perturbation branch. Verify the regular series and residual synchronous gauge mode removal in the pinned CAMB initial routine (`equations.f90:61-64, 1891-2053`). Do not claim the article's exact IC choice was identified.

## Regression contract for a later implementation

These are required checks, not checks run in this source audit:

1. **Background sign/unit:** evaluate both signs of \(\beta\); verify the ODE conservation residuals, \(\rho_v(a)>0\), \(\rho_c(a)>0\), and \(H^2>0\) over the entire integration interval. Check independent quadrature identities and confirm positive \(\beta\) transfers energy CDM-to-vacuum. Compare `H0rde` to the target formula and show it differs from `Hrde` away from today.
2. **Zero-coupling recovery:** at \(\beta=0\), recover stock CAMB background, density outputs, transfers, and spectra to declared precision; no branch-specific state should remain.
3. **Perturbation source:** on a fixed metric, compare the code source directly to Eq. (19) for \(\beta=0\) and both signs of \(\beta\). Check `gQ/grhoc_t = a*Gamma*rho_v/rho_c`, no `clxq`, no PPF stress, and preservation of \(\theta_c=0\). Confirm the total dark-sector perturbative conservation residual for the chosen covariant split.
4. **Radiation hierarchy:** compare photon/neutrino density, velocity and shear evolution to the standard CAMB branch in the \(\beta\to0\) limit; ensure nonzero free-streaming-neutrino shear is not suppressed by IVS dispatch.
5. **Convergence/observables:** vary integration and accuracy controls, inspect equality and \(r_{\rm drag}\), and check transfers/growth with an explicitly stated baryon+CDM or CDM-only definition. For CMB use ordinary radiation-era species settings and test TT/TE/EE and lensing as applicable.

Key failure modes are the wrong `Hrde` law, a flipped `Q` sign, raw dimensionless \(\beta\) treated as a rate, double-counting an effective \(\delta\rho_v\) through PPF, omitting the contrast-normalization source \(+Q\delta_c/\rho_c\), applying the stock \(a^{-3}\) CDM law in any background/output path, treating dynamic vacuum as a constant background via CAMB's cosmological-constant shortcut, suppressing radiation shear, and accepting negative/zero CDM density or a nonregular early-time branch without flagging it.

## Independent algebra check and limits

I independently checked the substitutions without executing either archive: (i) \(\Gamma=\beta H_0\) makes the two background rate laws identical; (ii) the IDECAMB `H0rde` scaled source has \(gQ/grhoc_t=aQ/\rho_c\), yielding exactly the target's positive contrast coefficient; and (iii) `uc` makes CDM geodesic but its ordinary IDECAMB PPF stress and `clxq` feedback do not vanish identically for nonzero coupling. The theory10 early-IC result supplies an independent asymptotic check of the regular ordinary branch and the tuned edge ambiguity. This audit inspected pinned source text and archive members only; no source was modified, compiled, or run.

Input SHA-256 values and exact archive/member/line provenance are machine-readable in `ivs_idecamb_patch_contract.json`.

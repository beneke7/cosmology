# NGC 3198 resolved-source audit and gravity pilot

**Status: H I map closure passed; stellar source remains conditional; bounded
TG branch pilot and matched-source alternatives completed. No physical
detection.** The unit and mass checks support an axisymmetric baryon-source
pilot, not a unique three-dimensional source. The resulting TG conclusion is
limited to this declared source, boundary family and one galaxy.

## Data and conversions

- THINGS natural-weighted moment-0 FITS: 1024x1024, `JY/B*M/S`, 1.5 arcsec
  pixels. The actual file's AIPS HISTORY records a restoring beam
  11.43108x9.36252 arcsec, whereas Walter et al. (2008) Table 2 gives
  13.01x11.56 arcsec for the nominal natural-weighted cube. Applying the
  file-recorded beam and the `m/s -> km/s` factor gives 218.769 Jy km/s and
  M_HI=9.823e9 Msun by direct column integration (9.816e9 by the independent
  integrated-flux relation). This is 3.5% below the published 227 Jy km/s /
  1.017e10 Msun. Using only the table beam on this particular exported map
  gives 155.680 Jy km/s and 6.985e9 Msun, a >30% underestimate; keep it as a
  visible metadata sensitivity, not the baseline conversion.
- HALOGAS DR1 HR `coldens` FITS is already N_HI in cm^-2 according to its
  HISTORY; its `BUNIT` string is stale. Direct full-map mass is 9.968e9 Msun.
  The THINGS/HALOGAS HI masses within the declared main-disk aperture R<35 kpc
  are 9.216e9 and 9.378e9 Msun, respectively (1.7% apart). About 6% of each
  full-map mass lies outside that aperture, so companions/off-disk emission is
  not silently included in the pilot disk source.
- The S4G P5 product is the cleaned 3.6-µm old-stellar map. Its official table
  marks NGC 3198 as included, ICA iteration 2, quality flag 2 (acceptable),
  with analysis ellipse SMA=215.4 arcsec, ellipticity 0.621, PA=33.7 deg. At
  the adopted 13.8 Mpc this is 14.4 kpc along the major axis. We now restrict
  baseline stellar photometry to that published ellipse; pixels beyond it
  require an explicitly labeled outer-disk extrapolation. The P5 paper's
  Eq. (6), 9308.23 I_MJy/sr D_Mpc^2 (M/L), agrees with the zero-point/pixel-solid-
  angle implementation. For M/L=0.6, measured stellar mass inside the P5
  ellipse is 1.632e10 Msun. Strictly dropping all nonzero mask labels lowers it
  1.7%; counting all finite values inside the ellipse raises it 1.8% and is a
  mask sensitivity only.

The mask treatment follows the official P5 README plus Querejeta et al. (2015):
positive labels are inherited P4 exclusions, while negative labels flag the
recursive second ICA treatment. The P5 paper recommends interpolation across
small second-iteration masked/oversubtracted regions. Accordingly, use negative
labels in the main source and retain mask==0 as a strict sensitivity. A single
M/L=0.6 carries an approximately 0.1 dex stellar-mass uncertainty; the finite
photometric aperture and outer stellar profile remain separate uncertainties.

## Physical scope and next computation

Use the HALOGAS radial HI profile as the primary gas input and THINGS as an
independent source conversion/map sensitivity, not an extra rotation-curve
likelihood. Multiply by 1.36 for helium; then test thin HI and the HALOGAS-
motivated 85% thin (h=0.2 kpc) + 15% thick (h=3 kpc) vertical mixture. The S4G
map supplies the stellar surface profile only inside its official 14.4-kpc
photometric ellipse; for a full-curve source test the outer stellar profile
must be extrapolated and varied. H2, ionized baryons, distance/inclination,
HI opacity, and vertical profiles are unresolved source nuisances. Azimuthal
averaging is a declared geometry simplification, not a recovered 3D map.

The first K=0 force baseline is now complete in
[`newtonian_hankel_baseline.py`](newtonian_hankel_baseline.py) and
[`newtonian_k0_baseline.json`](newtonian_k0_baseline.json). Its axisymmetric
Hankel kernel reproduces an analytic Freeman exponential disk to a maximum
relative error of 2.12e-4 in V²; increasing the k integration range/resolution
changes the predicted speeds by at most 2.11e-4 fractionally. Following the
independent P5 aperture audit, only full stellar rings through 13.75 kpc are
treated as measured; the 14.25-kpc bin is retained in the CSV but excluded from
the baseline as partial coverage. The explicit stellar tail is fitted over
5–12.75 kpc and varied separately.

This is not a TG fit. With the HALOGAS source, K=0 baryons alone give the
descriptive diagonal SPARC score 39286.3 on 43 points; the SPARC Rotmod
baryonic curve itself scores 26377.9 under the same diagonal errors. Neither
is a competitive full-galaxy gravity model without a halo or other additional
component. The resolved-map stellar force differs from the SPARC Rotmod
stellar component by 18.39 km/s RMS over 27 points within complete P5 support,
while gas differs by 4.95 km/s RMS. This unresolved stellar source-definition
difference must be treated as a sensitivity, not calibrated away to improve a
TG fit. Source tail choices move speeds by up to 3.76 km/s in this baseline;
the original 14.25-kpc partial-annulus value is not used as measured support.

## Nonzero-K axisymmetric TG pilot (2026-09-25)

The bounded production scan is now complete for the frozen source. The
transformed stationary equation was solved as `(B-KD)w=D`, `u=1+K w`, on a
320x480 cylindrical cell-centered finite-volume grid (`dR=.25`, `dz=.125 kpc`,
`Rmax=80`, `Zmax=60 kpc`) with an exterior monopole Robin approximation that
does not use the observed curve. The `K=0` force is extracted directly from
`w`, avoiding subtraction/division cancellation. An unshifted finest-grid
smallest-eigenvalue attempt did not converge after >80 CPU seconds and was
stopped; the generalized shift-invert eigenproblem below converged quickly.

Across fixed stellar normalizations `Υ*=0.4766, 0.6, 0.7547, 0.762`, the
sampled positive-K diagonal SPARC score declined monotonically through the
initial `K=3.9e-5 (km/s)^-2` scan edge; even there it was poor (best sampled
`χ²diag=13693.2`, 35.8 km/s RMS at `Υ*=0.762`). This is not an interior fit
or evidence for TG.

To close the admissible branch rather than infer it from CG failure, the
positive finite-volume operator was analyzed through
`B h=Kcrit diag(D) h`, where `D=4πGρ R dR dz`. On the fine grid at `Υ*=0.762`,
`Kcrit=5.3714678825e-5 (km/s)^-2`; the independent assembly gives
`5.3713874645e-5`, a fractional difference of `1.50e-5`, and eigenpair
residuals are `7.2e-12` / `4.5e-12`. The principal eigenvector is positive over
the full grid. For positive `h,D,b`, the identity
`(Kcrit-K) hᵀ D u = hᵀ b > 0` rules out a strictly positive solution for
`K>=Kcrit` in this declared discrete BVP. This is a statement about this
source, boundary and discretization, not every TG completion.

The finite-normalization solution approaches a finite limiting curve because
`v²=-R/K dR(log u)` tends to `-R/Kcrit dR(log h_mid)`. That eigenmode limit is
not itself a finite-`u` solution at `Kcrit`; it is reported separately. At
`Υ*=0.762`, the limiting curve has descriptive `χ²diag=3585.85` and 17.20
km/s RMS over the same 43 SPARC points. A finite solve at `0.998 Kcrit` has
`χ²diag=3641.76`, 17.36 km/s RMS, `u_min=31.53`, `u_max=986.92`, and relative
linear residual `1.91e-9`. Thus the score improves toward a singular
normalization but remains a poor curve match; proximity to the pole alone is
not called instability or falsification.

Numerical/source sensitivities are explicit. The coarse `dR=.5,dz=.25` grid
has `Kcrit=5.47063e-5`, 1.85% above fine, and its edge curve differs by 1.46
km/s RMS (7.14 maximum; mostly the innermost datum). At that same coarse
spacing, enlarging the box from `(80,60)` to `(100,80) kpc` changes `Kcrit`
negligibly (`5.47063e-5` to `5.47056e-5`). At `Υ*=0.762`, replacing HALOGAS gas
with THINGS moves the edge speed by only 0.12 km/s RMS; varying the exponential
stellar-tail scale by ±25% moves it by 0.28–0.31 km/s RMS. The M/L envelope
matters more: the edge score ranges from 17758 at `Υ*=0.4766` to 3586 at
`Υ*=0.762`. These scores use pointwise SPARC errors diagonally without a
covariance and are descriptive screening statistics only.

The independent conductance-graph review reproduces the `K=3.9e-5`,
`Υ*=0.762` curve to 0.00067 km/s RMS on the same fine grid, with positive `u`
over all cells and a `3.88e-10` discrete residual; it also checks the Robin
flux. The fine-grid K=0 curve is within 0.49 km/s RMS of the Hankel baseline at
matched M/L; the 1.37 km/s maximum is at the innermost point. Root and
independent code/results are in `analyze_branch_edge.py`,
`axisymmetric_branch_edge_summary.json`, and
`../tg_axisymmetric_independent_review/`.

Interpretation is deliberately limited: this is a conditional negative for a
good NGC 3198 curve within the declared axisymmetrized S4G+HI+He source,
vertical profiles, isolated monopole exterior, and M/L range, but not a
detection and not a TG-wide falsification. The same-resolved-source NFW and
simple-MOND controls are now complete; the NFW curve is far closer under the
same diagonal point errors. Distance/inclination systematics and full radial
covariance remain absent. A predictive positive result would require a
predeclared cross-galaxy K relation and a whole-galaxy holdout.

## Reproduction and provenance

Command from repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python work/tg_axisymmetric_pilot/analyze_sources.py
```

- Source script SHA-256: `b27b6dfa7b2f59283e295a20899b8b7ed13c64f571ca8d4db2225ab85c6ccf86`
- Summary JSON SHA-256: `14ccc543dc175192c8935d54b137e10f69caa7e22acc636c989d345722e76f42`
- Profiles CSV SHA-256: `2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5`
- H I map hashes: THINGS `6d98f314fa8b03bc2c048c84bdc629eccfa3f97b1f06ca4d52ce05162386448e`; HALOGAS `9e8070c2a2586a3eb9a5933c2814b18ab6f17c71054ff4051da1db7982485d90`.
- S4G stellar map hash `db2db49b15858e3bda491cb7657e5cadeeb7be2dbce29a1bbe48c66dbfc82478`; ICA mask hash `5f2c6804ed76ea44709aeb7e6a4245334dfac2d81f6b526963b6a620b6d6fc38`.
- Key sources: [Walter et al. / THINGS](https://arxiv.org/abs/0810.2125),
  [HALOGAS DR1](https://zenodo.org/records/2552349),
  [official S4G NGC 3198 products](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/galaxies/NGC3198.html),
  [P5 README](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/docs/P5_README.html),
  [Querejeta et al. (2015)](https://arxiv.org/abs/1410.0009).

The command above reproduces the source extraction and unit/mass checks. The
separate K=0 forward baseline is reproduced by:

```sh
.venv/bin/python scripts/run_bounded.py --state /tmp/tg-no-campaign-deadline.json --seconds 120 -- env OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python work/tg_axisymmetric_pilot/newtonian_hankel_baseline.py
```

K=0 code SHA-256: `85159bc98ed8d415f087d9768315b5d2b50029e5a851c931bd187e28418fcd6f`; result JSON SHA-256: `75a915e717e49ed7bf1af09f82fb5a8d902d7db0e599d6b6aff0b0c1fae9f0f6`. Its source-profile input SHA-256 matches the current CSV above.

The TG solver, branch-edge calculation, independent eigenmode reproduction,
matched alternatives, and both thermodynamics papers are summarized with
commands, hashes, limits and ranked continuation in [`REPORT.md`](REPORT.md).
The machine-readable branch result is
[`axisymmetric_branch_edge_summary.json`](axisymmetric_branch_edge_summary.json);
the same-source fits are in
[`map_matched/REPORT.md`](../tg_alternatives_contract/map_matched/REPORT.md).

# Independent S4G P5 photometry audit: NGC 3198

**Finding — independently checked.** The current source construction uses the correct P5 mask-sign policy and now clips the stellar source to the published NGC 3198 P5 analysis ellipse. The supported default radial profile domain is 0–14.0 kpc in 0.5 kpc bins (last full-ring midpoint 13.75 kpc). The ellipse-integrated P5 old-stellar mass is (1.6318\times10^{10}\,M_\odot) for the adopted (D=13.8\) Mpc and \(\Upsilon_{3.6}=0.6\). This is a source-construction audit only; I ran no thermodynamic-gravity calculation.

## Official P5 evidence

IRSA's [P5 README](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/docs/P5_README.html) defines `.stellar.fits` as the cleaned old-star (s_1) light map. It says positive mask labels are inherited from P4; negative labels mark the recursive second ICA iteration. The [official P5 table](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/docs/P5_table.txt), NGC 3198 row, records:

| P5 column | NGC 3198 value | Meaning |
|---|---:|---|
| Excluded? | 0 | ICA was run |
| ICA iteration | 2 | Second ICA solution selected |
| Quality flag | 2 | Acceptable (not bad) |
| ([3.6]-[4.5]_{s1}) | (-0.122130\pm0.022555) mag | Old-star component color |
| ([3.6]-[4.5]_{s2}) | (+0.209442\pm0.021859) mag | Dust component color |
| (F(s1/total)) | 0.727999 | Old-star fraction in the analysis area |
| Analysis ellipse | (a=215.4''), (e=0.621), PA (=33.7^\circ) | P3 25.5 mag arcsec(^{-2}) region |

The paper defines P5 quality 1 as excellent, 2 as acceptable, and 3 as bad; Q=2 is usable depending on application. Its §4.6 recommends interpolation over oversubtracted pixels and over the small regions masked in the second iteration. Thus the baseline policy in the current script is defensible: exclude positive P4-mask pixels, keep negative second-iteration pixels in the released (s_1) map as the pipeline's interpolated estimates, and retain `mask == 0` as a strict sensitivity. Do not conflate the table's `Excluded?=0` with the quality flag.

Sources: [Querejeta et al. 2015, P5 paper (arXiv PDF)](https://arxiv.org/pdf/1410.0009), especially §§4.1, 4.6–4.8 and 5.1–5.2; [IRSA NGC 3198 P5 product directory](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/galaxies/NGC3198/P5/). The paper explains that the ICA solution area is limited to the P3 25.5-mag/arcsec² isophote, specifically to avoid external field stars and background galaxies. Pixels outside that documented region are not a guaranteed dust-corrected stellar source.

## Local FITS inputs and hashes

Both maps are 701 × 891, single-precision FITS images with 0.75-arcsec pixels. The stellar map header identifies IRAC channel 1 and `MJy/sr`; the files are the public P5 products [stellar](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/galaxies/NGC3198/P5/NGC3198.stellar.fits) and [ICA mask](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/galaxies/NGC3198/P5/NGC3198.ICAmask.fits).

| Local file | Bytes | SHA-256 |
|---|---:|---|
| `work/tg_axisymmetric_pilot/data/NGC3198.stellar.fits` | 2,580,480 | `db2db49b15858e3bda491cb7657e5cadeeb7be2dbce29a1bbe48c66dbfc82478` |
| `work/tg_axisymmetric_pilot/data/NGC3198.ICAmask.fits` | 2,580,480 | `5f2c6804ed76ea44709aeb7e6a4245334dfac2d81f6b526963b6a620b6d6fc38` |

Mask counts independently read from the FITS data: 596,974 zero pixels, 3,029 negative pixels (−100: 2,964; −200: 64; −300: 1), and 24,588 positive pixels. The signs are not interchangeable. At the P5 ellipse aperture, `mask <= 0` gives 0.2029171 Jy and (1.631835\times10^{10}\,M_\odot); the strict `mask == 0` case gives 0.1994436 Jy and (1.603902\times10^{10}\,M_\odot), 1.71% lower. Including all finite pixels, including positive excluded regions, gives 0.2067256 Jy and (1.662463\times10^{10}\,M_\odot), 1.88% higher inside the analysis ellipse. Positive-mask contamination becomes much more consequential in the wider cutout: in an earlier pre-ellipse-gate R<35-kpc diagnostic, mass changed from (1.8008\times10^{10}\,M_\odot) for `mask <= 0` to (2.2439\times10^{10}\,M_\odot) when all finite pixels were included. That wider-aperture result should not be used as P5-cleaned photometry.

## Conversion and independent calculation

The current script's conversion matches P5 §5.1. P5 gives (1\,\mathrm{MJy\,sr^{-1}}=704.04\,L_\odot\,\mathrm{pc^{-2}}), and for a 0.75-arcsec pixel gives (M_*/M_\odot=9308.23\,S_{3.6}(\mathrm{MJy\,sr^{-1}})\,D^2(\mathrm{Mpc})\,\Upsilon_{3.6}(M_\odot/L_\odot)). Equivalently, the source script integrates (F_\nu=I_\nu\Omega_{pix}), uses (F_{\nu,\odot}(10\,pc)=280.9\,Jy\times10^{-0.4(3.24)}), scales luminosity by (D^2), then multiplies by \(\Upsilon=0.6\). A 1-MJy/sr pixel at 13.8 Mpc corresponds to (1.063\times10^6\,M_\odot), matching both expressions.

I independently parsed the big-endian float32 primary arrays with Python's standard library, applied the P5 ellipse and mask directly, integrated the pixel solid angle, and converted flux to luminosity/mass. The ellipse sum reproduces the current summary's (0.2029170528\) Jy and (1.631834948\times10^{10}\,M_\odot) to rounding. This does not validate a target-specific P5 background correction or the M/L assumption.

## Aperture support and recommended domain

At 13.8 Mpc, (1''=0.0669043\) kpc, so the P5 ellipse semimajor axis is 14.411 kpc. Its projected semiminor axis is (0.379a=5.462\) kpc; deprojecting that semiminor with the source script's (i=72^\circ) yields 17.675 kpc. Because the P5 ellipse is less flattened than the adopted disk projection ((b/a=0.379) versus \(\cos72^\circ=0.309\)), the ellipse contains complete deprojected circular annuli only through its 14.411-kpc major-axis radius. For the script's 0.5-kpc bins, all points in bins ending at 14.0 kpc are supported; the 14.0–14.5-kpc bin is clipped by the ellipse. The current coverage threshold reports its 14.25-kpc midpoint at 0.920035 coverage, so it is a partial annulus, not a fully observed radial ring. This is a radial-support bookkeeping correction, not a photometric recalculation. Use 13.75 kpc as the largest default radial sample and report the 14.25-kpc point only as a partial-coverage sensitivity.

Coverage-value reconciliation: the earlier draft's 0.832 figure was not tied to a hash-stamped output and cannot be supported as a separate snapshot; treat it as an error. The current [source-profile CSV](../tg_axisymmetric_pilot/source_profiles.csv) has SHA-256 `2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5`; the current [analyzer](../tg_axisymmetric_pilot/analyze_sources.py) has SHA-256 `b27b6dfa7b2f59283e295a20899b8b7ed13c64f571ca8d4db2225ab85c6ccf86`. The CSV gives 0.9200353027941259 at 14.25 kpc and 0.9791439102161514 at 13.75 kpc. I independently recomputed the analyzer's center-count/area formula from the same current FITS inputs and geometry: the 14.0–14.5-kpc annulus contains 5,504 map pixel centers, 5,340 inside the P5 ellipse, and 5,055 finite pixels with `mask <= 0`; these yield 0.920035302794 coverage. Thus the 14.25-kpc bin is mostly supported but still partial because its outer edge crosses the finite P5 ellipse and it also excludes masked pixels.

For integrated stellar mass, use the exact P5 ellipse, not a deprojected circular aperture. In the current script the stellar source mask is clipped to that ellipse before the R apertures are summed; therefore the `R_lt_20`, `R_lt_25`, `R_lt_30`, and `R_lt_35` stellar masses all plateau at the same full-ellipse sum. Those labels should be read as “within the ellipse and R aperture,” not as independent total masses to 20–35 kpc. The source summary's R<15-kpc value is lower because it omits the low-surface-brightness part of the ellipse extending to deprojected radii near 17.7 kpc.

## S4G–SPARC stellar-component RMS is a source-definition systematic

The reported 18.181 km/s is a model-to-model component RMS, not an observed-rotation-curve residual and not a TG result. In the [Newtonian baseline record](../tg_axisymmetric_pilot/newtonian_k0_baseline.json), it is the unweighted RMS over 28 SPARC radii through 14.05 kpc of the resolved-map stellar speed minus `sqrt(0.6) * Vdisk` from `NGC3198_rotmod.dat`; the calculation is shown in the [baseline script](../tg_axisymmetric_pilot/newtonian_hankel_baseline.py). The `sqrt(0.6)` rescaling is the correct velocity scaling for SPARC's tabulated `Vdisk`, which is normalized to M/L=1. The result is descriptive; there are no pointwise uncertainties or covariance for this cross-model component difference.

The two curves do not use identical stellar-light definitions. P5's `s1` map is the ICA-separated old-stellar component, and the P5 table reports `s1/total = 0.727999` in its analysis region. By contrast, SPARC derives its disk model from the observed 3.6-μm surface-brightness profile, extrapolating the fitted disk profile unless a clear truncation is present; its paper treats all stars as disk light for NGC 3198's Hubble type T=5. SPARC reports total `L[3.6]=38.279×10^9 Lsun`; the independently integrated P5 `s1` luminosity is `27.197×10^9 Lsun` within its ellipse, or 0.710 of that catalog total. This near-match to the P5 old-light fraction is suggestive, but not an aperture-matched equality: one is a finite P5 ellipse and the other is SPARC's total-light estimate ([SPARC master paper](https://astroweb.case.edu/ssm/papers/AJv152n157.pdf); local [SPARC Table 1](../data6/SPARC_Lelli2016c.mrt)).

As a direct radial check, I linearly interpolated the P5 face-on `s1` mass profile in [the source-profile CSV](../tg_axisymmetric_pilot/source_profiles.csv) to the 28 Rotmod radii and compared it with `0.6 × SBdisk`. The P5-to-SPARC surface-density ratio has median 0.693 and mean 0.736 (range 0.542–1.125); the central part is particularly lower in P5. This supports treating the 18.18 km/s RMS first as a mismatch between an old-star-only ICA source and a total-3.6-μm disk model, plus differences in measured profile/extrapolation, rather than as evidence about gravity. Other non-identical source choices remain: SPARC adopts an exponential vertical scale height `z_d=0.196 R_d^0.633` (about 0.60 kpc for its NGC 3198 `R_d=5.84` kpc), while the resolved-map forward model uses 0.30 kpc; adopted inclinations are 73° versus 72°. Both models use 13.8 Mpc. These should be harmonized or explicitly varied before interpreting the component RMS physically.

For a fair closure comparison, use either matched old-star-only photometry on both sides or a total-light profile and M/L prescription matched to SPARC, with identical distance, inclination, vertical profile, and outer continuation. Until then, keep 18.18 km/s as a source-definition/systematic diagnostic only—not evidence for or against thermodynamic gravity.

## Background and foreground limits

The IRSA [P1–P3 README](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/docs/pipelines_readme.html) says the original P1 mosaics are in MJy/sr at 0.75 arcsec/pixel and are not background-subtracted; P3 estimates sky level and noise. The [S4G catalog field definitions](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/gator_docs/s4g_colDescriptions.html) distinguish the sky level (`sky1`) from local noise (`ssky1`) and large-scale sky uncertainty (`esky1`). An official [Gator query for NGC 3198](https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog=s4gcat&objstr=NGC3198&spatial=cone&radius=5&outfmt=1&selcols=object%2Csky1%2Cssky1%2Cesky1%2Csky2%2Cssky2%2Cesky2%2Csma1_25p5%2Cpa1_25p5%2Cellip1_25p5) reports `sky1=0.033598`, `ssky1=0.009622`, and `esky1=0.002753 MJy/sr` at 3.6 μm. The map values outside the P5 ellipse are not a blank-sky measurement (they include any real outer disk); the finite, unmasked values have median 0 and mean 0.00321 MJy/sr, comparable in scale to `esky1`. This is consistent with a near-zeroed P5 map plus a low-level outer signal, but it does not separate outer-disk light from residual background.

Do not subtract `sky1=0.033598` from the P5 map a second time based only on P1 metadata: the local P5 map is not simply a positive sky pedestal, and P5 uses P3 sky uncertainties in its propagated errors. For scale, if `esky1` acted coherently as an additive residual over all 96,193 valid ellipse pixels, its mass-equivalent would be about (2.82\times10^8\,M_\odot) (1.7% of the aperture sum); this is a sensitivity scale, not an asserted residual correction. P5's fixed \(\Upsilon=0.6\) carries an approximately 0.1-dex old-population systematic (about 26%), and the total mass scales as (D^2).

The S4G survey paper reports a typical 3.6-μm mosaic rms of 0.0072 MJy/sr and profiles reaching about 1 (M_\odot\,\mathrm{pc^{-2}}) at 1σ ([Sheth et al. 2010](https://arxiv.org/pdf/1010.1592)). The current source profile near 27.25 kpc is about 0.73 (M_\odot\,\mathrm{pc^{-2}}), equivalent under its adopted (i,\Upsilon) to roughly 0.0056 MJy/sr, below that survey-typical rms. This reinforces that a high geometric map-coverage number beyond the P5 ellipse is not evidence for a clean stellar tail. Positive P4 masks are especially important because the P5 method explicitly masks foreground stars, background galaxies, and artifacts before ICA; excluded objects are visible as holes in the current masked map. Their individual identities were not assessed here.

## Audit of current local source artifacts

The current [source script](../tg_axisymmetric_pilot/analyze_sources.py) uses `mask <= 0` and the P5 ellipse for the stellar source (lines 193–206), the pixel-solid-angle/solar-zero-point conversion (164–181), and 0.5-kpc radial bins with a \(\ge0.8\) coverage rule (122–145). The [current summary](../tg_axisymmetric_pilot/source_summary.json) records the P5 row, mask counts and sensitivities, aperture masses, and source hashes. The ellipse clipping fixes the main aperture-domain issue. I recommend tightening the default radial-ring cut so the 14.25-kpc partial bin is not treated as a complete azimuthal ring; alternatively preserve it but flag its 0.920035 support fraction. Keep any stellar tail outside the ellipse as an explicit extrapolation/nuisance component in later source construction.

### Sources

- [S4G Pipeline 5 README](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/docs/P5_README.html)
- [S4G Pipeline 5 table](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/docs/P5_table.txt)
- [Querejeta et al. (2015), Pipeline 5 paper](https://arxiv.org/pdf/1410.0009)
- [S4G P1–P3 product README](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/docs/pipelines_readme.html)
- [S4G catalog column definitions](https://irsa.ipac.caltech.edu/data/SPITZER/S4G/gator_docs/s4g_colDescriptions.html)
- [S4G survey paper, Sheth et al. (2010)](https://arxiv.org/pdf/1010.1592)
- [SPARC master paper, Lelli, McGaugh & Schombert (2016)](https://astroweb.case.edu/ssm/papers/AJv152n157.pdf)
- [S4G data DOI](https://doi.org/10.26131/IRSA425)

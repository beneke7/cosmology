# NGC 3198 alternative-model comparison contract

Status: literature/data contract only; no model was fit and no model preference is claimed. Prepared 2026-09-25. Scope is an eventual matched comparison against the 2024 thermodynamic-gravity (TG) NGC 3198 calculation.

## Frozen observation and data limits

Use every row of `context/data/Rotmod_LTG.zip::NGC3198_rotmod.dat`: 43 radii from 0.32 to 44.08 kpc, with `Vobs`, `errV`, `Vgas`, `Vdisk`, `Vbul`, `SBdisk`, `SBbul`. The exact local archive SHA-256 is `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`; the member has 43 non-comment rows. Its header says D=13.8 Mpc. The SPARC Table 1 record independently specifies NGC3198 D=13.80±1.40 Mpc (Cepheids), inclination i=73±3 degrees, Q=1, and 3.6 μm luminosity; source row and field definitions: [SPARC Table 1](https://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt), data release page [SPARC](https://astroweb.case.edu/SPARC/). Do not drop inner or outer points after looking at fit residuals.

SPARC supplies per-point velocity errors, not a 43×43 covariance matrix. The reported error prescription includes local rotation-curve deviations and side-to-side asymmetry, but the released file does not give covariance between radii; shared distance, inclination, geometry, and calibration systematics are also not represented by the diagonal `errV`. Thus a diagonal Gaussian chi-square is a reproducible screening likelihood only, not a complete uncertainty model or decisive model-selection statistic. State this limitation beside any score; do not manufacture off-diagonal entries.

The file encodes a negative gas contribution in the inner region. Preserve SPARC's signed convention in the Newtonian baryonic acceleration: `Vgas*abs(Vgas)` and `Vdisk^2` (and `Vbul^2` if nonzero). NGC 3198's bulge column is zero. The SPARC mass-model paper explicitly gives this signed-sum convention and says the gas component includes helium correction ([Lelli et al. 2016, arXiv:1606.09251](https://arxiv.org/abs/1606.09251), mass-model equation and data description). Do not replace negative gas entries by positive magnitudes before summing squared velocities.

For the primary matched comparison, fix D=13.8 Mpc and i=73° for every model, and use the same 43 observed velocity/error pairs. A secondary nuisance sensitivity may profile D and i jointly with shared external constraints D=13.8±1.4 Mpc and i=73±3°; every model must use the same likelihood/prior and geometrically consistent rescaling. Do not compare one model at floating D/i to another held fixed. Distance and inclination are systematic model inputs, not innocuous cosmetic shifts (the TG article itself notes this; SPARC's error discussion notes inclination uncertainty is absent from pointwise velocity errors).

## Shared stellar normalization and candidate alternatives

`Vdisk` is the SPARC 3.6 μm disk contribution at unit mass-to-light ratio. SPARC's reference stellar normalization is Υdisk=0.5 M☉/L☉ at 3.6 μm, motivated as a population-synthesis consistent value; it is a benchmark, not a per-galaxy measurement with a fully specified Gaussian error ([Lelli et al. 2016](https://arxiv.org/abs/1606.09251)). The TG article instead fits a single Υ1=M*/L parameter. To prevent different baryon assumptions driving the outcome, publish (a) a benchmark with Υdisk fixed to 0.5 for all, and (b) a nuisance-matched run where one common Υdisk is fitted/profiled in every family over the same explicitly declared physical domain. Do not impose a narrow Gaussian width and call it data-derived. Use no bulge M/L parameter for this galaxy. Gas normalization stays at SPARC's built-in convention; no independent gas multiplier in the core contract.

Two concrete alternatives are suitable for an exploratory matched comparison:

| Family | Definition | Free parameters at fixed D/i | Core contract |
|---|---|---:|---|
| Newtonian + NFW halo | Spherical NFW `rho= rho_s/[x(1+x)^2]`, x=r/r_s; halo circular speed from enclosed mass. | 2 halo (`rho_s`,`r_s`) + 1 shared stellar Υ if floated; 2 if Υ fixed. | Flat, stated bounds on log halo parameters, no hidden concentration relation. Add a concentration–mass relation only as a separate prior-informed model, reporting its source and effective freedom. NFW is a conventional primary halo family (Navarro, Frenk & White 1996, [arXiv:astro-ph/9611107](https://arxiv.org/abs/astro-ph/9611107)). |
| Isolated MOND, simple μ | Algebraic radial relation `g μ(g/a0)=gN`, `μ(x)=x/(1+x)`, fixed `a0=1.2×10^-10 m s^-2`; equivalently `g=[gN+sqrt(gN^2+4a0 gN)]/2`. Derive gN from the signed SPARC baryon terms. | 1 shared stellar Υ if floated; 0 if Υ fixed. | No EFE, no fitted a0, no interpolating-function shape parameter in core case. This is an explicit one-dimensional/isolated MOND prescription; it is not a full disk-field solution and does not share TG's pseudo-spherical source construction. MOND's scale/interpolating-function setup originates in Milgrom (1983, [ADS record](https://ui.adsabs.harvard.edu/abs/1983ApJ...270..365M/abstract)); a concrete simple analytical μ family is discussed by Famaey & Binney (2005, [arXiv:astro-ph/0506723](https://arxiv.org/abs/astro-ph/0506723)). |

If using an EFE MOND variant, name it separately and count its external-field amplitude and orientation parameters (or externally constrained field components) explicitly. The TG 2024 paper's comparison cited a MOND EFE calculation by Chae et al.; that is not the isolated one-parameter baseline above. A fit with free `a0`, free EFE, or different interpolating function is a different model and requires a separate row in the comparison.

## TG comparator and non-shared physical ingredients

Exact user-provided TG article: Pszota & Ván, *Physics of the Dark Universe* 46 (2024) 101660, DOI [10.1016/j.dark.2024.101660](https://doi.org/10.1016/j.dark.2024.101660), local PDF SHA-256 `55cbfbbb0a916ea2ce0b70523334ee301897695860fb4bdb5fd2278ab47b62e`. Local audited notes are in `context/reading_cards/pszota_van_2024.md`. Its galaxy method (article §§3.1–3.2, Eqs. 46–61) pseudo-sphericalizes the disk/gas velocity contributions to construct a density source; it fits Υ1 and K. Its reported reduced-χ² degrees of freedom count two observed endpoint velocities used as fixed derivative boundary values as well as the two fitted parameters: effectively four data-dependent degrees of freedom are consumed by its stated accounting. The endpoints are therefore not independent predictive information. TG's current published construction cannot be treated as a same-information predictive holdout unless a later implementation replaces those endpoint conditions with a physical boundary condition independent of the observed curve and evaluates all 43 points consistently.

The following inputs are not physically shared and must be identified rather than hidden in a single parameter-count number:

- NFW adds a non-baryonic halo density profile with two freely chosen scales in the core fit.
- Isolated MOND modifies the acceleration law and adds a fixed universal acceleration scale; its interpolating function is a modeling choice. EFE MOND additionally depends on an external field.
- TG solves its modified field equation using a pseudo-spherical source reconstructed from disk rotation components, plus K and Υ1, and uses observed endpoint derivative data in the article's fit. These are different source geometry and boundary information from either alternative.
- SPARC `Vgas`, `Vdisk`, and `Vbul` are Newtonian mass-model components in the disk plane, not a unique three-dimensional density field. They do not uniquely determine TG's spherical source density; pseudo-spherical reconstruction is an added modeling assumption.
- D/i can be treated as the same nuisance quantities and external priors, but their transformations through a disk mass model, a spherical halo model, MOND's baryonic field, and TG's reconstructed source are not literally identical.

## Admissibility recommendation

The data contract is admissible for an exploratory, explicitly diagonal-error comparison of NFW and isolated simple-MOND alternatives over the same 43 SPARC points, with common baryon and geometry treatment. A formal claim that any one theory predicts NGC 3198 better is not admissible against the published TG fit as-is: its endpoint boundary conditions ingest observed data, the source geometry is model-specific, and no full covariance is available. A later comparison may proceed as a clearly labeled diagonal-error sensitivity study if TG gets an independent physical outer boundary and the same observations are scored without endpoint reuse. Otherwise compare descriptive curve residuals only and make no inferential preference claim.

## Source/version and hash record

- TG journal PDF: local `context/user_provided/pszota_van_2024.pdf`; SHA-256 above; DOI above. The local exact journal PDF was used for contract facts, with equation/page audit summarized in the reading card.
- SPARC Rotmod archive: local `context/data/Rotmod_LTG.zip`; SHA-256 above; SPARC release page linked above. `NGC3198_rotmod.dat` header identifies the adopted distance; content row count verified locally.
- SPARC galaxy metadata: [Table 1 MRT](https://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt), live primary data product (record fetched 2026-09-25); its NGC3198 row gives D, method, errors, i and error, and Q.
- SPARC paper: Lelli, McGaugh & Schombert (2016), [arXiv:1606.09251](https://arxiv.org/abs/1606.09251); local searchable text `context/papers/sparc_master.txt` was used for M/L, signed baryonic sum, gas/helium, and error caveats.
- NFW: Navarro, Frenk & White (1996), [arXiv:astro-ph/9611107](https://arxiv.org/abs/astro-ph/9611107).
- MOND: Milgrom (1983), [ADS](https://ui.adsabs.harvard.edu/abs/1983ApJ...270..365M/abstract); Famaey & Binney (2005), [arXiv:astro-ph/0506723](https://arxiv.org/abs/astro-ph/0506723).


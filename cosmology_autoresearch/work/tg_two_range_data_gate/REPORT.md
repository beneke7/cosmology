# Can the linked two-range model be refit to released Eöt-Wash data?

## Finding

**Partly feasible, but a calibrated full-apparatus 95% contour is not reproducible from the currently accessible release.** John Lee’s 2020 UW dissertation publicly prints the experiment’s *processed* Appendix A.1 torque dataset: 95 runs, three in-phase harmonics per run, a separation and marginal uncertainty, and one marginal torque uncertainty per harmonic. This allows a model-specific table-level likelihood to be constructed. It is substantially better than transferring the published single-Yukawa `|alpha|=1` curve.

The released table does not give the science-run horizontal stage positions or the underlying capacitance/actuator records used to reconstruct each run’s geometry and gap. Nor does it publish a covariance matrix for torque, separation, systematic, or nuisance errors. The article’s APS supplement is marked subscription-required and was not accessible in this audit. Without those inputs, a table-only profile must assume centered geometry and the published diagonal error model; it is a conditional sensitivity calculation, not a reproducible full-apparatus 95% exclusion contour. No fit, exclusion, or detection was produced here.

## Sources and provenance

The theory input is M. Pszota and P. Ván, “Ghost-free higher-gradient Newtonian gravity from the Second Law of Thermodynamics,” arXiv:2609.00317v1 (31 August 2026), especially pp. 3–4, Eqs. 11–14. The exact local PDF is `context/papers/tg_higher_gradient_2026.pdf`, SHA-256 `c0992fc15699bcd2c7abd5a4b044ace25203c2793c0cca42cb775888945fe165`; its extracted text companion has SHA-256 `bf3eb495efa371fa3ecd078af96ffea018e2c25c076499a7b1766eab3244acfe`. [arXiv v1 record](https://arxiv.org/abs/2609.00317v1).

The primary experiment paper is J. G. Lee et al., *Physical Review Letters* 124, 101101 (2020), arXiv:2002.11761v1. The local exact v1 PDF in `work/tg_lab_bound_independent/primary/Lee_Adelberger_2020_arxiv.pdf` has SHA-256 `ffcddf6d2c1a758f07112a3125ae8583254a3089b1022ac2cf459a052652504c`; its extracted text has SHA-256 `4cd6c5d32c6762564e0a7c2392314d5c9ad6e5aeb0a3aafbcd80ad4e16244583`. The [APS article page](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.101101) supplies the publication record; its [accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevLett.124.101101) is publicly served through CHORUS. The APS page labels its supplemental material “Subscription Required”; direct official-page access returned authorization errors, and no supplement bytes are present in the supplied source set.

The analyzed torque rows are in J. G. Lee, *A Fourier-Bessel Test of the Gravitational Inverse-Square Law*, UW PhD dissertation (2020), especially §§11.1–11.4, Tables 11.1, 11.3, 11.4 and Appendices A.1–A.2, pp. 123–152. The official [UW ResearchWorks item](https://digital.lib.washington.edu/researchworks/items/971237d1-100a-41ae-9027-d1bbce8cf315) dates the dissertation 30 April 2020 and lists the PDF as its file; the [official bitstream](https://digital.lib.washington.edu/researchworks/bitstreams/6013936a-da00-42b7-ba0a-4dff0cb05bf8/download) matches the local copy `context/papers/lee_2020_eotwash_thesis.pdf`, SHA-256 `00d13a466e4f8c14c3f6067d49d90fd0c49a89e72a8cf93f3a79c18d6aef924a`. The supplied extracted text companion has SHA-256 `73d7b783f6825c68730099a30bb9f3df96020e33227885b84f7e8569f907a848`.

## What is actually released

Appendix A.1 labels a complete list of the 95 runs used in the Newtonian and Newtonian-plus-Yukawa fits. Each row gives `s [µm]`, `N120`, `N18`, and `N54`; each torque is printed as a value ± error in fN·m. The table has 95 rows, hence 285 in-phase torque entries. The tabulated separations span 46.2 ± 0.7 µm to 3036.3 ± 1.2 µm. The two 46.2 µm rows have very large torque errors. The chapter and accepted manuscript describe the experiment’s face-to-face coverage as 52 µm–3.0 mm; the difference between that headline and the table’s lowest printed `s` is not resolved here, so the table minimum should not be silently substituted for the stated physical minimum.

These are analyzed measurements, not the full raw run streams. Appendix A.1 says that the values have been rotated into one common in-phase component, scaled by the measured torque calibration, and corrected by subtracting systematic torques. It says the unrotated sine/cosine amplitudes and errors are in `.sum` files on the UW `\\TYCHO-BRAHE\SRDATA\BATDATA\` share; those files are not in the dissertation record or the supplied local source set. The PRL manuscript separately describes the raw streams as including twist angle, turntable angle, capacitance, three-dimensional displacement, and other run variables.

The release contains useful associated information:

- Table 11.1 supplies the common phase, `φ0 = 51.77146 ± 0.00077°`; the data table is already reduced to its in-phase component.
- Table 11.3 prints measured means and 1σ uncertainties for the 17 apparatus/model parameters. The final simplified fit floats `x0`, `y0`, `s0`, overcut `ε`, and common torque scale `γ`, fixing the remaining measured quantities at their measured values. The full fit first allowed the full set to vary. The likelihood uses separate Gaussian penalty terms for the listed nuisance measurements.
- Appendix A.2 gives some gravitational-centering scan runs with `x`/`y` micrometer positions and `N120`. These are dedicated centering scans, not the `xj,yj` values for each of the 95 science runs.
- The thesis documents the mass-pattern geometry and corrections, including the Fourier–Bessel torque method, dimensions/material parameters, tilt correction, and an empirical off-center correction. It does not include executable torque-generation code or the generated torque tables needed to directly evaluate a new kernel.

**Errors and correlations.** The table supplies marginal errors, not covariance matrices. Equation 11.1 explicitly uses a diagonal sum over run `j` and harmonic `m`, with denominator `δNm² + δzp,j² (∂Ñm/∂sj)²`, plus independent quadratic nuisance penalties. Thus the author’s published approximation is available as a stated likelihood contract. The release does not provide the joint covariance of the original sine/cosine fits, shared separation calibration, systematic subtractions, or apparatus priors; using the published diagonal form reproduces its assumption, not an independently validated covariance model. A shared global `γ` prior and global separation-offset nuisance `s0` are present, but the underlying calibration and run records needed to inspect their cross-correlations are not.

## Physical/statistical refit contract

The 2026 stationary potential is linear in its Green kernel, so calculate the linked correction on the actual torsion-balance geometry. Do not fit independent Yukawa amplitudes. With `a,b` the roots of `x² − ℓ1² x + ℓ2⁴ = 0`, the real-root sector has

`a,b = (ℓ1² ± √(ℓ1⁴ − 4ℓ2⁴))/2`, `λa,b = √(a,b)`,

`αa = −a/(a−b)`, `αb = b/(a−b)`, and `αa+αb=−1`.

For each run and harmonic, the real-root template is `ÑTG = ÑN + αa ÑY(λa) + αb ÑY(λb)`, where `ÑY(λ)` is the apparatus torque for unit-strength Yukawa kernel, recomputed with the same geometry corrections. The exact limits are part of the fit: when `ℓ2=0`, `λ=ℓ1` and `α=−1`; at repeated roots, `ℓ1=√2 λ`, `ℓ2=λ`, and the finite correction factor is `−(1+r/(2λ)) exp(−r/λ)`.

For `ℓ1⁴ < 4ℓ2⁴`, write `a=A+iB`, `A=ℓ1²/2`, `B=√(4ℓ2⁴−ℓ1⁴)/2`, choose `λc=√a` with positive real part, and `αc=−1/2+iA/(2B)`. The real template is `ÑTG=ÑN+2 Re[αc ÑY(λc)]`, with the complex Yukawa torque evaluated from the same geometry integral (or an exactly equivalent real damped-oscillatory kernel). Neither this branch nor the repeated-root limit follows from the published `|α|=1` curve.

A defensible table-level implementation would proceed as follows:

1. Ingest the 95 Appendix A.1 rows without binning; retain the three torques, their printed 1σ errors, each `s` and `δs`, and the common in-phase convention. Do not digitize Fig. 5 in place of these rows.
2. Implement the full linked-kernel torque for each pattern/harmonic, including the thesis geometry, separation offset, common scale, tilt, overcut, runout and radial-misalignment response. Confirm the table’s `s` convention against the chapter’s `sj = zP,j + zA + s0` definition before adding `s0`, to avoid double-counting an offset.
3. Reproduce the published Newtonian reference (`χ²=274.99`, 285 DOF, `P=.654`) with its five floated nuisance parameters, and then reproduce the sign-resolved single-Yukawa results in Table 11.4. Extend only after these checks pass. For the released table, the directly stated objective is Equation 11.1 with its diagonal torque/gap-error terms and Gaussian nuisance penalties. A reduced *centered* version must set science-run lateral offsets to nominal center explicitly, because those runwise offsets are absent; it can only be labeled conditional/exploratory.
4. For every `(ℓ1,ℓ2)`, profile the nuisance parameters and compare to the global best fit. Scan the full allowed domain, including `ℓ2=0`, the repeated-root locus, and both real-root and complex-root sectors. Calibrate the 95% likelihood-ratio threshold with pseudoexperiments under the chosen error model, especially near `ℓ2=0` and the repeated-root locus where the separate-root parameterization is singular; the paper’s single-Yukawa threshold/interval is not a two-parameter threshold. A nominal two-parameter Wilks contour may be a diagnostic, but is not a substitute for that calibration.

## Feasibility boundary and safest fallback

The missing science-run `xj,yj` positions prevent applying the authors’ runwise radial-misalignment correction; the missing dial-indicator, capacitance, and actuator readings prevent reconstructing the exact per-run separation/geometry record and its shared calibration errors. A.2 centering scans and the A.1 `s` uncertainties do not fill this gap. The covariance and systematic-error correlations are also unavailable, and the generated Fourier–Bessel/Yukawa tables or code cannot be checked against the gated APS supplement. One can still compute a centered, diagonal table profile, but it is not the full data analysis and cannot support a calibrated 95% two-range contour on the evidence in hand.

Safest bounded fallback: keep the published sign-resolved single-Yukawa constraint only on the exact single-length branch, where `ℓ2=0`, `α=−1`, and `λ=ℓ1`. Table 11.4’s negative-strength entries bracket the `α=−1` crossing between 25 µm (−1.27) and 28 µm (−0.886); the prior audit’s ≈27.1 µm linear read-off is indicative only and is not a new one-parameter recalibration. Do not carry that number into the linked two-range or complex-root sectors. The generic 38.6 µm `|α|=1` statement is not a full-model bound either.

For this campaign, leave the torsion-balance profile secondary to the Warden’s galaxy-envelope work. Revisit a full laboratory contour if the UW runwise stage/capacitance data and the APS supplement or equivalent torque tables become available. A centered table-only diagnostic can be done later, with its assumptions and lack of calibrated coverage stated plainly.

## Checked source locators

- Lee et al. primary paper: raw-data description and run coordinates, accepted manuscript pp. 1–2; single-Yukawa result and supplement reference, pp. 4–5. [APS accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevLett.124.101101).
- Lee dissertation: §11.1, printed pp. 123–125 (95 runs, phase); §11.2, pp. 126–130 (geometry, coordinates, model nuisance); §11.3, pp. 131–133 (χ², nuisance table/reference fit); §11.4, pp. 133–136 (single-Yukawa limits); Appendix A.1, pp. 146–150 (torque rows and `.sum` location); Appendix A.2, pp. 150–152 (centering scan rows).
- APS supplement status: [PRL article/supplement page](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.101101) labels Supplemental Material as subscription-required; no supplement file was supplied locally.

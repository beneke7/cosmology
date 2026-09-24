# DES-Dovekie supernova-only profile screen

Task status: complete for the bounded first pass. Scientific status: exploratory profile calculation, with input and algebra checks independently passed. This is not a run of the official CosmoSIS pipeline and does not reproduce its published posterior constraints.

## Result

The same 1820-row DES-Dovekie Hubble diagram and full STAT+SYS precision were fitted under flat LCDM, constant-w, and CPL backgrounds. Each model profiles one additive magnitude offset. The primary comparison uses the campaign's common search domains: Ωm∈[0.05,0.6], w or w0∈[−2,−0.3], wa∈[−3,3].

| Model | Profile best fit | Profile χ² | Δχ² from LCDM | Offset (mag) | Domain edge |
|---|---|---:|---:|---:|---|
| LCDM | Ωm=0.330317 | 1631.421 | 0 | +0.007449 | none |
| constant-w | Ωm=0.261042, w=−0.832680 | 1630.177 | −1.243 | +0.016628 | none |
| CPL | Ωm=0.409838, w0=−0.788386, wa=−3.000000 | 1626.631 | −4.790 | +0.029198 | wa lower edge |

No model used BAO, CMB, or another SN compilation. Δχ² here is a difference between bounded profile minima, not a posterior probability, evidence ratio, or discovery significance. The CPL score improves as the wa lower bound is widened from −3 to −5; the minimum remains at that edge for wa∈[−5,3]. It reaches an interior minimum at wa≈−6.7818 for wa∈[−10,3], and the result is numerically unchanged for wa∈[−20,3]: Ωm≈0.465666, w0≈−0.538859, χ²=1625.308492, Δχ²=−6.112044. All 12 starts in the last two domains converged to the same minimum to the reported precision. This sensitivity is recorded in `experiments/dovekie_screen/boundary_sensitivity.json`.

The wider-domain CPL fit allows a strongly evolving equation of state. Its modest profile improvement over a one-parameter LCDM shape does not include an extra-parameter penalty or account for model selection, priors, or posterior volume. Do not read it as evidence for evolving dark energy.

## Data and likelihood contract

Inputs are the pinned official DES-SN5YR repository commit `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`:

- `context/data/des_dovekie_hd.csv`, SHA-256 `2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf`.
- `context/data/des_dovekie_stat_sys.npz`, SHA-256 `ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b`.

The pinned release README states that the Hubble diagram is ordered for the covariance and warns against substituting the differently ordered metadata file. The official HD and STAT+SYS files were retrieved from the same immutable commit. The parser retains the source `SN:` order without sorting or joining, checks 1820 distinct CIDs and 1820 precision rows, and applies the official likelihood cut zHD>0.00; all 1820 rows remain. The packed precision archive contains `nsn=1820` and a float32 upper triangle with 1,657,110 entries. It is unpacked and reflected to form the symmetric float64 matrix W=Cov(STAT+SYS)⁻¹. W is finite, exactly symmetric, and Cholesky-positive. The paired covariance has no per-row IDs, so the pairing check rests on the pinned-release provenance and ordering rule rather than a separate CID-to-matrix index file.

The pinned likelihood uses zHD (CMB-frame redshift with peculiar-velocity correction) for the cosmological distance integral, zHEL (heliocentric redshift without that correction) in the luminosity-distance prefactor, and MU directly. MUERR and the diagnostic MUERR_SYS are not added again: the supplied precision is already the full STAT+SYS likelihood matrix. The release README describes MU as a distance modulus formed under H0=70. That normalization is not used to calibrate H0 here.

The official module has a data-reader inconsistency: at the same pinned commit it calls Astropy `Table.read(..., format='ascii.csv')`, but the pinned file contains `VARNAMES:` and `SN:` whitespace records and no commas. This implementation reads that native official format and follows the remaining release semantics. It therefore provides an independent evaluation of the release likelihood definition, not an execution of the released module.

The covariance algebra is:

\[
W=C_{\rm STAT+SYS}^{-1},\qquad
d(\theta)=\mu_{\rm base}(\theta)-\mu_{\rm obs},\qquad
\widehat{\mathcal M}=-\frac{\mathbf 1^T Wd}{\mathbf 1^T W\mathbf 1},
\]

\[
\chi^2_{\rm prof}=d^TWd-\frac{(\mathbf 1^TWd)^2}{\mathbf 1^TW\mathbf 1}
=[d+\widehat{\mathcal M}\mathbf 1]^T W[d+\widehat{\mathcal M}\mathbf 1].
\]

The calculation evaluates the last expression to avoid subtracting large quadratic terms. The release code's analytically integrated offset adds log[(1ᵀW1)/(2π)] and the fixed Gaussian covariance log determinant. These terms do not depend on the cosmological parameters for this fixed covariance, so their omission does not change profile minima or Δχ². The offset and MU are in magnitudes; redshifts and E(z) are dimensionless; the distance calculation uses Mpc and the standard c=299792.458 km/s. The distance modulus is

\[
\mu_{\rm base}=5\log_{10}\!\left[(1+z_{\rm HEL})\frac{c}{H_{0,\rm gauge}}
\int_0^{z_{\rm HD}}\frac{dz'}{E(z')}\right]+25,
\]

with an arbitrary fixed H0,gauge=70 km/s/Mpc. This follows the release factor `(1+zHD)(1+zHEL) D_A(zHD)` using flat `D_A(z)=(c/H0) integral/E /(1+zHD)`. The fitted additive offset absorbs any change in the distance normalization. An explicit check changing the gauge from 70 to 100 km/s/Mpc shifted the offset by 0.7745097999 mag, exactly 5 log10(100/70), while every profile χ² changed by less than 6×10⁻¹². H0 is not inferred.

The screened flat expansion laws are

\[
E^2_{\Lambda}=\Omega_m(1+z)^3+(1-\Omega_m),
\]

\[
E^2_{w}=\Omega_m(1+z)^3+(1-\Omega_m)(1+z)^{3(1+w)},
\]

\[
E^2_{\rm CPL}=\Omega_m(1+z)^3+(1-\Omega_m)(1+z)^{3(1+w_0+w_a)}
e^{-3w_a z/(1+z)}.
\]

Radiation, curvature, perturbations, posterior sampling, and external calibration are outside this late-time background screen. The optimization bounds are search domains, not priors.

## Numerical checks and reproducibility

- `work/data_audit/validate_data.py` independently verified both input hashes, row counts, the 1820×1820 precision shape, positive-definite factorization, and a precision 2-norm condition estimate κ≈5.95424×10⁷ (eigenvalue range 4.56550×10⁻⁶ to 271.84070).
- The direct profile quadratic agreed with the Schur-complement expression and a Cholesky-whitened evaluation to <9.1×10⁻¹³ in χ²; the offset normal-equation residual was below 8.3×10⁻¹³.
- 96-point Gauss-Legendre integration and an independent 160-point calculation differed by at most 7.2×10⁻¹⁵ mag at the minima.
- The three model minima used 12 deterministic L-BFGS-B starts each and had no LCDM/constant-w boundary hits. All 36 baseline starts succeeded; the largest spread in profile χ² among starts was 2.22×10⁻¹⁰. Every widened-bound start also succeeded; profile χ² spreads were ≤1.52×10⁻¹⁰ for wa∈[−5,3], ≤1.75×10⁻⁷ for [−10,3], and ≤8.43×10⁻⁸ for [−20,3]. Small last-digit parameter variation in the two widest cases reflects a shallow profile direction. CPL hits wa=−3 in the baseline domain; the separately recorded wider-domain search resolves that local box boundary.
- The first invocation completed fitting and numerical checks but failed during Matplotlib label rendering (`\\mathcal` was not accepted by MathText) before it wrote result.json. The text label was simplified and the exact full command was rerun successfully; the failure remains in the machine-readable run history.

Primary command:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python work/inference/fit_dovekie.py --starts 12
```

Boundary-sensitivity command:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python work/inference/check_cpl_boundary.py
```

The primary fit and checks took 7.04 s (7.26 s through plot output), and the three widened-bound runs took 16.13 s on the recorded CPU environment. Python was 3.12.3, NumPy 2.5.3, SciPy 1.18.1, and Matplotlib 3.11.2. The project preflight reported 20 visible logical CPUs and an RTX 5090; this task used one CPU BLAS thread and no GPU.

Artifacts:

- `experiments/dovekie_screen/result.json`: fit parameters, domains, offsets, algebra checks, data/code hashes, commands, runtime, and limits; SHA-256 `55f9ab7cd1fbc5417a8a23703b094f9c36d635ca21dfae505a0e5bff129cd0e5`.
- `experiments/dovekie_screen/boundary_sensitivity.json`: wa-bound sensitivity runs, their parameters, checks, hashes, and runtime; SHA-256 `dd964607945ae059cb6aedee6c93e570d3fdab69fec5655b801f7eb30c5fbd38`.
- `experiments/dovekie_screen/dovekie_hubble_residuals.png`: residual and profiled model-difference visualization; SHA-256 `1ef79c39e1c2199400ee3a80a2c87288ef2f4013d3587bcb5e426d8c8bc57e4d`. The display bin medians are descriptive only; inference uses the full correlated precision.
- `work/inference/fit_dovekie.py`: exact primary implementation; SHA-256 `f796afccedf91b61af5b0cfc18281a4638a96d9abcfa9618f3ea79b672b10ce1`.
- `work/inference/check_cpl_boundary.py`: exact boundary-sensitivity implementation; SHA-256 `61602b976c364c4c3a58391821723dd24342a161feda68f770c23cfb79faedd7`.

Official release sources:

- [Pinned DES-Dovekie release README](https://raw.githubusercontent.com/des-science/DES-SN5YR/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/README.md)
- [Pinned DES-Dovekie Hubble diagram](https://raw.githubusercontent.com/des-science/DES-SN5YR/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv)
- [Pinned STAT+SYS archive](https://raw.githubusercontent.com/des-science/DES-SN5YR/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/STAT+SYS.npz)
- [Pinned CosmoSIS likelihood module](https://raw.githubusercontent.com/des-science/DES-SN5YR/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/5_COSMOLOGY/Dovekie_cosmosis_likelihood.py)
- [DES-Dovekie paper](https://arxiv.org/abs/2511.07517)

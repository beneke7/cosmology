# Numerical scope and data contract

`background_bao.py` is a low-redshift distance calculator and optional compressed-BAO likelihood fit. It assumes a spatially flat matter + dark-energy background with the CPL form

\[
E^2(z)=\Omega_m(1+z)^3+(1-\Omega_m)(1+z)^{3(1+w_0+w_a)}
\exp\!\left[-3w_a\frac{z}{1+z}\right].
\]

Radiation is deliberately neglected; the numerical domain is limited to `0 <= z <= 3`. This excludes early-universe physics. The code uses the amplitude `alpha = c/(H0 rd)` and returns `DM/rd = alpha * integral(dz/E)`, `DH/rd = alpha/E`, and `DV/rd = (z DM^2 DH)^(1/3)`. It does not predict `rd`, nor infer `H0` and `rd` separately.

## Public text format

The official DESI DR1 BAO-results documentation says the likelihood inputs are publicly available in the CobayaSampler `bao_data` repository. That repository identifies DESI DR1 and DR2 BAO releases, and the DR2 `*_mean.txt` file documents a whitespace-separated `z value quantity` layout. Its observable labels are `DV_over_rs`, `DM_over_rs`, and `DH_over_rs`; `background_bao.py` accepts exactly those labels and rejects unknown names. Example direct references:

- [DESI DR1 BAO-results documentation](https://data.desi.lbl.gov/doc/releases/dr1/vac/bao-cosmo-params/)
- [Public BAO likelihood-data repository and provenance](https://github.com/CobayaSampler/bao_data/blob/master/README.md)
- [DESI DR2 example mean vector](https://github.com/CobayaSampler/bao_data/blob/master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt)
- [Matching DESI DR2 covariance matrix](https://github.com/CobayaSampler/bao_data/blob/master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt)

The covariance is a plain numeric matrix without row labels. Its row and column order must match the mean-file rows exactly. The loader preserves mean rows unchanged and checks dimensions, finiteness, symmetry, and positive definiteness; it cannot independently establish source-file pairing or recover a reordered covariance from an unlabeled matrix.

## Numerics and interpretation

Predictions use NumPy float64 and 96-point fixed Gauss-Legendre quadrature. `--self-check` compares the vector quadrature with independent `scipy.quad` calculations and verifies the flat-LCDM limit, the Einstein-de Sitter analytic distance, covariance checks, row-order preservation, and recovery of a synthetic CPL model at the nested Lambda point. Correlated Gaussian chi-square uses a Cholesky factor and triangular solve, with no explicit covariance inverse. Optional fits use deterministic multi-start bounded optimization: `0.05 < Omega_m < 0.6`, `-2 < w0 < -0.3`, `-3 < wa < 3`, and positive `alpha` (implemented with a small positive numerical lower bound). LCDM and constant-w CDM fix `w0, wa` as appropriate.

Input-data fits emit JSON marked as exploratory and background-only. Their chi-square and parameter estimates are not posterior constraints, evidence, or a significance calculation. No real observational fit is included in this starter. The CPU float64 path is the reference; the benchmark reports if optional PyTorch/CUDA is absent rather than requiring an accelerator.

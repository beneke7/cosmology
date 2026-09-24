# Data audit: DESI DR2 BAO and DES-Dovekie

As of 2026-09-24. Status: DESI DR2 inputs independently checked; DES-Dovekie data and source code are publicly accessible and independently checked; a direct run of the published CosmoSIS plug-in is blocked by an input-reader mismatch and missing dependencies in the current project environment. No cosmological fit was run and no data were synthesized.

## DESI DR2 BAO

The local 13-row mean and 13×13 covariance are byte-for-byte the files at the immutable `CobayaSampler/bao_data` commit `bb0c1c9009dc76d1391300e169e8df38fd1096db` (2025-06-26). The DESI DR2 paper’s Data Availability section directs readers to this repository for the DR2 BAO likelihood; the official DESI data documentation separately lists the released DR2 cosmology chains and posterior maximizations, and says the underlying spectra/redshifts were not yet released. This is therefore the public compressed BAO likelihood input, not a reanalysis of the raw catalogues.

Sources: [DESI DR2 Results II](https://arxiv.org/abs/2503.14738), [official DESI paper/data index](https://data.desi.lbl.gov/doc/papers/), [pinned BAO mean](https://raw.githubusercontent.com/CobayaSampler/bao_data/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt), [pinned BAO covariance](https://raw.githubusercontent.com/CobayaSampler/bao_data/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt).

The released coordinates are dimensionless distance ratios: (D_V/r_s), (D_M/r_s), and (D_H/r_s). The files label the sound-horizon denominator `rs`; the local background code calls it (r_d). No Mpc or (h^{-1}) factor is needed in this compressed vector. The exact row order and quoted means are:

| Row | (z\) | Observable | Mean | σ from covariance |
|---:|---:|---|---:|---:|
| 1 | 0.295 | (D_V/r_s) | 7.94167639 | 0.076092 |
| 2 | 0.510 | (D_M/r_s) | 13.58758434 | 0.168367 |
| 3 | 0.510 | (D_H/r_s) | 21.86294686 | 0.428868 |
| 4 | 0.706 | (D_M/r_s) | 17.35069094 | 0.179931 |
| 5 | 0.706 | (D_H/r_s) | 19.45534918 | 0.333870 |
| 6 | 0.934 | (D_M/r_s) | 21.57563956 | 0.161782 |
| 7 | 0.934 | (D_H/r_s) | 17.64149464 | 0.201043 |
| 8 | 1.321 | (D_M/r_s) | 27.60085612 | 0.324556 |
| 9 | 1.321 | (D_H/r_s) | 14.17602155 | 0.224551 |
| 10 | 1.484 | (D_M/r_s) | 30.51190063 | 0.763558 |
| 11 | 1.484 | (D_H/r_s) | 12.81699964 | 0.518012 |
| 12 | 2.330 | (D_H/r_s) | 8.6315456748 | 0.101062 |
| 13 | 2.330 | (D_M/r_s) | 38.9889739620 | 0.531682 |

In particular, the final high-redshift pair is (D_H) then (D_M). The covariance must be indexed in this exact order. Independent byte checks reproduce the source repository's Git blob IDs (`8aff444f…` mean, `fd8e5697…` covariance) as well as the SHA-256 digests below. Numerical checks give exact symmetry, successful Cholesky factorization, λmin = 0.00578998687, λmax = 0.67640843627, and 2-norm condition number κ₂ = 116.823829. The file is block diagonal by tracer/redshift block, but each of its six (D_M,D_H) pairs is correlated (largest absolute pair correlation 0.493552); retaining the full matrix is required.

The mean entries round to the central distance values in Table IV of the DESI paper. The paper describes Table IV uncertainties as standard deviations of marginalized BAO-parameter posteriors, while the repository files are explicitly named as a Gaussian BAO compression. The supplied covariance diagonal does not exactly reproduce those tabulated marginal errors: the largest relative difference is 6.4% for (D_M/r_d) at (z=0.934) (0.1618 from the matrix versus 0.152 in Table IV); the corresponding pair correlation is −0.3472 versus −0.416 in the table. The pinned files match the paper-cited data repository exactly, so this is not a local corruption. The audited sources do not explain the residual difference; preserve it as a compression-versus-published-posterior ambiguity and do not tune the covariance to Table IV.

The current code is [scripts/background_bao.py](../../scripts/background_bao.py), SHA-256 `34bdf1d81397b39b740faa3e2f4a79748d5fa9d3bfdd8303619d631dbf8aa2e2`. `BAOData.from_files` checks the 13 aligned rows, finite symmetric covariance and Cholesky factorization. `GaussianBAOLikelihood.chi2` computes (r^T C^{-1}r) through a Cholesky solve. Its predictors are

\[
E^2(z)=\Omega_m(1+z)^3+(1-\Omega_m)(1+z)^{3(1+w_0+w_a)}e^{-3w_a z/(1+z)},
\]
\[
D_M/r_d=\alpha\int_0^z dz'/E(z'),\quad D_H/r_d=\alpha/E(z),\quad
D_V/r_d=[z(D_M/r_d)^2(D_H/r_d)]^{1/3},
\]

with α = (c/(H_0r_d)) fitted independently. This avoids claiming an (H_0) measurement without a sound-horizon calibration. The local models are flat ΛCDM, flat constant-(w) CDM and flat CPL. They omit radiation and curvature and use 96-point Gauss-Legendre distance quadrature for (0<z\le3). The included optimizer minimizes profile χ²; its bounds are (\alpha\in[10^{-6},10^4]), Ωm ∈ [0.05,0.6], (w_0\in[-2,-0.3]), and (w_a\in[-3,3]) as applicable. Those are optimizer bounds, not posterior priors. The fixed-covariance Gaussian determinant is absent from the function because it is constant for profile optimization; the code produces neither normalized posterior samples nor evidence. This is a transparent screening likelihood, not the full DESI collaboration inference. The covariance provenance and current likelihood definition agree at the compressed-Gaussian level; the result should remain labeled as background-only until radiation/model precision effects are assessed.

## DES-Dovekie access and checks

The official `des-science/DES-SN5YR` repository exposes the paired data and likelihood at commit `c9a4fcafc4cbd19bd750dee47fc76194a45c181f` (2026-01-28 commit; no formal Dovekie tag was identified, so this audit pins the repository commit). The exact paths are [DES-Dovekie_HD.csv](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv), [STAT+SYS.npz](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/STAT%2BSYS.npz), [distance/covariance README](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/README.md), and [CosmoSIS likelihood](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/5_COSMOLOGY/Dovekie_cosmosis_likelihood.py). The repository tree reports sizes 148,002 bytes and 6,244,951 bytes for the two data files; local Git blob IDs match the pinned tree (`f80ec4e2…`, `42896664…`). The upstream likelihood code (7,674 bytes; blob `b7142093…`) and two README files are retained in this audit directory. The official Dovekie paper is [arXiv:2511.07517](https://arxiv.org/abs/2511.07517).

The Hubble diagram has 1,820 rows and 1,820 unique `CID`s. It contains `zHD` (CMB-frame redshift with peculiar-velocity correction), `zHEL`, `MU` (distance modulus in magnitudes), `MUERR`, and the survey/classification fields. All 1,820 rows have (zHD>0), so the official plug-in's `zHD > 0` selection retains every row without changing order. The official README warns that the metadata table has a different ordering and must not be used with this covariance. Keep the Hubble diagram in its release order.

`STAT+SYS.npz` is the inverse of the full statistical-plus-systematic covariance, not a covariance matrix. It stores 1,657,110 float32 values, the upper triangle of a 1,820×1,820 matrix. Reconstructing that precision matrix according to the official source code gives a successful Cholesky factorization, eigenvalue range (4.56550\times10^{-6}) to 271.84070, and κ₂ ≈ (5.95424\times10^7). Since all rows survive the release selection, the inverse precision condition number is also that of the corresponding covariance. Do not add `MUERR` to this full covariance a second time.

The checked-in likelihood defines the model distance modulus using interpolated (D_A(zHD)):

\[
\mu_{model}=5\log_{10}[(1+zHD)(1+zHEL)D_A(zHD)]+25.
\]

For residual vector (d=\mu_{model}-\mu_{obs}) and release precision (P=C^{-1}), it analytically marginalizes one constant magnitude offset (M), which is fully degenerate with (H_0): (d^TPd-(\mathbf1^TPd)^2/(\mathbf1^TP\mathbf1)+\log[(\mathbf1^TP\mathbf1)/(2\pi)]). The CosmoSIS Gaussian base class also supplies the fixed covariance normalization. The supernova sample does not independently measure (H_0). These semantics, especially the paired redshifts, order and full inverse STAT+SYS, must be preserved in an adapter.

There are two concrete blockers to a direct run of the downloaded plug-in as committed. First, it calls `Table.read(..., format='ascii.csv')`, but the paired `.csv` file is actually a whitespace SNANA table (`VARNAMES:` header and `SN:` rows) with no commas; [Astropy documents `ascii.csv` as comma-separated](https://docs.astropy.org/en/stable/io/ascii/index.html). Second, its defaults point to `DES-SN5YR_HD.csv` and `STAT+SYS.txt.gz`, not these Dovekie release paths. The correct data paths must be passed, and the Hubble diagram must be read with a parser for its actual format or converted without reordering. The project `.venv` currently has NumPy 2.5.3 and SciPy 1.18.1 but lacks Astropy and CosmoSIS, so the official plug-in has not been executed. A custom parser and a focused likelihood-initialization/known-value smoke check are required before calling it runnable. The format-level check and a dependency check are not a cosmology fit.

For the simple local background comparison, an honest Dovekie branch is feasible with these two public files and an adapter implementing the released distance-modulus and offset-marginalization formulas. A reproduction of the paper's published multi-probe posterior or evidence additionally needs the stated model and sampler priors, matching curvature/model choices and sampler, plus its CMB and other probe likelihoods. The paper lists evidence priors (h\in(0.55,0.91)), Ωm ∈ (0.1,0.5), Ωk ∈ (-0.15,0.15), (w_0\in(-3,-0.4)), (w_a\in(-3,2)), and (w_0+w_a<0), with Σmν fixed at 0.06 eV; it uses CosmoSIS with Nautilus for nominal constraints. The flat, no-radiation local BAO screening code is not that published setup.

Access is public and anonymous at the pinned raw GitHub URLs. The official repo metadata has no detected `LICENSE` (`license: null` in the GitHub API); the data-release documentation asks users to cite the DES-SN5YR papers but does not state an explicit reuse license. The BAO data repo likewise reports no machine-detected license. Record this as unspecified rather than assigning a license.

## Likelihood contract and overlap

| Factor | One likelihood object | Treatment |
|---|---|---|
| DESI DR2 BAO | 13 distance-ratio entries across seven tracer/redshift blocks | One correlated Gaussian using the full released (C); no separate BAO-chain factor. |
| DES-Dovekie SNe | 1,820 released Hubble-diagram entries | One Gaussian using full `STAT+SYS` precision, with the constant magnitude/H0 offset marginalized. |

Treat DESI BAO and Dovekie as conditionally independent given the background parameters only as a recorded approximation: they are different observations, and the audited releases supply no cross-probe covariance. Keep standalone predictions separate if a zero cross-covariance assumption is not accepted. DESI posterior chains or a DES-Dovekie chain already containing BAO are alternative summaries of those same inputs, never extra factors. Do not multiply Dovekie by the original DES-SN5YR likelihood: the paper reports 1,718 overlapping SNe between the 1,820-event Dovekie and 1,829-event DES-SN5YR samples. Do not multiply Pantheon+ either; the paper describes low- and high-redshift SNe common to Pantheon+ and DES-SN5YR. Keep Pantheon+ as a separate robustness branch.

## Reproduction record

Smoke command:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./.venv/bin/python work/data_audit/validate_data.py
```

It verified pinned SHA-256 values, parsed DESI through the current `BAOData.from_files`, checked expected row order, symmetry, block structure, SPD and conditioning, parsed Dovekie's SNANA rows without sorting, reconstructed its packed inverse covariance, and checked its SPD/conditioning. It intentionally does not call `fit_model`. Data-file hashes and the exact fetched lock entries are recorded in `context/manifest.lock.json`; added data totals 6,392,953 bytes. No bootstrap-budget expansion was used.

| Local input | Bytes | SHA-256 |
|---|---:|---|
| `context/data/desi_dr2_mean.txt` | 472 | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `context/data/desi_dr2_cov.txt` | 2,547 | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| `context/data/des_dovekie_hd.csv` | 148,002 | `2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf` |
| `context/data/des_dovekie_stat_sys.npz` | 6,244,951 | `ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b` |
| official Dovekie CosmoSIS likelihood source | 7,674 | `78526c0e6013619e3586a5509be504f0bc5f0461d1e55a3d5bec0154353645ff` |
| DESI DR2 paper PDF already in context | 12,175,089 | `1e82f26e4cc3901b16168cd147f252bfa804f9c3caad3f4f7e3532640d237841` |
| DES-Dovekie paper PDF already in context | 2,691,260 | `cd3dd615c36afe0d45bcdfabd9c31611f0ca511e325997f93c3eb1f5a9744eb6` |

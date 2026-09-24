# Independent DES-Dovekie profile audit

Status: independently checked. The requested DES-Dovekie profile points and the reported bounded minima reproduce with a separately written parser, distance evaluator, full-matrix profile calculation, and seeded bounded global optimizer. This is a background-only profile validation; it is not a posterior, evidence, or significance calculation.

## Inputs and ordering

The run used the paired files from DES-SN5YR repository commit `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`:

| Input | SHA-256 |
|---|---|
| `context/data/des_dovekie_hd.csv` | `2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf` |
| `context/data/des_dovekie_stat_sys.npz` | `ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b` |

My parser reads the `VARNAMES:` / `SN:` whitespace records in file order. It found 1,820 rows, 1,820 distinct CIDs, and every row passes the official `zHD > 0` selection. The first and last IDs are `Gaia16agf` and `1257587`; the SHA-256 of the ordered CID list joined with newlines is `9c7ebdc43bd32c00623c43d0155c1d4bcfe54391fd79576447cd5a4c58d18683`. The precision archive reports `nsn=1820` and contains 1,657,110 float32 values, exactly the upper triangle for an 1,820×1,820 matrix.

I reconstructed and reflected that upper triangle as the symmetric precision `P`, checked that it is finite and Cholesky-positive, and used the full dense matrix. `1ᵀP1 = 43938.33706138445`. I did not add `MUERR` again because the paired release matrix is STAT+SYS. This matches the release README's statement that the distributed matrices are inverse covariances and its order warning, and the release likelihood's unpacking code. The Hubble diagram and precision have no common row IDs inside the NPZ, so row-to-row pairing is established by the pinned release and its ordering statement; it cannot be proven from an independent ID column in the matrix.

The official code uses `zHD` for the distance interpolation and both redshift factors in `5 log10[(1+zHD)(1+zHEL) D_A(zHD)] + 25` (see [`Dovekie_cosmosis_likelihood.py`](../data_audit/Dovekie_cosmosis_likelihood.py), lines 108–115). For a flat model, `D_A(zHD) = (c/H0) integral[0,zHD](dz/E)/(1+zHD)`, leaving

`mu_base = 5 log10[(1+zHEL) (c/H0_gauge) integral[0,zHD](dz/E)] + 25`.

I evaluated this with 128-point Gauss-Legendre quadrature and an arbitrary `H0_gauge=70 km/s/Mpc`. The profiled constant absorbs that gauge; the reported offsets are in magnitudes and do not measure `H0`.

## Requested fixed points

For `d = mu_base - MU`, I independently profiled the additive magnitude offset as `Mhat = -(1ᵀPd)/(1ᵀP1)`, then evaluated `(d + Mhat*1)ᵀP(d + Mhat*1)`. A reference value was subtracted from `d` before matrix products to improve numerical conditioning; this does not change the profiled result. The direct residual quadratic, the Schur-complement expression, and Cholesky whitening agreed to at most `1.37e-12` in chi-square. The largest offset normal-equation residual was `1.10e-12`.

| Fixed parameter point | Profile χ² | Profile offset (mag) | Δχ² vs. ΛCDM |
|---|---:|---:|---:|
| ΛCDM, Ωm = 0.3303167905483107 | 1631.420535571240 | +0.007449125568 | 0 |
| Constant-w, Ωm = 0.2610419, w = −0.83268 | 1630.177096764439 | +0.016628720155 | −1.243438807 |
| CPL, Ωm = 0.409838, w0 = −0.788386, wa = −3 | 1626.630538385166 | +0.029198374354 | −4.789997186 |
| CPL, Ωm = 0.465666, w0 = −0.538859, wa = −6.7818 | 1625.308491819304 | +0.042686489949 | −6.112043752 |

These agree with the values in [`work/inference/REPORT.md`](../inference/REPORT.md) to its stated precision, including the widened-CPL score. The release likelihood also adds `log[(1ᵀP1)/(2π)]` and a parameter-independent Gaussian covariance normalization; these constants do not affect the reported profile minima or Δχ².

## Independent bounded search

I also ran seeded `scipy.optimize.differential_evolution` searches over Ωm ∈ [0.05, 0.6], w/w0 ∈ [−2, −0.3], wa ∈ [−3, 3] for the baseline CPL, and wa ∈ [−10, 3] for the widened CPL. The four searches all returned success. Their profile minima were:

| Model/domain | Independently searched parameters | Profile χ² |
|---|---|---:|
| ΛCDM | Ωm = 0.3303168052 | 1631.420535571239 |
| Constant-w | Ωm = 0.2610394933, w = −0.8326761519 | 1630.177096764740 |
| CPL, wa ≥ −3 | Ωm = 0.4098391748, w0 = −0.7883928176, wa = −2.9999999603 | 1626.630538419629 |
| CPL, wa ≥ −10 | Ωm = 0.4656656140, w0 = −0.5388875991, wa = −6.7816252104 | 1625.308491840818 |

The baseline CPL solution remains at the `wa = −3` boundary. The widened CPL search finds an interior solution near the reported `wa ≈ −6.78`; its slightly different last digits lie along a shallow profile direction, while chi-square agrees within `2.5e-8`. Across all four models, independent minimum chi-square agrees with the supplied result to better than `3.5e-8`. This refutes an implementation-specific optimizer or profile-algebra explanation for the reported bounded minima; it does not validate any posterior interpretation.

As a convention sensitivity check, substituting `zHD` for `zHEL` only in the luminosity-distance prefactor at each requested parameter point raises profile χ² by `0.35746`, `0.61273`, `0.75956`, and `0.87097` respectively. The fitted constant offset shifts by `+0.001075185 mag` in each case. The pinned likelihood's separate `zHD` integral and `zHEL` prefactor convention therefore matters at the profile-score precision shown here.

An auxiliary DE search of this altered prefactor convention used a 40-iteration cap; it converged for ΛCDM and constant-w but stopped at the cap for both CPL domains. Those two auxiliary fits are marked unsuccessful in `independent_result.json` and are not used in the comparison above. The correctly specified searches converged successfully.

## Reproduction record

Command:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python work/critic/independent_dovekie.py
```

Elapsed wall time was 21.252 s (21.218 s measured inside the script). Python 3.12.3, NumPy 2.5.3, and SciPy 1.18.1 were used; 20 CPUs were in affinity, but all numerical-library thread settings were explicitly set to one, with one DE worker and no GPU. The main DE searches used seed `90210 + model_index`, `maxiter=65`, `popsize=8`, `tol=2e-10`, and `atol=1e-9`, followed by bounded polishing. The alternative-prefactor probe used `maxiter=40`, `popsize=7`, and `tol=1e-8`.

Code SHA-256: `4aaf8a9921135b0b5bee1af0605156dca3755701af3c95a7b7c16944a50d15b5` (`independent_dovekie.py`). Machine-readable output SHA-256: `8c7a5a276b872eef48ee72cf04163ebf6ad15b3c2bf26cb43d00335c08a05118` (`independent_result.json`). The output includes the input hashes, exact fixed-point values, optimizer metadata, and the recorded auxiliary convergence failures.

No posterior samples, evidence, significance, or joint-probe score were calculated. Radiation and curvature remain outside this flat late-time profile screen.

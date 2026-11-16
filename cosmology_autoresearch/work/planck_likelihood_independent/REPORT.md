# One-point Planck PR3 likelihood smoke test

**Status: reproduced software/data smoke; not a parameter-baseline reproduction.** A pinned Cobaya/CAMB/clipy model initialized the three full PR3 likelihood components, passed each component's bundled check vector, and returned finite likelihoods at the published Table 2 TT,TE,EE+lowE cosmological marginal centers. This establishes a usable likelihood path; it does not reproduce Planck's posterior or best fit.

## Input and model contract

The input is [`planck2018_nolens_camb1166.yaml`](../planck_baseline/planck2018_nolens_camb1166.yaml), SHA-256 `c1d38c87d4035dbe7c9c2d323ddc5d1b58694d9ac4bec2f6ef3808ff95bac909`. It has exactly the full Plik HM TTTEEE, Commander TT, and SimAll EE PR3 products. There is no separate Planck lensing-reconstruction likelihood; physical lensing in the CAMB spectra remains on. The products were extracted from the validated official PR3 archive audited in [`planck_stage/REPORT.md`](../planck_stage/REPORT.md) (archive SHA-256 `0b73171e3acc671c28184466a45485a2d1c1d93676b832abdfe688c7b04024e6`).

| Component | Data bytes / files | Product-tree SHA-256 | Loglike |
|---|---:|---|---:|
| `planck_2018_highl_plik.TTTEEE` | 57,891,282 / 120 | `15a03d959bc896da327ed027161e6ddfcae31cde34b9b5663107a81f26e93f9a` | −1981.831938415346 |
| `planck_2018_lowl.TT_clik` | 6,509,031 / 7 | `6a4660592e3a194b6888f62385b56670fb992805ae527824e24cc62437cc083b` | −11.808515424172015 |
| `planck_2018_lowl.EE_clik` | 694,285 / 7 | `810f760e2e1f9fd5a817f70b885193bfae08564a279735d8ea8670e6cf504eb1` | −198.0775858760956 |

The CAMB helium contract is explicit in the config: `PArthENoPE_880.2_marcucci.dat`, PArthENoPE 1.1, neutron lifetime 880.2 s, the Planck VI case-(b) Marcucci et al. 2016 rate, and `nnu = standard_neutrino_neff = 3.046`. The staged table SHA-256 is `994c68f11ecbd33ecf6e631a0432ffc4b80068bfe972ab98fd235420aad67d8e`; the config manifest gives helium mass fraction Y_P = 0.2453946975 and 10^5 D/H = 2.45471022 at omega_b h^2 = 0.02236.

## Point evaluated and score

The six cosmological inputs are the Planck VI Table 2 fourth-column (`TT,TE,EE+lowE`) central marginalized values: omega_b h^2 = 0.02236, omega_c h^2 = 0.1202, 100 theta_MC = 1.04090, tau = 0.0544, ln(10^10 A_s) = 3.045, and n_s = 0.9649. The full 21-parameter sampled nuisance vector was generated reproducibly from the Cobaya reference/prior measure with seed `20260925`; the point was checked to have finite cosmological and external SZ prior densities. The remaining 26 fixed model/nuisance defaults are recorded in [`result.json`](result.json), alongside the complete sampled nuisance vector, every likelihood score, and exact dependency/data/config hashes.

| Quantity | Value |
|---|---:|
| Plik HM TTTEEE loglike | −1981.831938415346 |
| Commander TT loglike | −11.808515424172015 |
| SimAll EE loglike | −198.0775858760956 |
| **Summed loglike** | **−2191.718039715614** |
| Prior log-density parts (base, external SZ) | −28.333819946554264, −2.623848439487015 |
| Summed log-prior | −30.95766838604128 |

Independent `math.fsum` summation agrees with the reported loglike sum to 4.55e-13. Cobaya warns that the external SZ prior may not be normalized; do not interpret the absolute log-posterior density as a normalized model probability.

Crucially, the Table 2 entries are *marginal posterior centers*, not a nuisance-optimized point or published likelihood maximum. The nuisance values here are one seeded reference/prior draw. This score is not the Planck published best fit, must not be compared as such, and says nothing yet about recovering the six-parameter posterior intervals.

## Validation, runtime, and exact reproduction

On initialization, clipy 0.15 compared all bundled test vectors to stored expected values. Printed absolute differences were `4.00569e-6` for Plik (test loglike about −1172.47), `1.02746e-6` for Commander (about −11.6257), and `4.1778e-8` for SimAll (about −197.99). The evaluated point and all returned likelihoods were finite. CAMB supplied the pinned BBN table; Cobaya/CAMB/clipy model construction and the exact likelihood call completed without runtime error.

The local environment was Python 3.12.3, Cobaya 3.6.2, CAMB 1.6.6, clipy-like 0.15, NumPy 2.5.3, SciPy 1.18.1. Exact installed-distribution `METADATA` and `RECORD` hashes are in `result.json`. One CPU was available to this process (affinity pinned to core 19); OpenMP, OpenBLAS, MKL, and NumExpr threads were all set to 1. Initialization took 0.279 s, likelihood evaluation 2.495 s, total wall time 3.231 s. No GPU, sampler, MCMC, or chain was used; a single likelihood evaluation does not justify GPU dispatch.

From the project root, reproduce with:

```sh
timeout --signal=TERM 600s taskset -c 19 env \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  cosmology_autoresearch/.venv/bin/python \
  cosmology_autoresearch/work/planck_likelihood_independent/evaluate_point.py
```

The executable source SHA-256 for this run is recorded as `source_sha256` in `result.json` (also printed at run completion). The JSON artifact preserves the exact point, parameters, versions and scores.

## Limitations and next step

This is a likelihood-integration smoke and a one-point diagnostic only. The clipy stored-vector checks validate data/software compatibility but are not an independent compiled-clik comparison. We have not checked spectrum or likelihood convergence under CAMB accuracy changes, independently reproduced a likelihood term, found the nuisance maximum, run chains, or compared a posterior to Planck's public chains. A chain start at the Table 2 centers is not itself a reproduction. Next, run a short, checkpointed no-lensing TTTEEE+lowE Lambda-CDM chain with the verified contract, measure likelihood throughput/convergence, then extend until multiple chains meet the recorded convergence criterion and compare marginalized parameters and derived quantities to Table 2. Treat that as the first posterior-baseline reproduction; retain nuisance freedom and document any mismatch before using this baseline in extension tests.

Two setup-wrapper failures occurred before the successful evaluation and are preserved in [`attempts.json`](attempts.json): an initial incorrect project-root path, and then a dict/vector mismatch when checking Cobaya prior validity. The second occurred after all three bundled likelihood checks passed but before evaluating a scientific score. Both were corrected; the final invocation exited successfully.

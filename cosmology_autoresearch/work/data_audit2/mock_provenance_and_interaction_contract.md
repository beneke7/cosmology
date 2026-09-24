# DESI mock provenance audit

**Disposition: reproduced and independently checked.** The 1,200 fitted-flat-LCDM null vectors can be regenerated deterministically for the interacting-vacuum search. The original BAO/CPL result files retain model-fit records rather than raw 13-vectors; the later curvature file retains all 1,200 flat/curved fit records. This audit reconstructs six representative vectors and refits each with an independent adaptive-quadrature flat-LCDM implementation. It did not fit any interacting-vacuum mocks.

## Generator and source verification

The two inputs are the 200-realization result (`seed=20260924`) and the 1,000-realization extension (`seed=20260925`). Both record matching current hashes for `scripts/run_bao_robustness.py`, `scripts/background_bao.py`, and the ordered DESI DR2 mean/covariance. The curvature record's hashes of both input result files, its fitter, result JSON and complete JSONL also all match the files on disk.

The exact draw recipe is:

```text
mean = fitted flat-LCDM prediction from experiments/campaign_seed_bao/result.json
normals = np.random.default_rng(seed).standard_normal((N, 13))  # PCG64
L = np.linalg.cholesky(full_covariance)                         # lower triangular
mock[i] = mean + normals[i] @ L.T                              # i is zero-based
```

Seed runs are regenerated as separate complete arrays, not as slices from one shared generator. Their normal-array digests match the stored curvature audit (`a5265e…` for 20260924; `b9cf26…` for 20260925); there are no identical 13-normal rows across the combined finite arrays. These checks establish distinct initialized streams and no observed row collision; they are not a mathematical proof that two infinite pseudorandom sequences can never intersect.

The mean hash is `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585`; the full covariance hash is `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`. The covariance is symmetric and Cholesky-positive, with minimum eigenvalue 0.00578998687 and 2-norm condition number 116.824. Its unchanged row order includes the older Ly-alpha BAO rows; a newer Ly-alpha result cannot be appended without its joint covariance.

## Independent fit spot checks

For each seed, the checker regenerated the first, middle and last indexed vector and independently refit flat LCDM using adaptive scalar quadrature, Cholesky-whitened residuals, analytic GLS profiling of the dimensionless amplitude α, and a bounded one-dimensional profile in Ωm. All six stored flat fits report optimizer success.

| Seed | Mock index | Stored χ² | Independent Δχ² | Maximum prediction difference |
|---:|---:|---:|---:|---:|
| 20260924 | 0 | 6.5718206795 | +1.3×10⁻¹³ | 6.3×10⁻⁹ |
| 20260924 | 99 | 8.3674465335 | −2.0×10⁻¹⁴ | 2.0×10⁻⁸ |
| 20260924 | 199 | 13.6008461972 | +6.750155989720952e-14 | 1.5×10⁻⁸ |
| 20260925 | 0 | 4.8868224404 | −6.0×10⁻¹⁴ | 2.8×10⁻⁸ |
| 20260925 | 499 | 3.9540899301 | −8.9×10⁻¹⁶ | 1.8×10⁻⁸ |
| 20260925 | 999 | 17.9532563432 | −1.3×10⁻¹³ | 2.3×10⁻⁸ |

The independent flat scores agree to (1.4\times10^{-13}) or better. Predictions differ by at most (2.9\times10^{-8}), consistent with the independent optimizer locating Ωm within (4\times10^{-9}) and the original and checker using different quadrature/profile implementations. The independently recomputed null mean differs from its stored prediction by (2.9\times10^{-14}). Exact selected-vector hashes and full-precision comparisons are in [provenance_check.json](provenance_check.json).

Reproduce from the project root:

```bash
.venv/bin/python work/data_audit2/check_mock_provenance.py
```

The checker is CPU-only and writes only `work/data_audit2/provenance_check.json`.

## Matched interacting-vacuum null-search contract

Reuse these exact 1,200 full-covariance vectors under the fixed fitted flat-ΛCDM null. For every realization, compute

```text
Δχ² = χ²(flat LCDM, refit) - χ²(interacting vacuum, refit)
```

with the same 13 observations and covariance, profiled α=α(Ωm,Γ/H₀), optimizer starts, component-positivity rule through (z=2.33), success/failure handling, and retry policy as for the observed fit. Keep per-mock parameters, scores, optimizer/domain statuses and checkpoints. The calibration then describes the false-improvement frequency for that exact finite search; it is not a posterior probability or a complete correction across other extension families.

The domain question is resolved: freeze the primary mock search at Γ/H₀ ∈ [−3, 3]. This is the source paper's Table-I box and was declared in the locked compute contract before the observed optimization; the observed minimum is interior at −0.4667. Nested ±1 and ±0.5 profiles retain the same interior minimum. The ±0.2 box from the initial theory sketch is only a compact exploratory alternative; its independent profile hits the lower boundary, so it does not replace the predeclared primary domain.

The 1,200-mock calibration still applies only to this exact interacting search. If the final campaign claims a winner among CPL, curvature and interaction, repeat the full across-model selection rule on each null realization; otherwise label the interacting bootstrap as lane-specific.

Use 20 CPU worker processes (the available affinity is 20 CPUs) with BLAS/OpenMP/MKL/NumExpr fixed to one thread per worker, `chunksize=1`, and checkpoint every 25 mocks. Do a 40-mock timing tranche from the existing stream, persist its actual final fit records, project runtime from observed throughput and resume the remaining 1,160. Current timings do not support an honest precise wall estimate: the observed interaction profile including a 41×41 grid and 39 starts took 12.10 s on one CPU; the 1,000-mock CPL selection audit took 711.57 s on 20 CPUs; and the simpler 1,200-mock curvature audit took 6.36 s on 20 CPUs. The interaction mock timing depends on its ODE/viability rejects and convergence behavior. An adaptive, low-dimensional ODE fit has no demonstrated GPU advantage, so keep this tranche CPU-only unless an end-to-end benchmark shows otherwise.

The selected interaction is a **late-time background screen**. The current best-fit Γ is negative, meaning DE → CDM under the paper's (Q=\Gamma\rho_x) convention. BAO does not constrain the baryon/CDM split, (H_0) or (r_d) separately; positive component feasibility is not a measured baryon abundance. The null search will not validate perturbations, growth, radiation-era evolution or a full CMB-calibrated interacting cosmology.

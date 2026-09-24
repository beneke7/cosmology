# RTX 5090 DESI BAO profile-surface benchmark

Status: reproduced float64 CPU/GPU benchmark; exploratory profile surface. The calculation profiles a three-dimensional CPL BAO likelihood at 1,048,576 scrambled Sobol points. It uses the local 13-row DESI DR2 mean vector and its full covariance, rather than evaluating a prediction alone.

## Calculation

Each Sobol point has θ = (Ωm, w0, wa), mapped uniformly over these screening domains:

| Parameter | Domain |
| --- | --- |
| Ωm | [0.05, 0.6] |
| w0 | [-2.0, -0.3] |
| wa | [-3.0, 3.0] |
| α = c/(H0 rd), profiled | [1e-6, 1e4] |

The shape prediction g(θ) uses the repository's flat, late-time CPL expansion history and 96-point Gauss–Legendre integration to each of the 13 observed redshifts. Radiation and curvature are omitted as in the existing screening kernel. With C = LLᵀ, the code forms gw = L⁻¹g and yw = L⁻¹y, computes α̂ = clip(gwᵀyw / gwᵀgw) inside the stated interval, and evaluates χ² = ||yw − α̂gw||². This is an analytic least-squares profile over the distance-scale nuisance parameter with the full covariance. It is not a posterior, Bayesian evidence calculation, or calibrated significance test; no H0 or rd constraint follows from the free amplitude.

Seed 20260924 and the exact row order, bounds, hashes, code versions, command arguments, and per-batch outputs are recorded in [result.json](/home/v/proj/bene/cosmology/cosmology_autoresearch/experiments/gpu_bao_surface/result.json). The saved compressed surface has 1,048,576 parameter triples plus CPU/GPU α and χ² arrays in [surface.npz](/home/v/proj/bene/cosmology/cosmology_autoresearch/experiments/gpu_bao_surface/surface.npz).

## Correctness and timing

The double-precision NumPy path agrees with the repository's `predict_bao()` API to 1.43e-14 maximum absolute prediction difference. For all 1,048,576 points, maximum GPU-versus-CPU differences were 3.56e-14 in profiled α and 1.10e-11 in χ² (2.90e-14 maximum relative χ² difference). The 128-point scalar check integrates each distance with adaptive SciPy quadrature (`epsabs=epsrel=1e-13`) and computes residualᵀ C⁻¹ residual directly: its largest prediction difference from the 96-point kernel was 4.98e-14, and its largest χ² difference from the float64 whitened path was 4.55e-12. The covariance Cholesky succeeded; its condition number was 116.82. No sample hit an α bound.

Two fresh-process 2²⁰ runs, each charging CuPy import and CUDA context initialization to its GPU time, gave:

| Run | CPU evaluator | GPU end-to-end evaluator | Shared input/Sobol setup | Full-pipeline speedup | Peak CuPy pool |
| --- | ---: | ---: | ---: | ---: | ---: |
| Repeat | 16.837 s | 0.5149 s | 0.0554 s | 29.62× | 246,004,736 B |
| Final recorded run | 16.735 s | 0.5189 s | 0.0559 s | 29.21× | 246,004,736 B |

The full-pipeline ratio is (CPU evaluator + shared setup) / (GPU end-to-end evaluator + shared setup). GPU timing includes CUDA context/device setup, host-to-device transfer of the covariance transform and each parameter chunk, the same quadrature and full-covariance profile, output transfer, and explicit synchronization. The matching CPU/GPU chunks were 4,096 points. CPU process time on the final run was 16.731 s over 16.735 s wall time (about one CPU core; BLAS/OpenMP thread caps were set to one). GPU evaluator process CPU time was 0.478 s over 0.519 s wall time. Including shared setup, the final full-pipeline rates were 62.45k points/s on CPU and 1.824M points/s on GPU.

The RTX 5090 had 5,003 MiB free before the run, with the other existing GPU services left untouched. CuPy reported 4,499 MiB free after CUDA context creation and a 246,004,736-byte (about 235 MiB) high-water allocation from its memory pool; the process pool was capped at 2 GiB. The maximum measured process allocation fit well within the initially free memory.

Environment: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, CuPy 14.2.0; RTX 5090; NVIDIA driver 580.159.03; CuPy CUDA runtime 12.9 (runtime API 12090), with system CUDA toolkit 12.8. Compute worker model/effort: gpt-6-luna/max.

## Reproduction

From the project root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
.venv/bin/python scripts/run_bounded.py --seconds 900 -- \
.venv/bin/python experiments/gpu_bao_surface/benchmark_bao_surface.py \
  --powers 20 --chunk-size 4096 --adaptive-check-points 128 \
  --pool-limit-gib 2 \
  --output experiments/gpu_bao_surface/result.json \
  --surface experiments/gpu_bao_surface/surface.npz
```

The final source hash is `42165f484b3e261b95f82402d5be1c050e7653e5b342a52628449cbf0f747b37`; the input mean and covariance hashes are `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` and `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509`. The surface SHA-256 is `fba7da5ae29c4fbf362e124e18e6d6eba45629cd11a5bc9cec1544988c83a4af`.

The sampled minimum was χ² = 5.74155 at (Ωm, w0, wa, α) = (0.37108, -0.31614, -2.23586, 32.33043). This is only the lowest point among the declared Sobol samples. The benchmark supports a material end-to-end speedup for this dense profile-surface workload; it does not establish a change to cosmology or a better calibrated fit.

## Attempt history

The first pilot invocation stopped before the GPU evaluation because a benchmark self-check unpacked the three-value profile-kernel return in the wrong order and mistook profiled α for the prediction vector. The failure was caught by the repository-API assertion (`max_abs=874.57`), the unpacking was corrected, and the pilot then passed. The corrected pilot (`pilot.json`) and an independent 2²⁰ timing repeat (`repeat.json`, `repeat.npz`) remain in the experiment directory. No lower-precision GPU run was performed; the CPU float64 reference preceded all CUDA timing.

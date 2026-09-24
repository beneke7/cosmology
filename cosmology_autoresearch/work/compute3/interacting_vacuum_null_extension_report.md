# Interacting-vacuum null calibration: completed extension

## Disposition and interpretation

The authorized extension completed all 1,160 remaining realizations: seed `20260924`, indices 20–199 (180 rows), and seed `20260925`, indices 20–999 (980 rows). Together with the unchanged 40-mock timing tranche, the matched calibration now has 1,200 unique vectors and fits. This is a finite parametric-bootstrap calibration of the profile-search improvement under the saved fitted flat-ΛCDM null, not a posterior p-value, model evidence, or discovery significance.

The observed profile statistic is `Δχ² = 1.5171630696615157`. The pooled raw mock statistics exceed or equal it in 256/1,200 cases (`0.21333`; exact two-sided Clopper–Pearson 95% interval `[0.19046, 0.23762]`). The pooled 95th percentile is `3.78610` (NumPy linear quantile); median `0.47047`, range `[2.3343e-7, 7.48847]`. These characterize this specific compressed-BAO background profile procedure and search budget only.

## Frozen search and integrity checks

The runner uses the locked paper Table-I box `Ωm ∈ [0.05, 0.60]`, `g=Γ/H0 ∈ [-3, 3]`, and free profiled `α ∈ [1e-6, 1e4]`. Each mock receives 16 deterministic flat-ΛCDM starts and 39 interaction starts plus the 41×41 grid/refinement, using the exact 13-row mean/covariance and the same physical-density feasibility test through `z=2.33`. All SciPy results/statuses are retained separately from re-evaluated physical endpoint validity; the exact nested flat fit is available as a selected-score fallback without overwriting the raw interaction-search result.

Extension integrity checks found 1,160 expected, complete, unique rows, 1,160 attempt records, no latest failed records, and zero fit exceptions. All stored mock-vector shapes, recomputed little-endian float64 SHA-256 digests, and exact regenerated PCG64 seed/index values passed. The source normal-stream digests and input hashes matched the passing `work/data_audit2` provenance audit. The extension began without a prior checkpoint; `jsonl_ahead_of_prior_checkpoint_on_resume=false`. Pilot JSONL, result, and checkpoint hashes remained unchanged.

## Optimizer/status audit

Counts below classify raw SciPy status and endpoint validity independently. A raw optimizer success at the `1e80` penalty plateau is not a valid fit.

| Search component | Extension (1,160 mocks) | Pooled (1,200 mocks) |
|---|---:|---:|
| Flat starts | 18,560 total: 17,844 status 0/success; 716 status 2/failure; all endpoints physically valid; 0 invalid penalties | 19,200 total: 18,458 status 0/success; 742 status 2/failure; all endpoints physically valid; 0 invalid penalties |
| Interaction multistarts | 45,240 total: 45,155 status 0/success, 85 status 2/failure; 31,320 valid endpoints, 13,920 invalid-penalty endpoints | 46,800 total: 46,712 status 0/success, 88 status 2/failure; 32,400 valid endpoints, 14,400 invalid-penalty endpoints |
| Invalid interaction endpoints | All 13,920 were status 0/success at `1e80`; none accepted as a fit | All 14,400 were status 0/success at `1e80`; none accepted as a fit |
| 41×41 grid evaluations | 1,382,640 valid, 564,920 invalid (520,840 `E²` failures; 44,080 nonpositive-state failures) | 1,432,800 valid, 584,400 invalid (538,800 `E²` failures; 45,600 nonpositive-state failures) |
| Grid refinement | 1,082 status 0/success; 78 status 2/failure; all 1,160 endpoints physically valid | 1,120 status 0/success; 80 status 2/failure; all 1,200 endpoints physically valid |
| Nested-flat fallback | 0 | 0 |

The extension grid contains 1,947,560 evaluations in total (`1,382,640 + 564,920`). All raw fit outputs/statuses remain in the JSONL; failed optimizer statuses are retained rather than silently discarded.

## Runtime and reproduction

The extension ran CPU-only with 20 worker processes, 20 CPUs in affinity, and one thread each for BLAS/OpenMP/MKL/NumExpr. It completed normally in `892.4302 s` wall time (`17,490.2974` summed worker CPU-seconds; about `1.300` extension mocks/s). Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, Linux x86_64; no GPU used. Combined pilot-plus-extension measured wall time was `928.2580 s` and summed worker CPU time `18,060.4542 s`, about `1.293` mocks/s pooled.

Exact bounded command:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python scripts/run_bounded.py --seconds 3600 -- \
  .venv/bin/python experiments/interacting_vacuum_screen/interacting_vacuum_null_extension.py \
  --workers 20 --checkpoint-every 25
```

## SHA-256 record

| Artifact/input | SHA-256 |
|---|---|
| Frozen reviewed extension runner | `26cdb4216af9e8bfed219115d5b35de86f03a806c6a76601693c561070e00534` |
| Extension rows `null_extension_1160.jsonl` | `caf01a78e711ede965f98b5e389e28affe74e1e292e888c1483fc8fcc85feed1` |
| Extension summary `null_extension_1160_summary.json` | `8bcc93486846237f020b3178403518a8d7f668da5761c0bf7d4c4e83ec8c287a` |
| Extension checkpoint `null_extension_1160_checkpoint.json` | `55b616637a08c0cf47e2559969dc15beedb34359e4f81f50ce8b44dd3e9d2a1d` |
| Unchanged pilot rows | `a5138b807ff21baf776e389848177db617e08b6f59d9318b3ff1e421f9dab09e` |
| Unchanged pilot summary | `402ed8fe5c9864be6629aa0022f7e666d3a919674b8f75b3333951363cb27b97` |
| Unchanged pilot checkpoint | `3f283ee1083f8353ebb8c00b8ae63695702174704afef8be97946d6bee5961ca` |
| Verified provenance audit | `8c4fed3713b9670c37133bfef13ca1a315acbc7d7880e23547f35673fdd7f581` |
| DESI mean vector | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| DESI full covariance | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| Seed 20260924 input result | `0f477285029024fa359740e726b46f2abf8c87ea8704f488d8ad6658c3fc2a3e` |
| Seed 20260925 input result | `5063ae50f267477b41febd30eabe00c3167b0ccb2d826cd94f3fc966dc025693` |
| Fitted flat-null result | `92662f3ee42500d2be34d9f87b485fe56c82e9550fa674e42f4dec2bf2cab1f0` |
| Profile fitter | `ae99ae789aa02885667e27aa2cbbaa0b035323e69d1a4bf7bbd3848949fc5282` |
| Common timing-tranche helper | `35f8edf052a7d76eaa63ea1e2390d79123c1111e7c605f8b2b49973f5b79d5d4` |

Complete PCG64 normal-array digests, independently checked against the audit, are seed `20260924`: `a5265eba65345219facd2ce66e1edf653f7a2a872f51a79905db228fe1f3de21`; seed `20260925`: `b9cf266dee55dd1f67366c030273a27ae1fb107691baf8bfc3eab07fa407c985`.

This remains a late-time background-only, finite profile-search calibration. It does not include early-time calibration/perturbations, establish a posterior or evidence, or account for selection among other tested models. The DR2 Lyα full-shape data overlap this BAO vector and cannot be appended without a joint covariance.

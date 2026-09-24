# Interacting-vacuum null-search timing tranche

## Disposition

Completed the requested 40-mock timing tranche only: indices `0–19` from each of the already-audited PCG64 streams with seeds `20260924` and `20260925`. All 40 rows completed, with no fit-level exceptions. This is a throughput/status check, not a stable estimate of the null tail and not the remaining 1,160 realizations. They were not launched pending root review.

## Locked search and mock reconstruction

The null mean is the saved best-fit flat-LCDM prediction from `experiments/campaign_seed_bao/result.json` (`Omega_m=0.297461814179106`, `alpha=29.524633230788915`). For each seed, the script regenerated the complete `standard_normal((N,13))` draw array using `np.random.default_rng(seed)` (NumPy PCG64), then selected its first 20 rows and formed `y_i = mean + normal_i @ L.T` using the unchanged lower Cholesky factor of the full ordered DESI covariance. Both complete-array little-endian float64 normal digests match the passing `work/data_audit2/provenance_check.json`; the fitted-mean independent prediction recheck differs by at most `3.55e-14`. JSONL rows retain each 13-vector and its SHA-256.

The calibration search uses the pre-fit-locked paper Table-I box `Gamma/H0 ∈ [−3,3]`, `Omega_m ∈ [0.05,0.60]`, and free profiled `alpha ∈ [1e−6,1e4]`, with the exact 13-row full covariance and the interacting-background positive-density feasibility mask through `z=2.33`. The earlier theory sketch `Gamma/H0 ∈ [−0.2,0.2]` is a distinct exploratory, boundary-limited range and was not used here; it excludes the observed best-fit `g≈−0.4667` and cannot substitute for the wider observed-data search in this matched calibration.

For every mock, the flat reference was independently refit with 16 deterministic L-BFGS-B starts. The interaction refit used 39 deterministic starts plus the same 41×41 `(Omega_m,Gamma/H0)` grid and bounded L-BFGS-B refinement used for the observed fit. All optimizer coordinates were re-evaluated under the physical objective. Raw SciPy `success/status/fun` remain separate from finite physically valid endpoint flags and scores. The selected interaction score is the best accepted interaction endpoint/refinement; if it had exceeded the independently refitted flat score, the exact nested flat fit at `Gamma/H0=0` would have been used as a fallback while the raw interaction-search result remained separately saved. No fallback was needed in this tranche.

## Tranche results

| Check | Result across 40 mocks |
|---|---:|
| Flat starts | 640 total; 614 raw SciPy success, 26 status-2 failures; all 640 returned endpoints valid |
| Interaction multistarts | 1,560 total; 1,557 raw success, 3 status-2 failures (all 3 endpoints physically valid); 1,080 valid endpoints and 480 invalid-penalty endpoints |
| Invalid interaction endpoints | All 480 had raw `success=True`, `status=0` at the `1e80` invalid penalty; they are not accepted fits |
| Grid evaluations | 67,240 total; 47,760 valid and 19,480 invalid (449 `E^2` failures and 38 nonpositive-state failures per grid) |
| Grid refinement | 40/40 returned physically valid endpoints; 38 status-0 successes and 2 status-2 abnormal terminations, all retained |
| Exact flat fallback | 0 uses |
| Mock profile runs | 40 complete; 0 exceptions |

The selected profile improvement `Delta chi2 = chi2_flat − chi2_interaction` was positive in all 40; the nested-flat fallback was available but was not needed in any row. Six of the 40 selected deltas were at least the observed profile improvement `1.51716306966`; the tranche median was `0.196362`, maximum `3.143359`. The pilot's exact two-sided 95% Clopper–Pearson interval for 6/40 is `[0.0571, 0.2984]`; both the count and interval are timing-tranche diagnostics only, not the completed calibration's tail result. They must not be read as a stable tail estimate, posterior p-value, evidence, or discovery significance. Full per-realization predictions, fit parameters, starts, statuses, invalid endpoint errors, grid-refinement statuses, raw improvements and fallback decisions are in the JSONL.

## Runtime and bounded command

The 20-process CPU run took `35.83 s` wall and summed `570.16 CPU-s` across workers, with 20 CPUs in affinity. That is an observed throughput of about `1.116` mocks/s; a straight-line rate calculation for the remaining 1,160 would be about `1,039 s` (`17.3 min`) of wall time. This is only a planning projection from a small tranche, not a guarantee; no further fit was started. BLAS, OpenMP, MKL and NumExpr were each fixed to one thread per worker. No GPU was used. Runtime environment was Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, Linux x86_64.

Exact command:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  .venv/bin/python scripts/run_bounded.py --seconds 1800 -- \
  .venv/bin/python experiments/interacting_vacuum_screen/interacting_vacuum_null_timing.py \
  --workers 20 --checkpoint-every 25
```

The JSONL stream was flushed and `fsync`ed after every completed row; checkpoint snapshots were written at 25 and 40 completed rows. The final checkpoint reports 40/40. The bounded command finished normally within its 1,800-second ceiling.

## Reproducibility and hashes

| Artifact/input | SHA-256 |
|---|---|
| Tranche driver `experiments/interacting_vacuum_screen/interacting_vacuum_null_timing.py` | `35f8edf052a7d76eaa63ea1e2390d79123c1111e7c605f8b2b49973f5b79d5d4` |
| Existing fitter `experiments/interacting_vacuum_screen/interacting_vacuum_profile.py` | `ae99ae789aa02885667e27aa2cbbaa0b035323e69d1a4bf7bbd3848949fc5282` |
| Full row file `experiments/interacting_vacuum_screen/null_timing_40.jsonl` | `a5138b807ff21baf776e389848177db617e08b6f59d9318b3ff1e421f9dab09e` |
| Result summary `experiments/interacting_vacuum_screen/null_timing_40.json` | `402ed8fe5c9864be6629aa0022f7e666d3a919674b8f75b3333951363cb27b97` |
| Final checkpoint `experiments/interacting_vacuum_screen/null_timing_40_checkpoint.json` | `3f283ee1083f8353ebb8c00b8ae63695702174704afef8be97946d6bee5961ca` |
| Verified stream audit `work/data_audit2/provenance_check.json` | `8c4fed3713b9670c37133bfef13ca1a315acbc7d7880e23547f35673fdd7f581` |
| Mean vector `context/data/desi_dr2_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| Full covariance `context/data/desi_dr2_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| Seed-20260924 source result `experiments/bao_robustness/result.json` | `0f477285029024fa359740e726b46f2abf8c87ea8704f488d8ad6658c3fc2a3e` |
| Seed-20260925 source result `experiments/bao_robustness/null_extension_1000_seed20260925.json` | `5063ae50f267477b41febd30eabe00c3167b0ccb2d826cd94f3fc966dc025693` |
| Fitted-flat null result `experiments/campaign_seed_bao/result.json` | `92662f3ee42500d2be34d9f87b485fe56c82e9550fa674e42f4dec2bf2cab1f0` |

Complete-normal stream digests (both verified against the audit): seed `20260924`, `a5265eba65345219facd2ce66e1edf653f7a2a872f51a79905db228fe1f3de21`; seed `20260925`, `b9cf266dee55dd1f67366c030273a27ae1fb107691baf8bfc3eab07fa407c985`.

The calculation is a finite parametric-bootstrap search under a fitted flat-LCDM null using this specific compressed-BAO background model and optimizer/search budget. It does not validate perturbations, radiation-era evolution, baryon calibration, or the full interacting cosmology, and it is not a posterior/model-evidence calculation. Root review is required before any further mock fits.

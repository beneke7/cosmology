# Independent compute33: DESI DR2 BAO optimizer stress test

Status: **exploratory independent numerical cross-check**. This scorer reparses the official 13-row DESI BAO mean vector and full 13×13 covariance, and uses pinned CAMB 2.0.4 without importing compute31. It evaluates flat fluid ΛCDM with positive flat-Λ closure over the declared standard-BBN-motivated interval, fixing omega_b h² slice by slice.

The one-point reproduction at compute31's reported BAO minimum agrees exactly in chi². The search covers 25 fixed omega_b h² values: compute31's 15-point coarse grid, its seven upper-edge refinement nodes, and all three selected prior centers. At each slice, 24 deterministic Powell starts were used (600 total). Five of those slices also received a Sobol-initialized differential-evolution search, with 6,400 evaluations total. There were 47 transient invalid Powell probes and 6 transient invalid differential-evolution probes; all 600 Powell endpoints were valid positive-Λ fits with successful optimizer flags. No start ended in a distinct basin more than Δchi²=10⁻⁶ above that slice's best.

The best independent point is chi²=10.28146188607339 at omega_b h²=0.024, H0=69.993181535, and omega_c h²=0.120954994. Compute31 reported chi²=10.281461886073446 at effectively the same parameters. The difference, −5.51×10⁻¹⁴, is floating-point roundoff. The selected prior-center profiles also reproduce compute31: chi²=10.281958582 at 0.022017, 10.281950004 at 0.02205, and 10.281911076 at 0.0222. This search found no credible lower basin. The sampled profile keeps decreasing toward the upper baryon edge; these finite searches cannot prove globality or locate an unsampled continuous minimum.

Dense covariance solve and Cholesky-whitened chi² scores agree to at most 7.11×10⁻¹⁵ over returned endpoints. Every selected slice has positive CAMB Omega_Lambda; its implied 100 theta_MC ranges from 1.0371 to 1.0457, inside the broad [0.5, 10] support. H(z) and c/H(z) are in km/s/Mpc and Mpc, while D_M and the drag sound horizon are in Mpc; the 13 predicted D_V/r_drag, D_M/r_drag, and D_H/r_drag values and closure checks are recorded in the JSON.

The run took 50.74 seconds on 20 worker processes. Their aggregate CPU time was 677.14 seconds, corresponding to 13.35 average busy cores over the pool wall time. This measures realized parallel CPU use rather than a serial-versus-parallel benchmark. Numerical-library and CAMB threads were capped at one per process; no GPU or downloads were used. The search parameterization spans the full physically feasible part of the declared box; below the closure threshold in H0, the lower omega_b and omega_c bounds make positive flat Λ impossible.

Reproduce from `cosmology_autoresearch/`:

```sh
.venv/bin/python scripts/run_bounded.py --state work/compute33/no_deadline_state.json --seconds 800 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python -B work/compute33/independent_bbn_bao_globality.py --workers 20
```

The executed-run script SHA-256 is `2b3a68fdb333f8f01d3ec5373cfa45f26601f5e032136e1b968ba97237924869`. The delivered script SHA-256 is `c2188bf0bcffd3527c39f23086a9891836707f61a2a6898f0f43faaa1a7fff8c`; a post-run bookkeeping/report-format correction fixed the transient-invalid counters, with no numerical search rerun. Exact input hashes, seeds, software versions, all starts and failed evaluations, and machine-readable results are in [independent_globality.json](independent_globality.json). [profile_globality.png](profile_globality.png) shows the sampled profile.

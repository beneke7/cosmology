# Independent check of the compute26 CAMB BAO boundary candidate

## Result

An independently written one-point calculation reproduces the compute26 candidate's 13-row DESI DR2 full-covariance score exactly at reported precision: dense covariance solve and Cholesky whitening both give \(\chi^2=10.273886272857048\), with zero difference from the compute26 JSON value. Independently integrating \(c/H(z)\) agrees with CAMB's comoving radial distance to at most \(9.1\times10^{-13}\) Mpc.

This validates the score at that parameter point, not the optimizer's global optimality or the candidate's physical plausibility. The point remains pinned to the declared upper \(\omega_bh^2=0.1\) boundary, with \(H_0=133.61\) km s\(^{-1}\) Mpc\(^{-1}\), and comes from BAO alone without a CMB likelihood. It must remain labeled an incomplete boundary profile, not a cosmological constraint or accepted interior baseline.

The apparent edge-slope warning is resolved by the independent reprofile in [compute29](../compute29/boundary_audit_report.md): at fixed \((H_0,\omega_ch^2)\), the score has a minimum at 0.1; after reoptimizing those nuisance coordinates, the in-box profile decreases toward 0.1 and the outside-box profile continues to decrease. Three independent full-box starts recover the same bound-constrained point. Thus the score is supported as the minimum on the tested box, but it is not an accepted interior minimum and is sensitive to extending the broad density limit.

## Independent method and provenance

The standalone script reads the official mean/covariance directly, configures pinned CAMB 2.0.4 at the supplied \((H_0,\omega_bh^2,\omega_ch^2)\), reconstructs \(D_M/r_d\), \(D_H/r_d\), and \(D_V/r_d\), solves the dense covariance system, performs an independent Cholesky whitening check, and integrates \(c/H\). It does not import the compute26 objective/scorer.

Reproduce from the campaign root:

```sh
.venv/bin/python scripts/run_bounded.py --state work/compute28/task_state.json --seconds 600 -- env PYTHONPATH=work/compute24/site OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -B work/compute28/independent_candidate_check.py
```

The run took 0.012 s. Script SHA-256: `c223d66a50bf6c254af134bb8bb4e3181877baeff0918c4d5e16746e98487a4f`; result SHA-256: `7cd8e1dc5dd629aeab4cc0212a0f0b1ea0236801c0738ad228b394608855106d`. Input means/covariance hashes and the compute26 reference hash are preserved in the JSON result.

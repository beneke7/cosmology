# Preparation checks and limits

Checked on 2026-09-23 in the preparation environment. This records completed preparation work; it is not the overnight morning report.

## What ran

- Seventeen automated tests passed: 12 context-ingestion tests and five background/BAO tests. See `runs/tests.log`.
- TG synthetic equation/ODE checks passed, with results in `runs/tg_audit.json`. These check the stated sign convention and transformed equation, not galaxy data or a full gravitational theory.
- CPL reduces correctly to ΛCDM; the distance integral agrees with independent adaptive integration and the analytic Einstein–de Sitter limit; synthetic ΛCDM parameters are recovered by the nested CPL fit. Covariance symmetry, positive definiteness, observable labels and supplied row order are checked.
- Real core downloads succeeded and are recorded in `context/manifest.lock.json`: ten open papers, the paired DESI DR2 BAO mean/covariance, two SPARC archives, and six documentation pages. The supplied journal PDF is included separately. Paper text was extracted, but full-text reading is not asserted for every cached paper.
- The BAO source is pinned to CobayaSampler/bao_data commit `bb0c1c9009dc76d1391300e169e8df38fd1096db`. The whole source collection records content hashes and retrieval metadata. A repeated fetch verified the cache rather than silently replacing it.
- Project/role TOML parsed successfully. The bounded-command wrapper was checked for successful completion, timeout termination and refusal to start after an expired deadline.

## One real-data numerical screen

Replay from the project root:

```bash
python scripts/run_seed_screen.py
```

The 13-row compressed DESI DR2 BAO vector gives the following bounded optimization results with free amplitude c/(H0 r_d), flat late-time geometry, and radiation omitted in the fit:

| Family | Free fitted parameters | Minimum chi-squared in the stated domain | Qualification |
|---|---:|---:|---|
| ΛCDM | 2 | 10.271041 | BAO-only background screen |
| constant-w | 3 | 9.041047 | Same covariance and amplitude freedom |
| CPL | 4 | 5.700616 | Best fit hits the chosen upper w0 bound, −0.3 |

These values establish a working calculation. They do **not** reproduce the collaboration's joint posterior constraints, give Bayesian evidence, establish significance, or demonstrate new physics. The CPL result is explicitly sensitive to the allowed parameter domain; expanding/checking that domain is an obvious follow-up, not an optional cosmetic check. Optimizer success is not proof of a global optimum.

`experiments/seed_bao/result.json` contains parameters, bounds, boundary flags, predictions, hashes, runtime and independent numerical checks. A separately written scalar adaptive integrator agreed with the vectorized predictions to less than 3×10^-14 in this run. An illustrative fixed-parameter radiation sensitivity check with Ωr=9×10^-5 shifted predictions by at most about 0.041 marginal measurement standard deviations. It is not a radiation refit, a complete early-universe model or a guarantee across the whole search domain.

## What remains untested

There was no GPU or Codex executable in the preparation container. The optional CUDA path, actual workstation throughput, target-client parsing/model access, and eight-hour orchestrator lifecycle have not been run here. No job has been started on the user's workstation. The model/configuration syntax was checked against current official documentation, and the package is ready for target-machine startup validation.

There is no completed SN, CMB, shear, SPARC galaxy or combined-probe inference, no posterior/evidence calculation, no adaptive-search mock calibration, and no Lean proof in this preparation. Those belong to the campaign and must be reported according to their actual completion status.

Seed numerical dependencies in this environment were Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. `requirements.txt` permits compatible installation; it is not a reproducibility lock. The target run must record its resolved environment and rerun the relevant checks.

# IVS physicality check — quadrature diagnosis and validated rerun

## Diagnosis at `g=Gamma/H0=2.5`

The initial failure was a Simpson-resolution issue, not an ODE sign/unit defect. Keeping the same DOP853 dense solution and fiducial cosmology fixed, the absolute residual in `q_end - (q0 + g I)` was:

| Uniform nodes | Simpson residual |
| ---: | ---: |
| 2049 | `2.3772067647698236e-8` |
| 4097 | `1.4856749164238181e-9` |
| 8193 | `9.285372470912989e-11` |
| 16385 | `5.803579838925543e-12` |

The original acceptance limit was `2e-8`; the 2049-node residual exceeded it only slightly. Refining the grid reduces the error by approximately a factor of 16 per doubling, as expected for composite Simpson on this smooth integrand. Adaptive `scipy.integrate.quad` gives `I=0.11291496559134644` (reported error estimate `1.23e-15`), with `|q_end - (q0+g I)|=3.33e-16`; this differs from the ODE-integrated `I` by `8.33e-17`. The ODE state identity residual is `1.11e-16` at the endpoint and at most `4.44e-16` on 4097 dense samples. Thus the equation and unit normalization are internally consistent; only the original sampled quadrature resolution failed its strict residual gate.

## Full acceptance rerun

The failed original producer at `../ivs_physicality.py` remains unchanged, and the prior `../physicality_result.json` remains the STOPPED/no-result attempt record (SHA-256 below). `rerun_full_acceptance.py` reruns that producer's full checks while overriding only its default dense-output Simpson grid from 2049 to 8193 nodes. The threshold logic, physical model, and other checks were not changed. Its output is separately stored as `fix/physicality_result.json`.

At the requested fixed flat benchmark (`Omega_m0=0.315`, `Omega_b0=0.049`, `H0=67.4`, `Tcmb=2.7255 K`, `Neff=3.044`, massless neutrinos), `Omega_c0=0.266`, `Omega_gamma0=5.435021880971039e-5`, `Omega_nu0=3.7573103850318476e-5`, and flat closure gives `Omega_x0=0.68490807667734`. The paper source used was the local arXiv v1 PDF, SHA-256 `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`. Table I gives the paper's broad IVS uniform range `g in [-3,3]`; Eqs. (3)–(8) give `E^2=(Omega_gamma0+Omega_nu0)e^(4u)+Omega_b0 e^(3u)+m+x`, `dx/du=-g x/E`, `dm/du=3m+g x/E`, `q=m e^(-3u)`, and `dq/du=g x e^(-3u)/E` for `u=ln(1+z)`. These signs follow the paper's forward-time convention `Gamma>0` transfers CDM to vacuum.

Every existing producer acceptance check passed with the 8193-node default: g=0 LCDM limit, positive `x` and `E^2`, positive CDM for trajectories reaching `z=1e10`, comoving-CDM identity (ODE and Simpson), `q` monotonicity, and bracket ordering. The maximum Simpson residual over the full 25-point g grid was `1.4206891219004092e-10`, below the unchanged `2e-8` limit; the maximum ODE-state identity residual was `1.3322676295501878e-15`. Runtime was 2.39 s on one CPU core under the 600 s bounded runner.

The benchmark-specific first-zero boundary by `z=1e10` is bracketed by `g=-1.0225570347975008` on the crossing side and `g=-1.0225570346228778` on the side still physical through the endpoint (bracket width `1.7462298274040222e-10`). The lower bracket endpoint's first zero is at `z=835.5507`; this is the crossing redshift for that endpoint, not a claim that the boundary itself is crossed at that redshift. The upper endpoint has `q(z=1e10)=7.0566e-11`. On the 0.25-spaced grid, 8 sampled negative-g values cross before the endpoint.

This only establishes that, holding this illustrative late-time benchmark fixed, the paper's broad numerical `g` interval contains values whose background CDM reaches zero and would become negative if extrapolated farther. It is **not** a statement about the paper's prior mass after varying its sampled physical densities and `100 theta_MC`, and says nothing about its posterior, evidence, or published analysis. No trajectory is continued after its first zero.

## Reproduction commands and integrity

Both runs used `scripts/run_bounded.py --state work/compute19/RUN_STATE.json --seconds 600` with the task-local state file absent and thread caps `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`. The diagnostic and rerun commands were:

```text
python3 scripts/run_bounded.py --state work/compute19/RUN_STATE.json --seconds 600 -- env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python work/compute19/fix/g2p5_quadrature_diagnostic.py
python3 scripts/run_bounded.py --state work/compute19/RUN_STATE.json --seconds 600 -- env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python work/compute19/fix/rerun_full_acceptance.py
```

The initial system `python3` attempt in the parent failure record lacked SciPy and stopped before numerical work. The successful runs used `.venv/bin/python` 3.12.3. Hashes for the preserved failed attempt, diagnostic, rerun driver, and validated output are in `fix/hashes.sha256`.

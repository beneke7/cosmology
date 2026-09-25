# Independent Dovekie IVS profile reproduction

Status: **passed**. This standalone checker parsed the SNANA release order and unpacked the full STAT+SYS precision itself. It imports none of the production profile, SN helper, ODE, distance, or scoring functions.

The local ODE uses the source convention `Q=Gamma*rho_x>0` (CDM to vacuum), with `g=Gamma/H0`, `Gamma/H=g/E`, `dx/du=-g*x/E`, `dm/du=3*m+g*x/E`, and `dD/du=exp(u)/E`. Distances use `zHD` in the integral and `zHEL` in the luminosity prefactor at a fixed 70 km/s/Mpc gauge. One unbounded intercept is profiled by GLS over the full 1,820-row precision; no MUERR term is added.

Input CID-order fingerprint: `b7d5c6ad8dfdadf1e12443852b14006c0d9d0ceb4bcb46efb378f058ce5aa5a3`. Both pinned data files and all code/result identities are in the JSON report. Full precision Cholesky passed; the independent intercept normal equation and direct/whitened/Schur quadratic evaluations are recorded per point.

| Point | Ωm0 | g | Independent χ²profile | Δχ² from IVS minimum | Abs. χ² difference from production |
|---|---:|---:|---:|---:|---:|
| ivs_minimum | 0.39301742027 | -0.629387814517 | 1629.999216924071 | 0.000000000000 | 4.55e-12 |
| bao_ivs_reference | 0.386971007697 | -0.466664729917 | 1630.598199134153 | 0.598982210083 | 2.5e-12 |
| sn_lcdm_reference | 0.330316790548 | 0 | 1631.420535571243 | 1.421318647172 | 2.05e-12 |
| bao_lcdm_reference | 0.297461814954 | 0 | 1636.296075323322 | 6.296858399252 | — |

The analytic flat-LCDM distance limit at g=0 agrees with the independent IVS ODE to 7.11e-14 mag (tolerance 2e-08); the source-term direction and standard-vs-tight RK45 refinements pass. Root-reported score and delta discrepancies are below 0.0001.

Physicality was checked over the observed 0 ≤ zHD ≤ 1.14418 range using 4,001 trajectory samples, with positive m, x, E² and an interior constant baryon/CDM split witness. This is a dense numerical check, not a formal continuous positivity proof. The reported finite-search minimum is reproduced at its coordinate; this check does not establish global optimality or turn the fixed BAO references into a joint likelihood.

Exact command from `/home/v/proj/bene/cosmology/cosmology_autoresearch`: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python work/compute10/dovekie_ivs_profile_independent.py`

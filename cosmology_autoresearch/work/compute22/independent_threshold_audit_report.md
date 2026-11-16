# Independent audit of the compute19 physicality threshold

**Disposition: independently reproduced for the stated fixed benchmark.** An RK45 solve using the independent state `(x,q)` reproduces the selected coupling outcomes and brackets the reported threshold. It does not establish posterior support or reproduce the paper's sampled prior.

## Source and independent equations

The primary source is the local arXiv v1 PDF, SHA-256 `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`. Its flat Friedmann equation (3), conservation equations (4)–(5), and `Q=Gamma*rho_x` in Eq. (8) imply, for `w_x=-1` and `n=ln(a)`:

```text
E^2 = Omega_rad*exp(-4n) + Omega_b*exp(-3n) + rho_c + x
x_n = g*x/E
(rho_c)_n = -3*rho_c - g*x/E
g = Gamma/H0
```

I used `x=rho_x/rho_crit0` and the independent comoving state `q=exp(3n)*rho_c`. Then:

```text
E^2 = Omega_rad*exp(-4n) + Omega_b*exp(-3n) + q*exp(-3n) + x
x_n = g*x/E
q_n = -g*exp(3n)*x/E
q(n) = Omega_c0 + g*integral_n^0[exp(3v)*x(v)/E(v)]dv
```

This is the compute19 comoving identity expressed in `n`; with `u=-n=ln(1+z)`, it becomes `q(u)=Omega_c0+g*integral_0^u[x/E*exp(-3v)]dv`. The sign agrees with the paper: positive `Gamma` transfers CDM energy to vacuum forward in time. For `g<0`, `q` decreases when integrating backward. I stop integration at the terminal `q=0` event and do not continue into negative CDM.

## Benchmark and selected couplings

I preserved the fix report's stated benchmark convention: `Omega_gamma*h^2 = 2.469e-5*(Tcmb/2.7255 K)^4`, `h=0.674`, `Neff=3.044` with massless neutrinos, `Omega_m0=0.315`, `Omega_b0=0.049`, flat closure. This yields `Omega_c0=0.266`, `Omega_gamma0=5.435021880971039e-5`, `Omega_nu0=3.7573103850318476e-5`, and `Omega_x0=0.68490807667734`, matching the saved result. This is the benchmark's declared radiation convention; it is not a recalibration to a different photon-density constant.

The independent solver used RK45, `rtol=3e-13`, `atol=2e-15`, `max_step=0.02` in `n`, and a terminal `q=0` event. It sampled positivity and monotonicity at 4,097 points through the event or endpoint, and checked the comoving identity independently with adaptive quadrature of dense output.

| `g` | Stop | `z` at stop | `x` at stop | `q` at stop | Identity absolute residual |
| ---: | --- | ---: | ---: | ---: | ---: |
| −3 | `q=0` | 0.142358005541392 | 1.01129594945626 | 7.37e-18 | 2.66e-15 |
| −1.25 | `q=0` | 0.704881639003588 | 1.24638803449283 | 4.00e-17 | 5.95e-16 |
| −1 | `zmax` | 9,999,999,999.999987 | 2.51021129973316 | 0.00895702737249135 | 1.08e-15 |
| 0 | `zmax` | 9,999,999,999.999987 | 0.68490807667734 | 0.266 | 0 |
| +3 | `zmax` | 9,999,999,999.999987 | 0.0631989825282942 | 0.580538369758141 | 5.06e-14 |

All sampled `x`, `E^2`, and CDM density values are positive through each stop. The event cases stop at numerical zero (`q` residual below `4.0e-17`); no post-event negative branch is evaluated. The endpoint mapping is `n_min=-ln(1+1e10)=-23.025850930040455`, `u_max=23.025850930040455`, `z=expm1(u)=9,999,999,999.999987` at double precision. The largest scaled comoving-identity residual among these selected cases is `8.73e-14` at `g=3`.

## Threshold comparison and event limitation

Starting with an independent event/no-event bracket `[-1.25,-1]`, 35 bisection steps gave:

- Crossing side: `g=-1.0225570347975008`, identical to the saved crossing-side endpoint. My event solution reaches `q=0` at `u=6.735254590763463`, `z=840.5577096358431`.
- Still positive through `zmax`: `g=-1.0225570347902249`, with `q(zmax)=2.8355040002502004e-12`.
- Independent bracket width: `7.275957614183426e-12`; the bracket is contained in the saved `[-1.0225570347975008,-1.0225570346228778]` bracket of width `1.7462298274040222e-10`.

The saved upper endpoint has `q(zmax)=7.056604622148212e-11`, also positive. Thus the event classification and benchmark-specific `g` threshold agree. There is one numerical limitation: the event redshift at the shared crossing-side `g` differs from the saved `z=835.5507089762066` (`delta u=0.00596745`). Near that event, `|dq/du|=5.28e-13`; my independent comoving-integral residual is `4.39e-15`, which corresponds to a local `u` sensitivity scale of about `0.0083`. The event is nearly tangent to zero, so the parameter-side classification is reliable at the bracket scale, while the quoted crossing redshift is not resolved to the displayed precision by this double-precision event calculation.

This threshold says only that this fixed benchmark has negative-`g` values in the broad interval whose backward CDM density reaches zero by `z=1e10`. The paper samples physical densities and `100 theta_MC`; this fixed-`Omega_m0` benchmark is not an exact mapping of those priors. No inference about posterior support, prior mass, likelihood, or evidence follows. The paper fixes `sum m_nu=0.06 eV`, while this benchmark uses massless neutrinos.

## Reproduction and integrity

Successful command, from the project root, with one CPU and a 600-second cap:

```sh
taskset -c 0 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/run_bounded.py --state work/compute22/no_campaign_deadline.json --seconds 600 -- .venv/bin/python work/compute22/independent_threshold_audit.py
```

The run took `0.8694 s` inside the script. Environment: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1; CPU affinity `[0]`; no GPU. Runtime model metadata exposed GPT-6 family only; exact variant and effort tier were not exposed, so the inherited default was used.

- Independent implementation SHA-256: `31190abef863091f55d1075c028c00eeaa259071d27a7f933622ff7edd69baa2`
- Independent result JSON SHA-256: `674d411a2950d94c8bdf094db9f0bf48106fc7b4a83a87c3860db137fbf0e715`
- Compared compute19 result JSON SHA-256: `cc12b1aa1f1e3b0a5fe5352497f646eb19c2c03cd6e90cedf3dc38ae55c2c6bc`
- Compared compute19 fix report SHA-256: `71802cdf1bbcbe5b6bf20c934b8e0f4686473fc01ca6a6857cb307cdb046ce5f`
- Primary PDF SHA-256: `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`

The independent implementation and machine-readable per-case comparison are [independent_threshold_audit.py](independent_threshold_audit.py) and [independent_threshold_audit.json](independent_threshold_audit.json). No compute19 producer or driver was imported or called. No failed numerical attempt occurred.

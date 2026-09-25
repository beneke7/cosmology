# PRIMAT standard-BBN forward feasibility

**Status:** forward feasibility and repeatability reproduced. This is not a reproduction of Poulin et al. figures, an abundance-likelihood fit, or a BBN posterior.

## Preflight and scope

- Public paper: Poulin et al., arXiv:2607.20635v1, analysis §III.1; it identifies Python PRIMAT 0.3.2 and the generic `primat_tools` Cobaya wrapper. It states the paper-specific code, data and notebooks are to be released separately.
- Official PRIMAT upstream: `https://github.com/CyrilPitrou/primat`, tag `v0.3.2`, commit `21ff8f39fa18e3937e9fdf386cfa982361bfdfce` (verified with `git ls-remote`). GitHub tree reports 21,580,336 bytes for `primat/`; full tag tree is 40,306,780 bytes before compression. Repository metadata reports 73,934 KiB. GPL-3.0-or-later per upstream package metadata; license text to be retained with selected source.
- Official `primat_tools`: `https://github.com/CyrilPitrou/primat_tools`, public `master` HEAD `f51331aad86fa62ccc912d331d9a44dde149db35` at preflight. Not needed for a direct PRIMAT forward solve and not fetched.
- Existing toolchain: project Python 3.12.3 virtualenv has NumPy, SciPy and Matplotlib; GCC, G++, CMake and Make exist. `gfortran`, Clang, and Flang are absent. The Python PRIMAT route does not require Fortran; do not install dependencies or alter global/project environments.
- Bootstrap accounting before checkout: `context/manifest.lock.json` reported 103,043,461 / 209,715,200 bytes (101.73 MiB headroom). The sparse checkout's Git object pack is 5,983 KiB; after it, estimated cap headroom is about 95.8 MiB. The checkout occupies 30 MiB on disk, including Git metadata; the versioned `primat/` tree is 21,580,336 bytes. Filesystem free space was 278 GB; the whole project tree was about 1,023 MB. These are separate limits.
- One initial cone-mode sparse-checkout command rejected file patterns (`README.md` is not a directory); it transferred no extra data and was immediately corrected using `--no-cone`. No solver or build command failed.

## Run contract

## Reproduced forward grid

The public source ran without an install or build using the existing project Python environment. PRIMAT's pure-Python backend was forced explicitly. The standard background used no extra energy component, `DeltaNeff=0`, and no neutrino chemical potential. The code's standard neutrino-decoupling calculation returns `Neff=3.0439772985579183` (the upstream SM value is conventionally quoted as 3.044). The five 0.001-MeV endpoint solves used the upstream `small` network (8 light nuclides, 12 tabulated thermonuclear reactions plus n-p weak rates) and central PRIMAT v0.3.2 rates.

| `Omega_b h^2` | `10^5 D/H` | `Yp` | returned `Neff` | runtime (s) |
|---:|---:|---:|---:|---:|
| 0.02150 | 2.6086843 | 0.24659184 | 3.04397730 | 1.957 |
| 0.02175 | 2.5599604 | 0.24670388 | 3.04397730 | 1.967 |
| 0.02205 | 2.5032624 | 0.24683664 | 3.04397730 | 1.951 |
| 0.02225 | 2.4665075 | 0.24692373 | 3.04397730 | 1.960 |
| 0.02250 | 2.4216930 | 0.24703141 | 3.04397730 | 1.964 |

At the broad PDG SBBN prior center, `Omega_b h^2=0.02205 +/- 0.00043` (PDG Eq. 24.6), PRIMAT gives `10^5 D/H=2.5032624`. PDG 2025 reports `10^5 D/H=2.508 +/- 0.029` (Eq. 24.2), so the point prediction is `-0.00474` in those units, or `-0.163` observational standard deviations. This descriptive residual uses only the observational error; it omits PRIMAT rate/lifetime uncertainty and is not used to infer a posterior. The broad PDG omega prior itself comes from SBBN abundance data, so it must not be multiplied by those same abundance measurements as if independent.

Assumptions: `tau_n=878.4 s` fixed (upstream default; `std_tau_n=0.5 s` is not sampled); full upstream default plasma and weak-rate corrections are retained, including QED plasma, non-instantaneous neutrino decoupling/NEVO spectrum, radiative, finite-nucleon-mass, thermal and spectral-distortion corrections. Nuclear rate inputs are the versioned central PRIMAT tables; no nuclear-rate nuisance Monte Carlo, PRIMAT abundance uncertainty, covariance or observational likelihood was evaluated. Abundances are dimensionless ratios; the table reports `10^5 D/H` and the mass fraction `Yp`.

## Reproducibility and exact commands

The source checkout is clean at commit `21ff8f39fa18e3937e9fdf386cfa982361bfdfce`, Git tree `9ae225a109aacfe8534f4e5ff25bf81a40a08326`. Hashing all 460 files under `primat/` gives 21,580,336 bytes and aggregate SHA-256 `583f054b68c100b249126618f41620718ed800f9ddde9d9fa9d1c24aaee98413`; individual file hashes are retained in [source_data_hashes.json](source_data_hashes.json), whose SHA-256 is `95398b44954a81bf4803efcca0a673578d461d958424353689151df0409b91ee`. License text is retained at `repo/LICENCE` (GPL-3.0-or-later). Environment: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1. No compiler was needed; no dependencies were installed.

Initial bounded run, with all thread pools limited to one:

```sh
env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --state work/data9_primat/task_state.json --seconds 540 -- .venv/bin/python -B work/data9_primat/forward_grid.py
```

The five solves took 9.799 seconds total. A second bounded execution repeated all five parameter points under the same pinned source and environment; `D/H`, `Yp`, and returned `Neff` matched bit-for-bit at every point. The repeat record is [reproducibility_check.json](reproducibility_check.json), status `passed`.

```sh
env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python scripts/run_bounded.py --state work/data9_primat/task_state.json --seconds 90 -- .venv/bin/python -B work/data9_primat/forward_repro_check.py
```

## Likelihood-reproduction gate

A standard-SBBN forward model is feasible and reproduced. A generic abundance likelihood could be assembled from the public solver and documented abundance measurements, but the official current `primat_tools` `master` HEAD inspected at preflight (`f51331aad86fa62ccc912d331d9a44dde149db35`) is not pinned by the paper. Its current README describes a Cobaya theory/Gaussian He-4+D/H likelihood wrapper, requires Cobaya >=3.5 and a separate `PyPRIMAT` sibling checkout; Cobaya is absent from this environment and was not installed. The paper itself names only the repository URL. Thus the exact wrapper revision and its compatibility with the paper's PRIMAT v0.3.2 source remain unspecified.

An exact Poulin et al. likelihood/figure reproduction is not feasible from the checked public inputs. The paper says its only PRIMAT change is adding the vEDE energy density to `H(T)`, normalized at `T_D`, but that paper-specific patch is absent from the generic v0.3.2 tag. The linked restricted Zenodo record is not accessed. The exact missing inputs are the unrestricted paper bundle (or equivalent source patch), its `primat_tools` commit and run configuration/YAMLs, and the associated profiles/chains/figure data/notebooks plus the exact PRIMAT rate/lifetime uncertainty summary used in their abundance likelihood. The public article contains enough information for a separately implemented exploratory likelihood, but such an implementation would not be a numerical reproduction of the paper's pipeline.

The restricted Zenodo files for the paper are out of scope and will not be requested or accessed. No posterior, inferred likelihood, or exact vEDE reproduction will be claimed.

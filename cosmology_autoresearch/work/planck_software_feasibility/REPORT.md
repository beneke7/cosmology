# Planck PR3 standard ΛCDM software feasibility

Checked 2026-09-25. This is a static compatibility and data-availability assessment for this workspace. I read the project `AGENTS.md` and `BRIEF.md`, the staged archive report and `inspection.json`, and version-pinned upstream Cobaya/CAMB/clipy sources. I did not install packages, extract the archive, evaluate a likelihood, run a fit, or use the GPU.

## Verdict

The official Planck PR3 baseline archive already staged here contains all three required likelihood data trees. A standard, project-local full-likelihood route is supported by the declared upstream Python/package constraints with this version set:

| Component | Recommended pin | Workspace / upstream compatibility evidence |
|---|---|---|
| Python | `3.12.3` | Existing project environment, as inventoried in `work/data7_planck_gate/planck_2018_full_likelihood_gate.json`. |
| NumPy | `2.5.3` | Existing environment; satisfies SciPy 1.18.1's `>=2.0,<2.8`, Cobaya `>=1.17`, CAMB 1.6.6 `>=1.22`, and clipy-like 0.15 `>=1.23.5`. |
| SciPy | `1.18.1` | Existing environment; upstream metadata requires Python `>=3.12` and NumPy `>=2.0,<2.8`, which admits Python 3.12.3 and NumPy 2.5.3. It also satisfies Cobaya `>=1.5` and CAMB 1.6.6 `>=1.0`. |
| Cobaya | `3.6.2` (PyPI; source tag commit `899f30a49f85de610dac321e91a1af50018e56aa`) | Python `>=3.10`; declares NumPy `>=1.17`, SciPy `>=1.5`; includes these exact Planck likelihood components and runs `clik`-format data through clipy. |
| Py-BOBYQA | `1.5.0` | Cobaya requires `py-bobyqa>=1.4`; 1.5.0 is pure Python, supports Python `>=3.8`, and declares NumPy/SciPy/Pandas dependencies without upper bounds. Pinning this release supplies the Python 3.12-compatible `pyproject.toml` packaging path. |
| CAMB | `1.6.6` | Python `>=3.10`; declares NumPy `>=1.22`, SciPy `>=1.0`; accepted by Cobaya 3.6.2's minimum CAMB version (`1.5.0`). Use its normal PyPI wheel so no Fortran compiler is needed. |
| Planck likelihood implementation | Cobaya 3.6.2's auto-installed `clipy` release `0.15` (`clipy_0.15` upstream tag), or the PyPI `clipy-like==0.15` wheel (SHA-256 `9bf56a1c648cab345e9fcb7f83c9a79cd82f92e4f6fe4392900eec23b96f49d0`) | Pure Python; requires Python `>=3.9`, NumPy `>=1.23.5`, Astropy `>=2.0`; v0.15 specifically changed binary parsing to `frombuffer` for newer NumPy. Cobaya's wrapper sets `CLIPY_NOJAX=1`. |
| Astropy | `8.0.1` | Deliberate exact pin for the clipy dependency: Python `>=3.11`, NumPy `>=2.0`, so compatible with Python 3.12.3 and NumPy 2.5.3. |

This is a **static compatibility pass**, not a reproduced runtime result: no package installation or likelihood evaluation was done, so I cannot certify that the complete transitive environment imports cleanly or that the three likelihoods successfully evaluate together until a small post-install smoke check runs. There is no direct version-constraint conflict in the pinned core set. Generate and retain a resolved package lock after setup; Cobaya's broad dependency ranges do not by themselves freeze its full transitive stack.

Do not combine Cobaya 3.6.2 with CAMB 2.0.4 as the default path. The cached 2.0.4 wheel is Python-installable, but this Cobaya release predates its CAMB v2 compatibility update. Upstream's 2026-07-24 commit `497fa9aebc3dbce2174416645300d2ee7e07f6d0` is explicitly titled “CAMB v2 changes”; it raises Cobaya's minimum CAMB version from 1.5.0 to 2.0.0, changes lensing-accuracy handling, raises its clipy minimum to 0.16, and updates Planck test tolerances. That source commit is a second, non-PyPI option if CAMB 2.0.4 is required, but it is not the stable Cobaya 3.6.2 release. CAMB 2.0 also changed defaults and can change spectra relative to 1.x. Keep the first reproduction on the stable `Cobaya 3.6.2 + CAMB 1.6.6 + clipy 0.15` line; assess the v2 stack separately if needed.

## Likelihood contract

Cobaya 3.6.2's version-pinned component files map the requested combination exactly:

| Requested product | Cobaya component | Planck PR3 archive member |
|---|---|---|
| Plik HM TT+TE+EE | `planck_2018_highl_plik.TTTEEE` | `baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik` |
| Commander low-ℓ TT | `planck_2018_lowl.TT_clik` | `baseline/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik` |
| SimAll low-ℓ EE | `planck_2018_lowl.EE_clik` | `baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik` |

The staged archive is `work/planck_stage/COM_Likelihood_Data-baseline_R3.00.tar.gz`, 60,323,470 bytes, SHA-256 `0b73171e3acc671c28184466a45485a2d1c1d93676b832abdfe688c7b04024e6`. Its checked inventory has 362 members and 92,163,979 regular-file bytes; the three target trees are present. It has not been extracted. The inspected target trees total about 65.1 MB, and the entire regular-file payload is about 92.2 MB, so a project-local extraction under the selected package path fits the already documented cache cap.

Keep the full `TTTEEE` component: `TTTEEE_lite` is a nuisance-marginalized alternative, not the requested full Plik likelihood. Cobaya automatically adds Plik nuisance parameters and its component file includes the SZ prior. Keep the two low-ℓ likelihoods as the exact Commander/SimAll components above. This three-way combination does not name or require the separate Planck lensing likelihood. The code that evaluates the `.clik` trees should be `clipy` (the Python implementation), not the unrelated PyPI project named `clipy` and not an install of the original compiled `clik` code. Cobaya 3.6.2's Planck wrapper imports `clipy`, validates its `clik` interface, disables JAX, and calls the same likelihood files; JAX is unnecessary. Numerical equivalence to the original compiled Planck clik implementation was not measured here, so preserve the implementation choice in the provenance record and compare reference log-likelihood values if exact numerical fidelity is required.

The Planck archive and the likelihood implementation are distinct artifacts: stage/extract the existing Planck data archive, then let the pinned Cobaya installer install its `clipy_0.15` code into the local Cobaya packages path. Do not let `cobaya-install` silently select an arbitrary newer CAMB source build; install the CAMB wheel into the project `.venv`, set `theory.camb.path: global`, and use `--skip-global` for the Cobaya package installation.

## Proposed local commands

These are future commands only; none was run in this task. Work from `cosmology_autoresearch/`. The package and data targets stay within the project. First pin/install the core wheel stack into the existing `.venv`, then create the matching hash lock from the resolved environment:

```bash
.venv/bin/python -m pip install \
  'numpy==2.5.3' 'scipy==1.18.1' \
  'cobaya==3.6.2' 'camb==1.6.6' 'astropy==8.0.1' 'py-bobyqa==1.5.0'
.venv/bin/python -m pip check
.venv/bin/python -m pip freeze --all > work/planck_software_feasibility/resolved-requirements.txt
```

Create `work/planck_software_feasibility/planck_lcdm.yaml` with the three likelihood components above, `theory.camb.path: global`, project-local `packages_path`, and an explicit flat-ΛCDM cosmological parameter/prior block and sampler configuration. Use the six sampled base parameters and explicitly record the fixed neutrino, radiation, BBN helium, recombination, and nuisance-prior choices; do not rely on an undocumented default or import priors from a different paper. The three Cobaya Planck components automatically supply their own nuisance parameters.

Once the configuration is reviewed, extract the already-staged archive only into its package data target, install the Cobaya Planck code/data requirements, and run a short loader/known-point check before any sampling:

```bash
mkdir -p work/planck_software_feasibility/packages/data/planck_2018
tar -xzf work/planck_stage/COM_Likelihood_Data-baseline_R3.00.tar.gz \
  -C work/planck_software_feasibility/packages/data/planck_2018
.venv/bin/cobaya-install work/planck_software_feasibility/planck_lcdm.yaml \
  --packages-path work/planck_software_feasibility/packages --skip-global
CLIPY_NOJAX=1 .venv/bin/clipy_print \
  work/planck_software_feasibility/packages/data/planck_2018/baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik
.venv/bin/cobaya-run work/planck_software_feasibility/planck_lcdm.yaml
```

Before MCMC, verify imports/versions and all three component initializations, evaluate one documented reference parameter point, compare the aggregate/individual likelihood values against an upstream Cobaya or Planck reference, and capture the run YAML, `pip freeze`, data hash, code tag/commit, command, log-likelihood values, output hash and runtime. If any exact version cannot resolve, if a likelihood fails to initialize, or if the reference values disagree beyond a stated tolerance, stop and record the blocker; do not silently substitute Plik-lite, distance priors/chains, original `clik`, or CAMB 2 with the older Cobaya wrapper.

## Caveats for scientific fidelity

- A compatible software stack does not make a paper-specific reproduction. The earlier gate found that the paper shorthand is identifiable but the authors' full nuisance/configuration contract is not provided. Explicitly distinguish a standard PR3 ΛCDM baseline from exact reproduction of another analysis.
- CAMB version, accuracy settings and physical assumptions affect spectra. The chosen 1.6.6 pin is a compatible released CAMB v1 implementation, not an assertion that it is the exact version used by the Planck collaboration or a target paper. Retain explicit CAMB accuracy/recombination settings and report a sensitivity comparison if matching published numbers.
- Cobaya 3.6.2's stock Plik YAML adds automatic nuisance parameters and the SZ prior. Record these defaults, and verify the chosen low-ℓ nuisance calibration setup and priors against the intended baseline.
- `clipy` is a pure-Python reimplementation of `clik` using the official Planck `.clik` data format. It supports Plik, Commander and SimAll. The package README notes that its optional JAX path can introduce up to `1e-4` changes in Plik log-likelihoods on some devices when 64-bit computation is unavailable; Cobaya sets `CLIPY_NOJAX=1` for its wrapper. Preserve that CPU path for the reference run.
- ESA lists the staged data under CC BY-NC 3.0 IGO. Keep the attribution and non-commercial restriction in any downstream distribution/use.

## Primary sources

- [Cobaya 3.6.2 Planck likelihood documentation](https://github.com/CobayaSampler/cobaya/blob/v3.6.2/docs/likelihood_planck.rst): component families, `clipy` versus clik, nuisance parameters and data setup.
- Pinned Cobaya 3.6.2 component definitions: [Plik TTTEEE](https://github.com/CobayaSampler/cobaya/blob/v3.6.2/cobaya/likelihoods/planck_2018_highl_plik/TTTEEE.yaml), [Commander TT](https://github.com/CobayaSampler/cobaya/blob/v3.6.2/cobaya/likelihoods/planck_2018_lowl/TT_clik.yaml), [SimAll EE](https://github.com/CobayaSampler/cobaya/blob/v3.6.2/cobaya/likelihoods/planck_2018_lowl/EE_clik.yaml), [CAMB wrapper](https://github.com/CobayaSampler/cobaya/blob/v3.6.2/cobaya/theories/camb/camb.py), and [Planck `clipy` wrapper](https://github.com/CobayaSampler/cobaya/blob/v3.6.2/cobaya/likelihoods/base_classes/planck_clik.py).
- [Cobaya 3.6.2 installation of cosmology data](https://github.com/CobayaSampler/cobaya/blob/v3.6.2/docs/installation_cosmo.rst).
- [Cobaya v3.6.2 release and pinned commit](https://github.com/CobayaSampler/cobaya/releases/tag/v3.6.2), [CAMB 1.6.6 source metadata](https://github.com/cmbant/CAMB/blob/1.6.6/pyproject.toml), [CAMB 1.6.6 PyPI artifact list](https://pypi.org/project/camb/1.6.6/), [clipy-like 0.15 PyPI metadata/artifact](https://pypi.org/project/clipy-like/0.15/), and [Astropy 8.0.1 metadata](https://pypi.org/project/astropy/8.0.1/).
- [Py-BOBYQA 1.5.0 release/requirements](https://pypi.org/project/Py-BOBYQA/1.5.0/) and [upstream `pyproject.toml`](https://github.com/numericalalgorithmsgroup/pybobyqa/blob/v1.5.0/pyproject.toml).
- [clipy v0.15 README/history](https://github.com/benabed/clipy/tree/clipy_0.15) and [Cobaya's CAMB-v2 support commit, 2026-07-24](https://github.com/CobayaSampler/cobaya/commit/497fa9aebc3dbce2174416645300d2ee7e07f6d0); [CAMB 2.0 release compatibility notes](https://github.com/cmbant/CAMB/releases/tag/v2.0.0).
- [SciPy 1.18.1 source metadata](https://github.com/scipy/scipy/blob/v1.18.1/pyproject.toml) and [release notes](https://github.com/scipy/scipy/releases/tag/v1.18.1), [NumPy release notes](https://numpy.org/doc/2.5/release/2.5.3-notes.html) (2.5.3).
- [Official ESA Planck PR3 data/code and license page](https://esdcdoi.esac.esa.int/doi/html/data/astronomy/planck/Cosmology.html); [Planck likelihood code overview](https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/CMB_spectrum_%26_Likelihood_Code).

## Execution record

Actual work: source/document inspection and local metadata reading only. No package manager command, extraction, import/runtime probe, likelihood evaluation, fitting, or GPU activity. The project brief configures this worker role as GPT-6 Luna Max; the delegated runtime did not expose an independently verifiable model/effort field, so the exact reasoning-effort setting is recorded as unavailable rather than inferred.

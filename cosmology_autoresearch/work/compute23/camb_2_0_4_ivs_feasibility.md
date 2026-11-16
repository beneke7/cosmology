# CAMB 2.0.4 IVS source-feasibility audit

Checked 2026-09-25 01:49 UTC. This is a static provenance/source/docs audit only: no install, compile, CAMB run, fit, or score; the separate IDECAMB implementation was not executed. The official CAMB archive was inspected in-place via read-only tar listing/streaming, not extracted.

## Verdict

The pinned [CAMB 2.0.4 source snapshot](https://github.com/cmbant/CAMB/tree/a6de8cc59124c8bbe3924f1c546f96cef73dcabb) is present and a stock ΛCDM/Γ=0 baseline is supported by CAMB without source changes. It was **not runnable in this workspace at audit time**: `.venv` has no `camb` package, no `gfortran`, `ifort`, `flang`, or `flang-new` is on `PATH`, and the source archive lacks the populated `forutils` submodule. The official [CAMB installation documentation](https://camb.readthedocs.io/en/latest/) says the standard PyPI install includes precompiled code (no Fortran compiler required), whereas building edited source requires a Fortran compiler and recommends a recursive source checkout. The official [PyPI 2.0.4 package attestation](https://pypi.org/project/camb/) ties its artifacts to the same commit. Thus an official 2.0.4 wheel is a straightforward future baseline route, but is not suitable for implementing the interaction; the source build needs its compiler and submodule prerequisites restored.

The IVS is not a parameter-only `DarkEnergyFluid(w=-1)` case. Implementing the paper's background and perturbation law requires coordinated non-Λ vacuum and CDM changes. A custom vacuum background hook can preserve the paper's δρ_v=0 perturbation closure, but changing only that hook while keeping CAMB's standard ρ_c∝a⁻³ evolution would violate the coupled continuity equations.

## Exact source and local artifact

The manifest entry is `camb_source_2_0_4`, tagged `2.0.4`, pinned to upstream commit `a6de8cc59124c8bbe3924f1c546f96cef73dcabb`. The local file `context/code/CAMB_2.0.4_a6de8cc59124c8bbe3924f1c546f96cef73dcabb.tar.gz` is 15,959,613 bytes, SHA-256 `66808001e908982857c312cc84671c66a51312251a7b15e07507d66b4788c0c7`; its `fortran/config.f90` identifies version `2.0.4`. The GitHub-linked PyPI 2.0.4 provenance also identifies the same upstream commit. The manifest itself is SHA-256 `7dfe09299bb4d87b835cc5a2ca997d377b5fe027cf2eb7263d7326c50f16b4d4`.

The archive includes the core Fortran files and `LICENCE.txt`, but its `.gitmodules` declares `forutils` at `https://github.com/cmbant/forutils` and only an empty `forutils/` directory is present. `fortran/Makefile_main` requires the `libforutils` submake, so this tarball cannot build until the submodule contents are restored. `setup.py` can remove the empty directory and try to clone an unpinned `master` automatically, which is not a reproducible source-build route. The local artifact/manifest does not record the submodule gitlink SHA, so that exact dependency revision remains unresolved until recovered from a recursive checkout of the pinned CAMB commit.

## Model-to-CAMB mapping

The checked primary paper defines (Q=\Gamma\rho_v), with (w_v=-1), and its conservation convention implies

\[
\dot\rho_v=\Gamma\rho_v,\qquad
\dot\rho_c+3H\rho_c=-\Gamma\rho_v.
\]

It states zero momentum transfer in the CDM rest frame and says the vacuum has no perturbative contribution in IVS. Its synchronous-gauge equations (19)–(20), with δρ_v=0, give

\[
\delta'_c=-\left(\theta_c+\frac{h'}2\right)
 + a\Gamma\frac{\rho_v}{\rho_c}\delta_c,
\qquad \theta'_c=-\mathcal H\theta_c.
\]

The paper uses Γ>0 for forward-time CDM-to-vacuum transfer. This sign is consistent across its background and perturbation equations. These are paper specifications, not a claim that the authors' unpublished modified solver/configuration is recovered. The exact local paper PDF is `context/papers/interacting_de_desi_dr2_2026.pdf`, SHA-256 `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`.

CAMB 2.0.4 has a useful but incomplete extension point:

- `TDarkEnergyModel` defaults to a cosmological-constant shortcut with zero DE perturbation equations. `BackgroundDensityAndPressure` returns constant `grhov*a*a` under that flag; when non-Λ it calls the model's density and (w(a)) functions. The base stress-energy callback returns zero DE density perturbation and heat flux, and the base initial perturbation callback initializes its vector to zero. Thus the requested smooth-vacuum closure is representable, provided the dynamic vacuum is **not** sent through the cosmological-constant shortcut and the new class explicitly retains zero DE perturbative stress-energy/equations.
- `DarkEnergyFluid` sets (w=-1, w_a=0) to `is_cosmological_constant` and allocates zero perturbation equations; it therefore routes to constant Λ background evolution, not ρ_v evolving under (\dot\rho_v=\Gamma\rho_v). It is not sufficient by setting parameters alone.
- CAMB assumes uncoupled dust in several places. In `results.f90`, `grho_no_de` includes ((grhoc+grhob)a), and background species output uses `grhoc*a`; in `equations.f90`, the perturbation/source background uses `grhoc_t=State%grhoc/a`, while the native synchronous CDM continuity line is `clxcdot=-k*z`. The latter lacks the IVS interaction source.

Accordingly, the minimum source implementation is a **joint, H-consistent background history** for ρ_c(a) and ρ_v(a), solving the two continuity equations with Friedmann closure, and wiring those histories consistently into the expansion, perturbation-source density weighting, background-density outputs, and any remaining routines that assume ρ_c∝a⁻³. The native synchronous CDM density equation then needs the additional conformal-time term (+a\Gamma(\rho_v/\rho_c)\delta_c), using CAMB internal units; its CDM-comoving/geodesic gauge choice supplies the no-momentum-transfer condition (with the corresponding velocity initialization checked). A Γ=0 regression must recover stock 2.0.4 ΛCDM. For a CMB implementation, early-time regular-mode/initial-condition and numerical-stability tests are additional requirements; a background-only profile is not a substitute.

The core inspected paths are the [pinned `DarkEnergyInterface.f90`](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/DarkEnergyInterface.f90), [`DarkEnergyFluid.f90`](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/DarkEnergyFluid.f90), [`equations.f90`](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/equations.f90), and [`results.f90`](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/results.f90). CAMB's [dark-energy model documentation](https://camb.readthedocs.io/en/latest/dark_energy.html) and [source-modification guide](https://camb.readthedocs.io/en/latest/modifying_code.html) support using a custom model/source extension rather than treating a background-only parameter change as sufficient.

## Next-command recommendation and limitations

For a stock baseline only, prepare an isolated environment with the official `camb==2.0.4` precompiled wheel, record the package hash/environment, and run a declared stock ΛCDM case. For the interacting source path, first create a recursive checkout at the immutable SHA above and record `git rev-parse HEAD` plus `git submodule status --recursive`; resolve the exact `forutils` gitlink, provision a supported Fortran compiler, then implement the coupled history and perturbation source. No command was run to do any of these actions.

The CAMB path is now clear but source execution is blocked on local compiler/submodule prerequisites. This audit does not run an LCDM reference, verify a modified build, establish perturbation initial conditions, reproduce the paper's private setup, or make any likelihood/inference claim.

## Primary links checked

- [CAMB official 2.0.4 source commit](https://github.com/cmbant/CAMB/tree/a6de8cc59124c8bbe3924f1c546f96cef73dcabb); [version field in pinned `config.f90`](https://github.com/cmbant/CAMB/blob/a6de8cc59124c8bbe3924f1c546f96cef73dcabb/fortran/config.f90).
- [CAMB current installation/docs](https://camb.readthedocs.io/en/latest/) and [official PyPI package/provenance](https://pypi.org/project/camb/) (2.0.4 source/wheels, provenance names the same CAMB commit).
- [Primary paper, arXiv v1 record](https://arxiv.org/abs/2609.05660v1), local PDF hash recorded above.

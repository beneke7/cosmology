# IDECAMB compatibility card

**Read depth:** exact arXiv v2 PDF/text; focused reading of §§3.2–3.3, 4.1, 5.1–5.2 (printed pp. 5–7, 10), plus a read-only inspection of the pinned code archive's README, `test_ide.ini`, `modification_manual.f90`, and `camb/equations_ppfi.f90`. This is not a full audit of the paper's likelihood or all code paths, and the downloaded code was not executed.

## Provenance

- Yun-He Li & Xin Zhang, *IDECAMB: an implementation of interacting dark energy cosmology in CAMB*, arXiv:2306.01593v2 (2023-09-13), JCAP 09 (2023) 046. Local PDF: `context/papers/idecamb_implementation_2023_v2.pdf`, SHA-256 `346ad0b9aaa68bde2fff1a316fd980b861d9d3d6be6f8645ed5f1d53c04a4b15`; searchable text is adjacent.
- Read-only code snapshot: `context/code/idecamb_4f1093d9efe46f28cf7e2acb4d07ae116ad5e075.tar.gz`, SHA-256 `c167aa0e0eae369a5b94816343ac6bdc8e3bbee32c5cf66fd116ad71450a6e44`. Git commit `4f1093d9efe46f28cf7e2acb4d07ae116ad5e075`, dated 2023-06-20. Repository README says it is a CAMB/CosmoMC patch to be copied over the Planck-2018 package. GitHub's repository API returned no license metadata; do not redistribute or execute it without resolving licensing and compatibility.
- The official arXiv record lists v2 and JCAP 09 (2023) 046: <https://arxiv.org/abs/2306.01593>. The code repository is <https://github.com/liaocrane/IDECAMB>.

## Relevance to the target IVS law

The target paper uses `Q = Gamma*rho_x`, `w_x=-1`, and chooses no momentum transfer in the CDM rest frame. Its synchronous/CDM-comoving density equation is

```text
delta_c' = -theta_c - h'/2 + a*Gamma*(rho_x/rho_c)*delta_c
theta_c' = -Hcal*theta_c
```

IDECAMB's coupled-fluid parameterization explicitly supports `Q = beta*H0*rho_de` and `Q_mu = Q*u_mu,c` (Sec. 3.2/Table 1; Sec. 3.3). With `beta = Gamma/H0`, these match the target background rate and transfer-vector direction. Its general synchronous CDM equation (Eq. 5.8) reduces algebraically to the target equation if the DE density perturbation is set to zero: its mapped coefficients are `C1=Q`, `C2=C3=0`, and `D1=0`, `D2=Q`; set synchronous lapse/shift to zero, CDM velocity to zero, and `delta_de=0`.

The code confirms a branch `H0rde` implementing the H0-scaled DE-density interaction and a `uc` momentum-transfer branch. This is a valuable implementation reference, not proof that the whole solver matches the target. IDECAMB's paper adopts PPF evolution for coupled-fluid perturbations to address possible large-scale instabilities. Its generic PPF DE contribution has not been shown equivalent to the target paper's stated `delta_de=0` closure. The code's sample `test_ide.ini` instead selects `Q=beta*H0*rho_c`, `w0=-1.05`, and PPF; it is not the target setup. The 2023 patch also targets an older CosmoMC/CAMB tree, not current CAMB 2.0.4.

## Assessment and next use

This is a meaningful shortcut for background-rate and transfer-direction implementation, and a source for differential checks. It does **not** recover the target authors' private/unspecified modified-CAMB version, initial conditions, likelihood YAML, or results; it must not be called a reproduction. A faithful independent path should implement the target's explicit vacuum closure in a current, pinned solver and compare it against IDECAMB only where the equations truly coincide. First verify the radiation-era regular mode, then recover standard CAMB spectra at `Gamma=0`, and test the interacting CMB source spectra under tolerance/initial-time refinement. Keep the two closures as separate implementations until a mode-by-mode comparison establishes equivalence.

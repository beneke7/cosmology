# Independent principal-mode / Kcrit check

## Result

On the requested `dr=.25`, `dz=.125 kpc`, `Rmax=80`, `Zmax=60 kpc` grid, the
independent first generalized eigenvalue is
`Kcrit=5.3713874645e-5 (km/s)^-2`. This differs from the reported production
`5.3714678825e-5` by `-8.04e-10`, or `-1.50e-5` fractionally. The generalized
eigen residual is `4.53e-12`; the sign-chosen principal eigenvector is positive
throughout the full grid, with normalized minimum `0.03092` and maximum `1`.

The principal-mode limiting curve has speeds from `24.78` to `141.32 km/s`
over the 43 SPARC radii and a descriptive diagonal score `chi2=3586.15`. This
curve is calculated from the eigenmode limit
`v²=-R/Kcrit * d(log h_mid)/dR`. It describes the branch-limit shape as
`K -> Kcrit` from below; it is not a finite-`u` solution and is not asserted to
be a physical solution at the singular point.

## Method and reproducibility

[principal_mode_audit.py](principal_mode_audit.py) uses the task-owned
independent conductance-graph finite-volume assembly from
[audit_axisymmetric.py](audit_axisymmetric.py), not the production solver or
branch-edge analysis. It forms the smallest positive generalized eigenpair
`B h = Kcrit diag(D) h`, where `B` is the integrated cylindrical `-Laplacian`
with axis regularity, midplane Neumann symmetry, and the monopole Robin outer
boundary; `D=4 pi G rho dV`. Shift-invert settings were
`eigsh(B,k=1,M=diag(D),sigma=0,which='LM',tol=2e-8,maxiter=4000)`. The matrix
is exactly symmetric in the assembled representation. The source CSV hash
matches the frozen snapshot:
`2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5`.

The prior task-owned coarse-grid eigenvalue at the same source normalization
was `5.4704685295e-5`, 1.84% above this fine-grid result, showing the expected
grid dependence. The fine-grid value agrees with the production value to
`0.0015%`.

Full eigenpair, normalized midplane mode, 43-point limiting curve, diagonal
score, hashes, and runtime are in [principal_mode_result.json](principal_mode_result.json).
Code SHA-256 is
`18eb3766ce69b2b32b3bd6cfc4de9f85d0510544a433982910a7bda19a126371`;
result SHA-256 is
`b013f8887c7ca12bb440f40ee6cb9e527014d103978ce417047a79cfc337b158`.

Recorded usage was `0.967 CPU seconds` and `0.968 wall seconds`, single-thread
BLAS, no GPU. The prior coarse-grid eigenvalue is the requested resolution
variant; no additional solve was needed.

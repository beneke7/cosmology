# Independent NGC 3198 axisymmetric BVP audit

## Result

The requested independent audit of the `K=3.9e-5 (km/s)^-2`,
`Upsilon*=0.762` production point is complete. The production profile hash
matches the frozen input. An independently assembled finite-volume solve
reproduces the production fine-grid speeds closely, keeps `u` positive over
the full domain, and has a small discrete residual. A coarser grid shows a
material speed error at small radii, which is quantified below. The K-to-zero
Newtonian check is consistent with the hashed Hankel reference within the
independent grid/source quadrature error.

The requested near-upper-K follow-up is still pending: the production file
available for this audit ends at `3.9e-5`. The parent agent is extending the K
list; no additional K was inferred or substituted while those output records
were unavailable.

## Equation, units, and independent discretization

For `u=exp(-K phi)`, the stationary axisymmetric field equation is

\[
  \frac1R\partial_R(R\partial_Ru)+\partial_z^2u
  +4\pi G K\rho u=0.
\]

I used `w=(u-1)/K`, so

\[
  (-\nabla^2-4\pi G K\rho)w=4\pi G\rho,
  \qquad u=1+Kw.
\]

This gives a direct K=0 limit without subtracting nearly equal values of `u`.
Here `G=4.30091e-6 kpc (km/s)^2/Msun`, `rho` is in `Msun/kpc^3`, and K is
in `(km/s)^-2`. Thus `4 pi G K rho` is in `kpc^-2`. With cylindrical volume
weights omitting common `2 pi`, the integrated matrix `B`, source vector
`D=4 pi G rho dV`, and `K D` all have units `(km/s)^2 kpc`; `w` is in
`(km/s)^2`.

The script [audit_axisymmetric.py](audit_axisymmetric.py) independently
constructs the cylindrical conductance graph and source from the frozen
profile CSV. It does not import or reuse the production `build_operator` or
`solve_point`. The tested cell-centered domain is `0<=R<=80 kpc`,
`0<=z<=60 kpc`; the mandatory coarse grid is `dr=.5`, `dz=.25 kpc` (38,400
cells). A second resolution run uses `.25`, `.125 kpc` (153,600 cells).
The radial source interpolation is PCHIP; the stellar tail is independently
fit to the specified measured rings and the published source model is
reconstructed: `Upsilon*=.762`, stellar `h=.3 kpc`, HALOGAS HI+He with 85%
`h=.2 kpc` plus 15% `h=3 kpc`, and a 5-kpc gas taper. The resulting stellar
tail scale is `4.05533 kpc`.

Even reflection at the midplane imposes zero vertical flux at `z=0`; no lower
face link is assembled. No inner radial link is assembled at `R=0`, which is
regular zero radial flux. At the outer faces I impose the isolated monopole
condition on `w`,
`dn w=-(n.r)/r^2 w`. The boundary face value includes the half-cell Robin
resistance, `w_b=w_c/(1+kappa*Delta_n/2)`, and the corresponding flux is
integrated over that face. This is an explicit monopole exterior
approximation on the finite cylindrical box.

The midplane speed is computed as
`v^2=-R (dw/dR)/(1+K w)`, using a quadratic even-reflection estimate from the
first two z-centers and a PCHIP radial derivative. The K=0 reference uses
`v^2=-R dw/dR`. The linear system is solved with Jacobi-preconditioned CG;
the measured residual is `||A w-D||_2/||D||_2`.

## Hashes and provenance

- Frozen source CSV SHA-256:
  `2dffddbcd53854b2f8baf810c7500a753f9a2939f23d6b5deb11fdf9683f10b5`
  (matches both the requested snapshot and production summary).
- Production candidate JSONL SHA-256:
  `8d767c07637ecef505b47c87c289cbc40d9f9688894c08c5242a6e5aa922f364`.
- Production summary SHA-256:
  `e8bd891ab4409ed3c0be65df019c6f9a98d2399c8e9caa5e499a03cecf91d028`.
- Independent code SHA-256:
  `1e7b2e3460394ec182f365027c2e85a95ca90d7052b2e8c641e54b0abe3c28e8`.
- Independent coarse result SHA-256:
  `166392208c58a34dd16db17b41b610315d8b6d33d7fcb060f4ffd4fcc0fa45b7`.
- Independent fine result SHA-256:
  `7587009b63f47e525d5f11d4ace6146bd391380a8934b7525941dfe56fdd99e7`.

The compact numerical outputs are [coarse.json](coarse.json) and
[fine.json](fine.json). They include the 43-point production and independent
speed vectors, residuals, positivity extrema, input hashes, and runtime.

## Mandatory K=3.9e-5 point

| Check | Independent result |
|---|---:|
| Fine grid (`dr=.25`, `dz=.125 kpc`) full-grid `u` range | `[1.1804200, 5.8482205]` |
| Fine-grid integrated PDE residual | `3.88e-10` |
| Fine-grid Robin flux condition relative mismatch | `4.29e-13` |
| Fine-grid vs production speed max absolute difference | `0.00206 km/s` |
| Fine-grid vs production speed RMS difference | `0.000666 km/s` |
| Fine-grid max relative speed difference | `2.30e-5` |

The production solve reports its positive-u check only at interpolated
midplane points. The independent extrema above include every cell in the
upper-half `(R,z)` domain. The Robin mismatch is computed from the reconstructed
face value and normal derivative on both outer faces.

On the required coarser grid, `u` ranges from `1.16815` to `5.53877`; the
integrated linear residual is `1.19e-10`. Speeds differ from production by
`7.03 km/s` maximum and `1.90 km/s` RMS over all 43 radii. The largest
difference is at the innermost datum (`R=.32 kpc`); even excluding that datum,
the coarse-to-fine maximum is `2.33 km/s`. This is a real radial/vertical
resolution effect and means the coarser grid is useful as an independent
operator check, but not as a high-precision curve prediction.

The transformed operator has the expected sign: positive baryon density
produces positive `w`, and increasing K reduces the matrix `A=B-KD`. For the
independent coarse operator I also solved the generalized eigenproblem

\[
  B h=K_{\rm crit}\,\mathrm{diag}(D)h
\]

using `scipy.sparse.linalg.eigsh(B, k=1, M=diags(D), sigma=0,
which="LM", tol=2e-8, maxiter=4000)`. This shift-invert calculation targets
the smallest positive generalized eigenvalue directly, rather than inferring
an edge from CG failure. It gives `Kcrit=5.47046853e-5 (km/s)^-2`,
`K/Kcrit=0.71292`, with generalized eigen residual
`||Bh-Kcrit Dh||/||Bh||=1.46e-12`. This is the coarse-grid, finite-box
critical value for the declared operator, not a continuum edge estimate.

## K-to-zero comparison with the Hankel baseline

The frozen Hankel table is normalized to stellar `Upsilon=0.6`; I compared at
the requested `0.762` by rescaling its stellar `V^2` by `.762/.6` while
leaving the gas `V^2` unchanged, then adding components in quadrature. This
avoids attributing a normalization change to the PDE discretization.

At fine resolution the maximum and RMS differences from this adjusted Hankel
curve are `1.365` and `0.490 km/s` over the 43 SPARC radii. For `R>=1 kpc`
they are `0.809` and `0.454 km/s`; the full-range maximum is driven by the
innermost point, where radial differentiation/extrapolation is most sensitive.
At coarse resolution these become `5.473` and `1.492 km/s`, respectively.
The improvement with resolution explicitly demonstrates that the K=0
discrepancy contains finite-volume source quadrature and interpolation error.
The remaining fine-grid offset also includes finite box and vertical-source
quadrature differences from the Hankel integration, so it is not interpreted
as a failure of the K->0 equation limit.

## CPU budget and remaining check

The independent numerical solves and coarse generalized eigenvalue used
approximately `17.7 CPU seconds` total (about `18 seconds` wall time), with
single-thread BLAS and no GPU. This is below the 300 CPU-second and 1,800
wall-second caps. The near-upper-K production point remains pending until the
extended production candidate output is available; no extra solve was started
without its selected K and root record.

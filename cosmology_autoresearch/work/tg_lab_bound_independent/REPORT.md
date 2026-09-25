# Independent static higher-gradient / torsion-balance audit

## Result

**Independently checked:** for the paper's stationary linear sector, the operator and real-root Yukawa algebra in Eqs. (11)–(13) are consistent, including the signs, dimensions, sum rule, one-length limit, and repeated-root limit. The quoted Lee et al. 38.6 μm crossing is present in the primary single-Yukawa |α| curve; the independent HUST primary curve crosses |α|=1 at 48 μm, so it does not strengthen that particular range threshold. The latest overview located (arXiv v2, 8 Sep 2026) still points to the UW and HUST studies as the leading tabletop torsion-balance tests in this band.

As a primary-source status check dated 25 Sep 2026, the Eöt-Wash group's own “Current Published Results” page points to the Lee et al. 2020 curve; the latest primary torsion-balance curves located for the relevant band are therefore Lee (38.6 μm) and Tan et al. (48 μm). The 2024 lattice-atom-interferometer paper measures attraction from a miniature mass and rules out screened-fifth-force models over their natural parameter space, but it does not report a replacement gravitational-strength Yukawa |α|=1 exclusion curve at these ranges. This is a scoped literature check, not a proof that no unpublished or differently parameterized result exists. Source links are in [SOURCES.md](SOURCES.md).

Two qualifications change the simple Eq. (14) reading:

1. The single-length branch has α = −1, so its relevant Lee curve is the negative-sign column, not an unsigned generic fifth-force amplitude. Table 11.4 gives −α₉₅ = −1.27 at λ = 25 μm and −0.886 at 28 μm. Thus the published sign-resolved profile contour places the α = −1 crossing between 25 and 28 μm; straight interpolation of the rounded points gives about 27.1 μm. This is a read-off from the paper's two-parameter contour, not a newly calibrated one-parameter likelihood. It shows the cited 38.6 μm unsigned bound is conservative for the one-length sign, but the exact one-length 95% limit should not be quoted as 27.1 μm without re-profiling the original likelihood.
2. For two real roots, ℓ1² = λₐ² + λᵦ². A limit on the longest range λₐ alone therefore does not imply ℓ1 < 38.6 μm: in the degenerate-root limit λₐ = λᵦ = ℓ1/√2. If one provisionally transferred λₐ < 38.6 μm to the two-root case, the algebra would give only ℓ1 < √2·38.6 = 54.6 μm (real-root sector). That transfer itself still needs an experiment-specific response calculation.

## Independent derivation

With ∆ → −k², the stationary equation

`Δ(1 − ℓ₁² Δ + ℓ₂⁴ Δ²) φ = 4πGρ`

becomes

`φ̃(k) = −4πG ρ̃(k) / [k²(1 + ℓ₁²k² + ℓ₂⁴k⁴)]`.

The denominator is positive at real k for ℓ1² > 0 and ℓ2⁴ ≥ 0. Writing it as `(1 + a k²)(1 + b k²)` requires `a+b=ℓ₁²`, `ab=ℓ₂⁴`; `a,b` have units length². For `a>b>0`, inverse Fourier transformation gives

`φ(r) = −GM/r [1 + αₐ exp(−r/√a) + αᵦ exp(−r/√b)]`,

`αₐ = −a/(a−b) = −λₐ²/(λₐ²−λᵦ²) ≤ −1`, `αᵦ = b/(a−b) = −αₐ−1 ≥ 0`, and `αₐ+αᵦ=−1`. Here `λₐ=√a`, `λᵦ=√b` and all α are dimensionless. Since `[φ]=L²T⁻²` and `[G]=L³M⁻¹T⁻²`, every term in the paper's energy density, including `ℓ₁²(Δφ)²/G` and `ℓ₂⁴|∇³φ|²/G`, has energy-density units.

For one length, b = 0 (ℓ2 = 0), giving λₐ = ℓ1 and α = −1 exactly. At a = b = λ² (ℓ1 = √2λ, ℓ2 = λ), the separate α coefficients diverge but their combined limit is finite:

`αₐe^(−r/λₐ)+αᵦe^(−r/λᵦ) → −(1 + r/(2λ))e^(−r/λ)`.

The regularized point potential tends to −GM/(λₐ+λᵦ), or −GM/(2λ) in the degenerate case and −GM/ℓ1 in the one-length case. Real roots require ℓ1⁴ ≥ 4ℓ2⁴; below this discriminant the roots are complex and the correction is damped/oscillatory, so real-α Yukawa curves do not directly test that branch.

Opposite signs in Eq. (13) do not create a zero in the static Green function. The exact Fourier transfer is

`φ̃/φ̃_N = 1/[(1 + a k²)(1 + b k²)]`,

so the full response is positive and the deviation from Newton is negative at every k>0. In fact the suppression relative to Newton is greater than for the long-pole one-length reference `α=−1, λ=√a` at every mode. In point space, the Yukawa sum is also strictly negative for r≥0; the two terms partly offset, but do not flip sign or cancel. The total point-mass force remains attractive for r>0 while being reduced relative to Newton.

This mode-by-mode statement is not a full torsion-balance fit. Torsion experiments constrain geometry-weighted, differentiated Yukawa kernels. The single-range curve cannot by itself provide a statistically calibrated contour in `(a,b)`; that needs both range templates evaluated on the experimental geometry and a nuisance re-profile against the measured torques. The root's stored thesis includes analyzed torque rows, but this task did not reimplement that full apparatus model. No two-range experimental exclusion is claimed.

## Reproducible check and next test

Run from `cosmology_autoresearch/`:

```sh
python3 work/tg_lab_bound_independent/check_factorization.py
```

It checks the factorization identities, amplitude sum, Fourier transfer, real-space sign, force sign, one-length and repeated-root limits, and the displayed Table 11.4 crossing brackets. It passed with Python 3 (standard library only). The task was CPU-only; no GPU or external software was used. Source versions, URLs, hashes and access notes are in [SOURCES.md](SOURCES.md).

The sector audited here is the stationary higher-gradient model in arXiv:2609.00317v1; it is distinct from the 2024 nonlinear galaxy K model.

**Next discriminating check:** evaluate the exact two-range torque templates for the UW apparatus and jointly profile `(ℓ₁,ℓ₂)` plus the paper's nuisance parameters against the Appendix A.1 torque rows, then calibrate a 95% contour. That is the step needed to turn the static no-cancellation result into a laboratory bound on the full two-length model.

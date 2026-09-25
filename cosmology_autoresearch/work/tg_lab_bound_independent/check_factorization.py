#!/usr/bin/env python3
"""Independent algebra and numerical-limit audit for arXiv:2609.00317v1.

Run from the project root with:
    python work/tg_lab_bound_independent/check_factorization.py

Uses only Python's standard library. Lengths in the numerical checks are SI.
The four Lee et al. table entries at the bottom are transcribed from Table
11.4 of the first author's 2020 UW dissertation (source record in SOURCES.md).
"""

from __future__ import annotations

from math import exp, isclose, sqrt


def real_roots(ell1: float, ell2: float) -> tuple[float, float] | None:
    """Return a >= b >= 0 in (1+a k^2)(1+b k^2), or None if roots complex."""
    disc = ell1**4 - 4.0 * ell2**4
    scale = max(ell1**4, 4.0 * ell2**4)
    if disc < -1e-14 * scale:
        return None
    disc = max(disc, 0.0)
    return (ell1**2 + sqrt(disc)) / 2.0, (ell1**2 - sqrt(disc)) / 2.0


def amplitudes(a: float, b: float) -> tuple[float, float]:
    if b == 0.0:
        return -1.0, 0.0
    if isclose(a, b, rel_tol=1e-14, abs_tol=0.0):
        raise ValueError("degenerate root: use the repeated-pole limit")
    return -a / (a - b), b / (a - b)


def yukawa_sum(r: float, a: float, b: float) -> float:
    """P(r)=sum alpha_i exp(-r/lambda_i), excluding the Newtonian 1."""
    if b == 0.0:
        return -exp(-r / sqrt(a))
    if isclose(a, b, rel_tol=1e-14, abs_tol=0.0):
        lam = sqrt(a)
        return -(1.0 + r / (2.0 * lam)) * exp(-r / lam)
    aa, ab = amplitudes(a, b)
    return aa * exp(-r / sqrt(a)) + ab * exp(-r / sqrt(b))


def fourier_transfer(k: float, a: float, b: float) -> float:
    """Static potential / Newtonian potential for real roots, including b=0."""
    return 1.0 / ((1.0 + a * k * k) * (1.0 + b * k * k))


def force_multiplier(r: float, a: float, b: float) -> float:
    """Radial force magnitude divided by Newtonian, for real roots."""
    if b == 0.0:
        lam = sqrt(a)
        return 1.0 - (1.0 + r / lam) * exp(-r / lam)
    if isclose(a, b, rel_tol=1e-14, abs_tol=0.0):
        lam = sqrt(a)
        # Force is -d/dr[-GM/r (1 + P)] in the inward radial direction.
        return 1.0 - (1.0 + r / lam + r * r / (2.0 * lam * lam)) * exp(-r / lam)
    aa, ab = amplitudes(a, b)
    return 1.0 + aa * (1.0 + r / sqrt(a)) * exp(-r / sqrt(a)) + ab * (
        1.0 + r / sqrt(b)
    ) * exp(-r / sqrt(b))


def close(x: float, y: float, tol: float = 2e-12) -> None:
    assert abs(x - y) <= tol * max(1.0, abs(x), abs(y)), (x, y)


def main() -> None:
    # Generic real-root example and Eq. (12) coefficient check.
    ell1, ell2 = 50e-6, 20e-6
    roots = real_roots(ell1, ell2)
    assert roots is not None
    a, b = roots
    close(a + b, ell1**2)
    close(a * b, ell2**4)
    la, lb = sqrt(a), sqrt(b)
    aa, ab = amplitudes(a, b)
    assert la >= lb > 0.0 and aa <= -1.0 and ab >= 0.0
    close(aa + ab, -1.0)

    # Eq. (11)/(13): the real-space potential correction has fixed sign.
    for r in (0.0, 1e-7, 1e-5, 1e-3):
        p = yukawa_sum(r, a, b)
        assert p <= 1e-15
        if r > 0.0:
            # More precisely, |P_two-range| >= exp(-r/lambda_a).
            assert -p + 1e-14 >= exp(-r / la)

    # Exact Fourier-space factorization, including relative mode suppression.
    for k in (0.0, 1e2, 1e4, 1e6, 1e8):
        s_partial = 1.0 + aa * a * k * k / (1.0 + a * k * k) + ab * b * k * k / (
            1.0 + b * k * k
        )
        close(s_partial, fourier_transfer(k, a, b))
        assert 0.0 < fourier_transfer(k, a, b) <= 1.0
        if k > 0.0:
            s_one_long_pole = 1.0 / (1.0 + a * k * k)
            assert fourier_transfer(k, a, b) < s_one_long_pole

    # The total point-pair force stays attractive for r>0 and vanishes at r=0.
    close(force_multiplier(0.0, a, b), 0.0)
    assert all(force_multiplier(r, a, b) > 0.0 for r in (1e-8, 1e-6, 1e-4, 1e-2))

    # One-length limit b=0: lambda=ell1, alpha=-1, finite potential at r=0.
    one = real_roots(40e-6, 0.0)
    assert one is not None
    ao, bo = one
    close(ao, (40e-6) ** 2)
    close(bo, 0.0)
    close(amplitudes(ao, bo)[0], -1.0)
    close(-yukawa_sum(0.0, ao, bo), 1.0)

    # Repeated-root limit: a=b=ell1^2/2=ell2^2; divergent separate alphas
    # combine into -(1+r/(2 lambda))*exp(-r/lambda).
    lam = 31e-6
    aeq = lam**2
    ell1_eq = sqrt(2.0) * lam
    ell2_eq = lam
    eqroots = real_roots(ell1_eq, ell2_eq)
    assert eqroots is not None
    close(eqroots[0], aeq)
    close(eqroots[1], aeq)
    r = 17e-6
    target = -(1.0 + r / (2.0 * lam)) * exp(-r / lam)
    eps = 1e-6
    near = yukawa_sum(r, aeq, aeq * (1.0 - eps))
    close(near, target, tol=2e-6)
    close(fourier_transfer(1.2e5, aeq, aeq), 1.0 / (1.0 + aeq * (1.2e5) ** 2) ** 2)
    assert real_roots(40e-6, 0.75 * 40e-6) is None

    # Table 11.4 of Lee (2020 UW thesis): alpha bounds at displayed points.
    # Values shown are (lambda_um, abs-alpha95, negative-alpha95).
    lee = {
        25.0: (6.4, -1.27),
        28.0: (3.73, -0.886),
        35.0: (1.45, -0.360),
        40.0: (0.877, -0.226),
    }
    # Rounded points bracket the published 38.6 um |alpha|=1 crossover.
    assert lee[35.0][0] > 1.0 > lee[40.0][0]
    # The alpha=-1 branch is already outside the tabulated negative-side
    # 95%-profile contour at 28 um, while it remains inside at 25 um.
    assert abs(lee[25.0][1]) > 1.0 > abs(lee[28.0][1])
    # A linear interpolation of the displayed, rounded sign-specific profile
    # points estimates the alpha=-1 crossing; it is not a refitted likelihood.
    neg_cross_um = 25.0 + (1.27 - 1.0) / (1.27 - 0.886) * (28.0 - 25.0)

    print("factorization, dimensions, amplitudes, Fourier transfer, and limits: PASS")
    print(f"real roots for ell1=50 um, ell2=20 um: lambda_a={la*1e6:.6f} um, lambda_b={lb*1e6:.6f} um")
    print(f"amplitudes: alpha_a={aa:.8f}, alpha_b={ab:.8f}, sum={aa+ab:.8f}")
    print(f"degenerate limit: lambda={lam*1e6:.3f} um, ell1={ell1_eq*1e6:.3f} um")
    print("Lee Table 11.4 absolute-|alpha| crossing bracket: 35--40 um; published 38.6 um")
    print(f"Lee Table 11.4 alpha=-1 sign-specific profile bracket: 25--28 um; linear read-off={neg_cross_um:.2f} um (not a re-fit for the one-parameter theory)")
    print("conversion: ell1^2=lambda_a^2+lambda_b^2, so lambda_a<38.6 um alone implies only ell1<sqrt(2)*38.6=54.6 um in the real-root sector")


if __name__ == "__main__":
    main()

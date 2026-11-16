# Independent integrated stationary TG force-identity screen: NGC 3198

Classification: spherical-surrogate necessary-condition screen, not a physical NGC 3198 prediction.

## Derivation

Pszota & Van (2024), Eq. (8), gives `g' = -2g/r - 4 pi G rho - K g^2`. Under the paper's sphericalized source construction, define `b(r)=Upsilon Vdisk^2 + Vgas*abs(Vgas)`, so `4 pi G rho = (r b)'/r^2`. Multiplication by `r^2` and integration from `r_min` gives `r^2 g + r b + K integral[r_min,r] s^2 g(s)^2 ds = C`. With `g=-Vobs^2/r`, the checked identity is:

`r * (b - Vobs^2) + K * integral[r_min,r] Vobs^4 ds = C`.

The constant is free because the first observed point is at 0.32 kpc. No outer observed endpoint is imposed as a boundary.

## Frozen independent calculation

The exact local archive member has 43 rows over 0.32–44.08 kpc. I used the published fixed stellar normalization `Upsilon=0.762 +/- 0.007`, with the quoted uncertainty propagated as a shared nuisance, and the SPARC signed-gas term `Vgas*abs(Vgas)`. The bulge column is zero for all rows.

The cumulative `Vobs^4` integral uses the exact polynomial panel integral for linearly interpolated Vobs; an independent 8-point Gauss-Legendre calculation agrees to maximum relative difference 1.89e-16.

Primary residual is the RMS of `Y+KX-C`, with `Y=r*(b-Vobs^2)`, `X=integral Vobs^4`, `C` free, and K fit as the slope; the reported fractional residual divides by RMS(Y). GLS uses the full covariance from shared observed-velocity errors in both Y and every cumulative X and includes the shared Upsilon error. A seeded 100,000-draw Monte Carlo refits C and K per draw to reflect nonlinear measurement propagation.

| Range | N | OLS K (s²/km²) | C (kpc km²/s²) | RMS residual | RMS residual / RMS(Y) | GLS K | GLS reduced χ² |
|---|---:|---:|---:|---:|---:|---:|---:|
| full_0p32_to_44p08_kpc | 43 | 3.6004424e-05 | -341.77231 | 11625.22 | 0.038507326 | 3.5971677e-05 | 1.8217501 |
| observed_points_r_ge_5_kpc | 28 | 3.6229252e-05 | 2718.0051 | 13880.048 | 0.037107132 | 3.5773668e-05 | 0.93573358 |

Monte Carlo K and residual summaries:

- full_0p32_to_44p08_kpc: K median 3.5971725e-05, mean 3.5977426e-05, SD 5.8677188e-07, 16–84% [3.5397502e-05, 3.6557127e-05], P(K>0)=1; median residual RMS 16652.151.
- observed_points_r_ge_5_kpc: K median 3.6210639e-05, mean 3.6217043e-05, SD 6.8928735e-07, 16–84% [3.5533562e-05, 3.689859e-05], P(K>0)=1; median residual RMS 20022.716.

Published Table 1 comparison: K=3.4e-05 +/- 4e-07 s²/km² and Upsilon=0.762 +/- 0.007. The screen does not refit Upsilon.

## Numerical and physical caveats

- This checks a necessary integrated identity only after replacing the disk geometry by the paper's sphericalized source prescription; it is not a physical galaxy prediction.
- Piecewise-linear interpolation of Vobs is used only to evaluate the cumulative integral; source values are evaluated at the data radii.
- Per-row uncertainties for the gas/disk component curves are unavailable and those curves are held fixed. Distance, inclination, and calibration uncertainties are omitted.
- Observed-velocity errors are independent Gaussian errors in Vobs. The innermost errV exceeds Vobs, so powers Vobs² and Vobs⁴ are nonlinear; Monte Carlo should be favored over the linearized GLS uncertainty for that range.
- The integration constant is free from r_min=0.32 kpc. The outermost observed point is not a boundary condition.

## Provenance and reproducibility

- ZIP SHA-256: `0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588`
- NGC3198_rotmod.dat member SHA-256: `17b774ad11e7dd745a067073b76f1909d4de32fa87616df12b64da3f225cf953`
- Supplied PDF SHA-256: `55cbfbbb0a916ea2ce0b70523334ee301897695860fb4bdb5bfd2278ab47b62e`
- Monte Carlo seed: `20260925`; draws: `100000`.
- Derivation reference: Pszota & Van (2024), Eqs. (8), (47), (53), Table 1; local PDF supplied in the workspace.

## Post-freeze comparison with the root experiment

I opened the root experiment only after freezing this calculation. Its input archive hash, row count, radial range, fixed `Upsilon`, signed-gas convention, exact segment integration, free inner integration constant, and lack of an outer endpoint boundary all agree with this independent setup. The root writes `Y=r*(Vobs^2-b)` and fits `Y=C+KX`; this report uses `Y=r*(b-Vobs^2)` and fits `Y=C-KX`, which is the same identity after reversing Y's sign. The root includes `Upsilon*Vbul^2`; NGC 3198's bulge column is zero in all 43 rows.

| Range | This check: central OLS K | This check: linearized-covariance GLS K | Root Monte Carlo profile K | GLS K difference vs root |
|---|---:|---:|---:|---:|
| Full, 0.32–44.08 kpc | 3.6004424e-05 | 3.5971677e-05 | 3.5953591e-05 | +0.0503% |
| Observed points r ≥ 5 kpc | 3.6229252e-05 | 3.5773668e-05 | 3.5774154e-05 | −0.00136% |

For the full range, my first-order GLS statistic is χ²=74.6918 for 41 degrees of freedom, versus the root's Monte Carlo-covariance profile minimum χ²=74.6494; the difference is 0.0424. For `r ≥ 5 kpc`, the corresponding values are 24.3291 and 24.6249 for 26 degrees of freedom. My 100,000-draw OLS K intervals include each root optimum. The unweighted outer-range slope differs more because it does not weight the velocity uncertainties; the covariance-weighted estimate agrees to 0.0014%.

The root's profile optimum is 3.5953591e-05, compared with the paper value 3.40e-05; its stored profile Δχ² is 54.08. My screen's point estimates are also above the paper value, but its `published_K_difference_sigma` field only divides a point-estimate offset by the paper's quoted K uncertainty and is not a calibrated significance. My single-pass GLS is not an exact reproduction of the root's Monte Carlo covariance profile: it uses analytic first-order propagation, while the root uses 12,000 shared-error draws and a nonlinear bias correction.

My exact polynomial quadrature and independent 8-point Gauss-Legendre integration agree to relative 1.89e-16. The root's second quadrature is trapezoidal integration of Vobs^4; its record reports a best-K shift of +6.103e-9 s²/km² (about +0.017% of K), consistent with a small quadrature sensitivity.

The frozen pre-comparison artifacts are preserved as `independent_result_frozen.json` and `report_frozen.md`. Their hashes are in `frozen_before_comparison.sha256`. The root record SHA-256 is `56f7d4a8dd5a19f13a564cfc2084658b9f1af1b462cca2ca6475a24a2610d968`; root script SHA-256 is `05231c937649c7cc115fb417e53e4b03d69d45457090cf5edecad48ef22daf49`.

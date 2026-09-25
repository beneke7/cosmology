# NGC 3198 Eq. (47) source-contract audit

**Status:** reproduced local source diagnostic; no fit or TG boundary-value solve.

## Inputs

Visually inspected Pszota & Ván, *Field equation of thermodynamic gravity and galactic rotational curves*, Physics of the Dark Universe 46 (2024) 101660, printed pp. 4–5 and 7. User-provided PDF SHA-256: 55cbfbbb0a916ea2ce0b70523334ee301897695860fb4bdb5bfd2278ab47b62e.

Read NGC3198_rotmod.dat from SPARC Rotmod_LTG.zip. Archive SHA-256: 0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588. Member SHA-256: 17b774ad11e7dd745a067073b76f1909d4de32fa87616df12b64da3f225cf953; 43 rows, 2,075 bytes. Local SPARC master paper PDF SHA-256: d089215877213661e40965543ee7e05736619082ad16d95e65ec059029588c63.

## Eq. (47): sign, units, M/L and signed gas

Visual reading gives rho(r) = [Upsilon_1 d(r v_SD^2)/dr + d(r v_HI^2)/dr] / (4 pi G r^2). The derivative sign is plus, as follows from M(<r)=r v_c^2/G and rho=(dM/dr)/(4 pi r^2).

With r in kpc, velocities in km/s and G=4.30091e-6 kpc (km/s)^2/Msun, rho has units Msun/kpc^3. Upsilon multiplies only the stellar-disk derivative, linearly, because Rotmod disk velocities are tabulated for M/L=1; it does not multiply gas or appear squared. Table 1 p.7 reports Upsilon_1=0.762 +/- 0.007 and K=(3.40 +/- 0.04)e-5 s^2/km^2.

There is a source-sign ambiguity: six local Vgas rows (r=1.93–3.54 kpc) are negative. Pszota & Ván p.5 say a negative gas value contributes negatively to total v^2; their Eq.53 preserves this as v_HI*abs(v_HI). Eq.47 visibly prints v_HI^2 with no signed-square definition. I retained that literal form and separately tested the signed convention. This ambiguity does not drive the robust negative-density band, where Vgas>0.

## Interpolation and pseudo-density

The paper only says to use an “interpolated function”; it supplies no interpolator, smoothing or endpoint derivative rule. The attached audit.py interpolates the component velocities with PCHIP, natural cubic spline and Akima, then differentiates r*v^2. A separate finite difference of the mass proxy on 20,001 grid points matches the analytic derivative minimum within 0.5 Msun/kpc^3. Negative density is never clipped.

| Interpolator | Minimum rho (Msun/kpc^3) | Radius (kpc) | Negative intervals (kpc) |
|---|---:|---:|---|
| PCHIP | -156,352 | 15.218 | 14.498–15.997; 43.058–44.080 |
| Natural cubic | -114,510 | 15.400 | 14.507–16.472; 42.912–44.080 |
| Akima | -142,449 | 15.198 | 14.468–15.993; 42.968–44.080 |

A common negative interval is 14.51–15.99 kpc. At the PCHIP minimum, the stellar term contributes -3.35e5 and gas +1.78e5 Msun/kpc^3. So the negative total comes from the spherical inversion of the disk curve, not the small negative-gas rows. The outermost negative interval is less robust to endpoint interpolation. SPARC I (2016), Table 2 note, gives Vdisk for M/L=1; its baryonic-velocity equation uses |Vgas|Vgas and confirms the signed-gas convention. Pszota & Ván §3.2 p.7 themselves say their spherical simplification makes the fitted M/L “not strictly physical.”

## Is the isolated-exterior test data-closed?

No, not as a physical isolated-galaxy prediction. The Rotmod data span 0.32–44.08 kpc (43 rows); the first observed velocity is 24.4 +/- 35.9 km/s and the last is 149.0 +/- 3.0 km/s. At the last radius Vgas=47.84 and Vdisk=61.63 km/s. Force-curve coverage does not specify an inner/outer 3D source or source taper. The local files do not reach r=0 and do not provide a resolved NGC 3198 HI surface-density profile or gas thickness; projected stellar profiles likewise do not uniquely specify 3D density.

A conditional spherical solve is mathematically possible only after declaring a center extrapolation and source taper/edge. For vacuum outside R and gauge u(infinity)=1, u=1+B/r and u'(R)=-(u(R)-1)/R. This must not be replaced by the observed endpoint velocity. The paper's Eqs.41–42 p.4 use observed inner/outer velocities as derivative boundary data, so its fit is not this isolated test.

The Warden's test can therefore stress-test a declared spherical surrogate with the signed gas convention and negative source retained. It cannot yield a unique physical source from these local inputs. A converged failure would limit/falsify that source prescription, not thermodynamic gravity generally. Never clip the negative density.

## Reproduction

From repository root: cosmology_autoresearch/.venv/bin/python cosmology_autoresearch/work/tg_sparc_source_audit/audit.py

The script deterministically regenerates audit.json. Script SHA-256: b8457b74730b1797e71ccf24441d35e04063ca2e77b999fd842662fec12bb1cc.

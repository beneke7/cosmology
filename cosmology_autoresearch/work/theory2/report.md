# Independent radiation-omission check

**Status: independently checked; exploratory fixed-point screen.** This check asks whether the radiation term changes the 13 supplied DESI DR2 BAO distance predictions by an amount comparable to their released full covariance. It is not a posterior analysis or a reproduction of DESI parameter constraints.

## SI derivation

I used (T_{\rm CMB}=2.7255\,\mathrm K), (h_{\rm P}=6.62607015\times10^{-34}\,\mathrm{J\,s}), (k_B=1.380649\times10^{-23}\,\mathrm{J\,K^{-1}}), (c=299792458\,\mathrm{m\,s^{-1}}), and (G=6.67430\times10^{-11}\,\mathrm{m^3\,kg^{-1}\,s^{-2}}). CODATA 2022 lists (c,h_{\rm P},k_B) as exact SI values and (G=6.67430(15)\times10^{-11}) with relative standard uncertainty (2.2\times10^{-5}) (Tables XXXII–XXXIII). Fixsen's literature combination is (2.72548\pm0.00057\,\mathrm K), consistent with the requested rounded input. The IAU exact conventions used for the conversion are (1\,\mathrm{au}=149597870700\,\mathrm m) and (1\,\mathrm{pc}=(648000/\pi)\,\mathrm{au}); thus (1\,\mathrm{Mpc}=3.085677581491367\times10^{22}\,\mathrm m). [CODATA 2022, Tables XXXII–XXXIII](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf), [Fixsen 2009, abstract](https://arxiv.org/abs/0911.1955), [IAU 2012 unit definition](https://iauarchive.eso.org/public/themes/measuring/), [IAU 2015 Resolution B2, note 4](https://iauarchive.eso.org/static/resolutions/IAU2015_English.pdf).

Integrating Planck's blackbody spectrum gives

\[
 u_\gamma=a_{\rm rad}T^4,\qquad
 a_{\rm rad}=\frac{8\pi^5 k_B^4}{15h_{\rm P}^3c^3}=\frac{4\sigma_{\rm SB}}{c}.
\]

This gives (a_{\rm rad}=7.56573325028\times10^{-16}\,\mathrm{J\,m^{-3}\,K^{-4}}), photon energy density (u_\gamma=4.17480091993\times10^{-14}\,\mathrm{J\,m^{-3}}), and mass density (u_\gamma/c^2=4.64509247757\times10^{-31}\,\mathrm{kg\,m^{-3}}). The (4\sigma_{\rm SB}/c) identity was checked against the CODATA expression for the Stefan–Boltzmann constant.

With (H_{100}=100\,\mathrm{km\,s^{-1}\,Mpc^{-1}}=3.24077928944\times10^{-18}\,\mathrm{s^{-1}}),

\[
 \omega_\gamma\equiv\Omega_\gamma h^2
 =\frac{u_\gamma/c^2}{3H_{100}^2/(8\pi G)}
 =\frac{8\pi G u_\gamma}{3c^2H_{100}^2}
 =2.472975328714088\times10^{-5}.
\]

For the requested massless-radiation reference, I used the effective-density relation

\[
 \frac{\rho_{\rm rel}}{\rho_\gamma}
 =\frac78 N_{\rm eff}\left(\frac4{11}\right)^{4/3},\qquad
 \omega_r=\omega_\gamma\left[1+\frac78 N_{\rm eff}\left(\frac4{11}\right)^{4/3}\right],
\]

with (N_{\rm eff}=3.046). The factor per effective species is (0.227107317660239), so (\omega_r=4.183702725849734\times10^{-5}). Planck 2018 gives this effective-density relation and the conventional Standard Model value near 3.046; PDG 2024 Eq. (26.1) explains its nominal temperature-ratio reference and notes refined neutrino-decoupling calculations can use (N_{\rm eff}=3.044). Here (3.046) is retained as requested and represents an effective energy-density normalization, not three ideal species each assigned an exact instantaneous-decoupling temperature. [Planck 2018 overview, Eq. for (\rho_{\rm rad})](https://doi.org/10.1051/0004-6361/201833880), [PDG 2024, Neutrinos in Cosmology, Eq. (26.1)](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-neutrinos-in-cosmology.pdf).

Since (h=H_0/(100\,\mathrm{km\,s^{-1}\,Mpc^{-1}})),

\[
 \Omega_r(H_0)=\frac{\omega_r}{h^2}.
\]

The requested (H_0\in[50,90]\,\mathrm{km\,s^{-1}\,Mpc^{-1}}) screen therefore spans (\Omega_r=1.67348109034\times10^{-4}\) down to (5.16506509364\times10^{-5}). This interval is only a sensitivity screen; it is not treated as an (H_0) prior or as an inference from BAO.

## Flat-background response

DESI DR2 writes the Friedmann background with separate photon and neutrino terms (its Eq. 6), the flat transverse distance as (D_M=(c/H_0)\int dz/E(z)) (Eq. 4), (D_H=c/H(z)) (Eq. 5), and CPL dark-energy evolution in Eqs. 8–10. Its Eq. 6 retains the redshift-dependent massive-neutrino density (\Omega_\nu\rho_\nu(z)/\rho_{\nu,0}). The present calculation instead adds a massless-radiation reference term, which is why this check is a low-redshift screening comparison rather than the DESI physical neutrino model. [DESI DR2 Results II, Eqs. (4)–(10)](https://doi.org/10.1103/tr6y-kpc6).

At fixed (\Omega_m,w_0,w_a), flat closure changes (\Omega_{\rm de}) from (1-\Omega_m) to (1-\Omega_m-\Omega_r). Put (x=1+z) and

\[
 f(z)=x^{3(1+w_0+w_a)}\exp\!\left(-\frac{3w_a z}{x}\right),\quad
 E^2(z)=\Omega_mx^3+\Omega_rx^4+(1-\Omega_m-\Omega_r)f(z).
\]

At the zero-radiation reference, (Q(z)\equiv\partial E^2/\partial\Omega_r=x^4-f(z)). For (I(z)=\int_0^z dz'/E(z')),

\[
 I'=-\frac12\int_0^z\frac{Q(z')}{E^3(z')}dz',\qquad
 (D_M/r_d)'=\alpha I',\qquad
 (D_H/r_d)'=-\frac{\alpha Q(z)}{2E^3(z)}.
\]

Here primes denote derivatives with respect to (\Omega_r), evaluated at zero radiation. With (D_V/r_d=[z(D_M/r_d)^2(D_H/r_d)]^{1/3}),

\[
 (D_V/r_d)'=\frac{D_V}{r_d}\left[\frac23\frac{I'}{I}-\frac{Q(z)}{6E^2(z)}\right],\qquad
 \Delta(D_i/r_d)=\Omega_r(D_i/r_d)' + O(\Omega_r^2).
\]

The BAO amplitude (\alpha=c/(H_0r_d)) is kept independent and fixed at each supplied seed fit point. It is not converted into an (H_0) or sound-horizon estimate.

## Numerical comparison

The reference parameters below are the reported exploratory fit points in `experiments/campaign_seed_bao/result.json`. The radiation contribution was evaluated at fixed seed parameters, with no radiation refit.

| Screen point | \(\alpha\) | \(\Omega_m\) | \(w_0\) | \(w_a\) |
|---|---:|---:|---:|---:|
| flat LCDM | 29.52463323 | 0.29746181 | -1 | 0 |
| flat wCDM | 30.05410706 | 0.29769348 | -0.91189922 | 0 |
| flat CPL | 32.40031269 | 0.37395105 | -0.3 | -2.30334177 |

The CPL seed point sits on its (w_0=-0.3) upper bound; the fixed-point conclusion therefore applies to that supplied screen point and does not establish a bound-wide CPL result.

| Model | Largest row shift / marginal \(\sigma_i\), over 41-point \(H_0\) screen | Largest full-covariance displacement \(\sqrt{\Delta d^T C^{-1}\Delta d}\) | Largest \(|\Delta\chi^2|\) at the same fit point |
|---|---:|---:|---:|
| LCDM | 0.07453 | 0.12761 | 0.03963 |
| wCDM | 0.07308 | 0.12258 | 0.00411 |
| CPL | 0.06237 | 0.10430 | 0.00843 |

The largest measurement-space shifts occur at the low endpoint (H_0=50\); (\Omega_r\propto H_0^{-2}). At (H_0=90), the corresponding maximum marginal shifts are 0.02302, 0.02257, and 0.01926 sigma, while full-covariance displacements are 0.03941, 0.03786, and 0.03221. Across the whole screen and all three fit points, the radiation-induced prediction movement remains below 0.075 of any row's marginal uncertainty and below 0.128 in the full-covariance metric. The fixed-point change in data (\chi^2) is below 0.04 in absolute value; its sign depends on alignment with the existing residuals.

The first-order prediction shifts agree with the exact finite-(\Omega_r) evaluations: at (H_0=50), the largest exact-minus-linear difference is below (9.8\times10^{-5}\) marginal sigma and relative (L_2) error in the shift vector is below (7.9\times10^{-4}). The independent zero-radiation predictions match the campaign seed vector within (7.2\times10^{-14}); changing Gauss–Legendre quadrature from 64 to 512 nodes changes the predictions by at most (5.0\times10^{-14}).

**Interpretation:** radiation neglect is acceptable at the supplied LCDM, wCDM, and CPL best-fit points for this 13-row, (z\le2.33), background-only covariance screen. This does not establish acceptability everywhere in the broad fit boxes, at earlier redshifts, for a sound-horizon calculation, or in precision inference using a physical massive-neutrino background. The conclusion is deliberately conditional on the specified massless-neutrino reference and fixed (\alpha\)/shape evaluation.

## Reproduction record

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 .venv/bin/python work/theory2/check_radiation.py
```

The calculation used Python 3.12.3, NumPy 2.5.3, CPU float64, one BLAS thread, and no GPU. The JSON output stores SI intermediate values, input SHA-256 hashes, both first-order and exact distance-shift vectors for all 13 rows, each row's marginal-standardized shift, the full-covariance metric, fit-point (\Delta\chi^2), and quadrature checks.

SHA-256 inputs:

| Input | SHA-256 |
|---|---|
| `context/data/desi_dr2_mean.txt` | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `context/data/desi_dr2_cov.txt` | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |
| `scripts/background_bao.py` | `34bdf1d81397b39b740faa3e2f4a79748d5fa9d3bfdd8303619d631dbf8aa2e2` |
| `experiments/campaign_seed_bao/result.json` | `92662f3ee42500d2be34d9f87b485fe56c82e9550fa674e42f4dec2bf2cab1f0` |

The standalone checker hash and complete raw result are recorded in [radiation_check.json](radiation_check.json); the executable derivation is [check_radiation.py](check_radiation.py).

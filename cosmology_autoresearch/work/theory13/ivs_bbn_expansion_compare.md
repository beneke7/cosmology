# IVS expansion-rate effect during standard BBN

**Finding.** For the archived IVS screen point, the largest absolute expansion-rate change relative to flat Lambda-CDM with identical present-day component densities is 7.2650 ppm over 1e8 <= z <= 1e10. The shift is negative for Gamma/H0=-0.4666647299: H_IVS is lower than H_Lambda. It is far below a 1% BBN-era expansion change, so this background point does not invalidate use of a standard-SBBN baryon prior on expansion-rate grounds. This is a homogeneous expansion comparison only; no abundance response or likelihood was computed.

## Setup and sign

This independently reimplements the theory5 screen point Omega_m0=0.3869710077, g=Gamma/H0=-0.4666647299, w_v=-1 over all 16 cases H0=[55.0, 65.0, 75.0, 85.0] km s^-1 Mpc^-1 and f_b=[0.1, 0.15, 0.2, 0.25]. For each case, IVS and Lambda-CDM share Omega_b0, Omega_c0, Omega_r0 and flat-closure Omega_v0. Photons use Omega_gamma h^2=2.4728e-05; neutrinos use N_eff=3.044 and the massless approximation at all redshifts. This extends the source paper's nonzero-mass setup with the explicit approximation requested here.

The source sign convention is dot(rho_v)=Gamma*rho_v and dot(rho_c)+3H*rho_c=-Gamma*rho_v. With u=ln(1+z), q_c=a^3 rho_c/rho_crit, v=rho_v/rho_crit, and E=H/H0, the direct equations are

    dq_c/du = g*a^3*v/E
    dv/du   = -g*v/E
    E^2     = (Omega_b0+q_c)*(1+z)^3 + Omega_r0*(1+z)^4 + v.

Thus negative g makes v increase backward and q_c decrease backward from Omega_c0. The Lambda reference is E_Lambda^2=Omega_m0*(1+z)^3+Omega_r0*(1+z)^4+Omega_v0. The reported delta-H uses an analytic difference of E^2 terms to avoid subtracting two radiation-dominated totals directly.

An independent second integration uses q_d=a^3(rho_c+rho_v)/rho_crit and v, with dq_d/du=-3*a^3*v and the same dv/du. It reconstructs q_c=q_d-a^3*v and H from total dark-sector density. The g=0 solution is the exact Lambda null in both formulations.

## Results

Ranges below are across all 16 H0/f_b cases. Delta-H is (H_IVS/H_Lambda)-1.

| z | min delta-H/H | max delta-H/H | largest magnitude (ppm) |
|---:|---:|---:|---:|
| 1e+08 | -7.2650021e-06 | -3.0412365e-06 | 7.26500 |
| 1e+09 | -7.2654156e-07 | -3.0413089e-07 | 0.72654 |
| 3e+09 | -2.4218154e-07 | -1.0137714e-07 | 0.24218 |
| 1e+10 | -7.2654569e-08 | -3.0413162e-08 | 0.07265 |
| 1.24748e+08 | -5.8238270e-06 | -2.4379216e-06 | 5.82383 |

On the 401-point logarithmic BBN bracket, delta-H/H ranges from -7.2650021e-06 to -3.0413162e-08; maximum |Omega_v/Omega_r| at requested target redshifts is 1.647e-28. This shows the vacuum component itself is negligible relative to radiation at BBN. The rate difference is dominated by the altered early comoving-CDM normalization, still strongly radiation-suppressed.

Independent direct-vs-total agreement is 8.216e-20 in delta-H/H. The maximum component-ratio versus stable-difference check is 1.531e-16. Across all 16 cases, the g=0 null has max |delta-H/H|=0.000e+00 (PASS); default-to-tight tolerance refinement changes target delta-H/H by at most 1.948e-20 (PASS). The full 401-point-per-case results, all checks, runtime, environment, and input hashes are in the companion JSON.

## Scope and disposition

The comparison supports using the standard-SBBN prior for this one screened IVS point with respect to the homogeneous H(T) effect, because |delta-H/H| stays below 7.2650 ppm over the declared BBN redshift bracket. It does not prove the prior valid for the full IVS posterior or different interaction points. It also does not account for electron-positron entropy transfer in the massless-radiation approximation, evolve neutrino masses, compute light-element abundances, or test perturbations. Reproducing PRIMAT with this modified expansion would be required to bound the actual abundance shift. No PRIMAT install, download, GPU, or shared-file edit was used.

Reproduce from the project root with:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
      .venv/bin/python scripts/run_bounded.py --seconds 300 -- \
      .venv/bin/python work/theory13/ivs_bbn_expansion_compare.py

Runtime: 3.284 s. Script SHA-256: 6ca7287e9c86fbfb090e9123f5cd286a8887706390d30f6e3e96d6cd50893431.

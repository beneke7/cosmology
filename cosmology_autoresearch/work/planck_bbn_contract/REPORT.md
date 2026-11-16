# Planck 2018 helium/BBN contract and CAMB table audit

Audit date: 2026-09-25. Status: reproduced local interpolation, independently checked interpolation, and literature/code audited. No root config, shared ledger, or campaign source catalog was edited. No Planck likelihood was evaluated.

## Conclusion

CAMB 1.6.6 contains a PArthENoPE table that reproduces Planck VI's published helium baseline to its quoted precision. The matching choice is `PArthENoPE_880.2_marcucci.dat`, evaluated at DeltaNeff = 0 for Planck's standard Neff = 3.046. CAMB 1.6.6 defaults instead to a PRIMAT 2021 table, so leaving `bbn_predictor` unspecified is not a PArthENoPE calculation.

At omega_b h2 = Omega_b h2 = 0.02236 and TCMB = 2.7255 K, the selected Marcucci table gives `Yp_BBN = 0.246724678`; CAMB's corresponding helium mass fraction is `YHe = 0.245394698`. These reproduce Planck's quoted 0.24672 and 0.2454. The table is a strong numerical match to Planck's published reference, but available metadata does **not** prove bitwise or exact executable-version equivalence to PArthENoPE 1.1: CAMB identifies it as a May 2017 PArthENoPE table, not explicitly as generated with version 1.1. No original Planck interpolation-grid hash was found for direct binary comparison.

Configuration implication for the root orchestrator: set CAMB `extra_args.bbn_predictor: PArthENoPE_880.2_marcucci.dat`, `nnu: 3.046`, and `standard_neutrino_neff: 3.046`. The last setting matters because the installed CAMB 1.6.6 default standard value is 3.044. If only `nnu: 3.046` is set, CAMB passes DeltaNeff = 0.002 into the BBN interpolator and its massive-neutrino density normalization also retains the 3.044 reference. This is a small numerical deviation, but it is not the stated Planck contract.

## Literature contract

Planck VI section 2.1 states that when the helium mass fraction is not sampled independently, it is derived by interpolating an updated BBN grid computed with PArthENoPE 1.1 (the paper says version 2.0 has identical physical inputs/results), with fixed neutron lifetime tau_n = 880.2 s and the uncertainty neglected. Its base model uses Neff = 3.046. At the Planck Lambda-CDM reference it reports Yp approximately 0.2454 as the helium **mass fraction** and Yp_BBN approximately 0.2467 as the **nucleon-number fraction**. These are different definitions and must not be interchanged.

Planck VI section 7.6.1 and Fig. 38 resolve the nuclear-rate choice. Its standard PArthENoPE prediction (green, case b) uses the theoretical d(p,gamma)3He rate of Marcucci et al. (2016). The alternate PArthENoPE case (a) uses the experimental Adelberger et al. rate. Consistently, Planck Eq. (72) gives y_DP = 10^5 D/H = 2.455 for case (b), and 2.587 for case (a). CAMB's Marcucci table yields 10^5 D/H = 2.454710, while its `standard` table yields 2.587454. This independent deuterium diagnostic identifies which of CAMB's near-identical-in-helium tables matches Planck's standard branch.

For uncertainty context, Planck VI section 7.6.1 quotes sigma(tau_n) = 1 s corresponding to sigma(Yp_BBN) = 3e-4. That uncertainty is not interpolated from CAMB's PArthENoPE files: their headers say errors are not given. Planck held tau_n fixed for baseline parameter inference, while separately discussing theoretical uncertainty in its abundance predictions.

Planck VI's 2021 corrigendum corrects only the 95th-percentile bounds on early-time optical depth in Eqs. (87a,b), not the baseline BBN prescription.

## Reproduced values and checks

Inputs: omega_b h2 = 0.02236, Neff = 3.046, CAMB reference `standard_neutrino_neff = 3.046`, and TCMB = 2.7255 K. The baryon density omega_b h2 is dimensionless and lies inside the table range 0.00499-0.03993. The table DeltaN coordinate is 0 for this Planck standard model. CAMB rescales the table density coordinate by `(2.7255 K / TCMB)^3`; the default temperature equals the table header, so no correction occurs.

| CAMB table | Yp_BBN (CAMB spline) | Yp table mass fraction | CAMB YHe from number fraction | 10^5 D/H | Interpretation |
|---|---:|---:|---:|---:|---|
| `PArthENoPE_880.2_marcucci.dat` | 0.246724678 | 0.245394694 | 0.245394698 | 2.454710 | Matches Planck standard case (b), Marcucci rate |
| `PArthENoPE_880.2_standard.dat` | 0.246721416 | 0.245391433 | 0.245391448 | 2.587454 | D/H matches Planck alternate case (a), experimental Adelberger rate |
| `PRIMAT_Yp_DH_ErrorMC_2021.dat` | 0.247182859 | 0.245855905 | 0.245851214 | 2.450054 | Installed CAMB default; not Planck's PArthENoPE baseline |

The Planck-relevant Marcucci output agrees with Planck Eq. (70a), Yp_BBN = 0.24672, to 4.7e-6, and its YHe rounds to Planck's 0.2454. The CAMB default PRIMAT table is higher than the selected PArthENoPE mass fraction by 4.565e-4, consistent in scale with Planck's stated approximately 5e-4 inter-code shift. Planck says the other parameter results remain consistent to well within 0.1 sigma. Do not silently call the CAMB default a Planck BBN reproduction.

An independently coded SciPy bilinear interpolator over the printed tables returns the same rounded Yp_BBN = 0.24672. At this off-grid point it differs from CAMB's bicubic spline by 2.68e-6 for the Marcucci table. This interpolation/table-rounding sensitivity is negligible compared with the approximately 3e-4 theoretical helium uncertainty quoted by Planck. Direct `camb.set_params` calls verify that the selected table name and `standard_neutrino_neff` reach CAMB's cosmology setter.

If Neff = 3.046 is passed but `standard_neutrino_neff` is omitted, CAMB 1.6.6's 3.044 default instead yields YHe = 0.245421469 for the Marcucci table (a +2.6771e-5 shift) and changes the massive-neutrino density normalization. Set the reference explicitly rather than relying on a package default that differs from the 2018 paper.

## Reproducibility and local provenance

The exact open paper used was downloaded from arXiv v4 and retained at [1807.06209v4.pdf](1807.06209v4.pdf), SHA-256 `8e172730faf07c9f4ff3fdcc7043f76ed67df6f76066d47df30d693025b6ce77` (9,317,371 bytes). The deterministic audit script is [reproduce_bbn_check.py](reproduce_bbn_check.py); it reports local product hashes, grid ranges, CAMB spline values, an independent bilinear check, and direct CAMB parameter-setting results.

Run from the project root:

```bash
source cosmology_autoresearch/.venv/bin/activate
python cosmology_autoresearch/work/planck_bbn_contract/reproduce_bbn_check.py
```

Environment used: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, CAMB 1.6.6, Cobaya 3.6.2. Key installed inputs and hashes:

| Artifact | SHA-256 |
|---|---|
| CAMB `bbn.py` | `879304e01c5df30eb38707e29ffc7cf7c6595f452a36382abcfbfca5189b78ec` |
| CAMB `model.py` | `5fafd2b1086369960f15ebc6f3dee37bf56ab707979e5de7d8c8e038a1a58fb8` |
| CAMB `constants.py` | `43e40ba181fc88a033e851f0ee86a345777b9bb2437bb09e84261d8acd28090d` |
| `PArthENoPE_880.2_marcucci.dat` | `994c68f11ecbd33ecf6e631a0432ffc4b80068bfe972ab98fd235420aad67d8e` |
| `PArthENoPE_880.2_standard.dat` | `99e82051a70b8a73eed112959a5bf42bc1ea9e0e08860d65ae75fe97ce09b187` |
| `PRIMAT_Yp_DH_ErrorMC_2021.dat` | `ea5adce061720b937d8abda3a04a384aedaab3168dbf17414ff600cc91a7160c` |

The exact local tables are under `/home/v/proj/bene/cosmology/cosmology_autoresearch/.venv/lib/python3.12/site-packages/camb/`. CAMB's local `bbn.py` lines 14-16 establish the PRIMAT default; lines 71-74 name the PArthENoPE 2017 and Marcucci tables. The bundled table headers say “Calculated from PArthENoPE by Ofelia Pisanti May 2017”, “neutron lifetime of 880.2 s and CMB temperature 2.7255 K”, identify `Yp^BBN` as a nucleon-number fraction and `Yp` as the mass fraction, and say errors are not given. `model.py` lines 605-625 and 684-690 expose and apply `bbn_predictor`, `standard_neutrino_neff`, DeltaN subtraction, and TCMB scaling. Cobaya's installed `theories/camb/camb.py` lines 949-960 passes `extra_args` to `camb.set_params`.

## Sources

1. Planck Collaboration, *Planck 2018 results. VI. Cosmological parameters*, arXiv v4, especially section 2.1, section 7.6.1, Fig. 38 and Eqs. (70a), (72): [arXiv HTML v4](https://arxiv.org/html/1807.06209v4), [exact PDF v4](https://arxiv.org/pdf/1807.06209v4), [journal DOI](https://doi.org/10.1051/0004-6361/201833910).
2. Planck Collaboration corrigendum, [A&A 652 C4 (2021), DOI](https://doi.org/10.1051/0004-6361/201833910e); the correction concerns the early-time optical-depth bounds in Eqs. (87a,b).
3. CAMB upstream pinned tag [1.6.6 `camb/bbn.py`](https://raw.githubusercontent.com/cmbant/CAMB/1.6.6/camb/bbn.py), [PArthENoPE Marcucci table](https://raw.githubusercontent.com/cmbant/CAMB/1.6.6/camb/PArthENoPE_880.2_marcucci.dat), [PArthENoPE standard table](https://raw.githubusercontent.com/cmbant/CAMB/1.6.6/camb/PArthENoPE_880.2_standard.dat), [CAMB `model.py`](https://raw.githubusercontent.com/cmbant/CAMB/1.6.6/camb/model.py), and [CAMB repository/tag](https://github.com/cmbant/CAMB/tree/1.6.6).
4. PArthENoPE official site: [parthenope.na.infn.it](https://parthenope.na.infn.it/).

## Limitations / disposition

- This reproduces the local deterministic interpolation and CAMB's number-to-mass-fraction mapping, not a fresh PArthENoPE nuclear-network run and not a complete Planck CMB likelihood or chain.
- The CAMB table's year, author, tau_n, and TCMB are explicit. Its metadata does not explicitly identify PArthENoPE executable version 1.1, and the paper provides no hash for its original grid. Therefore “matches the Planck published benchmark” is supported; “the exact same executable/table bits” is not.
- The independent bilinear check uses a less accurate, separately implemented interpolation over rounded grid values. It is a sanity check, not an independent BBN calculation.
- The project Planck configuration remains **not evaluated** by this task. It must explicitly select the Marcucci table and standard neutrino reference before the likelihood smoke test.

# Planck PR3 flat-\(\Lambda\)CDM reference for the no-reconstruction-lensing combination

Audit date: 2026-09-25. Status: literature/config audit only. No package was installed, no Planck likelihood was unpacked or evaluated, and this task did not download any file. All writes for this task are in this directory.

## Finding

The exact published reference is Planck Collaboration, *Planck 2018 results. VI. Cosmological parameters*, A&A 641, A6 (2020), Table 2, fourth data-combination column, labelled TT,TE,EE+lowE. The electronic paper numbers it as Table 2; it is on PDF p. 15. Its 68% intervals are the no-lensing-reconstruction result for full Plik high-\(\ell\) TT+TE+EE plus low-\(\ell\) Commander TT and SimAll EE. Planck §2.2.4 defines that likelihood label explicitly. The paper's §2.2.3 names Commander for low-\(\ell\) temperature and SimAll's EE-only likelihood for low-\(\ell\) polarization.

Here “no lensing” excludes the separate Planck CMB-lensing-reconstruction likelihood. The theoretical TT/TE/EE spectra remain physically lensed in standard \(\Lambda\)CDM. Do not switch off physical CMB lensing in CAMB or add the separate Planck lensing likelihood.

Table 2 reports marginalized 68% intervals, not a likelihood maximum. Do not substitute Table 1: its published Plik best fit and Plik posterior column are for TT,TE,EE+lowE+lensing. The fourth column of Table 2 is the exact match.

## Published result and locators

Primary source: [Planck 2018 results VI, A&A 641, A6 (2020), DOI 10.1051/0004-6361/201833910](https://doi.org/10.1051/0004-6361/201833910); [arXiv v4 HTML](https://arxiv.org/html/1807.06209v4), dated 2021-08-09.

- Table 2, fourth column, TT,TE,EE+lowE, 68% intervals (PDF p. 15; in arXiv v4 HTML, rows at lines 509–536). The rows for six sampled parameters are \(\Omega_bh^2\), \(\Omega_ch^2\), \(100\theta_{\rm MC}\), \(\tau\), \(\ln(10^{10}A_s)\), and \(n_s\).
- §2.1, “Theoretical model” (arXiv v4 HTML lines 259–274): base model and neutrinos, flat priors and their relation to Planck XIII Table 1, the 2018 \(A_s\) prior change, BBN helium convention, CAMB/CosmoMC versions.
- §2.2.3, “The low-\(\ell\) likelihood” (HTML lines 404–418): lowE uses SimAll EE; low-\(\ell\) TT uses Commander.
- §2.2.4, “Likelihood notation” (HTML lines 419–421): exact TT,TE,EE+lowE definition and exclusion of low-\(\ell\) TE.
- §2.2.5 on uncertainty construction (HTML line 423 onward); also §2.2.1 on Plik modelling and §3/Table 1 on high-\(\ell\) model sensitivity.
- Equation (2) in §2.2.3 is \(\tau=0.0506\pm0.0086\) for lowE alone, not the joint result. It is recorded here only to disambiguate it from the joint Table 2 result.
- Prior bounds carried forward from Planck 2013 results XVI, Table 1 (the 2018 paper explicitly says so in §2.1): [Planck 2013 XVI, A&A 571, A16, arXiv HTML](https://arxiv.org/html/1303.5076), Table 1 and surrounding text (HTML lines 271–287, 374–385). The 2018 paper replaces the earlier amplitude prior with its stated log-amplitude interval.

| Sampled base parameter | Planck Table 2: 68% marginal interval | Units / interpretation |
|---|---:|---|
| \(\omega_b=\Omega_bh^2\) | \(0.02236\pm0.00015\) | Dimensionless physical baryon density |
| \(\omega_c=\Omega_ch^2\) | \(0.1202\pm0.0014\) | Dimensionless physical cold-dark-matter density |
| \(100\theta_{\rm MC}\) | \(1.04090\pm0.00031\) | Dimensionless CosmoMC approximation to the acoustic angular scale |
| \(\tau\) | \(0.0544^{+0.0070}_{-0.0081}\) | Dimensionless Thomson optical depth |
| \(\ln(10^{10}A_s)\) | \(3.045\pm0.016\) | Dimensionless scalar amplitude at \(k_*=0.05\,{\rm Mpc}^{-1}\) |
| \(n_s\) | \(0.9649\pm0.0044\) | Dimensionless scalar tilt |

Useful derived checks from the same Table 2, same fourth column: \(H_0=67.27\pm0.60\ {\rm km\,s^{-1}\,Mpc^{-1}}\); \(\Omega_\Lambda=0.6834\pm0.0084\); \(\Omega_m=0.3166\pm0.0084\); \(\sigma_8=0.8120\pm0.0073\); \(S_8=\sigma_8(\Omega_m/0.3)^{1/2}=0.834\pm0.016\). Table 2 says \(\Omega_m\) includes one \(0.06\) eV neutrino.

## Model, priors, and published inference setup

Planck §2.1 defines the base model as spatially flat \(\Lambda\)CDM with adiabatic scalar perturbations and a power-law primordial scalar spectrum. There are three neutrino species, represented by two massless states and one massive eigenstate at \(m_\nu=0.06\) eV; the baseline fixes standard \(N_{\rm eff}=3.046\), \(\Omega_K=0\), and \(w=-1\). The helium fraction is predicted from BBN, using PArthENoPE 1.1 and fixed neutron lifetime \(880.2\) s; the Planck mean \(Y_P\) is about 0.2454.

The paper states flat priors on the six base parameters. The 2018 §2.1 directs readers to Planck 2013 XVI Table 1 for the inherited bounds, with one explicit 2018 change. In Planck notation:

| Parameter | Prior used for this base analysis | Comment |
|---|---|---|
| \(\omega_b\) | Uniform \([0.005,0.1]\) | Planck 2013 XVI Table 1 |
| \(\omega_c\) | Uniform \([0.001,0.99]\) | Planck 2013 XVI Table 1 |
| \(100\theta_{\rm MC}\) | Uniform \([0.5,10.0]\) | Planck 2013 XVI Table 1 |
| \(\tau\) | Uniform \([0.01,0.8]\) | Planck 2013 XVI Table 1 |
| \(n_s\) | Uniform \([0.9,1.1]\) | Pivot \(k_*=0.05\,{\rm Mpc}^{-1}\), Planck 2013 XVI Table 1 |
| \(\ln(10^{10}A_s)\) | Uniform \([1.61,3.91]\) | Explicit 2018 change, §2.1 |
| \(H_0\) | Hard bound \([20,100]\ {\rm km\,s^{-1}\,Mpc^{-1}}\) | Legacy hard bound described in Planck 2013 XVI §2; \(H_0\) is derived in the six-parameter base model |

The published parameter chains used the August 2017 CAMB and July 2018 CosmoMC versions (§2.1), with marginalized densities and intervals computed using the updated GetDist adaptive-kernel estimates. The paper says the reported errors do not include all likelihood-model systematics; its discussion estimates that roughly \(0.5\sigma\) total systematic uncertainty may be more realistic. In particular, its high-\(\ell\) polarization constraints respond to calibration and likelihood-construction details.

## Exact likelihood mapping: Cobaya 3.6.2

The following links pin the public Cobaya component definitions to the v3.6.2 Git tag. They reproduce the exact three likelihood components, without a fourth lensing-reconstruction component:

| Planck component | Cobaya v3.6.2 key and YAML | PR3 .clik file named by YAML |
|---|---|---|
| Full Plik high-\(\ell\) TT+TE+EE | planck_2018_highl_plik.TTTEEE, [TTTEEE.yaml](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/planck_2018_highl_plik/TTTEEE.yaml) | baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik |
| Commander low-\(\ell\) TT | planck_2018_lowl.TT_clik, [TT_clik.yaml](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/planck_2018_lowl/TT_clik.yaml) | baseline/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik |
| SimAll low-\(\ell\) EE | planck_2018_lowl.EE_clik, [EE_clik.yaml](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/planck_2018_lowl/EE_clik.yaml) | baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik |

All three YAMLs name Planck product_id 151902. The high-\(\ell\) component has alias plikHM_TTTEEE; the low-\(\ell\) components have aliases lowl and lowE. These exact filenames and aliases make the target unambiguous: this is full Plik, not TTTEEE_lite, and not an alternate lowE implementation.

The exact v3.6.2 defaults relevant to a later run are also version-pinned:

- [Planck calibration prior](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/base_classes/planck_calib.yaml): \(A_{\rm planck}\sim\mathcal N(1,0.0025)\).
- [Plik temperature calibration defaults](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/planck_2018_highl_plik/params_calib_temp.yaml): calib_100T \(\mathcal N(1.0002,0.0007)\), calib_217T \(\mathcal N(0.99805,0.00065)\).
- [Plik polarization calibration defaults](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/planck_2018_highl_plik/params_calib_pol.yaml): A_pol=1, calib_100P=1.021, calib_143P=0.966, calib_217P=1.040 are fixed by the YAML. This is consistent with the paper's baseline map-based polarization-efficiency choice.
- [Plik TT foreground defaults](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/planck_2018_highl_plik/params_TT.yaml), [TE foreground defaults](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/planck_2018_highl_plik/params_TE.yaml), and [EE foreground defaults](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/planck_2018_highl_plik/params_EE.yaml) define component-specific nuisance freedoms. Examples: CIB index -1.3 and TE/EE dust indices -2.4 are fixed; TT Galactic-dust and TE dust amplitudes have Gaussian priors; the four TT point-source amplitudes have uniform [0,400] priors; SZ amplitudes have uniform [0,10] priors; and EE dust amplitudes are fixed to their stated central values.
- [Plik SZ combination prior](https://raw.githubusercontent.com/CobayaSampler/cobaya/v3.6.2/cobaya/likelihoods/planck_2018_highl_plik/prior_SZ.yaml): a Gaussian prior on ksz_norm + 1.6*A_sz, mean 9.5, standard deviation 3.0. It is not two independent priors.

The YAML component definitions alone do not declare the six cosmological parameter priors. A later Cobaya input must specify those explicitly and match the Planck §2.1 ranges above. Full Plik samples nuisance parameters and includes multiple foreground/calibration terms. Holding them all fixed, using Plik-lite, or dropping the SZ combination prior changes the likelihood contract. The Planck paper further says low/high-multipole correlations are neglected in the hybrid baseline (§2.2.2); preserve that original likelihood definition.

## Valid initial point and chain-product caveat

For an initial standard-Cobaya evaluation, use the no-lensing Table 2 fourth-column central values as a posterior-centered cosmological initialization:

~~~yaml
ombh2: 0.02236
omch2: 0.1202
theta_MC_100: 1.04090
tau: 0.0544
logA: 3.045
ns: 0.9649
~~~

These names follow the usual Cobaya/CAMB parameterization; check the provider's spelling in the eventual pinned config. Set model constants to the Planck baseline (\(\Omega_K=0\), one \(0.06\) eV massive neutrino, standard radiation, BBN helium, \(w=-1\)). Initialize and sample the full Plik nuisance block using the pinned v3.6.2 definitions. This point lies near the exact combination's published posterior and inside the documented broad prior region. It is not the maximum-likelihood point: Table 2 gives marginalized intervals only. It supplies no official nuisance vector for scoring. Do not report a raw log-likelihood at this initial point as Planck's published best fit.

Planck's public PR3 archive lists COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip as baseline \(\Lambda\)CDM parameter chains based on the named combination (62.3 MB). These chains and their GetDist companions are useful for checking posterior locations or choosing a start, after verifying exact chain metadata. They are derived posterior samples, not the raw TT/TE/EE and low-\(\ell\) likelihood implementation and not a likelihood for rescoring arbitrary candidate parameters. Do not treat their samples or their stored likelihood columns as a new independent data factor or as a substitute for evaluating the PR3 likelihood. The larger full-grid and lensing-only chain archives are also analysis outputs; lensonly is a different observation and the full grid mixes many combinations.

The exact official 62.3 MB chain archive was not fetched or inspected in this task, to honor the 5 MiB download ceiling. The official [PR3 ancillary-data product index](https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/) and [ESA Planck cosmology product metadata](https://esdcdoi.esac.esa.int/doi/html/data/astronomy/planck/Cosmology.html) identify the products. The separate LAMBDA prose page contains an ambiguous summary phrase about “including lensing”; use the explicit product filename and the Planck Table 2 combination label when distinguishing the no-reconstruction-lensing chain set.

## Local artifacts and boundary

The shared campaign staging area already contained the official raw likelihood-data archive at audit time: [COM_Likelihood_Data-baseline_R3.00.tar.gz](../planck_stage/COM_Likelihood_Data-baseline_R3.00.tar.gz), 60,323,470 bytes, SHA-256 0b73171e3acc671c28184466a45485a2d1c1d93676b832abdfe688c7b04024e6. This task did not fetch, extract, or otherwise modify that file. It is the data archive, not a scored likelihood. No local Planck chains or local Cobaya component YAML files were found or created. Component YAML content was read from the pinned remote raw-file URLs above; no local hashes apply.

The raw-data route is described by the [Cobaya Planck likelihood documentation](https://cobaya.readthedocs.io/en/v3.6.2/likelihood_planck.html) and the versioned YAMLs above. The gate's environment audit says Cobaya 3.6.2, clipy-like, and the Planck likelihood data are not yet installed in the Python environment. That separate integration gate therefore remains; this reference task establishes the target and a valid starting cosmology, not a reproduced numerical likelihood or a new posterior.

## Source list

1. Planck Collaboration, 2020, *Planck 2018 results. VI. Cosmological parameters*, A&A 641 A6, DOI [10.1051/0004-6361/201833910](https://doi.org/10.1051/0004-6361/201833910), [arXiv v4 HTML](https://arxiv.org/html/1807.06209v4).
2. Planck Collaboration, 2014, *Planck 2013 results. XVI. Cosmological parameters*, A&A 571 A16, [arXiv HTML](https://arxiv.org/html/1303.5076), Table 1 and §2.
3. [Cobaya v3.6.2 Planck likelihood documentation](https://cobaya.readthedocs.io/en/v3.6.2/likelihood_planck.html); the three pinned component YAMLs above and their nuisance-default files.
4. [ESA Planck Legacy Archive](https://www.cosmos.esa.int/web/planck/pla), [ESA PR3 cosmology metadata](https://esdcdoi.esac.esa.int/doi/html/data/astronomy/planck/Cosmology.html), and [NASA/IPAC PR3 ancillary product index](https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/).

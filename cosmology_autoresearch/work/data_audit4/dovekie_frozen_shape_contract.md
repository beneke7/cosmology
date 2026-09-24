# DES-Dovekie frozen-shape contract audit

As of 2026-09-25. Status: **data contract verified; proceed with a frozen-shape, conditional predictive check through an adapter**. This audit checks the pinned Hubble diagram, packed precision, and likelihood source. It does not fit or score any cosmological shape.

## Pinned sources and local identity

The official `des-science/DES-SN5YR` repository is pinned here to commit [`c9a4fcafc4cbd19bd750dee47fc76194a45c181f`](https://github.com/des-science/DES-SN5YR/tree/c9a4fcafc4cbd19bd750dee47fc76194a45c181f), authored and committed 2026-01-28 10:43:40 UTC. This is a commit pin, not a release tag. I queried the official commit tree metadata (no data payload download) and compared its Git blob IDs with the local files:

| Contract item | Pinned upstream path | Upstream Git blob | Local SHA-256 |
|---|---|---|---|
| Ordered Hubble diagram | [`4_DISTANCES_COVMAT/DES-Dovekie_HD.csv`](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv) | `f80ec4e2795edcbf3442f460c539bea56226027a` | `2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf` |
| Full STAT+SYS inverse covariance | [`4_DISTANCES_COVMAT/STAT+SYS.npz`](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/STAT%2BSYS.npz) | `4289666487f427782ec81c327ae7f7741f0f5fe5` | `ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b` |
| Release distance/covariance README | [`4_DISTANCES_COVMAT/README.md`](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/README.md) | `d5900bc1a4a2ce345012c878c929204cef5b5c52` | Local snapshot `work/data_audit/Dovekie_distances_readme.md` has that same Git blob |
| CosmoSIS module | [`5_COSMOLOGY/Dovekie_cosmosis_likelihood.py`](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/5_COSMOLOGY/Dovekie_cosmosis_likelihood.py) | `b7142093d633bf62281a2253d85ccc54db48431a` | `78526c0e6013619e3586a5509be504f0bc5f0461d1e55a3d5bec0154353645ff` |

The upstream tree reports sizes of 148,002 bytes for the Hubble diagram and 6,244,951 bytes for the NPZ. The local copies’ Git blob IDs match those exact tree entries, establishing the row/covariance source mapping without re-downloading the data. The release README says its covariance files are already inverse covariance matrices and warns that the cosmology Hubble diagram ordering is the one paired with the covariance; its text has a filename typo (`DES-SN55YR_HD.csv`), while the actual paired release asset is `DES-Dovekie_HD.csv`. Do not substitute the metadata table or sort/join rows by any field.

The paper source is [arXiv:2511.07517v3](https://arxiv.org/abs/2511.07517v3), submitted 2025-11-10 and last revised 2026-03-27. The local PDF SHA-256 is `cd3dd615c36afe0d45bcdfabd9c31611f0ca511e325997f93c3eb1f5a9744eb6`. The paper reports the sample overlap with DES-SN5YR and describes its DESI DR2 BAO combination. The upstream README and likelihood links above are commit-pinned; the upstream README is dated by the commit, not by the mutable ReadTheDocs `latest` alias.

## Verified row and covariance contract

The local release Hubble diagram is a SNANA whitespace table despite its `.csv` suffix: `VARNAMES:` followed by `SN:` records, with no comma delimiters. Parsing in file order gives 1,820 records and 1,820 unique `CID`s; first/last IDs are `Gaia16agf` / `1257587`. The SHA-256 of original-order CIDs joined by LF, with a trailing LF, is `b7d5c6ad8dfdadf1e12443852b14006c0d9d0ceb4bcb46efb378f058ce5aa5a3`. `zHD` is nondecreasing in this file. All 1,820 rows pass the official module’s `zHD > 0` filter, so that filter preserves the full release ordering.

`zHD` is the Hubble-diagram/CMB-frame redshift with peculiar-velocity correction; `zHEL` is heliocentric redshift. Their ranges are 0.02509–1.14418 and 0.02385–1.14501. Redshifts are dimensionless. `MU` is a distance modulus in magnitudes (release README convention assumes H0 = 70); the observed range is 34.98002–44.62596 mag. `MUERR` is a magnitude-error diagnostic column (0.0409–468.0108 mag in the file), not an extra term to add to the released STAT+SYS covariance. `MUERR_SYS` is also diagnostic. The likelihood returns the vector `MU` and uses the packed STAT+SYS product for its Gaussian covariance.

`STAT+SYS.npz` has keys `nsn` (one int64), `cov` (1,657,110 float32 values), and scalar `allow_pickle` (bool). The packed length is (1820\cdot1821/2), the upper triangle of a 1,820×1,820 **precision** matrix (P=C^{-1}). I unpacked it in the official source order, reflected the upper triangle to the lower, and did not reorder or add diagonal errors. Checks with one BLAS/OMP/MKL thread:

- $P$ is exactly symmetric; Cholesky succeeds.
- Eigenvalue range of $P$: $4.5654988544\times10^{-6}$ to $271.8406992$; $\kappa_2(P)=5.95423869\times10^7$.
- Solving $PC=I$ gives a 1,820×1,820 covariance $C$; its maximum antisymmetric roundoff is $1.11\times10^{-16}$. After symmetrizing only that inversion roundoff, Cholesky succeeds. Thus the full released covariance is positive definite; use factorizations/solves for block operations because the condition number is large.

The official module reconstructs the packed upper-triangle `Covtot_inv`, reflects it, inverts it to covariance, applies the same `zHD > 0` mask, and returns it for the Gaussian parent to invert. This confirms both that the archive stores the inverse and that the CSV rows and matrix dimensions align. Since every row survives the mask, no covariance subsetting is needed for the official contract. The `MUERR` values must not be added on top of $C$.

## Likelihood semantics and source mismatch

The pinned module’s theory vector is

\[
g_i(\theta)=5\log_{10}\!\left[(1+z_{{\rm HD},i})(1+z_{{\rm HEL},i})D_A(z_{{\rm HD},i};\theta)\right]+25,
\]

with `D_A` interpolated at `zHD`. For residual $d=g-y$, data vector $y=\texttt{MU}$, and precision $P=C^{-1}$, it analytically integrates one additive magnitude offset $M$:

\[
-2\log L_{\rm marg}=d^TPd-\frac{(\mathbf1^TPd)^2}{\mathbf1^TP\mathbf1}
+\log\!\left(\frac{\mathbf1^TP\mathbf1}{2\pi}\right)+\log|C|,
\]

where the fixed Gaussian normalization is supplied by the CosmoSIS Gaussian base class. The source comments explicitly say that $M$ is fully degenerate with $H_0$. The `H0=70` convention in released `MU` is therefore an overall magnitude gauge: use a fixed common H0 convention and allow one free, integrated constant $M$; this SN score does not measure H0.

There are two implementation mismatches to account for before using the published module with this release:

1. Its defaults point to `DES-SN5YR_HD.csv` and `STAT+SYS.txt.gz`, not the paired `DES-Dovekie_HD.csv` and `STAT+SYS.npz` assets.
2. It calls `Table.read(..., format='ascii.csv')`, while the actual paired `.csv` is SNANA whitespace text. Parse that format explicitly (or make an order-preserving conversion) and pass the Dovekie paths. The source is exact upstream code, but the source/data pair is not directly runnable unchanged as a normal CSV ingestion.

Thus the release likelihood algebra and data are sufficiently specified for a faithful adapter; this audit did not run CosmoSIS or claim a reproduction of the paper’s sampler/posterior.

## Full-covariance conditional redshift-block score

For a BAO-frozen cosmological shape $\theta_{\rm BAO}$, let $g(\theta_{\rm BAO})$ be the distance-modulus shape with any fixed H0 gauge, and model the SN vector as

\[
y=g+M\mathbf1+\epsilon,\qquad \epsilon\sim\mathcal N(0,C),\qquad C=P^{-1}.
\]

Predeclare disjoint redshift index sets $T$ (training/calibration SNe) and $H$ (the held-out redshift block), preserving the release row map. Define $K_T=C_{TT}^{-1}$, $q_T=\mathbf1_T^TK_T\mathbf1_T>0$, and $r_T=y_T-g_T$,

\[
\widehat M_T=\frac{\mathbf1_T^TK_Tr_T}{q_T},\quad
A=C_{HT}K_T,\quad
a_H=\mathbf1_H-A\mathbf1_T,
\]

\[
S_H=C_{HH}-C_{HT}K_TC_{TH},\qquad
V_{H|T}=S_H+\frac{a_Ha_H^T}{q_T},
\]

\[
\bar y_{H|T}=g_H+A r_T+a_H\widehat M_T.
\]

The $a_Ha_H^T/q_T$ term is the rank-one predictive variance from the flat-prior posterior for the single intercept, $M|y_T,\theta\sim\mathcal N(\widehat M_T,1/q_T)$. This is the proper flat-prior-integrated predictive distribution for a nonempty training set. Its normalized held-out score is

\[
-2\log p(y_H|y_T,\theta)=e_H^TV_{H|T}^{-1}e_H+\log|V_{H|T}|+|H|\log(2\pi),\qquad e_H=y_H-\bar y_{H|T}.
\]

This is well-defined for the verified positive-definite $C$ and any nonempty $T$. A plug-in/profile score evaluates at $\widehat M_T$ but uses only the conditional noise covariance $S_H$; it drops $a_Ha_H^T/q_T$, so it understates predictive uncertainty. Fitting $M$ using held-out rows is leakage for a held-out prediction. A leave-one-block-out conditional score is valid per block but those scores are not additive. A sequential factorization may score blocks after a designated intercept-calibration anchor block; with the improper flat prior the first block alone has no normalized predictive density unless a proper prior for $M$ is supplied.

Equivalent computation directly from the released precision $P=C^{-1}$, after permuting rows into $(T,H)$ consistently, is:

\[
A=C_{HT}C_{TT}^{-1}=-P_{HH}^{-1}P_{HT},\quad
S_H=P_{HH}^{-1},\quad
K_T=C_{TT}^{-1}=P_{TT}-P_{TH}P_{HH}^{-1}P_{HT}.
\]

Use the Schur complement for the **marginal training precision** $K_T$; the block $P_{TT}$ alone is not $C_{TT}^{-1}$. This formulation permits Cholesky solves using the packed precision and avoids constructing a noisy explicit full inverse. The covariance-space vector governing intercept uncertainty remains $a_H=\mathbf1_H-C_{HT}C_{TT}^{-1}\mathbf1_T$.

If block scores are intended as an additive decomposition, fix a sequential redshift partition and an anchor/training rule in advance, then score only subsequent blocks conditionally. If the aim is simply to evaluate the fixed BAO shape on all SNe, the unpartitioned full-covariance, one-offset-marginalized likelihood is the direct contract. The conditional statistic is a held-out predictive diagnostic, not a second independent SN likelihood or a model-refitting procedure.

## BAO/SN independence and reuse caveats

The DES-Dovekie paper uses DESI DR2 BAO in its combined-probe analysis ([v3, §4.2](https://arxiv.org/html/2511.07517v3)); this supports evaluating that combination but is not a published DESI–SN cross-covariance. Here the BAO and SNe are different data products and no shared object rows are indicated. Their redshift ranges do overlap: local DESI BAO anchors span $z=0.295$ to $2.33$, while Dovekie spans $z_{HD}=0.02509$ to $1.14418$. Overlap in cosmic volume/sky can still induce small correlations. Neither public product supplies a cross-probe covariance, so multiplication of the two likelihoods assumes zero BAO–SN cross-covariance and should be labeled as that approximation.

Dovekie is not independent of original DES-SN5YR: the paper reports 1,718 common SNe (of 1,820 Dovekie and 1,829 DES-SN5YR rows). It is an updated/recalibrated analysis of substantially the same observations. Do not multiply Dovekie by the original DES-SN5YR likelihood, nor by an existing Dovekie chain/joint DESI+SN summary that already includes these factors.

## Recommendation

**Proceed** with the frozen-shape prediction using the exact pinned row order and full STAT+SYS precision, one common additive magnitude offset, no separate `MUERR` term, and the conditional formula above if scoring a held-out redshift block. Fit/freeze the cosmological shape using BAO only; calibrate/marginalize $M$ on the predeclared training set, and retain the rank-one predictive variance. An adapter is needed for the actual SNANA whitespace input and release-specific paths. Record the BAO–SN zero-cross-covariance assumption. No comparative shape score or refit was calculated in this audit.

Numerical checks were run with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`; the existing read-only parser/matrix check completed in under one second. The covariance verification used one thread and no model evaluations.

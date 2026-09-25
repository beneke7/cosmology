# 2026 short-range-gravity literature update

**Status: exact review PDF acquired; targeted full-text check; literature-only.**
This does not add a new experimental measurement or a model-specific bound.

## Exact source

J. Murata, T. Fujiie, and S. Suzuki, *Short-Range Tests of the Gravitational
Inverse-Square Law*, arXiv:2605.18212v2, revised 8 September 2026. The exact v2
PDF and extracted searchable text are local:

- PDF: [`short_range_tests_2026_v2.pdf`](short_range_tests_2026_v2.pdf), SHA-256
  `4b69e60159db4cc6ef3c0e84490f0f3b42bcaa5d743c6e1e8bdd878a56d47506`.
- Text: [`short_range_tests_2026_v2.txt`](short_range_tests_2026_v2.txt), SHA-256
  `d7830cd7dbc7169ef7ef86e12985442908222a3b0aaf632f8ec8d9cc413de37c`.
- Retrieval provenance: [`manifest.lock.json`](manifest.lock.json).

## Relevance to thermodynamic gravity

The review's Eq. (4) uses a Yukawa parameterization, but §3 explicitly limits
its plotted discussion to **attractive** interactions with `alpha > 0`; it says
repulsive `alpha < 0` constraints are omitted. Consequently, its displayed
curves cannot be applied directly to the single-length higher-gradient TG
branch (`alpha = -1`) or its correlated two-range/complex-root potential.
The reported HUST `Lambda < 11 micron` statement is for the ADD `n=2`
power-law transition scale, not a TG Yukawa-range limit.

The review identifies the University of Washington (Lee et al. 2020; minimum
surface gap 52 microns) and HUST (Tan et al. 2020; minimum gap 210 microns) as
the leading tabletop torsion-balance results in the micrometre-to-millimetre
range. It also summarizes the 2024 Berkeley lattice atom-interferometry result
but says its sensitivity does not yet surpass the existing bounds in the
review's Fig. 3. This is consistent with the independent primary-source audit
in [`work/tg_lab_bound_independent/REPORT.md`](../tg_lab_bound_independent/REPORT.md):
no newer primary gravitational-strength Yukawa contour was located for this
band in the bounded search.

For TG, retain the primary-data result rather than transferring the review's
extra-dimensional parameter conversion: the 2020 UW paper reports Newtonian
gravity as an excellent fit and a 95% single-Yukawa `|alpha|=1` limit of
`lambda < 38.6 microns`. The sign-resolved UW thesis table brackets the
`alpha=-1` one-range crossing between 25 and 28 microns; its interpolated
27.1-micron value is not a new calibrated limit. The linked two-range model
still requires evaluating both potential templates through the apparatus
geometry and jointly profiling the measured torque and nuisance parameters.

## Interpretation

There is no claimed physical detection here. The review is useful for updating
the literature map and preventing an attractive-only contour or ADD power-law
limit from being mistaken for a test of TG's sign-correlated Yukawa sector.
The complete two-range torque re-fit remains a valid future test, but the
current campaign's higher-priority local discriminator is the bounded
axisymmetric NGC 3198 calculation using resolved baryon sources.

Sources: [review v2](https://arxiv.org/abs/2605.18212v2), [UW primary
measurement](https://arxiv.org/abs/2002.11761), [Tan et al. primary
measurement](https://doi.org/10.1103/PhysRevLett.124.051301), [Panda et al.
atom interferometry](https://doi.org/10.1038/s41586-024-07561-3).

# TG galaxy prediction feasibility audit

Status: **sample split reproducible; a defensible 3D TG prediction is blocked by the currently staged SPARC inputs and boundary implementation.** No TG fit or spherical inversion was run.

## What can be predeclared

The official SPARC Table 1 lists 175 galaxies and defines rotation-curve quality as Q=1 (high), 2 (medium), or 3 (low). Its published dynamical sample excludes Q=3 and face-on systems with inclination below 30 degrees, giving 153 galaxies. The downloaded official machine-readable table reproduces the catalog counts (99/64/12) and this 153-object sample (93 Q=1, 60 Q=2). All 153 have a rotation-model file in the staged `Rotmod_LTG.zip`.

For a concrete future test, hold out the 31 galaxies whose SHA-256 of `tg-sparc-v1|<SPARC name>` sorts first; the other 122 are training galaxies. The split uses only the quality/inclination eligibility table and galaxy identifiers. It does not use any velocity values or residuals. The exact held-out IDs and a regeneration command are in `feasibility.json`. The simplest preregistered theory target would be one shared positive K estimated using only the 122 training galaxies; held-out velocities would be scored only after predictions are frozen.

## What the baryon data support

SPARC supplies useful stellar *surface* data. Its 3.6-micron radial surface-brightness profiles, distances, inclinations, and stellar disk/bulge decompositions can define an axisymmetric stellar source after fixing an M/L prescription. The small official `BulgeDiskDec_LTG.zip` acquired here contains radial `SBdisk` and `SBbulge` profiles in Lsun/pc², not volumetric density maps. The raw `sfb_LTG.zip` profiles are on-sky magnitudes without inclination or extinction corrections. The official page warns that the decompositions are uncorrected and that bulge light is linearly extrapolated at large radii.

One explicit 3D stellar prescription is available in the SPARC master paper: an exponential vertical disk with `z_d = 0.196 R_d^0.633`, and spherical bulges. With fixed stellar M/L values and stated radial extrapolations, this turns the 1D light profiles into a reproducible *assumed* stellar density. It does not make the geometry unique: scale heights, bars/lenses, non-axisymmetry, bulge flattening, warps, and outer extrapolations are not determined by the radial profiles.

The missing piece is HI source density. Table 1 gives total HI mass and an HI radius, while each Rotmod row gives `Vgas`, `Vdisk`, `Vbul`, observed velocity/error, and stellar surface-brightness columns. `Vgas` is a Newtonian force contribution computed from an HI distribution; it is not `Σ_HI(R)` or a 3D source. The SPARC paper says gas models use HI surface-density profiles where available (with a thin-disk calculation and a 1.33 helium factor); for D512-2, D564-8, and D631-7 it substitutes an exponential profile of scale length `2 R_d` because HI profiles were unavailable. The basic SPARC downloads staged here do not include the HI surface-density profiles or maps. A force curve sampled over a finite radial interval does not uniquely supply the missing mass profile without extra geometry, edge, and regularization assumptions.

For the 153-galaxy test, obtain a versioned per-galaxy HI surface-density profile (or map) with distance, inclination/warp handling, uncertainties, and an outer taper. A 3D source then still needs an explicit gas thickness/flaring rule; molecular gas is also absent from the basic files. The official basic-data page does not expose one complete HI-profile table for this sample, so this is a source-by-source recovery task from the referenced observations, not a single small SPARC download. Do not infer it by a spherical inversion of `Vgas`.

## Boundary condition

The supplied TG galaxy paper uses a spherical radial solver. It sets the inner and outer derivative boundaries at the first and last rotation-curve radii from the *observed* velocities (its Eqs. 41–42), then fits K and stellar M/L; it says the last radius is the furthest point in the velocity curve. That procedure cannot be used on held-out galaxies because it feeds the target velocities into the prediction.

The stationary 3D TG equation does admit a physical isolated-galaxy choice: regularity at the center and `phi -> 0` at infinity. Under `u = exp(-K phi)`, this is `u -> 1`, and the vacuum exterior obeys `nabla²u = 0`. A finite computational domain can therefore use an exterior harmonic/multipole match, then demonstrate convergence as the domain and baryon taper move outward. This is a boundary prescription that can be declared without target velocities; it has not been implemented or validated in the current TG work. Whether isolation is adequate for each SPARC environment, and how omitted matter beyond the HI profile is treated, must also be declared. For an axisymmetric source, the prediction to compare with a circular speed is `v_c²(R) = R * partial_R phi(R,0)`.

Thus a predeclared whole-galaxy split and a no-leakage physical boundary can be specified now, but the proposed out-of-sample TG calculation cannot yet be evaluated from the staged SPARC data: the resolved gas source and a 3D exterior-matching solver are absent. The 2024 conditional NGC 3198 exercise is not an out-of-sample validation under this contract.

## Exact sample split

Eligibility follows the SPARC paper's 153-object science sample: Q in {1,2}, inclination >= 30 degrees. There is no cut on Vflat or on the shape of a rotation curve. Select the 31 holdouts by ascending SHA-256 of `tg-sparc-v1|<Galaxy>`; train on the remaining 122.

```sh
awk 'NF==19 && $18 ~ /^[12]$/ && $6+0>=30 {print $1}' work/data6/SPARC_Lelli2016c.mrt |
while IFS= read -r galaxy; do
  digest=$(printf 'tg-sparc-v1|%s' "$galaxy" | sha256sum)
  printf '%s %s\n' "${digest%% *}" "$galaxy"
done | sort -k1,1 | awk 'NR<=31 {print $2}'
```

The machine-readable record preserves the IDs, criteria, inputs, hashes, evidence status, and outstanding source/boundary work. The official table, stellar-decomposition ZIP, and bulge-luminosity table fetched for this audit are small; no HI cubes, survey products, or images were acquired.

## Sources checked

- Official [SPARC data page](https://astroweb.case.edu/SPARC/) and its [machine-readable Galaxy Sample table](https://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt).
- Lelli, McGaugh & Schombert, [SPARC: Mass Models for 175 Disk Galaxies](https://astroweb.case.edu/ssm/papers/AJv152n157.pdf), especially §§3.1–3.3 and Table 1.
- Pszota & Ván, local extracted text `context/papers/tg_preprint.txt`, especially Eqs. (41)–(42), (46)–(49), and the NGC 3198 discussion. The source says the NGC 3198 data are based on SPARC.
- Existing equation audit `runs/tg_audit.json` checks synthetic spherical cases only; it is not a SPARC prediction.

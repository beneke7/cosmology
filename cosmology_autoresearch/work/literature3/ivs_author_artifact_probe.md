# Yang et al. IVS public-artifact probe

Checked 2026-09-24 UTC (2026-09-25 in the project timezone). This is a bounded follow-up to [`ivs_current_source.md`](../literature2/ivs_current_source.md) and [`solver_path.md`](../perturbation_audit/solver_path.md). It checks arXiv's live record and v1 source bundle, the APS accepted/version-of-record routes, and the code repository explicitly cited by the paper. No code was cloned, installed, or run; no fit or author contact was made. The arXiv PDF and source archive were streamed only for byte counts/hashes and text/source inspection; neither was newly saved.

## Verdict

The exact, versioned arXiv v1 paper PDF is publicly available and byte-matches the cached project PDF. The remaining implementation artifacts are still not publicly verifiable in the primary sources checked: no modified-CAMB repository or patch/base commit, no complete paper-run Cobaya configuration, and no explicit scalar perturbation initial-condition choice. APS still identifies the paper as accepted; its accepted-paper page says supplements appear upon version-of-record publication, and the current issue listing does not include this title. A faithful source-paper replication therefore remains blocked. A separately pinned independent implementation could test the stated equations, but could not establish equivalence to the authors' chains or numerical results.

## Artifact findings

| Requested artifact | Finding on 2026-09-24 UTC |
|---|---|
| Exact/versioned paper PDF | **Found:** official [arXiv v1 PDF](https://arxiv.org/pdf/2609.05660v1), 4,524,530 bytes, SHA-256 `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`. The direct streamed hash matches the existing `context/papers/interacting_de_desi_dr2_2026.pdf`; that local manifest records arXiv v1 and the same digest. The arXiv history lists only v1, submitted 2026-09-04 18:45:40 UTC; both prospective v2 record and PDF routes returned HTTP 404. |
| arXiv source bundle | **Found, but not implementation artifacts:** [v1 source bundle](https://arxiv.org/src/2609.05660v1), 3,662,580 bytes, streamed SHA-256 `bfa9fad7bd0c8c6b159f65c99b508d6ab4d67e3661e727c66ac743f78ca1354f`. Its archive listing contains `00README.json`, `clean.tex`, and plot PDFs only. No CAMB source, patch, Cobaya YAML/INI, chain, or likelihood/config lockfile is present. `clean.tex` is the accepted-paper text carrying the arXiv v1 date and APS DOI. |
| Modified CAMB repository, patch, and base commit | **Not found in checked public primary links.** The paper says it used a “modified version” of CAMB but gives no repository, diff, release, or commit. Its only explicit software repository link is `wgcosmo`, cited for the Cobaya/MCEvidence interface (see below), not as provenance for the CAMB modification. |
| Full Cobaya configuration | **Not found.** The arXiv source bundle has no run configuration. The linked `wgcosmo` default-branch tree has general Planck and DESI DR2 likelihood examples and standard model YAMLs, but no IVS/interacting-DM-DE configuration. These examples do not identify the paper's full parameter wiring, likelihood/data versions, sampler controls, or solver options. |
| Explicit perturbation initial conditions | **Not found in v1 text/source.** The paper states synchronous gauge and its dark-sector perturbation closure/equations, but does not specify adiabatic/isocurvature mode, regular-mode initialization, or another scalar IC prescription. Searches of `clean.tex` for “initial condition(s),” “adiabatic,” “isocurvature,” and “primordial” found no matches. Without the missing run config, CAMB's actual IC selection cannot be inferred. |
| APS version of record / supplement | **Not identified as of the check.** APS's [accepted-paper page](https://journals.aps.org/prd/accepted/10.1103/6kdf-vzq4) says accepted 2026-09-03 and states author supplements, if provided, become available upon publication of the version of record. The APS [Volume 114, Issue 6 listing](https://journals.aps.org/prd/issues/114/6), current to 2026-09-24, is marked partial and does not list this title or DOI. The direct APS abstract route returned an internal browser error; PDF and supplement routes returned a Cloudflare 403. Thus the exact arXiv v1 PDF is available; an APS VOR PDF or supplement was not verifiable. |

## Repository check

The article's TeX cites [William Giarè's `wgcosmo`](https://github.com/williamgiare/wgcosmo) specifically for the Cobaya interface used to calculate MCEvidence. I checked GitHub's official repository metadata and recursive default-branch tree without fetching repository source files. At the check, `main` pointed to commit [`80fa9742f10616bcaf55d5888a15e8f07e74cd25`](https://github.com/williamgiare/wgcosmo/commit/80fa9742f10616bcaf55d5888a15e8f07e74cd25), dated 2026-08-05; `main` was the only listed branch and no tags/releases were listed. The current tree contains generic CMB/Planck and BAO/DESI DR2 examples and model YAMLs, but no matching IVS config, interaction-specific CAMB patch, or CAMB fork. It does include an introductory CAMB notebook, which is not evidence for the paper's modified solver. The main-branch timestamp also predates the paper's 2026-09-04 submission. This finding is scoped to the linked public repository/tree; it does not prove that no unlinked, private, or later-posted author artifact exists.

## Checked primary URLs

All checks below were made 2026-09-24 UTC.

- [arXiv record and submission history](https://arxiv.org/abs/2609.05660): one v1 entry; accepted-for-publication comment and APS DOI `10.1103/6kdf-vzq4`.
- [arXiv v1 PDF](https://arxiv.org/pdf/2609.05660v1): HTTP 200; direct body hash and byte count recorded above.
- [arXiv v1 source bundle](https://arxiv.org/src/2609.05660v1): HTTP 200; archive streamed and inspected; hash and byte count recorded above.
- [Prospective arXiv v2 record](https://arxiv.org/abs/2609.05660v2) and [v2 PDF](https://arxiv.org/pdf/2609.05660v2): both HTTP 404.
- [arXiv v1 HTML article](https://arxiv.org/html/2609.05660): checked article methods, perturbation discussion, code citations, and IC terms.
- [APS accepted-paper page](https://journals.aps.org/prd/accepted/10.1103/6kdf-vzq4): accepted status/date; supplement availability is deferred to VOR publication.
- [APS Volume 114, Issue 6](https://journals.aps.org/prd/issues/114/6): searched current issue listing for title and DOI; neither appears.
- APS direct [abstract](https://journals.aps.org/prd/abstract/10.1103/6kdf-vzq4), [PDF](https://journals.aps.org/prd/pdf/10.1103/6kdf-vzq4), and [supplement](https://journals.aps.org/prd/supplemental/10.1103/6kdf-vzq4) routes: the abstract browser route returned an internal error; PDF and supplement returned HTTP 403 Cloudflare challenges. These failures alone do not establish nonexistence.
- [Paper-cited GitHub repository](https://github.com/williamgiare/wgcosmo), its [official repository API metadata](https://api.github.com/repos/williamgiare/wgcosmo), [recursive main tree](https://api.github.com/repos/williamgiare/wgcosmo/git/trees/main?recursive=1), and [main commit](https://github.com/williamgiare/wgcosmo/commit/80fa9742f10616bcaf55d5888a15e8f07e74cd25): metadata/tree-only inspection; details above.

## Remaining blockers

To change the conclusion from “independent implementation only” to a faithful source replication, the minimum missing provenance is the authors' exact modified-CAMB patch/repository and its base commit/version, the complete Cobaya input for the claimed run (including likelihood/data releases and solver/sampler settings), and the scalar initial-condition choice. A matching chain/output archive and exact APS version of record/supplement would further tie the implementation to the published results. No available public artifact in the checked sources closes the first three gaps.

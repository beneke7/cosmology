# Primary-source provenance check: Yang et al. DESI DR2 interacting vacuum

Checked 2026-09-24 (UTC date). Scope was restricted to arXiv's record/version history and article text, the APS accepted/final record, and links explicitly present in the article. No files were downloaded; no code was installed, cloned, or run; no parameter fit was performed.

## Findings

- The exact article is Weiqiang Yang et al., “Do DESI-DR2 BAO data imply a coupling of dark matter and dark energy?”, [arXiv:2609.05660](https://arxiv.org/abs/2609.05660), DOI [10.1103/6kdf-vzq4](https://doi.org/10.1103/6kdf-vzq4). The arXiv record lists submission on 2026-09-04, “Accepted for publication in Physical Review D,” and only one history entry: v1, 2026-09-04 18:45:40 UTC. The versioned HTML labels itself `2609.05660v1`. A direct check of the prospective v2 record returned an error, consistent with no v2 listed in the history on the audit date. No correction or replacement is recorded there.
- APS's official [accepted-paper page](https://journals.aps.org/prd/accepted/10.1103/6kdf-vzq4) says accepted 2026-09-03. The DOI-linked [APS abstract route](https://journals.aps.org/prd/abstract/10.1103/6kdf-vzq4) was inaccessible to the browser tool during this check; the accepted page is live. I found no version-of-record/final APS article page or article-linked supplemental artifacts to inspect as of this date. Thus “final” solver/config claims cannot be checked against an APS version of record here.
- The arXiv v1 text describes the model as `Q=Gamma rho_x`; for IVS, `w_x=-1`. It says “modified version” of CAMB with Cobaya for MCMC (Section III), and reports `Gamma/H0` as the dimensionless coupling. It does not give a CAMB tag/commit, modified-source URL or patch, full Cobaya YAML/INI, dependency/likelihood build lock, chain/config archive, or downloadable spectra/products. It names the Planck `plikTTTEEE+lowl+lowE` likelihood combination, which is not enough to identify the authors' complete runtime environment.
- The article's CAMB citation [153] is Lewis & Bridle's 2002 CAMB paper, not a release or source snapshot. The only explicit code link in the article is its reference [157], [William Giarè's `wgcosmo` repository](https://github.com/williamgiare/wgcosmo), used for the Cobaya interface to the MCEvidence calculation (Section IV). The repository's browsable README describes a general cosmology tools collection and says `yamls/` contains collected Cobaya configurations. It does not identify an IVS configuration or the CAMB modification. The GitHub `yamls/` directory listing could not be fetched in this browser session; so I make no claim that every file in that directory was exhaustively ruled out. Crucially, the paper points to `wgcosmo` for the MCEvidence interface, not explicitly as the source of its modified CAMB.
- The paper contains no link to CLASS and makes no claim that CLASS was used. The verifiable solver is CAMB only, and the precise implementation/configuration remains unverified from the allowed public primary sources. This independently confirms the existing `solver_path.md` caution; it narrows its phrasing: `wgcosmo` is explicitly provenance for the evidence interface, while the paper merely cites the 2002 CAMB paper and provides no identifiable modified CAMB artifact.

## What was checked / reproducibility limits

Opened the arXiv abstract record and history, the arXiv v1 HTML, the APS accepted-paper page, and the article-linked `wgcosmo` repository README. Searched the article HTML for `CAMB`, `CLASS`, `github`, `modified`, `code`, and `Cobaya`; inspected the methods paragraph, relevant repository citations, and references [153] and [157]. Checked the APS DOI abstract URL and prospective arXiv v2 URL; neither returned a usable page. Search was deliberately not expanded to secondary commentary or unrelated implementations.

These sources establish the named solver, sampler, likelihood label, model equations, and the general repository's stated scope. They do not establish which CAMB fork/patch generated the paper's chains, whether the same patch served all three DE equations-of-state, exact solver options/initial conditions, full likelihood releases, or a config-to-result chain. Absence of a link in the checked article/page is not proof no private or unlinked author artifact exists. The arXiv history can change after 2026-09-24.

If the user supplies the authors' modified-CAMB source/patch plus its base commit and the matching complete run config (including likelihood identifiers/releases), the next step changes from asking authors / treating a new solver as independent to a provenance review: verify patch applies to the declared base, compare the implemented background and IVS perturbation closure with the published equations, check unit/sign conventions and config wiring, then identify remaining reproduction gaps. User-provided artifacts would need to be inspected, not presumed equivalent or run automatically.

## Official URLs consulted

- https://arxiv.org/abs/2609.05660
- https://arxiv.org/html/2609.05660
- https://arxiv.org/abs/2609.05660v2 (unavailable / no v2 entry observed)
- https://doi.org/10.1103/6kdf-vzq4
- https://journals.aps.org/prd/accepted/10.1103/6kdf-vzq4
- https://journals.aps.org/prd/abstract/10.1103/6kdf-vzq4 (unavailable in this check)
- https://github.com/williamgiare/wgcosmo
- https://github.com/williamgiare/wgcosmo/tree/main/yamls (directory listing fetch failed)

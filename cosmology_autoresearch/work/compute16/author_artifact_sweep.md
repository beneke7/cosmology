# Final author-link and supplement sweep — arXiv:2609.05660

Read-only follow-up, accessed 2026-09-25 00:40–00:42 UTC. This supplements [the cited-repository audit](../compute14/primary_source_artifact_audit.md). No code was run; no repository was cloned, no paper source archive or dataset was downloaded, and no authors were contacted.

## Author-profile boundary

The official [arXiv abstract](https://arxiv.org/abs/2609.05660) and [HTML v1](https://arxiv.org/html/2609.05660) list six authors. Clicking each author-name link on arXiv leads to a generic arXiv author-search page—not an ORCID, institutional profile, or GitHub profile:

| Author | arXiv author-search URL |
| --- | --- |
| Weiqiang Yang | <https://arxiv.org/search/astro-ph?query=Yang%2C+W&searchtype=author> |
| Sibo Zhang | <https://arxiv.org/search/astro-ph?query=Zhang%2C+S&searchtype=author> |
| Dinorah Barbosa | <https://arxiv.org/search/astro-ph?query=Barbosa%2C+D&searchtype=author> |
| Marco Antonio Alvarez | <https://arxiv.org/search/astro-ph?query=Alvarez%2C+M+A&searchtype=author> |
| Supriya Pan | <https://arxiv.org/search/astro-ph?query=Pan%2C+S&searchtype=author> |
| Leila Graef | <https://arxiv.org/search/astro-ph?query=Graef%2C+L&searchtype=author> |

The article HTML has no `orcid.org` author links. Its only paper-content GitHub citation is reference [157], [`williamgiare/wgcosmo`](https://github.com/williamgiare/wgcosmo), cited for the Cobaya interface to MCEvidence. The cited repository is maintained by William Giarè, who is not among the paper's six listed authors; it was checked at pinned HEAD `80fa9742f10616bcaf55d5888a15e8f07e74cd25`. A fresh read-only `git ls-remote` at 00:41 UTC still returned that same SHA for `HEAD` and `refs/heads/main`, with no tag refs.

The APS [accepted-paper page](https://journals.aps.org/prd/accepted/10.1103/6kdf-vzq4) exposes the author names and affiliations but no author-profile or repository links. No arbitrary GitHub-user-name searches or unlinked author personal/institutional pages were followed; those are outside this bounded public-link sweep.

## Paper-linked supplements and artifact result

- The arXiv v1 page exposes PDF, HTML, and TeX-source navigation only. I did not download the TeX archive; its HTML rendering was searched for author profile/code links. It shows the cited `wgcosmo` link, but no CAMB patch, IVS Cobaya YAML, Dovekie likelihood/version, or explicit perturbation initial conditions.
- The APS page remains an “Accepted Paper” page (accepted 2026-09-03). It states that author supplemental materials, if provided, will become available with the version of record; it presents no supplement/code/data attachment. The two plausible APS supplemental URL forms were not retrievable by the browsing tool, so I treat the page's lack of a linked attachment—not those failed direct opens—as the evidence.
- The official arXiv page still displayed only v1 at access time. No post-v1 source or supplement was linked from that page.

**Conclusion:** this author-link/supplement sweep yielded no additional public artifact identifying a modified-CAMB patch/version, complete paper-specific Cobaya config, DES-Dovekie likelihood version, or IVS perturbation initial conditions. It is bounded to links surfaced by the primary article/APS pages and their listed author-name routes; it does not establish that unlinked or private code does not exist.

## Exact checked identities / URLs

- arXiv abstract/version history: <https://arxiv.org/abs/2609.05660>
- arXiv HTML v1: <https://arxiv.org/html/2609.05660>
- arXiv v1 permalink: <https://arxiv.org/abs/2609.05660v1>
- APS accepted paper / DOI: <https://journals.aps.org/prd/accepted/10.1103/6kdf-vzq4>
- Sole paper-cited repository: <https://github.com/williamgiare/wgcosmo>
- Current observed repository HEAD: `80fa9742f10616bcaf55d5888a15e8f07e74cd25` (`main`); direct refs check: <https://github.com/williamgiare/wgcosmo/commit/80fa9742f10616bcaf55d5888a15e8f07e74cd25>
- APS supplement paths attempted (not retrievable in browser): <https://journals.aps.org/prd/supplemental/10.1103/6kdf-vzq4> and <https://journals.aps.org/prd/supplemental/6kdf-vzq4>.

# DES-Dovekie — updated DES supernova calibration

- Source/version: arXiv:2511.07517 v3; local full PDF/text read at `context/papers/des_dovekie.pdf` (SHA-256 `cd3dd615c36afe0d45bcdfabd9c31611f0ca511e325997f93c3eb1f5a9744eb6`).
- Question: how does recalibration of DES-SN5YR and historical samples affect Type Ia supernova distances and cosmology?
- Likelihood: Eqs. (8)–(9), printed pp. 4–5, give the correlated distance-modulus statistic and systematic covariance construction. Tables 10–11 report calibrated cosmology/model comparisons.
- Data: the pinned release is repository commit `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`; its Hubble diagram has 1,820 rows and the STAT+SYS archive stores a packed inverse covariance. The release README states that row ordering matters.
- Overlap/limits: 1,718 objects overlap the original DES-SN5YR in the paper's comparison. The recalibration is an updated branch, not an independent likelihood to multiply by original DES. Uncalibrated SN distances leave a magnitude/H0 degeneracy.
- Reproduction caveat: the pinned module requests CSV although the paired `.csv` file contains whitespace SNANA records. Campaign implementations parsed that native format and followed release semantics; they did not run the full CosmoSIS pipeline or reproduce its published posterior.
- Experiment implication: preserve row order and full STAT+SYS precision, profile one magnitude offset, test CPL bounds, and never infer H0 from the SN-only screen.
- Local calculation: `work/inference/REPORT.md`, `work/critic/REPORT.md`, `experiments/dovekie_screen/`.

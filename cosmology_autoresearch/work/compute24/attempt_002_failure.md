# Attempt 002 — downstream spectrum-array alignment

Status: **fit and distance stages reached the CMB-accuracy check, but their values were not checkpointed before the subsequent error, so this run contributes no saved profile result.** `get_cmb_power_spectra` returned different array lengths for `AccuracyBoost=1` and `2` (2551 and 2527), and the check tried to combine unequal shapes. This is a comparison-indexing defect, not a physical result.

- Command was identical to attempt 001, with the objective now catching `CAMBError` / `CAMBFortranError` and mapping points without an H0 bracket to the invalid-domain penalty.
- No JSON, plot, or profile checkpoint was written before the CMB-stage exception; do not quote fit coordinates or scores from this attempt.
- The follow-up uses the multipole intersection for the accuracy comparison and writes a durable fit/background checkpoint before attempting spectra.
- Host observation during this run: 16 unrelated `interIsoFoam` processes each consumed about one CPU. They were left untouched; the next CAMB run is capped at four worker processes to use the available headroom without competing for those CPUs.

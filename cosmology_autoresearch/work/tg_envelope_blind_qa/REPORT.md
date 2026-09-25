# Independent blind method QA: NGC 3198 nuisance envelope

**Scope:** static review of `nuisance_envelope.py`, the frozen resolved-source audit, the SPARC input catalog, and the prior TG pilot report. No nuisance-envelope progress, results, grid, prediction, or score files were opened. No full run or numerical scan was performed. Findings below evaluate code and stated assumptions only.

## Severity-ranked findings

### High — source geometry transform and mass Jacobian are internally correct, subject to the declared thin-disk/axisymmetry assumptions

At `nuisance_envelope.py:120-123`, sky offsets are rotated to disk major/minor axes, the minor coordinate is divided by `cos(i)`, and physical radii scale with distance. At lines 131-151 and 185-204, each LOS column is multiplied by `cos(i)` before annular averaging; the deprojected annular area used for coverage is the face-on area times `cos(i)`, i.e. the projected annulus. Consequently the integrated mass cancels the `cos(i)` factors and scales as `D²`, as the module docstring states (lines 4-10). The map audit independently records the adopted THINGS/HALOGAS geometries, LOS-to-face-on conversion, and helium factor (`work/tg_resolved_source_audit/REPORT.md:15-21`). The sky-WCS helper assumes diagonal CD/CDELT terms (`:96-110`); that assumption should remain tied to the audited map headers.

This is a sound 2D thin-disk deprojection followed by axisymmetrization, not a unique reconstruction of a 3D source. That limitation is correctly stated in the frozen source audit (`work/tg_resolved_source_audit/REPORT.md:23-31`) and in the previous TG report (`work/tg_axisymmetric_pilot/REPORT.md:46-54`).

### Medium — the allowed inclination grid exceeds the HI map's reported geometry robustness

The rotation-curve inclination prior is correctly centered at 73° with σ=3°: the local SPARC Table 1 record lists NGC 3198 as `D=13.80±1.40 Mpc`, `Inc=73.0±3.0°` (`work/data6/SPARC_Lelli2016c.mrt:10-16,152`). The code uses that center for the five-node ±2σ grid (`nuisance_envelope.py:43-45,60-62,79-80,373-383`) and applies the expected LOS deprojection factor to both velocity and error (`:291-293,463-465`). Holding the LOS measurement fixed makes this rescaling internally consistent.

However, the HI source audit says HALOGAS models become marginally inconsistent for inclination shifts of about 2° and recommends geometry near 71° with at least a ±2° check (`work/tg_resolved_source_audit/REPORT.md:15-21`). The envelope re-deprojects its maps across 67°–79° around the SPARC 73° center (`nuisance_envelope.py:373-383`), considerably beyond that stated map-specific robustness. The transformed profiles are a mathematical sensitivity envelope there, but endpoints should not be described as equally validated physical HI sources without checking the map/model convention. Treat this as a source-validity limitation on those grid nodes, not an error in the `sin(i)` correction.

### Medium — nuisance scoring is a grid MAP screen, not posterior profiling/integration

The data objective is diagonal `Σ[(Vobs−Vmodel)/errV]²`; distance, inclination, and log-M/L Gaussian penalties are added once as squared standardized offsets (`nuisance_envelope.py:286-293,482-488,490-503,512-529,538-560,679-727`). The same transformed HALOGAS+stellar source and sampled nuisance point feed TG edge, NFW, and MOND. NFW alone profiles `log M200` and `log c200` at each nuisance node (`:296-320,490-503`); MOND is fixed in `a0` and TG is an edge-limit curve (`:512-529,538-553`). These choices are coherent for a declared descriptive, penalized grid comparison. The output documentation correctly calls it diagonal/descriptive and records no covariance and unequal model complexity (`:679-727`; prior report `work/tg_axisymmetric_pilot/REPORT.md:82-99`).

The sampled minimum is not a continuous nuisance optimum or a marginalized comparison: D, i, and M/L are discrete nodes, and each vertical case is selected separately. The source and code declare this, but any later interpretation must preserve it. The NFW bound interval has no concentration–mass relation or halo prior (`:298-320,699`), so its fit flexibility is intentionally greater than MOND/TG; the prior report also flags this (`work/tg_axisymmetric_pilot/REPORT.md:95-99`).

### Medium — NFW optimizer convergence is recorded but not enforced

Three bounded least-squares starts are run and the returned fit with the lowest residual is retained (`nuisance_envelope.py:296-320`). `optimizer_success` and its message are returned, but lines 491-503 do not reject or separately qualify a failed best fit. A failed/budget-exhausted solve could therefore enter a best-model comparison as though profiled. This is a code-level QA item to check in the already-produced records before using NFW scores; I did not inspect those records.

### Low — finite-u examples and the TG edge are properly separated

The principal-mode curve uses `Kcrit` and the logarithmic derivative of the positive eigenmode (`nuisance_envelope.py:330-358`), and is labeled as a limiting curve with no finite normalization or K prior (`:538-553`). Separate finite-u solves use fixed `K/Kcrit=0.95,0.99`, check positivity and finite positive speeds, and make no K-profile/prior claim (`:566-622`). This matches the prior report's distinction between the finite solution and singular eigenmode limit (`work/tg_axisymmetric_pilot/REPORT.md:63-73,109-117`).

### Low — source tail and clipping conventions are declared but remain physical model choices

The source is clipped to nonnegative H I per pixel (`nuisance_envelope.py:168-177`) and again after azimuthal averaging (`:147-150`); this avoids negative physical mass but can bias low-surface-brightness annuli upward if negative map noise is present. The code extrapolates the stellar profile exponentially beyond the P5 support and tapers gas linearly over 5 kpc (`:217-250`), then samples the radial source on a finite 0–50 kpc grid (`:64-67,372`). These are transparent assumptions, not data-determined outer 3D structure. The source audit explicitly calls for outer stellar extrapolation and notes unresolved opacity, missing baryons, and vertical uncertainty (`work/tg_resolved_source_audit/REPORT.md:19-31`); the previous report also identifies the stellar-tail source as a sensitivity (`work/tg_axisymmetric_pilot/REPORT.md:65-79`).

## Checks performed

- Read the nuisance-envelope source, focusing on geometry/profile conversion, NFW and MOND definitions, prior construction, edge and finite-u calculations, and result metadata.
- Read the frozen resolved-source audit and previous NGC 3198 pilot report for declared map conventions, geometry validity, source limits, and TG edge interpretation.
- Checked the SPARC Table 1 machine-readable source: it supports the code's 13.8 Mpc, 1.4 Mpc and 73°, 3° reference values.
- Analytically checked the cosine area/surface-density cancellation; the SPARC `sin(i)` LOS deprojection; `ρcrit=3H0²/(8πG)` and the standard `M200,c200` NFW expression with H0 converted from 73 km/s/Mpc to km/s/kpc; and the exponential vertical profile normalization on a reflected half-domain.
- No execution, fit, numerical scan, or score/output inspection was performed as part of these checks.

## Audit integrity note

One broad text search accidentally printed a line from a frozen, earlier pilot candidate JSONL while locating inclination provenance. It was not a current nuisance-envelope output and no score value was used in this audit. The current nuisance-envelope `envelope_progress.jsonl`, `results.json`, `envelope_grid.csv`, and live scores were not opened. This report is limited to method/code findings and was written before any output inspection.

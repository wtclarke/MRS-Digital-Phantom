# LCModel conformance check: saved legacy BASIS-water voxel

## Result

Three actual LCModel 6.3-1N fits completed with exit status 0, valid numerical
TABLE/COORD/PRINT/PS output and generated PDFs. No fatal errors were reported.
The original acquisition, BASIS and executable hashes are unchanged. No simulator,
calibration, signal amplitudes or noise realizations were changed/regenerated.

The result supports approximate recovery for prominent metabolites/groups, **not
agreement for every individual metabolite or scanner-level absolute calibration**.
The saved table values themselves are explicitly assumed nominal mM.

Input directory:
`outputs/legacy_basis_water_notebook/20261101T182122316149Z/voxel_0000_100_108_80/`.
Both the actual `total.nii.gz` and the saved `metabs.nii.gz` were fitted against the
matching supplied GE PRESS TE35 BASIS, using the separate `water_reference.nii.gz`.
All 18 metabolites with saved ground truth were included; a five-metabolite fit
would incorrectly omit actual simulated contributions. H2O and unused controls
were excluded from metabolite fitting.

## Reused protocol and necessary acquisition changes

The earlier Phase 7b protocol was retained: DOWS=T, DOECC=F, ATTH2O=ATTMET=1,
NAA scaling reference WSPPM=2.01/N1HMET=3, PPMCEN=4.68, 0.2–4.2 ppm fitting,
FWHMBA=1/127.731 ppm, NSIMUL=0, NRATIO=0, and original baseline/line-shape defaults.
No phantom-mode, modified-linewidth or post-hoc concentration corrections were used.
Necessary changes were matching the GE basis, TE35 (rather than TE30), 18 present
metabolites, input/output paths/title and the separately identified WCONC cases.
BASIS/source SHA256 is
`f0a42cf0bef1b5427d4f3f3f4e5dd35dd4406581230a7ead1055bd7fbd700a79`.

| Case | Acquisition | WCONC (mM) | Phased-real relative residual |
| --- | --- | ---: | ---: |
| total_previous_water | total | 39016.50000000 | 1.80624% |
| total_voxel_water | total | 45636.98676816 | 1.80623% |
| metabs_voxel_water | metabs | 45636.98676816 | 0.72117% |

`total_previous_water` retains the earlier pure-WM water setting. `total_voxel_water`
uses this mixed voxel's saved effective water concentration. `metabs_voxel_water`
uses exactly the existing metabolite-only file as a diagnostic, with the same water;
it still contains the simulated tissue-dependent metabolite linewidth and spatial
shim effects. No clean signal was regenerated.

The voxel fractions are approximately WM=0.194, GM=0.806, CSF=0. Nominal voxel and
WM/GM-only concentrations therefore differ only by legacy float32 fraction-sum
roundoff. No CSF or tissue-only concentration correction was applied.

## Concentration comparison (mM)

| Metabolite | Nominal voxel | Total, old water | Total, voxel water | Total error | Metabolites only, voxel water |
| --- | ---: | ---: | ---: | ---: | ---: |
| Asc | 1.073478 | 0.357 | 0.418 | -61.06% | 0.881 |
| Asp | 2.021706 | 1.584 | 1.853 | -8.34% | 1.693 |
| Cr | 5.471112 | 4.297 | 5.026 | -8.14% | 5.053 |
| GABA | 1.730572 | 1.348 | 1.576 | -8.93% | 1.338 |
| Gln | 2.409640 | 1.877 | 2.195 | -8.91% | 2.036 |
| Glu | 8.440927 | 7.616 | 8.909 | +5.55% | 8.175 |
| Gly | 4.722277 | 3.514 | 4.111 | -12.94% | 4.476 |
| GPC | 0.727706 | 0.312 | 0.365 | -49.84% | 0.641 |
| GSH | 1.597329 | 1.259 | 1.473 | -7.78% | 1.479 |
| Lac | 1.444271 | 0.838 | 0.980 | -32.15% | 1.342 |
| mI | 4.450162 | 3.377 | 3.950 | -11.24% | 4.068 |
| NAA | 8.176772 | 6.794 | 7.947 | -2.81% | 7.894 |
| NAAG | 0.960166 | 0.614 | 0.718 | -25.22% | 0.852 |
| PCh | 0.527690 | 0.724 | 0.847 | +60.51% | 0.543 |
| PCr | 2.666186 | 2.556 | 2.990 | +12.15% | 2.733 |
| PE | 1.826958 | 1.798 | 2.103 | +15.11% | 1.761 |
| sI | 0.304593 | 0.201 | 0.236 | -22.52% | 0.264 |
| Tau | 1.975306 | 1.274 | 1.491 | -24.52% | 1.683 |

Grouped values below sum the actual individual TABLE concentrations. No grouped
uncertainty is inferred by assuming independence.

| Group | Nominal mM | Total, voxel water | Error | Metabolites only | Error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Cr+PCr | 8.137297 | 8.016 | -1.49% | 7.786 | -4.32% |
| GPC+PCh | 1.255396 | 1.212 | -3.46% | 1.184 | -5.69% |
| NAA+NAAG | 9.136938 | 8.665 | -5.17% | 8.746 | -4.28% |

Individual GPC/PCh total-fit errors are -49.84%/+60.51%, despite their combined
error being only -3.46%. The metabolite-only fit reduces those individual errors
to -11.91%/+2.90%. This is consistent with overlap/model/background sensitivity;
the comparison does not uniquely identify its cause. Low %SD does not establish
unbiased recovery. In particular, a printed 0%SD in the clean diagnostic is rounded,
not proof of zero uncertainty. There is only one realization here.

## Water normalization and residual interpretation

Changing WCONC multiplies the quantification scale by **1.1696842815**. With the same
unchanged total and water samples, NAA moves from 6.794 to 7.947 mM. The earlier
39.0165 M WM assumption is inappropriate for this mixed voxel. It explains a
substantial common scale difference; it does not explain every individual error.

From the printed LCModel normalized NAA reference area 0.0014601 and water area
15931700, independently calculate:

```
FCALIB = normalized_basis_area * (2 * WCONC * ATTH2O) / water_area
```

The reference basis area already includes N1HMET=3 and ATTMET=1 normalization.
This reconstructs FCALIB=7.1515270e-6 (printed 7.15155e-6) for the old water setting,
and 8.3650288e-6 (printed 8.36505e-6) for voxel water. Differences are compatible
with rounding of the printed areas. These are quantification factors, not raw
metabolite coefficients or new simulator gain corrections.

The total fit has FWHM=0.046 ppm, S/N=97, data shift=-0.008 ppm, phase=-6 degrees
and first-order phase=-2.4 degrees/ppm in its TABLE. Reported diagnostics are
STARTV24, RFALSI11 and FINOUT9 informational messages, retained in the outputs.
No different optimizer/baseline configuration was tried to improve concentrations.
A low residual is not evidence of accurate individual coefficients.

The metabolite-only residual is 0.72117%, versus 1.80623% for the total fit. The
latter contains empirical MM/lipids, residual water and receiver noise not explicitly
represented by the earlier spin-system-only protocol. Remaining clean-fit errors
can include shared-response/integration conventions, individual tissue linewidths,
legacy attenuation and non-Voigt shim dephasing. Their contributions were not
uniquely decomposed here. No amplitude or concentration correction was introduced.
COORD arrays use LCModel's report scaling, so y-axis levels can differ between the
WCONC cases despite their input complex arrays being exactly identical.

## Conversion and input integrity

RAW/H2O conversion uses the validated conjugation convention and `(2E25.16)` text.
Independent parsing reproduces every input FID sample **exactly**, maximum complex
error **0**, without resampling, first-point replacement or amplitude normalization.
The same original spectra are retained in the LCModel-compatible BASIS; only
format/reference headers change, following the earlier documented ISHIFT correction.

The saved metabolite NIfTI-1 dwell is float32 0.0005000000237487257 s; the water
NIfTI-2 dwell is float64 0.0005 s. These agree with the declared 2000-Hz grid within
4.75e-8 relative float32 rounding. CONTROL uses the simulator's declared 0.0005 s.
This does not change any samples or source headers. All input hashes are recorded
and rechecked after fitting, including all source files, BASIS and executable.

Initial setup attempts discovered two conversion-helper limitations: interleaved
NMUSED blocks were parsed as spectrum numbers, and existing FMTBAS could be
duplicated. A failing small regression reproduced both; the helper now respects
all namelist boundaries and writes a single correct format. The original BASIS is
untouched. Failed setup artifacts are preserved separately; they ran no fits.

## Reproducibility and reports

All Python execution used Pixi. Exact run command (licensed local binary supplied
through LCMODEL_BINARY; no binary was installed or modified):

```bash
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m research.validation.lcmodel_legacy_water.run \
  --source outputs/legacy_basis_water_notebook/20261101T182122316149Z/voxel_0000_100_108_80 \
  --output outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits \
  --binary "$LCMODEL_BINARY"
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m research.validation.lcmodel_legacy_water.summarize \
  outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits
# Historical adapter/diagnostic suite: 25 passed, 1 Dash warning, 2.13 s.
# Quantitative-only diagnostic/noise tests have since been removed.
# Retained LCModel adapter checks: tests/test_lcmodel_quantitative_io.py.
```

The actual binary argument used the existing local `test_lcm/lcmodel` installation.
Each fit took approximately 2.2–2.5 seconds, exited 0, with stdout/stderr captured.
Complete actual command lines and per-case controls are saved in `commands.txt`
and case directories. No LCModel settings were tuned after inspecting results.

Outputs under `outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits/`:

- `comparison.csv`: all 18 metabolites, all three cases, nominal voxel and tissue-only values, errors and %SD.
- `grouped_comparison.csv`, `summary.json`, `analysis_summary.json`, `input_hashes.json`.
- `concentration_comparison.png`, `fit_comparison.png`.
- Each case: converted RAW/H2O, CONTROL, execution/results JSON, and
  `lcmodel_outputs/fit.{TABLE,PRINT,COORD,PS,pdf}`.

All three PDFs are nonempty and were generated from actual LCModel PostScript
using the already-installed Ghostscript. Existing earlier fits were not overwritten.

Production simulator code/calibration is unchanged. Modified source: the research
conversion helper and its regression test; added research runner, summary script
and this report. Nothing was staged or committed. The immediate conclusion is
approximate recovery of strong peaks/groups with voxel-matched water scaling,
not comprehensive individual-metabolite concentration conformance.

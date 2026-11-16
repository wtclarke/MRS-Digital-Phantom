# Phase 5b.2: optional water alignment and CLI reports

## Implementation and interpretation

`simulation/water_reference_preprocessing.py::align_water_reference` returns a
separate complex FID and JSON-serializable correction metadata. It requires
explicit known frequency and phase (the generator's values in this experiment):

\[
w_{aligned}[n]=w[n]\exp[-i(2\pi f_w n\Delta t+\phi_w)].
\]

The time origin is zero, frequency is in Hz, phase is in radians, and correction
operates on interpreted FID samples before storage conjugation. The original
array is not modified. Each multiplier has unit magnitude, so the initial
amplitude and every sample magnitude are preserved to floating-point precision.
There is no amplitude recalibration, inferred empirical factor, or fitted
concentration used in preprocessing. Invalid/nonfinite signals or corrections
are rejected. This function does not estimate parameters from noisy data.

Original acquisition/calibration metadata remain provenance. The aligned
dataset adds `water_preprocessing` to simulation_config.json and writes
preprocessing_metadata.json with input frequency/phase, applied signed
corrections, zero output frequency/phase, dwell, sample count, time origin,
unit magnitude scaling, original physical water position and analysis frame.

Both exports retain SpecFreqChemShift=4.68, EchoTime=.030 s, dwell=.0005 s,
127.731 MHz and 1H metadata. **Physical water remains 4.70 ppm**. Original
samples have +2.55462 Hz in the 4.68 ppm basis frame; aligned samples are
digitally represented at DC, which plots at 4.68 ppm. That is explicitly an
analysis transformation, not a claim that water's physical chemical shift has
changed. The metabolite FID and frequency metadata are identical between the
two pairs. FSL-MRS uses the metabolite MRS axes when handling its water
reference, so zero-offset alignment in that frame is intentional. Re-labeling
only the water NIfTI center as 4.70 would not be the same analysis convention.
All custom information is in JSON; no existing exporter was changed.

## CLI version and report option

Installed FSL-MRS is `2.5.1.dev78+ga311e38f2`. `pixi run fsl_mrs --help` and
the installed parser expose **`--report` (singular)**. They do not expose
`--reports`. The requested plural spelling is unavailable in this version;
the documented singular flag was explicitly used to generate actual HTML,
rather than omitting reports or altering the installed package. CLI help is
saved in the output directory.

The original/aligned driver is reproducible using a fresh directory:

```
pixi run python -m examples.fsl_mrs_water_preprocessing outputs/fsl_mrs_water_preprocessing
```

That command was executed and ran these two real CLI commands:

```
pixi run fsl_mrs --data outputs/fsl_mrs_water_preprocessing/original/metabolites.nii.gz --h2o outputs/fsl_mrs_water_preprocessing/original/water_reference.nii.gz --basis data/basissets/LCModel_Universal_Philips_UnEdited_PRESS_TE30.BASIS --output outputs/fsl_mrs_water_preprocessing/original/fsl_mrs_results --keep NAA --algo Newton --lorentzian --baseline off --ppmlim 0.2 4.2 --no_conj_basis --TE 30 --TR 2 --tissue_frac 1 0 0 --internal_ref NAA --wref_metabolite NAA --ref_protons 3 --ref_int_limits 1.8 2.2 --report
pixi run fsl_mrs --data outputs/fsl_mrs_water_preprocessing/aligned/metabolites.nii.gz --h2o outputs/fsl_mrs_water_preprocessing/aligned/water_reference.nii.gz --basis data/basissets/LCModel_Universal_Philips_UnEdited_PRESS_TE30.BASIS --output outputs/fsl_mrs_water_preprocessing/aligned/fsl_mrs_results --keep NAA --algo Newton --lorentzian --baseline off --ppmlim 0.2 4.2 --no_conj_basis --TE 30 --TR 2 --tissue_frac 1 0 0 --internal_ref NAA --wref_metabolite NAA --ref_protons 3 --ref_int_limits 1.8 2.2 --report
```

All settings match Phase 5b: 10 mM NAA, 40 M water, shared G=1, PRESS TE30,
2000 Hz, 2048 points, zero linewidth, no noise or additional relaxation. TR=2 s
and pure WM fractions are the explicit software quantification assumptions
used previously, not attenuation added to the simulated input. See Phase 5b
for the standard T1/T2, tissue-water densities and integration/proton model.

## Numerical results

| Quantity | Original | Aligned |
|---|---:|---:|
| Expected initial water amplitude | 61.73993435001698 | 61.73993435001698 |
| Fitted original-scale water amplitude | 61.73993439719796 | 61.739934352118276 |
| Fitted frequency Hz | 2.5546199998466563 | −2.314906889218109e-9 |
| Fitted phase rad | 3.2565282020350634e-10 | −2.6285071601920303e-8 |
| Fitted Lorentzian linewidth Hz | 0 | 0 |
| Actual complex water-fit relative residual | 8.326188710382894e-10 | 3.40014252030215e-8 |
| Integrated-model versus input complex relative difference | 1.4422720554676098 | 3.403460097795722e-11 |
| Unscaled source NAA coefficient (CLI) | 152.1269729860657 | 152.1269729860657 |
| Internal normalized NAA coefficient | 1.0000000905278943 | 1.0000000905278943 |
| Standard FSL-MRS molarity mM (CLI) | 6.003318300954079 | 6.003318305337433 |
| Synthetic water-referenced recovery mM | 10.000000897637063 | 10.00000090493861 |
| Metabolite relative complex residual (CLI fitted FID) | 1.8546135081783173e-5 | 1.8546135081783173e-5 |
| Water optimizer success | true | true |
| Metabolite optimizer success | true | true |
| CLI exit code | 0 | 0 |

CLI concentrations are read directly from each `concentrations.csv`, and CLI
predicted metabolite FIDs are reloaded from `fit.nii.gz`. The CLI does not
serialize water optimizer parameters or its uncorrected full-complex water
model. For those diagnostics the driver **separately replays the exact installed
WaterRef handler with its default initialization**, and records that provenance
explicitly. Likewise supplementary optimizer success/messages come from a
matching API replay. No optimization initialization or bounds are patched for
these two comparisons. The API predicted metabolite FID and CLI fitted file
are array-identical (max difference 0), and the CLI copies of the water inputs
match the supplied references exactly. The actual HTML traces are inspected
separately below, so an API fit alone is not treated as report verification.

Maximum magnitude change due to preprocessing: **1.4210854715202004e-14 SRU**.
Metabolite FID change: **0**. Both water export/reload errors: **0** at default
complex128 precision. Maximum complex difference from the independent ideal
DC reference is **2.2218607601071443e-13 SRU**; remodulating gives a maximum
complex reconstruction difference **2.2469334198890886e-13 SRU**. The independent analytical preprocessing tests compare
the complete aligned complex arrays (including a decaying envelope) at
rtol/atol 1e-13, while testing magnitude preservation more tightly. Ordinary
floating-point phase arithmetic, not amplitude rescaling, limits exact equality.

## HTML observations and comparison figures

Actual generated reports:

- [Original report](../outputs/fsl_mrs_water_preprocessing/original/fsl_mrs_results/report.html)
- [Aligned report](../outputs/fsl_mrs_water_preprocessing/aligned/fsl_mrs_results/report.html)

Both CLI runs exited successfully; both reports are nonempty and contain
Quantification information and Water referencing sections. Their original
sizes were 887447 and 887405 bytes respectively. Plotly traces were decoded
from the actual report HTML (including binary-encoded arrays), and the
overlapping real-spectrum traces compared numerically:

- Original: water-data peak bin 4.702936385842121 ppm; integrated-model peak
  4.68 ppm. Relative real-spectrum overlay difference 3.084510192820103.
- Aligned: both peak at 4.68 ppm. Relative real-spectrum overlay difference
  3.403451542701252e-11. No extra display transformation was required.

The original plot discrepancy remains the reference-frame mismatch explained
in Phase 5b.1, not an actual failed fit. The additional figure shows the true
fit before integration correction in the data's own frame, alongside the report
overlay and the unchanged suppressed metabolite fit:

![Original and aligned overlays, actual water fits, and metabolite fits](../outputs/fsl_mrs_water_preprocessing/water_fit_comparison.png)

HTML interactive traces use the installed report template's Plotly CDN;
the standalone PNG permits inspection without that browser dependency.

## Controlled tests and robustness

Tests exercise zero/zero, frequency-only (+2.55462 Hz), phase-only (0.6 rad),
and combined frequency/phase, as well as repeated deterministic preprocessing.
They independently calculate the expected complex demodulation and magnitude
envelope, assert input immutability, test correction metadata, export/reload
each aligned case and run the installed default WaterRef fitting. All aligned
cases converge with amplitude within 2e-6 relative, frequency and phase within
1e-6, and full-complex relative residual below 1e-6.

The Phase 5b.1 regression continues to preserve the unprocessed combined case's
known zero-amplitude failure as an expected failure. Preprocessing avoids that
initialization failure on these deterministic cases by removing rotation before
the default wrapped-phase regression. It does not fix the installed algorithm.
For the already accurate original case, preprocessing improves report
interpretability rather than numerical accuracy: its residual is slightly
larger while still small. Both amplitude and synthetic recovery are preserved
within floating-point/optimizer tolerances, and the standard FSL concentration
remains ~6.0033 mM because its quantification assumptions did not change.

Commands executed for testing:

```
pixi run python -m pytest -q tests/test_water_reference_preprocessing.py tests/integration/test_water_preprocessing_cli.py
pixi run python -m pytest -q tests
```

The focused preprocessing/CLI run passed **13 tests**, including two real CLI
runs with HTML reports in temporary directories. The full regression run
completed with **67 passed, 1 xfailed**, zero unexpected failures, in 46.86 s.
The expected failure is the documented installed WaterRef combined-phase/
frequency defect on the unprocessed signal. Five existing dependency,
deprecation and optional-display warnings were emitted. CLI processes completed
with exit code zero; stdout/stderr are preserved in each dataset's `cli.log`.

## Files and limitations

The dedicated output directory contains original and aligned input pairs,
unchanged calibration and ground-truth provenance, aligned preprocessing JSON,
actual `fsl_mrs_results/report.html` files and CLI-generated artifacts, diagnostic
arrays, `comparison.csv`, `validation_results.json`, `cli_commands.txt`, CLI
help/logs, the comparison PNG and `html_overlay_inspection.json`. Report files
were not renamed or replaced with links.

This is optional preprocessing using known synthetic parameters. It does not
estimate phase in noise, resolve multi-component water, remove eddy currents,
or establish physical scanner calibration. Alignment is performed only on the
water reference; it is not a physical shift of the metabolite signal. No legacy
model, calibration constant, basis, GUI, installed FSL package or existing
exporter was modified. Future LCModel use requires its own explicit frequency,
phase and water-scaling validation; that work has not started.

# Phase 5b: FSL-MRS fitting and water scaling

## Outcome and scope

FSL-MRS recovers the declared synthetic NAA coefficient and water amplitude.
At 10 mM, the synthetic water-referenced interpretation gives
10.000000897637065 mM. Standard FSL-MRS tissue/relaxation/integration scaling
gives 6.003318300954079 mM. The input signals and Phase 4b calibration were
not changed to resolve this difference.

This is a matching-basis numerical development test. It does not independently
validate the basis's molar response, physical scanner calibration, or LCModel.

## Installed implementation and settings

Pixi contains FSL-MRS `2.5.1.dev78+ga311e38f2`, Python 3.14.8, NumPy 2.5.4.
The installed package's reported version is recorded directly rather than
inferred from the manifest constraint/build name. Both `fsl_mrs` CLI and
Python `fitting.fit_FSLModel` are available. Inspection used the actual installed
`core/mrs.py`, `core/basis.py`, `utils/fitting.py`, `utils/results.py`,
`utils/quantify.py`, `utils/misc.py`, `utils/constants.py` and
`models/model_lorentzian.py`. Absolute source paths are recorded in the output
validation configuration.

The Python API uses Newton/TNC fitting, a Lorentzian model, NAA only, ppm range
0.2–4.2, one metabolite group, and baseline disabled (`baseline_order=-1`).
Consequently the recorded baseline is exactly zero by construction; this
experiment does not test flexible baseline recovery. Frequency, phase and
additional Lorentzian linewidth remain fitted. Optimizer initialization comes
from FSL-MRS, without a ground-truth parameter vector supplied as x0.

Data are the unchanged Phase 5a PRESS TE30 pair at 2000 Hz / 2048 complex
samples: G=1, water 40 M, zero water linewidth, default +2.55462 Hz offset,
zero phase, no noise, residual water, additional metabolite relaxation or
T1/TR attenuation. Four input NAA values are 5, 10, 15, 20 mM. Each case is
exported and reloaded before fitting. Both datasets retain dwell time,
127.731 MHz frequency, 1H nucleus, and SpecFreqChemShift=4.68 ppm.

NIfTI-MRS EchoTime is 0.030 s. The synthetic dataset has no TR-dependent
attenuation or physical repetition time. **TR=2 s is an explicit quantifier
assumption**, not recovered from the synthetic input. TE/TR and a homogeneous
pure-WM fraction dictionary are provided directly to QuantificationInfo. CLI
users must supply equivalent parameters or appropriate acquisition metadata.
Simulator-specific metadata remain unmodified in JSON sidecars.

## Expected fitted coefficient

Let b be the formatted source NAA basis, k the Phase 4a scalar, C the
concentration in mol/L, and d=(C/k)b the simulated FID. The API's
`rescaleForFitting()` sets the data norm to 100 with multiplier s_D, and the
selected basis's column-mean norm to 100 with multiplier s_B. This uses
`misc.rescale_FID`, not the simulator's mean-absolute-sample normalization.
It also scales H2O by s_D. With no independent basis scaling, the optimizer
therefore sees a coefficient

\[
a_{expected}=(C/k)\,s_D/s_B.
\]

The unscaled source coefficient is a*s_B/s_D; `getConc(scaling='raw')`
uses that same scaling ratio. The synthetic concentration is
1000*k*a*s_B/s_D mM. Because s_D changes with input concentration, the
internal coefficient remains approximately 1 at all concentrations. It is
the unscaled source coefficient that must be linear in C.

At 10 mM, k=6.573456836083183e-05, s_D=117.8317751259257 and
s_B=17925.38964873409. Expected internal coefficient=0.9999999999999988;
fitted=1.0000000905278943 (relative error 9.05278956e-08).
Expected source coefficient=152.12695921433226 approximately;
recovered=152.1269729860657.

Recorded 10 mM fit nuisance parameters, in model units:

- gamma=5.312649275675941e-07 s^-1 (additional FWHM gamma/pi Hz).
- eps=−1.2169327966303667e-06 rad/s, with exp(−i eps t).
- Phi0=5.7547802108371384e-05 rad and Phi1=1.6900500716563533e-07
  rad/Hz, applied as exp(−i(Phi0+Phi1*frequency)).
- Full-FID relative residual=1.8546135081785328e-05; baseline=0.

TNC reported successful convergence by objective-change tolerance. The small
nonzero nuisance parameters explain why the noiseless matching-basis residual
is not identically zero. Acceptance checks require coefficient recovery within
2e-6 relative and full-FID residual below 5e-5; parameters/calibration were not
tuned to force exact residual zero.

## Water fitting and amplitude

`quantify.WaterRef` actually fits a single Voigt FID in a water spectral window:
A exp(−gamma*t−sigma*t²−i*omega*t) exp(i*phi). It temporarily rescales the
water spectrum to norm 100, minimizes the complex spectral residual, then
undoes that temporary scaling and returns a phase/frequency-corrected FID.
Its integral is the trapezoidal sum of abs(real(spectrum)) over its limits.
This is not a peak-height amplitude measurement.

The returned water FID is still scaled by s_D from MRS, so the original-scale
amplitude is abs(wref.fid[0])/s_D. The script captures the water optimizer's
parameters without replacing its algorithm. At 10 mM:

- Expected original-scale amplitude=61.73993435001698.
- Fitted original-scale amplitude=61.739934397197956.
- Absolute error=4.7180976e-08; relative error approximately 7.64e-10.
- Frequency=−omega/(2*pi)=2.5546199998466563 Hz.
- Phase=3.2565282020350634e-10 rad; gamma=sigma=0.
- Temporary water-fit multiplier=0.0003038300392470479, after s_D scaling.
- Water optimizer succeeds with the projected-gradient convergence criterion.

For synthetic water referencing, invert the source coefficient convention and
then multiply the recovered synthetic concentration by (2q G C_w)/A_w_fit.
Here q is the fixed Phase 4b reference response, G=1 and C_w=40 M. This uses
the measured water amplitude and predeclared constants, never a target fitting
concentration to set an amplitude factor.

## Standard FSL-MRS concentration interpretation

QuantificationInfo selects stored 3T T1/T2 values at this field strength. We
explicitly specify WM=1, GM=CSF=0, TE=.030 s, TR=2 s, flip angle default 90°,
NAA as reference with 3 protons and integration range 1.8–2.2 ppm. The package
uses water molarity 55010 mM and water densities WM=.65, GM=.78, CSF=.97.
Thus even without relaxation its WM water assumption is 35.7565 M rather than
the simulation's explicit 40 M. There is no physical tissue density or
relaxation correction in the synthetic signals.

At 10 mM the defaults are T1_WM=.97 s, T1_GM=1.5 s, T1_CSF=4.47 s,
T1_metabolite=1.29 s; T2_WM=.073 s, T2_GM=.088 s, T2_CSF=2.03 s,
T2_metabolite=.194 s. These produce relaxation-corrected water concentration
20691.10854642174 mM, metabolite correction 1.4815764577216368, CSF correction
1 and additional correction 1.

The code quantifies an NAA reference spectral integral I_M from the fitted
phase-corrected model, relative to the fitted water integral I_W. The reported
NAA molarity is independently reconstructed as

\[
C_{FSL}=(I_M/I_W)(2/3)\,
C_{water,relaxed}\,R_{met}^{-1}\,F_{CSF}\,F_{additional}.
\]

At 10 mM I_M=115.17287561964316 and I_W=392080.6488209709 in FSL's scaled
spectral-integration units, giving 6.003318300954079 mM. For other metabolites
FSL divides its concentration scaling by the fitted reference coefficient
before multiplying each target coefficient. The script records both scalings
and reconstructs this expression, matching the package to roundoff.

The difference combines water/tissue/relaxation assumptions and a metabolite
singlet-integral/proton convention that is distinct from Phase 4b's full frozen
basis coefficient and Ref0ppm initial magnitude convention. No single physical
conversion follows from this synthetic input. The standard scaling interface
does not directly encode the complete Phase 4b response model; post-fit
interpretation using its declared q and k is appropriate for this experiment.
No extra `add_corr` or water-amplitude adjustment was introduced to force a
10 mM standard FSL result.

## Concentration series and determinism

| Truth mM | Synthetic water-referenced mM | Absolute error mM | Relative error | FSL standard molarity mM |
|---:|---:|---:|---:|---:|
| 5 | 5.000000448819 | 4.48819e-07 | 8.97637e-08 | 3.001659150477 |
| 10 | 10.000000897637 | 8.97637e-07 | 8.97637e-08 | 6.003318300954 |
| 15 | 15.000001150285 | 1.15028e-06 | 7.66856e-08 | 9.004977331965 |
| 20 | 20.000001795274 | 1.79527e-06 | 8.97637e-08 | 12.006636601908 |

Unscaled source coefficients grow linearly with concentration; relative
deviation is below 1e-7. A repeated independent 10 mM generation/export/read/
fit/quantification run produced maximum numeric comparison difference 0.
No fitting/quantification warnings were captured and both optimizers succeeded.
The full-FID residual is approximately 1.85e-5 throughout the series.

![Recovery under synthetic and FSL conventions](../outputs/fsl_mrs_water_scaling_phase5b/concentration_recovery.png)

## Reproduce and inspect outputs

Executed commands:

```
pixi run fsl_mrs --help
pixi run python -m examples.fsl_mrs_water_scaling_validation /private/tmp/phase5b_validation_20261030
pixi run python -m examples.fsl_mrs_water_scaling_validation /private/tmp/phase5b_validation_final_20261030
pixi run python -m examples.fsl_mrs_water_scaling_validation outputs/fsl_mrs_water_scaling_phase5b
pixi run python -m pytest -q tests/integration/test_fsl_mrs_water_scaling.py tests/test_mrscloud_quantitative_pair.py tests/test_mrscloud_water_calibration.py tests/test_frozen_basis_normalization.py tests/test_analytical_voxel_reference.py tests/test_analytical_water_reference.py tests/test_nifti_mrs_roundtrip.py
pixi run python -m pytest -q tests/integration/test_fsl_mrs_water_scaling.py
```

Combined run: **51 passed** (4 new integration tests plus 47 regressions),
5 existing dependency warnings. The final integration-only rerun after adding
optimizer/source-path reporting passed **4 tests**. The CLI inspection emitted
Matplotlib temporary-cache and package string-escape warnings; these did not
prevent execution. Use a fresh output directory when rerunning the script.

`outputs/fsl_mrs_water_scaling_phase5b` contains five case directories
(four concentrations plus repeated 10 mM), each with both input NIfTI-MRS files,
untouched synthetic configuration/calibration and ground-truth sidecars,
`fitting_results.json`, and `fit_arrays.npz`. The root contains
`comparison.csv`, `validation_config.json` (including software/source paths),
and the recovery figure. Expensive optimizer tests live under
`tests/integration/`, separately runnable from fast signal tests.

## Limits and LCModel implications

Matching generation and fitting bases can cancel common modeling errors. This
test isolates serialization, numerical fitting, coefficient normalization,
water amplitude estimation and concentration-convention differences; it does
not validate physical concentration or general fitting accuracy. The pure-WM
quantifier setting is explicitly an interpretation experiment, not tissue
mixing or a modification of the simulator. There is no noise, flexible
baseline, model mismatch, additional relaxation, or spatial effect here.

LCModel needs its own validation of BASIS amplitude conventions, receiver/water
scaling, spectral integration and tissue/relaxation assumptions. The FSL
concentration ratio observed here must not be transferred as an LCModel
calibration factor. Production simulators/calibration constants were unchanged.

## CLI HTML fitting reports

FSL-MRS was also executed with `--report` for the exact saved 5, 10, 15 and
20 mM inputs. Reports are under each case's `cli_report/report.html`:

- [5 mM](../outputs/fsl_mrs_water_scaling_phase5b/NAA_05mM/cli_report/report.html)
- [10 mM](../outputs/fsl_mrs_water_scaling_phase5b/NAA_10mM/cli_report/report.html)
- [15 mM](../outputs/fsl_mrs_water_scaling_phase5b/NAA_15mM/cli_report/report.html)
- [20 mM](../outputs/fsl_mrs_water_scaling_phase5b/NAA_20mM/cli_report/report.html)

Open **Quantification information** for tissue fractions, water densities,
T1/T2 values and concentration scaling factors; **Water referencing** shows
the fitted water and reference-metabolite integrals. At 10 mM, raw-to-molarity
scaling is approximately 6.00331775749 per internal fit coefficient, and the
internal coefficient is 1.00000009053. Their product is 6.003318300954079 mM.
The report's displayed relaxation-corrected water concentration comes from
`relax_corr_water_molar`; this installed report template labels it mmol/kg,
although the source property is molarity in mM (mmol/L). The underlying value
and the calculations above use the molar definition.

The following command was run separately for each concentration tag
(`05`, `10`, `15`, `20`), substituting that tag for `10` in the paths:

```
pixi run fsl_mrs --data outputs/fsl_mrs_water_scaling_phase5b/NAA_10mM/metabolites.nii.gz --h2o outputs/fsl_mrs_water_scaling_phase5b/NAA_10mM/water_reference.nii.gz --basis data/basissets/LCModel_Universal_Philips_UnEdited_PRESS_TE30.BASIS --output outputs/fsl_mrs_water_scaling_phase5b/NAA_10mM/cli_report --keep NAA --algo Newton --lorentzian --baseline off --ppmlim 0.2 4.2 --no_conj_basis --TE 30 --TR 2 --tissue_frac 1 0 0 --internal_ref NAA --wref_metabolite NAA --ref_protons 3 --ref_int_limits 1.8 2.2 --report
```

TE is in milliseconds on this CLI; TR is in seconds. Explicitly disabling
basis conjugation reproduces the API's already verified basis orientation.
CLI `--baseline off` reproduces `baseline_order=-1`. Each CLI directory
contains `options.txt` with the actual arguments, `concentrations.csv`,
`quantification_info.csv`, fitted/data/reference NIfTI files and fit plots.

All four CLI processes exited successfully. A Pixi Python numerical check
confirmed exact equality of saved CLI versus API molarity values (maximum
difference 0 mM), presence of the two quantification report sections, and
array-exact preservation of CLI copies of both input datasets.
`cli_report_comparison.csv` records these comparisons. No simulator source,
calibration constants, datasets or previous fitting results were changed.

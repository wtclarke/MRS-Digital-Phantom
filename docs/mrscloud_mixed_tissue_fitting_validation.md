# Phase 6c: mixed-tissue FSL-MRS fitting validation

## Result and scope

Actual FSL-MRS fits recover the Phase 6b voxel-average coefficients for all
five selected metabolites, rather than their WM/GM-only concentrations. Water
amplitudes recover the effective voxel water concentrations. Standard
water-referenced concentrations follow the independently predicted FSL-MRS
convention, including CSF correction; they do not directly equal voxel-average
ground truth. This is internal synthetic validation, not physical calibration.

Only a standalone validation script and tests were added. The Phase 6a bundle,
Phase 6b source datasets, production simulation, GUI, MRSCloud source and
calibration constants were preserved. Source and bundle SHA256 hashes before
and after execution are recorded in `results.json`; all match. The four input
files in each baseline output directory are byte-identical copies of Phase 6b.

## Reproduce and inspect

Run from the repository root with an unused output directory:

```sh
pixi run python --version
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.mrscloud_mixed_tissue_fitting_validation outputs/mrscloud_mixed_tissue_fitting_validation/final
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_mixed_tissue_fitting_predictions.py tests/integration/test_mixed_tissue_fitting_validation.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q
```

The writable Matplotlib cache avoids a sandbox warning and repeated font-cache
creation; no dependencies or Pixi configuration were changed. Outputs refuse
to overwrite an existing directory. The first script attempt failed because
receiver gain was looked up in ground truth rather than simulation config;
the script was corrected without touching either source. An initial fast-test
run had two fixture failures (a read-only MRS property and an overly precise
Phase 5g numerical anchor); both were corrected before the successful run.
No fit expectation or signal was adjusted to match a fitted concentration.

Installed FSL-MRS: `2.5.1.dev78+ga311e38f2`; Python 3.14.8. The matching complete 29-column
Phase 6a BASIS is loaded and the five known metabolite columns are selected:
NAA, Cr, GPC, PCh, Glu. The known-component test deliberately does not fit
27 unconstrained metabolites to five-component data. Receiver gain is 1.
Sampling is 2048 complex points, 2000 Hz bandwidth, dwell 0.0005 s,
spectrometer frequency 127.731 MHz, acquisition reference 4.68 ppm.
Water is unchanged at +2.55462 Hz (4.70 ppm), phase approximately zero and
1 Hz readout linewidth. No noise, additional relaxation, B0 or residual water.
TE=30 ms and TR=2 s are explicit quantification inputs, not altered metadata.

Newton fitting uses a Lorentzian model, baseline off, ppm 0.2–4.2, no basis
conjugation, one shared metabolite linewidth/frequency group, and fitted
zero/first-order phase. NAA is the water-scaling reference with 3 protons,
1.8–2.2 ppm; water uses 2 protons and corrected-model integration 1.65–7.65 ppm.

## Ground truth

All metabolite entries below are nominal mM per total voxel volume. The tissue
profiles are copied unchanged from the Phase 6b example. Water is mol/L.

| Voxel | WM/GM/CSF | NAA | Cr | GPC | PCh | Glu | Water M |
|---|---|---:|---:|---:|---:|---:|---:|
| WM | 1/0/0 | 10 | 8 | 2 | 1 | 8 | 39.0165 |
| GM | 0/1/0 | 14 | 9 | 1.5 | 1.2 | 12 | 47.2305 |
| CSF | 0/0/1 | 0 | 0 | 0 | 0 | 0 | 55.5 |
| Mixed | .2/.3/.5 | 6.2 | 4.3 | .85 | .56 | 5.2 | 49.72245 |

Mixed WM/GM-only concentrations are respectively 12.4, 8.6, 1.7, 1.12,
10.4 mM. These equal voxel values divided by .5 **only because CSF metabolite
concentrations are zero**. For pure CSF tissue-only values are null.

## Coefficient interpretation and recovery

With unit response B_m, data scaling s_D and common fitting basis scaling s_B:

\[
s(t)=G\sum_m10^{-3}C_{m,voxel}B_m(t),\qquad
 a_{m,expected}=G10^{-3}C_{m,voxel}s_D/s_B.
\]

FSL `getConc()` without a scaling returns internal coefficient a_m.
Its `raw` convention is a_m*s_B/s_D, in nominal mol/L times receiver gain
for this synthetic source model. Thus 1000*a_m*s_B/(s_D*G) recovers nominal mM.
No proton multiplier is applied to the generated basis or water responses.

| Voxel | Metabolite | Internal fitted coefficient | Expected internal coefficient | Recovered voxel mM | Relative coefficient error |
|---|---|---:|---:|---:|---:|
| pure_WM | NAA | 0.685054712076 | 0.685054711088 | 10.0000000144 | 1.442e-09 |
| pure_WM | Cr | 0.548043776844 | 0.548043768871 | 8.0000001164 | 1.455e-08 |
| pure_WM | GPC | 0.137010910455 | 0.137010942218 | 1.99999953634 | -2.318e-07 |
| pure_WM | PCh | 0.0685055024821 | 0.0685054711088 | 1.00000045797 | 4.58e-07 |
| pure_WM | Glu | 0.548043724607 | 0.548043768871 | 7.99999935386 | -8.077e-08 |
| pure_GM | NAA | 0.795789576191 | 0.795789594895 | 13.999999671 | -2.35e-08 |
| pure_GM | Cr | 0.511579066468 | 0.51157902529 | 9.00000072444 | 8.049e-08 |
| pure_GM | GPC | 0.0852632610714 | 0.0852631708816 | 1.50000158667 | 1.058e-06 |
| pure_GM | PCh | 0.0682104487527 | 0.0682105367053 | 1.19999845269 | -1.289e-06 |
| pure_GM | Glu | 0.682105372957 | 0.682105367053 | 12.0000001039 | 8.656e-09 |
| mixed | NAA | 0.758975799555 | 0.758975806279 | 6.19999994507 | -8.86e-09 |
| mixed | Cr | 0.526386443785 | 0.52638644629 | 4.29999997954 | -4.759e-09 |
| mixed | GPC | 0.104053145559 | 0.104053134732 | 0.850000088444 | 1.041e-07 |
| mixed | PCh | 0.0685526333537 | 0.0685526534704 | 0.559999835668 | -2.934e-07 |
| mixed | Glu | 0.636560332158 | 0.636560353654 | 5.1999998244 | -3.377e-08 |

Common s_B=15685.073904966044. s_D is 1074513.3772364885 (WM),
891572.7577666875 (GM), 1920095.4215438305 (mixed). Optimizer messages and all
nuisance parameters, including zero baseline coefficients, are in `results.json`.

| Voxel | Complex metabolite residual | gamma (s^-1) | epsilon (rad/s, model sign) | Phi0 (rad) | Phi1 (rad per Hz) | GPC/PCh fitted-parameter correlation |
|---|---:|---:|---:|---:|---:|---:|
| pure_WM | 5.3972e-08 | 3.8552e-19 | -1.1918e-07 | -2.8589e-09 | 7.137e-11 | -0.819264 |
| pure_GM | 1.4586e-07 | 0 | 7.6394e-08 | 6.6962e-08 | 1.8802e-11 | -0.821479 |
| mixed | 1.0636e-07 | 1.5489e-09 | -1.1129e-07 | 2.0056e-07 | 6.2171e-10 | -0.820825 |

All nonzero-data optimizers report success. gamma is **additional** decay; the
basis retains its original line shape. This table does not call gamma/pi the
measured total metabolite FWHM. Full-FID normalized complex design condition
number is 4.319608, GPC/PCh complex correlation magnitude .896531; the stacked
real/imaginary design in the fitted window has condition number 3.139449.
GPC/PCh parameter correlations are about -.82. The small residual alone is
therefore not sufficient evidence: each coefficient is separately compared
against known truth. Noise-free exact-basis recovery does not establish
robust individual GPC/PCh recovery under noise or model mismatch. Correlations
are diagnostics; basis amplitudes were not normalized for fitting by these
calculations.

## Water fitting

The installed `quantify.WaterRef` routine is exercised on each original water
FID, before metabolite rescaling, with its default optimizer initialization.
The diagnostic captures the optimizer without changing installed source or
initialization. Amplitude is the fitted initial complex-FID magnitude in input
units, not spectral peak height. Integral values used for quantification are
also recorded separately after the common fitting-data scaling.

| Voxel | Expected amplitude | Fitted amplitude | Relative amplitude error | Frequency Hz | Phase rad | Lorentzian Hz | Complex residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| pure_WM | 0.00788790829157 | 0.00788790812062 | -2.1672e-08 | 2.5546201458 | -1.4309e-07 | 0.9999999377 | 1.4428e-07 |
| pure_GM | 0.00954852056348 | 0.009548520364 | -2.0891e-08 | 2.5546201454 | -1.4239e-07 | 0.9999999382 | 1.4388e-07 |
| pure_CSF | 0.0112203531886 | 0.0112203529449 | -2.1719e-08 | 2.5546201454 | -1.4299e-07 | 0.9999999380 | 1.4393e-07 |
| mixed | 0.0100523144216 | 0.0100523141879 | -2.3254e-08 | 2.5546201464 | -1.4443e-07 | 0.9999999368 | 1.4486e-07 |

All water optimizers report a local minimum. Gaussian linewidth parameter is
zero. Water amplitude is proportional to the known effective water molarity;
the largest relative amplitude error is 2.33e-8. The previously documented
report overlay frame issue remains: WaterRef removes frequency and phase in
its integrated model while the original data remain unchanged. Its displayed
model/data difference (~1.316 relative here) is **not** the complex fit residual
above, which uses the model restored to the data frame. No display mismatch
was repaired by changing the exported water.

## Independent quantification prediction

Define J as the existing independent first-point-half-weighted orthonormal FFT,
absolute real spectrum and trapezoidal **bin** integration functional, with
nearest inclusive ppm endpoints. It reproduces installed `FIDIntegrator`.
Apply J to the unfitted formatted NAA unit response and to the independently
known corrected unit water line abs(W(0))*exp(-pi*t). No fitted amplitudes,
frequency shifts, final concentration functions or empirical ratios enter the
prediction. For this bundle:

- J_NAA = 0.006425124587755786, window 1.8–2.2 ppm.
- J_W = 0.004767430133932897, window 1.65–7.65 ppm.
- alpha = (2/3)*J_NAA/J_W = 0.898474916569691.

Alpha differs from the Phase 5g value 0.8984749391468103 by 2.26e-8 absolute:
the Phase 6a cross-calibration already measured a 1.58e-7 relative NAA array
difference from six-digit LCModel spectral serialization. Water cross-scale
error was 1.79e-16. The fast test allows the previously established 3e-7 basis
serialization tolerance against Phase 5g; predictions in this phase use the
actual Phase 6a arrays. No subset-dependent correction is introduced.

For target metabolite m, with the same NAA reference for every metabolite:

\[
C^{FSL}_{m,mM}=\frac{10^{-3}C_{m,voxel}J_{NAA}}{C_{w,voxel}J_W}
\frac{2}{3}\frac{55010\sum_j f_j d_jR_{w,j}}{(1-f_{CSF})R_M}.
\]

R(T1,T2)=(1-exp(-TR/T1))*exp(-TE/T2), excitation 90 degrees. The target-to-NAA
coefficient ratio uses **known source concentrations** in this prediction.
Integrating each target metabolite independently would incorrectly describe
FSL's NAA-reference quantification. Shared receiver gain and fitting-data and
basis scalings cancel. The same source NAA integral is used for all targets;
this does not apply per-metabolite amplitude calibration.

**Defaults:** d_WM=.65, d_GM=.78, d_CSF=.97, pure-water concentration 55010 mM;
installed 3T T1 values WM=.97, GM=1.50, CSF=4.47, METAB=1.29 seconds;
T2 values WM=.073, GM=.088, CSF=2.030, METAB=.194 seconds. R_M=.6749567292246237.
These corrections are intentionally retained in the actual baseline CLI run,
although synthetic signals have no corresponding attenuation.

**Matched API:** public `set_densities` uses each original tissue water
concentration times 1000/55010; public T1/T2 constructor overrides use T1=1e-6
seconds and T2=infinity to represent exactly zero attenuation in floating-point.
TE/TR stay .03/2 seconds. Fitted parameters are asserted unchanged between
quantification configurations. The FSL CSF correction is retained in both.

| Voxel | Default relaxed water sum mM | Matched water sum mM | CSF factor |
|---|---:|---:|---:|
| pure_WM | 20691.1085464 | 39016.5 | 1 |
| pure_GM | 22469.6743485 | 47230.5 | 1 |
| mixed | 20362.1439764 | 49722.45 | 2 |

Predictions below precede any fitted-result interpretation. Reported defaults
are both CLI and API results (they agree); matched results are Python API only,
not modified HTML reports. All concentrations in the following table are mM.

| Voxel | Metabolite | Voxel truth | Tissue-only truth | Default prediction | Default reported | Matched prediction | Matched reported |
|---|---|---:|---:|---:|---:|---:|---:|
| pure_WM | NAA | 10 | 10 | 7.059362383 | 7.059362483 | 8.984749166 | 8.984749293 |
| pure_WM | Cr | 8 | 8 | 5.647489907 | 5.647490061 | 7.187799333 | 7.187799529 |
| pure_WM | GPC | 2 | 2 | 1.411872477 | 1.411872167 | 1.796949833 | 1.796949439 |
| pure_WM | PCh | 1 | 1 | 0.705936238 | 0.705936571 | 0.898474917 | 0.898475339 |
| pure_WM | Glu | 8 | 8 | 5.647489907 | 5.647489522 | 7.187799333 | 7.187798843 |
| pure_GM | NAA | 14 | 14 | 8.866093281 | 8.866093169 | 12.578648832 | 12.578648673 |
| pure_GM | Cr | 9 | 9 | 5.699631395 | 5.699631916 | 8.086274249 | 8.086274988 |
| pure_GM | GPC | 1.5 | 1.5 | 0.949938566 | 0.949939581 | 1.347712375 | 1.347713815 |
| pure_GM | PCh | 1.2 | 1.2 | 0.759950853 | 0.759949881 | 1.078169900 | 1.078168521 |
| pure_GM | Glu | 12 | 12 | 7.599508527 | 7.599508675 | 10.781698999 | 10.781699209 |
| mixed | NAA | 6.2 | 12.4 | 6.759626423 | 6.759626453 | 11.141088965 | 11.141089016 |
| mixed | Cr | 4.3 | 8.6 | 4.688128003 | 4.688128043 | 7.726884282 | 7.726884349 |
| mixed | GPC | 0.85 | 1.7 | 0.926722977 | 0.926723086 | 1.527407358 | 1.527407538 |
| mixed | PCh | 0.56 | 1.12 | 0.610546903 | 0.610546732 | 1.006291907 | 1.006291625 |
| mixed | Glu | 5.2 | 10.4 | 5.669364096 | 5.669363981 | 9.344139132 | 9.344138942 |

The full `comparison.csv` additionally records differences from both ground-truth
conventions, prediction errors and the voxel-volume interpretation
C_reported*(1-f_CSF). Largest absolute prediction/report discrepancy is
1.441e-6 mM, largest relative discrepancy about 1.3e-6, consistent with measured
individual coefficient errors. The mixed matched NAA result 11.141089016 mM
has voxel-volume interpretation 5.570544508 mM, compared with voxel truth
6.2 mM and tissue-only truth 12.4 mM. Its remaining response bias is alpha,
not an unexplained fitted concentration scale or failed water fit.

### Pure CSF

The original all-zero complex metabolite array remains zero. An actual attempt
to rescale it in FSL-MRS produces infinite data scaling and non-finite samples
with divide-by-zero/invalid warnings. The invalid normalized spectrum is not
passed to an optimizer. `QuantificationInfo.csf_corr_molar` with f_CSF=1 raises
`ZeroDivisionError: division by zero`; the metabolite integration reference is
also zero. Therefore no fitted metabolite coefficients or water-scaled
concentrations are reported for CSF. Nominal voxel truth is zero, tissue-only
truth and water-scaled estimates are undefined. Its nonzero water acquisition
is independently fitted successfully. This preserves the unsupported case
instead of inventing a tiny reference signal or reporting fabricated zero fits.

## Linearity, repeatability and preservation

Additional mixed-voxel controls double all metabolite concentrations, double
only water concentrations, or double shared receiver gain. Actual fits verify:
metabolite doubling doubles raw coefficients and standard concentrations;
water doubling doubles its fitted amplitude and halves default water-scaled
concentrations; receiver doubling doubles both raw amplitudes while leaving
water-scaled concentrations unchanged. In matched quantification the doubled
actual water concentrations are also declared through densities, so matched
concentrations remain unchanged for the water-doubling control. The complete
repeat fit yields exactly equal parameter dictionaries and water diagnostics.

All baseline input NIfTI-MRS/JSON files are compared byte-for-byte with Phase 6b;
all source and calibration hashes are unchanged. CLI fitted complex arrays
match API fitted arrays exactly (maximum difference zero), and reported
concentrations match to 1e-12 relative. Calibration identity remains
`mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81`,
shared scalar 8857.353925964408. No simulation amplitudes were changed.

## Actual CLI reports and commands

- [WM report](../outputs/mrscloud_mixed_tissue_fitting_validation/final/pure_WM/fsl_mrs_results/report.html)
- [GM report](../outputs/mrscloud_mixed_tissue_fitting_validation/final/pure_GM/fsl_mrs_results/report.html)
- [Mixed report](../outputs/mrscloud_mixed_tissue_fitting_validation/final/mixed/fsl_mrs_results/report.html)

All three files exist, are nonempty and contain water referencing sections.
Complete commands are saved in
[cli_commands.txt](../outputs/mrscloud_mixed_tissue_fitting_validation/final/cli_commands.txt).
Each uses `--report` (singular), verified against this installed CLI. For example:

```sh
pixi run fsl_mrs --data outputs/mrscloud_mixed_tissue_fitting_validation/final/mixed/metabolites.nii.gz --h2o outputs/mrscloud_mixed_tissue_fitting_validation/final/mixed/water_reference.nii.gz --basis outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative/metabolite_basis.BASIS --output outputs/mrscloud_mixed_tissue_fitting_validation/final/mixed/fsl_mrs_results --keep NAA Cr GPC PCh Glu --algo Newton --lorentzian --baseline off --ppmlim 0.2 4.2 --no_conj_basis --TE 30 --TR 2 --tissue_frac 0.2 0.3 0.5 --internal_ref NAA --wref_metabolite NAA --ref_protons 3 --ref_int_limits 1.8 2.2 --report
```

Per-voxel `fsl_mrs_results/` contains actual CLI concentrations, fit arrays,
parameter diagnostics and HTML. `fit_arrays.npz` contains the separately
calculated API complex fitted signal and input. `results.json` contains all
optimizer messages, parameters, scaling, water diagnostics, prediction factors,
source hashes and software version. `comparison.csv` contains 70 rows covering
three baseline nonzero voxels and four controls, five metabolites and two
quantification configurations. CLI logs are empty for the final runs: no
fitting warnings or errors were emitted.

## Source locations and tests

Installed package source under `.pixi/envs/default/lib/python3.14/site-packages/`:

- `fsl_mrs/core/mrs.py:578`, `rescaleForFitting`.
- `fsl_mrs/utils/misc.py:491`, `rescale_FID` (zero-input behavior).
- `fsl_mrs/utils/results.py:258`, raw coefficient scaling; `:707`, `getConc`.
- `fsl_mrs/utils/quantify.py:69`, `FIDIntegrator`; `:92`, `WaterRef`;
  `:212`, `RefIntegral`; `:590` onward, relaxation and CSF properties;
  `:748`, `quantifyWater`, including the molarity equation and reference ratio.
- `fsl_mrs/utils/constants.py`, actual tissue densities, water molarity,
  proton counts and 3T relaxation values.

Fast tests forbid final FSL quantification calls during independent prediction
and verify public overrides without mutating defaults. Separate integration
tests execute all fits and three actual report-generating CLI invocations,
checking coefficients, water amplitudes, independent concentrations, linearity,
repeatability, calibration and immutable inputs. Tolerances are 2e-5 relative
for coefficients and 3e-5 for prediction/report agreement, exceeding measured
1.3e-6 fitting error; residual tolerance 5e-5, water amplitude/residual 1e-6.
Linearity and gain tests are much tighter (1e-9 to 1e-12) because normalized
optimization is identical for these scalar controls. Numerical precision of
fitting is distinct from exact preservation of the original exported samples.

Focused run: **5 passed** in 14.38 seconds. Full regression run: **190 passed,
1 xfailed** in 118.98 seconds; no failures. The existing expected failure
retains the documented installed-package combined phase/frequency water-fit
initialization defect. Warnings are external Dash/SWIG deprecations and a
Jupyter progress-widget warning. They do not indicate fitting failures.

## Limitations and next step

This is a noise-free, matching-basis validation of the known five-component
set, not a test of model selection, noisy identifiability or physical scanner
calibration. GPC/PCh overlap remains relevant. Metabolite DC/formatting and
water raw-acquisition processing remain asymmetric by the validated Phase 6a
convention. Finite-grid integrals, native sequence response and processing
retain the previously established alpha bias. Standard CLI water/relaxation
assumptions do not describe the synthetic no-attenuation tissue profiles;
matched API results and ground-truth volume conventions must remain explicit.
Pure-CSF quantification remains unsupported. No LCModel validation, spatial
simulation, relaxation extension or GUI integration was attempted. Review these
results before extending beyond Phase 6c.

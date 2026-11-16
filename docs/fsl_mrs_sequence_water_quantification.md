# Phase 5h: matched sequence-water quantification assumptions

Matching FSL-MRS's effective water concentration and relaxation assumptions
recovers **8.984750461126422 mM**, rather than 10 mM. The baseline CLI result
is reproduced exactly: **6.885791303512254 mM**. The remaining discrepancy is
predicted by the unit-response spectral integrals, without a fitted empirical
multiplier. No simulator, calibration constant, existing dataset, or installed
FSL-MRS source was changed.

## Inputs, fit, and source trace

FSL-MRS version: `2.5.1.dev78+ga311e38f2`. Inputs are the validated Phase 5g
sequence-mode pair and its matching `common_response.BASIS`: PRESS TE30,
2000 Hz bandwidth, 2048 points, 127.731 MHz, origin 4.68 ppm, NAA 10 mM,
water 40 mol/L, receiver factor 1. Water retains its original +2.55462 Hz
offset, phase, and 1 Hz readout decay. The fixed common normalization scalar
is 2779.198340948381, under the separate sequence calibration identity.

The 10 mM NIfTI-MRS files and JSON sidecars are copied byte-for-byte from
`outputs/mrscloud_sequence_water_integration/validated/mrscloud_sequence`.
Other concentrations use the unchanged generator and reference bundle. Input
and bundle SHA256 hashes are recorded before and after execution; all match.

Fitting uses Newton, Lorentzian, baseline off, NAA only, ppm 0.2–4.2,
and no basis conjugation. The metabolite model is fitted once per concentration.
All four quantification cases use exactly the same fitted parameters.
FSL's public `calculateConcScaling` refits water internally; assertions verify
the resulting full integration-model FID is identical between cases, and
metabolite parameters remain identical. No refit result is substituted by a
desired amplitude.

Installed source locations, relative to
`.pixi/envs/default/lib/python3.14/site-packages/`:

| Source | Function/line | Role |
|---|---|---|
| `fsl_mrs/utils/results.py` | `calculateConcScaling`, 230 | Raw/internal/water scaling; calls water quantification |
| `fsl_mrs/utils/quantify.py` | `FIDIntegrator._calculate_area`, 80 | Half-first-point FFT and abs-real spectral integration |
| same | `WaterRef._fit_w_ref`, 100 | Voigt water model; integration model removes fitted frequency/phase |
| same | `RefIntegral`, 212 | Fitted reference without fitted phase and baseline |
| same | `QuantificationInfo`, 240 onward | Public T1/T2 dictionary overrides |
| same | `set_densities`, 478 | Effective tissue-water density override |
| same | `_calc_relax`, 592 | Saturation and TE attenuation |
| same | `relax_corr_water_molar`, 643; `relax_corr_metab`, 660 | Water concentration and metabolite relaxation factors |
| same | `quantifyWater`, 748 | Integral ratio, proton count, water concentration, corrections |

Full source paths are also recorded in `results.json`. Basis normalization and
coefficient interpretation remain those documented in Phase 5g and the Phase
5d.1 scaling derivation.

## Public quantification configuration

All cases retain **TE=.030 s, TR=2 s**, WM/GM/CSF=1/0/0, 90-degree excitation,
water/reference protons 2/3, water window 1.65–7.65 ppm, NAA window 1.8–2.2
ppm, and additional correction exactly 1. Case A denotes the existing
pure-WM standard configuration, not FSL's no-tissue-input default.

- A: package WM water density .65, default 3T relaxation.
- B: only `set_densities` changes WM density to
  `40000/55010 = 0.7271405199054717`, representing effective water 40 M.
- C: only public T1/T2 dictionaries change to the synthetic no-attenuation
  limit: T1=1e-6 s, T2=infinity for water and metabolites.
- D: both overrides.

T1=1e-6 is a computational representation of instantaneous longitudinal
recovery: at the unchanged TR, exp(-TR/T1) underflows exactly to zero. T2=infinity
makes exp(-TE/T2)=1. These are explicitly synthetic quantification assumptions,
not measured tissue relaxation times. Acquisition metadata are not falsified;
TE/TR are not set to zero/infinity to suppress corrections. The existing 1 Hz
acquisition-time decay and TE-dependent sequence evolution are retained; neither
is an additional tissue TE/T2 attenuation model.

The installed CLI has T1/T2 overrides, but no direct tissue-water-density flag.
The comparison therefore uses `QuantificationInfo.set_densities` and T1/T2
overrides through the public Python API. The CLI HTML report contains **case A
only**. No matched-case concentrations are claimed to appear in that report.

## Independent equation and factor decomposition

For this pure-WM, NAA-self-reference experiment:

\[
R_x=\frac{1-e^{-TR/T_{1,x}}}{1-\cos(FA)e^{-TR/T_{1,x}}}
e^{-TE/T_{2,x}},\qquad
\widehat C_{NAA}[\mathrm{mM}]=\frac{I_m}{I_w}\frac{2}{3}
\underbrace{55010d_{WM}}_{C_{w,eff}[\mathrm{mM}]}
\frac{R_w}{R_m}.
\]

CSF and additional corrections equal 1. The reference coefficient divides out
of FSL's scaling and multiplies back when reporting NAA, so it cancels here.
Integrals are trapezoidal sums of abs(real spectrum), in bin units, using an
orthonormal shifted FFT with sample zero half-weighted and inclusive nearest-bin
window endpoints. They are not peak heights or initial FID amplitudes.

The independent implementation reconstructs metabolite and corrected water
model FIDs from the source basis, fitted nuisance parameters, fitted coefficient,
and captured water optimizer parameters. It independently applies NumPy FFT,
window integration, relaxation and final scalar equations. It does not call
FSL's final reporting function to calculate its prediction. Unit tests replace
`quantifyWater`/`getConc` with exceptions during this independent calculation.

At 10 mM, the fitted, data-scaled integrals are identical in A–D:

- Water: **341832.1888333719**.
- NAA reference: **115.17288435931849**.
- Data fitting multiplier: **562450.3298955784**; basis multiplier:
  **5624.503298955776**. These common receiver/data scales cancel in the ratio.

| Case | Effective water M | R_water | R_metabolite | Concentration mM | Difference from 10 mM |
|---|---:|---:|---:|---:|---:|
| A: defaults | 35.7565 | .5786670548409867 | .6749567292246237 | 6.885791303512254 | -3.114208696487746 |
| B: water only | 40 | .5786670548409867 | .6749567292246237 | 7.702981335994580 | -2.297018664005420 |
| C: relaxation only | 35.7565 | 1 | 1 | 8.031580746581675 | -1.968419253418325 |
| D: both matched | 40 | 1 | 1 | 8.984750461126422 | -1.015249538873578 |

Independent reconstruction gives **6.8857913035122555 mM** for A and
**8.984750461126424 mM** for D. The maximum reconstruction error across all
16 concentration/configuration combinations is **1.0658141036401503e-14 mM**.

The multiplicative prediction from the unfitted unit responses is:

| Contribution | Value |
|---|---:|
| Unit-response/integration factor | .8984749391468103 |
| Fitted source coefficient / true coefficient | 1.000000076436493 |
| Fitted nuisance-response versus unfitted response | 1.0000000426161648 |
| Default effective water / 40 M | .8939125 |
| Default water/metabolite relaxation | .8573394853115213 |

Multiplying 10 mM by these five factors reconstructs A. Omitting only the last
two gives D. Correcting density alone multiplies A by 40000/35756.5; correcting
relaxation alone multiplies A by .6749567292246237/.5786670548409867.

## Remaining response and integration mismatch

On the same target grid, let B be the unit-concentration NAA basis and W be
the unit-concentration raw-water envelope in the corrected integration frame.
The **unfitted**, response-derived factor is

\[
\alpha=\frac{2}{3}\frac{J[B;(1.8,2.2)]}{J[W;(1.65,7.65)]}
=\frac{2}{3}\frac{.020476985425512653}{.015193883571167494}
=.8984749391468103.
\]

This is calculated before using the fitted concentration. Both unit responses
contain the same division by S, which cancels; there is no unknown independent
water normalization in this factor. The existing two-proton H2O response and
FSL's integer 2/3 convention are retained unchanged.

Three independent response controls decompose alpha:

| Control/contribution | Factor |
|---|---:|
| Ideal 3H singlet at 2.01 ppm versus 2H water, same 1 Hz decay and windows | .9577941910199874 |
| Native exact-time NAA basis response versus that ideal 3H control | .9353584023864411 |
| Formatted versus native exact-time NAA integral | 1.002895580183198 |
| Product | .8984749391468103 |

The ideal control uses amplitude 1.5 times the water unit amplitude, by its
explicit 3H/2H definition, solely as a diagnostic. It does not rescale any
simulated dataset. Its non-unit integral ratio demonstrates finite-grid,
frequency-position and unequal-window effects even for perfectly defined
proton responses. Formatting raises the NAA integral by about 0.29%, rather
than explaining the entire deficit. Native exact-time NAA integral is
.020417863863526193; the formatted integral is .020476985425512653.

The native comparison still includes MRSCloud sequence response, intrinsic
phase, other NAA resonances and prior basis DC processing. Those contributions
cannot be separated into a unique physical excitation-versus-proton correction
from these stored responses: the bundle lacks a pre-DC isolated NAA methyl
capture. Therefore the .9353584 factor is explicitly an aggregate native-basis
response difference, not proof of a physical proton-count error. A further
isolated methyl/pre-DC experiment would be needed to attribute that portion.
The overall quantitative mismatch, however, is fully predicted numerically.

Original acquisition-frame input integrals at 10 mM are .00020476985425512656
for NAA and .5850274618881485 for water. They differ from fitted integration
values because the latter use fitting-scale units and remove fitted water
frequency/phase. This is the previously documented reporting convention, not
evidence of an amplitude loss.

No inverse-alpha correction is applied. Matching water content and relaxation
alone is insufficient for native standard quantification to recover 10 mM.
A future response-convention decision must address reference integration and
metabolite representation explicitly rather than changing water amplitudes
against the desired reported concentration.

## Amplitude recovery and linearity

At 10 mM the raw NAA source coefficient is **.01000000076436493 mol/L** and
the internal fitted coefficient is **1.0000000764364945**. Water initial
amplitude is **.02577258571845567**, fitted amplitude **.02577258514572865**.
Water complex residual is **1.4441623108503362e-7**, metabolite residual
**1.857590496917146e-5**. Both optimizers report success.

| Truth mM | Fully matched recovered mM | Signed error mM | Relative error |
|---:|---:|---:|---:|
| 5 | 4.492375230563211 | -.507624769436789 | -.101524953887358 |
| 10 | 8.984750461126422 | -1.015249538873578 | -.101524953887358 |
| 15 | 13.477126002855030 | -1.522873997144970 | -.101524933142998 |
| 20 | 17.969500922252845 | -2.030499077747155 | -.101524953887358 |

Maximum relative deviation of recovered/truth from the 10 mM slope is
**2.3088409406568644e-8**. This is accurate linearity with a systematic
response-convention offset, not accurate absolute concentration recovery.

## Reproduction, reports, and tests

Files added: `examples/fsl_mrs_sequence_water_quantification.py`,
`tests/test_sequence_water_quantification.py`,
`tests/integration/test_sequence_water_quantification_e2e.py`, and this document.
Production code and dependencies are unchanged.

Executed final experiment and regression commands:

```bash
pixi run python --version
pixi run fsl_mrs --help
pixi run python -m examples.fsl_mrs_sequence_water_quantification outputs/fsl_mrs_sequence_water_quantification/final
pixi run python -m pytest -q tests/test_sequence_water_quantification.py tests/integration/test_sequence_water_quantification_e2e.py
pixi run python -m pytest -q tests
```

The baseline command actually executed by the example is:

```bash
pixi run fsl_mrs --data outputs/fsl_mrs_sequence_water_quantification/final/NAA_10mM/metabolites.nii.gz --h2o outputs/fsl_mrs_sequence_water_quantification/final/NAA_10mM/water_reference.nii.gz --basis outputs/mrscloud_sequence_water_integration/validated/reference/common_response.BASIS --output outputs/fsl_mrs_sequence_water_quantification/final/NAA_10mM/fsl_mrs_results --keep NAA --algo Newton --lorentzian --baseline off --ppmlim 0.2 4.2 --no_conj_basis --TE 30 --TR 2 --tissue_frac 1 0 0 --internal_ref NAA --wref_metabolite NAA --ref_protons 3 --ref_int_limits 1.8 2.2 --report
```

[Actual baseline HTML report](../outputs/fsl_mrs_sequence_water_quantification/final/NAA_10mM/fsl_mrs_results/report.html)
is nonempty (887839 bytes) and reports **A only**. API comparison results are
separate in [comparison.csv](../outputs/fsl_mrs_sequence_water_quantification/final/comparison.csv)
and [results.json](../outputs/fsl_mrs_sequence_water_quantification/final/results.json).
Each concentration directory contains fitted parameters and model arrays;
the root contains source hashes, exact CLI command, and configuration.

![API quantification comparison and response controls](../outputs/fsl_mrs_sequence_water_quantification/final/quantification_comparison.png)

Focused tests executed **5 passed**. An initial focused invocation failed
collection because unit and integration files had the same basename; renaming
the integration file fixed collection without changing scientific expectations.
The full regression run completed with **112 passed, 1 xfailed, 5 dependency
warnings in 70.24 s**. The expected failure is the previously documented
installed-package combined phase/frequency water-initialization defect.
The suite includes analytical, Phase 3, frozen-basis, calibration, paired
export, processing, CLI-report, and quantification regressions. New tests cover
independent factor reconstruction, unchanged fitted amplitudes and source
files, public override isolation, concentration linearity, actual baseline
HTML generation, and exactly deterministic repeat fitting/quantification.

## Scientific limits

This establishes the numerical interpretation of these fixed synthetic signals,
not physical scanner calibration or universal cross-metabolite concentration
accuracy. Remaining issues include basis-specific DC processing, intrinsic
phase and sequence response, reference-window dependence, ideal excitation,
and physical acquisition/filter modeling. LCModel validation, simulator
corrections, tissue mixing, and GUI integration are outside this phase.

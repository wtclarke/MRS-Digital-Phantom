# Phase 6g: FSL-MRS fitting under spatial B0 inhomogeneity

## Findings

The zero-field signal is accurately fitted. A constant +37 Hz shift preserves
all true amplitudes but the default optimizer lands in a wrong local minimum;
known-frequency initialization recovers the exact metabolite signal. Standard
water quantification still falls sharply because its fixed NAA integration
window misses the shifted reference peak. Symmetric and tissue-dependent
frequency distributions are poorly represented by one shared metabolite
Lorentzian and one water Voigt. Their optimizer success flags do not establish
accurate amplitudes or concentrations.

No simulated amplitude, production module, original input, calibration constant,
raw water reference, MRSCloud source or installed package was modified.
Dedicated fitting-basis alignment changes frequency only, not amplitudes.
All Phase 6a and primary Phase 6f input hashes are unchanged.

## Source datasets and fit settings

The four actual Phase 6f input files per case are copied byte-for-byte from
`outputs/mrscloud_quantitative_spatial/validated`. Independent signals are also
calculated directly from cell labels, weights, nominal tissue profiles and
input T1/T2 equations, without using the simulator's returned modulation
functions or attenuation factors. Complex maximum differences from actual
stored/interpreted signals are below 2e-21 for metabolites and 2e-18 for water.

| Case | Tissue / normalized weights | Cell B0 Hz | Relaxation |
|---|---|---|---|
| A, zero_B0 | WM/WM, .5/.5 | 0/0 | none |
| B, constant_B0 | WM/WM, .5/.5 | +37/+37 | none |
| C, symmetric_B0 | WM/WM, .5/.5 | -5/+5 | none |
| D, mixed_relaxation_B0 | WM/GM/CSF, .2/.3/.5 | -3/+7/+1 | Phase 6e profiles |

Metabolites are NAA, Cr, GPC, PCh and Glu. Pure WM nominal values are
10/8/2/1/8 mM and water 39.0165 mol/L. Mixed nominal voxel concentrations are
6.2/4.3/.85/.56/5.2 mM, WM/GM-only concentrations 12.4/8.6/1.7/1.12/10.4 mM,
and nominal water 49.72245 M. Mixed relaxation-weighted effective coefficients
are 4.0101017306/2.8461097902/.5875894594/.3908850361/3.4308773164 mM-equivalent;
water is 21.7989112951 M-equivalent. These remain distinct from nominal truth.

The shared receiver is G=1, no noise. Sequence/basis is the unchanged full
PRESS TE30 Phase 6a bundle, S=8857.353925964408, identity
`mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81`.
Sampling is 2048 complex points, 2000 Hz, dwell .0005 s, 127.731 MHz,
acquisition reference 4.68 ppm, intrinsic water +2.55462 Hz.

Installed FSL-MRS version: **2.5.1.dev78+ga311e38f2**. All metabolite fits use Newton/TNC,
Lorentzian, baseline off, one shared frequency/linewidth group, fitted zero-
and first-order phase, no conjugation, fixed ppm 0.2–4.2. WaterRef separately
fits a single Voigt even when the metabolite CLI uses `--lorentzian`.
Quantification uses actual fractions, TE30/TR2, installed default 3T relaxation
and water densities, NAA 3-proton reference window 1.8–2.2 ppm, water 2 protons
and corrected-model window 1.65–7.65 ppm. No water gain correction is used.
These default quantification assumptions are intentionally not described as
fully matching the arbitrary metabolite-specific relaxation in case D.

## Constant-shift initialization and actual CLI adaptation

Default case B returns optimizer success with full-FID residual
0.96997109155, NAA coefficient zero, frequency about +0.4300 Hz rather than
+37 Hz, and large phase terms. Its complete original parameters and status are
preserved in `default_fit` in `results.json`. This is a local-minimum failure,
not attenuation or an amplitude-calibration defect. Frequency bounds are
unrestricted in the installed Lorentzian Newton model; no bounds were changed.

The Python retry uses public `fit_FSLModel(x0=...)`, setting the frequency from
the known +37 Hz source configuration and initial amplitudes from independent
nonnegative least squares on the shifted basis. It does not insert ground-truth
concentrations as initial amplitudes. Gamma and phase start at zero. The
original matching basis and data are retained. Residual becomes 2.958e-15,
with amplitude recovery to about 2e-15 relative.

The installed CLI has no x0/initial-frequency option. Its case B report therefore
uses a separate **target-grid fitting basis**, with every original formatted
basis column multiplied by exp(+i*2*pi*37*t). All 29 columns are retained;
there is no independent column scaling. Native-to-target formatting is done
before shifting, because shifting before Fourier formatting would introduce a
different boundary operation. Public FSL `Basis.save` writes the target-grid
complex columns without amplitude normalization. Reload error is exactly zero;
maximum sample-magnitude change from floating-point rotation is 2.169e-19.
The original Phase 6a bundle and input FIDs are untouched.

This adaptation is explicit in the command's `--basis` path and in result
metadata (`CLI_basis_frequency_alignment_hz=37`). CLI residual is 1.199e-7;
maximum complex difference from the accurately initialized API fit is 2.050e-12
in original signal units. Other cases use the original BASIS and CLI/API
fits match exactly. The aligned CLI's additional frequency parameter is near
zero **relative to the already +37 Hz fitting basis**; it must not be read as
zero absolute shift. The primary API parameter is exactly +37 Hz in the source
basis frame. This is a controlled known-shift validation, not blind B0 recovery.

## Metabolite coefficients and optimizer diagnostics

Let a be the internal fitted coefficient, s_D the data scale and s_B the
software's common fitting-basis scale. The expected initial coefficient is
G*1e-3*C_m_effective*s_D/s_B. The table converts fitted values back through
1000*a*s_B/(s_D*G), to compare against pre-acquisition-attenuated effective mM.
Neither internal a nor this conversion is automatically a nominal concentration
under mismatch. All internal values, expected values and scalings are saved.

| Case | Metabolite | True effective mM | Fitted effective mM | Relative coefficient error % | Fixed fitted-nuisance projection mM |
|---|---|---:|---:|---:|---:|
| zero_B0 | NAA | 10.000000000 | 9.999999416 | -5.839333e-06 | 10.000000088 |
| zero_B0 | Cr | 8.000000000 | 7.999999598 | -5.028455e-06 | 7.999999999 |
| zero_B0 | GPC | 2.000000000 | 1.999999582 | -2.08763e-05 | 1.999999955 |
| zero_B0 | PCh | 1.000000000 | 1.000000464 | 4.640896e-05 | 1.000000028 |
| zero_B0 | Glu | 8.000000000 | 8.000001965 | 2.456692e-05 | 7.999999682 |
| constant_B0 | NAA | 10.000000000 | 10.000000000 | -5.551115e-14 | 10.000000000 |
| constant_B0 | Cr | 8.000000000 | 8.000000000 | -4.440892e-14 | 8.000000000 |
| constant_B0 | GPC | 2.000000000 | 2.000000000 | 2.220446e-14 | 2.000000000 |
| constant_B0 | PCh | 1.000000000 | 1.000000000 | 1.776357e-13 | 1.000000000 |
| constant_B0 | Glu | 8.000000000 | 8.000000000 | 1.332268e-13 | 8.000000000 |
| symmetric_B0 | NAA | 10.000000000 | 5.379505675 | -46.20494 | 5.379505768 |
| symmetric_B0 | Cr | 8.000000000 | 4.273959871 | -46.5755 | 4.273958997 |
| symmetric_B0 | GPC | 2.000000000 | 1.272059534 | -36.39702 | 1.272062769 |
| symmetric_B0 | PCh | 1.000000000 | 0.374271497 | -62.57285 | 0.374268187 |
| symmetric_B0 | Glu | 8.000000000 | 6.668960578 | -16.63799 | 6.668961531 |
| mixed_relaxation_B0 | NAA | 4.010101731 | 2.751936627 | -31.37489 | 2.751936811 |
| mixed_relaxation_B0 | Cr | 2.846109790 | 1.805755625 | -36.55355 | 1.805755922 |
| mixed_relaxation_B0 | GPC | 0.587589459 | 0.359882784 | -38.75268 | 0.359882752 |
| mixed_relaxation_B0 | PCh | 0.390885036 | 0.226129882 | -42.14926 | 0.226129996 |
| mixed_relaxation_B0 | Glu | 3.430877316 | 3.048000771 | -11.15973 | 3.048000642 |

| Case | Metabolite full complex residual | Frequency Hz | Additional gamma/pi Hz | Phi0 (FSL spectral sign, rad) | Phi1 rad/Hz | Optimizer success |
|---|---:|---:|---:|---:|---:|---|
| zero_B0 | 1.14928379e-07 | -4.80638532e-08 | 9.24938748e-20 | -1.40263909e-08 | 3.40463884e-10 | True |
| constant_B0 | 2.9583282e-15 | 37 | 0 | 0 | 0 | True |
| symmetric_B0 | 0.674673313 | 5.09403448 | 0.143387514 | 0.16658938 | -7.517542e-05 | True |
| mixed_relaxation_B0 | 0.48459849 | 7.05149436 | 0.0644503622 | 0.0969817729 | -2.08381151e-05 | True |

Gamma is additional acquisition-time broadening on the processed basis, not
TE/T2 attenuation or a measured total spectral FWHM. Phi0/Phi1 enter the
spectral model as exp(-i*(Phi0+Phi1*frequency)); they are not directly the water
FID's positive phase parameter. All selected metabolite optimizers report
convergence, including the inadequate dispersed fits. Reported solutions are
selected local fits, not a demonstrated global nonlinear optimum.

Case C's model locks near the +5 Hz component. Case D locks near +7 Hz, where
the GM metabolite contribution is strong. It cannot represent the other cells
with one shared offset and monotonic line shape. Each individual coefficient
is checked; a fit summary alone would hide the GPC/PCh/Glu differences.

## Water amplitudes, phase, frequency and line shape

| Case | True initial amplitude | Fitted amplitude | Amplitude error % | Frequency Hz | Phase rad | Lorentzian gamma/pi Hz | Gaussian sigma s^-2 | Full complex residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| zero_B0 | 0.00788790829157 | 0.00788790812234 | -2.145412e-06 | 2.55462014 | -1.42262798e-07 | 0.999999938 | 0 | 1.43238671e-07 |
| constant_B0 | 0.00788790829157 | 0.00788789877029 | -0.0001207073 | 39.55462 | -1.10050855e-08 | 0.999995021 | 2.62390028e-05 | 1.01323861e-06 |
| symmetric_B0 | 0.00788790829157 | 0.00418093704268 | -46.99562 | -2.54899083 | 0.197520687 | 1.08255611 | 0 | 0.695736719 |
| mixed_relaxation_B0 | 0.00440705376321 | 0.00237438531283 | -46.12307 | 3.52992921 | 0.0356055373 | 1.28237171 | 0 | 0.638895304 |

All water optimizers report a local minimum. The Lorentzian column is only one
model component; Gaussian FWHM, if required, is 2*sqrt(log(2)*sigma)/pi Hz.
It is not the measured width of the original multiple-peak data.

After removing the **known** carrier, case C has envelope exp(-pi*t)*cos(2*pi*5*t).
It has zeros at 50 ms and sign reversals. A finite-amplitude single Voigt's
positive exponential magnitude cannot represent these zeros/revivals. FSL
instead selects approximately one lobe (frequency -2.548991 Hz versus actual
components -2.44538 and +7.55462 Hz), recovering about 53% of the true initial
amplitude with a 69.6% residual. This is actual line-shape mismatch, beyond the
previously documented HTML frequency/phase overlay issue.

Case D has water components at -.44538, +3.55462 and +9.55462 Hz. Their
water/relaxation-weighted true centroid is 4.768350 Hz, but a single center is
not a complete description. FSL selects 3.529929 Hz near the strongest CSF
component, underestimates initial amplitude by 46.1%, and leaves a 63.9%
residual. Metabolite and water contributions weight the same cell offsets
differently; different selected peaks are expected under this inadequate model.

Water plots and residuals use the fitted model restored to the **original
complex data frame**, unlike the standard HTML's integrated-model overlay.
No input phase/frequency correction is applied to improve the fit. The known
installed wrapped-phase initialization defect remains in existing regressions;
no package code was modified or water initialization changed in these cases.

## Independent least-squares diagnostics

The validation independently implements the half-first-point orthonormal FFT,
then stacks real and imaginary spectral samples within the fitting window.
SciPy nonnegative least squares solves coefficients on fixed designs:

1. **Known-offset, zero additional broadening:** uses known +37 Hz for B and
   spatial-weighted mean offsets for dispersed cases. Zero and constant controls
   recover truth; dispersed residuals are about .984 and .982.
2. **FSL fitted-nuisance design:** applies the fitted gamma, frequency and phase
   to the unit responses, then independently solves amplitudes. This is a
   conditional diagnostic, explicitly using fitted nuisance parameters; it is
   not an independent prediction of an unbiased concentration. It reproduces
   the biased coefficients, showing they are close to the least-squares optimum
   within those restricted line shapes. Remaining TNC/NNLS coefficient
   differences are below 1e-5 relative in these cases.
3. **Known spatial-response design (oracle control):** exact independently
   calculated spatial column shapes use declared tissue profiles and relaxation.
   Coefficients recover effective truth within about 1e-15 relative, with
   residuals below 8e-16 for metabolites and 3e-15 for water. This deliberately
   uses known within-metabolite tissue composition, so it is not a blind
   concentration-estimation pipeline or new calibration.

Water fixed-nuisance projection independently solves the amplitude of the
fitted Voigt shape. It reproduces the biased optimizer amplitudes; replacing
that shape with the known spatial water response recovers true initial
amplitude. Thus amplitude loss is not present in the input data.

## Water-referenced concentrations and predictions

FSL's NAA-reference concentration for target m is

\[
C_m^{FSL}=rac{a_m}{a_{NAA}}rac{J[s_{NAA,model}]}{J[s_{water,model,corrected}]}
rac{2}{3}rac{55010\sum_j f_jd_jR^Q_{w,j}}{(1-f_{CSF})R^Q_M}.
\]

J is nearest-inclusive-window trapezoidal integration of abs(real(FFT)), with
first-point half weighting and orthonormal transform, in spectral-bin units.
The reference NAA model removes fitted global phase, but **retains frequency
shift and broadening**. WaterRef removes its fitted frequency and phase.
Default tissue densities are .65/.78/.97; TE/TR and standard 3T T1/T2 values
are the previously documented quantification assumptions. Mixed CSF correction
is 2; nominal tissue-only values must not be called voxel averages.

Two different calculations are retained and labeled:

- **Unfitted source-response reference:** compute the NAA component directly
  from known cells and integrate it in the unchanged reference window. Align
  true water by its analytically computed water/relaxation-weighted centroid
  and initial phase, then numerically integrate its full spatial line shape.
  Use known effective metabolite ratios and declared default quantification
  factors. This is an ideal-model prediction for A/B. For C/D it is an explicit
  numerical counterfactual/reference, not an exact prediction of what a
  misspecified one-Voigt fitting algorithm will return. Its water-centroid
  frame is declared and not chosen from fitted concentrations.
- **Independent software reconstruction:** use the selected fitted coefficients
  and nuisance parameters to reconstruct model FIDs directly, then perform J
  and the equation above without calling final concentration reporting. This
  conditional reconstruction matches every reported value within 5.33e-15 mM.
  It checks software interpretation, not scientific accuracy of the biased fit.

No zero-B0 integral ratio is reused for dispersed data or a shifted NAA peak.

| Case | True NAA integral | True water integral, centroid frame | Fitted NAA integral | Fitted corrected water integral |
|---|---:|---:|---:|---:|
| zero_B0 | 6.42512458776e-05 | 0.186008437821 | 6.42512419838e-05 | 0.186008435618 |
| constant_B0 | 8.91904942167e-07 | 0.186008437821 | 8.91904942167e-07 | 0.186008119783 |
| symmetric_B0 | 6.47198057029e-05 | 0.183968912142 | 3.59743054941e-05 | 0.0974891521216 |
| mixed_relaxation_B0 | 2.61113816854e-05 | 0.105042209933 | 1.84818741924e-05 | 0.0544342184817 |

| Case | Metabolite | Nominal voxel mM | Tissue-only mM | Effective mM | Source-response reference mM | API standard FSL mM | Actual CLI mM |
|---|---|---:|---:|---:|---:|---:|---:|
| zero_B0 | NAA | 10 | 10 | 10.000000000 | 7.059362383 | 7.059362039 | 7.059362039 |
| zero_B0 | Cr | 8 | 8 | 8.000000000 | 5.647489907 | 5.647489677 | 5.647489677 |
| zero_B0 | GPC | 2 | 2 | 2.000000000 | 1.411872477 | 1.411872196 | 1.411872196 |
| zero_B0 | PCh | 1 | 1 | 1.000000000 | 0.705936238 | 0.705936573 | 0.705936573 |
| zero_B0 | Glu | 8 | 8 | 8.000000000 | 5.647489907 | 5.647491348 | 5.647491348 |
| constant_B0 | NAA | 10 | 10 | 10.000000000 | 0.097994679 | 0.097994846 | 0.097994824 |
| constant_B0 | Cr | 8 | 8 | 8.000000000 | 0.078395743 | 0.078395877 | 0.078395862 |
| constant_B0 | GPC | 2 | 2 | 2.000000000 | 0.019598936 | 0.019598969 | 0.019598954 |
| constant_B0 | PCh | 1 | 1 | 1.000000000 | 0.009799468 | 0.009799485 | 0.009799493 |
| constant_B0 | Glu | 8 | 8 | 8.000000000 | 0.078395743 | 0.078395877 | 0.078395857 |
| symmetric_B0 | NAA | 10 | 10 | 10.000000000 | 7.189676223 | 7.541412449 | 7.541412449 |
| symmetric_B0 | Cr | 8 | 8 | 8.000000000 | 5.751740978 | 5.991571740 | 5.991571740 |
| symmetric_B0 | GPC | 2 | 2 | 2.000000000 | 1.437935245 | 1.783272699 | 1.783272699 |
| symmetric_B0 | PCh | 1 | 1 | 1.000000000 | 0.718967622 | 0.524683102 | 0.524683102 |
| symmetric_B0 | Glu | 8 | 8 | 8.000000000 | 5.751740978 | 9.349071338 | 9.349071338 |
| mixed_relaxation_B0 | NAA | 6.2 | 12.4 | 4.010101731 | 9.998902004 | 13.657157771 | 13.657157771 |
| mixed_relaxation_B0 | Cr | 4.3 | 8.6 | 2.846109790 | 7.096571309 | 8.961503411 | 8.961503411 |
| mixed_relaxation_B0 | GPC | 0.85 | 1.7 | 0.587589459 | 1.465112313 | 1.786006230 | 1.786006230 |
| mixed_relaxation_B0 | PCh | 0.56 | 1.12 | 0.390885036 | 0.974643895 | 1.122224779 | 1.122224779 |
| mixed_relaxation_B0 | Glu | 5.2 | 10.4 | 3.430877316 | 8.554647333 | 15.126448412 | 15.126448412 |

### Separating three causes

- A reproduces the prior response/integration and default water/relaxation
  convention: NAA about 7.059362 mM for true 10 mM. Its largest individual
  coefficient error is only 4.65e-7 relative. That concentration difference
  is not B0 model mismatch.
- B has correctly recovered initial coefficients, yet NAA standard concentration
  is .097994846 mM (source prediction .097994679). The shifted methyl peak is
  near 2.30 ppm, outside fixed 1.8–2.2 ppm. The true NAA reference integral drops
  from 6.4251e-5 to 8.9190e-7. This analytically explains the additional loss.
  No reference-window correction, data demodulation or calibration adjustment
  was implemented. CLI/API differences are small optimizer roundoff (maximum
  relative concentration difference below 1e-6 in B).
- C/D add real model mismatch. NAA becomes 7.541412/13.657158 mM versus numerical
  source-response references 7.189676/9.998902. Water and metabolite amplitude
  biases can partly cancel or amplify in a ratio; plausible NAA concentration
  does not demonstrate accurate individual coefficients. D's result must also
  be compared with its 12.4 mM tissue-only and 6.2 mM voxel truth separately.

A globally optimal nonlinear Lorentzian/Voigt solution was not established;
these are reproducible actual optimizer solutions and conditional projections.
Their large residuals and exact known-spatial controls suffice to show the
selected model does not represent the data. No empirical correction is applied.

## Reports, plots and reproduction

Actual nonempty CLI HTML reports, all generated with `--report`:

- [A: zero B0](../outputs/mrscloud_spatial_b0_fitting_validation/validated/zero_B0/fsl_mrs_results/report.html)
- [B: constant B0, aligned fitting basis](../outputs/mrscloud_spatial_b0_fitting_validation/validated/constant_B0/fsl_mrs_results/report.html)
- [C: symmetric B0](../outputs/mrscloud_spatial_b0_fitting_validation/validated/symmetric_B0/fsl_mrs_results/report.html)
- [D: mixed B0](../outputs/mrscloud_spatial_b0_fitting_validation/validated/mixed_relaxation_B0/fsl_mrs_results/report.html)

All report sections are checked for water referencing. Their default integrated
water overlay may still use a different frame; the following diagnostic figure
uses original-frame complex models and was inspected directly:

![Complex input and fitted spectra](../outputs/mrscloud_spatial_b0_fitting_validation/validated/complex_spectral_fits.png)

![Coefficient, water amplitude, concentration errors and residuals](../outputs/mrscloud_spatial_b0_fitting_validation/validated/bias_and_residuals.png)

![Complex residuals and known-centroid water envelopes](../outputs/mrscloud_spatial_b0_fitting_validation/validated/residuals_and_water_envelopes.png)

Water time traces divide by the same true initial amplitude **for display only**;
this normalization is not used in numerical fitting or quantification.
[comparison.csv](../outputs/mrscloud_spatial_b0_fitting_validation/validated/comparison.csv)
contains all 20 metabolite rows. `results.json` retains versions, optimizer
messages, original failed constant-fit parameters, nuisance terms, projection
results, integration factors, source hashes and fitting-basis provenance.
Each case has original NIfTI-MRS inputs/sidecars, actual CLI output and an NPZ
with original, independently reconstructed and fitted complex signals.

Commands executed using Pixi, with an unused output directory:

```sh
pixi run python --version
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run fsl_mrs --help
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.mrscloud_spatial_b0_fitting_validation outputs/mrscloud_spatial_b0_fitting_validation/validated
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_spatial_b0_fitting_diagnostics.py tests/integration/test_spatial_b0_fitting_validation.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q
```

All four exact CLI commands are in
[cli_commands.txt](../outputs/mrscloud_spatial_b0_fitting_validation/validated/cli_commands.txt).
For case B the complete actual command was:

```sh
pixi run fsl_mrs --data outputs/mrscloud_spatial_b0_fitting_validation/validated/constant_B0/metabolites.nii.gz --h2o outputs/mrscloud_spatial_b0_fitting_validation/validated/constant_B0/water_reference.nii.gz --basis outputs/mrscloud_spatial_b0_fitting_validation/validated/constant_shifted_fitting_basis --output outputs/mrscloud_spatial_b0_fitting_validation/validated/constant_B0/fsl_mrs_results --keep NAA Cr GPC PCh Glu --algo Newton --lorentzian --baseline off --ppmlim .2 4.2 --no_conj_basis --TE 30 --TR 2 --tissue_frac 1.0 0.0 0.0 --internal_ref NAA --wref_metabolite NAA --ref_protons 3 --ref_int_limits 1.8 2.2 --report
```

Source locations in installed FSL-MRS under `.pixi/envs/default/lib/python3.14/site-packages/`:
`utils/fitting.py:24` (public x0), `models/model_lorentzian.py:43` (unrestricted
frequency bounds), `utils/quantify.py:92` (WaterRef), `:212` (RefIntegral),
`models/__init__.py:282` (no_phase zeros phase only), `utils/quantify.py:748`
(water quantification), `core/basis.py:236` (public fitting-basis save).

Initial exploratory API results are retained in the output root's `pilot.json`.
An initial script import used an unexported `fsl_mrs.core.Basis` symbol; using
its actual `fsl_mrs.core.basis.Basis` location fixed it. A fast linear-algebra
test initially assumed sub-ulp coefficient equality; its measured 1.7e-15
relative NNLS roundoff is now tested with 5e-15 tolerance. No data, amplitudes,
expected physical concentrations or genuine mismatch errors were altered.

Focused numerical/CLI run: **7 passed** in 18.56 s. Tests verify source-derived
signals, exact oracle projections, zero/constant amplitude recovery, retained
failed default initialization, large dispersed-fit residuals and biases,
conditional least-squares projections, independent quantification reconstruction,
all four reports, calibration/source immutability and deterministic refitting.
Oracle/projection precision and optimizer tolerances are kept separate;
large model mismatch is explicitly asserted rather than suppressed or xfailed.

Full regression: **260 passed, 1 xfailed**, 4 dependency warnings in
178.41 seconds. No failures. The existing expected failure remains the
installed WaterRef combined phase/frequency initialization defect. Warnings
are external Dash/SWIG deprecations. This full run also executes the final
residual/water-envelope plot assertion and all preceding phase regressions.

## Remaining limitations

This uses five known metabolites, noise-free matching sequence responses and
an acquisition-only spatial B0 model. It does not model B0-dependent RF or
localization efficiency. Individual GPC/PCh overlap and nonlinear local minima
remain relevant. Known-frequency initialization and truth-informed oracle
shapes are controlled diagnostics, not an unknown-data clinical pipeline.
The fixed quantification windows, asymmetric basis/water processing, tissue
water defaults, scalar metabolite relaxation correction and native response
bias remain. Standard reported mM under these conditions is not a demonstrated
physical concentration calibration. No production integration, GUI changes,
noise, empirical gain correction or LCModel validation was attempted. Stop at
Phase 6g and review model adequacy and reference-window behavior before further
pipeline work.

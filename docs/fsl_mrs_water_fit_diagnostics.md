# Phase 5b.1: water frequency, phase and integration diagnostics

## Findings

There are two distinct issues in installed FSL-MRS
`2.5.1.dev78+ga311e38f2`:

1. For the actual Phase 5b water (zero phase, +2.55462 Hz), the full complex
   fitted model matches the data with relative residual **8.32619e-10**. The
   report overlays the original data with the fitted model **after deliberately
   removing its fitted frequency and phase** for integration. These curves are
   in different reference frames. The apparent shifted absorptive fit in the
   HTML is therefore not the actual fitted curve.
2. A controlled nonzero-phase variant exposes a genuine initialization failure.
   Phase 0.6 rad plus the same frequency can converge to amplitude zero and
   residual 1 while the optimizer reports success. Unwrapped phase regression
   with the frequency sign required by the model recovers it without changing
   any signal, calibration constant, bounds or optimizer settings.

No production simulator, installed FSL-MRS code or calibration was modified.
Initialization replacement is scoped to diagnostic calls through a temporary
`unittest.mock.patch` of the optimizer entry point.

## Direct sample measurements and export trace

The original Phase 5b 10 mM dataset was loaded from
`outputs/fsl_mrs_water_scaling_phase5b/NAA_10mM`. Generator output was independently
recreated with the same Phase 5a configuration. Frequency was measured by a
linear fit to unwrapped complex phase, linewidth by
−slope(log(abs(FID)))/pi, initial amplitude by abs(s[0]), and phase by arg(s[0]).
There are 2048 samples, dt=0.0005 s, t[0]=0.

| Stage | Initial amplitude | Phase rad | Frequency Hz | Effective linewidth Hz |
|---|---:|---:|---:|---:|
| Generator | 61.73993435001698 | 0 | +2.5546200000000545 | 4.02e-16 |
| Raw NIfTI storage | 61.73993435001698 | 0 | −2.5546200000000545 | 4.02e-16 |
| FSL interpreted FID | 61.73993435001698 | 0 | +2.5546200000000545 | 4.02e-16 |
| Water-fit input after MRS scaling | 7274.926060620617 | 0 | +2.5546200000000545 | −3.28e-15 |

Near-zero linewidth measurements are numerical roundoff, consistent with the
configured zero linewidth. MRS applies a positive real multiplier
117.8317751259257 to both metabolite and water FIDs. Full-array comparisons:

- Stored versus generator: max complex difference **123.47985096532967**,
  because NIfTI-MRS stores the conjugate by the validated library convention.
- Stored versus conjugated generator: **0**.
- Interpreted versus generator: **0**.
- Fitting input divided by the MRS scaling versus generator: **1.00486e-14**.

The fitter subsequently applies its own temporary water scaling, which is
explicitly undone when reconstructing its fitted amplitude and actual model.
The trace therefore tests phases/frequencies at each boundary, not merely a
successful file reload. Full arrays are saved in `signal_trace.npz`.

## Chemical shift and FFT convention

The checked-in PRESS TE30 basis uses 4.68 ppm as its origin. For interpreted
complex FIDs the simulator uses exp(+i*2*pi*f*t). Both the repository FFT and
FSL's FIDToSpec use the forward FFT with negative exponent and fftshift. FSL
additionally halves the first FID sample and uses orthonormal FFT scaling.
That changes spectral amplitude conventions, not the rotation frequency.

Installed `nifti_mrs.axes.Axes` uses `fftshift(fftfreq(N,dt))`, with chemical
shift ppm = 4.68 + f/127.731. The numerical axis increases with array index;
plots reverse their horizontal axis to show high ppm on the left. Thus

\[
4.68 + 2.55462/127.731 = 4.70\ \mathrm{ppm}.
\]

The largest sampled magnitude bin is at 4.702936385842121 ppm (+2.9296875 Hz),
within one bin of 4.70 ppm; the zero-frequency bin is at 4.68 ppm. No ppm or
frequency sign correction to the generated water is warranted.

The strong dispersive appearance of the **real** spectrum is reproducible
despite zero initial phase: zero linewidth means the FID never decays before
the finite record ends. At this noninteger FFT-bin frequency, its rectangularly
truncated exponential has a complex Dirichlet spectrum. The phase of that
finite sum varies with frequency, so the real part has positive and negative
lobes. At zero frequency the same undamped FID falls exactly on a bin and
appears absorptive. This is a consequence of the deliberately idealized finite
signal, not evidence of an extra initial-phase error. Artificial broadening or
phase offsets would change the experiment and were not used as a remedy.

## Actual installed fitting and plotting code

Relevant installed source locations under
`.pixi/envs/default/lib/python3.14/site-packages/fsl_mrs/`:

- `utils/quantify.py::WaterRef._fit_w_ref`, lines 100–196: Voigt fitting.
- `utils/quantify.py::FIDIntegrator._calculate_area`, lines 80–89: integration.
- `utils/plotting.py::plot_references`, around lines 1223–1260: report overlay.
- `utils/misc.py::FIDToSpec`, lines 42–62: FFT convention.
- `core/mrs.py::H2O`: positive real MRS FID scaling.
- `utils/report.py::svs_plot_refs`: embedding the reference plot in HTML.

The water model is A*exp(−gamma*t−sigma*t²−i*omega*t)*exp(i*phi), with no
baseline term. Physical frequency in Hz is **−omega/(2*pi)**. The objective
uses the full complex spectral residual over the water window 1.65–7.65 ppm.
Amplitude initialization is from log-magnitude regression; gamma and sigma
are clipped to at least 0.001 initially. Bounds are A,gamma,sigma>=0 with
unbounded omega and phase. There is no extra water baseline correction.

Default frequency/phase initialization regresses `np.angle(rotating)` directly,
without unwrapping, onto [1,t]. It then passes the positive phase slope as
omega, although the model uses **−i*omega*t**. Thus initialization contains
both phase-wrap corruption and a sign inconsistency. The optimizer may escape
this initialization (case B) but is not guaranteed to do so (case D).

After fitting, the code deliberately calls
`fid_func(t, *pout.x[:3], 0, 0) / fid_scaling` and stores that frequency/phase
corrected signal as `WaterRef.fid`. `plot_references` plots original `mrs.H2O`
alongside this corrected signal, labeled “Fitted, integrated water”. It does
not reconstruct or plot the actual fitted frequency/phase model. In particular,
the corrected integration peak is at **4.68 ppm**, while the input is at
**4.70 ppm**. This is the shift visible in the original report.

## Controlled experiments

All cases use identical amplitude 61.73993435001698, zero linewidth, unchanged
sampling and receiver/calibration. Measured phase and frequency are used for
comparison. Relative residuals below refer to the **actual fitted complex
model before frequency/phase removal**, not the integrated model.

| Case | Expected f Hz / phase rad | Default fitted f Hz | Default phase rad | Default amplitude | Relative complex residual |
|---|---|---:|---:|---:|---:|
| A | 0 / 0 | −2.31491e-9 | −2.62851e-8 | 61.7399343520 | 3.40015e-8 |
| B (original Phase 5b) | 2.55462 / 0 | 2.55461999985 | 3.25653e-10 | 61.7399343972 | 8.32619e-10 |
| C | 0 / 0.6 | 1.85391e-12 | 0.599999999994 | 61.7399343500 | 3.44935e-12 |
| D | 2.55462 / 0.6 | 0.114948308119 | 0.534580707714 | **0** | **1** |

For B, actual optimizer parameters (internally scaled amplitude, gamma, sigma,
omega, phi) are approximately (2.21034107221, 0, 0, −16.0511508485,
3.25653e-10). D's failed parameters are (0, 0.248554223063,
0.173469542220, −0.722241520661, 0.534580707714), despite `success=True`.
All exact parameters, initializations and convergence messages are in JSON.

Only after verifying sample and frequency conventions, the diagnostic replaces
the frequency initialization by **minus** the unwrapped phase slope, and phase
by the measured initial phase. Amplitude/decay initialization, bounds and TNC
settings are retained. No expected target concentration is used.

| Case | Revised fitted f Hz | Revised phase rad | Revised amplitude | Relative complex residual |
|---|---:|---:|---:|---:|
| A | −2.31491e-9 | −2.62851e-8 | 61.7399343520 | 3.40015e-8 |
| B | 2.55462000660 | −2.82237e-8 | 61.7399343664 | 1.41196e-8 |
| C | 1.85803e-12 | 0.599999999994 | 61.7399343500 | 3.45697e-12 |
| D | 2.55462000337 | 0.599999985660 | 61.7399344352 | 7.34696e-9 |

The revised initialization resolves D. It does not improve the already accurate
B fit; this matters when distinguishing the actual Phase 5b issue from the
additional optimizer weakness.

![Original-frame fit versus report integration overlay](../outputs/fsl_mrs_water_fit_diagnostics/water_fit_before_after.png)

![Case D: initialization and consistent integration frame](../outputs/fsl_mrs_water_fit_diagnostics/initialization_and_frame_case_D.png)

## Recommendations and Phase 5b implications

Recommended upstream FSL-MRS changes are to retain the uncorrected fitted
model for plotting against raw data, and show the normalized integration model
in a separate panel or transform the data into the same frame. The integrated
model's normalization is intentional and should not be removed just to align
the report's curves. Frequency initialization should unwrap phase and apply
the model's minus sign. Optimizer success alone should not certify a fit;
check complex residual, nonzero amplitude/integral, and parameter plausibility.
Phase unwrapping itself needs noise-aware handling before a general upstream
fix; these tests concern deterministic noiseless signals.

The original Phase 5b dataset is B. Its amplitude fit is accurate and its
reported concentration mismatch is **not explained by this visual overlay**.
Replacing initialization changes B's integration area by only −4.99509e-10
relative; the standard ~6.0033 mM versus synthetic 10 mM discrepancy remains
the distinct quantification convention issue documented in Phase 5b. Case D
would instead produce a zero integral and invalidate water scaling; no valid
concentration should be inferred from that fit. LCModel validation remains
deferred.

## Reproduction and outputs

```
pixi run python -m examples.fsl_mrs_water_fit_diagnostics outputs/fsl_mrs_water_fit_diagnostics
pixi run python -m pytest -q tests/integration/test_fsl_mrs_water_fit_diagnostics.py tests/integration/test_fsl_mrs_water_scaling.py tests/test_mrscloud_quantitative_pair.py
```

The diagnostic script uses the actual saved Phase 5b 10 mM dataset and writes
the full trace, four case arrays, exact measured/fitted parameters, and two
standalone figures. Its integration tests generate equivalent inputs in a
temporary directory rather than requiring previously generated output files.
The default case-D defect is preserved as a specific expected failure, with
successful recovery under sample-derived initialization checked independently.
Actual combined result: **23 passed, 1 xfailed**, zero unexpected failures.
Two warnings concerned the installed pytest plugin deprecation and Octave's
forkpty use in a multithreaded process. The diagnostic's expected failure is
the installed default case-D recovery, not a waived signal-convention or
export error. No LCModel work was performed.

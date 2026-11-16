# Phase 7b.1: LCModel concentration and residual diagnosis

The original 10.639 mM NAA result contains two distinct effects: a
**normalization-only prediction of approximately 9.526 mM**, followed by a
**model-dependent increase of approximately 1.113 mM** in the nonlinear fit.
Water-reference arithmetic does not explain that upward fitting contribution.
There is a real acquisition-duration mismatch between LCModel's retained basis
and the finite synthetic acquisition, and its baseline/broadening/line-shape
trade-off changes the reported concentrations. These effects are measured below;
they are not reduced to a unique, fitted correction factor.

No simulator, calibration, BASIS or input FID was changed. The Phase 7b baseline
is immutable. All new executions use the same installed **LCModel 6.3-1N** binary
`$LCMODEL_BINARY` and exactly the same RAW/H2O
and LCModel-compatible BASIS files as the original fit.

## Immutable baseline and reproducibility

Original files remain in `outputs/lcmodel_validation/noise_free/`, including
`lcmodel.CONTROL`, `metabolites.RAW`, `water_reference.H2O`,
`lcmodel_compatible.BASIS`, and `lcmodel_outputs/fit.{COORD,PRINT,TABLE,PS}`.
The original Phase 6a BASIS remains in its calibration bundle.

`outputs/lcmodel_fit_bias_diagnostics/immutable_baseline_sha256.json` records
absolute paths and SHA256 hashes of all original fit files, including the
original failed-BASIS attempt, source acquisition files, calibration bundle and
executable. The runner checks this manifest before accepting cached results and
again after execution. The baseline normalization, concentration results,
residual, line-shape parameters and diagnostics are recorded in
`baseline_metrics.json`. A fresh **replay** fit, with only output paths changed,
reproduces every baseline COORD array exactly and all TABLE concentrations.

Original NAA10.639 mM; FCALIB121933; phased-real residual0.018163864406490208;
water area0.000933848; normalized NAA basis area0.00145922. Original diagnostics
MYBASI8, STARTV20, FINOUT7 and FINOUT9 remain preserved.

The analysis uses the supplied `$LCMODEL_MANUAL`
and local same-version source `$LCMODEL_SOURCE`.
Their hashes, and the Pixi lock hash, are recorded in `independent_analysis.json`.
No package, binary or local installation was modified. All Python ran in Pixi.

## 1. Independent normalization reconstruction

Manual §§10.2.1–10.2.2 (pp128–130) and source `WATER_SCALE`, lines4761–4802,
implement

\[
 A_{b,norm}={A_{NAA}\over N_{1H,NAA}\,ATTMET},\qquad
 A_{w,norm}={A_w\over2\,WCONC\,ATTH2O},\qquad
 F={A_{b,norm}\over A_{w,norm}}.
\]

Here the serialized basis is already scaled by
`TRAMP/(VOLUME*CONC)` in MYBASI, line3675. For every input component,
CONC=TRAMP=VOLUME=1. RAW and H2O also have TRAMP=VOLUME=1. WCONC is
**39016.5 mM**, matching the actual WM acquisition, not the earlier 40 M example.
ATTH2O=ATTMET=1, and N1HMET=3 for NAA's methyl reference. The water divisor2
belongs to LCModel quantification, not an additional multiplier on sysH2O.

AREABA integrates the NAA reference before basis bandwidth conversion, using
symmetric local-minimum bounds and bordering-window baseline subtraction.
The search is centered at WSPPM2.01, with reference bounds depending on
`RFWBAS*FWHMBA/2`; default RFWBAS10. AREAW2 integrates the water on its own
2048-point spectrum. Its transform is divided by sqrt2 to match the zero-filled
4096-point data normalization; its bordering-window count uses the **data** ppm
step, not the water-spectrum step. No spectral peak height is substituted for
an integrated area.

An independent NumPy implementation reconstructs these operations from the
unchanged serialized samples. It also emulates CSFT_R's sequential
single-precision complex-exponential recurrence. Known zero water phase is
used; GETPHA is not generally reimplemented. Source locations:
AREAW2 lines4885–4980, GETPHA immediately after, AREABA lines5174–5330,
CSFT_R lines12688–12721.

| Quantity | Independent calculation | Actual original PRINT |
|---|---:|---:|
| NAA normalized reference area | 0.0014592148729401777 | 0.00145922 |
| Water area | 0.0009338442469015718 | 0.000933848 |
| FCALIB | **121933.51788367616** | **121933** |

FCALIB agrees to **4.25e-6 relative**; the remaining difference is consistent
with single-precision complex exponentials/accumulation, known-zero versus
computed phase and rounded five-significant-digit PRINT output. The original
Phase 7b double-precision FFT estimate121935.4658 was also close, but omitted
the installed direct-transform rounding. No empirical multiplier is used.

The native BASIS has 8192 points at4000Hz. LCModel crops to4096 frequency
points at2000Hz, preserving native spectral spacing, then applies its unitary
inverse FFT. The native spectra use an unnormalized FFT. The derived common
basis-time scale is therefore

\[
 K={8192\over\sqrt{4096}}=128.
\]

This is an interpretation of reader transforms, not a modification of either
signal. For nominal concentration C in mM, the expected unscaled coefficient
is `C*1e-3/K`. Multiplication by the printed F gives the normalization-only
prediction. For NAA:

\[
 10\times{121933\over128\times1000}=9.526015625\ {\rm mM}.
\]

The independent F gives9.5260561mM; the earlier double-precision calculation
9.5262083mM differs only in the last few0.0001mM. These normalization errors
are far smaller than the final +0.639mM bias.

The arithmetic decomposition using printed F is
`10 - 0.473984375 + 1.112984375 = 10.639 mM`.
The last term is the observed nonlinear-model contribution relative to a
matching fixed-basis prediction, **not a transferable calibration correction**.

## 2. Basis/data compatibility

LCModel source MYBASI lines3662–3880 and FTDATA lines6038–6058 show:

1. The Phase 7b input adapter preserves all native spectral numbers exactly.
2. ISHIFT8 cancels the implicit `nint((4.65-4.68)/native_ppm_step)=-8` origin
   compensation. Net shift is zero; no frequency or phase correction is fitted
   independently for each metabolite in this audit.
3. LCModel copies the first2048 and last2048 native spectral points into a
   4096-point array, then uses a unitary inverse transform.
4. It retains this complete4096-sample basis time response.
5. Data contain2048 acquired samples and are padded with2048 zeros by FTDATA.

Thus matching bandwidth and field strength do **not** imply matching effective
acquisition duration. The simulator used FSL formatting to produce2048 samples
from the native8192-point basis. The model used by LCModel still contains basis
signal after the acquisition ends. This includes late response and effects of
periodic Fourier formatting; it is not uniquely a physical exponential tail.

After the single analytically derived factor128, acquired-sample differences
and unacquired-tail norms are:

| Component | Acquired full-complex relative difference | Max absolute difference | Retained tail/acquired norm |
|---|---:|---:|---:|
| NAA | 6.50696e-4 | 8.02130e-8 | 4.14800% |
| Cr | 5.07908e-4 | 7.45424e-8 | 4.11553% |
| GPC | 6.24498e-4 | 2.15515e-7 | 2.50697% |
| PCh | 4.23679e-4 | 1.52470e-7 | 4.06939% |
| Glu | 1.34542e-3 | 5.78867e-8 | 2.92875% |

The discrepancy persists when sample zero is excluded; it is not simply an
incorrect first-point factor. In this non-magnitude mode ABSVAL=False, LCModel's
special half-first-point operations in DATAIN/MYBASI are inactive. No factor-of-two
repair is justified. Native basis DC processing remains exactly as supplied;
this phase does not change metabolite or water DC correction.

### Independent projections

Expected coefficients come from the original `ground_truth.json`, not a fitting
result. The matching FSL-formatted basis, LCModel bandwidth-converted basis,
and unchanged acquired FID are compared without nuisance optimization.
Complex projection allows complex coefficients (their imaginary parts are
recorded); real projection uses real coefficients. It is an unconstrained
diagnostic, not LCModel's nonnegative regularized estimator.

| Design | Relative projection residual | NAA source-equivalent mM | Normalization-scaled NAA mM |
|---|---:|---:|---:|
| Matching FSL-formatted acquired complex FID | 2.388e-15 | 10.000000000 | — |
| LCModel acquired-sample complex FID | 0.00118893 | 9.999989663 | 9.526005778 |
| LCModel full-duration spectrum, complex | 0.03631848 | 9.982794678 | 9.509625808 |
| LCModel full-duration spectrum, real | 0.03464625 | 9.891878350 | 9.423018772 |
| Same real design + cubic baseline | 0.03435379 | 9.884618402 | 9.416102935 |
| LCModel diagnostic duration-matched real design | **8.88915e-6** | 10.000007175 | 9.526022460 |

The last row truncates the **independent analysis design only** after2048 time
samples; it does not alter the BASIS, RAW or any LCModel fit input. This
intervention isolates a substantial representation mismatch. A cubic baseline
alone does little to absorb it. These residuals and coefficient shifts are
conditional on fixed phase/frequency/broadening and cannot be equated to the
nonlinear LCModel residual or used to assign additive fractions of its bias.

## 3. Controlled LCModel executions

Every case references the original FILRAW, FILH2O and FILBAS paths; only its
output paths and listed documented options change. Each directory retains its
CONTROL, stdout/stderr/exit status, PRINT, TABLE, COORD, PS, PDF and metrics.
One baseline replay plus seven diagnostic controls were executed. Seven
produced finite final analyses; one produced the retained fatal diagnostic.

| Case | Changed options | Manual reference |
|---|---|---|
| replay | None | Immutable baseline replay |
| preliminary_matched | NUSE1=5; CHUSE1=NAA,Cr,GPC,Glu,PCh | §11.5 p140 |
| half_broadening_prior | DEEXT2=1, half the default2s^-1 | §11.13 pp148–149 |
| wider_convolution | RFWHM=3, from default1.8 | §11.12 p148 |
| phantom_mode | VITRO=True | §5.3.8 p40 |
| known_phase_prior | DEGZER=DEGPPM=0; SDDEGZ=3, SDDEGP=1 | §9.7 pp115–116 |
| default_basis_width | FWHMBA=.013 rather than1/127.731ppm | §8.6.4.1 p92 |
| zero_broadening_prior | DEEXT2=0 | Boundary control; SETUP14 rejects it |

No hidden regularization, minimum-support or initial-linewidth parameters were
changed. The manual warns against modifying DKNTMN, ALPB*, ALPS*, FWHMST and
NSIDMN. VITRO is the public baseline/statistical mode: source1083–1105 switches
knot spacing from0.15 to0.6ppm and rescales baseline penalties; source7459 also
changes a regularization-search branch. It is **not a pure single-parameter
baseline intervention**. Results must not be attributed solely to knot spacing.

DEEXT2=0 was rejected with SETUP14 despite process exit0. The source check at
6409 requires a strictly positive starting broadening. Manual SETUP14 explicitly
describes this boundary. That run is preserved as failed, not counted as a valid
concentration estimate. DEEXT2=1 is a separate, documented positive control,
not a silently altered replacement for the failed run.

### Actual concentration and residual results

All concentrations are reported mM. Baseline norm is `||baseline||/||data||`
in the phased-real output window, not a concentration or a signal correction.

| Case | NAA | Cr | GPC | PCh | Glu | Real residual | Baseline norm | FCALIB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| replay | 10.639 | 8.446 | 2.125 | 1.051 | 8.376 | 1.816386% | 4.051746% | 121933 |
| preliminary_matched | 10.639 | 8.446 | 2.125 | 1.051 | 8.376 | 1.816319% | 4.052361% | 121933 |
| half_broadening_prior | 10.424 | 8.262 | 2.076 | 1.040 | 8.222 | 1.917266% | 3.228370% | 121933 |
| wider_convolution | 10.639 | 8.446 | 2.125 | 1.051 | 8.376 | 1.816386% | 4.051746% | 121933 |
| phantom_mode | 9.959 | 7.890 | 1.974 | 1.008 | 7.795 | 2.233648% | 1.148688% | 121933 |
| known_phase_prior | 10.639 | 8.446 | 2.125 | 1.051 | 8.376 | 1.816335% | 4.051465% | 121933 |
| default_basis_width | 10.952 | 8.695 | 2.188 | 1.082 | 8.623 | 1.816390% | 4.051627% | 125523 |
| zero_broadening_prior | — | — | — | — | — | Fatal SETUP14 | — | No final fit |

Matching preliminary names removes MYBASI8 but does not correct the bias.
Phase priors and larger convolution support also leave concentrations unchanged
at reported precision. Increasing RFWHM does not change the final minimum
support of5 side points here. NSIDMN5 limits **support**, not the minimum
physical linewidth; a delta-like kernel is still representable.

The reduced broadening prior changes the NAA result from10.639 to10.424mM,
with more residual and less baseline. Phantom mode further reduces baseline and
returns NAA9.959mM, again with more residual. Closer agreement with nominal10
is **not proof** that either configuration has resolved all response conventions.

### Line shape and optimizer diagnostics

| Case | Baseline knots | alphaB | alphaS | NAA additional 1/T2, s^-1 | Report FWHM, ppm |
|---|---:|---:|---:|---:|---:|
| replay | 30 | .014 | .2 | 3.47 | .004 |
| preliminary_matched | 30 | .014 | .2 | 3.47 | .004 |
| half_broadening_prior | 30 | .031 | .2 | 2.69 | .004 |
| wider_convolution | 30 | .014 | .2 | 3.47 | .004 |
| phantom_mode | 10 | 6.3 | .4 | 1.46 | .004 |
| known_phase_prior | 30 | .014 | .2 | 3.47 | .004 |
| default_basis_width | 30 | .014 | .2 | 3.47 | .011 |

Metrics retain every component's shift, broadening, prior deviation in SD units,
and last printed off-center line-shape coefficients. PRINT's signed `SDs` fields
are deviations from prior expectations, **not estimated parameter standard
errors**. The convolution center is constrained to1 minus the sum of the
side coefficients (source9480 onward); large alternating side coefficients
compensate for added broadening and representation differences. Broadening and
convolution cannot be interpreted independently from their final combination.

Default fits retain FINOUT7 (flexible baseline) and FINOUT9 (unsmooth line shape).
Phantom mode retains MYBASI8 and STARTV20 but no FINOUT7/9. The half-prior case
adds RFALSI13 (numerical instability), RFALSI5 (step/gradient stopping issue)
and TWOREG7 (a refinement analysis failed to converge properly). It produced a
finite final result, but is not certified as fully converged. LCModel supplies
no generic optimizer-success boolean in these outputs. All process exits were0,
including the failed zero-prior run. Displayed FWHM is not an independent
measurement of the physical peak width, especially with FINOUT9 present.

## 4. Reference-integration sensitivity

Changing only FWHMBA to the manual default.013 changes the **NAA integration
window**: independently reconstructed normalized basis area0.0015021776890031043,
FCALIB125523.53563766646; actual PRINT gives0.00150218 and125523. Water area
is unchanged. The NAA concentration ratio10.952/10.639 closely matches
125523/121933, while residual and baseline are nearly identical. Thus this
control identifies concentration **normalization**, separate from fitting.
FWHMBA is not merely a display label.

The manual also recommends4times physical FWHM for Lorentzian basis singlets.
For this1Hz basis, that would be.0313158ppm, outside its stated acceptable range
(<.025ppm). No unreviewed out-of-range setting was tried. The effective-width
and integration policy for these unusually narrow, Lorentzian, finite-duration
synthetic responses remains a scientific compatibility question, not a reason
to change the water amplitude.

## 5. Inspection outputs

Root: `outputs/lcmodel_fit_bias_diagnostics/`. Actual PS/PDF reports:

- [Replay](../outputs/lcmodel_fit_bias_diagnostics/replay/lcmodel_outputs/fit.pdf)
- [Matched preliminary set](../outputs/lcmodel_fit_bias_diagnostics/preliminary_matched/lcmodel_outputs/fit.pdf)
- [Half broadening prior](../outputs/lcmodel_fit_bias_diagnostics/half_broadening_prior/lcmodel_outputs/fit.pdf)
- [Wider convolution](../outputs/lcmodel_fit_bias_diagnostics/wider_convolution/lcmodel_outputs/fit.pdf)
- [Phantom mode](../outputs/lcmodel_fit_bias_diagnostics/phantom_mode/lcmodel_outputs/fit.pdf)
- [Known phase prior](../outputs/lcmodel_fit_bias_diagnostics/known_phase_prior/lcmodel_outputs/fit.pdf)
- [Default basis width](../outputs/lcmodel_fit_bias_diagnostics/default_basis_width/lcmodel_outputs/fit.pdf)
- [Zero-prior fatal report](../outputs/lcmodel_fit_bias_diagnostics/zero_broadening_prior/lcmodel_outputs/fit.pdf)

`comparison.csv`, `baseline_metrics.json`, `independent_analysis.json` and
`results.json` contain unrounded numerical results; each case's `metrics.json`
contains its diagnostics and configuration. `commands.txt` records every exact
executable invocation and case CONTROL retains all exact options.

![Input, model, residual and baseline](../outputs/lcmodel_fit_bias_diagnostics/spectral_comparisons.png)

![Concentration error versus fitting residual](../outputs/lcmodel_fit_bias_diagnostics/concentration_vs_residual.png)

Plots use the original LCModel output scale without quantitative display
renormalization. COORD contains phased-real spectra; no full-complex LCModel
model residual is claimed. Independent complex projections are labeled separately.

## 6. Execution and tests

```sh
pixi run python --version
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.lcmodel_fit_bias_diagnostics
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_lcmodel_fit_bias_diagnostics.py tests/integration/test_lcmodel_bias_diagnostics.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q
```

The script invokes the local executable with CONTROL on stdin, a case-specific
working directory, captured stdout/stderr and a300-second timeout, then converts
PS to PDF using installed Ghostscript. Completed cases are reused only when
their CONTROL and immutable input hashes agree; no binary download is needed.

Focused tests: **10 passed**. They cover independent geometric-series transform
validation, duration-mismatch controls, coefficient projection, normalization,
actual replay equality, shared input paths, report existence, warning/fatal
retention, immutable source hashes and restart without new external executions.
An initial transform-test bound of3e-5 failed at3.23e-5 on the256-point
single-precision recurrence. The final analytical rounding bound is4*N*epsilon32,
reflecting accumulated complex multiply/add operations; the signal, fixture and
normalization checks were unchanged. Normalization reconstruction independently
retains its6e-6 relative tolerance against actual PRINT values.

Final full regression: **298 passed, 1 xfailed**, no failures or skips,
in **195.14 seconds**. The pre-existing expected failure is the installed
FSL-MRS WaterRef combined phase/frequency initialization defect. Four warnings
are external Dash/SWIG deprecations. The initial focused run was9passed/1failed
on the overly tight float32 transform bound described above; the corrected
focused run was10passed. The first script invocation stopped when it encountered
the retained SETUP14 failure; subsequent invocations reused the completed fits
and continued through the remaining controls without altering that failure.

Files added: `examples/lcmodel_fit_bias_diagnostics.py`,
`tests/test_lcmodel_fit_bias_diagnostics.py`,
`tests/integration/test_lcmodel_bias_diagnostics.py`, this document, and dedicated
output artifacts. No existing simulator, Phase7b utility or source data changed.

## Interpretation and remaining uncertainty

The mismatch is not a failed RAW conversion, a water-frequency error, or an
unexplained water gain. It combines a documented non-unit reference-integration
scale with an acquisition-duration/basis representation mismatch and nonlinear
regularized baseline/line-shape compensation. The controlled fits show genuine
model dependence: changing baseline/statistical mode or broadening priors changes
concentrations while water scaling stays fixed. FWHMBA independently changes
reference integration without materially changing the residual.

The retained basis tail is isolated as a substantial source of mismatch in an
independent fixed-nuisance design. These experiments do not uniquely decompose
the nonlinear +1.113mM contribution into baseline, convolution, broadening and
finite-duration terms. Rounded PRINT parameters and the standard COORD output
are insufficient to reconstruct all internal coefficients uniquely. Hidden
regularizers were deliberately left untouched. More direct internal-design
capture or a separately reviewed acquisition-matched basis representation would
be needed for a unique decomposition. No correction was implemented, and no
noise Monte Carlo fitting or physical calibration claim follows from this phase.

# Phase 5j.1: observed water and NAA linewidths

Raw sequence-matched water and isolated three-spin NAA methyl have **equal
absorption-mode FWHM to numerical precision in all twelve conditions**.
Their observed FWHM is not always the nominal decay linewidth: the shortest
1 Hz record measures **1.066865 Hz** for both. Processed methyl's zero-level
half-height width differs from water by at most **.021887 Hz**; full NAA's
acetyl peak differs by at most **.011167 Hz**. DC background and Fourier
formatting explain these small apparent changes; they do not establish
different intrinsic water/metabolite relaxation times.

This is a read-only analysis of Phase 5j. No signal was regenerated or refitted,
and no calibration, fitting pipeline, simulator or existing output was changed.
FSL-MRS version is `2.5.1.dev78+ga311e38f2`. Source-file SHA256 hashes, including
the saved reports and fit metadata, were verified unchanged after analysis.

## Inputs and reference conventions

Inputs are `outputs/mrscloud_sampling_sensitivity/validated/master_responses.npz`,
the twelve simulation configurations, matching BASIS files and saved fitted
parameters. Raw water, raw isolated methyl, DC-only methyl/full NAA and
formatted methyl/full NAA all retain their original common normalization.
Each condition uses the saved acquisition prefix and the original common
additional Lorentzian factor exp[-pi*(L-1)*t], with dt=.0005 s, N=2048/4096/8192,
nominal L=1/3/5/8 Hz. Quantitative integrals retain the original shared scale.

Spectral origin is 4.68 ppm and metadata frequency is 127.731 MHz. The water
center is +2.55462 Hz. The actual methyl center follows the original spin and
Hamiltonian/export conventions:

\[
f_m=(2.008-3)127.728+(3-4.68)127.731=-341.294256\ \mathrm{Hz}.
\]

Before analysis, temporary copies are demodulated to their own target resonance
and phased using the **raw isolated methyl/water phase at t=0**, both effectively
zero. These phases are known from the simulated single resonances. Processed
sample-zero phase is not used: DC and formatting change it without rotating
the acetyl resonance by the same phase. No phase is optimized to force widths
to agree. Full NAA is examined only within +/-20 Hz of the acetyl resonance;
the entire multi-resonance NAA spectrum is never treated as one peak.

For fitted water, stored gamma/sigma reconstruct its Voigt shape after removing
the known fitted frequency/phase. A unit model amplitude is sufficient because
FWHM is amplitude-invariant; no acquisition amplitude or fit coefficient changes.
For fitted NAA, the saved matching source basis and gamma/eps reconstruct the
reference model with the fitted Phi0/Phi1 removed, as FSL's no-phase reference
model does. The fitted center includes -eps/(2*pi). Intrinsic basis phase and
other NAA resonances remain. Fitting is not rerun.

## Width estimator and baseline conventions

The estimator half-weights sample zero **before** zero-filling, evaluates a
shifted FFT on **128*N points**, and takes the real, phase-corrected spectrum.
It finds the local positive central maximum and the first half-height crossing
on either side, with linear interpolation between spectral samples. The
interpolated grid spacing is .0076294, .0038147 or .00190735 Hz depending on N.
Doubling zero-filling changes processed-methyl FWHM by at most
**8.2804e-6 Hz**. Zero-filling evaluates the existing finite-record transform;
it does not improve acquired resolution or extend its decay.

The primary table reports **zero-level absorption half-height widths**: half
the peak height relative to zero. This is an unambiguous representation-level
measurement but requires care when DC processing introduces a pedestal.
It must not automatically be interpreted as background-corrected intrinsic
linewidth. Background corrections below are reported separately, rather than
quietly choosing a baseline that makes the widths equal.

No magnitude-spectrum widths are mixed into the table. Unit tests explicitly
distinguish absorption and magnitude Lorentzian widths (the latter approaches
sqrt(3) times the absorption width for a long record).

## Measured widths for every condition

All entries are Hz. “Full acetyl” means the local acetyl peak in the complete
formatted NAA response, including its local tails/background. “Fit NAA” is the
same local peak in the saved fitted NAA model.

| Points | Nominal | Raw water | Raw methyl | Processed methyl | Full acetyl | Fit water | Fit NAA |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2048 | 1 | 1.066865 | 1.066865 | 1.066563 | 1.066456 | 1.066865 | 1.066456 |
| 2048 | 3 | 3.000354 | 3.000354 | 2.997279 | 2.997464 | 3.000354 | 2.997464 |
| 2048 | 5 | 5.000032 | 5.000032 | 4.991481 | 4.993443 | 5.000032 | 4.993442 |
| 2048 | 8 | 8.000108 | 8.000108 | 7.978221 | 7.988941 | 8.000626 | 7.988943 |
| 4096 | 1 | 1.000265 | 1.000265 | .999917 | .999881 | 1.000265 | .999882 |
| 4096 | 3 | 3.000007 | 3.000007 | 2.996930 | 2.997116 | 3.000007 | 2.997117 |
| 4096 | 5 | 5.000027 | 5.000027 | 4.991479 | 4.993439 | 5.000026 | 4.993436 |
| 4096 | 8 | 8.000106 | 8.000106 | 7.978220 | 7.988941 | 8.000609 | 7.988941 |
| 8192 | 1 | 1.000002 | 1.000002 | .999640 | .999576 | 1.000002 | .999579 |
| 8192 | 3 | 3.000006 | 3.000006 | 2.996929 | 2.997115 | 3.000006 | 2.997113 |
| 8192 | 5 | 5.000026 | 5.000026 | 4.991478 | 4.993438 | 5.000028 | 4.993433 |
| 8192 | 8 | 8.000105 | 8.000105 | 7.978220 | 7.988941 | 8.000781 | 7.988933 |

The raw water/methyl differences are at roundoff, below 1e-11 Hz. Nominal
versus measured differences are largest for 2048/1 Hz (+.066865 Hz), falling
to +.000265 Hz at 4096 and approximately +.000002 Hz at 8192. Wider lines
have much smaller finite-record broadening; small residual differences include
the discrete-time Lorentzian and crossing interpolation.

| Points | Nominal | Processed methyl minus water | Full acetyl minus water |
|---:|---:|---:|---:|
| 2048 | 1 | -.000302 | -.000409 |
| 2048 | 3 | -.003074 | -.002889 |
| 2048 | 5 | -.008550 | -.006589 |
| 2048 | 8 | -.021887 | -.011167 |
| 4096 | 1 | -.000348 | -.000384 |
| 4096 | 3 | -.003077 | -.002891 |
| 4096 | 5 | -.008548 | -.006588 |
| 4096 | 8 | -.021886 | -.011165 |
| 8192 | 1 | -.000362 | -.000426 |
| 8192 | 3 | -.003077 | -.002891 |
| 8192 | 5 | -.008548 | -.006588 |
| 8192 | 8 | -.021886 | -.011165 |

## Isolating DC, formatting and overlap effects

The DC-only master differs from raw methyl only at sample zero. Under the
same half-first-point FFT, that difference is a known complex constant spectral
background. Subtracting its correctly phased real pedestal before calculating
half-height gives the **raw methyl FWHM exactly** (tested to 1e-11 Hz in all
conditions). DC subtraction itself therefore does not change the underlying
exponential decay; it changes the zero-level half-height measurement.

For the formatted response, the same known native sample-zero DC delta is
propagated through the original Fourier formatter. Its resulting background
is subtracted only in a separate diagnostic spectrum. This uses linearity and
the measured processing delta, not a fitted width or empirical multiplier.
The saved FIDs remain untouched. Some format-related peak-shape difference
persists after removing this DC contribution:

| Points | Nominal | Raw methyl | Processed methyl, DC background removed |
|---:|---:|---:|---:|
| 2048 | 1 | 1.066865 | 1.066737 |
| 2048 | 3 | 3.000354 | 2.999045 |
| 2048 | 5 | 5.000032 | 4.996389 |
| 2048 | 8 | 8.000108 | 7.990770 |
| 4096 | 1 | 1.000265 | 1.000117 |
| 4096 | 3 | 3.000007 | 2.998696 |
| 4096 | 5 | 5.000027 | 4.996383 |
| 4096 | 8 | 8.000106 | 7.990767 |
| 8192 | 1 | 1.000002 | .999836 |
| 8192 | 3 | 3.000006 | 2.998696 |
| 8192 | 5 | 5.000026 | 4.996383 |
| 8192 | 8 | 8.000105 | 7.990767 |

Thus DC explains part of the apparent narrowing, while Fourier filtering,
causal/periodic boundary effects and finite prefixes explain a residual shape
change. These are representation effects, not evidence that methyl has a
different native T2 from water. Full NAA's other resonances contribute local
tails and complex background: for example, at 8192/8 Hz raw full acetyl measures
8.030458 Hz versus isolated methyl 8.000105 Hz. The full-versus-isolated control
reveals overlap/background sensitivity without treating other resonances as
part of the methyl linewidth parameter.

Frequency-bin alignment is removed as a linewidth-comparison confound by
demodulating to each true center before interpolating. That does **not** alter
the original coarse-bin integrals used in Phase 5j quantification.

## Fitted parameters are not observed FWHM

The water model is exp(-gamma*t-sigma*t^2); gamma/pi is only its Lorentzian
parameter. NAA gamma/pi is **additional** broadening of its supplied processed
basis, not the complete acetyl peak width. Saved values are:

| Points | Nominal | Water gamma/pi Hz | Water sigma s^-2 | NAA additional gamma/pi Hz |
|---:|---:|---:|---:|---:|
| 2048 | 1 | .999999937 | 0 | .000000229 |
| 2048 | 3 | 3.000000003 | 0 | 1.999999276 |
| 2048 | 5 | 4.999998134 | .000045856 | 3.999998724 |
| 2048 | 8 | 7.998952587 | .041183260 | 7.000001263 |
| 4096 | 1 | 1.000000054 | 0 | .000001016 |
| 4096 | 3 | 2.999999950 | 0 | 2.000000594 |
| 4096 | 5 | 4.999999446 | 0 | 3.999997261 |
| 4096 | 8 | 7.998981276 | .040064553 | 7.000000201 |
| 8192 | 1 | 1.000000042 | 0 | .000002230 |
| 8192 | 3 | 3.000000037 | 0 | 1.999997550 |
| 8192 | 5 | 5.000001542 | 0 | 3.999995296 |
| 8192 | 8 | 7.998634582 | .053711330 | 6.999992488 |

The fitted NAA local width differs from generated full acetyl by less than
2e-5 Hz. Water's 8 Hz Voigt-model width differs from its generated peak by less
than .001 Hz, consistent with the previously reported finite optimizer error.
At 2048/1 Hz the fitted gamma/pi is approximately 1 Hz while both generated
and fitted finite-record water peaks have observed FWHM approximately 1.067 Hz.
This explicitly distinguishes parameter recovery from observed peak width.

## Quantification implications

Equal raw absorption FWHMs do **not** imply equal FSL reference integrals.
At 2048/1 Hz, raw water and methyl both measure 1.066865 Hz, but the independently
defined ideal 3H/2H response factor remains **.9098704421**. Original frequencies,
FFT-bin alignment, finite record, half-first-point weighting, abs-real operation
and unequal fixed windows matter to the integral ratio independently of FWHM.

The script also calculates a conditional width-only sensitivity using the
unchanged Phase 5j integration functional. It holds analytical control initial
amplitudes, water linewidth, frequencies and windows fixed, and changes only
the methyl exponential decay by the measured methyl-minus-water width delta.
This is a **hypothetical Lorentzian surrogate**, not a correction to data:
processed peak widths include pedestal/filtering effects and cannot generally
be identified with a pure decay parameter.

- Using the raw methyl/water width delta changes predicted concentration by
  less than 1e-10 mM (numerically zero).
- Using the apparent processed-methyl width delta changes that surrogate's
  prediction by at most **.002696 mM** across the twelve conditions.
- Actual Phase 5j concentration deficits are approximately .22–1.02 mM.

The surrogate therefore does not explain the measured bias under its declared
Lorentzian assumption. It is not a universal bound or a causal allocation.
The direct equal-width counterexample establishes that most of this bias can
exist without differential raw linewidth. Known DC backgrounds and other
peak-shape/integration effects change areas without a proportional FWHM change;
their effects must be examined through full complex responses and the actual
integration functional. No concentration estimate or calibration was changed.

## Plots and reproducibility

All overlaid spectra are in absorption mode and aligned to their own peak
frequency. Peak-height normalization is **only for display**. It is never used
for FWHM calculations, processing-background estimates or quantitative integrals.

![Short narrow-line acquisition](../outputs/mrscloud_measured_linewidth_validation/final/N2048_LW1_absorption.png)

![Short broad-line acquisition](../outputs/mrscloud_measured_linewidth_validation/final/N2048_LW8_absorption.png)

![Long narrow-line acquisition](../outputs/mrscloud_measured_linewidth_validation/final/N8192_LW1_absorption.png)

Each figure compares water/raw methyl, raw/DC/formatted methyl including the
DC-background control, and saved fitted/generated shapes. No new HTML report
or fitting run is required or implied in this phase.

Files added: `examples/mrscloud_measured_linewidth_validation.py`,
`tests/test_measured_linewidth_validation.py`,
`tests/integration/test_measured_linewidth_analysis.py`, and this document.
Outputs are in `outputs/mrscloud_measured_linewidth_validation/final`:
`linewidths.csv`, `results.json` (all crossing positions, fit parameters,
source hashes and baseline controls), and three PNG figures.

Executed commands:

```bash
pixi run python --version
pixi run python -m examples.mrscloud_measured_linewidth_validation outputs/mrscloud_measured_linewidth_validation/final
pixi run python -m pytest -q tests/test_measured_linewidth_validation.py tests/integration/test_measured_linewidth_analysis.py
pixi run python -m pytest -q tests
```

Use a fresh output directory when reproducing. Focused tests initially passed
**10 tests**. They check analytical Lorentzian half-height crossings, known
backgrounds, phase/frequency invariance, finite-record widths, absorption versus
magnitude conventions, interpolation precision, and source immutability.
The integration fixture replaces fitting and final quantification functions
with exceptions, demonstrating that the analysis uses saved results only.
The full regression run completed with **135 passed, 1 xfailed, 5 dependency
warnings in 87.46 s**. The expected failure remains the documented installed
FSL-MRS combined phase/frequency initialization defect. All existing analytical,
voxel, frozen-basis, calibration, export, processing, CLI and quantification
regression suites executed successfully.

## Limits

Observed width depends on truncation, phase, baseline convention, neighboring
resonances and filtering. The primary zero-level width and explicitly
background-corrected controls must not be conflated. Neither is a new physical
calibration. Nonparametric interpolation cannot uniquely attribute all effects
in the complete NAA peak; the isolated subsystem and processing controls provide
the necessary context. The original PRESS ideal-excitation/localization and
gyromagnetic conventions remain, and no physical scanner or LCModel validation
is claimed. Review this measurement convention before interpreting fitted
linewidth parameters as measured water/metabolite FWHM.

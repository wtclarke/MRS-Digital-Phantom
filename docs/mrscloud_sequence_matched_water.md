# Phase 5e: sequence-matched water feasibility

## Decision

**The bundled MRSCloud/Octave PRESS simulator can generate water with its RF
carrier at 4.7 ppm and a matching NAA acquisition at 3.0 ppm.** Both responses
were simulated, exported on one shared scale and fitted with the installed
FSL-MRS CLI. Proton-response scaling and deterministic water generation pass
numerical tests.

**Do not replace the analytical reference yet.** The preserved MRSCloud DC
processing changes the first FID sample, and subsequent Fourier resampling
introduces early-time ringing. The sequence-water FID therefore has a 2.85%
complex residual against FSL's single-Voigt water model and a fitted amplitude
nearly twice its processed initial sample. This is a processing-convention
issue requiring review before treating either quantity as a physical water
response. No amplitudes were tuned to make the concentration equal 10 mM.

## Code and sequence inspected

Relevant files under `simulation/basissets/MRSCloud-main/functions/`:

- `load_parameters.m`: waveform loading, spatial grid/gradients, timing,
  spin-system shifts, Hamiltonian and refocusing propagators.
- `sim_signals.m`: actual unedited PRESS evolution, coherence filters,
  spatial sums, ADC readout, spatial amplitude factors and DC correction.
- `FID-A/simulationTools/sim_lcmrawbasis.m`: an alternative wrapper selecting
  `sim_press`, the ideal-pulse PRESS implementation. It is **not** used here
  because it would omit the shaped refocusing of the MRSCloud path.
- `FID-A/simulationTools/metabolites/spinSystems.mat`: sysH2O is one spin,
  shift 4.7 ppm, scaleFactor 2; sysRef0ppm is one spin, scaleFactor 1;
  sysNAA includes its explicit spin/proton response.
- `experim_funct/sim_Hamiltonian_mgs.m`: spin Hamiltonian and scaleFactor
  applied to the initial density matrix.
- `FID-A/simulationTools/sim_readout.m`: density-matrix readout, subsystem
  spin-dimension factor, Lorentzian acquisition linewidth, receiver phase.
- `osprey_functions/fit_makeBasis_MRSCloud.m`: second DC correction and
  one common maximum-real-spectrum normalization.
- `osprey_functions/io_writelcmBASIS.m`: acquisition reference shift,
  conjugation, FFT and LCModel spectral serialization.

The wrapper `run_simulations_cloud.m` reads a single centreFreq for an entire
run and deletes its output folder contents. The standalone experiment bypasses
that wrapper and calls the actual parameter loader and sequence function in
isolated folders. No source JSON, pulse file or simulator function is edited.
The repository's configured Pixi Octave executable was used directly rather
than mutating MatlabRunner configuration.

## RF carrier versus acquisition reference

`load_parameters.m` subtracts `centreFreq` from each spin shift **before**
building the Hamiltonian. Shaped refocusing propagators are calculated with
that Hamiltonian, so changing centreFreq changes RF off-resonance and spatial
localization response. It is not an output frequency shift.

The high-level options do not expose two independent RF/ADC carrier inputs:
readout is initially expressed in the Hamiltonian/carrier frame. The exporter
separately multiplies FIDs by

\[
\exp[-i2\pi(4.68-f_{RF,ppm})\,\nu_{MHz}\,t]
\]

to express every acquisition in the 4.68 ppm analysis frame. This operation
can independently set the acquisition representation **after** simulating RF
response. It cannot reproduce the RF/localization consequences of changing
centreFreq. The experiment preserves that distinction explicitly.

The Hamiltonian uses 42.576 MHz/T, while readout/export metadata use 42.577
MHz/T. Consequently water simulated at RF=3 ppm and re-referenced to 4.68 ppm
has approximately 2.54952 Hz instead of 2.55462 Hz. On-carrier RF=4.7 ppm water
has zero spin offset during evolution, then exactly +2.55462 Hz from the
acquisition transformation. No compensating empirical frequency offset was
introduced. This small gyromagnetic-constant inconsistency is recorded as a
remaining source convention issue.

## Actual configuration

| Parameter | Value |
|---|---|
| Software | GNU Octave 10.3.0; FSL-MRS 2.5.1.dev78+ga311e38f2 |
| Sequence | Universal_Philips, UnEdited PRESS |
| TE / TE1 / TE2 | 30 / 14 / 16 ms |
| Refocusing waveform | univ_eddenrefo.pta, 200 rows, TBW 8.662192932992442 |
| Refocusing duration / nominal angle | 7 ms / 180° |
| Excitation | Ideal `sim_excite(...,'x')` |
| Spatial projection | 101 points in each refocusing direction |
| Refocusing FOV / thickness | 4.5 cm / 3 cm in each direction |
| Native ADC | 8192 points, 4000 Hz, linewidth 1 Hz |
| Final fitting grid | 2048 points, 2000 Hz, 127.731 MHz |
| NAA / Ref0ppm RF | 3.0 ppm |
| Main water RF | 4.7 ppm |
| Control water RF | 3.0 ppm |
| Final analysis frequency origin | 4.68 ppm |
| Concentrations / receiver multiplier | NAA .010 mol/L, water 40 mol/L, G=1 |

The loader reads an excitation waveform, but the actual unedited PRESS branch
uses ideal excitation; a shaped third spatial localization direction is not
applied there. The experiment preserves this bundled behavior. Thus matching
this simulator path and TE does not prove full physical 3D PRESS excitation
or identical sensitivity for different metabolites/carriers.

RF generation is also shared: `io_loadRFwaveform` reads the PTA pulse, and
`experim_funct/calc_shapedRF_propagator_refoc.m` derives w1max from pulse tw1,
duration and nominal flip angle (line 72), scales the waveform envelope to B1max
(lines 128–130), and derives RF time steps from the waveform weights (line
131). This is pulse/flip-angle calibration common to all spin systems, not
independent normalization of water and metabolite signal amplitudes.

## Normalization and export trace

1. Density-matrix scaleFactor already encodes two-proton water. **No additional
   factor of two is applied to simulated H2O.** A separate one-proton control
   halves its initial density matrix solely for the linearity test.
2. `sim_signals` sums spatial responses, divides by nX*nY and applies inverse
   voxel/FOV area ratio (4.5/3)^2. Both acquisitions use the identical rule.
3. `op_dccorr_FID_A` subtracts an edge-estimated complex spectral offset. It
   operates on the readout IFFT spectrum and changes the first FID sample.
4. Assembly recomputes the forward FFT and applies its additional edge DC
   correction (`op_dccorr`). The Python feasibility adapter explicitly
   reproduces this small processing step for each acquisition.
5. A **fixed reference set** NAA, Ref0ppm and off-carrier H2O, all acquired at
   RF=3 ppm, defines the maximum real assembled spectrum. Its recorded scalar
   is **2779.198340948381**. The on-carrier water does not determine that scalar.
   Every output FID, including both water acquisitions and controls, is divided
   by this same scalar. Water/metabolites are never independently normalized.
6. Acquisition-frame shifting is applied per carrier. LCModel entries are
   serialized as FFT(conj(FID)) with six-digit exponent precision, matching
   the bundled export operation; the FSL reader reverses that convention.
   This small writer permits diagnostic names/carriers in one master response
   file rather than using the bundled metabolite-filtering exporter on mixed
   acquisition carriers. It is isolated under examples, not a production export
   change. Round-trip relative errors are tested at 3e-7 for this serialization.
7. All columns are formatted to 2000 Hz / 2048 points with the same FSL-MRS
   formatter. No further frozen/selection-dependent normalization is applied.
   Datasets and their matching fit basis come from those exact common columns.

The bundled basis builder's `BASIS.scale=max(real(buffer.specs))` is verified
at lines 561–563. Its factor can be stored and shared across new acquisitions,
as done here. The old checked-in full-basis file does not store its original
scale. The current experiment uses a fixed three-entry reference instead of
claiming the historic full-29-metabolite maximum.

As a diagnostic comparison only, one complex least-squares scalar aligns the
new NAA column to checked-in frozen NAA: 4773.33324813688+0.0000841457355852i.
Full-array relative difference after that match is **1.43302606e-7**. Applying
the *same* comparison scalar to Ref0ppm yields **1.91316972e-7**. This supports
matching sequence/processing shapes up to export precision, without separate
metabolite corrections. It implies an approximate historical assembly scale
8857.35335631589 using the recorded Phase 4a scalar. This is an inferred global
scale estimate, not recoverable exact historical provenance; **the comparison
scalar is not applied to any fitted water data**. Exact historical recovery
would require the original full raw reference assembly/configuration.

## Response measurements and controls

On the final common response scale, unit-concentration initial samples are:

| Response | Initial complex FID | Magnitude | Initial phase rad |
|---|---|---:|---:|
| NAA, RF=3 | .0009230786561072083 − .0001119923229041747i | — | — |
| Ref0ppm, RF=3 | .00014397095663072873 − .00007357009311645507i | .00016167929661626762 | −.47241398013519376 |
| Water, RF=4.7 | .0003223024224859375 + .0000007387551635742228i | .00032230326914180603 | .002292113922451513 |
| Water, RF=3 | .0003214924550457697 + .00002823104803466795i | .00032272959381577016 | .08758780489431697 |
| Analytical 2*abs(Ref0ppm[0]) | .00032335859323253524 + 0i | .00032335859323253524 | 0 |

Water-on / water-off initial magnitude ratio is **0.9986790034687445** even
after placing both in the same acquisition frame; their phase/sensitivity
differences therefore cannot be reproduced by a post-acquisition frequency
shift alone. Sequence-water / same-scale analytical initial amplitude ratio
is **0.9967363660257814**, with relative phase **.002292113922451513 rad**.
Initial samples here are explicitly *processed*, not a proof of unattenuated
physical sensitivity.

Two water reruns have maximum complex-array difference **0**. The native
two-proton water equals exactly twice the native one-proton control at every
complex sample, maximum difference **0**. Six-digit export introduces ordinary
rounding in the control columns; no extra proton correction is applied.

## First-point and line-shape limitation

Before the second assembly correction, native on-carrier water starts at
0.8964257146226227−0.0001715488248337854i. Its remaining samples are a single
1 Hz exponential whose extrapolated initial response is
1.7906781867670487+approximately 0i. Thus processed sample zero has only
**0.5006068302287714** of the tail-extrapolated magnitude. Constant spectral
DC subtraction has changed that sample; it does not mean a missing water
proton, and multiplying water by two would not repair the array.

Fourier resampling of this discontinuity creates ringing. For example, the
final unit water response begins approximately .0003223, .0006873 and .0006179
in real part, rather than a monotonic exponential. FSL WaterRef assumes one
Voigt FID starting at its fitted amplitude, so its extrapolated water amplitude
and the processed first-sample magnitude disagree. The discrepancy is preserved
and measured, not concealed by amplitude normalization.

The lightweight 256-point test has a first-point/tail ratio approximately
0.9089, rather than 0.5006, because spectral DC estimation also depends on
record duration and the remaining tail. An initial test incorrectly assumed
the long-record ratio would apply to that short fixture. The corrected test
independently calculates the complete bundled DC operation and compares the
full complex array, rather than imposing a universal half-point rule.

## FSL-MRS validation

Three datasets use identical newly simulated NAA at 10 mM and water 40 M,
with the same shared normalization and fitting settings as Phase 5b:
Newton, Lorentzian, baseline off, NAA-only, .2–4.2 ppm; water reference NAA
3 protons, 1.8–2.2 ppm; TE=.03 s, TR=2 s, pure WM quantification defaults.
Actual CLI `--report` runs and supplementary parameter diagnostics were executed.
CLI and API molarity values agree to numerical precision.

| Water reference | Input initial magnitude | Fitted water amplitude | Water relative complex residual | Standard FSL mM |
|---|---:|---:|---:|---:|
| Sequence water, native/processed 1 Hz | .012892130765672243 | .02574996208670833 | .028477089309657616 | **6.891917170702524** |
| Same-scale Ref0ppm analytical, 1 Hz | .01293434372930141 | .012934343456327081 | 1.4370420115447315e-7 | **13.720421395387675** |
| Same-scale Ref0ppm analytical, 0 Hz | .01293434372930141 | .012934343726485146 | 2.6159561521890127e-10 | **6.003318896984144** |

The zero-linewidth analytical control reproduces the old checked-in-basis
Phase 5b result ~6.003318301 mM to export/fit precision despite a different
arbitrary global basis scale. This is the appropriate comparison to the old
analytical convention. The 1 Hz analytical control isolates linewidth effects:
it must not be conflated with the old zero-linewidth experiment. Abs(real)
spectral integration and first-point weighting give line-shape/window-dependent
areas, so changing linewidth can change the reported concentration even at
fixed initial FID amplitude.

All three cases fit the same NAA internal coefficient 1.0000000764364945 and
raw source coefficient .01000000076436493 mol/L. NAA relative residual is
1.857590496917146e-5. Raw source coefficients cannot be compared numerically
to the old value 152.12697 without accounting for the different source-basis
normalization; the underlying concentration interpretation is explicit.

Sequence water fits frequency 2.554581058538584 Hz, phase
3.555375893881878e-5 rad and linewidth .9982039705073966 Hz. Optimizer success
does not erase its 2.85% residual. Its fitted water integral per mol/L is
.015180378328683311, versus the NAA reference integral per mol/L
.020476985425512653. Using the unchanged Phase 5d.1 equation,

\[
C_{FSL}=\frac{.010J_{NAA}}{40J_w}\frac{2}{3}
\times30655.45930654141,
\]

predicts **6.8919164358192875 mM**, within 7.35e-7 mM of the fit result before
including its tiny metabolite coefficient/line-shape errors. The remaining
discrepancy from 10 mM is therefore consistent with the standard tissue,
relaxation and response-integral model plus this processing mismatch. It is
not grounds to tune water amplitude. Physical water response remains unresolved.

![Sequence versus same-scale analytical responses](../outputs/mrscloud_sequence_matched_water/validated/response_comparison.png)

Reports:

- [Sequence water](../outputs/mrscloud_sequence_matched_water/validated/sequence_water/fsl_mrs_results/report.html)
- [Analytical 1 Hz control](../outputs/mrscloud_sequence_matched_water/validated/analytical_reference/fsl_mrs_results/report.html)
- [Analytical zero-linewidth control](../outputs/mrscloud_sequence_matched_water/validated/analytical_reference_zero_lw/fsl_mrs_results/report.html)

## Reproduction, artifacts and tests

Files added: `examples/sequence_matched_water.m`,
`examples/mrscloud_sequence_matched_water.py`, the Octave integration test file,
and this document. Source digests for the pulse, spin systems, simulator,
Hamiltonian, basis builder and exporter are saved in results.json.

Commands actually executed included:

```
pixi run octave --version
pixi run octave --no-gui --quiet --eval "addpath('examples'); sequence_matched_water('/private/tmp/mrscloud_water_smoke',3,256);"
pixi run octave --no-gui --quiet --eval "addpath('examples'); sequence_matched_water('outputs/mrscloud_sequence_matched_water/octave',101,8192);"
pixi run python -m examples.mrscloud_sequence_matched_water outputs/mrscloud_sequence_matched_water --reuse-raw
pixi run python -m examples.mrscloud_sequence_matched_water outputs/mrscloud_sequence_matched_water/final --reuse-raw --raw-dir outputs/mrscloud_sequence_matched_water/octave
pixi run python -m examples.mrscloud_sequence_matched_water outputs/mrscloud_sequence_matched_water/review --reuse-raw --raw-dir outputs/mrscloud_sequence_matched_water/octave
pixi run python -m examples.mrscloud_sequence_matched_water outputs/mrscloud_sequence_matched_water/validated --reuse-raw --raw-dir outputs/mrscloud_sequence_matched_water/octave
pixi run python -m pytest -q tests/integration/test_mrscloud_sequence_matched_water.py
pixi run python -m pytest -q tests
```

For a fresh one-command reproduction, omit `--reuse-raw` and use a new output
folder. Default is 101 spatial points/axis and 8192 native ADC points. The
authoritative comparison after adding both line-shape controls is under
`outputs/mrscloud_sequence_matched_water/validated`; earlier subfolders retain
intermediate comparisons. Native MAT data and the Octave log are in `octave`.
Each dataset includes the true CLI HTML report, input pair, quantification
CSV, fitting output and command log. `cli_commands.txt` records exact fitting
commands with the newly generated common_response.BASIS as the fit basis.

Focused tests: **4 passed**. They execute the real bundled Octave simulation
with a small 3-point-per-axis grid/256 native points, verify deterministic
response, exact proton scaling, acquisition frequency, carrier effects,
common-normalization amplitude preservation, serialization precision and the
DC first-point limitation by independent full-array calculation. The final full
suite completed with **76 passed, 1 xfailed**, zero unexpected failures, in
56.04 s. The expected failure remains the previously documented installed
FSL-MRS combined-phase/frequency initialization defect. Five existing
dependency/deprecation/optional-display warnings were emitted.

Octave emitted expected bundled-code shadowing, obsolete-function and
load-path warnings. A shutdown `ignoring const execution_exception` message
appeared in the full native run despite exit code zero and complete MAT data;
the log is preserved. No dependency changes or production edits were made.

## Required review before replacement

Resolve a single documented physical/first-point processing convention for
metabolite basis and water acquisition, validate a water model compatible with
that processing, and retain a fixed source reference/common normalization.
Do not repair the sample discontinuity with an empirical scale or independently
normalize water. RF carrier, projection/localization, ideal excitation,
sequence sensitivity and the Hamiltonian/export gyromagnetic constants must
remain explicit. Existing analytical calibration is preserved and has not been
replaced. Physical scanner calibration, production integration and LCModel
validation are not claimed or performed in this feasibility phase.

# Phase 5i: raw NAA methyl and response-factor decomposition

The remaining **alpha = 0.8984749391468103** is reproduced exactly from
controlled response/integration ratios. There is no missing NAA proton factor.
At the established RF carrier, raw methyl response per proton differs from
water by only **2.68 parts per million**. The largest contribution is the
finite-record, discrete-grid abs-real spectral integral, followed by basis DC
processing. Full-NAA spectral tails and Fourier formatting partly offset it.
No signal, shared normalization or calibration constant was changed.

## Spin definitions and capture

Definitions are loaded from the bundled
`simulation/basissets/MRSCloud-main/functions/FID-A/simulationTools/metabolites/spinSystems.mat`:

| System | Actual spins | Shifts ppm | scaleFactor |
|---|---:|---|---:|
| sysNAA_Acetyl_only | 3 | 2.008, 2.008, 2.008 | 1 |
| sysNAA_Aspartate_only | 4 | 7.8205, 4.3817, 2.6727, 2.4863 | 1 |
| sysNAA | acetyl and aspartyl subsystems | the same shifts | 1 each |
| sysH2O | 1 represented spin | 4.7 | 2 |

The acetyl J matrix is zero. The aspartyl upper-triangle nonzero couplings are
6.4, 3.861, 9.821 and -15.592 Hz. Complete definitions are retained in
`results.json`. The acetyl simulation evolves the actual three-spin Hamiltonian;
it is not a one-spin FID multiplied by three. `sim_Hamiltonian_mgs` applies
scaleFactor to equilibrium density matrices, and `sim_readout` applies the
documented spin-dimension trace factor `2^(2-nspins)`. No extra proton factor
is applied to any simulated response. Analytical 3H/2H signals below are
explicitly separate diagnostic controls.

`examples/capture_naa_response.m` instruments a private copy of `sim_signals`
by saving `outA` immediately before the first `op_dccorr_FID_A` call. The
original source is untouched. This capture occurs **after** spatial integration,
division by nX*nY and voxel/FOV area adjustment. Source and existing reference
hashes are checked before/after analysis and retained in the output.

Actual configuration: Universal_Philips unedited PRESS TE30, TE1/TE2=14/16 ms,
7 ms universal refocusing waveform, 180-degree refocusing, ideal x excitation,
101 samples per refocusing spatial direction, FOV/thickness 4.5/3 cm, B0=3 T.
Native sampling is 8192 points at 4000 Hz with 1 Hz readout decay. Target
sampling is 2048 points at 2000 Hz. Full NAA and aspartyl use RF=3.0 ppm;
methyl controls use RF=3.0, 4.7 and its own resonance 2.008 ppm. A repeated
RF=3.0 methyl simulation verifies determinism.

All responses use the **unchanged fixed S=2779.198340948381** and existing
sequence calibration identity. The scalar is loaded from the validated
reference manifest, not recalculated from these new selections. Water is read
from the existing reference; it is neither resimulated nor rescaled separately.
Raw and post-DC captures, spectral arrays, carrier-relative spin systems,
waveform, timings and spatial/sample configurations are in `octave/naa_trace.mat`.

## Raw amplitudes and carrier sensitivity

The unnormalized, spatially integrated raw first samples are:

- Full NAA: **5.142476333409348 - 0.055353102715160515i**.
- Methyl, RF3: **2.686010077544293 + 1.9197571189341456e-15i**.
- Methyl, RF4.7: **2.5284093068710094 + 3.677931051966181e-15i**.
- Methyl, RF2.008: **2.6860172801505833 + 3.3947633990010587e-16i**.
- Existing raw water: **1.7906781867670487 + approximately 0i**.

After the shared division, nominal unit-concentration responses (SRU per
mol/L) on the target acquisition grid are:

| Response | Initial complex sample | Initial phase rad |
|---|---|---:|
| Water, RF4.7 | .0006443146429613918 + 3.86158e-20i | 5.99e-17 |
| Methyl, RF3 | .0009664693728291849 + 6.90759e-19i | 7.15e-16 |
| Methyl, RF4.7 | .0009097620956438137 + 1.32338e-18i | approximately 0 |
| Methyl, RF2.008 | .0009664719644420914 + 1.22149e-19i | approximately 0 |
| Full NAA, RF3 | .0018503452084152137 - .000019916931404137092i | -.0107634846540 |

Full NAA contains other spins, so its initial amplitude must not be interpreted
as the three-proton methyl amplitude. The complete native full FID equals
methyl plus aspartyl to maximum complex error **8.881784197001252e-16**.
Repeated methyl simulations differ by **0** at every sample. On the normalized
target grid, full-minus-methyl has maximum magnitude .0008841002074931652
and relative complex norm .2496808662448617 with respect to full NAA.

Water/methyl RF3 initial complex ratio is
**.6666684543507712 - 4.36529e-16i**, versus the ideal 2/3. Methyl response per
proton relative to water is **.999997318481034**. Methyl on its own RF resonance
has per-proton ratio **1.0000000000000038** against water. The common sequence
response therefore preserves 3H versus 2H to numerical precision on carrier.

Moving methyl RF3 to RF4.7 changes magnitude by a factor
**.9413253241337906** (5.8675% reduction); initial relative phase is
7.39918e-16 rad. RF3 versus on-resonance methyl is **.9999973184810301**.
These are localization/RF response comparisons, not proton-count corrections.

Acquisition representation is applied only after sequence evolution:

\[
x_{analysis}(t)=x_{carrier}(t)S^{-1}
e^{i2\pi(f_{RF,ppm}-4.68)127.731t}.
\]

Because the Hamiltonian uses 127.728 MHz while export metadata use 127.731 MHz,
methyl frequencies are respectively **-341.294256, -341.289156 and -341.297232
Hz** for RF3, RF4.7 and RF2.008. No compensating shift was introduced. Each
raw methyl FID matches its independently specified single-resonance exponential
at the correct frequency, with relative norm errors 4.26e-14, 1.05e-13 and
4.35e-14. The genuine frequency difference is kept separate from the magnitude
comparison.

## Incremental processing

The first correction operates on `fftshift(ifft(raw))` and returns
`fft(ifftshift(spec-offset))`; subtracting its complex edge offset changes only
sample zero, by native_points*offset. Assembly's second correction uses the
forward-FFT representation and subtracts its offset from sample zero without
that native-length multiplier. Both operations are independently reproduced
and compared against captures. Neither is an acquisition-amplitude correction.

Corrections are applied in the native carrier frame. Only then are the shared
division and unit-magnitude analysis-frame transformation applied. Direct
target-time sampling selects native indices 0,2,4,... and preserves those
samples exactly. Fourier formatting is a separate operation and changes the
causal record's initial boundary. No first-point repair or extra gain is used.

All integral values below use the **same target grid** and FSL's 1.8–2.2 ppm
NAA reference functional. Full complex arrays are retained in `responses.npz`;
CSV includes sample-zero phase, maximum and relative complex differences.

| Response/stage | First sample | J per nominal mol/L | Relative complex difference from raw |
|---|---|---:|---:|
| Methyl raw | .000966469373 + approximately 0i | .020736642736401665 | 0 |
| Methyl first DC | .000482269765 - .0000484890962i | .020455778689281792 | .0282218922 |
| Methyl second DC | .000482269768 - .0000483014664i | .020455778691340156 | .0282208097 |
| Methyl DC + Fourier | .000482967889 - .000110596611i | .020490407252205828 | .0293582358 |
| Methyl raw + Fourier only | .000725067691 - .0000864458778i | .020631868043540280 | .0159878653 |
| Full NAA raw | .001850345208 - .0000199169314i | .020927846443463932 | 0 |
| Full NAA first DC | .000923466606 - .00000752118158i | .020417861758646230 | .0519515179 |
| Full NAA second DC | .000923470421 - .00000715366657i | .020417863360434822 | .0519515836 |
| Full NAA DC + Fourier | .000923078704 - .000111992334i | .020476984963170750 | .0532853464 |

Methyl's raw initial phase is zero to roundoff, whereas DC/Fourier processing
changes its sample-zero phase to approximately -.2251 rad. A processed
sample-zero phase is therefore not evidence of an intrinsic methyl sequence
phase. Other NAA resonances affect the full response and its edge-estimated
DC offset. Since J uses abs(real), it is nonlinear: the other-resonance
contribution is an ordered comparison, not an additive peak-area attribution.

The captured full-NAA post-first-DC FID is **exactly equal** to the preserved
Phase 5e NAA response. Including assembly correction and Fourier formatting
matches the existing serialized basis to better than the tested 3e-7 relative
complex precision. The final tiny serialization-integral contribution is
measured separately below; the existing production basis remains unchanged.

## Finite-grid integration and alpha

Let J denote the independent NumPy implementation of FSL's orthonormal FFT,
half-first-point weighting, inclusive nearest-bin window selection, and
trapezoidal abs-real integral. Water is integrated in its fitted
frequency/phase-corrected frame over 1.65–7.65 ppm. The unit response integral
is **.015193883571167494**; integrating the original acquisition-frame water
instead gives .014625686547203712. The reference-frame distinction is retained.

The ideal analytical methyl control has three-proton amplitude 1.5 times
unit-water amplitude, 1 Hz decay, zero phase, and the **actual 2.008 ppm spin's
Hamiltonian/export frequency**. Its baseline ratio is **.9098704420596527**.
This differs from Phase 5h's **.9577941910199874** analytical control at
approximately 2.01 ppm. That earlier rounded-frequency control was explicitly
analytical; its difference cannot be assigned to native sequence sensitivity.
Finite-record bin alignment matters strongly for this narrow line.

Analytical-only record-length controls, holding dwell and frequency fixed:

| Points at 2000 Hz | Ideal 3H/2H window ratio |
|---:|---:|
| 2048 | .9098704420596526 |
| 4096 | .9887567195834266 |
| 8192 | .9882280042205992 |
| 16384 | .9882596738614721 |
| 32768 | .9882461104085961 |

At the end of the baseline 1.024 s record the 1 Hz FID is still about 4% of
its initial amplitude. Frequency-bin alignment, truncation and abs-real
integration explain why its ratio differs substantially from a long record.
The remaining long-record deficit is consistent with unequal finite reference
windows and discrete endpoint sampling; it should not be called RF loss.
These are diagnostic analytical records, not changes to the exported datasets.

The following ordered factors are derived from their respective response
functions, not selected to fit an observed concentration:

| Contribution | Factor |
|---|---:|
| Ideal finite-grid/window ratio at actual spin frequency | .9098704420596527 |
| Native methyl sequence amplitude per proton | .999997318481034 |
| Methyl line shape versus analytical exponential | .9999999999999984 |
| Intrinsic methyl phase | 1 |
| Other NAA resonances, raw full/raw methyl J | 1.0092205719842307 |
| Full NAA first DC J change | .9756312869460596 |
| Full NAA second DC J change | 1.0000000784503593 |
| Fourier formatting J change | 1.0028955822503198 |
| Existing BASIS serialization J change | 1.000000022578612 |
| **Product** | **.8984749391468103** |

Direct calculation from the unchanged final basis gives exactly the same:

\[
\alpha=\frac23\frac{J[B_{NAA}]}{J[W]}
=\frac23\frac{.020476985425512653}{.015193883571167494}
=.8984749391468103.
\]

The last serialization factor compares captured and existing basis response
integrals directly, not alpha or a fitted concentration. All multiplicative
contributions are conditional on this declared processing order and nonlinear
integration functional. They are not transferable physical correction factors.

## Actual CLI fits and reports

Two optional controls were fitted with the installed FSL-MRS
`2.5.1.dev78+ga311e38f2` CLI and **`--report`**, Newton/Lorentzian, baseline off,
ppm .2–4.2, NAA only, TE30/TR2, pure WM, reference protons 3 and reference
window 1.8–2.2. No matched-density/relaxation overrides were introduced here.
Both use exactly the unchanged Phase 5g water-reference FID.

| Control | Raw fitted coefficient mol/L | NAA complex residual | Water complex residual | Standard concentration mM |
|---|---:|---:|---:|---:|
| Existing formatted full NAA | .01000000076436493 | 1.85759050e-5 | 1.44416231e-7 | 6.885791303512254 |
| Raw isolated methyl | .009999999417418219 | 5.61978691e-7 | 1.44517852e-7 | 6.973105721414660 |

All optimizers report success. Actual CLI CSV coefficients and concentrations
equal the supplemental Python diagnostics exactly. Water fitted amplitudes are respectively
.02577258514572865 and .0257725851413698 versus expected .02577258571845567.
Water/NAA fitted integration-model areas (in each fit's internally scaled
units) are 341832.1888333719/115.17288435931849 for full NAA and
352471.20548623084/120.26335774483883 for methyl. Different data norm scalings
account for the different absolute water integrals; original water samples
are asserted equal, not inferred equal from fitted concentrations.

An initial methyl CLI run failed in the installed LCModel reader because it
squeezes a one-column BASIS to one dimension, then indexes `data.shape[1]`.
The successful control BASIS includes the existing shared-scale Ref0ppm column
and still fits only NAA (`--keep NAA`). No package or production exporter was
changed. The first Python attempt also exposed an experiment import typo,
corrected before these results. Failed logs remain in the earlier output root.

Actual nonempty reports:

- [Full NAA report](../outputs/mrscloud_naa_response_decomposition/validated/full_formatted/fsl_mrs_results/report.html), 887855 bytes.
- [Raw methyl report](../outputs/mrscloud_naa_response_decomposition/validated/methyl_raw/fsl_mrs_results/report.html), 893230 bytes.

Per-control sidecars explicitly mark diagnostic-only data and identify the
actual metabolite fitting basis; the original reference calibration remains
recorded. Fitted nuisance parameters, amplitudes, integrals, and residuals
are retained in `fitting_results.json`. The familiar water HTML frame mismatch
must not be interpreted as a failed numerical water fit.

![Raw and processed response comparisons](../outputs/mrscloud_naa_response_decomposition/validated/response_processing.png)

![Carrier controls and ordered factors](../outputs/mrscloud_naa_response_decomposition/validated/carrier_and_factors.png)

## Reproduce and tests

Files added: `examples/capture_naa_response.m`,
`examples/mrscloud_naa_response_decomposition.py`,
`tests/test_naa_response_decomposition.py`,
`tests/integration/test_naa_response_capture.py`, and this document.
No production simulation, MRSCloud source, GUI, reference or dependency file
was modified.

Actual capture and final analysis commands:

```bash
pixi run octave --no-gui --quiet --eval "addpath('examples'); capture_naa_response('outputs/mrscloud_naa_response_decomposition/octave',101,8192);"
pixi run python -m examples.mrscloud_naa_response_decomposition outputs/mrscloud_naa_response_decomposition/validated --reuse-capture
pixi run python -m pytest -q tests/test_naa_response_decomposition.py tests/integration/test_naa_response_capture.py
pixi run python -m pytest -q tests
```

The validated directory reuses a copy of the original `octave/naa_trace.mat`
and log. For a fully automatic fresh reproduction, run the Python module with
a new output directory and omit `--reuse-capture`: it invokes the configured
Pixi/Octave capture itself. Complete actual FSL CLI commands are retained in
`validated/cli_commands.txt`. Octave's known bundled load-path/deprecation and
shutdown messages remain in its log; the process exited successfully.

Focused tests: **7 passed**. They execute a small 3-spatial-point/256-readout
Octave fixture, independently verify frequency/exponential/DC behavior,
superposition, three-spin metadata and deterministic repeats, and inspect the
actual full-resolution experiment and CLI reports. Analytical tests verify
shared-scalar linearity, exact target-time selection, finite-grid controls,
and the response decomposition. Full regression execution completed with
**119 passed, 1 xfailed, 5 dependency warnings in 81.10 s**. The expected failure
remains the documented installed FSL-MRS combined phase/frequency initialization
defect. Existing analytical, voxel-pair, frozen-basis, water-calibration,
processing, paired export, CLI and quantification regressions all executed.

## Interpretation and remaining limits

The Phase 5h 10.15% concentration difference is numerically explained without
changing water gain: native methyl proton response is consistent, while the
specific FSL finite-record integration and basis-processing representations
are not a universal per-proton molar-area convention. Correct density and
relaxation settings alone therefore do not guarantee exact native FSL scaling.

Carrier response can depend on off-resonance localization, as the RF4.7 methyl
control demonstrates. This does not establish physical PRESS sensitivity:
excitation is ideal, localization follows the bundled two-refocusing-direction
model, and the Hamiltonian/export gamma mismatch remains. Other-resonance
effects include intrinsic J evolution, phase and spectral tails; their ordered
integral contribution does not uniquely assign each spin's physical peak area.
Generalization to other sequences, metabolites, windows, ADC filters and fitting
software remains unverified. Review these conventions before choosing any
future quantification adapter. No empirical correction, physical calibration,
production integration or LCModel validation is implemented in this phase.

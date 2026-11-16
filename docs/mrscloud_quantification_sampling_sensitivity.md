# Phase 5j: sampling and linewidth sensitivity

Increasing acquisition duration reduces the large 1 Hz finite-record bias but
does **not** make processed-basis FSL-MRS water referencing converge to 10 mM.
At 8192 points, the predicted concentration is 9.74561 mM for 1 Hz, 9.55491 mM
for 3 Hz, 9.34058 mM for 5 Hz and 9.02111 mM for 8 Hz. Wider lines substantially
reduce sampling sensitivity while increasing losses outside the fixed NAA
integration window. DC correction and Fourier boundary effects remain distinct
from physical sequence response. No calibration correction was applied.

## Conditions and acquisition extension

FSL-MRS version: `2.5.1.dev78+ga311e38f2`; Pixi Python 3.14.8 and configured
Octave. All twelve conditions retain:

- PRESS TE30, TE1/TE2=14/16 ms, universal 7 ms refocusing waveforms, ideal
  excitation, 101 samples per refocusing direction.
- NAA RF carrier 3.0 ppm; water RF carrier 4.7 ppm; analysis origin 4.68 ppm,
  frequency metadata 127.731 MHz, positive `exp(+2j*pi*f*t)` convention.
- NAA 10 mM per nominal homogeneous tissue volume; water 40 mol/L; receiver G=1.
- Existing sequence-response normalization **S=2779.198340948381**, from the
  original fixed NAA/Ref0ppm/off-carrier-water reference. S is never recomputed.
- Water/reference proton convention 2/3, NAA integration 1.8–2.2 ppm,
  corrected-water integration 1.65–7.65 ppm.
- Effective quantification water 40 M; no T1/TR or TE/T2 attenuation correction.
  Actual TE=.030 s and TR=2 s are retained. This does not remove sequence
  evolution or acquisition-time linewidth.
- No noise, residual water, spatial weighting changes or empirical gain factors.

| Points | Dwell s | Record length N*dt s | Last sample time s | Bin spacing Hz |
|---:|---:|---:|---:|---:|
| 2048 | .0005 | 1.024 | 1.0235 | .9765625 |
| 4096 | .0005 | 2.048 | 2.0475 | .48828125 |
| 8192 | .0005 | 4.096 | 4.0955 | .244140625 |

The previous 8192-point native responses at 4000 Hz covered only 2.048 s.
The experiment therefore captured **16384 genuine native readout samples** at
4000 Hz using the existing private NAA and water capture wrappers. It does not
zero-fill, repeat or analytically extrapolate full NAA. The extended raw full
NAA and water prefixes equal their previous 8192 native samples **exactly**.
Pulse/localization evolution and raw amplitudes have not changed.

## Fixed processing, shared scale and broadening

One longest-record metabolite master is prepared before the sweep:

1. Capture pre-DC full NAA and actual three-spin methyl responses.
2. Apply the captured first DC and independently reproduced second DC correction
   in the native carrier frame.
3. Divide every response by the same existing S and shift the acquisition frame.
4. Fourier-format metabolites once from 16384/4000 Hz to 8192/2000 Hz.
5. Use prefixes of this fixed master for 2048, 4096 and 8192 conditions.

Raw water bypasses DC and Fourier formatting; its target master uses exact
native samples 0,2,4,... with the same S and acquisition-frame transformation.
The main metabolite fitting basis is unbroadened; its coefficient-weighted
model fits the additional linewidth. A shared-scale methyl control column is
included solely to avoid the installed reader's one-column BASIS squeeze bug;
the fit selects NAA only.

The swept linewidth L=1,3,5,8 Hz means **total Lorentzian readout linewidth**.
Native simulation already contains 1 Hz. Both target FIDs receive the identical
additional factor

\[
E_L(t)=e^{-\pi(L-1)t},\quad
s_m(t)=.010\,B(t)E_L(t),\quad s_w(t)=40\,W(t)E_L(t).
\]

E_L(0)=1, so initial complex amplitudes are preserved exactly. It is not a
TE/T2 correction. There is no extra water factor of two: sysH2O already encodes
its two-proton response. Basis-formatting artifacts remain in B, so “linewidth”
describes the underlying resonances plus this common envelope, not a claim
that the complete processed NAA FID is a single exponential.

Every exported condition has the identical initial samples:

- Metabolites: **9.231164979192688e-6 - 1.1227377262236918e-6i**.
- Water: **.02577258571845567 + approximately 0i**.

The extended, once-processed master differs from Phase 5i's shorter processed
basis by relative norm **.00018491090152511795** in the first 2048 samples.
Its unit-response initial difference is 3.784181206042091e-8 -
2.814497181944666e-7i. This is the documented change in DC/Fourier representation
for the longer native record, not a new normalization. The 2048/1 Hz fitted
matched result is 8.9847896765 mM versus Phase 5h's 8.9847504611 mM.
Processing is fixed across all twelve new conditions, so it is not silently
recomputed or renormalized when N or L changes.

This standalone experiment does not bypass the production bundle's sampling
compatibility checks. Sidecars identify the original amplitude-reference
manifest, fixed scalar, extended native grid and current target grid, and
explicitly do not claim original production-grid compatibility.

## Independent predictions

J is the independently implemented FSL abs-real, half-first-point,
orthonormal shifted-FFT trapezoidal integral over inclusive nearest-bin limits.
For each grid and linewidth, the prediction is calculated **before fitting**:

\[
J_m=J[B E_L;(1.8,2.2)],\quad
J_w=J[|W(0)|e^{-\pi Lt};(1.65,7.65)],\quad
\alpha_{N,L}=\frac23\frac{J_m}{J_w},\quad
\widehat C_{pred}=10\alpha_{N,L}\ \mathrm{mM}.
\]

Water is in the corrected integration frame, not its physical acquisition
frequency. No optimizer parameter or fitted concentration enters these
predictions. Integrals are per nominal mol/L, in FSL spectral-bin units;
absolute numbers depend on FFT length and bin spacing, while their same-grid
ratio cancels the common scale. They must not be interpreted as physical
spectral-area changes across different grids without those factors.

| Points | Total Hz | J_metabolite | J_water | Response factor |
|---:|---:|---:|---:|---:|
| 2048 | 1 | .0204770738034 | .0151938835712 | .898478816932 |
| 2048 | 3 | .0208596200361 | .0145478786105 | .955906612214 |
| 2048 | 5 | .0203724992421 | .0145260464384 | .934987108778 |
| 2048 | 8 | .0196462015852 | .0144941692728 | .903637005812 |
| 4096 | 1 | .0302512515552 | .0206235106086 | .977888848294 |
| 4096 | 3 | .0294831088945 | .0205730594473 | .955395379090 |
| 4096 | 5 | .0287760839930 | .0205430540440 | .933846347982 |
| 4096 | 8 | .0277255803323 | .0204980483273 | .901730736817 |
| 8192 | 1 | .0425939217619 | .0291371724528 | .974560859921 |
| 8192 | 3 | .0416995607304 | .0290946739304 | .955491277204 |
| 8192 | 5 | .0407046746009 | .0290522228347 | .934057606826 |
| 8192 | 8 | .0392263303330 | .0289885494597 | .902110915382 |

## Fits and numerical errors

Each condition uses Newton, Lorentzian, baseline off, NAA only, ppm .2–4.2,
matching grid/basis, and identical public quantification settings. The Python
API sets effective WM density to 40000/55010 and uses the same synthetic
no-attenuation limit as Phase 5h: T1=1e-6 s and T2=infinity. Both attenuation
factors evaluate exactly to 1. Fitting-time data/basis norm scaling is recorded;
the comparable source coefficient is `internal * basis_scale / data_scale`.

| Points | Hz | Predicted mM | Fitted mM | Fit minus prediction mM | NAA relative residual | Water relative residual |
|---:|---:|---:|---:|---:|---:|---:|
| 2048 | 1 | 8.984788169 | 8.984789677 | 1.51e-6 | 1.86e-5 | 1.46e-7 |
| 2048 | 3 | 9.559066122 | 9.559065189 | -9.33e-7 | 1.41e-6 | 1.07e-9 |
| 2048 | 5 | 9.349871088 | 9.349871026 | -6.21e-8 | 1.45e-6 | 1.04e-7 |
| 2048 | 8 | 9.036370058 | 9.036663373 | 2.93e-4 | 1.58e-6 | 3.26e-5 |
| 4096 | 1 | 9.778888483 | 9.778895118 | 6.63e-6 | 1.65e-5 | 9.20e-8 |
| 4096 | 3 | 9.553953791 | 9.553955701 | 1.91e-6 | 4.42e-7 | 2.21e-7 |
| 4096 | 5 | 9.338463480 | 9.338464283 | 8.03e-7 | 1.12e-6 | 5.99e-7 |
| 4096 | 8 | 9.017307368 | 9.017592037 | 2.85e-4 | 1.35e-6 | 3.17e-5 |
| 8192 | 1 | 9.745608599 | 9.745615437 | 6.84e-6 | 1.86e-5 | 3.47e-7 |
| 8192 | 3 | 9.554912772 | 9.554907930 | -4.84e-6 | 6.96e-7 | 1.83e-7 |
| 8192 | 5 | 9.340576068 | 9.340570468 | -5.60e-6 | 3.34e-6 | 1.70e-7 |
| 8192 | 8 | 9.021109154 | 9.021485472 | 3.76e-4 | 4.28e-6 | 4.31e-5 |

All optimizers report success; success alone is not the accuracy criterion.
The maximum fitted-versus-unfitted prediction difference is
**.0003763183733695996 mM**, at 8192/8 Hz. The 8 Hz water fits systematically
retain larger numerical errors under the unchanged installed initialization
and optimizer settings. Their relative residuals are explicitly retained,
not called machine-precision fits. No per-condition initialization/bounds were
tuned. The largest source-coefficient relative error is **7.43458419e-7**.
Full coefficients, water amplitudes, parameters, optimizer messages, residuals
and export errors for all twelve conditions are in `comparison.csv` and
per-condition `fitting_results.json`.

Complex export/reload errors are **zero for every condition**, and direct
comparisons verify identical initial complex samples across N and L.

## Separating sampling, linewidth and processing

Unfitted response controls retain raw full NAA, isolated methyl, both-DC
responses before Fourier formatting, and analytical 3H/2H controls. All share
the same S, linewidth factor and target grid. At 1 Hz:

| Points | Ideal methyl | Raw methyl | Raw full NAA | Full NAA after DC | Formatted full NAA |
|---:|---:|---:|---:|---:|---:|
| 2048 | .9098704421 | .9098680022 | .9182575056 | .8958974543 | .8984788169 |
| 4096 | .9887567196 | .9887540682 | .9984962400 | .9750600411 | .9778888483 |
| 8192 | .9882280042 | .9882253543 | .9979865001 | .9745156898 | .9745608599 |

Thus most 2048/1 Hz bias already exists in the ideal finite-record control.
The near-unity raw full-NAA long-record result includes other-resonance tails;
it is not proof that every component or window has physical molar calibration.
DC processing leaves a deficit in the full processed basis even as truncation
of the underlying exponential becomes small.

The Fourier-formatted master also has periodic boundary ringing: its final
sample is approximately -.000136665 + .0000173369i, whereas raw full NAA's
final sample is 2.85e-10 + 2.36e-9i. The 8192 condition includes that master
boundary; shorter prefixes do not. This is a **basis-processing interaction
with acquisition length**, not RF attenuation or a change in receiver gain.
It contributes to the non-monotonic 4096→8192 formatted response. Wider-line
envelopes suppress this late artifact. The DC-versus-formatted control columns
expose it rather than assigning all point-count dependence to FFT-bin spacing.

The factor range across the three N values is .07941003 at 1 Hz, .00051123
at 3 Hz, .00114076 at 5 Hz and .00190627 at 8 Hz. Wider lines reduce finite-grid
and boundary sensitivity, but increasingly fall outside the fixed 1.8–2.2 ppm
reference window. At 8192 points, the ideal methyl factors are .9882280,
.9646688, .9411572 and .9061413 for 1/3/5/8 Hz: bias remains even without
NAA DC processing. Full-NAA other-resonance contributions partially offset
that window loss; their abs-real integral contribution is not additive or a
universal proton multiplier.

There is consequently no evidence that longer acquisition alone drives the
declared processed-reference factor to one. Three durations do not establish
a strict asymptotic limit, especially for a Fourier-formatted finite master.
Likewise, realistic linewidth stabilizes sampling sensitivity, not absolute
concentration recovery under a narrow fixed integration window.

## Representative CLI HTML reports

Three actual CLI runs used `--report`. Unlike earlier baseline-only reports,
these runs represent the **matched** assumptions. The public CLI lacks the
Python density setter, so its explicitly recorded configuration uses
`--h2o_scale 1.118677722931495`, derived solely as **40000/(55010*.65)**.
It is a tissue-water concentration configuration adapter, not an empirical
signal or fitted-concentration correction. CLI T1 JSON uses 1e-6 s and T2 JSON
uses finite 1e20 s; attenuation factors are exactly one in floating point,
matching the API infinity limit without nonstandard JSON Infinity literals.
TE/TR, proton counts, integration windows and all fitting options stay fixed.
API quantification uses density override directly and additional correction 1.
CLI and API reported molarity agree within 1e-12 relative tolerance.

- [2048 points, 1 Hz](../outputs/mrscloud_sampling_sensitivity/validated/N2048_LW1/fsl_mrs_results/report.html), 887814 bytes.
- [8192 points, 1 Hz](../outputs/mrscloud_sampling_sensitivity/validated/N8192_LW1/fsl_mrs_results/report.html), 3314288 bytes.
- [2048 points, 8 Hz](../outputs/mrscloud_sampling_sensitivity/validated/N2048_LW8/fsl_mrs_results/report.html), 880095 bytes.

The known corrected-water-model versus original-data HTML overlay convention
still applies. Numerical residuals above use consistent full complex frames.

![Sampling/linewidth bias and independent processing controls](../outputs/mrscloud_sampling_sensitivity/validated/sampling_linewidth_bias.png)

## Execution and files

Files added: `examples/mrscloud_quantification_sampling_sensitivity.py`,
`tests/test_quantification_sampling_sensitivity.py`,
`tests/integration/test_quantification_sampling_sensitivity_e2e.py`, and this
document. Existing capture helpers were reused without edits. No production
simulator, calibration, GUI, exporter, package, or dependency change was made.

Exact experiment/test commands:

```bash
pixi run python --version
pixi run octave --no-gui --quiet --eval "addpath('examples'); capture_naa_response('outputs/mrscloud_sampling_sensitivity/octave/naa',101,16384); capture_water_processing('outputs/mrscloud_sampling_sensitivity/octave/water',101,16384);"
pixi run python -m examples.mrscloud_quantification_sampling_sensitivity outputs/mrscloud_sampling_sensitivity/validated
pixi run python -m pytest -q tests/test_quantification_sampling_sensitivity.py tests/integration/test_quantification_sampling_sensitivity_e2e.py
pixi run python -m pytest -q tests
```

Use a new output directory when reproducing. `--capture-dir` can point to
another compatible extended raw capture. The original capture log, including
bundled Octave warnings, is preserved under `outputs/mrscloud_sampling_sensitivity/octave`.
The capture exited successfully. Complete actual CLI commands, including
`--report` and matched-assumption inputs, are in
`outputs/mrscloud_sampling_sensitivity/validated/cli_commands.txt`.

The output contains 12 pairs, matching fitting bases, simulation and ground-truth
sidecars, per-case fit diagnostics, three real HTML report folders,
`master_responses.npz`, `comparison.csv`, `results.json`, source/capture hashes
and the bias/control figure. Source reference hashes remain unchanged.

Focused tests completed **6 passed**. They test exact initial-amplitude
preservation, common-gain cancellation, deterministic generation/refitting,
independent NumPy integral predictions without final FSL reporting calls,
12-condition fitting error budgets, actual CLI/API agreement and complex export.
The complete regression run finished with **125 passed, 1 xfailed, 5 dependency
warnings in 87.08 s**. The expected failure is the previously documented
installed FSL-MRS combined phase/frequency water-initialization defect.
Existing analytical, voxel, frozen-basis, calibration, processing, export,
CLI-report and quantification regression suites all executed.

## Remaining scientific limits

This is a synthetic signal/quantifier sensitivity study. Fixed first-point/DC
and Fourier processing are fitting representations, not demonstrated physical
ADC behavior. Their artifacts, unequal windows, intrinsic metabolite phase and
J-evolved tails must not be absorbed into empirical water gains. The 8 Hz
optimizer accuracy limit is measured and retained. Ideal excitation,
two-direction localization, the Hamiltonian/export gyromagnetic mismatch,
and physical scanner calibration remain unresolved. Review the integration
and basis representation conventions before changing the simulator; LCModel,
production integration and calibration changes are outside this phase.

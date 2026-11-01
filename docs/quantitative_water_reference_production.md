# Phase 8: opt-in production quantitative water referencing

Phase 9 distribution and current commands are documented in
[the production user guide](quantitative_simulator_user_guide.md). The numerical
results below remain the historical Phase 8 record. Production lookup now uses
the packaged bundle and checks archived construction-source bytes, rather than
requiring live research files at their old locations. Signal equations are unchanged.

## Integration interface recorded before implementation

The production boundary is `simulation.headless.simulate_batch`, before legacy
basis/phantom construction. `settings['simulation_mode']` is `legacy` by default;
`sequence_matched_quantitative` dispatches to a separate production adapter.
The legacy `SignalModel`, normalization, volume-summed equations, residual-water
model, RNG behavior and existing exporters remain unchanged.

The new adapter accepts a `quantitative` configuration with a complete Phase 6a
bundle, explicit tissue metabolite profiles in mM, optional tissue water profiles
in mol/L, shared receiver factor, relaxation, optional normalized spatial cells,
declared B0 units/sign, independent additional linewidth and acquisition noise.
It returns the established paired-FID dictionary with ground truth and metadata.
A `DigitalPhantom` explicitly constructed in the new mode loads anatomy but
bypasses legacy basis, signal-model and lipid initialization; its quantitative
`simulate_data` returns that pair. Existing mode retains its five-element tuple.

Headless batches can use explicit fractions/cells without anatomy or use an
existing phantom/loaded skeleton and selected voxel centers/size. New-mode paired
exports use `gen_nifti_mrs` with the validated complex convention; legacy
component selectors/exporters are rejected for quantitative results rather
than silently calling them. The CLI reads a JSON configuration. No GUI changes
or fitting experiments are part of this integration.

## Audit of existing production interfaces

- `DigitalPhantom.__init__` loads a TorchIO subject, index-to-name tissue labels,
  spacing and affine, then creates `Basis`, `SignalModel` and lipid masks.
- `simulate_data` crops exclusive-stop six-index ROI bounds, one-hot labels,
  lipid mask and optional B0, then calls the legacy signal model.
- Legacy tissue metabolite CSV concentrations are sampled with Torch RNG;
  their definitions are not changed or inferred for the quantitative mode.
- Legacy signals multiply cell spectra by anatomical cell volume and sum cells.
  The new mode uses normalized represented-volume weights, with no volume gain.
- Legacy shim multiplies FIDs by exp(-i*2*pi*map*t). Quantitative positive Hz
  multiplies by exp(+i*2*pi*offset*t). B0 loading requires a unit/sign declaration;
  T1/T2/B0 map names alone do not establish units or justify implicit conversion.
- Existing headless batch settings combine `simulation_params`, basis settings
  and voxel size; seeds temporarily control/restore NumPy and Torch RNG states.
- The legacy headless spectrum exporter follows its existing FFT behavior.
  The quantitative branch exports paired FIDs directly and never exports the
  legacy residual-water component as an unsuppressed reference.

Implementation, validation results and executable configuration follow below.

## Implemented entry points

- `simulation.quantitative.simulate_quantitative(config, phantom=None,
  voi_bounds=None)` returns complex metabolite/water FIDs, clean per-metabolite
  components, time grid, ground truth and simulation configuration.
- `simulation.quantitative.export_quantitative(pair, output_dir)` exports paired
  FIDs using the validated `gen_nifti_mrs` convention and refuses existing targets.
- `simulation.headless.simulate_batch(..., settings, phantom=None)` dispatches on
  the explicit `simulation_mode`. Legacy remains the default. Quantitative
  results are pair dictionaries, not legacy component-spectrum dictionaries.
- `simulation.headless.run_headless(config, phantom=None)` and
  `pixi run python -m simulation.headless --config CONFIG.json` provide JSON
  execution. `--output-dir` can select a fresh destination without editing JSON.
- `DigitalPhantom(..., path2metabs=None, path2basis=None,
  simulation_params={'simulation_mode':'sequence_matched_quantitative',
  'quantitative':config})` loads the selected skeleton but does not initialize
  legacy `Basis`, `SignalModel` or lipid masks. `simulate_data(bounds)` returns
  a pair in this explicit mode. Legacy construction and tuple return remain.

An already loaded legacy phantom can also be supplied to quantitative headless
simulation: only its anatomy/getters are used, not its signal model. The
production adapter delegates to Phase 6b/6e/6f implementations; it does not
copy their signal equations. Existing analytical reference code is unchanged.
Legacy GUI-style updates and legacy component-export functions reject
quantitative objects so a reference cannot be silently confused with residual
water. No GUI controls were added.

## Configuration schema

A minimal executable JSON configuration is:

```json
{
  "simulation_mode": "sequence_matched_quantitative",
  "quantitative": {
    "calibration_bundle": "calibration/mrscloud_press_te30_quantitative",
    "tissue_fractions": {"WM": 0.2, "GM": 0.3, "CSF": 0.5},
    "metabolite_concentrations_mM": {
      "WM": {"NAA": 10},
      "GM": {"NAA": 14},
      "CSF": {"NAA": 0}
    },
    "receiver_factor": 1,
    "relaxation": {"mode": "none", "TE_s": 0.030, "TR_s": 2.0},
    "additional_linewidth_hz": {"metabolites": 0, "water": 0}
  },
  "output_dir": "outputs/new_quantitative_voxel"
}
```

`quantitative` accepts only the following keys; unknown/legacy settings are
rejected rather than silently ignored:

| Key | Meaning/default |
|---|---|
| `calibration_bundle` | Optional explicit bundle directory; default packaged calibration resource |
| `metabolite_concentrations_mM` | Required WM/GM/CSF mappings with identical selected names and explicit zeros |
| `water_concentrations_mol_l` | Optional complete WM/GM/CSF profile; defaults below |
| `receiver_factor` | Shared nonnegative finite gain, default 1 |
| `bandwidth_hz`, `points` | Fixed supported grid: 2000 Hz, 2048 complex points |
| `expected_calibration_identity` | Optional explicit identity check; only the pinned Phase 6a identity is accepted |
| `concentration_convention` | Only `nominal_per_total_tissue_volume` is accepted |
| `tissue_fractions` | Finite nonnegative volume fractions summing to 1, tolerance 1e-8 |
| `spatial` | Inline anatomical-cell selection and optional declared B0/weights/volumes |
| `phantom_spatial` | Options for spatial integration over a supplied/loaded phantom ROI |
| `relaxation` | Omitted/None: no attenuation; otherwise explicit Phase 6e configuration |
| `additional_linewidth_hz` | Nonnegative additional Lorentzian widths, separately for `metabolites` and `water` |
| `noise` | Omitted/None: disabled; otherwise explicit seed and independent sigma/SNR options |
| `voxel_affine_mm` | Optional nonsingular homogeneous 4x4 voxel affine for non-anatomical cases |

Select fractions, inline cells or phantom spatial options exclusively. Without
fractions or inline cells, supply a phantom or a `skeleton` in the outer headless
settings. Anatomical batches require outer `positions` (integer cell-center
indices) and `voxel_size_mm`. Voxel counts are rounded from size/spacing as in
the existing headless API; metadata/export use the actual covered extent.
Bounds are `(x0,x1,y0,y1,z0,z1)` with exclusive stops. Unknown/background/fat/skull
labels in the quantitative ROI are rejected, not dropped and renormalized away.

For standalone fractions/cells, the default position `(0,0,0)` is an output
identifier. With no affine supplied, identity geometry is explicitly labeled
as lacking anatomical registration. It does not imply a measured anatomical
volume. Existing output directories are not overwritten.

## Calibration and approximately physical relative response

The pinned identity is
`mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81`.
The scalar remains **8857.353925964408**, calculated from the fixed complete
reference set in Phase 6a. Selection never recomputes an applied normalization.

The loader verifies required file hashes, full reference arrays/order, recorded
normalization statistic, water identity, sequence/sampling and FSL/NumPy
formatting-version compatibility. Production additionally verifies the actual
**403 simulation/source files and 5 processing-source files** against recorded
hashes. These source paths resolve within the repository; missing or modified
sources, different calibration families or unsupported grids are rejected.
Source verification is read-only, with no global cache or RNG mutation.

PRESS TE/TE1/TE2 are 30/14/16 ms. The original sequence uses ideal excitation,
7 ms refocusing and the same spatial localization configuration. Metabolites
retain RF carrier 3.0 ppm; unsuppressed water retains RF carrier 4.7 ppm.
Acquisition origin is 4.68 ppm at 127.731 MHz, giving intrinsic water offset
+2.55462 Hz. Water uses raw pre-DC sequence response, one shared normalization,
acquisition-frame transformation and exact acquisition-time samples. Neither
basis DC operation nor Fourier basis resampling is applied to water.

Native spin-system response already includes proton multiplicity, including
sysH2O's two-proton scale. No additional proton multiplier, fitting-derived
water gain or independent water normalization is applied. This common sequence
and tissue-density framework makes relative responses more physically motivated
than an empirical residual-water amplitude. Receiver units remain arbitrary;
scanner-level absolute calibration is **not** established.

## Tissue water and concentration conventions

Defaults reuse the validated `H2O_CONCENTRATIONS` convention:

| Tissue | Water mol/L per total tissue-compartment volume |
|---|---:|
| WM | 39.0165 = 0.703 × 55.5 |
| GM | 47.2305 = 0.851 × 55.5 |
| CSF | 55.5 |

Explicit overrides use mol/L with the same convention. There is no implicit
conversion to concentration per tissue-water volume, nor reapplication of
these water fractions. Metabolites are nominal mM per total volume of each
tissue compartment. The legacy literature/CSV concentration definitions are
not reinterpreted; quantitative tissue profiles must be supplied explicitly.

For normalized spatial weights w_i and tissue j(i), production uses

\[
s_m(t)=G E_m(t)\sum_i w_i\sum_m10^{-3}C_{m,j(i)}R_{m,j(i)}B_m(t)e^{i2\pi\Delta f_i t},
\qquad
s_w(t)=G E_w(t)\sum_i w_i C_{w,j(i)}R_{w,j(i)}W(t)e^{i2\pi\Delta f_i t}.
\]

Fractions use f_j in place of the cell sum. B_m and W are common-normalized
unit-concentration responses. No absolute anatomical-volume multiplier is
used: these are represented-volume averages. Cell volume changes weighting,
not a hidden gain. Ground truth preserves nominal voxel averages, WM/GM-only
concentrations (null for pure CSF), effective water, tissue profiles and
relaxation-weighted signal concentrations separately.

## Relaxation and linewidth

Optional pre-acquisition attenuation is
`R=(1-exp(-TR/T1))*exp(-TE/T2)` with all times in seconds. TE must match .030 s;
TR must exceed TE, and all T1/T2 values must be finite and positive. Active
profiles explicitly specify every selected metabolite and water for WM/GM/CSF,
including zero-volume tissues. No literature relaxation defaults are inserted.
The example configuration supplies illustrative Phase 6e values, not a new
parameter database.

An active configuration uses `mode='tissue_specific'`, `TE_s`, `TR_s`,
`water={tissue:{T1_s,T2_s}}`, and
`metabolites={tissue:{name:{T1_s,T2_s}}}`. A no-attenuation control uses
`mode='none'` with TE/TR and no profiles. Omitted relaxation records R=1,
TE=.030 and unknown TR as null in the sidecar; it does not invent a TR header.
All applied factors and active input times are recorded.

Three mechanisms stay distinct: coherent PRESS evolution at TE30 already in
the basis, scalar pre-acquisition TE/TR attenuation, and acquisition-time line
shape. The frozen responses retain 1 Hz native readout decay and metabolite
processing effects. Optional additional broadening is
`E_a(t)=exp(-pi*additional_linewidth_hz[a]*t)`, independently for metabolites
and water; default 0 preserves the validated arrays exactly. It begins at 1,
preserves initial amplitudes and precedes noise injection. It is **additional**
linewidth; total measured peak FWHM is not inferred or forced equal. No tissue
T2 is silently reused as acquisition decay or applied twice.

## Spatial cells, geometry and B0

Inline spatial example (instead of `tissue_fractions`):

```json
"spatial": {
  "cell_labels": ["WM", "GM", "CSF"],
  "cell_volumes_mm3": [2, 3, 5],
  "b0_offsets_hz": [-3, 7, 1],
  "b0_units": "Hz",
  "b0_convention": "positive_frequency"
}
```

Equal-volume cells are default. `spatial_weights` must be finite, nonnegative
and already sum to one; alternatively positive cell volumes are normalized by
their total. These options are mutually exclusive. All arrays must match the
1D cell list or 3D ROI shape. Explicit weights describe represented-volume
quadrature; they are not silently interpreted as receive-coil sensitivity.

Positive-frequency offsets multiply both signals by exp(+i*2*pi*Hz*t).
`legacy_shim` declarations negate the supplied Hz map, matching the old negative
shim convention. Units/sign are never inferred from a map name or display axis.
Absent B0 means zero; no random shim map is generated in this mode.

For an actual phantom ROI use, for example,
`phantom_spatial={'use_phantom_b0':True,'b0_units':'Hz',
'b0_convention':'legacy_shim'}`. B0 loading is disabled by default, including
when the phantom has a B0 image. Label/B0 shapes, affines and spacing must agree;
no implicit registration/resampling is performed. Cell volume is
abs(det(anatomical_affine[:3,:3])). Export affine uses the geometric midpoint
of selected cell centers and actual covered extents. Both acquisitions share
that affine, spectral grid and B0 modulation. Detailed geometry/weights/B0
statistics remain in the JSON sidecars.

## Noise

`noise={'seed':70000,'metabolite_snr':50,'water_snr':1000}` uses the existing
Phase 7a acquisition noise implementation. Alternatively specify independent
`metabolite_sigma` and `water_sigma` in per-real/imaginary-component FID units.
Sigma and SNR are mutually exclusive for each acquisition. A missing arm remains
noiseless. Target metabolite SNR requires a nonzero clean NAA component.

The local PCG64 streams use separate SeedSequence children for metabolite and
water. Noise is added once **after** tissue/spatial integration, relaxation,
B0 modulation and optional broadening. Clean FIDs and components are preserved.
Target/realized SNR, sigma, child streams and seeds are recorded. In headless
batches the per-voxel seed is the configured root seed plus voxel index, with
no global NumPy/Torch RNG change. Keep batch ordering fixed for repeatability.
The legacy `seed=` option is rejected for quantitative mode; use `noise.seed`.

## Paired export

Each output voxel directory contains

```
metabolites.nii.gz
water_reference.nii.gz
ground_truth.json
simulation_config.json
calibration_provenance.json
```

When noise is configured, `metabolites_clean.nii.gz` and
`water_reference_clean.nii.gz` are also written. Files use the validated
`gen_nifti_mrs` conjugation convention, preserving full complex FIDs on reload;
there is no spectral-magnitude-only conversion or gain normalization. Standard
headers record dwell time, frequency, chemical-shift reference, TE and explicit
TR when provided. Calibration/custom truth stay in JSON sidecars. The water
file is the simulated unsuppressed acquisition, never the residual-water component.

## Examples and executed validation

Executed with Pixi, with no fitting:

```sh
pixi run python --version
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.quantitative_water_reference_production outputs/quantitative_water_reference_production
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m simulation.headless --config outputs/quantitative_water_reference_production/pure_WM.json --output-dir outputs/quantitative_water_reference_production_cli
```

The example writes runnable JSON configurations, four paired datasets and an
independent validation_results.json. Each case is beneath its named folder and
`voxel_0000_0_0_0/`:

| Case | Clean initial water amplitude | Independent metabolite max error | Independent water max error | Export/reload max error |
|---|---:|---:|---:|---:|
| pure_WM | .007887908291571334 | 0 | 0 | 0 |
| mixed (.2/.3/.5) | .010052314421648305 | 9.47e-22 | 1.79e-18 | 0 |
| mixed + relaxation/B0 | .004407053763208692 | 9.47e-22 | 1.23e-18 | 0 |
| noisy WM (clean truth) | .007887908291571334 | 0 | 0 | 0 |

Unit tests also cover pure GM and CSF, independent tissue/cell complex sums,
Phase 6b/6e/6f exact adapter regressions, Phase 7a seeded noise, unequal additional
linewidths, concentration/gain linearity, water selection independence, metadata
rejection, actual source hashes, no RNG mutation, CLI execution and paired export.
Fitting entry points are forbidden in the production tests. Tests involving
LCModel/FSL fitting experiments are intentionally not executed in Phase 8.

Legacy regression baseline is the exact pre-Phase-8 revision
`c85821a833209b45384c0a8835e5d2ead8051daf`; production files were clean and matched
that revision before editing. A small synthetic anatomy/two-resonance fixture
executes its original headless and phantom code against the current default and
explicit legacy modes, using controlled NumPy/Torch seeds. It runs real legacy
metabolite concentration sampling, T2 behavior, shim, residual water and noise.
Only external MM/lipid providers are substituted with fixed nonzero spectra.
Complex total and all component arrays, sampled/voxel concentrations, geometry
and metadata, time/ppm and restored RNG states match **exactly (zero error)**.
This is a meaningful small-path regression, not a claim of external MATLAB or
full anatomical dataset equivalence. The actual SignalModel source is also
verified byte-identical to that revision.

Validation completed with **187 passed** in the selected non-fitting regression
suite (107.24 s), **25 passed** in the production-focused run (35.18 s), and
**4 passed** in the final targeted checks after strengthening saved Phase 7a
noise, gain/linearity, anatomical export and legacy-mode verification (12.50 s).
There were no failures or skips. Warnings are existing Dash/SWIG/tqdm and the
unchanged legacy NumPy/Torch shim deprecations.

Exact test commands:

```sh
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_quantitative_production.py tests/test_quantitative_production_legacy.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q \
  tests/test_quantitative_production.py tests/test_quantitative_production_legacy.py \
  tests/test_analytical_water_reference.py tests/test_analytical_voxel_reference.py \
  tests/test_frozen_basis_normalization.py tests/test_mrscloud_water_calibration.py \
  tests/test_mrscloud_full_basis.py tests/test_mrscloud_quantitative_pair.py \
  tests/test_mrscloud_quantitative_voxel.py tests/test_mrscloud_quantitative_spatial.py \
  tests/test_quantitative_tissue_relaxation.py tests/test_quantitative_acquisition_noise.py \
  tests/test_nifti_mrs_roundtrip.py tests/integration/test_quantitative_spatial_bridge.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q \
  tests/test_quantitative_production.py::test_independent_linewidth_and_water_selection_linearity \
  tests/test_quantitative_production.py::test_noise_phase7a_rng_and_clean_export \
  tests/test_quantitative_production_legacy.py
```

The final noise check additionally reproduces both saved Phase 7a SNR50,
realization00 FIDs exactly, using its actual seed70100, without refitting.
The full unfiltered repository suite was deliberately not run because it
executes LCModel/FSL-MRS fitting experiments, which are excluded in this phase.
Reading/formatting a basis and reloading NIfTI-MRS use installed FSL-MRS I/O
only. No fitting outputs or concentration estimates are validation targets.

Files modified: `simulation/headless.py`, `simulation/mrs_phantom.py`.
Files added: `simulation/quantitative.py`,
`examples/quantitative_water_reference_production.py`,
`tests/test_quantitative_production.py`,
`tests/test_quantitative_production_legacy.py`, this document and dedicated
production example outputs. Legacy `SignalModel`, GUI, residual water,
analytical reference and calibration files were not edited. Existing Pixi
manifest/lock changes were not modified by this phase.

## Limitations and review boundary

This first production mode supports paired single-voxel/ROI acquisitions and
batches of such acquisitions. It does not add a full spatial imaging pipeline,
GUI controls, quantitative MM/lipids or residual-water contamination. Existing
legacy functionality remains available through the default mode.

Common sequence response, proton scale, tissue density and declared attenuation
provide an approximately physical relative signal framework, but excitation and
localization approximations, asymmetric metabolite basis processing, relaxation
assumptions, real acquisition filtering/coil sensitivity and scanner receiver
calibration remain unresolved. Tissue T1/T2 values must be chosen deliberately.
No fitting-derived amplitude correction or toolbox concentration target enters
this implementation. Different toolboxes may retain their documented integration
and quantification biases. No new LCModel/FSL-MRS fits were performed.

Bundle loading is deliberately strict about source and formatting-version
identity. Moving the bundle without the recorded repository sources, changing
FSL/NumPy formatting versions or choosing a new grid/sequence requires a newly
reviewed calibration rather than silent reuse. Bundle validation currently
occurs per simulation; large batches may benefit from a separately reviewed
validated in-memory response cache. Legacy module imports can emit existing
Octave/tqdm startup messages; no sequence regeneration is requested by this mode.

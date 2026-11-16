# Phase 6f: normalized spatial integration and declared B0 offsets

## Scope and equations

`simulation/mrscloud_quantitative_spatial.py` is a separate, opt-in adapter.
The existing Phase 6b/6e tissue-fraction pathway is unchanged. It integrates
selected anatomical cells with normalized nonnegative weights, the Phase 6a
unit responses and shared receiver factor G:

\[
s_m(t)=G\sum_i w_i\sum_m10^{-3}C_{m,j(i)}R_{m,j(i)}B_m(t)
                 e^{i2\pi\Delta f_i t},
\]
\[
s_w(t)=G\sum_i w_iC_{w,j(i)}R_{w,j(i)}W(t)e^{i2\pi\Delta f_i t}.
\]

Metabolites are nominal mM and water mol/L. R is exactly one by default, or
the Phase 6e pre-acquisition saturation/echo factor with explicit second-valued
T1/T2 and actual sequence TE. B0 offsets are declared Hz; time is n/bandwidth,
starting at zero. No additional proton, water gain or amplitude calibration
factor is applied. The same cell frequency modulation is used in both signals.

This is volume-average integration, not an anatomical-volume-summed signal.
The common scalar remains **8857.353925964408**, with calibration identity:

```
mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81
```

No legacy SignalModel, GUI, MRSCloud source, calibration bundle, raw reference
or earlier quantitative implementation was changed. Only a new module,
example, tests and this document were added.

## Audit of existing interfaces and signs

| Interface | Audited behavior |
|---|---|
| `DigitalPhantom.get_label_data()` (`mrs_phantom.py:101`) | CPU tensor, 3D integer label map; one-hot form is optional |
| `label_names` | Explicit index-to-name mapping, typically Background/WM/GM/CSF plus other classes |
| `simulate_data` (`:130`) | Six start/stop indices; stop bounds exclusive |
| `spacing` / `affine` | Project geometry convention is mm, voxel-index-to-world affine |
| `get_B0_data()` (`:117`) | Returns `subject['B0'].data.numpy()[0]`; no unit or sign validation |
| `utils/loading.py:107`, `load_custom_files` | Loads `B0_brain.nii.gz` as a ScalarImage; it does not establish Hz provenance or convert units |
| `SignalModel.simulate_mrs_data` (`signalModel.py:50`) | Uses label means for fractions, multiplies per-cell spectra by product(spacing), then sums cells |
| `SignalModel.apply_shim_map` (`:383`) | Converts spectra to FIDs and multiplies by **exp(-i*2*pi*map*t)** |
| `utils.auxillary.fid2spec` | Forward FFT followed by fftshift; positive exp(+i*2*pi*f*t) falls at positive FFT frequency |

The new sign is deliberately **positive**. A legacy shim value +37 therefore
corresponds to new offset -37 when reproducing the legacy modulation. An
executed test calls the actual legacy shim method and verifies its negative
sign. Its float32 shim-map intermediate introduces up to about 4e-7 relative
complex error against an ideal double-precision negative modulation in that
control; the legacy test allows 5e-7. No legacy arithmetic was changed.

A plot axis labeled Hz does not establish a source map's physical units.
Automatic loading is therefore disabled. The caller must declare both units
and the source sign convention, based on its actual provenance.

## API and weights

```python
import numpy as np
from examples.mrscloud_quantitative_voxel import BUNDLE, METABOLITES, WATER
from simulation.mrscloud_quantitative_spatial import simulate_spatial
from simulation.mrscloud_quantitative_pair import export_pair

pair = simulate_spatial(
    np.array(['WM', 'GM', 'CSF']),
    spatial_weights=np.array([.2, .3, .5]),
    b0_offsets_hz=np.array([-3., 7., 1.]),
    calibration_bundle=BUNDLE,
    metabolite_concentrations_mM=METABOLITES,
    water_concentrations_mol_l=WATER,
    receiver_factor=1.)
export_pair(pair, 'new_spatial_voxel')
```

Inputs are a nonempty 1D cell list or 3D selected ROI, with explicit named
WM/GM/CSF labels. Unknown labels are rejected even when their weight is zero.
B0/weights/volume arrays must exactly match the label array shape; scalar
broadcasting and implicit resampling are not provided. B0 values must be finite.
Missing B0 means exactly zero, with no random shim generation.

Three weighting modes:

- **Equal cells:** default 1/N. Absolute anatomical volume is unknown unless
  supplied by the phantom bridge. This assumes equal represented cell volume.
- **Explicit normalized weights:** finite, nonnegative and sum to one with
  absolute tolerance 1e-8. Accepted summation roundoff is divided out and
  recorded; arbitrary unnormalized weights are rejected.
- **Cell volumes:** finite strictly positive `cell_volumes_mm3`, divided by their
  finite total. This mode is mutually exclusive with explicit weights.

The normalized weights define represented-volume quadrature, not an additional
coil sensitivity model. Tissue fractions are sum(weights for that tissue), and
nominal voxel-average concentrations remain sum(f*C) under that declared
quadrature. Counts and count-based fractions are recorded separately: they
need not equal weighted tissue fractions. Explicit weights are not automatically
claimed to equal geometrical volume fractions. Use physical cell-volume weights
when the ground truth must represent known geometric volumes. In the phantom
bridge, geometric volume remains separately recorded when explicit weights
replace uniform volume quadrature.

Scaling every cell volume by ten leaves both signals unchanged, while changing
recorded absolute volume. No volume multiplier is introduced. WM/GM-only
concentrations retain the existing exclusion of CSF and are null for pure CSF.
The inherited nominal concentration fields are not replaced by attenuated
values. Relaxation-weighted effective concentrations are recorded separately,
including the R=1 case.

`relaxation=...` accepts exactly the Phase 6e configuration, and all existing
basis/sampling compatibility checks are reused. TE remains the real PRESS TE30;
2048 points, 2000 Hz and 127.731 MHz are inherited from the validated bundle.
The water unit response retains its raw pre-DC processing convention.

Identical offsets are grouped per tissue; the implementation sums three tissue
frequency envelopes without allocating an anatomical-cells-by-time-points
matrix. It loads the validated unit responses separately from the existing
voxel adapter rather than changing that established API. Cost still grows
with the number of distinct cell offsets; this is an ROI implementation.

## Phantom bridge and geometry

```python
from simulation.mrscloud_quantitative_spatial import simulate_phantom_spatial

pair = simulate_phantom_spatial(
    phantom, (x0, x1, y0, y1, z0, z1),
    use_phantom_b0=True,
    b0_units='Hz',
    b0_convention='positive_frequency',
    calibration_bundle=BUNDLE,
    metabolite_concentrations_mM=METABOLITES,
    water_concentrations_mol_l=WATER)
```

The bridge uses actual label/B0 getters and never invokes `simulate_data`,
legacy SignalModel, shim generation or external anatomical-data constructors.
Its default does not read the B0 map, even if one exists. Loading requires
`use_phantom_b0=True`, `b0_units='Hz'`, and one of:

- `positive_frequency`: values are the offsets in exp(+i*2*pi*offset*t).
- `legacy_shim`: negate the values to reproduce the legacy negative modulation.

No ppm/Tesla conversion or unknown sign inference is attempted. The unit/sign
declaration and applied conversion are recorded. Cropping copies the map;
negation does not modify the original B0 array. A missing map, wrong dimensions,
nonfinite selected offsets or mismatched B0/label affines is rejected.

Bounds must be integral, nonempty and inside the 3D map. Background and unknown
labels are rejected. Spacing must be finite positive 3D values matching affine
column lengths; the affine must be finite, homogeneous and have nonzero volume.
The phantom's affine must also match its label image. Cell volume is
abs(det(affine[:3,:3])) mm^3, which agrees with product(spacing) for orthogonal
cells and remains well-defined for a nonsingular sheared grid. No implicit
B0 registration/resampling occurs. ROI shape, C-order cell ordering, original
bounds, affine, spacing and absolute geometric volume are preserved.

An actual TorchIO/DigitalPhantom test ROI contains 4 WM, 3 GM and 1 CSF cells,
fractions .5/.375/.125, with 2x3x4 mm cells and volume 192 mm^3. It verifies
copied B0 offsets, explicit legacy-sign conversion, geometry and rejection of
unregistered maps. No large external dataset is needed.

## B0 interpretation and metadata

Metadata records cell labels/counts, normalized weights and divisor, weighting
mode, optional volumes, B0 offsets, Hz units, sign formula, weighted mean/std,
range, calibration, receiver gain, nominal and relaxed effective concentrations.
All per-cell lists follow flattened C-order ROI indices. B0 distribution
statistics use spatial weights, not metabolite- or water-concentration weights.

A homogeneous tissue with common concentration/relaxation has one shared
whole-voxel modulation envelope for metabolites and water. In a heterogeneous
ROI, different concentrations can weight the same field distribution differently.
Their normalized whole-voxel envelopes can therefore differ without violating
the shared **cell** modulation equation. A test explicitly covers this case.

B0 does not change nominal tissue concentrations or the t=0 attenuation-weighted
coefficients. It changes later complex samples. A symmetric two-point field
distribution produces cos(2*pi*d*t), including sign reversals and revivals;
it is not a monotonic Lorentzian decay or a newly specified linewidth.

The acquisition reference stays 4.68 ppm. The water config's +2.55462 Hz is
explicitly marked as its intrinsic pre-B0 offset; the cell distribution supplies
the additional shifts. Constant +37 Hz moves the entire water FID to
+39.55462 Hz relative to that origin. No metadata relabeling or phase alignment
is used to conceal the shift. Generic optional aligned-water export removes
its declared intrinsic offset only; it does not undo cell-wise B0 dispersion.

## Executed numerical controls

Outputs: [validated directory](../outputs/mrscloud_quantitative_spatial/validated/).
Each case contains two actual NIfTI-MRS files, both sidecars and an NPZ with
actual/independently expected complex arrays and time axis. The root contains
bundle hashes and [validation_results.json](../outputs/mrscloud_quantitative_spatial/validated/validation_results.json).
All source calibration hashes remain unchanged.

| Control | Tissue / offsets | Metabolite cell-sum error | Water cell-sum error | Water magnitude at 50 ms |
|---|---|---:|---:|---:|
| Zero B0 | WM; 0/0 Hz | 1.895e-21 | 0 | .006741290384 |
| Constant B0 | WM; +37/+37 Hz | 1.747e-21 | 1.940e-18 | .006741290384 |
| Symmetric B0 | WM; -5/+5 Hz | 1.747e-21 | 1.749e-18 | 4.128e-19 |
| Mixed + relaxation | .2 WM/.3 GM/.5 CSF; -3/+7/+1 Hz | 9.471e-22 | 1.227e-18 | .001620981665 |

Expected signals are independently summed per cell, not calculated from returned
tissue envelopes or simulated outputs. The mixed example retains Phase 6b
nominal NAA/Cr/GPC/PCh/Glu = 6.2/4.3/.85/.56/5.2 mM and water 49.72245 M;
it uses the unchanged Phase 6e illustrative relaxation profiles.

Across pure WM/GM/CSF and mixed voxels, with and without relaxation (eight
zero-field comparisons), maximum metabolite difference from Phase 6b/6e is
**1.895e-21**, water difference **zero**. The difference is floating-point
summation/multiplication order, not a fitted scale. Constant +37 Hz agrees with
baseline*exp(+i*2*pi*37*t) within **1.940e-18** full-array error; symmetric
+/-5 Hz agrees with baseline*cos(2*pi*5*t) within **1.735e-18**. Its first
null is at 50 ms. Initial complex sample differences are exactly zero in all
four exported controls, with original water amplitude .007887908291571334
for the three homogeneous WM cases and .004407053763208692 for relaxed mixed.

Repeated simulations have exactly equal complex arrays and metadata. All paired
complex export/reloads have **zero** error. Tests inspect stored conjugate
samples and FSL-interpreted FIDs, dwell .0005 s, spectrometer frequency
127.731 MHz, acquisition origin and JSON sidecars. Separate repeated exports
have identical interpreted samples and byte-identical JSON sidecars; binary
NIfTI header timestamps are not used as a determinism assertion.

## Tests, execution and scope limits

Executed from the repository root, using Pixi only:

```sh
pixi run python --version
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.mrscloud_quantitative_spatial outputs/mrscloud_quantitative_spatial/validated
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_mrscloud_quantitative_spatial.py tests/integration/test_quantitative_spatial_bridge.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q
```

Focused: **37 passed** in 25.81 seconds. Initial test runs exposed overly strict
control assertions: one used multiplication rather than the declared n/bandwidth
time grid, producing phase roundoff up to 5.7e-14; the other assumed the legacy
float32 shim arithmetic was double precision. The time grid was corrected and
the legacy sign test now uses its justified float32 tolerance and an explicit
negative-frequency assertion. Production signals were not changed to satisfy
either expectation. New double-precision analytical comparisons use 2–3e-15
relative tolerance plus 1e-20 absolute for summation/cancellation roundoff.

Full regression: **253 passed, 1 xfailed**, 4 dependency warnings in
161.40 seconds. No failures. The expected failure remains the previously
documented installed-package combined phase/frequency WaterRef initialization
defect. Warnings concern external Dash/SWIG deprecations. The full run includes
the final repeated-export assertions, earlier analytical/export/relaxation tests
and previous actual FSL-MRS CLI fitting regressions.

Optional new FSL-MRS spatial fits and CLI reports were **not performed**. This
phase establishes analytical spatial formation and export, not fitting recovery
under arbitrary field-distribution model mismatch. Earlier FSL fitting/report
regressions run in the full suite do not establish recovery for these new B0
cases. No spatial fitted coefficients or water-referenced concentration claim
is made here. Such fits must compare source-derived predictions with actual
model mismatch rather than assuming a Voigt can represent a cosine envelope.

B0 modulation acts during acquisition only: it does not rerun localization/RF
pulse response for each shifted cell, introduce TE phase, change the frozen
basis, or model susceptibility-dependent sequence efficiency. The intrinsic
readout line shape and phenomenological pre-acquisition relaxation remain
separate. Cells are piecewise homogeneous, no soft labels, coil sensitivity,
noise, B0 interpolation, dynamic field drift or full-brain performance claim is
included. Caller-declared map units/sign and sensible sampled frequency ranges
remain essential; frequencies outside the sampled band can alias. Existing
sequence-response/integration calibration limitations remain. This is internal
synthetic consistency, not physical scanner calibration. Stop after Phase 6f;
no LCModel, GUI or subsequent production integration was attempted.

# Approximate LCModel BASIS-derived water reference

**This mode provides an approximate, basis-derived water reference. It is not
equivalent to a separately simulated, water-centred unsuppressed acquisition and
does not establish scanner-level absolute calibration.**

## Supported production interface

Use the existing headless API with the explicit `water_reference_mode: basis_h2o`.
Legacy remains the default; `sequence_matched_quantitative` retains its frozen
PRESS TE30 calibration and behavior. The two modes cannot be combined.

```python
from simulation.headless import run_headless

config = {
    'water_reference_mode': 'basis_h2o',
    'path2basis': 'path/to/supplied.BASIS',
    'skeleton': 'MRiLab',
    'positions': [[100, 108, 80]],
    'voxel_size_mm': [10, 10, 10],
    'quantitative': {
        'basis_amplitude_convention': 'assume_common_scale_protons_encoded',
        'concentration_units': 'mM',
        'frequency_origin_ppm': 4.68,
        'metabolite_concentrations_mM': {
            'WM': {'NAA': 10, 'Cr': 8},
            'GM': {'NAA': 14, 'Cr': 9},
            'CSF': {'NAA': 0, 'Cr': 0}},
        'bandwidth_hz': 2000, 'points': 2048, 'receiver_factor': 1},
    'output_dir': 'outputs/my_basis_h2o_voxel'}
pair = run_headless(config)[0]
```

A supplied phantom may be passed as `run_headless(config, phantom=phantom)`.
Without one, the existing anatomy loader reads `skeleton`; DigitalPhantom's new
explicit mode skips legacy basis normalization, stochastic concentrations, lipids
and residual-water initialization. Geometry uses the existing exclusive-stop ROI
bridge and normalized represented-volume weighting, not absolute volume gain.
Positions are zero-based anatomical cell-center indices. Background/unsupported
tissue within the selected ROI is rejected rather than renormalized away.

For anatomy-free simulations, omit `skeleton`/`voxel_size_mm` and supply
`tissue_fractions` or `spatial` inside `quantitative`. The direct adapter is
`simulation.basis_h2o.simulate_basis_h2o(config, phantom=None, voi_bounds=None)`;
it takes `basis_path` instead of the headless top-level `path2basis`.

A runnable five-metabolite example, with explicit illustrative tissue profiles:

```bash
export MRS_BASIS_H2O=/path/to/GE_PRESS_TE35.BASIS
pixi run python -m examples.basis_h2o_headless \
  --basis "$MRS_BASIS_H2O" --position 100 108 80 --voxel-size 10 10 10 \
  --output outputs/my_basis_h2o_roi
# Small smoke run, independent of anatomical datasets:
pixi run python -m examples.basis_h2o_headless \
  --basis "$MRS_BASIS_H2O" --fractions 1 0 0 --output outputs/my_basis_h2o_WM
```

JSON configurations also run through the existing CLI:
`pixi run python -m simulation.headless --config CONFIG.json --output-dir NEW_OUTPUT`.
No Phase 6a calibration bundle, Octave, LCModel executable or fitting routine is
required. Output voxel directories must be new; existing acquisitions are preserved.

## Compatibility and explicit assumptions

The loader requires an unencrypted ASCII LCModel BASIS readable by FSL-MRS, unique
component names, exactly one SEQPAR/BASIS1 header, positive consistent HZPPPM,
BADELT, NDATAB and ECHOT, and explicit component CONC/TRAMP/VOLUME=1 plus integer
ISHIFT. The exact water component defaults to `H2O`; `water_component_name` can
name another explicitly identified water entry. It is excluded from metabolites,
as is Ref0ppm. Names are not translated or silently merged. Metabolite profiles
must explicitly select exact names from the file. Control/correction entries such
as CrCH2 should not be assigned physiological metabolite concentrations.

Missing/zero/nonfinite water, nonfinite samples, duplicate names, inconsistent
sampling, unknown components and grids extending beyond native bandwidth or
acquisition duration are rejected. Nonunit amplitude headers are deliberately
unsupported because the existing reader does not establish their molar meaning.
Known independently normalized bases without a recoverable common scale must not
be used. `basis_amplitude_convention` must explicitly equal
`assume_common_scale_protons_encoded`; declarations of independent normalization
or missing assumptions are rejected. A BASIS spectrum alone cannot prove or rule
out hidden normalization performed before export. For unknown provenance this
option is an explicit approximation, not a detector or repair of that history.

`frequency_origin_ppm` is required: the supplied file does not establish the
acquisition reference by itself. It labels the existing frame and performs no
frequency shift. `expected_basis_sha256` and `expected_calibration_identity` are
optional assertions. Every load computes and records the actual SHA256 and validates
sampling. No mutable global calibration cache or automatic regeneration is used.

## Normalization and signal equations

The fixed reference is **every component in the supplied BASIS**, including water,
in reader-sorted order. Existing global scaling is preserved:

```
F_all = FSL_BASIS.get_formatted_basis(bandwidth, points, scale_factor=None)
a_ref = 1
B_m = F_all[:, m] / a_ref
B_H2O = F_all[:, H2O] / a_ref
```

The statistic is `identity_preserve_supplied_common_scale`: no additional estimated
normalizer or amplitude calibration is imposed. All components are formatted before
selection. No selection-dependent mean-magnitude normalization, per-component norm,
water-only normalization or proton multiplier is applied. The file hash, full fixed
reference list, grid, declared reference, assumptions and formatter versions identify
the `basis-h2o-approximate-v1` calibration. This identity differs from Phase 6a.

For tissue `j(i)` and normalized spatial weights:

```
s_m(t) = G * sum_i w_i * sum_m 0.001*C_m,j(i)*R_m,j(i)*B_m(t)*exp(+i*2*pi*df_i*t)
s_w(t) = G * sum_i w_i * C_w,j(i)*R_w,j(i)*B_H2O(t)*exp(+i*2*pi*df_i*t)
```

Metabolites are nominal mM per total volume of their tissue compartment; water is
mol/L per total tissue volume. The 0.001 converts mM to mol/L. Water defaults remain
WM 39.0165, GM 47.2305 and CSF 55.5 mol/L, with explicit overrides supported.
Legacy CSV Conc_mean may mix IU/mM and is **not** implicitly reused as molar data.
Explicit `concentration_units: mM` and complete WM/GM/CSF profiles are mandatory.
There is no undocumented IU conversion or sampling of legacy concentration SDs.

Ground truth retains nominal voxel averages, WM/GM-only values (null for pure CSF),
effective water concentration and relaxation-weighted signal concentrations. These
are distinct conventions; no fitting output is used as ground truth.

## Processing, relaxation, B0 and noise

The actual installed FSL-MRS reader (`mrs_io/lcm_io.py::readLCModelBasis`) reads
spectra, applies stored ISHIFT, IFFTs and conjugates into its FID convention. Native
BASIS phase/area/first-point processing is retained. The entire basis is formatted
using the same FSL Fourier resampling for the requested grid, without scaling.
Resampling can change early samples and apparent linewidth, including water;
this is recorded as a formatting limitation, not repaired with empirical factors.
No first-point replacement, DC removal, extra factor of two or zero-phase/frequency
alignment is provided. No RF carrier is changed by this mode.

Existing quantitative utilities are shared through a private response-loader seam:
voxel/ROI tissue weighting, spatial quadrature, relaxation, receiver scaling,
additional linewidth, noise and paired export retain their equations. Legacy
`Basis.reformat()` and `SignalModel` are unmodified.

Optional `relaxation` uses the actual BASIS TE and explicit TR/T1/T2 in seconds:
`R=(1-exp(-TR/T1))*exp(-TE/T2)`. Active tissue relaxation additionally requires
`basis_relaxation_assumption: no_tissue_attenuation_encoded`, since arbitrary BASIS
provenance cannot exclude pre-encoded attenuation. TE mismatch is rejected. The
default has no additional attenuation. Sequence evolution and existing readout
linewidth are not removed. Independent `additional_linewidth_hz` values multiply
metabolites/water by their respective `exp(-pi*Hz*t)` envelopes, preserving t=0.

Spatial B0 uses the established explicit Hz/sign declarations and shared
`exp(+i*2*pi*df*t)` modulation. `phantom_spatial` is supported with explicit B0 units;
no map units/sign are inferred. Acquisition-level noise uses independent local
PCG64 child streams after deterministic integration/linewidth; clean acquisitions
remain available and NumPy/Torch global RNG state is preserved. The supported
noise/SNR schema is unchanged from the production user guide.

## GE PRESS TE35 inspection and provenance

The externally supplied `GE_PRESS_TE35.BASIS` was read without copying/modifying it.
SHA256: `f0a42cf0bef1b5427d4f3f3f4e5dd35dd4406581230a7ead1055bd7fbd700a79`.
It contains 25 entries, exact names:

```
Ala Asc Asp Cr CrCH2 Cystat EA GABA GPC GSH Gln Glu Gly H2O
Lac NAA NAAG PCh PCr PE Ser Tau bHG mI sI
```

Headers report PRESS, TE=35 ms, IDBASI='GE 35 ms', 127.731 MHz, native 4000 Hz,
8192 complex points, CONC=TRAMP=VOLUME=1 and ISHIFT=0 throughout. FWHMBA=0.007829 ppm
suggests approximately 1 Hz nominal basis broadening, not a measured water FWHM.
The H2O native FFT maximum is +2.44140625 Hz, or 4.69911365 ppm using origin 4.68.
A later-sample phase-increment estimate is about +2.55 Hz, consistent with 4.70 ppm;
it is not an RF-carrier estimate. NAA's acetyl maximum is near 2.008 ppm in that frame.

A read-only companion MRSCloud assembly/configuration in the supplied example data
records GE unedited PRESS TE35, RF carrier 3 ppm, 101 spatial points and one shared
assembly scalar **8964.19101422405**. The provided BASIS FIDs are exactly equal to
those read from the original companion LCModel export. All 25 columns match the
companion's already-normalized FIDs after the **single known** acquisition-frame
shift `(3-4.68)*127.731 = -214.58808 Hz`, with maximum relative error
**1.8991608464e-7**, consistent with the six-digit text spectra. No amplitude fit or
per-metabolite correction was applied to make this comparison.

This supports preservation of a shared assembly/export scale; it does not prove
absolute molar response or the absence of every upstream amplitude operation.
The supplied file itself carries no spin-system definitions. Bundled MRSCloud's
water spin system encodes two protons, but the exact external generating source
version cannot be authenticated from CONC=1 or this file alone. Encoded proton
response is therefore explicitly assumed for this mode; it is never multiplied
again. Classification: **compatible under the declared approximation**.

## Actual numerical validation (2000 Hz, 2048 points)

| Quantity | Value |
| --- | --- |
| Native H2O initial response | `9.959463211048405e-5 + 1.7063686777824097e-5j` |
| Native H2O initial phase | 0.1696838656350743 rad |
| Formatted H2O initial response | `9.964329578704734e-5 + 8.74990424560547e-6j` |
| Formatted NAA initial response | `2.5209703671603886e-4 - 3.712029818701174e-5j` |
| Pure WM water initial signal, 39.0165 M | `0.0038877326500753326 + 0.00034139063899866584j` |
| Pure WM NAA(10 mM)/water complex initial ratio | `0.0006351602323884909 - 0.00015125544691463568j` |
| Independent ratio error | 0 |
| Ratio difference after shared gain x3 | 1.0842021724855044e-19 |
| NAA selection-change full complex error | 0 |

Expected arrays were independently calculated from the full reader-formatted unit
responses and explicit tissue/relaxation/B0/linewidth sums. Across pure WM/GM/CSF,
mixed, relaxation+B0 and noisy cases, maximum clean-signal differences were
**1.1978854818e-21** for metabolites and **1.7881120307e-18** for water (summation
roundoff). All actual exported FIDs reloaded with **exact complex equality**;
quotient-based phase diagnostics were <=4.41e-17 rad, floating-point division
roundoff on identical arrays. Dwell time, origin, spectrometer frequency and TE
are preserved. Water/metabolite concentration linearity and independent receiver
scaling pass. No ratio is forced to the sequence-matched or ideal-proton result.

The 20/30/50% mixture gives NAA/Cr/GPC/PCh/Glu voxel-average concentrations
6.2/4.3/0.85/0.56/5.2 mM and water 49.72245 M. With the declared example relaxation,
B0 and additional linewidth, water t=0 is
`0.002099073870268963 + 0.00018432444676015419j`. Nominal concentrations remain
unchanged; attenuation-weighted concentrations are separate metadata.

The actual MRiLab ROI at (100,108,80), size 10 mm per axis, produced WM/GM/CSF
fractions 0.194/0.806/0, NAA 13.224 mM and water 45.636984 M. This exercises the
anatomical headless workflow rather than only hand-entered fractions.

## Export and notebook

Each voxel produces `metabolites.nii.gz`, `water_reference.nii.gz`, `ground_truth.json`,
`simulation_config.json` and `calibration_provenance.json`; noisy runs also retain
clean acquisitions. The validated `gen_nifti_mrs` convention preserves the actual
water phase/frequency and relative scale. Custom provenance stays in JSON.
Metadata includes BASIS SHA256/name, exact H2O name, scalar/full reference set,
approximate status, sequence/grid, units, tissue water and RF/processing limitations.
The basis-file hash is verified; unknown construction-source provenance is not
misrepresented as authenticated source hashes.

`notebooks/headless_simulation.ipynb` now has the section **Approximate water
referencing from an LCModel BASIS H2O component**. Set `MRS_BASIS_H2O` in the Pixi
kernel environment. Its tagged cells run independently; legacy cells and their
existing outputs are retained. The separate sequence-matched notebook is unchanged.

```bash
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m examples.execute_quantitative_notebook \
  --notebook notebooks/headless_simulation.ipynb --tag basis_h2o_demo \
  --output outputs/basis_h2o_validation/notebook_execution
```

Actual new code cells 9–12 all passed, including the GE/real-anatomy run, spectra,
metadata and exact complex export/reload. Legacy notebook cells were not executed.
Generated acquisitions/plots are under ignored `outputs/basis_h2o_notebook/<timestamp>/`.

## Tests, commands and organization

Production tests create a small independent analytical LCModel BASIS fixture; no
user BASIS, external phantom dataset, proprietary executable or research output is
required. Optional GE audits live only under `research/validation/basis_h2o/`.
Nothing is staged or committed, and neither external BASIS nor MAT is redistributed.

```bash
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m pytest -q tests/test_basis_h2o.py tests/test_quantitative_production.py tests/test_quantitative_production_legacy.py tests/test_calibration_resources.py
# Initial focused run: 61 passed in 54.77 s (before two final hardening tests).
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m pytest -q tests/test_basis_h2o.py
# 27 passed in 7.96 s.
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m pytest -q tests --ignore=tests/integration
# Initial broad run: 294 passed in 113.61 s.
# Final run including both hardening tests: 296 passed, 16 warnings, 112.39 s.
# No failures, skips or expected failures.
```

Actual GE validation commands used the supplied file via `--basis`; to reproduce:

```bash
pixi run python -m research.validation.basis_h2o.validate_basis_h2o \
  --basis "$MRS_BASIS_H2O" --output outputs/basis_h2o_validation/new_audit \
  --source-mat "$MRS_BASIS_SOURCE_MAT"
```

The optional MAT argument is not needed for the signal tests. Results are preserved
in `outputs/basis_h2o_validation/ge_press_te35_with_provenance/validation_results.json`.
An initial direct script-path launch failed with `ModuleNotFoundError: examples`;
module execution from the repository root (shown above) resolved it, without any
dependency/environment changes. No fitting or Octave experiments ran. Full integration
folders containing fitting/regeneration tests were intentionally excluded.

See the Phase 10 production/research path manifests for proposed separate staging.
The scientific model and Phase 6a bundle are unchanged; the remaining uncertainty
is the supplied BASIS's upstream molar/proton/RF-response provenance. Do not repair
that uncertainty with empirical water or metabolite factors.

## Proposed Git deliverables (not executed)

Production paths are listed explicitly in `docs/commit_manifests/phase10_production.txt`
(20 files, including the prior requested notebook split and shared documentation).
The two optional read-only audit files are listed in
`docs/commit_manifests/phase10_research.txt`. These lists are disjoint, contain no
calibration assets/external BASIS/generated outputs, and leave older optional
research work untouched. The complete Phase 9 manifests were updated for ongoing
production/research classification.

```bash
git add --pathspec-from-file=docs/commit_manifests/phase10_production.txt
# Proposed message: feat: support approximate H2O water reference from LCModel BASIS

git add --pathspec-from-file=docs/commit_manifests/phase10_research.txt
# Optional message: research: validate supplied BASIS H2O relative response
```

No commands above have been executed. Legacy exact fixed-seed regression, frozen
calibration hashes, analytical reference behavior, prior saved deterministic
quantitative arrays, units, tissue mixing, relaxation, noise RNG isolation and
complex export all pass in the 296-test non-fitting suite. Installed Dash/Swig and
legacy Torch/NumPy deprecation warnings remain visible. The fitted-concentration
experiments/integration suite were deliberately not run. Phase 10 ends here.

# Quantitative simulator user guide

The supported opt-in production mode is `sequence_matched_quantitative`.
Legacy simulation remains the default. No fitting executable, manual, research
output, Octave regeneration or experimental notebook is needed for this mode.
FSL-MRS libraries provide basis formatting and NIfTI-MRS I/O; fitting is not invoked.

## Quick start

From the repository root, use the existing locked Pixi environment:

```bash
pixi install
pixi run python --version
pixi run python -m simulation.headless \
  --config examples/configs/quantitative/pure_WM.json \
  --output-dir outputs/my_quantitative_example
```

The output directory must be new. The CLI creates one `voxel_0000_0_0_0/`
subdirectory for the example position. Existing acquisitions are never overwritten.
To repeat, choose another output directory. Other runnable configurations in
`examples/configs/quantitative/` are `mixed_tissue.json`, `relaxation.json`,
`spatial_b0.json` and `noisy.json`. Their `output_dir` fields declare their defaults.

Python uses exactly the same dispatch:

```python
import json
from simulation.headless import run_headless

config = json.load(open('examples/configs/quantitative/mixed_tissue.json'))
config['output_dir'] = 'outputs/my_mixed_example'
pair = run_headless(config)[0]
```

Omit `output_dir` for in-memory generation. Do not pass legacy component selectors,
residual-water settings or the legacy batch seed to quantitative mode.

## Notebook

Open `notebooks/headless_simulation.ipynb` using a kernel backed by the repository
Pixi Python. The original legacy demonstration precedes a standalone section titled
**Sequence-matched quantitative simulation with an unsuppressed water reference**.
That section can run in a fresh kernel; it does not depend on the preceding cells.
It demonstrates pure WM, 20/30/50% WM/GM/CSF, and relaxation/B0/noise. All examples
call `run_headless`. Spectra retain their relative units with separate axes.
New figures and acquisitions go into a unique timestamped directory below
`outputs/notebook_quantitative_demo/`; new embedded outputs are not committed.

For reproducible section execution without adding notebook dependencies:

```bash
pixi run python -m examples.execute_quantitative_notebook \
  --output outputs/my_notebook_execution
```

This uses installed `jupyter_client` and `ipykernel`, runs only the new section and
records cell status. It does not claim full execution of the anatomical legacy
example, which needs the MRiLab data and legacy lipid simulation dependencies.

## Configuration schema

Top-level fields are `simulation_mode`, `quantitative`, optional `positions`
(default `[[0,0,0]]`) and `output_dir`. Anatomical ROIs additionally support the
existing `skeleton` and `voxel_size_mm` interface. Inline fractions/cells do not
use skeleton or anatomical voxel size; an explicit `voxel_affine_mm` may be supplied.

The `quantitative` object accepts:

| Field | Units / behavior |
| --- | --- |
| `metabolite_concentrations_mM` | Required WM, GM, CSF dictionaries; same metabolite names, explicit zeros; mM per total volume of that tissue |
| `water_concentrations_mol_l` | Optional WM, GM, CSF dictionary; mol/L per total tissue volume |
| `tissue_fractions` | Finite nonnegative volume fractions, sum to 1; omitted tissues mean zero |
| `spatial` | Alternative cell labels, normalized `spatial_weights` or `cell_volumes_mm3`, optional explicit B0 |
| `phantom_spatial` | Alternative anatomical ROI bridge; explicitly declare B0 units/sign when used |
| `receiver_factor` | Shared finite synthetic receiver gain, default 1 |
| `calibration_bundle` | Optional explicit bundle directory; otherwise packaged resource |
| `bandwidth_hz`, `points` | Fixed to 2000 Hz and 2048 complex samples in this production family |
| `relaxation` | Omitted means no pre-acquisition attenuation; otherwise explicit second-valued configuration below |
| `additional_linewidth_hz` | Both `metabolites` and `water`, nonnegative additional Lorentzian widths; defaults zero |
| `noise` | Optional explicit root `seed`; independent metabolite/water sigma or target SNR |
| `expected_calibration_identity` | Optional assertion; production family is always checked |
| `concentration_convention` | Only `nominal_per_total_tissue_volume` |

Unknown or incompatible keys are rejected. See the JSON examples for full shapes.

## Frozen calibration and provenance

Default lookup is `simulation.calibration_resources.default_calibration_bundle()`:
`calibration/mrscloud_press_te30_quantitative/`, included as Python package data.
Explicit overrides remain supported and undergo the same scientific checks.
Missing assets produce an actionable error, with no download, regeneration or
analytical fallback. There is no mutable global calibration cache.

The frozen convention is `mrscloud-press-te30-full-reference-max-v1`, scalar
**8857.353925964408**, identity
`mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81`.
All seven original Phase 6a files retain their original bytes. `provenance_sources.zip`
adds 408 immutable construction-source snapshots checked against the original
metadata hashes. They are never imported or executed. An explicit historical bundle
without this archive can use the packaged snapshot only if its identity and complete
source hash maps match. This shares provenance, never response arrays or signal gain.
Checks include reference arrays/statistic, sequence, sampling, basis/water file
hashes, processing/proton convention and recorded NumPy/FSL-MRS formatting versions.

These checks intentionally reject incompatible library versions. Use the locked
Pixi environment; changing it may require a separately reviewed calibration family.
Ordinary Git stores this roughly 24 MiB asset set. No LFS or network fetch is needed.
One inert author-specific comparison path in the frozen manifest is preserved,
because editing it would change the frozen file hash; it is never used for lookup.

## Signal, water and ground truth

For normalized cell weights `w_i` and tissue `j(i)`:

```
s_m(t) = G * sum_i w_i * sum_m 0.001*C_m,j(i)*R_m,j(i)*B_m(t)*exp(+i*2*pi*df_i*t)
s_w(t) = G * sum_i w_i * C_w,j(i)*R_w,j(i)*W(t)*exp(+i*2*pi*df_i*t)
```

Metabolite input is mM; water is mol/L. Unit responses share the frozen sequence
scale and already encode proton response. No further proton or water gain factor
is applied. Water uses raw pre-DC PRESS TE30 response, water RF carrier 4.7 ppm;
metabolites use carrier 3 ppm and their validated basis DC/formatting pathway.
The acquisition reference is 4.68 ppm at 127.731 MHz. The unsuppressed water
reference is separate from legacy residual-water contamination, which is absent here.

Water defaults retain the existing tissue-volume convention:
WM `0.703*55.5 = 39.0165`, GM `0.851*55.5 = 47.2305`, CSF `55.5` mol/L.
They are not concentrations per volume of tissue water. Explicit overrides must
use the same convention. Mixed 20/30/50% effective water is **49.72245 mol/L**.

`ground_truth.json` distinguishes nominal voxel averages from WM/GM-only values
and relaxation-weighted signal concentrations. WM/GM-only divides WM/GM signal
contributions by `f_WM+f_GM`, excludes CSF contributions and is null for pure CSF.
Weights describe represented-volume averages; anatomical voxel volume is recorded
but is not an additional absolute signal multiplier.

## Relaxation and linewidth

Explicit `relaxation` uses `mode: tissue_specific`, actual `TE_s: 0.03`, a declared
`TR_s`, and WM/GM/CSF water `T1_s,T2_s` plus metabolite-specific values for every
selected metabolite. Times must be positive seconds. Pre-acquisition attenuation is
`R = (1-exp(-TR/T1))*exp(-TE/T2)`. `mode: none` preserves no attenuation.
The frozen PRESS sequence already includes its sequence evolution and 1 Hz
readout decay; these are not an additional tissue TE relaxation correction.

Additional acquisition linewidth multiplies by `exp(-pi*additional_Hz*t)` after
mixing. Independent water/metabolite values do not imply equal observed peak FWHM.
No first-point gain adjustment is made. Nominal concentrations remain distinct
from attenuated effective concentrations. Example T1/T2 values are explicit
illustrative assumptions, not universal constants.

## Spatial B0 and noise

Inline `spatial` accepts WM/GM/CSF cell labels with equal weights, explicit
normalized weights or positive cell-volume-derived weights. B0 must declare
`b0_units: Hz` and `b0_convention: positive_frequency` (positive means
`exp(+i*2*pi*df*t)`) or `legacy_shim` (converted with the opposite sign).
Unknown/background cells or ambiguous units are rejected. Metabolites and water
share cell modulation, but tissue concentrations and relaxation can weight it
differently. No automatic shim-map unit inference is made.

Noise is added once after deterministic integration and additional linewidth:
independent real/imaginary Gaussian FID components with SD `sigma`. Use
`metabolite_sigma`/`water_sigma`, or `metabolite_snr`/`water_snr` independently.
Target SNR uses the clean isolated NAA magnitude peak at 1.8–2.2 ppm (water at
4.1–5.3 ppm), divided by one-component spectral noise SD, with the orthonormal
FFT and no first-point weighting. This differs from toolbox-reported SNR.
A local PCG64 generator with separate SeedSequence child streams leaves NumPy
and Torch global RNG state unchanged. Batch seeds are root seed + voxel index.
Finite-sample realized SNR fluctuates; noise is never rescaled to force its target.

## Output interpretation

Each voxel directory contains:

- `metabolites.nii.gz`: water-suppressed complex acquisition.
- `water_reference.nii.gz`: separate unsuppressed complex water acquisition.
- `ground_truth.json`: nominal and effective concentrations, fractions, calibration.
- `simulation_config.json`: complete simulation, grid, relaxation, B0/noise settings.
- `calibration_provenance.json`: scientific identity, normalization and source hashes.

Noise-enabled output additionally retains `metabolites_clean.nii.gz` and
`water_reference_clean.nii.gz`. Export uses the validated 1D `gen_nifti_mrs`
complex convention. Both share dwell time 0.0005 s, frequency 127.731 MHz,
`SpecFreqChemShift=4.68`, and applicable actual TE/TR. The primary water retains
its phase/frequency. Custom scientific metadata remain in JSON sidecars.

## Production verification and limitations

```bash
pixi run python -m pytest -q tests/test_quantitative_production.py \
  tests/test_quantitative_production_legacy.py tests/test_calibration_resources.py
```

These tests require neither fitting binaries nor historical research outputs.
They cover the frozen bytes, relocation, tissue/relaxation/spatial/noise equations,
legacy fixed-seed regression, CLI, examples and exact complex export/reload.
Scientific validation history remains indexed in `research/README.md`; production
does not import it. See `repository_artifact_classification.md` for commit manifests.

This model preserves shared sequence response and explicit tissue assumptions;
it is not scanner-level absolute calibration. RF/sequence sensitivity, basis DC
processing, finite-grid integration and fitting-package quantification conventions
remain scientifically relevant. Generated ground truth is not a promise of unbiased
fitted concentration. No new fitting or scientific calibration occurs in Phase 9.

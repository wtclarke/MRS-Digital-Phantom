# Phase 6b: deterministic quantitative single-voxel adapter

`simulation/mrscloud_quantitative_voxel.py` provides an opt-in adapter for the
validated Phase 6a full basis and raw sequence-matched water. It accepts explicit
tissue concentrations and fractions, and can read a selected ROI from an
existing DigitalPhantom label map. It does not modify or invoke legacy
`SignalModel` simulation, and does not introduce spatial signal simulation,
noise, B0 effects, tissue relaxation or empirical amplitude factors.

## Existing voxel/concentration interfaces and units

| Existing interface | Representation and behavior |
|---|---|
| `DigitalPhantom.subject['labels']` / `get_label_data()` | TorchIO label map; a 3D integer-coded tensor after removing the image channel |
| `DigitalPhantom.label_names` | Index-to-tissue names; standard labels are Background/WM/GM/CSF, sometimes also Fat/Skull |
| `DigitalPhantom.spacing`, `affine` | Anatomical spacing in mm and 4x4 index-to-world affine |
| `DigitalPhantom.simulate_data(voi_coords)` | Six bounds `(x0,x1,y0,y1,z0,z1)`, with exclusive stops; crops one-hot labels and calls legacy simulation |
| `headless._position_and_bounds` / `simulate_batch` | Converts voxel-center indices and rounded mm/spacing counts to bounds; the adapter accepts bounds directly |
| `SignalModel.simulate_mrs_data` | Fractions are `one_hot.reshape(-1,n_labels).mean(0).real`; includes background if present |
| `SignalModel.simulate_metabs` | Samples `Conc_mean`/`Conc_std` by tissue using Torch normal draws, with basis amplitudes and legacy T2 processing |
| `sampled_concs_per_tissue` | Per-tissue sampled values, kept in the legacy table's existing units |
| `final_voxel_concentrations` | Sum of fraction times sampled tissue values; no cell-volume multiplier in concentration metadata |
| Legacy signal volume weighting | Multiplies cell spectra by `prod(voxel_spacing)` in mm^3 and sums anatomical cells; distinct from concentration averages |
| `utils.definitions.H2O_CONCENTRATIONS` | WM=.703*55.5=39.0165, GM=.851*55.5=47.2305, CSF=55.5, numerically mol/L; legacy residual water uses these with its own empirical/random amplitude model |
| Shim/B0 | Stored B0 map is passed through, or a random shim map is synthesized; `apply_shim_map` interprets values as Hz via exp(-i*2*pi*map*t). Loader provenance does not establish every input map's physical units. |

The legacy concentration preprocessing combines literature IU and mM columns
(`utils/preprocess_df.py`, around lines 294–297), then produces `Conc_mean`
through `calculate_metab_conc`. Those values cannot automatically be declared
explicit molarity. This adapter therefore does **not** consume or reinterpret
the legacy dataframe. The caller must provide deliberate nominal mM tissue
profiles. Legacy definitions, sampling, water defaults and volume weighting
remain unchanged.

There is no standalone reusable legacy fraction-only function: its mean is
embedded in a simulator that also draws random values and applies shim/noise.
The bridge reuses the actual label getter and the identical label-indicator
mean operation, with an independent comparison against the existing one-hot
getter in the tests. It never invokes that simulator to obtain fractions.

## API and tissue mixing

```python
from simulation.mrscloud_quantitative_voxel import simulate_voxel
from simulation.mrscloud_quantitative_pair import export_pair

bundle = 'outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative'
metabolites = {
    'WM': {'NAA': 10, 'Cr': 8},
    'GM': {'NAA': 14, 'Cr': 9},
    'CSF': {'NAA': 0, 'Cr': 0},
}
water = {'WM': 39.0165, 'GM': 47.2305, 'CSF': 55.5}
pair = simulate_voxel(
    {'WM': .2, 'GM': .3, 'CSF': .5}, calibration_bundle=bundle,
    metabolite_concentrations_mM=metabolites,
    water_concentrations_mol_l=water, receiver_factor=1)
export_pair(pair, 'new_quantitative_voxel')
```

All three tissue profiles are required, with identical nonempty metabolite
name sets and explicit zeros. Fractions may omit zero-volume tissues. Water
profiles are required explicitly; the example values are declared inputs,
not automatically imported assumptions. Unknown basis names, reference controls
as metabolites, missing profiles, negative/nonfinite values and incompatible
calibration metadata are rejected.

Fractions must be finite and nonnegative, with positive sum. Default sum
tolerance is absolute 1e-8 with no relative tolerance; accepted roundoff is
normalized and recorded. Arbitrary weights/counts require the explicit
`normalize_fractions=True` option. Input weights, input sum, divisor and policy
are preserved. No background tissue is silently removed and renormalized.

For unit responses B_m and W from the same Phase 6a calibration:

\[
\bar C_m=\sum_{j\in\{WM,GM,CSF\}} f_j C_{m,j}\ [\mathrm{mM}],\quad
\bar C_w=\sum_j f_j C_{w,j}\ [\mathrm{mol/L}],
\]

\[
s_m(t)=G\sum_m10^{-3}\bar C_m B_m(t),\qquad
s_w(t)=G\bar C_w W(t).
\]

The implementation averages explicit concentrations and delegates to the
unchanged Phase 6a synthesizer. Because responses are common across tissues
and no tissue attenuation is introduced, this is algebraically identical to
summing tissue-weighted complex FIDs. Tests independently calculate that
second formulation sample by sample. The receiver G is shared, no additional
proton multiplier is applied, and selection does not renormalize the basis.

There is **no anatomical volume multiplier**. These are volume-average signals
under the specified equations, not legacy volume-summed amplitudes. Anatomical
geometry is provenance only. A physical dependence on changed acquisition
voxel size/localization is not modeled by this phase.

## Ground-truth definitions

Each pair contains `ground_truth`, `simulation_config`, `components`, both
complex FIDs and a zero-origin `time_s` axis. The exported sidecars distinguish:

- `tissue_fractions` and `input_tissue_fractions` with normalization metadata.
- `tissue_specific_metabolite_concentrations_mM`: nominal concentration per
  total volume of that tissue compartment, not intracellular concentration.
- `tissue_specific_water_concentrations_mol_l`: the analogous bulk water input.
- `voxel_average_metabolite_concentrations_mM`: sum over WM, GM **and CSF**, per
  full selected voxel volume. Inherited `metabolite_concentrations_mM` is an
  alias of this voxel-average quantity for the existing pair export contract.
- `effective_voxel_water_concentration_mol_l`: weighted bulk water concentration;
  inherited `water_concentration_mol_l` has the same meaning.
- `tissue_only_metabolite_concentrations_mM`: specifically the WM/GM average
  \( (f_{WM}C_{m,WM}+f_{GM}C_{m,GM})/(f_{WM}+f_{GM}) \). It excludes CSF's
  contribution and is JSON `null` in pure CSF.
- Calibration identity, shared scalar, receiver factor, full sequence and
  sampling metadata, and `physical_calibration=false`.

The tissue-only definition is not automatically voxel-average/(1-f_CSF) when
CSF metabolites are nonzero. A test with CSF NAA=2 mM produces 7.2 mM voxel
average while preserving the correct 12.4 mM WM/GM-only average. Neither value
is a fitting estimate or an empirical CSF correction.

## Phantom label-map bridge

```python
from simulation.mrscloud_quantitative_voxel import simulate_phantom_voxel

# phantom is an existing DigitalPhantom; bounds use its usual index convention.
pair = simulate_phantom_voxel(
    phantom, (x0, x1, y0, y1, z0, z1), calibration_bundle=bundle,
    metabolite_concentrations_mM=metabolites,
    water_concentrations_mol_l=water)
```

The bridge reads `get_label_data()` on CPU, checks integral label codes and
nonempty in-bounds geometry, and computes equal-cell-volume WM/GM/CSF means.
Positive background/Fat/Skull or unknown labels are rejected. It records bounds,
anatomical spacing, cell count, volume estimate `count*prod(spacing)` and the
original anatomical affine. This assumes the phantom's uniform rectilinear
cell-volume convention; geometry is not used as an amplitude factor.

The integer getter is intentional: the existing one-hot getter chooses its
class count from the number of unique labels, which may be insufficient when
some indexed classes are absent. The bridge uses declared `label_names` and
indicator means instead, without changing that legacy method. It supports
hard-label maps, not soft-label images or fractional boundary resampling.

Tests construct a real TorchIO Subject and a DigitalPhantom instance without
calling its external-data constructor. A selected eight-cell ROI has four WM,
three GM and one CSF cell, giving **.5/.375/.125** exactly, matching the executed
legacy one-hot getter and mean. At 2x3x4 mm spacing its volume is **192 mm^3**.
Changing geometry to 1 mm cells changes metadata volume to 8 mm^3 but leaves
both concentration-based signals unchanged. Legacy simulation and B0 getters
are replaced by exceptions in the bridge test to prove they are not invoked.

## Calibration and paired export

The adapter calls the existing Phase 6a loader/synthesizer. It retains:

```
mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81
shared scalar = 8857.353925964408
```

Basis/source/reference/sampling compatibility is validated there; no frozen
scalar or unit response is recomputed from voxel fractions or concentrations.
`expected_calibration` can enforce exact metadata reuse, and mismatched point
counts/bandwidth are rejected. Phase 6a files are hash-checked before/after
the validation example and remain unchanged.

The unchanged `export_pair` produces:

```
quantitative_voxel/
    metabolites.nii.gz
    water_reference.nii.gz
    ground_truth.json
    simulation_config.json
```

It uses validated `gen_nifti_mrs` storage: raw file samples are conjugates of
the FSL-interpreted FID, and reader-side interpretation restores the original
complex signal exactly. Tests inspect both stored and interpreted arrays.
Headers retain dwell .0005 s, 2048 points, 127.731 MHz, 1H, TE .030 s and
`SpecFreqChemShift=4.68`. The primary water peak remains at 4.70 ppm with its
original phase and +2.55462 Hz offset. No alignment, independent gain or custom
calibration header field is introduced. Geometry stays in sidecars; these are
paired 1D datasets, not spatially registered MRS images.

## Executed analytical results

The reproducible example uses five metabolites with tissue profiles:

| Tissue | NAA mM | Cr mM | GPC mM | PCh mM | Glu mM | Water mol/L |
|---|---:|---:|---:|---:|---:|---:|
| WM | 10 | 8 | 2 | 1 | 8 | 39.0165 |
| GM | 14 | 9 | 1.5 | 1.2 | 12 | 47.2305 |
| CSF | 0 | 0 | 0 | 0 | 0 | 55.5 |

G=1, no noise or additional relaxation/B0 effects:

| Voxel | Effective water M | Voxel NAA mM | WM/GM-only NAA mM | Water initial amplitude |
|---|---:|---:|---:|---:|
| Pure WM | 39.0165 | 10 | 10 | .007887908291571334 |
| Pure GM | 47.2305 | 14 | 14 | .009548520563481090 |
| Pure CSF | 55.5 | 0 | null | .011220353188579424 |
| .2 WM / .3 GM / .5 CSF | 49.72245 | 6.2 | 12.4 | .010052314421648305 |

Mixed voxel-average concentrations are NAA 6.2, Cr 4.3, GPC .85, PCh .56 and
Glu 5.2 mM. WM/GM-only concentrations are NAA 12.4, Cr 8.6, GPC 1.7, PCh 1.12
and Glu 10.4 mM. The mixed NAA component begins
**1.7957503828180437e-6 - 2.1786916129965824e-7i**. Its water initial imaginary
part is approximately 6.02e-19, preserved from the raw sequence response.

Independent full complex tissue-sum errors across the four examples:

- Maximum metabolite absolute error: **3.4144989710325696e-21**.
- Maximum water absolute error: **1.788112030672749e-18**.
- Interpreted complex export/reload error: **0** for every signal.
- Wrapped phase differences: **0** for the nonzero exported signals; pure-CSF
  zero metabolite data have no physically defined phase.

Small independent-sum differences reflect floating-point operation ordering
between concentration averaging and summing FIDs. Deterministic repeats,
factor-of-two concentration/water/gain tests and subset comparisons pass exact
array assertions. No arbitrary amplitude factor is used to improve agreement.

## Regression evidence, commands and outputs

Executed commands:

```bash
pixi run python --version
pixi run python -m examples.mrscloud_quantitative_voxel outputs/mrscloud_quantitative_voxel/final
pixi run python -m pytest -q tests/test_mrscloud_quantitative_voxel.py tests/integration/test_quantitative_voxel_example.py
pixi run python -m pytest -q tests
```

Focused tests completed **30 passed**. The initial run had two erroneous test
assumptions: direct stored samples were compared against unconjugated FIDs,
and multiplying identical complex arrays introduced a tiny roundoff term in
the phase-error calculation. The checks were corrected to the established
storage convention and wrapped phase-angle difference. No simulation amplitude
or exporter was changed. Full regression execution completed with **185 passed,
1 xfailed, 5 dependency warnings in 106.87 s**. The expected failure remains
the installed FSL-MRS combined phase/frequency initialization defect. All
existing analytical, full-basis, processing, export, fitting and linewidth
regression suites executed.

The new fixed-seed legacy residual-water test executes the real legacy water
routine before/after the new adapter and verifies identical complex arrays and
subsequent NumPy random draws, with unchanged Torch RNG state. Existing legacy
entry points have no import/call of the new adapter. This is direct isolation
and residual-water regression evidence, not a claim of numerical equivalence
against an unavailable historical full-phantom baseline (see Phase 3 summary).

[Validation results](../outputs/mrscloud_quantitative_voxel/final/validation_results.json)
and [comparison CSV](../outputs/mrscloud_quantitative_voxel/final/comparison.csv)
record independent errors and hashes. Four subdirectories (`pure_WM`,
`pure_GM`, `pure_CSF`, `mixed`) contain actual paired NIfTI-MRS files and both
JSON sidecars. No external anatomical phantom was needed for the tests.

Files added: `simulation/mrscloud_quantitative_voxel.py`,
`examples/mrscloud_quantitative_voxel.py`,
`tests/test_mrscloud_quantitative_voxel.py`,
`tests/integration/test_quantitative_voxel_example.py`, and this document.
No legacy simulator, GUI, MRSCloud source, existing exporter, calibration bundle
or dependency file was changed.

## Limits and review boundary

The optional new FSL-MRS mixed-tissue fit was not performed. These results
validate tissue mixing and complex export, not fitted concentration recovery;
the response/integration, tissue-water and attenuation assumptions identified
in Phase 5h still require explicit treatment in any future quantification run.

Outputs are synthetic and have no demonstrated physical absolute calibration.
The bundle's fixed 1 Hz readout decay, ideal excitation, localization and
basis-processing conventions remain; no tissue-specific relaxation is added.
The bridge supports one hard-label ROI with uniform cell volumes and requires
explicit tissue concentration profiles. It does not sample literature values,
model physical voxel-size sensitivity, read shim/B0 fields, introduce noise,
perform full spatial simulation, integrate legacy SignalModel or change GUI.
Review these concentration and geometry conventions before extending the
pathway. No LCModel validation is performed in Phase 6b.

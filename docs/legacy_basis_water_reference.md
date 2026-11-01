# Semi-quantitative water in the legacy batch workflow

Use ordinary legacy settings and `simulate_batch(positions, settings, seed=..., save_dir=...)`.
Add:

```python
settings['water_reference_mode'] = 'basis_h2o'
settings['water_reference'] = {
    'concentration_convention': 'assume_table_values_nominal_mM',
    'basis_amplitude_convention': 'assume_common_scale_protons_encoded',
    'frequency_origin_ppm': 4.68,
}
```

`simulation_params` and `basis_set_settings` select this integrated legacy workflow.
The existing `quantitative` schema still selects the standalone adapter. Mixing
these schemas is rejected. Legacy defaults and the frozen sequence-matched mode
are unchanged. An optional `path2metabs` selects another concentration CSV; otherwise
the original table and its existing sampling/missing-value behavior are used.
No IU conversion is inferred: nominal mM is an explicit approximation.

All formatted BASIS columns, including H2O and controls, define one fixed
`mean(abs(full_formatted_basis))` scalar before metabolite selection. Both species
use that divisor. Enabled-mode signal levels can differ from the old selected-basis
normalization; changing selection cannot change the retained column amplitudes.
The default legacy normalization remains unchanged. H2O phase, RF attenuation,
DC processing and readout are retained, without an extra proton multiplier.

Legacy metabolites numerically use table values as nominal mM. Therefore water
is `1000*C_water_mol_l*R_water*B_H2O/scalar`. The factor 1000 is the mol/L-to-mM
conversion, not an empirical gain. Water uses the actual already-generated legacy
shim map, its negative-frequency sign, its dwell-shifted time vector, and the same
cell-volume sum. No shim or concentration draw is repeated. The existing
`shim_amplitude_hz`, correlation and boundary settings now work with this reference.
This preserves the legacy time-origin approximation rather than repairing it.

Optional `water_reference` fields:

- `water_concentrations_mol_l`: WM/GM/CSF, defaults 39.0165/47.2305/55.5.
- `additional_linewidth_hz`: default 0; additional Lorentzian decay with acquisition time starting at zero.
- `noise_sigma`: default 0; SD of each real/imaginary FID component after integration.
- `relaxation`: complete WM/GM/CSF dictionaries of positive `T1_s,T2_s`.
  Requires `basis_relaxation_assumption: no_tissue_attenuation_encoded`.
  Applies `(1-exp(-TR/T1))*exp(-TE/T2)` with BASIS TE and configured legacy TR.
  Legacy metabolite attenuation is retained, so water/metabolite corrections are
  explicitly asymmetric. No extra water relaxation is applied by default.

Water noise uses a local SeedSequence with the batch seed and voxel index. It
never changes legacy draws or concentration sampling. With seed=None its entropy
comes from the OS and is recorded. Residual contamination remains `components['Water']`;
the unsuppressed reference is never included in the total spectrum or component sums.

Results retain all legacy fields and add `water_reference_fid`,
`water_reference_spectrum`, `clean_water_reference_fid` and
`water_reference_provenance`. Ground truth additionally labels nominal voxel-average
and WM/GM-only metabolite concentrations and nominal water/factors. Pure-CSF
tissue-only concentration is null.

Saving preserves old filenames and adds `water_reference.nii.gz`, optional
`water_reference_clean.nii.gz`, `calibration_provenance.json`, `ground_truth.json`
and `simulation_config.json`. Both acquisitions carry the explicitly declared
frequency reference; setting this metadata does not shift the signal. The converter
also returns a separate `water_reference` NIfTI-MRS object, independently of the
legacy component selectors. Clean/noisy FID precision is retained on export.

```bash
pixi run python -m examples.legacy_basis_water --basis "$MRS_BASIS_H2O" \
  --output outputs/my_legacy_water --position 100 108 80 --seed 123
```

The legacy lipid/MM dependencies remain necessary; this is the full old simulation
workflow. The notebook's `legacy_basis_h2o_demo` tagged section uses the actual GE
PRESS TE35 file through the environment and actual MRiLab anatomy. It was executed
without fitting. Source BASIS/calibration files are not copied or modified.

Validation uses independent complex water sums, the real legacy signal model with
controlled nonzero MM/lipid providers, exact seeded repeats, independent water noise,
actual BASIS selection, relaxation factors and exact export/reload. A full seeded
legacy-default regression remains in the non-fitting suite. The notebook uses the
actual legacy lipid pathway, rather than those test doubles.

No scanner-level calibration is claimed. BASIS molar/proton provenance, ambiguous
table units and the retained legacy attenuation/time conventions remain limitations.
No fitting, calibration tuning, staging or commits are performed.

## Execution record

```bash
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m pytest -q tests --ignore=tests/integration
# 299 passed, 33 warnings, 137.24 s
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m pytest -q tests/test_legacy_basis_water.py
# Final targeted verification: 3 passed, 21 warnings, 3.85 s
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m examples.execute_quantitative_notebook --notebook notebooks/headless_simulation.ipynb --tag legacy_basis_h2o_demo --output outputs/legacy_basis_water_validation/notebook_final
# MRS_BASIS_H2O supplied externally; new code cell 14 passed.
```

The first test runs exposed fixture-only background-label and integer-concentration
errors, then an incorrect test assumption about the returned NIfTI-MRS object's
API. These were corrected using the existing legacy fixture conventions and an
actual save/read roundtrip. No legacy equations were changed to satisfy tests.
Exported water samples match exactly; analytical volume/shim expectations are
checked at rtol=5e-6 for the legacy complex64 arithmetic. No fitter ran.

Proposed staging only:

```bash
git add --pathspec-from-file=docs/commit_manifests/legacy_water_production.txt
# feat: integrate BASIS water referencing into legacy batch simulation
```

The manifest contains production source, tests, example, notebook and documentation;
no generated data or external BASIS. Existing prior-phase edits in shared files
should be reviewed before staging. No staging or commit was performed.

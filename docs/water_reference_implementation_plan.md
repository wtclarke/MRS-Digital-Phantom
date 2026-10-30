# Implementation Plan: Quantitative Water-Reference Simulation for MRS Digital Brain Phantom

## 1. Project Overview

**Repository:** https://github.com/dennisvds/MRS-Digital-Phantom

**Associated publication:** https://doi.org/10.1002/mrm.70138

**Objective:** Extend the MRS Digital Brain Phantom simulator to generate a separate, quantitatively calibrated, unsuppressed water reference suitable for validating metabolite concentration estimation in MRS fitting pipelines.

The implementation must preserve the existing residual-water simulation while introducing a distinct water-reference simulation pathway.

### Primary success criterion

When simulated metabolite data and their corresponding water reference are passed through a compatible fitting and quantification pipeline, the recovered metabolite concentrations should agree with the phantom's known ground-truth concentrations within expected numerical and fitting uncertainty.

The implementation must allow spectral fitting errors to be distinguished from errors introduced by water referencing, tissue corrections, and relaxation corrections.

### Scope

The new functionality should:

- Generate a complex unsuppressed water reference corresponding to the same voxel as the simulated metabolite spectrum.
- Support tissue-specific water concentrations and partial-volume mixing.
- Use a physically consistent signal-amplitude convention for metabolites and water.
- Support relaxation attenuation and relevant acquisition parameters.
- Preserve the existing residual-water contamination model.
- Export paired metabolite and water-reference spectra.
- Export ground-truth concentrations and simulation metadata.
- Support reproducible, automated validation.
- Remain backward compatible with existing simulations.

The implementation should initially prioritize quantitative correctness over realistic water line-shape complexity.

---

## 2. Existing Implementation

The repository already contains several components relevant to water-reference simulation.

### 2.1 Water signal generation

File: `simulation/signalModel.py`

The existing `simulate_water()` method generates a water signal using five components around 4.7 ppm.

It includes:

- Configurable water amplitude.
- Randomized phase and damping.
- Lorentzian and Gaussian decay.
- Tissue-specific water concentration scaling.
- Conversion between FID and spectral representations.

The current implementation contains:

```python
WF = self.simulation_params['water_amp_factor'] * 10e-3
```

This factor is used in the residual-water model and should not be reused as the quantitative calibration of an unsuppressed water reference.

### 2.2 Existing water storage

The simulator contains:

```python
if self.simulation_params.get("save_unsuppressed_water", True):
    self.unsuppressed_water = water_specs.clone()
```

This stores tissue-specific water spectra internally.

However, this variable is populated before the full voxel integration and shim processing applied to the final simulated spectrum.

It should therefore not be assumed to represent a complete, calibrated, voxel-level unsuppressed water acquisition.

### 2.3 Voxel-integrated water signal

The simulator already computes a voxel-integrated water component:

```python
specs_water = self.simulate_water()
spec_water = torch.matmul(voi_labels_flat, specs_water)
spec_water = spec_water.view(*voi_labels.shape[:3], -1)
spec_water = self.apply_shim_map(spec_water, self.shim_map)
spec_water *= voxel_volume
spec_water = spec_water.sum((0, 1, 2))
```

This provides a useful architectural reference for integrating the new water-reference pathway.

However, the existing signal remains tied to the residual-water model and its amplitude scaling.

### 2.4 Tissue water concentrations

File: `utils/definitions.py`

The existing simulator defines:

```python
H2O_CONCENTRATIONS = {
    'Background': 0.0,
    'WM': 0.703 * 55.5,
    'GM': 0.851 * 55.5,
    'CSF': 55.5,
    'Fat': 0.0,
    'Skull': 0.0,
}
```

These correspond approximately to:

| Tissue | Water concentration |
|---|---:|
| White matter | 39.02 M |
| Gray matter | 47.23 M |
| CSF | 55.50 M |

The values are based on the tissue-water assumptions described by Gasparovic et al. (2009).

Reference:

Gasparovic et al., *Quantitative spectroscopic imaging with in situ measurements of tissue water T1, T2, and density.*

https://doi.org/10.1002/mrm.22060

These values provide a useful starting point for quantitative water-reference simulation.

### 2.5 Basis-set normalization

File: `simulation/basis.py`

The basis-set handling includes normalization of the metabolite basis signals.

The coding agent must establish precisely how this normalization affects the relationship between metabolite concentration and simulated signal amplitude.

This is a critical prerequisite for quantitative water referencing.

---

## 3. Proposed Architecture

Introduce two distinct signal-generation pathways sharing a common quantitative framework.

```text
Phantom anatomy and tissue labels
                |
                v
    Tissue fractions within VOI
                |
                v
    Ground-truth concentrations
                |
       +--------+--------+
       |                 |
       v                 v
 Metabolite arm     Water-reference arm
       |                 |
 Metabolite basis    Water FID model
 Concentrations      Tissue water content
 Relaxation          Relaxation
       |                 |
       v                 v
 Shared acquisition and spatial model
       |                 |
       v                 v
 Metabolite spectrum   Unsuppressed water
       |                 |
       +--------+--------+
                |
                v
       Export paired datasets
                |
                v
       Fitting / quantification
                |
                v
    Compare recovered concentrations
          with ground truth
```

### Architectural principles

1. The metabolite and water-reference arms must share a documented amplitude convention.
2. Both arms must use the same voxel geometry and tissue-fraction definitions.
3. Relevant spatial effects should be applied consistently.
4. Water-reference generation must be independent of residual-water contamination.
5. Noise and acquisition-specific effects should be independently configurable.
6. Ground-truth concentrations must be retained separately from estimated concentrations.

---

## 4. Phase 1: Audit Existing Signal Normalization

### Objective

Determine exactly how known metabolite concentrations are converted into simulated spectral amplitudes.

### Files to inspect

- `simulation/signalModel.py`
- `simulation/basis.py`
- `simulation/water/water_sim.py`
- `utils/definitions.py`
- `gui/defaults.py`
- `gui/widgets.py`
- Existing export and serialization utilities.
- Existing simulation notebooks and tests.

### Required investigations

Trace the complete metabolite simulation pathway, including:

- Input concentration units.
- Tissue-specific metabolite concentrations.
- Voxel-weighted concentration calculation.
- Basis-set normalization.
- Metabolite proton-count conventions.
- T1 and T2 attenuation.
- Global amplitude scaling.
- Voxel-volume weighting.
- Shim and spatial weighting.
- FFT and FID normalization.
- Export-time amplitude transformations.

Determine whether the simulator currently preserves a consistent relationship between concentration and signal amplitude.

In particular, investigate the basis normalization in `simulation/basis.py` and whether its global scaling convention is suitable for quantitative comparison with water.

### Deliverable

Produce a short technical document:

`docs/quantitative_signal_model.md`

It should explain:

1. The existing metabolite amplitude convention.
2. The existing residual-water amplitude convention.
3. Whether the two conventions are quantitatively compatible.
4. Any required normalization changes.
5. Any assumptions that cannot be established from the code.

**Do not proceed with the quantitative water implementation until the amplitude conventions have been established.**

---

## 5. Phase 2: Define a Shared Quantitative Signal Model

### Objective

Introduce a consistent mathematical convention relating metabolite and water concentrations to simulated signal amplitudes.

### Simplified model

For a metabolite:

\[
S_m = K C_m n_m R_m
\]

For water:

\[
S_w = K C_w n_w R_w
\]

where:

- \(S_m\): Metabolite signal amplitude.
- \(S_w\): Water signal amplitude.
- \(C_m\): Metabolite concentration.
- \(C_w\): Water concentration.
- \(n_m\): Number of contributing metabolite protons, subject to the basis-set convention.
- \(n_w\): Number of water protons, normally 2.
- \(R_m\): Metabolite relaxation attenuation.
- \(R_w\): Water relaxation attenuation.
- \(K\): Shared signal calibration factor.

This model is a simplified reference convention, not a requirement to multiply every metabolite basis function by an explicit proton count.

If the basis already incorporates proton-number information, avoid applying it twice.

### Relaxation model

For a simplified acquisition model:

\[
R(T_1,T_2,TR,TE)
=
(1-e^{-TR/T_1})e^{-TE/T_2}
\]

This expression is appropriate only under the relevant excitation and steady-state assumptions.

The coding agent must verify whether the actual simulated sequence requires a more specific signal model.

### Mixed tissue model

For a voxel containing multiple tissue classes:

\[
S_w =
K n_w
\sum_t f_t C_{w,t} R_{w,t}
\]

where:

- \(f_t\): Tissue volume fraction.
- \(C_{w,t}\): Tissue-specific water concentration.
- \(R_{w,t}\): Tissue-specific relaxation attenuation.

For metabolites:

\[
S_m =
K \sum_t f_t C_{m,t} B_{m,t} R_{m,t}
\]

where \(B_{m,t}\) represents the relevant normalized basis response, including any modeled acquisition-dependent effects.

The implementation must not assume that a single averaged relaxation correction is equivalent to summing individually attenuated tissue signals.

### Concentration recovery

Under the simplified single-component model:

\[
C_m =
\frac{S_m}{S_w}
\frac{n_w}{n_m}
\frac{R_w}{R_m}
C_w
\]

For heterogeneous tissue, use the appropriate tissue-weighted signal model and define the target concentration convention explicitly.

### Requirements

- Document all units.
- Use consistent proton-count conventions.
- Ensure water and metabolites share a common amplitude calibration.
- Avoid arbitrary scaling introduced solely to force agreement with a fitting pipeline.
- Support analytically predictable amplitudes in idealized simulations.

### Deliverable

A documented quantitative signal convention, implemented through reusable code where appropriate.

---

## 6. Phase 3: Implement a Dedicated Water-Reference Model

### Objective

Implement a new water-reference generator independent of the existing residual-water model.

### Proposed interface

Introduce a method provisionally named:

```python
simulate_water_reference(...)
```

The precise signature should follow existing repository conventions.

### Required functionality

The method should:

1. Generate a complex unsuppressed water FID.
2. Use tissue-specific water concentrations.
3. Support partial-volume mixing.
4. Use the same spectral bandwidth and sampling as the metabolite simulation.
5. Use the same frequency-reference convention.
6. Support configurable water linewidth.
7. Support configurable frequency offset.
8. Support configurable phase.
9. Support T1 and T2 attenuation.
10. Use the shared quantitative amplitude calibration.
11. Support deterministic simulation.
12. Return a signal suitable for voxel integration and export.

### Initial water line shape

Begin with a single water resonance.

A simple model is:

\[
s_w(t)
=
A_w
e^{-t/T_2^*}
e^{i(2\pi \Delta f t+\phi)}
\]

where:

- \(A_w\) is the calibrated water amplitude.
- \(T_2^*\) controls effective linewidth.
- \(\Delta f\) is the frequency offset.
- \(\phi\) is the initial phase.

The initial implementation may use Lorentzian broadening.

Gaussian or Voigt broadening may be added later.

### Important distinction

The existing five-component residual-water model is intended to reproduce residual water contamination.

The new water-reference model is intended to represent a separate unsuppressed water acquisition.

These should remain separate components.

### Suggested code organization

Consider adding:

`simulation/water/water_reference.py`

Keep the existing:

`simulation/water/water_sim.py`

for residual-water generation unless the repository architecture suggests a cleaner separation.

### Deliverable

A tested quantitative water-reference generator with predictable amplitude and configurable spectral characteristics.

---

## 7. Phase 4: Integrate Water Reference with Voxel Simulation

### Objective

Ensure that the water reference corresponds to the same voxel and spatial simulation assumptions as the metabolite spectrum.

### Requirements

Use the same:

- Voxel position.
- Voxel dimensions.
- Tissue label maps.
- Tissue fractions.
- Spatial sampling.
- Voxel-volume weighting.
- Shim or B0 field map where applicable.
- Receive sensitivity assumptions.
- Spectral bandwidth.
- Sampling interval.
- Number of spectral points.
- Frequency reference.

The implementation should distinguish effects that are shared between acquisitions from those that are acquisition-specific.

### Recommended processing order

1. Obtain tissue fractions for the selected voxel.
2. Generate tissue-specific water signals.
3. Apply appropriate tissue-specific relaxation.
4. Apply relevant spatial and frequency effects.
5. Integrate signals across the voxel.
6. Add optional noise.
7. Return the final complex water-reference signal.
8. Export the signal and associated metadata.

### Important consideration

The current `self.unsuppressed_water` variable stores tissue-specific spectra before full voxel processing.

The new reference should represent the completed voxel-level signal, not merely a collection of intermediate tissue spectra.

### Deliverable

A paired metabolite and water-reference simulation generated from a common voxel definition.

---

## 8. Phase 5: Export Paired Spectra and Ground Truth

### Objective

Allow downstream fitting pipelines to consume the simulated spectra without custom in-memory integration.

### Required outputs

For each simulation, produce:

```text
simulation_001/
    metabolites.nii.gz
    water_reference.nii.gz
    ground_truth.json
    simulation_config.json
```

These filenames are proposed conventions and may be adapted to existing repository patterns.

### Metabolite dataset

Contains the normal simulated MRS signal, including configured contributions from:

- Metabolites.
- Macromolecules.
- Lipids.
- Residual water.
- Noise.

### Water-reference dataset

Contains the independently simulated unsuppressed water signal.

It should not inherit the residual-water amplitude factor.

### Ground-truth metadata

Include:

- Metabolite names.
- Ground-truth metabolite concentrations.
- Tissue-specific metabolite concentrations where available.
- Voxel-weighted concentrations.
- Tissue fractions.
- Tissue water concentrations.
- Effective water concentration.
- Metabolite and water relaxation parameters.
- Acquisition parameters.
- Signal normalization convention.
- Simulation seed.
- Version or commit identifier where available.

### NIfTI-MRS compatibility

Prefer NIfTI-MRS where practical.

Ensure that:

- Complex data are preserved.
- Spectral sampling is correct.
- Frequency metadata are correct.
- Dimension conventions are valid.
- Required acquisition metadata are present.
- Export does not unintentionally renormalize amplitudes.

Use existing repository export conventions where possible.

Custom simulator metadata may be stored in a separate JSON file rather than relying on downstream fitting tools to interpret custom NIfTI-MRS fields.

### Deliverable

Paired spectra and machine-readable ground truth that can be loaded by an external fitting pipeline.

---

## 9. Phase 6: Configuration and Backward Compatibility

### Objective

Allow users to enable the quantitative water-reference pathway without changing existing simulation behavior.

### Proposed configuration

```json
{
  "simulate_water_reference": true,
  "water_reference": {
    "add_noise": false,
    "apply_relaxation": true,
    "linewidth_hz": 5.0,
    "frequency_offset_hz": 0.0,
    "phase_deg": 0.0
  }
}
```

These are proposed parameters, not existing repository options.

The agent should follow established configuration patterns and avoid introducing incompatible schemas.

### Requirements

- Water-reference simulation must be optional.
- Existing residual-water settings must remain independent.
- Existing simulations must behave as before when the feature is disabled.
- Random-number generation must be reproducible.
- Enabling the water reference should not unintentionally alter metabolite noise or other stochastic simulation components.
- Default settings should be documented.

### GUI integration

After the core implementation is tested, expose the functionality through the GUI.

Suggested controls:

- Generate unsuppressed water reference.
- Water linewidth.
- Water frequency offset.
- Water phase.
- Apply relaxation.
- Add water-reference noise.

Avoid unnecessary GUI changes before the quantitative model is validated.

### Deliverable

Configuration support and optional GUI controls with preserved backward compatibility.

---

## 10. Phase 7: Automated Unit Tests

### Objective

Verify the quantitative correctness of the water-reference model independently of external fitting software.

### Test 1: Pure white matter voxel

Generate a voxel containing only white matter.

Expected result:

The water amplitude follows the configured white matter concentration and relaxation model.

### Test 2: Pure gray matter voxel

Generate a voxel containing only gray matter.

Expected result:

The water amplitude follows the configured gray matter concentration and relaxation model.

### Test 3: Pure CSF voxel

Generate a voxel containing only CSF.

Expected result:

The water amplitude follows the configured CSF concentration and relaxation model.

### Test 4: Mixed tissue voxel

Generate a voxel containing known WM, GM, and CSF fractions.

Expected result:

The simulated water signal agrees with the analytical tissue-weighted sum.

### Test 5: Water concentration linearity

Double the water concentration while holding other parameters constant.

Expected result:

The water signal amplitude doubles.

### Test 6: Residual-water independence

Change `water_amp_factor` while keeping water-reference settings constant.

Expected result:

The unsuppressed water-reference amplitude remains unchanged.

### Test 7: Noise-free deterministic simulation

Disable noise and stochastic spectral variation.

Expected result:

The generated signal agrees with the analytical model within numerical precision.

### Test 8: Reproducibility

Generate the same simulation twice using the same random seed.

Expected result:

The outputs are identical within the established deterministic numerical tolerance.

### Test 9: Backward compatibility

Run an existing simulation with the water-reference feature disabled.

Expected result:

The output remains unchanged relative to the baseline implementation.

### Test 10: Export and reload

Export a water reference and reload it.

Expected result:

Complex values, amplitude, spectral sampling, and required metadata are preserved.

### Test 11: Shared amplitude calibration

Generate a simplified metabolite signal and water reference with known concentrations and no relaxation attenuation.

Expected result:

The signal ratio agrees with the documented shared normalization convention.

### Deliverable

Automated tests integrated with the repository's existing testing framework.

---

## 11. Phase 8: End-to-End Fitting Pipeline Validation

### Objective

Demonstrate that the generated water reference can be used to recover known metabolite concentrations.

### Stage A: Idealized fitting validation

Generate a simple phantom with:

- One tissue class.
- One metabolite.
- Known metabolite concentration.
- Known water concentration.
- No noise.
- No macromolecules.
- No lipids.
- No residual water.
- No B0 inhomogeneity.
- No relaxation attenuation, or analytically matched attenuation.

Fit the metabolite spectrum and water reference using a compatible pipeline.

Compare recovered concentration against ground truth.

### Stage B: Multiple metabolites

Extend the simulation to several metabolites with known concentrations.

Verify that each concentration can be recovered independently.

Investigate any metabolite-specific bias arising from basis normalization or fitting-model mismatch.

### Stage C: Tissue mixing

Introduce mixed WM/GM/CSF voxels.

Verify that the quantification procedure recovers the explicitly defined voxel-average concentration.

Ensure that the quantification target is clear, particularly when distinguishing total voxel concentration from tissue-only concentration.

### Stage D: Relaxation

Enable tissue-specific water and metabolite relaxation.

Verify that the fitting pipeline recovers the correct concentrations when appropriate corrections are supplied.

Test the expected bias when corrections are intentionally omitted.

### Stage E: Realistic effects

Progressively introduce:

- Noise.
- Macromolecules.
- Lipids.
- Residual water.
- B0 inhomogeneity.
- Frequency offsets.
- Line broadening.
- Phase variation.
- Partial-volume effects.

Measure the effect of each addition on fitting and quantification accuracy.

### Validation metrics

For each metabolite, calculate:

\[
\text{Absolute error}
=
|\hat{C}_m-C_m|
\]

\[
\text{Relative error}
=
\frac{\hat{C}_m-C_m}{C_m}\times100\%
\]

Also report:

- Mean bias.
- Root mean squared error.
- Standard deviation across repeated simulations.
- Failure rate.
- Estimated fitting uncertainty, if available.

### Initial acceptance target

For an idealized, noise-free simulation with an exactly matched fitting model, target concentration recovery within 1%.

This is a proposed engineering target rather than a guaranteed tolerance for realistic spectra.

For noisy or model-mismatched simulations, use empirically justified tolerances and report uncertainty.

### Deliverable

A reproducible validation notebook or script that:

1. Generates paired spectra.
2. Exports the data.
3. Runs or prepares the fitting pipeline.
4. Loads estimated concentrations.
5. Compares estimates against ground truth.
6. Produces a summary table and diagnostic plots.

The implementation should avoid tightly coupling the simulator to one fitting package unless explicitly required.

---

## 12. Acceptance Criteria

The implementation is complete only when all the following conditions are met.

### Functional requirements

- [ ] A separate unsuppressed water reference can be generated.
- [ ] Water-reference generation is independent of residual-water contamination.
- [ ] Tissue-specific water concentrations are supported.
- [ ] Mixed tissue voxels are supported.
- [ ] Relaxation effects can be configured.
- [ ] The water reference uses the correct voxel geometry.
- [ ] Metabolite and water signals share a documented quantitative amplitude convention.
- [ ] Paired metabolite and water-reference datasets can be exported.
- [ ] Ground-truth concentrations and simulation assumptions are exported.

### Quantitative requirements

- [ ] Pure tissue water amplitudes agree with analytical expectations.
- [ ] Mixed tissue water signals agree with the weighted signal model.
- [ ] Water amplitude is linear with concentration.
- [ ] Metabolite-to-water amplitude ratios follow the documented calibration.
- [ ] Noise-free concentration recovery meets the defined numerical tolerance.
- [ ] Relaxation correction tests recover expected concentrations.
- [ ] No arbitrary post-hoc calibration is required to obtain agreement.

### Software requirements

- [ ] Existing simulations remain unchanged when the new feature is disabled.
- [ ] Random seeds support reproducibility.
- [ ] Unit tests pass.
- [ ] Integration tests pass.
- [ ] Export/reload tests pass.
- [ ] The implementation follows repository coding conventions.
- [ ] New configuration options are documented.
- [ ] Limitations and assumptions are explicitly documented.

---

## 13. Recommended Implementation Order

The coding agent should work in the following sequence.

### Step 1: Audit normalization

Establish the current relationship between metabolite concentrations and simulated amplitudes.

Do not modify the water model before this is understood.

### Step 2: Implement minimal water-reference FID

Start with:

- One tissue.
- One water resonance.
- No noise.
- No relaxation.
- No spatial effects.

Verify the amplitude analytically.

### Step 3: Establish shared calibration

Ensure the metabolite and water amplitudes follow a common quantitative convention.

Test known signal ratios.

### Step 4: Add tissue mixing

Introduce WM, GM, and CSF contributions.

Validate partial-volume behavior.

### Step 5: Add relaxation

Implement tissue-specific attenuation and verify analytically.

### Step 6: Integrate with the voxel model

Apply the existing spatial weighting and voxel integration architecture.

### Step 7: Implement export

Write paired spectra and ground-truth metadata.

Verify amplitude preservation.

### Step 8: Run end-to-end fitting validation

Begin with idealized data and progressively introduce realistic effects.

### Step 9: Add GUI and documentation

Expose the feature only after the core implementation is validated.

---

## 14. Critical Engineering Constraints

The coding agent must follow these constraints throughout development.

1. **Do not use arbitrary calibration factors.** Water amplitudes must follow the documented quantitative model, not empirical adjustments introduced solely to match fitted concentrations.

2. **Do not conflate residual water and unsuppressed water.** These represent different acquisition conditions and serve different purposes.

3. **Do not assume basis-set normalization is physically calibrated.** Verify the basis normalization and proton-count conventions.

4. **Do not apply corrections twice.** Relaxation, tissue fractions, and proton-number factors must be accounted for exactly once under the chosen convention.

5. **Do not silently modify existing simulation outputs.** Any change to metabolite amplitude conventions must be deliberate, documented, tested, and backward compatible where feasible.

6. **Do not use fitted concentrations as ground truth.** Ground truth must originate from the phantom's concentration model.

7. **Do not confuse voxel-average and tissue-only concentrations.** Define the quantification target and use it consistently.

8. **Do not overcomplicate the initial water model.** Start with a single physically interpretable resonance.

9. **Do not assume all fitting pipelines use the same scaling convention.** Keep the simulator's physical model independent of pipeline-specific quantification interfaces.

10. **Do not claim quantitative validation solely from visual spectral agreement.** Validate recovered numerical amplitudes and concentrations.

---

## 15. Expected Deliverables

The coding agent should produce:

### Code

- A dedicated water-reference simulation implementation.
- Any necessary shared amplitude-normalization utilities.
- Integration with the existing signal model.
- Configuration support.
- Optional GUI support.
- Export functionality.

### Tests

- Analytical unit tests.
- Tissue mixing tests.
- Relaxation tests.
- Backward-compatibility tests.
- Export/reload tests.
- End-to-end concentration recovery tests.

### Documentation

- Quantitative signal model.
- Description of amplitude conventions.
- Configuration instructions.
- Output format description.
- Ground-truth metadata specification.
- Validation example.
- Known limitations.

### Example dataset

At least one small, reproducible example containing:

```text
example_water_reference/
    metabolites.nii.gz
    water_reference.nii.gz
    ground_truth.json
    simulation_config.json
    validation_results.csv
```

The exact directory structure may be adapted to repository conventions.

---

## 16. Instructions for the Coding Agent

Implement a quantitatively calibrated unsuppressed water-reference simulation pathway in the MRS Digital Brain Phantom repository.

**Repository:** https://github.com/dennisvds/MRS-Digital-Phantom

**Primary goal:** Enable end-to-end validation of water-referenced metabolite concentration quantification against known phantom ground truth.

### Required workflow

1. Inspect the existing codebase and identify the complete metabolite amplitude-generation pathway.
2. Document the current basis normalization and signal scaling.
3. Identify whether the current metabolite simulation supports physically consistent water referencing.
4. Propose the smallest appropriate architectural changes.
5. Implement a dedicated water-reference generator independent of residual-water contamination.
6. Establish compatible quantitative amplitude conventions for water and metabolites.
7. Integrate the reference with the existing voxel model.
8. Export paired metabolite and water-reference datasets.
9. Export machine-readable ground-truth metadata.
10. Implement analytical and integration tests.
11. Demonstrate end-to-end concentration recovery.
12. Document assumptions, limitations, and usage.

### Development approach

Work in small, independently testable changes.

Before implementing the new water model, produce an amplitude-normalization audit and identify any ambiguities.

Where the current implementation is ambiguous, document the ambiguity instead of silently selecting a calibration convention.

Prefer reusing existing simulation infrastructure rather than duplicating voxel-processing logic.

Keep the implementation independent of any particular fitting pipeline unless integration is explicitly requested.

### Final report

At completion, provide:

- A summary of modified files.
- A description of the quantitative amplitude convention.
- An explanation of how the water reference is generated.
- Example usage.
- Automated test results.
- End-to-end validation results.
- Remaining limitations.
- Any decisions requiring scientific review.

---

## 17. Final Recommendation

The highest-priority task is the **metabolite amplitude-normalization audit**.

The simulator already contains much of the infrastructure needed to generate water signals, including tissue water concentrations, complex FID generation, voxel integration, and spatial effects.

However, generating a realistic water peak is not equivalent to generating a quantitatively calibrated water reference.

The essential scientific requirement is that the simulated water and metabolite signals share a defensible, documented relationship between concentration and measured amplitude.

The implementation should therefore prioritize:

**Quantitative signal consistency → analytical validation → voxel integration → export → end-to-end fitting validation → realistic spectral effects.**

This ordering will make it possible to identify whether concentration errors originate in the simulator, the fitting algorithm, or the subsequent water-referenced quantification procedure.
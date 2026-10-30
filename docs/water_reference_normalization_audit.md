## Phase 1 audit

Files modified: None.
Tests executed: None; this was a read-only audit.

### Verified metabolite signal path

1. Basis.reformat() in simulation/basis.py (simulation/basis.py:108) obtains formatted FIDs, then divides the entire basis array by mean(abs(all FID samples)). This is a global average-magnitude normalization across metabolites and time samples—not a per-metabolite peak or area normalization. The original basis amplitudes are therefore not preserved. There is no explicit proton-count normalization here.
2. SignalModel.simulate_metabs() in simulation/signalModel.py (simulation/signalModel.py:146) samples Conc_mean/Conc_std for each tissue and multiplies each basis FID by its concentration. Concentration units are not declared in this function. The code does not explicitly multiply by metabolite proton count.
3. It applies a T2-like decay using gamma = 1/(π T2) after converting the dataframe T2 from milliseconds to seconds, then exp(-gamma*t). No T1 or TR attenuation is applied in this metabolite path. The T2 units and decay equation are consistent with a Lorentzian linewidth if the input T2 is in milliseconds.
4. The basis FID is transformed with an unnormalized FFT. Tissue signals are mixed using the VOI labels, shim-modulated, multiplied by voxel volume, and summed in simulate_mrs_data() (simulation/signalModel.py (simulation/signalModel.py:80)). The ground-truth voxel concentration is separately computed as tissue-fraction-weighted sampled concentrations, without relaxation or basis-response weighting.

### Verified water path

simulate_water() (simulation/signalModel.py (simulation/signalModel.py:282)) creates five damped water components. Their amplitudes are determined by water_amp_factor, a fixed component-amplitude array, and generate_sum_fid()’s scaling argument (simulation/water/water_sim.py (simulation/water/water_sim.py:10)). It multiplies this signal by tissue water concentrations from utils/definitions.py (utils/definitions.py:172), then applies the same VOI label mixing, shim, voxel-volume weighting, summation, and unnormalized FFT convention as the metabolite path.

The water concentrations are specified as tissue fractions times 55.5 (WM, GM, CSF), consistent with the plan’s approximate molar values. The residual-water signal has no explicit water-proton factor, T1 saturation, or TE/T2 attenuation. save_unsuppressed_water stores the tissue-specific intermediate before voxel integration and shim; it is not the final voxel-level reference.

### FFT and export

utils/auxillary.py (utils/auxillary.py:21) uses unnormalized FFT and inverse FFT, so spectral values scale with point count relative to FID values. The simulation uses the same convention for both arms, so this alone does not establish or destroy their relative scale.

save_nifti_mrs() in utils/nii_processing.py (utils/nii_processing.py:36) converts spectra back to FIDs. Its one-dimensional path conjugates the FID for the NIfTI-MRS phase convention. results_to_nifti_mrs() uses a separate conversion path in simulation/headless.py (simulation/headless.py:220): inverse FFT, complex64 conversion, then NIfTI-MRS creation, without that explicit conjugation. Neither path visibly applies amplitude normalization, but they differ in phase handling. Reload preservation of complex amplitudes was not tested.

### Can the amplitudes be quantitatively compared?

Not established; the current code does not support a justified quantitative comparison. Both arms share voxel-volume weighting, shim application, and FFT conventions, but that is insufficient:

- The metabolite basis is arbitrarily rescaled by a global mean-FID magnitude, erasing the source basis calibration.
- Metabolites are concentration-scaled with no explicit proton-count convention; water uses molar water concentration and a separate empirical water_amp_factor plus component amplitudes/scaling.
- Metabolite T2 attenuation is applied, but no sequence T1/TR attenuation; water has randomized linewidth/phase and no analogous metabolite relaxation model.
- The metabolite ground-truth concentration mixture is not a signal-equivalent mixture when tissue-specific basis responses and T2 attenuation differ.
- Basis TE, sequence, and original basis amplitude conventions have not been reconciled with the concentration dataframe’s units.

The plan’s shared factor model S=K C n R cannot be confirmed from current code: there is no demonstrated shared K, and neither arm’s proton-count convention is defined well enough to derive one.

### Assumptions and uncertainties

Verified: code-level scaling, transforms, tissue weighting, and the different export phase handling described above.

Not verified: concentration units in the metabolite dataframe; proton content encoded in the source basis; whether source basis FIDs are already TE-attenuated or calibrated; intended meaning/units of water_amp_factor; whether the two export paths preserve the desired phase convention after reload. Those require inspecting the input basis provenance, concentration data provenance, and round-trip behavior.

### Proposed changes after review

Document the existing conventions first, including the basis normalization and residual-water model. Then decide whether to preserve source basis amplitudes or define a calibrated per-metabolite response, and define metabolite concentration units and proton-count handling without double counting basis-encoded protons. Specify sequence-specific T1/T2 attenuation and the target concentration convention. Only then define a water-reference amplitude on a demonstrably shared scale. Finally, reconcile or explicitly document the export phase convention and validate complex amplitude round trips.

Recommended next step: review these findings and resolve the basis and concentration provenance questions before Phase 2. No implementation was started.
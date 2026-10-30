# Phase 1b: Basis, concentration, and calibration provenance audit

**Scope:** provenance and amplitude audit only. No simulation code changed and no water-reference pathway implemented.

## Findings at a glance

- Basis files are LCModel `.BASIS` files with MRSCloud identifiers and sequence/TE metadata. MRSCloud generation code simulates spin systems with sequence pulse timings and emits metabolite-specific signals. Thus the resonance pattern and TE-dependent J evolution are represented. The checked-in BASIS header does not state a molar/proton normalization or an absolute receiver calibration.
- LCModel records `CONC = 1`, `TRAMP = 1`, and `VOLUME = 1` for each basis entry. These are basis metadata, not evidence that the sampled FID is absolutely calibrated to one mole or one proton.
- `Basis.reformat()` further divides all formatted metabolite FIDs by one shared mean absolute sample magnitude. This preserves relative amplitudes in that input array under that common scalar, but changes its absolute scale; it is not per-metabolite normalization. Any input relative amplitudes survive only insofar as the reader/reformatter preserves the raw basis consistently.
- `Conc_mean`/`Conc_std` are populated from the MRS Database workbook’s `mM_u`/`mM_std`, then combined-metabolite rows are split by fixed ratios and mM/IU fields are collapsed. The resulting values are nominally millimolar concentration summaries. The source represents reported in-vivo concentrations by tissue group; it does not establish intracellular, tissue-water, or wet-tissue-mass concentration semantics consistently across literature sources.
- The current decay implementation likely has a π error if `T2` is a physical relaxation time: it computes linewidth `1/(πT2)` in s⁻¹, then applies `exp(-linewidth*t)`. The conventional Lorentzian relation is `linewidth_FWHM = 1/(πT2)` Hz and the corresponding FID is `exp(-π linewidth_FWHM t) = exp(-t/T2)`. The current code therefore applies `exp(-t/(πT2))`, a much slower decay. It may have been intended as additional broadening, but the implementation uses the column named T2 and has no separate interpretation documented.
- No absolute or defensible relative metabolite-to-water calibration follows from current files alone. The basis response has no known concentration/proton calibration and the water arm has empirical scaling. A new calibrated mode needs an explicit convention and reference response; preserve current simulation behavior.

## 1. Basis provenance and amplitude

### Checked-in basis files

Example: `data/basissets/LCModel_Universal_Philips_UnEdited_PRESS_TE35.BASIS` declares `ECHOT = 35.00`, `SEQ = 'UnEdited'`, and `IDBASI = 'Universal_Philips UnEdited 35 ms Osprey'`. Entries include `ID`, `METABO`, and unit-valued `CONC`, `TRAMP`, and `VOLUME`, followed by complex FID sample pairs. Other checked-in files encode PRESS or sLASER and TE in their filenames and headers.

These headers verify sequence and TE labels, but do not document absolute signal units, a receiver gain, sample concentration, coil sensitivity, or a calibrated response per proton. `CONC=1` is LCModel basis bookkeeping; treating it as 1 M or one mole is unsupported.

### Generation code

`simulation/basissets/MRSCloud-main/functions/sim_signals.m` loops over metabolites and reads TE from `MRS_opt(ii).TEs{1}`; output names include vendor, sequence, localization, TE, and metabolite. The code evolves spin-system density matrices through excitation/refocusing/editing pulses and delays, using sequence timings derived from TE, then calculates the signal. This supports the conclusion that sequence-dependent resonance evolution, including TE-dependent coupling evolution, is modeled. A scalar relaxation attenuation corresponding to in-vivo metabolite T2 is not established by the inspected generation path. A universal numerical amplitude normalization by proton count or concentration is also not stated there.

**Proton counts:** A spin-system simulation sums magnetization from modeled spins, so proton multiplicity and relative resonance amplitudes are expected to arise from its spin system. That is a code-structure inference, not a verified per-metabolite calibration. The BASIS file alone does not encode explicit proton counts, and sequence pulse profiles, simulated metabolite structures/spin systems, and export normalization would need auditing to establish the exact convention. Therefore, adding an explicit proton count to basis amplitudes risks double counting.

### Two distinct normalizations

1. Original BASIS generation/export: the file stores the complex FID values produced by the simulation/export pipeline. Their absolute scale has no documented physical unit in the checked-in headers. Its relative resonance pattern is meaningful for fitting; cross-metabolite absolute response is not proven by metadata.
2. Runtime normalization: `Basis.reformat()` in `simulation/basis.py` calls `get_formatted_basis(bw, points, ...)` and divides the resulting whole array by `mean(abs(_raw_fids))`. This is one scalar shared across all metabolites and points, not individual unit-peak normalization. It destroys any absolute scale in the file while preserving within-array ratios up to formatting/interpolation/selection effects. Since the factor depends on the selected/formatted basis array, it can also vary with basis selection or requested sampling parameters.

## 2. Concentration provenance and units

`utils/preprocess_df.py` loads `MRS_Database.xlsx` sheets `Refs_Concs` and `Values_Concs`; it merges `mM_u`, `mM_std`, `tCr_u`, `tCr_std`, `IU_u`, and `IU_std`. `calculate_metab_conc()` combines entries by metabolite and tissue using inverse-variance/random-effects logic after selecting `IU_u`/`IU_std`. The preprocessing explicitly collapses mM and IU fields, and splits combined totals such as tCr, tNAA, Glx, and tCho using prescribed ratios. `create_metab_df()` writes those results to `Conc_mean` and `Conc_std`; checked-in `data/metabolites/metab_df.csv` contains typical values around 1–10, consistent with mM, but no column-unit declaration.

The workbook field names and the preprocessor establish the intended nominal unit as mM. However, the database combines heterogeneous studies and uses IU and mM measures, and `IU` is a fitting/reference convention whose conversion to absolute molar concentration varies by source. The pipeline adds the fields numerically, so it does not provide a single traceable physical conversion for every row. Treat units as nominal mM pending confirmation of the database’s definitions and row-level provenance.

Tissue labels GM/WM describe the tissue group to which concentration estimates are assigned. They do not prove that values are per total tissue volume versus per tissue-water volume, nor correct them for tissue water fraction. The simulation then treats concentration as a local tissue signal coefficient and multiplies it by voxel tissue fractions. This is a model convention, not demonstrated semantics of all source measurements.

## 3. T2/linewidth equation

In `SignalModel.simulate_metabs()` (`simulation/signalModel.py`), dataframe `T2` values are read in ms, and code calculates:

\[
T_2[s] = T_2[ms] \cdot 10^{-3},\quad \gamma = \frac{1}{\pi T_2},\quad s(t)=s_0(t)e^{-\gamma t}.
\]

Dimensionally, `gamma` is s⁻¹ (equivalently numerically Hz as cycles/s, but not angular frequency). A Lorentzian FID with FWHM linewidth \(\Delta f\) in Hz is:

\[
e^{-\pi\Delta f t},\qquad T_2=\frac{1}{\pi\Delta f}.
\]

Thus if `gamma` is the FWHM linewidth, the code should mathematically apply `exp(-π gamma t)` to represent T2, yielding `exp(-t/T2)`. As written it applies `exp(-t/(πT2))`, which is neither the conventional T2 decay nor conventional linewidth application. This is verified algebraically from the code; whether it was an intentional broadening surrogate is unresolved. The pre-existing audit called the equation consistent with Lorentzian linewidth, but that interpretation omitted the required π in the exponent and should be corrected in the audit record. Do not silently change runtime behavior; first decide whether to preserve legacy output or fix behind a calibration mode.

The basis already has a TE-specific sequence simulation, while this runtime multiplier adds an FID decay across acquisition time. It does not apply a TE echo attenuation factor `exp(-TE/T2)`. These are distinct effects: echo-time relaxation and acquisition-time FID decay. Avoid applying both to the same effect without clarifying whether the basis generation already includes T2 relaxation.

## 4. Calibration conclusion and explicit convention

Let \(b_m(t)\) be the raw per-metabolite basis response for unit reference amount, \(q_m\) its molar concentration, and \(n_m\) the effective contributing proton/spin response already present in that basis. An absolute model would require:

\[
S_m(t)=K\,q_m\,b_m(t)\,R_m(TR,TE,T_2),
\quad
S_w(t)=K\,C_w\,n_w\,b_w(t)\,R_w(TR,TE,T_2),
\]

where `K` is a common acquisition/receive calibration, `C_w` is molar water content per stated tissue-volume convention, and `n_w=2` only if water basis response is defined per molecule rather than already per proton. The simulator has no documented `K`, no calibrated `b_m` per molar amount, no guaranteed interpretation of `n_m`, and no paired water response. Therefore the ratio cannot be derived from existing data. It is not enough to use `Conc_mean` in mM and water concentration in M because the basis normalization and relative response are unknown.

Recommended explicit calibration convention for a separate mode:

- Define metabolite concentrations in mol/L of a clearly specified volume (recommend total tissue volume for direct volume-fraction mixing; record literature provenance and limitations).
- Define each basis as the complex response to a stated reference amount (e.g. 1 mol/L of metabolite molecules) under its sequence, with spin-system proton response already included exactly once.
- Define water as response to 1 mol/L of water molecules, with the two-proton response either in the simulated water basis or explicit, but never both.
- Apply a documented sequence-specific relaxation model once, and state whether the basis already includes sequence evolution or any relaxation.
- Keep an explicit shared arbitrary unitless signal convention only if the paired arms are calibrated to the same reference; do not add a fitted empirical multiplier.
- Preserve the current legacy mode and select the calibrated mode separately until its end-to-end numerical validation passes.

A more conservative initial validation can define a relative reference empirically from a known simulated spin-system / reference concentration pair, but that is calibration data and must be identified and justified, not an arbitrary amplitude factor.

## 5. NIfTI-MRS export paths and minimal round-trip test design

Two paths exist:

- `utils/nii_processing.py::save_nifti_mrs`: takes spectrum to FID through `spec2fid`; for a 1-D signal it reshapes and conjugates the FID, with a comment that this is the NIfTI-MRS phase convention. It saves via nibabel and writes a header extension.
- `simulation/headless.py::_as_nifti_mrs` / `results_to_nifti_mrs`: takes spectrum to FID using `spec2fid`, converts to `complex64`, reshapes, then calls `gen_nifti_mrs`; no explicit conjugation appears in this path.

Both `fid2spec`/`spec2fid` use unnormalized FFT/IFFT, so the inverse pair reconstructs samples (subject to dtype precision). The export paths have different explicit phase treatment; equivalence must not be assumed.

Minimal test: construct a deterministic complex FID with nonzero real and imaginary parts at multiple time points, compute its spectrum using the project’s `fid2spec`, export through each path to a temporary NIfTI-MRS file, reload each using its matching loader, and compare both complex FID samples and the reconstructed spectrum against the source. Assert expected phase/conjugation explicitly, peak magnitude, and relative complex error; use `complex64` tolerance and avoid a peak-only check that can miss phase inversion. Run separately for a 1-D result and the supported image-shaped input. This is a proposed test only; it was not executed.

## 6. Least disruptive strategy

Use a separate quantitative calibration mode and retain the existing simulation path unchanged. The mode needs a new calibrated basis response convention and a corresponding water reference model; changing runtime basis handling globally would silently alter existing outputs and may break saved simulations. For the new mode, either preserve raw basis amplitudes and attach an explicit calibration record, or generate/derive calibrated per-concentration basis responses. The current global `Basis.reformat()` normalization cannot be used as an absolute response without recording and undoing its exact factor. No current data supports a universal undo factor or proton multiplier.

Before implementation, decide:

1. Meaning/unit of each concentration field after IU/mM merging and combined-metabolite splitting.
2. Whether basis spin systems encode proton multiplicity and what concentration amount their amplitude represents.
3. Whether basis sequence simulations include any relaxation beyond sequence evolution.
4. Tissue-volume versus tissue-water concentration target.
5. T2 formula and legacy-compatibility policy.
6. Which export path and phase convention is canonical.

## References inspected

- `data/basissets/LCModel_Universal_Philips_UnEdited_PRESS_TE35.BASIS` and sibling PRESS/sLASER TE variants.
- `simulation/basis.py::Basis.reformat` and basis reader/reformat calls.
- `simulation/basissets/MRSCloud-main/functions/sim_signals.m`, `simMRS.json`, and repository basis README.
- `data/metabolites/MRS_Database.xlsx`, `data/metabolites/metab_df.csv`, `data/metabolites/concs_df.csv`, `utils/preprocess_df.py`.
- `simulation/signalModel.py::simulate_metabs`, `simulate_water`, and spatial integration.
- `utils/auxillary.py::fid2spec/spec2fid`, `utils/nii_processing.py::save_nifti_mrs`, and `simulation/headless.py::results_to_nifti_mrs`.

**Verification status:** code and file metadata reviewed. No simulations, tests, basis parsing, export/reload round trips, or external literature verification were performed. Proton calibration and concentration semantics remain unresolved beyond the workbook/preprocessing evidence above.

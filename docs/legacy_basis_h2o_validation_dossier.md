# Legacy BASIS H2O validation dossier

**Status:** Technical evidence synthesis for submission authors, reviewed 11 November 2026. No device-specific adequacy determination or FDA acceptance is asserted.

The legacy batch extension generates a separate BASIS-derived water reference alongside simulated MRS spectra and saved nominal voxel-average concentrations. Existing verification and a single-voxel LCModel experiment support its use for controlled numerical testing of fitting software under declared assumptions. The common synthesis mechanism motivates extension to other metabolites and compatible supplied bases. Observed concentration-recovery errors remain specific to the tested responses, mixture, acquisition and fitter configuration.

This dossier concerns the section **Batch simulation with a separate BASIS-H2O reference** in [headless_simulation.ipynb](../notebooks/headless_simulation.ipynb). It assembles existing evidence; it adds no simulation, fitting study, calibration, performance threshold or scientific code change.

## 1. Intended use of the evidence

A proposed context of use is to supply paired synthetic metabolite and water acquisitions, with known nominal concentration inputs, for numerical evaluation of MRS fitting and water-referenced quantification software. The simulator contributes supporting test evidence alongside other evidence selected by the device sponsor. This context does not establish patient-level accuracy, diagnostic performance or physical scanner calibration.

FDA's [Assessing the Credibility of Computational Modeling and Simulation in Medical Device Submissions](https://www.fda.gov/media/154985/download), issued 17 November 2023, recommends relating model credibility to its context of use, influence on the decision, and consequences of error. Its reporting framework can help organize the sponsor's assessment. This dossier does not determine the adequacy of evidence for a particular device or qualify the toolkit as a medical device development tool.

| Sponsor field | Information to complete for a submission |
| --- | --- |
| Product and version | Fitting software, configuration, and version tested |
| Question of interest | Specific numerical performance question addressed by simulation |
| Intended clinical use | Claimed outputs, metabolites, population and interpretation |
| Acquisition scope | Applicable vendors, field strengths, sequences, TE/TR, sampling and preprocessing |
| Model influence and consequence | Reliance on simulated results, complementary evidence, and consequences of incorrect estimates |
| Acceptance and adequacy | Sponsor-defined criteria, their justification, and assessment against the complete evidence |
| Controlled records | Simulator source snapshot, dependencies, basis identity, dataset identities, review owner and approval status |

The present record has no prospectively approved concentration-recovery acceptance limits. Its results characterize the exercised case; proximity to a chosen percentage cannot be retroactively treated as a prespecified pass criterion.

## 2. Model and concentration conventions

### Pathway identity

The relevant interface is [simulate_batch](../simulation/headless.py) with legacy `simulation_params`, `basis_set_settings`, and `water_reference_mode='basis_h2o'`. The example configuration is [legacy_basis_water.py](../examples/legacy_basis_water.py); reference preparation, synthesis and export are in [legacy_water_reference.py](../simulation/legacy_water_reference.py). The [legacy workflow description](legacy_basis_water_reference.md) documents the options.

The similarly named [standalone BASIS-H2O pathway](basis_h2o_water_reference.md) uses the `quantitative` schema and different spatial and normalization conventions. Its pure-tissue, relaxation and noise results are supporting component context, not direct validation of this legacy integration. Frozen sequence-matched quantitative simulations are also a separate pathway.

### Shared signal formation

Let `B_m` and `B_w` denote formatted metabolite and water FIDs from the supplied basis. A single divisor is calculated before metabolite selection:

```text
q = mean(abs(all formatted BASIS columns and samples))
```

Metabolites and water use this same divisor, including H2O and control columns in its definition. The selection-independence test demonstrates that changing selected metabolites preserves retained responses for its fixture. Changing the basis or sampling can change `q`; the new value and basis identity must remain attached to that acquisition.

In schematic FID form, the implemented linear sums are:

```text
S_m(t) = sum_j volume_j * sum_tissue fraction_j,tissue
         * C_m,tissue * B_m(t)/q * D_m,tissue(t) * P_j(t)

S_w(t) = sum_j volume_j * sum_tissue fraction_j,tissue
         * 1000*C_w,tissue*R_w,tissue * B_w(t)/q * D_w(t) * P_j(t)
```

Here `C_m` is nominal mM, `C_w` is mol/L, `volume_j` is in mm³, `D` represents added acquisition-time decay and `P` the actual legacy shim modulation. The factor 1000 converts mol/L to mM. No additional proton-count multiplier is introduced: the explicit assumption is that the supplied responses already encode proton contributions on a common scale. A shared runtime divisor preserves that scale assumption; it cannot establish or repair missing upstream molar/proton calibration.

The legacy shim has a negative frequency sign and a time vector starting at one dwell time. Legacy metabolite decay is `exp(-t/(pi*T2))`, with table T2 converted from milliseconds. It is not conventional physical `exp(-t/T2)` decay and includes no separate metabolite T1/TR saturation or echo-time T2 attenuation. These retained conventions affect interpretation, rather than supplying physiological relaxation calibration. The later [provenance audit](water_reference_provenance_audit.md) corrects the earlier normalization audit's interpretation of this decay.

Extra water relaxation defaults to unity. If explicitly enabled with the required basis assumption, it is `(1-exp(-TR/T1))*exp(-TE/T2)`, with times in seconds. It therefore does not make the legacy metabolite and water relaxation models symmetric. The basis H2O response also retains the supplied response's RF attenuation, phase and processing; it is not an independently generated water-centered acquisition.

Residual contamination remains the legacy `Water` component. The separate reference is not added to `total`. Water synthesis reuses the existing shim and tissue representation; optional water noise has a separate random stream. This separation is visible in the implementation; the available direct legacy tests do not include a sweep over residual-water settings.

### Ground truth and tissue correction

Read `ground_truth.json` from the same voxel directory as the fitted acquisition:

```python
import json
from pathlib import Path

voxel_dir = Path("outputs/legacy_basis_water_notebook/<run>/<voxel>")
truth = json.loads((voxel_dir / "ground_truth.json").read_text())
concentrations_mM = truth["nominal_voxel_average_metabolite_concentrations_mM"]
wconc_mM = 1000 * truth["water_reference"]["nominal_water_concentration_mol_l"]
```

`concentrations_mM` is a copy of `voxel_weighted`: the actual sampled tissue concentrations weighted by voxel tissue fractions, not the concentration-table means and not fitted estimates. The underlying literature table pools mM and institutional-unit sources. Its interpretation as nominal mM is an explicit extension assumption, not a demonstrated IU-to-molar conversion.

With default tissue water contents:

```text
WCONC [mM] = 1000*(39.0165*f_WM + 47.2305*f_GM + 55.5*f_CSF)
```

Use the saved water concentration when defaults are overridden. LCModel handles water's two protons; do not double `WCONC`. Fractions are relative to the represented voxel, including any background; do not silently renormalize them to WM+GM.

Compare LCModel's concentration column against voxel-average truth, matching metabolite names. Its ratio-to-reference column, %SD and COORD plotting levels are different quantities. No further CSF correction is needed for this comparison. If a WM+GM-only target is selected, use the saved `nominal_tissue_only_metabolite_concentrations_mM`. Dividing voxel-average estimates by `f_WM+f_GM` is equivalent only when other tissues contribute no metabolite. The saved tissue-only target is undefined for zero WM+GM. Group comparisons sum the same named constituents on both sides; individual %SD values cannot be combined as independent errors.

## 3. Evidence and applicability

| Claim or question | Evidence | Supported interpretation and boundary |
| --- | --- | --- |
| Anatomical and spectral model foundations | Original paper, §§2.2–2.3 | Literature profiles, basis responses, tissue integration and empirical residual-water/MM/lipid models; no separate water-reference validation |
| Realism relative to acquired spectra | Original paper, §§2.4.2 and 3.3 | Comparison of 104 acquired and 480 simulated spectra; not independent absolute-concentration validation |
| Shared normalization and selection behavior | [Direct legacy tests](../tests/test_legacy_basis_water.py), `test_fixed_normalizer_and_validation` and `test_actual_basis_selection_independence` | One full-basis scalar and retained-column equality for a synthetic basis fixture |
| Tissue/volume/shim water formation | Same test file, `test_real_legacy_metabolites_shim_noise_and_paired_export` | Independent complex water expectation with volume and dwell-shifted shim; `rtol=5e-6`, `atol=0.01` in fixture signal units; not a concentration-range sweep |
| Repetition, water noise and export | Same test | Exact seeded repeat and water save/reload; water noise leaves spectrum and clean reference unchanged; NumPy/Torch global RNG restoration checked |
| Optional water relaxation | Same test file, `test_fixed_normalizer_and_validation` | Required declaration and calculated factor checked; not a full mixed-tissue relaxed-signal recovery study |
| Actual legacy notebook execution | [Execution record](../outputs/legacy_basis_water_validation/notebook_final/notebook_execution.json) and [workflow report](legacy_basis_water_reference.md) | Historical GE/MRiLab execution without fitting; does not establish execution of the current notebook source |
| Water-scaled concentration recovery | [LCModel conformance report](lcmodel_legacy_water_conformance.md) and saved artifacts below | One voxel, one seed, one supplied basis, one fitter version; individual and grouped errors reported below |
| Broader nominal quantitative behavior | [Standalone BASIS-H2O report](basis_h2o_water_reference.md) and related quantitative reports | Different pathway; no automatic transfer of its numerical tolerances to the legacy example |

The [original paper](Magnetic%20Resonance%20in%20Med%20-%202025%20-%20Sande%20-%20A%20Digital%20Phantom%20for%20MR%20Spectroscopy%20Data%20Simulation.pdf), DOI [10.1002/mrm.70138](https://doi.org/10.1002/mrm.70138), used Osprey estimates from acquired spectra as simulation inputs and adjusted remaining parameters for resemblance (PDF p8). Its comparison reported incomplete distributional overlap and significant differences for some fitted components (PDF pp10–11). Similar Cr linewidths therefore support a realism argument, not a known-concentration recovery bound. Residual water in that paper is distinct from this extension. Its reproducibility review excluded scientific accuracy and efficiency (PDF p13).

The legacy report records a historical non-fitting suite of 299 passes and a subsequent focused run of three passes. The currently inspected test file contains four tests. These historical counts are not current execution results or evidence that every current source change passed. The model/specification audits also describe earlier implementation states; where they conflict, use the identified current pathway and dated execution artifacts.

## 4. Observed LCModel experiment

### Acquisition and fitting protocol

The saved acquisition is [voxel_0000_100_108_80](../outputs/legacy_basis_water_notebook/20261101T182122316149Z/voxel_0000_100_108_80/). It uses a 10×10×10 mm voxel, GE PRESS TE35 basis, TR 2000 ms, 2048 points, 2000 Hz bandwidth, 127.731 MHz, and seed 123. Tissue fractions are approximately WM=0.194, GM=0.806, CSF=0. The saved nominal water concentration is 45.63698676815629 mol/L; float32 fraction rounding is retained.

The saved source used shim amplitude 3 Hz and residual-water amplitude factor 5. The current notebook uses 5 Hz and 0, respectively. Thus these results demonstrate the pathway under a related saved configuration, not the exact current example. `individual_components=["Metabs"]` saves a separate metabolite acquisition; it does not remove MM, lipid, residual-water or noise contributions from the default total spectrum.

The [runner](../research/validation/lcmodel_legacy_water/run.py) fitted saved total and metabolite-only acquisitions against the matching supplied basis with LCModel 6.3-1N. All 18 present metabolites were fitted. The [saved CONTROL](../outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits/total_voxel_water/lcmodel.CONTROL) specifies water scaling, no eddy-current correction, `ATTH2O=ATTMET=1`, NAA scaling at 2.01 ppm with three protons, `PPMCEN=4.68`, window 0.2–4.2 ppm, `NSIMUL=0`, `NRATIO=0`, and default baseline/line-shape settings. H2O and six unused entries, including bHG, are omitted. Total-spectrum MM/lipid contributions are not explicitly represented by this spin-system-only fitting protocol.

`FWHMBA=1/127.731=0.007828953... ppm` was inherited from an earlier 1 Hz synthetic-basis experiment. That conversion alone does not establish the supplied GE basis's appropriate effective width. The [integration-sensitivity investigation](lcmodel_fit_bias_diagnostics.md#4-reference-integration-sensitivity) shows that FWHMBA can change concentration scaling. The setting is recorded for reproducibility, not endorsed as a universal basis-width choice. LCModel parameter definitions are in the [manual, §§8.6.4.1 and 10.2](https://www.lcmodel.com/pub/LCModel/manual/manual.pdf).

| Case | Acquisition | WCONC, mM | Phased-real relative residual |
| --- | --- | ---: | ---: |
| total_previous_water | Total, historical pure-WM scaling | 39016.50000000 | 1.80624% |
| total_voxel_water | Same total, voxel-specific scaling | 45636.98676816 | 1.80623% |
| metabs_voxel_water | Saved metabolites only, voxel-specific scaling | 45636.98676816 | 0.72117% |

Residual is `norm(data-fit)/norm(data)` on LCModel's phased-real COORD arrays. It is not concentration error. The near-identical total residuals accompany a WCONC scale increase of 1.1696842815: NAA changes from 6.794 to 7.947 mM without changing input signals. This supports the expected role of voxel-specific water scaling, not a new signal-calibration factor.

### Individual and grouped concentration recovery

The following values reproduce the saved [comparison CSV](../outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits/comparison.csv), rounded for display. Relative error is `100*(estimate-truth)/truth`. Truth is nominal voxel-average mM; the LCModel values are its reported concentrations under the declared scaling assumptions.

| Metabolite | Truth | Total estimate | Total error | Metabolites-only estimate |
| --- | ---: | ---: | ---: | ---: |
| Asc | 1.073478 | 0.418 | -61.06% | 0.881 |
| Asp | 2.021706 | 1.853 | -8.34% | 1.693 |
| Cr | 5.471112 | 5.026 | -8.14% | 5.053 |
| GABA | 1.730572 | 1.576 | -8.93% | 1.338 |
| Gln | 2.409640 | 2.195 | -8.91% | 2.036 |
| Glu | 8.440927 | 8.909 | +5.55% | 8.175 |
| Gly | 4.722277 | 4.111 | -12.94% | 4.476 |
| GPC | 0.727706 | 0.365 | -49.84% | 0.641 |
| GSH | 1.597329 | 1.473 | -7.78% | 1.479 |
| Lac | 1.444271 | 0.980 | -32.15% | 1.342 |
| mI | 4.450162 | 3.950 | -11.24% | 4.068 |
| NAA | 8.176772 | 7.947 | -2.81% | 7.894 |
| NAAG | 0.960166 | 0.718 | -25.22% | 0.852 |
| PCh | 0.527690 | 0.847 | +60.51% | 0.543 |
| PCr | 2.666186 | 2.990 | +12.15% | 2.733 |
| PE | 1.826958 | 2.103 | +15.11% | 1.761 |
| sI | 0.304593 | 0.236 | -22.52% | 0.264 |
| Tau | 1.975306 | 1.491 | -24.52% | 1.683 |

| Group | Truth | Total estimate | Total error | Metabolites-only estimate | Metabolites-only error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Cr+PCr | 8.137297 | 8.016 | -1.49% | 7.786 | -4.32% |
| GPC+PCh | 1.255396 | 1.212 | -3.46% | 1.184 | -5.69% |
| NAA+NAAG | 9.136938 | 8.665 | -5.17% | 8.746 | -4.28% |

The [grouped CSV](../outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits/grouped_comparison.csv) sums individual TABLE estimates. Good grouped agreement does not establish individual recovery: GPC and PCh total-fit errors have opposite signs and large magnitudes. Low residual or low reported %SD does not establish unbiased quantification; there is only one realization, so these results provide no empirical precision or coverage estimate. The metabolite-only fit is a diagnostic of the same saved simulation, not an independent experimental reference.

## 5. Rationale for extension to other metabolites and bases

The 18-metabolite experiment exercises a general synthesis-and-comparison chain rather than a metabolite-specific amplitude correction. Metabolite responses enter through the same concentration multiplication, tissue weighting, common normalization, shim, volume sum and export operations. Their linear form motivates applying that implementation to additional basis entries, such as bHG, and to other compatible supplied bases.

This is a mechanistic applicability argument. A different entry can use the same implementation properties without inheriting a measured recovery error. It must have a supplied response, consistent identity/alias mapping, concentration-table rows and the declared common response scale. A different basis must also satisfy the loader's sampling, header and H2O requirements and match the acquisition's sequence/TE and frequency convention. Independently normalized responses without a recoverable common scale do not satisfy that assumption. Preserving a common divisor cannot fix such a basis.

Transfer is strongest for the shared arithmetic and file-handling mechanisms. Basis-specific proton/RF response quality, spectral overlap, concentration, line shape, nuisance signals and fitting identifiability remain separate dependencies. Therefore this case motivates broader use of the toolkit for controlled fitting-software tests, but supplies no observed bHG accuracy, universal recovery tolerance, vendor/sequence coverage or scanner-level molar calibration. Other bases' scalar and identities must be recomputed and recorded rather than copied from the GE case.

## 6. Artifact traceability and reproducibility limits

The [fit artifact directory](../outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits/) contains per-case TABLE, PRINT, COORD, PS and PDF reports, CONTROL/RAW/H2O files, [summary.json](../outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits/summary.json), [input_hashes.json](../outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits/input_hashes.json), and [commands.txt](../outputs/lcmodel_legacy_water_conformance/20261101T182122316149Z_fits/commands.txt). The [conformance report](lcmodel_legacy_water_conformance.md#reproducibility-and-reports) gives the historical run and analysis commands. Paths in those records describe that execution environment; availability and redistribution are separate issues.

On 11 November 2026, read-only inspection found all three cases' reports present and nonempty, with saved exit codes of zero. The historical summary records unchanged inputs and zero complex RAW/H2O conversion error. Nine of the ten recorded input files were available and their SHA256 values matched the manifest, including the saved acquisition/sidecar files and local LCModel executable. The original GE BASIS was absent from its recorded external path, so its recorded source SHA256 (`f0a42cf0bef1b5427d4f3f3f4e5dd35dd4406581230a7ead1055bd7fbd700a79`) could not be independently rechecked there. An adapted BASIS remains with the fit artifacts; it is not a byte-identical substitute for the absent source.

This inspection is anchored to checkout HEAD `01036a30f7b827894a27d2c7ec7fbee7b80027c5`, with modified and untracked sources present. HEAD alone does not identify the reviewed working tree or the historical executed simulator. A submission copy needs a controlled source/artifact snapshot; historical commands and input hashes do not establish full source-version traceability.

The current notebook also accesses `original_wref_settings` before initializing it in the following cell. That section cannot run top-to-bottom in a fresh kernel as written, and its configured source-basis path is unavailable locally. These limitations are recorded rather than repaired by this documentation task. The historical notebook record cannot certify the current cell order or personal path.

### Documentation acceptance checks

The dossier's checks concern documentation integrity, not new validation execution:

- Confirm that all 18 displayed rows and three grouped rows match the saved CSVs to their displayed precision, and that each case's TABLE concentration agrees with its summary.
- Recalculate displayed errors, group sums, WCONC, and residual percentages from saved values; confirm consistent nominal-mM and mol/L conversion.
- Check local references and report availability, and compare available input hashes against the saved manifest while identifying missing external files.
- Keep pathway, configuration, historical execution status and evidence limits explicit. No simulation, fitter or pytest rerun is represented by these checks.

**Inspection outcome, 11 November 2026:** all 18 individual rows and three grouped rows matched the saved CSVs at displayed precision; all three TABLE/summary concentration comparisons, displayed residuals and nominal water calculation passed. All 24 local references resolved, all 15 expected report files were nonempty, and all nine available recorded input hashes matched. The missing original basis remains the tenth manifest entry. These checks establish consistency of this dossier with the available records, not a new simulator or fitting validation run.

The supported conclusion is that the extended legacy framework supplies paired signals and nominal ground truth for controlled numerical evaluation, with direct implementation checks and an informative single-voxel LCModel demonstration. Its shared mechanism motivates additional metabolites and compatible bases under the same assumptions. Absolute calibration, unmeasured recovery bounds and device-specific sufficiency remain outside the demonstrated evidence.

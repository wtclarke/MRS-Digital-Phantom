# Headless water reference characterization

This study measures concentration sensitivity of the existing headless simulator and its approximate BASIS-H2O water reference using FSL-MRS and LCModel. It characterizes bias, tissue mixing, cross-talk, noise and poorer shimming without tuning signal amplitudes or imposing fit-accuracy pass/fail thresholds.

## Run and reproduce

From the repository root:

```bash
MPLCONFIGDIR=/tmp/headless-validation-mpl OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  pixi run python -m research.validation.headless_water.run \
  --output outputs/headless_water_validation
```

`--bases` selects the external BASIS directory. The default is the supplied VoxelSystems `LCM_basis_2` directory. `--limit 2` in a new output directory runs a smoke check; `--stages clean noise` selects stages. Existing matching case records are reused. A changed driver, adapter, simulator, input table, BASIS, or selected design requires a new output directory. No source acquisitions, external BASIS or prior notebook outputs are overwritten. LCModel uses the existing licensed local executable specified by `BINARY` in its adapter.

The driver generates 1,312 acquisitions and sends each to both fitters. Four isolated worker processes generate acquisitions and four fit them. Scoped loader replacements supply controlled anatomy; production code is unchanged. Cached BASIS readers return independent object copies. Zero-amplitude MM/lipid computations are bypassed only in clean research cases; realism uses the actual legacy providers.

## Methods

Detailed bases are GE PRESS TE35 and SH PRESS TE97. The other ten supplied vendor/TE combinations receive smaller NAA/bHG transfer checks. Sampling is 2048 points at 2000 Hz, TR 2000 ms, 20-mm cubes, actual BASIS TE and spectrometer frequency. Source amplitudes, intrinsic sequence evolution, shared full-BASIS divisor, and encoded-proton assumptions are preserved.

Controlled WM/GM/CSF percentages are 100/0/0, 0/100/0, 50/50/0, 40/40/20, 25/25/50, 10/10/80 and 0/0/100. Exact cell counts are shuffled deterministically within the same 1,000-cell volume. Additional label classes exist only outside the simulated region. The existing nominal concentration table supplies baseline values; missing values and CSF metabolites are zero. Concentration SD is zero. bHG is explicitly 4 mM in WM and 3 mM in GM. These are nominal demonstration concentrations, not newly asserted physiological defaults.

NAA, Cr, bHG, Glu, GABA and GSH are independently scaled by 0, 0.5, 0.75, 1, 1.25, 1.5 and 2. GPC/PCh are jointly scaled with their baseline ratio preserved. All applicable BASIS metabolites remain in the fit. tNAA, tCr and tCho are sums of individual estimates; their uncertainty is not obtained by assuming independence. The zero-CSF-metabolite condition characterizes boundary behavior and is excluded from meaningful sensitivity slopes.

Clean cases use NaN metabolite T2 to disable extra legacy decay, and no additional water attenuation. Both preserve evolution already encoded in the BASIS. Realism checks restore the original metabolite T2 expression and then MM, lipid and residual-water signals. Those checks deliberately contain model mismatch and do not imply matched physical relaxation corrections.

Noise uses the toolkit's native complex spectral generator, with `noise_level` set from baseline isolated NAA peak / desired per-component noise SD / 20000. Nominal SNR levels are 100, 30 and 10; noise amplitude stays fixed during concentration changes. Ten independent seeds are paired across concentration conditions. Reference noise uses the native independent FID-noise generator at nominal water SNR 200 and 50, with the Fourier-domain square-root-of-N conversion recorded. The authoritative `acquisition_checks.csv` derives achieved SNR from saved signal/noise pairs and preserves the native applied noise SD. Its water SNR supersedes the generation-time `achieved_snr.json` value, which used the already-noisy water peak. Metabolite achieved SNR uses the total metabolite peak and is distinct from the nominal isolated-NAA reference.

Poorer shim means increasing B0 dispersion, not increasing shimming quality. A fixed smooth spatial field is scaled to SD 0, 2, 5 and 10 Hz and applied to both acquisitions through the existing B0 path. Anatomical confirmation selects fully WM/GM/CSF voxels closest to pure WM, pure GM and 40/40/20 on a deterministic grid. Anatomical cases retain the notebook's random shim defaults unless supplied anatomy B0 overrides them. `applied_field.json`, rather than the controlled-design `b0_sd` column, records their actual field statistics.

Each controlled case independently reconstructs complex water from the tissue-specific molar concentrations, formatted H2O response, full-BASIS normalization, applied field and cell-volume sum. The downstream water-scaling controls deliberately multiply saved reference amplitudes and declared water concentrations together; an explicitly mismatched case scales only the signal. These perturbations isolate the fitters’ scaling response. Native generator concentration linearity is checked separately by the existing doubled-water-profile regression test, not inferred from this multiplication. Existing tests cover residual-water independence and export/reload preservation.

## Fitting and concentration conventions

Both fitters use all applicable supplied metabolite responses, 0.2–4.2 ppm, disabled ECC and no concentration-ratio priors or added empirical MM signals. Individual optimizer behavior is retained. LCModel uses actual TE, source FWHMBA, unit water/metabolite attenuation and voxel-weighted WCONC. Converted RAW/H2O and adapted BASIS preserve complex arrays and amplitudes; original inputs are hashed.

FSL-MRS uses native deterministic Newton/Voigt fitting with baseline off. Its fixed PCr methyl reference (3 protons, 2.9–3.15 ppm) has nonzero ground-truth concentration when NAA or Cr is swept to zero; under stress its fitted coefficient can still reach zero and prevent quantification. Quantification is configured in the unrelaxed TE=0/TR=infinity limit while acquisition metadata retain actual timing. Water densities explicitly match the simulator. Native tissue-only molarity is retained and multiplied by `(1-f_CSF)` to yield voxel-average estimates. This is a declared concentration-convention conversion, not an empirical amplitude correction. LCModel's NAA reference uses its BASIS area, independently of the fitted NAA concentration.

Primary truth is nominal voxel-average mM. Pure CSF has no tissue-only target or fitted reference metabolite; quantification failures are recorded, not discarded. The study reports raw bias, free-intercept slope/intercept, R², baseline-relative delta/fold error, zero-input leakage, cross-talk, empirical SD/RMSE, and failed/missing estimates. Ten repeats support screening-level variability estimates.

## Outputs and interpretation

The output root includes `acquisition_checks.csv` with independent water errors, actual applied field statistics, measured magnitude-spectrum water FWHM, and noise diagnostics. Magnitude FWHM is not an absorption-spectrum linewidth.

The output root contains cases and estimates CSVs, stratified accuracy and failure tables, concentration-response and relative-change tables, noise/shim summaries, a generated Markdown results report, and PDF/PNG figure pairs. Each case retains source acquisition/ground-truth/config/provenance, independent water checks, SNR and applied-field diagnostics, and native fitter records. `manifest.json` links the study to source and input hashes.

Compare response slopes and relative changes separately from additive offsets and multiplicative gain differences. A low residual does not establish accurate individual metabolite separation. Increasing B0 dispersion need not worsen every fit because the fitters model linewidth. Failures and optimizer warnings remain auditable. Matched-BASIS results characterize internal consistency rather than independent scanner-level calibration.

References: [FSL-MRS quantitation](https://open.oxcin.ox.ac.uk/pages/fsl/fsl_mrs/quantitation.html), [LCModel manual](https://www.lcmodel.ca/pub/LCModel/manual/manual.pdf), and the repository's [legacy reference conventions](../../../docs/legacy_basis_water_reference.md).

## Completed study

The [results report](../../../docs/headless_water_validation.md) gives the numerical conclusions and limitations. Compact CSV tables and PDF/PNG figures are retained in [results/](results/). The full native fitting output remains under `outputs/headless_water_validation/`; it is approximately 14 GiB and is ignored by Git.

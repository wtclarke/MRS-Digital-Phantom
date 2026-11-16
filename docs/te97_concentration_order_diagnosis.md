# Siemens PRESS TE97: collaborator concentration-order diagnosis

## Finding

The simulator applies the concentration and T2 vectors in CSV row order to BASIS
FID columns sorted by metabolite name. Adding bHG at the bottom of the supplied
CSV changes the intended name-to-column mapping. This is a simulator indexing bug,
not evidence for an NAA spin-response correction factor.

The PDF reports NAA **2.534**, while saved nominal voxel NAA is
**10.008234539734982**. The NAA FID column actually receives the **NAAG** concentration
**2.6149676046421693**. Thus the apparent factor is 10.00823454/2.61496760 =
**3.82728816**. LCModel's NAA is only **3.0963% below the wrongly assigned input**.
Multiplying NAA by four would hide this symptom while leaving other metabolites
and their linewidth assignments wrong.

## Evidence and reproduction

Read-only inputs: the supplied `97ms_basis` directory, including its report PDF,
three BASIS files, CSV and JSON metadata. Exact saved PDF table values were
extracted with installed Ghostscript through Pixi. No LCModel fits were rerun.

- Supplied TE97 SHA256 matches saved acquisition provenance:
  `87b6b1400b2d6819af75a669523621b47a35b6fa96ee53a03915eb10d0a7ff91`.
- The CSV/seed=1 replay reproduces **every saved concentration sample exactly**,
  maximum difference 0, after the same filtering and tissue iteration as production.
- Supplied basis columns preserve their correctly named original responses after
  a common normalization; checked relative complex column-shape errors are below
  4e-15 for the inspected TE97 columns. The error occurs when applying parameters.
- The concentration table order is identical for WM, GM and CSF in this example.

Actual selected orders:

```
BASIS: Asc Asp bHG Cr GABA Gln Glu Gly GPC GSH Lac mI NAA NAAG PCh PCr PE sI Tau
CSV:   Asc Asp Cr GABA Gln Glu Gly GPC GSH Lac mI NAA NAAG PCh PCr PE sI Tau bHG
```

The bug is visible in `simulation/basis.py::rename_metabolites` (name sorting)
and `simulation/signalModel.py::simulate_metabs`: `df` retains CSV order;
`concs` and `t2s` come from it, while `fids` follows `basis.names`.
`specs = concs * spec` performs positional multiplication with no name alignment.
Saved sampled concentration dictionaries are correctly labelled according to the
CSV, so saved ground truth alone conceals the synthesis mismatch.

## Complete TE97 mapping

All values below are nominal concentration units treated as mM by this mode.
The assigned column includes the same misassigned tissue T2 values; it is not
merely a labelling error in the exported ground-truth JSON.

| Signal fitted by LCModel | Concentration label actually applied | Saved nominal | Actually assigned | PDF LCModel | Fit / assigned |
| --- | --- | ---: | ---: | ---: | ---: |
| Asc | Asc | 0.049973 | 0.049973 | 0.049300 | 0.9865 |
| Asp | Asp | 1.429047 | 1.429047 | 1.319000 | 0.9230 |
| bHG | Cr | 3.966743 | 4.873432 | 4.683000 | 0.9609 |
| Cr | GABA | 4.873432 | 1.391023 | 1.402000 | 1.0079 |
| GABA | Gln | 1.391023 | 2.498749 | 2.391000 | 0.9569 |
| Gln | Glu | 2.498749 | 7.915198 | 7.671000 | 0.9691 |
| Glu | Gly | 7.915198 | 0.308245 | 0.216000 | 0.7007 |
| Gly | GPC | 0.308245 | 0.810108 | 0.820000 | 1.0122 |
| GPC | GSH | 0.810108 | 1.086856 | 0.958000 | 0.8814 |
| GSH | Lac | 1.086856 | 0.053800 | 0.042000 | 0.7807 |
| Lac | mI | 0.053800 | 5.227995 | 5.063000 | 0.9684 |
| mI | NAA | 5.227995 | 10.008235 | 9.509000 | 0.9501 |
| NAA | NAAG | 10.008235 | 2.614968 | 2.534000 | 0.9690 |
| NAAG | PCh | 2.614968 | 0.787307 | 0.767000 | 0.9742 |
| PCh | PCr | 0.787307 | 2.138553 | 2.171000 | 1.0152 |
| PCr | PE | 2.138553 | 0.082819 | 0.028800 | 0.3477 |
| PE | sI | 0.082819 | 0.573339 | 0.678000 | 1.1825 |
| sI | Tau | 0.573339 | 1.615352 | 1.641000 | 1.0159 |
| Tau | bHG | 1.615352 | 3.966743 | 4.244000 | 1.0699 |

NAA's nominal input is being applied to **mI**, NAAG's input to **NAA**, and Cr's
input to **bHG**. The bHG signal therefore carries approximately **4.8734**, not
the recorded bHG concentration **3.9667**; its reported LCModel value is **4.683**.
The requested bHG standard deviation of 0.4 is also attached to the wrong basis
column (Tau) in this case. The issue affects the full mixture.

## Independent minimal control

A standalone control invokes the actual `SignalModel.simulate_metabs()` method,
using two known complex FIDs, zero concentration SD and no decay:

```
BASIS order:                  NAA, Cr
Table order / values:         Cr=2, NAA=8
Recorded concentration dict:  Cr=2, NAA=8
Actual synthesized weights:  NAA=1.999999996, Cr=7.999999980
Name-aligned table control:   NAA=8.000000220, Cr=1.999999987
```

The coefficients are recovered independently by complex least squares, without
LCModel or an empirical gain. Assertions pass within complex64 arithmetic tolerance.
Only the control's in-memory table is reordered; no production source is patched.
Torch RNG state is restored after both the control and seed replay.

## Why this is not established as TE97-specific

All three supplied bases put bHG in the same alphabetically sorted position.
Using the supplied bHG04.csv therefore produces this mismatch for TE35 and TE144
as well. Removing bHG from both the selected basis and this table makes their
remaining orders agree exactly. A difference in metabolite selection/table or
software path between the collaborator's TE experiments is the leading explanation
for the reported TE specificity; their other acquisition outputs are not supplied.
No claim is made that the actual TE35/144 acquisitions were reproduced.

The 20-mm cube contributes 8000 mm^3 to legacy volume summation. Metabolites and
the separate reference share that convention, so volume alone would not produce
the observed metabolite permutation or an NAA-specific factor of four.

Saved voxel water is 39.33278931 M; CONTROL uses WCONC=39330 mM. This is a
-0.00709% water-concentration difference, not the factor of four. The config's
vendor label says GE while the supplied basis is Siemens, but the loaded basis
hash matches; a metadata label mismatch does not explain these swapped inputs.
Do not change spin-system proton counts, basis amplitudes or water calibration
on the strength of the reported NAA error.

## Recommended correction

Fix the name-based join in synthesis:

1. Sample concentration means/SDs in the existing table order and preserve the
   existing name-labelled sampled dictionaries and seed behavior.
2. Align the sampled concentrations **and T2 values** to `basis.names` before
   forming FID/spectral tensors. Reject duplicate metabolite rows; handle missing
   per-tissue entries explicitly rather than relying on positional broadcasting.
3. Add row-permutation and appended-bHG regression tests that assert spectra and
   labelled ground truth agree. Include deliberately shuffled tissue rows.
4. Regenerate affected acquisitions and ground truth from the corrected model;
   rerun LCModel with the matched basis. Do not repair existing results by scaling
   NAA or relabelling alone: tissue-specific decay was also misassigned.

No production fix was applied during this investigation. Correcting this bug
intentionally changes spectra previously generated with mismatched row orders;
that change should be recorded, not disguised as calibration.

## Reproducible commands

```bash
mkdir -p outputs/te97_basis_diagnosis
pixi run gs -q -dBATCH -dNOPAUSE -sDEVICE=txtwrite \
  -sOutputFile=outputs/te97_basis_diagnosis/collaborator_report.txt \
  "$TE97_SOURCE/0001.pdf"
MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m research.validation.te97_concentration_order.diagnose \
  --source "$TE97_SOURCE" \
  --report-text outputs/te97_basis_diagnosis/collaborator_report.txt \
  --output outputs/te97_basis_diagnosis/verified_replay
```

`TE97_SOURCE` was the supplied external `97ms_basis` directory. Choose a new output
path when repeating. The numerical script was executed successfully; all built-in
assertions passed. Machine-readable source hashes, all three mappings, coefficient
control and exact seed replay are saved in `verified_replay/diagnosis.json`.
All external input hashes are unchanged. No fitting, arbitrary correction factor,
production edits, staging or commits occurred. Only the research diagnostic and
this report were added.

Limitations: the supplied directory contains no actual output.RAW, water.RAW,
NIfTI-MRS acquisitions, or machine-readable LCModel TABLE/PRINT/COORD files. The
observed PDF/metadata pattern is explained and the synthesis defect is reproduced,
but exact original fit replay and post-fix concentration recovery remain unverified.

## Name-keyed correction implemented

Legacy Basis now exposes `responses_by_name`; `names` and `fids` are derived
compatibility views, rebuilt on reset. Synthesis samples in existing CSV order,
then builds named concentration/T2 records and assembles all numerical vectors
from the same name sequence. Duplicate tissue/metabolite records and basis alias
collisions are rejected. Missing records contribute zero, with explicit metadata.
New exports record `metabolite_mapping.version = name-keyed-v1`. Historical
acquisitions are not repaired or overwritten. The old minimal diagnostic artifacts
remain evidence; running the diagnostic now verifies corrected production weights.

Regression coverage includes shuffled CSV rows with distinct T2 values, reversed
BASIS order, appended bHG, missing records, duplicate rejection before RNG draws,
fixed-seed CSV-order sample preservation, reset/alias handling and complex exports.
The obsolete AST-equality check was replaced with an actual independently executed
historical SignalModel comparison for aligned inputs. Arrays match exactly there;
new mapping metadata and missing-entry zeros are tested explicitly. The existing
notebook-import test was updated to locate the user-edited import cell by content
rather than a removed tag; the notebook itself was not modified.

A new controlled 20-mm cube uses the supplied saved tissue fractions and exactly
replayed seed-1 samples, with no B0, noise, MM/lipids or residual contamination.
Existing tissue-dependent readout decay is retained. Tissue signals are integrated
with the same 8000-mm^3 gain without allocating an anatomical cell-by-time array.
Independent named complex sums agree with corrected production to relative error
1.2810147432341715e-7, consistent with its complex64 arithmetic. Original inputs,
normalization and sequence calibration are unchanged.

Using the earlier LCModel settings, adapted only to actual TE97, present metabolite
names, supplied basis/grid and water concentration, the matched-water results are:

| Metabolite | Nominal mM | Historical positional control | Corrected control |
| --- | ---: | ---: | ---: |
| NAA | 10.008235 | 2.561 | 9.705 |
| NAAG | 2.614968 | 0.762 | 2.506 |
| bHG | 3.966743 | 4.729 | 3.927 |
| Cr | 4.873432 | 1.375 | 4.611 |
| Glu | 7.915198 | 0.265 | 7.695 |
| mI | 5.227995 | 9.813 | 5.059 |

Corrected NAA remains about 3.03% low and bHG about 1.00% low; no gain was tuned to
remove these residual fitting/quantification differences. Phased-real residuals
are 1.20336% historical versus 0.80589% corrected. Six actual LCModel fits completed
(three water/input cases per control), with PDFs and complete outputs under
`outputs/te97_name_fix_validated/fit_corrected/` and `fit_historical/`. These are
controlled new datasets, not claims of original RAW replay. The first fitting
attempts exposed a research-runner formatting error (`ECHOT=97.0.`); formatting
floating TE as valid Fortran input fixed that setup failure without changing fit
settings. Failed attempts remain preserved.

Commands executed with Pixi:

```bash
pixi run python -m pytest -q tests/test_named_metabolites.py
pixi run python -m pytest -q tests --ignore=tests/integration
pixi run python -m research.validation.te97_concentration_order.validate_fix --source "$TE97_SOURCE" --output outputs/te97_name_fix_validated
pixi run python -m research.validation.lcmodel_legacy_water.run --source outputs/te97_name_fix_validated/corrected --output outputs/te97_name_fix_validated/fit_corrected --binary "$LCMODEL_BINARY"
pixi run python -m research.validation.lcmodel_legacy_water.run --source outputs/te97_name_fix_validated/historical --output outputs/te97_name_fix_validated/fit_historical --binary "$LCMODEL_BINARY"
```

All Python commands used the existing Pixi environment (MPLCONFIGDIR directed to a
writable temporary cache). No Git staging/commits or frozen calibration changes.

Final verification: the non-fitting suite passed **305 tests** (33 existing dependency/legacy warnings, 127.39 s); the focused name-keying suite passed **5 tests**. No expected failures or skips were used.

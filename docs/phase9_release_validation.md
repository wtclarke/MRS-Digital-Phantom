# Phase 9 release verification

Scientific signal generation is unchanged. Only resource lookup, packaging/import
boundaries, examples, notebook, tests and repository policy were hardened. No
LCModel/FSL-MRS fitting, basis regeneration, staging or commits were performed.

## Executed commands and results

All Python execution used the repository Pixi environment. No environment manifest
solve or lockfile modification was performed in Phase 9. The already-present lock
and dependency changes are explicitly classified for the production commit.

```bash
pixi run python --version
# Python 3.14.8

pixi run python -m pytest -q tests/test_quantitative_production.py tests/test_quantitative_production_legacy.py tests/test_calibration_resources.py
# Initial focused run: 34 passed, 16 warnings, 59.15 s.

MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m simulation.headless --config examples/configs/quantitative/pure_WM.json --output-dir outputs/phase9_release_validation/cli_pure_WM
# Generated 1 voxel; exit 0.

MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m examples.execute_quantitative_notebook --output outputs/phase9_release_validation/notebook_execution_permitted
# Cells 9, 10, 12, 14, 16: all passed.

MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m examples.verify_production_release --output outputs/phase9_release_validation/independent_production
# Production-only tree: 36 passed, 16 warnings, 61.20 s; CLI and wheel build succeeded.

MPLCONFIGDIR=/tmp/mrs-phase9-mpl pixi run python -m pytest -q tests --ignore=tests/integration
# First broad run before the final two new tests: 267 passed, 16 warnings, 106.26 s.
# Final broad run including those tests: 269 passed, 16 warnings, 116.09 s.
```

No failures, skips or expected failures occurred in these pytest runs. The two new
controls test forbidden research/Octave imports and independent noise/SNR statistics.
Warnings come from installed Dash/Swig deprecations and the pre-existing legacy
Torch/NumPy array-wrap interaction; they are not hidden or suppressed.

The unrestricted full suite was deliberately not run: `tests/integration/` includes
fitting and Octave experiments prohibited by this phase. Fast saved-output numerical
analyses are included in the broad non-fitting command; none launches a fitter.

The first attempted notebook command used output
`outputs/phase9_release_validation/notebook_execution` and failed before executing
cells with `PermissionError: [Errno 1] Operation not permitted` while Jupyter bound
its local sockets. It was rerun with approved execution outside the sandbox at the
`notebook_execution_permitted` destination. No dependency was added. The kernel also
reported its installed TCP transport warning; local execution completed normally.

## Production-only independence

The release verifier exports the existing tracked HEAD with `git archive`, excluding
historical generated `outputs/` files and external anatomical datasets. It overlays
only the production manifest and asserts that **all 101 optional research paths
are absent**. Its 36 focused tests, actual CLI and wheel build all succeeded using
the active Pixi Python. The wheel contains all eight bundle assets, including exact
construction-source provenance, and excludes optional experimental calibration
modules. Both planning manifests belong to the production commit so its verifier
has no dependency on the optional research commit.

The tool internally ran the active Pixi Python with:

```text
-m pytest -q tests/test_quantitative_production.py tests/test_quantitative_production_legacy.py tests/test_calibration_resources.py
-m simulation.headless --config examples/configs/quantitative/pure_WM.json --output-dir outputs/independent_cli
-c "import setuptools.build_meta; setuptools.build_meta.build_wheel('dist')"
```

Complete subprocess command lines, stdout/stderr and return codes are saved in
`outputs/phase9_release_validation/independent_production/verification.json`.
The final tree has all 48 production files byte-identical to the proposed manifest
and no research-manifest files; see `final_manifest_validation.json` alongside it.
Final documentation/manifest-only updates were synchronized after execution; tested
runtime, tests, configuration and notebook code are unchanged.

The standalone quantitative notebook section was additionally executed in that tree
using the active Pixi interpreter, via a Pixi Python `subprocess.run` with cwd
`outputs/phase9_release_validation/independent_production/production_checkout`:

```text
-m examples.execute_quantitative_notebook --output outputs/independent_notebook_execution
```

All five code cells passed. The captured output is
`outputs/phase9_release_validation/independent_notebook_execution.log`; cell records
and actual acquisition files are under the production tree's ignored `outputs/`.
This is a production-only source-tree verification, not a fresh dependency download
or a clean-platform installation test. Only the available macOS arm64 environment
was executed. Linux is present in the existing lock but was not exercised.

## Calibration relocation: exact preservation

The seven original files were copied from
`outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative/` to
`calibration/mrscloud_press_te30_quantitative/` without editing any byte.
`calibration/frozen_file_hashes.json` records all original hashes. The original
files remain on disk. No response was regenerated or rescaled.

| Check | Actual result |
| --- | --- |
| Seven original file SHA256 comparisons | All exactly identical |
| Identity | `mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81` |
| Scalar | `8857.353925964408`, unchanged |
| 27 loaded metabolite unit-response arrays | Exact array equality |
| Loaded water unit response and calibration metadata | Exact equality |
| Pure WM, mixed, relaxation, spatial B0 and noisy generated pairs | Exact complex equality to original-bundle outputs |
| Ground-truth dictionaries for all five cases | Exact equality |
| Full complex export/reload error for both acquisitions in all five cases | **0** |
| Amplitude/phase changes due to relocation/export | **0**, phase considered where sample amplitude is nonzero |
| Dwell time / spectrometer frequency | 0.0005 s / 127.731 MHz, preserved |

The direct original-bundle comparison was executed with Pixi Python using
`load_library`, `simulate_quantitative`, `export_quantitative` and FSL-MRS's NIfTI-MRS
reader, with independent original/stable paths and exact `array_equal` assertions.
Its machine-readable results are
`outputs/phase9_release_validation/original_bundle_comparison/relocation_validation.json`.
The committed tests repeat hash checks, arbitrary-location resource loading and
exact paired exports without depending on the original `outputs/` directory.

`provenance_sources.zip` contains exactly 408 historical source files. Its members
are individually SHA256-checked against the unchanged original normalization record,
never imported/executed. The archive is 13,885,482 bytes; the stable asset set is
roughly 24 MiB, suitable for ordinary Git. No LFS or download mechanism was added.

## Notebook coverage and outputs

Cells 8–17 were added; code cells **9, 10, 12, 14, 16** were executed in a fresh
Pixi-backed kernel, both in the working tree and the production-only projection.
Each uses `simulation.headless.run_headless`, the same dispatcher used by the CLI.
It loads the default packaged calibration and exports pairs, verifies full complex
reload equality, dwell time, spectrometer frequency, ground-truth sidecar equality,
identity/scalar and fixed-seed noise repeatability. Three spectrum/tissue/ground-truth
figures were saved, without embedding new binary outputs in the committed notebook.

Working-tree outputs:

`outputs/notebook_quantitative_demo/20261031T175757564217Z/`

This contains `pure_WM/`, `mixed_tissue/`, `noisy/` with paired acquisitions/JSON,
plus `pure_WM.png`, `mixed_tissue.png`, `relaxation_B0_noise.png`. Noise-enabled
output also contains clean acquisitions. Each subsequent execution chooses a new
timestamped output root, preserving all prior runs.

| Case | NAA/Cr/GPC/PCh/Glu nominal voxel-average mM | Effective water mol/L |
| --- | --- | --- |
| Pure WM | 10 / 8 / 2 / 1 / 8 | 39.0165 |
| 20/30/50% WM/GM/CSF | 6.2 / 4.3 / 0.85 / 0.56 / 5.2 | 49.72245 |
| Mixed relaxation/B0/noise | Same nominal values; separate relaxation-weighted metadata | Same nominal water; separate attenuation metadata |

Mixed WM/GM-only metabolite concentrations are 12.4 / 8.6 / 1.7 / 1.12 / 10.4 mM.
The signals are never independently normalized for display or export. Water and
metabolite plots have separate axes with actual relative signal units.

Original legacy cells 0–7 and all their existing outputs were retained. Stable cell
IDs were supplied to satisfy the existing nbformat 4.5 declaration. The broken
personal Siemens TE97 basis path was replaced by the checked-in Philips PRESS TE30
basis and matching declared metadata. Historical legacy outputs are labeled as such.
**Those cells were not executed**: they retain the anatomical MRiLab/legacy lipid
requirements. Full-notebook validation is not claimed. A separate executed legacy
regression supplies the backward-compatibility evidence below.

## Legacy and scientific invariance

`tests/fixtures/legacy_pre_quantitative/` contains exact files from
`c85821a833209b45384c0a8835e5d2ead8051daf`, before Phase 8 integration. Testing no longer
requires that revision to be present in Git history. It executes the pre-feature
phantom/headless logic with fixed seed 123, controlled small anatomy and explicit
basis data. External MM/lipid providers are controlled nonzero fixtures; actual
metabolite, residual-water, shim and acquisition-noise equations execute.

Both omitted mode (default legacy) and explicit `legacy` match the baseline with
**zero difference** in total complex arrays, every nonzero component array,
concentrations, axes and relevant metadata. NumPy/Torch RNG state is restored exactly.
SignalModel AST comparison removes imports only; scientific code is unchanged.
The sole Phase 9 SignalModel edit delays MatlabRunner import until lipid generation,
so importing production quantitative/headless simulation does not initialize Octave.
This demonstrates the exercised legacy configuration, not every external anatomy,
lipid engine or platform.

The shared unit-response scaling function was moved verbatim into
`simulation/quantitative_signals.py`, with its previous research API re-exported.
All analytical/frozen-calibration non-fitting regressions pass. No scientific
normalization, proton factor, tissue weighting, relaxation equation, B0 sign,
noise stream or concentration convention was changed.

## Git hygiene and remaining limits

Explicit manifests contain **48 production** and **101 optional research** paths,
with no overlap and no new generated output paths. The research manifest includes
no production calibration assets. Both manifests and the classification document
belong to the first commit, so review/verification planning is self-contained.
No staging or commit command was executed; `git diff --cached --exit-code` passes.
`git diff --check` also passes.

New generated artifacts are ignored and preserved on disk. **4,772 old tracked
output files** remain tracked; removing them requires a separately reviewed migration,
not an automatic cleanup. No proprietary LCModel binary, manual or source is included.
Research executable/manual/source locations now use explicit environment variables.

Frozen metadata contains one historical absolute comparison path; six archived
upstream text sources retain historical example paths/comments. Existing legacy
notebook outputs also retain old user-directory warnings. These inert provenance
exceptions are documented rather than rewritten against byte-preservation and
output-preservation requirements. New runtime/configuration/notebook source paths
are repository-relative or resolved package resources.

Optional research scripts are logically separated by explicit manifests and indexes;
physical movement is deferred to preserve cross-imports and capture signatures.
No scientific model or empirical calibration correction was added. Relative sequence
response remains distinct from physical scanner calibration and fitted concentration.
Recommended next step: review the two manifests and proposed production commit;
do not launch further fitting as part of Phase 9.

# Repository artifact classification and proposed commits

Phase 9 separates executable production from optional scientific experiments without
changing the signal model. **No staging, committing, history rewriting or deletion
of existing output files is performed.** Explicit non-overlapping path manifests:

- [Production commit: every path](commit_manifests/production.txt).
- [Optional research commit: every path](commit_manifests/research.txt).

The manifests cover current Phase 2–9 additions/modifications. Existing tracked
baseline files remain in the repository; a manifest does not re-add unchanged files.
The production file list includes packaging/environment changes already present
before Phase 9, because the frozen formatting-version checks require that environment.
Phase 9 does not solve dependencies or alter the Pixi lockfile.

## Classification

| File or directory | Category | Runtime requirement | Git policy | Current location | Proposed location | Migration action | Reason |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `simulation/quantitative.py`, `headless.py`, `mrs_phantom.py` | A production | Dispatch, config, batch/ROI, export | Commit 1 | `simulation/` | Same | Keep | Supported public headless and phantom entry points |
| `simulation/mrscloud_full_basis.py`, `quantitative_signals.py` | A production | Load/validate frozen responses; concentration scaling | Commit 1 | `simulation/` | Same | Keep; pure scaling moved verbatim into production module | Full-library runtime no longer imports experimental calibration adapter |
| `simulation/mrscloud_quantitative_voxel.py`, `mrscloud_quantitative_spatial.py`, `tissue_relaxation.py`, `mrscloud_quantitative_noise.py` | A production | Validated tissue, relaxation, B0 and noise equations | Commit 1 | `simulation/` | Same | Keep | Direct runtime import closure |
| `simulation/calibration_resources.py` | A production | Stable lookup, historical source-byte checks | Commit 1 | `simulation/` | Same | Add | Removes dependence on experimental output/source locations without dropping hash checks |
| `simulation/signalModel.py` | A legacy production | Existing legacy signal equations | Commit 1 import-only change | `simulation/` | Same | Delay MatlabRunner import until lipid generation | Avoid Octave startup for quantitative imports; scientific AST unchanged after imports |
| `gui/defaults.py`, `utils/definitions.py`, loading/auxiliary/NIfTI helpers, `simulation/basis.py`, legacy signal components | A existing production | Shared transitive dependencies | Already tracked baseline | Existing packages | Same | No changes | Reused by legacy imports, anatomy and I/O; no research dependency |
| Seven original Phase 6a calibration files | B required asset | Complete basis, raw water, reference arrays/metadata | Commit 1 ordinary Git | `outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative/` (original retained) | `calibration/mrscloud_press_te30_quantitative/` | Byte-copy, never regenerate | Stable production location; frozen identity/scalar/array bytes unchanged |
| `calibration/.../provenance_sources.zip` | B required provenance | Verify 408 historical source hashes | Commit 1 ordinary Git | New archive | Stable calibration bundle | Archive exact original bytes; never execute | Replaces dependence on live research sources and ignored spin-system captures |
| `calibration/__init__.py`, README, `frozen_file_hashes.json` | B/C resources and integrity | Package resource discovery/test record | Commit 1 | `calibration/` | Same | Add | Package data and original-seven-file hash record |
| `tests/test_quantitative_production.py`, `test_quantitative_production_legacy.py`, `test_calibration_resources.py`, `tests/fixtures/legacy_pre_quantitative/` | C production tests | Release checks only | Commit 1 | `tests/` | Same | Remove historical-output/Git-history dependencies; retain fixed baseline sources | Independent production coverage without fitting binaries, external anatomy or research captures |
| `tests/test_analytical_water_reference.py`, `test_analytical_voxel_reference.py`, `test_nifti_mrs_roundtrip.py` and analytical runtime helpers | C existing regression | Existing supported analytical/legacy paths | Already tracked baseline | `tests/`, `simulation/` | Same | Preserve | No analytical behavior changes |
| `notebooks/headless_simulation.ipynb` | C demonstration | Optional interactive use | Commit 1 | `notebooks/` | Same | Append cells 8–17; retain legacy cells and outputs | New standalone section calls production API, no fitting/regeneration |
| `examples/configs/quantitative/*.json`, `quantitative_water_reference_production.py`, `execute_quantitative_notebook.py`, `verify_production_release.py` | C usage/verification | Production demos/release tooling | Commit 1 | `examples/` | Same | Add supported examples | No dependency on research examples or their output directories |
| User guide, Phase 8 integration report, Phase 9 release report/classification/manifests, AGENTS.md | C production documentation | Usage and release policy | Commit 1 | `docs/`, root | Same | Add/update | Core installation/use independent of scientific archive |
| `pyproject.toml`, `pixi.lock`, `.gitignore` | A/C packaging/config | Reproducible environment, resources, artifact policy | Commit 1 | Root | Same | Package calibration; retain validated existing lock; ignore new outputs | No new dependency needed for notebook execution |
| `simulation/frozen_basis_normalization.py`, `mrscloud_water_calibration.py`, `mrscloud_quantitative_pair.py`, `mrscloud_sequence_reference.py`, `water_reference_preprocessing.py` | D research/optional analytical comparison | Not imported by production full-library mode | Commit 2 | `simulation/` | Retain for now | Index and classify; defer physical relocation | Existing experimental cross-imports/Phase 4–5 APIs must stay reproducible; calibration adapter re-exports pure production scaling |
| `examples/sequence_matched_water.m`, `capture_water_processing.m`, `capture_naa_response.m`, `full_basis_responses.m`, associated Python experiment drivers | D generation | No production runtime requirement | Commit 2 | `examples/` | Indexed by `research/calibration_generation/` | Keep paths; documented later migration | Frozen source signatures/cross-imports refer to these historical paths |
| Remaining scientific `examples/mrscloud_*.py`, `fsl_mrs_*.py`, `lcmodel_*.py` | D validation | No production runtime requirement | Commit 2 | `examples/` | Indexed by `research/validation/` | Preserve cross-imports; replace personal LCModel defaults with explicit environment configuration | Optional research outputs and fitters are prerequisites only for rerunning research |
| `utils/lcmodel_quantitative_io.py`, `lcmodel_noise_analysis.py` | D research helpers | No production imports | Commit 2 | `utils/` | Retain for now | Document dependency boundary | Independent fitting/conversion analysis; not needed by headless simulation |
| Other new `tests/test_*.py`, `tests/integration/*` listed in research manifest | D optional tests | May need historical captures, fitting or Octave | Commit 2 | `tests/` | Retain; indexed by research | Keep tests with source workflows; run explicitly | Production release suite excludes these external prerequisites |
| Scientific Phase 4–7 reports in research manifest | D scientific documentation | Optional provenance | Commit 2 | `docs/` | Same; research index | Keep discoverable; remove personal tool paths from prose | Does not impose research dependencies on normal use |
| `research/README.md`, `calibration_generation/README.md`, `validation/README.md` | D organization | None | Commit 2 | `research/` | Same | Add indexes and migration instructions | Clear separation without breaking historical imports |
| All generated fitting/simulation outputs, HTML/PDF/PS/PRINT/COORD/TABLE/RAW/H2O, NIfTI acquisitions, NPZ/MAT captures, plots and logs | E generated | Not needed except designated frozen calibration | Exclude new files | `outputs/` | Same | Preserve on disk; ignore | Regenerable/local artifacts; calibration is the explicit immutable exception copied out |
| `.pixi/`, caches, notebook checkpoints, build/egg-info | E generated | Local environment/build | Ignore | Existing | Same | Preserve/ignore | No source or scientific asset requirement |
| Licensed LCModel binary/manual/source | E external proprietary | Never production dependencies | Never include | User-managed external installation | Remains external | No copying | Research uses explicit environment variables, no download or redistribution |

## Migration, integrity and ambiguities

The seven original files total about 10 MiB. The added provenance archive is
13,885,482 bytes; total stable assets are about 24 MiB. Ordinary Git is appropriate;
Git LFS is neither configured nor necessary. Packaging includes all bundle files.
The archive retains original MRSCloud/FID-A source attribution and is inert data,
not a research runtime dependency. It contains no proprietary LCModel material.

**Immutable provenance exception:** the frozen `reference_manifest.json` retains
one historical absolute comparison path. It is not used at runtime. Removing it
would change the manifest hash returned in scientific metadata, conflicting with
the explicit byte-preservation requirement. All new code/configuration/notebook
source paths are portable. Six archived upstream source files retain historical
example paths/comments, and preserved legacy notebook output text contains old
user-directory warnings. These are inert historical content, not new execution
defaults. These byte-preservation exceptions remain visible during review.

There are **4,772 historically tracked generated files under `outputs/`**. Ignore
rules cannot untrack them. They were not deleted or silently removed from the
index. A later separately authorized cleanup can inspect them, archive any unique
provenance, and remove tracked artifacts in a dedicated migration commit. Neither
proposed commit stages new output files or changes those existing tracked files.

Physical research relocation is intentionally deferred, not implied completed.
Several experiments import each other under `examples`; the old generation source
paths also participate in frozen signatures. The research indexes and explicit
manifest provide logical separation now. A future move must update the entire
cross-import group and its optional tests, and record a new research configuration
without changing production calibration metadata.

The runtime import graph now terminates in production modules, FSL-MRS basis/I/O
libraries, NumPy/Torch/TorchIO, nibabel/nifti_mrs and existing anatomy/utilities.
Importing the quantitative headless path does not initialize oct2py. Octave remains
in Pixi for legacy lipid simulation/research; it is not removed speculatively.
Existing FSL-MRS development-channel/pytest lock changes predate Phase 9 and remain
part of Commit 1 because the calibration enforces exact formatter versions.
Installation on the second locked platform has not been executed in this session.

## Proposed staging and commit messages — not executed

Commit 1:

```bash
git add --pathspec-from-file=docs/commit_manifests/production.txt
# Proposed message: feat: add production sequence-matched quantitative MRS simulation
```

Commit 2, after reviewing Commit 1:

```bash
git add --pathspec-from-file=docs/commit_manifests/research.txt
# Proposed message: research: archive MRS quantitative calibration and validation workflows
```

No `git add`, `git commit`, push or history rewrite was run. Both planning manifests belong to Commit 1 so its verification tool is self-contained.
Each manifest lists individual paths, not a broad directory glob. They are disjoint, and the research
manifest contains no `calibration/` paths. Commit 2 depends on production helpers
from Commit 1; Commit 1 has no reverse dependency.

## Verification that Commit 1 stands alone

`examples/verify_production_release.py` exports the tracked baseline, excluding
historically generated outputs and external anatomical data, overlays **only** the
production manifest, and asserts that every research-manifest path is absent.
It runs the production suite and the actual headless CLI using the active Pixi
Python, then builds a wheel and checks calibration inclusion and research exclusion.
It writes no Git index/history data. This verifies a proposed tree rather than
creating commits or claiming a clean-install/network test.

```bash
pixi run python -m examples.verify_production_release \
  --output outputs/phase9_release_validation/independent_production
```

The same environment's existing dependencies are reused; no installation/download
is performed. The release report records actual counts, numerical comparisons and
notebook coverage. A reviewer can additionally create a clean checkout after
manually committing Commit 1, run `pixi install`, the documented CLI, the notebook
section and the production regression command without Commit 2.

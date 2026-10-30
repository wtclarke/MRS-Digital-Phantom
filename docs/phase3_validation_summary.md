# Phase 3 validation summary

## Static isolation assessment

The Phase 3 production addition is `simulation/analytical_voxel_reference.py`.
It imports only NumPy, the NIfTI-MRS generator, and the concentration-unit
converter from `simulation/analytical_water_reference.py`. Search of the
legacy entry points (`DigitalPhantom`, `SignalModel`, `simulate_batch`, and
`simulate_data`) found no import, call, or initialization of either analytical
module. The analytical functions do not set seeds, draw random values, or
modify global state. Their import therefore adds no code path or RNG action to
legacy simulation. No legacy source, defaults, configuration, GUI, MRSCloud
basis code, or exporter implementation was changed specifically for this
feature. The Pixi manifest/lockfile are already modified in the working tree;
this assessment did not change them.

This is static isolation evidence, not numerical legacy equivalence.

## Historical baseline and remaining verification

The available parent revision `0adf4f8` is not a clean pre-feature baseline:
diff inspection shows unrelated prior changes in `simulation/mrs_phantom.py`
and `simulation/signalModel.py` (including logging/verbosity, label caching,
and concentration reporting). No verified legacy output fixture was found.
Per task scope, historical baseline search stopped here. No fixed-seed legacy
simulation comparison was run, so complex signals, components, concentrations,
metadata, and random-stream equivalence have not been numerically verified.

No lightweight legacy simulator smoke test was added because constructing a
valid run requires the repository's external phantom/basis data and a clean
baseline reference; fabricating expected outputs would not be meaningful.

## Pixi test results

- `pixi run python -m pytest -q tests/test_analytical_voxel_reference.py` — **6 passed**.
- `pixi run python -m pytest -q tests/test_analytical_water_reference.py` — **8 passed**.
- `pixi run python -m pytest -q tests/test_nifti_mrs_roundtrip.py` — **4 passed**.

The Phase 3 suite covers pure WM/GM/CSF water responses, independent mixed
tissue sums and ratio, deterministic repeated output, and paired NIfTI-MRS
reload. The Phase 2b suite confirms both validated 1D exporter round trips and
records current image-shaped behavior. Test warnings are dependency deprecation
and optional progress-display warnings; all commands exited successfully.

## Conclusion and limitations

The analytical feature is statically isolated from legacy simulation and the
requested analytical/export suites pass. This supports the claim that the new
path is opt-in and does not invoke or alter legacy behavior when unused.
Numerical backward compatibility is **not established** because a clean
pre-feature baseline is unavailable. The analytical outputs remain synthetic,
not physically calibrated; MRSCloud calibration, relaxation, and spatial
integration remain outside this validation.

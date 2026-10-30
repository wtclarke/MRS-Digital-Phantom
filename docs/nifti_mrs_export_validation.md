# NIfTI-MRS export round-trip validation

## Method

Pixi Python 3.14.8 was used. A 256-sample complex FID was generated as
\(s(t)=e^{-t/0.085}e^{i(2\pi 37t+35^\circ)}\), with dwell time 0.5 ms and
spectrometer frequency 123.2 MHz. Its spectrum was generated with repository
`fid2spec`. Both exporter paths wrote NIfTI-MRS; raw stored samples were read
with nibabel, and interpreted data were read with `mrs_io.read_FID(...).mrs().FID`.

## Results

Command `pixi run python -m pytest -q -s tests/test_nifti_mrs_roundtrip.py`:
**4 passed**. Both 1D paths preserve complex samples through the NIfTI-MRS
reader with identical measured errors:

- Relative complex-array error: `2.4367957221621032e-08`.
- Maximum absolute amplitude error: `1.1920928955078125e-07`.
- Maximum phase error: `4.213424249058947e-08` radians.
- Relative error against the conjugated interpreted FID: `1.4480114054527462`.
- Stored sample order and offset phase advance are preserved.
- Dwell time is preserved as 0.0005 s; spectrometer frequency is `[123.2]` MHz.

The raw stored arrays differ by exporter: `save_nifti_mrs` explicitly
conjugates the FID before writing, while `results_to_nifti_mrs` passes it to
`gen_nifti_mrs`, whose default convention also results in conjugated stored
samples. The NIfTI-MRS reader interprets both as the original FID. Therefore
the effective external convention is conjugated storage with reader-side
interpretation restoring the signal phase.

Image-shaped behavior was also tested. `save_nifti_mrs` accepts a `(2,2,2,N)`
spectrum and writes that same shape, without adding a dynamic dimension.
`results_to_nifti_mrs` adds the singleton dynamic axis to its supplied
image-shaped spectrum and the current helper reaches an 8D array; the installed
NIfTI-MRS generator rejects this because its supported maximum is 7D. The test
records this existing failure as an expected `ValueError`, rather than
asserting unsupported headless image export.

The required non-verbose command `pixi run python -m pytest -q
tests/test_nifti_mrs_roundtrip.py` passed: **4 passed**. The analytical
regression command `pixi run python -m pytest -q
tests/test_analytical_water_reference.py` passed: **8 passed**.

## Recommendation

For 1D spectra, either exporter currently preserves the interpreted complex
FID to complex64 precision. Prefer `gen_nifti_mrs`/`results_to_nifti_mrs` as
the canonical convention because it delegates storage and metadata creation to
the NIfTI-MRS library; explicitly retain its default conjugation behavior and
keep the round-trip test as the contract. Do not remove conjugation based on
raw-array inspection alone. Before using `results_to_nifti_mrs` with spatial
data, its dimensionality handling needs a separate production fix and tests.

These tests establish serialization consistency, not absolute signal
calibration. Exporter implementations were not modified.

# Validation workflows (optional)

Existing `examples/fsl_mrs_*.py`, `examples/lcmodel_*.py` and scientific
`examples/mrscloud_*.py` are retained in place to preserve cross-imports and historical
output references. Their optional tests remain in `tests/` and `tests/integration/`;
the explicit research commit manifest identifies them individually.

Scientific reports remain in `docs/` for discoverability: frozen normalization,
synthetic water calibration, sequence response and processing, NAA decomposition,
sampling/linewidth, tissue/relaxation/spatial validation, and fitting/noise benchmarks.
Production users should start with `docs/quantitative_simulator_user_guide.md`.

Fast numerical tests and saved-output analyses do not necessarily launch a fitter.
Integration tests may launch FSL-MRS, licensed LCModel or Octave, and need external
research outputs; they are excluded from production verification. Keep those
prerequisites explicit. Do not substitute regenerated noise or missing captures.
Generated HTML/PDF/PS/RAW/NIfTI/capture/Monte Carlo artifacts stay under ignored
`outputs/`, never in the optional research commit.

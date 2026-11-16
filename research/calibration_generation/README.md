# Calibration generation workflows (optional)

The historical entry points remain under `examples/`:

- `mrscloud_sequence_matched_water.py` / `sequence_matched_water.m` (Phase 5e).
- `mrscloud_water_processing_validation.py` / `capture_water_processing.m` (5f).
- `mrscloud_naa_response_decomposition.py` / `capture_naa_response.m` (5i).
- `mrscloud_full_basis_water_reference.py` / `full_basis_responses.m` (6a).

They use bundled MRSCloud/FID-A source and the configured Octave environment.
Do not use them to regenerate or modify the production bundle during normal use.
Other phases reference their captures and may need historical ignored output caches.
See the corresponding `docs/mrscloud_*.md` scientific reports for recorded commands.
A future physical relocation must update cross-imports and capture source signatures
as a new research configuration, without rebinding existing frozen calibration.

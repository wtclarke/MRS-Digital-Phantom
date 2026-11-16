# Optional scientific research archive

Production simulation is documented in [the user guide](../docs/quantitative_simulator_user_guide.md).
It imports no research modules and requires no generated `outputs/` data.

The optional research commit retains existing paths deliberately: scientific scripts
import each other under `examples`, and frozen provenance records their historical
paths. Moving them wholesale would create avoidable reproducibility risk. The two
indexes below separate their purpose without breaking those imports. A later migration
can move groups together, update imports and fixture paths, and run the corresponding
optional integrations. Archived frozen source bytes must never be rewritten.

- [Calibration generation](calibration_generation/README.md): MRSCloud/Octave work.
- [Validation](validation/README.md): analytical, FSL-MRS and LCModel experiments.
- [File classification and explicit commit manifests](../docs/repository_artifact_classification.md).

Existing generated outputs are local research prerequisites, not committed assets.
No fitting runs are part of Phase 9. For optional LCModel research, configure your
licensed executable using `LCMODEL_BINARY`; normalization diagnostics additionally
require explicit `LCMODEL_MANUAL` and `LCMODEL_SOURCE`. None is bundled or downloaded.

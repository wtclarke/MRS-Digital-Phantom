# Optional BASIS-H2O read-only validation

Run from the repository root using Pixi and an explicitly supplied external BASIS:

```bash
pixi run python -m research.validation.basis_h2o.validate_basis_h2o \
  --basis "$MRS_BASIS_H2O" --output outputs/basis_h2o_validation/new_audit
```

Optional `--source-mat "$MRS_BASIS_SOURCE_MAT"` compares all supplied BASIS columns
with a companion MRSCloud assembly after one declared 3-to-4.68 ppm frame change,
using its recorded RF carrier. No amplitude factors are estimated or applied.
The validation does not run fitting, Octave or basis generation, and never copies
or modifies the supplied data. Output directories must be new.

Production runtime and tests do not import this research module or need its outputs.
See `docs/basis_h2o_water_reference.md` for the supported production interface,
recorded GE results and scientific limitations.

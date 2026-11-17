# Optional BASIS/MAT provenance comparison

Compare all supplied BASIS columns with a companion MRSCloud assembly after the declared 3-to-4.68 ppm frame change, using its recorded RF carrier. No amplitude factor is estimated or applied. This preserves the independent provenance comparison; the standalone quantitative simulation audit has been removed.

Run from the repository root with external input files and an existing output parent directory:

```bash
pixi run python -m research.validation.basis_h2o.validate_basis_h2o \
  --basis "$MRS_BASIS_H2O" --source-mat "$MRS_BASIS_SOURCE_MAT" \
  --output outputs/basis_h2o_validation/new_comparison.json
```

The output must be a new JSON file. The comparison does not modify its inputs or run fitting, Octave, basis generation or simulation. The historical `outputs/basis_h2o_validation/mat_common_scale_comparison.json` is retained as provenance evidence; matching source responses does not prove molar/proton calibration or scanner accuracy.

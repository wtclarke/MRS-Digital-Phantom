# Optional notebook validation

- [BASIS/MAT provenance comparison](basis_h2o/README.md): compare supplied responses after a declared frequency-frame change, without amplitude adjustment.
- [Legacy-water LCModel conformance](../../docs/lcmodel_legacy_water_conformance.md): fit and summarize saved notebook acquisitions with `lcmodel_legacy_water/run.py` and `summarize.py`.
- [TE97 concentration-order diagnosis](../../docs/te97_concentration_order_diagnosis.md): replay and validate the metabolite-name mapping fix with `te97_concentration_order/diagnose.py` and `validate_fix.py`.

These workflows require the identified external inputs or saved notebook outputs. LCModel runs require a licensed executable supplied through `--binary`. Historical reports describe the runs that produced them; they do not certify current source execution. Production simulation does not import this research code.

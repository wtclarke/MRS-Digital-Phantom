# Phase 4b: synthetic MRSCloud-to-water response

## Source findings

The checked-in MRSCloud FID-A `simulationTools/metabolites/spinSystems.mat`
contains `sysRef0ppm` with one shift at 0 ppm, zero coupling and scaleFactor 1;
`sysH2O` has one shift at 4.7 ppm, zero coupling and scaleFactor 2.
`sysNAA_Acetyl_only` explicitly contains three equal shifts at 2.008 ppm and
scaleFactor 1. These structures were inspected using SciPy `loadmat` in Pixi.

Verified code behavior:

- `functions/FID-A/simulationTools/sim_Hamiltonian.m` multiplies the initial
  density matrix by each spin system's scaleFactor (lines 102–103).
- `sim_readout.m` calculates the magnetization trace with Fx+iFy, scales each
  subsystem by 2^(2-nspins), and sums subsystem FIDs. Proton multiplicity can
  therefore be encoded through explicit spins and scaleFactor. It must not be
  added a second time to metabolite basis columns.
- `fit_makeBasis_MRSCloud.m` copies metabolite FIDs, applies processing including
  DC correction, and divides the entire assembled basis by a common maximum
  real spectral value (lines 561–563). No separate unit-peak normalization of
  each metabolite was found in this inspected assembly/export path. Separate
  MM/lipid synthesis and scaling are present; this conclusion does not apply
  indiscriminately to those components.
- `io_writelcmBASIS.m` frequency-shifts and conjugates FIDs, applies an FFT,
  and writes each entry. These operations preserve relative scaling under
  common sampling, up to processing and output precision. `CONC=1`, `TRAMP=1`,
  and `VOLUME=1` are bookkeeping and do not demonstrate molar calibration.

The code supports proton response being encoded, but does not establish that
the checked-in file is a traceable calibrated simulation of unit molarity.
The exact historical generator configuration and sequence-specific response
of each spin have not been independently reproduced. A reference outside the
excitation passband may respond differently from NAA or water. Thus the
reference below is deliberately a synthetic convention.

## Definition

Use the unchanged Phase 4a frozen columns B_m[n] at a specified sampling grid.
Let q = abs(B_Ref0ppm[0]) for the **same full reference file and frozen scalar**.
Define an analytical water singlet w[n] = 2q for all acquired samples, with
zero offset, zero phase, no relaxation and no noise. This sets one nominal
mol/L of the analytical single-proton reference to q SRU and one nominal
mol/L of water molecules to 2q SRU. No additional proton multiplier is applied
to any metabolite column.

For concentration C_m and water concentration C_w in nominal mol/L per total
tissue volume and shared nonnegative receiver multiplier G:

\[
s_m[n]=G C_m B_m[n],\qquad s_w[n]=G C_w 2q.
\]

Metabolite responses are exactly the Phase 4a arrays. Water is constructed
analytically; the existing MRSCloud H2O column is not substituted for it. The
reference's initial magnitude defines amplitude only: it deliberately does
not transfer the reference frequency, phase, decay, or localization response.
There are no per-metabolite correction factors and no fitting results used
to set constants. Reference identity and grid are independent of selection.

For a known basis coefficient, recovery is

\[
C_m=\frac{s_m[0]}{s_w[0]} C_w\frac{2q}{B_m[0]}.
\]

This is a coefficient-consistency check, not independent validation of a
spectral fitting pipeline or absolute concentration accuracy. Future fitting
validation must use held-out responses and document mismatch sensitivity.

## Implementation and metadata

`simulation/mrscloud_water_calibration.py` provides `calibrated_responses`
and `concentration_signals`. Both are opt-in standalone functions. The former
loads the frozen basis, requires Ref0ppm, and returns selected unit-mol/L
metabolite arrays, water array, and JSON-serializable calibration metadata.
The latter applies explicitly supplied mol/L concentrations and one receiver
factor. No exporter or production SignalModel integration was added.

Metadata include convention version, physical-calibration flag (false), source
file identity/hash, frozen scalar/statistic, grid, reference name, initial
real/imaginary reference samples, q, assumed reference/water proton counts,
water frequency and phase, concentration units, and relaxation exclusion.
Concentrations and G are supplied at signal generation and must be retained by
callers alongside this metadata when exporting a simulation configuration.

Phase 4a compatibility rejects mismatched source identity, bandwidth, or sample
count. Optional expected synthetic metadata must match the regenerated record
exactly, including the convention and constants; selection is deliberately
excluded. Unknown or duplicate selections are rejected. A zero/nonfinite
reference amplitude is rejected. No constants are obtained from residual-water
settings or a fitting pipeline.

## Actual numerical results

Reference: `LCModel_Universal_Philips_UnEdited_PRESS_TE30.BASIS`, SHA-256
`bac0042c2e9a688ef63abd88867a79fd24c659ba92da5c6a333e71f96ffa5abe`.
Sampling: 2000 Hz, 2048 samples, 0.0005 s dwell, PRESS TE30 as labeled by the
file. Frozen scalar: 6.573456836083183e-05.

- q = 0.77174917937521237 SRU per nominal mol/L of one-proton reference.
- NAA unit-concentration initial response = 4.406162020098117 −
  0.5345769841144177i SRU per nominal mol/L.
- At C_NAA = 0.010 mol/L, C_w = 40 mol/L and G = 1, water initial amplitude
  = 61.73993435001698 + 0i SRU.
- Independently calculated NAA/water initial complex ratio =
  0.0007136648372702561 − 0.00008658528547889047i.
- Recovered concentration = 10 − 1.0112437029688739e-16i mM.
- NAA full FID and water array differences across NAA-only, NAA+Cr,
  NAA+Cr+GPC+PCh and all-entry selections: exactly zero, including metadata.
- Doubling either concentration produces exact doubling; tripling G triples
  both FIDs exactly and preserves their ratio within relative tolerance 2e-15.

Expected ratios are calculated from an independent direct FSL-MRS reader call,
the full-reference mean magnitude, and the Ref0ppm and NAA columns, bypassing
the calibration utility. Full complex-array ratios are compared at rtol 2e-15.

Commands actually executed:

```
pixi run python -m pytest -q -s tests/test_mrscloud_water_calibration.py
pixi run python -m pytest -q tests/test_mrscloud_water_calibration.py tests/test_frozen_basis_normalization.py tests/test_analytical_voxel_reference.py tests/test_analytical_water_reference.py tests/test_nifti_mrs_roundtrip.py
```

Results: Phase 4b 8 passed; combined suites 31 passed (8 + 5 + 6 + 8 + 4),
zero failures. Dependency deprecation and optional display warnings remain.

The fixed-seed compatibility test runs the actual legacy `simulate_water`
method on a small in-memory basis/time grid. It verifies exact complex-array
equality and identical subsequent NumPy random draws with and without an
intervening calibration call. Changing its residual-water amplitude setting
does not change the new water response. Legacy production sources are
unchanged. This demonstrates local residual-water/RNG isolation; a historical
end-to-end legacy phantom equivalence claim is still outside the available
baseline evidence described in Phase 3.

## Scientific limits

This convention does not establish physical scanner calibration or a universal
water-to-metabolite sensitivity. Shared receiver gain, per-molecule basis
response, excitation/localization, spin-system accuracy, and source-file
processing remain physical uncertainties. Frozen normalization cannot recover
missing acquisition calibration. No extra TE attenuation, T1/T2 modeling,
production integration or GUI behavior was introduced. Phase 4b stops here.

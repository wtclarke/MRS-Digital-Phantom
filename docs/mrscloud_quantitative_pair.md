# Phase 5a: homogeneous synthetic MRSCloud/water pair

`simulation/mrscloud_quantitative_pair.py` provides the opt-in `generate_pair`
and `export_pair` functions. The supported reference is the checked-in
`LCModel_Universal_Philips_UnEdited_PRESS_TE30.BASIS` (SHA-256
`bac0042c2e9a688ef63abd88867a79fd24c659ba92da5c6a333e71f96ffa5abe`).
Other file identities are rejected because their frequency origins have not
been verified. Optional Phase 4a calibration and Phase 4b expected metadata are
checked by the existing utilities; different metabolite selections use the
same complete-reference scalar. No legacy simulator integration is involved.

## Signal convention

For frozen columns B_m[n], input C_m in mM, C_w in mol/L, and common receiver
factor G, the returned metabolite signal is

\[
s_m[n]=G\sum_m(10^{-3}C_m)B_m[n].
\]

Let q = abs(B_Ref0ppm[0]) as in Phase 4b and t_n = n/BW, starting at zero.
Water is

\[
s_w[n]=2qGC_w\exp(-\pi\Delta f_wt_n)
\exp(i(2\pi f_wt_n+\phi_w)).
\]

Linewidth is nonnegative FWHM in Hz, frequency is in Hz, phase is in radians.
Phase rotates the initial complex sample; its magnitude remains 2qGC_w.
There is no added metabolite attenuation, T1/TR correction, noise, residual
water, tissue mixing, spatial weighting, or additional water amplitude factor.
Both FIDs use the same complex grid and repository `fid2spec` FFT convention.

## Frequency reference evidence

The checked-in MRSCloud exporter shifts to 4.68 ppm before conjugating and
writing the BASIS spectrum (`io_writelcmBASIS.m::shift_centerFreq`). FSL-MRS
reads the spectral entries, IFFTs, and conjugates them (`readLCModelBasis`).
The resulting FIDs follow positive complex phase advance for frequencies
defined by f_hz = (chemical_shift_ppm − 4.68) × cf_MHz.

At BW=2000 Hz, N=2048, cf=127.731 MHz, measured FFT maxima are −597.65625 Hz
for Ref0ppm, +2.9296875 Hz for H2O and −340.8203125 Hz for NAA's main peak.
Ref0ppm and H2O match their expected frequencies within one FFT bin (0.9765625
Hz). Thus the default water at 4.70 ppm is +2.554620000000059 Hz. An explicit
water offset overrides this value relative to the same 4.68 ppm origin.

Both exported datasets include standard `SpecFreqChemShift=4.68` and EchoTime
0.030 s (basis sequence label, with no extra scalar echo attenuation). Required
dwell time, frequency, nucleus, and dimension fields come from `gen_nifti_mrs`.
The library's default conjugated storage convention is retained, and FSL-MRS's
interpreted FID recovers the original array. Custom provenance remains in JSON.

## Reproducible example

From the repository root, choose a new output directory and run:

```
pixi run python -m examples.mrscloud_quantitative_pair /private/tmp/phase5a_quantitative_pair_20261030
```

The example was executed at that path and writes:

```
metabolites.nii.gz
water_reference.nii.gz
ground_truth.json
simulation_config.json
```

It uses NAA 10 mM, water 40 M, G=1, BW=2000 Hz, 2048 points, water linewidth
5 Hz, default water frequency and phase 0. Existing outputs are refused.
Measured results:

- Frozen NAA initial response: 4.406162020098117 − 0.5345769841144177i.
- Water initial amplitude: 61.73993435001698 SRU.
- Initial NAA/water ratio: 0.0007136648372702563 − 0.00008658528547889048i.
- Recovery using ratio × C_w × 2q / B_NAA[0]:
  10 + 9.287204115879089e-17i mM.
- Maximum absolute complex reload error: 0 for both example FIDs at the
  generator/library default complex128 storage precision.

Ground truth preserves the explicit concentration values/units and synthetic
status. Simulation configuration preserves the source hash, frozen scalar,
Phase 4b convention/constants, receiver factor, grid, frequency origin and
sign convention, linewidth/frequency/phase, and excluded effects. Nested
Phase 4b water phase/frequency describe its unit response before this generator
applies the water envelope; the top-level configuration records the final
acquisition's phase/frequency.

## Validation

Phase 5a tests independently read and normalize all columns through FSL-MRS,
then check single/multiple metabolite signals, selection independence,
linearity, gain cancellation, water envelope and initial amplitude, measured
frequency origin, both raw/interpreted complex exports, sidecars, incompatible
metadata, invalid inputs, and legacy residual-water/RNG isolation.

The independent time construction multiplies by dwell time while production
divides by bandwidth. These equivalent forms accumulate slightly different
floating-point phases: the initial strict 4e-15 comparison failed at a maximum
absolute difference 9.58451107e-14. The final complex signal check uses
rtol=1e-13, atol=1e-14 and the separate magnitude-envelope check uses rtol=1e-14.
No signal formula was changed to satisfy the test. Export comparisons use
rtol/atol 2e-7 to allow complex64 precision even though default complex128
storage was exact. A test with 0.6 rad initial water phase measured maximum
phase differences around 5.1e-17 radians through reload.

Commands:

```
pixi run python -m pytest -q -s tests/test_mrscloud_quantitative_pair.py
pixi run python -m pytest -q tests/test_mrscloud_quantitative_pair.py tests/test_mrscloud_water_calibration.py tests/test_frozen_basis_normalization.py tests/test_analytical_voxel_reference.py tests/test_analytical_water_reference.py tests/test_nifti_mrs_roundtrip.py
```

Final combined result: **47 passed, zero failures** (16 Phase 5a + 8 Phase 4b +
5 Phase 4a + 6 Phase 3 + 8 analytical + 4 export). Five existing dependency
deprecation/optional-display warnings were emitted. The first standalone run
and intermediate combined run exposed the floating-point tolerance issues
described above; the final combined command executed every relevant suite
successfully after the test tolerances were corrected.

## Limits

This is an internally consistent synthetic convention. Ratio-based recovery
is an arithmetic consistency check and does not independently validate a
fitting pipeline or physical molar sensitivity. Proton response, localization,
sequence-specific effects and scanner calibration remain unresolved. Only the
verified PRESS TE30 source identity is currently accepted; no production,
tissue, relaxation, or GUI integration is included. Existing local
residual-water/RNG isolation tests do not establish historical full-phantom
numerical equivalence.

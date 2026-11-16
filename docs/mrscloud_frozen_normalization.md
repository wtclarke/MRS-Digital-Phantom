# Frozen MRSCloud basis normalization (Phase 4a)

This is an opt-in **synthetic** response convention. It does not change
`Basis.reformat()` and does not establish a physical water-to-metabolite ratio.

## Definition and reference

For the full reference basis formatted at bandwidth (BW) and point count (N),
the frozen scalar is

\[
k_{ref}=\frac{1}{N M}\sum_{n=0}^{N-1}\sum_{m=1}^{M}|B_{ref}[n,m]|.
\]

Every selected basis column is (B_m[n]/k_{ref}). Selection occurs after
formatting and scalar calculation, so adding or removing other metabolites
cannot change NAA's response. Concentration scaling is separate: for a nominal
concentration (C) mol/L the prototype response is (C B_m/k_{ref}).

Reference identity is the basis filename plus SHA-256 of the source file.
Calibration metadata also records statistic, scalar, bandwidth, and point
count. Reuse rejects different identity, bandwidth, or point count. Metabolite
selection is intentionally not part of compatibility metadata because it must
not affect the frozen normalization.

Reference used: `LCModel_Universal_Philips_UnEdited_PRESS_TE30.BASIS`, SHA-256
`bac0042c2e9a688ef63abd88867a79fd24c659ba92da5c6a333e71f96ffa5abe`, formatted
at 2000 Hz and 2048 samples. FSL-MRS formatting performs resampling and
selection; this utility calculates the scalar on the entire formatted
reference basis before selecting columns.

## Selection experiment

Pixi experiment compares the legacy rule (format selected columns, then divide
by their own mean absolute sample magnitude) with a frozen scalar computed
once from the complete formatted reference. Legacy NAA initial samples are
reported to demonstrate the effect; the frozen NAA vector is exactly the same
for each selection.

| Selection | Legacy scalar | Legacy NAA initial FID |
|---|---:|---:|
| NAA | 9.202857787877271e-05 | 3.1472523556820953 − 0.38183994707259217i |
| NAA + Cr | 9.732599189914159e-05 | 2.9759486943548445 − 0.36105655458217256i |
| NAA + Cr + Cho proxy (GPC + PCh) | 1.7290755315395803e-04 | 1.6750983588388746 − 0.20323107154901118i |
| Complete reference basis (29 entries) | 6.573456836083183e-05 | 4.406162020098117 − 0.5345769841144177i |

“Cho proxy” uses the available GPC and PCh basis entries; there is no single
Cho entry in this reference. The frozen scalar is `6.573456836083183e-05` for
all four selections. NAA normalized vector max absolute and relative difference
across selections: `0` and `0`, respectively, because all selections index the
same full-reference formatted NAA column before division by the same scalar.

Concentration linearity is checked by multiplying the frozen NAA vector by
0.010 and 0.020 mol/L: the 20 mM vector equals twice the 10 mM vector exactly
in the test arithmetic (maximum absolute residual `0`). Repeated formatting
returns array-identical results and calibration metadata.

## Compatibility and legacy behavior

Calibration reuse rejects a different source identity, bandwidth, or point
count with `ValueError`. Selection differences are compatible by design. The
legacy `Basis.reformat()` implementation remains unchanged: its scalar is
`mean(abs(formatted_selected_array))`, hence selection-dependent. This utility
does not call or modify that method.

## Limitations

The frozen scalar stabilizes relative synthetic amplitudes against selection
for a fixed source file and sampling grid. The MRSCloud headers do not establish
absolute signal per mole, receiver calibration, or a validated proton-count
response. This convention therefore does **not** establish physical
concentration calibration or physical water-to-metabolite calibration. It is
not yet integrated with the analytical water-reference path.

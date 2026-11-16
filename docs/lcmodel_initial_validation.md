# Phase 7b: first noise-free LCModel validation

LCModel **6.3-1N** executed successfully on the exact Phase 7a noise-free WM
acquisition. Conversion preserved both complex FIDs exactly. The first successful
fit reports NAA **10.639 mM**, compared with nominal 10 mM and standard FSL-MRS
7.059366427 mM. This establishes an executable conversion/fitting baseline,
**not unbiased concentration recovery or physical scanner calibration**.

## Inputs, installation and execution

- Executable: `$LCMODEL_BINARY`, native arm64,
  executable; version appears in actual TABLE/PS output.
- Manual supplied by the user: `$LCMODEL_MANUAL`.
- Matching local source: `$LCMODEL_SOURCE`.
- Input directory, unchanged:
  `outputs/mrscloud_noise_fitting_validation/full/SNR_infinite/realization_00/`.
- Original BASIS, unchanged:
  `outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative/metabolite_basis.BASIS`.
- 2048 acquired complex points, dwell time 0.0005 s, 127.731 MHz, PRESS TE30,
  receiver factor 1, acquisition reference 4.68 ppm, water 4.70 ppm.
- Existing WM water concentration is **39.0165 M**, not the 40 M used in earlier
  single-metabolite examples. No signals or concentrations were regenerated.
- Phase 6a normalization: `mrscloud-press-te30-full-reference-max-v1`, scalar
  **8857.353925964408**, identity suffix
  `d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81`.

A non-destructive missing-RAW probe ran and returned exit 0 with fatal ZEROVX3 in
its PostScript report. Therefore process status alone is insufficient. Probe
artifacts are retained under `outputs/lcmodel_validation/probe/`. No additional
runtime files or environment setup were needed; the key in the supplied local
example CONTROL was used without modifying the installation.

The actual successful executable invocation was, with CONTROL on stdin:

```sh
cd outputs/lcmodel_validation/noise_free
$LCMODEL_BINARY < lcmodel.CONTROL
```

Python conversion/execution used Pixi. Reproducible commands from repository root:

```sh
pixi run python --version
pixi run python -m examples.lcmodel_initial_validation --reuse-fit
# One fit in a new directory, preserving this first-fit baseline:
pixi run python -m examples.lcmodel_initial_validation --output outputs/lcmodel_validation_reproduction
```

The second command was executed; the fresh-directory example is provided for
reproduction, not claimed as another executed fit. The runner invokes the exact
local binary with `subprocess.run`, CONTROL text on stdin, working directory
`noise_free`, captured stdout/stderr, and a 300-second timeout. A cached run checks
conversion, CONTROL, outputs and source hashes without refitting. Original
`execution.json` records exit **0**, stdout **empty**, stderr **empty**.

## RAW and H2O conversion

The writer stores

\[
 d_{\rm RAW}[n]=\overline{s_{\rm FSL}[n]},
\]

as real, imaginary pairs in unchanged sample order, `(2E25.16)` fixed-width text.
This follows the LCModel spectral convention; FSL's LCModel reader also
conjugates on input. No first-point weighting, gain, or amplitude normalization is
applied. An independent ASCII parser reads the stored pairs without conjugation
and compares their conjugates against the FSL-interpreted source FIDs.

| Metric | Metabolites | Water |
|---|---:|---:|
| Maximum absolute complex error | 0 | 0 |
| Relative complex-array error | 0 | 0 |
| Initial magnitude error | 0 | 0 |
| Maximum phase comparison error, rad | 5.07e-17 | 4.69e-17 |
| Sample order | Exact | Exact |

Initial complex metabolite FSL sample is
`9.24027864591584e-6 - 8.126164138134769e-7 i`;
water is `0.007887908291571334 + 4.72747582344067e-19 i`.
The recovered initial complex metabolite/water ratio has **zero error**.
Water frequency is +2.5546200000000585 Hz in the positive-frequency FSL FID and
-2.5546200000000585 Hz in stored RAW. Thus
`4.68 - f_RAW/127.731 = 4.70 ppm`. The conjugation changes the frequency sign;
it preserves the intended chemical shift.

RAW `SEQPAR` records HZPPPM, ECHOT and SEQ. Dwell time and acquired point count
are CONTROL `DELTAT` and `NUNFIL`; they are not invented RAW namelist fields.
TRAMP and VOLUME are both 1 for data and reference. No gain difference exists.

## BASIS incompatibility and input-only adapter

The original Phase 6a BASIS failed in actual LCModel with **FATAL MYBASI31**,
although FSL-MRS could read it. Its BASIS1 namelist lacks required `FMTBAS`, and
its numeric lines are not in the corresponding declared fixed-width format.
The failed CONTROL, execution log and original-basis PRINT/TABLE/COORD/PS files
are preserved. Manual §8.6.3 and MYBASI31 explain this requirement.

The derived `noise_free/lcmodel_compatible.BASIS`:

1. Declares `FMTBAS='(2E25.16)'` and serializes all 29 spectra at full precision.
2. Sets the existing per-spectrum ISHIFT from 0 to **8**, compensating for the
   installed LCModel reader's implicit 4.65 ppm origin. Local source line 3709
   computes `jshift = ishift + nint((4.65-ppmcen)/ppminc_basis)`.
   With 8192 native points at 4000 Hz and 127.731 MHz, the native step is
   0.0038227309736868892 ppm; the automatic contribution is -8 at PPMCEN4.68.
   Added ISHIFT8 therefore produces **net zero LCModel shift**.

All original complex spectral arrays, names, CONC, TRAMP, VOLUME and sampling
are retained; maximum spectral-array error is **0**. No per-metabolite or water
amplitude correction was made. This is a reader-format/reference adapter, not a
new calibration. Manual §8.6.7 describes ISHIFT. Do not feed this adapted copy
back to FSL-MRS: its reader applies ISHIFT without LCModel's compensating 4.65
term. Continue using the original basis with FSL-MRS.

Original BASIS SHA256:
`21461662797fafdd6aa8c13f6eb1d3084da1db5707a46a253d15fdc4b38e298b`.
Source NIfTI-MRS hashes match Phase 7a; source-sidecar, bundle and executable
hashes are saved and checked before/after conversion validation.

## CONTROL and quantification assumptions

The complete actual CONTROL is retained alongside RAW/H2O. Key settings:

```text
NUNFIL=2048, DELTAT=0.0005, HZPPPM=127.731, ECHOT=30.
PPMCEN=4.68, PPMH2O=4.70, PPMST=4.2, PPMEND=.2
DOWS=T, DOECC=F, WCONC=39016.5, ATTH2O=1., ATTMET=1.
WSMET='NAA', WSPPM=2.01, N1HMET=3
FWHMBA=0.007828953034110749
NSIMUL=0, NRATIO=0, NOMIT=24
```

Explicit CHOMIT entries omit all basis entries except NAA, Cr, GPC, PCh and Glu,
including H2O and Ref0ppm controls. NSIMUL0 disables simulated MM/lipid components;
NRATIO0 removes default concentration-ratio priors. FWHMBA is the known basis
readout linewidth **1 Hz / 127.731 MHz**, not a fitted gain. LCModel retains its
regularized baseline and line-shape model; this is not FSL's baseline-off Voigt
configuration. DOECCFalse avoids changing the already-correct input phase and
frequency frame through eddy-current correction.

Manual §§10.2.1–10.2.2 (pp128–130) define:

\[
 A_{b,\rm norm}=\frac{A_{\rm reference\ basis}}
 {N_{1H,\rm met}\,C_{\rm basis}\,ATTMET},\qquad
 A_{w,\rm norm}=\frac{A_{\rm water}}{2\,WCONC\,ATTH2O},\qquad
 f_{\rm scale}=\frac{A_{b,\rm norm}}{A_{w,\rm norm}}.
\]

LCModel multiplies data by this factor before basis fitting. WCONC is NMR-visible
voxel water in **mM** (default 35880); ATTH2O is water attenuation, also capable of
representing known acquisition-gain differences (default .7). Here WCONC39016.5
and ATTH2O1 describe this exact simulation. ATTMET refers to attenuation of the
reference singlet **in the basis**, not a generic tissue-specific in-vivo
relaxation correction; it is 1 here. No CSF correction is needed for pure WM.
NAA's three-proton methyl reference is used only by the quantification formula;
no extra proton multiplier was applied to either simulated FID.

Relevant installed Fortran source locations:

| Operation | `LCModel.f` line / routine |
|---|---|
| RAW namelists and complex data | 2777 onward, MYDATA |
| BASIS format, CONC/TRAMP/VOLUME, origin shift | 3662–3710, MYBASI |
| BASIS bandwidth conversion | 3770–3880 |
| Water scaling | 4761–4802, WATER_SCALE |
| Water peak/phase/integral | 4900 onward, AREAW2; GETPHA; INTEGRATE |
| Reference basis peak/integral | 5174–5330, AREABA |
| Phased-real COORD arrays | 10775 onward |
| Unitary Fourier transforms | 12847,12910, CFFT/CFFTIN |

These are local source line numbers for this installation, not portable APIs.

## Independent normalization calculation

The validation integrates serialized NAA and water spectra independently using
the source-defined symmetric minima, bordering-mean subtraction and fixed ppm
windows. It does not call LCModel's final reporting function or fit a factor to
returned concentrations. The pure water control has known zero phase; this
limited independent calculation does not implement GETPHA for arbitrary data.

| Factor | Independent NumPy calculation | Actual PRINT |
|---|---:|---:|
| Normalized basis reference area | 0.0014592148729401777 | 0.00145922 |
| Water area | 0.0009338293285514907 | 0.000933848 |
| Water normalized area | 1.1967107871688781e-8 | Derived from printed area |
| FCALIB | 121935.46582732153 | 121933 |

Differences are below 3e-5 relative, accounting for installed single-precision
transforms/integration and five-significant-digit PRINT output. LCModel uses
unitary transforms; the stored BASIS has an unnormalized native 8192-point FFT.
After bandwidth cropping to 4096 frequency points, the ideal LCModel/FSL basis
amplitude conversion is `8192/sqrt(4096) = 128`, derived from transform
normalizations. This is **not an applied signal/calibration factor**.
For an exact matching basis fit, the normalization-only prediction is
`C_mM * FCALIB/(128*1000)`: NAA **9.526208268 mM**.

A direct comparison of LCModel-style bandwidth cropping/inverse transform with
FSL basis formatting, after only this analytically derived factor128, gives
relative full-FID differences Cr0.0005079, GPC0.0006245, Glu0.0013454,
NAA0.0006507, PCh0.0004237. These are residual formatting differences, not tuned
alignment factors. They do not uniquely explain the final concentration bias.
The actual nonlinear LCModel NAA result is another ~11.7% above the
normalization-only prediction. Reference integration and fit/line-shape/baseline
contributions must be distinguished. The remaining aggregate is unresolved;
this first fit does not uniquely attribute it to one processing operation.

## Actual fitting results

| Metabolite | Nominal voxel mM | FSL recovered source coefficient, mM | FSL standard mM | LCModel mM | LCModel error | LCModel %SD |
|---|---:|---:|---:|---:|---:|---:|
| NAA | 10 | 10.000002738 | 7.059366427 | 10.639 | +6.390% | 0 |
| Cr | 8 | 8.000002322 | 5.647493234 | 8.446 | +5.575% | 0 |
| GPC | 2 | 1.999997410 | 1.411871071 | 2.125 | +6.250% | 1 |
| PCh | 1 | 1.000002171 | 0.705937982 | 1.051 | +5.100% | 2 |
| Glu | 8 | 7.999997790 | 5.647490035 | 8.376 | +4.700% | 1 |

LCModel GPC+PCh is 3.176 mM (nominal3); its TABLE rounds %SD to integer percentages.
Reported 0%SD is rounded, not proof of zero uncertainty. FSL raw internal
coefficients, scaling and uncertainty remain in the original fit_result.json;
its standard concentrations include its previously documented water assumptions.
Neither package's standard concentration is identical to its raw fitted amplitude.

LCModel COORD phased-real spectral residual is **0.01816386440649** relative to
phased-real data over its output window. COORD contains real spectra, not the
full complex fit; no complex residual is claimed. Existing FSL complex residual
is 6.0416953e-7, with fitted water initial magnitude0.00788790812234321 versus
expected0.007887908291571334. LCModel reports an integrated water area, not a
Voigt-fitted initial water amplitude; those quantities cannot be equated.
Unbiased per-metabolite initial FID coefficients cannot be uniquely extracted
from these standard LCModel outputs without additional modeling assumptions.

Actual diagnostics, retained without suppressing them:

- MYBASI8 warning: one default preliminary-analysis metabolite (Ins) is absent
  from the explicitly selected five-component spectrum.
- STARTV20 info: narrow preliminary line shape reaches a minimum convolution
  grid extent.
- FINOUT7 info: baseline required substantial freedom and may have a wild shape.
- FINOUT9 info: line shape is not very smooth; displayed FWHM may be unreliable.

TABLE reports FWHM0.004 ppm, S/N393, data shift0ppm, phase0degrees and
0.2degrees/ppm. PRINT retains iterations and metabolite-specific line-shape
parameters. A successful completed analysis is established; the diagnostics
preclude claiming an unbiased ideal fit merely from exit0 or small %SD.

## Inspection outputs and tests

Actual outputs under `outputs/lcmodel_validation/` include comparison.csv,
execution_log.txt, RAW/H2O/CONTROL, basis_compatibility.json,
conversion_validation.json, fitting_results.json and original failure artifacts.
The original LCModel report is `noise_free/lcmodel_outputs/fit.PS`; installed
Ghostscript generated **fit.pdf** and **fit.png**, and the PNG was visually
inspected. The overlay shows aligned peaks but visible residual/line-shape and
baseline differences. No input amplitudes were adjusted after inspection.

```sh
pixi run gs -q -dNOPAUSE -dBATCH -sDEVICE=txtwrite -sOutputFile=outputs/lcmodel_validation/manual.txt $LCMODEL_MANUAL
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_lcmodel_quantitative_io.py tests/integration/test_lcmodel_initial_validation.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q
```

Focused tests: **13 passed** (2.77 s on the final focused run).
Full regression suite: **288 passed, 1 xfailed**, no failures or skips,
in **196.12 s**. The existing expected failure is retained. Warnings were the
installed Dash pytest-plugin deprecation and three SWIG type deprecations.
Fast tests cover complex precision, conjugation, phase/frequency, ratio,
deterministic text conversion, invalid inputs, BASIS preservation/reference
compatibility, CONTROL and output parsing. Separate integration tests verify
actual first-fit artifacts and independent area calculations; pytest does not
run another proprietary fit. If binary/artifacts are absent these tests skip,
which does not count as successful LCModel execution on that machine.

No noise Monte Carlo LCModel experiment was performed. Production simulation,
calibration, legacy behavior, original source signals, MRSCloud and installed
LCModel were unchanged. The next scientific task is to resolve the remaining
LCModel basis-formatting/model/integration aggregate before interpreting its
noise precision or claiming concentration equivalence.

Files created for this phase: `utils/lcmodel_quantitative_io.py`,
`examples/lcmodel_initial_validation.py`, `tests/test_lcmodel_quantitative_io.py`,
`tests/integration/test_lcmodel_initial_validation.py`, this document, and the
dedicated `outputs/lcmodel_validation/` artifacts. Pre-existing Pixi manifest and
lockfile changes were not modified by this phase.

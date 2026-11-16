# Phase 6a: complete PRESS TE30 shared-response library

This standalone library preserves the bundled sequence response and the
validated asymmetric processing conventions: metabolites follow MRSCloud basis
DC/assembly/export processing, while unsuppressed water uses raw pre-DC exact
acquisition-time samples. All responses share one complete-reference scalar.
It remains a synthetic convention, not physical scanner calibration.
`SignalModel`, GUI, analytical mode and Phase 5g calibration are unchanged.

## Supported library and explicit controls

The checked-in PRESS TE30 BASIS contains 29 entries. Every name has an exact
`sys<name>` definition in the bundled spinSystems.mat. No aliases, mergers or
per-metabolite correction factors are introduced. The 27 metabolites are:

```
Ala Asc Asp Cit Cr EA EtOH GABA GPC GSH Glc Gln Glu Gly Ins Lac
NAA NAAG PCh PCr PE Phenyl Scyllo Ser Tau Thr Tyros
```

Two remaining entries are explicitly labeled controls, not metabolites:
`H2O` is the original **RF=3 ppm off-carrier water basis control**; `Ref0ppm`
is the synthetic reference spin. Both remain in the fixed reference assembly
and BASIS file to reproduce the checked-in library convention. The separate
unsuppressed water response is sysH2O acquired with RF=4.7 ppm. Synthesis refuses
H2O/Ref0ppm as metabolite concentrations to prevent mixing those roles.

There are no MM/lipid entries in this checked-in target and no unsupported
entries to exclude. The manifest records all available source spin-system
identifiers, including non-target definitions, and embeds each selected spin
definition with its original scale factors. Exact `sysIns`, `sysGSH`, etc. are
used; other similarly named available definitions are not substituted.

## Sequence and processing

Production uses Universal_Philips unedited PRESS, TE/TE1/TE2=30/14/16 ms,
7 ms universal refocusing waveforms, ideal x excitation, refocusing angle 180
degrees, 101 spatial samples per axis, refocusing FOV 4.5 cm and thickness
3 cm. Native readout is 8192 points at 4000 Hz with 1 Hz Lorentzian decay;
the validated target grid is 2048 points at 2000 Hz. Spectrometer metadata
are 127.731 MHz, analysis origin 4.68 ppm. Metabolite RF is 3 ppm and the
unsuppressed water RF is 4.7 ppm. Source hashes bind the actual pulse files,
spin definitions and bundled implementation; the private wrapper preserves
the spatial amplitude rule (4.5/3)^2/(nX*nY).

`examples/full_basis_responses.m` instruments a private copy of `sim_signals`
to capture each spatially integrated response immediately before first DC.
It records both pre-DC and post-first-DC arrays and the actual sequence and
carrier-relative spin metadata. Bundled source files are not modified.
Repeated NAA, Cr, Glu and water acquisitions are independent sequence runs.

Metabolites receive:

1. Original first `op_dccorr_FID_A` correction, independently checked against
   the captured raw response and actual ppm edge windows.
2. The established second assembly DC correction.
3. One common division by S, followed by the original acquisition-frame shift.
4. LCModel conjugated-FFT spectral serialization at six-digit precision and
   FSL-MRS formatting to the target grid, without subset normalization.

Water receives only:

\[
W(t_n)=\frac{x_{H2O,raw}(t_n)}{S}
e^{i2\pi(4.7-4.68)127.731t_n},\qquad t_n=n/2000.
\]

Native indices 0,2,4,... implement those exact acquisition times. The exported
water is a **unit nominal mol/L response**, not a 40 M acquisition. sysH2O
already encodes its two-proton response. There is no additional factor of two,
water DC correction, Fourier resampling or independent normalization.

DC correction is metabolite-specific **processing**, not a per-metabolite
gain calibration. It changes sample zero by a complex, signal-dependent
offset and can change reference integrals. Fourier formatting also changes
the causal boundary. These asymmetric transformations are retained rather
than silently removing them to make water referencing return a chosen value;
their physical limitations are documented in Phases 5f–5j.1.

## Complete-reference normalization and identity

For the exact checked-in order

```
Ala Asc Asp Cit Cr EA EtOH GABA GPC GSH Glc Gln Glu Gly H2O Ins Lac
NAA NAAG PCh PCr PE Phenyl Ref0ppm Scyllo Ser Tau Thr Tyros
```

let x_m be the native response after both DC operations and before the
acquisition-frame transformation. The frozen normalization is

\[
S_{full}=\max_{m\in\mathcal R,k}
\operatorname{Re}\{\operatorname{fftshift}(\operatorname{FFT}(x_m))[k]\}.
\]

The reference set and ordering are fixed to all 29 entries above. S is
calculated once during bundle assembly and applied unchanged to every basis
column and the raw on-carrier water. It is never calculated from selected
metabolites, fitting coefficients, returned concentrations or a water target.

Convention: **`mrscloud-press-te30-full-reference-max-v1`**, explicitly separate
from Phase 5g's limited-reference identity and scalar. A content-derived
calibration identity hashes the normalization record, including S, complete
reference list, sequence-config hash, raw-source/capture signature, processing
source hashes, package versions and reference-array identity.

The bundle retains the complete unnormalized native processed reference
matrix. Loading independently recomputes S, validates native BASIS amplitudes
against that matrix and the common frame transformation (relative norm
tolerance 3e-7 for LCModel serialization), and validates the water against its
stored raw native response and exact sample indices. File identity, software
versions, reference order, sampling, carrier, sequence and water-proton policy
must agree. Optional expected-calibration metadata provides strict reuse
identity checking. A smoke or Phase 5g calibration cannot be silently reused.

The executed production scalar is **8857.353925964408**, determined by **PCh**.
Its complete identity is:

```
mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81
```

For comparison, the executed smoke scalar is **405.1177739168092** and its
independent identity is
`mrscloud-press-te30-full-reference-max-v1:9072f3ea90368eea4ca16148475c378d057c1f010bec28508d4ed9459fa2db4a`.
That coarse family is not used for production cross-calibration.

## Bundle and API

```
mrscloud_press_te30_quantitative/
    metabolite_basis.BASIS
    water_reference.npz
    reference_manifest.json
    normalization.json
    sequence_config.json
    validation_results.json
    normalization_reference.npz
```

The extra reference-array file enables independent normalization and
amplitude-compatibility checks. The manifest identifies file hashes, raw-cache
hashes, source spin systems, control roles, source BASIS, unsupported entries
and processing source hashes. No external phantom dataset is required.

```python
from simulation.mrscloud_full_basis import load_library, synthesize
from simulation.mrscloud_quantitative_pair import export_pair

bundle = 'outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative'
responses, unit_water, calibration = load_library(
    bundle, selected_names=['NAA', 'Cr', 'GPC', 'PCh', 'Glu'])
pair = synthesize(bundle, {'NAA': 10, 'Cr': 8, 'GPC': 2, 'PCh': 1, 'Glu': 12},
                  water_concentration_mol_l=40, receiver_factor=1)
export_pair(pair, 'new_mixture_directory')
```

Metabolite inputs are mM, water is mol/L, and both use the same receiver gain:

\[
s_m=G\sum_m10^{-3}C_{m,mM}B_m,\qquad s_w=G C_{w,mol/L}W.
\]

These remain nominal concentrations per total tissue volume. There is no
phantom voxel/tissue integration, relaxation correction or residual-water
contamination in this API. The example mixture is intentionally **not fitted**.

## Reproducible workflow and performance

```bash
pixi run python -m examples.mrscloud_full_basis_water_reference outputs/new_full_library
pixi run python -m examples.mrscloud_full_basis_water_reference outputs/new_smoke_library --smoke
pixi run python -m examples.mrscloud_full_basis_water_reference outputs/existing_raw_cache --reuse-raw
```

Use a fresh bundle destination. `--reuse-raw` assembles from already generated
per-response caches without invoking Octave; it refuses to overwrite an existing
bundle. If generation is interrupted before assembly, rerunning without that
flag calls Octave but skips caches with matching source/configuration signatures.
Completed per-response MAT files are never independently renormalized.

Production uses 101 spatial points and 8192 native samples. Smoke uses 3 spatial
points and 256 native samples, exporting 128 target samples. Both use the
complete ordered reference but have distinct identities and scalars. The
coarse spatial grid changes off-resonance sensitivity, so smoke results cannot
stand in for physical or Phase 5g cross-calibration evidence. Coupled systems,
especially seven-spin glucose, dominate full-resolution runtime.

Cache signatures bind all bundled MRSCloud MATLAB/spin/pulse source hashes,
the private capture wrapper, spatial resolution and native point count.
Cached metabolite name, carrier, signature and sequence/grid must match.
Missing cache files and changed configurations fail explicitly. Existing
Phase 5e–5j.1 outputs are read-only comparison inputs.

## Numerical validation

The complete authoritative production numbers are in
[validation_results.json](../outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative/validation_results.json)
and [checked_basis_comparison.csv](../outputs/mrscloud_full_basis/checked_basis_comparison.csv).
The validation includes independently repeated sequence arrays; NAA, Cr and
Glu subset checks; concentration/water/gain linearity; exact mixture summation;
complex NIfTI-MRS round trips; and per-column BASIS serialization checks.

Phase 5g comparison uses only the derived common multiplier
S_full/S_Phase5g, applied to **both** new NAA and water. No separate amplitude
alignment is fitted for either signal. The checked-in basis comparison uses
one explicitly declared complex least-squares global scalar anchored to NAA,
applied unchanged to every column. This scalar is a **comparison only** and
never changes a generated signal or normalization. Per-column correlations,
complex errors, sample-zero ratios, full-band abs-real integral ratios, phase
differences and dominant-peak positions are recorded, including discrepancies.

The on-carrier acetyl three-spin control checks the 3:2 relationship against
raw water. A second methyl control at the actual metabolite carrier retains
the known off-resonance localization difference, rather than treating it as
a missing proton multiplier.

Executed production results:

| Check | Result |
|---|---:|
| Independently repeated NAA, Cr and Glu raw/post-DC maximum complex difference | 0 |
| Independently repeated raw water maximum complex difference | 0 |
| NAA/Cr/Glu subset maximum complex differences | 0 |
| Mixed NAA/Cr/GPC/PCh/Glu superposition maximum difference | 0 |
| Paired mixture NIfTI-MRS export/reload maximum difference | 0 |
| Largest BASIS serialization relative complex error, all 29 columns | 1.9141829392079292e-7 |
| On-resonance methyl per-proton response / water per-proton response | 1.0000000000000038 |
| Metabolite-RF3 methyl per-proton response / water per-proton response | .999997318481034 |

The mixture is NAA 10, Cr 8, GPC 2, PCh 1 and Glu 12 mM with water 40 mol/L,
G=1. Tests independently double individual metabolite concentrations, water
concentration and receiver factor; corresponding arrays scale exactly, without
changing other responses. Metadata and generated arrays are deterministic.
Reassembling a new bundle from copied raw caches without Octave reproduces
the complete calibration metadata and arrays **exactly**.

The new unit NAA initial sample is
**.0002896371585190393 - .000035140187306396487i**; unit water starts at
**.00020216852592034996 + 1.21166066239685e-20i**. These share the new scale and
are not directly comparable in amplitude to Phase 5g without the declared
normalization change.

### Phase 5g cross-calibration

The one derived multiplier from new responses to the old scale is
**S_full/S_Phase5g = 3.1870175638281006**. Applied unchanged to both responses:

| Response | Maximum complex difference | Relative complex norm error |
|---|---:|---:|
| NAA | 1.3729124268490591e-10 | 1.5805897446725846e-7 |
| Water | 2.1852793414608445e-19 | 1.7853983543054977e-16 |

NAA differences are within the justified six-digit LCModel serialization
precision; water differs only by floating-point division/multiplication order.
Neither value is fitted using a concentration result. There is no separate
NAA/water correction.

### Checked-in library comparison

The shared NAA-anchored comparison scalar is
**1.000000000000001 + 3.066771406749467e-17i**. The newly generated complete
library therefore reproduces the checked-in source scale directly to numerical
precision, rather than needing distinct metabolite gains. Across all 29 entries:

- Maximum relative complex-array error: **3.709428105613202e-12**, for GPC.
- Minimum complex-correlation magnitude: **.9999999999999942**.
- Maximum initial complex-amplitude ratio deviation from 1: **8.5843e-13**.
- Maximum full-band abs-real integral-ratio deviation from 1: **3.1586e-13**.
- Maximum correlation phase difference: **3.4256e-14 rad**.
- Dominant chemical-shift bins agree; the CSV also handles near-degenerate
  peak maxima within the export rounding tolerance without inventing shifts.

| Entry | Relative complex error | Dominant peak ppm |
|---|---:|---:|
| Ala | 2.74e-13 | 1.438324 |
| Asc | 2.01e-13 | 3.747254 |
| Asp | 1.12e-15 | 2.661598 |
| Cit | 2.78e-13 | 2.615725 |
| Cr | 8.56e-13 | 3.028580 |
| EA | 6.50e-14 | 3.777835 |
| EtOH | 1.14e-15 | 1.170733 |
| GABA | 2.31e-12 | 2.286970 |
| GPC | 3.71e-12 | 3.212071 |
| GSH | 1.55e-15 | 3.770190 |
| Glc | 1.14e-15 | 3.418499 |
| Gln | 1.14e-15 | 6.813084 |
| Glu | 1.15e-15 | 2.348134 |
| Gly | 1.61e-13 | 3.548472 |
| H2O control | 2.18e-13 | 4.702936 |
| Ins | 1.10e-15 | 3.563763 |
| Lac | 1.53e-13 | 1.285415 |
| NAA | 9.97e-15 | 2.011734 |
| NAAG | 1.15e-15 | 2.042316 |
| PCh | 3.59e-13 | 3.204426 |
| PCr | 1.64e-14 | 3.028580 |
| PE | 3.64e-12 | 3.257944 |
| Phenyl | 1.17e-15 | 7.386494 |
| Ref0ppm control | 9.76e-14 | .000977 |
| Scyllo | 5.16e-13 | 3.342044 |
| Ser | 7.31e-14 | 3.938390 |
| Tau | 2.26e-13 | 3.242653 |
| Thr | 1.12e-15 | 1.293060 |
| Tyros | 1.15e-15 | 6.920120 |

The linked CSV contains every entry's full complex correlation, initial
amplitude ratio, integral ratio and phase difference, all under that one global
alignment. Peak positions describe the strongest observed response, which need
not be the resonance normally chosen for fitting (for example Gln's strongest
full-band peak is outside the usual .2–4.2 ppm fitting window).

### Processing effects

Every entry's raw, first-DC, second-DC and formatted initial samples and
complex-change norms are recorded under `processing_effects`. First-DC
sample-zero magnitude ratios span **.49824236973095026–.5229702199576024**,
demonstrating signal-dependent, non-multiplicative processing changes. These
ratios are diagnostics only and are not applied as calibration factors. Water
bypasses them entirely. Agreement with the checked-in processed library and
the Phase 5g pair establishes representation consistency, not equality of
raw and processed physical acquisition signals.

No optional mixture fit or new HTML report is generated: this phase validates
basis and signal generation and follows the instruction to leave the mixture
unfitted. Existing concentration-bias results are not used as acceptance targets.

## Executed commands and tests

The full-resolution production run used the generated, hash-bound command file:

```bash
pixi run python --version
pixi run octave --version
pixi run octave --no-gui --quiet outputs/mrscloud_full_basis/run_capture.m
pixi run python -m examples.mrscloud_full_basis_water_reference outputs/mrscloud_full_basis --reuse-raw
pixi run python -m examples.mrscloud_full_basis_water_reference outputs/mrscloud_full_basis_smoke --smoke
pixi run python -m examples.mrscloud_full_basis_water_reference outputs/mrscloud_full_basis_smoke_final --smoke --reuse-raw
pixi run python -m pytest -q tests/test_mrscloud_full_basis.py tests/integration/test_full_basis_sequence_validation.py
pixi run python -m pytest -q tests
```

`run_capture.m` records the exact `full_basis_responses` call, entry file, 101/
8192 configuration and source signature. For new runs, the Python entry point
generates these inputs and calls the same Octave function automatically.
Full generation completed all 36 reference/control/repeat responses. The final
smoke bundle was rebuilt from copied smoke raw caches to exercise reuse; its
source processing hashes identify that assembly version. Runtime was GNU
Octave **10.3.0**, FSL-MRS **2.5.1.dev78+ga311e38f2**. Package/lockfile and source
hashes are retained with the bundle.

An initial launch preceded creation of the generated command file and failed
with “no such file”; rerunning after creation completed successfully. No partial
response was substituted. Bundled Octave load-path/deprecation/shutdown
messages remain in `outputs/mrscloud_full_basis/octave.log`.

Focused suite: **20 passed** in **12.90 s**. It covers manifest integrity,
shared-scalar validation, three-metabolite subset invariance, concentration/
water/gain linearity, complex paired export, metadata tampering, missing-cache
errors, smoke identity separation, actual independent repeats, proton controls,
all-column comparisons and exact cached rebuilding without Octave.
Full regression suite: **155 passed, 1 xfailed, 5 dependency warnings in
95.16 s**. The expected failure is the existing installed FSL-MRS combined
phase/frequency water-initialization defect. All prior analytical, voxel,
frozen-basis, calibration, processing, export, CLI, quantification and linewidth
regressions executed. No required Phase 6a verification remains pending.

Files created: `examples/full_basis_responses.m`,
`examples/mrscloud_full_basis_water_reference.py`,
`simulation/mrscloud_full_basis.py`, `tests/test_mrscloud_full_basis.py`,
`tests/integration/test_full_basis_sequence_validation.py`, and this document.

## Limitations and review boundary

This does not establish physical molar/scanner calibration for all 27
metabolites. Original spin-system amplitudes and couplings are preserved;
they have not been independently validated as a complete physical proton
response library. Ideal excitation, two-direction localization, finite RF
sensitivity, the Hamiltonian/export gyromagnetic mismatch, DC backgrounds,
Fourier boundary artifacts and fitting integration conventions remain relevant.
Subset invariance and numerical agreement do not remove those limitations.
Review the complete common-normalized library before any phantom integration.
LCModel, GUI and `SignalModel` integration are outside Phase 6a.

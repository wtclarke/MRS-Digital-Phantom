# Phase 5g: sequence-matched quantitative pairs

The standalone generator now accepts `water_reference_mode='analytical'`
(the default) or `'mrscloud_sequence'`. Both remain synthetic conventions,
not physical scanner calibrations. Legacy `SignalModel`, residual-water
simulation, GUI, MRSCloud source, and existing exporters are unchanged.

## API and reproducibility

```python
from examples.mrscloud_sequence_water_integration import build_reference_bundle
from simulation.mrscloud_quantitative_pair import generate_pair, export_pair

# Package the existing validated Phase 5e/5f simulations; no new Octave run.
reference = build_reference_bundle('new_output/reference')
pair = generate_pair(
    {'NAA': 10}, water_concentration_mol_l=40, receiver_factor=1,
    water_reference_mode='mrscloud_sequence', sequence_reference=reference)
export_pair(pair, 'new_output/pair', aligned_water=True)
```

`sequence_reference` is an explicit directory containing `reference.json`,
`common_response.BASIS`, and `raw_water.npz`. The packaging helper reads the
preserved Phase 5e native simulation and Phase 5f pre-DC capture, verifies the
post-first-DC response against Phase 5e, and recomputes the shared scalar.
It refuses a nonempty destination. Generation itself needs no Octave runtime.
The bundle, not the original Phase 4a basis, supplies sequence-mode metabolites.

Currently this bundle contains NAA and the explicitly named Ref0ppm reference
control. Other metabolite names are rejected; this phase does not generate a
new full metabolite basis. Selecting NAA with/without the reference control
does not change its response or normalization.

## Distinct normalization identities

Analytical mode retains `synthetic-frozen-Ref0ppm-initial-magnitude-v1`, the
Phase 4a frozen scalar `6.573456836083183e-05`, and the original basis SHA256
`bac0042c2e9a688ef63abd88867a79fd24c659ba92da5c6a333e71f96ffa5abe`.
Its signal arrays match the saved pre-feature Phase 5b dataset exactly.
Configuration adds only the explicit mode field; prior metadata values remain
unchanged.

Sequence mode uses `mrscloud-press-te30-raw-water-common-max-v1`, not the Phase
4a scalar. Its fixed scalar is

\[
S=\max_{m\in\{\mathrm{NAA,Ref0ppm,water\_off}\},k}
\operatorname{Re}\{\mathrm{FFT}_{\mathrm{shift}}(x_m)[k]\}
=2779.198340948381.
\]

Here the fixed reference arrays are the Phase 5e native, basis-processed
responses before shared division and acquisition-frame shifting. They are
stored in the bundle so the scalar can be checked independently on loading.
The common basis SHA256 is
`1ab8a9706a8cfd5d8183fda7a7d2ba917b04c92c33d45a5ed69d07df3d82463e`.
The validated manifest SHA256 is
`f05598c72643b424ea0159b391e7167219504f38d22147fe0ab470e06932b598`;
the complete calibration identity appends this hash to the convention name.

The loader verifies file hashes, formatting software version, grid, sequence
configuration, independently recomputed scalar, and native NAA/Ref0ppm basis
responses against the stored normalization anchors. The latter check allows
relative norm error `3e-7` for LCModel text serialization, not selection-dependent
rescaling. Exact expected metadata can also be supplied. Original Phase 4a
calibration objects/bases, incompatible sampling, and extra water-linewidth,
phase, or frequency overrides are rejected. Selection never recomputes S.

## Sequence and signal equations

PRESS TE30 uses TE1/TE2 = 14/16 ms, 7 ms universal refocusing waveforms,
101 spatial samples per refocusing axis, RF carrier 3.0 ppm for NAA and
4.7 ppm for water. Excitation is the existing ideal `sim_excite` operation;
this is not a complete shaped-excitation physical acquisition model.

For nominal metabolite concentration in mM and water concentration in mol/L,

\[
s_m(t_n)=G\sum_m(10^{-3} C_{m,\mathrm{mM}}) B^{\mathrm{common}}_m(t_n),
\qquad
s_w(t_n)=G C_{w,\mathrm{mol/L}}\,
\frac{x^{\mathrm{raw}}_w(t_n)}{S}
e^{i2\pi(4.7-4.68)\nu_0t_n}.
\]

Signals are in arbitrary shared receiver units. `sysH2O` already encodes two
protons: there is no extra multiplier of two. The water response bypasses both
DC corrections and Fourier-crop resampling. Native acquisition is 8192 points
at 4000 Hz; target acquisition is 2048 points at 2000 Hz, using native indices
0, 2, 4, ... exactly. The existing 1 Hz simulated readout decay is retained
once, with no additional relaxation. Metabolites use the unchanged matching
Phase 5e formatted basis, including its basis-processing convention.

The common spectral origin is 4.68 ppm at 127.731 MHz. Positive complex
frequency follows `exp(+2j*pi*f*t)`; water offset is +2.55462 Hz and its
physical peak is 4.70 ppm. No RF-carrier change is implemented as a substitute
for the original sequence simulation.

## Export and metadata

Primary files are `metabolites.nii.gz`, `water_reference.nii.gz`,
`ground_truth.json`, and `simulation_config.json`, using the validated
`gen_nifti_mrs` convention. Headers retain 0.0005 s dwell time, 127.731 MHz,
`SpecFreqChemShift=4.68`, nucleus 1H, and TE 0.030 s. Custom calibration fields
stay in JSON. Sidecars include nominal concentration units, shared receiver
factor, basis/calibration identities, S, sequence/carriers, grid, source
hashes, water processing, and `physical_calibration=false`.

`aligned_water=True` additionally exports `water_reference_aligned.nii.gz`
and `water_preprocessing.json`. It demodulates the declared water frequency
and phase without modifying either primary FID. Its analysis peak is 4.68 ppm;
the correction sidecar retains source position 4.70 ppm and explicitly states
that the metabolite frame and physical water position were not changed.
Both CLI comparisons below use the original primary water, not this aligned
alternative.

## Executed numerical validation

FSL-MRS `2.5.1.dev78+ga311e38f2`; NAA 10 mM, water 40 M, G=1:

| Quantity | Analytical | Sequence-matched |
|---|---:|---:|
| Initial NAA complex sample | 0.04406162020098117 - 0.005345769841144177i | 9.230786561072083e-6 - 1.119923229041747e-6i |
| Initial water complex sample | 61.73993435001698 | 0.02577258571845567 + 1.5446335249833405e-18i |
| Fitted water amplitude | 61.73993439719796 | 0.02577258514572865 |
| Fitted water frequency, Hz | 2.5546199998466563 | 2.554620145877357 |
| Fitted water phase, rad | 3.2565282020350634e-10 | -1.4340570100684988e-7 |
| Fitted water linewidth, Hz | 0 | 0.9999999372792765 |
| Relative complex water residual | 8.326188710382894e-10 | 1.4441623108503362e-7 |
| Internal fitted NAA coefficient | 1.0000000905278943 | 1.0000000764364945 |
| Raw NAA coefficient | 152.1269729860657 | 0.01000000076436493 |
| Relative complex metabolite residual | 1.8546135081785328e-5 | 1.857590496917146e-5 |
| Standard FSL-MRS concentration, mM | 6.003318300954079 | 6.885791303512254 |
| Maximum complex export/reload difference | 0 | 0 |

The sequence-mode full water and NAA arrays match Phase 5f exactly (maximum
absolute error zero). Water optimizers report success. The two raw coefficient
scales reflect different basis normalization identities; they are not directly
comparable across modes. Neither FSL-MRS concentration was used to tune signals.
The standard result still includes FSL-MRS proton/integral, water-content and
relaxation assumptions described in the quantification-scaling derivation.

Reproduce into a new directory:

```bash
pixi run python -m examples.mrscloud_sequence_water_integration outputs/mrscloud_sequence_water_integration/validated
pixi run python -m pytest -q tests
```

The example runs the actual CLI twice with `--report`, Newton, Lorentzian,
baseline off, ppm limits 0.2–4.2, NAA only, `--no_conj_basis`, TE 30 ms, TR 2 s,
WM/GM/CSF fractions 1/0/0, reference protons 3 and integration 1.8–2.2 ppm.
Complete commands are in `outputs/mrscloud_sequence_water_integration/validated/cli_commands.txt`.
CLI concentrations are read from its result CSV; supplemental Python replay
records actual complex water-model residuals and NAA nuisance fitting behavior.

The complete `tests` run executed successfully: **107 passed, 1 xfailed,
5 dependency warnings in 70.25 s**. The expected failure retains the installed
FSL-MRS combined phase/frequency water-initialization defect. Coverage includes
the existing analytical, export, frozen-normalization, calibration, Phase 3,
and fitting suites, plus 23 new integration-unit cases and two real-CLI tests.
Exact-array tests establish selection independence, deterministic output,
concentration/receiver linearity, Phase 5f equivalence, and preservation of the
saved analytical signals. The fixed-seed legacy residual-water comparison
also verifies that generation consumes neither NumPy nor Torch RNG state.

Actual nonempty HTML reports are:

- `outputs/mrscloud_sequence_water_integration/validated/analytical/fsl_mrs_results/report.html` (887487 bytes).
- `outputs/mrscloud_sequence_water_integration/validated/mrscloud_sequence/fsl_mrs_results/report.html` (887861 bytes).

`results.json`, `comparison.csv`, per-mode `cli.log`, and
`integration_validation.json` retain numerical evidence. The known HTML water
overlay may compare corrected model and uncorrected data; the residuals above
use the same complex reference frame and do not infer accuracy from that plot.

## Limits

No new physical calibration, tissue mixing, relaxation correction, LCModel
validation, or production integration is claimed. The reference bundle relies
on preserved Phase 5e/5f experiment artifacts; tests require those small local
simulation artifacts, not a full phantom dataset. Rebuilding the underlying
sequence simulation is outside this integration phase. A full sequence-matched
metabolite library remains future work. Fixed-seed residual-water/RNG tests
and analytical saved-output comparisons cover isolation; they do not establish
a historical full-phantom numerical baseline.

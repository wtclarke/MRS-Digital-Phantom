# Phase 6e: optional tissue-specific pre-acquisition relaxation

## Scope and scientific convention

`simulate_voxel(..., relaxation=None)` retains the Phase 6b no-attenuation
behavior exactly, including existing metadata. The opt-in configuration applies
only scalar saturation and echo attenuation, assuming 90-degree excitation:

\[
R_{m,j}=(1-e^{-TR/T_{1,m,j}})e^{-TE/T_{2,m,j}},\qquad
R_{w,j}=(1-e^{-TR/T_{1,w,j}})e^{-TE/T_{2,w,j}}.
\]

All times are seconds. With nominal tissue concentrations C_m,j in mM, water
C_w,j in mol/L, tissue volume fractions f_j and shared receiver G:

\[
s_m(t)=G\sum_j f_j10^{-3}C_{m,j}R_{m,j}B_m(t),\qquad
s_w(t)=G\sum_j f_j C_{w,j}R_{w,j}W(t).
\]

The implementation computes relaxation-weighted effective concentrations and
passes them to the unchanged Phase 6a synthesizer. It preserves nominal ground
truth separately. B_m and W keep their original common scalar
8857.353925964408 and calibration identity
`mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81`.
There are no new gain factors, proton multipliers or changes to unit water.

**Three distinct mechanisms remain separate:** PRESS coherent TE-dependent
sequence evolution is already in the basis; the new TE/T2 attenuation is a
single pre-acquisition scalar; the original 1 Hz acquisition-time decay stays
inside both response shapes. The new T2 values do not replace that readout
linewidth, and no additional exp(-t/T2) is applied during acquisition.

Bundled `sim_signals.m:368` calls `sim_readout` only after spatial integration
and sequence propagators; `functions/FID-A/simulationTools/sim_readout.m:62`
applies its linewidth-derived decay on acquisition sample times. Sequence
propagation is Hamiltonian evolution. No tissue T1/TR saturation or TE/T2
attenuation is present in this frozen simulation. Applying these declared
factors once therefore does not double-count those effects. This is a
phenomenological synthetic extension, not a full relaxation-aware Bloch or
localization simulation.

## API and validation rules

```python
from examples.mrscloud_tissue_relaxation_validation import configuration
from examples.mrscloud_quantitative_voxel import BUNDLE, METABOLITES, WATER
from simulation.mrscloud_quantitative_voxel import simulate_voxel
from simulation.mrscloud_quantitative_pair import export_pair

pair = simulate_voxel(
    {'WM': .2, 'GM': .3, 'CSF': .5}, calibration_bundle=BUNDLE,
    metabolite_concentrations_mM=METABOLITES,
    water_concentrations_mol_l=WATER,
    relaxation=configuration(), receiver_factor=1.)
export_pair(pair, 'new_relaxed_voxel')
```

`relaxation` is either None, or an explicit dictionary:

```python
{
    'mode': 'tissue_specific',
    'TE_s': .030,
    'TR_s': 2.,
    'water': {
        'WM': {'T1_s': .97, 'T2_s': .073},
        'GM': {'T1_s': 1.5, 'T2_s': .088},
        'CSF': {'T1_s': 4.47, 'T2_s': 2.03},
    },
    'metabolites': {
        # All selected names explicitly required in each of WM, GM, CSF:
        'WM': {'NAA': {'T1_s': 1.29, 'T2_s': .194}},
        'GM': {'NAA': {'T1_s': 1.45, 'T2_s': .180}},
        'CSF': {'NAA': {'T1_s': 1.5, 'T2_s': .25}},
    },
}
```

The abbreviated dictionary above selects NAA only; the executable example
uses complete five-metabolite profiles. No default tissue relaxation times
are inserted by the adapter. All profiles must match selected names, including
zero concentrations and zero-volume tissues. Missing, additional, misspelled
or millisecond-valued keys are rejected rather than ignored. All supplied
times must be finite, T1/T2 positive, TE nonnegative and TR greater than TE.
Units are explicit `_s` fields, not inferred from numerical magnitude.

An explicit no-attenuation acquisition is:
`{'mode':'none', 'TE_s':.030, 'TR_s':2.}`. It requires no T1/T2 profiles and
records factors of exactly one. Arbitrary contradictory profile keys are
rejected in this mode. None preserves the historical metadata without adding
TR. There is no use of infinite T1/T2 as a hidden no-attenuation switch.

The adapter checks TE against the loaded calibration sequence. This bundle
is PRESS TE30; TE other than .030 s is rejected within 1e-12 s tolerance.
Changing sequence TE requires a new coherent sequence simulation/basis, not
just changing the attenuation exponent. The independent scalar helper tests
TE dependence at 30 and 60 ms, but does not claim to simulate PRESS TE60 with
the TE30 basis. TR may vary while retaining the same sequence unit responses.

## Metadata and export

Existing nominal fields, including `metabolite_concentrations_mM`,
`voxel_average_metabolite_concentrations_mM`, tissue-only concentrations and
`effective_voxel_water_concentration_mol_l`, retain their original meanings.
Additional fields when relaxation is explicitly configured are:

- `relaxation`: mode, TE/TR, seconds unit, tissue-specific T1/T2 values,
  water and metabolite factors, equation, application point and basis assumption.
- `relaxation_weighted_metabolite_concentrations_mM`: sum(f*C*R) per metabolite;
  these are signal coefficients in mM-equivalent units, not nominal truth.
- `relaxation_weighted_water_concentration_mol_l`: sum(f*Cw*Rw), likewise
  distinguished from nominal effective water concentration.
- Simulation config includes explicit `TE_s`, `TR_s`, the relaxation metadata,
  and `tissue_specific_relaxation` mode flag.

Both JSON sidecars are reproducible and keep the complete tissue profiles and
unchanged calibration identity. NIfTI-MRS uses the existing validated 1D
`gen_nifti_mrs` convention, adding standard EchoTime and RepetitionTime only
from explicit acquisition config. Default EchoTime remains .030 and default
exports still have no newly invented TR. Optional aligned-water exports receive
the same acquisition timing metadata. Primary water frequency and phase are
unchanged. Tests compare full complex interpreted arrays exactly; no magnitude
or phase information is discarded.

## Executed illustrative configuration

Nominal concentrations, tissue fractions and water profiles are unchanged from
Phase 6b/6c: WM/GM/CSF water 39.0165/47.2305/55.5 M; WM NAA/Cr/GPC/PCh/Glu
10/8/2/1/8 mM; GM 14/9/1.5/1.2/12 mM; CSF all zero. Mixed fractions .2/.3/.5.
TE=.030 s, TR=2 s, G=1, no noise/B0/residual water. These relaxation values
are **declared illustrative inputs**, not a newly validated literature database.
Water times coincide with installed FSL 3T defaults for comparison.

| Species | WM T1/T2 s | GM T1/T2 s | CSF T1/T2 s |
|---|---|---|---|
| water | 0.97/0.073 | 1.5/0.088 | 4.47/2.03 |
| NAA | 1.29/0.194 | 1.45/0.18 | 1.5/0.25 |
| Cr | 1.1/0.15 | 1.35/0.17 | 1.5/0.25 |
| GPC | 1.2/0.2 | 1.3/0.22 | 1.5/0.25 |
| PCh | 1.1/0.18 | 1.25/0.21 | 1.5/0.25 |
| Glu | 1.15/0.16 | 1.4/0.19 | 1.5/0.25 |

| Species | R WM | R GM | R CSF |
|---|---:|---:|---:|
| water | 0.578667054841 | 0.523673419482 | 0.355437529169 |
| NAA | 0.674956729225 | 0.633378160042 | 0.653130747862 |
| Cr | 0.685833876845 | 0.647694661951 | 0.653130747862 |
| GPC | 0.698141238511 | 0.685184364340 | 0.653130747862 |
| PCh | 0.709080293952 | 0.691858270264 | 0.653130747862 |
| Glu | 0.683391008445 | 0.649292139698 | 0.653130747862 |

## Analytical and numerical fitting results

Independent expected FIDs are direct tissue sums using NumPy's separate
(1-exp(-TR/T1))*exp(-TE/T2) implementation and loaded unit responses, rather
than the adapter's factor metadata or returned effective concentrations.
Maximum metabolite complex-array error is 2.396e-21; water error 1.940e-18.
All four paired export/reloads have exactly zero complex-array error.

FSL-MRS version `2.5.1.dev78+ga311e38f2`, Newton, Lorentzian, baseline off, 0.2–4.2 ppm,
matching complete Phase 6a BASIS with NAA/Cr/GPC/PCh/Glu selected, one shared
linewidth/frequency group, no conjugation. NAA water-scaling reference uses
3 protons and 1.8–2.2 ppm; corrected water uses 2 protons, 1.65–7.65 ppm.
Sampling remains 2048 points at 2000 Hz, 127.731 MHz, reference 4.68 ppm.

Expected internal coefficient a_m = G*1e-3*sum(f*C_m*R_m)*s_D/s_B.
Reported recovered effective mM = 1000*a_m*s_B/(s_D*G); this is not nominal mM.
The complete internal coefficients, expected coefficients, nuisance parameters,
scaling and optimizer messages are saved in `results.json` and `comparison.csv`.

| Voxel | Metabolite | Nominal voxel mM | Weighted expected mM | Fitted effective mM | Relative coefficient error |
|---|---|---:|---:|---:|---:|
| pure_WM | NAA | 10 | 6.7495672922 | 6.7495677926 | 7.4138e-08 |
| pure_WM | Cr | 8 | 5.4866710148 | 5.4866716881 | 1.2272e-07 |
| pure_WM | GPC | 2 | 1.3962824770 | 1.3962808934 | -1.1342e-06 |
| pure_WM | PCh | 1 | 0.7090802940 | 0.7090820833 | 2.5235e-06 |
| pure_WM | Glu | 8 | 5.4671280676 | 5.4671164400 | -2.1268e-06 |
| pure_GM | NAA | 14 | 8.8672942406 | 8.8672941881 | -5.9206e-09 |
| pure_GM | Cr | 9 | 5.8292519576 | 5.8292522228 | 4.5498e-08 |
| pure_GM | GPC | 1.5 | 1.0277765465 | 1.0277765489 | 2.3614e-09 |
| pure_GM | PCh | 1.2 | 0.8302299243 | 0.8302299322 | 9.5471e-09 |
| pure_GM | Glu | 12 | 7.7915056764 | 7.7915050198 | -8.4265e-08 |
| mixed | NAA | 6.2 | 4.0101017306 | 4.0101016502 | -2.0052e-08 |
| mixed | Cr | 4.3 | 2.8461097902 | 2.8461097625 | -9.7424e-09 |
| mixed | GPC | 0.85 | 0.5875894594 | 0.5875895188 | 1.0115e-07 |
| mixed | PCh | 0.56 | 0.3908850361 | 0.3908849560 | -2.0501e-07 |
| mixed | Glu | 5.2 | 3.4308773164 | 3.4308772453 | -2.0745e-08 |

All three metabolite optimizers converge, with maximum coefficient error
2.524e-6 relative; individual coefficients are tested, not inferred from a
small residual alone. GPC/PCh remain overlapping components. The pure WM
optimizer stops with a 2.078e-6 complex residual; GM and mixed residuals are
7.975e-8 and 4.297e-8. No amplitudes or convergence criteria were adjusted.
The saved parameters show only tiny additional linewidth/frequency/phase terms.

| Voxel | Nominal water M | Weighted water M | Expected initial amplitude | Fitted amplitude | Relative amplitude error | Complex water residual |
|---|---:|---:|---:|---:|---:|---:|
| pure_WM | 39.0165 | 22.5775631452 | 0.00456447265994 | 0.0045644725572 | -2.2507e-08 | 1.4531e-07 |
| pure_GM | 47.2305 | 24.7333574389 | 0.00500030641447 | 0.00500030631008 | -2.0877e-08 | 1.4467e-07 |
| pure_CSF | 55.5 | 19.7267828689 | 0.00398813461376 | 0.00398813452252 | -2.2877e-08 | 1.4594e-07 |
| mixed | 49.7225 | 21.7989112951 | 0.00440705376321 | 0.00440705366783 | -2.1643e-08 | 1.4266e-07 |

All four water optimizers converge; the largest relative amplitude error is
2.288e-8. Linewidth and frequency remain the same unit-response values, about
1 Hz and +2.55462 Hz, phase about zero. Scalar echo attenuation does not
broaden the water peak. Residuals compare the model and data in the same
complex frame; the known HTML integrated-water overlay frame issue remains.
Pure CSF stays exactly zero in metabolites and has no defined metabolite
reference or tissue-only quantification. Only its water fit is performed.

## Quantification equations and model limits

Let C_m^R=sum(f*C_m*R_m) in mM and C_w^R=sum(f*Cw*Rw) in M. Independent,
unfitted unit-response integrals are J_NAA=.006425124587755786 and
J_W=.004767430133932897; alpha=(2/3)*J_NAA/J_W=.898474916569691.
They use the existing first-point-half-weighted, orthonormal FFT and absolute
real-spectrum bin integration with nearest inclusive ppm endpoints. The water
integral uses the known unit amplitude and exp(-pi*t) envelope after removing
its known frequency/phase, rather than a fitted amplitude. NAA is the reference
for every metabolite; targets are not independently integrated for scaling.

For the actual FSL quantification inputs Q:

\[
\widehat C_m=\frac{10^{-3}C_m^R J_{NAA}}{C_w^R J_W}\frac{2}{3}
\frac{55010\sum_j f_jd_jR^{Q}_{w,j}}{(1-f_{CSF})R^{Q}_{M}}.
\]

This prediction uses source concentrations and declared relaxation parameters,
not fitted results. The fitting data scale, basis scale and receiver gain cancel.
The integration/sequence factor alpha is unchanged by pre-acquisition relaxation.

Two separate quantification configurations are recorded:

1. **Actual CLI defaults:** TE30/TR2 and correct tissue fractions, installed
   water densities .65/.78/.97 and 3T relaxation constants. Its single METAB
   T1/T2 is 1.29/.194 s in every voxel. These are the actual HTML concentrations.
2. **API declared-water/NAA-scalar:** public `QuantificationInfo` T1/T2 inputs
   match every tissue's water parameters; `set_densities` uses the declared
   tissue water molarity*1000/55010. The single METAB pair equals the actual
   NAA pair for pure WM/GM. For mixed tissue it is the explicit WM/GM volume
   mean NAA T1=1.386, T2=.1856 s. This is an approximation, not an effective
   physical T1/T2 or an assertion that averaging exponentials commutes.

All fitted parameters are asserted unchanged across these configurations.
CLI help also supports `--t1-values`/`--t2-values` JSON overrides, still with
only one METAB pair. Those switches cannot represent the full simulator's
metabolite/tissue tensor; no `--h2o_scale` gain is introduced. API results use
explicit density overrides and are kept separate from baseline HTML output.

With exact water settings the equation simplifies to
alpha*C_m^R/[(1-f_CSF)*R_M^Q]. For pure-tissue NAA the correct R_M^Q cancels
its attenuation, leaving alpha times nominal concentration. Other metabolites
have different R_m, so one NAA scalar does not exactly correct them. In the
mixed voxel the true concentration-weighted NAA factor is .6467906017140265,
whereas the documented scalar approximation is .6497871646937979; their ratio
is approximately .995388. This small extra bias is predicted, not repaired.
There is no single FSL scalar correction that removes all independently chosen
metabolite/tissue attenuation factors simultaneously.

| Voxel | Metabolite | Nominal voxel mM | Default predicted | Default reported | API predicted | API reported |
|---|---|---:|---:|---:|---:|---:|
| pure_WM | NAA | 10 | 8.234033897 | 8.234035461 | 8.984749166 | 8.984750872 |
| pure_WM | Cr | 8 | 6.693382429 | 6.693384025 | 7.303633061 | 7.303634803 |
| pure_WM | GPC | 2 | 1.703373972 | 1.703372237 | 1.858674383 | 1.858672490 |
| pure_WM | PCh | 1 | 0.865031923 | 0.865034206 | 0.943898817 | 0.943901308 |
| pure_WM | Glu | 8 | 6.669541302 | 6.669527890 | 7.277618285 | 7.277603650 |
| pure_GM | NAA | 14 | 10.723457866 | 10.723457911 | 12.578648832 | 12.578648885 |
| pure_GM | Cr | 9 | 7.049471469 | 7.049471861 | 8.269051566 | 8.269052026 |
| pure_GM | GPC | 1.5 | 1.242917872 | 1.242917888 | 1.457946461 | 1.457946479 |
| pure_GM | PCh | 1.2 | 1.004019419 | 1.004019439 | 1.177717845 | 1.177717869 |
| pure_GM | Glu | 12 | 9.422477767 | 9.422477068 | 11.052595202 | 11.052594382 |
| mixed | NAA | 6.2 | 9.972501520 | 9.972501426 | 11.089710643 | 11.089710538 |
| mixed | Cr | 4.3 | 7.077833959 | 7.077833965 | 7.870756442 | 7.870756449 |
| mixed | GPC | 0.85 | 1.461243921 | 1.461244084 | 1.624945580 | 1.624945761 |
| mixed | PCh | 0.56 | 0.972070506 | 0.972070317 | 1.080970568 | 1.080970358 |
| mixed | Glu | 5.2 | 8.532060170 | 8.532060084 | 9.487898125 | 9.487898029 |

Predictions agree within 1.464e-5 mM absolute and 2.64e-6 relative, consistent
with the measured coefficient/nuisance fitting error. These are successful
predictions of the standard quantification convention, not proof that every
reported value equals nominal ground truth. Mixed values include FSL's CSF
correction of 2: nominal mixed NAA is 6.2 mM per full voxel, 12.4 mM per WM/GM
volume. The API NAA result is 11.089710538 mM on that corrected convention,
with voxel-volume interpretation 5.544855269 mM. The factor alpha and the
explicit scalar relaxation approximation remain. Pure CSF quantification is
undefined, not a reported zero concentration.

## Actual reports and reproducible execution

Three nonempty actual `--report` HTML files:

- [WM](../outputs/mrscloud_tissue_relaxation_validation/final/pure_WM/fsl_mrs_results/report.html)
- [GM](../outputs/mrscloud_tissue_relaxation_validation/final/pure_GM/fsl_mrs_results/report.html)
- [Mixed](../outputs/mrscloud_tissue_relaxation_validation/final/mixed/fsl_mrs_results/report.html)

These contain CLI-default quantification only; API comparisons are in
[comparison.csv](../outputs/mrscloud_tissue_relaxation_validation/final/comparison.csv)
and [results.json](../outputs/mrscloud_tissue_relaxation_validation/final/results.json).
Primary datasets preserve original complex water phase/frequency. The original
Phase 6a bundle is hash-checked before/after, unchanged. Phase 6b/6c outputs
are not overwritten. No simulator normalization or physical water scale was
changed. All final CLI logs are empty, with no emitted fitting errors/warnings.

Commands executed from the repository root (output directory must be new):

```sh
pixi run python --version
pixi run fsl_mrs --help
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.mrscloud_tissue_relaxation_validation outputs/mrscloud_tissue_relaxation_validation/final
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_quantitative_tissue_relaxation.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_quantitative_tissue_relaxation.py tests/integration/test_tissue_relaxation_fitting.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q
```

The writable Matplotlib cache avoids repeated sandbox font-cache warnings.
The complete three executed CLI commands are saved in
[cli_commands.txt](../outputs/mrscloud_tissue_relaxation_validation/final/cli_commands.txt).
The mixed command is:

```sh
pixi run fsl_mrs --data outputs/mrscloud_tissue_relaxation_validation/final/mixed/metabolites.nii.gz --h2o outputs/mrscloud_tissue_relaxation_validation/final/mixed/water_reference.nii.gz --basis outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative/metabolite_basis.BASIS --output outputs/mrscloud_tissue_relaxation_validation/final/mixed/fsl_mrs_results --keep NAA Cr GPC PCh Glu --algo Newton --lorentzian --baseline off --ppmlim 0.2 4.2 --no_conj_basis --TE 30 --TR 2 --tissue_frac 0.2 0.3 0.5 --internal_ref NAA --wref_metabolite NAA --ref_protons 3 --ref_int_limits 1.8 2.2 --report
```

Tests cover exact saved Phase 6b default regression for all four voxels,
explicit no-attenuation, pure/mixed tissue sums, water/metabolite attenuation
ratios, TR dependence, scalar TE dependence and incompatible sequence TE
rejection, invalid profiles including unused tissues, concentration and water
linearity, receiver cancellation, exact determinism, unchanged nominal metadata,
standard acquisition headers and full-complex paired reload. Fitting tests
exercise independent predictions, all three actual reports and repeat fitting.

Initial analytical run: **23 passed** in 11.10 s. Focused analytical/fitting
run: **26 passed** in 22.30 s. Coefficient and prediction comparison tolerance
is 5e-6 relative; water tolerance is 1e-6. Analytical FID comparisons use
3e-15 relative plus at most 1e-20 absolute; exported samples compare exactly.

Full regression run: **216 passed, 1 xfailed**, 5 dependency warnings in
138.79 seconds. No failures. The expected failure remains the previously
documented installed FSL-MRS combined phase/frequency water-initialization
defect. Dependency warnings concern Dash/SWIG deprecations and a missing
Jupyter progress widget, not numerical fitting failures.

## Files and limitations

Production changes are limited to the opt-in quantitative adapter,
`simulation/tissue_relaxation.py` and optional standard acquisition timing
headers in `simulation/mrscloud_quantitative_pair.py`. Added example,
analytical and integration tests, this document and dedicated output artifacts.
Legacy SignalModel, GUI, residual-water model, MRSCloud source, raw water
reference, basis normalization and calibration files are unchanged.

This assumes 90-degree excitation and the prescribed saturation/echo product,
without exchange, magnetization transfer, inversion, sequence-specific steady
state or relaxation during shaped pulses. T2 is an amplitude attenuation
parameter here, not a request to change observed readout linewidth. Tissue
profiles are explicit synthetic inputs; their physical accuracy is not newly
established. PRESS TE is fixed by the bundle; no TE60 acquisition is simulated.
FSL's scalar metabolite correction cannot invert arbitrary heterogeneous
relaxation; the existing response/integration bias also remains. No empirical
inverse-alpha or inverse-fit gain was applied. No spatial B0, noise, GUI,
LCModel validation or physical absolute calibration was attempted. Stop at
Phase 6e and review the relaxation assumptions and model limitations.

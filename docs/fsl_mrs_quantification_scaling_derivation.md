# Phase 5d.1: derivation of FSL-MRS water-referenced concentration

## Result

The reported **6.003318300954079 mM** is independently reconstructed as
**6.003318300954086 mM**, a difference of **7.105427357601002e-15 mM**.
There is no missing empirical multiplier. Three different conventions account
for the systematic discrepancy, apart from tiny fitting errors:

1. NAA's selected real-spectrum integral is not three times the Ref0ppm-based
   synthetic single-proton response: response factor **0.78332772**.
2. FSL's pure-WM water model is **35.7565 M**, while the simulation specifies
   **40 M**: concentration factor **0.8939125**.
3. FSL corrects for T1/TR and TE/T2 attenuation absent from the synthetic input:
   water/metabolite relaxation factor **0.8573394853**.

No signal, basis, calibration constant, fitted dataset or FSL configuration
was modified. Proposed compatibility changes below are calculations only.

## Installed source trace

Installed version: **2.5.1.dev78+ga311e38f2**. Locations below are relative to
`.pixi/envs/default/lib/python3.14/site-packages/`:

| Source | Function / line | Role |
|---|---|---|
| fsl_mrs/core/mrs.py | rescaleForFitting, 578; scaling, 402; H2O, 192 | Data norm scaling; matching water scaling; coefficient interpretation |
| fsl_mrs/core/basis.py | get_formatted_basis, 264; _rescale_basis, 459 | Grid formatting, selection and fitting-time basis norm scaling |
| fsl_mrs/utils/misc.py | rescale_FID, 491; FIDToSpec, 42 | Norm-to-100 scaling; first-sample half-weighting and orthonormal FFT |
| fsl_mrs/models/model_lorentzian.py | forward, around 187 | Basis-weighted signal, fitted gamma/eps and phase |
| fsl_mrs/models/__init__.py | getFittedModel, 260 | Removes fitted Phi0/Phi1 for reference integration; retains gamma/eps |
| fsl_mrs/utils/results.py | calculateConcScaling, 230; getConc, 707 | Raw/internal/water scaling and multiplication of fitted coefficients |
| fsl_mrs/utils/quantify.py | WaterRef._fit_w_ref, 100 | Voigt water fit; temporary scaling restored; frequency/phase removed |
| fsl_mrs/utils/quantify.py | FIDIntegrator._calculate_area, 80; RefIntegral, 212 | Integrates abs(real spectrum) of fitted reference |
| fsl_mrs/utils/quantify.py | _calc_relax, 592; relax_corr_water_molar, 643; relax_corr_metab, 660; csf_corr_molar, 665 | Water-density, tissue and relaxation assumptions |
| fsl_mrs/utils/quantify.py | quantifyWater, 748, especially 812–837 | Molar scaling, proton ratio and division by reference coefficient |
| fsl_mrs/utils/constants.py | TISSUE_WATER_DENSITY, 73; H2O_MOLARITY, 82 | Tissue densities and pure-water molarity |
| nifti_mrs/axes.py | axis_indices / ppmShiftIndices | Nearest-bin, inclusive integration windows |
| fsl_mrs/scripts/fsl_mrs.py | h2o_scale, 164; assignment, 399–400 | Documented optional additional concentration scaling |

These are the actual functions used in Phase 5b, not an assumed mapping from
LCModel metadata. `CONC=1` is not used as evidence of physical concentration.

## Signal, fitting and integral equations

Let b_m be the source basis column and k the frozen full-reference mean
magnitude. The simulator uses B_m=b_m/k, q=abs(B_Ref0ppm[0]),
s_m=G*C_m*B_m, and s_w=G*C_w*2q*L_w, where C_m and C_w are mol/L.
The fit applies data multiplier s_D and basis multiplier s_B, giving fitted
coefficient a_m approximately G*C_m*s_D/(k*s_B). Therefore

\[
\widehat{G C_m}=a_m\frac{s_B}{s_D}k.
\]

At 10 mM, s_D=117.8317751259257, s_B=17925.38964873409,
k=6.573456836083183e-05 and a_NAA=1.0000000905278943. The resulting synthetic
coefficient is 0.010000000905278954 mol/L. FSL's unscaled source coefficient
a_m*s_B/s_D is 152.1269729860657; it is not itself mM.

Define the spectral functional over a specified inclusive-bin ppm window:

\[
J[x;W]=\operatorname{trapz}_{bins\in W}
\left|\Re\left[\operatorname{fftshift}
\operatorname{FFT}_{ortho}(h x)\right]\right|,
\quad h[0]=\tfrac12,\ h[n>0]=1.
\]

FSL integrates in spectral-bin units without a Hz-bin-width multiplier; the
common multiplier would cancel in the ratio. It uses **abs(real)**, not complex
magnitude, signed area, a peak height, or initial FID amplitude. The first point
is halved during the FFT, and the nearest window endpoints are included.

For the 10 mM fit, with fitted phase removed but intrinsic basis phase retained,

\[
\widetilde B_m=B_m e^{-(\gamma_m+i\epsilon_m)t},\quad
\widetilde L_w=e^{-\gamma_w t-\sigma_wt^2},
\]

\[
I_m=s_D\widehat{G C_m}J[\widetilde B_m;(1.8,2.2)],\quad
I_w=s_D\widehat A_w J[\widetilde L_w;(1.65,7.65)].
\]

The independent script implements these equations directly in NumPy. It
reconstructs the unit FIDs from the saved optimizer parameters and source
basis, not the saved integral scalars. The fitted water amplitude is reconstructed
from the saved optimizer amplitude divided by its temporary scaling and s_D.
It never calls `quantifyWater`, `calculateConcScaling` or `getConc`.

Numerical values:

| Term | Value |
|---|---:|
| q (independently calculated; roundoff may affect last digit) | 0.7717491793752124 |
| Fitted water amplitude, original scale | 61.73993439719796 |
| J[metabolite reference] | 97.74346950998739 |
| J[unit fitted water envelope] | 53.894794978562125 |
| Reconstructed I_m | 115.1728756196432 |
| Reconstructed I_w | 392080.6488209705 |
| I_m/I_w | approximately 0.00029375 |

The reconstructed areas differ from saved FSL integrals by 2.84e-14 and
−4.07e-10 respectively, consistent with floating-point operation ordering.

## Complete FSL concentration equation

For tissue fractions f_t and package tissue-water densities d_t:

\[
R_x=\frac{1-e^{-TR/T_{1,x}}}{1-\cos(FA)e^{-TR/T_{1,x}}}e^{-TE/T_{2,x}},
\qquad
Q_{eff}=55010\left(\sum_t f_t d_t R_{w,t}\right)
\frac{1}{R_m}\frac{1}{1-f_{CSF}}F_{additional}.
\]

Q_eff is in mM. In this actual experiment FA=90°, TE=.030 s, TR=2 s,
f_WM=1, f_GM=f_CSF=0, d_WM=.65, F_additional=1. Package pure-water molarity
is 55010 mM. T1_WM=.97 s, T2_WM=.073 s, T1_metabolite=1.29 s and
T2_metabolite=.194 s give:

- R_water,WM = 0.5786670548409867.
- R_metabolite = 0.6749567292246237.
- Relaxed tissue water = 20691.10854642174 mM.
- Q_eff = **30655.45930654141 mM**.

The water has n_w=2 protons and the selected NAA reference has n_ref=3.
For NAA referencing itself, final reported molarity is

\[
\boxed{\widehat C_{FSL}[mM]=\frac{I_m}{I_w}\frac{2}{3}Q_{eff}.}
\]

For another target j referenced to NAA, multiply this by a_j/a_NAA, in the
package's internal coefficient convention. `quantifyWater` divides by the
reference coefficient, then `getConc` multiplies by the target coefficient.
For NAA itself those coefficients cancel. There is no further factor of 1000
in this final equation because Q_eff already uses mM.

## Factor decomposition of the mismatch

With the declared water amplitude A_w=2qGC_w and known C_m, the exact fitted
result can be factored as

\[
\frac{\widehat C_{FSL}}{1000C_m}=
\underbrace{\frac{\widehat{GC_m}}{GC_m}}_{\text{coefficient fit}}
\underbrace{\frac{2qGC_w}{\widehat A_w}}_{\text{water fit}}
\underbrace{\frac{J[\widetilde B_m]}{3qJ[\widetilde L_w]}}_{\text{reference response}}
\underbrace{\frac{55010\sum f_td_t}{1000C_w}}_{\text{tissue water}}
\underbrace{\frac{\sum f_td_tR_{w,t}}{\sum f_td_t}\frac{1}{R_m}}_{\text{relaxation}}
\frac{F_{additional}}{1-f_{CSF}}.
\]

Here G=1. This factorization is derived from the signal and quantification
equations, not fit from the observed 0.6003318 ratio.

| Factor | Value | Cumulative prediction from 10 mM |
|---|---:|---:|
| Coefficient recovery | 1.0000000905278954 | 10.000000905278954 |
| Water-amplitude recovery | 0.9999999992358111 | 10.000000897637063 |
| Reference integral versus Ref0ppm | 0.7833277201353844 | 7.833277904497838 |
| Tissue water versus 40 M | 0.8939125 | 7.002265034804424 |
| Water/metabolite relaxation | 0.8573394853115213 | 6.003318300954086 |
| CSF correction | 1 | 6.003318300954086 |
| Additional correction | 1 | 6.003318300954086 |

The product is 0.6003318300954087. Thus neither inaccurate fitting, data gain,
nor water phase/frequency explains the discrepancy. The frozen normalization
k scales B_m and q identically and cancels in the reference-response factor.
The explicit water factor 2 cancels the matching synthetic water proton factor;
there is no missing extra water-proton multiplier. The reference factor includes
FSL's 3-proton NAA singlet convention and measured frozen response in its
window, relative to the Ref0ppm initial magnitude. It cannot by itself establish
that the basis has a missing proton multiplier. Source response, intrinsic
phase, sequence evolution and real-spectrum integration are combined in it.

Correcting tissue-water/relaxation assumptions alone would still leave the
ideal reference-response factor **0.7833276969013767**. Conversely, changing
only Ref0ppm calibration would not align tissue and relaxation assumptions.

## Proposed compatible convention, derived without changing signals

Let J_m^0=J[B_NAA;(1.8,2.2)] and J_w^0=J[L_w;(1.65,7.65)] be calculated
from the **unfitted frozen reference basis and specified water envelope**. With
no noise and no additional fitted nuisance terms, the amplitude that would
make the standard FSL equation return 1000*C_m is

\[
\boxed{A_{w,required}=G\frac{J_m^0}{J_w^0}\frac{2}{3}\frac{Q_{eff}}{1000}.}
\]

This derives from substituting I_m=GC_m J_m^0 and I_w=A_w J_w^0 into the
FSL equation and cancelling C_m. It uses neither fitted concentration nor
the observed 0.6003318 ratio. For the current grid, basis, integration windows,
zero water linewidth and default Phase 5b quantifier settings it is
**37.06444335191874 SRU**. This number is a *proposed FSL-model compatibility
amplitude*, not a new physical water calibration or a recommended unconditional
change to every simulated water reference.

Under the declared unattenuated synthetic assumptions and a quantifier whose
water concentration is actually 40 M, R_water=R_metabolite=1 and Q_eff=40000
mM. The corresponding required amplitude would instead be
**48.362600581241004 SRU**. The distinction demonstrates why compatibility
must declare acquisition/quantification assumptions explicitly. The original
Ref0ppm convention gives **61.73993435001698 SRU**, which remains unchanged.

## Can configuration alone resolve this?

Yes, for this explicitly synthetic NAA/reference/grid experiment, a documented
configuration adapter can leave the existing signals intact. Define from
the unfitted reference responses and declared parameters:

\[
\alpha=\frac{J_m^0}{3qJ_w^0},\qquad
F_{adapter}=\frac{1000C_w}{Q_{eff}\alpha}
=\mathbf{1.6657456248245763}.
\]

The existing CLI `--h2o_scale` maps directly to QuantificationInfo.add_corr;
the analogous Python property is `add_corr`. This is a place to declare a
derived response-convention adapter, rather than silently modifying the
simulator. The factor above was calculated with the current additional
correction set to 1 and would be a replacement value, not recursively applied.
That option affects both molarity and molality outputs; this derivation targets
molarity in mM and does not identify molality in mmol/kg as the same quantity.
The script evaluates predicted consequences only; **no corrected CLI run,
configuration file or water signal has been produced**.

This is not the empirical operation 10/reported_concentration. Tests deliberately
alter saved fitted coefficients, fitted water amplitudes and reported
concentrations in temporary copies, and verify the proposed constants remain
unchanged. Only source-basis response, known water concentration, grid, windows,
and the explicitly declared quantifier assumptions determine the adapter.

To remove physically inappropriate attenuation assumptions independently,
FSL offers T1/T2 inputs; its Python API also exposes `set_densities`. A
synthetic no-attenuation limit and effective water density 40/55.01 would
match the declared unattenuated 40 M water. That still leaves alpha unequal
to one. Selecting false tissue fractions, changing the integer proton count
to absorb alpha, or changing TE metadata solely to force agreement would not
be justified. A response adapter is more explicit.

**Recommendation:** retain the mathematically consistent simulator in this
phase. For a future development comparison, prefer a recorded, derived
FSL configuration adapter, with acquisition assumptions stated separately.
If a future opt-in simulator aims to reproduce FSL's native assumptions with
additional scaling exactly 1, use the response-integral amplitude formula
above and label its reference, grid, windows and quantifier convention. Neither
alternative establishes general cross-metabolite or physical calibration.

## Independent four-concentration verification

The script reconstructs all four saved datasets read-only. It uses saved fit
parameters and real source basis arrays, independently implements the FFT/bin
integration and relaxation equations, and reads the reported concentration
only after calculating its prediction, for comparison.

| Truth mM | Saved standard FSL mM | Reconstructed mM | Absolute difference mM | Relative difference |
|---:|---:|---:|---:|---:|
| 5 | 3.0016591504770394 | 3.001659150477043 | 3.55e-15 | 1.11e-15 |
| 10 | 6.003318300954079 | 6.003318300954086 | 7.11e-15 | 1.11e-15 |
| 15 | 9.0049773319647 | 9.004977331964712 | 1.24e-14 | 1.33e-15 |
| 20 | 12.006636601908157 | 12.006636601908172 | 1.42e-14 | 1.11e-15 |

Using unfitted ideal responses instead of saved nuisance parameters gives
6.0033175840117385 mM at 10 mM, within 7.17e-7 mM of the saved result. This
small difference is predicted by coefficient/line-shape fit errors, not an
unresolved scale. The proposed adapter predicts 10.000001194243577 mM from
the actual fitted result, without correcting the known small optimizer error.

## Reproduce, outputs and tests

```
pixi run python -m examples.fsl_mrs_quantification_scaling_derivation outputs/fsl_mrs_quantification_scaling_derivation
pixi run python -m pytest -q tests/test_fsl_mrs_quantification_scaling_derivation.py
pixi run python -m pytest -q tests
```

The input is the existing `outputs/fsl_mrs_water_scaling_phase5b` series; an
alternative read-only input directory can be supplied with `--source-dir`.
Output contains `derivation.json`, `comparison.csv`, and `factors_10mM.csv`.
No simulator generation or FSL fitting is performed by the derivation script.
Its final-quantification independence is tested by replacing FSL's
quantifyWater/getConc functions with exceptions. Other tests check all four
reconstructions, linearity, response-derived proposals, first-point/integral
units, deterministic calculations and unchanged hashes of existing inputs.

Focused tests: **5 passed**. The full suite completed with **72 passed,
1 xfailed**, zero unexpected failures, in 45.74 s. The expected failure remains
the installed FSL-MRS nonzero-phase/nonzero-frequency initialization defect.
Five existing dependency deprecation/optional-display warnings were emitted.

## Remaining uncertainties

The MRSCloud source response has no demonstrated scanner/molar calibration.
The response adapter is specific to NAA, its integration limits, reference
basis identity, spectral grid, water line shape, software integration and
quantification assumptions. A water-reference area may change with linewidth
and integration limits. This experiment cannot prove one response constant
will quantify every metabolite correctly; additional independent references
and basis provenance would be needed. No correction is implemented in this
phase, and no LCModel validation, tissue mixing or new relaxation model was
introduced.

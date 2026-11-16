# Phase 7a: controlled complex noise and Voigt fitting

## Completed experiment

A noise-free Voigt baseline and **120 actual noisy fits** were executed: 20
independent realizations at each simulator SNR 100, 50, 25, 15, 10 and 5.
Water target SNR is 1000 for every finite-SNR acquisition, with a separate
independent realization. Five actual CLI HTML reports were generated. A full
restart reused all 121 completed fits and verified regenerated acquisition
arrays against the exported inputs exactly. Original clean acquisition and
Phase 6a calibration hashes remain unchanged.

All optimizers and water quantifications returned success, but precision and
individual recovery deteriorate substantially at low SNR. The established
water-scaling convention bias is preserved, not corrected using these fits.
This establishes a reproducible synthetic noise benchmark, not physical
scanner calibration or validation against LCModel.

## Noise model and implementation boundary

For each acquisition after deterministic signal generation:

\[
s_{noisy}[n]=s_{clean}[n]+\eta_R[n]+i\eta_I[n],\qquad
\eta_R,\eta_I\stackrel{independent}{\sim}\mathcal N(0,\sigma^2).
\]

Sigma is the SD of **one real or imaginary FID component**, in the same
arbitrary receiver signal units as the clean FID. Thus E[abs(eta)^2]=2*sigma^2;
complex magnitude SD is not sigma. All samples use the same white receiver
noise model. No apodization, noise renormalization, additional gain or proton
factor is applied.

`simulation/mrscloud_quantitative_noise.py` provides:

- `complex_gaussian_noise(points, sigma=..., seed=...)`: reusable complex draws.
- `reference_peak(fid, config, window_ppm)`: declared clean spectral reference.
- `add_pair_noise(pair, seed=..., metabolite_sigma=... OR metabolite_snr=...,
  water_sigma=... OR water_snr=...)`: separate output with clean signals retained.

Sigma must be finite and nonnegative; target SNR positive or infinite. Missing
specifications leave that acquisition noiseless. Infinite SNR gives sigma zero
and is recorded as a string instead of nonstandard JSON Infinity. Target
metabolite SNR requires a nonzero clean NAA component. Explicit sigma supports
other metabolite selections. Seeds must be explicit nonnegative integers.
Applying the wrapper twice to an already marked acquisition is rejected.

Noise uses local NumPy PCG64 generators, with `SeedSequence(root_seed).spawn(2)`:
metabolite stream [0], water stream [1]. Real and imaginary draws are independent
normal vectors within each stream. No global NumPy seed/state is changed and no
Torch random API is used. Generator, NumPy version, root seed, child spawn keys,
sigma, target SNR, reference peaks/windows and realized SNR are recorded.

The pair, voxel and spatial APIs accept optional `noise={...}`. It is applied
only at the outer acquisition boundary: after concentration summation,
relaxation and B0 modulation. The spatial adapter calls its internal voxel
adapter without noise, then adds one noise vector to each integrated FID;
it does not accumulate independent receiver noise per anatomical cell.

```python
from simulation.mrscloud_quantitative_voxel import simulate_voxel
from examples.mrscloud_quantitative_voxel import BUNDLE, METABOLITES, WATER

pair = simulate_voxel(
    {'WM': 1}, calibration_bundle=BUNDLE,
    metabolite_concentrations_mM=METABOLITES,
    water_concentrations_mol_l=WATER,
    noise={'seed': 70000, 'metabolite_snr': 50, 'water_snr': 1000})
```

The same keyword works with `generate_pair` in analytical and sequence modes,
and `simulate_spatial` with optional Phase 6e relaxation. For a standalone
Phase 6a `synthesize(...)` result, use `add_pair_noise` on the returned pair;
the frozen library implementation is unchanged. Direct SD configuration uses
`metabolite_sigma`/`water_sigma` instead of the corresponding SNR fields.

Default `noise=None` preserves old arrays, metadata and export behavior.
Explicit zero-noise settings preserve samples but record their provenance.
The wrapper does not mutate its input. Components remain **clean ground-truth
responses**, and `clean_metabolite_fid` / `clean_water_reference_fid` are retained
alongside the noisy primary FIDs. Existing spectrum fields are updated to the
noisy FIDs using the repository's unnormalized FFT convention; the SNR
calculation below deliberately uses a separately specified orthonormal FFT.

## Simulator SNR convention

For N samples, define unweighted unitary Fourier samples

\[
X[k]=\operatorname{fftshift}(\operatorname{FFT}(s))[k]/\sqrt N,
\qquad P_{NAA}=\max_{ppm\in[1.8,2.2]}|X_{NAA,clean}[ppm]|,
\]
\[
SNR_{sim}=P_{NAA}/\sigma.
\]

The reference is the **isolated clean NAA component** at its actual acquisition
concentration, not the noisy total spectrum, a fitted peak or a first FID sample.
The magnitude spectrum is used, so intrinsic basis phase is not silently
removed. There is no first-point half weighting, zero filling or matched filter.
The fixed window is the NAA acetyl region; the controlled sweep has zero B0.
Peak height depends on sampling/bin placement, which is part of this convention.

For circular white complex noise, an orthonormal transform preserves sigma^2
in each real/imaginary component, so spectral and FID per-component sigma are
identical in expectation. With an unnormalized FFT both peak and SD acquire
sqrt(N); their ratio is unchanged if conventions are applied consistently.
Water target SNR independently uses its own clean magnitude peak at 4.1–5.3 ppm.
There is no water-to-metabolite calibration multiplier in either SNR definition.

In this acquisition P_NAA=3.2090203050950436e-5. At target 100, sigma is
3.2090203050950437e-7; at target 5, it is 6.418040610190087e-6.
P_water=.09152146403516252, giving fixed sigma_water=9.152146403516253e-5
for water SNR 1000. All clean amplitudes are identical between conditions.

Achieved SNR is independently estimated using the actual difference between
exported noisy data and clean data: compute a unitary FFT, estimate
sqrt((var(real)+var(imag))/2) with ddof=1, and divide the independently calculated
clean peak by it. This fluctuates statistically; draws are not rescaled to force
exact realized SNR. A direct DFT of the clean reference bin independently checks
the production FFT reference peak. Metadata reports both target and realized SNR.

### Distinction from FSL-MRS SNR

Installed FSL-MRS `utils/qc.py` uses fitted-component real peak heights,
matched filtering/apodization and a detrended identified noise region.
`noiseSD` multiplies real-spectrum SD by sqrt(2). Its residual SNR uses
max abs(real(data-baseline)) / (2*sqrt(MSE)). FSL's spectrum transform also
half-weights sample zero. These are different definitions. Mean reported
FSL NAA peak SNR is about **52.28** at simulator SNR 100 and **2.88** at
simulator SNR 5. Both are saved separately; equality is not an acceptance test.

## Clean configuration, calibration and Voigt baseline

The exact Phase 6b/6c pure WM acquisition is used and compared sample-for-sample
against `outputs/mrscloud_quantitative_voxel/final/pure_WM`. Nominal NAA/Cr/GPC/
PCh/Glu are 10/8/2/1/8 mM per total tissue volume; water is **39.0165 M**, not
40 M. Receiver G=1, no added relaxation, no B0, no noise in the baseline.
PRESS TE30 uses 2048 points at 2000 Hz, .0005 s dwell, 127.731 MHz and origin
4.68 ppm. The common response scalar remains 8857.353925964408 and identity
`mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81`.
Basis proton responses and raw pre-DC water processing are unchanged.

Runtime FSL-MRS version: `2.5.1.dev78+ga311e38f2`; NumPy `2.5.4`, Pixi Python
3.14.8. Python fitting explicitly uses model='voigt', Newton/TNC, baseline off,
one shared frequency/linewidth group, ppm .2–4.2 and the matching five selected
columns of the complete basis. CLI omits `--lorentzian` and `--free_shift`:
installed `scripts/fsl_mrs.py:336–344` explicitly selects Voigt in that branch.
No initialization changes, warm starts across realizations, data alignment or
basis frequency shifts are used. The clean baseline is checked before any
finite-SNR sweep; failure to recover coefficients within 5e-6 relative stops it.

Voigt metabolite attenuation is exp(-gamma*t-sigma_model^2*t^2). The model's
sigma_model is not the acquisition-noise sigma. It is also different from
WaterRef's Gaussian parameter, which enters exp(-sigma_water_model*t^2).

| Metabolite | True coefficient mM-equivalent | Fitted coefficient mM-equivalent | Relative error | Standard FSL mM |
|---|---:|---:|---:|---:|
| NAA | 10 | 10.000002737823 | 2.73782e-07 | 7.059366426934 |
| Cr | 8 | 8.000002321679 | 2.9021e-07 | 5.647493234321 |
| GPC | 2 | 1.999997410311 | -1.29484e-06 | 1.411871070685 |
| PCh | 1 | 1.000002170525 | 2.17052e-06 | 0.705937981673 |
| Glu | 8 | 7.999997789934 | -2.76258e-07 | 5.647490035197 |
| GPC+PCh | 3 | 2.999999580836 | -1.39721e-07 | 2.117809052358 |

Maximum relative coefficient error is **2.171e-6**, full complex residual
**6.042e-7**. Additional gamma=0 s^-1, sigma_model=.0025052457010807704 s^-1,
epsilon=-8.420409222451153e-7 rad/s, Phi0=2.0470794598796493e-7 rad and
Phi1=5.515286048169447e-10 rad/Hz. The original readout linewidth remains in
the basis; these parameters are additional fitted broadening.

Expected clean initial water amplitude is .007887908291571334; fitted amplitude
.00788790812234321, relative error -2.145e-8 and residual 1.432e-7. Both
optimizers converge. The non-unit standard concentration is the established
quantification convention, not a failed clean coefficient fit.

## Monte Carlo method and statistics

Full mode: one noise-free baseline, then 20 draws at each level
[100,50,25,15,10,5]. For zero-based level index l and realization r,
root_seed=70000+100*l+r. Metabolite and water use distinct child streams.
Seeds are unique across the 120 noisy acquisitions; noise is not shared between
levels. Water SNR is fixed at 1000, with independent water noise per acquisition.

Per fit, results contain expected/internal fitted coefficients and normalization,
signed/relative errors, all nuisance parameters, noisy-data residual, fitted
complex error to clean data, optimizer diagnostics, warnings, water fit,
standard concentrations, independent predictions and FSL uncertainty estimates.
A coefficient is recovered from internal amplitude a through
1000*a*s_basis/(s_data*G). Expected internal amplitude is C_mM*1e-3*G*s_data/s_basis.

Bias is mean(fitted-truth), SD uses n-1, RMSE uses mean squared error, and
median absolute error is median(abs(fitted-truth)). Every returned numerical
fit is included, including non-success exits. Exceptions remain recorded and
count as failures; they do not silently disappear from denominators. Failure
rate is optimizer exception/non-success; quantification failures are separately
counted. A descriptive poor-recovery fraction counts abs(relative error)>20%,
without excluding those fits from any statistics.

Coefficient uncertainty uses the FSL covariance in internal parameter units,
converted to mM-equivalent units. GPC+PCh variance includes both variances and
**twice the covariance**, not independent-error quadrature. Reported covariance
matrices and name ordering are retained. Approximate coefficient Wald intervals
are fitted +/-1.96*reported SD. Coverage uses available finite positive SDs;
no-noise coverage is inapplicable and not fabricated. These are not fully
propagated water-scaled concentration intervals.

### Full Monte Carlo numerical results

All entries except coverage/poor recovery are mM-equivalent. Each finite-SNR
row contains 20 returned fits, with 20 valid coefficient intervals.

| SNR | Metabolite | Bias | SD | RMSE | Median absolute error | 95% interval coverage | Error >20% |
|---:|---|---:|---:|---:|---:|---:|---:|
| 100 | NAA | 0.010260 | 0.076119 | 0.074898 | 0.051529 | 95% | 0% |
| 100 | Cr | 0.007092 | 0.039833 | 0.039467 | 0.032030 | 100% | 0% |
| 100 | GPC | -0.009086 | 0.031268 | 0.031802 | 0.020884 | 100% | 0% |
| 100 | PCh | 0.011607 | 0.036225 | 0.037166 | 0.033667 | 95% | 0% |
| 100 | Glu | -0.018120 | 0.167674 | 0.164429 | 0.134412 | 95% | 0% |
| 100 | GPC+PCh | 0.002520 | 0.019378 | 0.019054 | 0.013913 | 100% | 0% |
| 50 | NAA | 0.005969 | 0.140739 | 0.137305 | 0.066722 | 95% | 0% |
| 50 | Cr | 0.040910 | 0.103169 | 0.108560 | 0.073592 | 100% | 0% |
| 50 | GPC | 0.007314 | 0.072550 | 0.071091 | 0.048116 | 100% | 0% |
| 50 | PCh | 0.005230 | 0.072321 | 0.070683 | 0.057972 | 100% | 0% |
| 50 | Glu | -0.056826 | 0.350531 | 0.346349 | 0.250451 | 100% | 0% |
| 50 | GPC+PCh | 0.012544 | 0.028696 | 0.030654 | 0.029042 | 100% | 0% |
| 25 | NAA | 0.007803 | 0.235694 | 0.229858 | 0.211376 | 100% | 0% |
| 25 | Cr | 0.024519 | 0.223963 | 0.219665 | 0.178068 | 100% | 0% |
| 25 | GPC | 0.045985 | 0.218696 | 0.218063 | 0.114919 | 85% | 10% |
| 25 | PCh | -0.003305 | 0.186731 | 0.182032 | 0.099617 | 85% | 35% |
| 25 | Glu | 0.180508 | 0.590450 | 0.603144 | 0.255693 | 100% | 0% |
| 25 | GPC+PCh | 0.042680 | 0.088721 | 0.096433 | 0.047397 | 85% | 0% |
| 15 | NAA | -0.029395 | 0.374968 | 0.366654 | 0.217975 | 100% | 0% |
| 15 | Cr | 0.006511 | 0.344312 | 0.335657 | 0.175826 | 95% | 0% |
| 15 | GPC | -0.068560 | 0.225515 | 0.230249 | 0.117023 | 95% | 10% |
| 15 | PCh | 0.035279 | 0.252645 | 0.248762 | 0.153521 | 95% | 35% |
| 15 | Glu | 0.392416 | 1.325162 | 1.349905 | 0.950664 | 90% | 30% |
| 15 | GPC+PCh | -0.033281 | 0.122381 | 0.123838 | 0.102138 | 100% | 0% |
| 10 | NAA | 0.086349 | 0.759322 | 0.745116 | 0.727969 | 95% | 0% |
| 10 | Cr | 0.165458 | 0.555499 | 0.566150 | 0.390885 | 100% | 0% |
| 10 | GPC | -0.025569 | 0.385384 | 0.376495 | 0.248120 | 100% | 35% |
| 10 | PCh | 0.101432 | 0.416063 | 0.418020 | 0.329735 | 95% | 75% |
| 10 | Glu | 0.292398 | 1.652992 | 1.637455 | 1.158451 | 100% | 40% |
| 10 | GPC+PCh | 0.075863 | 0.160995 | 0.174295 | 0.105209 | 100% | 0% |
| 5 | NAA | 0.396038 | 1.239756 | 1.271610 | 0.694177 | 95% | 5% |
| 5 | Cr | 0.245338 | 1.603058 | 1.581612 | 0.720196 | 90% | 15% |
| 5 | GPC | -0.212008 | 0.722292 | 0.735233 | 0.607708 | 100% | 65% |
| 5 | PCh | 0.477106 | 0.592658 | 0.749207 | 0.680592 | 100% | 75% |
| 5 | Glu | 0.975562 | 4.456976 | 4.452317 | 3.333927 | 90% | 80% |
| 5 | GPC+PCh | 0.265098 | 0.519066 | 0.571170 | 0.479409 | 100% | 40% |

All 120 noisy metabolite optimizers and all 120 quantifications returned success:
**0% optimizer/exception failure, 0% quantification failure** at every level.
This does not mean accurate recovery. At SNR 5, 75% of PCh and 80% of Glu fits
have more than 20% coefficient error. Grouped GPC+PCh is substantially more
stable than the individual overlapping components; at SNR 10 its SD=.160995 mM
compared with .385384 for GPC and .416063 for PCh. Noise-related broadening,
phase uncertainty, parameter overlap and nonnegative constraints remain part
of the actual nonlinear estimator; their biases are not removed.

Coverage ranges 85–100% in this seed set. Twenty draws provide limited precision:
at nominal 95% coverage, binomial standard error is approximately 4.9 percentage
points. These observations do not establish exact interval calibration, nor
justify rescaling reported uncertainties. Baseline covariance reflects numerical
residuals rather than stochastic variance and is excluded from coverage.

## Water amplitudes and quantification

The clean theoretical FSL prediction uses the same independent Phase 6c
unit-response integration and declared quantification assumptions. It is
unchanged across draws. In particular, no concentrations returned by fitting
are used to calibrate the water amplitude or target noise.

Standard FSL quantification uses NAA's 3-proton reference window 1.8–2.2 ppm,
water's 2-proton corrected-model window 1.65–7.65 ppm, TE=.03 s, TR=2 s,
pure-WM fractions, installed water density .65 and standard 3T relaxation.
The synthetic signals have no corresponding T1/TR or TE/T2 attenuation; this
is the previously established convention mismatch. The independent clean NAA
prediction is **7.059362383312279 mM**, compared with baseline fitted
**7.059366426933524 mM** and nominal 10 mM.

For each noisy fit, an additional independent numerical reconstruction uses
its returned coefficient and nuisance parameters to rebuild NAA and corrected
water model FIDs and integrate them directly. It does not call final reporting.
This conditional reconstruction checks software interpretation, not unbiased
recovery. It matches reported values within **1.422e-14 mM**. Actual concentrations,
differences from nominal, differences from independent clean prediction, and
differences from the fitted noise-free reference are all saved per realization.
No zero-B0 calibration factor is tuned using the Monte Carlo output.

Across all 120 noisy water fits, mean initial-amplitude relative error is
**-0.03464%**, sample SD **0.10210%**, maximum absolute error **0.34556%**.
The input clean water amplitude is exactly the same for every draw. These
small but nonzero water-fitting errors contribute to empirical quantification
variability; the experiment does not claim perfectly noiseless water.

Mean standard FSL concentration below is followed by its sample SD in brackets,
all mM. The noiseless row contains one deterministic fit and has no sample SD.

| Simulator SNR | NAA | Cr | GPC | PCh | Glu | GPC+PCh |
|---|---:|---:|---:|---:|---:|---:|
| infinite | 7.059366 | 5.647493 | 1.411871 | 0.705938 | 5.647490 | 2.117809 |
| 100 | 7.074544 (0.056742) | 5.658833 (0.029047) | 1.407043 (0.022794) | 0.714918 (0.025267) | 5.640970 (0.116354) | 2.121962 (0.013524) |
| 50 | 7.077354 (0.105067) | 5.687468 (0.079804) | 1.419756 (0.050487) | 0.711035 (0.051575) | 5.618254 (0.248030) | 2.130791 (0.020413) |
| 25 | 7.092034 (0.182938) | 5.686546 (0.168262) | 1.449732 (0.154000) | 0.706361 (0.132698) | 5.796766 (0.418027) | 2.156093 (0.063091) |
| 15 | 7.062166 (0.282810) | 5.670778 (0.250824) | 1.367983 (0.159892) | 0.733182 (0.178328) | 5.944084 (0.940230) | 2.101165 (0.087199) |
| 10 | 7.163691 (0.553177) | 5.799598 (0.411446) | 1.401793 (0.271459) | 0.782540 (0.296179) | 5.889201 (1.172965) | 2.184333 (0.115450) |
| 5 | 7.461687 (0.923483) | 5.918384 (1.172777) | 1.282121 (0.518051) | 1.060578 (0.425681) | 6.442443 (3.207736) | 2.342698 (0.376067) |

Thus NAA coefficient mean bias grows to about +.396 mM at SNR 5 while its
standard concentration mean grows from 7.05937 to 7.46169 mM. The latter still
falls below nominal 10 mM because the pre-existing convention bias remains.
Fitting variability, estimator bias and response/quantification convention
are separate effects; a seemingly closer reported concentration at low SNR
is not improved calibration.

## Reproducibility, performance and outputs

`examples/mrscloud_noise_fitting_validation.py OUTPUT_DIR` runs full mode;
`--smoke` runs two draws each at SNR 50/10 plus the clean baseline (5 fits).
The complete basis is read once for fitting and copied in memory per fit;
the clean pair and unit responses are generated/loaded once, not regenerated
for each realization. No Octave/MRSCloud simulation is rerun.

An existing output directory is restartable only if configuration, implementation
hashes, software versions and original source hashes match exactly. Each
completed fit is saved atomically; unfinished fits are retried on restart.
Input NIfTI hashes are checked before accepting cached fits, and seeded
regenerated noisy arrays are compared exactly against the saved acquisitions.
An incompatible configuration is rejected rather than silently reusing fits.
An incomplete/corrupt acquisition export is rejected and requires a fresh case
or output directory; it is not guessed from partial files.

The full run completed 121 fits. Running the same command again reused **121
cached fits**, with unchanged records/statistics and no additional fitting.
The smoke integration test also forbids the optimizer on restart to verify
that completed fits are actually skipped. A separately rerun noisy fit returns
exactly equal numerical results. Existing clean datasets and calibration bundle
hashes are checked before/after and remain unchanged.

Per-realization datasets contain:

```
metabolites.nii.gz            # noisy primary, clean for infinite baseline
water_reference.nii.gz       # independent noisy water, clean for baseline
metabolites_clean.nii.gz      # preserved noiseless acquisition when noise requested
water_reference_clean.nii.gz
ground_truth.json
simulation_config.json       # complete noise parameters and seeds
fit_result.json
fitted_arrays.npz
```

Default noise-off export still writes the original two NIfTI-MRS files and two
sidecars. No custom NIfTI calibration fields are introduced. Export tests compare
all four clean/noisy complex FIDs exactly on reload; the full restart also
compares every regenerated primary acquisition against exported data exactly.
Root `validation_config.json`, `results.json`, `coefficients.csv`, `summary.csv`,
figures and CLI commands provide machine-readable provenance and comparisons.

## Actual CLI reports and plots

The installed CLI selected Voigt by omitting `--lorentzian`; other settings
match Python exactly. All five actual HTML reports exist, are nonempty and
contain water referencing. CLI/API standard concentrations agree exactly for
all reported metabolites in these representative runs. Realization 00 was
used for all five reports; the general runner selects the first complete,
optimizer-successful quantified fit and retains every realization in statistics.

- [SNR infinite, realization 0](../outputs/mrscloud_noise_fitting_validation/full/SNR_infinite/realization_00/fsl_mrs_results/report.html)
- [SNR 50, realization 0](../outputs/mrscloud_noise_fitting_validation/full/SNR_50/realization_00/fsl_mrs_results/report.html)
- [SNR 25, realization 0](../outputs/mrscloud_noise_fitting_validation/full/SNR_25/realization_00/fsl_mrs_results/report.html)
- [SNR 10, realization 0](../outputs/mrscloud_noise_fitting_validation/full/SNR_10/realization_00/fsl_mrs_results/report.html)
- [SNR 5, realization 0](../outputs/mrscloud_noise_fitting_validation/full/SNR_5/realization_00/fsl_mrs_results/report.html)


![Bias, precision, RMSE, coverage and water-scaled concentration](../outputs/mrscloud_noise_fitting_validation/full/noise_bias_precision.png)

![Water amplitude and fitted-signal errors, realized SNR](../outputs/mrscloud_noise_fitting_validation/full/noise_water_and_signal_errors.png)

Plots were inspected; display scaling does not change numerical calculations.
Exact CLI invocations are saved in
[cli_commands.txt](../outputs/mrscloud_noise_fitting_validation/full/cli_commands.txt).
For example, the actual SNR 25 invocation was:

```sh
pixi run fsl_mrs --data outputs/mrscloud_noise_fitting_validation/full/SNR_25/realization_00/metabolites.nii.gz --h2o outputs/mrscloud_noise_fitting_validation/full/SNR_25/realization_00/water_reference.nii.gz --basis outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative/metabolite_basis.BASIS --output outputs/mrscloud_noise_fitting_validation/full/SNR_25/realization_00/fsl_mrs_results --keep NAA Cr GPC PCh Glu --algo Newton --baseline off --ppmlim .2 4.2 --no_conj_basis --TE 30 --TR 2 --tissue_frac 1 0 0 --internal_ref NAA --wref_metabolite NAA --ref_protons 3 --ref_int_limits 1.8 2.2 --report
```

Executed with the repository's Pixi environment:

```sh
pixi run python --version
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.mrscloud_noise_fitting_validation outputs/mrscloud_noise_fitting_validation/smoke_validated --smoke
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.mrscloud_noise_fitting_validation outputs/mrscloud_noise_fitting_validation/full
# Same full command repeated to verify restart (121 cache hits).
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_quantitative_acquisition_noise.py tests/test_noise_monte_carlo_summary.py tests/integration/test_noise_voigt_validation.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q
```

The exploratory smoke folder predates addition of covariance metadata and the
final cache signature; it is retained as exploratory output, and the final
full output above is authoritative. Final smoke tests use fresh temporary
directories. Reusing a folder from a changed implementation correctly fails
its signature check rather than accepting stale results.

## Tests and source changes

Focused run: **15 passed** in 20.04 s. Tests cover 200,000-sample noise means,
real/imaginary variance/correlation and complex power; independent channel
streams; deterministic seeds; unchanged global NumPy/Torch states; direct-DFT
SNR calibration; target SNR versus realized statistics; default source regression;
analytical/sequence pair modes; voxel/spatial post-integration noise; optional
relaxation compatibility; all clean/noisy complex exports; invalid noise inputs;
hand-calculated summary statistics; exception/failure denominators; actual
Voigt fits and three smoke CLI reports; covariance-correct group SDs;
independent quantification; deterministic refitting and restart.

Full regression: **275 passed, 1 xfailed**, 4 dependency warnings in
196.86 seconds. No failures. The expected failure remains the documented
installed WaterRef combined phase/frequency initialization defect. Warnings
are external Dash/SWIG deprecations. The full run includes final post-relaxation
spatial noise and global-state assertions plus all prior phase regressions.

Changed existing files are `simulation/mrscloud_quantitative_pair.py`,
`simulation/mrscloud_quantitative_voxel.py`, and
`simulation/mrscloud_quantitative_spatial.py`. Added
`simulation/mrscloud_quantitative_noise.py`,
`examples/mrscloud_noise_fitting_validation.py`,
`tests/test_quantitative_acquisition_noise.py`,
`tests/test_noise_monte_carlo_summary.py`,
`tests/integration/test_noise_voigt_validation.py`, and this document.
The frozen library, normalization/calibration files, raw water response,
legacy SignalModel, GUI and environment dependencies were not modified.

## Limitations and review boundary

Only one homogeneous WM, zero-B0, no-additional-relaxation acquisition was swept.
The optional moderate spatial B0 experiment was not performed; prior severe
B0 mismatch results are not recast as noise errors. Twenty realizations per
level provide limited bias/coverage precision. Noise is ideal circular white
Gaussian receiver noise, not correlated coil-combination noise, colored
filter noise, drift or physiological variability. Target SNR is sampling- and
reference-window-specific and differs from FSL/LCModel conventions. All true
coefficients and water amplitudes retain their original synthetic calibration;
physical concentration calibration is not established. No uncertainty, noise
or water gain was adjusted against desired fit outcomes. No LCModel validation,
GUI integration or further spatial model development was attempted. Stop at
Phase 7a and review these statistical results before extending the benchmark.

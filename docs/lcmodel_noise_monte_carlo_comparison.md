# Phase 7c: LCModel versus FSL-MRS noise Monte Carlo

All **121 LCModel acquisitions were executed**: one clean acquisition and20
realizations at each simulator SNR100,50,25,15,10,5. All121 produced valid final
numerical outputs; **zero fatal/software/incomplete failures**. There were121
actual attempts in the full experiment, no retries or changes of parameters.
A subsequent restart reused **121 cached fits** without invoking LCModel.
A separate five-acquisition smoke run is not included in Monte Carlo counts.

The protocol was frozen before fitting and was not tuned. Signal, noise,
calibration, basis duration, water gain and original FSL fits were unchanged.
The clean LCModel results reproduce Phase7b exactly: NAA10.639,Cr8.446,
GPC2.125,PCh1.051,Glu8.376mM. These are toolbox results, not ground truth.

The high-SNR noisy mean does not simply preserve the clean offset. At SNR100,
NAA mean is9.96095mM, SD0.29690mM, a **-0.67805mM change from clean LCModel**.
Automatic regularization/baseline behavior and water integration differ between
clean and noisy data under the same settings. Nearer nominal concentration is
not proof of improved calibration. Individual GPC/PCh recovery deteriorates
strongly at low SNR; their sum is more stable. There is no physical calibration
claim and no inference that a successful optimizer/output guarantees accuracy.

## Frozen environment, signals and protocol

LCModel runtime version **6.3-1N**, executable
`$LCMODEL_BINARY`. Existing FSL-MRS runtime
**2.5.1.dev78+ga311e38f2**, NumPy2.5.4; Python3.14.8 in
Pixi. The supplied LCModel manual and Phase7b/7b.1 source traces were used;
LCModel was neither installed nor modified.

Exact original CONTROL SHA256:
`c944f45c493f00265307270fe7995a34142f0c0cf28c27e97616d13b22310be1`.
Path-independent CONTROL SHA256:
`6df398958a4283357fcf1bf871c0de52c2f69a5f3a6fddcaba51a93a1c2141ae`.
Unchanged compatible BASIS SHA256:
`9ec8bf2992ac6471aeeaf13d2af8ea3fb070791844830f51228dae4ddcbad629`.

The complete reference CONTROL is copied to `frozen_protocol.CONTROL` and
all acquisition CONTROL files differ **only** in FILRAW,FILH2O and output paths.
FILBAS still points to the identical Phase7b compatible BASIS. Settings include
WCONC39016.5mM,ATTH2O1,ATTMET1,NAA reference,N1HMET3,WSPPM2.01,
FWHMBA1/127.731ppm,NSIMUL0,NRATIO0,original preliminary list,default baseline
and line-shape regularization. No VITRO,DEEXT2 or experimental FWHMBA override
was introduced. LCModel's native bandwidth conversion,4096-point basis retention
and intentional data zero-filling were preserved exactly.

The inherited CONTROL TITLE and RAW ID refer to the Phase7b clean protocol;
they were preserved to honor path-only changes. The SNR directory and per-case
JSON identify the actual noisy acquisition. They are not evidence of clean data.

Inputs are the exact saved acquisitions under
`outputs/mrscloud_noise_fitting_validation/full/`. Their original hash records,
individual FSL results and aggregate Phase7a records agree. The inventory rejects
missing/duplicate/unexpected realizations, mismatching seeds, metadata or hashes.
Ground truth is WM,NAA10,Cr8,GPC2,PCh1,Glu8mM; water39.0165mol/L; gain1;
PRESS TE30,2048complex points,2000Hz bandwidth,127.731MHz,origin4.68ppm,
water4.70ppm. There is no added relaxation,B0 or residual water.

The same Phase6a calibration identity is retained:
`mrscloud-press-te30-full-reference-max-v1:d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81`.
The frozen provenance protects **1118 files**, including
all121 acquisitions, FSL result/report artifacts, earlier LCModel fits and
calibration/binary files. Hashes agree before/after execution and restart.

## Noise, conversion and SNR

No benchmark noise was regenerated. The saved experiment uses independent circular complex
Gaussian receiver noise, SD sigma per real/imaginary component; complex power
is2sigma². The noise is added after deterministic acquisition formation.

Simulator SNR is the maximum magnitude of the **clean isolated NAA** unitary
spectrum in1.8–2.2ppm divided by the per-component spectral noise SD. The FFT is
`fftshift(fft(fid,norm='ortho'))`, without first-point weighting. It is different
from either toolbox's reported SNR. Root seeds are
`70000+100*level_index+realization`, with separate child streams0/1 for metabolite
and water. The stored water target SNR is1000 for all finite-SNR cases.

Actual realized SNR is independently checked from saved noisy-minus-clean FIDs,
using both real and imaginary unitary spectral variances. Target,realized SNR,
seed,sigma and full calibration metadata are retained in per-case records.
Nothing is scaled to achieve a post-hoc target.

Conversion reuses Phase7b: store conjugated FSL-interpreted complex FIDs as
real/imaginary pairs, `(2E25.16)`, unchanged time order and no gain normalization.
Independent parsing and comparison over all242 converted arrays gives:

- Maximum absolute complex error **0**.
- Maximum relative complex-array error **0**.
- Initial amplitude and complex metabolite/water ratio errors **0**.
- Maximum phase-comparison rounding error **5.52e-17rad**.

The sign reversal is the declared LCModel conjugation convention, not a new
frequency offset. Sampling/reference metadata stay fixed in the CONTROL.

## Execution, restart and failure definitions

`examples/lcmodel_noise_comparison.py` uses the installed executable with CONTROL
on stdin, captured stdout/stderr, a **180-second timeout**, per-attempt output
directories and atomic completion records. Full-fit executable time summed to
**102.716s**, maximum1.747s per case; this excludes conversion/analysis I/O.

The cache validates frozen protocol,basis,binary,execution/conversion source
hashes, exact input hashes, converted samples and retained output hashes. It
also reparses completed outputs to reject inconsistent cached result records.
Completed failures remain failures; they are not silently retried or treated
as successful. Interrupted/incomplete attempts retry with the same settings in
new directories, preserving every earlier attempt and execution record.
Configuration mismatch or corrupted completed artifacts is rejected.

Valid returned fits require finite final five-metabolite results and readable
COORD/PRINT outputs without fatal diagnostics or nonzero execution. All actual
PS/TABLE/COORD/PRINT outputs are retained and checked. A numerical/convergence
warning is recorded separately, not automatically equated to a fatal/software
failure. No favorable-fit or outlier filtering is used. Failure-rate denominators
are all20 acquisitions per finite-SNR level. The timeout and interrupted/fatal
cache branches are tested with controlled subprocess substitutes, not claimed
as failures observed in these121 fits.

## Comparison references and statistics

Three concentration representations are retained separately:

1. **LCModel standard concentration**, with its Phase7b clean result as its
   noise reference. It is not interpreted as a source-equivalent coefficient.
2. **FSL source-equivalent coefficient**, with its saved clean fitted coefficient
   as its noise reference. Nominal truth is also reported.
3. **FSL standard water-referenced concentration**, with its own clean result as
   noise reference. Its established quantification bias is not corrected.

For each representation and SNR, summary.csv reports mean,nominal bias,
clean-reference bias,sample SD(ddof1),RMSE against nominal and clean reference,
median absolute errors,median uncertainty and >20% error fractions. The clean
row has no sample SD. Concentrations are never rescaled to force agreement.
The additional SD/clean-baseline percentage is a dimensionless statistic for
comparison across different toolbox scales, not a corrected concentration.

Grouped concentration is **the sum of actual individual GPC and PCh results**.
Its empirical SD is computed from those sums, preserving covariance. LCModel's
separate reported group concentration and group %SD are also preserved; small
rounding differences from the sum are not overwritten. No independent-error
quadrature is used. The saved FSL grouped uncertainty already includes its
GPC/PCh cross-covariance.

## LCModel concentration results

Each row has20/20 valid fits,zero software failure rate. All values except the
last column are mM. Complete median errors,uncertainties and baseline-relative
poor fractions are in summary.csv; no cases are omitted.

| SNR | Metabolite | Mean | Bias nominal | Bias clean | SD | RMSE nominal | RMSE clean | Error >20% nominal |
|---|---|---|---|---|---|---|---|---|
| 100 | NAA | 9.9610 | -0.0391 | -0.6780 | 0.2969 | 0.2920 | 0.7372 | 0% |
| 100 | Cr | 7.9493 | -0.0507 | -0.4967 | 0.2162 | 0.2168 | 0.5396 | 0% |
| 100 | GPC | 1.9836 | -0.0164 | -0.1414 | 0.1111 | 0.1096 | 0.1781 | 0% |
| 100 | PCh | 1.0220 | 0.0220 | -0.0290 | 0.1027 | 0.1025 | 0.1042 | 0% |
| 100 | Glu | 7.8937 | -0.1063 | -0.4823 | 0.3020 | 0.3130 | 0.5650 | 0% |
| 100 | GPC+PCh | 3.0056 | 0.0056 | -0.1704 | 0.0838 | 0.0819 | 0.1890 | 0% |
| 50 | NAA | 9.9831 | -0.0168 | -0.6558 | 0.2747 | 0.2683 | 0.7084 | 0% |
| 50 | Cr | 8.0075 | 0.0075 | -0.4384 | 0.2691 | 0.2624 | 0.5109 | 0% |
| 50 | GPC | 2.0653 | 0.0653 | -0.0597 | 0.2210 | 0.2251 | 0.2236 | 5% |
| 50 | PCh | 0.9951 | -0.0049 | -0.0559 | 0.2163 | 0.2109 | 0.2181 | 35% |
| 50 | Glu | 8.0472 | 0.0472 | -0.3287 | 0.4412 | 0.4326 | 0.5413 | 0% |
| 50 | GPC+PCh | 3.0604 | 0.0604 | -0.1156 | 0.0933 | 0.1092 | 0.1471 | 0% |
| 25 | NAA | 9.8945 | -0.1054 | -0.7444 | 0.4110 | 0.4142 | 0.8454 | 0% |
| 25 | Cr | 8.0186 | 0.0185 | -0.4274 | 0.3228 | 0.3152 | 0.5308 | 0% |
| 25 | GPC | 2.1044 | 0.1044 | -0.0206 | 0.3882 | 0.3925 | 0.3789 | 40% |
| 25 | PCh | 1.0690 | 0.0690 | 0.0180 | 0.4111 | 0.4066 | 0.4011 | 65% |
| 25 | Glu | 8.5449 | 0.5449 | 0.1689 | 0.6757 | 0.8548 | 0.6799 | 5% |
| 25 | GPC+PCh | 3.1734 | 0.1734 | -0.0027 | 0.1653 | 0.2367 | 0.1612 | 0% |
| 15 | NAA | 9.8589 | -0.1411 | -0.7801 | 0.6506 | 0.6496 | 1.0053 | 0% |
| 15 | Cr | 8.1306 | 0.1305 | -0.3154 | 0.4613 | 0.4682 | 0.5492 | 0% |
| 15 | GPC | 1.7129 | -0.2871 | -0.4121 | 0.6471 | 0.6930 | 0.7534 | 75% |
| 15 | PCh | 1.3500 | 0.3500 | 0.2990 | 0.6685 | 0.7397 | 0.7169 | 80% |
| 15 | Glu | 8.7824 | 0.7824 | 0.4064 | 1.5051 | 1.6626 | 1.5223 | 20% |
| 15 | GPC+PCh | 3.0629 | 0.0629 | -0.1131 | 0.1592 | 0.1674 | 0.1920 | 0% |
| 10 | NAA | 10.4026 | 0.4026 | -0.2364 | 1.1003 | 1.1456 | 1.0982 | 10% |
| 10 | Cr | 8.3853 | 0.3853 | -0.0607 | 0.6848 | 0.7706 | 0.6702 | 0% |
| 10 | GPC | 1.8568 | -0.1432 | -0.2682 | 1.0060 | 0.9909 | 1.0165 | 80% |
| 10 | PCh | 1.3826 | 0.3826 | 0.3316 | 1.0566 | 1.0986 | 1.0819 | 90% |
| 10 | Glu | 8.5825 | 0.5825 | 0.2065 | 1.8108 | 1.8586 | 1.7770 | 50% |
| 10 | GPC+PCh | 3.2394 | 0.2394 | 0.0634 | 0.1962 | 0.3064 | 0.2015 | 0% |
| 5 | NAA | 11.0213 | 1.0213 | 0.3824 | 1.8739 | 2.0927 | 1.8661 | 25% |
| 5 | Cr | 9.0616 | 1.0615 | 0.6156 | 2.1810 | 2.3761 | 2.2131 | 45% |
| 5 | GPC | 1.0786 | -0.9214 | -1.0464 | 1.3332 | 1.5930 | 1.6684 | 80% |
| 5 | PCh | 2.4178 | 1.4178 | 1.3668 | 1.3511 | 1.9350 | 1.8980 | 100% |
| 5 | Glu | 9.9808 | 1.9808 | 1.6048 | 5.4567 | 5.6754 | 5.5553 | 75% |
| 5 | GPC+PCh | 3.4964 | 0.4963 | 0.3203 | 0.6496 | 0.8045 | 0.7095 | 50% |

At SNR5, poor-recovery fractions are NAA25%,Cr45%,GPC80%,PCh100%,Glu75%,
GPC+PCh50%. These are scientifically poor estimates despite zero fatal failures.
GPC hits zero in9/20 cases, and PCh hits zero in2/20 cases. Glu ranges from0
to20.595mM. Those boundary/extreme results remain in all concentration statistics.

## Precision compared with FSL-MRS

Cells contain **LCModel standard SD / FSL standard SD / FSL source-equivalent
SD**, in mM. These are different reported scales; absolute SD alone does not
establish a superior estimator. The dimensionless precision plot and paired
clean-reference errors address the different clean baselines explicitly.

| SNR | NAA | Cr | GPC | PCh | Glu | GPC+PCh |
|---|---|---|---|---|---|---|
| 100 | 0.297 / 0.057 / 0.076 | 0.216 / 0.029 / 0.040 | 0.111 / 0.023 / 0.031 | 0.103 / 0.025 / 0.036 | 0.302 / 0.116 / 0.168 | 0.084 / 0.014 / 0.019 |
| 50 | 0.275 / 0.105 / 0.141 | 0.269 / 0.080 / 0.103 | 0.221 / 0.050 / 0.073 | 0.216 / 0.052 / 0.072 | 0.441 / 0.248 / 0.351 | 0.093 / 0.020 / 0.029 |
| 25 | 0.411 / 0.183 / 0.236 | 0.323 / 0.168 / 0.224 | 0.388 / 0.154 / 0.219 | 0.411 / 0.133 / 0.187 | 0.676 / 0.418 / 0.590 | 0.165 / 0.063 / 0.089 |
| 15 | 0.651 / 0.283 / 0.375 | 0.461 / 0.251 / 0.344 | 0.647 / 0.160 / 0.226 | 0.669 / 0.178 / 0.253 | 1.505 / 0.940 / 1.325 | 0.159 / 0.087 / 0.122 |
| 10 | 1.100 / 0.553 / 0.759 | 0.685 / 0.411 / 0.555 | 1.006 / 0.271 / 0.385 | 1.057 / 0.296 / 0.416 | 1.811 / 1.173 / 1.653 | 0.196 / 0.115 / 0.161 |
| 5 | 1.874 / 0.923 / 1.240 | 2.181 / 1.173 / 1.603 | 1.333 / 0.518 / 0.722 | 1.351 / 0.426 / 0.593 | 5.457 / 3.208 / 4.457 | 0.650 / 0.376 / 0.519 |

For example, at SNR100 LCModel NAA SD0.297mM exceeds FSL source-equivalent
SD0.076mM. At SNR5, individual PCh remains much less stable than its group in
both packages. The comparison includes different water-estimation and baseline
models, not just different metabolite optimizers.

## Uncertainty behavior

LCModel %SD is a displayed,integer-rounded uncertainty estimate and is not a
fully calibrated frequentist interval. A displayed0% is **not zero uncertainty**;
999% is a capped/undefined marker. All displayed products `C*%SD/100` are saved,
but zero/capped or zero-concentration cases are marked unavailable for the
positive uncertainty summaries. Their concentrations are still included.
The valid uncertainty count and capped count are explicit below and in CSV.

LCModel clean NAA is10.639mM with displayed0%SD despite+6.39% nominal bias.
Thus low %SD can coexist with systematic error. Median reported uncertainty
usually increases as SNR decreases, but positivity boundaries/capping can make
individual-metabolite medians misleading at the lowest SNR.

| SNR | Metabolite | Median %SD (incl. cap) | Empirical SD mM | Median valid reported SD mM | Valid uncertainties | Capped |
|---|---|---|---|---|---|---|
| 100 | NAA | 1.0 | 0.297 | 0.099 | 20/20 | 0 |
| 100 | GPC | 6.0 | 0.111 | 0.118 | 20/20 | 0 |
| 100 | PCh | 12.0 | 0.103 | 0.118 | 20/20 | 0 |
| 100 | Glu | 3.0 | 0.302 | 0.239 | 20/20 | 0 |
| 100 | GPC+PCh | 1.0 | 0.084 | 0.030 | 20/20 | 0 |
| 50 | NAA | 2.0 | 0.275 | 0.199 | 20/20 | 0 |
| 50 | GPC | 10.0 | 0.221 | 0.212 | 20/20 | 0 |
| 50 | PCh | 23.0 | 0.216 | 0.217 | 20/20 | 0 |
| 50 | Glu | 5.0 | 0.441 | 0.439 | 20/20 | 0 |
| 50 | GPC+PCh | 2.0 | 0.093 | 0.061 | 20/20 | 0 |
| 25 | NAA | 4.0 | 0.411 | 0.399 | 20/20 | 0 |
| 25 | GPC | 20.0 | 0.388 | 0.421 | 20/20 | 0 |
| 25 | PCh | 36.5 | 0.411 | 0.429 | 20/20 | 0 |
| 25 | Glu | 10.0 | 0.676 | 0.819 | 20/20 | 0 |
| 25 | GPC+PCh | 4.0 | 0.165 | 0.128 | 20/20 | 0 |
| 15 | NAA | 6.0 | 0.651 | 0.610 | 20/20 | 0 |
| 15 | GPC | 43.5 | 0.647 | 0.683 | 20/20 | 0 |
| 15 | PCh | 45.5 | 0.669 | 0.684 | 19/20 | 1 |
| 15 | Glu | 14.5 | 1.505 | 1.311 | 20/20 | 0 |
| 15 | GPC+PCh | 7.0 | 0.159 | 0.209 | 20/20 | 0 |
| 10 | NAA | 9.0 | 1.100 | 0.906 | 20/20 | 0 |
| 10 | GPC | 54.5 | 1.006 | 0.943 | 19/20 | 1 |
| 10 | PCh | 69.5 | 1.057 | 1.002 | 16/20 | 4 |
| 10 | Glu | 22.0 | 1.811 | 1.935 | 20/20 | 0 |
| 10 | GPC+PCh | 9.0 | 0.196 | 0.301 | 20/20 | 0 |
| 5 | NAA | 16.0 | 1.874 | 1.734 | 20/20 | 0 |
| 5 | GPC | 223.5 | 1.333 | 1.861 | 11/20 | 9 |
| 5 | PCh | 55.0 | 1.351 | 0.639 | 17/20 | 3 |
| 5 | Glu | 37.0 | 5.457 | 3.959 | 19/20 | 1 |
| 5 | GPC+PCh | 16.5 | 0.650 | 0.588 | 20/20 | 0 |

The grouped LCModel uncertainty uses its **directly reported group %SD**;
no group covariance is invented. FSL source SD is its saved coefficient
covariance result. Its standard-concentration uncertainty is a **conditional
coefficient-only conversion** by the per-fit standard/source ratio; it omits
water-scaling uncertainty. This conversion is applied to uncertainty only,
not to concentration estimates or simulated amplitudes. No LCModel %SD-derived
95% confidence-coverage claim is made. Positive-only uncertainty medians can
be selection-biased when many estimates are capped, so counts must accompany them.

## Paired comparisons

All20 paired acquisitions are used at every SNR for every metabolite. Let
`dL=C_L-C_L,clean` and `dF=C_F-C_F,clean`. paired_comparison.csv includes means
of `abs(dL)-abs(dF)`, `dL²-dF²`, `dL-dF`, and analogous nominal-truth errors,
for **both** FSL standard and source-equivalent comparators. Fractional
clean-baseline absolute/squared error differences are also saved as dimensionless
statistics; original concentrations remain unchanged.

Intervals are paired percentile-bootstrap95% intervals for **mean differences**,
10000 draws,resampling the paired acquisition indices together. Local PCG64
seeds are `77000+100*level_index+metabolite_index`. No global RNG state changes,
independent unpaired resampling,p-value or statistical-significance assertion
is used. With20 realizations,small effects and multiple comparisons require
caution; these are conditional summaries of one fixed synthetic experiment.

NAA baseline-centered absolute-error differences versus FSL standard (mM):

| SNR | Mean paired MAE difference | 95% paired bootstrap interval | Fractional MAE difference, percentage points |
|---|---|---|---|
| 100 | 0.6312 | [0.4928, 0.7641] | 5.710 |
| 50 | 0.5751 | [0.4414, 0.7027] | 5.020 |
| 25 | 0.5941 | [0.4082, 0.7930] | 4.868 |
| 15 | 0.6464 | [0.4162, 0.8808] | 5.111 |
| 10 | 0.4090 | [0.1673, 0.6632] | 1.541 |
| 5 | 0.7609 | [0.4569, 1.0780] | 3.686 |

Positive values mean larger deviation from LCModel's own clean result. They do
not mean LCModel has larger nominal-truth error: FSL standard concentrations
have a different systematic water-scaling bias. The LCModel high-SNR shift away
from its clean result contributes to these paired errors. All72 paired
comparisons and intervals are in the CSV, not only this illustrative NAA subset.

## Water scaling and residuals

Even at water target SNR1000, LCModel's directly integrated water area varies
more than FSL's fitted water initial amplitude. They are **different estimands**;
only their relative variation is compared here. No fitted amplitude is equated
to an integrated area. LCModel's FCALIB contains the inverse water area; it is
expected to vary with the noisy reference under the frozen settings.

| SNR | LCModel water-area CV % | FSL water-amplitude CV % | NAA–FCALIB correlation |
|---|---|---|---|
| 100 | 2.286 | 0.143 | 0.919 |
| 50 | 2.205 | 0.115 | 0.819 |
| 25 | 2.471 | 0.078 | 0.691 |
| 15 | 2.061 | 0.075 | 0.453 |
| 10 | 2.264 | 0.095 | 0.302 |
| 5 | 2.465 | 0.092 | 0.129 |

At SNR100 the correlation is0.919, consistent with a substantial shared
water-normalization contribution to concentration variability. It does not
uniquely partition optimizer and water noise. Mean water scaling stays near
the clean value, while NAA's mean decreases by0.678mM; the systematic mean
change cannot be explained by mean FCALIB alone. Adaptive fitting/baseline
behavior remains relevant. LCModel %SD should not be assumed to propagate all
of this independent reference variability.

LCModel residuals are phased-real spectra over its analysis window; saved
FSL residuals are full complex FIDs. They are recorded and plotted with separate
labels/domains, **not interpreted as interchangeable fit-quality scores**.
Increasing noise raises both residuals. The clean LCModel1.816% versus FSL
6.04e-7 residual discrepancy remains documented; no basis-duration alteration
was made to address it in this benchmark.

## Warning and poor-recovery cases

| Diagnostic | Acquisitions affected |
|---|---|
| warning:MYBASI8 | 121 |
| info:STARTV20 | 121 |
| info:FINOUT7 | 1 |
| info:FINOUT9 | 103 |
| info:RFALSI4 | 97 |
| info:RFALSI11 | 32 |
| info:TWOREG8 | 8 |
| info:RFALSI13 | 1 |
| info:RFALSI5 | 1 |
| info:TWOREG7 | 1 |
| info:RFALSI12 | 1 |

MYBASI8 is the preserved default preliminary-list mismatch (Ins is absent).
STARTV20 is the narrow-line minimum convolution-support diagnostic, not a
minimum physical linewidth or fatal error. FINOUT9 reports an unsmooth line
shape; FINOUT7 appears only in the clean case. RFALSI informational diagnostics
are retained without reclassifying every message as optimizer failure.

SNR10 realization17 has RFALSI13,RFALSI5 and TWOREG7: a numerical/refinement
convergence issue, despite a finite final result. It remains included and is
flagged. LCModel provides no universal final optimizer-success boolean here.
Zero fatal/software failures must not be read as proof that every refinement
converged or every concentration is reliable.

The additional poor-recovery PDF is SNR5 realization15, selected by the largest
individual nominal relative concentration error. It reports PCh4.565mM
(+356.5%),GPC0mM and Glu16.863mM. Its uncertainty/diagnostics are preserved;
it is not excluded from summaries. Original FSL results for the same case
remain available, without creating a new FSL report.

## Reports and plots

All121 PS/PRINT/TABLE/COORD outputs are retained. Six representative PDFs were
converted using the installed Ghostscript workflow. Original FSL HTML reports
were read only, never regenerated. PDF/HTML index:

- SNRinfinite,realization00: [LCModel PDF](../outputs/lcmodel_noise_comparison/SNR_infinite/realization_00/attempt_001/lcmodel_outputs/fit.pdf); [existing FSL HTML](../outputs/mrscloud_noise_fitting_validation/full/SNR_infinite/realization_00/fsl_mrs_results/report.html).
- SNR50,realization00: [LCModel PDF](../outputs/lcmodel_noise_comparison/SNR_50/realization_00/attempt_001/lcmodel_outputs/fit.pdf); [existing FSL HTML](../outputs/mrscloud_noise_fitting_validation/full/SNR_50/realization_00/fsl_mrs_results/report.html).
- SNR25,realization00: [LCModel PDF](../outputs/lcmodel_noise_comparison/SNR_25/realization_00/attempt_001/lcmodel_outputs/fit.pdf); [existing FSL HTML](../outputs/mrscloud_noise_fitting_validation/full/SNR_25/realization_00/fsl_mrs_results/report.html).
- SNR10,realization00: [LCModel PDF](../outputs/lcmodel_noise_comparison/SNR_10/realization_00/attempt_001/lcmodel_outputs/fit.pdf); [existing FSL HTML](../outputs/mrscloud_noise_fitting_validation/full/SNR_10/realization_00/fsl_mrs_results/report.html).
- SNR5,realization00: [LCModel PDF](../outputs/lcmodel_noise_comparison/SNR_5/realization_00/attempt_001/lcmodel_outputs/fit.pdf); [existing FSL HTML](../outputs/mrscloud_noise_fitting_validation/full/SNR_5/realization_00/fsl_mrs_results/report.html).
- SNR5,realization15: [LCModel PDF](../outputs/lcmodel_noise_comparison/SNR_5/realization_15/attempt_001/lcmodel_outputs/fit.pdf); no Phase7a HTML for this realization; original numerical FSL results retained.

![Bias relative to nominal concentration](../outputs/lcmodel_noise_comparison/nominal_bias_mM.png)

![Bias relative to each toolbox clean result](../outputs/lcmodel_noise_comparison/noise_bias_mM.png)

![Empirical precision, all six metabolite/group estimates](../outputs/lcmodel_noise_comparison/SD_mM.png)

![Dimensionless precision relative to clean toolbox results](../outputs/lcmodel_noise_comparison/SD_percent_clean_baseline.png)

![Nominal concentration RMSE](../outputs/lcmodel_noise_comparison/nominal_RMSE_mM.png)

![Noise-reference RMSE](../outputs/lcmodel_noise_comparison/noise_RMSE_mM.png)

![Poor recovery fraction](../outputs/lcmodel_noise_comparison/nominal_poor_fraction.png)

![Reported uncertainty and empirical variability](../outputs/lcmodel_noise_comparison/uncertainty_comparison.png)

![Paired error intervals](../outputs/lcmodel_noise_comparison/paired_error_comparison.png)

![Paired per-acquisition deviations](../outputs/lcmodel_noise_comparison/paired_deviations.png)

![Different residual domains](../outputs/lcmodel_noise_comparison/residuals.png)

## Reproduction, output files and tests

Executed from repository root:

```sh
pixi run python --version
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.lcmodel_noise_comparison --smoke --output outputs/lcmodel_noise_comparison_smoke
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m examples.lcmodel_noise_comparison
# Repeated full command verifies121 cache hits; no extra fitting.
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q tests/test_lcmodel_noise_analysis.py tests/integration/test_lcmodel_noise_comparison.py
MPLCONFIGDIR=/tmp/mrs-phase6c-matplotlib pixi run python -m pytest -q
```

Each actual command is recorded in commands.txt and the per-attempt execution
JSON, including timestamps,timeout,elapsed time,stdout/stderr,exit status and
CONTROL. A representative invocation is:

```sh
cd outputs/lcmodel_noise_comparison/SNR_100/realization_00/attempt_001
$LCMODEL_BINARY < lcmodel.CONTROL
```

This executable was invoked by Python launched with Pixi, using the exact local
binary; there was no download or alternate installation. Frozen configuration
also records Python/NumPy/FSL versions and analysis source hashes.

Root outputs include frozen_configuration.json,frozen_protocol.CONTROL,
per_realization.csv,summary.csv,paired_comparison.csv,water_scaling_comparison.csv,
results.json,reports.json,failure_summary.json,warning_summary.json,commands.txt
and diagnostic plots. Per-acquisition directories contain exact RAW/H2O,
conversion checks,identity/signature,atomic completion and separate attempts.
Raw printed concentration tokens and report values are retained at their actual
serialization precision; unavailable internal unrounded coefficients are not
invented. This implementation does not derive LCModel source-equivalent coefficients.

Focused tests: **16 passed**. Tests cover input/CONTROL identity,missing inventory,
exact noisy RAW round-trip,fatal detection,grouped uncertainty,hand-calculated
statistics,failure denominators,paired bootstrap and RNG isolation,timeout retry,
cached fatal preservation,corruption rejection,all121 actual outputs and source
hashes,report existence and cache reuse without invoking the executable.
Proprietary-execution artifacts are checked in separate integration tests and
skip if this local experiment is absent; no proprietary fit is run by pytest.

The first full-suite run was **313passed,1failed,1xfailed**: an older Phase7b.1
cache test redrew two old PNGs under another test's global Matplotlib style.
The scientific input and numerical-result hashes remained unchanged. The PNGs
were restored to their **original recorded hashes**, and that cache test now
suppresses savefig while validating cached numerical results, preserving earlier
reports. No provenance check was weakened and no LCModel fits were repeated.
The corrected full suite completed with **314passed,1xfailed**,no failures or
skips,in **203.22seconds**. The retained expected failure is the installed
FSL-MRS combined water-phase/frequency initialization defect. Four warnings
are external Dash/SWIG deprecations. Final source checks confirm all1118
protected files match their frozen hashes. The uncertainty/bias plots and the
SNR5 realization15 rendered report were inspected; no fitting changes followed.

Files created: utils/lcmodel_noise_analysis.py,examples/lcmodel_noise_comparison.py,
tests/test_lcmodel_noise_analysis.py,tests/integration/test_lcmodel_noise_comparison.py,
this document and outputs/lcmodel_noise_comparison/. The existing Phase7b.1 cache
test was adjusted only to preserve its reports. Simulator,conversion conventions,
calibration,compatible BASIS and fitting parameters were unchanged.

## Limitations and review boundary

There are only20 realizations per SNR and one homogeneous WM acquisition,
with ideal Gaussian receiver noise. Bootstrap intervals are conditional and
exploratory,not general proof of population performance or interval coverage.
LCModel and FSL use different basis formatting,water estimation,baseline and
concentration conventions; the experiment compares the frozen complete
pipelines,not isolated optimizers with identical model families. The clean
LCModel fit has known representation/regularization bias,which may respond
nonlinearly to added noise. Standard concentrations and source coefficients
must remain distinct. Reported %SD may omit independent water-normalization
uncertainty,and boundary/capped values limit its interpretation. No signals,
noise amplitudes,water gains or fitting parameters were changed to improve
agreement. Stop at Phase7c; no new simulator physics or calibration correction
is proposed by these results.

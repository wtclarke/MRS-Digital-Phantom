# Proposed synthetic quantitative signal model

**Status:** an isolated analytical prototype now implements the initial response convention in `simulation/analytical_water_reference.py`. The integrated water-reference mode remains a design only. This prototype does not establish absolute quantification.

## Decisions and scope

- The existing simulator remains unchanged. The proposed path is explicitly opt-in and must not modify legacy output when disabled.
- Ground-truth concentrations are nominal **mM per total tissue volume**. These are simulation labels for the requested model, not guaranteed intracellular or tissue-water concentrations. Their provenance combines heterogeneous literature/database values and conversions; it does not establish a traceable absolute standard.
- The checked-in MRSCloud basis is not treated as absolutely calibrated. Its complex response is useful as a fitting template, but its headers do not establish signal per molar concentration or per proton.
- Metabolite proton response may already be represented in spin-system basis amplitudes. No additional proton multiplier is applied unless a separate provenance study verifies that it is required and is not double-counted.
- Water and metabolite branches use a single declared **synthetic response convention** (below). This supports internally consistent simulation and controlled recovery experiments only. It is not physical calibration without independently established basis response per concentration and a common acquisition/receiver scale.
- Correct the conventional T2 FID decay only inside the new mode. Preserve the legacy equation and outputs.
- Adopt one complex NIfTI-MRS convention only after the proposed complex round-trip tests identify and validate it.

## Units and signal equations

Let (t_n=n\Delta t) seconds, (n=0,\ldots,N-1), where \(\Delta t=1/BW\) seconds and (BW\) is the complex sampling bandwidth in Hz. Signals below are complex FID samples in **synthetic response units (SRU)**. Concentrations are in mol/L of total tissue volume; input nominal mM values are converted as \(C[\mathrm{mol/L}]=10^{-3}C[\mathrm{mM}]\). Tissue fractions (f_j\) are dimensionless and sum to one over included tissue.

Represent the formatted metabolite basis for metabolite (m), tissue/sequence condition (j), and sampling definition \(\mathcal{P}=(BW,N,\text{sequence},TE,\text{selected basis set})\) as a complex vector \(b_{mj}[n;\mathcal{P}]\). Its units are deliberately **SRU per nominal mol/L under this exact formatting convention**, not a claimed physical unit. The basis is sequence-evolved: its frequency pattern and coupling evolution at its specified TE are retained. It is not assumed to encode calibrated concentration or a separately verified T2 attenuation.

For an acquisition with repetition time (TR\) seconds and echo time (TE\) seconds, define the separately modeled relaxation factor under a fully spoiled, single-pulse, steady-state approximation:

\[
R_{mj}(TR,TE)=\left(1-e^{-TR/T_{1,mj}}\right)e^{-TE/T_{2,mj}}.
\]

This factor is dimensionless. It is an optional model in the new mode; its sequence assumptions and parameters must be stored. Do not apply the TE factor if the supplied response is later demonstrated to already include that relaxation attenuation. Sequence-dependent evolution (including J evolution) remains part of (b\), distinct from this scalar TE attenuation.

Acquisition-time transverse decay is a separate FID envelope. For physical (T_2\) in seconds, use

\[
L_{mj}(t_n)=e^{-t_n/T_{2,mj}}=e^{-\pi\Delta f_{mj}t_n},\qquad
\Delta f_{mj}=1/(\pi T_{2,mj})\ \mathrm{Hz}.
\]

An independently specified additional linewidth \(\Delta f_{\mathrm{extra}}\) Hz may contribute \(e^{-\pi\Delta f_{\mathrm{extra}}t_n}\); it is not a substitute for TE relaxation or sequence evolution. The erroneous legacy expression \(e^{-t/(\pi T_2)}\) is not reused in the new path.

The proposed synthetic metabolite signal is

\[
s_m[n]=K_{\mathrm{syn}}\sum_j f_j C_{mj}\,R_{mj}(TR,TE)
 b_{mj}[n;\mathcal{P}]\,L_{mj}(t_n)\,e^{i(2\pi\Delta f_{0,mj}t_n+\phi_{mj})}.
\]

Here \(C_{mj}\) is mol/L, (\Delta f_{0,mj}\) is an offset in Hz, \(\phi\) is radians, and \(K_{\mathrm{syn}}\) is the shared synthetic scale defined below. If the basis response itself is normalized at runtime, that normalization must be captured in its response metadata and (b\) must be interpreted using that exact factor. Tissue fractions and relaxation are applied tissue by tissue before summing; a single averaged relaxation factor is generally not equivalent.

For water, let (C_{w,j}\) be tissue water concentration in mol/L (the repository's approximate WM/GM/CSF values can be input as such, with provenance recorded). Define a unit-concentration water response (b_{w,j}[n;\mathcal{P}]\) in the same SRU per mol/L convention. Then

\[
s_w[n]=K_{\mathrm{syn}}\sum_j f_j C_{w,j}\,R_{w,j}(TR,TE)
 b_{w,j}[n;\mathcal{P}]\,e^{-\pi\Delta f_{w,j}t_n}
 e^{i(2\pi\Delta f_{0,w,j}t_n+\phi_{w,j})}.
\]

Water proton response is part of (b_w\) by definition if the water response is derived per water molecule; do not also multiply by two. Likewise, no explicit metabolite proton multiplier is currently proposed. The model's unit-concentration water and metabolite responses must share (K_{syn}\), but available provenance cannot derive their ratio from first principles. Until independent calibration exists, both response functions are declared synthetic and are only mutually comparable under documented assumptions.

For the simplest synthetic water response, use a unit initial complex response per mol/L with a stated phase and line shape, so its amplitude is linear in water concentration. This defines a reproducible convention, not an experimentally calibrated water sensitivity. No value is assigned to (K_{syn}\) to force a desired fitting result. It may be fixed to a unitless identity convention (one SRU per chosen response unit), provided the same convention is applied to both arms and the output is labeled synthetic.

Residual water remains its current independent contamination model and settings. It is not the unsuppressed water reference and must not influence (s_w\).

The FFT convention must be explicit. Prefer preserving the complex FID as the canonical internal signal; any spectral conversion uses the existing unnormalized transform pair and must not rescale amplitudes. Export metadata must state sampling and phase convention.

### Implemented analytical prototype

The isolated prototype defines response per mol/L as the number of equivalent contributing protons: (A=K_{syn} C n_H) SRU. The single metabolite resonance uses (n_H=3\); water uses (n_H=2\). For nominal metabolite concentration (10\,\mathrm{mM}=0.010\,\mathrm{mol/L}) and water concentration (40\,\mathrm{mol/L}), the initial amplitudes are (0.030K_{syn}\) and (80K_{syn}\) SRU, respectively, giving (0.030/80=0.000375\). A shared arbitrary receiver factor (K_{syn}\) cancels exactly in this ratio. Both FIDs are constant across samples in this zero-offset, zero-phase, no-relaxation prototype. This proton-count convention is an explicit synthetic definition for the prototype; it is not evidence that MRSCloud basis amplitudes use the same convention. The prototype tests validate internal synthetic consistency only, not physical scanner calibration.

## Basis normalization and provenance

`Basis.reformat()` currently divides the formatted basis array by one shared mean absolute sample magnitude. This is not per-metabolite normalization, but it removes the input's absolute scale. The factor can depend on selected metabolites and on bandwidth/point sampling because those choices change the formatted array over which the mean is computed. Therefore the quantitative mode must record at least: source basis identifier/hash, selected metabolite names, original and formatted sampling, normalization statistic and numeric value, sequence, TE, and any interpolation/truncation applied. Do not silently reuse a normalization factor from another basis selection or sampling grid.

The legacy normalizer remains untouched. The new mode may either use the same normalized response while explicitly declaring its synthetic meaning, or preserve source amplitudes where possible; neither option turns MRSCloud amplitudes into absolute calibration. A future physical calibration requires traceable basis response per amount and a water response on the same scale.

## Avoiding circular fitting validation

Do not validate only by simulating with the same basis and recovering concentrations with a fitter that uses that exact basis normalization and the same assumptions. Such a round trip can pass because shared errors cancel. Use two layers:

1. **Independent analytical validation:** construct signals from simple analytic single-resonance or orthogonal-component responses with declared unit response and known concentration. Check predicted complex FID samples, tissue mixing, attenuation, and export numerically without invoking the production fitting pipeline.
2. **Cross-model recovery:** generate a held-out synthetic dataset from independently specified analytic responses (or a separately generated, provenance-documented basis); quantify it with the target fitting pipeline and compare to the known nominal concentrations. Separately perturb basis selection/sampling, tissue fractions, relaxation, and noise to expose sensitivity. If the target fitter cannot consume the independent response, report that limitation and restrict the claim to internal consistency.

Never tune (K_{syn}\) against the same fitter's returned concentrations and then cite that recovery as independent validation. A physical-calibration claim additionally requires an external calibrated reference or independently justified response ratio.

## Proposed interfaces and files

Names are provisional and should follow existing project conventions. Keep this as a separate path so default simulation behavior is unaffected.

- `simulation/signalModel.py`: add an opt-in branch or separate method such as `simulate_quantitative_pair(...)` that receives explicit mode/configuration and returns metabolite and water-reference complex FIDs plus metadata. Avoid routing the legacy path through changed normalization or relaxation logic.
- `simulation/basis.py`: expose/record the formatted basis normalization provenance needed by the new path (selected entries, grid, shared scalar, source identity). No proton multiplier by default.
- `simulation/water/water_reference.py` (new): implement only the independent unsuppressed water response and documented attenuation/line-shape terms. Do not call the residual-water generator.
- `utils/nii_processing.py` and `simulation/headless.py`: converge export on one complex FID-to-NIfTI-MRS convention after round-trip evidence; share one conversion route if feasible and retain the complex phase exactly.
- `utils/definitions.py` / concentration loading: document water units and nominal metabolite concentration semantics at the boundary; do not alter legacy concentration values as part of this design.
- `tests/`: add independent analytical signal tests and complex export/reload round-trip tests. Proposed tests are not yet implemented.
- Configuration/export metadata: add an opt-in flag and record concentration convention, synthetic calibration label, basis normalization provenance, tissue fractions, sequence/TE/TR, T1/T2, linewidth, residual-water independence, seed, and export convention. GUI changes should follow validation, not precede it.

## Minimal independent analytical tests

1. **Concentration linearity:** for fixed analytic unit response, doubling one metabolite or tissue water concentration doubles each complex FID sample; zero concentration gives zero contribution.
2. **Tissue weighting:** two-tissue signal equals the explicit sum \(\sum_j f_j C_j R_j b_j\); pure tissue limits and fraction scaling match closed-form expectations.
3. **Relaxation separation:** verify independently that sequence/TE response is not duplicated by (e^{-TE/T2}\), that (R(TR,TE)\) matches the declared expression, and that acquisition-time decay equals \(e^{-t/T2}\). Verify legacy output remains byte/value-identical for a fixed seeded case once implemented.
4. **Residual-water independence:** changing legacy residual-water settings leaves the separate reference FID unchanged for the same seed/configuration.
5. **Basis-normalization provenance:** selection or sampling changes that alter the shared normalization scalar are reflected in metadata; a fixed recorded basis response is deterministic.
6. **Determinism:** fixed configuration and seed produce identical complex outputs without changing the random stream/output of legacy simulation.
7. **Complex export round trip:** use a multi-sample FID with nonzero real and imaginary parts and verify export/reload preserves complex samples and reconstructed spectrum within complex64 tolerance, with the expected phase convention. Cover each supported shape.

These tests establish model arithmetic and serialization, not absolute quantification. End-to-end quantification remains a separate cross-model validation gate.

## Staged implementation plan

1. **Design review:** review this model and resolve sequence/relaxation and export convention questions below. No implementation yet.
2. **Analytical core:** implement a minimal isolated synthetic pair generator with explicit inputs and deterministic response functions; preserve legacy code path exactly.
3. **Independent tests:** complete tests 1–6 against closed-form signals; demonstrate legacy parity and residual-water independence.
4. **Export convention:** select one export path from complex round-trip evidence; add paired export/reload tests and metadata.
5. **Voxel integration:** apply tissue fractions, voxel volume, and shared spatial effects consistently; verify analytic tissue mixtures.
6. **Cross-model validation:** quantify held-out independent synthetic data; report recovery and sensitivity without claiming physical calibration.
7. **Optional user configuration:** add configuration and later GUI controls only after numerical validation; keep disabled default.
8. **Physical calibration gate:** do not label output absolutely calibrated until independent metabolite and water response provenance, concentration semantics, and external/analytical validation support that claim.

## Remaining scientific questions

- Do MRSCloud basis spin systems and export retain a known signal per molecule or per proton response? The checked-in metadata do not establish this; additional source-generation provenance or a calibrated reference is needed.
- Do basis responses include any relaxation beyond TE-dependent sequence evolution? This determines whether scalar TE/T2 attenuation can be applied without duplication.
- What exact sequence and steady-state assumptions justify (1-e^{-TR/T1}\) for the intended acquisitions, including flip angle and spoiling?
- Are nominal database mM values consistently per total tissue volume across source studies? For now they are the declared synthetic ground-truth convention, with known provenance limitations.
- What water line shape, proton response convention, and tissue-specific relaxation assumptions should the synthetic response use?
- Which existing NIfTI-MRS export path has the intended conjugation/phase convention? Resolve by implementing and running complex round-trip tests before selecting the canonical path.
- Should synthetic outputs use the runtime-normalized basis or preserve raw file amplitude where available? Either choice remains synthetic; record the choice and normalization factor.
- Which independent fitting pipeline/data format will be used for cross-model recovery, and can it accept an independently generated analytic basis?

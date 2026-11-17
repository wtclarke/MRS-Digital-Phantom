# Codex Instructions — MRS Digital Brain Phantom

## Project objective

Preserve legacy MRS simulation and the headless notebook workflow, including its approximate BASIS-H2O water reference.

For changes to the water reference, read `docs/legacy_basis_water_reference.md` before starting work.

## Scientific correctness

- Prioritize physical correctness over plausible-looking spectra.
- Do not introduce arbitrary amplitude calibration factors.
- Do not conflate residual water with an unsuppressed water reference.
- Verify metabolite basis normalization and proton-count conventions.
- Distinguish ground-truth concentrations from fitted estimates.
- Explicitly document concentration units, signal normalization, and relaxation assumptions.
- Preserve existing simulation behavior unless a change is explicitly justified.
- Do not silently change existing concentration or amplitude conventions.

## Development workflow

- Work on one implementation phase at a time.
- Begin by inspecting the existing code.
- Make small, reviewable changes.
- Add tests for changes to concentration, signal amplitude, and water-reference behavior.
- Run relevant tests after making changes.
- Report any scientific ambiguities requiring a decision.
- Do not modify unrelated code.
- Do not create commits unless explicitly requested.

## Validation requirements

- Verify water amplitude linearity with concentration.
- Verify tissue-weighted signal formation.
- Verify independence from residual-water settings.
- Verify compatibility of metabolite and water amplitude conventions.
- Verify deterministic simulation.
- Verify export and reload preserve complex signal amplitudes.
- Do not claim successful absolute quantification without an end-to-end numerical validation.

## Reporting

After each task, report:

1. Files modified.
2. Scientific assumptions.
3. Tests executed and their results.
4. Remaining uncertainties.
5. Recommended next step.

## Scientific scope

The BASIS-H2O reference assumes a common basis scale and nominal mM table values. Preserve these explicit assumptions; introducing a calibrated quantitative pathway requires a separate reviewed normalization audit and implementation plan.

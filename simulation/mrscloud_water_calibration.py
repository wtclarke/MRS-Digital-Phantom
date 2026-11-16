"""Synthetic water response anchored to frozen MRSCloud Ref0ppm, without fitting."""

import numpy as np

from simulation.frozen_basis_normalization import format_frozen_basis


def calibrated_responses(path, *, bandwidth_hz, points, selected_names,
                         calibration=None, expected_metadata=None):
    """Return unit-mol/L metabolite responses, water response, and provenance.

    Frozen metabolite FIDs are unchanged. Water is an analytical zero-offset,
    zero-phase singlet with two times the frozen Ref0ppm initial magnitude.
    This is a synthetic definition, not a physical per-proton calibration.
    """
    if not np.isfinite(bandwidth_hz) or bandwidth_hz <= 0 or points <= 0:
        raise ValueError("sampling must be positive and finite")
    selected_names = list(selected_names)
    if len(set(selected_names)) != len(selected_names) or not selected_names:
        raise ValueError("select distinct metabolite names")
    requested = list(dict.fromkeys([*selected_names, "Ref0ppm"]))
    formatted, names, frozen = format_frozen_basis(
        path, bandwidth_hz=bandwidth_hz, points=points,
        selected_names=requested, calibration=calibration,
    )
    proton_response = float(abs(formatted[0, names.index("Ref0ppm")]))
    if not np.isfinite(proton_response) or proton_response <= 0:
        raise ValueError("reference singlet must have a nonzero finite initial amplitude")
    metadata = {
        "convention": "synthetic-frozen-Ref0ppm-initial-magnitude-v1",
        "physical_calibration": False,
        "frozen_basis": frozen.to_dict(),
        "reference_metabolite": "Ref0ppm",
        "reference_protons_assumed": 1,
        "reference_initial_real": float(formatted[0, names.index("Ref0ppm")].real),
        "reference_initial_imag": float(formatted[0, names.index("Ref0ppm")].imag),
        "proton_response_sru_per_mol_l": proton_response,
        "water_protons": 2,
        "water_frequency_offset_hz": 0.0,
        "water_phase_rad": 0.0,
        "concentration_units": "nominal mol/L per total tissue volume",
        "relaxation": False,
        "metabolite_proton_multiplier": "none; encoded source response retained",
    }
    if expected_metadata is not None and expected_metadata != metadata:
        raise ValueError("synthetic calibration metadata mismatch")
    responses = {name: formatted[:, names.index(name)].copy() for name in selected_names}
    water = np.full(points, 2 * proton_response, dtype=np.complex128)
    return responses, water, metadata


from simulation.quantitative_signals import concentration_signals

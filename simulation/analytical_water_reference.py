"""Isolated synthetic single-resonance signals for analytical validation.

Signals are complex FID samples in arbitrary synthetic response units (SRU).
Concentrations are mol/L; convenience inputs in mM are converted explicitly.
This module does not depend on the legacy simulator or its residual-water model.
"""

import numpy as np


def concentration_mol_per_litre(value: float, unit: str) -> float:
    """Convert a concentration in ``mol/L`` or ``mM`` to ``mol/L``."""
    if unit == "mol/L":
        return float(value)
    if unit == "mM":
        return float(value) * 1e-3
    raise ValueError("unit must be 'mol/L' or 'mM'")


def single_resonance_fid(
    concentration: float,
    *,
    unit: str,
    equivalent_protons: int,
    n_points: int,
    dwell_time_s: float,
    receiver_factor: float = 1.0,
) -> np.ndarray:
    """Return an unattenuated, zero-offset, zero-phase complex FID.

    The initial response is ``receiver_factor * concentration[mol/L] *
    equivalent_protons`` SRU. No relaxation or line broadening is applied.
    """
    if equivalent_protons <= 0 or n_points <= 0 or dwell_time_s <= 0:
        raise ValueError("proton count, n_points and dwell_time_s must be positive")
    concentration_mol_l = concentration_mol_per_litre(concentration, unit)
    initial_amplitude = receiver_factor * concentration_mol_l * equivalent_protons
    return np.full(n_points, initial_amplitude, dtype=np.complex128)


def synthetic_pair(
    *,
    metabolite_concentration_mM: float = 10.0,
    water_concentration_mol_l: float = 40.0,
    n_points: int = 8,
    dwell_time_s: float = 1e-3,
    receiver_factor: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (metabolite, water) FIDs under the shared synthetic convention."""
    metabolite = single_resonance_fid(
        metabolite_concentration_mM,
        unit="mM",
        equivalent_protons=3,
        n_points=n_points,
        dwell_time_s=dwell_time_s,
        receiver_factor=receiver_factor,
    )
    water = single_resonance_fid(
        water_concentration_mol_l,
        unit="mol/L",
        equivalent_protons=2,
        n_points=n_points,
        dwell_time_s=dwell_time_s,
        receiver_factor=receiver_factor,
    )
    return metabolite, water


def recover_metabolite_concentration_mM(
    metabolite_amplitude: complex,
    water_amplitude: complex,
    *,
    water_concentration_mol_l: float,
    metabolite_equivalent_protons: int = 3,
    water_equivalent_protons: int = 2,
) -> float:
    """Recover nominal metabolite concentration in mM from initial amplitudes."""
    if water_amplitude == 0:
        raise ValueError("water amplitude must be nonzero")
    concentration_mol_l = (
        metabolite_amplitude / water_amplitude
        * water_concentration_mol_l
        * water_equivalent_protons / metabolite_equivalent_protons
    )
    return float(np.real(concentration_mol_l) * 1e3)

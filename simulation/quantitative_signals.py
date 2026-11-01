"""Unit-response concentration scaling shared by production and research."""
import numpy as np


def concentration_signals(responses, water_response, *, concentrations_mol_l,
                          water_concentration_mol_l, receiver_factor=1.0):
    """Scale unit responses by explicit concentrations and one shared receiver."""
    values = [*concentrations_mol_l.values(), water_concentration_mol_l, receiver_factor]
    if not all(np.isfinite(x) and x >= 0 for x in values):
        raise ValueError("concentrations and receiver factor must be finite and nonnegative")
    if set(responses) != set(concentrations_mol_l):
        raise ValueError("concentrations must match selected metabolites")
    return ({name: response * concentrations_mol_l[name] * receiver_factor
             for name, response in responses.items()},
            water_response * water_concentration_mol_l * receiver_factor)

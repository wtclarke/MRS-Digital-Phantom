import numpy as np
import pytest

from simulation.analytical_water_reference import (
    recover_metabolite_concentration_mM,
    single_resonance_fid,
    synthetic_pair,
)


def test_expected_complex_fid_samples_and_amplitude_ratio():
    metabolite, water = synthetic_pair(
        n_points=4, dwell_time_s=2e-3, receiver_factor=7.25
    )

    np.testing.assert_allclose(metabolite, np.full(4, 0.2175 + 0j), rtol=1e-15, atol=0)
    np.testing.assert_array_equal(water, np.full(4, 580.0 + 0j))
    assert metabolite[0] / water[0] == pytest.approx(0.000375, abs=1e-15)


def test_fid_is_deterministic_and_initial_amplitude_linear_in_concentration():
    one = single_resonance_fid(
        10, unit="mM", equivalent_protons=3, n_points=3, dwell_time_s=1e-3
    )
    doubled = single_resonance_fid(
        20, unit="mM", equivalent_protons=3, n_points=3, dwell_time_s=1e-3
    )
    repeated, _ = synthetic_pair(n_points=3, dwell_time_s=1e-3)

    np.testing.assert_array_equal(doubled, 2 * one)
    np.testing.assert_array_equal(repeated, one)


def test_water_reference_amplitude_is_linear_in_molar_concentration():
    water_40 = single_resonance_fid(
        40, unit="mol/L", equivalent_protons=2, n_points=2, dwell_time_s=1e-3
    )
    water_20 = single_resonance_fid(
        20, unit="mol/L", equivalent_protons=2, n_points=2, dwell_time_s=1e-3
    )
    np.testing.assert_array_equal(water_40, 2 * water_20)


def test_water_reference_ratio_recovers_nominal_concentration():
    metabolite, water = synthetic_pair(
        metabolite_concentration_mM=10,
        water_concentration_mol_l=40,
        receiver_factor=3.7,
    )
    recovered = recover_metabolite_concentration_mM(
        metabolite[0], water[0], water_concentration_mol_l=40
    )
    assert recovered == pytest.approx(10.0, abs=1e-12)


@pytest.mark.parametrize("legacy_water_amp_factor", [0.0, 5.0, 1234.5])
def test_reference_is_independent_of_legacy_residual_water_setting(
    legacy_water_amp_factor,
):
    # The legacy setting is intentionally not an input to this isolated reference.
    del legacy_water_amp_factor
    metabolite, water = synthetic_pair(n_points=5, dwell_time_s=1e-3)
    np.testing.assert_array_equal(water, np.full(5, 80.0 + 0j))
    np.testing.assert_array_equal(metabolite, np.full(5, 0.03 + 0j))


def test_concentration_units_are_explicit():
    with pytest.raises(ValueError, match="unit must be"):
        single_resonance_fid(
            10, unit="molar", equivalent_protons=3, n_points=1, dwell_time_s=1e-3
        )

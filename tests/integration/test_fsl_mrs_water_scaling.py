"""Optimizer-based integration checks, separate from the fast signal tests."""

import json

import numpy as np
import pytest

from examples.fsl_mrs_water_scaling_validation import run_validation
from simulation.mrscloud_quantitative_pair import generate_pair


@pytest.fixture(scope='module')
def validation(tmp_path_factory):
    path = tmp_path_factory.mktemp('fsl_mrs_validation')
    return run_validation(path), path


def test_coefficient_and_water_amplitude_recovery(validation):
    result, _ = validation
    for case in result['cases']:
        row = case['comparison']
        assert row['fitted_coefficient'] == pytest.approx(row['expected_coefficient'], rel=2e-6)
        assert row['water_amplitude'] == pytest.approx(row['expected_water_amplitude'], rel=2e-6)
        assert row['water_frequency_hz'] == pytest.approx(2.55462, abs=1e-5)
        assert abs(row['water_phase_rad']) < 1e-5
        assert abs(row['water_gamma_per_s']) < 1e-4
        assert case['water_optimizer']['success']


def test_concentration_linearity_and_distinct_standard_scaling(validation):
    result, _ = validation
    rows = [case['comparison'] for case in result['cases']]
    c = np.array([row['ground_truth_mM'] for row in rows])
    source_coeffs = np.array([row['unscaled_source_coefficient'] for row in rows])
    np.testing.assert_allclose(source_coeffs / source_coeffs[1], c / 10, rtol=2e-6)
    for row in rows:
        assert row['synthetic_water_referenced_mM'] == pytest.approx(row['ground_truth_mM'], rel=2e-6)
        assert row['fsl_molarity_mM'] == pytest.approx(row['independent_fsl_molarity_mM'], rel=2e-14)
        # Standard tissue/relaxation/integration scaling is a different model.
        assert row['fsl_molarity_mM'] / row['ground_truth_mM'] == pytest.approx(.6003318, rel=1e-5)


def test_metadata_preservation(validation):
    _, path = validation
    original = generate_pair({'NAA': 10})
    saved = json.loads((path / 'NAA_10mM' / 'simulation_config.json').read_text())
    assert saved == original['simulation_config']
    ground_truth = json.loads((path / 'NAA_10mM' / 'ground_truth.json').read_text())
    assert ground_truth == original['ground_truth']


def test_deterministic_fitting_and_residual(validation):
    result, path = validation
    assert result['configuration']['repeat_10mM_max_numeric_difference'] < 1e-10
    for case in result['cases']:
        assert case['comparison']['fit_relative_residual'] < 5e-5
        assert case['comparison']['baseline_max_abs'] == 0
    a = np.load(path / 'NAA_10mM' / 'fit_arrays.npz')
    b = np.load(path / 'repeat_10mM' / 'fit_arrays.npz')
    np.testing.assert_allclose(a['predicted'], b['predicted'], rtol=1e-12, atol=1e-12)

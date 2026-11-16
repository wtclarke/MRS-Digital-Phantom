import numpy as np
import pytest

from examples.fsl_mrs_water_fit_diagnostics import diagnose
from simulation.mrscloud_quantitative_pair import generate_pair, export_pair


@pytest.fixture(scope='module')
def diagnostics(tmp_path_factory):
    root = tmp_path_factory.mktemp('water_fit_diagnostics')
    export_pair(generate_pair({'NAA': 10}), root / 'input')
    return diagnose(root, source_dir=root / 'input')


def test_trace_and_chemical_shift(diagnostics):
    d = diagnostics
    assert d['max_complex_differences']['stored_vs_conjugate'] == 0
    assert d['max_complex_differences']['interpreted_vs_generator'] == 0
    assert d['max_complex_differences']['fitting_input_unscaled_vs_generator'] < 2e-13
    for name in ['generator', 'interpreted', 'fitting_input']:
        assert d['trace'][name]['frequency_hz'] == pytest.approx(2.55462, abs=1e-10)
        assert d['trace'][name]['initial_phase_rad'] == 0
        assert abs(d['trace'][name]['linewidth_hz']) < 1e-10
    assert d['trace']['stored']['frequency_hz'] == pytest.approx(-2.55462, abs=1e-10)
    assert abs(d['ppm_peak_bin'] - 4.70) < (2000 / 2048) / 127.731
    assert d['ppm_center_bin'] == 4.68


def test_original_phase5b_fit_and_integration_frame_are_distinct(diagnostics):
    case = diagnostics['cases'][1]['default_fit']
    assert case['relative_complex_residual'] < 1e-6
    assert case['integrated_model_relative_difference'] > .5
    assert case['fitted_amplitude'] == pytest.approx(61.73993435001698, rel=2e-6)


def test_unwrapped_sign_correct_initialization_recovers_all_cases(diagnostics):
    for case in diagnostics['cases']:
        fitted = case['unwrapped_initialization_fit']
        expected = case['expected']
        assert fitted['relative_complex_residual'] < 1e-6
        assert fitted['fitted_amplitude'] == pytest.approx(expected['initial_amplitude'], rel=2e-6)
        assert fitted['fitted_frequency_hz'] == pytest.approx(expected['frequency_hz'], abs=1e-6)
        assert abs(np.angle(np.exp(1j * (fitted['fitted_phase_rad'] - expected['initial_phase_rad'])))) < 1e-6


@pytest.mark.xfail(reason='Installed WaterRef regresses wrapped phase with the wrong initial frequency sign; case D converges to zero amplitude', strict=False)
def test_default_nonzero_phase_nonzero_frequency_recovery(diagnostics):
    assert diagnostics['cases'][3]['default_fit']['relative_complex_residual'] < 1e-6

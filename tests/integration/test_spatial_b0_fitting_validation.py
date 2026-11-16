"""Actual four-case CLI fitting, model-mismatch diagnostics and immutable inputs."""
import numpy as np
import pytest

from examples.mrscloud_spatial_b0_fitting_validation import run_validation, fit_case, SOURCE, CASES, NAMES, INPUTS
from examples.mrscloud_quantitative_voxel import BUNDLE
from simulation.mrscloud_full_basis import load_library


@pytest.fixture(scope='module')
def validation(tmp_path_factory):
    out=tmp_path_factory.mktemp('spatial_b0_fitting')/'results'
    return run_validation(out),out


def test_zero_and_constant_initial_amplitudes_and_explicit_initialization(validation):
    result,_=validation
    zero=result['cases']['zero_B0'];constant=result['cases']['constant_B0']
    assert zero['relative_complex_residual']<1e-6
    assert constant['default_fit']['relative_complex_residual']>.9
    assert constant['default_fit']['success'] # Success is explicitly not an accurate-fit assertion.
    assert constant['default_fit']['parameters']['NAA']==0
    assert constant['relative_complex_residual']<1e-11
    assert constant['parameters']['eps_0']==pytest.approx(-2*np.pi*37,abs=1e-10)
    for case in [zero,constant]:
        assert case['optimizer_success'] and case['water_fit']['success']
        assert abs(case['water_fit']['relative_amplitude_error'])<5e-6
        assert case['water_fit']['relative_complex_residual']<5e-6
        for c in case['coefficients'].values():assert abs(c['relative_coefficient_error'])<2e-6
    assert constant['water_fit']['fitted_frequency_hz']==pytest.approx(39.55462,abs=1e-5)
    assert constant['true_response_integrals']['NAA']<zero['true_response_integrals']['NAA']/50
    for n in NAMES:
        c=constant['coefficients'][n]
        assert c['standard_FSL_mM']==pytest.approx(c['source_response_reference_mM'],rel=5e-6)
    assert constant['coefficients']['NAA']['standard_FSL_mM']<.1


def test_non_lorentzian_envelopes_and_conditional_projections(validation):
    result,_=validation
    for name in ['symmetric_B0','mixed_relaxation_B0']:
        d=result['cases'][name]
        assert d['optimizer_success'] and d['water_fit']['success']
        assert d['relative_complex_residual']>.4
        assert d['water_fit']['relative_complex_residual']>.6
        assert d['water_fit']['relative_amplitude_error']<-.4
        assert not d['source_response_reference_is_exact_model_prediction']
        assert d['fixed_fitted_nuisance_projection_residual']>.4
        assert d['water_fixed_nuisance_projection_residual']>.6
        assert d['water_fixed_nuisance_projection_amplitude']==pytest.approx(d['water_fit']['fitted_amplitude'],rel=2e-7)
        for c in d['coefficients'].values():
            assert abs(c['relative_coefficient_error'])>.05
            assert c['fixed_fitted_nuisance_projection_mM']==pytest.approx(c['fitted_effective_mM'],rel=2e-5)
    for d in result['cases'].values():
        assert d['oracle_spatial_projection_residual']<1e-12
        assert d['water_oracle_residual']<1e-12
        assert d['water_oracle_amplitude']==pytest.approx(d['water_fit']['expected_initial_amplitude'],rel=1e-12)
        assert d['independent_metabolite_max_error']<1e-18 and d['independent_water_max_error']<1e-16
        for c in d['coefficients'].values():
            assert c['oracle_spatial_projection_mM']==pytest.approx(c['true_effective_mM'],rel=1e-12)
            assert c['independently_reconstructed_FSL_mM']==pytest.approx(c['standard_FSL_mM'],rel=3e-12)


def test_reports_calibration_and_original_arrays_preserved(validation):
    result,out=validation
    assert result['source_unchanged']
    assert result['constant_CLI_basis_alignment']['export_reload_max_complex_error']==0
    assert result['constant_CLI_basis_alignment']['magnitude_max_error']<1e-18
    for name in CASES:
        report=out/name/'fsl_mrs_results/report.html'
        assert report.stat().st_size>1000 and 'Water referencing' in report.read_text()
        for filename in INPUTS:assert (out/name/filename).read_bytes()==(SOURCE/name/filename).read_bytes()
        d=result['cases'][name]
        for n in NAMES:assert d['CLI_mM'][n]==pytest.approx(d['coefficients'][n]['standard_FSL_mM'],rel=2e-5)
        if name!='constant_B0':assert d['CLI_API_max_complex_fit_difference']<1e-18
    for filename in ['complex_spectral_fits.png','bias_and_residuals.png','residuals_and_water_envelopes.png','comparison.csv','results.json']:
        assert (out/filename).stat().st_size>1000
    assert (out/'cli_commands.txt').read_text().count('--report')==4


def test_deterministic_repeated_mismatched_fit(validation):
    result,out=validation
    responses,water,_=load_library(BUNDLE,selected_names=NAMES)
    original={k:v for k,v in result['cases']['symmetric_B0'].items() if not k.startswith('CLI')}
    assert fit_case(out/'symmetric_B0',responses,water)==original

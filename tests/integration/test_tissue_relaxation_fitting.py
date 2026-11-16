"""Actual Phase 6e fits and CLI reports; separate from analytical unit tests."""
import numpy as np
import pytest

from examples.mrscloud_quantitative_voxel import BUNDLE
from examples.mrscloud_tissue_relaxation_validation import run_validation, fit_case, NAMES
from simulation.mrscloud_full_basis import load_library


@pytest.fixture(scope='module')
def validation(tmp_path_factory):
    out=tmp_path_factory.mktemp('tissue_relaxation')/'validation'
    return run_validation(out),out


def test_fitted_attenuated_coefficients_water_and_independent_quantification(validation):
    result,_=validation
    assert result['bundle_unchanged']
    for data in result['cases'].values():
        w=data['water_fit']
        assert w['success']
        assert abs(w['relative_amplitude_error'])<1e-6
        assert w['relative_complex_residual']<1e-6
        assert w['fitted_frequency_hz']==pytest.approx(2.55462,abs=1e-5)
        if 'coefficients' not in data:
            assert 'undefined' in data['status']
            assert all(c==0 for c in data['truth']['relaxation_weighted_metabolite_concentrations_mM'].values())
            continue
        assert data['optimizer_success'],data['optimizer_message']
        assert data['relative_complex_residual']<5e-6
        for name,c in data['coefficients'].items():
            assert abs(c['relative_error'])<5e-6
            assert c['fitted_effective_mM']==pytest.approx(data['truth']['relaxation_weighted_metabolite_concentrations_mM'][name],rel=5e-6)
        for q in data['quantification'].values():
            for name in NAMES:
                assert q['reported_mM'][name]==pytest.approx(q['prediction']['predicted_mM'][name],rel=5e-6)
            np.testing.assert_allclose(list(q['prediction']['effective_metabolite_mM'].values()),
                [data['truth']['relaxation_weighted_metabolite_concentrations_mM'][n] for n in NAMES],rtol=5e-16)


def test_actual_reports_nominal_vs_weighted_metadata_and_exact_reload(validation):
    result,out=validation
    for row in result['analytical']:
        assert row['metabolite_max_error']<1e-19 and row['water_max_error']<1e-17
        assert row['export_max_error']==0
    for name in ['pure_WM','pure_GM','mixed']:
        data=result['cases'][name]
        report=out/name/'fsl_mrs_results/report.html'
        assert report.stat().st_size>1000 and 'Water referencing' in report.read_text()
        for n in NAMES:
            assert data['cli_mM'][n]==pytest.approx(data['quantification']['CLI_defaults']['reported_mM'][n],rel=1e-12)
        gt=data['truth']
        assert gt['calibration_identity']==result['calibration_identity']
        assert gt['metabolite_concentrations_mM']==gt['voxel_average_metabolite_concentrations_mM']
        assert gt['relaxation_weighted_water_concentration_mol_l']<gt['effective_voxel_water_concentration_mol_l']
    assert (out/'cli_commands.txt').read_text().count('--report')==3
    assert (out/'comparison.csv').stat().st_size>0


def test_scalar_fsl_model_limitations_and_deterministic_fitting(validation):
    result,out=validation
    responses,water,_=load_library(BUNDLE,selected_names=NAMES)
    # CLI values are added by the driver, not the Python fitting function.
    expected={k:v for k,v in result['cases']['mixed'].items() if k!='cli_mM'}
    assert fit_case(out/'mixed',responses,water)==expected
    for tissue in ['pure_WM','pure_GM']:
        d=result['cases'][tissue];q=d['quantification']['API_declared_water_NAA_scalar']
        assert q['prediction']['predicted_mM']['NAA']==pytest.approx(
            d['truth']['voxel_average_metabolite_concentrations_mM']['NAA']*q['prediction']['response_factor'],rel=2e-15)
    mixed=result['cases']['mixed'];p=mixed['quantification']['API_declared_water_NAA_scalar']['prediction']
    nominal=mixed['truth']['voxel_average_metabolite_concentrations_mM']['NAA']
    actual_effective_factor=p['effective_metabolite_mM']['NAA']/nominal
    assert abs(actual_effective_factor/p['quantification_metabolite_factor']-1)>.001

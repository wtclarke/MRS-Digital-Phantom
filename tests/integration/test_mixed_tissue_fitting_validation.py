"""Actual Newton fits and three CLI HTML reports; kept separate from fast predictions."""
import json

import numpy as np
import pytest

from examples.mrscloud_mixed_tissue_fitting_validation import run_validation, SOURCE, INPUTS, NAMES


@pytest.fixture(scope='module')
def validation(tmp_path_factory):
    out=tmp_path_factory.mktemp('mixed_tissue_fitting')/'results'
    return run_validation(out),out


def test_metabolite_and_water_recovery_and_independent_quantification(validation):
    result,_=validation
    for label,data in result['cases'].items():
        assert data['water_fit']['success'], (label,data['water_fit'])
        assert abs(data['water_fit']['relative_amplitude_error'])<1e-6
        assert data['water_fit']['relative_complex_residual']<1e-6
        assert data['water_fit']['fitted_frequency_hz']==pytest.approx(2.55462,abs=1e-5)
        assert abs(data['water_fit']['fitted_phase_rad'])<1e-5
        assert data['water_fit']['fitted_linewidth_hz']==pytest.approx(1,abs=1e-5)
        if label=='pure_CSF':
            assert 'unsupported' in data['status']
            assert 'coefficients' not in data
            assert 'csf_correction_error' in data
            continue
        assert data['optimizer_success'], (label,data['optimizer_message'])
        assert data['relative_complex_residual']<5e-5
        for coef in data['coefficients'].values():
            assert abs(coef['relative_error'])<2e-5
        for setting in ['defaults','matched']:
            for n in NAMES:
                assert data['quantification'][setting]['reported_mM'][n]==pytest.approx(
                    data['predictions'][setting]['predicted_mM'][n],rel=3e-5)


def test_linearity_gain_cancellation_and_deterministic_fits(validation):
    result,_=validation
    base=result['cases']['mixed']
    repeat=result['cases']['repeat']
    assert base['parameters']==repeat['parameters']
    assert base['water_fit']==repeat['water_fit']
    for label,cm,cw,g in [('double_metabolites',2,1,1),('double_water',1,2,1),('double_receiver',1,1,2)]:
        data=result['cases'][label]
        assert data['water_fit']['fitted_amplitude']/base['water_fit']['fitted_amplitude']==pytest.approx(cw*g,rel=1e-12)
        for n in NAMES:
            assert data['coefficients'][n]['raw_mol_l']/base['coefficients'][n]['raw_mol_l']==pytest.approx(cm*g,rel=1e-10)
            assert data['quantification']['defaults']['reported_mM'][n]/base['quantification']['defaults']['reported_mM'][n]==pytest.approx(cm/cw,rel=1e-9)


def test_cli_reports_provenance_and_unchanged_phase6b_arrays(validation):
    result,out=validation
    assert result['source_unchanged']
    assert (out/'comparison.csv').stat().st_size>0
    for label in ['pure_WM','pure_GM','mixed']:
        data=result['cases'][label]
        report=out/label/'fsl_mrs_results/report.html'
        assert report.stat().st_size>1000 and 'Water referencing' in report.read_text()
        assert data['cli_api_max_complex_difference']<1e-18
        for n in NAMES:
            assert data['cli_mM'][n]==pytest.approx(data['quantification']['defaults']['reported_mM'][n],rel=1e-12)
    for label in ['pure_WM','pure_GM','pure_CSF','mixed']:
        for filename in INPUTS:
            assert (out/label/filename).read_bytes()==(SOURCE/label/filename).read_bytes()
        truth=json.loads((out/label/'ground_truth.json').read_text())
        assert truth['calibration_identity']==result['calibration_identity']
    assert (out/'cli_commands.txt').read_text().count('--report')==3

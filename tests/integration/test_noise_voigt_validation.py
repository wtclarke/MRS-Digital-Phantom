"""Small actual Voigt/CLI run; the full 120-realization experiment is separate."""
import json
from pathlib import Path

import numpy as np
import pytest
from fsl_mrs.utils import fitting,mrs_io

from examples.mrscloud_noise_fitting_validation import run_validation,fit_one,NAMES,SOURCE
from examples.mrscloud_quantitative_voxel import BUNDLE
from simulation.mrscloud_full_basis import load_library


@pytest.fixture(scope='module')
def validation(tmp_path_factory):
    out=tmp_path_factory.mktemp('noise_voigt')/'smoke'
    return run_validation(out,smoke=True),out


def test_noise_free_voigt_baseline_and_actual_cli_reports(validation):
    result,out=validation
    base=result['records'][0]['fit']
    assert result['source_unchanged'] and len(result['records'])==5
    assert base['optimizer_success'] and base['relative_complex_residual']<2e-6
    assert 'sigma_0' in base['parameters']
    for c in base['coefficients'].values():
        assert abs(c['relative_error'])<5e-6
        assert c['covered_95'] is None
    for report in result['reports']:
        path=Path(report['path']);assert path.stat().st_size>1000 and 'Water referencing' in path.read_text()
        assert all(abs(v)<1e-12 for v in report['CLI_API_concentration_difference_mM'].values())
    assert '--lorentzian' not in (out/'cli_commands.txt').read_text()
    for filename in ['noise_bias_precision.png','noise_water_and_signal_errors.png','summary.csv','coefficients.csv']:
        assert (out/filename).stat().st_size>100


def test_noisy_fits_uncertainty_covariance_quantification_and_clean_preservation(validation):
    result,_=validation
    for record in result['records'][1:]:
        fitted=record['fit'];assert fitted['optimizer_success'] and fitted['quantification_success']
        assert fitted['independently_measured_SNR']==pytest.approx(float(record['snr']),rel=.06)
        assert fitted['water_fit']['expected_initial_amplitude']==result['records'][0]['fit']['water_fit']['expected_initial_amplitude']
        for n,c in fitted['coefficients'].items():
            assert c['reported_SD_mM']>0
            assert c['standard_FSL_mM']==pytest.approx(c['reconstructed_FSL_mM'],rel=2e-12)
            assert c['independent_clean_FSL_prediction_mM']==result['records'][0]['fit']['coefficients'][n]['independent_clean_FSL_prediction_mM']
        cov=np.array(fitted['coefficient_covariance_mM2']);names=fitted['coefficient_covariance_names']
        i,j=names.index('GPC'),names.index('PCh')
        assert fitted['coefficients']['GPC+PCh']['reported_SD_mM']**2==pytest.approx(cov[i,i]+cov[j,j]+2*cov[i,j],rel=2e-14)
        directory=Path(record['directory'])
        for filename in ['metabolites_clean','water_reference_clean']:
            clean=mrs_io.read_FID(str(directory/f'{filename}.nii.gz')).mrs().FID
            source=mrs_io.read_FID(str(SOURCE/f'{filename.removesuffix("_clean")}.nii.gz')).mrs().FID
            np.testing.assert_array_equal(clean,source)
        cfg=json.loads((directory/'simulation_config.json').read_text())
        assert cfg['noise_parameters']['root_seed']==record['seed']
        assert cfg['synthetic_calibration']['calibration_identity']==result['calibration_identity']


def test_deterministic_refitting_and_restart_cache(validation,monkeypatch):
    result,out=validation
    record=result['records'][1]
    responses,_,_=load_library(BUNDLE,selected_names=NAMES)
    expected=fit_one(Path(record['directory']),mrs_io.read_basis(str(BUNDLE/'metabolite_basis.BASIS')),
        responses,result['independent_clean_prediction'],result['records'][0]['fit'])
    assert expected==record['fit']
    def forbidden(*args,**kwargs):raise AssertionError('completed fit rerun instead of using cache')
    monkeypatch.setattr(fitting,'fit_FSLModel',forbidden)
    resumed=run_validation(out,smoke=True)
    assert resumed['cache_hits']==5 and resumed['records']==result['records']
    assert resumed['summary']==result['summary']

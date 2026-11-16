import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from examples.fsl_mrs_quantification_scaling_derivation import derive_case, run_derivation, spectral_integral

SOURCE=Path('outputs/fsl_mrs_water_scaling_phase5b')


@pytest.fixture(scope='module')
def derived(tmp_path_factory):
    return run_derivation(SOURCE,tmp_path_factory.mktemp('scaling_derivation'))


def test_four_concentration_reconstruction_and_factors(derived):
    for row in derived['cases']:
        assert abs(row['absolute_error_mM'])<5e-12
        assert abs(row['relative_error'])<5e-13
        assert abs(row['saved_integral_differences']['metabolite'])<1e-10
        assert abs(row['saved_integral_differences']['water'])<2e-8
        assert row['factor_product']*row['ground_truth_mM']==pytest.approx(row['reconstructed_fsl_mM'],rel=2e-15)
        assert row['ideal_prediction_mM']==pytest.approx(row['recorded_fsl_mM'],rel=2e-6)
        assert row['proposal_applied'] is False


def test_proposal_is_concentration_independent(derived):
    reference=derived['cases'][0]
    for row in derived['cases']:
        assert row['proposed_water_amplitude_default_FSL']==reference['proposed_water_amplitude_default_FSL']
        assert row['proposed_configuration_additional_scale']==reference['proposed_configuration_additional_scale']
        assert row['configuration_only_predicted_mM']==pytest.approx(row['ground_truth_mM'],rel=2e-6)
    # Removing relaxation/density assumptions still leaves the independently
    # computed Ref0ppm-versus-integral response mismatch.
    assert reference['ideal_response_factor'] != pytest.approx(1,rel=.01)


def test_does_not_call_final_quantification_or_change_saved_inputs(monkeypatch):
    from fsl_mrs.utils import quantify
    from fsl_mrs.utils.results import FitRes
    def forbidden(*args,**kwargs):
        raise AssertionError('derivation must not call final concentration reporting')
    monkeypatch.setattr(quantify,'quantifyWater',forbidden)
    monkeypatch.setattr(FitRes,'getConc',forbidden)
    directory=SOURCE/'NAA_10mM'
    files=list(directory.glob('*.json'))+list(directory.glob('*.nii.gz'))
    before={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    a=derive_case(directory)
    b=derive_case(directory)
    assert a==b
    assert before=={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def test_proposed_constants_do_not_depend_on_reported_or_fitted_concentrations(tmp_path):
    original=derive_case(SOURCE/'NAA_10mM')
    for name in ['simulation_config.json','ground_truth.json','fitting_results.json']:
        (tmp_path/name).write_bytes((SOURCE/'NAA_10mM'/name).read_bytes())
    fit=json.loads((tmp_path/'fitting_results.json').read_text())
    fit['comparison']['fsl_molarity_mM']=999
    fit['fit_parameters']['NAA']*=.5
    fit['water_optimizer']['params'][0]*=3
    (tmp_path/'fitting_results.json').write_text(json.dumps(fit))
    changed=derive_case(tmp_path)
    for key in ['proposed_water_amplitude_default_FSL','proposed_water_amplitude_no_relaxation_40M','proposed_configuration_additional_scale']:
        assert original[key]==changed[key]


def test_first_point_fft_and_bin_integral_units():
    impulse=np.zeros(8,dtype=complex);impulse[0]=1
    # Half the first sample, orthonormal FFT: each real spectral bin is
    # .5/sqrt(8); trapezoid across all eight bins spans seven intervals.
    assert spectral_integral(impulse,np.arange(8),(0,7))==pytest.approx(7*.5/np.sqrt(8),rel=1e-15)
    assert impulse[0]==1

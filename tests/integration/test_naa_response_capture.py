import hashlib
import json
import csv
from pathlib import Path

import numpy as np
import pytest

from examples.mrscloud_naa_response_decomposition import capture, first_dc, decompose


@pytest.fixture(scope='module')
def captured(tmp_path_factory):
    source=Path('simulation/basissets/MRSCloud-main/functions/sim_signals.m')
    before=hashlib.sha256(source.read_bytes()).hexdigest()
    trace=capture(tmp_path_factory.mktemp('naa_capture'),spatial_points=3,points=256)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==before
    return trace


def test_actual_three_spin_capture_frequency_superposition_and_determinism(captured):
    raw=captured['raw'];time=np.arange(256)/4000
    np.testing.assert_allclose(raw['full'],raw['methyl']+raw['aspartyl'],rtol=1e-12,atol=1e-12)
    np.testing.assert_array_equal(raw['methyl'],raw['methyl_repeat'])
    for label,carrier in [('methyl',3.),('methyl_water_carrier',4.7),('methyl_on_carrier',2.008)]:
        expected=raw[label][0]*np.exp(-np.pi*time)*np.exp(2j*np.pi*(2.008-carrier)*127.728*time)
        np.testing.assert_allclose(raw[label],expected,rtol=1e-11,atol=1e-12)
        sys=captured['configs'][label]['spin_system_carrier_relative']
        assert len(sys['shifts'])==3 and sys['scaleFactor']==1
        np.testing.assert_array_equal(sys['J'],np.zeros((3,3)))
        np.testing.assert_allclose(captured['first_dc'][label],first_dc(raw[label],captured['ppm'][label]),rtol=1e-11,atol=1e-12)


@pytest.fixture(scope='module')
def experiment():
    root=Path('outputs/mrscloud_naa_response_decomposition/validated')
    return json.loads((root/'results.json').read_text()),np.load(root/'responses.npz'),root


def test_full_resolution_alpha_and_proton_response(experiment):
    result,arrays,_=experiment
    assert result['source_hashes_unchanged']
    assert result['methyl_repeat_max_error']==0
    assert result['raw_superposition_max_error']<2e-14
    assert result['carrier_sensitivity']['on_resonance_methyl_per_proton_vs_water']==pytest.approx(1,rel=1e-12)
    d=decompose(arrays,arrays['raw_water_target'])
    assert d['product']==pytest.approx(result['phase5h_alpha'],rel=5e-14)
    assert d['direct_alpha']==pytest.approx(d['product'],rel=2e-15)
    assert result['formatted_vs_existing_relative_error']<3e-7
    assert .94<result['carrier_sensitivity']['methyl_RF4p7_vs_RF3_magnitude']<.95


def test_actual_cli_reports_and_unmodified_water(experiment):
    result,_,root=experiment
    from fsl_mrs.utils import mrs_io
    source=mrs_io.read_FID('outputs/mrscloud_sequence_water_integration/validated/mrscloud_sequence/water_reference.nii.gz').mrs().FID
    for label,fit in result['fits'].items():
        assert fit['success'] and fit['water']['success']
        assert fit['metabolite_relative_residual']<5e-5
        assert fit['water']['relative_complex_residual']<1e-6
        assert fit['raw_coefficient_mol_l']==pytest.approx(.010,rel=2e-6)
        assert (root/label/'fsl_mrs_results/report.html').stat().st_size>1000
        with (root/label/'fsl_mrs_results/concentrations.csv').open() as stream:
            rows=list(csv.reader(stream))
        cli=dict(zip(rows[1][1:5],map(float,rows[2][1:5])))
        assert cli['molarity']==fit['standard_mM']
        assert cli['raw']==fit['raw_coefficient_mol_l']
        np.testing.assert_array_equal(mrs_io.read_FID(str(root/label/'water_reference.nii.gz')).mrs().FID,source)
    assert result['fits']['full_formatted']['standard_mM']==pytest.approx(6.885791303512254,rel=1e-12)

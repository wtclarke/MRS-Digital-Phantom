"""Fast, independent numerical tests; no LCModel executable needed."""
import json
from pathlib import Path

import numpy as np
import pytest

from examples.lcmodel_noise_comparison import BASE,canonical,control_for_case,execute_case
from utils.lcmodel_noise_analysis import descriptive,paired_bootstrap,parse_fit
from utils.lcmodel_quantitative_io import write_raw,parse_raw


def fixture_outputs(tmp_path):
    table=' LCModel (Version 6.3-1N)\n'
    table+='\n'.join(f' {c} {sd}% 1.0 {n}' for n,c,sd in [('NAA',10.6,0),('Cr',8.4,1),('GPC',2.1,8),('PCh',1.0,12),('Glu',8.3,6),('GPC+PCh',3.1,3)])
    (tmp_path/'fit.TABLE').write_text(table+'\n 1 warning MYBASI 8\n')
    (tmp_path/'fit.PRINT').write_text('FCALIB = 1.21933E+05\n FWHM = .004 ppm S/N = 42\nData shift = .001 ppm\n')
    (tmp_path/'fit.COORD').write_text('3 points on ppm-axis = NY\n 3 2 1\nNY phased data points follow\n 4 5 6\nNY points of the fit to the data follow\n 3.9 4.9 5.9\nNY background values follow\n .1 .2 .3\n')
    return tmp_path


def test_parse_grouped_uncertainty_is_reported_not_quadrature(tmp_path):
    result=parse_fit(fixture_outputs(tmp_path))
    assert result['valid']
    group=result['values']['GPC+PCh']
    assert group['concentration']==pytest.approx(3.1)
    assert group['reported_group_percent_SD']==3
    assert group['uncertainty_mM']==pytest.approx(3.1*.03)
    assert result['values']['NAA']['zero_percent_SD']
    assert result['values']['NAA']['uncertainty_mM'] is None
    assert result['values']['NAA']['displayed_SD_product_mM']==0
    assert result['values']['GPC']['concentration_text']=='2.1'
    assert result['reported_SNR']==42 and result['FCALIB']==121933


def test_fatal_diagnostic_overrides_apparent_concentrations(tmp_path):
    fixture_outputs(tmp_path)
    with (tmp_path/'fit.TABLE').open('a') as stream:stream.write(' 1 FATAL TWOREG 4\n')
    result=parse_fit(tmp_path)
    assert not result['valid'] and result['fatal']


def test_missing_concentration_is_incomplete(tmp_path):
    (tmp_path/'fit.TABLE').write_text('LCModel (Version 6.3-1N)\n 10.6 0% 1.0 NAA\n')
    assert not parse_fit(tmp_path)['valid']


def test_statistics_hand_calculated_and_failure_denominator():
    result=descriptive([9,10,13],10,11,[.4,None,.6],4,1)
    assert result['count']==3 and result['failure_rate']==.25
    assert result['mean_mM']==pytest.approx(32/3)
    assert result['nominal_bias_mM']==pytest.approx(2/3)
    assert result['noise_bias_mM']==pytest.approx(-1/3)
    assert result['SD_mM']==pytest.approx(np.sqrt(13/3))
    assert result['nominal_RMSE_mM']==pytest.approx(np.sqrt(10/3))
    assert result['noise_RMSE_mM']==pytest.approx(np.sqrt(3))
    assert result['median_reported_uncertainty_mM']==.5
    assert result['nominal_poor_fraction']==pytest.approx(1/3)
    assert descriptive([],10,11,[],2,2)==dict(count=0,total=2,failure_rate=1)


def test_paired_bootstrap_baselines_and_determinism():
    lcm=np.array([11.,12.,10.]);fsl=np.array([7.,8.,6.])
    before=np.random.get_state()
    result=paired_bootstrap(lcm,fsl,11.,7.,10.,seed=123,draws=2000)
    assert result['baseline_centered_estimate_mean']==0
    assert result['noise_absolute_error_mean']==0
    assert result['noise_absolute_error_CI_low']==result['noise_absolute_error_CI_high']==0
    assert result['nominal_absolute_error_mean']==pytest.approx(-2)
    assert result==paired_bootstrap(lcm,fsl,11.,7.,10.,seed=123,draws=2000)
    after=np.random.get_state();assert before[0]==after[0];np.testing.assert_array_equal(before[1],after[1])
    with pytest.raises(ValueError):paired_bootstrap([1],[1,2],1,1,1)


def test_frozen_control_changes_only_paths(tmp_path):
    if not (BASE/'lcmodel.CONTROL').exists():pytest.skip('Phase 7b protocol unavailable')
    original=(BASE/'lcmodel.CONTROL').read_text()
    changed=control_for_case(original,tmp_path,tmp_path/'attempt_001')
    assert canonical(changed)==canonical(original)
    for setting in ['WCONC=39016.5','ATTH2O=1.','ATTMET=1.','FWHMBA=0.007828953034110749']:
        assert setting in changed
    assert 'VITRO' not in changed and 'DEEXT2=' not in changed


def test_complex_noisy_raw_preservation(tmp_path):
    rng=np.random.default_rng(77);t=np.arange(2048)*.0005
    x=np.exp(-np.pi*t+1j*(2*np.pi*37*t+.3))+.01*(rng.normal(size=2048)+1j*rng.normal(size=2048))
    path=tmp_path/'signal.RAW';write_raw(path,x,frequency_mhz=127.731)
    stored,_=parse_raw(path);np.testing.assert_array_equal(stored.conj(),x)


def case_setup(tmp_path,monkeypatch):
    directory=tmp_path/'case';directory.mkdir()
    for name in ['metabolites.RAW','water_reference.H2O']:(directory/name).write_text('fixed converted samples')
    monkeypatch.setattr('examples.lcmodel_noise_comparison.prepare_input',lambda *args:{'conversion_error':0.})
    protocol=' $LCMODL\n'+''.join(f" {field}='old'\n" for field in ['FILRAW','FILH2O','FILBAS','FILPS','FILPRI','FILTAB','FILCOO'])+' $END\n'
    return directory,protocol,dict(hashes={})


def test_completed_cache_integrity_and_no_unnecessary_retry(tmp_path,monkeypatch):
    from types import SimpleNamespace
    directory,protocol,record=case_setup(tmp_path,monkeypatch)
    def fake_run(*args,**kwargs):
        fixture_outputs(kwargs['cwd']/'lcmodel_outputs')
        return SimpleNamespace(returncode=0,stdout='captured',stderr='')
    monkeypatch.setattr('examples.lcmodel_noise_comparison.subprocess.run',fake_run)
    first,cached=execute_case(record,directory,protocol,'frozen',30)
    assert first['status']=='complete' and not cached
    def forbidden(*args,**kwargs):raise AssertionError('completed execution retried')
    monkeypatch.setattr('examples.lcmodel_noise_comparison.subprocess.run',forbidden)
    second,cached=execute_case(record,directory,protocol,'frozen',30)
    assert cached and first==second
    with pytest.raises(ValueError,match='identity'):execute_case(record,directory,protocol,'different',30)
    (Path(first['attempt_directory'])/'lcmodel_outputs/fit.COORD').write_text('corruption')
    with pytest.raises(ValueError,match='artifact'):execute_case(record,directory,protocol,'frozen',30)


def test_timed_out_attempt_preserved_and_retried_unchanged(tmp_path,monkeypatch):
    import subprocess
    from types import SimpleNamespace
    directory,protocol,record=case_setup(tmp_path,monkeypatch)
    def timeout(*args,**kwargs):raise subprocess.TimeoutExpired(args[0],30,output=b'partial')
    monkeypatch.setattr('examples.lcmodel_noise_comparison.subprocess.run',timeout)
    first,_=execute_case(record,directory,protocol,'frozen',30)
    assert first['status']=='interrupted' and first['execution']['stdout']=='partial'
    def complete(*args,**kwargs):
        fixture_outputs(kwargs['cwd']/'lcmodel_outputs')
        return SimpleNamespace(returncode=0,stdout='',stderr='')
    monkeypatch.setattr('examples.lcmodel_noise_comparison.subprocess.run',complete)
    second,cached=execute_case(record,directory,protocol,'frozen',30)
    assert second['status']=='complete' and second['attempt']==2 and not cached
    assert (directory/'attempt_001/execution.json').exists()
    assert canonical((directory/'attempt_001/lcmodel.CONTROL').read_text())==canonical((directory/'attempt_002/lcmodel.CONTROL').read_text())


def test_cached_fatal_remains_failure_without_parameter_retry(tmp_path,monkeypatch):
    from types import SimpleNamespace
    directory,protocol,record=case_setup(tmp_path,monkeypatch)
    def fatal(*args,**kwargs):
        (kwargs['cwd']/'lcmodel_outputs/fit.TABLE').write_text('LCModel (Version 6.3-1N)\n 1 FATAL TWOREG 4\n')
        return SimpleNamespace(returncode=0,stdout='',stderr='')
    monkeypatch.setattr('examples.lcmodel_noise_comparison.subprocess.run',fatal)
    first,_=execute_case(record,directory,protocol,'frozen',30)
    assert first['status']=='failed' and not first['fit']['valid']
    def forbidden(*args,**kwargs):raise AssertionError('terminal fatal fit retried')
    monkeypatch.setattr('examples.lcmodel_noise_comparison.subprocess.run',forbidden)
    second,cached=execute_case(record,directory,protocol,'frozen',30)
    assert cached and second['status']=='failed'


def test_inventory_rejects_missing_realizations(tmp_path,monkeypatch):
    import examples.lcmodel_noise_comparison as comparison
    (tmp_path/'results.json').write_text(json.dumps({'records':[{'snr':'infinite','realization':0}]}))
    monkeypatch.setattr(comparison,'SOURCE',tmp_path)
    with pytest.raises(ValueError,match='incomplete'):comparison.inventory()

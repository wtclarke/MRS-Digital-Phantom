"""Validate actual local Monte Carlo artifacts without executing proprietary fits."""
import json
from pathlib import Path

import numpy as np
import pytest

from examples.lcmodel_noise_comparison import BASE,OUTPUT,canonical,digest,execute_case,inventory
from utils.lcmodel_quantitative_io import parse_raw


@pytest.fixture(scope='module')
def comparison():
    path=OUTPUT/'results.json'
    if not path.exists():pytest.skip('requires completed local Phase 7c experiment')
    return json.loads(path.read_text())


def test_all_121_actual_fits_and_protocol_identity(comparison):
    assert comparison['fits']==comparison['attempts']==comparison['valid_fits']==121
    assert comparison['failures']==0 and comparison['source_files_unchanged']
    protocol=(BASE/'lcmodel.CONTROL').read_text()
    frozen=json.loads((OUTPUT/'frozen_configuration.json').read_text())
    assert frozen['basis_sha256']==digest(BASE/'lcmodel_compatible.BASIS')
    assert frozen['CONTROL_sha256']==digest(BASE/'lcmodel.CONTROL')
    assert len(comparison['records'])==121
    for record in comparison['records']:
        result=record['result']
        assert result['status']=='complete' and result['execution']['returncode']==0
        control=Path(result['attempt_directory'])/'lcmodel.CONTROL'
        assert canonical(control.read_text())==canonical(protocol)
        for path,expected in result['artifact_hashes'].items():assert digest(path)==expected


def test_all_source_identities_and_complex_conversions(comparison):
    from fsl_mrs.utils import mrs_io
    frozen=json.loads((OUTPUT/'frozen_configuration.json').read_text())
    assert frozen['source_hashes']=={p:digest(p) for p in frozen['source_hashes']}
    for record in comparison['records']:
        result=record['result'];directory=Path(result['attempt_directory']).parent
        for name,filename in [('metabolites','metabolites.RAW'),('water_reference','water_reference.H2O')]:
            source=mrs_io.read_FID(str(Path(record['source'])/f'{name}.nii.gz')).mrs().FID
            raw,_=parse_raw(directory/filename)
            np.testing.assert_array_equal(raw.conj(),source)
            assert result['conversion'][name]['max_absolute_complex_error']==0
        assert result['conversion']['initial_ratio_error']==0
        if record['noise']:
            assert record['noise']['root_seed']==record['seed']
            assert result['conversion']['metabolites']['realized_snr']==pytest.approx(record['noise']['metabolite']['realized_snr'],rel=1e-11)


def test_grouped_results_and_paired_statistics(comparison):
    for record in comparison['records']:
        values=record['result']['fit']['values']
        assert values['GPC+PCh']['concentration']==values['GPC']['concentration']+values['PCh']['concentration']
        assert values['GPC+PCh']['reported_group_percent_SD'] is not None
    assert len(comparison['paired'])==72
    assert all(p['pairs']==20 for p in comparison['paired'])
    for level in ['100','50','25','15','10','5']:
        rows=[s for s in comparison['summary'] if s['snr']==level]
        assert len(rows)==18 and all(s['count']==s['total']==20 and s['failure_rate']==0 for s in rows)
    assert len(comparison['summary'])==126


def test_cached_execution_and_corrupt_source_rejection(comparison,monkeypatch):
    import examples.lcmodel_noise_comparison as experiment
    def forbidden(*args,**kwargs):raise AssertionError('completed proprietary fit rerun')
    monkeypatch.setattr(experiment.subprocess,'run',forbidden)
    for record in comparison['records']:
        result=record['result'];directory=Path(result['attempt_directory']).parent
        resumed,cached=execute_case(record,directory,(BASE/'lcmodel.CONTROL').read_text(),result['signature']['config_hash'],180)
        assert cached and resumed==result
    original=experiment.digest
    def wrong(path):return 'incorrect' if Path(path).name=='metabolites.nii.gz' else original(path)
    monkeypatch.setattr(experiment,'digest',wrong)
    with pytest.raises(ValueError,match='input hash mismatch'):inventory()


def test_reports_and_summaries_exist(comparison):
    assert len(comparison['reports'])>=6
    for report in comparison['reports']:
        assert Path(report['PDF']).stat().st_size>1000
        if report['FSL_HTML']:assert Path(report['FSL_HTML']).stat().st_size>1000
    for filename in ['per_realization.csv','summary.csv','paired_comparison.csv','warning_summary.json','failure_summary.json','paired_error_comparison.png','paired_deviations.png','uncertainty_comparison.png']:
        assert (OUTPUT/filename).stat().st_size>100

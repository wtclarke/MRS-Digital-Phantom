import json
from pathlib import Path
import shutil

import numpy as np
import pytest
from fsl_mrs.utils import mrs_io

from simulation.mrscloud_full_basis import load_library, synthesize, sha256, normalization_identity
from simulation.mrscloud_quantitative_pair import export_pair

ROOT=Path('outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative')
SMOKE=Path('outputs/mrscloud_full_basis_smoke_final/mrscloud_press_te30_quantitative')


@pytest.fixture(scope='module')
def library():
    responses,water,metadata=load_library(ROOT)
    return responses,water,metadata,json.loads((ROOT/'validation_results.json').read_text())


def test_complete_manifest_reference_scalar_and_export_precision(library):
    responses,water,metadata,validation=library
    manifest=json.loads((ROOT/'reference_manifest.json').read_text())
    norm=metadata['normalization']
    assert len(responses)==27 and len(manifest['entries'])==29
    assert manifest['unsupported_entries']==manifest['excluded_MM_lipid_entries']==[]
    assert all(e['source_spin_system']=='sys'+e['name'] for e in manifest['entries'])
    assert norm['calibration_identity']==normalization_identity(norm)
    with np.load(ROOT/'normalization_reference.npz') as data:
        expected=np.max(np.fft.fftshift(np.fft.fft(data['processed_reference'],axis=0),axes=0).real)
    assert norm['scalar']==expected
    assert all(x.shape==(2048,) and np.iscomplexobj(x) and np.all(np.isfinite(x)) for x in [*responses.values(),water])
    assert max(validation['basis_serialization_relative_errors'].values())<3e-7
    assert metadata['sequence']['water_proton_multiplier_applied'] is False
    assert 'full-reference' in metadata['calibration_identity']


@pytest.mark.parametrize('name',['NAA','Cr','Glu'])
def test_selection_independence_and_deterministic_metadata(library,name):
    responses,water,metadata,_=library
    selected,w,m=load_library(ROOT,selected_names=[name],expected_calibration=metadata)
    together,w2,m2=load_library(ROOT,selected_names=[name,'GPC','PCh'])
    np.testing.assert_array_equal(selected[name],responses[name])
    np.testing.assert_array_equal(together[name],responses[name])
    np.testing.assert_array_equal(w,water);np.testing.assert_array_equal(w2,water)
    assert m==m2==metadata


def test_explicit_concentrations_linearity_and_mixture_sum(library):
    responses,water,_,_=library
    concentrations={'NAA':10.,'Cr':8.,'GPC':2.,'PCh':1.,'Glu':12.}
    a=synthesize(ROOT,concentrations)
    b=synthesize(ROOT,{n:2*c for n,c in concentrations.items()},water_concentration_mol_l=80)
    g=synthesize(ROOT,concentrations,receiver_factor=2)
    expected=sum(responses[n]*(c/1000) for n,c in concentrations.items())
    np.testing.assert_array_equal(a['metabolite_fid'],expected)
    np.testing.assert_array_equal(a['water_reference_fid'],water*40)
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_array_equal(b[key],2*a[key])
        np.testing.assert_array_equal(g[key],2*a[key])
    for n in concentrations:np.testing.assert_array_equal(b['components'][n],2*a['components'][n])
    only_water=synthesize(ROOT,concentrations,water_concentration_mol_l=80)
    np.testing.assert_array_equal(only_water['metabolite_fid'],a['metabolite_fid'])
    np.testing.assert_array_equal(only_water['water_reference_fid'],2*a['water_reference_fid'])
    for kwargs in [{'concentrations_mM':{'H2O':1}},{'concentrations_mM':{'NAA':-1}}]:
        with pytest.raises(ValueError):synthesize(ROOT,**kwargs)


def test_paired_complex_export_and_metadata(library,tmp_path):
    pair=synthesize(ROOT,{'NAA':10,'Cr':8,'GPC':2,'PCh':1,'Glu':12})
    paths=export_pair(pair,tmp_path)
    for key,path in paths.items():
        data=mrs_io.read_FID(str(path))
        np.testing.assert_array_equal(data.mrs().FID,pair[key])
        assert data.dwelltime==.0005 and data.SpecFreqChemShift==4.68
        assert data.spectrometer_frequency==[127.731]
    assert json.loads((tmp_path/'ground_truth.json').read_text())==pair['ground_truth']


@pytest.mark.parametrize('change',['scalar','statistic','order','missing','source','capture','water','reference'])
def test_rejects_rebound_but_incompatible_calibration_metadata(tmp_path,change):
    root=tmp_path/'bundle';shutil.copytree(ROOT,root)
    norm=json.loads((root/'normalization.json').read_text())
    manifest=json.loads((root/'reference_manifest.json').read_text())
    if change=='scalar':norm['scalar']*=2
    elif change=='statistic':norm['statistic']='mean_absolute_selected_FID'
    elif change=='order':norm['reference_names'].reverse()
    elif change=='source':norm['source_hashes']={}
    elif change=='capture':norm['capture_configuration_signature']='0'*64
    elif change=='missing':del norm['scalar']
    elif change=='water':
        with np.load(root/'water_reference.npz') as a:raw=a['raw_native'];fid=a['fid']
        np.savez(root/'water_reference.npz',raw_native=raw,fid=2*fid)
        manifest['files']['water_reference.npz']=sha256(root/'water_reference.npz')
    else:
        with np.load(root/'normalization_reference.npz') as a:reference=a['processed_reference']
        reference[:,0]*=2
        np.savez(root/'normalization_reference.npz',processed_reference=reference)
        norm['reference_array_sha256']=sha256(root/'normalization_reference.npz')
        manifest['files']['normalization_reference.npz']=norm['reference_array_sha256']
    norm['calibration_identity']=normalization_identity(norm)
    (root/'normalization.json').write_text(json.dumps(norm))
    manifest['files']['normalization.json']=sha256(root/'normalization.json')
    (root/'reference_manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError):load_library(root)


def test_rejects_other_sampling_selection_and_expected_identity(library):
    for kwargs in [dict(points=1024),dict(bandwidth_hz=4000),dict(selected_names=['Cho']),
                   dict(selected_names=['NAA','NAA']),dict(expected_calibration={})]:
        with pytest.raises(ValueError):load_library(ROOT,**kwargs)
    _,_,smoke_metadata=load_library(SMOKE,points=128)
    assert smoke_metadata['calibration_identity']!=library[2]['calibration_identity']
    assert smoke_metadata['sequence']['resolution_family']=='smoke'
    with pytest.raises(ValueError):load_library(ROOT,expected_calibration=smoke_metadata)


def test_missing_raw_cache_and_changed_capture_configuration_are_explicit(tmp_path):
    from examples.mrscloud_full_basis_water_reference import build, prepare_capture
    with pytest.raises(FileNotFoundError,match='run without --reuse-raw'):
        build(tmp_path/'missing',smoke=True,reuse_raw=True)
    path=tmp_path/'missing/capture_configuration.json'
    cfg=json.loads(path.read_text());cfg['native_points']=8192;path.write_text(json.dumps(cfg))
    with pytest.raises(ValueError,match='different source/configuration'):
        prepare_capture(tmp_path/'missing',smoke=True)

"""Production resource integrity and relocation, with no research outputs."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

import numpy as np
import pytest
from fsl_mrs.utils import mrs_io

from simulation.calibration_resources import calibration_bundle,default_calibration_bundle,verify_source_provenance
from simulation.mrscloud_full_basis import load_library
from simulation.headless import run_headless
from simulation.quantitative import simulate_quantitative,export_quantitative,CALIBRATION_ID

ROOT=Path(__file__).resolve().parents[1]


def config():
    return json.loads((ROOT/'examples/configs/quantitative/mixed_tissue.json').read_text())


def test_original_frozen_bytes_and_sources():
    root=default_calibration_bundle()
    hashes=json.loads((ROOT/'calibration/frozen_file_hashes.json').read_text())
    for name,expected in hashes.items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected
    norm=json.loads((root/'normalization.json').read_text())
    assert norm['scalar']==8857.353925964408 and norm['calibration_identity']==CALIBRATION_ID
    verify_source_provenance(root,norm)


def test_relocation_preserves_arrays_and_exports_exactly(tmp_path):
    root=default_calibration_bundle();moved=tmp_path/'arbitrary_location'
    shutil.copytree(root,moved)
    before=load_library(root);after=load_library(moved)
    assert before[2]==after[2]
    for name in before[0]:np.testing.assert_array_equal(before[0][name],after[0][name])
    np.testing.assert_array_equal(before[1],after[1])
    cfg=config()['quantitative'];original=simulate_quantitative(cfg)
    relocated=simulate_quantitative(dict(cfg,calibration_bundle=str(moved)))
    assert original['ground_truth']==relocated['ground_truth']
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(original[key],relocated[key])
    paths=[export_quantitative(p,tmp_path/n) for n,p in [('original',original),('moved',relocated)]]
    for label,key in [('metabolites','metabolite_fid'),('water_reference','water_reference_fid')]:
        a,b=[mrs_io.read_FID(p[label]).mrs().FID for p in paths]
        np.testing.assert_array_equal(a,b);np.testing.assert_array_equal(a,original[key])


def test_provenance_corruption_rejected(tmp_path):
    norm=json.loads((default_calibration_bundle()/'normalization.json').read_text())
    with zipfile.ZipFile(default_calibration_bundle()/'provenance_sources.zip') as source:
        with zipfile.ZipFile(tmp_path/'provenance_sources.zip','w') as target:
            for i,name in enumerate(source.namelist()):target.writestr(name,b'changed' if i==0 else source.read(name))
    with pytest.raises(ValueError,match='provenance identity mismatch'):verify_source_provenance(tmp_path,norm)


def test_missing_calibration_actionable(monkeypatch,tmp_path):
    import simulation.calibration_resources as resources
    monkeypatch.setattr(resources,'default_calibration_bundle',lambda:tmp_path/'absent')
    with pytest.raises(FileNotFoundError,match='No download or regeneration'):calibration_bundle()


def test_incompatible_calibration_cannot_borrow_source_snapshot(tmp_path):
    norm=json.loads((default_calibration_bundle()/'normalization.json').read_text())
    norm['calibration_identity']='different'
    with pytest.raises(ValueError,match='no matching frozen source provenance'):verify_source_provenance(tmp_path,norm)


@pytest.mark.parametrize('name',['pure_WM','mixed_tissue','relaxation','spatial_b0','noisy'])
def test_supported_headless_configurations(name,tmp_path):
    cfg=json.loads((ROOT/f'examples/configs/quantitative/{name}.json').read_text())
    cfg['output_dir']=str(tmp_path/name)
    result=run_headless(cfg)[0]
    assert result['simulation_config']['synthetic_calibration']['calibration_identity']==CALIBRATION_ID
    for label,key in [('metabolites','metabolite_fid'),('water_reference','water_reference_fid')]:
        reader=mrs_io.read_FID(result['saved_paths'][label])
        np.testing.assert_array_equal(reader.mrs().FID,result[key])
        assert reader.dwelltime==.0005
        assert reader.hdr_ext['SpectrometerFrequency']==[127.731]


def test_runtime_does_not_import_research_or_start_octave():
    import subprocess
    import sys
    code='''
import importlib.abc,sys,json
class NoResearch(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname.split('.')[0] in {'examples','research','oct2py'} or fullname=='simulation.run_matlab':
            raise AssertionError('unexpected runtime dependency: '+fullname)
sys.meta_path.insert(0,NoResearch())
from simulation.headless import run_headless
cfg=json.load(open('examples/configs/quantitative/pure_WM.json'))
cfg.pop('output_dir')
assert len(run_headless(cfg)[0]['metabolite_fid'])==2048
'''
    result=subprocess.run([sys.executable,'-c',code],cwd=ROOT,capture_output=True,text=True)
    assert result.returncode==0,result.stderr


def test_noise_statistics_and_direct_snr_reference():
    from simulation.mrscloud_quantitative_noise import complex_gaussian_noise,reference_peak
    state=np.random.get_state()
    noise=complex_gaussian_noise(200000,sigma=2.,seed=8671)
    assert abs(noise.real.mean())<.02 and abs(noise.imag.mean())<.02
    assert abs(noise.real.std()/2-1)<.008 and abs(noise.imag.std()/2-1)<.008
    assert abs(np.corrcoef(noise.real,noise.imag)[0,1])<.008
    np.testing.assert_array_equal(noise,complex_gaussian_noise(200000,sigma=2.,seed=8671))
    cfg=config()['quantitative'];cfg['noise']={'seed':18,'metabolite_snr':50,'water_snr':1000}
    pair=simulate_quantitative(cfg);c=pair['simulation_config'];fid=pair['components']['NAA']
    freq=np.fft.fftshift(np.fft.fftfreq(2048,.0005));ppm=4.68+freq/127.731
    bins=np.where((ppm>=1.8)&(ppm<=2.2))[0]
    independent=max(abs(np.exp(-2j*np.pi*freq[bins,None]*pair['time_s'])@fid))/np.sqrt(2048)
    np.testing.assert_allclose(independent,reference_peak(fid,c,(1.8,2.2)),rtol=2e-13)
    assert c['noise_parameters']['metabolite']['sigma_per_component']==pytest.approx(independent/50,rel=2e-13)
    assert not np.array_equal(pair['metabolite_fid']-pair['clean_metabolite_fid'],pair['water_reference_fid']-pair['clean_water_reference_fid'])
    for actual,old in zip(np.random.get_state(),state):np.testing.assert_array_equal(actual,old)

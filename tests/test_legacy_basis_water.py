from copy import deepcopy
import numpy as np
import pandas as pd
import pytest
import torch
from fsl_mrs.utils import mrs_io
from test_basis_h2o import basis
from test_quantitative_production_legacy import anatomy
from simulation import mrs_phantom
from simulation.signalModel import SignalModel
from simulation.headless import simulate_batch,results_to_nifti_mrs
from simulation.legacy_water_reference import prepare


def settings(path):
    return dict(path2basis=str(path),basis_set_settings={'TE':35,'vendor':'test','localization':'PRESS'},voxel_size_mm=[4,4,4],
        simulation_params=dict(bandwidth=2000,spectral_points=256,TR=2000,lipid_sigma=1,mm_level=.5,lipid_amp_factor=.1,
        noise_level=.001,water_amp_factor=.3,water_damping_min=1,water_damping_max=4,water_phase_min=-30,water_phase_max=30),
        water_reference_mode='basis_h2o',water_reference=dict(concentration_convention='assume_table_values_nominal_mM',
        basis_amplitude_convention='assume_common_scale_protons_encoded',frequency_origin_ppm=4.68))


def test_real_legacy_metabolites_shim_noise_and_paired_export(basis,monkeypatch,tmp_path):
    path,raw=basis;cfg=settings(path)
    df=pd.DataFrame([dict(Tissue=t,Metabolite=m,Conc_mean=float(c),Conc_std=.1,T2=200.) for t in ['WM','GM','CSF'] for m,c in [('NAA',10),('Cr',8)]])
    monkeypatch.setattr(pd,'read_csv',lambda *a,**kw:df.copy())
    subject=anatomy()
    monkeypatch.setattr(mrs_phantom,'load_skeleton_as_subject',lambda *a:(subject,['Background','WM','GM','CSF']))
    monkeypatch.setattr(mrs_phantom.DigitalPhantom,'init_lipid_labels',lambda self:torch.ones(2,2,2))
    monkeypatch.setattr(SignalModel,'simulate_mms',lambda self:torch.ones(4,256,dtype=torch.complex64)*.001)
    monkeypatch.setattr(SignalModel,'simulate_lipids',lambda self,*a:torch.ones(256,dtype=torch.complex64)*.002)
    state=np.random.get_state();ts=torch.random.get_rng_state().clone()
    a=simulate_batch([(1,1,1)],cfg,seed=123,save_dir=tmp_path/'saved')[0]
    b=simulate_batch([(1,1,1)],cfg,seed=123)[0]
    for key in ['spectrum','water_reference_fid']:np.testing.assert_array_equal(a[key],b[key])
    ref=mrs_io.read_FID(str(tmp_path/'saved/voxel_0000_1_1_1/water_reference.nii.gz'))
    np.testing.assert_array_equal(ref.mrs().FID,a['water_reference_fid'])
    scalar=prepare(cfg)['normalization']['scalar']
    cw=sum(a['ground_truth']['tissue_fractions'][t]*v for t,v in {'WM':39.0165,'GM':47.2305,'CSF':55.5}.items())
    # 8 cells * 8 mm^3, with the real legacy exp(-i*2*pi*shim*t) and shifted time origin.
    expected=1000*cw*64*raw['H2O']/scalar*np.exp(-2j*np.pi*.5*np.arange(1,257)/2000)
    np.testing.assert_allclose(a['water_reference_fid'],expected,rtol=5e-6,atol=.01)
    cfg['water_reference']['noise_sigma']=1
    noisy=simulate_batch([(1,1,1)],cfg,seed=123)[0]
    np.testing.assert_array_equal(noisy['spectrum'],a['spectrum'])
    np.testing.assert_array_equal(noisy['clean_water_reference_fid'],a['clean_water_reference_fid'])
    assert not np.array_equal(noisy['water_reference_fid'],a['water_reference_fid'])
    converted=results_to_nifti_mrs([a]);converted[0]['water_reference'].save(str(tmp_path/'converted.nii.gz'))
    np.testing.assert_array_equal(mrs_io.read_FID(str(tmp_path/'converted.nii.gz')).mrs().FID,a['water_reference_fid'])
    assert torch.equal(ts,torch.random.get_rng_state())
    for x,y in zip(state,np.random.get_state()):np.testing.assert_array_equal(x,y)


def test_fixed_normalizer_and_validation(basis):
    cfg=settings(basis[0]);a=prepare(cfg)
    direct=mrs_io.read_basis(str(basis[0])).get_formatted_basis(2000,256)
    assert a['normalization']['scalar']==float(np.mean(abs(direct)))
    cfg['water_reference']['relaxation']={t:{'T1_s':1,'T2_s':.1} for t in ['WM','GM','CSF']}
    with pytest.raises(ValueError):prepare(cfg)
    cfg['water_reference']['basis_relaxation_assumption']='no_tissue_attenuation_encoded'
    b=prepare(cfg)
    assert b['water_factors']['WM']==pytest.approx((1-np.exp(-2))*np.exp(-.35))


def test_actual_basis_selection_independence(basis):
    from simulation.basis import Basis
    cfg=settings(basis[0]);meta=prepare(cfg)
    a=Basis(str(basis[0]),bw=2000,points=256,metabs=['NAA'],water_reference=meta)
    b=Basis(str(basis[0]),bw=2000,points=256,metabs=['NAA','Cr'],water_reference=meta)
    np.testing.assert_array_equal(a.water_fid,b.water_fid)
    np.testing.assert_array_equal(a.fids[:,a.names.index('NAA')],b.fids[:,b.names.index('NAA')])


def test_notebook_legacy_reference_imports_from_notebooks_directory():
    import json,os,subprocess,sys
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    notebook=json.loads((root/'notebooks/headless_simulation.ipynb').read_text())
    cell=next(c for c in notebook['cells'] if c['cell_type']=='code' and 'legacy_basis_h2o_demo' in c.get('metadata',{}).get('tags',[]))
    imports=''.join(cell['source']).split('legacy_reference_settings =')[0]
    result=subprocess.run([sys.executable,'-c',imports],cwd=root/'notebooks',
        env=dict(os.environ,PYTHONPATH=''),capture_output=True,text=True)
    assert result.returncode==0,result.stderr

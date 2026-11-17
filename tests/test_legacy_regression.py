"""Fixed-seed execution against exact pre-Phase-8 production files, without fitting."""
from copy import deepcopy
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from legacy_support import anatomy

from simulation import headless,mrs_phantom
from simulation.signalModel import SignalModel

C={'WM':{'NAA':10.,'Cr':8.},'GM':{'NAA':14.,'Cr':9.},'CSF':{'NAA':0.,'Cr':0.}}


def test_executed_legacy_default_and_explicit_mode_equal_pre_phase8(monkeypatch,basis):
    def source(path):
        return (Path(__file__).parent/'fixtures/legacy_pre_quantitative'/Path(path).name).read_text()
    old_phantom={'__name__':'phase8_baseline_phantom'};exec(compile(source('simulation/mrs_phantom.py'),'baseline_phantom.py','exec'),old_phantom)
    old_headless={'__name__':'phase8_baseline_headless'};exec(compile(source('simulation/headless.py'),'baseline_headless.py','exec'),old_headless)
    old_model={'__name__':'baseline_signal_model'}
    exec(compile(source('simulation/signalModel.py'),'baseline_signal_model.py','exec'),old_model)
    baseline_class=old_model['SignalModel']
    old_phantom['SignalModel']=baseline_class
    basis_path=basis[0]
    t=np.arange(64)/2000
    basis=SimpleNamespace(n=64,B0=3.,TE=30.,cf=127.731,dwelltime=.0005,t=t,
        f=np.fft.fftfreq(64,.0005),ppm=4.68+np.fft.fftshift(np.fft.fftfreq(64,.0005))/127.731,
        names=['NAA','Cr'],fids=np.column_stack([np.exp(-np.pi*t+2j*np.pi*f*t+.2j) for f in [-340,-210]]))
    df=pd.DataFrame([dict(Tissue=tissue,Metabolite=m,Conc_mean=c,Conc_std=.1,T2=200.)
        for tissue,profile in C.items() for m,c in profile.items() if m in ['NAA','Cr']])
    monkeypatch.setattr(pd,'read_csv',lambda *args,**kwargs:df.copy())
    def load(*args):return anatomy(),['Background','WM','GM','CSF']
    monkeypatch.setattr(mrs_phantom,'load_skeleton_as_subject',load);old_phantom['load_skeleton_as_subject']=load
    monkeypatch.setattr(mrs_phantom,'Basis',lambda *a,**k:deepcopy(basis));old_phantom['Basis']=lambda *a,**k:deepcopy(basis)
    old_headless['DigitalPhantom']=old_phantom['DigitalPhantom']
    for cls in [old_phantom['DigitalPhantom'],mrs_phantom.DigitalPhantom]:
        monkeypatch.setattr(cls,'init_lipid_labels',lambda self:torch.ones(2,2,2))
    for cls in [SignalModel,baseline_class]:
        monkeypatch.setattr(cls,'simulate_mms',lambda self:torch.tensor([0,.001,.002,0])[:,None].repeat(1,64).to(torch.complex64))
        monkeypatch.setattr(cls,'simulate_lipids',lambda self,*args:torch.tensor(np.exp(2j*np.pi*40*t))*.01)
    settings=dict(path2basis=str(basis_path),basis_set_settings={'TE':30.,'vendor':'test','localization':'PRESS'},
        voxel_size_mm=(4.,4.,4.),simulation_params=dict(bandwidth=2000,spectral_points=64,TR=2000.,lipid_sigma=1.,
            mm_level=.5,lipid_amp_factor=.1,water_amp_factor=.3,water_damping_min=1.,water_damping_max=4.,water_phase_min=-30.,water_phase_max=30.,noise_level=.001))
    npstate=np.random.get_state();torchstate=torch.random.get_rng_state().clone()
    reference=old_headless['simulate_batch']([(1,1,1)],deepcopy(settings),seed=123)[0]
    for mode in [None,'legacy']:
        config=deepcopy(settings)
        if mode:config['simulation_mode']=mode
        actual=headless.simulate_batch([(1,1,1)],config,seed=123)[0]
        torch.testing.assert_close(actual['spectrum'],reference['spectrum'],rtol=0,atol=0)
        for name in reference['components']:
            torch.testing.assert_close(actual['components'][name],reference['components'][name],rtol=0,atol=0)
            assert torch.linalg.vector_norm(actual['components'][name])>0
        for key in ['voxel_weighted','tissue_fractions']:assert actual['ground_truth'][key]==reference['ground_truth'][key]
        for tissue,values in reference['ground_truth']['sampled_per_tissue'].items():
            for name,value in values.items():assert actual['ground_truth']['sampled_per_tissue'][tissue][name]==value
        assert actual['ground_truth']['metabolite_mapping']['version']=='name-keyed-v1'
        metadata=deepcopy(actual['nifti_metadata']);metadata['sim_params'].pop('metabolite_mapping')
        assert metadata==reference['nifti_metadata']
        np.testing.assert_array_equal(actual['time_axis'],reference['time_axis']);np.testing.assert_array_equal(actual['ppm_axis'],reference['ppm_axis'])
    np.testing.assert_array_equal(torch.random.get_rng_state(),torchstate)
    after=np.random.get_state();np.testing.assert_array_equal(after[1],npstate[1]);assert after[2:]==npstate[2:]

"""Fixed-seed execution against exact pre-Phase-8 production files, without fitting."""
from copy import deepcopy
from types import SimpleNamespace
import ast
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch
import torchio as tio

from examples.quantitative_water_reference_production import BUNDLE,C
from simulation import headless,mrs_phantom
from simulation.signalModel import SignalModel
from simulation.quantitative import MODE,simulate_quantitative

BASELINE='c85821a833209b45384c0a8835e5d2ead8051daf'


def anatomy():
    affine=np.diag([2.,2.,2.,1.])
    labels=torch.tensor([0,1,1,2,2,2,3,3],dtype=torch.int16).reshape(1,2,2,2)
    return tio.Subject(labels=tio.LabelMap(tensor=labels,affine=affine),
        B0=tio.ScalarImage(tensor=torch.full((1,2,2,2),.5),affine=affine),
        brain_mask=tio.LabelMap(tensor=torch.ones(1,2,2,2),affine=affine))


def test_executed_legacy_default_and_explicit_mode_equal_pre_phase8(monkeypatch):
    def source(path):
        return (Path(__file__).parent/'fixtures/legacy_pre_quantitative'/Path(path).name).read_text()
    old_phantom={'__name__':'phase8_baseline_phantom'};exec(compile(source('simulation/mrs_phantom.py'),'baseline_phantom.py','exec'),old_phantom)
    old_headless={'__name__':'phase8_baseline_headless'};exec(compile(source('simulation/headless.py'),'baseline_headless.py','exec'),old_headless)
    # The actual signal model is unchanged; only external MM/lipid providers are replaced.
    class WithoutImports(ast.NodeTransformer):
        def visit_Import(self,node):return None
        def visit_ImportFrom(self,node):return None
    def equations(text):return ast.dump(WithoutImports().visit(ast.parse(text)))
    assert equations(source('simulation/signalModel.py'))==equations(Path('simulation/signalModel.py').read_text())
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
    monkeypatch.setattr(SignalModel,'simulate_mms',lambda self:torch.tensor([0,.001,.002,0])[:,None].repeat(1,64).to(torch.complex64))
    monkeypatch.setattr(SignalModel,'simulate_lipids',lambda self,*args:torch.tensor(np.exp(2j*np.pi*40*t))*.01)
    settings=dict(path2basis=str(BUNDLE/'metabolite_basis.BASIS'),basis_set_settings={'TE':30.,'vendor':'test','localization':'PRESS'},
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
        assert actual['ground_truth']==reference['ground_truth']
        assert actual['nifti_metadata']==reference['nifti_metadata']
        np.testing.assert_array_equal(actual['time_axis'],reference['time_axis']);np.testing.assert_array_equal(actual['ppm_axis'],reference['ppm_axis'])
    np.testing.assert_array_equal(torch.random.get_rng_state(),torchstate)
    after=np.random.get_state();np.testing.assert_array_equal(after[1],npstate[1]);assert after[2:]==npstate[2:]


def test_quantitative_digital_phantom_skips_legacy_initialization_and_roi(monkeypatch,tmp_path):
    from simulation.mrs_phantom import DigitalPhantom
    from simulation.mrscloud_quantitative_spatial import simulate_phantom_spatial
    def forbidden(*args,**kwargs):raise AssertionError('legacy basis/signal/lipid initialization called')
    monkeypatch.setattr(mrs_phantom,'load_skeleton_as_subject',lambda *args:(anatomy(),['Background','WM','GM','CSF']))
    monkeypatch.setattr(DigitalPhantom,'_create_basis_and_signal_model',forbidden)
    monkeypatch.setattr(DigitalPhantom,'init_lipid_labels',forbidden)
    q=dict(calibration_bundle=str(BUNDLE),metabolite_concentrations_mM=C,
        phantom_spatial={'use_phantom_b0':True,'b0_units':'Hz','b0_convention':'legacy_shim'})
    phantom=DigitalPhantom('test',None,None,{'simulation_mode':MODE,'quantitative':q},verbose=False)
    # Select a pure GM cell, excluding background.
    bounds=(1,2,0,1,0,1);pair=phantom.simulate_data(bounds)
    assert pair['ground_truth']['tissue_fractions']=={'WM':0.,'GM':1.,'CSF':0.}
    assert pair['ground_truth']['spatial_integration']['b0_offsets_hz']==[-.5]
    from simulation.quantitative import DEFAULT_WATER
    reference=simulate_phantom_spatial(phantom,bounds,calibration_bundle=BUNDLE,metabolite_concentrations_mM=C,
        water_concentrations_mol_l=DEFAULT_WATER,use_phantom_b0=True,b0_units='Hz',b0_convention='legacy_shim')
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(pair[key],reference[key])
    np.testing.assert_array_equal(pair['simulation_config']['voxel_affine_mm'],[[2,0,0,2],[0,2,0,0],[0,0,2,0],[0,0,0,1]])
    from simulation.quantitative import export_quantitative
    from fsl_mrs.utils import mrs_io
    import nibabel as nib
    paths=export_quantitative(pair,tmp_path/'anatomical')
    for name,key in [('metabolites','metabolite_fid'),('water_reference','water_reference_fid')]:
        np.testing.assert_array_equal(mrs_io.read_FID(paths[name]).mrs().FID,pair[key])
        np.testing.assert_array_equal(nib.load(paths[name]).affine,pair['simulation_config']['voxel_affine_mm'])
    result=headless.simulate_batch([(1,0,0)],{'simulation_mode':MODE,'quantitative':q,'voxel_size_mm':[2,2,2]},phantom=phantom)[0]
    np.testing.assert_array_equal(result['water_reference_fid'],pair['water_reference_fid'])
    with pytest.raises(ValueError):phantom.update_phantom({})

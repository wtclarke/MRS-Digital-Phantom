"""Production signals/exports only: explicitly forbid fitting."""
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import torch
from fsl_mrs.utils import mrs_io

from examples.quantitative_water_reference_production import BUNDLE,C,relaxation_config,validate
from simulation.headless import simulate_batch,run_headless,results_to_nifti_mrs
from simulation.mrscloud_full_basis import load_library,sha256
from simulation.mrscloud_quantitative_voxel import simulate_voxel
from simulation.mrscloud_quantitative_spatial import simulate_spatial
from simulation.quantitative import MODE,DEFAULT_WATER,simulate_quantitative,export_quantitative


@pytest.fixture(autouse=True)
def no_fitting(monkeypatch):
    from fsl_mrs.utils import fitting,quantify
    def forbidden(*args,**kwargs):raise AssertionError('Phase 8 must not run fitting')
    monkeypatch.setattr(fitting,'fit_FSLModel',forbidden)
    monkeypatch.setattr(quantify.WaterRef,'_fit_w_ref',forbidden)


def configuration(**kwargs):
    return dict(calibration_bundle=str(BUNDLE),metabolite_concentrations_mM=deepcopy(C),**kwargs)


@pytest.fixture(scope='module')
def units():return load_library(BUNDLE,selected_names=list(C['WM']))


@pytest.mark.parametrize('fractions',[{'WM':1},{'GM':1},{'CSF':1},{'WM':.2,'GM':.3,'CSF':.5}])
@pytest.mark.parametrize('relaxed',[False,True])
def test_independent_tissue_signals_and_phase6_regression(units,fractions,relaxed):
    responses,water,_=units;relax=relaxation_config() if relaxed else None
    pair=simulate_quantitative(configuration(tissue_fractions=fractions,relaxation=relax))
    em=np.zeros(2048,complex);ew=np.zeros(2048,complex)
    for tissue in C:
        f=fractions.get(tissue,0.)
        def R(p):return (1-np.exp(-2/p['T1_s']))*np.exp(-.03/p['T2_s'])
        for m,c in C[tissue].items():em+=f*.001*c*(R(relax['metabolites'][tissue][m]) if relax else 1)*responses[m]
        ew+=f*DEFAULT_WATER[tissue]*(R(relax['water'][tissue]) if relax else 1)*water
    np.testing.assert_allclose(pair['metabolite_fid'],em,rtol=3e-15,atol=1e-20)
    np.testing.assert_allclose(pair['water_reference_fid'],ew,rtol=3e-15,atol=1e-20)
    reference=simulate_voxel(fractions,calibration_bundle=BUNDLE,metabolite_concentrations_mM=C,water_concentrations_mol_l=DEFAULT_WATER,relaxation=relax)
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(pair[key],reference[key])


def test_spatial_b0_relaxation_independent_cells(units):
    responses,water,_=units;relax=relaxation_config();weights=[.2,.3,.5];offsets=[-3,7,1];labels=['WM','GM','CSF']
    cfg=configuration(spatial={'cell_labels':labels,'cell_volumes_mm3':[2,3,5],'b0_offsets_hz':offsets,'b0_units':'Hz','b0_convention':'positive_frequency'},relaxation=relax)
    pair=simulate_quantitative(cfg)
    reference=simulate_spatial(labels,calibration_bundle=BUNDLE,metabolite_concentrations_mM=C,water_concentrations_mol_l=DEFAULT_WATER,cell_volumes_mm3=[2,3,5],b0_offsets_hz=offsets,relaxation=relax)
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(pair[key],reference[key])
    t=np.arange(2048)/2000;em=np.zeros(2048,complex);ew=em.copy()
    for tissue,f,df in zip(labels,weights,offsets):
        def R(p):return (1-np.exp(-2/p['T1_s']))*np.exp(-.03/p['T2_s'])
        for m,c in C[tissue].items():em+=f*.001*c*R(relax['metabolites'][tissue][m])*responses[m]*np.exp(2j*np.pi*df*t)
        ew+=f*DEFAULT_WATER[tissue]*R(relax['water'][tissue])*water*np.exp(2j*np.pi*df*t)
    np.testing.assert_allclose(pair['metabolite_fid'],em,rtol=3e-15,atol=1e-20)
    np.testing.assert_allclose(pair['water_reference_fid'],ew,rtol=3e-15,atol=1e-20)


def test_independent_linewidth_and_water_selection_linearity():
    cfg=configuration(tissue_fractions={'WM':1});base=simulate_quantitative(cfg)
    broad=simulate_quantitative(dict(cfg,additional_linewidth_hz={'metabolites':2.,'water':7.}))
    t=np.arange(2048)/2000
    np.testing.assert_array_equal(broad['metabolite_fid'],sum(base['components'][m]*np.exp(-2*np.pi*t) for m in base['components']))
    np.testing.assert_array_equal(broad['water_reference_fid'],base['water_reference_fid']*np.exp(-7*np.pi*t))
    assert broad['metabolite_fid'][0]==base['metabolite_fid'][0] and broad['water_reference_fid'][0]==base['water_reference_fid'][0]
    subset={t:{'NAA':p['NAA']} for t,p in C.items()}
    selected=simulate_quantitative(dict(cfg,metabolite_concentrations_mM=subset))
    np.testing.assert_array_equal(selected['water_reference_fid'],base['water_reference_fid'])
    twice=simulate_quantitative(dict(cfg,water_concentrations_mol_l={t:2*c for t,c in DEFAULT_WATER.items()}))
    np.testing.assert_array_equal(twice['water_reference_fid'],2*base['water_reference_fid'])
    gain=simulate_quantitative(dict(cfg,receiver_factor=2))
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(gain[key],2*base[key])
    np.testing.assert_array_equal(gain['metabolite_fid']/gain['water_reference_fid'],base['metabolite_fid']/base['water_reference_fid'])
    doubled=simulate_quantitative(dict(cfg,metabolite_concentrations_mM={t:{m:2*c for m,c in p.items()} for t,p in C.items()}))
    np.testing.assert_array_equal(doubled['metabolite_fid'],2*base['metabolite_fid'])


def test_noise_phase7a_rng_and_clean_export(tmp_path):
    npstate=np.random.get_state();torchstate=torch.random.get_rng_state().clone()
    cfg=configuration(tissue_fractions={'WM':1},noise={'seed':70100,'metabolite_snr':50,'water_snr':1000})
    pair=simulate_quantitative(cfg);repeat=simulate_quantitative(cfg)
    reference=simulate_voxel({'WM':1},calibration_bundle=BUNDLE,metabolite_concentrations_mM=C,water_concentrations_mol_l=DEFAULT_WATER,noise=cfg['noise'])
    for key in ['metabolite_fid','water_reference_fid','clean_metabolite_fid','clean_water_reference_fid']:
        np.testing.assert_array_equal(pair[key],repeat[key]);np.testing.assert_array_equal(pair[key],reference[key])
    saved=export_quantitative(pair,tmp_path/'pair')
    for name,key in [('metabolites','metabolite_fid'),('water_reference','water_reference_fid'),('metabolites_clean','clean_metabolite_fid'),('water_reference_clean','clean_water_reference_fid')]:
        reader=mrs_io.read_FID(saved[name]);np.testing.assert_array_equal(reader.mrs().FID,pair[key])
        assert reader.dwelltime==.0005 and reader.hdr_ext['SpectrometerFrequency']==[127.731]
    np.testing.assert_array_equal(torch.random.get_rng_state(),torchstate)
    after=np.random.get_state();assert after[0]==npstate[0] and after[2:]==npstate[2:];np.testing.assert_array_equal(after[1],npstate[1])
    with pytest.raises(FileExistsError):export_quantitative(pair,tmp_path/'pair')


@pytest.mark.parametrize('extra',[
    {'water_amp_factor':5},{'expected_calibration_identity':'wrong'}, {'points':4096},{'bandwidth_hz':1000},
    {'concentration_convention':'per_tissue_water_volume'}, {'noise':{'seed':True}},
    {'additional_linewidth_hz':{'water':-1,'metabolites':0}},
    {'spatial':{'cell_labels':['WM'],'b0_offsets_hz':[37]}},
    {'relaxation':{'mode':'none','TE_s':.06,'TR_s':2}},
])
def test_incompatible_configuration_rejected(extra):
    cfg=configuration(tissue_fractions={'WM':1})
    if 'spatial' in extra:cfg.pop('tissue_fractions')
    with pytest.raises((ValueError,TypeError)):simulate_quantitative(dict(cfg,**extra))


def test_headless_opt_in_batch_and_cli(tmp_path):
    cfg={'simulation_mode':MODE,'quantitative':configuration(tissue_fractions={'WM':1},noise={'seed':10,'metabolite_sigma':1e-8,'water_sigma':1e-6})}
    result=simulate_batch([(0,0,0),(1,0,0)],cfg,save_dir=tmp_path/'batch')
    assert [r['simulation_config']['noise_parameters']['root_seed'] for r in result]==[10,11]
    for r in result:
        assert Path(r['saved_path'],'calibration_provenance.json').is_file()
    with pytest.raises(ValueError):results_to_nifti_mrs(result)
    with pytest.raises(ValueError):simulate_batch([(0,0,0)],cfg,total_components=['Water'])
    with pytest.raises(ValueError):simulate_batch([],{'simulation_mode':'typo'})
    cli=dict(cfg,output_dir=str(tmp_path/'cli'));path=tmp_path/'cli.json';path.write_text(json.dumps(cli))
    completed=subprocess.run([sys.executable,'-m','simulation.headless','--config',str(path)],capture_output=True,text=True)
    assert completed.returncode==0,completed.stderr
    assert (tmp_path/'cli/voxel_0000_0_0_0/water_reference.nii.gz').is_file()


def test_example_bundle_immutability(tmp_path):
    before={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}
    result=validate(tmp_path/'examples')
    assert result['bundle_unchanged'] and result['fitting']=='not performed'
    for case in result['cases']:
        assert case['independent_metabolite_error']<1e-19 and case['independent_water_error']<1e-17
        assert case['export_max_complex_error']==0
    assert before=={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}



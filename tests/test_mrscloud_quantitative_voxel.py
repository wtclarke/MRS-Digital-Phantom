from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from fsl_mrs.utils import mrs_io

from simulation.mrscloud_full_basis import load_library
from simulation.mrscloud_quantitative_voxel import simulate_voxel, simulate_phantom_voxel
from simulation.mrscloud_quantitative_pair import export_pair

BUNDLE=Path('outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative')
TISSUES=['WM','GM','CSF']
C={'WM':{'NAA':10.,'Cr':8.},'GM':{'NAA':14.,'Cr':9.},'CSF':{'NAA':0.,'Cr':0.}}
W={'WM':39.0165,'GM':47.2305,'CSF':55.5}


def simulate(fractions,**kwargs):
    args=dict(calibration_bundle=BUNDLE,metabolite_concentrations_mM=C,water_concentrations_mol_l=W)
    args.update(kwargs)
    return simulate_voxel(fractions,**args)


@pytest.fixture(scope='module')
def unit():
    return load_library(BUNDLE,selected_names=['NAA','Cr'])


@pytest.mark.parametrize('fractions,expected_c,expected_w',[
    ({'WM':1},[10.,8.],39.0165),({'GM':1},[14.,9.],47.2305),
    ({'CSF':1},[0.,0.],55.5),({'WM':.2,'GM':.3,'CSF':.5},[6.2,4.3],49.72245)])
def test_pure_and_mixed_independent_full_complex_fids(unit,fractions,expected_c,expected_w):
    responses,water,_=unit
    pair=simulate(fractions)
    expected=.001*expected_c[0]*responses['NAA']+.001*expected_c[1]*responses['Cr']
    np.testing.assert_allclose(pair['metabolite_fid'],expected,rtol=3e-15,atol=1e-20)
    np.testing.assert_allclose(pair['water_reference_fid'],expected_w*water,rtol=3e-15,atol=1e-20)
    for name,value in zip(['NAA','Cr'],expected_c):
        assert pair['ground_truth']['voxel_average_metabolite_concentrations_mM'][name]==pytest.approx(value)
    assert pair['ground_truth']['effective_voxel_water_concentration_mol_l']==pytest.approx(expected_w)
    assert pair['time_s'][0]==0 and pair['time_s'][1]==.0005


def test_fraction_and_tissue_only_conventions():
    pair=simulate({'WM':.2,'GM':.3,'CSF':.5})
    gt=pair['ground_truth']
    assert gt['tissue_only_metabolite_concentrations_mM']==pytest.approx({'NAA':12.4,'Cr':8.6})
    csf=simulate({'CSF':1})
    assert csf['ground_truth']['tissue_only_metabolite_concentrations_mM']=={'NAA':None,'Cr':None}
    nonzero=deepcopy(C);nonzero['CSF']['NAA']=2
    extra=simulate({'WM':.2,'GM':.3,'CSF':.5},metabolite_concentrations_mM=nonzero)
    assert extra['ground_truth']['voxel_average_metabolite_concentrations_mM']['NAA']==pytest.approx(7.2)
    assert extra['ground_truth']['tissue_only_metabolite_concentrations_mM']['NAA']==pytest.approx(12.4)
    normalized=simulate({'WM':2,'GM':3,'CSF':5},normalize_fractions=True)
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(normalized[key],pair[key])
    assert normalized['ground_truth']['fraction_normalization']['applied_divisor']==10
    assert normalized['ground_truth']['input_tissue_fractions']=={'WM':2,'GM':3,'CSF':5}


@pytest.mark.parametrize('fractions',[{'WM':2,'GM':3,'CSF':5},{'WM':0},{'WM':-.1,'GM':1.1},
                                    {'WM':np.nan},{'WM':np.inf},{'Fat':1}])
def test_invalid_fractions_are_rejected(fractions):
    with pytest.raises(ValueError):simulate(fractions)


@pytest.mark.parametrize('change',['missing_tissue','missing_water','different_metabolites','negative_metabolite',
                                  'negative_water','nan_zero_fraction','unknown_metabolite','control','receiver'])
def test_invalid_profiles_do_not_hide_behind_tissue_weighting(change):
    c=deepcopy(C);w=dict(W);args={}
    if change=='missing_tissue':del c['CSF']
    elif change=='missing_water':del w['GM']
    elif change=='different_metabolites':del c['GM']['Cr']
    elif change=='negative_metabolite':c['WM']['NAA']=-1
    elif change=='negative_water':w['WM']=-1
    elif change=='nan_zero_fraction':c['CSF']['NAA']=np.nan
    elif change=='receiver':args['receiver_factor']=-1
    else:
        for t in TISSUES:c[t]['Cho' if change=='unknown_metabolite' else 'H2O']=0
    with pytest.raises(ValueError):simulate({'WM':1},metabolite_concentrations_mM=c,water_concentrations_mol_l=w,**args)


def test_linearity_gain_ratios_selection_and_determinism(unit):
    fractions={'WM':.2,'GM':.3,'CSF':.5}
    a=simulate(fractions);repeat=simulate(fractions)
    doubled={t:{m:2*v for m,v in values.items()} for t,values in C.items()}
    m=simulate(fractions,metabolite_concentrations_mM=doubled)
    w=simulate(fractions,water_concentrations_mol_l={t:2*v for t,v in W.items()})
    g=simulate(fractions,receiver_factor=2)
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_array_equal(a[key],repeat[key]);np.testing.assert_array_equal(g[key],2*a[key])
    np.testing.assert_array_equal(m['metabolite_fid'],2*a['metabolite_fid'])
    np.testing.assert_array_equal(m['water_reference_fid'],a['water_reference_fid'])
    np.testing.assert_array_equal(w['water_reference_fid'],2*a['water_reference_fid'])
    np.testing.assert_array_equal(w['metabolite_fid'],a['metabolite_fid'])
    np.testing.assert_array_equal(a['metabolite_fid']/a['water_reference_fid'],g['metabolite_fid']/g['water_reference_fid'])
    alone=simulate(fractions,metabolite_concentrations_mM={t:{'NAA':C[t]['NAA']} for t in TISSUES})
    np.testing.assert_array_equal(alone['components']['NAA'],a['components']['NAA'])
    assert repeat['ground_truth']==a['ground_truth']
    assert repeat['simulation_config']==a['simulation_config']
    valid=simulate(fractions,expected_calibration=unit[2])
    np.testing.assert_array_equal(valid['metabolite_fid'],a['metabolite_fid'])
    for kwargs in [dict(points=1024),dict(bandwidth_hz=4000),dict(expected_calibration={})]:
        with pytest.raises(ValueError):simulate(fractions,**kwargs)


def test_complex_export_phase_frequency_and_json_sidecars(unit,tmp_path):
    import nibabel as nib
    pair=simulate({'WM':.2,'GM':.3,'CSF':.5});paths=export_pair(pair,tmp_path)
    for key,path in paths.items():
        obj=mrs_io.read_FID(str(path))
        np.testing.assert_array_equal(obj.mrs().FID,pair[key])
        # Validated gen_nifti_mrs convention: stored conjugate, interpreted FID.
        np.testing.assert_array_equal(np.asarray(nib.load(path).dataobj).squeeze(),pair[key].conj())
        assert obj.dwelltime==.0005 and obj.spectrometer_frequency==[127.731] and obj.SpecFreqChemShift==4.68
    gt=json.loads((tmp_path/'ground_truth.json').read_text())
    assert gt==pair['ground_truth']
    assert gt['shared_normalization_scalar']==8857.353925964408
    assert gt['calibration_identity']==unit[2]['calibration_identity']
    assert pair['simulation_config']['residual_water'] is False
    np.testing.assert_array_equal(pair['water_reference_fid'],49.72245*unit[1])


def phantom_fixture():
    import torch
    import torchio as tio
    from simulation.mrs_phantom import DigitalPhantom
    labels=torch.zeros((1,3,2,2),dtype=torch.int16)
    labels[:,1:]=torch.tensor([1,1,1,1,2,2,2,3]).reshape(1,2,2,2)
    affine=np.diag([2.,3.,4.,1.])
    phantom=DigitalPhantom.__new__(DigitalPhantom)
    phantom.subject=tio.Subject(labels=tio.LabelMap(tensor=labels,affine=affine))
    phantom.label_names=['Background','WM','GM','CSF'];phantom.spacing=(2.,3.,4.);phantom.affine=affine
    return phantom


def test_actual_phantom_label_getter_and_legacy_one_hot_mean_agree(monkeypatch):
    phantom=phantom_fixture()
    expected=phantom.get_label_data(one_hot=True)[1:3].reshape(-1,4).float().mean(dim=0).numpy()
    def forbidden(*args,**kwargs):raise AssertionError('legacy simulation/shim must not run')
    monkeypatch.setattr(phantom,'simulate_data',forbidden)
    monkeypatch.setattr(phantom,'get_B0_data',forbidden)
    pair=simulate_phantom_voxel(phantom,(1,3,0,2,0,2),calibration_bundle=BUNDLE,
        metabolite_concentrations_mM=C,water_concentrations_mol_l=W)
    f=pair['ground_truth']['tissue_fractions']
    np.testing.assert_array_equal(list(f.values()),expected[1:])
    assert f=={'WM':.5,'GM':.375,'CSF':.125}
    assert pair['ground_truth']['voxel_geometry']['voxel_volume_mm3']==192
    original=pair['metabolite_fid'].copy()
    phantom.spacing=(1.,1.,1.)
    phantom.affine=np.eye(4)
    resized=simulate_phantom_voxel(phantom,(1,3,0,2,0,2),calibration_bundle=BUNDLE,
        metabolite_concentrations_mM=C,water_concentrations_mol_l=W)
    np.testing.assert_array_equal(resized['metabolite_fid'],original)
    assert resized['ground_truth']['voxel_geometry']['voxel_volume_mm3']==8


@pytest.mark.parametrize('bounds',[(0,3,0,2,0,2),(-1,2,0,2,0,2),(1,4,0,2,0,2),(1,1,0,2,0,2),(1,2,3)])
def test_background_and_bad_voxel_geometry_are_not_silently_renormalized(bounds):
    with pytest.raises(ValueError):simulate_phantom_voxel(phantom_fixture(),bounds,calibration_bundle=BUNDLE,
        metabolite_concentrations_mM=C,water_concentrations_mol_l=W)


def test_legacy_residual_water_and_rng_are_unchanged():
    import torch
    from simulation.signalModel import SignalModel
    model=SignalModel(None,SimpleNamespace(B0=3,TE=30,t=np.arange(32)/2000,cf=127.731),
        labels=['WM'],simulation_params=dict(water_amp_factor=5,water_damping_min=1,
            water_damping_max=8,water_phase_min=-180,water_phase_max=180),verbose=False)
    state=np.random.get_state();torch_state=torch.random.get_rng_state().clone()
    try:
        np.random.seed(61);baseline=model.simulate_water().numpy().copy();next_values=np.random.random(8)
        np.random.seed(61);simulate({'WM':.2,'GM':.3,'CSF':.5})
        np.testing.assert_array_equal(model.simulate_water().numpy(),baseline)
        np.testing.assert_array_equal(np.random.random(8),next_values)
        assert torch.equal(torch.random.get_rng_state(),torch_state)
    finally:np.random.set_state(state)

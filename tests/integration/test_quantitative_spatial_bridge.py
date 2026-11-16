"""Actual phantom getters, geometry, B0 declarations and executed legacy sign."""
from types import SimpleNamespace

import numpy as np
import pytest
import torch
import torchio as tio

from examples.mrscloud_quantitative_voxel import BUNDLE, METABOLITES, WATER
from simulation.mrscloud_quantitative_spatial import simulate_phantom_spatial, simulate_spatial


def phantom_fixture():
    from simulation.mrs_phantom import DigitalPhantom
    labels=torch.zeros((1,3,2,2),dtype=torch.int16)
    labels[:,1:]=torch.tensor([1,1,1,1,2,2,2,3]).reshape(1,2,2,2)
    affine=np.diag([2.,3.,4.,1.])
    phantom=DigitalPhantom.__new__(DigitalPhantom)
    phantom.subject=tio.Subject(labels=tio.LabelMap(tensor=labels,affine=affine),
        B0=tio.ScalarImage(tensor=torch.arange(12,dtype=torch.float64).reshape(1,3,2,2),affine=affine))
    phantom.label_names=['Background','WM','GM','CSF'];phantom.spacing=(2.,3.,4.);phantom.affine=affine
    return phantom


def bridge(phantom,**kwargs):
    return simulate_phantom_spatial(phantom,(1,3,0,2,0,2),calibration_bundle=BUNDLE,
        metabolite_concentrations_mM=METABOLITES,water_concentrations_mol_l=WATER,**kwargs)


def test_phantom_geometry_fractions_and_explicit_b0_sign(monkeypatch):
    phantom=phantom_fixture()
    def forbidden(*args,**kwargs):raise AssertionError('legacy simulation must not run')
    monkeypatch.setattr(phantom,'simulate_data',forbidden)
    original_b0=phantom.get_B0_data().copy()
    pair=bridge(phantom,use_phantom_b0=True,b0_units='Hz',b0_convention='legacy_shim')
    spatial=pair['ground_truth']['spatial_integration']
    assert pair['ground_truth']['tissue_fractions']=={'WM':.5,'GM':.375,'CSF':.125}
    assert spatial['normalized_weights']==[.125]*8 and spatial['b0_sign_conversion']==-1
    assert spatial['b0_offsets_hz']==list(-np.arange(4.,12.))
    np.testing.assert_array_equal(phantom.get_B0_data(),original_b0)
    assert pair['ground_truth']['voxel_geometry']['voxel_volume_mm3']==pytest.approx(192)
    labels=np.array(['WM']*4+['GM']*3+['CSF']).reshape(2,2,2)
    direct=simulate_spatial(labels,b0_offsets_hz=-np.arange(4.,12.).reshape(2,2,2),
        calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,water_concentrations_mol_l=WATER)
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_array_equal(pair[key],direct[key])
    monkeypatch.setattr(phantom,'get_B0_data',forbidden)
    no_map=bridge(phantom)
    assert no_map['ground_truth']['spatial_integration']['b0_offsets_hz']==[0.]*8


@pytest.mark.parametrize('change',['units','sign','affine','shape','spacing','background','label_affine'])
def test_unknown_units_sign_or_inconsistent_geometry_rejected(change):
    phantom=phantom_fixture();kwargs=dict(use_phantom_b0=True,b0_units='Hz',b0_convention='positive_frequency')
    if change=='units':kwargs['b0_units']='ppm'
    elif change=='sign':kwargs['b0_convention']=None
    elif change=='affine':phantom.subject['B0'].affine[0,3]=1.
    elif change=='shape':phantom.subject['B0']=tio.ScalarImage(tensor=torch.zeros(1,2,2,2),affine=phantom.affine)
    elif change=='spacing':phantom.spacing=(1.,1.,1.)
    elif change=='background':phantom.subject['labels'].data[0,1,0,0]=0
    elif change=='label_affine':phantom.subject['labels'].affine[0,3]=1.
    with pytest.raises(ValueError):bridge(phantom,**kwargs)


def test_executed_legacy_sign_is_negative_new_sign_is_positive():
    from simulation.signalModel import SignalModel
    from utils.auxillary import fid2spec, spec2fid
    t=np.arange(256)/2000
    fid=torch.from_numpy(np.exp(-np.pi*t)*np.exp(1j*(.3+2*np.pi*7*t)))
    model=SignalModel.__new__(SignalModel);model.basis=SimpleNamespace(t=torch.from_numpy(t))
    result=model.apply_shim_map(fid2spec(fid.reshape(1,1,1,-1)),torch.tensor([[[37.]]]))
    interpreted=spec2fid(result).numpy().reshape(-1)
    # Legacy multiplies 2*pi by its float32 shim map before the float64 time axis.
    np.testing.assert_allclose(interpreted,fid.numpy()*np.exp(-2j*np.pi*37*t),rtol=5e-7,atol=2e-15)
    ratio=interpreted/fid.numpy()
    assert np.angle(ratio[1]/ratio[0])/(2*np.pi*.0005)==pytest.approx(-37,abs=1e-5)


def test_spatial_example_exports_independent_controls(tmp_path):
    from examples.mrscloud_quantitative_spatial import validate
    result=validate(tmp_path/'spatial_controls')
    assert result['bundle_unchanged']
    for row in result['cases']:
        assert row['independent_metabolite_max_error']<1e-19
        assert row['independent_water_max_error']<1e-17
        assert row['export_max_error']==row['repeat_max_error']==0
        assert row['initial_amplitude_max_difference']<1e-17
    assert result['cases'][0]['zero_reference_max_difference']<1e-17
    assert result['cases'][2]['water_magnitude_at_50ms']<1e-16

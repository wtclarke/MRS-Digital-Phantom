"""Independent cell sums, frequency-sign controls and zero-field regressions."""
from copy import deepcopy
import json

import nibabel as nib
import numpy as np
import pytest
from fsl_mrs.utils import mrs_io

from examples.mrscloud_quantitative_voxel import BUNDLE, METABOLITES, WATER
from examples.mrscloud_tissue_relaxation_validation import configuration
from simulation.mrscloud_full_basis import load_library
from simulation.mrscloud_quantitative_pair import export_pair
from simulation.mrscloud_quantitative_spatial import simulate_spatial
from simulation.mrscloud_quantitative_voxel import simulate_voxel


def spatial(labels, **kwargs):
    args=dict(calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,water_concentrations_mol_l=WATER)
    args.update(kwargs)
    return simulate_spatial(labels,**args)


def voxel(fractions,**kwargs):
    return simulate_voxel(fractions,calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,
        water_concentrations_mol_l=WATER,**kwargs)


@pytest.mark.parametrize('labels',[['WM']*5,['GM']*3,['CSF']*4,['WM']*2+['GM']*3+['CSF']*5])
@pytest.mark.parametrize('relaxed',[False,True])
def test_zero_b0_matches_volume_average_and_relaxation(labels,relaxed):
    cfg=configuration() if relaxed else None
    fractions={t:labels.count(t)/len(labels) for t in WATER}
    actual=spatial(labels,relaxation=cfg);expected=voxel(fractions,relaxation=cfg)
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_allclose(actual[key],expected[key],rtol=3e-15,atol=1e-20)
    for key in ['voxel_average_metabolite_concentrations_mM','tissue_only_metabolite_concentrations_mM']:
        for n,v in actual['ground_truth'][key].items():
            if v is None:assert expected['ground_truth'][key][n] is None
            else:assert v==pytest.approx(expected['ground_truth'][key][n],rel=5e-16)


def test_constant_positive_frequency_shift_and_symmetric_dephasing():
    baseline=spatial(['WM','WM'])
    shifted=spatial(['WM','WM'],b0_offsets_hz=[37.,37.])
    split=spatial(['WM','WM'],b0_offsets_hz=[-5.,5.])
    t=np.arange(2048)/2000
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_allclose(shifted[key],baseline[key]*np.exp(2j*np.pi*37*t),rtol=2e-15,atol=1e-20)
        np.testing.assert_allclose(split[key],baseline[key]*np.cos(2*np.pi*5*t),rtol=2e-15,atol=1e-20)
        assert shifted[key][0]==baseline[key][0] and split[key][0]==baseline[key][0]
    ratio=shifted['water_reference_fid']/baseline['water_reference_fid']
    assert np.angle(ratio[1]/ratio[0])/(2*np.pi*.0005)==pytest.approx(37,abs=1e-12)
    assert abs(split['water_reference_fid'][100])<1e-16 # first cosine null at 50 ms


def test_weights_volumes_and_tissue_correlated_field_independent_cell_sum():
    labels=['WM','GM','CSF'];weights=[.2,.3,.5];offsets=[-3.,7.,1.]
    cfg=configuration()
    pair=spatial(labels,spatial_weights=weights,b0_offsets_hz=offsets,relaxation=cfg)
    responses,water,_=load_library(BUNDLE,selected_names=list(METABOLITES['WM']))
    t=np.arange(2048)/2000;em=np.zeros(2048,complex);ew=np.zeros(2048,complex)
    for tissue,w,f in zip(labels,weights,offsets):
        modulation=np.exp(2j*np.pi*f*t)
        for n,c in METABOLITES[tissue].items():
            p=cfg['metabolites'][tissue][n]
            attenuation=(1-np.exp(-2/p['T1_s']))*np.exp(-.03/p['T2_s'])
            em+=w*c*.001*attenuation*responses[n]*modulation
        p=cfg['water'][tissue];attenuation=(1-np.exp(-2/p['T1_s']))*np.exp(-.03/p['T2_s'])
        ew+=w*WATER[tissue]*attenuation*water*modulation
    np.testing.assert_allclose(pair['metabolite_fid'],em,rtol=3e-15,atol=1e-20)
    np.testing.assert_allclose(pair['water_reference_fid'],ew,rtol=3e-15,atol=1e-20)
    volume=spatial(labels,cell_volumes_mm3=[2.,3.,5.],b0_offsets_hz=offsets,relaxation=cfg)
    bigger=spatial(labels,cell_volumes_mm3=[20.,30.,50.],b0_offsets_hz=offsets,relaxation=cfg)
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_array_equal(volume[key],pair[key]);np.testing.assert_array_equal(bigger[key],pair[key])
    assert volume['ground_truth']['spatial_integration']['total_anatomical_volume_mm3']==10
    # Shared cell modulation need not produce a shared whole-voxel envelope when tissue concentrations differ.
    zero=spatial(labels,spatial_weights=weights,relaxation=cfg)
    assert not np.allclose(pair['components']['NAA']/zero['components']['NAA'],pair['water_reference_fid']/zero['water_reference_fid'])


@pytest.mark.parametrize('kwargs',[
    {'spatial_weights':[1.,1.]},{'spatial_weights':[-.1,1.1]}, {'spatial_weights':[np.nan,0.]},
    {'spatial_weights':[0.,0.]},{'spatial_weights':[1.]},{'cell_volumes_mm3':[1.,0.]},
    {'cell_volumes_mm3':[1.,np.inf]}, {'cell_volumes_mm3':[1.,1.],'spatial_weights':[.5,.5]},
    {'b0_offsets_hz':[np.nan,0.]},{'b0_offsets_hz':[0.]},{'b0_offsets_hz':37.}])
def test_invalid_spatial_inputs_rejected(kwargs):
    with pytest.raises(ValueError):spatial(['WM','GM'],**kwargs)


@pytest.mark.parametrize('labels',[[],['WM','Background'],['Fat'],[['WM','GM']]])
def test_invalid_tissues_and_dimensions_rejected(labels):
    with pytest.raises(ValueError):spatial(labels)


def test_linearity_gain_selection_and_determinism():
    args=dict(spatial_weights=[.2,.3,.5],b0_offsets_hz=[-3.,7.,1.],relaxation=configuration())
    original=deepcopy(args);a=spatial(['WM','GM','CSF'],**args);b=spatial(['WM','GM','CSF'],**args)
    gain=spatial(['WM','GM','CSF'],receiver_factor=2,**args)
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_array_equal(a[key],b[key]);np.testing.assert_array_equal(gain[key],2*a[key])
    np.testing.assert_array_equal(a['metabolite_fid']/a['water_reference_fid'],gain['metabolite_fid']/gain['water_reference_fid'])
    assert a['ground_truth']==b['ground_truth'] and a['simulation_config']==b['simulation_config'] and args==original
    doubled={t:{n:2*c for n,c in p.items()} for t,p in METABOLITES.items()}
    m=spatial(['WM','GM','CSF'],metabolite_concentrations_mM=doubled,**args)
    w=spatial(['WM','GM','CSF'],water_concentrations_mol_l={t:2*c for t,c in WATER.items()},**args)
    np.testing.assert_array_equal(m['metabolite_fid'],2*a['metabolite_fid'])
    np.testing.assert_array_equal(w['water_reference_fid'],2*a['water_reference_fid'])
    selected=deepcopy(args);selected['relaxation']['metabolites']={t:{'NAA':p['NAA']} for t,p in args['relaxation']['metabolites'].items()}
    alone=spatial(['WM','GM','CSF'],metabolite_concentrations_mM={t:{'NAA':p['NAA']} for t,p in METABOLITES.items()},**selected)
    np.testing.assert_array_equal(alone['components']['NAA'],a['components']['NAA'])


def test_spatial_complex_export_and_sidecars(tmp_path):
    pair=spatial(np.array(['WM','GM','CSF','WM']).reshape(2,1,2),
        b0_offsets_hz=np.array([-3.,7.,1.,-3.]).reshape(2,1,2),relaxation=configuration())
    paths=export_pair(pair,tmp_path)
    repeated=export_pair(pair,tmp_path/'repeat')
    for key,path in paths.items():
        obj=mrs_io.read_FID(str(path))
        np.testing.assert_array_equal(obj.mrs().FID,pair[key])
        np.testing.assert_array_equal(mrs_io.read_FID(str(repeated[key])).mrs().FID,obj.mrs().FID)
        np.testing.assert_array_equal(np.asarray(nib.load(path).dataobj).squeeze(),pair[key].conj())
        assert obj.dwelltime==.0005 and obj.spectrometer_frequency==[127.731] and obj.SpecFreqChemShift==4.68
    for key in ['ground_truth','simulation_config']:
        assert json.loads((tmp_path/f'{key}.json').read_text())==pair[key]
        assert (tmp_path/f'{key}.json').read_bytes()==(tmp_path/'repeat'/f'{key}.json').read_bytes()

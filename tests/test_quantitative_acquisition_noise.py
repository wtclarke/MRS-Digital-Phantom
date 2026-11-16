"""Statistical noise controls and exact default-path regressions."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from examples.mrscloud_quantitative_voxel import BUNDLE,METABOLITES,WATER
from simulation.mrscloud_quantitative_noise import add_pair_noise,complex_gaussian_noise,reference_peak
from simulation.mrscloud_quantitative_pair import generate_pair,export_pair
from simulation.mrscloud_quantitative_spatial import simulate_spatial
from simulation.mrscloud_quantitative_voxel import simulate_voxel


def voxel(**kwargs):
    return simulate_voxel({'WM':1},calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,
                          water_concentrations_mol_l=WATER,**kwargs)


def test_complex_component_statistics_and_global_rng_isolation():
    import torch
    state=np.random.get_state();ts=torch.random.get_rng_state().clone()
    noise=complex_gaussian_noise(200000,sigma=2.,seed=8671)
    assert abs(noise.real.mean())<.02 and abs(noise.imag.mean())<.02
    assert noise.real.std()==pytest.approx(2,rel=.008) and noise.imag.std()==pytest.approx(2,rel=.008)
    assert abs(np.corrcoef(noise.real,noise.imag)[0,1])<.008
    assert np.mean(abs(noise)**2)==pytest.approx(8,rel=.008)
    for actual,old in zip(np.random.get_state(),state):np.testing.assert_array_equal(actual,old)
    assert torch.equal(ts,torch.random.get_rng_state())
    np.testing.assert_array_equal(noise,complex_gaussian_noise(200000,sigma=2.,seed=8671))


def test_snr_calculation_and_independent_streams():
    clean=voxel();original=deepcopy(clean)
    state=np.random.get_state()
    noisy=add_pair_noise(clean,seed=18,metabolite_snr=50,water_snr=1000)
    for actual,old in zip(np.random.get_state(),state):np.testing.assert_array_equal(actual,old)
    cfg=clean['simulation_config'];t=np.arange(2048)/2000
    ppm=np.fft.fftshift(np.fft.fftfreq(2048,.0005))/127.731+4.68
    # Direct DFT checks the chosen clean peak independently of reference_peak's FFT.
    spec=np.fft.fftshift(np.fft.fft(clean['components']['NAA'],norm='ortho'))
    indices=np.where((ppm>=1.8)&(ppm<=2.2))[0];peak_index=indices[np.argmax(abs(spec[indices]))]
    frequency=(ppm[peak_index]-4.68)*127.731
    peak=abs(np.sum(clean['components']['NAA']*np.exp(-2j*np.pi*frequency*t))/np.sqrt(2048))
    assert peak==pytest.approx(reference_peak(clean['components']['NAA'],cfg,(1.8,2.2)),rel=1e-12)
    meta=noisy['simulation_config']['noise_parameters']
    assert meta['metabolite']['sigma_per_component']==pytest.approx(peak/50,rel=1e-12)
    a=(noisy['metabolite_fid']-clean['metabolite_fid'])/meta['metabolite']['sigma_per_component']
    b=(noisy['water_reference_fid']-clean['water_reference_fid'])/meta['water']['sigma_per_component']
    assert abs(np.corrcoef(a.real,b.real)[0,1])<.06
    assert not np.array_equal(a,b)
    assert meta['metabolite']['realized_snr']==pytest.approx(50,rel=.06)
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_array_equal(noisy['clean_'+key],clean[key]);np.testing.assert_array_equal(clean[key],original[key])
    assert clean['simulation_config']==original['simulation_config'] and noisy['ground_truth']==clean['ground_truth']


def test_default_regression_integration_order_and_pair_modes(tmp_path):
    from fsl_mrs.utils import mrs_io
    clean=voxel();source=Path('outputs/mrscloud_quantitative_voxel/final/pure_WM')
    for key,filename in [('metabolite_fid','metabolites'),('water_reference_fid','water_reference')]:
        np.testing.assert_array_equal(clean[key],mrs_io.read_FID(str(source/f'{filename}.nii.gz')).mrs().FID)
    assert clean['simulation_config']==json.loads((source/'simulation_config.json').read_text())
    settings=dict(seed=7,metabolite_snr=25,water_snr=1000)
    direct=voxel(noise=settings);wrapped=add_pair_noise(clean,**settings)
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(direct[key],wrapped[key])
    spatial_args=dict(calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,water_concentrations_mol_l=WATER,b0_offsets_hz=[-3.,7.])
    from examples.mrscloud_tissue_relaxation_validation import configuration
    spatial_args['relaxation']=configuration()
    a=simulate_spatial(['WM','GM'],**spatial_args);b=simulate_spatial(['WM','GM'],noise=settings,**spatial_args)
    expected=add_pair_noise(a,**settings)
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(b[key],expected[key])
    for mode,args in [('analytical',{}),('mrscloud_sequence',dict(sequence_reference=Path('outputs/mrscloud_sequence_water_integration/validated/reference')))]:
        pair=generate_pair({'NAA':10},water_reference_mode=mode,**args)
        noisy=generate_pair({'NAA':10},water_reference_mode=mode,noise=settings,**args)
        np.testing.assert_array_equal(noisy['clean_metabolite_fid'],pair['metabolite_fid'])
    paths=export_pair(direct,tmp_path)
    assert len(paths)==4
    for key,path in paths.items():np.testing.assert_array_equal(mrs_io.read_FID(str(path)).mrs().FID,direct[key])
    off=add_pair_noise(clean,seed=7,metabolite_sigma=0,water_snr=np.inf)
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(off[key],clean[key])
    assert off['simulation_config']['noise'] is False


@pytest.mark.parametrize('kwargs',[dict(seed=-1),dict(seed=True),dict(seed=1,metabolite_sigma=-1),
    dict(seed=1,metabolite_snr=0),dict(seed=1,water_snr=np.nan),
    dict(seed=1,metabolite_sigma=1,metabolite_snr=10),dict(seed=1,water_sigma=np.inf)])
def test_invalid_noise_configurations(kwargs):
    with pytest.raises(ValueError):add_pair_noise(voxel(),**kwargs)

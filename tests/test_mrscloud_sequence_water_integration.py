from copy import deepcopy
import json
from pathlib import Path
import shutil

import numpy as np
import pytest
from fsl_mrs.utils import mrs_io

from examples.mrscloud_sequence_water_integration import build_reference_bundle
from simulation.mrscloud_quantitative_pair import generate_pair, export_pair, REFERENCE_BASIS
from simulation.mrscloud_sequence_reference import file_sha256
from simulation.frozen_basis_normalization import format_frozen_basis


@pytest.fixture(scope='module')
def reference(tmp_path_factory):
    return build_reference_bundle(tmp_path_factory.mktemp('sequence_reference'))


def pair(reference, **kwargs):
    return generate_pair({'NAA':10}, water_reference_mode='mrscloud_sequence', sequence_reference=reference, **kwargs)


def test_sequence_matches_full_phase5f_arrays(reference):
    generated=pair(reference)
    previous=np.load('outputs/mrscloud_water_processing/final/E_raw_acquisition_grid.npz')
    np.testing.assert_array_equal(generated['water_reference_fid'],40*previous['water'])
    np.testing.assert_array_equal(generated['metabolite_fid'],.010*previous['fixed_NAA'])
    metadata=generated['simulation_config']['synthetic_calibration']
    assert metadata['normalization_scalar']==2779.198340948381
    assert metadata['water_proton_multiplier_applied'] is False
    assert 'raw-water-common-max' in metadata['calibration_identity']
    assert generated['simulation_config']['water_linewidth_hz']==1
    matching=pair(reference,expected_metadata=metadata)
    np.testing.assert_array_equal(matching['water_reference_fid'],generated['water_reference_fid'])


def test_analytical_mode_preserves_saved_pre_feature_signals_and_metadata():
    generated=generate_pair({'NAA':10})
    explicit=generate_pair({'NAA':10},water_reference_mode='analytical')
    source=Path('outputs/fsl_mrs_water_scaling_phase5b/NAA_10mM')
    for key,filename in [('metabolite_fid','metabolites.nii.gz'),('water_reference_fid','water_reference.nii.gz')]:
        old=mrs_io.read_FID(str(source/filename)).mrs().FID
        np.testing.assert_array_equal(generated[key],old)
        np.testing.assert_array_equal(generated[key],explicit[key])
    old_metadata=json.loads((source/'simulation_config.json').read_text())
    for key,value in old_metadata.items():
        assert generated['simulation_config'][key]==value
    assert generated['simulation_config']['water_reference_mode']=='analytical'


def test_selection_linearity_determinism_and_metadata(reference):
    first=pair(reference)
    repeat=pair(reference)
    # Ref0ppm is an explicit reference-signal control, not a newly simulated Cr.
    multiple=generate_pair({'NAA':10,'Ref0ppm':0},water_reference_mode='mrscloud_sequence',sequence_reference=reference)
    doubled=generate_pair({'NAA':20},water_reference_mode='mrscloud_sequence',sequence_reference=reference,water_concentration_mol_l=80,receiver_factor=2)
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_array_equal(first[key],repeat[key])
        np.testing.assert_array_equal(first[key],multiple[key])
        np.testing.assert_array_equal(doubled[key],4*first[key])
    assert first['simulation_config']['synthetic_calibration']==multiple['simulation_config']['synthetic_calibration']
    assert first['simulation_config']==repeat['simulation_config']


@pytest.mark.parametrize('mode',['analytical','mrscloud_sequence'])
def test_both_modes_export_primary_and_optional_aligned_water(reference,tmp_path,mode):
    kwargs={'sequence_reference':reference} if mode=='mrscloud_sequence' else {'water_phase_rad':.6}
    generated=generate_pair({'NAA':10},water_reference_mode=mode,**kwargs)
    saved=generated['water_reference_fid'].copy()
    paths=export_pair(generated,tmp_path,aligned_water=True)
    for key in ['metabolite_fid','water_reference_fid']:
        loaded=mrs_io.read_FID(str(paths[key]))
        np.testing.assert_array_equal(loaded.mrs().FID,generated[key])
        assert loaded.dwelltime==.0005
        assert loaded.SpecFreqChemShift==4.68
    np.testing.assert_array_equal(generated['water_reference_fid'],saved)
    aligned=mrs_io.read_FID(str(paths['aligned_water_reference_fid'])).mrs().FID
    np.testing.assert_allclose(abs(aligned),abs(saved),rtol=1e-15,atol=1e-14)
    assert np.max(abs(np.angle(aligned)))<1e-12
    preprocessing=json.loads((tmp_path/'water_preprocessing.json').read_text())
    assert preprocessing['source_water_position_ppm']==pytest.approx(4.7)
    assert preprocessing['aligned_analysis_position_ppm']==4.68
    assert preprocessing['metabolite_frame_changed'] is False
    assert json.loads((tmp_path/'simulation_config.json').read_text())==generated['simulation_config']


@pytest.mark.parametrize('kwargs', [dict(points=1024),dict(bandwidth_hz=4000),dict(water_linewidth_hz=0),
                                  dict(water_frequency_offset_hz=37),dict(water_phase_rad=.6),
                                  dict(basis_path=REFERENCE_BASIS),dict(expected_metadata={'convention':'analytical'})])
def test_sequence_rejects_incompatible_configuration(reference,kwargs):
    with pytest.raises(ValueError):pair(reference,**kwargs)


def test_rejects_frozen_analytical_calibration_object(reference):
    _,_,calibration=format_frozen_basis(REFERENCE_BASIS,bandwidth_hz=2000,points=2048)
    with pytest.raises(ValueError):pair(reference,calibration=calibration)
    with pytest.raises(ValueError):generate_pair({'NAA':10},sequence_reference=reference)
    with pytest.raises(ValueError):generate_pair({'NAA':10},water_reference_mode='mrscloud_sequence')


@pytest.mark.parametrize('change',['basis_hash','scalar','scalar_both_fields','RF','version','replace_basis_and_hash','water_hash','statistic','missing'])
def test_reference_validation(reference,tmp_path,change):
    bad=tmp_path/'reference';shutil.copytree(reference,bad)
    path=bad/'reference.json';manifest=json.loads(path.read_text())
    if change=='basis_hash':manifest['basis_sha256']='0'*64
    elif change=='scalar':manifest['normalization_scalar']*=2
    elif change=='scalar_both_fields':
        manifest['normalization_scalar']*=2
        manifest['sequence_configuration']['shared_max_real_spectrum_scalar']*=2
    elif change=='RF':manifest['sequence_configuration']['RF_carriers_ppm']['water_on']=3
    elif change=='version':manifest['fsl_mrs_version']='incompatible'
    elif change=='water_hash':manifest['raw_water_sha256']='0'*64
    elif change=='statistic':manifest['normalization_statistic']='mean_absolute_sample'
    elif change=='missing':del manifest['sequence_configuration']
    else:
        shutil.copyfile(REFERENCE_BASIS,bad/'common_response.BASIS')
        manifest['basis_sha256']=file_sha256(bad/'common_response.BASIS')
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):pair(bad)


def test_legacy_residual_water_and_random_stream_unchanged(reference):
    from types import SimpleNamespace
    import torch
    from simulation.signalModel import SignalModel
    model=SignalModel(None,SimpleNamespace(B0=3,TE=30,t=np.arange(32)/2000,cf=127.731),
        labels=['WM'],simulation_params=dict(water_amp_factor=5,water_damping_min=1,
            water_damping_max=8,water_phase_min=-180,water_phase_max=180),verbose=False)
    saved=np.random.get_state();torch_state=torch.random.get_rng_state().clone()
    try:
        np.random.seed(17);baseline=model.simulate_water().numpy().copy();expected=np.random.random(8)
        np.random.seed(17);pair(reference)
        np.testing.assert_array_equal(model.simulate_water().numpy(),baseline)
        np.testing.assert_array_equal(np.random.random(8),expected)
        assert torch.equal(torch.random.get_rng_state(),torch_state)
    finally:np.random.set_state(saved)

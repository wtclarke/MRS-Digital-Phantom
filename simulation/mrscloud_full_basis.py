"""Standalone full-library synthetic responses; no legacy simulator integration."""

import hashlib
import json
from pathlib import Path

import fsl_mrs
import numpy as np
from fsl_mrs.utils import mrs_io

from simulation.quantitative_signals import concentration_signals

CONVENTION='mrscloud-press-te30-full-reference-max-v1'
STATISTIC='maximum_real_native_fft_after_both_DC_before_acquisition_shift'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalization_identity(normalization):
    content={k:v for k,v in normalization.items() if k!='calibration_identity'}
    return CONVENTION+':'+hashlib.sha256(json.dumps(content,sort_keys=True,allow_nan=False).encode()).hexdigest()


def load_library(directory, *, bandwidth_hz=2000, points=2048, selected_names=None, expected_calibration=None):
    """Check the complete fixed reference, then select unchanged response columns."""
    root=Path(directory)
    manifest=json.loads((root/'reference_manifest.json').read_text())
    required={'metabolite_basis.BASIS','water_reference.npz','normalization_reference.npz',
              'normalization.json','sequence_config.json'}
    if not required<=manifest.get('files',{}).keys():
        raise ValueError('incomplete full-library manifest')
    for name in required:
        if sha256(root/name)!=manifest['files'][name]:raise ValueError('bundle file identity mismatch: '+name)
    norm=json.loads((root/'normalization.json').read_text())
    cfg=json.loads((root/'sequence_config.json').read_text())
    if not {'convention','statistic','calibration_identity','sequence_config_sha256','source_hashes',
            'reference_names','fsl_mrs_version','scalar','reference_array_sha256','capture_configuration_signature',
            'processing_source_hashes','numpy_version'}<=norm.keys():
        raise ValueError('incomplete normalization metadata')
    if (norm['convention']!=CONVENTION or norm['statistic']!=STATISTIC
            or norm['calibration_identity']!=normalization_identity(norm)
            or norm['sequence_config_sha256']!=sha256(root/'sequence_config.json')
            or norm['source_hashes']!=manifest['source_hashes']
            or norm['processing_source_hashes']!=manifest.get('processing_source_hashes')
            or norm['numpy_version']!=np.__version__
            or norm['reference_names']!=[e['name'] for e in manifest['entries']]
            or norm['reference_array_sha256']!=manifest['files']['normalization_reference.npz']
            or norm['fsl_mrs_version']!=fsl_mrs.__version__):
        raise ValueError('incompatible normalization/reference metadata')
    if (cfg['target_bandwidth_hz'],cfg['target_points'])!=(bandwidth_hz,points):
        raise ValueError('incompatible full-library spectral grid')
    if ((cfg['TE_ms'],cfg['TE1_ms'],cfg['TE2_ms'])!=(30,14,16)
            or cfg['metabolite_RF_carrier_ppm']!=3. or cfg['water_RF_carrier_ppm']!=4.7
            or cfg['frequency_origin_ppm']!=4.68 or cfg['water_proton_multiplier_applied'] is not False
            or cfg['readout_linewidth_hz']!=1 or cfg['refocusing_duration_ms']!=7
            or cfg['refocusing_waveform']!='univ_eddenrefo.pta'
            or cfg['water_processing']!='raw_pre_DC_shared_scale_frame_shift_exact_samples'):
        raise ValueError('incompatible sequence/water processing convention')
    capture_config={'spatial_points':cfg['spatial_points_per_axis'],'native_points':cfg['native_points'],
                    'source_hashes':norm['source_hashes']}
    if hashlib.sha256(json.dumps(capture_config,sort_keys=True).encode()).hexdigest()!=norm['capture_configuration_signature']:
        raise ValueError('sequence configuration does not match raw capture identity')
    scalar=norm['scalar']
    with np.load(root/'normalization_reference.npz',allow_pickle=False) as archive:
        reference=archive['processed_reference']
    names=norm['reference_names']
    if reference.shape!=(cfg['native_points'],len(names)) or not np.all(np.isfinite(reference)):
        raise ValueError('invalid normalization reference arrays')
    measured=float(np.max(np.fft.fftshift(np.fft.fft(reference,axis=0),axes=0).real))
    if not np.isfinite(scalar) or scalar<=0 or not np.isclose(scalar,measured,rtol=1e-12,atol=0):
        raise ValueError('normalization scalar does not match fixed reference')
    basis=mrs_io.read_basis(str(root/'metabolite_basis.BASIS'))
    if (list(basis.names)!=names or basis.original_points!=cfg['native_points']
            or not np.isclose(basis.original_bw,cfg['native_bandwidth_hz'],rtol=1e-12)
            or not np.isclose(basis.cf,cfg['spectrometer_frequency_mhz'],rtol=1e-12)):
        raise ValueError('basis order/sampling/frequency mismatch')
    t=np.arange(cfg['native_points'])/cfg['native_bandwidth_hz']
    expected=reference/scalar*np.exp(2j*np.pi*(3.-4.68)*basis.cf*t)[:,None]
    errors=np.linalg.norm(basis.original_basis_array-expected,axis=0)/np.linalg.norm(expected,axis=0)
    if not np.all(np.isfinite(errors)) or np.max(errors)>3e-7:
        raise ValueError('basis does not retain shared reference amplitudes')
    with np.load(root/'water_reference.npz',allow_pickle=False) as archive:
        water=archive['fid'].copy();raw=archive['raw_native']
    if raw.shape!=(cfg['native_points'],) or not np.iscomplexobj(raw) or not np.all(np.isfinite(raw)):
        raise ValueError('invalid native water response')
    step=cfg['native_bandwidth_hz']/bandwidth_hz
    if step<1 or not np.isclose(step,round(step)) or (points-1)*round(step)>=len(raw):
        raise ValueError('water sampling requires covered exact native acquisition times')
    expected_water=(raw*np.exp(2j*np.pi*(4.7-4.68)*basis.cf*t)/scalar)[np.arange(points)*int(round(step))]
    if water.shape!=(points,) or not np.all(np.isfinite(water)) or not np.allclose(water,expected_water,rtol=2e-15,atol=0):
        raise ValueError('water reference does not share raw sequence scale and sampling')
    metadata={'calibration_identity':norm['calibration_identity'],'normalization':norm,'sequence':cfg,
              'reference_manifest_sha256':sha256(root/'reference_manifest.json'),'physical_calibration':False}
    if expected_calibration is not None and expected_calibration!=metadata:
        raise ValueError('full-library calibration metadata mismatch')
    selected=list(selected_names) if selected_names is not None else [e['name'] for e in manifest['entries'] if e['role']=='metabolite']
    if not selected or len(set(selected))!=len(selected) or not set(selected)<=set(names):
        raise ValueError('unknown, repeated or empty response selection')
    formatted=basis.get_formatted_basis(bandwidth_hz,points)
    return {name:formatted[:,names.index(name)].copy() for name in selected},water,metadata


def synthesize(directory, concentrations_mM, *, water_concentration_mol_l=40., receiver_factor=1.,
               bandwidth_hz=2000, points=2048, expected_calibration=None):
    if {'H2O','Ref0ppm'} & set(concentrations_mM):
        raise ValueError('reference controls are not metabolite concentrations; water is supplied separately')
    responses,water,metadata=load_library(directory,bandwidth_hz=bandwidth_hz,points=points,
        selected_names=list(concentrations_mM),expected_calibration=expected_calibration)
    components,water_signal=concentration_signals(responses,water,
        concentrations_mol_l={k:float(v)*.001 for k,v in concentrations_mM.items()},
        water_concentration_mol_l=water_concentration_mol_l,receiver_factor=receiver_factor)
    cfg=metadata['sequence']
    return {'metabolite_fid':sum(components.values()),'water_reference_fid':water_signal,'components':components,
        'ground_truth':{'metabolite_concentrations_mM':{k:float(v) for k,v in concentrations_mM.items()},
            'water_concentration_mol_l':float(water_concentration_mol_l),
            'concentration_convention':'nominal concentration per total tissue volume','physical_calibration':False},
        'simulation_config':{'synthetic_calibration':metadata,'receiver_factor':float(receiver_factor),
            'points':points,'bandwidth_hz':bandwidth_hz,'dwell_time_s':1/bandwidth_hz,
            'frequency_origin_ppm':4.68,'spectrometer_frequency_mhz':cfg['spectrometer_frequency_mhz'],
            'water_frequency_offset_hz':(4.7-4.68)*cfg['spectrometer_frequency_mhz'],
            'water_phase_rad':0.,'residual_water':False,'physical_calibration':False}}

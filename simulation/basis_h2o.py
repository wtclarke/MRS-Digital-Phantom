"""Explicitly approximate water referencing from a supplied, unchanged LCModel BASIS."""
from copy import deepcopy
import hashlib
import json
from numbers import Integral
from pathlib import Path
import re

import fsl_mrs
import numpy as np
from fsl_mrs.utils import mrs_io

MODE='basis_h2o'
ASSUMPTION='assume_common_scale_protons_encoded'
EXTRA={'basis_path','concentration_units','basis_amplitude_convention','frequency_origin_ppm',
       'water_component_name','expected_basis_sha256','basis_relaxation_assumption'}


def _field(header,name, *, string=False):
    pattern=rf"\b{name}\s*=\s*'([^']*)'" if string else rf'\b{name}\s*=\s*([-+\d.eEdD]+)'
    match=re.search(pattern,header,re.I)
    if match is None:raise ValueError(f'BASIS requires explicit {name}')
    return match[1].strip() if string else float(match[1].replace('D','E').replace('d','e'))


def load_responses(path, *, frequency_origin_ppm, bandwidth_hz=2000,points=2048,
                   water_component_name='H2O',expected_basis_sha256=None):
    """Format the entire fixed reference without any additional amplitude scaling."""
    path=Path(path);content=path.read_bytes();digest=hashlib.sha256(content).hexdigest()
    if expected_basis_sha256 is not None and digest!=expected_basis_sha256:
        raise ValueError('BASIS SHA256 mismatch')
    if (isinstance(points,bool) or not isinstance(points,Integral) or points<2
            or isinstance(bandwidth_hz,bool) or not np.isfinite(bandwidth_hz) or bandwidth_hz<=0
            or not np.isfinite(frequency_origin_ppm) or not 0<=frequency_origin_ppm<=10):
        raise ValueError('finite Hz bandwidth, integer points >=2 and declared ppm origin are required')
    text=content.decode('ascii')
    blocks=re.findall(r'\$(SEQPAR|BASIS1|BASIS)\b(.*?)\$END',text,re.I|re.S)
    headers={kind.upper():header for kind,header in blocks if kind.upper()!='BASIS'}
    entries=[header for kind,header in blocks if kind.upper()=='BASIS']
    if any(sum(kind.upper()==key for kind,_ in blocks)!=1 for key in ['SEQPAR','BASIS1']) or not entries:
        raise ValueError('missing LCModel BASIS headers/components')
    native_points=_field(headers['BASIS1'],'NDATAB');dt=_field(headers['BASIS1'],'BADELT')
    cf=_field(headers['SEQPAR'],'HZPPPM');te=_field(headers['SEQPAR'],'ECHOT')
    if not all(np.isfinite(v) and v>0 for v in [native_points,dt,cf,te]) or native_points!=int(native_points):
        raise ValueError('inconsistent BASIS sampling/sequence metadata')
    component_headers={}
    for header in entries:
        name=_field(header,'METABO',string=True)
        if name in component_headers:raise ValueError('duplicate BASIS component names')
        values={key:_field(header,key) for key in ['CONC','TRAMP','VOLUME','ISHIFT']}
        if any(values[key]!=1. for key in ['CONC','TRAMP','VOLUME']):
            raise ValueError('unsupported CONC/TRAMP/VOLUME: common scale cannot be established')
        if not np.isfinite(values['ISHIFT']) or values['ISHIFT']!=int(values['ISHIFT']):
            raise ValueError('invalid component frequency shift')
        component_headers[name]=values
    basis=mrs_io.read_basis(str(path))
    names=list(map(str,basis.names));native=basis.original_basis_array
    if (len(set(names))!=len(names) or set(names)!=set(component_headers)
            or basis.original_points!=native_points or native.shape!=(int(native_points),len(names))
            or not np.isclose(basis.original_dwell,dt,rtol=1e-12,atol=0)
            or not np.isclose(basis.cf,cf,rtol=1e-12,atol=0)
            or not np.iscomplexobj(native) or not np.all(np.isfinite(native))):
        raise ValueError('inconsistent or nonfinite BASIS complex samples/sampling')
    if water_component_name not in names:raise ValueError(f'BASIS has no {water_component_name} component')
    if bandwidth_hz>1/dt or (points-1)/bandwidth_hz>(native_points-1)*dt:
        raise ValueError('requested sampling exceeds native BASIS bandwidth/time coverage')
    # Full reference always includes water; no legacy mean-magnitude or per-column scaling.
    formatted=basis.get_formatted_basis(bandwidth_hz,points,scale_factor=None)
    water=formatted[:,names.index(water_component_name)].copy()
    if not np.all(np.isfinite(formatted)) or not np.any(abs(water)>0):
        raise ValueError('H2O response must be finite and nonzero')
    frequencies=np.fft.fftshift(np.fft.fftfreq(points,1/bandwidth_hz))
    peak=float(frequencies[np.argmax(abs(np.fft.fftshift(np.fft.fft(water))))])
    descriptor={'basis_sha256':digest,'bandwidth_hz':float(bandwidth_hz),'points':int(points),
        'frequency_origin_ppm':float(frequency_origin_ppm),'water_component_name':water_component_name,
        'reference_names':names,'normalization':'identity_preserve_supplied_common_scale',
        'amplitude_assumption':ASSUMPTION,'fsl_mrs_version':fsl_mrs.__version__,'numpy_version':np.__version__}
    identity='basis-h2o-approximate-v1:'+hashlib.sha256(json.dumps(descriptor,sort_keys=True).encode()).hexdigest()
    metadata={**descriptor,'calibration_identity':identity,'basis_filename':path.name,
        'physical_calibration':False,'approximate_calibration':True,
        'normalization':{'scalar':1.,'operation':'divide all formatted responses by 1; no extra normalization',
            'statistic':'identity_preserve_supplied_common_scale','reference_names':names,
            'selection_independent':True},'component_headers':component_headers,
        'sequence':{'sequence':_field(headers['SEQPAR'],'SEQ',string=True),'TE_ms':te,
            'spectrometer_frequency_mhz':cf,'frequency_origin_ppm':float(frequency_origin_ppm),
            'native_points':int(native_points),'native_bandwidth_hz':1/dt,
            'water_frequency_offset_hz':peak,'water_phase_rad':float(np.angle(water[0])),
            'water_frequency_estimate':'magnitude FFT maximum bin, not a fitted carrier',
            'water_proton_multiplier_applied':False},
        'limitations':['Common molar response and encoded proton counts are explicitly assumed, not proven by CONC=1.',
            'RF carrier/refocusing sensitivity, intrinsic phase, DC correction, first point and readout are retained.',
            'No off-resonance compensation or per-metabolite correction; not a water-centred acquisition.']}
    return {n:formatted[:,i].copy() for i,n in enumerate(names) if n not in {water_component_name,'Ref0ppm'}},water,metadata


def validate_config(config):
    from simulation.quantitative import _validate_common_config
    if not isinstance(config,dict):raise ValueError('basis_h2o configuration must be an object')
    if config.get('basis_amplitude_convention')!=ASSUMPTION:
        raise ValueError(f'approximate mode requires basis_amplitude_convention={ASSUMPTION}; independently normalized bases are unsupported')
    if config.get('concentration_units')!='mM':raise ValueError('explicit concentration_units=mM required; IU is not converted')
    if 'basis_path' not in config or 'frequency_origin_ppm' not in config:
        raise ValueError('basis_path and explicit frequency_origin_ppm are required')
    if 'calibration_bundle' in config:raise ValueError('basis_h2o takes a BASIS file, not a Phase 6a calibration bundle')
    common={k:v for k,v in config.items() if k not in EXTRA}
    common['calibration_bundle']=str(config['basis_path'])
    cfg=_validate_common_config(common)
    cfg.setdefault('bandwidth_hz',2000);cfg.setdefault('points',2048)
    _,_,meta=load_responses(config['basis_path'],frequency_origin_ppm=config['frequency_origin_ppm'],
        bandwidth_hz=cfg['bandwidth_hz'],points=cfg['points'],water_component_name=config.get('water_component_name','H2O'),
        expected_basis_sha256=config.get('expected_basis_sha256'))
    if cfg.get('expected_calibration_identity',meta['calibration_identity'])!=meta['calibration_identity']:
        raise ValueError('expected calibration identity mismatch')
    if (cfg.get('relaxation') or {}).get('mode')=='tissue_specific' and config.get('basis_relaxation_assumption')!='no_tissue_attenuation_encoded':
        raise ValueError('active relaxation requires explicit basis_relaxation_assumption=no_tissue_attenuation_encoded')
    cfg.update({k:deepcopy(v) for k,v in config.items() if k in EXTRA})
    return cfg,meta


def simulate_basis_h2o(config, *, phantom=None,voi_bounds=None):
    from simulation.quantitative import _simulate_with_config
    cfg,metadata=validate_config(config)
    responses,water,meta=load_responses(cfg['basis_path'],frequency_origin_ppm=cfg['frequency_origin_ppm'],
        bandwidth_hz=cfg['bandwidth_hz'],points=cfg['points'],water_component_name=cfg.get('water_component_name','H2O'),
        expected_basis_sha256=metadata['basis_sha256'])
    def loader(path, *, bandwidth_hz,points,selected_names,expected_calibration):
        if (str(path)!=str(cfg['basis_path']) or bandwidth_hz!=cfg['bandwidth_hz'] or points!=cfg['points']
                or expected_calibration!=meta):raise ValueError('incompatible BASIS identity/sampling metadata')
        if not selected_names or len(set(selected_names))!=len(selected_names) or not set(selected_names)<=set(responses):
            raise ValueError('unknown/repeated metabolite selection or H2O/reference selected as metabolite')
        return {n:responses[n].copy() for n in selected_names},water.copy(),deepcopy(meta)
    pair=_simulate_with_config(config,cfg,metadata,phantom=phantom,voi_bounds=voi_bounds,
        response_loader=loader,simulation_mode=MODE)
    pair['simulation_config'].update(water_reference_mode=MODE,approximate_calibration=True,
        source_hashes_verified=False,basis_file_sha256_verified=True,basis_sha256=meta['basis_sha256'],concentration_units='mM',
        pathway='approximate_basis_h2o_spatial_voxel' if pair['simulation_config']['spatial_simulation'] else 'approximate_basis_h2o_voxel',
        basis_amplitude_convention=ASSUMPTION,TE_s=meta['sequence']['TE_ms']/1000)
    pair['ground_truth'].update(water_reference_kind='approximate supplied BASIS H2O response; not a separate water-centred sequence',
        approximate_calibration=True)
    pair['ground_truth']['relaxation']['basis_assumption']='supplied sequence/readout preserved; additional tissue attenuation explicitly assumed absent in BASIS when enabled'
    pair['simulation_config']['relaxation']=deepcopy(pair['ground_truth']['relaxation'])
    return pair


def simulate_batch(positions,settings, *, phantom=None,save_dir=None):
    from simulation.quantitative import simulate_quantitative_batch
    if settings.get('simulation_mode','legacy') not in ('legacy',MODE):
        raise ValueError('basis_h2o cannot be combined with sequence_matched_quantitative')
    settings=deepcopy(settings);settings.pop('water_reference_mode')
    if 'path2basis' not in settings:raise ValueError('basis_h2o requires path2basis')
    config=settings['quantitative']
    if 'basis_path' in config:raise ValueError('headless basis_h2o uses top-level path2basis only')
    config['basis_path']=settings.pop('path2basis')
    settings['simulation_mode']=MODE
    return simulate_quantitative_batch(positions,settings,phantom=phantom,save_dir=save_dir,
        _simulate=simulate_basis_h2o,_validate=validate_config,_mode=MODE)

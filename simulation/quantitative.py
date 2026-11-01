"""Production dispatch for the opt-in Phase 6a sequence-response family."""
from copy import deepcopy
import json
from numbers import Integral
from pathlib import Path

import numpy as np
from nifti_mrs.create_nmrs import gen_nifti_mrs

from simulation.mrscloud_full_basis import CONVENTION,load_library
from simulation.mrscloud_quantitative_noise import add_pair_noise
from simulation.mrscloud_quantitative_voxel import simulate_voxel,simulate_phantom_voxel
from simulation.mrscloud_quantitative_spatial import simulate_spatial,simulate_phantom_spatial
from utils.definitions import H2O_CONCENTRATIONS
from simulation.calibration_resources import calibration_bundle,verify_source_provenance

MODE='sequence_matched_quantitative'
CALIBRATION_ID=CONVENTION+':d56c1adc96044bd550a3b43923f5b036b8fa480e234a21ea422b79accda91a81'
DEFAULT_WATER={t:float(H2O_CONCENTRATIONS[t]) for t in ('WM','GM','CSF')}
KEYS={'calibration_bundle','metabolite_concentrations_mM','water_concentrations_mol_l',
      'receiver_factor','bandwidth_hz','points','relaxation','noise','tissue_fractions',
      'spatial','phantom_spatial','additional_linewidth_hz','voxel_affine_mm',
      'expected_calibration_identity','concentration_convention'}


def _validate_common_config(config):
    """Shared units, profiles and acquisition controls; no calibration-family decision."""
    if not isinstance(config,dict) or set(config)-KEYS:
        raise ValueError(f'unknown quantitative configuration keys: {set(config)-KEYS if isinstance(config,dict) else config}')
    cfg=deepcopy(config)
    if 'metabolite_concentrations_mM' not in cfg:
        raise ValueError('quantitative mode requires explicit tissue metabolite_concentrations_mM')
    profiles=cfg['metabolite_concentrations_mM']
    if not isinstance(profiles,dict) or set(profiles)!={'WM','GM','CSF'} or any(not isinstance(v,dict) or not v for v in profiles.values()):
        raise ValueError('explicit nonempty WM, GM and CSF metabolite profiles are required')
    if cfg.get('noise') is not None:
        noise=cfg['noise']
        if (not isinstance(noise,dict) or set(noise)-{'seed','metabolite_sigma','metabolite_snr','water_sigma','water_snr'}
                or isinstance(noise.get('seed'),bool) or not isinstance(noise.get('seed'),Integral) or noise['seed']<0):
            raise ValueError('noise requires a nonnegative integer seed and independent metabolite/water sigma or SNR')
    if cfg.get('concentration_convention','nominal_per_total_tissue_volume')!='nominal_per_total_tissue_volume':
        raise ValueError('only nominal concentrations per total tissue volume are supported')
    if sum(k in cfg for k in ['tissue_fractions','spatial','phantom_spatial'])>1:
        raise ValueError('choose tissue_fractions, spatial cells or phantom_spatial, not a mixture')
    linewidth=cfg.get('additional_linewidth_hz',{'metabolites':0.,'water':0.})
    if not isinstance(linewidth,dict) or set(linewidth)!={'metabolites','water'}:
        raise ValueError('additional_linewidth_hz requires metabolites and water in Hz')
    if any(isinstance(v,bool) for v in linewidth.values()):raise ValueError('linewidths require numeric Hz, not booleans')
    try:linewidth={k:float(v) for k,v in linewidth.items()}
    except (TypeError,ValueError) as error:raise ValueError('linewidths require numeric Hz') from error
    if any(not np.isfinite(v) or v<0 for v in linewidth.values()):
        raise ValueError('additional linewidths must be finite and nonnegative Hz')
    cfg['additional_linewidth_hz']={k:float(v) for k,v in linewidth.items()}
    cfg.setdefault('water_concentrations_mol_l',deepcopy(DEFAULT_WATER))
    return cfg


def validate_config(config):
    """Validate the frozen sequence family and its construction-source hashes."""
    cfg=_validate_common_config(config)
    cfg['calibration_bundle']=str(calibration_bundle(cfg.get('calibration_bundle')))
    bw=cfg.get('bandwidth_hz',2000);points=cfg.get('points',2048)
    if isinstance(points,bool) or not isinstance(points,int) or points!=2048 or bw!=2000:
        raise ValueError('production calibration requires 2048 complex points at 2000 Hz')
    _,_,meta=load_library(cfg['calibration_bundle'],bandwidth_hz=bw,points=points,
                         selected_names=list(cfg['metabolite_concentrations_mM'].get('WM',{})))
    sequence=meta['sequence']
    if (meta['normalization']['convention']!=CONVENTION or sequence['resolution_family']!='production'
            or sequence['sequence']!='Universal_Philips_UnEdited_PRESS'
            or sequence['spectrometer_frequency_mhz']!=127.731):
        raise ValueError('requires the full production PRESS TE30 calibration family at 127.731 MHz')
    if meta['calibration_identity']!=CALIBRATION_ID or cfg.get('expected_calibration_identity',CALIBRATION_ID)!=meta['calibration_identity']:
        raise ValueError('expected calibration identity mismatch')
    verify_source_provenance(cfg['calibration_bundle'],meta['normalization'])
    return cfg,meta


def simulate_quantitative(config, *, phantom=None, voi_bounds=None):
    """Return paired complex FIDs; no legacy signal generation is called."""
    cfg,calibration=validate_config(config)
    return _simulate_with_config(config,cfg,calibration,phantom=phantom,voi_bounds=voi_bounds)


def _simulate_with_config(config,cfg,calibration, *, phantom=None,voi_bounds=None,
                          response_loader=None,simulation_mode=MODE):
    kwargs={k:cfg[k] for k in ['calibration_bundle','metabolite_concentrations_mM','water_concentrations_mol_l',
                             'receiver_factor','bandwidth_hz','points','relaxation'] if k in cfg}
    kwargs['expected_calibration']=calibration
    if response_loader is not None:kwargs['_response_loader']=response_loader
    if phantom is not None:
        if 'tissue_fractions' in cfg or 'spatial' in cfg:
            raise ValueError('phantom ROIs cannot also specify fractions or inline cells')
        if voi_bounds is None:raise ValueError('phantom simulation requires exclusive-stop ROI bounds')
        if 'voxel_affine_mm' in cfg:raise ValueError('phantom affine is derived from its ROI, not overridden')
        if 'phantom_spatial' in cfg:
            options=cfg['phantom_spatial']
            if not isinstance(options,dict) or set(options)-{'use_phantom_b0','b0_units','b0_convention','spatial_weights'}:
                raise ValueError('unknown phantom_spatial option')
            if not isinstance(options.get('use_phantom_b0',False),bool):raise ValueError('use_phantom_b0 must be boolean')
            pair=simulate_phantom_spatial(phantom,voi_bounds,**options,**kwargs)
        else:pair=simulate_phantom_voxel(phantom,voi_bounds,**kwargs)
        original_affine=np.asarray(phantom.affine,float)
        if (not np.allclose(np.linalg.norm(original_affine[:3,:3],axis=0),phantom.spacing,rtol=1e-6,atol=1e-8)
                or not np.allclose(phantom.subject['labels'].affine,original_affine,rtol=0,atol=1e-6)):
            raise ValueError('anatomical spacing, phantom affine and label affine must agree')
        bounds=tuple(voi_bounds);counts=np.array(bounds[1::2])-np.array(bounds[::2])
        center=(np.array(bounds[::2])+np.array(bounds[1::2])-1)/2
        affine=np.array(phantom.affine,float,copy=True)
        affine[:3,3]=(affine@np.r_[center,1])[:3]
        affine[:3,:3]=affine[:3,:3]@np.diag(counts)
        pair['ground_truth']['voxel_geometry'].update(cell_volume_mm3=float(abs(np.linalg.det(original_affine[:3,:3]))),
            voxel_volume_mm3=float(abs(np.linalg.det(original_affine[:3,:3]))*np.prod(counts)),
            cell_volume_definition='abs(det(anatomical_affine[:3,:3]))')
        geometry_scope='anatomical ROI; geometric cell-center midpoint and actual covered extent'
    elif 'spatial' in cfg:
        options=deepcopy(cfg['spatial'])
        allowed={'cell_labels','spatial_weights','cell_volumes_mm3','b0_offsets_hz','b0_units','b0_convention'}
        if not isinstance(options,dict) or set(options)-allowed or 'cell_labels' not in options:raise ValueError('spatial requires cell_labels and recognized spatial options')
        if 'b0_offsets_hz' in options:
            if options.pop('b0_units',None)!='Hz' or options.get('b0_convention') not in ['positive_frequency','legacy_shim']:
                raise ValueError('explicit B0 requires Hz units and positive_frequency or legacy_shim sign')
            sign=-1 if options.pop('b0_convention')=='legacy_shim' else 1
            options['b0_offsets_hz']=np.asarray(options['b0_offsets_hz'],float)*sign
        elif 'b0_units' in options or 'b0_convention' in options:
            raise ValueError('B0 declarations require explicit offsets')
        pair=simulate_spatial(**options,**kwargs)
        pair['simulation_config']['input_spatial_configuration']=deepcopy(cfg['spatial'])
        affine=np.asarray(cfg.get('voxel_affine_mm',np.eye(4)),float)
        geometry_scope='explicit voxel affine' if 'voxel_affine_mm' in cfg else 'identity; no anatomical registration supplied'
    else:
        if 'phantom_spatial' in cfg or voi_bounds is not None:raise ValueError('phantom options require an anatomical phantom')
        if 'tissue_fractions' not in cfg:raise ValueError('supply tissue_fractions, spatial cells or a phantom ROI')
        pair=simulate_voxel(cfg['tissue_fractions'],**kwargs)
        affine=np.asarray(cfg.get('voxel_affine_mm',np.eye(4)),float)
        geometry_scope='explicit voxel affine' if 'voxel_affine_mm' in cfg else 'identity; no anatomical registration supplied'
    if (affine.shape!=(4,4) or not np.all(np.isfinite(affine))
            or not np.allclose(affine[3],[0,0,0,1],rtol=0,atol=1e-12)
            or abs(np.linalg.det(affine[:3,:3]))<=0):
        raise ValueError('voxel affine requires finite, nonsingular homogeneous geometry in mm')
    line=cfg['additional_linewidth_hz'];t=pair['time_s']
    if line['metabolites']:
        envelope=np.exp(-np.pi*line['metabolites']*t)
        pair['components']={m:fid*envelope for m,fid in pair['components'].items()}
        pair['metabolite_fid']=sum(pair['components'].values())
    if line['water']:pair['water_reference_fid']=pair['water_reference_fid']*np.exp(-np.pi*line['water']*t)
    pair['simulation_config'].update(simulation_mode=simulation_mode,production_quantitative=True,
        water_concentration_source='explicit mol/L override' if 'water_concentrations_mol_l' in config else 'H2O_CONCENTRATIONS validated tissue-volume defaults',
        additional_linewidth_hz=line,linewidth_convention='additional Lorentzian exp(-pi*Hz*t); frozen readout retained; observed FWHM not inferred',
        voxel_affine_mm=affine.tolist(),geometry_scope=geometry_scope,
        source_hashes_verified=True,concentration_convention='nominal_per_total_tissue_volume')
    pair['ground_truth'].update(signal_volume_convention='normalized represented-volume average; no absolute anatomical volume gain',
        water_reference_kind='separate raw pre-DC sequence-matched unsuppressed acquisition',
        residual_water_contamination=False)
    truth=pair['ground_truth']
    truth.setdefault('relaxation',{'mode':'none','TE_s':calibration['sequence']['TE_ms']/1000,'TR_s':None,'application':'pre_acquisition_amplitude_only',
        'times_unit':'s','water_factors':dict.fromkeys(DEFAULT_WATER,1.),
        'metabolite_factors':{t:dict.fromkeys(pair['components'],1.) for t in DEFAULT_WATER}})
    truth.setdefault('relaxation_weighted_metabolite_concentrations_mM',dict(truth['voxel_average_metabolite_concentrations_mM']))
    truth.setdefault('relaxation_weighted_water_concentration_mol_l',truth['effective_voxel_water_concentration_mol_l'])
    pair['simulation_config']['relaxation']=deepcopy(truth['relaxation'])
    return add_pair_noise(pair,**cfg['noise']) if cfg.get('noise') is not None else pair


def export_quantitative(pair,output_dir):
    """Validated complex 1D export with anatomical affine and JSON-only provenance."""
    output=Path(output_dir)
    if output.exists():raise FileExistsError(f'quantitative output already exists: {output}')
    output.mkdir(parents=True)
    cfg=pair['simulation_config'];paths={}
    signals={'metabolites':pair['metabolite_fid'],'water_reference':pair['water_reference_fid']}
    if 'clean_metabolite_fid' in pair:
        signals.update(metabolites_clean=pair['clean_metabolite_fid'],water_reference_clean=pair['clean_water_reference_fid'])
    for name,fid in signals.items():
        obj=gen_nifti_mrs(fid[None,None,None,:,None],cfg['dwell_time_s'],cfg['spectrometer_frequency_mhz'],
                         nucleus='1H',dim_tags=['DIM_DYN',None,None],affine=np.asarray(cfg['voxel_affine_mm']))
        obj.add_hdr_field('SpecFreqChemShift',cfg['frequency_origin_ppm'])
        obj.add_hdr_field('EchoTime',cfg.get('TE_s',.030))
        if 'TR_s' in cfg:obj.add_hdr_field('RepetitionTime',cfg['TR_s'])
        path=output/f'{name}.nii.gz';obj.save(str(path));paths[name]=str(path)
    for name in ['ground_truth','simulation_config']:
        (output/f'{name}.json').write_text(json.dumps(pair[name],indent=2,allow_nan=False)+'\n')
    (output/'calibration_provenance.json').write_text(json.dumps(cfg['synthetic_calibration'],indent=2,allow_nan=False)+'\n')
    return paths


def simulate_quantitative_batch(positions,settings,*,phantom=None,save_dir=None,
                                _simulate=None,_validate=None,_mode=MODE):
    from simulation.headless import _position_and_bounds
    allowed={'simulation_mode','quantitative','positions','voxel_size_mm','skeleton','output_dir'}
    if set(settings)-allowed:raise ValueError('legacy or unknown settings supplied to quantitative mode')
    simulate=simulate_quantitative if _simulate is None else _simulate
    validate=validate_config if _validate is None else _validate
    config=settings['quantitative'];validate(config)
    if phantom is None and not {'tissue_fractions','spatial'}&set(config):
        from simulation.mrs_phantom import DigitalPhantom
        if 'skeleton' not in settings:raise ValueError('anatomical mode requires a skeleton or supplied phantom')
        phantom=DigitalPhantom(settings['skeleton'],None,None,{'simulation_mode':_mode,'quantitative':config},verbose=False)
    positions=list(positions);bounds=[]
    if not positions and save_dir is not None:raise ValueError('at least one quantitative voxel is required when saving')
    if phantom is None and ('skeleton' in settings or 'voxel_size_mm' in settings):
        raise ValueError('fractions/cells use voxel_affine_mm; skeleton/voxel_size_mm require an anatomical ROI')
    for position in positions:
        try:position=tuple(position)
        except TypeError as error:raise ValueError('voxel positions must contain three integer indices') from error
        if any(isinstance(v,bool) for v in position):raise ValueError('voxel positions must be integer indices, not booleans')
        if phantom is not None:
            size=np.asarray(settings.get('voxel_size_mm'),float)
            if size.shape!=(3,) or np.any(size<=0) or not np.all(np.isfinite(size)):raise ValueError('anatomical batch requires positive voxel_size_mm')
            counts=tuple(round(a/b) for a,b in zip(size,phantom.spacing))
            if min(counts)<1:raise ValueError('voxel must cover at least one anatomical cell per axis')
            bounds.append(_position_and_bounds(position,counts,phantom.get_label_data().shape))
        else:bounds.append(_position_and_bounds(position,(1,1,1),(2**31,)*3))
    outputs=[Path(save_dir)/f'voxel_{i:04d}_{p[0]}_{p[1]}_{p[2]}' for i,(p,_) in enumerate(bounds)] if save_dir is not None else []
    if any(p.exists() for p in outputs):raise FileExistsError('one or more quantitative voxel outputs already exist')
    results=[]
    for i,(position,roi) in enumerate(bounds):
        voxel_config=deepcopy(config)
        if voxel_config.get('noise') is not None:
            voxel_config['noise']['seed']+=i
        pair=simulate(voxel_config,phantom=phantom,voi_bounds=roi if phantom is not None else None)
        pair['position']=position
        pair['simulation_config']['batch']={'index':i,'position_index':list(position),
            'noise_seed_policy':'configured root seed + voxel index; separate metabolite/water child streams'}
        if save_dir is not None:
            pair['saved_paths']=export_quantitative(pair,outputs[i]);pair['saved_path']=str(outputs[i])
        results.append(pair)
    return results

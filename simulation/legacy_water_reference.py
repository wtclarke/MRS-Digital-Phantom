"""Opt-in BASIS water on the legacy concentration, shim and volume convention."""
import json
from numbers import Integral
from utils.definitions import H2O_CONCENTRATIONS
import hashlib
from pathlib import Path
import numpy as np
import torch
from fsl_mrs.utils import mrs_io
from nifti_mrs.create_nmrs import gen_nifti_mrs
from simulation.basis_h2o import ASSUMPTION,load_responses
from utils.auxillary import fid2spec,spec2fid


DEFAULT_WATER={t:float(H2O_CONCENTRATIONS[t]) for t in ('WM','GM','CSF')}

def attenuation(*, TE_s, TR_s, T1_s, T2_s):
    """90-degree excitation, independent saturation and echo attenuation."""
    values = np.asarray([TE_s, TR_s, T1_s, T2_s], dtype=float)
    TE_s, TR_s, T1_s, T2_s = values
    if not np.all(np.isfinite(values)) or TE_s < 0 or TR_s <= TE_s or T1_s <= 0 or T2_s <= 0:
        raise ValueError('times must be finite seconds: TE >= 0, TR > TE, T1 > 0, T2 > 0')
    return float(-np.expm1(-TR_s/T1_s)*np.exp(-TE_s/T2_s))

def complex_gaussian_noise(points, *, sigma, seed):
    """sigma is the SD of EACH real/imaginary FID component, in signal units."""
    if isinstance(points,bool) or not isinstance(points,Integral) or points<=0:
        raise ValueError('points must be a positive integer')
    if not np.isfinite(sigma) or sigma<0:raise ValueError('sigma must be finite and nonnegative')
    if isinstance(seed,bool) or not isinstance(seed,(Integral,np.random.SeedSequence)) or (isinstance(seed,Integral) and seed<0):
        raise ValueError('seed must be a nonnegative integer or explicit SeedSequence')
    if sigma==0:return np.zeros(points,complex)
    values=np.random.Generator(np.random.PCG64(seed)).normal(0,sigma,size=(2,points))
    return values[0]+1j*values[1]


def prepare(settings):
    options=dict(settings.get('water_reference',{}))
    allowed={'concentration_convention','basis_amplitude_convention','frequency_origin_ppm',
             'water_concentrations_mol_l','additional_linewidth_hz','noise_sigma','relaxation','basis_relaxation_assumption'}
    if set(options)-allowed:raise ValueError('unknown legacy water-reference settings')
    if options.get('concentration_convention')!='assume_table_values_nominal_mM' or options.get('basis_amplitude_convention')!=ASSUMPTION:
        raise ValueError('explicit nominal-mM table and common-scale/proton assumptions required')
    if 'frequency_origin_ppm' not in options:raise ValueError('declare frequency_origin_ppm')
    sim=settings['simulation_params'];bw=sim['bandwidth'];points=sim['spectral_points']
    _,_,meta=load_responses(settings['path2basis'],frequency_origin_ppm=options['frequency_origin_ppm'],bandwidth_hz=bw,points=points)
    if not np.isclose(settings['basis_set_settings']['TE'],meta['sequence']['TE_ms']):raise ValueError('BASIS and configured TE mismatch')
    basis=mrs_io.read_basis(settings['path2basis']);all_fids=basis.get_formatted_basis(bw,points,scale_factor=None)
    scalar=float(np.mean(abs(all_fids)))
    if not np.isfinite(scalar) or scalar<=0:raise ValueError('invalid fixed reference normalization')
    water=options.setdefault('water_concentrations_mol_l',dict(DEFAULT_WATER))
    if set(water)!=set(DEFAULT_WATER) or any(not np.isfinite(v) or v<0 for v in water.values()):raise ValueError('water profiles require nonnegative WM/GM/CSF mol/L')
    for key in ['noise_sigma','additional_linewidth_hz']:
        value=options.setdefault(key,0.)
        if isinstance(value,bool) or not np.isfinite(value) or value<0:raise ValueError(f'{key} requires finite nonnegative values')
    factors=dict.fromkeys(DEFAULT_WATER,1.)
    if options.get('relaxation') is not None:
        if options.get('basis_relaxation_assumption')!='no_tissue_attenuation_encoded':raise ValueError('declare no_tissue_attenuation_encoded before adding water relaxation')
        if set(options['relaxation'])!=set(DEFAULT_WATER):raise ValueError('relaxation requires WM/GM/CSF T1_s,T2_s')
        for tissue,p in options['relaxation'].items():
            if set(p)!={'T1_s','T2_s'}:raise ValueError('relaxation requires T1_s,T2_s seconds')
            factors[tissue]=attenuation(TE_s=meta['sequence']['TE_ms']/1000,TR_s=sim['TR']/1000,**p)
    meta['normalization']={'scalar':scalar,'statistic':'mean_absolute_full_formatted_basis',
        'reference_names':list(map(str,basis.names)),'bandwidth_hz':bw,'points':points}
    meta['calibration_identity']='legacy-basis-h2o-full-mean-v1:'+hashlib.sha256((meta['calibration_identity']+json.dumps(meta['normalization'],sort_keys=True)).encode()).hexdigest()
    meta['concentration_convention']=options['concentration_convention']
    meta['signal_volume_convention']='legacy cell-volume sum in mm^3'
    meta['shim_convention']='legacy exp(-2j*pi*shim_Hz*basis.t); basis.t starts at dwell time'
    meta['water_units_conversion']='mol/L * 1000 -> nominal mM; no proton multiplier'
    meta['water_factors']=factors;meta['water_options']=options
    return meta


def generate(model,labels,metadata,seed,index):
    options=metadata['water_options'];basis=model.basis
    t=np.arange(basis.n)*basis.dwelltime
    fid=basis.water_fid*np.exp(-np.pi*options['additional_linewidth_hz']*t)
    concentrations=[1000*options['water_concentrations_mol_l'].get(name,0)*metadata['water_factors'].get(name,1) for name in model.labels]
    spectra=torch.tensor(fid2spec(fid),dtype=torch.complex64)[None,:]*torch.tensor(concentrations)[:,None]
    cells=torch.matmul(labels.to(torch.complex64).reshape(-1,labels.shape[-1]),spectra)
    cells=cells.reshape(*labels.shape[:3],basis.n)
    cells=model.apply_shim_map(cells,model.shim_map)
    clean=spec2fid((cells*torch.prod(torch.tensor(model.simulation_params['voxel_spacing']))).sum((0,1,2)))
    clean=clean.detach().cpu().numpy()
    stream=np.random.SeedSequence(seed,spawn_key=(index,))
    noise=complex_gaussian_noise(basis.n,sigma=options['noise_sigma'],seed=stream)
    noisy=clean if options['noise_sigma']==0 else clean+noise
    provenance={**metadata,'water_noise':{'sigma_per_FID_component':options['noise_sigma'],'entropy':stream.entropy,'spawn_key':list(stream.spawn_key)}}
    return dict(water_reference_fid=noisy,clean_water_reference_fid=clean,
                water_reference_spectrum=torch.as_tensor(fid2spec(noisy)),water_reference_provenance=provenance)


def export_reference(result,output_dir=None):
    meta=result['nifti_metadata'];provenance=result['water_reference_provenance']
    obj=gen_nifti_mrs(result['water_reference_fid'][None,None,None,:,None],meta['dwelltime'],meta['spec_freq'],affine=np.asarray(meta['affine']),dim_tags=['DIM_DYN',None,None])
    obj.add_hdr_field('SpecFreqChemShift',provenance['frequency_origin_ppm'])
    obj.add_hdr_field('EchoTime',provenance['sequence']['TE_ms']/1000)
    obj.add_hdr_field('RepetitionTime',meta['sim_params']['simulation_params']['TR']/1000)
    if output_dir is not None:
        output=Path(output_dir)
        for existing in output.glob('*.nii.gz'):
            if existing.name.startswith('water_reference'):continue
            acquisition=mrs_io.read_FID(str(existing))
            acquisition.add_hdr_field('SpecFreqChemShift',provenance['frequency_origin_ppm'])
            acquisition.save(str(existing))
        obj.save(str(output/'water_reference.nii.gz'))
        (output/'ground_truth.json').write_text(json.dumps(result['ground_truth'],indent=2)+'\n')
        (output/'simulation_config.json').write_text(json.dumps(meta['sim_params'],indent=2)+'\n')
        (output/'calibration_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
        if provenance['water_options']['noise_sigma']:
            clean=gen_nifti_mrs(result['clean_water_reference_fid'][None,None,None,:,None],meta['dwelltime'],meta['spec_freq'],affine=np.asarray(meta['affine']),dim_tags=['DIM_DYN',None,None])
            for name in ['SpecFreqChemShift','EchoTime','RepetitionTime']:clean.add_hdr_field(name,obj.hdr_ext[name])
            clean.save(str(output/'water_reference_clean.nii.gz'))
    return obj

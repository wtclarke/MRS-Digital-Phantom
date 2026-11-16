"""Generate/cache the complete PRESS TE30 library and its shared-scale water."""

import argparse
import csv
import json
from pathlib import Path
import shlex
import subprocess

import fsl_mrs
import numpy as np
from scipy.io import loadmat
from fsl_mrs.utils import mrs_io

from examples.mrscloud_sequence_matched_water import write_basis
from examples.mrscloud_water_processing_validation import second_dc
from examples.mrscloud_naa_response_decomposition import first_dc
from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral
from simulation.mrscloud_quantitative_pair import REFERENCE_BASIS, export_pair
from simulation.mrscloud_full_basis import CONVENTION, STATISTIC, sha256, normalization_identity, load_library, synthesize

MRS_SOURCE=Path('simulation/basissets/MRSCloud-main')
PHASE5G=Path('outputs/mrscloud_sequence_water_integration/validated/reference')


def json_write(path,content):
    Path(path).write_text(json.dumps(content,indent=2,sort_keys=True,
        default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item(),allow_nan=False)+'\n')


def prepare_capture(directory, smoke=False):
    root=Path(directory);root.mkdir(parents=True,exist_ok=True)
    basis=mrs_io.read_basis(str(REFERENCE_BASIS));spins=loadmat(MRS_SOURCE/'functions/FID-A/simulationTools/metabolites/spinSystems.mat',simplify_cells=True)
    names=list(map(str,basis.names))
    missing=[n for n in names if 'sys'+n not in spins]
    if missing:raise ValueError('unmapped checked-in spin systems: '+str(missing))
    entries=[(n,n,3.) for n in names]+[('water_reference','H2O',4.7),
        ('methyl_control','NAA_Acetyl_only',3.),('methyl_on_carrier','NAA_Acetyl_only',2.008)]
    entries += [('repeat_'+n,n,3.) for n in ['NAA','Cr','Glu']]+[('water_repeat','H2O',4.7)]
    sources={str(p):sha256(p) for p in sorted(MRS_SOURCE.rglob('*')) if p.is_file() and p.suffix in ['.m','.mat','.pta']}
    sources['examples/full_basis_responses.m']=sha256('examples/full_basis_responses.m')
    cfg={'spatial_points':3 if smoke else 101,'native_points':256 if smoke else 8192,'source_hashes':sources}
    import hashlib
    signature=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest()
    cfg['signature']=signature
    config_path=root/'capture_configuration.json'
    if config_path.exists() and json.loads(config_path.read_text())!=cfg:
        raise ValueError('capture directory belongs to different source/configuration; use a new directory')
    json_write(config_path,cfg)
    text=''.join(f'{a} {b} {c}\n' for a,b,c in entries)
    entry_file=root/'entries.txt'
    if entry_file.exists() and entry_file.read_text()!=text:raise ValueError('capture entry list changed')
    entry_file.write_text(text)
    command=['pixi','run','octave','--no-gui','--quiet','--eval',
        "addpath('examples'); full_basis_responses('%s','%s',%d,%d,'%s');" %
        (str((root/'raw').resolve()).replace("'","''"),str(entry_file.resolve()).replace("'","''"),
         cfg['spatial_points'],cfg['native_points'],signature)]
    return cfg,names,entries,spins,command


def build(directory, smoke=False, reuse_raw=False):
    root=Path(directory)
    cfg,names,entries,spins,command=prepare_capture(root,smoke)
    if not reuse_raw:
        try:result=subprocess.run(command,capture_output=True,text=True)
        except FileNotFoundError as exc:raise RuntimeError('Pixi/Octave is required to generate missing raw responses') from exc
        (root/'octave.log').write_text(result.stdout+result.stderr)
        (root/'octave_command.txt').write_text(shlex.join(command)+'\n')
        result.check_returncode()
    raw={}
    for label,metabolite,carrier in entries:
        path=root/'raw'/label/'response.mat'
        if not path.exists():raise FileNotFoundError(f'missing raw response {label}; run without --reuse-raw')
        record=loadmat(path,simplify_cells=True)
        if (record['signature']!=cfg['signature'] or record['metabolite']!=metabolite or record['carrier']!=carrier
                or len(record['raw_fid'])!=cfg['native_points']):
            raise ValueError('raw cache source/configuration mismatch: '+label)
        seq=record['sequence']
        if ((seq['TE_ms'],seq['TE1_ms'],seq['TE2_ms'],seq['sw'],seq['readout_linewidth_hz'])!=(30,14,16,4000,1)
                or seq['spatial_points']!=cfg['spatial_points'] or seq['refocusing_duration_ms']!=7):
            raise ValueError('raw sequence mismatch: '+label)
        np.testing.assert_allclose(record['after_first_dc'],first_dc(record['raw_fid'],record['raw_ppm']),rtol=1e-11,atol=1e-12)
        raw[label]=record
    out=root/'mrscloud_press_te30_quantitative'
    if out.exists():raise FileExistsError('bundle already exists; validate/load it without regeneration')
    out.mkdir()
    native_points=cfg['native_points'];target_points=128 if smoke else 2048
    processed=np.column_stack([second_dc(raw[name]['after_first_dc'],4000,127.731,3.)[0] for name in names])
    spectra=np.fft.fftshift(np.fft.fft(processed,axis=0),axes=0).real
    scalar=float(np.max(spectra))
    if not np.isfinite(scalar) or scalar<=0:raise ValueError('invalid complete-reference normalization')
    t=np.arange(native_points)/4000
    unit=processed/scalar*np.exp(2j*np.pi*(3.-4.68)*127.731*t)[:,None]
    write_basis(out/'metabolite_basis.BASIS',dict(zip(names,unit.T)),{'Bfield':3.,'sw':4000.,'points':native_points})
    water=(raw['water_reference']['raw_fid']*np.exp(2j*np.pi*(4.7-4.68)*127.731*t)/scalar)[2*np.arange(target_points)]
    np.savez(out/'water_reference.npz',fid=water,raw_native=raw['water_reference']['raw_fid'])
    np.savez(out/'normalization_reference.npz',processed_reference=processed)
    seq={'sequence':'Universal_Philips_UnEdited_PRESS','TE_ms':30,'TE1_ms':14,'TE2_ms':16,
        'metabolite_RF_carrier_ppm':3.,'water_RF_carrier_ppm':4.7,'frequency_origin_ppm':4.68,
        'spectrometer_frequency_mhz':127.731,'native_points':native_points,'native_bandwidth_hz':4000,
        'target_points':target_points,'target_bandwidth_hz':2000,'dwell_time_s':.0005,
        'readout_linewidth_hz':1,'spatial_points_per_axis':cfg['spatial_points'],
        'refocusing_duration_ms':7,'refocusing_waveform':'univ_eddenrefo.pta',
        'refocusing_flip_angle_deg':180,'excitation_flip_angle_deg':90,
        'refocusing_FOV_cm':[4.5,4.5],'refocusing_thickness_cm':[3.,3.],
        'spatial_signal_scale':(4.5/3.)**2/cfg['spatial_points']**2,
        'unit_response_concentration_convention':'nominal mol/L per total tissue volume',
        'excitation':'ideal sim_excite x; bundled two-direction refocusing localization',
        'metabolite_processing':'first DC; second DC; shared scalar; frame shift; LCModel serialization; FSL formatting',
        'water_processing':'raw_pre_DC_shared_scale_frame_shift_exact_samples',
        'water_proton_multiplier_applied':False,'resolution_family':'smoke' if smoke else 'production',
        'physical_calibration':False}
    json_write(out/'sequence_config.json',seq)
    processing_sources={p:sha256(p) for p in ['examples/mrscloud_full_basis_water_reference.py',
        'examples/mrscloud_sequence_matched_water.py','examples/mrscloud_water_processing_validation.py',
        'examples/mrscloud_naa_response_decomposition.py','simulation/mrscloud_full_basis.py']}
    norm={'convention':CONVENTION,'statistic':STATISTIC,'scalar':scalar,'reference_names':names,
        'sequence_config_sha256':sha256(out/'sequence_config.json'),'source_hashes':cfg['source_hashes'],
        'processing_source_hashes':processing_sources,'numpy_version':np.__version__,
        'pixi_lock_sha256':sha256('pixi.lock'),
        'capture_configuration_signature':cfg['signature'],'fsl_mrs_version':fsl_mrs.__version__,
        'reference_array_sha256':sha256(out/'normalization_reference.npz')}
    norm['calibration_identity']=normalization_identity(norm)
    json_write(out/'normalization.json',norm)
    manifest={'entries':[{'name':n,'source_spin_system':'sys'+n,
        'role':'off_carrier_water_control' if n=='H2O' else 'reference_control' if n=='Ref0ppm' else 'metabolite',
        'RF_carrier_ppm':3.,'spin_system_definition':spins['sys'+n]} for n in names],
        'unsupported_entries':[],'excluded_MM_lipid_entries':[],
        'available_source_spin_systems':sorted(k for k in spins if k.startswith('sys')),
        'non_target_source_spin_systems':sorted(k for k in spins if k.startswith('sys') and k not in ['sys'+n for n in names]),
        'water_source_spin_system':'sysH2O','water_RF_carrier_ppm':4.7,
        'source_checked_basis':{'path':str(REFERENCE_BASIS),'sha256':sha256(REFERENCE_BASIS)},
        'source_hashes':cfg['source_hashes'],'processing_source_hashes':processing_sources,
        'raw_response_hashes':{label:sha256(root/'raw'/label/'response.mat') for label,_,_ in entries},
        'files':{n:sha256(out/n) for n in ['metabolite_basis.BASIS','water_reference.npz','normalization_reference.npz',
            'normalization.json','sequence_config.json']}}
    json_write(out/'reference_manifest.json',manifest)
    responses,loaded_water,metadata=load_library(out,points=target_points,selected_names=names)
    basis=mrs_io.read_basis(str(out/'metabolite_basis.BASIS'))
    formatted=basis.get_formatted_basis(2000,target_points)
    from fsl_mrs.utils.misc import ts_to_ts_ft
    direct=np.column_stack([ts_to_ts_ft(unit[:,i,None],.00025,.0005,target_points).ravel() for i in range(len(names))])
    serialization={n:float(np.linalg.norm(formatted[:,i]-direct[:,i])/np.linalg.norm(direct[:,i])) for i,n in enumerate(names)}
    processing_effects={}
    for i,n in enumerate(names):
        r=raw[n]['raw_fid'];a=raw[n]['after_first_dc'];b=processed[:,i]
        processing_effects[n]={'raw_initial_complex':str(r[0]),'first_DC_initial_complex':str(a[0]),
            'both_DC_initial_complex':str(b[0]),'formatted_unit_initial_complex':str(formatted[0,i]),
            'first_DC_initial_magnitude_ratio':float(abs(a[0])/abs(r[0])),
            'second_DC_initial_magnitude_ratio':float(abs(b[0])/abs(a[0])),
            'first_DC_relative_complex_change':float(np.linalg.norm(a-r)/np.linalg.norm(r)),
            'second_DC_relative_complex_change':float(np.linalg.norm(b-a)/np.linalg.norm(a))}
    determinism={n:{'raw_max_error':float(np.max(abs(raw[n]['raw_fid']-raw['repeat_'+n]['raw_fid']))),
        'post_DC_max_error':float(np.max(abs(raw[n]['after_first_dc']-raw['repeat_'+n]['after_first_dc'])))} for n in ['NAA','Cr','Glu']}
    determinism['water']={'raw_max_error':float(np.max(abs(raw['water_reference']['raw_fid']-raw['water_repeat']['raw_fid'])))}
    proton={'RF3_methyl_per_proton_vs_water':float(abs(raw['methyl_control']['raw_fid'][0])/3/(abs(raw['water_reference']['raw_fid'][0])/2)),
        'on_resonance_methyl_per_proton_vs_water':float(abs(raw['methyl_on_carrier']['raw_fid'][0])/3/(abs(raw['water_reference']['raw_fid'][0])/2))}
    comparisons=[];cross=None
    if not smoke:
        old=mrs_io.read_basis(str(REFERENCE_BASIS));old_array=old.get_formatted_basis(2000,2048)
        index=names.index('NAA');common=np.vdot(formatted[:,index],old_array[:,index])/np.vdot(formatted[:,index],formatted[:,index])
        ppm=np.fft.fftshift(np.fft.fftfreq(2048,.0005))/127.731+4.68
        for i,n in enumerate(names):
            new=common*formatted[:,i];historic=old_array[:,i]
            corr=np.vdot(historic,new)/(np.linalg.norm(historic)*np.linalg.norm(new))
            jold=spectral_integral(historic,ppm,(float(min(ppm)),float(max(ppm))))
            jnew=spectral_integral(new,ppm,(float(min(ppm)),float(max(ppm))))
            peak=float(ppm[np.argmax(abs(np.fft.fftshift(np.fft.fft(new))))])
            oldspec=abs(np.fft.fftshift(np.fft.fft(historic)))
            oldpeak=float(ppm[np.argmax(oldspec)])
            # Near-degenerate peaks can exchange argmax after text rounding.
            candidates=ppm[oldspec>=np.max(oldspec)*(1-1e-6)]
            comparisons.append({'name':n,'correlation_real':float(corr.real),'correlation_imag':float(corr.imag),
                'relative_complex_error':float(np.linalg.norm(new-historic)/np.linalg.norm(historic)),
                'initial_complex_ratio':str(new[0]/historic[0]),'spectral_integral_ratio':jnew/jold,
                'phase_difference_rad':float(np.angle(corr)),'dominant_peak_ppm':peak,'checked_peak_ppm':oldpeak,
                'dominant_peak_difference_ppm':peak-oldpeak,
                'distance_to_checked_near_max_ppm':float(np.min(abs(candidates-peak)))})
        oldnorm=json.loads((PHASE5G/'reference.json').read_text())['normalization_scalar']
        oldbasis=mrs_io.read_basis(str(PHASE5G/'common_response.BASIS'))
        oldN=oldbasis.get_formatted_basis(2000,2048)[:,oldbasis.names.index('NAA')]
        oldraw=np.load(PHASE5G/'raw_water.npz')['raw_water']
        oldW=(oldraw*np.exp(2j*np.pi*(4.7-4.68)*127.731*t)/oldnorm)[2*np.arange(2048)]
        ratio=scalar/oldnorm
        cross={'new_to_Phase5g_amplitude_multiplier':ratio,'NAA_max_complex_error':float(np.max(abs(responses['NAA']*ratio-oldN))),
            'NAA_relative_error':float(np.linalg.norm(responses['NAA']*ratio-oldN)/np.linalg.norm(oldN)),
            'water_max_complex_error':float(np.max(abs(loaded_water*ratio-oldW))),
            'water_relative_error':float(np.linalg.norm(loaded_water*ratio-oldW)/np.linalg.norm(oldW)),
            'shared_checked_basis_comparison_scalar':str(common)}
    validation={'calibration_identity':norm['calibration_identity'],'scalar':scalar,
        'maximum_reference_entry':names[np.unravel_index(np.argmax(spectra),spectra.shape)[1]],
        'determinism':determinism,'proton_controls':proton,'cross_calibration':cross,
        'basis_serialization_relative_errors':serialization,'checked_basis_comparison':comparisons,
        'processing_effects':processing_effects}
    # Exercise independent selections for three metabolites and a real mixture.
    selection={}
    for n in ['NAA','Cr','Glu']:
        alone,_,_=load_library(out,points=target_points,selected_names=[n])
        selection[n]=float(np.max(abs(alone[n]-responses[n])))
    mixture={'NAA':10.,'Cr':8.,'GPC':2.,'PCh':1.,'Glu':12.}
    pair=synthesize(out,mixture,points=target_points)
    expected=sum(responses[n]*(c*.001) for n,c in mixture.items())
    validation['selection_max_errors']=selection
    validation['mixture_concentrations_mM']=mixture
    validation['mixture_superposition_max_error']=float(np.max(abs(expected-pair['metabolite_fid'])))
    paths=export_pair(pair,root/'mixture_example')
    validation['mixture_export_max_error']=max(float(np.max(abs(mrs_io.read_FID(str(p)).mrs().FID-pair[k]))) for k,p in paths.items())
    json_write(out/'validation_results.json',validation)
    if comparisons:
        with (root/'checked_basis_comparison.csv').open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=comparisons[0]);writer.writeheader();writer.writerows(comparisons)
    if cfg['source_hashes']!={str(p):sha256(p) for p in cfg['source_hashes']}:
        raise RuntimeError('source changed during experiment')
    print(json.dumps({k:v for k,v in validation.items() if k not in ['checked_basis_comparison','basis_serialization_relative_errors','processing_effects']},indent=2))
    return out,validation


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir');parser.add_argument('--smoke',action='store_true');parser.add_argument('--reuse-raw',action='store_true')
    args=parser.parse_args();build(args.output_dir,args.smoke,args.reuse_raw)

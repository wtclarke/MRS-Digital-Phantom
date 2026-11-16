"""Isolate DC, formatting, and acquisition sampling of sequence-matched water."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

import numpy as np
from scipy.io import loadmat
from fsl_mrs.core import MRS
from fsl_mrs.utils import mrs_io
from fsl_mrs.utils.misc import FIDToSpec, ts_to_ts_ft

from examples.fsl_mrs_water_fit_diagnostics import fit_water
from examples.fsl_mrs_water_preprocessing import cli_command
from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral
from simulation.mrscloud_quantitative_pair import export_pair

ROOT=Path(__file__).resolve().parents[1]


def capture_raw(directory, spatial_points=101, points=8192):
    directory=Path(directory).resolve();directory.mkdir(parents=True,exist_ok=True)
    quote=lambda value:str(value).replace("'","''")
    command=['pixi','run','octave','--no-gui','--quiet','--eval',
             f"addpath('{quote(ROOT/'examples')}'); capture_water_processing('{quote(directory)}',{int(spatial_points)},{int(points)});"]
    result=subprocess.run(command,capture_output=True,text=True)
    (directory/'octave.log').write_text(result.stdout+result.stderr)
    (directory/'command.txt').write_text(shlex.join(command)+'\n')
    result.check_returncode()
    return loadmat(directory/'processing_trace.mat',simplify_cells=True)


def second_dc(fid, sw, cf, carrier):
    spectrum=np.fft.fftshift(np.fft.fft(fid))
    ppm=(np.arange(len(fid))*sw/len(fid)-sw/2+sw/(2*len(fid)))/cf+carrier
    a=(ppm>max(ppm)-.5)&(ppm<max(ppm));b=(ppm>min(ppm))&(ppm<min(ppm)+.5)
    offset=(np.mean(spectrum[a])+np.mean(spectrum[b]))/2
    return np.fft.ifft(np.fft.ifftshift(spectrum-offset)),complex(offset)


def sample_on_acquisition_grid(fid, old_dwell, new_dwell, points):
    ratio=new_dwell/old_dwell
    step=int(round(ratio))
    if step<1 or not np.isclose(ratio,step,rtol=0,atol=1e-12) or (points-1)*step>=len(fid):
        raise ValueError('requires covered acquisition times on an integer native-grid stride')
    return np.asarray(fid)[np.arange(points)*step].copy()


def response_measurements(water, naa, raw_reference, dwell, cf, linewidth):
    t=np.arange(len(water))*dwell
    selected=(t>=.05)&(t<=.5)
    if sum(selected)<3:selected=t>0
    log_fit=np.polyfit(t[selected],np.log(abs(water[selected])),1)
    phase_fit=np.polyfit(t[selected],np.unwrap(np.angle(water))[selected],1)
    ppm=np.fft.fftshift(np.fft.fftfreq(len(water),dwell))/cf+4.68
    return {'initial_complex':str(water[0]),'initial_magnitude':float(abs(water[0])),
            'tail_extrapolated_initial_magnitude':float(np.exp(log_fit[1])),
            'initial_phase_rad':float(np.angle(water[0])),
            'tail_frequency_hz':float(phase_fit[0]/(2*np.pi)),
            'tail_linewidth_hz':float(-log_fit[0]/np.pi),
            'max_complex_difference_from_raw_grid':float(np.max(abs(water-raw_reference))),
            'relative_complex_difference_from_raw_grid':float(np.linalg.norm(water-raw_reference)/np.linalg.norm(raw_reference)),
            'spectral_integral_per_mol_l':spectral_integral(water,ppm,(1.65,7.65)),
            'water_to_NAA_initial_complex_ratio':str(water[0]/naa[0]),
            'water_to_NAA_initial_magnitude_ratio':float(abs(water[0]/naa[0])),
            'points':len(water),'dwell_time_s':dwell,'configured_linewidth_hz':linewidth}


def run_validation(output_dir,source_dir='outputs/mrscloud_sequence_matched_water',reuse_capture=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out,source=Path(output_dir),Path(source_dir)
    out.mkdir(parents=True,exist_ok=True)
    trace=(loadmat(out/'octave/processing_trace.mat',simplify_cells=True) if reuse_capture else capture_raw(out/'octave'))
    previous=loadmat(source/'octave/raw_sequence.mat',simplify_cells=True)
    # Independently captured samples must reproduce Phase 5e before varying processing.
    np.testing.assert_array_equal(trace['after_first_dc'],previous['signals']['water_on'])
    saved=json.loads((source/'validated/results.json').read_text())
    scalar=saved['configuration']['shared_max_real_spectrum_scalar']
    cf=127.731;sw=float(trace['sw']);dt=1/sw
    mix=np.exp(2j*np.pi*(4.7-4.68)*cf*np.arange(int(trace['points']))*dt)
    a=np.asarray(trace['raw_fid'])*mix/scalar
    b=np.asarray(trace['after_first_dc'])*mix/scalar
    c_unshifted,offset=second_dc(trace['after_first_dc'],sw,cf,4.7)
    c=c_unshifted*mix/scalar
    basis_path=source/'validated/common_response.BASIS'
    basis=mrs_io.read_basis(str(basis_path))
    native_naa=basis.original_basis_array[:,basis.names.index('NAA')]
    formatted=basis.get_formatted_basis(2000,2048)
    final_naa=formatted[:,basis.names.index('NAA')]
    d=formatted[:,basis.names.index('water_on')]
    e=sample_on_acquisition_grid(a,dt,.0005,2048)
    f=ts_to_ts_ft(a[:,None],dt,.0005,2048).ravel()
    d_direct=ts_to_ts_ft(c[:,None],dt,.0005,2048).ravel()
    variants={'A_raw':(a,native_naa,a,dt),'B_first_DC':(b,native_naa,a,dt),
              'C_both_DC':(c,native_naa,a,dt),'D_DC_and_Fourier_resampling':(d,final_naa,e,.0005),
              'E_raw_acquisition_grid':(e,final_naa,e,.0005),
              'F_raw_Fourier_control':(f,final_naa,e,.0005)}
    results={};fig,axes=plt.subplots(3,2,figsize=(12,11));commands=[]
    for ax,(name,(water,naa,reference,dwell)) in zip(axes.flat,variants.items()):
        mrs=MRS(FID=.010*naa,H2O=40*water,cf=cf,bw=1/dwell,nucleus='1H',chemShift=4.68)
        mrs.rescaleForFitting()
        fit_info,integrator,model=fit_water(mrs)
        measured=response_measurements(water,naa,reference,dwell,cf,float(trace['lw']))
        measured['fitted_water_integral_per_mol_l']=float(integrator.integral/mrs.scaling['FID']/40)
        results[name]={'response':measured,'water_fit':fit_info}
        np.savez(out/f'{name}.npz',water=water,raw_reference=reference,fixed_NAA=naa,
                 water_fit=model/mrs.scaling['FID']/40,fitted_integration_model=integrator.fid/mrs.scaling['FID']/40)
        t=np.arange(len(water))*dwell
        select=t<.012
        ax.plot(t[select],water.real[select],'.-',label='real data')
        ax.plot(t[select],model.real[select]/mrs.scaling['FID']/40,'--',label='actual water fit')
        ax.set_title(f'{name}; residual {fit_info["relative_complex_residual"]:.3g}',fontsize=9)
        ax.set_xlabel('Time (s)');ax.set_ylabel('SRU per mol/L');ax.legend(fontsize=8)
        if name in ['D_DC_and_Fourier_resampling','E_raw_acquisition_grid']:
            directory=out/name
            config={**saved['configuration'],'dwell_time_s':dwell,
                    'water_processing_variant':name,'source_response_scale':scalar}
            pair={'metabolite_fid':.010*naa,'water_reference_fid':40*water,
                  'simulation_config':config,
                  'ground_truth':{'NAA_mM':10,'water_mol_l':40,'physical_calibration':False}}
            export_pair(pair,directory)
            command=cli_command(directory);command[command.index('--basis')+1]=str(basis_path)
            commands.append(shlex.join(command))
            process=subprocess.run(command,capture_output=True,text=True)
            (directory/'cli.log').write_text(process.stdout+process.stderr);process.check_returncode()
            report=directory/'fsl_mrs_results/report.html'
            if not report.is_file() or report.stat().st_size==0:raise RuntimeError('CLI HTML missing')
            with (directory/'fsl_mrs_results/concentrations.csv').open() as stream:records=list(csv.reader(stream))
            values=dict(zip(records[1][1:5],map(float,records[2][1:5])))
            results[name]['CLI_concentrations']=values
            reloaded=mrs_io.read_FID(str(directory/'water_reference.nii.gz')).mrs().FID
            results[name]['export_max_complex_error']=float(np.max(abs(reloaded-40*water)))
    fig.tight_layout();fig.savefig(out/'early_time_processing.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    ppm=np.fft.fftshift(np.fft.fftfreq(2048,.0005))/cf+4.68;select=(ppm>4.5)&(ppm<4.9)
    for name in ['D_DC_and_Fourier_resampling','E_raw_acquisition_grid']:
        arrays=np.load(out/f'{name}.npz')
        axes[0].plot(ppm[select],FIDToSpec(arrays['water'].copy()).real[select],label=name)
        axes[1].plot(np.arange(8)*.0005,arrays['water'][:8].real,'o-',label=name)
    axes[0].invert_xaxis();axes[0].set_xlabel('ppm');axes[1].set_xlabel('Time (s)')
    for ax in axes:ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(out/'preferred_vs_basis_processing.png',dpi=160);plt.close(fig)
    (out/'cli_commands.txt').write_text('\n'.join(commands)+'\n')
    with (out/'comparison.csv').open('w') as stream:
        rows=[{'variant':name,**value['response'],**{f'fit_{k}':v for k,v in value['water_fit'].items() if isinstance(v,(int,float,bool))}} for name,value in results.items()]
        writer=csv.DictWriter(stream,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    summary={'source_scalar':scalar,'source_normalization_changed':False,'variants':results,
             'capture_vs_phase5e_first_DC_max_error':float(np.max(abs(trace['after_first_dc']-previous['signals']['water_on']))),
             'DC2_complex_offset':str(offset),
             'fourier_roundtrip_no_DC_max_error':float(np.max(abs(np.fft.fft(np.fft.ifft(trace['raw_fid']))-trace['raw_fid']))),
             'LCModel_serialization_then_formatting_relative_difference':float(np.linalg.norm(d-d_direct)/np.linalg.norm(d_direct)),
             'preferred_initial_vs_raw_initial_error':float(abs(e[0]-a[0])),
             'native_vs_final_raw_common_time_max_error':float(np.max(abs(e-a[np.arange(2048)*2]))),
             'metabolite_array_difference_D_vs_E':float(np.max(abs(variants['D_DC_and_Fourier_resampling'][1]-variants['E_raw_acquisition_grid'][1]))),
             'instrumented_original_source_sha256':hashlib.sha256((ROOT/'simulation/basissets/MRSCloud-main/functions/sim_signals.m').read_bytes()).hexdigest()}
    (out/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir');parser.add_argument('--reuse-capture',action='store_true')
    parser.add_argument('--source-dir',default='outputs/mrscloud_sequence_matched_water')
    args=parser.parse_args()
    run_validation(args.output_dir,args.source_dir,args.reuse_capture)

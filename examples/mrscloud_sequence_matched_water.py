"""Isolated MRSCloud/Octave water-carrier feasibility experiment."""

import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

import numpy as np
import fsl_mrs
from scipy.io import loadmat
from fsl_mrs.utils import mrs_io, fitting, quantify
from fsl_mrs.utils.misc import FIDToSpec

from examples.fsl_mrs_water_fit_diagnostics import fit_water
from examples.fsl_mrs_water_preprocessing import cli_command
from simulation.mrscloud_quantitative_pair import REFERENCE_BASIS, export_pair
from simulation.frozen_basis_normalization import format_frozen_basis

ROOT=Path(__file__).resolve().parents[1]


def simulate_raw(output_dir, spatial_points=101, points=8192):
    directory=Path(output_dir).resolve()
    directory.mkdir(parents=True,exist_ok=True)
    quote=lambda path: str(path).replace("'","''")
    command=['pixi','run','octave','--no-gui','--quiet','--eval',
             f"addpath('{quote(ROOT/'examples')}'); sequence_matched_water('{quote(directory)}',{int(spatial_points)},{int(points)});"]
    completed=subprocess.run(command,capture_output=True,text=True)
    (directory/'octave.log').write_text(completed.stdout+completed.stderr)
    (directory/'command.txt').write_text(shlex.join(command)+'\n')
    completed.check_returncode()
    return loadmat(directory/'raw_sequence.mat',simplify_cells=True)


def assemble_common(raw):
    """Reproduce assembly's second DC correction and one shared maximum scale."""
    points=int(raw['points']);sw=float(raw['sw']);cf=float(raw['Bfield'])*42.577
    t=np.arange(points)/sw
    shifted={};processed={}
    carriers=dict(zip(raw['signals'],raw['carriers']))
    native_f=np.arange(points)*sw/points-sw/2+sw/(2*points)
    for name,fid in raw['signals'].items():
        fid=np.asarray(fid,dtype=complex).reshape(-1)
        spec=np.fft.fftshift(np.fft.fft(fid))
        ppm=native_f/cf+carriers[name]
        a=(ppm>max(ppm)-.5)&(ppm<max(ppm))
        b=(ppm>min(ppm))&(ppm<min(ppm)+.5)
        dc=(np.mean(spec[a])+np.mean(spec[b]))/2
        processed[name]=np.fft.ifft(np.fft.ifftshift(spec-dc))
    reference_names=['NAA','Ref0ppm','water_off']
    scalar=float(max(np.max(np.fft.fftshift(np.fft.fft(processed[n])).real) for n in reference_names))
    if not np.isfinite(scalar) or scalar<=0:
        raise ValueError('reference maximum must be finite and positive')
    for name,fid in processed.items():
        # Acquisition frequency shift AFTER the carrier-sensitive sequence.
        shifted[name]=fid/scalar*np.exp(-2j*np.pi*(4.68-carriers[name])*cf*t)
    return shifted,scalar,carriers,processed


def write_basis(path,signals,raw):
    """LCModel spectral serialization equivalent to bundled export, common scale."""
    cf=float(raw['Bfield'])*42.577
    with Path(path).open('w') as stream:
        stream.write(f" $SEQPAR\n HZPPPM = {cf:.9f},\n ECHOT = 30.,\n SEQ = 'UnEdited'\n $END\n")
        stream.write(f" $BASIS1\n IDBASI = 'sequence matched feasibility',\n BADELT = {1/float(raw['sw']):.12f},\n NDATAB = {raw['points']}\n $END\n")
        for name,fid in signals.items():
            stream.write(f" $BASIS\n ID = '{name}',\n METABO = '{name}',\n CONC = 1.,\n TRAMP = 1.,\n VOLUME = 1.,\n ISHIFT = 0\n $END\n")
            # Bundled exporter FFTs conjugated FID; FSL reader reverses this.
            values=np.fft.fft(fid.conj())
            np.savetxt(stream,np.column_stack([values.real,values.imag]),fmt='%.6e')


def run_experiment(output_dir,reuse_raw=False,spatial_points=101,raw_dir=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    raw=(loadmat(Path(raw_dir or out/'octave')/'raw_sequence.mat',simplify_cells=True) if reuse_raw
         else simulate_raw(out/'octave',spatial_points=spatial_points))
    native,scalar,carriers,processed=assemble_common(raw)
    basis_path=out/'common_response.BASIS'
    write_basis(basis_path,native,raw)
    basis=mrs_io.read_basis(str(basis_path))
    formatted=basis.get_formatted_basis(2000,2048)
    signals={str(name):formatted[:,i] for i,name in enumerate(basis.names)}
    t=np.arange(2048)/2000;cf=float(raw['Bfield'])*42.577
    q=abs(signals['Ref0ppm'][0])
    analytical=2*q*np.exp(-np.pi*float(raw['lw'])*t)*np.exp(2j*np.pi*(4.7-4.68)*cf*t)
    common={'bandwidth_hz':2000,'points':2048,'dwell_time_s':.0005,
            'spectrometer_frequency_mhz':cf,'frequency_origin_ppm':4.68,
            'TE_ms':30,'TE1_ms':float(raw['TE1']),'TE2_ms':float(raw['TE2']),
            'refocusing_duration_ms':float(raw['refTp']),'readout_linewidth_hz':float(raw['lw']),
            'RF_carriers_ppm':carriers,'shared_max_real_spectrum_scalar':scalar,
            'normalization_reference_names':['NAA','Ref0ppm','water_off'],
            'spatial_points_per_axis':int(raw['spatial_points']),
            'native_points':int(raw['points']),'native_bandwidth_hz':float(raw['sw']),
            'refocusing_waveform':'univ_eddenrefo.pta','refocusing_flip_angle_deg':180,
            'excitation':'ideal sim_excite x; loaded excitation waveform not applied in PRESS',
            'refocusing_FOV_cm':[4.5,4.5],'refocusing_thickness_cm':[3.,3.],
            'spatial_signal_scale':1/(int(raw['spatial_points'])**2)*(4.5/3.)**2,
            'normalization':'one common assembly maximum; no per-acquisition normalization',
            'physical_calibration':False}
    results={};commands=[]
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    analytical_zero_lw=2*q*np.exp(2j*np.pi*(4.7-4.68)*cf*t)
    for name,water_response in [('sequence_water',signals['water_on']),('analytical_reference',analytical),
                               ('analytical_reference_zero_lw',analytical_zero_lw)]:
        case=out/name
        pair={'metabolite_fid':.010*signals['NAA'],'water_reference_fid':40*water_response,
              'simulation_config':{**common,'water_response':name,
                                   'water_response_linewidth_hz':0. if name=='analytical_reference_zero_lw' else 1.},
              'ground_truth':{'NAA_mM':10,'water_mol_l':40,'physical_calibration':False}}
        export_pair(pair,case)
        command=cli_command(case)
        command[command.index('--basis')+1]=str(basis_path)
        commands.append(shlex.join(command))
        completed=subprocess.run(command,capture_output=True,text=True)
        (case/'cli.log').write_text(completed.stdout+completed.stderr)
        completed.check_returncode()
        if not (case/'fsl_mrs_results/report.html').is_file():
            raise RuntimeError('CLI report missing')
        mrs=mrs_io.read_FID(str(case/'metabolites.nii.gz')).mrs(basis=basis,
            ref_data=str(case/'water_reference.nii.gz'))
        mrs.keep=['NAA'];mrs.rescaleForFitting()
        fit=fitting.fit_FSLModel(mrs,method='Newton',model='lorentzian',baseline_order=-1,
                                ppmlim=(.2,4.2),metab_groups=[0],capture_minimize_output=True)
        water_info,reference,water_model=fit_water(mrs)
        quant=quantify.QuantificationInfo(.03,2.,['NAA'],cf,water_ref_metab='NAA',
                                        water_ref_metab_protons=3,water_ref_metab_limits=(1.8,2.2))
        quant.set_fractions({'WM':1.,'GM':0.,'CSF':0.})
        fit.calculateConcScaling(mrs,quant_info=quant,internal_reference=['NAA'])
        c=float(fit.getConc(scaling='molarity',metab='NAA'))
        import csv
        with (case/'fsl_mrs_results/concentrations.csv').open() as stream:
            records=list(csv.reader(stream))
        cli_concentrations=dict(zip(records[1][1:5],map(float,records[2][1:5])))
        if not np.isclose(cli_concentrations['molarity'],c,rtol=1e-12,atol=1e-12):
            raise ValueError('CLI/API quantification mismatch')
        qeff=quant.relax_corr_water_molar*quant.relax_corr_metab
        # Predict standard scaling from the measured input response integrals.
        ppm=mrs.getAxes()
        from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral
        j_m=spectral_integral(signals['NAA'],ppm,(1.8,2.2))
        water_unit_fit=reference.fid/mrs.scaling['FID']/40
        j_w=spectral_integral(water_unit_fit,ppm,(1.65,7.65))
        predicted=.010*j_m/(40*j_w)*(2/3)*qeff
        results[name]={'expected_initial_water':str(pair['water_reference_fid'][0]),
                       'input_initial_magnitude':float(abs(pair['water_reference_fid'][0])),
                       'fitted_water':water_info,'internal_NAA_coefficient':float(fit.getConc(metab='NAA')),
                       'raw_NAA_coefficient_mol_l':float(fit.getConc(scaling='raw',metab='NAA')),
                       'fsl_molarity_mM':c,'integral_predicted_mM':predicted,
                       'cli_molarity_mM':cli_concentrations['molarity'],
                       'metabolite_integral_per_mol_l':j_m,'water_integral_per_mol_l':j_w,
                       'effective_water_mM':qeff,
                       'NAA_relative_residual':float(np.linalg.norm(mrs.FID-fit.predictedFID(mrs))/np.linalg.norm(mrs.FID))}
        ppm=mrs.getAxes();selection=(ppm>4.4)&(ppm<5.)
        axes[0].plot(ppm[selection],FIDToSpec(pair['water_reference_fid'].copy()).real[selection],label=name)
        axes[1].plot(t,abs(water_response),label=name)
    axes[0].invert_xaxis();axes[0].set_xlabel('ppm');axes[1].set_xlabel('Time (s)')
    for ax in axes:ax.legend()
    fig.tight_layout();fig.savefig(out/'response_comparison.png',dpi=160);plt.close(fig)
    original,_,old_calibration=format_frozen_basis(REFERENCE_BASIS,bandwidth_hz=2000,points=2048,selected_names=['NAA','Ref0ppm'])
    alpha=np.vdot(signals['NAA'],original[:,0])/np.vdot(signals['NAA'],signals['NAA'])
    source_files=[ROOT/'examples/sequence_matched_water.m',ROOT/'simulation/basissets/MRSCloud-main/functions/sim_signals.m',
                  ROOT/'simulation/basissets/MRSCloud-main/functions/load_parameters.m',
                  ROOT/'simulation/basissets/MRSCloud-main/functions/osprey_functions/fit_makeBasis_MRSCloud.m',
                  ROOT/'simulation/basissets/MRSCloud-main/functions/osprey_functions/io_writelcmBASIS.m',
                  ROOT/'simulation/basissets/MRSCloud-main/functions/FID-A/simulationTools/sim_readout.m',
                  ROOT/'simulation/basissets/MRSCloud-main/functions/experim_funct/sim_Hamiltonian_mgs.m',
                  ROOT/'simulation/basissets/MRSCloud-main/functions/FID-A/simulationTools/metabolites/spinSystems.mat',
                  ROOT/'simulation/basissets/MRSCloud-main/pulses_universal/univ_eddenrefo.pta']
    measurements={name:{'initial_complex':str(signal[0]),'initial_magnitude':float(abs(signal[0])),
                         'initial_phase_rad':float(np.angle(signal[0])),
                         'frequency_hz_from_tail':float(np.polyfit(t[1:],np.unwrap(np.angle(signal[1:])),1)[0]/(2*np.pi))}
                  for name,signal in signals.items() if name!='NAA'}
    measurements['NAA']={'initial_complex':str(signals['NAA'][0])}
    tail_extrapolated_water = raw['signals']['water_on'][1]*np.exp(np.pi*float(raw['lw'])/float(raw['sw']))
    result={'FSL_MRS_version':fsl_mrs.__version__,
            'configuration':common,'fitting':results,'responses':measurements,
            'raw_water_before_assembly_initial':str(raw['signals']['water_on'][0]),
            'raw_water_tail_extrapolated_initial':str(tail_extrapolated_water),
            'raw_water_first_sample_vs_tail_extrapolation':float(abs(raw['signals']['water_on'][0]/tail_extrapolated_water)),
            'water_on_off_initial_amplitude_ratio':float(abs(signals['water_on'][0]/signals['water_off'][0])),
            'water_vs_same_scale_analytical_initial_ratio':float(abs(signals['water_on'][0])/(2*q)),
            'water_vs_same_scale_analytical_phase_difference_rad':float(np.angle(signals['water_on'][0])),
            'water_determinism_max_error':float(np.max(abs(raw['signals']['water_on']-raw['signals']['water_repeat']))),
            'two_proton_max_error':float(np.max(abs(raw['signals']['water_on']-2*raw['signals']['water_single_proton']))),
            'checked_basis_NAA_relative_difference_after_global_complex_match':float(np.linalg.norm(original[:,0]-alpha*signals['NAA'])/np.linalg.norm(original[:,0])),
            'checked_basis_Ref0ppm_relative_difference_after_same_global_match':float(np.linalg.norm(original[:,1]-alpha*signals['Ref0ppm'])/np.linalg.norm(original[:,1])),
            'checked_basis_global_complex_match':str(alpha),
            'estimated_historical_assembly_scale':float(scalar/(alpha.real*old_calibration.normalization_scalar)),
            'comparison_global_match_applied_to_water':False,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}}
    np.savez(out/'common_responses.npz',**signals,analytical=analytical,analytical_zero_lw=analytical_zero_lw)
    (out/'sequence_configuration.json').write_text(json.dumps(common,indent=2)+'\n')
    (out/'cli_commands.txt').write_text('\n'.join(commands)+'\n')
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir');parser.add_argument('--reuse-raw',action='store_true')
    parser.add_argument('--spatial-points',type=int,default=101)
    parser.add_argument('--raw-dir',type=Path)
    args=parser.parse_args()
    run_experiment(args.output_dir,args.reuse_raw,args.spatial_points,args.raw_dir)

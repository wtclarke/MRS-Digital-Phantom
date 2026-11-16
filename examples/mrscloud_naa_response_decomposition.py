"""Capture actual NAA spin responses and decompose processing/integration effects."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

import fsl_mrs
import numpy as np
from scipy.io import loadmat
from fsl_mrs.utils import mrs_io, fitting, quantify
from fsl_mrs.utils.misc import ts_to_ts_ft, FIDToSpec

from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral
from examples.mrscloud_water_processing_validation import second_dc, sample_on_acquisition_grid
from examples.fsl_mrs_water_fit_diagnostics import fit_water
from examples.fsl_mrs_water_preprocessing import cli_command
from examples.mrscloud_sequence_matched_water import write_basis
from simulation.mrscloud_quantitative_pair import generate_pair, export_pair

REFERENCE = Path('outputs/mrscloud_sequence_water_integration/validated/reference')


def capture(directory, spatial_points=101, points=8192):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    command = ['pixi','run','octave','--no-gui','--quiet','--eval',
        "addpath('examples'); capture_naa_response('%s',%d,%d);" %
        (str(directory).replace("'","''"), spatial_points, points)]
    result = subprocess.run(command, capture_output=True, text=True)
    (directory/'command.txt').write_text(shlex.join(command)+'\n')
    (directory/'octave.log').write_text(result.stdout+result.stderr)
    result.check_returncode()
    return loadmat(directory/'naa_trace.mat', simplify_cells=True)


def first_dc(raw, ppm):
    spectrum = np.fft.fftshift(np.fft.ifft(raw))
    edge1 = (ppm>max(ppm)-.5)&(ppm<max(ppm))
    edge2 = (ppm>min(ppm))&(ppm<min(ppm)+.5)
    offset = (np.mean(spectrum[edge1])+np.mean(spectrum[edge2]))/2
    return np.fft.fft(np.fft.ifftshift(spectrum-offset))


def stages(raw, after_first, *, carrier, scalar, sw=4000., points=2048):
    """DC in native carrier frame, then shared division and analysis-frame mix."""
    after_second, _ = second_dc(after_first, sw, 127.731, carrier)
    mix = np.exp(-2j*np.pi*(4.68-carrier)*127.731*np.arange(len(raw))/sw)
    native = {key: value/scalar*mix for key,value in
              [('raw',raw),('first_dc',after_first),('second_dc',after_second)]}
    target = {key+'_target': sample_on_acquisition_grid(value,1/sw,.0005,points)
              for key,value in native.items()}
    target['formatted'] = ts_to_ts_ft(native['second_dc'][:,None],1/sw,.0005,points).ravel()
    target['raw_fourier'] = ts_to_ts_ft(native['raw'][:,None],1/sw,.0005,points).ravel()
    return native, target


def measurements(fid, reference, sw, protons=None):
    ppm = np.fft.fftshift(np.fft.fftfreq(len(fid),1/sw))/127.731+4.68
    jm = spectral_integral(fid,ppm,(1.8,2.2))
    return {'initial_real':float(fid[0].real),'initial_imag':float(fid[0].imag),
            'initial_magnitude':float(abs(fid[0])),'initial_phase_rad':float(np.angle(fid[0])),
            'spectral_integral':jm,
            'initial_magnitude_per_proton':None if protons is None else float(abs(fid[0])/protons),
            'max_complex_difference_from_raw':float(np.max(abs(fid-reference))),
            'relative_complex_difference_from_raw':float(np.linalg.norm(fid-reference)/np.linalg.norm(reference)),
            'points':len(fid),'bandwidth_hz':sw}


def decompose(arrays, water):
    """Ordered controlled-response ratios; no factor is fitted to alpha."""
    time = np.arange(len(water))*.0005
    ppm = np.fft.fftshift(np.fft.fftfreq(len(water),.0005))/127.731+4.68
    jw = spectral_integral(abs(water[0])*np.exp(-np.pi*time),ppm,(1.65,7.65))
    # Hamiltonian and analysis-frame frequencies use different bundled gamma values.
    frequency = (2.008-3.)*127.728+(3.-4.68)*127.731
    ideal = 1.5*abs(water[0])*np.exp(-np.pi*time)*np.exp(2j*np.pi*frequency*time)
    raw = arrays['methyl_raw_target']
    amplitude_ratio = float(abs(raw[0])/(1.5*abs(water[0])))
    phase_zero = raw*np.exp(-1j*np.angle(raw[0]))
    area = lambda x: spectral_integral(x,ppm,(1.8,2.2))
    factor = {
        'ideal_finite_grid_window':area(ideal)/jw*(2/3),
        'methyl_sequence_amplitude_per_proton':amplitude_ratio,
        'methyl_line_shape_vs_analytical':area(phase_zero)/(amplitude_ratio*area(ideal)),
        'methyl_intrinsic_phase':area(raw)/area(phase_zero),
        'other_NAA_resonances':area(arrays['full_raw_target'])/area(raw),
        'first_DC':area(arrays['full_first_dc_target'])/area(arrays['full_raw_target']),
        'second_DC':area(arrays['full_second_dc_target'])/area(arrays['full_first_dc_target']),
        'Fourier_formatting':area(arrays['full_formatted'])/area(arrays['full_second_dc_target']),
    }
    final = arrays['full_formatted']
    if 'full_existing_formatted' in arrays:
        final = arrays['full_existing_formatted']
        factor['BASIS_serialization'] = area(final)/area(arrays['full_formatted'])
    product = float(np.prod(list(factor.values())))
    return {'factors':factor,'product':product,'direct_alpha':area(final)/jw*(2/3),
            'unit_water_integral':jw,'unit_formatted_NAA_integral':area(final),
            'methyl_frequency_hz':frequency,'water_to_methyl_initial_complex':str(water[0]/raw[0]),
            'ideal_3H_target':ideal,'water_integration_target':abs(water[0])*np.exp(-np.pi*time)}


def ideal_grid_controls():
    """Analytical controls only: distinguish finite record, windows and spin shift."""
    rows=[]
    frequency=(2.008-3.)*127.728+(3.-4.68)*127.731
    for points in [2048,4096,8192,16384,32768]:
        time=np.arange(points)/2000
        ppm=np.fft.fftshift(np.fft.fftfreq(points,.0005))/127.731+4.68
        water=2*np.exp(-np.pi*time)
        methyl=3*np.exp(-np.pi*time)*np.exp(2j*np.pi*frequency*time)
        ratio=spectral_integral(methyl,ppm,(1.8,2.2))/spectral_integral(water,ppm,(1.65,7.65))*2/3
        rows.append({'points':points,'dwell_time_s':.0005,'ideal_response_factor':ratio})
    return rows


def fit_control(directory, response, water_pair, basis_path):
    config={**water_pair['simulation_config'], 'diagnostic_metabolite_basis':str(basis_path),
            'diagnostic_basis_sha256':hashlib.sha256(Path(basis_path).read_bytes()).hexdigest(),
            'diagnostic_only':True}
    pair = {**water_pair, 'metabolite_fid':.010*response, 'simulation_config':config}
    export_pair(pair,directory)
    command = cli_command(directory)
    command[command.index('--basis')+1] = str(basis_path)
    completed = subprocess.run(command,capture_output=True,text=True)
    (directory/'cli.log').write_text(completed.stdout+completed.stderr)
    completed.check_returncode()
    if (directory/'fsl_mrs_results/report.html').stat().st_size<1000:
        raise RuntimeError('CLI report missing')
    mrs = mrs_io.read_FID(str(directory/'metabolites.nii.gz')).mrs(
        basis=mrs_io.read_basis(str(basis_path)),ref_data=str(directory/'water_reference.nii.gz'))
    mrs.keep=['NAA'];mrs.rescaleForFitting()
    fit=fitting.fit_FSLModel(mrs,method='Newton',model='lorentzian',baseline_order=-1,
        ppmlim=(.2,4.2),metab_groups=[0],capture_minimize_output=True)
    info,_,_=fit_water(mrs)
    q=quantify.QuantificationInfo(.03,2.,['NAA'],127.731,water_ref_metab='NAA',
        water_ref_metab_protons=3,water_ref_metab_limits=(1.8,2.2))
    q.set_fractions({'WM':1.,'GM':0.,'CSF':0.})
    fit.calculateConcScaling(mrs,quant_info=q,internal_reference=['NAA'])
    refs=fit.concScalings['ref_info']
    result={'raw_coefficient_mol_l':float(fit.getConc(scaling='raw',metab='NAA')),
        'standard_mM':float(fit.getConc(scaling='molarity',metab='NAA')),
        'water':info, 'water_integral':float(refs['water_ref'].integral),
        'metabolite_integral':float(refs['metab_ref'].integral),
        'metabolite_relative_residual':float(np.linalg.norm(mrs.FID-fit.predictedFID(mrs))/np.linalg.norm(mrs.FID)),
        'fit_parameters':dict(zip(fit.fitResults.columns,fit.fitResults.iloc[0].tolist())),
        'success':bool(fit.scipy_minimize_result.success)}
    (directory/'fitting_results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result,shlex.join(command)


def run(output_dir, reuse_capture=False):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    if (out/'results.json').exists():
        raise FileExistsError('experiment results already exist')
    source=Path('simulation/basissets/MRSCloud-main/functions')
    tracked=[source/'sim_signals.m',source/'load_parameters.m',
        source/'FID-A/simulationTools/metabolites/spinSystems.mat',
        source/'FID-A/processingTools/op_dccorr_FID_A.m',
        source/'FID-A/processingTools/op_dccorr.m',
        source/'FID-A/simulationTools/sim_readout.m',
        source/'experim_funct/sim_Hamiltonian_mgs.m',
        Path('examples/capture_naa_response.m'),
        REFERENCE/'common_response.BASIS',REFERENCE/'raw_water.npz',REFERENCE/'reference.json']
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in tracked}
    trace=loadmat(out/'octave/naa_trace.mat',simplify_cells=True) if reuse_capture else capture(out/'octave')
    manifest=json.loads((REFERENCE/'reference.json').read_text());scalar=manifest['normalization_scalar']
    old=loadmat('outputs/mrscloud_sequence_matched_water/octave/raw_sequence.mat',simplify_cells=True)
    np.testing.assert_array_equal(trace['first_dc']['full'],old['signals']['NAA'])
    pair=generate_pair({'NAA':10},water_reference_mode='mrscloud_sequence',sequence_reference=REFERENCE)
    water=pair['water_reference_fid']/40
    arrays={};rows=[]
    for label in ['full','methyl','methyl_water_carrier','methyl_on_carrier','aspartyl','methyl_repeat']:
        raw,first=trace['raw'][label],trace['first_dc'][label]
        np.testing.assert_allclose(first,first_dc(raw,trace['ppm'][label]),rtol=1e-11,atol=1e-12)
        config=trace['configs'][label]
        native,target=stages(raw,first,carrier=config['carrier_ppm'],scalar=scalar)
        for stage,fid in {**{k+'_native':v for k,v in native.items()},**target}.items():
            reference=native['raw'] if stage.endswith('_native') else target['raw_target']
            row={'response':label,'stage':stage,'carrier_ppm':config['carrier_ppm'],
                 **measurements(fid,reference,4000 if stage.endswith('_native') else 2000,
                                None if label in ['full','aspartyl'] else 3)}
            rows.append(row);arrays[label+'_'+stage]=fid
    basis=mrs_io.read_basis(str(REFERENCE/'common_response.BASIS'))
    formatted=basis.get_formatted_basis(2000,2048)[:,basis.names.index('NAA')]
    arrays['full_existing_formatted']=formatted
    serial_error=float(np.linalg.norm(arrays['full_formatted']-formatted)/np.linalg.norm(formatted))
    if serial_error>3e-7:raise ValueError('capture does not match existing formatted NAA')
    decomposition=decompose(arrays,water)
    arrays['ideal_3H_target']=decomposition.pop('ideal_3H_target')
    arrays['water_integration_target']=decomposition.pop('water_integration_target')
    arrays['raw_water_target']=water
    for key,protons in [('ideal_3H_target',3),('water_integration_target',2),('raw_water_target',2)]:
        row={'response':key,'stage':'analytical_control' if key!='raw_water_target' else 'raw_target',
             'carrier_ppm':None,**measurements(arrays[key],arrays[key],2000,protons)}
        row['spectral_integral']=spectral_integral(arrays[key],np.fft.fftshift(np.fft.fftfreq(2048,.0005))/127.731+4.68,
                                                 (1.8,2.2) if protons==3 else (1.65,7.65))
        rows.append(row)
    np.savez(out/'responses.npz',**arrays)
    with (out/'response_comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    methyl=arrays['methyl_raw_target'];on=arrays['methyl_on_carrier_raw_target'];off=arrays['methyl_water_carrier_raw_target']
    sensitivity={'methyl_RF3_vs_on_resonance_magnitude':float(abs(methyl[0]/on[0])),
        'methyl_RF4p7_vs_RF3_magnitude':float(abs(off[0]/methyl[0])),
        'methyl_RF4p7_vs_RF3_phase_rad':float(np.angle(off[0]/methyl[0])),
        'on_resonance_methyl_per_proton_vs_water':float(abs(on[0])/3/(abs(water[0])/2))}
    for label,carrier in [('methyl',3.),('methyl_water_carrier',4.7),('methyl_on_carrier',2.008)]:
        response=arrays[label+'_raw_target'];time=np.arange(len(response))*.0005
        frequency=(2.008-carrier)*127.728+(carrier-4.68)*127.731
        sensitivity[label+'_frequency_hz']=frequency
        expected=response[0]*np.exp(-np.pi*time)*np.exp(2j*np.pi*frequency*time)
        sensitivity[label+'_analytical_shape_relative_error']=float(np.linalg.norm(response-expected)/np.linalg.norm(expected))
    spin=loadmat(source/'FID-A/simulationTools/metabolites/spinSystems.mat',simplify_cells=True)
    to_json=lambda obj: obj.tolist() if isinstance(obj,np.ndarray) else obj.item() if isinstance(obj,np.generic) else str(obj)
    fits={};commands=[]
    for label,response,path in [('full_formatted',formatted,REFERENCE/'common_response.BASIS'),
        ('methyl_raw',methyl,out/'methyl_raw.BASIS')]:
        if label=='methyl_raw':
            # Installed LCModel reader squeezes a one-column BASIS incorrectly.
            ref=basis.get_formatted_basis(2000,2048)[:,basis.names.index('Ref0ppm')]
            write_basis(path,{'NAA':response,'Ref0ppm':ref},{'Bfield':3.,'sw':2000.,'points':2048})
        fits[label],command=fit_control(out/label,response,pair,path);commands.append(command)
    (out/'cli_commands.txt').write_text('\n'.join(commands)+'\n')
    after={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in tracked}
    if hashes!=after:raise RuntimeError('source or existing reference changed')
    result={'FSL_MRS_version':fsl_mrs.__version__,'normalization_scalar':scalar,
        'source_hashes':hashes,'source_hashes_unchanged':True,'sequence_configuration':manifest['sequence_configuration'],
        'spin_systems':{key:spin[key] for key in ['sysNAA','sysNAA_Acetyl_only','sysNAA_Aspartate_only','sysH2O']},
        'response_comparison':rows,'decomposition':decomposition,'carrier_sensitivity':sensitivity,
        'ideal_grid_controls':ideal_grid_controls(),
        'raw_superposition_max_error':float(np.max(abs(trace['raw']['full']-trace['raw']['methyl']-trace['raw']['aspartyl']))),
        'methyl_repeat_max_error':float(np.max(abs(trace['raw']['methyl']-trace['raw']['methyl_repeat']))),
        'formatted_vs_existing_relative_error':serial_error,
        'phase5h_alpha':.8984749391468103,'fits':fits}
    (out/'results.json').write_text(json.dumps(result,indent=2,default=to_json)+'\n')
    plot(out,arrays)
    print(json.dumps({'decomposition':decomposition,'carrier_sensitivity':sensitivity,
                      'superposition_max_error':result['raw_superposition_max_error']},indent=2))
    return result


def plot(out,arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t=np.arange(2048)*.0005;ppm=np.fft.fftshift(np.fft.fftfreq(2048,.0005))/127.731+4.68
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for key in ['ideal_3H_target','methyl_raw_target','full_raw_target','full_first_dc_target','full_formatted']:
        axes[0].plot(ppm,FIDToSpec(arrays[key]).real,label=key)
        axes[1].plot(t[:12]*1000,abs(arrays[key][:12]),'.-',label=key)
    axes[0].set(xlim=(2.2,1.8),xlabel='ppm',ylabel='Real spectrum per mol/L')
    axes[1].set(xlabel='Time (ms)',ylabel='FID magnitude per mol/L')
    for ax in axes:ax.legend(fontsize=6)
    fig.tight_layout();fig.savefig(out/'response_processing.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for key in ['methyl_raw_target','methyl_water_carrier_raw_target','methyl_on_carrier_raw_target']:
        axes[0].plot(t[:30]*1000,arrays[key][:30].real,label=key)
    axes[0].set(xlabel='Time (ms)',ylabel='Real FID per mol/L');axes[0].legend(fontsize=6)
    factors=decompose(arrays,arrays['raw_water_target'])['factors']
    axes[1].barh(list(factors),list(factors.values()));axes[1].axvline(1,color='k',ls='--')
    axes[1].set_xlabel('Ordered contribution to alpha')
    fig.tight_layout();fig.savefig(out/'carrier_and_factors.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir');parser.add_argument('--reuse-capture',action='store_true')
    args=parser.parse_args();run(args.output_dir,args.reuse_capture)

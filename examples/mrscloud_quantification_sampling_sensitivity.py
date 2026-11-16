"""Fixed-response PRESS sampling/linewidth sweep; no production calibration changes."""

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
from fsl_mrs.utils import fitting, mrs_io
from fsl_mrs.utils.misc import ts_to_ts_ft

from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral
from examples.fsl_mrs_sequence_water_quantification import quantification_info
from examples.fsl_mrs_water_fit_diagnostics import fit_water
from examples.fsl_mrs_water_preprocessing import cli_command
from examples.mrscloud_naa_response_decomposition import first_dc, REFERENCE
from examples.mrscloud_water_processing_validation import second_dc
from examples.mrscloud_sequence_matched_water import write_basis
from simulation.mrscloud_quantitative_pair import export_pair

POINTS = [2048, 4096, 8192]
LINEWIDTHS = [1, 3, 5, 8]


def broaden(response, points, linewidth_hz):
    """Total linewidth = original 1 Hz + nonnegative additional Lorentzian width."""
    if points not in POINTS or linewidth_hz not in LINEWIDTHS or len(response) < points:
        raise ValueError('unsupported or uncovered sampling/linewidth condition')
    return response[:points]*np.exp(-np.pi*(linewidth_hz-1)*np.arange(points)/2000)


def load_master(capture_dir):
    directory=Path(capture_dir)
    naa=loadmat(directory/'naa/naa_trace.mat',simplify_cells=True)
    water=loadmat(directory/'water/processing_trace.mat',simplify_cells=True)
    manifest=json.loads((REFERENCE/'reference.json').read_text())
    scalar=manifest['normalization_scalar']
    if len(naa['raw']['full'])!=16384 or len(water['raw_fid'])!=16384:
        raise ValueError('requires genuine 16384-point native readouts, not zero filling')
    if naa['configs']['full']['carrier_ppm']!=3 or water['sw']!=4000 or water['spatial_points']!=101:
        raise ValueError('incompatible sequence capture')
    old=loadmat('outputs/mrscloud_naa_response_decomposition/validated/octave/naa_trace.mat',simplify_cells=True)
    old_water=np.load(REFERENCE/'raw_water.npz')['raw_water']
    np.testing.assert_array_equal(naa['raw']['full'][:8192],old['raw']['full'])
    np.testing.assert_array_equal(water['raw_fid'][:8192],old_water)
    t=np.arange(16384)/4000
    mix_m=np.exp(2j*np.pi*(3.-4.68)*127.731*t)
    mix_w=np.exp(2j*np.pi*(4.7-4.68)*127.731*t)
    responses={}
    for name in ['full','methyl']:
        raw=naa['raw'][name]
        first=naa['first_dc'][name]
        np.testing.assert_allclose(first,first_dc(raw,naa['ppm'][name]),rtol=1e-11,atol=1e-12)
        second,_=second_dc(first,4000,127.731,3.)
        for label,values in [('raw',raw),('DC',second)]:
            native=values/scalar*mix_m
            responses[name+'_'+label]=native[::2].copy()
            if label=='DC':
                responses[name+'_formatted']=ts_to_ts_ft(native[:,None],.00025,.0005,8192).ravel()
    responses['water_raw']=water['raw_fid']*mix_w/scalar
    responses['water_raw']=responses['water_raw'][::2].copy()
    responses['ideal_methyl']=1.5*abs(responses['water_raw'][0])*np.exp(-np.pi*np.arange(8192)/2000)*np.exp(
        2j*np.pi*((2.008-3)*127.728+(3-4.68)*127.731)*np.arange(8192)/2000)
    return responses,manifest


def predict(master, points, linewidth_hz):
    time=np.arange(points)/2000
    ppm=np.fft.fftshift(np.fft.fftfreq(points,.0005))/127.731+4.68
    jw=spectral_integral(abs(master['water_raw'][0])*np.exp(-np.pi*linewidth_hz*time),ppm,(1.65,7.65))
    integrals={name:spectral_integral(broaden(fid,points,linewidth_hz),ppm,(1.8,2.2))
               for name,fid in master.items() if name!='water_raw'}
    factors={name:jm/jw*(2/3) for name,jm in integrals.items()}
    return {'unit_metabolite_integral':integrals['full_formatted'],'unit_water_integral':jw,
            'response_factor':factors['full_formatted'],'predicted_mM':10*factors['full_formatted'],
            'control_factors':factors, 'control_integrals':integrals}


def fit_condition(directory, master, manifest, points, linewidth_hz):
    prediction=predict(master,points,linewidth_hz)
    metab=broaden(master['full_formatted'],points,linewidth_hz)
    water=broaden(master['water_raw'],points,linewidth_hz)
    pair={'metabolite_fid':.010*metab,'water_reference_fid':40*water,
        'ground_truth':{'metabolite_concentrations_mM':{'NAA':10.},'water_concentration_mol_l':40.,
                        'physical_calibration':False},
        'simulation_config':{'dwell_time_s':.0005,'spectrometer_frequency_mhz':127.731,'frequency_origin_ppm':4.68,
            'points':points,'bandwidth_hz':2000,'receiver_factor':1.,'physical_calibration':False,
            'total_linewidth_hz':linewidth_hz,'additional_linewidth_hz':linewidth_hz-1,
            'common_normalization_scalar':manifest['normalization_scalar'],
            'amplitude_reference_manifest_sha256':hashlib.sha256((REFERENCE/'reference.json').read_bytes()).hexdigest(),
            'processing':'fixed longest-record DC/format master; prefixes; equal additional decay',
            'sequence_configuration':manifest['sequence_configuration'],
            'native_readout_points':16384,'native_bandwidth_hz':4000,
            'original_production_grid_compatibility_claimed':False}}
    export_pair(pair,directory)
    basis_path=directory/'matching.BASIS'
    write_basis(basis_path,{'NAA':master['full_formatted'][:points],
        'MethylControl':master['methyl_formatted'][:points]}, {'Bfield':3.,'sw':2000.,'points':points})
    mrs=mrs_io.read_FID(str(directory/'metabolites.nii.gz')).mrs(basis=mrs_io.read_basis(str(basis_path)),
        ref_data=str(directory/'water_reference.nii.gz'))
    mrs.keep=['NAA'];mrs.rescaleForFitting()
    fit=fitting.fit_FSLModel(mrs,method='Newton',model='lorentzian',baseline_order=-1,
        ppmlim=(.2,4.2),metab_groups=[0],capture_minimize_output=True)
    water_info,_,_=fit_water(mrs)
    q=quantification_info(127.731,match_water=True,match_relaxation=True)
    fit.calculateConcScaling(mrs,quant_info=q,internal_reference=['NAA'])
    reported=float(fit.getConc(scaling='molarity',metab='NAA'))
    row={'points':points,'record_duration_s':points/2000,'bin_width_hz':2000/points,
        'total_linewidth_hz':linewidth_hz,'additional_linewidth_hz':linewidth_hz-1,
        **{k:v for k,v in prediction.items() if k not in ['control_factors','control_integrals']},
        'reported_mM':reported,'difference_from_prediction_mM':reported-prediction['predicted_mM'],
        'bias_from_truth_mM':reported-10,
        'raw_NAA_coefficient_mol_l':float(fit.getConc(scaling='raw',metab='NAA')),
        'internal_NAA_coefficient':float(fit.getConc(metab='NAA')),
        'fitted_water_amplitude':water_info['fitted_amplitude'],
        'expected_water_amplitude':float(abs(pair['water_reference_fid'][0])),
        'metabolite_relative_residual':float(np.linalg.norm(mrs.FID-fit.predictedFID(mrs))/np.linalg.norm(mrs.FID)),
        'water_relative_residual':water_info['relative_complex_residual'],
        'metabolite_success':bool(fit.scipy_minimize_result.success),'water_success':water_info['success'],
        'initial_metabolite_real':float(pair['metabolite_fid'][0].real),
        'initial_metabolite_imag':float(pair['metabolite_fid'][0].imag),
        'max_export_error':max(float(np.max(abs(mrs_io.read_FID(str(directory/f)).mrs().FID-pair[k]))) for k,f in
             [('metabolite_fid','metabolites.nii.gz'),('water_reference_fid','water_reference.nii.gz')])}
    details={'row':row,'prediction':prediction,'water_fit':water_info,
        'fit_parameters':dict(zip(fit.fitResults.columns,fit.fitResults.iloc[0].tolist())),
        'scaling':mrs.scaling,'quantification':{'TE_s':.030,'TR_s':2.,'effective_water_M':40.,
            'water_relaxation_factor':q.R_H2O_WM,'metabolite_relaxation_factor':q.R_M,'additional_correction':1.}}
    (directory/'fitting_results.json').write_text(json.dumps(details,indent=2)+'\n')
    return details


def cli_report(directory):
    command=cli_command(directory)
    command[command.index('--basis')+1]=str(directory/'matching.BASIS')
    # Public CLI lacks set_densities: this is the declared 40M/WM-water ratio,
    # independent of any fitted result. It reproduces the API density override.
    density_ratio=40000/(55010*.65)
    t1=directory/'synthetic_t1.json';t2=directory/'synthetic_t2.json'
    keys=['H2O_WM','H2O_GM','H2O_CSF','METAB']
    t1.write_text(json.dumps(dict.fromkeys(keys,1e-6)))
    t2.write_text(json.dumps(dict.fromkeys(keys,1e20)))
    command+=['--t1-values',str(t1),'--t2-values',str(t2),'--h2o_scale',str(density_ratio)]
    completed=subprocess.run(command,capture_output=True,text=True)
    (directory/'cli.log').write_text(completed.stdout+completed.stderr)
    completed.check_returncode()
    if (directory/'fsl_mrs_results/report.html').stat().st_size<1000:
        raise RuntimeError('CLI HTML report missing')
    with (directory/'fsl_mrs_results/concentrations.csv').open() as stream:
        records=list(csv.reader(stream))
    cli=dict(zip(records[1][1:5],map(float,records[2][1:5])))
    return {'command':shlex.join(command),'molarity_mM':cli['molarity'],
            'density_configuration_ratio':density_ratio,'T1_s':1e-6,'T2_s':1e20}


def run(output_dir, capture_dir='outputs/mrscloud_sampling_sensitivity/octave'):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=False)
    paths=list(REFERENCE.iterdir())+[Path(capture_dir)/'naa/naa_trace.mat',Path(capture_dir)/'water/processing_trace.mat']
    paths+=list(Path('simulation/basissets/MRSCloud-main/functions').rglob('sim_signals.m'))
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    master,manifest=load_master(capture_dir)
    np.savez(out/'master_responses.npz',**master)
    old=np.load('outputs/mrscloud_naa_response_decomposition/validated/responses.npz')
    baseline={'relative_master_vs_phase5i_metabolite':float(np.linalg.norm(master['full_formatted'][:2048]-old['full_existing_formatted'])/np.linalg.norm(old['full_existing_formatted'])),
              'initial_master_vs_phase5i_metabolite':str(master['full_formatted'][0]-old['full_existing_formatted'][0])}
    cases=[];commands=[]
    for points in POINTS:
        for linewidth in LINEWIDTHS:
            directory=out/f'N{points}_LW{linewidth}'
            case=fit_condition(directory,master,manifest,points,linewidth)
            if (points,linewidth) in [(2048,1),(8192,1),(2048,8)]:
                cli=cli_report(directory)
                if not np.isclose(cli['molarity_mM'],case['row']['reported_mM'],rtol=1e-12,atol=1e-12):
                    raise ValueError('CLI matched-assumption quantification differs from API')
                case['CLI']=cli;commands.append(cli['command'])
            cases.append(case)
    after={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    if hashes!=after:raise RuntimeError('source calibration or capture changed')
    result={'FSL_MRS_version':fsl_mrs.__version__,'normalization_scalar':manifest['normalization_scalar'],
        'source_hashes':hashes,'source_hashes_unchanged':True,'configuration':{
            'points':POINTS,'total_linewidths_hz':LINEWIDTHS,'bandwidth_hz':2000,
            'NAA_mM':10,'water_M':40,'receiver_factor':1,'sequence':manifest['sequence_configuration'],
            'native_points':16384,'master_target_points':8192,
            'processing_fixed_before_sweep':True},'baseline_comparison':baseline,'cases':cases}
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'cli_commands.txt').write_text('\n'.join(commands)+'\n')
    with (out/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=cases[0]['row'].keys());writer.writeheader()
        writer.writerows(case['row'] for case in cases)
    plot(out,cases)
    print(json.dumps([case['row'] for case in cases],indent=2))
    return result


def plot(out,cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for points in POINTS:
        selected=[c for c in cases if c['row']['points']==points]
        axes[0].plot(LINEWIDTHS,[c['row']['reported_mM']-10 for c in selected],'o-',label=str(points)+' points')
    axes[0].axhline(0,color='k',ls='--');axes[0].set(xlabel='Total linewidth (Hz)',ylabel='Concentration bias (mM)');axes[0].legend()
    for name in ['ideal_methyl','methyl_raw','full_raw','full_DC','full_formatted']:
        selected=[c for c in cases if c['row']['total_linewidth_hz']==1]
        axes[1].plot([p/2000 for p in POINTS],[c['prediction']['control_factors'][name] for c in selected],'o-',label=name)
    axes[1].axhline(1,color='k',ls='--');axes[1].set(xlabel='Record duration (s)',ylabel='Unfitted response factor (1 Hz)');axes[1].legend(fontsize=7)
    fig.tight_layout();fig.savefig(out/'sampling_linewidth_bias.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir');parser.add_argument('--capture-dir',default='outputs/mrscloud_sampling_sensitivity/octave')
    args=parser.parse_args();run(args.output_dir,args.capture_dir)

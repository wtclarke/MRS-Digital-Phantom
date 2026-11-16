"""Restartable Phase 7a Voigt Monte Carlo validation on one cached WM acquisition."""
import argparse
from copy import deepcopy
import csv
import json
from pathlib import Path
import shlex
import subprocess
import warnings

import fsl_mrs
import numpy as np
from fsl_mrs.utils import fitting,mrs_io

from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral
from examples.fsl_mrs_water_fit_diagnostics import fit_water
from examples.mrscloud_mixed_tissue_fitting_validation import prediction,quant_info
from examples.mrscloud_quantitative_voxel import BUNDLE,METABOLITES,WATER
from simulation.mrscloud_full_basis import load_library,sha256
from simulation.mrscloud_quantitative_noise import add_pair_noise
from simulation.mrscloud_quantitative_pair import export_pair
from simulation.mrscloud_quantitative_voxel import simulate_voxel

NAMES=['NAA','Cr','GPC','PCh','Glu']
LEVELS=[100,50,25,15,10,5]
SOURCE=Path('outputs/mrscloud_quantitative_voxel/final/pure_WM')


def finite(value):
    value=float(value)
    return value if np.isfinite(value) else None


def save_json(path,value):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def fit_one(directory,basis,responses,expected,baseline=None):
    gt=json.loads((directory/'ground_truth.json').read_text());config=json.loads((directory/'simulation_config.json').read_text())
    mrs=mrs_io.read_FID(str(directory/'metabolites.nii.gz')).mrs(basis=deepcopy(basis),ref_data=str(directory/'water_reference.nii.gz'))
    mrs.keep=NAMES
    clean=mrs_io.read_FID(str(directory/('metabolites_clean.nii.gz' if config['noise'] else 'metabolites.nii.gz'))).mrs().FID
    mrs.rescaleForFitting()
    with warnings.catch_warnings(record=True) as caught:
        fit=fitting.fit_FSLModel(mrs,method='Newton',model='voigt',baseline_order=-1,
            ppmlim=(.2,4.2),metab_groups=[0]*5,capture_minimize_output=True)
        water,_,water_model=fit_water(mrs)
    sd,sb=mrs.scaling['FID'],mrs.scaling['basis'];gain=config['receiver_factor'];conversion=1000*sb/sd/gain
    params=dict(zip(fit.fitResults.columns,fit.fitResults.iloc[0].tolist()))
    output=dict(optimizer_success=bool(fit.scipy_minimize_result.success),optimizer_message=str(fit.scipy_minimize_result.message),
        parameters=params,scaling=mrs.scaling,warnings=[str(w.message) for w in caught],water_fit=water,
        relative_complex_residual=float(np.linalg.norm(mrs.FID-fit.predictedFID(mrs))/np.linalg.norm(mrs.FID)),
        relative_complex_error_to_clean=float(np.linalg.norm(clean-fit.predictedFID(mrs)/sd)/np.linalg.norm(clean)),
        FSL_reported_SNR={k:finite(v) for k,v in fit.SNR.peaks.mean().items()},
        FSL_residual_SNR=finite(fit.SNR.residual),coefficients={})
    output['coefficient_covariance_names']=list(fit.metabs)
    output['coefficient_covariance_mM2']=(fit.cov[:5,:5]*conversion**2).tolist()
    clean_water=mrs_io.read_FID(str(directory/('water_reference_clean.nii.gz' if config['noise'] else 'water_reference.nii.gz'))).mrs().FID
    water.update(expected_initial_amplitude=float(abs(clean_water[0])),
                 relative_amplitude_error=float(water['fitted_amplitude']/abs(clean_water[0])-1))
    q=quant_info(mrs,gt,False)
    original=fit.params.copy()
    try:
        fit.calculateConcScaling(mrs,quant_info=q,internal_reference=['NAA'])
        output['quantification_success']=True
    except (ValueError,ZeroDivisionError) as error:
        output.update(quantification_success=False,quantification_error=f'{type(error).__name__}: {error}')
    np.testing.assert_array_equal(fit.params,original)
    t=np.arange(2048)/2000;ppm=np.fft.fftshift(np.fft.fftfreq(2048,.0005))/127.731+4.68
    naa=float(fit.getConc(metab='NAA')*sb/sd)
    naa_model=naa*responses['NAA']*np.exp(-(params['gamma_0']+1j*params['eps_0'])*t-params['sigma_0']**2*t*t)
    _,gamma,sigma,_,_=water['parameters']
    water_corrected=water['fitted_amplitude']*np.exp(-gamma*t-sigma*t*t)
    jm=spectral_integral(naa_model,ppm,(1.8,2.2));jw=spectral_integral(water_corrected,ppm,(1.65,7.65))
    factor=q.relax_corr_water_molar*q.csf_corr_molar*q.relax_corr_metab*2/3
    for names in [[n] for n in NAMES]+[['GPC','PCh']]:
        name='+'.join(names);indices=[fit.metabs.index(n) for n in names]
        estimate=float(sum(fit.getConc(metab=n) for n in names)*conversion)
        truth=sum(gt['voxel_average_metabolite_concentrations_mM'][n] for n in names)
        variance=float(np.sum(fit.cov[np.ix_(indices,indices)]))*conversion**2
        uncertainty=finite(np.sqrt(variance)) if variance>=0 else None
        reported=float(sum(fit.getConc(scaling='molarity',metab=n) for n in names)) if output['quantification_success'] else None
        reconstruction=(sum(float(fit.getConc(metab=n))*sb/sd for n in names)/naa*jm/jw*factor) if naa>0 and jw>0 else None
        reference=sum(expected['predicted_mM'][n] for n in names)
        baseline_value=baseline['coefficients'][name]['standard_FSL_mM'] if baseline else reported
        output['coefficients'][name]=dict(truth_mM=truth,estimate_mM=estimate,absolute_error_mM=estimate-truth,
            relative_error=estimate/truth-1,reported_SD_mM=uncertainty,
            covered_95=bool(abs(estimate-truth)<=1.96*uncertainty) if config['noise'] and uncertainty is not None and uncertainty>0 else None,
            fitted_internal=float(sum(fit.getConc(metab=n) for n in names)),expected_internal=truth/1000*gain*sd/sb,
            standard_FSL_mM=reported,independent_clean_FSL_prediction_mM=reference,
            reconstructed_FSL_mM=finite(reconstruction) if reconstruction is not None else None,
            difference_from_nominal_mM=reported-truth if reported is not None else None,
            difference_from_clean_prediction_mM=reported-reference if reported is not None else None,
            difference_from_noise_free_fit_mM=reported-baseline_value if reported is not None and baseline_value is not None else None)
    noise=config.get('noise_parameters');output['noise_parameters']=noise
    if noise:
        # Independent estimate from the actual acquisition difference, not utility metadata.
        spectral_noise=np.fft.fft(mrs.FID/sd-clean,norm='ortho')
        estimate=float(np.sqrt((spectral_noise.real.var(ddof=1)+spectral_noise.imag.var(ddof=1))/2))
        clean_naa=responses['NAA']*gt['voxel_average_metabolite_concentrations_mM']['NAA']/1000*gain
        mask=(ppm>=1.8)&(ppm<=2.2)
        peak=float(abs(np.fft.fftshift(np.fft.fft(clean_naa,norm='ortho'))[mask]).max())
        output['independently_measured_SNR']=peak/estimate
    np.savez(directory/'fitted_arrays.npz',data=mrs.FID/sd,clean=clean,fitted=fit.predictedFID(mrs)/sd,
             water=mrs.H2O/sd,water_fitted=water_model/sd)
    output['input_hashes']={p.name:sha256(p) for p in directory.glob('*.nii.gz')}
    return output


def summarize(records):
    """All returned numerical fits, including non-success exits; no favorable-fit filtering."""
    rows=[]
    for level in dict.fromkeys(record['snr'] for record in records):
        selected=[r for r in records if r['snr']==level]
        for name in [*NAMES,'GPC+PCh']:
            values=[r['fit']['coefficients'][name] for r in selected if 'coefficients' in r['fit']]
            errors=np.array([v['absolute_error_mM'] for v in values]);estimates=np.array([v['estimate_mM'] for v in values])
            intervals=[v['covered_95'] for v in values if v['covered_95'] is not None]
            concentrations=[v['standard_FSL_mM'] for v in values if v['standard_FSL_mM'] is not None]
            rows.append(dict(snr=level,metabolite=name,realizations=len(selected),returned_fits=len(values),
                mean_bias_mM=finite(errors.mean()) if len(errors) else None,
                SD_mM=finite(estimates.std(ddof=1)) if len(estimates)>1 else None,
                RMSE_mM=finite(np.sqrt(np.mean(errors**2))) if len(errors) else None,
                median_absolute_error_mM=finite(np.median(abs(errors))) if len(errors) else None,
                coverage_95=float(np.mean(intervals)) if intervals else None,coverage_valid_count=len(intervals),
                failure_rate=float(np.mean([not r['fit'].get('optimizer_success',False) for r in selected])),
                quantification_failure_rate=float(np.mean([not r['fit'].get('quantification_success',False) for r in selected])),
                poor_recovery_fraction=float(np.mean([abs(v['relative_error'])>.2 for v in values])) if values else None,
                mean_standard_FSL_mM=finite(np.mean(concentrations)) if concentrations else None,
                SD_standard_FSL_mM=finite(np.std(concentrations,ddof=1)) if len(concentrations)>1 else None))
    return rows


def write_csv(path,rows):
    with path.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)


def run_validation(output_dir,*,smoke=False):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    levels=[50,10] if smoke else LEVELS;realizations=2 if smoke else 20
    tracked=[p for p in BUNDLE.iterdir() if p.is_file()]+list(SOURCE.glob('*.nii.gz'))+list(SOURCE.glob('*.json'))
    source_hashes={str(p):sha256(p) for p in tracked}
    config=dict(smoke=smoke,levels=levels,realizations=realizations,water_SNR=1000,seed_recipe='70000+100*level_index+realization; spawn 0 metabolite,1 water',
        FSL_MRS_version=fsl_mrs.__version__,numpy_version=np.__version__,source_hashes=source_hashes,
        implementation_hashes={str(p):sha256(p) for p in [Path(__file__),Path('simulation/mrscloud_quantitative_noise.py')]},
        fitting=dict(method='Newton',model='voigt',baseline='off',ppm=[.2,4.2],shared_group=True,TE_s=.03,TR_s=2.))
    path=out/'validation_config.json'
    if path.exists() and json.loads(path.read_text())!=config:raise ValueError('incompatible cached validation configuration; select a new output directory')
    save_json(path,config)
    clean=simulate_voxel({'WM':1},calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,water_concentrations_mol_l=WATER)
    responses,water,calibration=load_library(BUNDLE,selected_names=NAMES)
    basis=mrs_io.read_basis(str(BUNDLE/'metabolite_basis.BASIS'))
    expected=prediction(clean['ground_truth'],responses,water)
    for key,filename in [('metabolite_fid','metabolites'),('water_reference_fid','water_reference')]:
        np.testing.assert_array_equal(clean[key],mrs_io.read_FID(str(SOURCE/f'{filename}.nii.gz')).mrs().FID)
    records=[];cache_hits=0;baseline=None
    for level,realization in [('infinite',0)]+[(level,r) for level in levels for r in range(realizations)]:
        directory=out/f'SNR_{level}'/f'realization_{realization:02d}'
        seed=None if level=='infinite' else 70000+100*levels.index(level)+realization
        pair=clean if seed is None else add_pair_noise(clean,seed=seed,metabolite_snr=level,water_snr=1000)
        if not directory.exists():export_pair(pair,directory)
        else:
            for key,filename in [('metabolite_fid','metabolites'),('water_reference_fid','water_reference')]:
                np.testing.assert_array_equal(mrs_io.read_FID(str(directory/f'{filename}.nii.gz')).mrs().FID,pair[key])
            if json.loads((directory/'simulation_config.json').read_text())!=pair['simulation_config']:
                raise ValueError('cached simulation config mismatch')
        path=directory/'fit_result.json'
        if path.exists():
            fitted=json.loads(path.read_text());cache_hits+=1
            if fitted['input_hashes']!={p.name:sha256(p) for p in directory.glob('*.nii.gz')}:
                raise ValueError('cached fit inputs changed')
        else:
            try:fitted=fit_one(directory,basis,responses,expected,baseline)
            except Exception as error:
                fitted=dict(optimizer_success=False,quantification_success=False,error=f'{type(error).__name__}: {error}',
                    input_hashes={p.name:sha256(p) for p in directory.glob('*.nii.gz')})
            save_json(path,fitted)
        if level=='infinite':
            baseline=fitted
            if not fitted.get('optimizer_success') or max(abs(c['relative_error']) for c in fitted['coefficients'].values())>5e-6:
                raise RuntimeError('noise-free Voigt baseline failed; do not run SNR sweep')
        records.append(dict(snr=str(level),realization=realization,seed=seed,directory=str(directory),fit=fitted))
        print(f'SNR {level}, realization {realization}: optimizer={fitted.get("optimizer_success")}',flush=True)
    summary=summarize(records);write_csv(out/'summary.csv',summary)
    rows=[dict(snr=r['snr'],realization=r['realization'],seed=r['seed'],metabolite=n,**v)
          for r in records for n,v in r['fit'].get('coefficients',{}).items()]
    write_csv(out/'coefficients.csv',rows)
    commands=[];reports=[]
    for level in ['infinite',*([50,10] if smoke else [50,25,10,5])]:
        record=next(r for r in records if r['snr']==str(level) and r['fit'].get('optimizer_success') and r['fit'].get('quantification_success'))
        directory=Path(record['directory']);result_dir=directory/'fsl_mrs_results';report=result_dir/'report.html'
        command=['pixi','run','fsl_mrs','--data',str(directory/'metabolites.nii.gz'),'--h2o',str(directory/'water_reference.nii.gz'),
            '--basis',str(BUNDLE/'metabolite_basis.BASIS'),'--output',str(result_dir),'--keep',*NAMES,'--algo','Newton',
            '--baseline','off','--ppmlim','.2','4.2','--no_conj_basis','--TE','30','--TR','2','--tissue_frac','1','0','0',
            '--internal_ref','NAA','--wref_metabolite','NAA','--ref_protons','3','--ref_int_limits','1.8','2.2','--report']
        commands.append(shlex.join(command))
        if not report.exists():
            completed=subprocess.run(command,capture_output=True,text=True)
            (directory/'cli.log').write_text(completed.stdout+completed.stderr);completed.check_returncode()
        if report.stat().st_size<1000 or 'Water referencing' not in report.read_text():raise RuntimeError('CLI report invalid')
        with (result_dir/'concentrations.csv').open() as stream:data=list(csv.reader(stream))
        comparisons={row[0]:float(row[4])-record['fit']['coefficients'][row[0]]['standard_FSL_mM'] for row in data[2:] if row[0] in NAMES}
        reports.append(dict(snr=str(level),realization=record['realization'],path=str(report),CLI_API_concentration_difference_mM=comparisons))
    (out/'cli_commands.txt').write_text('\n'.join(commands)+'\n')
    if source_hashes!={str(p):sha256(p) for p in tracked}:raise RuntimeError('source calibration or clean acquisition changed')
    result=dict(config=config,records=records,summary=summary,reports=reports,cache_hits=cache_hits,
        source_unchanged=True,calibration_identity=calibration['calibration_identity'],independent_clean_prediction=expected,
        failure_definition='exception or optimizer success false; all returned fits included in coefficient statistics',
        poor_recovery_definition='absolute relative coefficient error >20%; descriptive, never an exclusion',
        coverage_definition='approximate coefficient Wald +/-1.96 reported SD; GPC+PCh uses covariance including cross term')
    save_json(out/'results.json',result);plots(out,result)
    return result


def plots(out,result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for name in [*NAMES,'GPC+PCh']:
        rows=[r for r in result['summary'] if r['metabolite']==name and r['snr']!='infinite']
        x=[int(r['snr']) for r in rows]
        for ax,key in zip(axes.flat,['mean_bias_mM','SD_mM','RMSE_mM','coverage_95','mean_standard_FSL_mM','poor_recovery_fraction']):
            ax.plot(x,[r[key] for r in rows],'-o',label=name);ax.set_xscale('log');ax.set_xlabel('Simulator SNR');ax.set_title(key);ax.grid(alpha=.3)
    for ax in axes.flat:ax.invert_xaxis()
    axes[0,0].legend(fontsize=8);fig.tight_layout();fig.savefig(out/'noise_bias_precision.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    levels=list(dict.fromkeys(r['snr'] for r in result['records'] if r['snr']!='infinite'))
    for level in levels:
        records=[r for r in result['records'] if r['snr']==level and 'water_fit' in r['fit']]
        axes[0].scatter([int(level)]*len(records),[100*r['fit']['water_fit']['relative_amplitude_error'] for r in records],alpha=.5)
        axes[1].scatter([int(level)]*len(records),[r['fit']['relative_complex_error_to_clean'] for r in records],alpha=.5)
        axes[2].scatter([int(level)]*len(records),[r['fit']['independently_measured_SNR'] for r in records],alpha=.5)
    for ax,title in zip(axes,['Water amplitude error %','Fitted complex error to clean FID','Independently measured simulator SNR']):
        ax.set_title(title);ax.set_xlabel('Nominal simulator SNR');ax.set_xscale('log');ax.invert_xaxis();ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(out/'noise_water_and_signal_errors.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output_dir');parser.add_argument('--smoke',action='store_true')
    args=parser.parse_args();result=run_validation(args.output_dir,smoke=args.smoke)
    print(f'Completed {len(result["records"])} fits; {len(result["reports"])} HTML reports; {result["cache_hits"]} cached fits')

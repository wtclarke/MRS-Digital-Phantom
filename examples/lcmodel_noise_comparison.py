"""Frozen-protocol LCModel fits of existing Phase 7a acquisitions; no simulation."""
import argparse
import csv
import hashlib
import inspect
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime,timezone

import numpy as np
from fsl_mrs.utils import mrs_io

from examples.lcmodel_initial_validation import BINARY,digest,array_errors
from utils.lcmodel_quantitative_io import NAMES,parse_raw,write_raw
from utils.lcmodel_noise_analysis import METABOLITES,parse_fit,descriptive,paired_bootstrap

SOURCE=Path('outputs/mrscloud_noise_fitting_validation/full')
BASE=Path('outputs/lcmodel_validation/noise_free')
OUTPUT=Path('outputs/lcmodel_noise_comparison')
LEVELS=[100,50,25,15,10,5]
PATH_FIELDS=['FILRAW','FILH2O','FILBAS','FILPS','FILPRI','FILTAB','FILCOO']


def atomic_json(path,data):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');temp.replace(path)


def text_hash(text):return hashlib.sha256(text.encode()).hexdigest()


def canonical(control):
    for field in PATH_FIELDS:control=re.sub(rf"({field}=)'[^']+'",rf"\1'<{field}>'",control)
    return control


def control_for_case(protocol,directory,attempt):
    mapping=dict(FILRAW=directory/'metabolites.RAW',FILH2O=directory/'water_reference.H2O',
        FILBAS=(BASE/'lcmodel_compatible.BASIS').resolve(),
        **{f:attempt/'lcmodel_outputs'/f'fit.{s}' for f,s in [('FILPS','PS'),('FILPRI','PRINT'),('FILTAB','TABLE'),('FILCOO','COORD')]})
    result=protocol
    for field,path in mapping.items():
        path=str(path.resolve())
        if "'" in path:raise ValueError('single quotes unsupported in LCModel paths')
        result,count=re.subn(rf"({field}=)'[^']+'",lambda m:m[1]+f"'{path}'",result)
        if count!=1:raise ValueError('missing or duplicate path in frozen protocol')
    if canonical(result)!=canonical(protocol):raise ValueError('non-path CONTROL settings changed')
    return result


def inventory(smoke=False):
    root=json.loads((SOURCE/'results.json').read_text())
    expected=[('infinite',0)]+[(str(s),r) for s in LEVELS for r in range(20)]
    indexed={(r['snr'],r['realization']):r for r in root['records']}
    if len(indexed)!=121 or set(indexed)!=set(expected):raise ValueError('Phase 7a acquisition inventory incomplete or duplicated')
    actual={p.relative_to(SOURCE).as_posix() for p in SOURCE.glob('SNR_*/realization_*')}
    if actual!={f'SNR_{s}/realization_{r:02d}' for s,r in expected}:raise ValueError('unexpected Phase 7a acquisition directories')
    records=[];hashes={str((SOURCE/'results.json').resolve()):digest(SOURCE/'results.json')}
    truth=json.loads((SOURCE/'SNR_infinite/realization_00/ground_truth.json').read_text())
    for level,r in expected:
        record=indexed[level,r];directory=SOURCE/f'SNR_{level}'/f'realization_{r:02d}'
        config=json.loads((directory/'simulation_config.json').read_text());gt=json.loads((directory/'ground_truth.json').read_text())
        fit=json.loads((directory/'fit_result.json').read_text())
        if fit!=record['fit'] or gt!=truth:raise ValueError(f'Phase 7a result/ground truth mismatch: {directory}')
        if config['synthetic_calibration']['calibration_identity']!=root['calibration_identity']:raise ValueError('calibration mismatch')
        if (config['points'],config['dwell_time_s'],config['spectrometer_frequency_mhz'],config['frequency_origin_ppm'])!=(2048,.0005,127.731,4.68):raise ValueError('acquisition grid mismatch')
        if config['receiver_factor']!=1 or gt['effective_voxel_water_concentration_mol_l']!=39.0165:raise ValueError('unexpected gain or water concentration')
        for filename,expected_hash in fit['input_hashes'].items():
            if digest(directory/filename)!=expected_hash:raise ValueError(f'input hash mismatch: {directory/filename}')
        for path in directory.rglob('*'):
            if path.is_file():hashes[str(path.resolve())]=digest(path)
        noise=config.get('noise_parameters');seed=None if level=='infinite' else 70000+100*LEVELS.index(int(level))+r
        if record['seed']!=seed or (noise and noise['root_seed']!=seed):raise ValueError('noise seed mismatch')
        if level=='infinite':
            if config['noise'] or noise:raise ValueError('baseline unexpectedly noisy')
        elif not noise or noise!=fit['noise_parameters'] or noise['metabolite']['requested_snr']!=int(level) or noise['water']['requested_snr']!=1000:
            raise ValueError('noise metadata mismatch')
        chosen=not smoke or (level=='infinite' or (level in ['50','10'] and r<2))
        if chosen:records.append(dict(snr=level,realization=r,seed=seed,source=str(directory.resolve()),fsl=fit,noise=noise,
            hashes={str(p.resolve()):hashes[str(p.resolve())] for p in directory.iterdir() if p.is_file()},
            calibration_identity=root['calibration_identity'],truth=gt['voxel_average_metabolite_concentrations_mM']))
    # Protect all earlier LCModel diagnostics, both basis inputs and calibration metadata.
    for folder in [BASE,Path('outputs/lcmodel_fit_bias_diagnostics'),Path('outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative')]:
        for path in folder.rglob('*'):
            if path.is_file():hashes[str(path.resolve())]=digest(path)
    hashes[str(BINARY)]=digest(BINARY)
    return records,hashes,root


def prepare_input(record,directory):
    checks={};signals={};source=Path(record['source'])
    for name,target in [('metabolites','metabolites.RAW'),('water_reference','water_reference.H2O')]:
        reader=mrs_io.read_FID(str(source/f'{name}.nii.gz'));mrs=reader.mrs();fid=mrs.FID
        if len(fid)!=2048 or not np.isclose(mrs.dwellTime,.0005) or not np.isclose(mrs.centralFrequency/1e6,127.731):raise ValueError('NIfTI-MRS grid mismatch')
        path=directory/target
        if not path.exists():
            temporary=directory/(target+'.tmp')
            if temporary.exists():temporary.unlink()  # Interrupted conversion of this output only.
            write_raw(temporary,fid,frequency_mhz=127.731);temporary.replace(path)
        stored,header=parse_raw(path);checks[name]=array_errors(stored.conj(),fid)
        if not np.array_equal(stored.conj(),fid) or header['frequency_mhz']!=127.731:raise ValueError('converted RAW identity mismatch')
        signals[name]=fid
        if record['noise']:
            clean=mrs_io.read_FID(str(source/f'{name}_clean.nii.gz')).mrs().FID
            spectrum=np.fft.fft(fid-clean,norm='ortho')
            sd=float(np.sqrt((spectrum.real.var(ddof=1)+spectrum.imag.var(ddof=1))/2))
            metadata=record['noise']['metabolite' if name=='metabolites' else 'water']
            measured=metadata['reference_peak']/sd
            if not np.isclose(measured,metadata['realized_snr'],rtol=1e-11):raise ValueError('saved realized SNR inconsistent with actual noise')
            checks[name].update(realized_snr=measured,target_snr=metadata['requested_snr'],sigma_metadata=metadata['sigma_per_component'])
    raw_met,_=parse_raw(directory/'metabolites.RAW');raw_water,_=parse_raw(directory/'water_reference.H2O')
    checks['initial_ratio_error']=float(abs((raw_met[0]/raw_water[0]).conjugate()-signals['metabolites'][0]/signals['water_reference'][0]))
    atomic_json(directory/'conversion_validation.json',checks)
    return checks


def execute_case(record,directory,protocol,config_hash,timeout):
    directory.mkdir(parents=True,exist_ok=True)
    conversion=prepare_input(record,directory)
    signature=dict(config_hash=config_hash,source_hashes=record['hashes'],converted_hashes={name:digest(directory/name) for name in ['metabolites.RAW','water_reference.H2O']})
    identity=directory/'identity.json'
    if identity.exists() and json.loads(identity.read_text())!=signature:raise ValueError('incompatible case cache identity')
    atomic_json(identity,signature)
    completion=directory/'completion.json'
    if completion.exists():
        previous=json.loads(completion.read_text())
        if previous['signature']!=signature:raise ValueError('incompatible completion signature')
        for path,expected in previous['artifact_hashes'].items():
            if not Path(path).exists() or digest(path)!=expected:raise ValueError('cached artifact changed')
        if previous['status']=='complete' and parse_fit(Path(previous['attempt_directory'])/'lcmodel_outputs')!=previous['fit']:
            raise ValueError('cached parsed fit disagrees with retained outputs')
        if previous['status'] in ['complete','failed']:return previous,True
    # Every actual attempt has a separate directory, including interrupted attempts.
    number=max([int(p.name.split('_')[1]) for p in directory.glob('attempt_*')]+[0])+1
    attempt=directory/f'attempt_{number:03d}';attempt.mkdir();(attempt/'lcmodel_outputs').mkdir()
    control=control_for_case(protocol,directory,attempt);(attempt/'lcmodel.CONTROL').write_text(control)
    command=f'cd {shlex.quote(str(attempt))} && {shlex.quote(str(BINARY))} < lcmodel.CONTROL'
    started=datetime.now(timezone.utc).isoformat();clock=time.monotonic()
    execution=dict(started_utc=started,command=command,timeout_s=timeout,finished=False)
    atomic_json(attempt/'execution.json',execution)
    try:
        p=subprocess.run([str(BINARY)],input=control,text=True,cwd=attempt,capture_output=True,timeout=timeout)
        execution.update(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,finished=True,timed_out=False)
    except subprocess.TimeoutExpired as error:
        def decode(v):return v.decode(errors='replace') if isinstance(v,bytes) else v or ''
        execution.update(returncode=None,stdout=decode(error.stdout),stderr=decode(error.stderr),finished=True,timed_out=True)
    except OSError as error:execution.update(returncode=None,stdout='',stderr=str(error),finished=True,timed_out=False)
    execution['elapsed_s']=time.monotonic()-clock;atomic_json(attempt/'execution.json',execution)
    (attempt/'stdout.txt').write_text(execution['stdout']);(attempt/'stderr.txt').write_text(execution['stderr'])
    fit=dict(valid=False,failure='incomplete outputs',diagnostics=[])
    try:
        fit=parse_fit(attempt/'lcmodel_outputs')
        if execution['returncode']!=0:fit.update(valid=False,failure='nonzero execution or timeout')
    except (OSError,ValueError,IndexError,KeyError,ZeroDivisionError) as error:
        fit.update(failure=f'{type(error).__name__}: {error}')
    status='complete' if fit['valid'] else ('interrupted' if execution.get('timed_out') else 'failed')
    if not fit['valid'] and execution['returncode']==0 and not fit.get('fatal'):status='incomplete'
    artifacts={str(p):digest(p) for p in attempt.rglob('*') if p.is_file()}
    result=dict(signature=signature,status=status,fit=fit,execution=execution,attempt=number,attempt_directory=str(attempt),conversion=conversion,artifact_hashes=artifacts)
    atomic_json(completion,result)
    return result,False


def write_csv(path,rows):
    if not rows:raise ValueError('no rows to write')
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def summarize(records,clean_fsl,clean_lcm):
    rows=[];summaries=[];paired=[]
    for record in records:
        fit=record['result']['fit'];fsl=record['fsl'];noise=record['noise']
        for name in METABOLITES:
            truth=sum(record['truth'][n] for n in name.split('+'));base_lcm=sum(clean_lcm['concentrations'][n]['concentration'] for n in name.split('+'))
            v=fsl['coefficients'][name];lcm=fit['values'][name] if fit['valid'] else None
            source_sd=v.get('reported_SD_mM');source=v['estimate_mM'];standard=v.get('standard_FSL_mM')
            rows.append(dict(snr=record['snr'],realization=record['realization'],seed=record['seed'],metabolite=name,truth_mM=truth,
                LCModel_valid=fit['valid'],LCModel_failure=fit.get('failure'),LCModel_mM=lcm['concentration'] if lcm else None,
                LCModel_clean_mM=base_lcm,LCModel_percent_SD=(lcm.get('reported_group_percent_SD') if name=='GPC+PCh' else lcm['percent_SD']) if lcm else None,
                LCModel_uncertainty_mM=lcm['uncertainty_mM'] if lcm else None,
                LCModel_displayed_SD_product_mM=lcm['displayed_SD_product_mM'] if lcm else None,
                LCModel_zero_percent_SD=lcm['zero_percent_SD'] if lcm else None,
                FSL_source_mM=source,FSL_standard_mM=standard,FSL_source_SD_mM=source_sd,
                FSL_standard_conditional_SD_mM=source_sd*abs(standard/source) if source_sd is not None and standard is not None and source!=0 else None,
                FSL_source_clean_mM=clean_fsl['coefficients'][name]['estimate_mM'],FSL_standard_clean_mM=clean_fsl['coefficients'][name]['standard_FSL_mM'],
                FSL_optimizer_success=fsl['optimizer_success'],FSL_quantification_success=fsl['quantification_success'],
                LCModel_real_residual=fit.get('relative_real_residual'),FSL_complex_residual=fsl['relative_complex_residual'],
                FSL_water_fitted_amplitude=fsl['water_fit']['fitted_amplitude'],FSL_water_complex_residual=fsl['water_fit']['relative_complex_residual'],
                LCModel_reported_SNR=fit.get('reported_SNR'),LCModel_FWHM_ppm=fit.get('FWHM_ppm'),LCModel_shift_ppm=fit.get('data_shift_ppm'),
                LCModel_phase_zero_deg=fit.get('phase_zero_deg'),LCModel_phase_first_deg_ppm=fit.get('phase_first_deg_ppm'),
                LCModel_baseline_norm=fit.get('baseline_relative_norm'),LCModel_convergence_diagnostic=fit.get('convergence_diagnostic'),
                LCModel_FCALIB=fit.get('FCALIB'),LCModel_water_area=fit.get('water_area'),
                target_SNR=noise['metabolite']['requested_snr'] if noise else None,
                realized_SNR=record['result']['conversion']['metabolites'].get('realized_snr'),
                realized_water_SNR=record['result']['conversion']['water_reference'].get('realized_snr'),
                diagnostics=';'.join(d['routine']+str(d['code']) for d in fit['diagnostics'])))
    modes=[('LCModel_standard','LCModel_mM','LCModel_clean_mM','LCModel_uncertainty_mM'),
        ('FSL_source_equivalent','FSL_source_mM','FSL_source_clean_mM','FSL_source_SD_mM'),
        ('FSL_standard','FSL_standard_mM','FSL_standard_clean_mM','FSL_standard_conditional_SD_mM')]
    for level in dict.fromkeys(r['snr'] for r in rows):
        for index,name in enumerate(METABOLITES):
            selected=[r for r in rows if r['snr']==level and r['metabolite']==name]
            for toolbox,key,bkey,ukey in modes:
                returned=[r for r in selected if r[key] is not None]
                failures=sum(not r['LCModel_valid'] if toolbox=='LCModel_standard' else (not r['FSL_optimizer_success'] or (toolbox=='FSL_standard' and not r['FSL_quantification_success'])) for r in selected)
                stats=descriptive([r[key] for r in returned],selected[0]['truth_mM'],selected[0][bkey],[r[ukey] for r in returned],len(selected),failures)
                stats.update(snr=level,metabolite=name,toolbox=toolbox)
                baseline=selected[0][bkey]
                stats['SD_percent_clean_baseline']=100*stats['SD_mM']/baseline if stats.get('SD_mM') is not None else None
                stats['noise_RMSE_percent_clean_baseline']=100*stats['noise_RMSE_mM']/baseline if stats.get('noise_RMSE_mM') is not None else None
                if name=='GPC+PCh':
                    individual=[[r for r in rows if r['snr']==level and r['metabolite']==n] for n in ['GPC','PCh']]
                    lookup=[{r['realization']:r[key] for r in group if r[key] is not None} for group in individual]
                    common=sorted(set(lookup[0])&set(lookup[1]))
                    stats['empirical_GPC_PCh_correlation']=float(np.corrcoef([[group[i] for i in common] for group in lookup])[0,1]) if len(common)>1 else None
                if toolbox=='LCModel_standard':
                    stats['rounded_zero_SD_fraction']=float(np.mean([r['LCModel_zero_percent_SD'] for r in returned])) if returned else None
                    percentages=[r['LCModel_percent_SD'] for r in returned if r['LCModel_percent_SD'] is not None]
                    stats['median_reported_percent_SD']=float(np.median(percentages)) if percentages else None
                    stats['capped_percent_SD_count']=sum(p>=999 for p in percentages)
                    stats['median_displayed_SD_product_mM']=float(np.median([r['LCModel_displayed_SD_product_mM'] for r in returned if r['LCModel_displayed_SD_product_mM'] is not None])) if returned else None
                summaries.append(stats)
            if level!='infinite':
                for comparator,key,bkey in [('FSL_standard','FSL_standard_mM','FSL_standard_clean_mM'),('FSL_source_equivalent','FSL_source_mM','FSL_source_clean_mM')]:
                    pairs=[r for r in selected if r['LCModel_mM'] is not None and r[key] is not None]
                    stats=paired_bootstrap([r['LCModel_mM'] for r in pairs],[r[key] for r in pairs],selected[0]['LCModel_clean_mM'],selected[0][bkey],selected[0]['truth_mM'],seed=77000+100*LEVELS.index(int(level))+index) if pairs else dict(pairs=0)
                    paired.append(dict(snr=level,metabolite=name,comparator=comparator,**stats))
    return rows,summaries,paired


def reports(out,records):
    requested=[r for r in records if r['realization']==0 and r['snr'] in ['infinite','50','25','10','5']]
    valid=[r for r in records if r['snr']!='infinite' and r['result']['fit']['valid']]
    if valid:
        worst=max(valid,key=lambda r:max(abs(r['result']['fit']['values'][n]['concentration']/sum(r['truth'][k] for k in n.split('+'))-1) for n in METABOLITES))
        if worst not in requested:requested.append(worst)
    output=[]
    for r in requested:
        directory=Path(r['result']['attempt_directory'])/'lcmodel_outputs';ps=directory/'fit.PS';pdf=directory/'fit.pdf'
        if ps.exists() and not pdf.exists():
            subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-sDEVICE=pdfwrite',f'-sOutputFile={pdf}',str(ps)],check=True,capture_output=True)
        html=Path(r['source'])/'fsl_mrs_results/report.html'
        if (not pdf.exists() or pdf.stat().st_size==0) and r['result']['fit']['valid']:raise RuntimeError('valid fit representative report missing')
        output.append(dict(snr=r['snr'],realization=r['realization'],PDF=str(pdf) if pdf.exists() else None,FSL_HTML=str(html) if html.exists() else None))
    atomic_json(out/'reports.json',output)
    return output


def plots(out,rows,summary,paired):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for metric,title in [('nominal_bias_mM','Bias relative to nominal truth (mM)'),('SD_mM','Empirical SD (mM)'),('nominal_RMSE_mM','RMSE relative to nominal truth (mM)'),('noise_RMSE_mM','RMSE relative to each toolbox clean baseline (mM)'),('noise_bias_mM','Noise-induced bias (mM)'),('nominal_poor_fraction','Fraction with >20% nominal error'),('SD_percent_clean_baseline','Empirical SD / clean toolbox baseline (%)')]:
        fig,axes=plt.subplots(2,3,figsize=(13,7))
        for name,ax in zip(METABOLITES,axes.flat):
            for toolbox in ['LCModel_standard','FSL_standard','FSL_source_equivalent']:
                selected=[s for s in summary if s['snr']!='infinite' and s['metabolite']==name and s['toolbox']==toolbox]
                ax.plot([int(s['snr']) for s in selected],[s.get(metric) for s in selected],'-o',label=toolbox,ms=3)
            ax.set_title(name);ax.set_xscale('log');ax.invert_xaxis();ax.grid(alpha=.3);ax.set_xlabel('Simulator SNR')
        axes.flat[0].legend(fontsize=7);fig.suptitle(title);fig.tight_layout();fig.savefig(out/f'{metric}.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,7))
    for name,ax in zip(METABOLITES,axes.flat):
        selected=[s for s in summary if s['snr']!='infinite' and s['metabolite']==name and s['toolbox']=='LCModel_standard']
        x=[int(s['snr']) for s in selected]
        ax.plot(x,[s.get('SD_mM') for s in selected],'-o',label='Empirical SD')
        ax.plot(x,[s.get('median_reported_uncertainty_mM') for s in selected],'-o',label='Median positive reported SD')
        ax.set_title(name);ax.set_xscale('log');ax.invert_xaxis();ax.grid(alpha=.3)
    axes.flat[0].legend(fontsize=7);fig.suptitle('LCModel rounded %SD products versus empirical variability');fig.tight_layout();fig.savefig(out/'uncertainty_comparison.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,7))
    for name,ax in zip(METABOLITES,axes.flat):
        selected=[p for p in paired if p['metabolite']==name and p['comparator']=='FSL_standard']
        x=[int(p['snr']) for p in selected];y=[p['noise_absolute_error_mean'] for p in selected]
        ax.errorbar(x,y,yerr=[[v-p['noise_absolute_error_CI_low'] for v,p in zip(y,selected)],[p['noise_absolute_error_CI_high']-v for v,p in zip(y,selected)]],fmt='o-')
        ax.axhline(0,c='k',lw=.7);ax.set_title(name);ax.set_xscale('log');ax.invert_xaxis();ax.grid(alpha=.3)
    fig.suptitle('Paired noise absolute error: LCModel minus FSL standard, 95% bootstrap CI (mM)');fig.tight_layout();fig.savefig(out/'paired_error_comparison.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,7))
    for name,ax in zip(METABOLITES,axes.flat):
        for level in LEVELS:
            selected=[r for r in rows if r['metabolite']==name and r['snr']==str(level) and r['LCModel_valid']]
            ax.scatter([r['FSL_standard_mM']-r['FSL_standard_clean_mM'] for r in selected],
                       [r['LCModel_mM']-r['LCModel_clean_mM'] for r in selected],s=12,label=str(level),alpha=.7)
        ax.set_title(name);ax.set_xlabel('FSL standard minus clean (mM)');ax.set_ylabel('LCModel minus clean (mM)');ax.grid(alpha=.3)
    axes.flat[0].legend(title='SNR',fontsize=7);fig.tight_layout();fig.savefig(out/'paired_deviations.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,key,title in zip(axes,['LCModel_real_residual','FSL_complex_residual'],['LCModel phased-real spectral window','FSL-MRS full complex FID']):
        for level in LEVELS:
            selected=[r[key] for r in rows if r['metabolite']=='NAA' and r['snr']==str(level) and r[key] is not None]
            ax.scatter([level]*len(selected),selected,s=12,alpha=.6)
        ax.set_title(title);ax.set_xlabel('Simulator SNR');ax.set_ylabel('Relative residual (different domains)');ax.set_xscale('log');ax.invert_xaxis();ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(out/'residuals.png',dpi=140);plt.close(fig)


def run(out=OUTPUT,*,smoke=False,timeout=180):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    records,source_hashes,root=inventory(smoke)
    protocol=(BASE/'lcmodel.CONTROL').read_text()
    config=dict(smoke=smoke,timeout_s=timeout,CONTROL_sha256=text_hash(protocol),canonical_CONTROL_sha256=text_hash(canonical(protocol)),
        python_version=sys.version,numpy_version=np.__version__,FSL_MRS_version=root['config']['FSL_MRS_version'],
        LCModel_version=json.loads((BASE/'fitting_results.json').read_text())['version'],
        basis_sha256=digest(BASE/'lcmodel_compatible.BASIS'),binary_sha256=digest(BINARY),source_hashes=source_hashes,
        conversion_sha256=digest(Path('utils/lcmodel_quantitative_io.py')),
        execution_sha256=text_hash(inspect.getsource(prepare_input)+inspect.getsource(execute_case)+inspect.getsource(control_for_case)+inspect.getsource(parse_fit)),
        calibration_identity=root['calibration_identity'],bootstrap_seed_recipe='77000+100*level_index+metabolite_index',bootstrap_draws=10000)
    config_path=out/'frozen_configuration.json'
    if config_path.exists() and json.loads(config_path.read_text())!=config:raise ValueError('incompatible frozen experiment configuration')
    atomic_json(config_path,config);(out/'frozen_protocol.CONTROL').write_text(protocol)
    results=[];hits=0;commands=[]
    for record in records:
        directory=out/f'SNR_{record["snr"]}'/f'realization_{record["realization"]:02d}'
        result,hit=execute_case(record,directory,protocol,text_hash(json.dumps(config,sort_keys=True)),timeout)
        record=dict(record,result=result);results.append(record);hits+=int(hit)
        commands.append(result['execution']['command'])
        print(f'SNR {record["snr"]} realization {record["realization"]:02d}: {result["status"]}; cache={hit}',flush=True)
    clean_fsl=root['records'][0]['fit'];clean_lcm=json.loads((BASE/'fitting_results.json').read_text())
    baseline=results[0]['result']['fit']
    if not baseline['valid'] or baseline['concentrations']!=clean_lcm['concentrations']:raise RuntimeError('frozen-protocol clean baseline did not reproduce Phase 7b')
    rows,summary,paired=summarize(results,clean_fsl,clean_lcm)
    write_csv(out/'per_realization.csv',rows);write_csv(out/'summary.csv',summary);write_csv(out/'paired_comparison.csv',paired)
    water_rows=[]
    for level in dict.fromkeys(r['snr'] for r in results):
        selected=[r for r in results if r['snr']==level and r['result']['fit']['valid']]
        areas=np.array([r['result']['fit']['water_area'] for r in selected])
        factors=np.array([r['result']['fit']['FCALIB'] for r in selected])
        naa=np.array([r['result']['fit']['values']['NAA']['concentration'] for r in selected])
        amplitudes=np.array([r['fsl']['water_fit']['fitted_amplitude'] for r in selected])
        water_rows.append(dict(snr=level,count=len(selected),LCModel_mean_water_area=float(areas.mean()) if len(areas) else None,
            LCModel_water_area_CV_percent=float(100*areas.std(ddof=1)/areas.mean()) if len(areas)>1 else None,
            LCModel_FCALIB_mean=float(factors.mean()) if len(factors) else None,
            NAA_FCALIB_correlation=float(np.corrcoef(naa,factors)[0,1]) if len(naa)>1 else None,
            FSL_mean_water_amplitude=float(amplitudes.mean()) if len(amplitudes) else None,
            FSL_water_amplitude_CV_percent=float(100*amplitudes.std(ddof=1)/amplitudes.mean()) if len(amplitudes)>1 else None))
    write_csv(out/'water_scaling_comparison.csv',water_rows)
    report_list=reports(out,results);plots(out,rows,summary,paired)
    warnings={}
    for record in results:
        for diagnostic in record['result']['fit']['diagnostics']:
            key=diagnostic['severity']+':'+diagnostic['routine']+str(diagnostic['code']);warnings.setdefault(key,[]).append(dict(snr=record['snr'],realization=record['realization'],count=diagnostic['count']))
    atomic_json(out/'warning_summary.json',warnings)
    failure_summary=[dict(snr=level,total=sum(r['snr']==level for r in results),
        failed=[r['realization'] for r in results if r['snr']==level and not r['result']['fit']['valid']],
        numerical_diagnostics=[r['realization'] for r in results if r['snr']==level and r['result']['fit'].get('convergence_diagnostic')])
        for level in dict.fromkeys(r['snr'] for r in results)]
    atomic_json(out/'failure_summary.json',failure_summary)
    (out/'commands.txt').write_text('\n'.join(commands)+'\n')
    if source_hashes!={p:digest(p) for p in source_hashes}:raise RuntimeError('protected Phase 7a/7b/calibration data changed')
    attempts=sum(len(list((out/f'SNR_{r["snr"]}'/f'realization_{r["realization"]:02d}').glob('attempt_*'))) for r in results)
    result=dict(records=results,summary=summary,paired=paired,reports=report_list,cache_hits=hits,
        fits=len(results),attempts=attempts,valid_fits=sum(r['result']['fit']['valid'] for r in results),
        failures=sum(not r['result']['fit']['valid'] for r in results),source_files_unchanged=True,
        analysis_implementation_sha256=digest(Path(__file__)),statistics_implementation_sha256=digest(Path('utils/lcmodel_noise_analysis.py')))
    atomic_json(out/'results.json',result)
    print(f'{len(results)} acquisitions; {result["valid_fits"]} valid; {result["failures"]} failures; {hits} cache hits',flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--smoke',action='store_true');p.add_argument('--timeout',type=int,default=180)
    a=p.parse_args();run(a.output,smoke=a.smoke,timeout=a.timeout)

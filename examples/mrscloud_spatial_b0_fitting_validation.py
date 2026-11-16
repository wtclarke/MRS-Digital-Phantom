"""Phase 6g: immutable spatial inputs, actual CLI reports and mismatch diagnostics."""
import argparse
import csv
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import warnings

import fsl_mrs
import numpy as np
from scipy.optimize import nnls
from fsl_mrs.core.basis import Basis
from fsl_mrs.utils import fitting, mrs_io

from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral, relaxation
from examples.fsl_mrs_water_fit_diagnostics import fit_water
from examples.mrscloud_mixed_tissue_fitting_validation import quant_info
from examples.mrscloud_quantitative_voxel import BUNDLE
from simulation.mrscloud_full_basis import load_library, sha256

SOURCE=Path('outputs/mrscloud_quantitative_spatial/validated')
CASES=['zero_B0','constant_B0','symmetric_B0','mixed_relaxation_B0']
NAMES=['NAA','Cr','GPC','PCh','Glu']
INPUTS=['metabolites.nii.gz','water_reference.nii.gz','ground_truth.json','simulation_config.json']


def spectrum(fid):
    """Independent FSL-compatible first-point weighting and orthonormal FFT."""
    weighted=np.array(fid,complex,copy=True);weighted[0]*=.5
    return np.fft.fftshift(np.fft.fft(weighted,axis=0,norm='ortho'),axes=0)


def window_indices(n, limits):
    ppm=np.fft.fftshift(np.fft.fftfreq(n,.0005))/127.731+4.68
    lo,hi=sorted(np.argmin(abs(ppm-p)) for p in limits)
    return slice(lo,hi+1)


def projection(columns, data, limits):
    """Real nonnegative coefficients on a fixed complex spectral design."""
    selection=window_indices(len(data),limits)
    matrix=columns[selection];observed=spectrum(data)[selection]
    coef,_=nnls(np.concatenate([matrix.real,matrix.imag]),np.r_[observed.real,observed.imag])
    return coef,float(np.linalg.norm(matrix@coef-observed)/np.linalg.norm(observed))


def independent_signals(truth, config, responses, water):
    """Cell equations use declared T1/T2, not simulator envelopes or factors."""
    spatial=truth['spatial_integration'];t=np.arange(len(water))/2000
    components={n:np.zeros(len(t),complex) for n in NAMES};water_signal=np.zeros(len(t),complex)
    effective=dict.fromkeys(NAMES,0.);cw=0.;water_moment=0.
    cfg=config.get('relaxation');gain=config['receiver_factor']
    def r(tissue,name=None):
        if not cfg or cfg['mode']=='none':return 1.
        p=cfg['metabolites'][tissue][name] if name else cfg['water'][tissue]
        return relaxation(cfg['TE_s'],cfg['TR_s'],p['T1_s'],p['T2_s'])
    for tissue,w,offset in zip(spatial['cell_labels'],spatial['normalized_weights'],spatial['b0_offsets_hz']):
        mod=np.exp(2j*np.pi*offset*t)
        for n in NAMES:
            c=w*truth['tissue_specific_metabolite_concentrations_mM'][tissue][n]*r(tissue,n)
            effective[n]+=c;components[n]+=gain*.001*c*responses[n]*mod
        c=w*truth['tissue_specific_water_concentrations_mol_l'][tissue]*r(tissue)
        cw+=c;water_moment+=c*offset;water_signal+=gain*c*water*mod
    water_centroid=config['water_frequency_offset_hz']+water_moment/cw
    return components,water_signal,effective,cw,water_centroid


def nuisance_design(responses, params, *, frequency_hz=0.):
    t=np.arange(len(next(iter(responses.values()))))/2000
    frequency=np.fft.fftshift(np.fft.fftfreq(len(t),.0005))
    decay=np.exp(-(params.get('gamma_0',0)+1j*params.get('eps_0',0))*t+2j*np.pi*frequency_hz*t)
    phase=np.exp(-1j*(params.get('Phi0',0)+params.get('Phi1',0)*frequency))
    return spectrum(np.column_stack([responses[n]*decay for n in NAMES]))*phase[:,None]


def aligned_cli_basis(directory):
    """Known +37 Hz nuisance alignment only; preserve every target basis sample magnitude."""
    original=mrs_io.read_basis(str(BUNDLE/'metabolite_basis.BASIS'))
    unshifted=original.get_formatted_basis(2000,2048)
    shifted=unshifted*np.exp(2j*np.pi*37*np.arange(2048)/2000)[:,None]
    headers=[dict(dwelltime=.0005,bandwidth=2000,centralFrequency=127.731e6,
                  fwhm=width,nucleus='1H',centralShift=4.68) for width in original.basis_fwhm]
    Basis(shifted,list(original.names),headers).save(directory,info_str='Phase 6g known +37 Hz target-grid fitting alignment; no amplitude normalization')
    reloaded=mrs_io.read_basis(str(directory))
    error=max(float(abs(reloaded.original_basis_array[:,reloaded.names.index(n)]-shifted[:,original.names.index(n)]).max()) for n in original.names)
    return dict(applied_frequency_hz=37.,export_reload_max_complex_error=error,
                magnitude_max_error=float(abs(abs(shifted)-abs(unshifted)).max()),
                source_basis_sha256=sha256(BUNDLE/'metabolite_basis.BASIS'),
                purpose='CLI has no x0 option; fitting-coordinate alignment, not simulation or calibration change')


def fit_case(directory, responses, water):
    truth=json.loads((directory/'ground_truth.json').read_text());config=json.loads((directory/'simulation_config.json').read_text())
    actual=mrs_io.read_FID(str(directory/'metabolites.nii.gz')).mrs(
        basis=mrs_io.read_basis(str(BUNDLE/'metabolite_basis.BASIS')),ref_data=str(directory/'water_reference.nii.gz'))
    actual.keep=NAMES;actual.rescaleForFitting()
    components,known_water,effective,cw,water_centroid=independent_signals(truth,config,responses,water)
    original_data=actual.FID/actual.scaling['FID'];original_water=actual.H2O/actual.scaling['FID']
    kwargs=dict(method='Newton',model='lorentzian',baseline_order=-1,ppmlim=(.2,4.2),
                metab_groups=[0]*5,capture_minimize_output=True)
    with warnings.catch_warnings(record=True) as caught:
        fit=fitting.fit_FSLModel(actual,**kwargs)
    default=dict(parameters=dict(zip(fit.fitResults.columns,fit.fitResults.iloc[0].tolist())),
        relative_complex_residual=float(np.linalg.norm(actual.FID-fit.predictedFID(actual))/np.linalg.norm(actual.FID)),
        success=bool(fit.scipy_minimize_result.success),message=str(fit.scipy_minimize_result.message))
    offsets=truth['spatial_integration']['b0_offsets_hz']
    known_shift=offsets[0] if len(set(offsets))==1 else truth['spatial_integration']['b0_weighted_mean_hz']
    fixed_coef,fixed_resid=projection(nuisance_design(responses,{},frequency_hz=known_shift),original_data,(.2,4.2))
    if known_shift==37. and len(set(offsets))==1:
        # Public x0: known acquisition shift and data-derived linear amplitudes.
        initial=fit.params.copy();order=list(fit.fitResults.columns)
        for n,c in zip(NAMES,fixed_coef):initial[order.index(n)]=c*actual.scaling['FID']/actual.scaling['basis']
        for n,v in [('gamma_0',0.),('eps_0',-2*np.pi*37),('Phi0',0.),('Phi1',0.)]:initial[order.index(n)]=v
        fit=fitting.fit_FSLModel(actual,x0=initial,**kwargs)
    params=dict(zip(fit.fitResults.columns,fit.fitResults.iloc[0].tolist()))
    w,_,water_model=fit_water(actual)
    sd,sb=actual.scaling['FID'],actual.scaling['basis'];gain=config['receiver_factor']
    w.update(expected_initial_amplitude=float(abs(known_water[0])),
             relative_amplitude_error=float(w['fitted_amplitude']/abs(known_water[0])-1))
    projected,projected_resid=projection(nuisance_design(responses,params),original_data,(.2,4.2))
    t=np.arange(len(water))/2000;ppm=np.fft.fftshift(np.fft.fftfreq(len(t),.0005))/127.731+4.68
    amp,gamma,sigma,omega,phi=w['parameters']
    water_unit=np.exp(-t*(gamma+t*sigma+1j*omega)+1j*phi)
    water_projection,water_projection_resid=projection(spectrum(water_unit)[:,None],original_water,(1.65,7.65))
    # Correctly specified spatial design is an oracle diagnostic, not a standard fitting basis.
    oracle_design=spectrum(np.column_stack([components[n]/(gain*.001*effective[n]) for n in NAMES]))
    oracle_coef,oracle_resid=projection(oracle_design,original_data,(.2,4.2))
    water_oracle,water_oracle_resid=projection(spectrum(known_water/abs(known_water[0]))[:,None],original_water,(1.65,7.65))
    q=quant_info(actual,truth,False)
    before=fit.params.copy();fit.calculateConcScaling(actual,quant_info=q,internal_reference=['NAA'])
    np.testing.assert_array_equal(before,fit.params)
    correction=q.relax_corr_water_molar*q.csf_corr_molar*q.relax_corr_metab*(2/3)
    true_water_aligned=known_water*np.exp(-2j*np.pi*water_centroid*t-1j*np.angle(known_water[0]))
    j_true_water=spectral_integral(true_water_aligned,ppm,(1.65,7.65))
    j_true_naa=spectral_integral(components['NAA'],ppm,(1.8,2.2))
    raw={n:float(fit.getConc(metab=n)*sb/sd) for n in NAMES}
    corrected_water=w['fitted_amplitude']*np.exp(-gamma*t-sigma*t*t)
    naa_model=raw['NAA']*responses['NAA']*np.exp(-(params['gamma_0']+1j*params['eps_0'])*t)
    jm=spectral_integral(naa_model,ppm,(1.8,2.2));jw=spectral_integral(corrected_water,ppm,(1.65,7.65))
    # Reconstruction and projection are explicitly conditional on fitted nuisance parameters.
    rows={}
    for i,n in enumerate(NAMES):
        reported=float(fit.getConc(scaling='molarity',metab=n))
        rows[n]=dict(nominal_voxel_mM=truth['voxel_average_metabolite_concentrations_mM'][n],
            tissue_only_mM=truth['tissue_only_metabolite_concentrations_mM'][n],true_effective_mM=effective[n],
            expected_zero_B0_internal=effective[n]*.001*gain*sd/sb,fitted_internal=float(fit.getConc(metab=n)),
            fitted_effective_mM=raw[n]*1000/gain,relative_coefficient_error=raw[n]/(effective[n]*.001*gain)-1,
            fixed_known_shift_projection_mM=float(fixed_coef[i]*1000/gain),
            fixed_fitted_nuisance_projection_mM=float(projected[i]*1000/gain),
            oracle_spatial_projection_mM=float(oracle_coef[i]*1000/gain),
            standard_FSL_mM=reported,source_response_reference_mM=effective[n]/effective['NAA']*j_true_naa/j_true_water*correction,
            independently_reconstructed_FSL_mM=raw[n]/raw['NAA']*jm/jw*correction)
    fitted=fit.predictedFID(actual)/sd
    np.savez(directory/'diagnostic_arrays.npz',time_s=t,ppm=ppm,data=original_data,fitted=fitted,
        independent_data=sum(components.values()),water=original_water,water_fitted=water_model/sd,
        water_independent=known_water,water_true_aligned=true_water_aligned,water_integrated_model=corrected_water)
    result=dict(default_fit=default,parameters=params,scaling=actual.scaling,coefficients=rows,
        optimizer_success=bool(fit.scipy_minimize_result.success),optimizer_message=str(fit.scipy_minimize_result.message),
        relative_complex_residual=float(np.linalg.norm(original_data-fitted)/np.linalg.norm(original_data)),
        water_fit=w,nominal_effective_water_M=truth['effective_voxel_water_concentration_mol_l'],
        relaxation_weighted_water_M=cw,true_water_centroid_hz=water_centroid,
        independent_metabolite_max_error=float(abs(sum(components.values())-original_data).max()),
        independent_water_max_error=float(abs(known_water-original_water).max()),
        fixed_known_shift_projection_residual=fixed_resid,fixed_fitted_nuisance_projection_residual=projected_resid,
        oracle_spatial_projection_residual=oracle_resid,
        water_fixed_nuisance_projection_amplitude=float(water_projection[0]),water_fixed_nuisance_projection_residual=water_projection_resid,
        water_oracle_amplitude=float(water_oracle[0]),water_oracle_residual=water_oracle_resid,
        source_response_reference_is_exact_model_prediction=len(set(offsets))==1,
        true_response_integrals=dict(NAA=j_true_naa,water_centroid_aligned=j_true_water),
        fitted_response_integrals=dict(NAA=jm,water_corrected=jw),quantification_factor=correction,
        quantification_water_relaxation_mM=q.relax_corr_water_molar,quantification_R_M=q.R_M,
        csf_correction=q.csf_corr_molar,warnings=[str(x.message) for x in caught])
    return result


def run_validation(output_dir):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=False)
    tracked=[p for p in BUNDLE.iterdir() if p.is_file()]+[SOURCE/c/n for c in CASES for n in INPUTS]
    hashes={str(p):sha256(p) for p in tracked}
    responses,water,calibration=load_library(BUNDLE,selected_names=NAMES)
    help_result=subprocess.run(['pixi','run','fsl_mrs','--help'],capture_output=True,text=True,check=True)
    (out/'cli_help.txt').write_text(help_result.stdout)
    basis_alignment=aligned_cli_basis(out/'constant_shifted_fitting_basis')
    cases={};commands=[]
    for name in CASES:
        directory=out/name;directory.mkdir()
        for filename in INPUTS:shutil.copyfile(SOURCE/name/filename,directory/filename)
        cases[name]=fit_case(directory,responses,water)
        truth=json.loads((directory/'ground_truth.json').read_text());f=truth['tissue_fractions']
        basis=out/'constant_shifted_fitting_basis' if name=='constant_B0' else BUNDLE/'metabolite_basis.BASIS'
        command=['pixi','run','fsl_mrs','--data',str(directory/'metabolites.nii.gz'),'--h2o',str(directory/'water_reference.nii.gz'),
            '--basis',str(basis),'--output',str(directory/'fsl_mrs_results'),'--keep',*NAMES,'--algo','Newton',
            '--lorentzian','--baseline','off','--ppmlim','.2','4.2','--no_conj_basis','--TE','30','--TR','2',
            '--tissue_frac',*[str(f[t]) for t in ['WM','GM','CSF']],'--internal_ref','NAA',
            '--wref_metabolite','NAA','--ref_protons','3','--ref_int_limits','1.8','2.2','--report']
        commands.append(shlex.join(command));(out/'cli_commands.txt').write_text('\n'.join(commands)+'\n')
        completed=subprocess.run(command,capture_output=True,text=True)
        (directory/'cli.log').write_text(completed.stdout+completed.stderr);completed.check_returncode()
        report=directory/'fsl_mrs_results/report.html'
        if report.stat().st_size<1000 or 'Water referencing' not in report.read_text():raise RuntimeError('CLI HTML report missing or empty')
        with (directory/'fsl_mrs_results/concentrations.csv').open() as stream:records=list(csv.reader(stream))
        cases[name]['CLI_mM']={row[0]:float(row[4]) for row in records[2:] if row[0] in NAMES}
        cases[name]['CLI_raw_coefficient_mol_l']={row[0]:float(row[1]) for row in records[2:] if row[0] in NAMES}
        cli_fit=mrs_io.read_FID(str(directory/'fsl_mrs_results/fit.nii.gz')).mrs().FID
        arrays=np.load(directory/'diagnostic_arrays.npz')
        cases[name]['CLI_fit_relative_complex_residual']=float(np.linalg.norm(cli_fit-arrays['data'])/np.linalg.norm(arrays['data']))
        cases[name]['CLI_API_max_complex_fit_difference']=float(abs(cli_fit-arrays['fitted']).max())
        cases[name]['CLI_basis_frequency_alignment_hz']=37. if name=='constant_B0' else 0.
    if hashes!={str(p):sha256(p) for p in tracked}:raise RuntimeError('Phase 6a/6f inputs changed')
    result=dict(FSL_MRS_version=fsl_mrs.__version__,cases=cases,calibration_identity=calibration['calibration_identity'],
        source_hashes=hashes,source_unchanged=True,constant_CLI_basis_alignment=basis_alignment,
        settings=dict(method='Newton',model='lorentzian',baseline='off',ppm=[.2,4.2],one_shared_group=True,
            TE_s=.03,TR_s=2.,water_reference='NAA 3 protons 1.8-2.2 ppm',water_window_ppm=[1.65,7.65],
            quantification='installed defaults with actual fractions; no gain correction',
            constant_API_initialization='public x0: known +37 Hz, NNLS amplitudes; original basis and unchanged bounds',
            constant_CLI_initialization='target-grid basis with common known +37 Hz; original input signals'))
    (out/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (out/'comparison.csv').open('w') as stream:
        first=next(iter(cases.values()))['coefficients'][NAMES[0]]
        writer=csv.DictWriter(stream,fieldnames=['case','metabolite',*first]);writer.writeheader()
        writer.writerows(dict(case=c,metabolite=n,**d) for c,v in cases.items() for n,d in v['coefficients'].items())
    plots(out,result)
    return result


def plots(out,result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(4,3,figsize=(16,13))
    for row,name in enumerate(CASES):
        a=np.load(out/name/'diagnostic_arrays.npz');ppm=a['ppm']
        for col,part in [(0,'real'),(1,'imag')]:
            ax=axes[row,col];ax.plot(ppm,getattr(spectrum(a['data']),part),label='data')
            ax.plot(ppm,getattr(spectrum(a['fitted']),part),label='fitted');ax.set_xlim(4.2,.2)
            ax.set_title(f'{name}: metabolite {part}');ax.legend()
        ax=axes[row,2]
        for key,style in [('water','-'),('water_fitted','--')]:
            spec=spectrum(a[key]);ax.plot(ppm,spec.real,style,label=key+' real');ax.plot(ppm,spec.imag,style,alpha=.6,label=key+' imag')
        ax.set_xlim(5.3,4.1);ax.set_title('Water: same original complex frame');ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(out/'complex_spectral_fits.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,8));x=np.arange(4)
    for n in NAMES:
        d=[result['cases'][c]['coefficients'][n] for c in CASES]
        axes[0,0].plot(x,[100*z['relative_coefficient_error'] for z in d],'-o',label=n)
        axes[1,0].plot(x,[z['standard_FSL_mM']-z['nominal_voxel_mM'] for z in d],'-o',label=n)
        axes[1,1].plot(x,[z['standard_FSL_mM']-z['source_response_reference_mM'] for z in d],'-o',label=n)
    data=[result['cases'][c] for c in CASES]
    axes[0,1].plot(x,[100*d['water_fit']['relative_amplitude_error'] for d in data],'-o')
    axes[0,2].plot(x,[d['relative_complex_residual'] for d in data],'-o',label='metabolites')
    axes[0,2].plot(x,[d['water_fit']['relative_complex_residual'] for d in data],'-o',label='water')
    axes[1,2].plot(x,[d['fixed_fitted_nuisance_projection_residual'] for d in data],'-o',label='fixed-nuisance projection')
    for ax,title in zip(axes.flat,['Coefficient error %','Water amplitude error %','Full-FID relative residual',
                                 'FSL minus nominal voxel mM','FSL minus source-response reference mM','Projection residual in fit window']):
        ax.set_xticks(x,CASES,rotation=25,ha='right');ax.set_title(title);ax.grid(alpha=.3)
    axes[0,0].legend();axes[0,2].legend();axes[1,2].legend();fig.tight_layout()
    fig.savefig(out/'bias_and_residuals.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(4,3,figsize=(15,12))
    for row,name in enumerate(CASES):
        a=np.load(out/name/'diagnostic_arrays.npz');d=result['cases'][name]
        for col,data,model,limits,title in [(0,'data','fitted',(4.2,.2),'Metabolite residual'),
                                           (1,'water','water_fitted',(5.3,4.1),'Water residual')]:
            residual=spectrum(a[data]-a[model]);ax=axes[row,col]
            ax.plot(a['ppm'],residual.real,label='real');ax.plot(a['ppm'],residual.imag,label='imag')
            ax.set_xlim(*limits);ax.set_title(name+': '+title);ax.legend(fontsize=7)
        ax=axes[row,2];time=a['time_s'];initial=abs(a['water'][0])
        demod=np.exp(-2j*np.pi*d['true_water_centroid_hz']*time)/initial
        for key,style in [('water','-'),('water_fitted','--')]:
            value=a[key]*demod
            ax.plot(time,value.real,style,label=key+' real');ax.plot(time,value.imag,style,label=key+' imag',alpha=.6)
        ax.set_xlim(0,.4);ax.set_title('Known-centroid frame; display / true initial amplitude')
        ax.set_xlabel('Time (s)');ax.legend(fontsize=6)
    fig.tight_layout();fig.savefig(out/'residuals_and_water_envelopes.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output_dir')
    result=run_validation(parser.parse_args().output_dir)
    print(json.dumps({c:dict(residual=d['relative_complex_residual'],water_residual=d['water_fit']['relative_complex_residual']) for c,d in result['cases'].items()},indent=2))

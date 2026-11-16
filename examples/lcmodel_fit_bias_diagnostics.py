"""Phase 7b.1: immutable-input LCModel controls and independent basis projections."""
import argparse
import os
import csv
import json
from pathlib import Path
import re
import shlex
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from fsl_mrs.utils import mrs_io

from examples.lcmodel_initial_validation import BINARY,BUNDLE,SOURCE,digest,save
from utils.lcmodel_quantitative_io import NAMES,independent_water_scaling,local_peak_area,parse_basis,parse_coord,parse_raw,parse_table

BASE=Path('outputs/lcmodel_validation/noise_free')
OUTPUT=Path('outputs/lcmodel_fit_bias_diagnostics')
# Public controls: manual §§5.3.8,8.6.4.1,9.7,11.5,11.12,11.13.
CASES={
    'replay':('', 'Unchanged fitting settings; only output paths differ'),
    'preliminary_matched':(" NUSE1=5, CHUSE1(1)='NAA', CHUSE1(2)='Cr', CHUSE1(3)='GPC', CHUSE1(4)='Glu', CHUSE1(5)='PCh'\n",'Use exactly the simulated five components in preliminary analysis'),
    'zero_broadening_prior':(' DEEXT2=0.\n','Expected additional 1/T2=0 rather than default2; SD prior unchanged'),
    'half_broadening_prior':(' DEEXT2=1.\n','Half the default positive expected additional 1/T2; SD prior unchanged'),
    'wider_convolution':(' RFWHM=3.\n','Manual §11.12 documented wider convolution support'),
    'phantom_mode':(' VITRO=T\n','Bundled phantom statistical/baseline mode; not an isolated baseline parameter'),
    'known_phase_prior':(' DEGZER=0., DEGPPM=0., SDDEGZ=3., SDDEGP=1.\n','Known near-zero phase with recommended finite prior uncertainties'),
    'default_basis_width':(' FWHMBA=.013\n','Manual default effective basis width assumption; actual FID is unchanged'),
}


def metrics(directory):
    fits=directory/'lcmodel_outputs';table=parse_table(fits/'fit.TABLE')
    execution=json.loads((directory/'execution.json').read_text())
    if execution['returncode'] or table['fatal'] or set(table['concentrations'])!=set(NAMES):
        return dict(table,execution=execution,optimizer_status='failed; retained diagnostic, no valid final fit'),None
    arrays=parse_coord(fits/'fit.COORD');printed=(fits/'fit.PRINT').read_text()
    def number(pattern):return float(re.findall(pattern,printed)[-1])
    result=dict(table,execution=execution,
        relative_real_residual=float(np.linalg.norm(arrays['data']-arrays['fit'])/np.linalg.norm(arrays['data'])),
        baseline_relative_norm=float(np.linalg.norm(arrays['baseline'])/np.linalg.norm(arrays['data'])),
        baseline_signed_area=float(np.trapezoid(arrays['baseline'][::-1],arrays['ppm'][::-1])),
        normalized_basis_area=number(r'Normalized area of reference Basis singlet =\s*([\d.E+-]+)'),
        water_area=number(r'Area of unsuppressed water peak =\s*([\d.E+-]+)'),
        FCALIB=number(r'FCALIB =\s*([\d.E+-]+)'),
        FWHM_ppm=number(r'FWHM =\s*([\d.E+-]+)'),
        alpha_baseline_and_lineshape=[float(x) for x in re.findall(r'alphaB,S =\s*([\d.Ee+-]+),\s*([\d.Ee+-]+)',printed)[-1]],
        final_parameters={})
    for name,shift,shift_sd,rate,rate_sd,concentration,sd in re.findall(
        r'\d+ met=\s*(\w+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([\d.E+-]+)\s+(\d+)',printed):
        result['final_parameters'][name]=dict(shift_millippm=float(shift),delta_inverse_T2_s=float(rate),rate_prior_deviation_SDs=float(rate_sd))
    line_shapes=re.findall(r' Lineshape =([^\n]+)',printed)
    result['last_printed_lineshape_coefficients']=[float(x) for x in line_shapes[-1].split()] if line_shapes else []
    shape=re.findall(r'(\d+) spline knots\.\s+Ns =\s*(\d+)\((\d+)\)',printed)[-1]
    result['baseline_knots'],result['convolution_sides'],result['convolution_stride']=map(int,shape)
    result['optimizer_status']='completed finite final analysis; no generic optimizer-success flag in TABLE'
    return result,arrays


def converted_basis():
    """MYBASI positive/negative-half bandwidth crop and unitary inverse FFT."""
    entries,meta=parse_basis(BASE/'lcmodel_compatible.BASIS')
    dt=.0005;n=4096
    columns=[]
    for name in NAMES:
        _,header,spectrum=next(e for e in entries if e[0]==name)
        fields={key:float(re.search(rf'\b{key}\s*=\s*([\d.Ee+-]+)',header)[1]) for key in ['CONC','TRAMP','VOLUME','ISHIFT']}
        step=1/(meta['points']*meta['dwell_time_s']*meta['frequency_mhz'])
        shift=int(fields['ISHIFT'])+int(np.rint((4.65-4.68)/step))
        spectrum=np.roll(spectrum,-shift)*fields['TRAMP']/fields['VOLUME']/fields['CONC']
        crop=np.r_[spectrum[:n//2],spectrum[-n//2:]]
        columns.append(np.fft.ifft(crop,norm='ortho'))
    return np.array(columns).T,meta


def project(design,data,*,real=False):
    if real:design=design.real;data=data.real
    coeff=np.linalg.lstsq(design,data,rcond=None)[0]
    return dict(coefficients_real=coeff.real.tolist(),coefficients_imag=coeff.imag.tolist() if np.iscomplexobj(coeff) else [0.]*len(coeff),
        relative_residual=float(np.linalg.norm(design@coeff-data)/np.linalg.norm(data)))


def water_recurrence_spectrum(fid):
    """Independent CSFT_R single-precision recurrence; known zero-phase control."""
    values=np.asarray(fid,dtype=np.complex64);n=len(values)
    angle=-np.float32(6.28318530717959)*np.arange(n,dtype=np.float32)/np.float32(n)
    factor=np.exp(np.complex64(1j)*angle)
    term=np.ones(n,dtype=np.complex64);spectrum=np.zeros(n,dtype=np.complex64)
    for value in values:
        spectrum+=value*term
        term*=factor
    return np.fft.fftshift(spectrum*np.float32(1/np.sqrt(np.float32(n))))*np.sqrt(np.float32(.5))


def normalization_for_width(width_ppm):
    """AREABA window effect separated from unchanged AREAW2 water area."""
    entries,meta=parse_basis(BASE/'lcmodel_compatible.BASIS')
    _,header,basis=next(e for e in entries if e[0]=='NAA')
    step=1/(meta['points']*meta['dwell_time_s']*meta['frequency_mhz'])
    center=round((4.68-2.01)/step);search=round(20*.004/step)
    center=center-search+int(np.argmax(basis.real[center-search:center+search+1]))
    half=round(.5*10*width_ppm/step)
    area=local_peak_area(basis.real,center,half,round(.1/step),step)/3
    water,_=parse_raw(BASE/'water_reference.H2O');ws=water_recurrence_spectrum(water)
    water_step=1/(len(water)*.0005*meta['frequency_mhz'])
    ppm=4.68-np.fft.fftshift(np.fft.fftfreq(len(water),.0005))/meta['frequency_mhz']
    peak=int(np.argmax(np.where(abs(ppm-4.7)<=1,abs(ws),0)))
    water_area=local_peak_area(ws.real,peak,round(2/water_step),round(.2/(water_step/2)),water_step)
    return dict(basis_CONC=1.,basis_TRAMP=1.,basis_VOLUME=1.,N1HMET=3,
        WCONC_mM=39016.5,ATTH2O=1.,ATTMET=1.,FWHMBA_ppm=width_ppm,
        normalized_basis_area=area,water_area=water_area,
        FCALIB=area*2*39016.5/water_area,
        scope='independent AREABA and single-precision CSFT_R/AREAW2 zero-phase calculation; no fitted gain')


def independent_analysis():
    basis,meta=converted_basis()
    data,_=parse_raw(BASE/'metabolites.RAW');water,_=parse_raw(BASE/'water_reference.H2O')
    source_basis=mrs_io.read_basis(str(BUNDLE/'metabolite_basis.BASIS'))
    source=source_basis.get_formatted_basis(2000,2048)
    source=np.array([source[:,source_basis.names.index(n)].conj() for n in NAMES]).T
    ground_truth=json.loads((SOURCE/'ground_truth.json').read_text())
    truth=np.array([ground_truth['metabolite_concentrations_mM'][n] for n in NAMES])
    # Derived from the native FFT length and LCModel unitary bandwidth transform.
    fft_factor=meta['points']/np.sqrt(len(basis));expected=truth*.001/fft_factor
    fixed=basis[:len(data)]
    normalized=independent_water_scaling(BUNDLE/'metabolite_basis.BASIS',water)
    original=metrics(BASE)[0];cal=original['FCALIB']
    zero=np.pad(data,(0,len(basis)-len(data)))
    frequencies=np.fft.fftfreq(len(basis),.0005);ppm=4.68-frequencies/meta['frequency_mhz']
    mask=(ppm>=.2)&(ppm<=4.2)
    spectral=np.fft.fft(basis,axis=0,norm='ortho')[mask]
    data_spectrum=np.fft.fft(zero,norm='ortho')[mask]
    truncated=basis.copy();truncated[len(data):]=0
    truncated_spectrum=np.fft.fft(truncated,axis=0,norm='ortho')[mask]
    projections={
        'matching_FSL_complex_acquired':project(source,data),
        'LCModel_complex_acquired':project(fixed,data),
        'LCModel_complex_excluding_first':project(fixed[1:],data[1:]),
        'LCModel_full_spectral_complex':project(spectral,data_spectrum),
        'LCModel_full_spectral_real':project(spectral,data_spectrum,real=True),
        'LCModel_truncated_spectral_real':project(truncated_spectrum,data_spectrum,real=True),
    }
    x=ppm[mask];x=(x-x.mean())/np.ptp(x)
    polynomial=np.column_stack([x**order for order in range(4)])
    projections['LCModel_real_with_cubic_baseline']=project(np.column_stack([spectral.real,polynomial]),data_spectrum.real,real=True)
    for name,p in projections.items():
        scale=1000 if name.startswith('matching') else fft_factor*1000
        p['source_equivalent_mM']=(np.array(p['coefficients_real'][:5])*scale).tolist()
        if not name.startswith('matching'):
            p['normalization_scaled_mM']=(np.array(p['coefficients_real'][:5])*cal).tolist()
    differences={}
    for i,name in enumerate(NAMES):
        a=fixed[:,i]/fft_factor;b=source[:,i]
        differences[name]=dict(relative_acquired_FID_error=float(np.linalg.norm(a-b)/np.linalg.norm(b)),
            relative_FID_error_excluding_first=float(np.linalg.norm(a[1:]-b[1:])/np.linalg.norm(b[1:])),
            max_absolute_FID_error=float(abs(a-b).max()),
            first_complex_difference=[float((a[0]-b[0]).real),float((a[0]-b[0]).imag)],
            LCModel_unacquired_tail_relative_norm=float(np.linalg.norm(basis[2048:,i])/np.linalg.norm(basis[:2048,i])))
    output=dict(names=list(NAMES),truth_mM=truth.tolist(),expected_unscaled_LCModel_coefficients=expected.tolist(),
        unitary_FFT_factor=float(fft_factor),independent_normalization=normalized,printed_normalization=cal,
        normalization_only_mM=(expected*cal).tolist(),basis_data_differences=differences,projections=projections,
        normalization_checks={key:normalization_for_width(width) for key,width in [('baseline',1/127.731),('manual_default',.013)]},
        provenance={str(p):digest(p) for p in [Path(os.environ['LCMODEL_MANUAL']),Path(os.environ['LCMODEL_SOURCE']),Path('pixi.lock')]},
        scope='Fixed frequency, phase and linewidth; polynomial projection is diagnostic, not LCModel spline/convolution model')
    return output


def run(output=OUTPUT):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    protected=list(BASE.rglob('*'))+list(BUNDLE.iterdir())+list(SOURCE.iterdir())+[BINARY]
    protected=[p for p in protected if p.is_file()]
    hashes={str(p.resolve()):digest(p) for p in protected}
    manifest=output/'immutable_baseline_sha256.json'
    if manifest.exists() and json.loads(manifest.read_text())!=hashes:raise ValueError('immutable baseline/source changed')
    if not manifest.exists():save(manifest,hashes)
    baseline=metrics(BASE)[0]
    save(output/'baseline_metrics.json',baseline)
    analysis=independent_analysis();save(output/'independent_analysis.json',analysis)
    rows=[];commands=[];all_results={};all_arrays={}
    for name,(options,description) in CASES.items():
        directory=output/name;directory.mkdir(exist_ok=True);fits=directory/'lcmodel_outputs';fits.mkdir(exist_ok=True)
        control=(BASE/'lcmodel.CONTROL').read_text().replace(str((BASE/'lcmodel_outputs').resolve()),str(fits))
        control=control.rsplit(' $END',1)[0]+options+' $END\n'
        path=directory/'lcmodel.CONTROL'
        if path.exists() and path.read_text()!=control:raise ValueError(f'cached control differs: {name}')
        if not path.exists():path.write_text(control)
        command=f'cd {shlex.quote(str(directory))} && {shlex.quote(str(BINARY))} < lcmodel.CONTROL'
        commands.append(command)
        if not (directory/'execution.json').exists():
            p=subprocess.run([str(BINARY)],input=control,cwd=directory,text=True,capture_output=True,timeout=300)
            save(directory/'execution.json',dict(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr))
        result,arrays=metrics(directory)
        for device,extension in [('pdfwrite','pdf')]:
            target=fits/f'fit.{extension}'
            if not target.exists():subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE',f'-sDEVICE={device}',f'-sOutputFile={target}',str(fits/'fit.PS')],check=True,capture_output=True)
        result['options']=options;result['description']=description
        save(directory/'metrics.json',result);all_results[name]=result
        if arrays is None:
            print(json.dumps(dict(case=name,fatal=result['fatal'],diagnostics=result['diagnostics'])),flush=True)
            continue
        all_arrays[name]=arrays
        row=dict(case=name,**{n:result['concentrations'][n]['concentration'] for n in NAMES},
            FCALIB=result['FCALIB'],water_area=result['water_area'],basis_area=result['normalized_basis_area'],
            residual=result['relative_real_residual'],baseline_norm=result['baseline_relative_norm'],
            FWHM_ppm=result['FWHM_ppm'],baseline_alpha=result['alpha_baseline_and_lineshape'][0],
            lineshape_alpha=result['alpha_baseline_and_lineshape'][1],knots=result['baseline_knots'],
            diagnostics=';'.join(d['routine']+str(d['code']) for d in result['diagnostics']))
        rows.append(row);print(json.dumps(row),flush=True)
    (output/'commands.txt').write_text('\n'.join(commands)+'\n')
    with (output/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    figure,axes=plt.subplots(len(all_arrays),1,figsize=(12,3*len(all_arrays)),sharex=True)
    for ax,(name,arr) in zip(axes,all_arrays.items()):
        ax.plot(arr['ppm'],arr['data'],label='Input',lw=1);ax.plot(arr['ppm'],arr['fit'],label='LCModel fit',lw=.8)
        ax.plot(arr['ppm'],arr['data']-arr['fit'],label='Residual',lw=.8);ax.plot(arr['ppm'],arr['baseline'],label='Baseline',lw=.8)
        ax.set_title(name);ax.legend(ncol=4,fontsize=8)
    axes[-1].set_xlim(4.2,.2);axes[-1].set_xlabel('ppm');figure.tight_layout();figure.savefig(output/'spectral_comparisons.png',dpi=130);plt.close(figure)
    fig,axes=plt.subplots(1,2,figsize=(14,4))
    axes[0].bar([r['case'] for r in rows],[r['NAA']-10 for r in rows]);axes[0].tick_params(axis='x',rotation=65);axes[0].set_ylabel('NAA error (mM)')
    for row in rows:
        axes[1].scatter(row['residual']*100,row['NAA']-10,label=row['case'])
    axes[1].legend(loc='upper left',bbox_to_anchor=(1.02,1),fontsize=7)
    axes[1].set_xlabel('Phased-real residual (%)');axes[1].set_ylabel('NAA error (mM)');fig.tight_layout();fig.savefig(output/'concentration_vs_residual.png',dpi=150);plt.close(fig)
    if hashes!={str(p.resolve()):digest(p) for p in protected}:raise RuntimeError('immutable input changed')
    save(output/'results.json',dict(baseline=baseline,analysis=analysis,cases=all_results,immutable_files_unchanged=True))
    return all_results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=OUTPUT)
    run(parser.parse_args().output)

"""Measure saved Phase 5j peaks; no regeneration, fitting or calibration changes."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import fsl_mrs
import numpy as np
from fsl_mrs.utils import mrs_io
from fsl_mrs.utils.misc import ts_to_ts_ft

from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral

SOURCE=Path('outputs/mrscloud_sampling_sensitivity/validated')
METHYL_FREQUENCY=(2.008-3)*127.728+(3-4.68)*127.731
ZERO_FILL=128


def absorption_spectrum(fid, *, dwell_time_s, center_hz, phase_rad, zero_fill=ZERO_FILL):
    """Demodulate known peak reference, half-weight sample zero, then zero fill."""
    fid=np.asarray(fid)
    if fid.ndim!=1 or len(fid)<4 or not np.iscomplexobj(fid) or not np.all(np.isfinite(fid)):
        raise ValueError('requires a finite one-dimensional complex FID')
    if dwell_time_s<=0 or not isinstance(zero_fill,int) or zero_fill<1:
        raise ValueError('positive dwell and integer zero-filling factor required')
    time=np.arange(len(fid))*dwell_time_s
    weighted=fid*np.exp(-1j*(2*np.pi*center_hz*time+phase_rad))
    weighted[0]*=.5
    size=len(fid)*zero_fill
    freq=np.fft.fftshift(np.fft.fftfreq(size,dwell_time_s))
    spec=np.fft.fftshift(np.fft.fft(weighted,n=size))/np.sqrt(len(fid))
    return freq,spec


def fwhm(frequency, absorption, *, search_hz=20., baseline=0.):
    """Local central-lobe half-height crossings with linear interpolation."""
    selected=(frequency>=-search_hz)&(frequency<=search_hz)
    x=np.asarray(frequency)[selected];y=np.asarray(absorption)[selected]-baseline
    peak=int(np.argmax(y))
    if len(x)<5 or peak==0 or peak==len(x)-1 or y[peak]<=0:
        raise ValueError('positive interior target peak required')
    half=y[peak]/2
    left=np.where(y[:peak]<=half)[0]
    right=np.where(y[peak+1:]<=half)[0]
    if not len(left) or not len(right):raise ValueError('half-height crossings not bracketed')
    lo=int(left[-1]);hi=int(peak+1+right[0])
    xl=x[lo]+(half-y[lo])*(x[lo+1]-x[lo])/(y[lo+1]-y[lo])
    xr=x[hi-1]+(half-y[hi-1])*(x[hi]-x[hi-1])/(y[hi]-y[hi-1])
    return {'fwhm_hz':float(xr-xl),'left_half_hz':float(xl),'right_half_hz':float(xr),
            'peak_offset_hz':float(x[peak]),'half_height':float(half+baseline),
            'baseline':float(baseline)}


def measure(fid, center, phase, zero_fill=ZERO_FILL):
    freq,spec=absorption_spectrum(fid,dwell_time_s=.0005,center_hz=center,
                                   phase_rad=phase,zero_fill=zero_fill)
    result=fwhm(freq,spec.real)
    return result,freq,spec.real


def run(output_dir, source=SOURCE):
    source,out=Path(source),Path(output_dir)
    out.mkdir(parents=True,exist_ok=False)
    files=sorted(p for p in source.rglob('*') if p.is_file())
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    metadata=json.loads((source/'results.json').read_text())
    master=np.load(source/'master_responses.npz')
    phase_m=float(np.angle(master['methyl_raw'][0]));phase_w=float(np.angle(master['water_raw'][0]))
    dc_kernels={}
    for name in ['methyl','full']:
        impulse=np.zeros(16384,complex)
        impulse[0]=master[name+'_DC'][0]-master[name+'_raw'][0]
        dc_kernels[name]=ts_to_ts_ft(impulse[:,None],.00025,.0005,8192).ravel()
    rows=[];details=[];plot_cases={}
    for case in metadata['cases']:
        original=case['row'];n=original['points'];width=original['total_linewidth_hz']
        directory=source/f'N{n}_LW{width}'
        config=json.loads((directory/'simulation_config.json').read_text())
        if config['dwell_time_s']!=.0005 or config['frequency_origin_ppm']!=4.68:
            raise ValueError('unexpected saved spectral grid/reference')
        t=np.arange(n)*.0005;envelope=np.exp(-np.pi*(width-1)*t)
        generated={name:master[name][:n]*envelope for name in
                   ['water_raw','methyl_raw','methyl_DC','methyl_formatted','full_raw','full_DC','full_formatted']}
        params=case['fit_parameters'];wp=case['water_fit']['parameters']
        b=mrs_io.read_basis(str(directory/'matching.BASIS'))
        basis=b.get_formatted_basis(2000,n)[:,b.names.index('NAA')]
        # FSL's phased reference model removes fitted Phi0/Phi1; intrinsic
        # basis phase is retained. No saved optimizer parameter is changed.
        generated['fitted_NAA']=basis*np.exp(-(params['gamma_0']+1j*params['eps_0'])*t)
        generated['fitted_water']=np.exp(-wp[1]*t-wp[2]*t*t).astype(complex)
        row={'points':n,'nominal_linewidth_hz':width,'additional_linewidth_hz':width-1}
        measured={};spectra={}
        for name,fid in generated.items():
            if name in ['water_raw','fitted_water']:
                center=2.55462 if name=='water_raw' else 0.
                phase=phase_w if name=='water_raw' else 0.
            else:
                center=METHYL_FREQUENCY-(params['eps_0']/(2*np.pi) if name=='fitted_NAA' else 0.)
                phase=phase_m
            m,freq,spec=measure(fid,center,phase)
            measured[name]=m
            local=abs(freq)<25
            spectra[name]=(freq[local],spec[local])
            row[name+'_fwhm_hz']=m['fwhm_hz']
            if name=='methyl_DC':
                delta=master['methyl_DC'][0]-master['methyl_raw'][0]
                pedestal=float((.5*delta*np.exp(-1j*phase_m)/np.sqrt(n)).real)
                corrected=fwhm(freq,spec,baseline=pedestal)
                row['methyl_DC_background_removed_fwhm_hz']=corrected['fwhm_hz']
                measured['methyl_DC_background_removed']=corrected
            if name in ['methyl_formatted','full_formatted']:
                component='methyl' if name=='methyl_formatted' else 'full'
                _,pedestal=absorption_spectrum(dc_kernels[component][:n]*envelope,
                    dwell_time_s=.0005,center_hz=center,phase_rad=phase)
                corrected=fwhm(freq,spec-pedestal.real)
                row[name+'_DC_background_removed_fwhm_hz']=corrected['fwhm_hz']
                measured[name+'_DC_background_removed']=corrected
                spectra[name+'_DC_background_removed']=(freq[local],(spec-pedestal.real)[local])
        row['raw_methyl_minus_water_hz']=row['methyl_raw_fwhm_hz']-row['water_raw_fwhm_hz']
        row['processed_methyl_minus_water_hz']=row['methyl_formatted_fwhm_hz']-row['water_raw_fwhm_hz']
        row['full_NAA_minus_water_hz']=row['full_formatted_fwhm_hz']-row['water_raw_fwhm_hz']
        # Check that half-height interpolation is converged beyond the smallest
        # processed/raw width difference, without altering acquired samples.
        finer,_,_=measure(generated['methyl_formatted'],METHYL_FREQUENCY,phase_m,zero_fill=2*ZERO_FILL)
        row['methyl_zero_fill_doubling_difference_hz']=finer['fwhm_hz']-row['methyl_formatted_fwhm_hz']
        row['fitted_water_Lorentzian_parameter_hz']=wp[1]/np.pi
        row['fitted_water_sigma_per_s2']=wp[2]
        row['fitted_NAA_additional_Lorentzian_parameter_hz']=params['gamma_0']/np.pi
        row['fitted_NAA_nominal_plus_additional_hz']=1+params['gamma_0']/np.pi
        # Peak widths alone do not specify area. A fixed-amplitude, hypothetical
        # Lorentzian-width perturbation isolates their possible area effect;
        # it is not applied to any saved signal or quantification result.
        ppm=np.fft.fftshift(np.fft.fftfreq(n,.0005))/127.731+4.68
        amp=abs(master['methyl_raw'][0])
        simple=amp*np.exp(-np.pi*width*t)*np.exp(2j*np.pi*METHYL_FREQUENCY*t)
        jw=spectral_integral(abs(master['water_raw'][0])*np.exp(-np.pi*width*t),ppm,(1.65,7.65))
        base=10*(2/3)*spectral_integral(simple,ppm,(1.8,2.2))/jw
        for kind,delta in [('raw',row['raw_methyl_minus_water_hz']),
                           ('processed',row['processed_methyl_minus_water_hz'])]:
            control=simple*np.exp(-np.pi*delta*t)
            row[kind+'_width_only_surrogate_delta_mM']=10*(2/3)*spectral_integral(control,ppm,(1.8,2.2))/jw-base
        row['original_predicted_mM']=original['predicted_mM']
        row['original_reported_mM']=original['reported_mM']
        row['original_ideal_methyl_factor']=case['prediction']['control_factors']['ideal_methyl']
        details.append({'row':row,'measurements':measured})
        rows.append(row)
        if (n,width) in [(2048,1),(2048,8),(8192,1)]:plot_cases[f'N{n}_LW{width}']=spectra
    if hashes!={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}:
        raise RuntimeError('existing Phase 5j dataset changed')
    result={'FSL_MRS_version':fsl_mrs.__version__,'source':str(source),'source_hashes':hashes,
        'source_hashes_unchanged':True,'zero_fill_factor':ZERO_FILL,'interpolation':'linear half-height crossings',
        'measurement_background':'zero absorption level; DC/overlap contributions retained and diagnosed separately',
        'methyl_reference_frequency_hz':METHYL_FREQUENCY,'methyl_phase_rad':phase_m,
        'water_phase_rad':phase_w,'fit_reexecuted':False,'cases':details}
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    with (out/'linewidths.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    plot(out,plot_cases)
    print(json.dumps(rows,indent=2))
    return result


def plot(out,cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for label,spectra in cases.items():
        fig,axes=plt.subplots(1,3,figsize=(13,4))
        groups=[['water_raw','methyl_raw'],['methyl_raw','methyl_DC','methyl_formatted','methyl_formatted_DC_background_removed'],
                ['water_raw','fitted_water','full_formatted','fitted_NAA']]
        for ax,names in zip(axes,groups):
            for name in names:
                freq,spec=spectra[name];sel=abs(freq)<15
                ax.plot(freq[sel],spec[sel]/max(spec[sel]),label=name)
            ax.axhline(.5,color='k',ls=':',lw=.7)
            ax.set(xlabel='Frequency relative to target resonance (Hz)',ylabel='Display-normalized absorption')
            ax.legend(fontsize=7)
        fig.tight_layout();fig.savefig(out/f'{label}_absorption.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output_dir')
    args=parser.parse_args();run(args.output_dir)

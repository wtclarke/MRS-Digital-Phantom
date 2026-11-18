"""Aggregate measured acquisition properties separately from requested settings."""
import json
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from fsl_mrs.utils import mrs_io


def linewidth(fid,dwell,frequency):
    """Magnitude-peak half-height width, interpolated on an eightfold padded grid."""
    spectrum=np.abs(np.fft.fftshift(np.fft.fft(fid,n=8*len(fid))))
    hz=np.fft.fftshift(np.fft.fftfreq(len(spectrum),dwell))
    select=np.flatnonzero(np.abs(hz)<frequency)
    peak=select[np.argmax(spectrum[select])];half=spectrum[peak]/2
    left=np.flatnonzero(spectrum[:peak]<half);right=np.flatnonzero(spectrum[peak+1:]<half)
    if not len(left) or not len(right):return float('nan')
    a=left[-1];b=peak+1+right[0]
    lo=np.interp(half,spectrum[a:a+2],hz[a:a+2])
    hi=np.interp(half,spectrum[b-1:b+1][::-1],hz[b-1:b+1][::-1])
    return float(hi-lo)


def audit(root):
    rows=[];normalizers={};noise_hashes={};fit_status=[]
    for case in pd.read_csv(root/'cases.csv').to_dict('records'):
        source=case.get('source')
        if not isinstance(source,str) or not source:continue
        path=Path(source);directory=path.parent.parent
        row={k:case.get(k) for k in ['case_id','stage','basis_name','tissue','snr','b0_sd','water_snr']}
        for name in ['applied_field','achieved_snr','water_check']:
            p=directory/f'{name}.json'
            if p.exists():row.update({f'{name}_{k}':v for k,v in json.loads(p.read_text()).items()})
        ref=mrs_io.read_FID(str(path/'water_reference.nii.gz'))
        row['water_magnitude_fwhm_Hz']=linewidth(ref.mrs().FID,ref.dwelltime,ref.spectrometer_frequency[0])
        truth=json.loads((path/'ground_truth.json').read_text())
        table=pd.read_csv(directory/'concentrations.csv')
        fractions=truth['tissue_fractions']
        errors=[]
        for name,value in truth['voxel_weighted'].items():
            profile=table.loc[table.Metabolite.eq(name)]
            expected=sum(fractions.get(r.Tissue,0)*r.Conc_mean for r in profile.itertuples())
            errors.append(abs(value-expected))
        row['concentration_truth_max_abs_error_mM']=max(errors,default=0.)
        provenance=json.loads((path/'calibration_provenance.json').read_text())
        normalizers.setdefault(case['basis_name'],set()).add(provenance['normalization']['scalar'])
        expected_water=sum(fractions.get(t,0)*v for t,v in provenance['water_options']['water_concentrations_mol_l'].items())
        row['water_truth_abs_error_mol_l']=abs(expected_water-truth['water_reference']['nominal_water_concentration_mol_l'])
        sigma=provenance['water_options'].get('noise_sigma',0.)
        row['native_water_FID_sigma']=sigma
        clean_path=path/'water_reference_clean.nii.gz'
        if sigma and clean_path.exists():
            clean=mrs_io.read_FID(str(clean_path)).mrs().FID
            noise=np.fft.fft(ref.mrs().FID-clean)
            measured=np.sqrt((np.var(noise.real)+np.var(noise.imag))/2)
            peak=float(np.max(np.abs(np.fft.fft(clean))))
            row['achieved_snr_water']=peak/measured
            row['water_snr_from_applied_sigma']=peak/(sigma*np.sqrt(len(clean)))
            row['water_noise_measured_spectral_sigma']=float(measured)
            row['water_snr_note']='Derived from native applied sigma and clean saved reference; supersedes generation-time achieved_snr.json water value'
        total=mrs_io.read_FID(str(path/'total.nii.gz')).mrs().FID
        metabs=mrs_io.read_FID(str(path/'metabs.nii.gz')).mrs().FID
        if np.isfinite(case.get('snr',float('inf'))) and case['stage'] in ['noise','combined','water_noise']:
            noise=np.fft.fft(total-metabs)
            measured=np.sqrt((np.var(noise.real)+np.var(noise.imag))/2)
            row['metabolite_noise_measured_spectral_sigma']=float(measured)
            row['achieved_snr_metabolite']=float(np.max(np.abs(np.fft.fft(metabs)))/measured)

        if case['stage']=='noise':
            key=tuple(case[k] for k in ['basis_name','tissue','target','scale','snr','b0_sd'])
            entry=noise_hashes.setdefault(key,{'clean':set(),'noise':set(),'count':0})
            entry['clean'].add(hashlib.sha256(np.asarray(metabs).tobytes()).hexdigest())
            entry['noise'].add(hashlib.sha256(np.asarray(total-metabs).tobytes()).hexdigest())
            entry['count']+=1
        for fitter in ['fsl','lcmodel']:
            result_path=directory/fitter/'study_result.json'
            if result_path.exists():
                result=json.loads(result_path.read_text())
                fit_status.append(dict(case_id=case['case_id'],stage=case['stage'],basis_name=case['basis_name'],fitter=fitter,
                    status=result['status'],error=result.get('error',''),tissue=case['tissue'],snr=case['snr'],b0_sd=case['b0_sd']))
        rows.append(row)
    pd.DataFrame(rows).to_csv(root/'acquisition_checks.csv',index=False)

    data=pd.DataFrame(rows)
    summary=dict(cases=len(data),concentration_truth_max_abs_error_mM=float(data.concentration_truth_max_abs_error_mM.max()),
        water_truth_max_abs_error_mol_l=float(data.water_truth_abs_error_mol_l.max()),
        fixed_normalization_per_basis=all(len(v)==1 for v in normalizers.values()),
        controlled_water_checks=int(data.water_check_relative_complex_error.notna().sum()) if 'water_check_relative_complex_error' in data else 0,
        max_water_relative_complex_error=float(data.water_check_relative_complex_error.max()) if 'water_check_relative_complex_error' in data else None,
        noise_groups=len(noise_hashes),noise_groups_with_identical_clean_signals=sum(len(v['clean'])==1 for v in noise_hashes.values()),
        noise_groups_with_distinct_draws=sum(len(v['noise'])==v['count'] for v in noise_hashes.values()))
    (root/'acquisition_audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')

    pd.DataFrame(fit_status).to_csv(root/'fitter_status.csv',index=False)

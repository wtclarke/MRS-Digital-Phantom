"""Read-only GE/external BASIS numerical audit; no fitting or calibration changes."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
from fsl_mrs.utils import mrs_io

from examples.basis_h2o_headless import configuration
from examples.quantitative_water_reference_production import relaxation_config
from simulation.basis_h2o import load_responses,simulate_basis_h2o
from simulation.quantitative import export_quantitative,DEFAULT_WATER


def compare_source_mat(basis_path,source_mat):
    """One declared acquisition-frame shift, no amplitude alignment/fitting."""
    from scipy.io import loadmat
    mat=loadmat(source_mat,simplify_cells=True)['BASIS']
    basis=mrs_io.read_basis(str(basis_path))
    columns=np.column_stack([mat['fids'][:,list(mat['name']).index(n)] for n in basis.names])
    time=np.arange(basis.original_points)*basis.original_dwell
    shift=(float(mat['centerFreq'])-4.68)*basis.cf
    transformed=columns*np.exp(2j*np.pi*shift*time)[:,None]
    errors=np.linalg.norm(transformed-basis.original_basis_array,axis=0)/np.linalg.norm(transformed,axis=0)
    return {'source_mat_sha256':hashlib.sha256(source_mat.read_bytes()).hexdigest(),
        'source_global_scale':float(mat['scale']),'rf_carrier_ppm':float(mat['centerFreq']),
        'acquisition_origin_ppm':4.68,'frame_shift_hz':shift,'amplitude_factor_applied':1.,
        'relative_errors':dict(zip(map(str,basis.names),map(float,errors))),
        'max_relative_error':float(errors.max()),
        'normalized_source_max_real_spectrum':float(np.max(mat['specs'].real)),
        'individual_max_real_spectrum_range':[float(np.max(mat['specs'].real,axis=0).min()),float(np.max(mat['specs'].real,axis=0).max())]}


def validate(basis_path,output,source_mat=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    before=hashlib.sha256(basis_path.read_bytes()).hexdigest()
    b=mrs_io.read_basis(str(basis_path));units=b.get_formatted_basis(2000,2048,scale_factor=None)
    names=list(map(str,b.names));water=units[:,names.index('H2O')]
    cfg=configuration(basis_path,fractions=[1,0,0])['quantitative'];cfg['basis_path']=str(basis_path)
    cases={'pure_WM':dict(cfg,tissue_fractions={'WM':1}), 'pure_GM':dict(cfg,tissue_fractions={'GM':1}),
        'pure_CSF':dict(cfg,tissue_fractions={'CSF':1}), 'mixed':dict(cfg,tissue_fractions={'WM':.2,'GM':.3,'CSF':.5})}
    relax=relaxation_config();relax['TE_s']=.035
    cases['relaxation_B0']=dict(cfg,basis_relaxation_assumption='no_tissue_attenuation_encoded',relaxation=relax,
        additional_linewidth_hz={'metabolites':2.,'water':7.},spatial={'cell_labels':['WM','GM','CSF'],
        'spatial_weights':[.2,.3,.5],'b0_offsets_hz':[-3,7,1],'b0_units':'Hz','b0_convention':'positive_frequency'})
    cases['relaxation_B0'].pop('tissue_fractions')
    cases['noisy']=dict(cases['relaxation_B0'],noise={'seed':10010,'metabolite_snr':50,'water_snr':1000})
    records=[];t=np.arange(2048)/2000
    for name,q in cases.items():
        pair=simulate_basis_h2o(q);repeat=simulate_basis_h2o(q)
        em=np.zeros(2048,complex);ew=em.copy()
        labels=q.get('spatial',{}).get('cell_labels',['WM','GM','CSF'])
        weights=q.get('spatial',{}).get('spatial_weights',[q.get('tissue_fractions',{}).get(k,0) for k in labels])
        offsets=q.get('spatial',{}).get('b0_offsets_hz',[0]*len(labels))
        for tissue,f,df in zip(labels,weights,offsets):
            def r(p):return (1-np.exp(-2/p['T1_s']))*np.exp(-.035/p['T2_s'])
            for m,c in q['metabolite_concentrations_mM'][tissue].items():
                R=r(q['relaxation']['metabolites'][tissue][m]) if q.get('relaxation') else 1.
                em+=f*.001*c*R*units[:,names.index(m)]*np.exp(2j*np.pi*df*t)
            R=r(q['relaxation']['water'][tissue]) if q.get('relaxation') else 1.
            ew+=f*DEFAULT_WATER[tissue]*R*water*np.exp(2j*np.pi*df*t)
        line=q.get('additional_linewidth_hz',{'metabolites':0,'water':0})
        em*=np.exp(-np.pi*line['metabolites']*t);ew*=np.exp(-np.pi*line['water']*t)
        cm=pair.get('clean_metabolite_fid',pair['metabolite_fid']);cw=pair.get('clean_water_reference_fid',pair['water_reference_fid'])
        np.testing.assert_allclose(cm,em,rtol=1e-13,atol=1e-20)
        np.testing.assert_allclose(cw,ew,rtol=1e-13,atol=1e-17)
        paths=export_quantitative(pair,output/name);export_error=0.;phase_error=0.
        for label,key in [('metabolites','metabolite_fid'),('water_reference','water_reference_fid')]:
            np.testing.assert_array_equal(pair[key],repeat[key])
            x=mrs_io.read_FID(paths[label]).mrs().FID;np.testing.assert_array_equal(x,pair[key])
            export_error=max(export_error,float(abs(x-pair[key]).max()))
            mask=abs(x)>0
            if mask.any():phase_error=max(phase_error,float(abs(np.angle(x[mask]/pair[key][mask])).max()))
        records.append({'case':name,'independent_metabolite_error':float(abs(cm-em).max()),
            'independent_water_error':float(abs(cw-ew).max()),'export_max_complex_error':export_error,
            'export_max_phase_error_rad':phase_error,'water_initial':[float(cw[0].real),float(cw[0].imag)],
            'voxel_concentrations_mM':pair['ground_truth']['voxel_average_metabolite_concentrations_mM'],
            'water_concentration_mol_l':pair['ground_truth']['effective_voxel_water_concentration_mol_l']})
    q=cases['pure_WM'];a=simulate_basis_h2o(q);gain=simulate_basis_h2o(dict(q,receiver_factor=3))
    one=deepcopy(q)
    for p in one['metabolite_concentrations_mM'].values():
        for n in list(p):
            if n!='NAA':del p[n]
    selected=simulate_basis_h2o(one)
    assert np.array_equal(a['components']['NAA'],selected['components']['NAA'])
    assert np.array_equal(a['water_reference_fid'],selected['water_reference_fid'])
    ratio=a['components']['NAA'][0]/a['water_reference_fid'][0]
    expected=.01*units[0,names.index('NAA')]/(39.0165*water[0])
    np.testing.assert_allclose(ratio,expected,rtol=5e-16)
    gain_ratio=gain['components']['NAA'][0]/gain['water_reference_fid'][0]
    assert hashlib.sha256(basis_path.read_bytes()).hexdigest()==before
    native=b.original_basis_array[:,names.index('H2O')]
    result={'basis_sha256':before,'native_bandwidth_hz':b.original_bw,'native_points':b.original_points,
        'spectrometer_frequency_mhz':b.cf,'components':names,'normalization':a['simulation_config']['synthetic_calibration']['normalization'],
        'native_H2O_initial':[float(native[0].real),float(native[0].imag)],'native_H2O_phase_rad':float(np.angle(native[0])),
        'formatted_H2O_initial':[float(water[0].real),float(water[0].imag)],
        'formatted_NAA_initial':[float(units[0,names.index('NAA')].real),float(units[0,names.index('NAA')].imag)],
        'pure_WM_NAA_water_ratio':[float(ratio.real),float(ratio.imag)],'ratio_independent_error':float(abs(ratio-expected)),
        'receiver_ratio_error':float(abs(ratio-gain_ratio)),'selection_max_complex_error':float(abs(a['components']['NAA']-selected['components']['NAA']).max()),
        'cases':records,'source_unchanged':True,'fitting':'not performed'}
    if source_mat is not None:result['companion_mat_comparison']=compare_source_mat(basis_path,source_mat)
    (output/'validation_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--basis',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source-mat',type=Path,help='Optional original MRSCloud assembly for read-only provenance comparison')
    args=p.parse_args();validate(args.basis,args.output,args.source_mat)

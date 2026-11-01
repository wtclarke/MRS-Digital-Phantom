"""Generate four production headless pairs; no fitting or MRSCloud regeneration."""
import argparse
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
from fsl_mrs.utils import mrs_io

from simulation.headless import run_headless
from simulation.mrscloud_full_basis import load_library,sha256
from simulation.quantitative import MODE,DEFAULT_WATER
from simulation.calibration_resources import default_calibration_bundle

BUNDLE=default_calibration_bundle()
C={'WM':{'NAA':10.,'Cr':8.,'GPC':2.,'PCh':1.,'Glu':8.},
   'GM':{'NAA':14.,'Cr':9.,'GPC':1.5,'PCh':1.2,'Glu':12.},
   'CSF':{'NAA':0.,'Cr':0.,'GPC':0.,'PCh':0.,'Glu':0.}}


def relaxation_config():
    water={'WM':(.97,.073),'GM':(1.5,.088),'CSF':(4.47,2.03)}
    profiles={'WM':{'NAA':(1.29,.194),'Cr':(1.1,.15),'GPC':(1.2,.2),'PCh':(1.1,.18),'Glu':(1.15,.16)},
        'GM':{'NAA':(1.45,.18),'Cr':(1.35,.17),'GPC':(1.3,.22),'PCh':(1.25,.21),'Glu':(1.4,.19)},
        'CSF':{m:(1.5,.25) for m in C['CSF']}}
    return dict(mode='tissue_specific',TE_s=.03,TR_s=2.,
        water={t:dict(T1_s=a,T2_s=b) for t,(a,b) in water.items()},
        metabolites={t:{m:dict(T1_s=a,T2_s=b) for m,(a,b) in p.items()} for t,p in profiles.items()})


def validate(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    before={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}
    unit,water,_=load_library(BUNDLE,selected_names=list(C['WM']))
    base=dict(calibration_bundle=str(BUNDLE.resolve()),metabolite_concentrations_mM=C)
    configurations={
        'pure_WM':dict(base,tissue_fractions={'WM':1}),
        'mixed':dict(base,tissue_fractions={'WM':.2,'GM':.3,'CSF':.5}),
        'mixed_relaxation_B0':dict(base,relaxation=relaxation_config(),spatial={
            'cell_labels':['WM','GM','CSF'],'spatial_weights':[.2,.3,.5],
            'b0_offsets_hz':[-3.,7.,1.],'b0_units':'Hz','b0_convention':'positive_frequency'}),
        'noisy':dict(base,tissue_fractions={'WM':1},noise={'seed':70000,'metabolite_snr':50,'water_snr':1000})}
    records=[];t=np.arange(2048)/2000
    for name,q in configurations.items():
        config={'simulation_mode':MODE,'quantitative':deepcopy(q),'output_dir':str((output/name).resolve())}
        path=output/f'{name}.json';path.write_text(json.dumps(config,indent=2)+'\n')
        result=run_headless(config)[0]
        labels=q.get('spatial',{}).get('cell_labels',['WM','GM','CSF'])
        weights=q.get('spatial',{}).get('spatial_weights',[q.get('tissue_fractions',{}).get(k,0.) for k in labels])
        offsets=q.get('spatial',{}).get('b0_offsets_hz',[0.]*len(labels))
        expected_m=np.zeros(2048,complex);expected_w=np.zeros(2048,complex)
        for tissue,fraction,offset in zip(labels,weights,offsets):
            modulation=np.exp(2j*np.pi*offset*t)
            def R(profile):return (1-np.exp(-2/profile['T1_s']))*np.exp(-.03/profile['T2_s'])
            relax=q.get('relaxation')
            for m,c in C[tissue].items():
                expected_m+=fraction*.001*c*(R(relax['metabolites'][tissue][m]) if relax else 1)*unit[m]*modulation
            expected_w+=fraction*DEFAULT_WATER[tissue]*(R(relax['water'][tissue]) if relax else 1)*water*modulation
        clean_m=result.get('clean_metabolite_fid',result['metabolite_fid']);clean_w=result.get('clean_water_reference_fid',result['water_reference_fid'])
        errors=[]
        for label,key in [('metabolites','metabolite_fid'),('water_reference','water_reference_fid')]:
            loaded=mrs_io.read_FID(result['saved_paths'][label]).mrs().FID
            errors.append(float(abs(loaded-result[key]).max()))
        records.append(dict(case=name,output=result['saved_path'],configuration=str(path),
            water_initial_complex=[float(clean_w[0].real),float(clean_w[0].imag)],
            independent_metabolite_error=float(abs(clean_m-expected_m).max()),
            independent_water_error=float(abs(clean_w-expected_w).max()),export_max_complex_error=max(errors)))
    if before!={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}:raise RuntimeError('calibration bundle changed')
    validation=dict(cases=records,bundle_unchanged=True,source_bundle_hashes=before,fitting='not performed')
    (output/'validation_results.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation['cases'],indent=2));return validation


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    validate(parser.parse_args().output)

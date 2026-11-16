"""Export and independently validate pure/mixed tissue quantitative voxels."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from fsl_mrs.utils import mrs_io

from simulation.mrscloud_full_basis import load_library, sha256
from simulation.mrscloud_quantitative_voxel import simulate_voxel
from simulation.mrscloud_quantitative_pair import export_pair

BUNDLE=Path('outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative')
METABOLITES={
    'WM':{'NAA':10.,'Cr':8.,'GPC':2.,'PCh':1.,'Glu':8.},
    'GM':{'NAA':14.,'Cr':9.,'GPC':1.5,'PCh':1.2,'Glu':12.},
    'CSF':{'NAA':0.,'Cr':0.,'GPC':0.,'PCh':0.,'Glu':0.}}
WATER={'WM':39.0165,'GM':47.2305,'CSF':55.5}


def validate(output_dir):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=False)
    hashes={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}
    responses,water,calibration=load_library(BUNDLE,selected_names=list(METABOLITES['WM']))
    cases=[('pure_WM',[1.,0.,0.]),('pure_GM',[0.,1.,0.]),('pure_CSF',[0.,0.,1.]),('mixed',[.2,.3,.5])]
    rows=[]
    for label,weights in cases:
        fractions=dict(zip(['WM','GM','CSF'],weights))
        pair=simulate_voxel(fractions,calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,
            water_concentrations_mol_l=WATER,expected_calibration=calibration)
        expected_m=np.zeros(2048,complex);expected_w=np.zeros(2048,complex)
        for tissue,fraction in fractions.items():
            for name,concentration in METABOLITES[tissue].items():
                expected_m+=fraction*(concentration/1000)*responses[name]
            expected_w+=fraction*WATER[tissue]*water
        paths=export_pair(pair,out/label)
        errors={key:float(np.max(abs(mrs_io.read_FID(str(path)).mrs().FID-pair[key]))) for key,path in paths.items()}
        phase_errors={key:float(np.max(abs(np.angle(np.exp(1j*(
            np.angle(mrs_io.read_FID(str(path)).mrs().FID)-np.angle(pair[key])))))))
                      for key,path in paths.items()}
        gt=pair['ground_truth']
        rows.append({'case':label,'effective_water_M':gt['effective_voxel_water_concentration_mol_l'],
            'voxel_NAA_mM':gt['voxel_average_metabolite_concentrations_mM']['NAA'],
            'tissue_only_NAA_mM':gt['tissue_only_metabolite_concentrations_mM']['NAA'],
            'water_initial_real':float(pair['water_reference_fid'][0].real),
            'water_initial_imag':float(pair['water_reference_fid'][0].imag),
            'NAA_initial_complex':str(pair['components']['NAA'][0]),
            'independent_metabolite_max_error':float(np.max(abs(pair['metabolite_fid']-expected_m))),
            'independent_water_max_error':float(np.max(abs(pair['water_reference_fid']-expected_w))),
            'export_max_error':max(errors.values()),'export_max_phase_error_rad':max(phase_errors.values())})
    if hashes!={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}:
        raise RuntimeError('Phase 6a bundle changed')
    result={'calibration_identity':calibration['calibration_identity'],'shared_normalization_scalar':calibration['normalization']['scalar'],
            'source_bundle_hashes':hashes,'source_bundle_unchanged':True,'cases':rows,'FSL_MRS_fitting':'not performed'}
    (out/'validation_results.json').write_text(json.dumps(result,indent=2)+'\n')
    with (out/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output_dir')
    args=parser.parse_args();validate(args.output_dir)

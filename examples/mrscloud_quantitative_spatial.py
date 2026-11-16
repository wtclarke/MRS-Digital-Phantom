"""Phase 6f deterministic spatial controls; no production or fitting changes."""
import argparse
import json
from pathlib import Path

import numpy as np
from fsl_mrs.utils import mrs_io

from examples.mrscloud_quantitative_voxel import BUNDLE, METABOLITES, WATER
from examples.mrscloud_tissue_relaxation_validation import configuration
from simulation.mrscloud_full_basis import load_library, sha256
from simulation.mrscloud_quantitative_pair import export_pair
from simulation.mrscloud_quantitative_spatial import simulate_spatial
from simulation.mrscloud_quantitative_voxel import simulate_voxel


def validate(output_dir):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=False)
    before={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}
    responses,water,calibration=load_library(BUNDLE,selected_names=list(METABOLITES['WM']))
    t=np.arange(2048)/2000;results=[]
    for label,labels,offsets,relaxation in [
        ('zero_B0',['WM','WM'],[0.,0.],None),
        ('constant_B0',['WM','WM'],[37.,37.],None),
        ('symmetric_B0',['WM','WM'],[-5.,5.],None),
        ('mixed_relaxation_B0',['WM','GM','CSF'],[-3.,7.,1.],configuration())]:
        weights=[.2,.3,.5] if len(labels)==3 else [.5,.5]
        args=dict(calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,
            water_concentrations_mol_l=WATER,relaxation=relaxation)
        pair=simulate_spatial(labels,spatial_weights=weights,b0_offsets_hz=offsets,**args)
        repeat=simulate_spatial(labels,spatial_weights=weights,b0_offsets_hz=offsets,**args)
        fractions={tissue:sum(w for tissue_i,w in zip(labels,weights) if tissue_i==tissue) for tissue in WATER}
        baseline=simulate_voxel(fractions,**args)
        expected_m=np.zeros(2048,complex);expected_w=np.zeros(2048,complex)
        for tissue,w,f in zip(labels,weights,offsets):
            modulation=np.exp(2j*np.pi*f*t)
            for n,c in METABOLITES[tissue].items():
                p=relaxation['metabolites'][tissue][n] if relaxation else None
                r=(1-np.exp(-2/p['T1_s']))*np.exp(-.03/p['T2_s']) if p else 1.
                expected_m+=w*c*.001*r*responses[n]*modulation
            p=relaxation['water'][tissue] if relaxation else None
            r=(1-np.exp(-2/p['T1_s']))*np.exp(-.03/p['T2_s']) if p else 1.
            expected_w+=w*WATER[tissue]*r*water*modulation
        directory=out/label;paths=export_pair(pair,directory)
        np.savez(directory/'validation_arrays.npz',metabolites=pair['metabolite_fid'],water=pair['water_reference_fid'],
            expected_metabolites=expected_m,expected_water=expected_w,time_s=t)
        results.append(dict(case=label,
            independent_metabolite_max_error=float(np.max(abs(pair['metabolite_fid']-expected_m))),
            independent_water_max_error=float(np.max(abs(pair['water_reference_fid']-expected_w))),
            export_max_error=max(float(np.max(abs(mrs_io.read_FID(str(path)).mrs().FID-pair[key]))) for key,path in paths.items()),
            repeat_max_error=max(float(np.max(abs(repeat[key]-pair[key]))) for key in paths),
            zero_reference_max_difference=max(float(np.max(abs(pair[key]-baseline[key]))) for key in paths) if label=='zero_B0' else None,
            initial_amplitude_max_difference=max(float(abs(pair[key][0]-baseline[key][0])) for key in paths),
            water_magnitude_at_50ms=float(abs(pair['water_reference_fid'][100])),
            nominal_metabolite_concentrations_mM=pair['ground_truth']['voxel_average_metabolite_concentrations_mM'],
            water_initial_complex=str(pair['water_reference_fid'][0])))
    assert before=={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}
    result=dict(calibration_identity=calibration['calibration_identity'],bundle_unchanged=True,bundle_hashes=before,
        cases=results,FSL_MRS_fitting='not performed: optional spatial-model-mismatch fits deferred; no fitted concentration claim')
    (out/'validation_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(results,indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output_dir')
    validate(parser.parse_args().output_dir)

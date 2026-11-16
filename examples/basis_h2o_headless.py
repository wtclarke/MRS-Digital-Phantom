"""Headless approximate BASIS-H2O example; no fitting, basis copying or regeneration."""
import argparse
from copy import deepcopy
import json
from pathlib import Path

from examples.quantitative_water_reference_production import C
from simulation.basis_h2o import ASSUMPTION
from simulation.headless import run_headless


def configuration(basis, *, position=(100,108,80),voxel_size_mm=(10,10,10),fractions=None):
    config={'water_reference_mode':'basis_h2o','path2basis':str(basis),'positions':[list(position)],
        'quantitative':{'basis_amplitude_convention':ASSUMPTION,'concentration_units':'mM',
            'frequency_origin_ppm':4.68,'metabolite_concentrations_mM':deepcopy(C),
            'bandwidth_hz':2000,'points':2048,'receiver_factor':1.}}
    if fractions is None:config.update(skeleton='MRiLab',voxel_size_mm=list(voxel_size_mm))
    else:config['quantitative']['tissue_fractions']=dict(zip(('WM','GM','CSF'),fractions))
    return config


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--basis',type=Path,required=True)
    p.add_argument('--position',nargs=3,type=int,default=[100,108,80])
    p.add_argument('--voxel-size',nargs=3,type=float,default=[10,10,10],metavar=('X','Y','Z'))
    p.add_argument('--fractions',nargs=3,type=float,metavar=('WM','GM','CSF'),help='Small anatomy-free smoke example')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    cfg=configuration(args.basis,position=args.position,voxel_size_mm=args.voxel_size,fractions=args.fractions)
    cfg['output_dir']=str(args.output)
    print(json.dumps(cfg,indent=2))
    pair=run_headless(cfg)[0]
    print('Output:',pair['saved_path'])
    print('Ground truth:',json.dumps(pair['ground_truth']['voxel_average_metabolite_concentrations_mM']))
    print('Effective water (mol/L):',pair['ground_truth']['effective_voxel_water_concentration_mol_l'])

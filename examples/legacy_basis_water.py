"""Legacy components/shim with a separate semi-quantitative BASIS water acquisition."""
import argparse
from copy import deepcopy
from simulation.headless import simulate_batch

DEFAULTS = {'path2basis': '', 'basis_set_settings': {'vendor': 'GE', 'localization': 'PRESS', 'TE': 35}, 'voxel_size_mm': (10.0, 10.0, 10.0), 'simulation_params': {'spectral_points': 2048, 'TR': 2000.0, 'bandwidth': 2000.0, 'noise_level': 1.0, 'mm_level': 25.0, 'mm_json_file': './data/macromolecules/mm_params.json', 'lipid_amp_factor': 40.0, 'lipid_sigma': 5.0, 'lipid_phase_min': -180.0, 'lipid_phase_max': 180.0, 'lipid_lw_min': 30.0, 'lipid_lw_max': 50.0, 'water_amp_factor': 5.0, 'water_phase_min': -180.0, 'water_phase_max': 180.0, 'water_damping_min': 0.0, 'water_damping_max': 50.0, 'shim_amplitude_hz': 3.0, 'shim_corr_length': 1.0, 'shim_boundary_amp_factor': 1.5, 'shim_boundary_smoothing': 2.0}}

def configuration(basis):
    settings=deepcopy(DEFAULTS)
    settings['path2basis']=str(basis)
    settings['water_reference_mode']='basis_h2o'
    settings['water_reference']={
        'concentration_convention':'assume_table_values_nominal_mM',
        'basis_amplitude_convention':'assume_common_scale_protons_encoded',
        'frequency_origin_ppm':4.68}
    return settings

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--basis',required=True);p.add_argument('--output',required=True)
    p.add_argument('--position',nargs=3,type=int,default=[100,108,80]);p.add_argument('--seed',type=int,default=123)
    a=p.parse_args()
    r=simulate_batch([a.position],configuration(a.basis),save_dir=a.output,seed=a.seed)[0]
    print(r['saved_path'])

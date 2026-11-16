"""Run with pixi run python -m examples.mrscloud_quantitative_pair OUTPUT_DIR."""

import argparse
import json

import numpy as np
from fsl_mrs.utils import mrs_io

from simulation.mrscloud_quantitative_pair import generate_pair, export_pair


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir')
    args = parser.parse_args()
    pair = generate_pair({'NAA': 10}, water_concentration_mol_l=40, water_linewidth_hz=5)
    paths = export_pair(pair, args.output_dir)
    naa = pair['metabolite_fid'][0] / .010
    amplitude = abs(pair['water_reference_fid'][0])
    ratio = pair['metabolite_fid'][0] / pair['water_reference_fid'][0]
    q = pair['simulation_config']['synthetic_calibration']['proton_response_sru_per_mol_l']
    recovered = ratio * 40 * 2 * q / naa * 1000
    errors = {key: float(np.max(abs(mrs_io.read_FID(str(path)).mrs().FID - pair[key])))
              for key, path in paths.items()}
    print(json.dumps({'frozen_NAA_initial': str(naa), 'water_initial_amplitude': amplitude,
                      'complex_ratio': str(ratio), 'recovered_mM': str(recovered),
                      'max_complex_reload_errors': errors,
                      'water_frequency_offset_hz': pair['simulation_config']['water_frequency_offset_hz']}, indent=2))


if __name__ == '__main__':
    main()

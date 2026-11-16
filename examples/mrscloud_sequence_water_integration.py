"""Package Phase 5e/5f responses and validate both integrated water modes."""

import argparse
import csv
import json
from pathlib import Path
import shlex
import shutil
import subprocess

import fsl_mrs
import numpy as np
from scipy.io import loadmat
from fsl_mrs.utils import mrs_io, fitting

from examples.mrscloud_sequence_matched_water import assemble_common
from examples.fsl_mrs_water_preprocessing import cli_command
from examples.fsl_mrs_water_fit_diagnostics import fit_water
from simulation.mrscloud_quantitative_pair import generate_pair, export_pair
from simulation.mrscloud_sequence_reference import CONVENTION, file_sha256


def build_reference_bundle(directory, source='outputs/mrscloud_sequence_matched_water',
                           capture='outputs/mrscloud_water_processing/final/octave/processing_trace.mat'):
    """Bind unchanged Phase 5e basis and pre-DC Phase 5f water to one identity."""
    directory, source = Path(directory), Path(source)
    directory.mkdir(parents=True, exist_ok=True)
    if any(directory.iterdir()):
        raise FileExistsError('reference bundle directory must be empty')
    result = json.loads((source/'validated/results.json').read_text())
    raw = loadmat(source/'octave/raw_sequence.mat', simplify_cells=True)
    water = loadmat(capture, simplify_cells=True)
    np.testing.assert_array_equal(water['after_first_dc'], raw['signals']['water_on'])
    _, scalar, _, processed = assemble_common(raw)
    config = result['configuration']
    if scalar != config['shared_max_real_spectrum_scalar'] or water['points'] != config['native_points']:
        raise ValueError('sequence capture and reference normalization mismatch')
    if water['sw'] != config['native_bandwidth_hz'] or water['spatial_points'] != config['spatial_points_per_axis']:
        raise ValueError('sequence capture acquisition mismatch')
    shutil.copyfile(source/'validated/common_response.BASIS', directory/'common_response.BASIS')
    np.savez(directory/'raw_water.npz', raw_water=np.asarray(water['raw_fid'], dtype=np.complex128),
             normalization_reference_fids=np.column_stack([processed[name] for name in config['normalization_reference_names']]))
    manifest = {
        'convention': CONVENTION, 'fsl_mrs_version': fsl_mrs.__version__,
        'basis_sha256': file_sha256(directory/'common_response.BASIS'),
        'raw_water_sha256': file_sha256(directory/'raw_water.npz'),
        'normalization_scalar': float(scalar),
        'normalization_statistic': 'max_real_native_spectrum_fixed_three_response_reference',
        'normalization_reference_names': config['normalization_reference_names'],
        'sequence_configuration': config, 'source_sha256': result['source_sha256'],
        'spectral_grid': {'bandwidth_hz': 2000, 'points': 2048},
        'response_names': ['NAA', 'Ref0ppm'],
        'reference_signal_names': ['Ref0ppm'],
        'water_source': 'captured before both DC corrections; sysH2O already encodes two protons',
    }
    (directory/'reference.json').write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    return directory


def validate(output_dir):
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    reference = build_reference_bundle(directory/'reference')
    rows, commands = [], []
    for mode in ['analytical', 'mrscloud_sequence']:
        kwargs = {'sequence_reference': reference} if mode == 'mrscloud_sequence' else {}
        pair = generate_pair({'NAA': 10}, water_reference_mode=mode, **kwargs)
        case = directory/mode
        export_pair(pair, case, aligned_water=True)
        command = cli_command(case)
        if mode == 'mrscloud_sequence':
            command[command.index('--basis')+1] = str(reference/'common_response.BASIS')
        commands.append(shlex.join(command))
        completed = subprocess.run(command, capture_output=True, text=True)
        (case/'cli.log').write_text(completed.stdout+completed.stderr)
        completed.check_returncode()
        report = case/'fsl_mrs_results/report.html'
        if not report.exists() or not report.stat().st_size:
            raise RuntimeError('CLI HTML report missing')
        with (case/'fsl_mrs_results/concentrations.csv').open() as stream:
            records = list(csv.reader(stream))
        concentrations = dict(zip(records[1][1:5], map(float, records[2][1:5])))
        basis_path = command[command.index('--basis')+1]
        mrs = mrs_io.read_FID(str(case/'metabolites.nii.gz')).mrs(
            basis=mrs_io.read_basis(basis_path), ref_data=str(case/'water_reference.nii.gz'))
        mrs.keep = ['NAA']; mrs.rescaleForFitting()
        water_info, _, _ = fit_water(mrs)
        fit = fitting.fit_FSLModel(mrs, method='Newton', model='lorentzian', baseline_order=-1,
                                  ppmlim=(.2, 4.2), metab_groups=[0], capture_minimize_output=True)
        max_export_error = max(float(np.max(abs(mrs_io.read_FID(str(case/filename)).mrs().FID-pair[key])))
                               for key, filename in [('metabolite_fid','metabolites.nii.gz'),
                                                     ('water_reference_fid','water_reference.nii.gz')])
        row = {'mode': mode, 'initial_water_complex': str(pair['water_reference_fid'][0]),
               'initial_water_amplitude': float(abs(pair['water_reference_fid'][0])),
               'NAA_initial_complex': str(pair['metabolite_fid'][0]),
               'water_fit_amplitude': water_info['fitted_amplitude'],
               'water_relative_residual': water_info['relative_complex_residual'],
               'water_fit_frequency_hz': water_info['fitted_frequency_hz'],
               'water_fit_phase_rad': water_info['fitted_phase_rad'],
               'water_fit_linewidth_hz': water_info['fitted_linewidth_hz'],
               'water_optimizer_success': water_info['success'],
               'NAA_internal_coefficient': float(fit.getConc(metab='NAA')),
               'NAA_raw_coefficient': concentrations['raw'],
               'NAA_relative_residual': float(np.linalg.norm(mrs.FID-fit.predictedFID(mrs))/np.linalg.norm(mrs.FID)),
               'standard_fsl_mM': concentrations['molarity'], 'max_export_error': max_export_error}
        if mode == 'mrscloud_sequence':
            previous = np.load('outputs/mrscloud_water_processing/final/E_raw_acquisition_grid.npz')
            row['phase5f_water_max_error'] = float(np.max(abs(pair['water_reference_fid']-40*previous['water'])))
            row['phase5f_metabolite_max_error'] = float(np.max(abs(pair['metabolite_fid']-.010*previous['fixed_NAA'])))
        rows.append(row)
        (case/'integration_validation.json').write_text(json.dumps({'comparison':row,'water_diagnostic':water_info}, indent=2)+'\n')
    (directory/'cli_commands.txt').write_text('\n'.join(commands)+'\n')
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with (directory/'comparison.csv').open('w') as stream:
        writer = csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    (directory/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows,indent=2))
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir')
    validate(parser.parse_args().output_dir)

"""Compare original/aligned water with real FSL-MRS CLI HTML reports."""

import argparse
from copy import deepcopy
import csv
import json
from pathlib import Path
import shlex
import subprocess

import fsl_mrs
import numpy as np
from fsl_mrs.utils import fitting, mrs_io
from fsl_mrs.utils.misc import FIDToSpec

from examples.fsl_mrs_water_fit_diagnostics import fit_water
from simulation.mrscloud_quantitative_pair import REFERENCE_BASIS, generate_pair, export_pair
from simulation.water_reference_preprocessing import align_water_reference
from utils.auxillary import fid2spec


def prepare_datasets(directory):
    directory = Path(directory)
    original = generate_pair({'NAA': 10})
    aligned = deepcopy(original)
    config = original['simulation_config']
    corrected, metadata = align_water_reference(original['water_reference_fid'],
        dwell_time_s=config['dwell_time_s'], frequency_offset_hz=config['water_frequency_offset_hz'],
        initial_phase_rad=config['water_phase_rad'])
    metadata.update(acquisition_water_chemical_shift_ppm=4.70,
                    analysis_frequency_origin_ppm=config['frequency_origin_ppm'],
                    aligned_water_analysis_position_ppm=config['frequency_origin_ppm'],
                    metabolite_frame_changed=False, nifti_reference_metadata_changed=False)
    aligned['water_reference_fid'] = corrected
    aligned['water_reference_spectrum'] = fid2spec(corrected)
    # Original acquisition parameters remain provenance; digital corrections
    # describe the separately represented processed water signal.
    aligned['simulation_config']['water_preprocessing'] = metadata
    for name, pair in [('original', original), ('aligned', aligned)]:
        export_pair(pair, directory / name)
    (directory / 'aligned' / 'preprocessing_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return original, aligned


def cli_command(dataset):
    return ['pixi', 'run', 'fsl_mrs', '--data', str(dataset/'metabolites.nii.gz'),
            '--h2o', str(dataset/'water_reference.nii.gz'), '--basis', str(REFERENCE_BASIS),
            '--output', str(dataset/'fsl_mrs_results'), '--keep', 'NAA', '--algo', 'Newton',
            '--lorentzian', '--baseline', 'off', '--ppmlim', '0.2', '4.2', '--no_conj_basis',
            '--TE', '30', '--TR', '2', '--tissue_frac', '1', '0', '0', '--internal_ref', 'NAA',
            '--wref_metabolite', 'NAA', '--ref_protons', '3', '--ref_int_limits', '1.8', '2.2', '--report']


def validate(output_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    help_result = subprocess.run(['pixi', 'run', 'fsl_mrs', '--help'], capture_output=True, text=True, check=True)
    (root/'cli_help.txt').write_text(help_result.stdout + help_result.stderr)
    if '--report ' not in help_result.stdout:
        raise RuntimeError('Installed CLI does not advertise --report; inspect cli_help.txt')
    pairs = prepare_datasets(root)
    commands = [cli_command(root/name) for name in ['original', 'aligned']]
    (root/'cli_commands.txt').write_text('\n'.join(shlex.join(c) for c in commands)+'\n')
    rows, details = [], {}
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for i, (name, pair, command) in enumerate(zip(['original', 'aligned'], pairs, commands)):
        dataset = root/name
        completed = subprocess.run(command, capture_output=True, text=True)
        (dataset/'cli.log').write_text(completed.stdout + completed.stderr)
        completed.check_returncode()
        result_dir = dataset/'fsl_mrs_results'
        report_path = result_dir/'report.html'
        if not report_path.is_file() or not report_path.stat().st_size:
            raise RuntimeError(f'CLI report missing: {report_path}')
        html = report_path.read_text()
        if 'Quantification information' not in html or 'Water referencing' not in html:
            raise RuntimeError('HTML report missing quantification/reference sections')
        with (result_dir/'concentrations.csv').open() as stream:
            records = list(csv.reader(stream))
        concentrations = dict(zip(records[1][1:5], map(float, records[2][1:5])))
        data = mrs_io.read_FID(str(dataset/'metabolites.nii.gz'))
        mrs = data.mrs(basis=mrs_io.read_basis(str(REFERENCE_BASIS)),
                       ref_data=str(dataset/'water_reference.nii.gz'))
        mrs.keep = ['NAA']
        mrs.rescaleForFitting()
        # CLI does not serialize water optimizer parameters. Repeat the exact
        # installed WaterRef handler for diagnostics, without changing p0.
        water_info, water_ref, model = fit_water(mrs)
        fit = fitting.fit_FSLModel(mrs, method='Newton', model='lorentzian', baseline_order=-1,
                                  ppmlim=(.2, 4.2), metab_groups=[0], capture_minimize_output=True)
        cli_fit = np.asarray(mrs_io.read_FID(str(result_dir/'fit.nii.gz')).mrs().FID)
        api_fit = fit.predictedFID(mrs) / mrs.scaling['FID']
        np.testing.assert_allclose(cli_fit, api_fit, rtol=1e-12, atol=1e-12)
        raw_cli_water = np.asarray(mrs_io.read_FID(str(result_dir/'h2o.nii.gz')).mrs().FID)
        np.testing.assert_allclose(raw_cli_water, pair['water_reference_fid'], rtol=1e-14, atol=1e-13)
        k = pair['simulation_config']['synthetic_calibration']['frozen_basis']['normalization_scalar']
        row = {'dataset': name, 'expected_water_amplitude': float(abs(pair['water_reference_fid'][0])),
               'fitted_water_amplitude': water_info['fitted_amplitude'],
               'water_frequency_hz': water_info['fitted_frequency_hz'],
               'water_phase_rad': water_info['fitted_phase_rad'],
               'water_linewidth_hz': water_info['fitted_linewidth_hz'],
               'water_relative_complex_residual': water_info['relative_complex_residual'],
               'water_report_overlay_relative_difference': water_info['integrated_model_relative_difference'],
               'raw_NAA_source_coefficient': concentrations['raw'],
               'internal_NAA_coefficient': float(fit.getConc(metab='NAA')),
               'fsl_standard_mM': concentrations['molarity'],
               'synthetic_water_referenced_mM': concentrations['raw']*k*1000
                   *abs(pair['water_reference_fid'][0])/water_info['fitted_amplitude'],
               'metabolite_relative_residual': float(np.linalg.norm(cli_fit-pair['metabolite_fid'])/np.linalg.norm(pair['metabolite_fid'])),
               'water_optimizer_success': water_info['success'],
               'metabolite_optimizer_success': bool(fit.scipy_minimize_result.success),
               'cli_exit_code': completed.returncode,
               'water_export_max_complex_error': float(np.max(abs(raw_cli_water-pair['water_reference_fid'])))}
        rows.append(row)
        details[name] = {'water_diagnostic': water_info, 'comparison': row,
                         'water_diagnostic_provenance': 'supplementary API replay of unmodified CLI handler',
                         'metabolite_optimizer_message': str(fit.scipy_minimize_result.message),
                         'cli_vs_api_predicted_fid_max_difference': float(np.max(abs(cli_fit-api_fit)))}
        np.savez(dataset/'diagnostic_arrays.npz', water_input=water_ref.original_fid,
                 water_fit=model, water_integration_model=water_ref.fid,
                 metabolite_input=pair['metabolite_fid'], metabolite_cli_fit=cli_fit)
        ppm = mrs.getAxes()
        selected = (ppm>4.50)&(ppm<4.86)
        for col, water_model, title in [(0, water_ref.fid, 'HTML integration overlay'),
                                         (1, model, 'Actual fit in the data frame')]:
            axes[i,col].plot(ppm[selected], FIDToSpec(water_ref.original_fid.copy()).real[selected], label='water data')
            axes[i,col].plot(ppm[selected], FIDToSpec(water_model.copy()).real[selected], '--', label='water model')
            axes[i,col].set_title(f'{name}: {title}', fontsize=9)
        axes[i,2].plot(ppm, FIDToSpec(pair['metabolite_fid'].copy()).real, label='metabolite data')
        axes[i,2].plot(ppm, FIDToSpec(cli_fit.copy()).real, '--', label='CLI fit')
        axes[i,2].set_xlim(4.2,.2)
        axes[i,2].set_title(f'{name}: metabolite CLI fit')
        for col in [0,1]:
            axes[i,col].invert_xaxis()
        for ax in axes[i]:
            ax.set_xlabel('Chemical shift (ppm)')
            ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(root/'water_fit_comparison.png', dpi=160)
    plt.close(fig)
    with (root/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    summary={'FSL_MRS_version': fsl_mrs.__version__, 'report_flag': '--report (singular; --reports unsupported)',
             'datasets': details, 'max_magnitude_change': float(np.max(abs(abs(pairs[0]['water_reference_fid'])-abs(pairs[1]['water_reference_fid'])))),
             'metabolite_array_difference': float(np.max(abs(pairs[0]['metabolite_fid']-pairs[1]['metabolite_fid'])))}
    (root/'validation_results.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return summary


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir')
    validate(parser.parse_args().output_dir)

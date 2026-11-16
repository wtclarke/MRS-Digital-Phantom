"""Execute synthetic coefficient and standard FSL-MRS water-scaling validation."""

import argparse
import csv
import json
import platform
from pathlib import Path
from unittest.mock import patch
import warnings

import fsl_mrs
import numpy as np
from fsl_mrs.utils import fitting, quantify, mrs_io
from fsl_mrs.utils.misc import rescale_FID

from simulation.mrscloud_quantitative_pair import REFERENCE_BASIS, generate_pair, export_pair


def fit_case(concentration, directory):
    pair = generate_pair({'NAA': concentration}, water_concentration_mol_l=40)
    paths = export_pair(pair, directory)
    basis = mrs_io.read_basis(str(REFERENCE_BASIS))
    mrs = mrs_io.read_FID(str(paths['metabolite_fid'])).mrs(
        basis=basis, ref_data=str(paths['water_reference_fid']))
    mrs.keep = ['NAA']
    mrs.rescaleForFitting()
    scaling = {key: float(value) for key, value in mrs.scaling.items()}
    k = pair['simulation_config']['synthetic_calibration']['frozen_basis']['normalization_scalar']
    expected_coefficient = concentration * .001 / k * scaling['FID'] / scaling['basis']
    water_optimizations = []
    original_minimize = quantify.minimize
    def capture(*args, **kwargs):
        result = original_minimize(*args, **kwargs)
        water_optimizations.append({'params': result.x.tolist(), 'success': bool(result.success),
                                    'message': str(result.message)})
        return result
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit = fitting.fit_FSLModel(mrs, method='Newton', model='lorentzian', baseline_order=-1,
                                  ppmlim=(.2, 4.2), metab_groups=[0], capture_minimize_output=True)
        quant = quantify.QuantificationInfo(te=.030, tr=2.0, metabolites=['NAA'], cf=mrs.centralFrequency,
            water_ref_metab='NAA', water_ref_metab_protons=3, water_ref_metab_limits=(1.8, 2.2))
        quant.set_fractions({'WM': 1.0, 'GM': 0.0, 'CSF': 0.0})
        with patch.object(quantify, 'minimize', side_effect=capture):
            fit.calculateConcScaling(mrs, quant_info=quant, internal_reference=['NAA'])
    refs = fit.concScalings['ref_info']
    water_fit = refs['water_ref']
    water_amplitude = float(abs(water_fit.fid[0]) / scaling['FID'])
    # Correct for the acquisition's arbitrary scaling, then invert the source
    # coefficient convention; no fitted concentration sets a calibration factor.
    coefficient = float(fit.getConc(metab='NAA'))
    synthetic_mM = coefficient * scaling['basis'] / scaling['FID'] * k * 1000
    independent_mM = synthetic_mM * abs(pair['water_reference_fid'][0]) / water_amplitude
    reported_mM = float(fit.getConc(scaling='molarity', metab='NAA'))
    area_m = float(refs['metab_ref'].integral)
    area_w = float(water_fit.integral)
    calculated_fsl = area_m / area_w * (2 / 3) * quant.add_corr * quant.relax_corr_water_molar * quant.csf_corr_molar * quant.relax_corr_metab
    fit_fid = np.asarray(fit.predictedFID(mrs))
    residual_relative = float(np.linalg.norm(mrs.FID - fit_fid) / np.linalg.norm(mrs.FID))
    final_water = water_optimizations[-1]
    row = {
        'ground_truth_mM': concentration, 'fitted_coefficient': coefficient,
        'expected_coefficient': expected_coefficient,
        'unscaled_source_coefficient': coefficient * scaling['basis'] / scaling['FID'],
        'expected_source_coefficient': concentration * .001 / k,
        'coefficient_relative_error': (coefficient / expected_coefficient - 1),
        'synthetic_coefficient_mM': synthetic_mM, 'synthetic_water_referenced_mM': independent_mM,
        'synthetic_absolute_error_mM': independent_mM - concentration,
        'synthetic_relative_error': independent_mM / concentration - 1,
        'water_amplitude': water_amplitude,
        'expected_water_amplitude': float(abs(pair['water_reference_fid'][0])),
        'water_frequency_hz': -final_water['params'][3] / (2 * np.pi),
        'water_phase_rad': final_water['params'][4],
        'water_gamma_per_s': final_water['params'][1], 'water_sigma_per_s2': final_water['params'][2],
        'fsl_molarity_mM': reported_mM, 'independent_fsl_molarity_mM': calculated_fsl,
        'fit_relative_residual': residual_relative,
        'baseline_max_abs': float(np.max(abs(fit.baseline))),
    }
    details = {'comparison': row, 'scaling': scaling, 'frozen_scalar': k,
               'metabolite_optimizer': {'success': bool(fit.scipy_minimize_result.success),
                                        'message': str(fit.scipy_minimize_result.message)},
               'water_internal_fit_scaling': float(rescale_FID(
                   water_fit.original_spec[water_fit.limits], scale=100.0)[1]),
               'fit_parameters': dict(zip(fit.fitResults.columns, fit.fitResults.iloc[0].tolist())),
               'water_optimizer': final_water, 'warnings': [str(w.message) for w in caught],
               'water_integral': area_w, 'metabolite_integral': area_m,
               'quantification': {'TE_s': .030, 'TR_s': 2.0, 'WM': 1., 'GM': 0., 'CSF': 0.,
                                 'T1_s': quant.t1, 'T2_s': quant.t2,
                                 'relax_corrected_water_mM': quant.relax_corr_water_molar,
                                 'metabolite_relaxation_factor': quant.relax_corr_metab,
                                 'csf_factor': quant.csf_corr_molar, 'additional_factor': quant.add_corr,
                                 'reference_protons': 3, 'reference_ppm_limits': [1.8, 2.2]},
               'concentration_scalings': {key: value for key, value in fit.concScalings.items()
                                          if isinstance(value, (float, str, int))}}
    (directory / 'fitting_results.json').write_text(json.dumps(details, indent=2) + '\n')
    np.savez(directory / 'fit_arrays.npz', fid=mrs.FID, predicted=fit_fid, residual=mrs.FID-fit_fid,
             baseline=fit.baseline, fitted_water=water_fit.fid / scaling['FID'])
    return details


def run_validation(output_dir):
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    cases = [fit_case(c, directory / f'NAA_{c:02d}mM') for c in [5, 10, 15, 20]]
    repeat = fit_case(10, directory / 'repeat_10mM')
    rows = [case['comparison'] for case in cases]
    determinism = max(abs(rows[1][key] - repeat['comparison'][key]) for key in rows[1])
    with open(directory / 'comparison.csv', 'w') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    config = {'FSL_MRS_version': fsl_mrs.__version__, 'python': platform.python_version(),
              'fitting_source': str(Path(fitting.__file__).resolve()),
              'quantification_source': str(Path(quantify.__file__).resolve()),
              'numpy': np.__version__, 'settings': 'Newton; lorentzian; baseline disabled; ppm .2–4.2',
              'concentrations_mM': [5, 10, 15, 20], 'water_linewidth_hz': 0,
              'water_concentration_mol_l': 40, 'receiver_factor': 1,
              'repeat_10mM_max_numeric_difference': determinism}
    (directory / 'validation_config.json').write_text(json.dumps(config, indent=2) + '\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    x = [row['ground_truth_mM'] for row in rows]
    axes[0].plot(x, [row['synthetic_water_referenced_mM'] for row in rows], 'o-', label='synthetic recovery')
    axes[0].plot(x, x, '--', label='ground truth')
    axes[1].plot(x, [row['fsl_molarity_mM'] for row in rows], 'o-', label='FSL standard scaling')
    for ax in axes:
        ax.set_xlabel('Ground truth (mM)')
        ax.set_ylabel('Recovered (mM)')
        ax.legend()
    fig.tight_layout()
    fig.savefig(directory / 'concentration_recovery.png', dpi=160)
    plt.close(fig)
    return {'cases': cases, 'configuration': config}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir')
    arguments = parser.parse_args()
    result = run_validation(arguments.output_dir)
    print(json.dumps({'configuration': result['configuration'],
                      'results': [case['comparison'] for case in result['cases']]}, indent=2))

"""Hold sequence-mode signals fixed and vary FSL-MRS quantification assumptions."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess

import fsl_mrs
import numpy as np
from fsl_mrs.utils import fitting, mrs_io, quantify
from fsl_mrs.utils.constants import H2O_MOLARITY

from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral, relaxation
from examples.fsl_mrs_water_fit_diagnostics import fit_water
from examples.fsl_mrs_water_preprocessing import cli_command
from simulation.mrscloud_quantitative_pair import generate_pair, export_pair

SOURCE = Path('outputs/mrscloud_sequence_water_integration/validated')


def quantification_info(cf, *, match_water=False, match_relaxation=False):
    """Use public overrides; keep the real TE/TR and integer proton convention."""
    keys = ['H2O_WM', 'H2O_GM', 'H2O_CSF', 'METAB']
    # Exact floating-point no-attenuation limit: exp(-TR/T1)=0 and exp(-TE/inf)=1.
    t1 = dict.fromkeys(keys, 1e-6) if match_relaxation else None
    t2 = dict.fromkeys(keys, np.inf) if match_relaxation else None
    info = quantify.QuantificationInfo(.030, 2., ['NAA'], cf, t1=t1, t2=t2,
        water_ref_metab='NAA', water_ref_metab_protons=3, water_ref_metab_limits=(1.8, 2.2))
    info.set_fractions({'WM': 1., 'GM': 0., 'CSF': 0.})
    if match_water:
        info.set_densities({'WM': 40000 / H2O_MOLARITY, 'GM': .78, 'CSF': .97})
    return info


def independent_concentration(metabolite, water, ppm, settings):
    """NumPy reconstruction: no FSL final quantification/reporting calls."""
    im = spectral_integral(metabolite, ppm, (1.8, 2.2))
    iw = spectral_integral(water, ppm, (1.65, 7.65))
    rw = relaxation(settings['TE_s'], settings['TR_s'], settings['water_T1_s'], settings['water_T2_s'])
    rm = relaxation(settings['TE_s'], settings['TR_s'], settings['metabolite_T1_s'], settings['metabolite_T2_s'])
    concentration = im / iw * (2 / 3) * settings['effective_water_mM'] * rw / rm
    return {'independent_mM': concentration, 'water_integral': iw, 'metabolite_integral': im,
            'water_relaxation_factor': rw, 'metabolite_relaxation_factor': rm}


def fit_case(directory, reference):
    config = json.loads((directory/'simulation_config.json').read_text())
    truth = json.loads((directory/'ground_truth.json').read_text())['metabolite_concentrations_mM']['NAA']
    basis = mrs_io.read_basis(str(reference/'common_response.BASIS'))
    unit = basis.get_formatted_basis(2000, 2048)[:, basis.names.index('NAA')]
    mrs = mrs_io.read_FID(str(directory/'metabolites.nii.gz')).mrs(
        basis=basis, ref_data=str(directory/'water_reference.nii.gz'))
    mrs.keep = ['NAA']; mrs.rescaleForFitting()
    fit = fitting.fit_FSLModel(mrs, method='Newton', model='lorentzian', baseline_order=-1,
                              ppmlim=(.2, 4.2), metab_groups=[0], capture_minimize_output=True)
    original_params = fit.params.copy()
    water_info, water_ref, raw_water_model = fit_water(mrs)
    params = dict(zip(fit.fitResults.columns, fit.fitResults.iloc[0].tolist()))
    scale_d, scale_b = mrs.scaling['FID'], mrs.scaling['basis']
    source_coefficient = float(fit.getConc(metab='NAA') * scale_b / scale_d)
    time = np.arange(2048) * .0005
    # Reconstruct reference models from basis/optimizer parameters, not integral scalars.
    metabolite = source_coefficient * scale_d * unit * np.exp(-(params['gamma_0'] + 1j*params['eps_0'])*time)
    amp, gamma, sigma, _, _ = water_info['parameters']
    water = amp / water_info['temporary_water_scaling'] * np.exp(-gamma*time-sigma*time*time)
    ppm = np.fft.fftshift(np.fft.fftfreq(2048, .0005))/config['spectrometer_frequency_mhz'] + 4.68
    rows = []
    for label, mw, mr in [('A_defaults', False, False), ('B_water_only', True, False),
                          ('C_relaxation_only', False, True), ('D_both_matched', True, True)]:
        info = quantification_info(mrs.centralFrequency, match_water=mw, match_relaxation=mr)
        fit.calculateConcScaling(mrs, quant_info=info, internal_reference=['NAA'])
        refs = fit.concScalings['ref_info']
        np.testing.assert_array_equal(fit.params, original_params)
        np.testing.assert_array_equal(refs['water_ref'].fid, water_ref.fid)
        settings = dict(TE_s=info.te, TR_s=info.tr, water_T1_s=info.t1['H2O_WM'],
                        water_T2_s=info.t2['H2O_WM'], metabolite_T1_s=info.t1['METAB'],
                        metabolite_T2_s=info.t2['METAB'], effective_water_mM=float(H2O_MOLARITY*info.d_WM))
        independent = independent_concentration(metabolite, water, ppm, settings)
        reported = float(fit.getConc(scaling='molarity', metab='NAA'))
        row = dict(case=label, ground_truth_mM=truth, reported_mM=reported,
                   difference_from_truth_mM=reported-truth, relative_error=reported/truth-1,
                   effective_water_M=settings['effective_water_mM']/1000,
                   **independent, reconstruction_error_mM=independent['independent_mM']-reported,
                   fsl_water_integral=float(refs['water_ref'].integral),
                   fsl_metabolite_integral=float(refs['metab_ref'].integral),
                   tissue_water_factor=settings['effective_water_mM']/40000,
                   relaxation_factor=independent['water_relaxation_factor']/independent['metabolite_relaxation_factor'],
                   additional_factor=info.add_corr)
        rows.append(row)
    j_m = spectral_integral(unit*np.exp(-(params['gamma_0']+1j*params['eps_0'])*time), ppm, (1.8, 2.2))
    water_unit = water/scale_d/40
    j_w = spectral_integral(water_unit, ppm, (1.65, 7.65))
    ideal_water = abs(mrs.H2O[0]/scale_d/40)*np.exp(-np.pi*time)
    j_w_ideal = spectral_integral(ideal_water, ppm, (1.65, 7.65))
    j_m_ideal = spectral_integral(unit, ppm, (1.8, 2.2))
    # Control processing separately: same native basis, exact target-time samples.
    native_exact = basis.original_basis_array[2*np.arange(2048), basis.names.index('NAA')]
    j_m_exact = spectral_integral(native_exact, ppm, (1.8, 2.2))
    ideal_three_protons = 1.5*ideal_water*np.exp(2j*np.pi*(2.01-4.68)*127.731*time)
    control = spectral_integral(ideal_three_protons, ppm, (1.8, 2.2))/j_w_ideal*(2/3)
    details = {'ground_truth_mM': truth, 'cases': rows, 'fit_parameters': params,
        'scaling': mrs.scaling, 'raw_NAA_coefficient_mol_l': source_coefficient,
        'water_fit': water_info, 'metabolite_optimizer_success': bool(fit.scipy_minimize_result.success),
        'metabolite_relative_residual': float(np.linalg.norm(mrs.FID-fit.predictedFID(mrs))/np.linalg.norm(mrs.FID)),
        'response_factor_fitted': j_m/j_w*(2/3),
        'response_factor_unfitted': j_m_ideal/j_w_ideal*(2/3),
        'unit_NAA_integral': j_m_ideal, 'unit_water_integral': j_w_ideal,
        'input_data_metabolite_integral': spectral_integral(mrs.FID/scale_d, ppm, (1.8,2.2)),
        'input_data_water_integral': spectral_integral(mrs.H2O/scale_d, ppm, (1.65,7.65)),
        'fitted_unit_NAA_integral': j_m, 'fitted_unit_water_integral': j_w,
        'ideal_3H_2H_window_control_factor': control,
        'native_exact_NAA_integral': j_m_exact,
        'formatting_integral_factor': j_m_ideal/j_m_exact,
        'native_response_vs_ideal_3H_factor': (j_m_exact/j_w_ideal*2/3)/control,
        'coefficient_fit_factor': source_coefficient/(truth*.001),
        'nuisance_fit_response_factor': (j_m/j_w)/(j_m_ideal/j_w_ideal),
        'input_water_initial_amplitude': float(abs(mrs.H2O[0])/scale_d),
        'calibration': config['synthetic_calibration']}
    np.savez(directory/'quantification_arrays.npz', metabolite_reference=metabolite,
             water_reference_model=water, water_original=mrs.H2O, water_raw_fit=raw_water_model,
             ppm=ppm, unit_NAA=unit, unit_water=ideal_water)
    (directory/'quantification_results.json').write_text(json.dumps(details, indent=2, allow_nan=False)+'\n')
    return details


def run_validation(output_dir, source=SOURCE):
    source, out = Path(source), Path(output_dir)
    out.mkdir(parents=True, exist_ok=False)
    inputs = sorted((source/'mrscloud_sequence').glob('*.nii.gz')) + sorted((source/'mrscloud_sequence').glob('*.json'))
    inputs += sorted((source/'reference').iterdir())
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    reference = source/'reference'
    cases = []
    for c in [5, 10, 15, 20]:
        directory = out/f'NAA_{c:02d}mM'
        if c == 10:
            directory.mkdir()
            for name in ['metabolites.nii.gz', 'water_reference.nii.gz', 'ground_truth.json', 'simulation_config.json']:
                shutil.copyfile(source/'mrscloud_sequence'/name, directory/name)
            command = cli_command(directory)
            command[command.index('--basis')+1] = str(reference/'common_response.BASIS')
            (out/'cli_commands.txt').write_text(shlex.join(command)+'\n')
            completed = subprocess.run(command, capture_output=True, text=True)
            (directory/'cli.log').write_text(completed.stdout+completed.stderr)
            completed.check_returncode()
            if (directory/'fsl_mrs_results/report.html').stat().st_size < 1000:
                raise RuntimeError('baseline CLI report missing or empty')
        else:
            export_pair(generate_pair({'NAA': c}, water_reference_mode='mrscloud_sequence',
                sequence_reference=reference), directory)
        cases.append(fit_case(directory, reference))
    with (out/'comparison.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=cases[0]['cases'][0].keys())
        writer.writeheader(); writer.writerows(row for case in cases for row in case['cases'])
    # Preserve the actual CLI result separately from API assumption changes.
    with (out/'NAA_10mM/fsl_mrs_results/concentrations.csv').open() as stream:
        records = list(csv.reader(stream))
    cli_mM = dict(zip(records[1][1:5], map(float, records[2][1:5])))['molarity']
    after = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    if after != hashes:
        raise RuntimeError('existing Phase 5g input changed')
    result = {'FSL_MRS_version': fsl_mrs.__version__, 'source': str(source),
        'source_hashes': hashes, 'source_hashes_unchanged': True, 'cases': cases,
        'baseline_cli_mM': cli_mM,
        'quantification_configuration': {'TE_s': .030, 'TR_s': 2., 'fractions': {'WM': 1, 'GM': 0, 'CSF': 0},
            'matched_WM_density': 40000/H2O_MOLARITY, 'matched_T1_s': 1e-6, 'matched_T2_s': 'infinity',
            'protons_water': 2, 'protons_NAA_reference': 3, 'metabolite_window_ppm': [1.8,2.2],
            'water_window_ppm': [1.65,7.65], 'additional_correction': 1.,
            'CLI_report_quantification': 'defaults only; matched cases are Python API results'},
        'source_locations': {'quantification': quantify.__file__, 'fitting': fitting.__file__}}
    (out/'results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10,4))
    for i, label in enumerate(['defaults','water only','relaxation only','both matched']):
        axes[0].plot([5,10,15,20], [c['cases'][i]['reported_mM'] for c in cases], 'o-', label=label)
    axes[0].plot([5,20],[5,20],'k--',label='ground truth')
    axes[0].set(xlabel='Ground truth (mM)',ylabel='FSL-MRS molarity (mM)');axes[0].legend(fontsize=7)
    control = cases[1]
    names = ['window-only 3H/2H control','native exact-time basis','formatted basis']
    vals = [control['ideal_3H_2H_window_control_factor'],
            control['native_exact_NAA_integral']/control['unit_water_integral']*2/3,
            control['response_factor_unfitted']]
    axes[1].barh(names,vals);axes[1].axvline(1,color='k',ls='--');axes[1].set_xlabel('Unit response/integration factor')
    fig.tight_layout();fig.savefig(out/'quantification_comparison.png',dpi=160);plt.close(fig)
    print(json.dumps({'baseline_cli_mM':cli_mM,'10mM_cases':cases[1]['cases'],
                      'response_factor_unfitted':control['response_factor_unfitted'],
                      'formatting_integral_factor':control['formatting_integral_factor']},indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir')
    args = parser.parse_args()
    run_validation(args.output_dir)

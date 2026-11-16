"""Phase 6c: fit immutable Phase 6b voxels; predict quantification from unit responses."""
import argparse
import csv
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import warnings

import fsl_mrs
import numpy as np
from fsl_mrs.utils import fitting, mrs_io, quantify
from fsl_mrs.utils.constants import H2O_MOLARITY, STANDARD_T1, STANDARD_T2, TISSUE_WATER_DENSITY

from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral, relaxation
from examples.fsl_mrs_water_fit_diagnostics import fit_water
from examples.mrscloud_quantitative_voxel import BUNDLE, METABOLITES, WATER
from simulation.mrscloud_full_basis import load_library, sha256
from simulation.mrscloud_quantitative_pair import export_pair
from simulation.mrscloud_quantitative_voxel import simulate_voxel

SOURCE = Path('outputs/mrscloud_quantitative_voxel/final')
NAMES = list(METABOLITES['WM'])
INPUTS = ['metabolites.nii.gz', 'water_reference.nii.gz', 'ground_truth.json', 'simulation_config.json']


def prediction(truth, responses, water, *, matched=False):
    """Independent equation A21: input truth and unfitted responses only."""
    t = np.arange(len(water)) / 2000
    ppm = np.fft.fftshift(np.fft.fftfreq(len(t), .0005))/127.731 + 4.68
    jm = spectral_integral(responses['NAA'], ppm, (1.8, 2.2))
    jw = spectral_integral(abs(water[0])*np.exp(-np.pi*t), ppm, (1.65, 7.65))
    f = truth['tissue_fractions']
    if f['CSF'] == 1:
        return {'status': 'undefined: no metabolite reference and CSF correction divides by zero'}
    rm = 1. if matched else relaxation(.03, 2., STANDARD_T1['3T']['METAB'], STANDARD_T2['3T']['METAB'])
    weighted_water = 0.
    factors = {}
    for tissue in f:
        rw = 1. if matched else relaxation(.03, 2., STANDARD_T1['3T']['H2O_'+tissue], STANDARD_T2['3T']['H2O_'+tissue])
        density = truth['tissue_specific_water_concentrations_mol_l'][tissue]*1000/H2O_MOLARITY if matched else TISSUE_WATER_DENSITY[tissue]
        weighted_water += f[tissue]*density*rw*H2O_MOLARITY
        factors[tissue] = dict(fraction=f[tissue], water_density=density, water_relaxation=rw)
    multiplier = jm/jw*(2/3)*weighted_water/rm/(1-f['CSF'])/truth['effective_voxel_water_concentration_mol_l']/1000
    return dict(unit_NAA_integral=jm, unit_water_integral=jw, response_factor=jm/jw*2/3,
                weighted_relaxed_water_mM=weighted_water, metabolite_relaxation=rm,
                csf_correction=1/(1-f['CSF']), tissue_factors=factors,
                predicted_mM={n:c*multiplier for n,c in truth['voxel_average_metabolite_concentrations_mM'].items()})


def quant_info(mrs, truth, matched):
    keys = ['H2O_WM', 'H2O_GM', 'H2O_CSF', 'METAB']
    q = quantify.QuantificationInfo(.03, 2., mrs.names, mrs.centralFrequency,
        water_ref_metab='NAA', water_ref_metab_protons=3, water_ref_metab_limits=(1.8,2.2),
        t1=dict.fromkeys(keys,1e-6) if matched else None,
        t2=dict.fromkeys(keys,np.inf) if matched else None)
    q.set_fractions(truth['tissue_fractions'])
    if matched:
        q.set_densities({k:v*1000/H2O_MOLARITY for k,v in truth['tissue_specific_water_concentrations_mol_l'].items()})
    return q


def fit_case(directory, responses, water):
    truth = json.loads((directory/'ground_truth.json').read_text())
    receiver = json.loads((directory/'simulation_config.json').read_text())['receiver_factor']
    mrs = mrs_io.read_FID(str(directory/'metabolites.nii.gz')).mrs(
        basis=mrs_io.read_basis(str(BUNDLE/'metabolite_basis.BASIS')), ref_data=str(directory/'water_reference.nii.gz'))
    mrs.keep = NAMES
    water_result, _, _ = fit_water(mrs)
    expected = abs(water[0])*truth['effective_voxel_water_concentration_mol_l']*receiver
    water_result.update(expected_amplitude=float(expected), relative_amplitude_error=water_result['fitted_amplitude']/expected-1)
    result = dict(truth=truth, receiver_factor=receiver, water_fit=water_result,
                  predictions={label:prediction(truth,responses,water,matched=matched)
                               for label,matched in [('defaults',False),('matched',True)]})
    if not np.any(mrs.FID):
        # Zero spectra cannot be normalized or provide an NAA integration reference.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:
                mrs.rescaleForFitting()
                result['zero_signal_normalization'] = str(mrs.scaling)
                result['zero_signal_finite'] = bool(np.all(np.isfinite(mrs.FID)))
            except Exception as error:
                result['zero_signal_normalization'] = f'{type(error).__name__}: {error}'
        result['warnings'] = [str(w.message) for w in caught]
        try:
            result['csf_correction'] = quant_info(mrs,truth,False).csf_corr_molar
        except (ValueError, ZeroDivisionError) as error:
            result['csf_correction_error'] = f'{type(error).__name__}: {error}'
        result['status'] = 'unsupported zero metabolite signal; quantification undefined'
        return result
    mrs.rescaleForFitting()
    with warnings.catch_warnings(record=True) as caught:
        fit = fitting.fit_FSLModel(mrs, method='Newton', model='lorentzian', baseline_order=-1,
                                  ppmlim=(.2,4.2), metab_groups=[0]*len(NAMES), capture_minimize_output=True)
    result['warnings'] = [str(w.message) for w in caught]
    sd,sb = mrs.scaling['FID'],mrs.scaling['basis']
    result.update(parameters=dict(zip(fit.fitResults.columns,fit.fitResults.iloc[0].tolist())),
                  scaling=mrs.scaling, optimizer_success=bool(fit.scipy_minimize_result.success),
                  optimizer_message=str(fit.scipy_minimize_result.message),
                  relative_complex_residual=float(np.linalg.norm(mrs.FID-fit.predictedFID(mrs))/np.linalg.norm(mrs.FID)))
    result['coefficients'] = {}
    result['GPC_PCh_parameter_correlation'] = float(fit.corr[fit.metabs.index('GPC'),fit.metabs.index('PCh')])
    for name in NAMES:
        internal = float(fit.getConc(metab=name))
        expected = receiver*truth['voxel_average_metabolite_concentrations_mM'][name]*.001*sd/sb
        raw = internal*sb/sd
        result['coefficients'][name] = dict(internal=internal, expected_internal=expected,
            raw_mol_l=raw, recovered_voxel_mM=raw*1000/receiver, relative_error=internal/expected-1)
    result['quantification'] = {}
    fitted_parameters = fit.params.copy()
    for label,matched in [('defaults',False),('matched',True)]:
        fit.calculateConcScaling(mrs,quant_info=quant_info(mrs,truth,matched),internal_reference=['NAA'])
        np.testing.assert_array_equal(fit.params,fitted_parameters)
        refs = fit.concScalings['ref_info']
        result['quantification'][label] = dict(reported_mM={n:float(fit.getConc(scaling='molarity',metab=n)) for n in NAMES},
            water_integral=float(refs['water_ref'].integral),NAA_integral=float(refs['metab_ref'].integral))
    np.savez(directory/'fit_arrays.npz', fitted=fit.predictedFID(mrs)/sd, input=mrs.FID/sd)
    return result


def run_validation(output_dir, *, cli_reports=True):
    out = Path(output_dir); out.mkdir(parents=True,exist_ok=False)
    tracked = [p for p in BUNDLE.iterdir() if p.is_file()]
    tracked += [SOURCE/c/n for c in ['pure_WM','pure_GM','pure_CSF','mixed'] for n in INPUTS]
    before = {str(p):sha256(p) for p in tracked}
    responses,water,calibration = load_library(BUNDLE,selected_names=NAMES)
    cases,commands = {},[]
    for label in ['pure_WM','pure_GM','pure_CSF','mixed']:
        directory = out/label; directory.mkdir()
        for name in INPUTS:
            shutil.copyfile(SOURCE/label/name,directory/name)
        cases[label] = fit_case(directory,responses,water)
        if cli_reports and label != 'pure_CSF':
            f = cases[label]['truth']['tissue_fractions']
            command = ['pixi','run','fsl_mrs','--data',str(directory/'metabolites.nii.gz'),
                '--h2o',str(directory/'water_reference.nii.gz'),'--basis',str(BUNDLE/'metabolite_basis.BASIS'),
                '--output',str(directory/'fsl_mrs_results'),'--keep',*NAMES,'--algo','Newton',
                '--lorentzian','--baseline','off','--ppmlim','0.2','4.2','--no_conj_basis',
                '--TE','30','--TR','2','--tissue_frac',*[str(f[t]) for t in ['WM','GM','CSF']],
                '--internal_ref','NAA','--wref_metabolite','NAA','--ref_protons','3',
                '--ref_int_limits','1.8','2.2','--report']
            commands.append(shlex.join(command))
            (out/'cli_commands.txt').write_text('\n'.join(commands)+'\n')
            completed = subprocess.run(command,capture_output=True,text=True)
            (directory/'cli.log').write_text(completed.stdout+completed.stderr)
            completed.check_returncode()
            report = directory/'fsl_mrs_results/report.html'
            if report.stat().st_size < 1000 or 'Water referencing' not in report.read_text():
                raise RuntimeError(f'Invalid HTML report: {report}')
            with (directory/'fsl_mrs_results/concentrations.csv').open() as stream:
                records = list(csv.reader(stream))
            cases[label]['cli_mM'] = {row[0]:float(row[4]) for row in records[2:] if row[0] in NAMES}
            cli_fid = mrs_io.read_FID(str(directory/'fsl_mrs_results/fit.nii.gz')).mrs().FID
            api_fid = np.load(directory/'fit_arrays.npz')['fitted']
            cases[label]['cli_api_max_complex_difference'] = float(np.max(abs(cli_fid-api_fid)))
    fractions = cases['mixed']['truth']['tissue_fractions']
    for label,cm,cw,g in [('repeat',1,1,1),('double_metabolites',2,1,1),('double_water',1,2,1),('double_receiver',1,1,2)]:
        pair = simulate_voxel(fractions,calibration_bundle=BUNDLE,
            metabolite_concentrations_mM={t:{n:cm*c for n,c in profile.items()} for t,profile in METABOLITES.items()},
            water_concentrations_mol_l={t:cw*c for t,c in WATER.items()},receiver_factor=g)
        directory = out/label; export_pair(pair,directory)
        cases[label] = fit_case(directory,responses,water)
    matrix = np.column_stack([responses[n] for n in NAMES])
    normalized = matrix/np.linalg.norm(matrix,axis=0)
    gram = normalized.conj().T@normalized
    weighted = matrix.copy(); weighted[0] *= .5
    spectrum = np.fft.fftshift(np.fft.fft(weighted,axis=0),axes=0)
    ppm = np.fft.fftshift(np.fft.fftfreq(2048,.0005))/127.731+4.68
    lo,hi = sorted(np.argmin(abs(ppm-p)) for p in [.2,4.2])
    window = spectrum[lo:hi+1]
    window /= np.linalg.norm(window,axis=0)
    real_design = np.concatenate([window.real,window.imag],axis=0)
    rows = []
    for case,data in cases.items():
        if 'coefficients' not in data:
            continue
        truth = data['truth']
        for n in NAMES:
            for setting in ['defaults','matched']:
                reported = data['quantification'][setting]['reported_mM'][n]
                predicted = data['predictions'][setting]['predicted_mM'][n]
                rows.append(dict(case=case,metabolite=n,configuration=setting,
                    voxel_mM=truth['voxel_average_metabolite_concentrations_mM'][n],
                    tissue_only_mM=truth['tissue_only_metabolite_concentrations_mM'][n],
                    **data['coefficients'][n],reported_mM=reported,predicted_mM=predicted,
                    prediction_error_mM=reported-predicted,
                    difference_from_voxel_mM=reported-truth['voxel_average_metabolite_concentrations_mM'][n],
                    difference_from_tissue_mM=reported-truth['tissue_only_metabolite_concentrations_mM'][n],
                    voxel_interpretation_mM=reported*(1-truth['tissue_fractions']['CSF'])))
    with (out/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    if before != {str(p):sha256(p) for p in tracked}:
        raise RuntimeError('Phase 6a/6b inputs changed')
    result = dict(version=fsl_mrs.__version__,cases=cases,source_hashes=before,source_unchanged=True,
        calibration_identity=calibration['calibration_identity'],
        full_FID_normalized_design_condition=float(np.linalg.cond(normalized)),
        GPC_PCh_complex_correlation=float(abs(gram[NAMES.index('GPC'),NAMES.index('PCh')])),
        fitted_window_real_design_condition=float(np.linalg.cond(real_design)),
        settings=dict(method='Newton',model='lorentzian',baseline='off',groups='one shared group',
            ppm=[.2,4.2],TE_s=.03,TR_s=2.,reference='NAA',reference_protons=3,
            metabolite_integration_ppm=[1.8,2.2],water_integration_ppm=[1.65,7.65],
            matched='public density and T1/T2 overrides; original fitted parameters unchanged'))
    (out/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output_dir')
    args=parser.parse_args(); result=run_validation(args.output_dir)
    print(json.dumps({k:{'water_residual':v['water_fit']['relative_complex_residual'],
        'metabolite_residual':v.get('relative_complex_residual'),'success':v.get('optimizer_success')}
        for k,v in result['cases'].items()},indent=2))

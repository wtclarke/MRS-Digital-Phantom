"""Trace synthetic water conventions and compare raw versus integrated FSL fits."""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import nibabel as nib
import numpy as np
import fsl_mrs
from fsl_mrs.utils import mrs_io, quantify
from fsl_mrs.utils.misc import FIDToSpec, rescale_FID

from simulation.mrscloud_quantitative_pair import generate_pair, REFERENCE_BASIS


def measure(fid, dwell):
    t = np.arange(len(fid)) * dwell
    phase = np.unwrap(np.angle(fid))
    return {'initial_amplitude': float(abs(fid[0])), 'initial_phase_rad': float(np.angle(fid[0])),
            'frequency_hz': float(np.polyfit(t, phase, 1)[0] / (2 * np.pi)),
            'linewidth_hz': float(-np.polyfit(t, np.log(abs(fid)), 1)[0] / np.pi),
            'dwell_time_s': dwell}


def fit_water(mrs, use_unwrapped_initialization=False):
    original_minimize = quantify.minimize
    captured = {}
    def capture(fun, p0, **kwargs):
        captured['default_initial_parameters'] = list(p0)
        if use_unwrapped_initialization:
            measured = measure(mrs.H2O, mrs.dwellTime)
            p0 = list(p0)
            p0[3] = -2 * np.pi * measured['frequency_hz']
            p0[4] = measured['initial_phase_rad']
        result = original_minimize(fun, p0, **kwargs)
        captured.update(parameters=result.x.tolist(), success=bool(result.success),
                        message=str(result.message), initial_parameters=list(p0))
        return result
    with patch.object(quantify, 'minimize', side_effect=capture):
        reference = quantify.WaterRef(mrs, limits=(1.65, 7.65))
    scaling = rescale_FID(reference.original_spec[reference.limits], scale=100)[1]
    amp, gamma, sigma, omega, phi = captured['parameters']
    t = reference.t_axis
    raw_model = amp * np.exp(-t * (gamma + t * sigma + 1j * omega)) * np.exp(1j * phi) / scaling
    original = reference.original_fid
    captured.update(fitted_amplitude=float(amp / scaling / mrs.scaling['FID']),
                    fitted_frequency_hz=float(-omega / (2 * np.pi)), fitted_phase_rad=float(phi),
                    fitted_linewidth_hz=float(gamma / np.pi),
                    relative_complex_residual=float(np.linalg.norm(raw_model-original) / np.linalg.norm(original)),
                    integrated_model_relative_difference=float(np.linalg.norm(reference.fid-original) / np.linalg.norm(original)),
                    integral=float(reference.integral), temporary_water_scaling=float(scaling))
    return captured, reference, raw_model


def diagnose(output_dir, source_dir=Path('outputs/fsl_mrs_water_scaling_phase5b/NAA_10mM')):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    pair = generate_pair({'NAA': 10})
    data = mrs_io.read_FID(str(source_dir / 'metabolites.nii.gz'))
    loaded = mrs_io.read_FID(str(source_dir / 'water_reference.nii.gz'))
    mrs = data.mrs(basis=mrs_io.read_basis(str(REFERENCE_BASIS)), ref_data=loaded)
    mrs.keep = ['NAA']
    mrs.rescaleForFitting()
    generator = pair['water_reference_fid']
    stored = np.asarray(nib.load(source_dir / 'water_reference.nii.gz').dataobj).squeeze()
    interpreted = np.asarray(loaded.mrs().FID)
    input_water = mrs.H2O
    dwell = mrs.dwellTime
    trace = {name: measure(array, dwell) for name, array in
             [('generator', generator), ('stored', stored), ('interpreted', interpreted), ('fitting_input', input_water)]}
    differences = {'stored_vs_generator': float(np.max(abs(stored-generator))),
                   'stored_vs_conjugate': float(np.max(abs(stored-generator.conj()))),
                   'interpreted_vs_generator': float(np.max(abs(interpreted-generator))),
                   'fitting_input_unscaled_vs_generator': float(np.max(abs(input_water/mrs.scaling['FID']-generator)))}
    np.savez(out / 'signal_trace.npz', generator=generator, stored=stored,
             interpreted=interpreted, fitting_input=input_water)
    cases = []
    fig, axes = plt.subplots(4, 2, figsize=(12, 13))
    frequency = pair['simulation_config']['water_frequency_offset_hz']
    for i, (label, phase, freq) in enumerate([('A', 0., 0.), ('B', 0., frequency),
                                             ('C', .6, 0.), ('D', .6, frequency)]):
        test_pair = generate_pair({'NAA': 10}, water_phase_rad=phase, water_frequency_offset_hz=freq)
        mrs.H2O = test_pair['water_reference_fid']
        before, reference, raw_model = fit_water(mrs)
        after, corrected, alternative_model = fit_water(mrs, use_unwrapped_initialization=True)
        case = {'case': label, 'expected': measure(test_pair['water_reference_fid'], dwell),
                'default_fit': before, 'unwrapped_initialization_fit': after}
        cases.append(case)
        ppm = mrs.getAxes()
        limits = (ppm > 4.50) & (ppm < 4.86)
        for column, model, title in [(0, reference.fid, 'Report overlay: integration model'),
                                     (1, raw_model, 'Actual fitted model before correction')]:
            ax = axes[i, column]
            ax.plot(ppm[limits], FIDToSpec(reference.original_fid.copy()).real[limits], label='original data')
            ax.plot(ppm[limits], FIDToSpec(model.copy()).real[limits], '--', label='model')
            ax.invert_xaxis()
            ax.set_title(f'{label}: phase={phase:.2f}, f={freq:.5f} Hz\n{title}', fontsize=9)
            ax.set_xlabel('Chemical shift (ppm)')
            ax.legend()
        np.savez(out / f'case_{label}.npz', original=reference.original_fid,
                 actual_fit=raw_model, integrated_fit=reference.fid, alternate_fit=alternative_model)
    fig.tight_layout()
    fig.savefig(out / 'water_fit_before_after.png', dpi=160)
    plt.close(fig)
    # Show initialization and integration-frame effects separately for case D.
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    ppm = mrs.getAxes()
    limits = (ppm > 4.50) & (ppm < 4.86)
    axes[0].plot(ppm[limits], FIDToSpec(reference.original_fid.copy()).real[limits], label='data')
    axes[0].plot(ppm[limits], FIDToSpec(raw_model.copy()).real[limits], '--', label='default fit')
    axes[1].plot(ppm[limits], FIDToSpec(reference.original_fid.copy()).real[limits], label='data')
    axes[1].plot(ppm[limits], FIDToSpec(alternative_model.copy()).real[limits], '--', label='unwrapped/sign-correct fit')
    amp, gamma, sigma, omega, phi = after['parameters']
    corrected_data = reference.original_fid * np.exp(1j * omega * reference.t_axis - 1j * phi)
    axes[2].plot(ppm[limits], FIDToSpec(corrected_data).real[limits], label='data corrected to integration frame')
    axes[2].plot(ppm[limits], FIDToSpec(corrected.fid.copy()).real[limits], '--', label='integrated model')
    for ax, title in zip(axes, ['D: default initialization', 'D: corrected initialization', 'D: consistent integration frame']):
        ax.set_title(title)
        ax.set_xlabel('Chemical shift (ppm)')
        ax.invert_xaxis()
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / 'initialization_and_frame_case_D.png', dpi=160)
    plt.close(fig)
    original_fit = cases[1]['default_fit']
    ppm = mrs.getAxes()
    measured_ppm_bin = float(ppm[np.argmax(abs(FIDToSpec(generator.copy())))])
    center_ppm_bin = float(ppm[np.argmax(abs(FIDToSpec(np.ones_like(generator))))])
    result = {'FSL_MRS_version': fsl_mrs.__version__,
              'trace': trace, 'max_complex_differences': differences, 'mrs_scaling': mrs.scaling,
              'configured': pair['simulation_config'], 'ppm_peak_bin': measured_ppm_bin,
              'ppm_center_bin': center_ppm_bin, 'cases': cases,
              'phase5b_water_amplitude_fit': original_fit['fitted_amplitude']}
    (out / 'diagnostics.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir')
    args = parser.parse_args()
    diagnose(args.output_dir)

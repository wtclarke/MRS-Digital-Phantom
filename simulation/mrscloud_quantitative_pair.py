"""Opt-in homogeneous synthetic MRSCloud/water pair; no legacy integration."""

import json
from pathlib import Path

import numpy as np
from fsl_mrs.utils.mrs_io.main import read_basis
from nifti_mrs.create_nmrs import gen_nifti_mrs

from simulation.mrscloud_water_calibration import calibrated_responses, concentration_signals
from simulation.mrscloud_sequence_reference import sequence_responses
from simulation.water_reference_preprocessing import align_water_reference
from simulation.mrscloud_quantitative_noise import add_pair_noise
from utils.auxillary import fid2spec

REFERENCE_BASIS = Path(__file__).resolve().parents[1] / 'data/basissets/LCModel_Universal_Philips_UnEdited_PRESS_TE30.BASIS'
REFERENCE_SHA256 = 'bac0042c2e9a688ef63abd88867a79fd24c659ba92da5c6a333e71f96ffa5abe'


def generate_pair(concentrations_mM, *, water_concentration_mol_l=40.0,
                  basis_path=None, bandwidth_hz=2000, points=2048,
                  receiver_factor=1.0, water_linewidth_hz=None,
                  water_frequency_offset_hz=None, water_phase_rad=0.0,
                  calibration=None, expected_metadata=None,
                  water_reference_mode='analytical', sequence_reference=None, noise=None):
    """Generate concentration-weighted FIDs and FFT spectra on one common grid.

    Analytical mode uses the verified PRESS TE30 frozen reference. Sequence
    mode requires an explicitly shared-scale bundle and retains its native
    1 Hz water decay. It rejects analytical calibration objects and overrides.
    Default water is 4.70 ppm
    relative to the basis' 4.68 ppm origin. Explicit offsets are in Hz relative
    to that origin, using exp(+2j*pi*f*t). Calibration metadata refer to the
    Phase 4b unit response before applying this acquisition's water envelope.
    """
    if water_reference_mode not in ['analytical', 'mrscloud_sequence']:
        raise ValueError('unknown water reference mode')
    if water_reference_mode == 'analytical':
        if sequence_reference is not None:
            raise ValueError('sequence reference supplied to analytical mode')
        basis_path = REFERENCE_BASIS if basis_path is None else basis_path
        water_linewidth_hz = 0.0 if water_linewidth_hz is None else water_linewidth_hz
    elif sequence_reference is None or calibration is not None:
        raise ValueError('sequence mode requires its own reference bundle, not a frozen analytical calibration')
    if isinstance(points, bool) or not isinstance(points, (int, np.integer)) or points <= 0:
        raise ValueError('points must be a positive integer')
    if water_linewidth_hz is not None and (not np.isfinite(water_linewidth_hz) or water_linewidth_hz < 0):
        raise ValueError('water linewidth must be finite and nonnegative')
    if not np.isfinite(water_phase_rad):
        raise ValueError('water phase must be finite')
    if water_reference_mode == 'analytical':
        responses, water_response, metadata = calibrated_responses(
            basis_path, bandwidth_hz=bandwidth_hz, points=points,
            selected_names=list(concentrations_mM), calibration=calibration,
            expected_metadata=expected_metadata)
        if metadata['frozen_basis']['reference_identity'].split(':')[-1] != REFERENCE_SHA256:
            raise ValueError('frequency convention verified only for the PRESS TE30 reference identity')
        cf = float(read_basis(str(basis_path)).cf)
    else:
        responses, water_response, metadata = sequence_responses(
            sequence_reference, bandwidth_hz=bandwidth_hz, points=points,
            selected_names=list(concentrations_mM), expected_metadata=expected_metadata,
            basis_path=basis_path)
        sequence_config = metadata['sequence_configuration']
        cf = sequence_config['spectrometer_frequency_mhz']
        if ((water_linewidth_hz is not None and water_linewidth_hz != sequence_config['readout_linewidth_hz'])
                or water_phase_rad != 0
                or (water_frequency_offset_hz is not None
                    and not np.isclose(water_frequency_offset_hz, (4.70 - 4.68) * cf, rtol=0, atol=1e-10))):
            raise ValueError('sequence response preserves its simulated linewidth, phase and frequency')
        water_linewidth_hz = sequence_config['readout_linewidth_hz']
        water_frequency_offset_hz = (4.70 - 4.68) * cf
    frequency = (4.70 - 4.68) * cf if water_frequency_offset_hz is None else float(water_frequency_offset_hz)
    if not np.isfinite(frequency) or abs(frequency) >= bandwidth_hz / 2:
        raise ValueError('water frequency must lie inside the sampled bandwidth')
    concentrations = {name: float(value) * 1e-3 for name, value in concentrations_mM.items()}
    components, water = concentration_signals(
        responses, water_response, concentrations_mol_l=concentrations,
        water_concentration_mol_l=water_concentration_mol_l, receiver_factor=receiver_factor)
    time = np.arange(points) / bandwidth_hz
    if water_reference_mode == 'analytical':
        water = water * np.exp(-np.pi * water_linewidth_hz * time) * np.exp(
            1j * (2 * np.pi * frequency * time + water_phase_rad))
    metabolites = sum(components.values())
    config = {
        'synthetic_calibration': metadata, 'bandwidth_hz': float(bandwidth_hz),
        'water_reference_mode': water_reference_mode,
        'points': int(points), 'dwell_time_s': 1.0 / bandwidth_hz,
        'spectrometer_frequency_mhz': cf, 'frequency_origin_ppm': 4.68,
        'frequency_convention': 'f_hz=(ppm-4.68)*spectrometer_frequency_mhz; exp(+2j*pi*f*t)',
        'water_frequency_offset_hz': frequency, 'water_linewidth_hz': float(water_linewidth_hz),
        'water_phase_rad': float(water_phase_rad), 'receiver_factor': float(receiver_factor),
        'noise': False, 'additional_metabolite_relaxation': False, 'T1_TR_attenuation': False,
        'residual_water': False, 'spatial_weighting': False,
    }
    ground_truth = {
        'metabolite_concentrations_mM': {name: float(value) for name, value in concentrations_mM.items()},
        'water_concentration_mol_l': float(water_concentration_mol_l),
        'concentration_convention': 'nominal concentration per total homogeneous tissue volume',
        'physical_calibration': False,
    }
    pair = {'metabolite_fid': metabolites, 'water_reference_fid': water,
            'metabolite_spectrum': fid2spec(metabolites), 'water_reference_spectrum': fid2spec(water),
            'components': components, 'time_s': time, 'ground_truth': ground_truth,
            'simulation_config': config}
    return add_pair_noise(pair, **noise) if noise is not None else pair


def export_pair(pair, output_dir, *, aligned_water=False):
    """Write two 1D NIfTI-MRS datasets and their JSON provenance sidecars."""
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    filenames = ['metabolites.nii.gz', 'water_reference.nii.gz', 'ground_truth.json', 'simulation_config.json']
    signals = [('metabolite_fid', 'metabolites.nii.gz'), ('water_reference_fid', 'water_reference.nii.gz')]
    if 'clean_metabolite_fid' in pair:
        signals += [('clean_metabolite_fid', 'metabolites_clean.nii.gz'),
                    ('clean_water_reference_fid', 'water_reference_clean.nii.gz')]
        filenames += ['metabolites_clean.nii.gz', 'water_reference_clean.nii.gz']
    if aligned_water:
        filenames += ['water_reference_aligned.nii.gz', 'water_preprocessing.json']
    if any((directory / name).exists() for name in filenames):
        raise FileExistsError('pair output already exists')
    config = pair['simulation_config']
    paths = {}
    for key, filename in signals:
        obj = gen_nifti_mrs(pair[key][None, None, None, :, None], config['dwell_time_s'],
                            config['spectrometer_frequency_mhz'], nucleus='1H',
                            dim_tags=['DIM_DYN', None, None])
        obj.add_hdr_field('SpecFreqChemShift', config['frequency_origin_ppm'])
        obj.add_hdr_field('EchoTime', config.get('TE_s', .030))
        if 'TR_s' in config:
            obj.add_hdr_field('RepetitionTime', config['TR_s'])
        obj.save(str(directory / filename))
        paths[key] = directory / filename
    for key in ['ground_truth', 'simulation_config']:
        (directory / f'{key}.json').write_text(json.dumps(pair[key], indent=2) + '\n')
    if aligned_water:
        aligned, preprocessing = align_water_reference(
            pair['water_reference_fid'], dwell_time_s=config['dwell_time_s'],
            frequency_offset_hz=config['water_frequency_offset_hz'], initial_phase_rad=config['water_phase_rad'])
        obj = gen_nifti_mrs(aligned[None, None, None, :, None], config['dwell_time_s'],
                            config['spectrometer_frequency_mhz'], nucleus='1H',
                            dim_tags=['DIM_DYN', None, None])
        obj.add_hdr_field('SpecFreqChemShift', config['frequency_origin_ppm'])
        obj.add_hdr_field('EchoTime', config.get('TE_s', .030))
        if 'TR_s' in config:
            obj.add_hdr_field('RepetitionTime', config['TR_s'])
        path = directory / 'water_reference_aligned.nii.gz'
        obj.save(str(path))
        paths['aligned_water_reference_fid'] = path
        preprocessing.update(analysis_frequency_origin_ppm=config['frequency_origin_ppm'],
                             aligned_analysis_position_ppm=config['frequency_origin_ppm'],
                             source_water_position_ppm=config['frequency_origin_ppm']+config['water_frequency_offset_hz']/config['spectrometer_frequency_mhz'],
                             physical_water_position_unchanged=True, metabolite_frame_changed=False)
        (directory / 'water_preprocessing.json').write_text(json.dumps(preprocessing, indent=2) + '\n')
    return paths

"""Load a validated, explicitly shared-scale sequence response bundle."""

import hashlib
import json
from pathlib import Path

import numpy as np
import fsl_mrs
from fsl_mrs.utils.mrs_io.main import read_basis

CONVENTION = 'mrscloud-press-te30-raw-water-common-max-v1'


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sequence_responses(reference_dir, *, bandwidth_hz, points, selected_names,
                       expected_metadata=None, basis_path=None):
    """Read paired response provenance; never renormalize selected columns."""
    directory = Path(reference_dir)
    manifest_path = directory / 'reference.json'
    manifest = json.loads(manifest_path.read_text())
    required = {'convention', 'fsl_mrs_version', 'spectral_grid', 'basis_sha256',
                'raw_water_sha256', 'normalization_scalar', 'sequence_configuration',
                'normalization_statistic', 'normalization_reference_names', 'response_names', 'source_sha256'}
    if not required <= manifest.keys():
        raise ValueError('incomplete sequence reference metadata')
    if manifest['convention'] != CONVENTION:
        raise ValueError('incompatible sequence calibration convention')
    if manifest['fsl_mrs_version'] != fsl_mrs.__version__:
        raise ValueError('sequence reference formatting software version mismatch')
    if manifest['spectral_grid'] != {'bandwidth_hz': bandwidth_hz, 'points': points}:
        raise ValueError('incompatible sequence sampling configuration')
    source_basis = directory / 'common_response.BASIS'
    water_path = directory / 'raw_water.npz'
    if (file_sha256(source_basis) != manifest['basis_sha256']
            or file_sha256(water_path) != manifest['raw_water_sha256']):
        raise ValueError('sequence response identity mismatch')
    if basis_path is not None and file_sha256(basis_path) != manifest['basis_sha256']:
        raise ValueError('basis is not on the sequence reference scale')
    scalar = manifest['normalization_scalar']
    config = manifest['sequence_configuration']
    if (not np.isfinite(scalar) or scalar <= 0
            or config['shared_max_real_spectrum_scalar'] != scalar
            or (config['TE_ms'], config['TE1_ms'], config['TE2_ms']) != (30, 14, 16)
            or config['RF_carriers_ppm']['water_on'] != 4.7
            or config['RF_carriers_ppm']['NAA'] != 3.0
            or config['frequency_origin_ppm'] != 4.68
            or config['readout_linewidth_hz'] != 1.0):
        raise ValueError('incompatible sequence or normalization metadata')
    if (manifest['normalization_reference_names'] != ['NAA', 'Ref0ppm', 'water_off']
            or manifest['normalization_statistic'] != 'max_real_native_spectrum_fixed_three_response_reference'):
        raise ValueError('incompatible sequence normalization convention')
    names = list(selected_names)
    if not names or len(names) != len(set(names)) or not set(names) <= set(manifest['response_names']):
        raise ValueError('select distinct responses present in the sequence reference')
    basis = read_basis(str(source_basis))
    if (basis.original_points != config['native_points']
            or not np.isclose(basis.original_bw, config['native_bandwidth_hz'], rtol=1e-12, atol=1e-12)
            or not np.isclose(basis.cf, config['spectrometer_frequency_mhz'], rtol=1e-12, atol=1e-12)):
        raise ValueError('sequence basis/native grid metadata mismatch')
    formatted = basis.get_formatted_basis(bandwidth_hz, points)
    with np.load(water_path, allow_pickle=False) as archive:
        raw = archive['raw_water'].copy()
        reference = archive['normalization_reference_fids'].copy()
    if raw.shape != (config['native_points'],) or not np.iscomplexobj(raw) or not np.all(np.isfinite(raw)):
        raise ValueError('invalid raw sequence water response')
    reference_names = manifest['normalization_reference_names']
    if (reference.shape != (config['native_points'], len(reference_names))) or not np.all(np.isfinite(reference)):
        raise ValueError('invalid normalization reference responses')
    measured_scalar = float(np.max(np.fft.fftshift(np.fft.fft(reference, axis=0), axes=0).real))
    if not np.isclose(measured_scalar, scalar, rtol=1e-12, atol=0):
        raise ValueError('normalization scalar does not match reference responses')
    dt = 1 / config['native_bandwidth_hz']
    native_time = np.arange(config['native_points']) * dt
    for name in ['NAA', 'Ref0ppm']:
        if name not in reference_names or name not in basis.names:
            raise ValueError('missing sequence normalization anchor')
        expected = reference[:, reference_names.index(name)] / scalar * np.exp(
            -2j * np.pi * (4.68 - config['RF_carriers_ppm'][name])
            * config['spectrometer_frequency_mhz'] * native_time)
        actual = basis.original_basis_array[:, basis.names.index(name)]
        if np.linalg.norm(actual - expected) > 3e-7 * np.linalg.norm(expected):
            raise ValueError('basis does not share the declared raw reference normalization')
    stride = config['native_bandwidth_hz'] / bandwidth_hz
    step = int(round(stride))
    if step < 1 or not np.isclose(stride, step) or (points - 1) * step >= len(raw):
        raise ValueError('sequence water requires covered exact acquisition times')
    # Match Phase 5f operation order: unit-magnitude frame change, shared scale,
    # then exact acquisition-time selection. sysH2O already represents 2H.
    native = raw * np.exp(2j * np.pi * (4.7 - config['frequency_origin_ppm'])
                         * config['spectrometer_frequency_mhz'] * np.arange(len(raw)) * dt) / scalar
    water = native[np.arange(points) * step].copy()
    metadata = {
        'convention': CONVENTION,
        'calibration_identity': CONVENTION + ':' + file_sha256(manifest_path),
        'reference_basis_identity': source_basis.name + ':' + manifest['basis_sha256'],
        'normalization_scalar': float(scalar),
        'normalization_statistic': manifest['normalization_statistic'],
        'normalization_reference_names': manifest['normalization_reference_names'],
        'sequence_configuration': config,
        'spectral_grid': manifest['spectral_grid'],
        'raw_water_sha256': manifest['raw_water_sha256'],
        'source_sha256': manifest['source_sha256'],
        'water_processing': 'raw pre-DC; acquisition-frame shift; shared scalar; exact time samples',
        'water_proton_multiplier_applied': False,
        'concentration_units': 'nominal mol/L per total tissue volume',
        'physical_calibration': False,
    }
    if expected_metadata is not None and expected_metadata != metadata:
        raise ValueError('sequence calibration metadata mismatch')
    return {name: formatted[:, basis.names.index(name)].copy() for name in names}, water, metadata

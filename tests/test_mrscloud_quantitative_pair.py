import json
from dataclasses import replace

import nibabel as nib
import numpy as np
import pytest
from fsl_mrs.utils.mrs_io.main import read_basis
from fsl_mrs.utils import mrs_io

from simulation.frozen_basis_normalization import format_frozen_basis
from simulation.mrscloud_quantitative_pair import REFERENCE_BASIS, generate_pair, export_pair


@pytest.fixture(scope='module')
def independent():
    b = read_basis(str(REFERENCE_BASIS))
    raw = b.get_formatted_basis(2000, 2048)
    frozen = raw / np.mean(abs(raw))
    return {str(name): frozen[:, i] for i, name in enumerate(b.names)}, float(b.cf)


def test_single_superposition_selection_and_determinism(independent):
    basis, _ = independent
    one = generate_pair({'NAA': 10})
    multiple = generate_pair({'NAA': 10, 'Cr': 8})
    np.testing.assert_allclose(one['metabolite_fid'], .010 * basis['NAA'], rtol=1e-15)
    np.testing.assert_allclose(multiple['metabolite_fid'], .010 * basis['NAA'] + .008 * basis['Cr'], rtol=1e-15)
    np.testing.assert_array_equal(one['metabolite_fid'], multiple['components']['NAA'])
    again = generate_pair({'NAA': 10})
    for key in ['metabolite_fid', 'water_reference_fid']:
        np.testing.assert_array_equal(one[key], again[key])
    assert one['simulation_config'] == again['simulation_config']


def test_concentration_linearity_receiver_and_ratio(independent):
    basis, _ = independent
    one = generate_pair({'NAA': 10})
    twice = generate_pair({'NAA': 20}, water_concentration_mol_l=80)
    receiver = generate_pair({'NAA': 10}, receiver_factor=3)
    for key in ['metabolite_fid', 'water_reference_fid']:
        np.testing.assert_allclose(twice[key], 2 * one[key], rtol=2e-15)
        np.testing.assert_allclose(receiver[key], 3 * one[key], rtol=2e-15)
    np.testing.assert_allclose(receiver['metabolite_fid'] / receiver['water_reference_fid'],
                               one['metabolite_fid'] / one['water_reference_fid'], rtol=2e-15)
    q = abs(basis['Ref0ppm'][0])
    ratio = one['metabolite_fid'][0] / one['water_reference_fid'][0]
    expected = .010 * basis['NAA'][0] / (80 * q)
    assert ratio == pytest.approx(expected, rel=2e-15)
    recovered = ratio * 40 * (2 * q) / basis['NAA'][0] * 1000
    assert recovered == pytest.approx(10, abs=1e-12)
    print(f'NAA unit response={basis["NAA"][0]!r}; water={one["water_reference_fid"][0]!r}; '
          f'ratio={ratio!r}; recovered_mM={recovered!r}')


def test_water_linewidth_frequency_phase_and_initial_amplitude(independent):
    basis, _ = independent
    pair = generate_pair({'NAA': 10}, water_linewidth_hz=5, water_frequency_offset_hz=37,
                         water_phase_rad=.6, receiver_factor=2)
    t = np.arange(2048) * .0005
    amplitude = 2 * abs(basis['Ref0ppm'][0]) * 2 * 40
    expected = amplitude * np.exp(-np.pi * 5 * t + 1j * (2 * np.pi * 37 * t + .6))
    # Independent multiply/divide time grids differ by roundoff in long phases.
    np.testing.assert_allclose(pair['water_reference_fid'], expected, rtol=1e-13, atol=1e-14)
    np.testing.assert_allclose(abs(pair['water_reference_fid']), amplitude * np.exp(-np.pi * 5 * t), rtol=1e-14)
    assert abs(pair['water_reference_fid'][0]) == pytest.approx(amplitude)
    assert pair['time_s'][0] == 0


def test_frequency_origin_against_real_basis(independent):
    basis, cf = independent
    frequencies = np.fft.fftfreq(2048, .0005)
    for name, ppm in [('Ref0ppm', 0), ('H2O', 4.7)]:
        measured = frequencies[np.argmax(abs(np.fft.fft(basis[name])))]
        assert abs(measured - (ppm - 4.68) * cf) < 2000 / 2048
    assert generate_pair({'NAA': 10})['simulation_config']['water_frequency_offset_hz'] == pytest.approx(.02 * cf)


def test_paired_export_example(tmp_path):
    pair = generate_pair({'NAA': 10}, water_linewidth_hz=5, water_phase_rad=.6)
    directory = tmp_path / 'quantitative_pair'
    paths = export_pair(pair, directory)
    errors = []
    for key, path in paths.items():
        reader = mrs_io.read_FID(str(path))
        actual = np.asarray(reader.mrs().FID)
        original = pair[key]
        error = np.max(abs(actual - original))
        phase_error = np.max(abs(np.angle(actual * original.conj())))
        errors.append(error)
        np.testing.assert_allclose(actual, original, rtol=2e-7, atol=2e-7)
        np.testing.assert_allclose(np.asarray(nib.load(str(path)).dataobj).squeeze(),
                                   original.conj(), rtol=2e-7, atol=2e-7)
        assert phase_error < 2e-7
        assert reader.dwelltime == pytest.approx(.0005)
        assert reader.spectrometer_frequency[0] == pytest.approx(127.731)
        assert reader.SpecFreqChemShift == pytest.approx(4.68)
        assert actual.shape == (2048,)
        print(f'{key}: max_complex_error={error:.17g}; max_phase_error_rad={phase_error:.17g}')
    assert json.loads((directory / 'ground_truth.json').read_text()) == pair['ground_truth']
    assert json.loads((directory / 'simulation_config.json').read_text()) == pair['simulation_config']
    with pytest.raises(FileExistsError):
        export_pair(pair, directory)


@pytest.mark.parametrize('field,value', [('reference_identity', 'wrong'), ('bandwidth_hz', 1000), ('points', 1024)])
def test_calibration_incompatibility(field, value):
    _, _, frozen = format_frozen_basis(REFERENCE_BASIS, bandwidth_hz=2000, points=2048)
    with pytest.raises(ValueError, match='metadata mismatch'):
        generate_pair({'NAA': 10}, calibration=replace(frozen, **{field: value}))


def test_synthetic_metadata_rejected():
    with pytest.raises(ValueError, match='metadata mismatch'):
        generate_pair({'NAA': 10}, expected_metadata={'convention': 'wrong'})


@pytest.mark.parametrize('kwargs', [{'water_linewidth_hz': -1}, {'water_phase_rad': np.nan},
                                  {'water_frequency_offset_hz': 1000}, {'points': 1.5},
                                  {'water_concentration_mol_l': -1}, {'receiver_factor': np.inf}])
def test_invalid_configuration(kwargs):
    with pytest.raises(ValueError):
        generate_pair({'NAA': 10}, **kwargs)


def test_no_legacy_residual_configuration_or_rng_dependency():
    import torch
    from types import SimpleNamespace
    from simulation.signalModel import SignalModel
    model = SignalModel(None, SimpleNamespace(B0=3, TE=30, t=np.arange(32) / 2000, cf=127.731),
                        labels=['WM'], simulation_params=dict(water_amp_factor=5,
                            water_damping_min=1, water_damping_max=8,
                            water_phase_min=-180, water_phase_max=180), verbose=False)
    saved = np.random.get_state()
    torch_state = torch.random.get_rng_state().clone()
    try:
        np.random.seed(17)
        baseline = model.simulate_water().clone()
        next_random = np.random.random(8)
        np.random.seed(17)
        one = generate_pair({'NAA': 10})
        np.testing.assert_array_equal(model.simulate_water().numpy(), baseline.numpy())
        np.testing.assert_array_equal(np.random.random(8), next_random)
        model.simulation_params['water_amp_factor'] = 100
        two = generate_pair({'NAA': 10})
        np.testing.assert_array_equal(one['water_reference_fid'], two['water_reference_fid'])
        assert torch.equal(torch.random.get_rng_state(), torch_state)
    finally:
        np.random.set_state(saved)

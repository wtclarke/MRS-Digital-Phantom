from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from fsl_mrs.utils.mrs_io.main import read_basis

from simulation.frozen_basis_normalization import format_frozen_basis
from simulation.mrscloud_water_calibration import calibrated_responses, concentration_signals

PATH = Path('data/basissets/LCModel_Universal_Philips_UnEdited_PRESS_TE30.BASIS')
GRID = dict(bandwidth_hz=2000, points=2048)


@pytest.fixture(scope='module')
def calibrated():
    return calibrated_responses(PATH, selected_names=['NAA'], **GRID)


def test_independent_ratio_and_recovery(calibrated):
    responses, water, metadata = calibrated
    raw = read_basis(str(PATH))
    all_fids = raw.get_formatted_basis(2000, 2048)
    scalar = np.mean(abs(all_fids))
    ref = abs(all_fids[0, raw.names.index('Ref0ppm')]) / scalar
    naa = all_fids[:, raw.names.index('NAA')] / scalar
    signals, water_signal = concentration_signals(responses, water,
        concentrations_mol_l={'NAA': .010}, water_concentration_mol_l=40)
    expected_ratio = .010 * naa / (80 * ref)
    np.testing.assert_allclose(signals['NAA'] / water_signal, expected_ratio, rtol=2e-15)
    recovered = (signals['NAA'][0] / water_signal[0]) * 40 * (2 * ref) / naa[0]
    assert recovered == pytest.approx(.010, abs=1e-17)
    print(f'proton_response={ref:.17g}; NAA_initial={naa[0]!r}; '
          f'water_40M={water_signal[0]!r}; ratio={expected_ratio[0]!r}; '
          f'recovery_mM={recovered * 1000!r}')
    assert metadata['physical_calibration'] is False


def test_reproducibility_selection_and_frozen_compatibility(calibrated):
    responses, water, metadata = calibrated
    _, _, frozen = format_frozen_basis(PATH, **GRID)
    for selection in [['NAA'], ['Cr', 'NAA'], ['NAA', 'Cr', 'GPC', 'PCh'], read_basis(str(PATH)).names]:
        other, other_water, other_meta = calibrated_responses(
            PATH, selected_names=selection, calibration=frozen, expected_metadata=metadata, **GRID)
        np.testing.assert_array_equal(other['NAA'], responses['NAA'])
        np.testing.assert_array_equal(other_water, water)
        assert other_meta == metadata
    frozen_fids, _, _ = format_frozen_basis(PATH, selected_names=['NAA'], **GRID)
    np.testing.assert_array_equal(responses['NAA'], frozen_fids[:, 0])


def test_linearity_and_shared_receiver(calibrated):
    responses, water, _ = calibrated
    def generate(cm, cw, receiver=1):
        return concentration_signals(responses, water, concentrations_mol_l={'NAA': cm},
                                     water_concentration_mol_l=cw, receiver_factor=receiver)
    a, w = generate(.010, 40)
    b, _ = generate(.020, 40)
    _, v = generate(.010, 80)
    c, z = generate(.010, 40, 3)
    np.testing.assert_array_equal(b['NAA'], 2 * a['NAA'])
    np.testing.assert_array_equal(v, 2 * w)
    np.testing.assert_array_equal(c['NAA'], 3 * a['NAA'])
    np.testing.assert_array_equal(z, 3 * w)
    np.testing.assert_allclose(c['NAA'] / z, a['NAA'] / w, rtol=2e-15)


@pytest.mark.parametrize('field,value', [('reference_identity', 'wrong'), ('bandwidth_hz', 1000), ('points', 1024)])
def test_incompatible_frozen_metadata(field, value):
    _, _, frozen = format_frozen_basis(PATH, **GRID)
    with pytest.raises(ValueError, match='metadata mismatch'):
        calibrated_responses(PATH, selected_names=['NAA'], calibration=replace(frozen, **{field: value}), **GRID)


def test_incompatible_synthetic_metadata(calibrated):
    with pytest.raises(ValueError, match='metadata mismatch'):
        calibrated_responses(PATH, selected_names=['NAA'], expected_metadata={**calibrated[2], 'water_protons': 3}, **GRID)


def test_legacy_residual_water_and_rng_unchanged(calibrated):
    from simulation.signalModel import SignalModel
    basis = SimpleNamespace(B0=3, TE=30, t=np.arange(32) / 2000, cf=123.2)
    params = dict(water_amp_factor=5, water_damping_min=1, water_damping_max=8,
                  water_phase_min=-180, water_phase_max=180)
    model = SignalModel(None, basis, labels=['WM', 'GM', 'CSF'], simulation_params=params, verbose=False)
    saved = np.random.get_state()
    try:
        np.random.seed(17)
        baseline = model.simulate_water().numpy().copy()
        next_baseline = np.random.random(8)
        np.random.seed(17)
        calibrated_responses(PATH, selected_names=['NAA'], **GRID)
        current = model.simulate_water().numpy()
        np.testing.assert_array_equal(current, baseline)
        np.testing.assert_array_equal(np.random.random(8), next_baseline)
        params['water_amp_factor'] = 100
        again = calibrated_responses(PATH, selected_names=['NAA'], **GRID)
        np.testing.assert_array_equal(again[1], calibrated[1])
    finally:
        np.random.set_state(saved)

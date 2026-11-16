from pathlib import Path
from dataclasses import replace

import numpy as np
import pytest
from fsl_mrs.utils.mrs_io.main import read_basis

from simulation.frozen_basis_normalization import format_frozen_basis

BASIS = Path("data/basissets/LCModel_Universal_Philips_UnEdited_PRESS_TE30.BASIS")
BW, POINTS = 2000, 2048


def _legacy(path, selected):
    b = read_basis(str(path))
    names = list(b.names)
    ignore = [name for name in names if name not in selected]
    arr = b.get_formatted_basis(BW, POINTS, ignore=ignore)
    scalar = float(np.mean(np.abs(arr)))
    return arr / scalar, [n for n in names if n in selected], scalar


def test_legacy_selection_experiment_and_frozen_selection_invariance():
    selections = [["NAA"], ["NAA", "Cr"], ["NAA", "Cr", "GPC", "PCh"], None]
    legacy_scalars = []
    frozen_fids = []
    calibration = None
    for selected in selections:
        legacy_names = read_basis(str(BASIS)).names
        requested = legacy_names if selected is None else selected
        legacy_arr, legacy_order, scalar = _legacy(BASIS, requested)
        legacy_scalars.append(scalar)
        legacy_naa = legacy_arr[:, legacy_order.index("NAA")]
        print({"selection": selected or "complete", "legacy_scalar": scalar,
               "legacy_NAA_initial": legacy_naa[0]})
        if calibration is None:
            frozen, _, calibration = format_frozen_basis(
                BASIS, bandwidth_hz=BW, points=POINTS, selected_names=["NAA"]
            )
        else:
            frozen, _, same = format_frozen_basis(
                BASIS, bandwidth_hz=BW, points=POINTS, selected_names=["NAA"], calibration=calibration
            )
            assert same == calibration
        frozen_fids.append(frozen[:, 0])
        assert legacy_naa.shape == frozen[:, 0].shape
    assert len(set(legacy_scalars)) > 1
    assert max(np.max(np.abs(x - frozen_fids[0])) for x in frozen_fids) == 0


def test_frozen_basis_concentration_linearity_and_reproducibility():
    basis, _, cal = format_frozen_basis(BASIS, bandwidth_hz=BW, points=POINTS, selected_names=["NAA"])
    signal_10mM = basis[:, 0] * 0.010
    signal_20mM = basis[:, 0] * 0.020
    np.testing.assert_array_equal(signal_20mM, 2 * signal_10mM)
    again, _, cal_again = format_frozen_basis(BASIS, bandwidth_hz=BW, points=POINTS, selected_names=["NAA"])
    np.testing.assert_array_equal(basis, again)
    assert cal == cal_again


@pytest.mark.parametrize("field,value", [("reference_identity", "wrong"), ("bandwidth_hz", 1000), ("points", 1024)])
def test_incompatible_frozen_metadata_rejected(field, value):
    _, _, calibration = format_frozen_basis(BASIS, bandwidth_hz=BW, points=POINTS, selected_names=["NAA"])
    incompatible = replace(calibration, **{field: value})
    with pytest.raises(ValueError, match="metadata mismatch"):
        format_frozen_basis(BASIS, bandwidth_hz=BW, points=POINTS,
                            selected_names=["NAA"], calibration=incompatible)

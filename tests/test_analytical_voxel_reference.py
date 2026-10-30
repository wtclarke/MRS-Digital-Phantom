import numpy as np
import pytest
from fsl_mrs.utils import mrs_io

from simulation.analytical_voxel_reference import simulate_analytical_voxel_pair, export_analytical_pair


@pytest.mark.parametrize("fractions,expected", [({"WM": 1}, 39.0165), ({"GM": 1}, 47.2305), ({"CSF": 1}, 55.5)])
def test_pure_tissue_water_amplitude(fractions, expected):
    pair = simulate_analytical_voxel_pair(fractions, metabolite_concentrations_mM={k: 10 for k in fractions})
    assert pair["water_reference_fid"][0] == pytest.approx(2 * expected)


def test_mixed_weighting_ratio_and_determinism():
    fractions = {"WM": .2, "GM": .3, "CSF": .5}
    kwargs = dict(tissue_fractions=fractions, metabolite_concentrations_mM={k: 10 for k in fractions})
    pair = simulate_analytical_voxel_pair(**kwargs)
    repeated = simulate_analytical_voxel_pair(**kwargs)
    expected_w = 2 * (.2 * 39.0165 + .3 * 47.2305 + .5 * 55.5)
    expected_m = 3 * .010
    np.testing.assert_array_equal(pair["water_reference_fid"], repeated["water_reference_fid"])
    assert pair["water_reference_fid"][0] == pytest.approx(expected_w)
    assert pair["metabolite_fid"][0] == pytest.approx(expected_m)
    assert pair["metabolite_fid"][0] / pair["water_reference_fid"][0] == pytest.approx(expected_m / expected_w)


def test_changed_tissue_fraction_changes_water_as_weighted_sum():
    wm = simulate_analytical_voxel_pair({"WM": 1}, metabolite_concentrations_mM={"WM": 10})
    gm = simulate_analytical_voxel_pair({"GM": 1}, metabolite_concentrations_mM={"GM": 10})
    assert wm["water_reference_fid"][0] != gm["water_reference_fid"][0]


def test_paired_nifti_export_roundtrip(tmp_path):
    pair = simulate_analytical_voxel_pair({"WM": .5, "CSF": .5}, metabolite_concentrations_mM={"WM": 10, "CSF": 10})
    paths = export_analytical_pair(pair, tmp_path, spectrometer_frequency_mhz=123.2)
    for key, path in (("metabolite_fid", paths["metabolite_fid"]), ("water_reference_fid", paths["water_reference_fid"])):
        loaded = mrs_io.read_FID(path)
        actual = np.asarray(loaded.mrs().FID)
        np.testing.assert_allclose(actual, pair[key], rtol=2e-7, atol=2e-7)
        assert loaded.dwelltime == pytest.approx(pair["metadata"]["DwellTime_s"])
        assert loaded.spectrometer_frequency[0] == pytest.approx(123.2)
    metadata = __import__("json").loads((tmp_path / "ground_truth.json").read_text())
    assert metadata["CalibrationIdentifier"] == "synthetic-protons-per-mol-v1"

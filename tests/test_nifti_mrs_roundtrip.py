import numpy as np
import nibabel as nib
import pytest
import torch
from fsl_mrs.utils import mrs_io

from simulation.headless import results_to_nifti_mrs
from utils.auxillary import fid2spec
from utils.nii_processing import save_nifti_mrs

N = 256
DWELL = 1 / 2000
SFREQ = 123.2
OFFSET = 37.0


def fid():
    t = np.arange(N) * DWELL
    return np.exp(-t / 0.085) * np.exp(1j * (np.deg2rad(35) + 2 * np.pi * OFFSET * t))


def save_utils(tmp_path, spectrum):
    save_nifti_mrs(str(tmp_path), ["signal"], [], {"signal": spectrum}, np.eye(4), DWELL,
                   SFREQ, base_name="utils")
    return tmp_path / "utils" / "total.nii.gz"


def save_headless(tmp_path, spectrum):
    metadata = {"dwelltime": DWELL, "spec_freq": SFREQ, "nucleus": "1H", "affine": np.eye(4),
                "sim_params": {"basis_set_settings": {"TE": 0, "vendor": "test"},
                               "simulation_params": {"TR": 0}}}
    result = {"position": (0, 0, 0), "components": {"signal": torch.as_tensor(spectrum)},
              "spectrum": torch.as_tensor(spectrum), "nifti_metadata": metadata}
    obj = results_to_nifti_mrs([result], total_components=["signal"])[0]["total"]
    path = tmp_path / "headless.nii.gz"
    obj.save(str(path))
    return path


def inspect_roundtrip(path, original, expected_storage):
    img = nib.load(str(path))
    stored = np.asarray(img.dataobj).reshape(-1)
    assert np.allclose(stored, expected_storage, rtol=2e-7, atol=2e-7)
    interpreted = np.asarray(mrs_io.read_FID(str(path)).mrs().FID).reshape(-1)
    error = np.linalg.norm(interpreted - original) / np.linalg.norm(original)
    amp_error = np.max(np.abs(np.abs(interpreted) - np.abs(original)))
    phase_error = np.max(np.abs(np.angle(interpreted * np.conj(original))))
    conjugate_error = np.linalg.norm(interpreted - np.conj(original)) / np.linalg.norm(original)
    assert error < 2e-7
    assert amp_error < 2e-7
    assert phase_error < 2e-7
    assert conjugate_error > 0.1
    assert img.header["pixdim"][4] == pytest.approx(DWELL)
    reader = mrs_io.read_FID(str(path))
    assert reader.dwelltime == pytest.approx(DWELL)
    assert reader.spectrometer_frequency[0] == pytest.approx(SFREQ)
    assert np.angle(interpreted[1] * np.conj(interpreted[0])) == pytest.approx(
        2 * np.pi * OFFSET * DWELL, abs=2e-6)
    return error, amp_error, phase_error, conjugate_error


@pytest.mark.parametrize("exporter", ["utils", "headless"])
def test_one_dimensional_export_roundtrip(tmp_path, exporter):
    original = fid()
    spectrum = fid2spec(original)
    if exporter == "utils":
        path = save_utils(tmp_path, spectrum)
        expected = np.conj(original).astype(np.complex64)
    else:
        path = save_headless(tmp_path, spectrum)
        expected = np.conj(original).astype(np.complex64)
    errors = inspect_roundtrip(path, original, expected)
    print(f"{exporter} errors relative/amplitude/phase/conjugate={errors}")


def test_headless_image_shape_not_supported_by_current_interface(tmp_path):
    # results_to_nifti_mrs adds a time dimension to this 3D spatial spectrum,
    # producing an unsupported 8D array in the current NIfTI-MRS generator.
    original = np.broadcast_to(fid(), (2, 2, 2, N)).copy()
    spectrum = fid2spec(original)
    with pytest.raises(ValueError, match="dimensions kwarg"):
        save_headless(tmp_path, spectrum)


def test_utils_image_shape_is_not_expanded_as_spatial_data(tmp_path):
    original = np.broadcast_to(fid(), (2, 2, 2, N)).copy()
    path = save_utils(tmp_path, fid2spec(original))
    stored = np.asarray(nib.load(str(path)).dataobj)
    assert stored.shape == (2, 2, 2, N)

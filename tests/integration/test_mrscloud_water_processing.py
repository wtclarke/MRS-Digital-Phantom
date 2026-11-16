import hashlib
from pathlib import Path

import numpy as np
import pytest
from fsl_mrs.core import MRS
from fsl_mrs.utils.misc import FIDToSpec, ts_to_ts_ft

from examples.mrscloud_water_processing_validation import capture_raw, second_dc, sample_on_acquisition_grid
from examples.fsl_mrs_water_fit_diagnostics import fit_water


@pytest.fixture(scope='module')
def captured(tmp_path_factory):
    source=Path('simulation/basissets/MRSCloud-main/functions/sim_signals.m')
    before=hashlib.sha256(source.read_bytes()).hexdigest()
    result=capture_raw(tmp_path_factory.mktemp('capture_processing'),spatial_points=3,points=256)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==before
    return result


def test_capture_raw_exponential_and_independent_first_dc(captured):
    raw=captured['raw_fid'];after=captured['after_first_dc']
    t=np.arange(len(raw))/float(captured['sw'])
    np.testing.assert_allclose(raw,raw[0]*np.exp(-np.pi*float(captured['lw'])*t),rtol=1e-12,atol=1e-13)
    spectrum=np.fft.fftshift(np.fft.ifft(raw))
    np.testing.assert_allclose(spectrum,captured['raw_spec'],rtol=1e-12,atol=1e-13)
    ppm=captured['ppm']
    a=(ppm>max(ppm)-.5)&(ppm<max(ppm));b=(ppm>min(ppm))&(ppm<min(ppm)+.5)
    offset=(np.mean(spectrum[a])+np.mean(spectrum[b]))/2
    np.testing.assert_allclose(after,np.fft.fft(np.fft.ifftshift(spectrum-offset)),rtol=1e-12,atol=1e-13)
    assert after[0]==pytest.approx(raw[0]-len(raw)*offset,rel=1e-12,abs=1e-13)
    np.testing.assert_allclose(raw[1:],after[1:],rtol=1e-12,atol=1e-13)
    assert abs(raw[0]-after[0])>.001*abs(raw[0])


def test_second_dc_only_changes_first_sample(captured):
    before=captured['after_first_dc']
    after,offset=second_dc(before,4000,127.731,4.7)
    assert after[0]==pytest.approx(before[0]-offset,rel=1e-12,abs=1e-13)
    np.testing.assert_allclose(after[1:],before[1:],rtol=1e-12,atol=1e-13)


def test_first_point_weighting_is_analysis_not_persistent_fid_edit(captured):
    raw=captured['raw_fid'].copy();saved=raw.copy()
    spectrum=FIDToSpec(raw)
    np.testing.assert_array_equal(raw,saved)
    expected=np.fft.fftshift(np.fft.fft(raw,norm='ortho'))-.5*raw[0]/np.sqrt(len(raw))
    np.testing.assert_allclose(spectrum,expected,rtol=1e-12,atol=1e-13)


def test_fourier_cropping_and_integer_acquisition_are_different(captured):
    raw=captured['raw_fid']
    sampled=sample_on_acquisition_grid(raw,.00025,.0005,128)
    np.testing.assert_array_equal(sampled,raw[np.arange(128)*2])
    assert sampled[0]==raw[0]
    actual=ts_to_ts_ft(raw[:,None],.00025,.0005,128).ravel()
    spec=np.fft.fftshift(np.fft.fft(raw))
    expected=np.fft.ifft(np.fft.ifftshift(spec[64:-64]))*.5
    np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-13)
    assert abs(actual[0]-raw[0])>1e-4*abs(raw[0])


def test_raw_acquisition_water_fit_is_consistent(captured):
    raw=captured['raw_fid'];native_t=np.arange(len(raw))/4000
    shifted=raw*np.exp(2j*np.pi*2.55462*native_t)
    sampled=sample_on_acquisition_grid(shifted,.00025,.0005,128)
    mrs=MRS(FID=np.ones(128,dtype=complex),H2O=40*sampled,cf=127.731,bw=2000,nucleus='1H',chemShift=4.68)
    mrs.rescaleForFitting()
    info,_,_=fit_water(mrs)
    assert info['success']
    assert info['relative_complex_residual']<1e-6
    assert info['fitted_amplitude']==pytest.approx(40*abs(raw[0]),rel=2e-6)


def test_sampling_rejects_fractional_or_uncovered_grid():
    with pytest.raises(ValueError):sample_on_acquisition_grid(np.ones(16),.001,.0015,8)
    with pytest.raises(ValueError):sample_on_acquisition_grid(np.ones(16),.001,.002,32)

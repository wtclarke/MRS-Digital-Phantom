import numpy as np
import pytest

from examples.mrscloud_measured_linewidth_validation import fwhm, absorption_spectrum, measure


@pytest.mark.parametrize('width',[1.,3.,5.,8.])
def test_analytical_lorentzian_and_known_constant_background(width):
    frequency=np.linspace(-20,20,80001)
    absorption=.13+3.7/(1+(2*frequency/width)**2)
    result=fwhm(frequency,absorption,baseline=.13)
    assert result['fwhm_hz']==pytest.approx(width,abs=1e-8)
    assert result['left_half_hz']==pytest.approx(-width/2,abs=1e-8)
    assert result['right_half_hz']==pytest.approx(width/2,abs=1e-8)


def test_long_complex_fids_measure_absorption_not_magnitude_width():
    time=np.arange(32768)*.0005
    for width in [1,3,5,8]:
        fid=(2*np.exp(.6j))*np.exp(-np.pi*width*time)*np.exp(2j*np.pi*37*time)
        saved=fid.copy()
        result,freq,real=measure(fid,37,.6,zero_fill=8)
        assert result['fwhm_hz']==pytest.approx(width,rel=5e-5)
        _,complex_spec=absorption_spectrum(fid,dwell_time_s=.0005,center_hz=37,phase_rad=.6,zero_fill=8)
        # Infinite-record Lorentzian magnitude width is sqrt(3) times the
        # absorption width; explicitly distinguish the two measurements.
        magnitude=fwhm(freq,abs(complex_spec))['fwhm_hz']
        assert magnitude/result['fwhm_hz']==pytest.approx(np.sqrt(3),rel=3e-4)
        np.testing.assert_array_equal(fid,saved)


def test_truncation_and_demodulation_are_distinct_from_inherent_linewidth():
    time=np.arange(2048)*.0005
    decay=np.exp(-np.pi*time).astype(complex)
    shifted=decay*np.exp(1j*(2*np.pi*37*time+.6))
    a,_,_=measure(decay,0,0)
    b,_,_=measure(shifted,37,.6)
    assert a['fwhm_hz']==pytest.approx(b['fwhm_hz'],abs=2e-13)
    assert 1.066<a['fwhm_hz']<1.068


def test_known_dc_pedestal_changes_half_zero_width_but_not_background_corrected_width():
    time=np.arange(8192)*.0005
    fid=np.exp(-np.pi*3*time).astype(complex)
    modified=fid.copy();modified[0]-=.5
    raw,freq,_=measure(fid,0,0)
    apparent,_,spec=measure(modified,0,0)
    assert apparent['fwhm_hz']<raw['fwhm_hz']
    corrected=fwhm(freq,spec,baseline=-.25/np.sqrt(len(fid)))
    assert corrected['fwhm_hz']==pytest.approx(raw['fwhm_hz'],abs=2e-13)
    with pytest.raises(ValueError):measure(np.ones(16),0,0)
    with pytest.raises(ValueError):fwhm(np.linspace(-1,1,10),np.zeros(10))

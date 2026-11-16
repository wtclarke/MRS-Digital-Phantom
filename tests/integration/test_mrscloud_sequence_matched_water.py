from copy import deepcopy

import numpy as np
import pytest

from examples.mrscloud_sequence_matched_water import simulate_raw, assemble_common, write_basis
from fsl_mrs.utils import mrs_io


@pytest.fixture(scope='module')
def raw(tmp_path_factory):
    return simulate_raw(tmp_path_factory.mktemp('sequence_water_octave'),spatial_points=3,points=256)


def test_deterministic_water_and_two_proton_response(raw):
    np.testing.assert_array_equal(raw['signals']['water_on'],raw['signals']['water_repeat'])
    np.testing.assert_array_equal(raw['signals']['water_on'],2*raw['signals']['water_single_proton'])


def test_sequence_dc_processing_alters_first_point(raw):
    water=raw['signals']['water_on']
    t=np.arange(int(raw['points']))/float(raw['sw'])
    extrapolated=water[1]*np.exp(np.pi*float(raw['lw'])/float(raw['sw']))
    np.testing.assert_allclose(water[1:],extrapolated*np.exp(-np.pi*float(raw['lw'])*t[1:]),rtol=1e-12,atol=1e-13)
    before_dc=extrapolated*np.exp(-np.pi*float(raw['lw'])*t)
    spec=np.fft.fftshift(np.fft.ifft(before_dc))
    freq=np.arange(int(raw['points']))*float(raw['sw'])/int(raw['points'])-float(raw['sw'])/2+float(raw['sw'])/(2*int(raw['points']))
    ppm=-freq/127.731+4.65-(4.68-4.7)
    a=(ppm>max(ppm)-.5)&(ppm<max(ppm))
    b=(ppm>min(ppm))&(ppm<min(ppm)+.5)
    dc=(np.mean(spec[a])+np.mean(spec[b]))/2
    expected=np.fft.fft(np.fft.ifftshift(spec-dc))
    np.testing.assert_allclose(water,expected,rtol=1e-12,atol=1e-13)
    # Preserve this processing limitation; do not count sample zero as an
    # undistorted physical water amplitude or compensate it by a new factor.
    assert abs(water[0]-extrapolated)>1e-3*abs(extrapolated)


def test_frequency_placement_and_carrier_effect(raw):
    signals,_,_,_=assemble_common(raw)
    water=signals['water_on']
    # Beyond sample zero (changed by DC correction), native on-carrier water
    # is an exact single exponential; post-acquisition shift sets its origin.
    f=np.angle(water[2]/water[1])*float(raw['sw'])/(2*np.pi)
    assert f==pytest.approx((4.7-4.68)*127.731,abs=1e-10)
    assert not np.allclose(abs(signals['water_on']),abs(signals['water_off']),rtol=1e-5,atol=0)
    assert raw['TE1']+raw['TE2']==30
    assert raw['refTp']==7


def test_one_common_normalization_preserves_relative_amplitudes(raw,tmp_path):
    signals,scalar,carriers,processed=assemble_common(raw)
    assert scalar>0
    t=np.arange(int(raw['points']))/float(raw['sw'])
    for name,signal in signals.items():
        expected=processed[name]/scalar*np.exp(-2j*np.pi*(4.68-carriers[name])*127.731*t)
        np.testing.assert_array_equal(signal,expected)
    changed=deepcopy(raw)
    changed['signals']['water_on']*=3
    again,same,_,_=assemble_common(changed)
    assert same==scalar
    np.testing.assert_array_equal(again['NAA'],signals['NAA'])
    np.testing.assert_allclose(again['water_on'],3*signals['water_on'],rtol=2e-13,atol=1e-16)
    path=tmp_path/'response.BASIS'
    write_basis(path,signals,raw)
    reloaded=mrs_io.read_basis(str(path))
    for name in signals:
        actual=reloaded.original_basis_array[:,reloaded.names.index(name)]
        assert np.linalg.norm(actual-signals[name])/np.linalg.norm(signals[name])<3e-7

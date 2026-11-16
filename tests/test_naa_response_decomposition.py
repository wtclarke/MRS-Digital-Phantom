import numpy as np
import pytest

from examples.mrscloud_naa_response_decomposition import first_dc, stages, decompose, ideal_grid_controls


def test_first_dc_changes_only_sample_zero_with_independent_edge_estimate():
    t=np.arange(256)/4000
    raw=(2+.3j)*np.exp(-np.pi*t)*np.exp(-2j*np.pi*126.727*t)
    ppm=np.linspace(20,-10,256)
    spectrum=np.fft.fftshift(np.fft.ifft(raw))
    edges=[spectrum[(ppm>max(ppm)-.5)&(ppm<max(ppm))],
           spectrum[(ppm>min(ppm))&(ppm<min(ppm)+.5)]]
    offset=sum(np.mean(x) for x in edges)/2
    after=first_dc(raw,ppm)
    assert after[0]==pytest.approx(raw[0]-256*offset,rel=1e-12)
    np.testing.assert_allclose(after[1:],raw[1:],rtol=1e-12,atol=1e-13)


def test_frame_change_exact_sampling_and_shared_scalar():
    raw=np.exp((-np.pi-2j*np.pi*126.727)*np.arange(256)/4000)
    first=raw.copy();first[0]*=.5
    native,target=stages(raw,first,carrier=3.,scalar=7.,points=128)
    _,doubled=stages(raw,first,carrier=3.,scalar=14.,points=128)
    expected=raw/7*np.exp(2j*np.pi*(3.-4.68)*127.731*np.arange(256)/4000)
    np.testing.assert_array_equal(native['raw'],expected)
    np.testing.assert_array_equal(target['raw_target'],expected[::2])
    for key in target:
        np.testing.assert_array_equal(doubled[key],target[key]/2)
    np.testing.assert_allclose(abs(native['raw']),abs(raw)/7,rtol=1e-15)
    assert abs(target['raw_fourier'][0]-target['raw_target'][0])>1e-4


def test_independently_specified_ideal_control_has_no_proton_or_phase_deficit():
    t=np.arange(2048)*.0005
    water=2*np.exp(-np.pi*t)*np.exp(2j*np.pi*2.55462*t)
    freq=(2.008-3.)*127.728+(3.-4.68)*127.731
    methyl=3*np.exp(-np.pi*t)*np.exp(2j*np.pi*freq*t)
    keys=['methyl_raw_target','full_raw_target','full_first_dc_target','full_second_dc_target','full_formatted']
    result=decompose(dict.fromkeys(keys,methyl),water)
    for key,value in result['factors'].items():
        if key!='ideal_finite_grid_window':assert value==pytest.approx(1,rel=2e-14)
    assert result['product']==pytest.approx(result['direct_alpha'],rel=2e-15)
    assert result['factors']['ideal_finite_grid_window']==pytest.approx(.9098704420596527,rel=2e-13)


def test_longer_ideal_records_separate_truncation_from_sequence_response():
    rows=ideal_grid_controls()
    assert rows==ideal_grid_controls()
    assert rows[0]['ideal_response_factor']==pytest.approx(.9098704420596527,rel=2e-13)
    assert rows[-1]['ideal_response_factor']==pytest.approx(.9882461104085961,rel=2e-13)
    assert abs(rows[-1]['ideal_response_factor']-rows[-2]['ideal_response_factor'])<2e-5

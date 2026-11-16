import numpy as np
import pytest

from simulation.water_reference_preprocessing import align_water_reference


@pytest.mark.parametrize('frequency,phase', [(0.,0.),(2.55462,0.),(0.,.6),(2.55462,.6)])
def test_alignment_magnitudes_metadata_and_determinism(frequency, phase):
    t = np.arange(2048)/2000
    amplitude = 61.73993435001698
    original = amplitude*np.exp(-np.pi*5*t+1j*(2*np.pi*frequency*t+phase))
    saved = original.copy()
    kwargs=dict(dwell_time_s=.0005,frequency_offset_hz=frequency,initial_phase_rad=phase)
    aligned, metadata=align_water_reference(original,**kwargs)
    again, repeated_metadata=align_water_reference(original,**kwargs)
    np.testing.assert_array_equal(original,saved)
    assert not np.shares_memory(original,aligned)
    np.testing.assert_allclose(aligned,amplitude*np.exp(-np.pi*5*t),rtol=1e-13,atol=1e-13)
    np.testing.assert_allclose(abs(aligned),abs(original),rtol=5e-16,atol=3e-14)
    assert abs(aligned[0]) == pytest.approx(amplitude,rel=2e-15)
    np.testing.assert_array_equal(aligned,again)
    assert metadata == repeated_metadata
    assert metadata['applied_frequency_correction_hz'] == -frequency
    assert metadata['applied_phase_correction_rad'] == -phase
    assert metadata['magnitude_scaling'] == 1


@pytest.mark.parametrize('kwargs', [dict(dwell_time_s=0),dict(initial_phase_rad=np.nan),
                                  dict(frequency_offset_hz=1000)])
def test_invalid_corrections(kwargs):
    defaults=dict(dwell_time_s=.0005,frequency_offset_hz=0,initial_phase_rad=0)
    defaults.update(kwargs)
    with pytest.raises(ValueError):
        align_water_reference(np.ones(8,dtype=complex),**defaults)

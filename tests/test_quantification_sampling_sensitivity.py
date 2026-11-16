import numpy as np
import pytest

from examples.mrscloud_quantification_sampling_sensitivity import broaden, predict, POINTS, LINEWIDTHS


def analytical_master():
    time=np.arange(8192)/2000
    methyl=3*np.exp(-np.pi*time)*np.exp(-2j*np.pi*341.294256*time)
    water=2*np.exp(-np.pi*time)*np.exp(2j*np.pi*2.55462*time)
    return {'full_formatted':methyl,'water_raw':water}


def test_common_broadening_preserves_initial_samples_and_input_arrays():
    master=analytical_master();saved={k:v.copy() for k,v in master.items()}
    for points in POINTS:
        for width in LINEWIDTHS:
            for key,raw in master.items():
                actual=broaden(raw,points,width)
                assert actual[0]==raw[0]
                expected=raw[:points]*np.exp(-np.pi*(width-1)*np.arange(points)/2000)
                np.testing.assert_array_equal(actual,expected)
                np.testing.assert_array_equal(actual,broaden(raw,points,width))
                np.testing.assert_array_equal(raw,saved[key])
    with pytest.raises(ValueError):broaden(master['water_raw'],16384,1)
    with pytest.raises(ValueError):broaden(master['water_raw'],2048,0)


def test_prediction_independent_of_final_fsl_quantification_and_common_gain(monkeypatch):
    from fsl_mrs.utils import quantify
    from fsl_mrs.utils.results import FitRes
    def forbidden(*args,**kwargs):
        raise AssertionError('independent prediction called fitting/reporting')
    monkeypatch.setattr(quantify,'quantifyWater',forbidden)
    monkeypatch.setattr(FitRes,'getConc',forbidden)
    master=analytical_master()
    for n in POINTS:
        for width in LINEWIDTHS:
            a=predict(master,n,width);b=predict({k:2*v for k,v in master.items()},n,width)
            assert a['response_factor']==b['response_factor']
            assert b['unit_metabolite_integral']==2*a['unit_metabolite_integral']
            assert b['unit_water_integral']==2*a['unit_water_integral']
            assert a==predict(master,n,width)
    assert predict(master,2048,1)['response_factor']==pytest.approx(.9098704420596527,rel=2e-13)

import numpy as np
import pytest

from examples.fsl_mrs_sequence_water_quantification import independent_concentration, quantification_info


def test_public_quantification_overrides_preserve_acquisition_and_defaults():
    for water, relax in [(False,False),(True,False),(False,True),(True,True)]:
        info = quantification_info(127.731, match_water=water, match_relaxation=relax)
        assert info.te == .030 and info.tr == 2.
        assert info.ref_protons == 3 and info.ref_limits == (1.8,2.2)
        assert info.h2o_limits == (1.65,7.65) and info.add_corr == 1.
        assert info.d_WM == pytest.approx(40000/55010 if water else .65)
        if relax:
            assert info.R_H2O_WM == 1. and info.R_M == 1.
        else:
            assert info.R_H2O_WM == pytest.approx(.5786670548409867)
            assert info.R_M == pytest.approx(.6749567292246237)
    # Overrides must not mutate package defaults for subsequent runs.
    default = quantification_info(127.731)
    assert default.d_WM == .65 and default.t1['METAB'] == 1.29


def test_independent_equation_calls_no_final_fsl_quantification(monkeypatch):
    from fsl_mrs.utils import quantify
    from fsl_mrs.utils.results import FitRes
    def forbidden(*args, **kwargs):
        raise AssertionError('final FSL reporting called by independent reconstruction')
    monkeypatch.setattr(quantify, 'quantifyWater', forbidden)
    monkeypatch.setattr(FitRes, 'getConc', forbidden)
    # Constant real spectrum from an impulse: known trapz interval counts.
    m = np.zeros(8,complex);m[0] = 3
    w = np.zeros(8,complex);w[0] = 2
    ppm = np.linspace(1.8,2.2,8)
    settings = dict(TE_s=.030, TR_s=2., water_T1_s=1e-6, water_T2_s=np.inf,
                    metabolite_T1_s=1e-6, metabolite_T2_s=np.inf, effective_water_mM=40000)
    result = independent_concentration(m,w,ppm,settings)
    assert result['independent_mM'] == pytest.approx(40000)
    assert result['metabolite_integral'] == pytest.approx(7*1.5/np.sqrt(8))
    assert result['water_integral'] == pytest.approx(7/np.sqrt(8))
    assert m[0] == 3 and w[0] == 2
    assert independent_concentration(m,w,ppm,settings) == result

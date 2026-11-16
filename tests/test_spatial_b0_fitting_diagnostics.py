"""Fast independent numerical controls; actual optimizers live in integration tests."""
import json

import numpy as np
import pytest

from examples.mrscloud_spatial_b0_fitting_validation import (
    SOURCE, NAMES, independent_signals, nuisance_design, projection, spectrum)


def test_fixed_complex_projection_recovers_known_real_coefficients():
    t=np.arange(2048)/2000
    columns=np.column_stack([np.exp(-np.pi*t)*np.exp(-2j*np.pi*f*t) for f in [300.,350.]])
    expected=np.array([.01,.004])
    actual,resid=projection(spectrum(columns),columns@expected,(.2,4.2))
    np.testing.assert_allclose(actual,expected,rtol=5e-15)
    assert resid<5e-15


def test_independent_spatial_response_uses_no_fitted_amplitudes(monkeypatch):
    from fsl_mrs.utils import fitting
    from fsl_mrs.utils.results import FitRes
    def forbidden(*args,**kwargs):raise AssertionError('fitting/reporting called by source prediction')
    monkeypatch.setattr(fitting,'fit_FSLModel',forbidden)
    monkeypatch.setattr(FitRes,'getConc',forbidden)
    truth=json.loads((SOURCE/'symmetric_B0/ground_truth.json').read_text())
    config=json.loads((SOURCE/'symmetric_B0/simulation_config.json').read_text())
    t=np.arange(2048)/2000
    responses={n:(i+1)*np.exp(-np.pi*t-2j*np.pi*(300+10*i)*t) for i,n in enumerate(NAMES)}
    water=np.exp(-np.pi*t+2j*np.pi*2.55462*t)
    components,w,effective,cw,centroid=independent_signals(truth,config,responses,water)
    for n in NAMES:
        c=truth['tissue_specific_metabolite_concentrations_mM']['WM'][n]
        np.testing.assert_allclose(components[n],.001*c*responses[n]*np.cos(2*np.pi*5*t),rtol=2e-15,atol=1e-17)
        assert effective[n]==c
    np.testing.assert_allclose(w,39.0165*water*np.cos(2*np.pi*5*t),rtol=2e-15,atol=1e-14)
    assert cw==39.0165 and centroid==pytest.approx(2.55462,abs=1e-12)
    # Removing a carrier cannot remove the sign-changing cosine envelope.
    aligned=w*np.exp(-2j*np.pi*centroid*t)
    assert aligned[0].real>0 and aligned[200].real<0


def test_nuisance_design_frequency_direction_and_phase():
    t=np.arange(2048)/2000;nu=np.fft.fftshift(np.fft.fftfreq(2048,.0005))
    response={n:np.exp(-np.pi*t-2j*np.pi*300*t) for n in NAMES}
    p=dict(gamma_0=.3,eps_0=-2*np.pi*37,Phi0=.2,Phi1=.001)
    expected=spectrum(response['NAA']*np.exp(-.3*t+2j*np.pi*37*t))*np.exp(-1j*(.2+.001*nu))
    np.testing.assert_allclose(nuisance_design(response,p)[:,0],expected,rtol=1e-13,atol=1e-14)

"""Numerical diagnostic controls, independent of proprietary fitting execution."""
import numpy as np
import pytest

from examples.lcmodel_fit_bias_diagnostics import project,water_recurrence_spectrum


def test_independent_projection_recovers_complex_mixture():
    t=np.arange(256)*.0005
    design=np.column_stack([np.exp(-np.pi*2*t+1j*(2*np.pi*f*t+.2)) for f in [37.,-119.]])
    truth=np.array([.01,.008]);result=project(design,design@truth)
    np.testing.assert_allclose(result['coefficients_real'],truth,rtol=1e-13)
    np.testing.assert_allclose(result['coefficients_imag'],0,atol=1e-17)
    assert result['relative_residual']<1e-14


def test_unacquired_basis_tail_is_a_real_model_mismatch():
    # A real spectral projection of full-duration versus finite-duration signals.
    t=np.arange(512)*.0005;basis=np.exp(-np.pi*t+1j*2*np.pi*37*t)
    acquired=basis.copy();acquired[256:]=0
    spectrum=np.fft.fft(basis)[:,None];data=np.fft.fft(acquired)
    full=project(spectrum,data,real=True)
    exact=project(np.fft.fft(acquired)[:,None],data,real=True)
    assert full['relative_residual']>.1
    assert exact['relative_residual']<1e-14
    assert abs(full['coefficients_real'][0]-1)>.01


def test_single_precision_water_transform_against_geometric_sum():
    n=256;dt=.0005
    z=np.exp(-np.pi*3*dt-1j*2*np.pi*37*dt)
    fid=2*np.exp(1j*.3)*z**np.arange(n)
    grid=z*np.exp(-1j*2*np.pi*np.arange(n)/n)
    expected=np.fft.fftshift(2*np.exp(1j*.3)*(1-grid**n)/(1-grid)/np.sqrt(2*n))
    actual=water_recurrence_spectrum(fid)
    # Complex multiply/add recurrence accumulates single-precision rounding over n steps.
    assert np.linalg.norm(actual-expected)/np.linalg.norm(expected)<4*n*np.finfo(np.float32).eps
    np.testing.assert_array_equal(actual,water_recurrence_spectrum(fid))


def test_normalization_and_fitted_coefficients_are_distinct():
    # Values in this arithmetic example are not used as a calibration.
    reference=.0015;water=.001;concentration_water_mM=40000
    scale=reference/(water/(2*concentration_water_mM))
    coefficient=.01/128
    predicted=scale*coefficient
    assert predicted==pytest.approx(9.375)
    assert predicted!=pytest.approx(10.)

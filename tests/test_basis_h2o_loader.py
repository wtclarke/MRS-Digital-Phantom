"""Supplied BASIS loader validation without standalone simulation adapters."""
import numpy as np
import pytest
from simulation.basis_h2o import load_responses


def load(path,**options):
    return load_responses(path,frequency_origin_ppm=4.68,bandwidth_hz=2000,points=256,**options)


def test_full_basis_complex_responses_and_identity(basis):
    path,expected=basis
    responses,water,metadata=load(path)
    for name in responses:np.testing.assert_allclose(responses[name],expected[name],rtol=2e-13,atol=1e-14)
    np.testing.assert_allclose(water,expected['H2O'],rtol=2e-13,atol=1e-14)
    assert metadata['physical_calibration'] is False
    assert metadata['normalization']['scalar']==1.
    with pytest.raises(ValueError,match='SHA256'):load(path,expected_basis_sha256='wrong')


@pytest.mark.parametrize('change',['missing_water','zero_water','nonfinite','sampling','scale','duplicate'])
def test_bad_basis_rejected(basis,change):
    path,_=basis;text=path.read_text()
    if change=='missing_water':text=text.replace("METABO='H2O'","METABO='other'")
    elif change=='scale':text=text.replace('TRAMP=1','TRAMP=2',1)
    elif change=='sampling':text=text.replace('NDATAB=256','NDATAB=255')
    elif change=='duplicate':text=text.replace("METABO='Cr'","METABO='NAA'")
    else:
        idx=text.index("METABO='H2O'");a=text.index('$END',idx)+4
        text=text[:a]+'\n'+(' nan 0\n' if change=='nonfinite' else ' 0 0\n'*256)
    path.write_text(text)
    with pytest.raises((ValueError,TypeError)):load(path)

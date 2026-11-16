"""Verify the actual first local fit; never run proprietary fitting during pytest."""
import json
from pathlib import Path

import numpy as np
import pytest

from examples.lcmodel_initial_validation import BINARY,OUTPUT,SOURCE,run
from utils.lcmodel_quantitative_io import NAMES,independent_water_scaling,parse_raw,parse_table


@pytest.fixture(scope='module')
def actual_fit():
    if not BINARY.exists() or not (OUTPUT/'noise_free/execution.json').exists():
        pytest.skip('requires local LCModel binary and actual Phase 7b first-fit artifacts')
    return json.loads((OUTPUT/'noise_free/fitting_results.json').read_text())


def test_actual_first_fit_outputs_and_warnings(actual_fit):
    assert actual_fit['version']=='6.3-1N' and not actual_fit['fatal']
    assert set(actual_fit['concentrations'])==set(NAMES)
    assert actual_fit['execution']['returncode']==0
    assert actual_fit['source_files_unchanged']
    assert 0<actual_fit['relative_phased_real_spectrum_residual']<.03
    assert not actual_fit['complex_residual_available']
    assert any(d['routine']=='FINOUT' and d['code']==7 for d in actual_fit['diagnostics'])
    for extension in ['TABLE','COORD','PRINT','PS','pdf','png']:
        assert (OUTPUT/f'noise_free/lcmodel_outputs/fit.{extension}').stat().st_size>100


def test_source_arrays_and_basis_areas_verified_independently(actual_fit):
    from fsl_mrs.utils import mrs_io
    for name,rawname in [('metabolites','metabolites.RAW'),('water_reference','water_reference.H2O')]:
        source=mrs_io.read_FID(str(SOURCE/f'{name}.nii.gz')).mrs().FID
        stored,_=parse_raw(OUTPUT/f'noise_free/{rawname}')
        np.testing.assert_array_equal(stored.conj(),source)
        assert actual_fit['validation'][name]['max_absolute_complex_error']==0
    stored,_=parse_raw(OUTPUT/'noise_free/water_reference.H2O')
    reconstructed=independent_water_scaling('outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative/metabolite_basis.BASIS',stored)
    for key in ['water_area','normalized_basis_area','FCALIB']:
        # LCModel single-precision transforms and five-significant-digit PRINT output.
        assert reconstructed[key]==pytest.approx(actual_fit['printed_scaling'][key],rel=3e-5)
    assert reconstructed['derived_basis_time_normalization']==128


def test_deterministic_conversion_revalidation_no_refitting(actual_fit,monkeypatch):
    import examples.lcmodel_initial_validation as validation
    def no_external_call(*args,**kwargs):raise AssertionError('cached validation must not rerun LCModel')
    monkeypatch.setattr(validation.subprocess,'run',no_external_call)
    repeated=run(reuse_fit=True)
    assert repeated==actual_fit


def test_original_basis_failure_preserved(actual_fit):
    original=parse_table(OUTPUT/'noise_free/lcmodel_outputs/original_basis.TABLE')
    assert original['fatal'] and any(d['routine']=='MYBASI' and d['code']==31 for d in original['diagnostics'])

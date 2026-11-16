"""Actual cached LCModel controls; fitting remains an explicit example command."""
import json
import re

import numpy as np
import pytest

from examples.lcmodel_fit_bias_diagnostics import BASE,CASES,OUTPUT,run
from examples.lcmodel_initial_validation import digest
from utils.lcmodel_quantitative_io import parse_coord


@pytest.fixture(scope='module')
def diagnostics():
    path=OUTPUT/'results.json'
    if not path.exists():pytest.skip('requires actual Phase 7b.1 local LCModel outputs')
    return json.loads(path.read_text())


def test_baseline_preservation_and_actual_replay(diagnostics):
    hashes=json.loads((OUTPUT/'immutable_baseline_sha256.json').read_text())
    assert hashes=={p:digest(p) for p in hashes}
    assert diagnostics['immutable_files_unchanged']
    original=parse_coord(BASE/'lcmodel_outputs/fit.COORD')
    replay=parse_coord(OUTPUT/'replay/lcmodel_outputs/fit.COORD')
    for name in original:np.testing.assert_array_equal(original[name],replay[name])
    assert diagnostics['cases']['replay']['concentrations']==diagnostics['baseline']['concentrations']


def test_independent_normalization_and_coefficient_projection(diagnostics):
    analysis=diagnostics['analysis'];cases=diagnostics['cases']
    for key,case in [('baseline','replay'),('manual_default','default_basis_width')]:
        independent=analysis['normalization_checks'][key]
        for factor in ['FCALIB','normalized_basis_area','water_area']:
            # Single-precision complex-exp recurrence and rounded five-digit PRINT.
            assert independent[factor]==pytest.approx(cases[case][factor],rel=6e-6)
    matching=analysis['projections']['matching_FSL_complex_acquired']
    assert matching['relative_residual']<1e-13
    np.testing.assert_allclose(matching['source_equivalent_mM'],analysis['truth_mM'],rtol=1e-13)
    projection=analysis['projections']['LCModel_complex_acquired']
    np.testing.assert_allclose(projection['source_equivalent_mM'],analysis['truth_mM'],rtol=2e-5)
    assert analysis['projections']['LCModel_full_spectral_real']['relative_residual']>.03
    assert analysis['projections']['LCModel_truncated_spectral_real']['relative_residual']<1e-5
    assert analysis['unitary_FFT_factor']==128


def test_controls_keep_identical_raw_water_and_basis(diagnostics):
    reference=(BASE/'lcmodel.CONTROL').read_text()
    for case in CASES:
        control=(OUTPUT/case/'lcmodel.CONTROL').read_text()
        for field in ['FILRAW','FILH2O','FILBAS']:
            expression=rf"{field}='([^']+)'"
            assert re.search(expression,control)[1]==re.search(expression,reference)[1]
        for output in ['PS','pdf']:
            assert (OUTPUT/case/f'lcmodel_outputs/fit.{output}').stat().st_size>100
    assert len(diagnostics['cases'])==8


def test_fitting_controls_separate_model_and_scaling_effects(diagnostics):
    cases=diagnostics['cases'];base=cases['replay']
    for name in ['preliminary_matched','half_broadening_prior','wider_convolution','phantom_mode','known_phase_prior']:
        assert cases[name]['FCALIB']==base['FCALIB']
    assert cases['preliminary_matched']['concentrations']==base['concentrations']
    assert cases['wider_convolution']['relative_real_residual']==base['relative_real_residual']
    assert cases['phantom_mode']['relative_real_residual']>base['relative_real_residual']
    assert abs(cases['phantom_mode']['concentrations']['NAA']['concentration']-10)<abs(base['concentrations']['NAA']['concentration']-10)
    assert cases['phantom_mode']['baseline_knots']<base['baseline_knots']
    assert cases['phantom_mode']['baseline_relative_norm']<base['baseline_relative_norm']
    assert any(d['routine']=='TWOREG' and d['code']==7 for d in cases['half_broadening_prior']['diagnostics'])
    width=cases['default_basis_width']
    assert width['FCALIB']!=base['FCALIB']
    assert width['concentrations']['NAA']['concentration']/base['concentrations']['NAA']['concentration']==pytest.approx(width['FCALIB']/base['FCALIB'],rel=5e-5)


def test_zero_broadening_failure_is_retained(diagnostics):
    failed=diagnostics['cases']['zero_broadening_prior']
    assert failed['fatal'] and not failed['concentrations']
    assert failed['execution']['returncode']==0
    assert any(d['severity']=='fatal' and d['routine']=='SETUP' and d['code']==14 for d in failed['diagnostics'])


def test_cached_analysis_reproducible_without_new_fits(diagnostics,monkeypatch):
    from matplotlib.figure import Figure
    def forbidden(*args,**kwargs):raise AssertionError('cached diagnostic reran an external program')
    monkeypatch.setattr('examples.lcmodel_fit_bias_diagnostics.subprocess.run',forbidden)
    # Cache verification must not redraw protected reports using another test's global style.
    monkeypatch.setattr(Figure,'savefig',lambda *args,**kwargs: None)
    repeated=run()
    assert repeated==diagnostics['cases']
    for filename in ['spectral_comparisons.png','concentration_vs_residual.png','comparison.csv','commands.txt']:
        assert (OUTPUT/filename).stat().st_size>100

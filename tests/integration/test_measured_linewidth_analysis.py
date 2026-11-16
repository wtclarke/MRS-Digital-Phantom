from pathlib import Path
from unittest.mock import patch

import pytest

from examples.mrscloud_measured_linewidth_validation import run


@pytest.fixture(scope='module')
def analysis(tmp_path_factory):
    root=tmp_path_factory.mktemp('measured_linewidths')/'analysis'
    def forbidden(*args,**kwargs):
        raise AssertionError('read-only linewidth validation must not refit or requantify')
    with patch('fsl_mrs.utils.fitting.fit_FSLModel',side_effect=forbidden), \
         patch('fsl_mrs.utils.quantify.quantifyWater',side_effect=forbidden):
        result=run(root)
    return result,root


def test_all_conditions_raw_equal_and_processed_differences_are_small(analysis):
    result,_=analysis
    assert result['source_hashes_unchanged'] and result['fit_reexecuted'] is False
    assert len(result['cases'])==12
    for case in result['cases']:
        row=case['row'];width=row['nominal_linewidth_hz']
        assert abs(row['raw_methyl_minus_water_hz'])<1e-11
        assert abs(row['processed_methyl_minus_water_hz'])<.023
        assert abs(row['full_NAA_minus_water_hz'])<.012
        assert row['methyl_DC_background_removed_fwhm_hz']==pytest.approx(row['methyl_raw_fwhm_hz'],abs=1e-11)
        assert abs(row['methyl_formatted_DC_background_removed_fwhm_hz']-row['water_raw_fwhm_hz'])<abs(row['processed_methyl_minus_water_hz'])
        assert abs(row['methyl_zero_fill_doubling_difference_hz'])<3e-5
        assert abs(row['fitted_water_fwhm_hz']-row['water_raw_fwhm_hz'])<.001
        assert abs(row['fitted_NAA_fwhm_hz']-row['full_formatted_fwhm_hz'])<2e-5
        assert row['fitted_NAA_additional_Lorentzian_parameter_hz']==pytest.approx(width-1,abs=1e-5)


def test_equal_linewidths_do_not_imply_equal_reference_integrals(analysis):
    result,_=analysis
    for case in result['cases']:
        row=case['row']
        assert abs(row['raw_width_only_surrogate_delta_mM'])<1e-10
        assert abs(row['processed_width_only_surrogate_delta_mM'])<.003
        assert abs(row['original_reported_mM']-10)>.2
    first=result['cases'][0]['row']
    assert first['original_ideal_methyl_factor']==pytest.approx(.9098704420596523,rel=1e-12)
    assert first['water_raw_fwhm_hz']==pytest.approx(first['methyl_raw_fwhm_hz'],abs=1e-12)


def test_actual_tables_phasing_and_plots(analysis):
    result,root=analysis
    assert result['zero_fill_factor']==128
    assert abs(result['methyl_phase_rad'])<1e-12
    assert abs(result['water_phase_rad'])<1e-12
    assert (root/'linewidths.csv').stat().st_size>1000
    for label in ['N2048_LW1','N2048_LW8','N8192_LW1']:
        assert (root/f'{label}_absorption.png').stat().st_size>1000

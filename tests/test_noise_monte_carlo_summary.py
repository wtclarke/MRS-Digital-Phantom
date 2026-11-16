"""Closed numerical checks of bias, sample SD, RMSE and coverage accounting."""
import math

import pytest

from examples.mrscloud_noise_fitting_validation import summarize,NAMES


def test_summary_statistics_against_hand_calculation():
    records=[]
    for index,error in enumerate([-1.,0.,1.]):
        coefficients={n:dict(absolute_error_mM=error,estimate_mM=10+error,
            covered_95=(error==0),relative_error=error/10,standard_FSL_mM=7+error) for n in [*NAMES,'GPC+PCh']}
        records.append(dict(snr='25',fit=dict(coefficients=coefficients,optimizer_success=index!=0,quantification_success=True)))
    for row in summarize(records):
        assert row['mean_bias_mM']==0 and row['SD_mM']==1
        assert row['RMSE_mM']==pytest.approx(math.sqrt(2/3))
        assert row['median_absolute_error_mM']==1
        assert row['coverage_95']==pytest.approx(1/3) and row['coverage_valid_count']==3
        assert row['failure_rate']==pytest.approx(1/3)
        assert row['poor_recovery_fraction']==0
        assert row['mean_standard_FSL_mM']==7 and row['SD_standard_FSL_mM']==1


def test_missing_fit_is_counted_as_failure_not_favorable_exclusion():
    good=dict(coefficients={n:dict(absolute_error_mM=0.,estimate_mM=10.,covered_95=True,
        relative_error=0.,standard_FSL_mM=7.) for n in [*NAMES,'GPC+PCh']},optimizer_success=True,quantification_success=True)
    rows=summarize([dict(snr='5',fit=good),dict(snr='5',fit=dict(optimizer_success=False,quantification_success=False,error='failure'))])
    for row in rows:
        assert row['realizations']==2 and row['returned_fits']==1
        assert row['failure_rate']==row['quantification_failure_rate']==.5
        assert row['coverage_valid_count']==1 and row['SD_mM'] is None

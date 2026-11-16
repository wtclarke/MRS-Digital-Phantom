"""Fast independent concentration predictions without fitting or report generation."""
import copy
import json

import numpy as np
import pytest

from examples.mrscloud_mixed_tissue_fitting_validation import SOURCE, BUNDLE, NAMES, prediction, quant_info
from simulation.mrscloud_full_basis import load_library


def test_prediction_uses_source_truth_not_fitted_results(monkeypatch):
    from fsl_mrs.utils import quantify
    from fsl_mrs.utils.results import FitRes
    def forbidden(*args, **kwargs):
        raise AssertionError('prediction called fitted concentration reporting')
    monkeypatch.setattr(quantify, 'quantifyWater', forbidden)
    monkeypatch.setattr(FitRes, 'getConc', forbidden)
    responses,water,_ = load_library(BUNDLE,selected_names=NAMES)
    truth = json.loads((SOURCE/'mixed/ground_truth.json').read_text())
    result = prediction(truth,responses,water,matched=True)
    # Matched quantification still corrects to WM/GM volume; finite-grid response remains.
    assert result['csf_correction'] == 2
    assert result['response_factor'] == pytest.approx(.8984749391468103,rel=3e-7)
    assert result['weighted_relaxed_water_mM'] == pytest.approx(49722.45)
    for name,c in truth['tissue_only_metabolite_concentrations_mM'].items():
        assert result['predicted_mM'][name] == pytest.approx(c*result['response_factor'],rel=2e-15)
    double = copy.deepcopy(truth)
    double['voxel_average_metabolite_concentrations_mM'] = {n:2*c for n,c in truth['voxel_average_metabolite_concentrations_mM'].items()}
    for name,value in prediction(double,responses,water,matched=True)['predicted_mM'].items():
        assert value == 2*result['predicted_mM'][name]
    csf = json.loads((SOURCE/'pure_CSF/ground_truth.json').read_text())
    assert 'undefined' in prediction(csf,responses,water)['status']


def test_public_assumption_overrides_and_defaults_are_isolated():
    from fsl_mrs.utils import mrs_io
    mrs = mrs_io.read_FID(str(SOURCE/'mixed/metabolites.nii.gz')).mrs(basis=mrs_io.read_basis(str(BUNDLE/'metabolite_basis.BASIS')))
    # QuantificationInfo needs only selected names and central frequency here.
    mrs.keep = NAMES
    truth=json.loads((SOURCE/'mixed/ground_truth.json').read_text())
    q=quant_info(mrs,truth,True)
    assert q.te == .03 and q.tr == 2
    assert q.R_M == q.R_H2O_WM == q.R_H2O_GM == q.R_H2O_CSF == 1
    assert q.relax_corr_water_molar == pytest.approx(49722.45)
    assert q.csf_corr_molar == 2
    default=quant_info(mrs,truth,False)
    assert default.d_WM == .65 and default.t1['METAB'] == 1.29

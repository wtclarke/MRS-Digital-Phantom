import numpy as np
import pytest

from examples.fsl_mrs_sequence_water_quantification import run_validation, fit_case, SOURCE


@pytest.fixture(scope='module')
def validation(tmp_path_factory):
    out = tmp_path_factory.mktemp('matched_quantification')/'results'
    return run_validation(out), out


def test_default_reproduction_and_matched_quantification(validation):
    result, _ = validation
    ten = result['cases'][1]
    assert result['baseline_cli_mM'] == pytest.approx(6.885791303512254, rel=1e-12)
    assert ten['cases'][0]['reported_mM'] == pytest.approx(result['baseline_cli_mM'], rel=1e-12)
    assert ten['cases'][3]['reported_mM'] == pytest.approx(8.98475046112642, rel=2e-7)
    for case in result['cases']:
        assert case['metabolite_optimizer_success'] and case['water_fit']['success']
        assert case['water_fit']['relative_complex_residual'] < 1e-6
        assert case['metabolite_relative_residual'] < 5e-5
        a,b,c,d = case['cases']
        assert b['reported_mM']/a['reported_mM'] == pytest.approx(40000/35756.5, rel=2e-15)
        assert c['reported_mM']/a['reported_mM'] == pytest.approx(.6749567292246237/.5786670548409867, rel=2e-15)
        assert d['reported_mM'] == pytest.approx(b['reported_mM']*c['reported_mM']/a['reported_mM'], rel=2e-15)
        for row in case['cases']:
            assert abs(row['reconstruction_error_mM']) < 2e-11
            assert row['water_integral'] == pytest.approx(row['fsl_water_integral'], rel=2e-13)
            assert row['metabolite_integral'] == pytest.approx(row['fsl_metabolite_integral'], rel=2e-13)
            assert row['fsl_water_integral'] == a['fsl_water_integral']
            assert row['fsl_metabolite_integral'] == a['fsl_metabolite_integral']
            assert row['additional_factor'] == 1.
        assert d['water_relaxation_factor'] == d['metabolite_relaxation_factor'] == 1.


def test_linearity_and_response_factor_without_empirical_correction(validation):
    result, _ = validation
    concentrations = np.array([5,10,15,20])
    recovered = np.array([case['cases'][3]['reported_mM'] for case in result['cases']])
    np.testing.assert_allclose(recovered/concentrations, recovered[1]/10, rtol=3e-7)
    for case in result['cases']:
        predicted = case['ground_truth_mM']*case['response_factor_unfitted']
        assert case['cases'][3]['reported_mM'] == pytest.approx(predicted, rel=3e-7)
        assert .89 < case['response_factor_unfitted'] < .91
        assert case['raw_NAA_coefficient_mol_l'] == pytest.approx(case['ground_truth_mM']*.001, rel=2e-7)
        product = case['ideal_3H_2H_window_control_factor']*case['native_response_vs_ideal_3H_factor']*case['formatting_integral_factor']
        assert product == pytest.approx(case['response_factor_unfitted'], rel=2e-15)


def test_reports_provenance_and_deterministic_refitting(validation):
    result, out = validation
    assert result['source_hashes_unchanged']
    report = out/'NAA_10mM/fsl_mrs_results/report.html'
    assert report.stat().st_size > 1000 and 'Water referencing' in report.read_text()
    assert '--report' in (out/'cli_commands.txt').read_text()
    for filename in ['comparison.csv','results.json','quantification_comparison.png']:
        assert (out/filename).stat().st_size > 0
    repeat = fit_case(out/'NAA_10mM', SOURCE/'reference')
    assert repeat == result['cases'][1]
    for filename in ['metabolites.nii.gz','water_reference.nii.gz','simulation_config.json','ground_truth.json']:
        assert (out/'NAA_10mM'/filename).read_bytes() == (SOURCE/'mrscloud_sequence'/filename).read_bytes()

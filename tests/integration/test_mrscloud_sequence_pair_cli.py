import pytest

from examples.mrscloud_sequence_water_integration import validate


@pytest.fixture(scope='module')
def validated(tmp_path_factory):
    root=tmp_path_factory.mktemp('integrated_sequence_cli')
    return validate(root),root


def test_both_modes_cli_quantification_and_phase5f_reference(validated):
    rows,_=validated
    analytical,sequence=rows
    assert analytical['standard_fsl_mM']==pytest.approx(6.003318300954079,rel=2e-6)
    assert sequence['standard_fsl_mM']==pytest.approx(6.885791303512254,rel=2e-6)
    assert sequence['phase5f_water_max_error']==0
    assert sequence['phase5f_metabolite_max_error']==0
    for row in rows:
        assert row['water_optimizer_success']
        assert row['water_relative_residual']<1e-6
        assert row['NAA_relative_residual']<5e-5
        assert row['max_export_error']==0
        assert row['water_fit_amplitude']==pytest.approx(row['initial_water_amplitude'],rel=2e-6)


def test_actual_reports_primary_and_aligned_outputs(validated):
    rows,root=validated
    for row in rows:
        case=root/row['mode']
        report=case/'fsl_mrs_results/report.html'
        assert report.stat().st_size>1000
        assert 'Water referencing' in report.read_text()
        assert (case/'water_reference_aligned.nii.gz').stat().st_size>0
        assert (case/'water_preprocessing.json').stat().st_size>0
    assert (root/'cli_commands.txt').read_text().count('--report')==2

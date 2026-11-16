from examples.mrscloud_quantitative_voxel import validate


def test_real_example_outputs_and_phase6a_bundle_preservation(tmp_path):
    result=validate(tmp_path/'voxel_validation')
    assert result['source_bundle_unchanged']
    assert len(result['cases'])==4
    for row in result['cases']:
        assert row['independent_metabolite_max_error']<1e-18
        assert row['independent_water_max_error']<1e-17
        assert row['export_max_error']==0
        assert row['export_max_phase_error_rad']==0

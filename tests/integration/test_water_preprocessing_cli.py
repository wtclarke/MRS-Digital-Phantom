import json
from copy import deepcopy

import numpy as np
import pytest
from fsl_mrs.utils import mrs_io

from examples.fsl_mrs_water_preprocessing import validate
from examples.fsl_mrs_water_fit_diagnostics import fit_water
from simulation.mrscloud_quantitative_pair import generate_pair, export_pair, REFERENCE_BASIS
from simulation.water_reference_preprocessing import align_water_reference


@pytest.mark.parametrize('frequency,phase',[(0.,0.),(2.55462,0.),(0.,.6),(2.55462,.6)])
def test_controlled_export_and_installed_water_fit(tmp_path,frequency,phase):
    pair=generate_pair({'NAA':10},water_frequency_offset_hz=frequency,water_phase_rad=phase)
    aligned,metadata=align_water_reference(pair['water_reference_fid'],dwell_time_s=.0005,
                                         frequency_offset_hz=frequency,initial_phase_rad=phase)
    modified=deepcopy(pair)
    modified['water_reference_fid']=aligned
    modified['simulation_config']['water_preprocessing']=metadata
    paths=export_pair(modified,tmp_path)
    reader=mrs_io.read_FID(str(paths['metabolite_fid']))
    water=mrs_io.read_FID(str(paths['water_reference_fid']))
    np.testing.assert_allclose(water.mrs().FID,aligned,rtol=1e-14,atol=1e-13)
    np.testing.assert_array_equal(reader.mrs().FID,pair['metabolite_fid'])
    assert water.SpecFreqChemShift==reader.SpecFreqChemShift==4.68
    assert json.loads((tmp_path/'simulation_config.json').read_text())['water_preprocessing']==metadata
    mrs=reader.mrs(basis=mrs_io.read_basis(str(REFERENCE_BASIS)),ref_data=water)
    mrs.keep=['NAA'];mrs.rescaleForFitting()
    info,_,_=fit_water(mrs)
    assert info['success']
    assert info['relative_complex_residual']<1e-6
    assert info['fitted_amplitude']==pytest.approx(abs(pair['water_reference_fid'][0]),rel=2e-6)
    assert abs(info['fitted_frequency_hz'])<1e-6
    assert abs(info['fitted_phase_rad'])<1e-6


@pytest.fixture(scope='module')
def cli_validation(tmp_path_factory):
    root=tmp_path_factory.mktemp('preprocessing_cli')
    return validate(root),root


def test_real_cli_reports_and_signal_preservation(cli_validation):
    result,root=cli_validation
    assert result['metabolite_array_difference']==0
    assert result['max_magnitude_change']<5e-14
    for name in ['original','aligned']:
        report=root/name/'fsl_mrs_results'/'report.html'
        assert report.stat().st_size>1000
        row=result['datasets'][name]['comparison']
        assert row['cli_exit_code']==0
        assert row['water_optimizer_success'] and row['metabolite_optimizer_success']
        assert row['water_relative_complex_residual']<1e-6
        assert row['water_export_max_complex_error']<1e-12
        assert row['synthetic_water_referenced_mM']==pytest.approx(10,rel=2e-6)
    assert (root/'comparison.csv').stat().st_size>0
    assert (root/'water_fit_comparison.png').stat().st_size>0
    commands=(root/'cli_commands.txt').read_text()
    assert commands.count('--report')==2


def test_alignment_preserves_metabolite_coefficients_and_fixes_overlay(cli_validation):
    result,_=cli_validation
    original=result['datasets']['original']['comparison']
    aligned=result['datasets']['aligned']['comparison']
    assert original['raw_NAA_source_coefficient']==aligned['raw_NAA_source_coefficient']
    assert original['internal_NAA_coefficient']==aligned['internal_NAA_coefficient']
    assert original['water_report_overlay_relative_difference']>.5
    assert aligned['water_report_overlay_relative_difference']<1e-6
    assert aligned['fsl_standard_mM']==pytest.approx(original['fsl_standard_mM'],rel=2e-6)
    assert aligned['fitted_water_amplitude']==pytest.approx(original['fitted_water_amplitude'],rel=2e-6)

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from examples.mrscloud_quantification_sampling_sensitivity import load_master, fit_condition

ROOT=Path('outputs/mrscloud_sampling_sensitivity/validated')


@pytest.fixture(scope='module')
def experiment():
    return json.loads((ROOT/'results.json').read_text()),np.load(ROOT/'master_responses.npz')


def area(fid, n, limits):
    # Independent full-array reproduction; no production integral helper.
    weighted=fid.copy();weighted[0]*=.5
    spec=np.fft.fftshift(np.fft.fft(weighted,norm='ortho'))
    ppm=np.fft.fftshift(np.fft.fftfreq(n,.0005))/127.731+4.68
    lo,hi=sorted(np.argmin(abs(ppm-x)) for x in limits)
    return float(np.trapezoid(abs(spec[lo:hi+1].real)))


def test_all_twelve_unfitted_predictions_and_fitting_errors(experiment):
    result,master=experiment
    assert result['source_hashes_unchanged'] and len(result['cases'])==12
    assert result['normalization_scalar']==2779.198340948381
    for case in result['cases']:
        row=case['row'];n=row['points'];width=row['total_linewidth_hz'];time=np.arange(n)/2000
        metabolite=master['full_formatted'][:n]*np.exp(-np.pi*(width-1)*time)
        water=abs(master['water_raw'][0])*np.exp(-np.pi*width*time)
        im=area(metabolite,n,(1.8,2.2));iw=area(water,n,(1.65,7.65))
        assert row['unit_metabolite_integral']==pytest.approx(im,rel=3e-15)
        assert row['unit_water_integral']==pytest.approx(iw,rel=3e-15)
        assert row['predicted_mM']==pytest.approx(10*im/iw*2/3,rel=3e-15)
        assert abs(row['difference_from_prediction_mM'])<5e-4
        assert row['raw_NAA_coefficient_mol_l']==pytest.approx(.010,rel=2e-6)
        assert row['metabolite_success'] and row['water_success']
        assert row['metabolite_relative_residual']<3e-5
        # Preserve the installed optimizer's measured 8Hz residual, not a false
        # assertion that every noiseless water fit attains machine precision.
        assert row['water_relative_residual']<(5e-5 if width==8 else 1e-6)
        assert row['max_export_error']==0
        assert case['quantification']['water_relaxation_factor']==1
        assert case['quantification']['metabolite_relaxation_factor']==1
        assert case['quantification']['effective_water_M']==40


def test_initial_amplitudes_fixed_and_actual_cli_matches_api(experiment):
    from fsl_mrs.utils import mrs_io
    result,master=experiment
    report_count=0
    for case in result['cases']:
        row=case['row'];directory=ROOT/f"N{row['points']}_LW{row['total_linewidth_hz']}"
        config=json.loads((directory/'simulation_config.json').read_text())
        assert config['common_normalization_scalar']==result['normalization_scalar']
        for file,key,scale in [('metabolites.nii.gz','full_formatted',.010),('water_reference.nii.gz','water_raw',40)]:
            loaded=mrs_io.read_FID(str(directory/file))
            assert loaded.mrs().FID[0]==scale*master[key][0]
            assert loaded.dwelltime==.0005
            assert loaded.SpecFreqChemShift==4.68
        if 'CLI' in case:
            report_count+=1
            assert (directory/'fsl_mrs_results/report.html').stat().st_size>1000
            with (directory/'fsl_mrs_results/concentrations.csv').open() as stream:
                rows=list(csv.reader(stream))
            reported=dict(zip(rows[1][1:5],map(float,rows[2][1:5])))['molarity']
            assert reported==pytest.approx(row['reported_mM'],rel=1e-12)
            assert '--report' in case['CLI']['command']
            assert case['CLI']['density_configuration_ratio']==40000/(55010*.65)
    assert report_count==3


def test_sampling_dependence_reduces_but_convergence_is_not_unity(experiment):
    result,_=experiment
    spreads={width:np.ptp([case['row']['response_factor'] for case in result['cases']
                          if case['row']['total_linewidth_hz']==width]) for width in [1,3,5,8]}
    assert spreads[1]>.07
    assert all(spreads[w]<.002 for w in [3,5,8])
    assert all(abs(case['row']['response_factor']-1)>.015 for case in result['cases'])


def test_repeat_master_and_fit_are_deterministic(experiment,tmp_path):
    result,archive=experiment
    master,metadata=load_master('outputs/mrscloud_sampling_sensitivity/octave')
    for name in master:np.testing.assert_array_equal(master[name],archive[name])
    repeat=fit_condition(tmp_path/'repeat',master,metadata,2048,3)
    original=next(c for c in result['cases'] if c['row']['points']==2048 and c['row']['total_linewidth_hz']==3)
    assert repeat==original

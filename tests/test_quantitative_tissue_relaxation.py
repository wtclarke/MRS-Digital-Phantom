"""Independent attenuation and unchanged-baseline tests for the opt-in adapter."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest
from fsl_mrs.utils import mrs_io

from examples.mrscloud_quantitative_voxel import BUNDLE, METABOLITES, WATER
from examples.mrscloud_mixed_tissue_fitting_validation import SOURCE
from simulation.mrscloud_full_basis import load_library
from simulation.mrscloud_quantitative_pair import export_pair
from simulation.mrscloud_quantitative_voxel import simulate_voxel
from simulation.tissue_relaxation import attenuation
from examples.mrscloud_tissue_relaxation_validation import configuration


def simulate(fractions, **kwargs):
    args=dict(calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,water_concentrations_mol_l=WATER)
    args.update(kwargs)
    return simulate_voxel(fractions,**args)


@pytest.fixture(scope='module')
def units():
    return load_library(BUNDLE,selected_names=list(METABOLITES['WM']))


@pytest.mark.parametrize('label',['pure_WM','pure_GM','pure_CSF','mixed'])
def test_default_no_relaxation_matches_saved_phase6b_exactly(label):
    gt=json.loads((SOURCE/label/'ground_truth.json').read_text())
    pair=simulate(gt['tissue_fractions'])
    for key,name in [('metabolite_fid','metabolites'),('water_reference_fid','water_reference')]:
        np.testing.assert_array_equal(pair[key],mrs_io.read_FID(str(SOURCE/label/f'{name}.nii.gz')).mrs().FID)
    assert pair['ground_truth']==gt
    assert pair['simulation_config']==json.loads((SOURCE/label/'simulation_config.json').read_text())
    explicit=simulate(gt['tissue_fractions'],relaxation=dict(mode='none',TE_s=.03,TR_s=2.))
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_array_equal(explicit[key],pair[key])
    assert explicit['ground_truth']['relaxation']['water_factors']==dict.fromkeys(WATER,1.)


@pytest.mark.parametrize('fractions',[{'WM':1},{'GM':1},{'CSF':1},{'WM':.2,'GM':.3,'CSF':.5}])
def test_independent_pure_and_mixed_relaxation_fids(units,fractions):
    responses,water,_=units
    cfg=configuration(); pair=simulate(fractions,relaxation=cfg)
    expected_m=np.zeros(2048,complex); expected_w=np.zeros(2048,complex)
    for t,f in fractions.items():
        rw=(1-math.exp(-2/cfg['water'][t]['T1_s']))*math.exp(-.03/cfg['water'][t]['T2_s'])
        expected_w+=f*WATER[t]*rw*water
        for n,c in METABOLITES[t].items():
            p=cfg['metabolites'][t][n]
            r=(1-math.exp(-2/p['T1_s']))*math.exp(-.03/p['T2_s'])
            expected_m+=f*c*.001*r*responses[n]
    np.testing.assert_allclose(pair['metabolite_fid'],expected_m,rtol=3e-15,atol=1e-20)
    np.testing.assert_allclose(pair['water_reference_fid'],expected_w,rtol=3e-15,atol=1e-20)
    no_relax=simulate(fractions)
    assert pair['ground_truth']['voxel_average_metabolite_concentrations_mM']==no_relax['ground_truth']['voxel_average_metabolite_concentrations_mM']
    assert pair['ground_truth']['effective_voxel_water_concentration_mol_l']==no_relax['ground_truth']['effective_voxel_water_concentration_mol_l']
    # Scalar tissue attenuation changes no FID shape, frequency, phase or linewidth.
    for n in responses:
        c=pair['ground_truth']['relaxation_weighted_metabolite_concentrations_mM'][n]
        np.testing.assert_allclose(pair['components'][n],c*.001*responses[n],rtol=2e-15,atol=0)
    assert pair['simulation_config']['tissue_specific_relaxation']


def test_tr_te_dependence_and_attenuation_ratios():
    cfg=configuration(); fractions={'WM':1}
    a=simulate(fractions,relaxation=cfg)
    cfg['TR_s']=4.;b=simulate(fractions,relaxation=cfg)
    assert b['ground_truth']['relaxation_weighted_water_concentration_mol_l']>a['ground_truth']['relaxation_weighted_water_concentration_mol_l']
    ratio=(1-math.exp(-4/.97))/(1-math.exp(-2/.97))
    np.testing.assert_allclose(b['water_reference_fid'],ratio*a['water_reference_fid'],rtol=5e-16)
    rm=(1-math.exp(-2/1.29))*math.exp(-.03/.194)
    rw=(1-math.exp(-2/.97))*math.exp(-.03/.073)
    base=simulate(fractions)
    assert (a['components']['NAA'][0]/a['water_reference_fid'][0])/(base['components']['NAA'][0]/base['water_reference_fid'][0])==pytest.approx(rm/rw)
    r1=attenuation(TE_s=.03,TR_s=2.,T1_s=1.29,T2_s=.194)
    r2=attenuation(TE_s=.06,TR_s=2.,T1_s=1.29,T2_s=.194)
    assert r2/r1==pytest.approx(math.exp(-.03/.194),rel=2e-15)
    cfg['TE_s']=.06
    with pytest.raises(ValueError,match='match the frozen sequence TE'):
        simulate(fractions,relaxation=cfg)


def test_linearity_receiver_cancellation_and_determinism():
    cfg=configuration();before=deepcopy(cfg);fractions={'WM':.2,'GM':.3,'CSF':.5}
    a=simulate(fractions,relaxation=cfg);repeat=simulate(fractions,relaxation=cfg)
    g=simulate(fractions,relaxation=cfg,receiver_factor=2)
    m=simulate(fractions,relaxation=cfg,metabolite_concentrations_mM={t:{n:2*c for n,c in p.items()} for t,p in METABOLITES.items()})
    w=simulate(fractions,relaxation=cfg,water_concentrations_mol_l={t:2*c for t,c in WATER.items()})
    for key in ['metabolite_fid','water_reference_fid']:
        np.testing.assert_array_equal(a[key],repeat[key]);np.testing.assert_array_equal(g[key],2*a[key])
    np.testing.assert_array_equal(m['metabolite_fid'],2*a['metabolite_fid'])
    np.testing.assert_array_equal(w['water_reference_fid'],2*a['water_reference_fid'])
    np.testing.assert_array_equal(a['metabolite_fid']/a['water_reference_fid'],g['metabolite_fid']/g['water_reference_fid'])
    assert a['ground_truth']==repeat['ground_truth'] and a['simulation_config']==repeat['simulation_config']
    assert cfg==before


@pytest.mark.parametrize('change',['mode','missing','unknown_key','units','negative','zero','nan','infinite','unused_tissue','metabolite','tr','te'])
def test_invalid_relaxation_rejected_even_for_unused_tissues(change):
    cfg=configuration()
    if change=='mode':cfg['mode']='unknown'
    elif change=='missing':del cfg['water']
    elif change=='unknown_key':cfg['TE_ms']=30
    elif change=='units':cfg['water']['WM']={'T1_ms':970,'T2_ms':73}
    elif change in ('negative','zero','nan','infinite'):
        cfg['water']['CSF']['T1_s']={'negative':-1,'zero':0,'nan':np.nan,'infinite':np.inf}[change]
    elif change=='unused_tissue':del cfg['metabolites']['CSF']
    elif change=='metabolite':del cfg['metabolites']['GM']['Cr']
    elif change=='tr':cfg['TR_s']=.01
    elif change=='te':cfg['TE_s']=-.03
    with pytest.raises(ValueError):simulate({'WM':1},relaxation=cfg)


def test_relaxation_export_complex_metadata_and_acquisition_times(tmp_path):
    pair=simulate({'WM':.2,'GM':.3,'CSF':.5},relaxation=configuration())
    paths=export_pair(pair,tmp_path,aligned_water=True)
    for key in ['metabolite_fid','water_reference_fid']:
        obj=mrs_io.read_FID(str(paths[key]))
        np.testing.assert_array_equal(obj.mrs().FID,pair[key])
        assert obj.hdr_ext['EchoTime']==.03 and obj.hdr_ext['RepetitionTime']==2.
        assert obj.dwelltime==.0005 and obj.SpecFreqChemShift==4.68
    aligned=mrs_io.read_FID(str(paths['aligned_water_reference_fid']))
    assert aligned.hdr_ext['RepetitionTime']==2.
    assert json.loads((tmp_path/'ground_truth.json').read_text())==pair['ground_truth']
    assert json.loads((tmp_path/'simulation_config.json').read_text())==pair['simulation_config']

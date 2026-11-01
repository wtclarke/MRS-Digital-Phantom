"""Independent basis-derived response tests; no external BASIS or fitting required."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import torch
import torchio as tio
from fsl_mrs.utils import mrs_io

from simulation.basis_h2o import ASSUMPTION,load_responses,simulate_basis_h2o
from simulation.headless import run_headless,simulate_batch,results_to_nifti_mrs
from simulation.quantitative import export_quantitative,DEFAULT_WATER


@pytest.fixture
def basis(tmp_path):
    path=tmp_path/'test.BASIS';n=256;cf=100.;bw=2000;t=np.arange(n)/bw
    signals={m:a*np.exp(-np.pi*3*t+1j*(2*np.pi*f*t+p)) for m,a,f,p in [
        ('NAA',3,-267,.2),('Cr',3,-165,-.1),('GPC',9,-147,.05),('PCh',9,-148,.07),
        ('Glu',2,-233,-.2),('H2O',2,2,.4)]}
    text=" $SEQPAR\n HZPPPM=100, ECHOT=35, SEQ='PRESS'\n $END\n $BASIS1\n BADELT=0.0005, NDATAB=256, FMTBAS='(2E25.16)'\n $END\n"
    for name,x in signals.items():
        text+=f" $BASIS\n METABO='{name}', CONC=1, TRAMP=1, VOLUME=1, ISHIFT=0\n $END\n"
        text+=''.join(f'{z.real:25.16E} {z.imag:25.16E}\n' for z in np.fft.fft(x.conj()))
    path.write_text(text)
    return path,signals


def config(path):
    return dict(basis_path=str(path),concentration_units='mM',basis_amplitude_convention=ASSUMPTION,
        frequency_origin_ppm=4.68,bandwidth_hz=2000,points=256,tissue_fractions={'WM':1},
        metabolite_concentrations_mM={'WM':{'NAA':10,'Cr':8},'GM':{'NAA':14,'Cr':9},'CSF':{'NAA':0,'Cr':0}})


@pytest.mark.parametrize('fractions',[{'WM':1},{'GM':1},{'CSF':1},{'WM':.2,'GM':.3,'CSF':.5}])
def test_independent_tissue_mixing_and_no_extra_protons(basis,fractions):
    path,raw=basis;q=config(path);q['tissue_fractions']=fractions;pair=simulate_basis_h2o(q)
    mean={m:sum(fractions.get(t,0)*p[m] for t,p in q['metabolite_concentrations_mM'].items()) for m in ['NAA','Cr']}
    cw=sum(fractions.get(t,0)*c for t,c in DEFAULT_WATER.items())
    np.testing.assert_allclose(pair['metabolite_fid'],sum(.001*mean[m]*raw[m] for m in mean),rtol=2e-13,atol=1e-15)
    np.testing.assert_allclose(pair['water_reference_fid'],cw*raw['H2O'],rtol=2e-13,atol=1e-13)
    assert pair['ground_truth']['voxel_average_metabolite_concentrations_mM']==mean
    assert pair['ground_truth']['effective_voxel_water_concentration_mol_l']==cw
    if fractions=={'CSF':1}:assert pair['ground_truth']['tissue_only_metabolite_concentrations_mM']['NAA'] is None


def test_selection_linearity_gain_and_ratio(basis):
    path,_=basis;q=config(path);a=simulate_basis_h2o(q)
    only=deepcopy(q)
    for p in only['metabolite_concentrations_mM'].values():del p['Cr']
    b=simulate_basis_h2o(only)
    np.testing.assert_array_equal(a['components']['NAA'],b['components']['NAA'])
    np.testing.assert_array_equal(a['water_reference_fid'],b['water_reference_fid'])
    for field in ['receiver_factor','water_concentrations_mol_l','metabolite_concentrations_mM']:
        c=deepcopy(q)
        if field=='receiver_factor':c[field]=2
        elif field=='water_concentrations_mol_l':c[field]={t:2*v for t,v in DEFAULT_WATER.items()}
        else:c[field]={t:{m:2*v for m,v in p.items()} for t,p in c[field].items()}
        result=simulate_basis_h2o(c)
        if field!='water_concentrations_mol_l':np.testing.assert_array_equal(result['metabolite_fid'],2*a['metabolite_fid'])
        if field!='metabolite_concentrations_mM':np.testing.assert_array_equal(result['water_reference_fid'],2*a['water_reference_fid'])
    gain=simulate_basis_h2o(dict(q,receiver_factor=3))
    np.testing.assert_allclose(gain['components']['NAA']/gain['water_reference_fid'],a['components']['NAA']/a['water_reference_fid'],rtol=5e-16)
    assert a['simulation_config']['synthetic_calibration']['normalization']['scalar']==1.


def test_spatial_relaxation_linewidth_noise_and_export(basis,tmp_path):
    path,raw=basis;q=config(path);q.pop('tissue_fractions')
    q.update(spatial={'cell_labels':['WM','GM','CSF'],'spatial_weights':[.2,.3,.5],
        'b0_offsets_hz':[-5,3,37],'b0_units':'Hz','b0_convention':'positive_frequency'},
        relaxation={'mode':'tissue_specific','TE_s':.035,'TR_s':2.,
          'water':{t:{'T1_s':1.+i,'T2_s':.1+i*.1} for i,t in enumerate(DEFAULT_WATER)},
          'metabolites':{t:{m:{'T1_s':1.3+i*.1,'T2_s':.2+i*.1} for m in ['NAA','Cr']} for i,t in enumerate(DEFAULT_WATER)}},
        basis_relaxation_assumption='no_tissue_attenuation_encoded',additional_linewidth_hz={'metabolites':2,'water':7})
    clean=simulate_basis_h2o(q);time=np.arange(256)/2000;em=np.zeros(256,complex);ew=em.copy()
    def r(p):return (1-np.exp(-2/p['T1_s']))*np.exp(-.035/p['T2_s'])
    for tissue,weight,df in zip(DEFAULT_WATER,[.2,.3,.5],[-5,3,37]):
        for m,c in q['metabolite_concentrations_mM'][tissue].items():em+=weight*.001*c*r(q['relaxation']['metabolites'][tissue][m])*raw[m]*np.exp(2j*np.pi*df*time)*np.exp(-2*np.pi*time)
        ew+=weight*DEFAULT_WATER[tissue]*r(q['relaxation']['water'][tissue])*raw['H2O']*np.exp(2j*np.pi*df*time)*np.exp(-7*np.pi*time)
    np.testing.assert_allclose(clean['metabolite_fid'],em,rtol=2e-13,atol=1e-15)
    np.testing.assert_allclose(clean['water_reference_fid'],ew,rtol=2e-13,atol=1e-13)
    state=np.random.get_state();ts=torch.random.get_rng_state().clone()
    q['noise']={'seed':23,'metabolite_snr':50,'water_snr':1000}
    a=simulate_basis_h2o(q);b=simulate_basis_h2o(q)
    for key in ['metabolite_fid','water_reference_fid','clean_metabolite_fid','clean_water_reference_fid']:np.testing.assert_array_equal(a[key],b[key])
    np.testing.assert_array_equal(a['clean_metabolite_fid'],clean['metabolite_fid'])
    for actual,expected in zip(np.random.get_state(),state):np.testing.assert_array_equal(actual,expected)
    assert torch.equal(torch.random.get_rng_state(),ts)
    for name,path in export_quantitative(a,tmp_path/'out').items():
        key={'metabolites':'metabolite_fid','water_reference':'water_reference_fid','metabolites_clean':'clean_metabolite_fid','water_reference_clean':'clean_water_reference_fid'}[name]
        reader=mrs_io.read_FID(path);np.testing.assert_array_equal(reader.mrs().FID,a[key]);assert reader.dwelltime==.0005
        assert reader.SpecFreqChemShift==4.68 and reader.hdr_ext['EchoTime']==.035


@pytest.mark.parametrize('change',[{'basis_amplitude_convention':'independently_normalized'}, {'concentration_units':'IU'},
    {'frequency_origin_ppm':None},{'expected_basis_sha256':'wrong'},{'water_component_name':'missing'},
    {'points':10000},{'bandwidth_hz':4000},{'points':True},{'receiver_factor':-1},
    {'water_amp_factor':3},{'noise':{'seed':True}},{'calibration_bundle':'anything'},
    {'spatial':{'cell_labels':['WM'],'b0_offsets_hz':[1]}}])
def test_rejected_configurations(basis,change):
    q=config(basis[0]);q.update(change)
    with pytest.raises((ValueError,TypeError)):simulate_basis_h2o(q)


@pytest.mark.parametrize('change',['missing_water','zero_water','nonfinite','sampling','scale'])
def test_bad_basis_rejected(basis,change):
    path,_=basis;text=path.read_text()
    if change=='missing_water':text=text.replace("METABO='H2O'","METABO='other'")
    elif change=='scale':text=text.replace('TRAMP=1','TRAMP=2',1)
    elif change=='sampling':text=text.replace('NDATAB=256','NDATAB=255')
    elif change=='nonfinite':text=text.replace(' $END\n', ' $END\n',1);idx=text.index("METABO='H2O'");a=text.index('$END',idx)+4;text=text[:a]+'\n nan 0\n'+text[a:]
    else:
        idx=text.index("METABO='H2O'");a=text.index('$END',idx)+4;text=text[:a]+'\n'+' 0 0\n'*256
    path.write_text(text)
    with pytest.raises((ValueError,TypeError)):simulate_basis_h2o(config(path))


def test_headless_voxel_position_size_anatomy_and_cli(basis,tmp_path):
    path,_=basis;labels=np.array([1,1,2,2,2,3,3,3],dtype=np.int16).reshape(2,2,2)
    image=tio.LabelMap(tensor=torch.tensor(labels)[None],affine=np.diag([2,2,2,1]))
    from simulation.mrs_phantom import DigitalPhantom
    phantom=DigitalPhantom.__new__(DigitalPhantom);phantom.subject=tio.Subject(labels=image)
    phantom.label_names=['Background','WM','GM','CSF'];phantom.spacing=(2,2,2);phantom.affine=image.affine
    q=config(path);q.pop('tissue_fractions');q.pop('basis_path')
    cfg={'water_reference_mode':'basis_h2o','path2basis':str(path),'quantitative':q,
         'positions':[[1,1,1]],'voxel_size_mm':[4,4,4]}
    a=run_headless(cfg,phantom=phantom)[0]
    assert a['ground_truth']['tissue_fractions']=={'WM':.25,'GM':.375,'CSF':.375}
    reference=simulate_basis_h2o(dict(q,basis_path=str(path),tissue_fractions={'WM':.25,'GM':.375,'CSF':.375}))
    for key in ['metabolite_fid','water_reference_fid']:np.testing.assert_array_equal(a[key],reference[key])
    with pytest.raises(ValueError):results_to_nifti_mrs([a])
    with pytest.raises(ValueError):run_headless(dict(cfg,simulation_mode='sequence_matched_quantitative'),phantom=phantom)
    simple=deepcopy(cfg);simple.pop('voxel_size_mm');simple['quantitative']['tissue_fractions']={'WM':1}
    simple['output_dir']=str(tmp_path/'cli');file=tmp_path/'config.json';file.write_text(json.dumps(simple))
    done=subprocess.run([sys.executable,'-m','simulation.headless','--config',str(file)],capture_output=True,text=True)
    assert done.returncode==0,done.stderr
    assert (tmp_path/'cli/voxel_0000_1_1_1/calibration_provenance.json').exists()


def test_no_frozen_bundle_or_octave_required(basis,monkeypatch):
    import simulation.quantitative as production
    import simulation.mrscloud_full_basis as frozen
    def forbidden(*a,**kw):raise AssertionError('frozen Phase 6a calibration loaded')
    monkeypatch.setattr(production,'validate_config',forbidden)
    monkeypatch.setattr(frozen,'load_library',forbidden)
    result=simulate_basis_h2o(config(basis[0]))
    assert result['simulation_config']['water_reference_mode']=='basis_h2o'
    assert result['simulation_config']['basis_file_sha256_verified']


def test_relaxation_requires_declaration_and_actual_te(basis):
    from examples.quantitative_water_reference_production import relaxation_config
    q=config(basis[0]);q['relaxation']=relaxation_config()
    for tissue in q['relaxation']['metabolites']:
        q['relaxation']['metabolites'][tissue]={n:q['relaxation']['metabolites'][tissue][n] for n in ['NAA','Cr']}
    with pytest.raises(ValueError,match='basis_relaxation_assumption'):simulate_basis_h2o(q)
    q['basis_relaxation_assumption']='no_tissue_attenuation_encoded'
    with pytest.raises(ValueError,match='TE_s must match'):simulate_basis_h2o(q)
    q['relaxation']['TE_s']=.035
    result=simulate_basis_h2o(q)
    assert result['simulation_config']['TE_s']==.035

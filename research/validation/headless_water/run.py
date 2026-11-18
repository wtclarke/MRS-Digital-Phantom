"""Reproducible paired-fitter characterization of the existing headless simulator."""
import argparse
from contextlib import ExitStack, redirect_stdout
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import subprocess
from unittest.mock import patch

import numpy as np
import pandas as pd
import torch
import torchio as tio
from scipy.ndimage import gaussian_filter
from examples.legacy_basis_water import configuration
from simulation import mrs_phantom
from simulation.headless import simulate_batch, results_to_nifti_mrs
from simulation.signalModel import SignalModel
from simulation.basis_h2o import load_responses
from utils.auxillary import fid2spec

BASES=Path('/Users/wclarke/Documents/VoxelSystems/voxel_systems_example_basis/LCM_basis_2')
MIXTURES=[(1,0,0),(0,1,0),(.5,.5,0),(.4,.4,.2),(.25,.25,.5),(.1,.1,.8),(0,0,1)]
TARGETS=['NAA','Cr','bHG','Glu','GABA','GSH','tCho']
GROUPS={'tNAA':['NAA','NAAG'],'tCr':['Cr','PCr'],'tCho':['GPC','PCh']}
SCALES=[0,.5,.75,1,1.25,1.5,2]


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,data):Path(path).write_text(json.dumps(data,indent=2,default=str)+'\n')
def groups(values):
    return {**values,**{g:sum(values.get(n,0) for n in ns) for g,ns in GROUPS.items()}}


def design(bases):
    cases=[]
    def add(b,stage,mix,target='baseline',scale=1,snr=float('inf'),b0_sd=0,water_snr=float('inf'),repeat=0,**extra):
        cases.append(dict(basis_path=str(b.resolve()),basis_name=b.stem,stage=stage,
            tissue='/'.join(str(round(f*100)) for f in mix),wm=mix[0],gm=mix[1],csf=mix[2],
            target=target,scale=scale,snr=snr,b0_sd=b0_sd,water_snr=water_snr,repeat=repeat,**extra))
    detailed=[bases/'GE_PRESS_TE35.BASIS',bases/'SH_PRESS_TE97.BASIS']
    changes=[('baseline',1)]+[(t,s) for t in ['NAA','bHG'] for s in [.5,1.5]]
    for b in detailed:
        for mix in MIXTURES:
            add(b,'clean',mix)
            if mix[-1]<1:
                for t in TARGETS:
                    for s in SCALES:
                        if s!=1:add(b,'clean',mix,t,s)
        for mix in [MIXTURES[0],MIXTURES[1],MIXTURES[3]]:
            selected=changes if mix==MIXTURES[3] else changes[:1]
            for t,s in selected:
                for snr in [100,30,10]:
                    for r in range(10):add(b,'noise',mix,t,s,snr=snr,repeat=r)
                for sd in [2,5,10]:add(b,'shim',mix,t,s,b0_sd=sd)
        for t,s in changes:
            for sd in [5,10]:
                for r in range(10):add(b,'combined',MIXTURES[3],t,s,snr=30,b0_sd=sd,repeat=r)
        for ws in [200,50]:
            for r in range(10):add(b,'water_noise',MIXTURES[3],snr=30,water_snr=ws,repeat=r)
        for factor in [.5,1,2]:add(b,'water_control',MIXTURES[3],water_factor=factor)
        add(b,'water_control',MIXTURES[3],water_factor=2,mismatched=True)
    for b in sorted(bases.glob('*.BASIS')):
        if b not in detailed:
            for t,s in changes:add(b,'transfer',MIXTURES[3],t,s)
    from utils.loading import load_skeleton_as_subject
    anatomy,_=load_skeleton_as_subject('data/skeletons/MRiLab','MRiLab')
    labels=anatomy.labels.data[0].numpy();counts=np.round(20/np.array(anatomy.spacing)).astype(int)
    candidates=[]
    for x in range(counts[0],labels.shape[0]-counts[0],10):
        for y in range(counts[1],labels.shape[1]-counts[1],10):
            for z in range(counts[2],labels.shape[2]-counts[2],10):
                pos=np.array([x,y,z]);start=pos-counts//2;end=start+counts
                cube=labels[tuple(slice(a,b) for a,b in zip(start,end))]
                f=np.array([(cube==i).mean() for i in [1,2,3]])
                if np.isclose(f.sum(),1,atol=1e-12,rtol=0):candidates.append((tuple(pos.tolist()),f))
    chosen=[]
    for desired in [(1,0,0),(0,1,0),(.4,.4,.2)]:
        pos,f=min(candidates,key=lambda p:np.linalg.norm(p[1]-desired))
        chosen.append((pos,f))
    for b in detailed:
        for pos,f in chosen:
            for t,scale in changes:add(b,'anatomy',tuple(f),t,scale,anatomical=True,position=pos)
        pos,f=chosen[-1]
        for realistic in [False,True]:add(b,'realism',tuple(f),anatomical=True,position=pos,relax=True,realistic=realistic)
    return cases


def concentration_table(names,target,scale,relax=False):
    source=pd.read_csv('data/metabolites/metab_df.csv')
    rows=[]
    for i,t in enumerate(['Background','WM','GM','CSF','Fat','Skull']):
        for name in names:
            alias={'Ins':'mI','Scyllo':'sI'}.get(name,name)
            match=source[(source.Tissue==t)&(source.Metabolite==alias)]
            c=float(match.iloc[0].Conc_mean) if len(match) else 0
            c=c if np.isfinite(c) and t in ['WM','GM'] else 0.
            if name=='bHG':c={'WM':4.,'GM':3.}.get(t,0.)
            if name==target or name in GROUPS.get(target,[]):c*=scale
            t2=float(match.iloc[0].T2) if relax and len(match) else float('nan')
            rows.append(dict(Metabolite=alias,Label=i,Tissue=t,Conc_mean=c,Conc_std=0.,T1=np.nan,T2=t2))
    return pd.DataFrame(rows)


def subject(case):
    labels=np.repeat([1,2,3],[round(1000*case[k]) for k in ['wm','gm','csf']])
    np.random.default_rng(715).shuffle(labels)
    cube=np.zeros((1,12,12,12),dtype=np.int16)
    cube[:,1:11,1:11,1:11]=labels.reshape(1,10,10,10)
    cube[0,0,0,:6]=np.arange(6)
    labels=torch.tensor(cube)
    field=gaussian_filter(np.random.default_rng(310).normal(size=(10,10,10)),1)
    field=(field-field.mean())/field.std()*case['b0_sd']
    padded=np.zeros((1,12,12,12));padded[:,1:11,1:11,1:11]=field
    affine=np.diag([2.,2.,2.,1.])
    return tio.Subject(labels=tio.LabelMap(tensor=labels,affine=affine),
        B0=tio.ScalarImage(tensor=torch.tensor(padded,dtype=torch.float32),affine=affine),
        brain_mask=tio.LabelMap(tensor=torch.ones(1,12,12,12),affine=affine))


def generate(case,directory):
    directory.mkdir(parents=True,exist_ok=True)
    responses,water,meta=load_responses(case['basis_path'],frequency_origin_ppm=4.68)
    table=concentration_table(list(responses),case['target'],case['scale'],case.get('relax',False))
    csv=directory/'concentrations.csv';table.to_csv(csv,index=False)
    cfg=configuration(case['basis_path']);cfg['path2metabs']=str(csv.resolve())
    cfg['voxel_size_mm']=(20.,20.,20.);cfg['basis_set_settings']['TE']=meta['sequence']['TE_ms']
    cfg['simulation_params'].update(noise_level=0.,water_amp_factor=0.,mm_level=0.,lipid_amp_factor=0.)
    if case.get('realistic'):cfg['simulation_params'].update(water_amp_factor=5.,mm_level=25.,lipid_amp_factor=40.)
    sub=subject(case) if not case.get('anatomical') else None
    applied={}
    original_apply=SignalModel.apply_shim_map
    def capture_apply(model,specs,field):
        applied['b0_sd_Hz']=float(torch.std(torch.as_tensor(field),unbiased=False))
        applied['b0_mean_Hz']=float(torch.as_tensor(field).mean())
        return original_apply(model,specs,field)
    with ExitStack() as stack:
        stack.enter_context(patch.object(SignalModel,'apply_shim_map',capture_apply))
        if sub is not None:
            stack.enter_context(patch.object(mrs_phantom,'load_skeleton_as_subject',lambda *a:(sub,['Background','WM','GM','CSF','Fat','Skull'])))
            stack.enter_context(patch.object(mrs_phantom.DigitalPhantom,'init_lipid_labels',lambda self:torch.zeros(12,12,12)))
        if not case.get('realistic'):
            stack.enter_context(patch.object(SignalModel,'simulate_mms',lambda self:torch.zeros(len(self.labels),self.basis.n,dtype=torch.complex64)))
            stack.enter_context(patch.object(SignalModel,'simulate_lipids',lambda self,*a:torch.zeros(self.basis.n,dtype=torch.complex64)))
        with redirect_stdout(io.StringIO()):result=simulate_batch([case.get('position',(6,6,6))],cfg,seed=123)[0]
        if np.isfinite(case['snr']):
            baseline=concentration_table(list(responses),'baseline',1)
            f=result['ground_truth']['tissue_fractions']
            concentration=sum(f.get(t,0)*float(baseline[(baseline.Tissue==t)&(baseline.Metabolite=='NAA')].Conc_mean.iloc[0]) for t in ['WM','GM'])
            divisor=result['water_reference_provenance']['normalization']['scalar']
            cfg['simulation_params']['noise_level']=float(np.max(abs(fid2spec(responses['NAA'])))*concentration*8000/divisor/case['snr']/2e4)
        if np.isfinite(case['water_snr']):
            cfg['water_reference']['noise_sigma']=float(np.max(abs(fid2spec(result['water_reference_fid'])))/(case['water_snr']*np.sqrt(2048)))
        if np.isfinite(case['snr']) or np.isfinite(case['water_snr']):
            with redirect_stdout(io.StringIO()):result=simulate_batch([case.get('position',(6,6,6))],cfg,seed=123+case['repeat'])[0]

    # SNR fixes the native spectral noise_level; concentration draws remain zero-SD.
    base=concentration_table(list(responses),'baseline',1)
    fractions=result['ground_truth']['tissue_fractions']
    naa=sum(fractions.get(t,0)*float(base[(base.Tissue==t)&(base.Metabolite=='NAA')].Conc_mean.iloc[0]) for t in ['WM','GM'])
    scalar=result['water_reference_provenance']['normalization']['scalar']
    peak=np.max(abs(fid2spec(responses['NAA'])))*naa*8000/scalar
    sigma=peak/case['snr']
    wp=np.max(abs(fid2spec(result['water_reference_fid'])))
    wsigma=wp/(case['water_snr']*np.sqrt(2048))
    factor=case.get('water_factor',1.)
    result['water_reference_fid']=factor*result['water_reference_fid']
    if not case.get('mismatched'):
        result['ground_truth']['water_reference']['nominal_water_concentration_mol_l']*=factor
        result['water_reference_provenance']['water_options']['water_concentrations_mol_l']={t:c*factor for t,c in result['water_reference_provenance']['water_options']['water_concentrations_mol_l'].items()}
    result['water_reference_provenance']['study_noise']={'metabolite_spectral_sigma':sigma,'water_FID_sigma':wsigma,'metabolite_peak_definition':'baseline isolated NAA magnitude spectral peak; achieved_metabolite_snr uses total metabolite peak','water_peak':float(wp)}
    if sub is not None:
        cw=sum(fractions[t]*c for t,c in {'WM':39.0165,'GM':47.2305,'CSF':55.5}.items())
        field=sub.B0.data.numpy()[:,1:11,1:11,1:11].reshape(-1);lab=sub.labels.data.numpy()[:,1:11,1:11,1:11].reshape(-1)
        c=np.array([0.,39016.5,47230.5,55500.])
        times=np.arange(1,2049)/2000
        expected=8*np.sum(c[lab,None]*np.exp(-2j*np.pi*field[:,None]*times),axis=0)*water/scalar
        err=float(np.linalg.norm(result['clean_water_reference_fid']-expected)/np.linalg.norm(expected))
        assert err<5e-6,err
        save(directory/'water_check.json',{'relative_complex_error':err,'applied_b0_sd_Hz':float(field.std())})
    result['water_reference_spectrum']=torch.as_tensor(fid2spec(result['water_reference_fid']))
    converted=results_to_nifti_mrs([result],output_dir=str(directory/'acquisitions'),individual_components=['Metabs'])
    source=Path(converted[0]['saved_paths']['total']).parent
    achieved_water_snr=float(wp/(wsigma*np.sqrt(2048))) if wsigma else float('inf')
    achieved_metabolite_snr=float(np.max(abs(np.asarray(result['components']['Metabs'])))/sigma) if sigma else float('inf')
    save(directory/'applied_field.json',applied)
    save(directory/'achieved_snr.json',dict(water=achieved_water_snr,metabolite=achieved_metabolite_snr))
    return source


def generate_job(case,root):
    torch.set_num_threads(1)
    from functools import lru_cache
    from fsl_mrs.utils import mrs_io
    from simulation import basis,legacy_water_reference
    global _read_basis_cached, _load_cached
    if '_read_basis_cached' not in globals():
        _read_basis_cached=lru_cache(maxsize=12)(mrs_io.read_basis)
        _load_cached=lru_cache(maxsize=32)(lambda path,kw:load_responses(path,**dict(kw)))
    def read_cached(path):return deepcopy(_read_basis_cached(str(path)))
    def load_cached(path,**kw):return deepcopy(_load_cached(str(path),tuple(sorted(kw.items()))))
    case_id=hashlib.sha256(json.dumps(case,sort_keys=True).encode()).hexdigest()[:16]
    directory=root/'cases'/case_id;marker=directory/'generated.json'
    try:
        if marker.exists():source=Path(json.loads(marker.read_text())['source'])
        else:
            with patch.object(basis,'read_basis',read_cached),patch.object(mrs_io,'read_basis',read_cached),patch.object(legacy_water_reference,'load_responses',load_cached):
                source=generate(case,directory)
            save(marker,{'source':str(source.resolve()),'case':case})
        return case,case_id,str(source.resolve()),''
    except Exception as exc:
        directory.mkdir(parents=True,exist_ok=True);save(directory/'error.json',{'error':repr(exc)})
        return case,case_id,'',repr(exc)


def fit_job(job):
    case_id,fitter,source,basis,out,truth=job
    output=Path(out);result_path=output/'study_result.json'
    try:
        if result_path.exists():res=json.loads(result_path.read_text())
        else:
            if fitter=='fsl':from .fsl_adapter import fit_case
            else:from .lcmodel_adapter import fit_case
            res=fit_case(Path(source),Path(basis),output)
            output.mkdir(parents=True,exist_ok=True);save(result_path,res)
        values=groups({{'Ins':'mI','Scyllo':'sI'}.get(n,n):v for n,v in res.get('estimates_mM',{}).items()}) if res.get('status') in ('ok','optimizer_warning') else {}
        return [dict(case_id=case_id,fitter=fitter,metabolite=name,truth_mM=gt,estimate_mM=values.get(name,np.nan),status='ok' if name in values else 'failed') for name,gt in truth.items()]
    except Exception as exc:
        output.mkdir(parents=True,exist_ok=True);save(output/'driver_error.json',{'error':repr(exc)})
        return [dict(case_id=case_id,fitter=fitter,metabolite=n,truth_mM=v,estimate_mM=np.nan,status='failed') for n,v in truth.items()]


def run(root,bases,limit=None,stages=None,fitters=('fsl','lcmodel')):
    torch.set_num_threads(1)
    root.mkdir(parents=True,exist_ok=True)
    cases=design(bases)
    if stages:cases=[c for c in cases if c['stage'] in stages]
    if limit:cases=cases[:limit]
    fingerprint={p.name:digest(p) for p in Path(__file__).parent.glob('*adapter.py')}
    fingerprint['run.py']=digest(__file__)
    fingerprint['simulator']={str(p):digest(p) for p in [Path('simulation/headless.py'),Path('simulation/signalModel.py'),Path('simulation/legacy_water_reference.py'),Path('simulation/basis.py'),Path('data/metabolites/metab_df.csv')]}
    fingerprint['cases']=hashlib.sha256(json.dumps(cases,sort_keys=True).encode()).hexdigest()
    manifest=dict(driver_sha256=fingerprint,commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),basis_hashes={str(p):digest(p) for p in bases.glob('*.BASIS')},cases=len(cases))
    old=root/'manifest.json'
    if old.exists() and json.loads(old.read_text())!=manifest:raise ValueError('Manifest changed; use a new output directory')
    save(old,manifest)
    rows=[];estimates=[];jobs=[]
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import multiprocessing
    with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(generate_job,c,root) for c in cases]
        for i,future in enumerate(as_completed(futures)):
            c,case_id,source,error=future.result()
            rows.append({**c,'case_id':case_id,'source':source,'generation_error':error})
            if source:
                truth=groups(json.loads((Path(source)/'ground_truth.json').read_text())['voxel_weighted'])
                jobs.extend((case_id,fitter,source,c['basis_path'],str(root/'cases'/case_id/fitter),truth) for fitter in fitters)
            pd.DataFrame(rows).to_csv(root/'cases.csv',index=False)
            print('GENERATE',i+1,'/',len(cases),error or case_id,flush=True)
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import multiprocessing
    with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn')) as pool:
        futures=[pool.submit(fit_job,job) for job in jobs]
        for j,future in enumerate(as_completed(futures)):
            estimates.extend(future.result())
            pd.DataFrame(estimates).to_csv(root/'estimates.csv',index=False)
            print('FIT',j+1,'/',len(jobs),flush=True)
    from .report import report
    report(root)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--bases',type=Path,default=BASES);p.add_argument('--limit',type=int);p.add_argument('--stages',nargs='+');p.add_argument('--fitters',nargs='+',choices=['fsl','lcmodel'],default=['fsl','lcmodel']);a=p.parse_args()
    run(a.output.resolve(),a.bases.resolve(),a.limit,a.stages,a.fitters)

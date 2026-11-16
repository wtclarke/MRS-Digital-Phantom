"""Controlled old/new TE97 mapping comparison; original acquisitions remain untouched."""
import ast,json
from pathlib import Path
from types import SimpleNamespace,MethodType
import numpy as np
import pandas as pd
import torch
from nifti_mrs.create_nmrs import gen_nifti_mrs
from simulation.basis import Basis
from simulation.signalModel import SignalModel


def run(source,output):
    output.mkdir(parents=True,exist_ok=False)
    config=json.loads((source/'simulation_config.json').read_text());pro=json.loads((source/'calibration_provenance.json').read_text())
    truth=json.loads((source/'ground_truth.json').read_text());df=pd.read_csv(source/'bHG04.csv')
    basis_path=source/'SH_PRESS_TE97.BASIS'
    basis=Basis(str(basis_path),df,bw=2000,points=2048,TE=97,water_reference=pro)
    code=ast.parse(Path('tests/fixtures/legacy_pre_quantitative/signalModel.py').read_text())
    cls=next(n for n in code.body if isinstance(n,ast.ClassDef) and n.name=='SignalModel')
    method=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='simulate_metabs')
    namespace={'torch':torch,'np':np};exec(compile(ast.Module(body=[method],type_ignores=[]),'historical_method','exec'),namespace)
    results={};state=torch.random.get_rng_state()
    try:
        for label in ['historical','corrected']:
            m=SignalModel.__new__(SignalModel);m.basis=basis;m.metab_df=df;m.labels=['Background','WM','GM','CSF'];m.gui=None;m.verbose=False
            if label=='historical':m.simulate_metabs=MethodType(namespace['simulate_metabs'],m)
            torch.manual_seed(1);spec=m.simulate_metabs().numpy()
            sampled=m.sampled_concs_per_tissue
            for tissue,values in truth['sampled_per_tissue'].items():
                if tissue in m.labels:
                    for name,value in values.items():assert sampled[tissue][name]==value
            fractions=truth['tissue_fractions'];fid=sum(fractions[t]*np.fft.ifft(np.fft.ifftshift(spec[i])) for i,t in enumerate(m.labels))*8000
            if label=='corrected':
                expected=np.zeros(2048,complex)
                for t in m.labels:
                    rows=df[df.Tissue==t].set_index('Metabolite')
                    for name,response in basis.responses_by_name.items():
                        if sampled[t][name]==0:continue
                        t2=float(rows.loc[name,'T2']);gamma=1/(np.pi*t2*.001) if np.isfinite(t2) else 0.
                        expected+=8000*fractions[t]*sampled[t][name]*response*np.exp(-gamma*basis.t)
                error=float(abs(expected-fid).max());relative=float(np.linalg.norm(expected-fid)/np.linalg.norm(expected))
                assert relative<2e-7
                results['independent_signal_relative_error']=relative
            water=8000*1000*truth['water_reference']['nominal_water_concentration_mol_l']*basis.water_fid
            dest=output/label;dest.mkdir()
            for filename,data in [('total',fid),('metabs',fid),('water_reference',water)]:
                obj=gen_nifti_mrs(data[None,None,None,:,None],.0005,basis.cf,dim_tags=['DIM_DYN',None,None])
                obj.add_hdr_field('SpecFreqChemShift',4.68);obj.add_hdr_field('EchoTime',.097);obj.save(str(dest/f'{filename}.nii.gz'))
            cfg=dict(config,path2basis=str(basis_path.resolve()))
            (dest/'simulation_config.json').write_text(json.dumps(cfg,indent=2))
            (dest/'ground_truth.json').write_text(json.dumps(truth,indent=2))
            (dest/'calibration_provenance.json').write_text(json.dumps(pro,indent=2))
        results.update(conditions='20 mm cube, saved tissue fractions/draws, zero B0/noise/MM/lipids/residual water; existing tissue linewidth retained',original_files_modified=False)
        (output/'validation.json').write_text(json.dumps(results,indent=2))
    finally:torch.random.set_rng_state(state)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.source,a.output)

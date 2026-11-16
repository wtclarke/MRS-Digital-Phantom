"""Read-only concentration/BASIS ordering diagnosis; no fitting or empirical scaling."""
import argparse,hashlib,json,re
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import torch
from fsl_mrs.utils import mrs_io
from simulation.basis import Basis
from simulation.signalModel import SignalModel


def probe():
    """Exercise the actual legacy synthesis with deliberately reversed table rows."""
    t=np.arange(256)/2000
    fids=np.column_stack([np.exp(-np.pi*t+2j*np.pi*f*t+.2j) for f in [-300,-200]])
    model=SignalModel.__new__(SignalModel);model.gui=None;model.verbose=False;model.labels=['WM']
    model.basis=SimpleNamespace(names=['NAA','Cr'],fids=fids,n=256,t=t,f=np.fft.fftfreq(256,.0005))
    model.metab_df=pd.DataFrame([dict(Tissue='WM',Metabolite=n,Conc_mean=c,Conc_std=0.,T2=np.inf) for n,c in [('Cr',2.),('NAA',8.)]])
    before=torch.random.get_rng_state().clone()
    try:
        a=model.simulate_metabs().numpy().reshape(-1)
        recorded=dict(model.sampled_concs_per_tissue['WM'])
        design=np.fft.fftshift(np.fft.fft(fids,axis=0),axes=0)
        coefficients=np.linalg.lstsq(design,a,rcond=None)[0]
        # Control: only align the table by name; all actual synthesis code is unchanged.
        model.metab_df=model.metab_df.set_index('Metabolite').loc[model.basis.names].reset_index()
        b=model.simulate_metabs().numpy().reshape(-1)
        corrected=np.linalg.lstsq(design,b,rcond=None)[0]
    finally:torch.random.set_rng_state(before)
    np.testing.assert_allclose(coefficients,[2,8],rtol=2e-7,atol=2e-7)
    np.testing.assert_allclose(corrected,[8,2],rtol=2e-7,atol=2e-7)
    assert recorded=={'Cr':2.,'NAA':8.}
    return {'metadata_concentrations':recorded,'basis_order':['NAA','Cr'],
        'actual_synthesized_coefficients':coefficients.real.tolist(),
        'name_aligned_control_coefficients':corrected.real.tolist(),
        'note':'Control changes the in-memory test table only, not production code or collaborator data.'}


def diagnose(source,report_text,output):
    output.mkdir(parents=True,exist_ok=False)
    files=[p for p in source.iterdir() if p.is_file()]
    before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    gt=json.loads((source/'ground_truth.json').read_text());pro=json.loads((source/'calibration_provenance.json').read_text())
    df=pd.read_csv(source/'bHG04.csv');fits={}
    for line in report_text.read_text().splitlines():
        m=re.search(r'([\d.E+-]+)\s+(\d+)%\s+[\d.E+-]+\s+(\S+)',line)
        if m and m[3] in gt['voxel_weighted']:fits[m[3]]=float(m[1])
    controls=probe();tables={}
    # Replay only the documented concentration draws, not spatial/noise acquisitions.
    state=torch.random.get_rng_state().clone();replay_error=0.
    try:
        torch.manual_seed(1)
        filtered=df[df.Metabolite.isin(gt['sampled_per_tissue']['WM'])]
        for tissue in ['Background','WM','GM','CSF']:
            rows=filtered[filtered.Tissue==tissue]
            values=torch.normal(torch.tensor(rows.Conc_mean.fillna(0).values),torch.tensor(rows.Conc_std.fillna(0).values))
            for name,value in zip(rows.Metabolite,values):
                replay_error=max(replay_error,abs(float(value)-gt['sampled_per_tissue'][tissue][name]))
    finally:torch.random.set_rng_state(state)
    assert replay_error==0.,'provided table/seed do not reproduce saved concentrations'

    for te in [35,97,144]:
        file=source/f'SH_PRESS_TE{te}.BASIS'
        b=Basis(str(file),df,bw=2000,points=2048,TE=te)
        orders={t:list(df[(df.Tissue==t)&df.Metabolite.isin(b.names)].Metabolite) for t in ['WM','GM','CSF']}
        assert all(v==orders['WM'] for v in orders.values())
        mapping=[]
        for signal,label in zip(b.names,orders['WM']):
            nominal=gt['voxel_weighted'][signal];assigned=gt['voxel_weighted'][label]
            mapping.append(dict(signal=signal,assigned_table_metabolite=label,nominal_mM=nominal,
                actual_assigned_mM=assigned,LCModel_reported_mM=fits.get(signal),
                fit_over_assigned=fits[signal]/assigned if signal in fits and assigned else None))
        no_bhg=[n for n in b.names if n!='bHG'];table_no_bhg=[n for n in orders['WM'] if n!='bHG']
        assert no_bhg==table_no_bhg
        original=mrs_io.read_basis(str(file));formatted=original.get_formatted_basis(2000,2048)
        errors={}
        for name in b.names:
            x=b.fids[:,b.names.index(name)];y=formatted[:,original.names.index(name)]
            c=np.vdot(y,x)/np.vdot(y,y)
            errors[name]=float(np.linalg.norm(x-c*y)/np.linalg.norm(x))
        tables[str(te)]={'basis_order':b.names,'table_order':orders['WM'],'mapping':mapping,
            'without_bhg_orders_match':True,'column_identity_errors_after_common_scale':errors,'basis_sha256':before[file.name]}
    assert tables['97']['basis_sha256']==pro['basis_sha256']
    naa=next(r for r in tables['97']['mapping'] if r['signal']=='NAA')
    assert naa['assigned_table_metabolite']=='NAAG'
    assert abs(naa['LCModel_reported_mM']/naa['actual_assigned_mM']-1)<.05
    assert abs(naa['LCModel_reported_mM']/naa['nominal_mM']-1)>.5
    assert before=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    result={'concentration_seed_replay_max_error':replay_error,'controls':controls,'bases':tables,'source_hashes':before,'source_files_unchanged':True,
        'NAA_nominal_to_assigned_ratio':naa['nominal_mM']/naa['actual_assigned_mM'],
        'LCModel_fits_rerun':False,'missing_inputs':['output.RAW','water.RAW','TABLE/PRINT/COORD','original NIfTI-MRS signals']}
    (output/'diagnosis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'controls':controls,'NAA':naa,'nominal_to_assigned':result['NAA_nominal_to_assigned_ratio'],
        'all_three_TEs_affected_with_supplied_csv':True,'source_files_unchanged':True},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--report-text',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();diagnose(a.source,a.report_text,a.output)

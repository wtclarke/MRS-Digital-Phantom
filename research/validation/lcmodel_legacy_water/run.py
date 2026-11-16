"""Fit immutable saved legacy acquisitions with the earlier LCModel protocol."""
import argparse,csv,hashlib,json,os,re,shlex,subprocess,time
from pathlib import Path
import numpy as np
from fsl_mrs.utils import mrs_io
from utils.lcmodel_quantitative_io import compatible_basis,control_text,parse_basis,parse_raw,write_raw,parse_coord


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')


def run(source,output,binary,basis=None):
    source=source.resolve();output=output.resolve();output.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((source/'simulation_config.json').read_text());truth=json.loads((source/'ground_truth.json').read_text())
    provenance=json.loads((source/'calibration_provenance.json').read_text())
    basis=Path(basis or cfg['path2basis']).resolve()
    if digest(basis)!=provenance['basis_sha256']:raise ValueError('supplied basis differs from acquisition')
    if not binary.is_file() or not os.access(binary,os.X_OK):raise RuntimeError('licensed LCModel binary missing/not executable')
    protected=[p for p in source.iterdir() if p.is_file()]+[basis,binary]
    before={str(p):digest(p) for p in protected};save(output/'input_hashes.json',before)
    adapted=output/'lcmodel_compatible.BASIS'
    conversion=compatible_basis(basis,adapted,data_origin_ppm=provenance['frequency_origin_ppm'])
    entries,_=parse_basis(basis);copies,_=parse_basis(adapted)
    for a,b in zip(entries,copies):np.testing.assert_array_equal(a[2],b[2])
    nominal=truth['nominal_voxel_average_metabolite_concentrations_mM'];names=list(nominal)
    if not set(names)<=set(e[0] for e in entries):raise ValueError('ground-truth names missing from basis')
    readers={n:mrs_io.read_FID(str(source/f'{n}.nii.gz')) for n in ['total','metabs','water_reference']}
    origin=provenance['frequency_origin_ppm'];cf=readers['total'].spectrometer_frequency[0];dt=1/cfg['simulation_params']['bandwidth']
    for reader in readers.values():
        assert np.isclose(reader.dwelltime,dt,rtol=1e-7,atol=0) and reader.spectrometer_frequency[0]==cf and reader.SpecFreqChemShift==origin
    water_mol_l=truth['water_reference']['nominal_water_concentration_mol_l']
    rows=[];records=[]
    for label,acquisition,cw in [('total_previous_water','total',39.0165),('total_voxel_water','total',water_mol_l),('metabs_voxel_water','metabs',water_mol_l)]:
        directory=output/label;directory.mkdir();fits=directory/'lcmodel_outputs';fits.mkdir()
        errors={}
        for name,filename in [(acquisition,'metabolites.RAW'),('water_reference','water_reference.H2O')]:
            fid=readers[name].mrs().FID
            write_raw(directory/filename,fid,frequency_mhz=cf,echo_time_ms=cfg['basis_set_settings']['TE'])
            stored,_=parse_raw(directory/filename);np.testing.assert_array_equal(stored.conj(),fid)
            errors[name]=float(abs(stored.conj()-fid).max())
        control=control_text(directory,adapted,names=[e[0] for e in entries],points=len(fid),dwell_time_s=dt,
            frequency_mhz=cf,origin_ppm=origin,water_mol_l=cw)
        control=control.replace('ECHOT=30.',f"ECHOT={float(cfg['basis_set_settings']['TE']):g}.").replace('Phase7b exact Phase7a noise-free acquisition',f"Legacy saved {acquisition}, TE{float(cfg['basis_set_settings']['TE']):g}")
        omit=[e[0] for e in entries if e[0] not in names]
        control=control[:control.index(' NSIMUL=')]+f' NSIMUL=0, NRATIO=0, NOMIT={len(omit)}\n'+''.join(f" CHOMIT({i})='{n}'\n" for i,n in enumerate(omit,1))+' $END\n'
        (directory/'lcmodel.CONTROL').write_text(control)
        started=time.monotonic()
        p=subprocess.run([str(binary)],input=control,text=True,cwd=directory,capture_output=True,timeout=300)
        execution={'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr,'elapsed_s':time.monotonic()-started}
        save(directory/'execution.json',execution)
        table=(fits/'fit.TABLE').read_text() if (fits/'fit.TABLE').exists() else ''
        if p.returncode or 'FATAL' in table or not table:raise RuntimeError(f'LCModel failure in {label}: {execution}; {table}')
        values={}
        for line in table.splitlines():
            match=re.match(r'^\s*([\d.Ee+-]+)\s+(\d+(?:\.\d+)?)%\s+([\d.Ee+-]+)\s+(\S+)\s*$',line)
            if match:values[match[4]]={'mM':float(match[1]),'percent_SD':float(match[2])}
        if not set(names)<=values.keys():raise ValueError('incomplete LCModel concentration table')
        coord=parse_coord(fits/'fit.COORD');residual=float(np.linalg.norm(coord['data']-coord['fit'])/np.linalg.norm(coord['data']))
        printed=(fits/'fit.PRINT').read_text()
        scaling={k:float(m[1]) for k,pattern in {'FCALIB':r'FCALIB =\s*([\d.E+-]+)','water_area':r'Area of unsuppressed water peak =\s*([\d.E+-]+)','reference_basis_area':r'Normalized area of reference Basis singlet =\s*([\d.E+-]+)'}.items() if (m:=re.search(pattern,printed))}
        diagnostics=[line.strip() for line in table.splitlines() if re.search(r'\b(fatal|warning|info|error)\b',line,re.I)]
        record={'case':label,'input_acquisition':acquisition,'water_mol_l':cw,'concentrations':values,
            'relative_phased_real_residual':residual,'scaling':scaling,'diagnostics':diagnostics,
            'conversion_max_complex_errors':errors,'execution':execution}
        save(directory/'fitting_results.json',record);records.append(record)
        for name,gt in nominal.items():
            v=values[name]
            rows.append(dict(case=label,metabolite=name,nominal_voxel_mM=gt,
                tissue_only_mM=truth['nominal_tissue_only_metabolite_concentrations_mM'][name],LCModel_mM=v['mM'],
                absolute_error_mM=v['mM']-gt,relative_error_percent=100*(v['mM']-gt)/gt if gt else None,percent_SD=v['percent_SD']))
        subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-sDEVICE=pdfwrite',f'-sOutputFile={fits/"fit.pdf"}',str(fits/'fit.PS')],check=True,capture_output=True)
    with (output/'comparison.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    assert before=={str(p):digest(p) for p in protected},'original inputs changed'
    summary={'sampling_metadata':{k:float(v.dwelltime) for k,v in readers.items()},'nominal_dwell_time_s':dt,'cases':records,'truth':truth,'basis_compatibility':conversion,'inputs_unchanged':True,
        'protocol_changes':['matching supplied GE basis','actual TE35','all 18 simulated metabolites','case-specific WCONC','paths/title'],
        'commands':[f'cd {shlex.quote(str(output/r["case"]))} && {shlex.quote(str(binary))} < lcmodel.CONTROL' for r in records]}
    save(output/'summary.json',summary)
    (output/'commands.txt').write_text('\n'.join(summary['commands'])+'\n')
    print(json.dumps({'cases':[{k:r[k] for k in ['case','relative_phased_real_residual','scaling','diagnostics']} for r in records],'comparison':rows},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--basis',type=Path)
    a=p.parse_args();run(a.source,a.output,a.binary,a.basis)

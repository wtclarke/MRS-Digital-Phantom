"""Phase 6e: analytical and FSL-MRS validation of pre-acquisition relaxation."""
import argparse
import csv
import json
from pathlib import Path
import shlex
import subprocess
import warnings

import fsl_mrs
import numpy as np
from fsl_mrs.utils import fitting, mrs_io, quantify
from fsl_mrs.utils.constants import H2O_MOLARITY

from examples.fsl_mrs_quantification_scaling_derivation import spectral_integral, relaxation as independent_R
from examples.fsl_mrs_water_fit_diagnostics import fit_water
from examples.mrscloud_quantitative_voxel import BUNDLE, METABOLITES, WATER
from simulation.mrscloud_full_basis import load_library, sha256
from simulation.mrscloud_quantitative_pair import export_pair
from simulation.mrscloud_quantitative_voxel import simulate_voxel

NAMES = list(METABOLITES['WM'])


def configuration():
    """Declared illustrative seconds, not a new literature parameter database."""
    water = {'WM':(.97,.073),'GM':(1.5,.088),'CSF':(4.47,2.03)}
    metabolites = {'WM':{'NAA':(1.29,.194),'Cr':(1.1,.15),'GPC':(1.2,.2),'PCh':(1.1,.18),'Glu':(1.15,.16)},
                   'GM':{'NAA':(1.45,.18),'Cr':(1.35,.17),'GPC':(1.3,.22),'PCh':(1.25,.21),'Glu':(1.4,.19)},
                   'CSF':{n:(1.5,.25) for n in METABOLITES['CSF']}}
    return dict(mode='tissue_specific',TE_s=.03,TR_s=2.,
                water={t:dict(T1_s=a,T2_s=b) for t,(a,b) in water.items()},
                metabolites={t:{n:dict(T1_s=a,T2_s=b) for n,(a,b) in profile.items()} for t,profile in metabolites.items()})


def quantification_info(mrs, truth, cfg, *, matched):
    fractions = truth['tissue_fractions']
    kwargs = {}
    if matched:
        # FSL has one METAB pair: use actual NAA times in pure tissue and an
        # explicitly approximate WM/GM volume mean in mixtures, never a fitted factor.
        tissue_sum = fractions['WM']+fractions['GM']
        times = {key:sum(fractions[t]*cfg['metabolites'][t]['NAA'][key] for t in ['WM','GM'])/tissue_sum
                 for key in ['T1_s','T2_s']}
        kwargs = dict(t1={**{'H2O_'+t:cfg['water'][t]['T1_s'] for t in fractions},'METAB':times['T1_s']},
                      t2={**{'H2O_'+t:cfg['water'][t]['T2_s'] for t in fractions},'METAB':times['T2_s']})
    q = quantify.QuantificationInfo(cfg['TE_s'],cfg['TR_s'],mrs.names,mrs.centralFrequency,
        water_ref_metab='NAA',water_ref_metab_protons=3,water_ref_metab_limits=(1.8,2.2),**kwargs)
    q.set_fractions(fractions)
    if matched:
        q.set_densities({t:c*1000/H2O_MOLARITY for t,c in truth['tissue_specific_water_concentrations_mol_l'].items()})
    return q


def predict(truth, cfg, responses, water, q):
    """Use nominal profiles and input relaxation times, never fitted amplitudes."""
    f = truth['tissue_fractions']
    effective = {}; effective_water = 0.
    for n in NAMES:
        effective[n] = sum(f[t]*truth['tissue_specific_metabolite_concentrations_mM'][t][n]*
            independent_R(cfg['TE_s'],cfg['TR_s'],cfg['metabolites'][t][n]['T1_s'],cfg['metabolites'][t][n]['T2_s']) for t in f)
    for t in f:
        effective_water += f[t]*truth['tissue_specific_water_concentrations_mol_l'][t]*independent_R(
            cfg['TE_s'],cfg['TR_s'],cfg['water'][t]['T1_s'],cfg['water'][t]['T2_s'])
    t = np.arange(len(water))/2000
    ppm = np.fft.fftshift(np.fft.fftfreq(len(t),.0005))/127.731+4.68
    jm = spectral_integral(responses['NAA'],ppm,(1.8,2.2))
    jw = spectral_integral(abs(water[0])*np.exp(-np.pi*t),ppm,(1.65,7.65))
    rm = independent_R(q.te,q.tr,q.t1['METAB'],q.t2['METAB'])
    density = dict(WM=q.d_WM,GM=q.d_GM,CSF=q.d_CSF)
    sum_water = H2O_MOLARITY*sum(f[t]*density[t]*independent_R(q.te,q.tr,q.t1['H2O_'+t],q.t2['H2O_'+t]) for t in f)
    multiplier = (jm/jw)*(2/3)*sum_water/rm/(1-f['CSF'])/effective_water/1000
    return dict(predicted_mM={n:effective[n]*multiplier for n in NAMES},
        effective_metabolite_mM=effective,effective_water_M=effective_water,
        NAA_unit_integral=jm,water_unit_integral=jw,response_factor=jm/jw*2/3,
        quantification_water_mM=sum_water,quantification_metabolite_factor=rm,csf_correction=1/(1-f['CSF']),
        quantification_T1_s=q.t1,quantification_T2_s=q.t2,quantification_densities=density)


def fit_case(directory, responses, water):
    truth=json.loads((directory/'ground_truth.json').read_text())
    config=json.loads((directory/'simulation_config.json').read_text())
    cfg=config['relaxation'];gain=config['receiver_factor']
    mrs=mrs_io.read_FID(str(directory/'metabolites.nii.gz')).mrs(
        basis=mrs_io.read_basis(str(BUNDLE/'metabolite_basis.BASIS')),ref_data=str(directory/'water_reference.nii.gz'))
    mrs.keep=NAMES
    water_fit,_,_=fit_water(mrs)
    expected_water=float(abs(water[0])*truth['relaxation_weighted_water_concentration_mol_l']*gain)
    water_fit.update(expected_amplitude=expected_water,relative_amplitude_error=water_fit['fitted_amplitude']/expected_water-1)
    result=dict(truth=truth,water_fit=water_fit)
    if not np.any(mrs.FID):
        result['status']='pure CSF: zero metabolite FID; tissue-only and water-referenced quantification undefined'
        return result
    mrs.rescaleForFitting()
    with warnings.catch_warnings(record=True) as caught:
        fit=fitting.fit_FSLModel(mrs,method='Newton',model='lorentzian',baseline_order=-1,
            ppmlim=(.2,4.2),metab_groups=[0]*len(NAMES),capture_minimize_output=True)
    sd,sb=mrs.scaling['FID'],mrs.scaling['basis']
    result.update(parameters=dict(zip(fit.fitResults.columns,fit.fitResults.iloc[0].tolist())),
        optimizer_success=bool(fit.scipy_minimize_result.success),optimizer_message=str(fit.scipy_minimize_result.message),
        relative_complex_residual=float(np.linalg.norm(mrs.FID-fit.predictedFID(mrs))/np.linalg.norm(mrs.FID)),
        scaling=mrs.scaling,warnings=[str(w.message) for w in caught],coefficients={},quantification={})
    for n in NAMES:
        expected=truth['relaxation_weighted_metabolite_concentrations_mM'][n]*.001*gain*sd/sb
        fitted=float(fit.getConc(metab=n))
        result['coefficients'][n]=dict(expected_internal=expected,fitted_internal=fitted,
            fitted_effective_mM=fitted*sb/sd/gain*1000,relative_error=fitted/expected-1)
    original=fit.params.copy()
    for label,matched in [('CLI_defaults',False),('API_declared_water_NAA_scalar',True)]:
        q=quantification_info(mrs,truth,cfg,matched=matched)
        expected=predict(truth,cfg,responses,water,q)
        fit.calculateConcScaling(mrs,quant_info=q,internal_reference=['NAA'])
        np.testing.assert_array_equal(fit.params,original)
        result['quantification'][label]=dict(prediction=expected,
            reported_mM={n:float(fit.getConc(scaling='molarity',metab=n)) for n in NAMES})
    np.savez(directory/'fit_arrays.npz',fitted=fit.predictedFID(mrs)/sd,input=mrs.FID/sd)
    return result


def run_validation(output_dir):
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=False)
    before={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}
    responses,water,metadata=load_library(BUNDLE,selected_names=NAMES)
    cfg=configuration();cases={};commands=[];analytical=[]
    for name,f in [('pure_WM',{'WM':1}),('pure_GM',{'GM':1}),('pure_CSF',{'CSF':1}),('mixed',{'WM':.2,'GM':.3,'CSF':.5})]:
        pair=simulate_voxel(f,calibration_bundle=BUNDLE,metabolite_concentrations_mM=METABOLITES,
            water_concentrations_mol_l=WATER,relaxation=cfg)
        # Independently sum the tissue signals using scalar math and unit arrays.
        expected_m=np.zeros(2048,complex);expected_w=np.zeros(2048,complex)
        for tissue,fraction in f.items():
            for n,c in METABOLITES[tissue].items():
                p=cfg['metabolites'][tissue][n]
                expected_m+=fraction*c*.001*independent_R(.03,2.,p['T1_s'],p['T2_s'])*responses[n]
            p=cfg['water'][tissue]
            expected_w+=fraction*WATER[tissue]*independent_R(.03,2.,p['T1_s'],p['T2_s'])*water
        directory=out/name;paths=export_pair(pair,directory)
        analytical.append(dict(case=name,metabolite_max_error=float(np.max(abs(pair['metabolite_fid']-expected_m))),
            water_max_error=float(np.max(abs(pair['water_reference_fid']-expected_w))),
            export_max_error=max(float(np.max(abs(mrs_io.read_FID(str(path)).mrs().FID-pair[key]))) for key,path in paths.items())))
        cases[name]=fit_case(directory,responses,water)
        if name=='pure_CSF':continue
        fractions=pair['ground_truth']['tissue_fractions']
        command=['pixi','run','fsl_mrs','--data',str(directory/'metabolites.nii.gz'),
            '--h2o',str(directory/'water_reference.nii.gz'),'--basis',str(BUNDLE/'metabolite_basis.BASIS'),
            '--output',str(directory/'fsl_mrs_results'),'--keep',*NAMES,'--algo','Newton','--lorentzian',
            '--baseline','off','--ppmlim','0.2','4.2','--no_conj_basis','--TE','30','--TR','2',
            '--tissue_frac',*[str(fractions[t]) for t in ['WM','GM','CSF']],
            '--internal_ref','NAA','--wref_metabolite','NAA','--ref_protons','3','--ref_int_limits','1.8','2.2','--report']
        commands.append(shlex.join(command));(out/'cli_commands.txt').write_text('\n'.join(commands)+'\n')
        completed=subprocess.run(command,capture_output=True,text=True)
        (directory/'cli.log').write_text(completed.stdout+completed.stderr);completed.check_returncode()
        if (directory/'fsl_mrs_results/report.html').stat().st_size<1000:
            raise RuntimeError('CLI HTML report missing or empty')
        with (directory/'fsl_mrs_results/concentrations.csv').open() as stream:
            records=list(csv.reader(stream))
        cases[name]['cli_mM']={row[0]:float(row[4]) for row in records[2:] if row[0] in NAMES}
    rows=[]
    for name,data in cases.items():
        if 'coefficients' not in data:continue
        gt=data['truth']
        for n,coef in data['coefficients'].items():
            for setting,q in data['quantification'].items():
                predicted=q['prediction']['predicted_mM'][n];reported=q['reported_mM'][n]
                rows.append(dict(case=name,metabolite=n,quantification=setting,
                    nominal_voxel_mM=gt['voxel_average_metabolite_concentrations_mM'][n],
                    nominal_tissue_only_mM=gt['tissue_only_metabolite_concentrations_mM'][n],
                    weighted_effective_mM=gt['relaxation_weighted_metabolite_concentrations_mM'][n],**coef,
                    reported_mM=reported,predicted_mM=predicted,prediction_error_mM=reported-predicted))
    with (out/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    if before!={p.name:sha256(p) for p in BUNDLE.iterdir() if p.is_file()}:raise RuntimeError('calibration bundle changed')
    result=dict(FSL_MRS_version=fsl_mrs.__version__,configuration=cfg,cases=cases,analytical=analytical,
        calibration_identity=metadata['calibration_identity'],bundle_hashes=before,bundle_unchanged=True,
        quantification_policy='CLI defaults; API exact water profiles and scalar NAA T1/T2 (WM/GM volume means for mixtures)',
        fitting_settings=dict(method='Newton',model='lorentzian',baseline='off',ppm=[.2,4.2],metabolites=NAMES,shared_group=True))
    (out/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output_dir')
    result=run_validation(parser.parse_args().output_dir)
    print(json.dumps(result['analytical'],indent=2))

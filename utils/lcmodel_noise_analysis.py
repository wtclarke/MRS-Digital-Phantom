"""Phase 7c statistics and output parsing; no fitting or signal generation."""
from pathlib import Path
import re

import numpy as np

from utils.lcmodel_quantitative_io import NAMES,parse_coord,parse_table

METABOLITES=(*NAMES,'GPC+PCh')


def parse_fit(directory):
    directory=Path(directory);table=directory/'fit.TABLE'
    result=parse_table(table)
    if result['fatal']:return dict(result,valid=False,failure='fatal diagnostic')
    if set(result['concentrations'])!=set(NAMES):
        return dict(result,valid=False,failure='missing final concentrations')
    text=table.read_text();raw={}
    for line in text.splitlines():
        m=re.match(r'^\s*([\d.Ee+-]+)\s+(\d+(?:\.\d+)?)%\s+([\d.Ee+-]+)\s+(\S+)\s*$',line)
        if m and m[4] in METABOLITES:
            raw[m[4]]=dict(concentration_text=m[1],percent_SD_text=m[2],concentration=float(m[1]),percent_SD=float(m[2]))
    if any(not np.isfinite(v['concentration']) or v['concentration']<0 for v in raw.values()):
        return dict(result,valid=False,failure='nonfinite or negative final concentrations')
    group=sum(raw[n]['concentration'] for n in ['GPC','PCh'])
    group_report=raw.get('GPC+PCh')
    # Use the actual reported group %SD, never independent-error quadrature.
    values={n:dict(v) for n,v in raw.items()}
    values['GPC+PCh']=dict(concentration=group,percent_SD=None,
        reported_group_concentration=group_report['concentration'] if group_report else None,
        reported_group_percent_SD=group_report['percent_SD'] if group_report else None)
    for name,v in values.items():
        c=group_report['concentration'] if name=='GPC+PCh' and group_report else v['concentration']
        p=group_report['percent_SD'] if name=='GPC+PCh' and group_report else v['percent_SD']
        v['displayed_SD_product_mM']=c*p/100 if p is not None else None
        v['uncertainty_mM']=c*p/100 if p is not None and c>0 and 0<p<999 else None
        v['zero_percent_SD']=p==0 if p is not None else None
    coords=parse_coord(directory/'fit.COORD');printed=(directory/'fit.PRINT').read_text()
    def last(pattern):
        found=re.findall(pattern,printed)
        return float(found[-1]) if found else None
    shape=re.findall(r'(\d+) spline knots\.\s+Ns =\s*(\d+)\((\d+)\)',text)
    phase=re.search(r'Ph:\s*([-\d.]+) deg\s+([-\d.]+) deg/ppm',text)
    result.update(valid=True,values=values,raw_concentration_tokens=raw,
        reported_SNR=last(r'S/N =\s*([-\d.]+)'),FWHM_ppm=last(r'FWHM =\s*([\d.E+-]+)'),
        data_shift_ppm=last(r'Data shift =\s*([-\d.E+]+) ppm'),
        phase_zero_deg=float(phase[1]) if phase else None,phase_first_deg_ppm=float(phase[2]) if phase else None,
        relative_real_residual=float(np.linalg.norm(coords['data']-coords['fit'])/np.linalg.norm(coords['data'])),
        baseline_relative_norm=float(np.linalg.norm(coords['baseline'])/np.linalg.norm(coords['data'])),
        FCALIB=last(r'FCALIB =\s*([\d.E+-]+)'),water_area=last(r'Area of unsuppressed water peak =\s*([\d.E+-]+)'),
        baseline_knots=int(shape[-1][0]) if shape else None,
        convergence_diagnostic=any(d['routine'] in ['RFALSI','SOLVE'] or (d['routine']=='TWOREG' and d['code'] in [6,7,10]) for d in result['diagnostics']))
    return result


def descriptive(values,truth,baseline,uncertainties,total,failures):
    values=np.asarray(values,dtype=float);unc=[u for u in uncertainties if u is not None and np.isfinite(u)]
    if not len(values):return dict(count=0,total=total,failure_rate=failures/total)
    error=values-truth;noise=values-baseline
    return dict(count=len(values),total=total,failure_rate=failures/total,
        mean_mM=float(values.mean()),nominal_bias_mM=float(error.mean()),noise_bias_mM=float(noise.mean()),
        SD_mM=float(values.std(ddof=1)) if len(values)>1 else None,
        nominal_RMSE_mM=float(np.sqrt(np.mean(error**2))),noise_RMSE_mM=float(np.sqrt(np.mean(noise**2))),
        nominal_median_absolute_error_mM=float(np.median(abs(error))),noise_median_absolute_error_mM=float(np.median(abs(noise))),
        median_reported_uncertainty_mM=float(np.median(unc)) if unc else None,uncertainty_valid_count=len(unc),
        nominal_poor_fraction=float(np.mean(abs(error)/truth>.2)),noise_poor_fraction=float(np.mean(abs(noise)/baseline>.2)))


def paired_bootstrap(lcm,fsl,lcm_baseline,fsl_baseline,truth,*,seed=77001,draws=10000):
    """Percentile intervals for paired mean differences; no significance test."""
    lcm=np.asarray(lcm,float);fsl=np.asarray(fsl,float)
    if len(lcm)!=len(fsl) or not len(lcm):raise ValueError('paired arrays must have equal nonzero length')
    a=lcm-lcm_baseline;b=fsl-fsl_baseline
    deltas=dict(nominal_absolute_error=abs(lcm-truth)-abs(fsl-truth),
        nominal_squared_error=(lcm-truth)**2-(fsl-truth)**2,
        noise_absolute_error=abs(a)-abs(b),noise_squared_error=a*a-b*b,
        baseline_centered_estimate=a-b,
        noise_fractional_absolute_error=abs(a/lcm_baseline)-abs(b/fsl_baseline),
        noise_fractional_squared_error=(a/lcm_baseline)**2-(b/fsl_baseline)**2)
    indices=np.random.default_rng(seed).integers(0,len(lcm),size=(draws,len(lcm)))
    output=dict(pairs=len(lcm),bootstrap_seed=seed,bootstrap_draws=draws)
    for name,delta in deltas.items():
        means=delta[indices].mean(axis=1);low,high=np.quantile(means,[.025,.975])
        output.update({name+'_mean':float(delta.mean()),name+'_CI_low':float(low),name+'_CI_high':float(high)})
    return output

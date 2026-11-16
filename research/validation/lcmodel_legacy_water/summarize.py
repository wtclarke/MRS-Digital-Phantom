"""Summarize saved LCModel results without rerunning fitting."""
import argparse,csv,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from utils.lcmodel_quantitative_io import parse_coord


def summarize(root):
    result=json.loads((root/'summary.json').read_text());truth=result['truth']['voxel_weighted'];names=list(truth)
    grouped=[]
    for case in result['cases']:
        values=case['concentrations']
        for group,members in [('tCr',['Cr','PCr']),('tCho',['GPC','PCh']),('tNAA',['NAA','NAAG'])]:
            expected=sum(truth[n] for n in members);actual=sum(values[n]['mM'] for n in members)
            grouped.append(dict(case=case['case'],group=group,nominal_mM=expected,LCModel_mM=actual,
                relative_error_percent=100*(actual/expected-1)))
        scale=case['scaling'];case['FCALIB_from_printed_areas']=scale['reference_basis_area']*2*case['water_mol_l']*1000/scale['water_area']
    with (root/'grouped_comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=grouped[0]);writer.writeheader();writer.writerows(grouped)
    fig,axes=plt.subplots(2,1,figsize=(12,8),layout='constrained')
    x=np.arange(len(names));width=.25
    axes[0].bar(x-width,[truth[n] for n in names],width,label='Nominal voxel truth')
    for i,case in enumerate(result['cases'][1:]):axes[0].bar(x+i*width,[case['concentrations'][n]['mM'] for n in names],width,label=case['case'])
    axes[0].set_xticks(x,names,rotation=45);axes[0].set_ylabel('Nominal / LCModel-reported mM');axes[0].legend()
    for case in result['cases'][1:]:axes[1].plot(x,[100*(case['concentrations'][n]['mM']/truth[n]-1) for n in names],'o-',label=case['case'])
    axes[1].axhline(0,color='black',lw=.8);axes[1].set_xticks(x,names,rotation=45);axes[1].set_ylabel('Difference from nominal truth (%)');axes[1].legend()
    fig.savefig(root/'concentration_comparison.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,9),layout='constrained')
    for ax,case in zip(axes,result['cases']):
        coord=parse_coord(root/case['case']/'lcmodel_outputs/fit.COORD')
        for key in ['data','fit','baseline']:ax.plot(coord['ppm'],coord[key],label=key)
        ax.plot(coord['ppm'],coord['data']-coord['fit'],label='residual')
        ax.set_xlim(4.2,.2);ax.set_title(case['case']);ax.set_ylabel('LCModel report units');ax.legend()
    axes[-1].set_xlabel('Chemical shift (ppm)');fig.savefig(root/'fit_comparison.png',dpi=150);plt.close(fig)
    (root/'analysis_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Comparisons and plots saved:',root)
    for case in result['cases']:print(case['case'],'FCALIB printed/reconstructed',case['scaling']['FCALIB'],case['FCALIB_from_printed_areas'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);summarize(p.parse_args().output)

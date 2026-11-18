"""Summarize paired headless-water fitting experiments without acceptance thresholds."""
from pathlib import Path
import json
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress


CASE_COLUMNS = ['stage', 'basis_name', 'tissue', 'target', 'scale', 'snr', 'b0_sd', 'water_snr']


def summarize(data, keys):
    """Count all attempts; compute errors only for successful finite estimates."""
    records = []
    for values, rows in data.groupby(keys, dropna=False, sort=True):
        if not isinstance(values, tuple):
            values = (values,)
        good = rows.loc[rows['valid']]
        error = good['estimate_mM'] - good['truth_mM']
        relative = 100 * error / good['truth_mM'].replace(0, np.nan)
        records.append(dict(zip(keys, values), attempted=len(rows), successful=len(good),
            failed=len(rows)-len(good), failure_rate=1-len(good)/len(rows),
            bias_mM=error.mean(), rmse_mM=np.sqrt((error**2).mean()),
            mean_mM=good['estimate_mM'].mean(), sd_mM=good['estimate_mM'].std(),
            truth_mean_mM=good['truth_mM'].mean(), bias_percent=relative.mean()))
    return pd.DataFrame(records, columns=keys + ['attempted', 'successful', 'failed',
        'failure_rate', 'bias_mM', 'rmse_mM', 'mean_mM', 'sd_mM', 'truth_mean_mM', 'bias_percent'])


def responses(data):
    clean = data.loc[data.stage.isin(['clean', 'transfer']) &
        np.isinf(data.snr) & data.b0_sd.eq(0)].copy()
    baseline = clean.loc[clean.target.eq('baseline')].copy()
    targets = clean.loc[~clean.target.eq('baseline'), ['stage', 'basis_name', 'tissue', 'target']].drop_duplicates()
    injected = baseline.drop(columns=['target', 'stage']).merge(targets, on=['basis_name', 'tissue'])
    clean = pd.concat([clean.loc[~clean.target.eq('baseline')], injected], ignore_index=True)
    target_truth = data[['case_id', 'metabolite', 'truth_mM']].drop_duplicates(['case_id', 'metabolite'])
    clean = clean.merge(target_truth.rename(columns={'metabolite': 'target', 'truth_mM': 'target_truth_mM'}),
        on=['case_id', 'target'], how='left', validate='many_to_one')
    baseline_keys = ['basis_name', 'tissue', 'fitter', 'metabolite']
    base_truth = baseline.groupby(baseline_keys, dropna=False).truth_mM.mean().rename('baseline_truth_mM')
    base_fit = baseline.loc[baseline.valid].groupby(baseline_keys, dropna=False).estimate_mM.mean().rename('baseline_estimate_mM')
    clean = clean.join(base_truth, on=baseline_keys).join(base_fit, on=baseline_keys)
    clean['truth_delta_mM'] = clean.truth_mM-clean.baseline_truth_mM
    clean['estimate_delta_mM'] = clean.estimate_mM-clean.baseline_estimate_mM
    clean['delta_error_mM'] = (clean.estimate_delta_mM-clean.truth_delta_mM).where(clean.valid)
    clean['truth_fold_change'] = clean.truth_mM/clean.baseline_truth_mM.replace(0, np.nan)
    clean['estimate_fold_change'] = clean.estimate_mM/clean.baseline_estimate_mM.replace(0, np.nan)
    clean['fold_change_error'] = (clean.estimate_fold_change-clean.truth_fold_change).where(clean.valid)
    keys = ['stage', 'basis_name', 'tissue', 'target', 'fitter', 'metabolite']
    records = []
    for values, rows in clean.groupby(keys, dropna=False):
        good = rows.loc[rows.valid & np.isfinite(rows.target_truth_mM)]
        record = dict(zip(keys, values), attempted=len(rows), successful=len(good),
            failed=len(rows)-len(good), slope=np.nan, intercept_mM=np.nan, r_squared=np.nan,
            delta_rmse_mM=np.sqrt((good.delta_error_mM**2).mean()),
            fold_change_rmse=np.sqrt((good.fold_change_error**2).mean()))
        if good.target_truth_mM.nunique() >= 2:
            fit = linregress(good.target_truth_mM, good.estimate_mM)
            record.update(slope=fit.slope, intercept_mM=fit.intercept,
                r_squared=fit.rvalue**2 if good.estimate_mM.nunique() > 1 else np.nan)
        records.append(record)
    columns = keys + ['attempted', 'successful', 'failed', 'slope', 'intercept_mM', 'r_squared', 'delta_rmse_mM', 'fold_change_rmse']
    return clean, pd.DataFrame(records, columns=columns)


def save_figure(fig, root, name):
    fig.savefig(root / f'{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(root / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


def empty_axes(ax, text='No successful finite estimates available'):
    ax.text(.5, .5, text, ha='center', va='center', transform=ax.transAxes)
    ax.set_axis_off()


def response_figure(clean, root, mode='absolute'):
    x, y, xlabel, ylabel, filename = {
        'absolute': ('target_truth_mM', 'estimate_mM', 'Target ground truth (nominal mM)', 'Fitted estimate (mM)', 'target_response'),
        'residual': ('target_truth_mM', 'signed_error_mM', 'Target ground truth (nominal mM)', 'Signed error: fitted − ground truth (mM)', 'target_residual'),
        'delta': ('truth_delta_mM', 'estimate_delta_mM', 'Ground-truth change from baseline (mM)', 'Fitted change from baseline (mM)', 'delta_response'),
        'fold': ('truth_fold_change', 'estimate_fold_change', 'Ground-truth fold change', 'Fitted fold change', 'fold_response'),
    }[mode]
    targeted = clean.loc[clean.metabolite.eq(clean.target)]
    if mode == 'residual':
        targeted = targeted.assign(signed_error_mM=targeted.estimate_mM-targeted.truth_mM)
    panels = targeted[['basis_name', 'target']].drop_duplicates().values.tolist()

    def draw_panel(ax, basis, target):
        rows = targeted.loc[targeted.basis_name.eq(basis) & targeted.target.eq(target) & targeted.valid & np.isfinite(targeted[x]) & np.isfinite(targeted[y])]
        for values, group in rows.groupby(['tissue', 'fitter', 'stage']):
            means = group.groupby(x)[y].mean().sort_index()
            if mode == 'residual':
                ax.plot(means.index, means, marker={'fsl': 'o', 'lcmodel': 's'}[values[1]],
                        linestyle='none', markersize=3, label='/'.join(map(str, values)))
            else:
                ax.plot(means.index, means, 'o-', markersize=3, label='/'.join(map(str, values)))
        if len(rows):
            limits = [rows[x].min(), rows[x].max()]
            ax.plot(limits, [0, 0] if mode == 'residual' else limits,
                    '--', color='.5', label='Zero error' if mode == 'residual' else 'Identity')
            ax.legend(fontsize=6, ncol=2)
        else:
            empty_axes(ax)
        ax.set(title=f'{basis}: {target}', xlabel=xlabel, ylabel=ylabel)

    # Each complete basis stays on a readable page; the shared filename is an overview.
    for basis in targeted.basis_name.drop_duplicates():
        targets = targeted.loc[targeted.basis_name.eq(basis), 'target'].drop_duplicates()
        nrows = max(1, (len(targets)+1)//2)
        fig, axes = plt.subplots(nrows, 2, figsize=(11, 3.2*nrows), squeeze=False)
        for ax, target in zip(axes.flat, targets):
            draw_panel(ax, basis, target)
        for ax in list(axes.flat)[len(targets):]:
            ax.set_axis_off()
        fig.tight_layout()
        save_figure(fig, root, filename+'_'+re.sub(r'[^A-Za-z0-9_.-]+', '_', str(basis)))

    bases = list(targeted.basis_name.drop_duplicates())
    detailed = [b for b in bases if b in ('GE_PRESS_TE35', 'SH_PRESS_TE97')]
    overview_bases = detailed or bases[:2]
    targets = [t for t in ['NAA', 'Cr', 'tCho', 'bHG'] if t in targeted.target.values]
    if not targets:
        targets = list(targeted.target.drop_duplicates())[:4]
    fig, axes = plt.subplots(max(1, len(targets)), max(1, len(overview_bases)),
        figsize=(5.5*max(1, len(overview_bases)), 3.2*max(1, len(targets))), squeeze=False)
    if not panels:
        empty_axes(axes[0, 0])
    for row, target in enumerate(targets):
        for col, basis in enumerate(overview_bases):
            draw_panel(axes[row, col], basis, target)
    fig.tight_layout()
    save_figure(fig, root, filename)


def bias_figure(data, root):
    rows = data.loc[data.stage.isin(['clean', 'transfer']) & data.valid & data.scale.eq(1)].copy()
    rows['relative_error'] = 100*(rows.estimate_mM-rows.truth_mM)/rows.truth_mM.replace(0, np.nan)
    rows['basis_tissue'] = rows.basis_name.astype(str) + '/' + rows.tissue.astype(str)
    fitters = sorted(rows.fitter.unique())
    grids = [rows.loc[rows.fitter.eq(fitter)].pivot_table(index='basis_tissue', columns='metabolite', values='relative_error') for fitter in fitters]
    heights = [max(3, 1.8+.28*len(grid)) for grid in grids] or [3]
    fig, axes = plt.subplots(len(heights), 1, figsize=(10, sum(heights)),
        gridspec_kw={'height_ratios': heights}, squeeze=False)
    if not fitters:
        empty_axes(axes[0, 0])
    for ax, fitter, grid in zip(axes[:, 0], fitters, grids):
        if grid.empty:
            empty_axes(ax)
            continue
        bound = max(1., float(np.nanmax(np.abs(grid.to_numpy()))))
        plotted = ax.imshow(grid, aspect='auto', cmap='RdBu_r', vmin=-bound, vmax=bound)
        ax.set_xticks(range(len(grid.columns)), grid.columns, rotation=45, ha='right')
        ax.set_yticks(range(len(grid.index)), grid.index, fontsize=8)
        ax.set_title(f'{fitter}: baseline concentration bias (scale=1; clean/transfer pooled)')
        fig.colorbar(plotted, ax=ax, label='Mean relative error (%)')
    fig.tight_layout()
    save_figure(fig, root, 'tissue_bias')


def stress_figure(data, root, stage, variable, label, filename):
    rows = data.loc[data.stage.eq(stage) & np.isfinite(data[variable])]
    keys = ['basis_name', 'fitter', 'metabolite', variable]
    pooled = summarize(rows, keys)
    pooled.to_csv(root / f'{filename}.csv', index=False)
    fitters = sorted(rows.fitter.unique())
    bases = list(pooled.basis_name.drop_duplicates())
    detailed = [b for b in bases if b in ('GE_PRESS_TE35', 'SH_PRESS_TE97')]
    overview = pooled.loc[pooled.basis_name.isin(detailed or bases[:2]) &
        pooled.metabolite.isin(['NAA', 'bHG', 'tCho', 'tCr'])]
    fig, axes = plt.subplots(max(1, len(fitters)), 1,
        figsize=(9, max(3.5, 3.5*len(fitters))), squeeze=False)
    if not fitters:
        empty_axes(axes[0, 0])
    for ax, fitter in zip(axes[:, 0], fitters):
        for (basis, metabolite), group in overview.loc[overview.fitter.eq(fitter)].groupby(['basis_name', 'metabolite']):
            group = group.sort_values(variable)
            ax.plot(group[variable], group.rmse_mM, 'o-', markersize=3, label=f'{basis}/{metabolite}')
        ax.set(title=f'{fitter}: {stage} overview', xlabel=label, ylabel='RMSE (mM)')
        if ax.lines:
            ax.legend(fontsize=7, ncol=2)
        else:
            empty_axes(ax)
    fig.tight_layout()
    save_figure(fig, root, filename)

    for basis in bases:
        subset = pooled.loc[pooled.basis_name.eq(basis)]
        metabolites = [m for m in ['NAA', 'Cr', 'bHG', 'Glu', 'GABA', 'GSH', 'tCho', 'tCr', 'tNAA']
            if m in subset.metabolite.values]
        if not metabolites:
            continue
        nrows = (len(metabolites)+1)//2
        fig, axes = plt.subplots(nrows, 2, figsize=(11, 3*nrows), squeeze=False)
        for ax, metabolite in zip(axes.flat, metabolites):
            for fitter, group in subset.loc[subset.metabolite.eq(metabolite)].groupby('fitter'):
                group = group.sort_values(variable)
                ax.plot(group[variable], group.rmse_mM, 'o-', markersize=3, label=fitter)
            ax.set(title=f'{basis}: {metabolite}', xlabel=label, ylabel='RMSE (mM)')
            ax.legend(fontsize=8)
        for ax in list(axes.flat)[len(metabolites):]:
            ax.set_axis_off()
        fig.suptitle(f'{stage}: tissues, target changes and repeats pooled', fontsize=11)
        fig.tight_layout()
        save_figure(fig, root, filename+'_'+re.sub(r'[^A-Za-z0-9_.-]+', '_', str(basis)))


def crosstalk_figure(slopes, root):
    rows = slopes.loc[slopes.metabolite.ne(slopes.target)].copy()
    rows['basis_tissue_fitter'] = rows.basis_name.astype(str)+'/'+rows.tissue.astype(str)+'/'+rows.fitter.astype(str)+'/'+rows.stage.astype(str)
    targets = sorted(rows.target.dropna().unique())

    def draw_grid(ax, grid, title):
        if grid.empty:
            empty_axes(ax)
            return
        bound = max(.01, float(np.nanmax(np.abs(grid.to_numpy()))))
        plotted = ax.imshow(grid, aspect='auto', cmap='RdBu_r', vmin=-bound, vmax=bound)
        ax.set_xticks(range(len(grid.columns)), grid.columns, rotation=45, ha='right')
        ax.set_yticks(range(len(grid.index)), grid.index, fontsize=8)
        ax.set_title(title)
        ax.figure.colorbar(plotted, ax=ax, label='Fitted non-target change / target truth change')

    for target in targets:
        grid = rows.loc[rows.target.eq(target)].pivot_table(index='basis_tissue_fitter', columns='metabolite', values='slope')
        fig, ax = plt.subplots(figsize=(11, max(3.5, 2+.28*len(grid))))
        draw_grid(ax, grid, f'Cross talk when varying {target}: all bases/tissues/fitters')
        fig.tight_layout()
        save_figure(fig, root, 'cross_talk_'+re.sub(r'[^A-Za-z0-9_.-]+', '_', str(target)))

    bases = list(rows.basis_name.drop_duplicates())
    detailed = [b for b in bases if b in ('GE_PRESS_TE35', 'SH_PRESS_TE97')]
    overview = rows.loc[rows.basis_name.isin(detailed or bases[:2])].copy()
    overview['basis_fitter'] = overview.basis_name.astype(str)+'/'+overview.fitter.astype(str)+'/'+overview.stage.astype(str)
    primary_targets = [t for t in ['NAA', 'bHG'] if t in targets] or targets[:2]
    grids = [overview.loc[overview.target.eq(t)].pivot_table(index='basis_fitter', columns='metabolite', values='slope')
        for t in primary_targets]
    heights = [max(3.5, 2+.28*len(grid)) for grid in grids] or [3.5]
    fig, axes = plt.subplots(len(heights), 1, figsize=(11, sum(heights)),
        gridspec_kw={'height_ratios': heights}, squeeze=False)
    if not primary_targets:
        empty_axes(axes[0, 0])
    for ax, target, grid in zip(axes[:, 0], primary_targets, grids):
        draw_grid(ax, grid, f'{target} cross talk: mean per-tissue regression slope')
    fig.tight_layout()
    save_figure(fig, root, 'cross_talk')


def robust_responses(data, clean_slopes):
    """Fit each stressed sweep against its exact same-condition baseline."""
    keys = ['basis_name', 'tissue', 'fitter', 'stage', 'snr', 'b0_sd', 'repeat', 'target']
    stressed = data.loc[data.stage.isin(['noise', 'shim', 'combined'])]
    perturbed = stressed.loc[stressed.target.isin(['NAA', 'bHG']) & stressed.metabolite.eq(stressed.target)]
    baseline = stressed.loc[stressed.target.eq('baseline') & stressed.metabolite.isin(['NAA', 'bHG'])].copy()
    baseline['target'] = baseline.metabolite
    anchors = baseline[keys+['case_id', 'truth_mM', 'estimate_mM', 'valid']].rename(columns={
        'case_id': 'baseline_case_id', 'truth_mM': 'baseline_truth_mM',
        'estimate_mM': 'baseline_estimate_mM', 'valid': 'baseline_valid'})
    # Only baselines with a perturbed sweep enter the response calculation.
    baseline = baseline.merge(perturbed[keys].drop_duplicates(), on=keys, how='inner')
    paired = pd.concat([perturbed, baseline], ignore_index=True).merge(anchors,
        on=keys, how='left', validate='many_to_one')
    records = []
    for values, rows in paired.groupby(keys, dropna=False):
        good = rows.loc[rows.valid & rows.baseline_valid.eq(True) &
            np.isfinite(rows.truth_mM) & np.isfinite(rows.estimate_mM)]
        record = dict(zip(keys, values), attempted=len(rows), successful=len(good),
            truth_levels=good.truth_mM.nunique(), status='insufficient_levels',
            slope=np.nan, intercept_mM=np.nan, r_squared=np.nan,
            baseline_bias_mM=np.nan, delta_rmse_mM=np.nan, fold_change_rmse=np.nan)
        if good.truth_mM.nunique() >= 3:
            fitted = linregress(good.truth_mM, good.estimate_mM)
            changes = good.loc[good.case_id.ne(good.baseline_case_id)]
            delta_error = (changes.estimate_mM-changes.baseline_estimate_mM)-(
                changes.truth_mM-changes.baseline_truth_mM)
            fold_error = changes.estimate_mM/changes.baseline_estimate_mM.replace(0, np.nan)-(
                changes.truth_mM/changes.baseline_truth_mM.replace(0, np.nan))
            record.update(status='ok', slope=fitted.slope, intercept_mM=fitted.intercept,
                r_squared=fitted.rvalue**2 if good.estimate_mM.nunique()>1 else np.nan,
                baseline_bias_mM=(good.baseline_estimate_mM-good.baseline_truth_mM).iloc[0],
                delta_rmse_mM=np.sqrt((delta_error**2).mean()),
                fold_change_rmse=np.sqrt((fold_error**2).mean()))
        records.append(record)
    columns = keys+['attempted', 'successful', 'truth_levels', 'status', 'slope',
        'intercept_mM', 'r_squared', 'baseline_bias_mM', 'delta_rmse_mM', 'fold_change_rmse']
    fits = pd.DataFrame(records, columns=columns)
    reference_keys = ['basis_name', 'tissue', 'fitter', 'target']
    reference = clean_slopes.loc[clean_slopes.metabolite.eq(clean_slopes.target),
        reference_keys+['slope', 'intercept_mM']].rename(columns={
            'slope': 'clean_slope', 'intercept_mM': 'clean_intercept_mM'})
    fits = fits.merge(reference, on=reference_keys, how='left', validate='many_to_one')
    summary_keys = [key for key in keys if key != 'repeat']
    stats = fits.groupby(summary_keys, dropna=False).agg(
        attempted_repeats=('status', 'size'), successful_repeats=('slope', 'count'),
        slope_mean=('slope', 'mean'), slope_sd=('slope', 'std'),
        intercept_mean_mM=('intercept_mM', 'mean'), intercept_sd_mM=('intercept_mM', 'std'),
        delta_rmse_mean_mM=('delta_rmse_mM', 'mean'), fold_change_rmse_mean=('fold_change_rmse', 'mean'),
        clean_slope=('clean_slope', 'first'), clean_intercept_mM=('clean_intercept_mM', 'first')).reset_index()
    stats['failed_repeats'] = stats.attempted_repeats-stats.successful_repeats
    return fits, stats


def robust_figure(stats, root):
    stats = stats.assign(snr=pd.to_numeric(stats.snr, errors='coerce'))
    rows = stats.loc[stats.stage.eq('noise') & stats.tissue.eq('40/40/20') & np.isfinite(stats.snr)]
    bases = list(rows.basis_name.drop_duplicates())
    targets = [target for target in ['NAA', 'bHG'] if target in rows.target.values]
    fig, axes = plt.subplots(max(1, len(targets)), max(1, len(bases)),
        figsize=(5.5*max(1, len(bases)), 3.6*max(1, len(targets))), squeeze=False)
    if not targets:
        empty_axes(axes[0, 0], 'No complete mixed-tissue stressed sweeps available')
    for row, target in enumerate(targets):
        for col, basis in enumerate(bases):
            ax = axes[row, col]
            for fitter, group in rows.loc[rows.target.eq(target) & rows.basis_name.eq(basis)].groupby('fitter'):
                group = group.sort_values('snr')
                ax.errorbar(group.snr, group.slope_mean, yerr=group.slope_sd.fillna(0),
                    fmt='o-', capsize=3, label=fitter)
            ax.axhline(1, color='.5', linestyle='--', label='Ideal slope')
            ax.set(title=f'{basis}: {target}; WM/GM/CSF 40/40/20',
                xlabel='Metabolite SNR', ylabel='Absolute response slope (mean ± SD)')
            ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, root, 'robust_sensitivity')


def markdown_table(frame):
    text = frame.copy()
    for column in text.select_dtypes(include='number'):
        text[column] = text[column].map(lambda v: f'{v:.4g}' if np.isfinite(v) else '—')
    return '\n'.join(['| '+' | '.join(map(str, text.columns))+' |',
        '| '+' | '.join(['---']*len(text.columns))+' |'] +
        ['| '+' | '.join(str(v).replace('|', '/') for v in row)+' |' for row in text.itertuples(index=False, name=None)])


def reference_noise_effect(data):
    """Compare water noise with the same metabolite-noise realization and truth."""
    keys=['basis_name','fitter','metabolite','repeat']
    controls=data.loc[data.stage.eq('noise') & data.target.eq('baseline') & data.tissue.eq('40/40/20') & data.snr.eq(30),keys+['estimate_mM','valid']]
    noisy=data.loc[data.stage.eq('water_noise')].merge(controls.rename(columns={'estimate_mM':'control_mM','valid':'control_valid'}),on=keys,how='left',validate='many_to_one')
    noisy=noisy.loc[noisy.valid & noisy.control_valid.fillna(False)].copy()
    noisy['delta_mM']=noisy.estimate_mM-noisy.control_mM
    return noisy.groupby(['basis_name','fitter','metabolite','water_snr']).delta_mM.agg(
        n='size',mean_delta_mM='mean',sd_delta_mM='std',rmse_delta_mM=lambda x:np.sqrt(np.mean(x*x))).reset_index()


def report(root: Path):
    root = Path(root)
    from .audit import audit
    audit(root)
    cases = pd.read_csv(root / 'cases.csv')
    try:
        estimates = pd.read_csv(root / 'estimates.csv')
    except pd.errors.EmptyDataError:
        estimates = pd.DataFrame(columns=['case_id', 'fitter', 'metabolite', 'truth_mM', 'estimate_mM', 'status'])
    for column in ['scale', 'snr', 'b0_sd', 'water_snr']:
        cases[column] = pd.to_numeric(cases[column], errors='coerce')
    # Metadata belongs to cases even when a fitter emitted its own copied columns.
    data = estimates.drop(columns=[c for c in cases if c != 'case_id' and c in estimates]).merge(
        cases, on='case_id', how='left', validate='many_to_one', indicator=True)
    if data['_merge'].ne('both').any():
        raise ValueError('estimates contain case IDs absent from cases.csv')
    data = data.drop(columns='_merge')
    for column in ['truth_mM', 'estimate_mM']:
        data[column] = pd.to_numeric(data[column], errors='coerce')
    data['valid'] = data.status.eq('ok') & np.isfinite(data.truth_mM) & np.isfinite(data.estimate_mM)
    fitters = sorted(data.fitter.dropna().unique()) or ['unreported']
    expected = cases[['case_id']].merge(pd.DataFrame({'fitter': fitters}), how='cross')
    missing = expected.merge(data[['case_id', 'fitter']].drop_duplicates(), how='left', indicator=True)
    missing = missing.loc[missing['_merge'].eq('left_only'), ['case_id', 'fitter']]
    missing.to_csv(root / 'missing_estimates.csv', index=False)
    data.loc[~data.valid].to_csv(root / 'failed_estimates.csv', index=False)
    data.to_csv(root / 'estimates_with_metadata.csv', index=False)
    conditions = [c for c in ['wm', 'gm', 'csf', 'water_factor', 'mismatched', 'anatomical', 'position', 'relax', 'realistic'] if c in data]
    summary = summarize(data, CASE_COLUMNS + conditions + ['fitter', 'metabolite'])
    summary.to_csv(root / 'summary.csv', index=False)
    clean, slopes = responses(data)
    robust_fits, robust_stats = robust_responses(data, slopes)
    robust_fits.to_csv(root / 'robust_response_slopes.csv', index=False)
    robust_stats.to_csv(root / 'robust_response_summary.csv', index=False)
    clean.to_csv(root / 'response_changes.csv', index=False)
    slopes.loc[slopes.target.eq(slopes.metabolite)].to_csv(root / 'response_slopes.csv', index=False)
    slopes.loc[slopes.target.ne(slopes.metabolite)].to_csv(root / 'cross_talk.csv', index=False)
    with plt.rc_context({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42}):
        response_figure(clean, root)
        response_figure(clean, root, 'residual')
        response_figure(clean, root, 'delta')
        response_figure(clean, root, 'fold')
        bias_figure(data, root)
        stress_figure(data, root, 'noise', 'snr', 'Metabolite SNR', 'noise_rmse')
        stress_figure(data, root, 'shim', 'b0_sd', 'B0 standard deviation (Hz)', 'shim_rmse')
        stress_figure(data, root, 'water_noise', 'water_snr', 'Water SNR', 'water_noise_rmse')
        crosstalk_figure(slopes, root)
        robust_figure(robust_stats, root)
    reference_effect=reference_noise_effect(data)
    reference_effect.to_csv(root / 'reference_noise_effect.csv',index=False)
    overview = summarize(data, ['stage', 'basis_name', 'fitter'])
    results = dict(cases=len(cases), estimates=len(data), successful=int(data.valid.sum()),
        failed=int((~data.valid).sum()), missing_case_fitter_pairs=len(missing),
        stage_summary=json.loads(overview.to_json(orient='records')),
        robust_response_summary=json.loads(robust_stats.to_json(orient='records')),
        response_slopes=json.loads(slopes.loc[slopes.target.eq(slopes.metabolite)].to_json(orient='records')))
    (root / 'summary.json').write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')
    (root / 'report.md').write_text(f'''# Headless BASIS-water fitting results

{len(cases)} cases; {len(data)} reported estimates; {int(data.valid.sum())} successful finite estimates;
{int((~data.valid).sum())} failed or nonfinite estimates; {len(missing)} missing case/fitter pairs.

## Error summaries

Bias is fitted minus ground-truth concentration. RMSE and mean/SD use only successful finite
estimates; failures remain in attempted counts. SD is the sample SD (undefined for one result).
Relative error is undefined for zero truth. This overview pools metabolites, tissues, sweep
levels, and repeats within each stage/basis/fitter; use [summary.csv](summary.csv) for stratified results.

{markdown_table(overview)}

## Figures and response estimates

- [Target response](target_response.pdf) ([PNG](target_response.png)): clean/transfer concentration sweeps, compared with identity. These combined figures show key targets for the detailed bases; `target_response_<basis>`, `delta_response_<basis>` and `fold_response_<basis>` PDF/PNG files contain every target for each basis.
- [Signed concentration error](target_residual.pdf) ([PNG](target_residual.png)): the same panels and curves as target response, with fitted minus ground truth on the y-axis and a zero-error reference. `target_residual_<basis>` PDF/PNG files pair with the per-basis response plots.
- [Response slopes](response_slopes.csv): estimate versus target ground truth, with intercept and R²; at least two distinct concentrations required. Baseline rows are included in each target response group.
- [Baseline-relative changes](response_changes.csv): truth/fitted deltas, fold changes, and errors. The shared baseline is matched by basis/tissue/fitter/metabolite; zero baseline concentrations make fold changes undefined. Delta and fold-change RMSE are also included with response slopes.
- [Stressed sensitivity](robust_sensitivity.pdf) ([PNG](robust_sensitivity.png)): mixed 40/40/20 NAA/bHG slopes versus SNR, with per-fitter mean ± sample SD and ideal slope=1. [Per-repeat fits](robust_response_slopes.csv) and [replicate summaries](robust_response_summary.csv) also cover shim and combined stages. Baselines match the exact basis/tissue/fitter/stage/SNR/B0/repeat/target; at least three successful distinct truth levels are required. Paired delta/fold RMSE excludes the baseline anchor itself. Single-repeat shim SD is undefined. Clean slope/intercept are included as context and never used to correct fitted concentrations.
- [Cross talk](cross_talk.pdf) ([CSV](cross_talk.csv)): overview for NAA/bHG in the detailed bases, using arithmetic means of per-tissue regression slopes. `cross_talk_<target>` PDF/PNG files retain all basis/tissue/fitter rows with height scaled to row count. The CSV is unpooled; slopes measure non-target estimate change per nominal-mM target concentration change.
- [Delta response](delta_response.pdf) and [fold response](fold_response.pdf): target fitted changes relative to the shared baseline, with baseline points and identity lines.
- [Tissue bias](tissue_bias.pdf): mean relative error at scale=1 in clean/transfer cases; zero-truth entries excluded.
- [Noise RMSE](noise_rmse.pdf), [shim RMSE](shim_rmse.pdf), [water-noise RMSE](water_noise_rmse.pdf): tissues, target changes and repeats pooled; finite stress levels only. Main figures show NAA/bHG/tCho/tCr for the two detailed bases; per-basis PDF/PNG files facet NAA/Cr/bHG/Glu/GABA/GSH/tCho/tCr/tNAA with separate fitter curves. CSVs retain all metabolites.
- [Failed estimates](failed_estimates.csv) and [missing case/fitter pairs](missing_estimates.csv) are explicit. An absent metabolite within a reported fit cannot be inferred from these input tables.

## Reference noise with paired metabolite controls

[Paired reference-noise effects](reference_noise_effect.csv) compare against the same metabolite-noise realization with clean water. Changes here reflect the water-reference pipeline rather than a different metabolite noise draw. Estimates retain their uncorrected scale.

{markdown_table(reference_effect.loc[reference_effect.metabolite.eq('NAA')])}

## Interpretation and limits

Concentrations and water referencing retain the simulator's nominal-mM and approximate BASIS
response assumptions. Ground truth is separate from fitted estimates. tCho, when supplied,
is the GPC+PCh truth and corresponding grouped fit estimate. No additional proton factor,
calibration, or response correction is inferred by this report.

Noise experiments with ten repetitions characterize only this small numerical sample;
SD and RMSE are descriptive, and no pass/fail threshold or population precision is asserted.
Failures are excluded from numerical accuracy metrics and must be assessed alongside them.
Missing case/fitter pairs are listed separately and are not included in reported-row failure rates.
R² is undefined for constant fitted responses; it does not establish unbiased recovery.
These results do not establish scanner calibration, clinical accuracy, or general absolute quantification.
''')
    return results


def self_check():
    rows = pd.DataFrame({'group': ['a']*3, 'truth_mM': [1., 3., 9.],
        'estimate_mM': [2., 2., np.nan], 'valid': [True, True, False]})
    result = summarize(rows, ['group']).iloc[0]
    assert result.attempted == 3 and result.failed == 1
    assert result.bias_mM == 0 and result.rmse_mM == 1
    assert result.mean_mM == 2 and result.sd_mM == 0
    assert summarize(rows.iloc[:0], ['group']).empty
    sweep = pd.DataFrame({'case_id': [0, 1, 2, 0, 1, 2], 'stage': ['clean']*6,
        'basis_name': ['test']*6, 'tissue': ['WM']*6, 'target': ['NAA']*6,
        'fitter': ['fit']*6, 'metabolite': ['NAA']*3+['Cr']*3, 'snr': [np.inf]*6,
        'b0_sd': [0]*6, 'truth_mM': [1., 2., 3., 8., 8., 8.],
        'estimate_mM': [3., 5., 7., 8.1, 8.2, 8.3], 'valid': [True]*6})
    _, fits = responses(sweep)
    targeted = fits.loc[fits.metabolite.eq('NAA')].iloc[0]
    other = fits.loc[fits.metabolite.eq('Cr')].iloc[0]
    assert np.isclose(targeted.slope, 2) and np.isclose(targeted.intercept_mM, 1)
    assert np.isclose(targeted.r_squared, 1) and np.isclose(other.slope, .1)
    anchored = sweep.copy()
    anchored.loc[anchored.case_id.eq(1), 'target'] = 'baseline'
    changes, fits = responses(anchored)
    assert len(changes) == 6
    targeted = fits.loc[fits.metabolite.eq('NAA')].iloc[0]
    assert np.isclose(targeted.slope, 2) and np.isclose(targeted.delta_rmse_mM, np.sqrt(2/3))
    assert changes.baseline_estimate_mM.notna().all()
    stressed = []
    for snr in [20, 100]:
        for repeat in [0, 1]:
            for index, truth in enumerate([5., 10., 15.]):
                noise = [.2, -.1, .3][index]*(1 if repeat == 0 else -1)
                stressed.append(dict(case_id=f'{snr}-{repeat}-{index}', basis_name='test',
                    tissue='40/40/20', fitter='fit', stage='noise', snr=snr, b0_sd=0,
                    repeat=repeat, target='baseline' if index == 1 else 'NAA',
                    metabolite='NAA', truth_mM=truth,
                    estimate_mM=1.5*truth+snr+repeat+noise, valid=True))
    stressed = pd.DataFrame(stressed)
    fits = fits.assign(tissue='40/40/20')
    robust, stats = robust_responses(stressed, fits)
    np.testing.assert_allclose(robust.loc[robust.repeat.eq(0), 'slope'], 1.51)
    np.testing.assert_allclose(robust.loc[robust.repeat.eq(1), 'slope'], 1.49)
    np.testing.assert_allclose(stats.slope_mean, 1.5)
    np.testing.assert_allclose(stats.slope_sd, np.sqrt(2)*.01)
    assert robust.status.eq('ok').all()
    assert robust.clean_slope.eq(2).all()
    for _, row in robust.iterrows():
        paired = stressed.loc[stressed.snr.eq(row.snr) & stressed.repeat.eq(row['repeat'])]
        base = paired.loc[paired.target.eq('baseline')].iloc[0]
        points = paired.loc[paired.target.ne('baseline')]
        error = (points.estimate_mM-base.estimate_mM)-(points.truth_mM-base.truth_mM)
        assert np.isclose(row.delta_rmse_mM, np.sqrt((error**2).mean()))
    incomplete = stressed.loc[~stressed.case_id.eq('20-0-1')]
    failed, _ = robust_responses(incomplete, fits)
    assert failed.status.eq('ok').sum() == 3
    empty, summary = robust_responses(stressed.iloc[:0], fits)
    assert empty.empty and summary.empty



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path, nargs='?')
    parser.add_argument('--self-check', action='store_true')
    args = parser.parse_args()
    if args.self_check:
        self_check()
    if args.root:
        report(args.root)
    elif not args.self_check:
        parser.error('root is required unless --self-check is used')

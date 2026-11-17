"""Recompute the five published baseline chromatin and methylation panels."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def calculate(inputs, output, projected=None, windows_dir=None, baseline=None):
    inputs, output = Path(inputs), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    def read(name):
        return pd.read_csv(inputs / name, sep='\t', float_precision='round_trip')
    def save(data, name):
        data.to_csv(output / name, sep='\t', index=False)
    projected = Path(projected) if projected else inputs
    windows_dir = Path(windows_dir) if windows_dir else inputs
    lads = pd.read_csv(projected / 'LAD_values.tsv', sep='\t', float_precision='round_trip')
    names = ['H3K9me2_H3mean', 'H3K9me2_H3rep1', 'H3K9me2_H3rep2',
             'DNMT3A_rep1', 'DNMT3A_rep2', 'DNMT3A_rep3']
    supported = np.isfinite(lads[names]).all(axis=1)
    values = lads.loc[supported].copy()
    ordered = values.sort_values('B1_TTX_mean', kind='stable')
    groups = []
    for i, ids in enumerate(np.array_split(np.arange(len(ordered)), 20), 1):
        group = ordered.iloc[ids]
        groups.append({'group': i, 'LADs': len(group), 'B1': group.B1_TTX_mean.mean(),
                       'DNMT3A_mean': group.DNMT3A_mean.mean()})
    groups = pd.DataFrame(groups)
    save(groups, 'ranked_LAD_groups.tsv')
    rho = spearmanr(values.B1_TTX_mean, values.DNMT3A_mean).correlation
    save(pd.DataFrame([{'measurement': 'DNMT3A_mean', 'LADs': len(values), 'Spearman': rho}]), 'B1_correlations.tsv')

    profiles = np.load(projected / 'body_profile.npz')
    np.testing.assert_array_equal(profiles['LAD_ids'], lads.LAD_id)
    paired = np.isfinite(lads[names + [name + '_outside' for name in names]]).all(axis=1).to_numpy()
    body = pd.DataFrame({'x': np.r_[np.linspace(-1, 0, 18)[:-1]+1/34, np.linspace(0, 1, 51)[:-1]+.01, np.linspace(1, 2, 18)[:-1]+1/34], 'LADs': np.isfinite(profiles['body'][paired, :, 0]).sum(axis=0),
                         'H3K9me2_H3mean': np.nanmean(profiles['body'][paired, :, 0], axis=0)})
    save(body, 'paired_body_profiles.tsv')

    windows = pd.read_csv(windows_dir / 'mCG_DNMT3A_windows.tsv.gz', sep='\t', float_precision='round_trip')
    windows = windows.loc[windows.eligible].copy()
    windows['DNMT3A_mean'] = windows[['DNMT3A_rep1', 'DNMT3A_rep2', 'DNMT3A_rep3']].mean(axis=1)
    ranked = windows.sort_values('DNMT3A_mean', kind='stable')
    rows = []
    for i, start in enumerate(range(0, len(ranked), 5000), 1):
        group = ranked.iloc[start:start + 5000]
        rows.append({'group': i, 'windows': len(group), 'x': group.DNMT3A_mean.mean(),
                     **{name: group[name].mean() for name in ['Baseline', 'A0', 'B0', 'C0', '1wk', '2wk']}})
    rows = pd.DataFrame(rows)
    save(rows, 'mean_DNMT3A_rank_groups.tsv')

    bins = read('B1_mCG_bins.tsv.gz')
    unrounded = np.load((Path(baseline) if baseline else inputs) / 'baseline_counts.npz')
    original_mcg = (100 * unrounded['meth'] / unrounded['depth']).mean(axis=1)
    original_b1 = unrounded['B1']
    trends, correlations = [], []
    for scope, keep in [('All eligible bins', bins.eligible), ('TTX LAD bins', bins.eligible & (bins.TTX_state >= 1))]:
        subset = bins.loc[keep].copy()
        quantile = pd.qcut(subset.TTX_mean, 20, duplicates='drop')
        trend = subset.assign(q=quantile).groupby('q', observed=True).agg(
            x=('TTX_mean', 'median'), median=('mCG_mean', 'median'),
            lo=('mCG_mean', lambda x: x.quantile(.25)), hi=('mCG_mean', lambda x: x.quantile(.75)),
            n=('mCG_mean', 'size')).reset_index(drop=True)
        trend['scope'] = scope
        trends.append(trend)
        correlations.append({'scope': scope, 'n_bins': len(subset),
                             'Spearman_rho': spearmanr(original_b1[keep], original_mcg[keep]).correlation})
    trends = pd.concat(trends, ignore_index=True)
    correlations = pd.DataFrame(correlations)
    save(trends, 'quantile_trends.tsv')
    save(correlations, 'association_summary.tsv')
    return correlations


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--projected', type=Path)
    parser.add_argument('--windows', type=Path)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    print(calculate(args.inputs, args.output, args.projected, args.windows, args.baseline).to_string(index=False))

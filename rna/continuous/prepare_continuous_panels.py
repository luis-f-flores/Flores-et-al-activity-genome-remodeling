"""Calculate regression curves and plotting summaries from continuous RNA results."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy import stats


def fit(x, y):
    f = stats.linregress(x, y)
    xx = np.linspace(min(x), max(x), 151)
    residual = y - (f.intercept + f.slope * x)
    se = np.sqrt(np.sum(residual**2) / (len(x) - 2))
    ci = stats.t.ppf(.975, len(x) - 2) * se * np.sqrt(
        1 / len(x) + (xx - np.mean(x))**2 / np.sum((x - np.mean(x))**2))
    yy = f.intercept + f.slope * xx
    return dict(x=xx.tolist(), y=yy.tolist(), lo=(yy-ci).tolist(), hi=(yy+ci).tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True, help='Output from continuous_rna.py.')
    parser.add_argument('--b1-groups', type=Path, help='pairwise_gene_RNA.tsv from the B1 LAD/RNA analysis.')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    read = lambda name: pd.read_csv(args.inputs/name, sep='\t', float_precision='round_trip')
    panels = {'M03': [], 'S06': [], 'S17': []}
    prom, body = read('promoter_points.tsv'), read('gene_body_points.tsv')
    for tp in (1, 6, 24):
        d = prom[prom.timepoint == tp].sort_values('gene_id')
        x, y = d.delta_mCG_percent.to_numpy(), d.RNA_change.to_numpy()
        result = stats.pearsonr(x, y)
        scaled = stats.linregress((x-x.min())/(x.max()-x.min())-.5, y)
        panels['M03'].append(dict(tp=tp, n=len(d), x=x.tolist(), y=y.tolist(),
                                 r=result.statistic, p=result.pvalue, R2=result.statistic**2,
                                 slope=scaled.slope, se=scaled.stderr, fit=fit(x, y)))
    for cohort in ('differentially_methylated', 'DEG_subset'):
        for tp in (1, 6, 24):
            d = body[(body.timepoint == tp) & (body.cohort == cohort)]
            x, y = d.delta_mCG_percent.to_numpy(), d.RNA_change.to_numpy()
            result = stats.pearsonr(x, y)
            panels['S06'].append(dict(tp=tp, cohort=cohort, n=len(d), x=x.tolist(), y=y.tolist(),
                                     r=result.statistic, p=result.pvalue, fit=fit(x, y)))
    regions, spatial = read('lad_points.tsv'), read('spatial_tests.tsv').set_index('contrast')
    contrasts = ('Bic_vs_TTX', 'Bic_Nimo_vs_TTX', 'Bic_Nimo_vs_Bic')
    for contrast in contrasts:
        d = regions[regions.contrast == contrast]
        panels['S17'].append(dict(contrast=contrast, n=len(d), x=d.B1_delta.tolist(),
                                 y=d.RNA_change.tolist(), rho=float(spatial.loc[contrast, 'rho']),
                                 p=float(spatial.loc[contrast, 'P'])))
    if args.b1_groups:
        means = read('Bic_mean_TPM.tsv').set_index('gene_id')
        groups = pd.read_csv(args.b1_groups, sep='\t', float_precision='round_trip')
        panels['S17boxes'] = []
        for contrast, target in zip(contrasts[:2], ('Bic', 'Bic_Nimo')):
            effects = np.log2((means[target]+1)/(means.TTX+1))
            d = groups[groups.contrast == contrast].copy()
            d['RNA_change'] = effects.loc[d.gene_id].to_numpy()
            values = [d.loc[d.B1_group == g, 'RNA_change'].to_numpy()
                      for g in ('Lower B1', 'Not detected', 'Higher B1')]
            summaries = []
            for group, v in zip(('Lower B1', 'Not detected', 'Higher B1'), values):
                q1, median, q3 = np.quantile(v, [.25, .5, .75]); iqr = q3-q1
                summaries.append(dict(group=group, n=len(v), q25=q1, median=median, q75=q3,
                                      whislo=float(min(v[v >= q1-1.5*iqr])),
                                      whishi=float(max(v[v <= q3+1.5*iqr]))))
            p = stats.mannwhitneyu(values[0], values[2], alternative='two-sided').pvalue
            panels['S17boxes'].append(dict(contrast=contrast, groups=summaries, p=p))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'panel_data.json').write_text(json.dumps(panels, allow_nan=False), encoding='utf-8')


if __name__ == '__main__':
    main()

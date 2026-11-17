"""Conditional DEG–LAD enrichment with expression, gene-length and Input matching."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

CONTRASTS = ['Bic_vs_TTX', 'Bic_Nimo_vs_Bic', 'Bic_Nimo_vs_TTX']


def cdf(values, foreground, weights):
    order = np.argsort(values, kind='stable')
    values = values[order]
    cumulative = np.cumsum((foreground / foreground.sum() - weights / weights.sum())[order])
    last = np.r_[values[1:] != values[:-1], True]
    return float(abs(cumulative[last]).max())


def strata(labels, foreground):
    weights = np.zeros(len(labels))
    cells = []
    for label in np.unique(labels):
        index = np.flatnonzero(labels == label)
        count = int(foreground[index].sum())
        weights[index] = count / len(index)
        cells.append((int(label), index, count))
    return weights, cells


def calculate(data, contrast):
    foreground = data.DEG.to_numpy(bool)
    lad = data.LAD.to_numpy(bool)
    included = ['log2_reference_FPKM_plus1', 'log2_RNA_counting_length']
    base_labels = data.stratum_5.to_numpy(int)
    selected = None
    if contrast == 'Bic_Nimo_vs_TTX':
        weights, cells = strata(base_labels, foreground)
        selected = (None, base_labels, weights, cells)
    else:
        reference = 'TTX' if contrast == 'Bic_vs_TTX' else 'Bic'
        target = 'Bic' if contrast == 'Bic_vs_TTX' else 'Bic_Nimo'
        added = [f'log2_Input_{reference}_CPM_per_clean_kb', f'log2_Input_{target}_CPM_per_clean_kb']
        included += added
        for quantiles in [5, 10]:
            bins = [pd.qcut(data[column], quantiles, labels=False, duplicates='drop').to_numpy(int) for column in added]
            labels = base_labels * quantiles * quantiles + bins[0] * quantiles + bins[1]
            weights, cells = strata(labels, foreground)
            adequate = True
            for column in included:
                values = data[column].to_numpy(float)
                smd = (values[foreground].mean() - np.average(values, weights=weights)) / values.std(ddof=1)
                adequate &= abs(smd) <= .1 and cdf(values, foreground, weights) <= .1
            fixed = sum(count for _, index, count in cells if len(index) == count)
            if adequate and fixed / foreground.sum() <= .1:
                selected = (quantiles, labels, weights, cells)
                break
        if selected is None:
            raise ValueError(f'No Input matching partition met covariate balance and support criteria: {contrast}')
    quantiles, labels, weights, cells = selected
    balance = []
    for column in included:
        values = data[column].to_numpy(float)
        balance.append(dict(contrast=contrast, Input_quantiles=quantiles, covariate=column,
                            foreground_mean=float(values[foreground].mean()),
                            expected_mean=float(np.average(values, weights=weights)),
                            SMD=float((values[foreground].mean() - np.average(values, weights=weights)) / values.std(ddof=1)),
                            CDF=cdf(values, foreground, weights)))
    distribution = np.array([1.])
    offset = 0
    expected = variance = 0.
    rows = []
    for label, index, count in cells:
        total, hits = len(index), int(lad[index].sum())
        low, high = max(0, count - (total - hits)), min(count, hits)
        distribution = np.convolve(distribution, hypergeom.pmf(np.arange(low, high + 1), total, hits, count))
        offset += low
        expected += count * hits / total
        variance += count * hits / total * (1 - hits / total) * (total - count) / (total - 1) if total > 1 else 0
        rows.append(dict(contrast=contrast, stratum=label, N=total, DEG=count, LAD=hits,
                         observed=int((foreground[index] & lad[index]).sum())))
    distribution /= distribution.sum()
    support = np.arange(len(distribution)) + offset
    observed = int((foreground & lad).sum())
    pvalue = min(1., 2 * min(float(distribution[support <= observed].sum()), float(distribution[support >= observed].sum())))
    lower, upper = support[np.searchsorted(np.cumsum(distribution), [.025, .975])]
    result = dict(contrast=contrast, eligible=len(data), DEG=int(foreground.sum()), LAD=int(lad.sum()),
                  observed=observed, Input_quantiles=quantiles, matched_expected=expected,
                  fold=observed / expected, P=pvalue, null_2p5=int(lower), null_97p5=int(upper), variance=variance)
    membership = data.copy()
    membership['coverage_match_stratum'] = labels
    membership['coverage_expected_inclusion_weight'] = weights
    return result, balance, rows, membership, dict(offset=offset, pmf=distribution.tolist(), P=pvalue)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output / 'tables'
    output.mkdir(parents=True, exist_ok=True)
    results, balance, cells = [], [], []
    for contrast in CONTRASTS:
        data = pd.read_csv(args.inputs / (contrast + '_membership.tsv'), sep='\t')
        result, b, c, membership, null = calculate(data, contrast)
        results.append(result)
        balance += b
        cells += c
        membership.to_csv(output / (contrast + '_membership.tsv'), sep='\t', index=False)
        (output / (contrast + '_exact_null.json')).write_text(json.dumps(null, indent=2) + '\n')
    pvalues = np.array([row['P'] for row in results])
    order = np.argsort(pvalues)
    qvalues = np.empty(3)
    qvalues[order] = np.minimum(1, np.minimum.accumulate((pvalues[order] * 3 / np.arange(1, 4))[::-1])[::-1])
    for row, qvalue in zip(results, qvalues):
        row['BH_q'] = float(qvalue)
    (output / 'RESULTS.json').write_text(json.dumps({'comparisons': results}, indent=2) + '\n')
    pd.DataFrame(results).to_csv(output / 'overlap_results.tsv', sep='\t', index=False)
    pd.DataFrame(balance).to_csv(output / 'balance.tsv', sep='\t', index=False)
    pd.DataFrame(cells).to_csv(output / 'strata.tsv', sep='\t', index=False)
    print(pd.DataFrame(results).to_string(index=False))


if __name__ == '__main__':
    main()

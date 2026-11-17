"""KCl time-course versus bicuculline RNA changes on pairwise-testable genes."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def tpm(counts, columns):
    rpk = counts[columns].div(counts.gene_length / 1000, axis=0)
    return rpk.div(rpk.sum(), axis=1) * 1e6


def calculate(counts, testability, output):
    counts, output = Path(counts), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    kcl = pd.read_csv(counts / 'counts_KCl.tsv', sep='\t').set_index('gene_id')
    bic = pd.read_csv(counts / 'counts_Bic.tsv', sep='\t').set_index('gene_id')
    tests = pd.read_csv(testability, sep='\t')
    bic_cols = [f'{condition}_{rep}' for condition in ['Basal', 'Bic', 'Bic_nimo'] for rep in 'ABC']
    bt = tpm(bic, bic_cols)
    bic_change = np.log2((bt[[f'Bic_{r}' for r in 'ABC']].mean(axis=1) + 1) /
                         (bt[[f'Basal_{r}' for r in 'ABC']].mean(axis=1) + 1))
    rows, tables = [], []
    for tp, label in [(1, '1 h'), (6, '6 h'), (24, 'Recovery')]:
        tested = tests.loc[tests.timepoint.eq(label)].set_index('gene_id')
        finite = np.isfinite(tested[['KCl_pvalue', 'KCl_log2FoldChange',
                                     'Bic_pvalue', 'Bic_log2FoldChange']]).all(axis=1)
        ids = tested.index[finite]
        assert ids.is_unique
        cols = [f'{r}{tp}{c}' for c in ['K', 'M'] for r in 'ABC']
        kt = tpm(kcl, cols)
        kchange = np.log2((kt.loc[ids, cols[:3]].mean(axis=1) + 1) /
                         (kt.loc[ids, cols[3:]].mean(axis=1) + 1))
        d = pd.DataFrame({'gene_id': ids, 'timepoint': label,
                          'KCl_TPM_change': kchange.to_numpy(),
                          'Bic_TPM_change': bic_change.loc[ids].to_numpy()})
        rows.append({'timepoint': label, 'n': len(d),
                     'spearman_rho': spearmanr(d.KCl_TPM_change, d.Bic_TPM_change).correlation})
        tables.append(d)
    result = pd.DataFrame(rows)
    pd.concat(tables).to_csv(output / 'gene_expression_changes.tsv.gz', sep='\t', index=False)
    result.to_csv(output / 'correlations.tsv', sep='\t', index=False)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', type=Path, required=True)
    parser.add_argument('--testability', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(calculate(args.counts, args.testability, args.output).to_string(index=False))

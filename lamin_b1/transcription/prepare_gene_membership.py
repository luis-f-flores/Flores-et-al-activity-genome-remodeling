"""Link RNA comparisons to assessed RefGene TSSs and shared-LAD states."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['genes', 'rna', 'signals', 'decoded', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    genes = pd.read_csv(args.genes, sep='\t')
    rna = pd.read_csv(args.rna, sep='\t')
    signals = dict(np.load(args.signals))
    states = np.load(args.decoded)['states']
    if not genes.gene_id.is_unique:
        raise ValueError('Gene coordinates require unique RNA gene IDs.')
    index = {(str(c), int(s)): i for i, (c, s) in enumerate(zip(signals['chromosome'], signals['start']))}
    groups = {'TTX': [0, 1], 'Bic': [2, 3, 4], 'Bic_Nimo': [5, 6, 7]}
    plans = [('Bic_vs_TTX', 'TTX', 'Bic'), ('Bic_Nimo_vs_Bic', 'Bic', 'Bic_Nimo'),
             ('Bic_Nimo_vs_TTX', 'TTX', 'Bic_Nimo')]
    rows = []
    for contrast, reference, target in plans:
        stats = rna[rna.contrast == contrast].copy()
        d = genes.merge(stats.drop(columns='gene_name'), on='gene_id', how='inner', validate='one_to_one')
        d = d[(d.annotation_status == 'unique') & np.isfinite(d.pvalue)].copy()
        d['TSS_bin_index'] = [index.get((str(c), int(t) // 30000 * 30000), -1)
                              for c, t in zip(d.chromosome, d.tss0)]
        d['baseline_LAD'] = [i >= 0 and bool((states[i, groups[reference]] >= 1).any()) for i in d.TSS_bin_index]
        d['target_LAD'] = [i >= 0 and bool((states[i, groups[target]] >= 1).any()) for i in d.TSS_bin_index]
        d['DEG'] = d.classification.isin(['up', 'down'])
        expected = np.isfinite(d.padj) & (d.padj <= .01) & (abs(d.log2FoldChange) >= .37)
        if not (d.DEG == expected).all():
            raise ValueError(f'RNA classification does not match the specified thresholds: {contrast}')
        d['reference_condition'] = reference
        d['target_condition'] = target
        rows.append(d)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()

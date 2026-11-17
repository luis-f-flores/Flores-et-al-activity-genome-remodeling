"""Prepare expression, counting-length and Input-coverage covariates for DEG–LAD overlap."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

CONTRASTS = ['Bic_vs_TTX', 'Bic_Nimo_vs_TTX', 'Bic_Nimo_vs_Bic']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['eligible', 'rna-counts', 'bins', 'normalization', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    counts = pd.read_csv(args.rna_counts, sep='\t').set_index('gene_id')
    bins = dict(np.load(args.bins))
    norm = json.loads(args.normalization.read_text())
    totals = np.array(norm['autosomal_Input_totals'])
    input_counts = bins['counts'][:, :6]
    args.output.mkdir(parents=True, exist_ok=True)
    for contrast in CONTRASTS:
        d = pd.read_csv(args.eligible / (contrast + '_eligible.tsv'), sep='\t')
        if not d.gene_id.is_unique:
            raise ValueError(f'Duplicate gene IDs in {contrast}.')
        length = counts.loc[d.gene_id, 'gene_length'].to_numpy(float)
        if not (np.isfinite(length).all() and (length > 0).all()):
            raise ValueError('RNA counting lengths must be finite and positive.')
        reference = 'Bic_mean_FPKM' if contrast == 'Bic_Nimo_vs_Bic' else 'TTX_mean_FPKM'
        d['RNA_counting_length'] = length
        d['reference_mean_FPKM'] = d[reference]
        d['log2_reference_FPKM_plus1'] = np.log2(d[reference] + 1)
        d['log2_RNA_counting_length'] = np.log2(length)
        labels = [pd.qcut(d[column], 5, labels=False, duplicates='drop').to_numpy(int)
                  for column in ['log2_reference_FPKM_plus1', 'log2_RNA_counting_length']]
        d['stratum_5'] = labels[0] * 5 + labels[1]
        index = d.TSS_bin_index.to_numpy(int)
        if (index < 0).any():
            raise ValueError('The eligible gene set contains an unassessed TSS.')
        clean = bins['clean_bp'][index]
        for name, pool in [('TTX', [0]), ('Bic', [1, 2]), ('Bic_Nimo', [3, 4, 5])]:
            values = input_counts[index][:, pool].sum(1) / totals[pool].sum() * 1e6 / (clean / 1000)
            d['log2_Input_' + name + '_CPM_per_clean_kb'] = np.log2(1 + values)
        d.to_csv(args.output / (contrast + '_membership.tsv'), sep='\t', index=False, float_format='%.17g')


if __name__ == '__main__':
    main()

"""Combine gene-span B1 measurements and RNA contrasts for nimodipine comparisons."""
from pathlib import Path
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--membership', type=Path, required=True)
    parser.add_argument('--replicates', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    genes = pd.read_csv(args.membership, sep='\t', low_memory=False)
    replicates = pd.read_csv(args.replicates, sep='\t', low_memory=False)
    replicates = replicates[replicates.scheme == 'condition_depth'].set_index('gene_id')
    values = genes[genes.contrast == 'Bic_vs_TTX'].set_index('gene_id').copy()
    if not values.index.is_unique or not replicates.index.is_unique:
        raise ValueError('Gene measurements must contain one row per gene.')
    values = values.loc[values.index.intersection(replicates.index)]
    groups = [['L9', 'L12'], ['L11', 'L14', 'L16'], ['L10', 'L13', 'L15']]
    for name, group in zip(['B1_TTX', 'B1_Bic', 'B1_Nimo'], groups):
        values[name] = replicates.loc[values.index, group].mean(axis=1)
    for contrast in ['Bic_vs_TTX', 'Bic_Nimo_vs_TTX', 'Bic_Nimo_vs_Bic']:
        result = genes[genes.contrast == contrast].set_index('gene_id')
        if not result.index.is_unique:
            raise ValueError(f'Duplicate genes in {contrast}.')
        for key, column in [('fc', 'log2FoldChange'), ('q', 'padj'), ('class', 'classification')]:
            values[key + '_' + contrast] = result[column].reindex(values.index)
        values['LAD_' + contrast] = (result.baseline_LAD | result.target_LAD).reindex(values.index)
    values['b'] = values.B1_Bic - values.B1_TTX
    values['bn'] = values.B1_Nimo - values.B1_TTX
    values['n'] = values.B1_Nimo - values.B1_Bic
    values['B1_toward'] = (values.b * values.n < 0) & (values.bn.abs() < values.b.abs())
    values['RNA_toward'] = (values.fc_Bic_vs_TTX * values.fc_Bic_Nimo_vs_Bic < 0) & (
        values.fc_Bic_Nimo_vs_TTX.abs() < values.fc_Bic_vs_TTX.abs())
    values['RNA_initial_DEG'] = (values.q_Bic_vs_TTX <= .01) & (values.fc_Bic_vs_TTX.abs() >= .37)
    values['RNA_direct_DEG'] = (values.q_Bic_Nimo_vs_Bic <= .01) & (values.fc_Bic_Nimo_vs_Bic.abs() >= .37)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    values.to_csv(args.output, sep='\t')
    print(f'Wrote B1 and RNA responses for {len(values):,} jointly assessable genes.')


if __name__ == '__main__':
    main()

"""Match LAD-associated DEGs to non-DEGs by baseline B1 and RNA abundance."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def build(measurements, output):
    d = pd.read_csv(measurements, sep='\t')
    matched, pairs = [], []
    for contrast in ['Bic_vs_TTX', 'Bic_Nimo_vs_TTX']:
        a = d[(d.contrast == contrast) & (d.scheme == 'condition_depth')].sort_values('gene_id').copy()
        features = np.column_stack([a.B1_TTX, np.log2(a.TTX_mean_FPKM + 1)])
        scale = features.std(0, ddof=1)
        z = features / scale
        controls = np.flatnonzero(a.RNA_classification == 'not_DEG')
        for i in np.flatnonzero(a.RNA_classification.isin(['up', 'down'])):
            distance = abs(z[controls] - z[i])
            eligible = controls[(distance <= .2).all(1) & ((features[controls, 1] == 0) == (features[i, 1] == 0))]
            squared_distance = ((z[eligible] - z[i]) ** 2).sum(axis=1)
            selected = eligible[np.lexsort((a.iloc[eligible].gene_id.to_numpy(), squared_distance))[:5]]
            row = a.iloc[i]
            record = dict(contrast=contrast, gene_id=row.gene_id, gene_name=row.gene_name,
                          RNA_direction=row.RNA_classification, RNA_log2FC=row.RNA_log2FC,
                          DEG_B1_delta=row.B1_delta_log2, controls=len(selected),
                          control_B1_delta=a.iloc[selected].B1_delta_log2.mean() if len(selected) else np.nan)
            record['excess_delta'] = record['DEG_B1_delta'] - record['control_B1_delta']
            matched.append(record)
            for j in selected:
                pairs.append(dict(contrast=contrast, gene_id=row.gene_id, control_gene_id=a.iloc[j].gene_id))
    output.mkdir(parents=True, exist_ok=True)
    matched = pd.DataFrame(matched)
    pairs = pd.DataFrame(pairs)
    matched.to_csv(output / 'matched_gene_scatter.tsv', sep='\t', index=False)
    pairs.to_csv(output / 'matched_gene_pairs.tsv', sep='\t', index=False)
    differences = matched[matched.contrast == 'Bic_vs_TTX'].sort_values('gene_id').copy()
    if (differences.controls == 0).any():
        raise ValueError('A DEG has no eligible matched control.')
    selected_pairs = pairs[pairs.contrast == 'Bic_vs_TTX'].copy()
    genes = d[(d.contrast == 'Bic_vs_TTX') & (d.scheme == 'condition_depth')].set_index('gene_id')
    observations = []
    for row in differences.itertuples():
        controls = selected_pairs.loc[selected_pairs.gene_id == row.gene_id, 'control_gene_id'].tolist()
        for gene_id, is_deg, weight in [(row.gene_id, 1, 1.0)] + [(c, 0, 1 / len(controls)) for c in controls]:
            g = genes.loc[gene_id]
            observations.append(dict(RNA_direction=row.RNA_direction, matched_set=row.gene_id,
                                     gene_id=gene_id, gene_name=g.gene_name, is_DEG=is_deg,
                                     B1_delta=g.B1_delta_log2, analysis_weight=weight,
                                     LAD_cluster=int(g.cluster_id), chromosome=g.chromosome))
    differences.to_csv(output / 'paired_differences.tsv', sep='\t', index=False)
    selected_pairs.to_csv(output / 'control_pairs.tsv', sep='\t', index=False)
    pd.DataFrame(observations).to_csv(output / 'fixed_observations.tsv', sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--measurements', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build(args.measurements, args.output)


if __name__ == '__main__':
    main()

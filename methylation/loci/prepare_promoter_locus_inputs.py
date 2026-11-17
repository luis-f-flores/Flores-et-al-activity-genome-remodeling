"""Prepare promoter-example methylation tracks and DESeq2 RNA measurements."""
from pathlib import Path
import argparse
import pandas as pd
from locus_mcg import extract


def prepare(bigwigs, rna, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for time, path in zip([1, 6, 24], rna):
        table = pd.read_csv(path, sep='\t')
        for gene in ['Ccdc25', 'Serpinb8']:
            found = table[table.gene_name == gene]
            if len(found) != 1:
                raise ValueError(f'Expected one RNA record for {gene}.')
            row = found.iloc[0]
            rows.append(dict(gene=gene, time_h=time, gene_id=row.gene_id,
                             log2FC=row.log2FoldChange, padj=row.padj))
    samples, tracks = extract(bigwigs, [('Ccdc25', 'chr14', 65835119, 65839660, 6),
                                        ('Serpinb8', 'chr1', 107588761, 107592291, 6)])
    samples.to_csv(output/'track_sample_values.tsv', sep='\t', index=False)
    tracks.to_csv(output/'track_plot_values.tsv', sep='\t', index=False)
    pd.DataFrame(rows).to_csv(output/'locus_RNA_values.tsv', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bigwigs', type=Path, required=True)
    parser.add_argument('--rna', type=Path, nargs=3, required=True, help='1 h, 6 h and recovery DE tables.')
    parser.add_argument('--output', type=Path, required=True)
    a = parser.parse_args()
    prepare(a.bigwigs, a.rna, a.output)

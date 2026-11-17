"""Link assessed gene TSSs to baseline B1 states and shared LAD segments."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

GROUPS = [[0, 1], [2, 3, 4], [5, 6, 7]]
CONTRASTS = [('Bic_vs_TTX', 0, 1), ('Bic_Nimo_vs_Bic', 1, 2),
             ('Bic_Nimo_vs_TTX', 0, 2)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['signals', 'decoded', 'genes', 'rna', 'segments', 'output']:
        p.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    signals = np.load(a.signals)
    states = np.load(a.decoded)['states']
    seg = np.load(a.segments / 'all_condition_regions.npz')
    meta = pd.read_csv(a.segments / 'all_condition_regions.tsv', sep='\t')
    genes = pd.read_csv(a.genes, sep='\t', low_memory=False)
    genes = genes[genes.annotation_status == 'unique'].copy()
    index = {(str(c), int(s)): i for i, (c, s) in
             enumerate(zip(signals['chromosome'], signals['start']))}
    genes['bin_index'] = [index.get((str(c), int(t) // 30000 * 30000), -1)
                          for c, t in zip(genes.chromosome, genes.tss0)]
    genes = genes[genes.bin_index >= 0].copy()
    genes['state'] = states[genes.bin_index.to_numpy(int)][:, GROUPS[0]].max(1)
    genes[['gene_id', 'state']].to_csv(a.output / 'baseline_genes_states.tsv', sep='\t', index=False)
    rna = pd.read_csv(a.rna, sep='\t')
    rna = rna[np.isfinite(rna.pvalue) & np.isfinite(rna.log2FoldChange)]
    rna = rna.merge(genes[['gene_id', 'bin_index']], on='gene_id', validate='many_to_one')
    common = set.intersection(*(set(rna.loc[rna.contrast == c, 'gene_id']) for c, _, _ in CONTRASTS))
    rna = rna[rna.gene_id.isin(common)].copy()
    rna['region_index'] = seg['owner'][rna.bin_index.to_numpy(int)]
    rna = rna[rna.region_index >= 0]
    links, regions = [], []
    for contrast, baseline, target in CONTRASTS:
        d = rna[rna.contrast == contrast]
        links.append(d[['contrast', 'gene_id', 'region_index']])
        delta = seg['signal'][:, GROUPS[target]].mean(1) - seg['signal'][:, GROUPS[baseline]].mean(1)
        reg = meta[meta.region_index.isin(d.region_index)].copy()
        reg['B1_delta'] = delta[reg.region_index.to_numpy(int)]
        reg['contrast'] = contrast
        regions.append(reg[['contrast', 'region_index', 'chromosome', 'start0', 'end0', 'B1_delta']].sort_values(['chromosome', 'start0']))
    pd.concat(links).to_csv(a.output / 'S17_links.tsv', sep='\t', index=False)
    pd.concat(regions).to_csv(a.output / 'S17_regions.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()

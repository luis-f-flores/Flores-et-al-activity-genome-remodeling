"""Prepare pairwise shared-LAD segments, standardized B1 values and gene links."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

CONTRASTS = ['Bic_vs_TTX', 'Bic_Nimo_vs_Bic', 'Bic_Nimo_vs_TTX']
PAIRS = [(0, 1), (1, 2), (0, 2)]
GROUPS = [[0, 1], [2, 3, 4], [5, 6, 7]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['signals', 'decoded', 'genes', 'rna', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    tables = args.output / 'tables'
    tables.mkdir(parents=True, exist_ok=True)
    a = dict(np.load(args.signals))
    states = np.load(args.decoded)['states']
    n = len(a['start'])
    if states.shape != a['condition_depth'].shape:
        raise ValueError('B1 values and decoded states must cover the same bins and libraries.')
    z = (a['condition_depth'] - a['condition_depth'].mean(0)) / a['condition_depth'].std(0, ddof=1)
    support = np.column_stack([(states[:, group] >= 1).sum(1) for group in GROUPS])
    lookup = {(str(c), int(s)): i for i, (c, s) in enumerate(zip(a['chromosome'], a['start']))}
    genes = pd.read_csv(args.genes, sep='\t', low_memory=False)
    genes = genes[genes.annotation_status == 'unique'].copy()
    genes['bin_index'] = [lookup.get((str(c), int(t) // 30000 * 30000), -1)
                          for c, t in zip(genes.chromosome, genes.tss0)]
    genes = genes[genes.bin_index >= 0].copy()
    if not genes.gene_id.is_unique:
        raise ValueError('Gene coordinates must be unique.')
    rna = pd.read_csv(args.rna, sep='\t', low_memory=False)
    rna = rna[np.isfinite(rna.pvalue) & np.isfinite(rna.log2FoldChange)]
    rna = rna.merge(genes[['gene_id', 'chromosome', 'tss0', 'bin_index']], on='gene_id', validate='many_to_one')
    np.savetxt(tables / '64_library_sizes.tsv', a['library_sizes'], delimiter='\t', fmt='%.17g')
    region_rows, rna_rows, go_rows = [], [], []
    for contrast, (baseline, target) in zip(CONTRASTS, PAIRS):
        signature = (support[:, baseline] > 0).astype(int) + 2 * (support[:, target] > 0)
        members, owner = [], np.full(n, -1, int)
        for i in np.flatnonzero(signature):
            if (members and i == members[-1][-1] + 1 and a['chromosome'][i] == a['chromosome'][i - 1]
                    and a['start'][i] == a['start'][i - 1] + 30000 and signature[i] == signature[i - 1]):
                members[-1].append(int(i))
            else:
                members.append([int(i)])
        rows = []
        for k, indices in enumerate(members):
            owner[indices] = k
            i = indices[0]
            row = dict(contrast=contrast, region_index=k, region_id=f'{contrast}_seg{k+1:05}',
                       chromosome=str(a['chromosome'][i]), start0=int(a['start'][i]),
                       end0=int(a['start'][indices[-1]] + 30000), bins=len(indices),
                       clean_bp=int(a['clean_bp'][indices].sum()), membership_signature=int(signature[i]),
                       baseline_support_min=int(support[indices, baseline].min()),
                       baseline_support_max=int(support[indices, baseline].max()),
                       target_support_min=int(support[indices, target].min()),
                       target_support_max=int(support[indices, target].max()))
            rows.append(row)
            region_rows.append(row)
        counts = np.array([a['counts'][indices].sum(0) for indices in members])
        values = np.array([np.average(z[indices], weights=a['clean_bp'][indices], axis=0) for indices in members])
        np.savetxt(tables / f'64_{contrast}_counts.tsv', counts, delimiter='\t', fmt='%.17g')
        np.savetxt(tables / f'64_{contrast}_condition_depth.tsv', values, delimiter='\t', fmt='%.17g')
        np.savez_compressed(tables / f'64_{contrast}_arrays.npz', region_of_bin=owner, signature=signature,
                            support=support, counts=counts, condition_depth=values)
        for gene in rna[rna.contrast == contrast].itertuples():
            bi = int(gene.bin_index)
            ri = int(owner[bi])
            if ri < 0:
                continue
            rna_rows.append(dict(contrast=contrast, gene_id=gene.gene_id, gene_name=gene.gene_name,
                                 chromosome=gene.chromosome, tss0=int(gene.tss0), bin_index=bi,
                                 region_id=f'{a["chromosome"][bi]}:{int(a["start"][bi])}-{int(a["start"][bi]) + 30000}',
                                 RNA_log2FC=gene.log2FoldChange, RNA_pvalue=gene.pvalue, RNA_padj=gene.padj,
                                 RNA_classification=gene.classification, pair_region_index=ri,
                                 pair_region_id=rows[ri]['region_id'], membership_signature=int(signature[bi]),
                                 baseline_support=int(support[bi, baseline]), target_support=int(support[bi, target])))
        for gene in genes.itertuples():
            bi = int(gene.bin_index)
            ri = int(owner[bi])
            if ri < 0:
                continue
            tss = int(gene.tss0)
            go_rows.append(dict(contrast=contrast, gene_id=gene.gene_id, gene_name=gene.gene_name,
                                gene_biotype=gene.gene_biotype, chromosome=gene.chromosome, tss0=tss,
                                bin_index=bi, pair_region_index=ri, pair_region_id=rows[ri]['region_id'],
                                gene_length=int(gene.end0) - int(gene.start0), membership_signature=int(signature[bi]),
                                baseline_support=int(support[bi, baseline]), target_support=int(support[bi, target]),
                                segment_edge_distance_bp=min(tss - rows[ri]['start0'], rows[ri]['end0'] - tss - 1)))
    for name, rows in [('64_regions.tsv', region_rows), ('64_gene_RNA_links.tsv', rna_rows),
                       ('64_gene_GO_links.tsv', go_rows)]:
        pd.DataFrame(rows).to_csv(tables / name, sep='\t', index=False)
    print(f'Prepared {len(region_rows):,} segments across three LAD comparisons.')


if __name__ == '__main__':
    main()

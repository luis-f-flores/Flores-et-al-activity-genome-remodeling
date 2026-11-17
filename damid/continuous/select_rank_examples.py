"""Select rank-band callouts and annotate their overlapping genes."""
import numpy as np
import pandas as pd

# Locus panels in Supplementary Figure 1.
LOCUS_WINDOWS = [('Bdnf', 'chr2', 108800000, 110800000),
                 ('Dlg2', 'chr7', 90800000, 92800000)]


def select(windows, annotation):
    genes = pd.read_csv(annotation, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'gene', 'score', 'strand'])
    genes = genes.groupby(['gene', 'chrom'], as_index=False).agg(start=('start', 'min'), end=('end', 'max'))
    labels = []
    anchors = pd.DataFrame(LOCUS_WINDOWS, columns=['gene', 'chrom', 'start0', 'end0'])
    anchors = anchors.merge(windows, on=['chrom', 'start0', 'end0'], validate='one_to_many')
    for time in ('6 h', 'Recovery'):
        data = windows[(windows.time == time) & windows.eligible].copy()
        data['chrom_number'] = data.chrom.str.replace('chr', '', regex=False).astype(int)
        for target in (1, 3, 5, 7, 9, 91, 93, 95, 97, 99):
            candidates = data[data.percentile.ge(target-1) & data.percentile.lt(target+1)].copy()
            candidates['tie_distance'] = (candidates.percentile-target).abs().round(12)
            anchor = anchors[(anchors.time == time) & anchors.percentile.ge(target-1) & anchors.percentile.lt(target+1)]
            if len(anchor):
                a = anchor.iloc[0]
                row = candidates[(candidates.chrom == a.chrom) & (candidates.start0 == a.start0) & (candidates.end0 == a.end0)].iloc[0].copy()
                row['anchor_gene'] = a.gene
            else:
                row = candidates.sort_values(['tie_distance', 'chrom_number', 'start0', 'end0']).iloc[0].copy()
                row['anchor_gene'] = ''
            row['target_percentile'] = target
            row['tail'] = 'Loss' if target < 10 else 'Gain'
            overlap = genes[(genes.chrom == row.chrom) & (genes.start < row.end0) & (genes.end > row.start0)].copy()
            overlap['center_distance'] = abs((overlap.start+overlap.end)/2-(row.start0+row.end0)/2)
            overlap = overlap.sort_values(['center_distance', 'gene'])
            row['example_genes'] = ' / '.join(overlap.gene.head(2)) if len(overlap) else 'No annotated gene'
            if row.anchor_gene: row['example_genes'] = row.anchor_gene
            labels.append(row)
    labels = pd.DataFrame(labels)
    context = labels[labels.anchor_gene.eq('')]
    loci = list(LOCUS_WINDOWS)
    for ascending in (True, False):
        row = context.sort_values(['mean_change', 'time', 'chrom_number', 'start0', 'end0'], ascending=[ascending, True, True, True, True]).iloc[0]
        loci.append((row.example_genes.split(' / ')[0], row.chrom, int(row.start0), int(row.end0)))
    shown = []
    for (time, tail), rows in labels.groupby(['time', 'tail'], sort=False):
        targets = [91, 95, 97] if time == '6 h' and tail == 'Gain' else [1, 5, 9] if tail == 'Loss' else [91, 95, 99]
        shown.append(rows[rows.target_percentile.isin(targets)])
    selected = pd.DataFrame(loci, columns=['gene', 'chrom', 'start0', 'end0'])
    selected = selected.merge(windows, on=['chrom', 'start0', 'end0'], validate='one_to_many')
    return pd.concat(shown, ignore_index=True), selected

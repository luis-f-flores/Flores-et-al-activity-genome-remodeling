"""Prepare the three DMR-example tracks from CpG bigWigs and merged DMR calls."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from locus_mcg import extract


def prepare(matrix, regions, dmrs, bigwigs, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    data = np.load(matrix)
    union = pd.read_csv(regions, sep='\t')
    calls = pd.concat([pd.read_csv(path, sep='\t').assign(source_file=Path(path).name)
                       for path in dmrs if Path(path).name.startswith('1M_')])
    ids = list(data['samples'])
    values, n = data['mCG_percent'], data['covered_CpG_records']
    im = [ids.index(rep+'1M') for rep in 'ABC']
    ik = [ids.index(rep+'1K') for rep in 'ABC']
    i0 = [ids.index(rep+'0') for rep in 'ABC']
    union['mean_M'] = values[:, im].mean(axis=1)
    union['mean_K'] = values[:, ik].mean(axis=1)
    union['delta'] = union.mean_K-union.mean_M
    union['mean0'] = values[:, i0].mean(axis=1)
    union['minCpG'] = n[:, im+ik+i0].min(axis=1)
    union['rep_separation'] = np.where(union.delta > 0,
        values[:, ik].min(axis=1)-values[:, im].max(axis=1),
        values[:, im].min(axis=1)-values[:, ik].max(axis=1))
    candidates = union.merge(calls, on=['chr', 'start', 'end'])
    candidates = candidates[(candidates.end-candidates.start+1).between(750, 1500)
        & candidates.delta.abs().between(15, 30) & (candidates.rep_separation > 3)
        & (candidates.minCpG >= 12) & candidates.chr.str.match(r'^chr([1-9]|1[0-9])$')].copy()
    candidates['median_distance'] = (candidates.delta.abs()-candidates.delta.abs().median()).abs()
    candidates = candidates.sort_values(['median_distance', 'chr', 'start'])
    selected = candidates.iloc[0]
    candidates.to_csv(output/'one_hour_dmr_candidates.tsv', sep='\t', index=False)
    windows = [('Adssl1', 'chr12', 112628750, 112630000, 1),
               ('Bdnf', 'chr2', 109692150, 109693150, 6),
               ('Styx', 'chr14', 45349099, 45353609, 24)]
    if not (selected.chr == windows[0][1] and selected.start >= windows[0][2]
            and selected.end <= windows[0][3]):
        raise ValueError('The selected one-hour DMR falls outside the configured example window.')
    samples, tracks = extract(bigwigs, windows)
    samples.to_csv(output/'locus_sample_values.tsv', sep='\t', index=False)
    tracks.to_csv(output/'three_example_plot_values.tsv', sep='\t', index=False)
    other = []
    for gene, chrom, lo, hi, time in windows[1:]:
        for path in dmrs:
            if Path(path).name.startswith(str(time)+'M_'):
                table = pd.read_csv(path, sep='\t')
                table = table[(table.chr == chrom) & (table.start <= hi) & (table.end >= lo)].copy()
                table['direction'] = 'hyper' if '_hyper.' in str(path) else 'hypo'
                table['gene'] = gene
                other.append(table)
    other = pd.concat(other, ignore_index=True)
    first = pd.DataFrame([dict(gene='Adssl1', chr=selected.chr, start=int(selected.start),
        end=int(selected.end), meth_diff=float(selected.meth_diff), pvalue=float(selected.pvalue),
        qvalue=float(selected.qvalue), direction='hyper', source_file=selected.source_file)])
    result = pd.concat([first.reindex(columns=list(other.columns)+['source_file']), other], ignore_index=True)
    result.to_csv(output/'three_example_DMRs.tsv', sep='\t', index=False)
    index = union.index[union.region_id == selected.region_id][0]
    pd.DataFrame({'sample': [ids[k] for k in i0+im+ik],
                  'mCG_percent': values[index, i0+im+ik],
                  'covered_CpG_records': n[index, i0+im+ik].astype(int)}).to_csv(
        output/'one_hour_region_replicates.tsv', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--regions', type=Path, required=True)
    parser.add_argument('--dmrs', type=Path, nargs=6, required=True)
    parser.add_argument('--bigwigs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    a = parser.parse_args()
    prepare(a.matrix, a.regions, a.dmrs, a.bigwigs, a.output)

"""Generate the published B1 locus tracks from fragment counts and RNA bigWigs."""
from pathlib import Path
import argparse
import shutil
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import rankdata

AUTOSOMES = ['chr' + str(i) for i in range(1, 20)]
CHROMOSOMES = AUTOSOMES + ['chrX', 'chrY']
IPS = ['L9', 'L12', 'L11', 'L14', 'L16', 'L10', 'L13', 'L15']
INPUTS = ['L4', 'L6', 'L8', 'L2', 'L5', 'L7']
POOLS = [[0], [0], [1, 2], [1, 2], [1, 2], [3, 4, 5], [3, 4, 5], [3, 4, 5]]
GROUPS = [[0, 1], [2, 3, 4], [5, 6, 7]]
WINDOWS = 3 * np.unique(2 * np.rint((np.geomspace(1, 33, 28) - 1) / 2).astype(int) + 1)


def fine_vectors(auto_path, sex_path, counts, sex_counts):
    auto, sex = dict(np.load(auto_path)), dict(np.load(sex_path))
    ids = INPUTS + IPS
    per = {s: dict(np.load(counts / (s + '.counts.npz'))) for s in ids}
    xy = {s: dict(np.load(sex_counts / (s + '.XY_counts.npz'))) for s in ids}
    totals = {s: sum(int(per[s][c].sum()) for c in AUTOSOMES) for s in ids}
    h, t = np.array([totals[s] for s in IPS]), np.array([totals[s] for s in INPUTS])
    cc, ss, xx, bb = [], [], [], []
    for chrom in CHROMOSOMES:
        ref, source = (auto, per) if chrom in AUTOSOMES else (sex, xy)
        length = int(ref['lengths'][list(ref['chromosomes']).index(chrom)])
        n = length // 10000
        x = np.column_stack([source[s][chrom][:n * 2].reshape(n, 2).sum(1) for s in ids])
        bp = ref[chrom + '_clean_bp'][:n * 2].reshape(n, 2).sum(1)
        keep = (bp > 0) & (x[:, :6].min(1) >= 10) & (x[:, 6:].min(1) > 0)
        cc.extend([chrom] * int(keep.sum()))
        ss.extend((np.flatnonzero(keep) * 10000).tolist())
        xx.append(x[keep]); bb.extend(bp[keep].tolist())
    x = np.concatenate(xx)
    denominator = np.column_stack([x[:, pool].sum(1) / t[pool].sum() for pool in POOLS])
    raw = np.log2((x[:, 6:] / h) / denominator)
    return np.array(cc), np.array(ss), raw, np.array(bb)


def domainograms(chrom, start, raw):
    means = np.column_stack([raw[:, group].mean(1) for group in GROUPS])
    segments = np.split(np.arange(len(start)), np.flatnonzero(
        (chrom[1:] != chrom[:-1]) | (np.diff(start) != 10000)) + 1)
    colors = np.full((3, len(WINDOWS), len(start)), np.nan)
    for wi, width in enumerate(WINDOWS):
        half = int(width) // 2
        local = np.full(means.shape, np.nan)
        for seg in segments:
            if len(seg) < width:
                continue
            centers = seg[half:len(seg) - half]
            cs = np.vstack([np.zeros((1, 3)), np.cumsum(means[seg], axis=0)])
            local[centers] = (cs[width:] - cs[:-width]) / width
        valid = np.isfinite(local).all(1)
        for ci, (baseline, target) in enumerate([(0, 1), (0, 2), (1, 2)]):
            qa, qb = np.quantile(local[valid][:, [baseline, target]], .3, axis=0)
            eligible = valid & ~((local[:, baseline] <= qa) & (local[:, target] <= qb))
            delta = local[:, target] - local[:, baseline]
            q = np.full(len(start), np.nan)
            q[eligible] = rankdata(delta[eligible], method='average') / eligible.sum()
            z = np.zeros(len(start)); z[~np.isfinite(q)] = np.nan
            low, high = (q < .1) & (delta < 0), (q > .9) & (delta > 0)
            z[low] = (q[low] - .1) / .1; z[high] = (q[high] - .9) / .1
            colors[ci, wi] = z
    return colors


def curves(chrom, starts, raw, lookup, gene_chrom, grid, width):
    out = np.full((len(grid), 3), np.nan)
    half = width // 2
    for j, start in enumerate(grid):
        ix = [lookup.get((str(gene_chrom), int(start) + 10000 * k)) for k in range(-half, half + 1)]
        if any(i is None for i in ix):
            continue
        ip = np.log2(np.mean(2 ** raw[ix], axis=0))
        out[j] = [ip[group].mean() for group in GROUPS]
    return out


def rna_values(bigwigs, chromosome, grid):
    values = np.full((len(grid), len(bigwigs)), np.nan)
    lo, hi = int(grid[0]), int(grid[-1]) + 50
    for k, bw in enumerate(bigwigs):
        for left, right, value in bw.intervals(chromosome, lo, hi) or []:
            start = max(0, (left - lo) // 50)
            stop = min(len(grid), (right - lo + 49) // 50)
            values[start:stop, k] = value
    return values


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['autosomal-reference', 'sex-reference', 'counts-dir', 'sex-counts-dir',
                'signals', 'decoded', 'gene-responses', 'gene-metadata', 'exons', 'rna-tracks', 'output']:
        p.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=True)
    ch, st, raw, bp = fine_vectors(a.autosomal_reference, a.sex_reference, a.counts_dir, a.sex_counts_dir)
    colors = domainograms(ch, st, raw)
    lookup = {(str(c), int(s)): i for i, (c, s) in enumerate(zip(ch, st))}
    signal = np.load(a.signals)
    states = np.load(a.decoded)['states']
    calls = np.column_stack([states[:, group].max(1) for group in GROUPS])
    genes = pd.read_csv(a.gene_responses, sep='\t', low_memory=False).set_index('gene_name')
    metadata = pd.read_csv(a.gene_metadata, sep='\t').set_index('gene_name')
    bws = [pyBigWig.open(str(a.rna_tracks / (name + '.meanCPM50.mm10.bw')))
           for name in ['TTX', 'Bic', 'Bic_Nimo']]
    arrays, rows, selected = {}, [], []
    for name in ['Ccbe1', 'Grin2a']:
        r = genes.loc[name]
        if isinstance(r, pd.DataFrame):
            raise ValueError('Locus name is not unique: ' + name)
        if name == 'Ccbe1':
            lo = int((r.start0 - 50000) // 10000 * 10000)
            hi = int((r.end0 + 50000 + 9999) // 10000 * 10000)
        else:
            center = (int(r.start0) + int(r.end0)) // 2
            lo = (center - 280000) // 50 * 50; hi = lo + 560000
        grid = np.arange(lo // 10000 * 10000, (hi + 9999) // 10000 * 10000, 10000)
        curve = curves(ch, st, raw, lookup, r.chromosome, grid, 3)
        z = np.full((2, len(WINDOWS), len(grid)), np.nan)
        for j, start in enumerate(grid):
            k = lookup.get((r.chromosome, int(start)))
            if k is not None:
                z[:, :, j] = colors[:2, :, k]
        rx = np.arange(lo // 50 * 50, (hi + 49) // 50 * 50, 50)
        rna = rna_values(bws, r.chromosome, rx)
        visible = (rx + 50 > r.start0) & (rx < r.end0)
        rmax = float(np.ceil(np.nanmax(rna[visible]) * 1.05 * 10) / 10)
        take = ((signal['chromosome'] == r.chromosome) & (signal['start'] < hi) &
                (signal['start'] + 30000 > lo))
        arrays.update({name + '_starts': grid, name + '_domainograms': z,
                       name + '_RNA_starts': rx, name + '_RNA_CPM': rna,
                       name + '_HMM_starts': signal['start'][take], name + '_HMM_states': calls[take],
                       name + '_delta': curve[:, 1:] - curve[:, [0]], name + '_RNA_display_mask': visible})
        selected.append(dict(r, gene_name=name))
        rows.append(dict(gene=name, chromosome=r.chromosome, view_start0=lo, view_end0=hi,
                         RNA_start0=lo, RNA_end0=hi, RNA_ymax=rmax, gene_start0=r.start0,
                         gene_end0=r.end0, Bic_RNA_log2FC=r.fc_Bic_vs_TTX, Bic_RNA_q=r.q_Bic_vs_TTX,
                         Bic_gene_B1_delta=r.b, combined_gene_B1_delta=r.bn))
    for bw in bws:
        bw.close()
    np.savez_compressed(a.output / 'main.npz', **arrays)
    pd.DataFrame(selected).to_csv(a.output / 'main_genes.tsv', sep='\t', index=False)
    pd.DataFrame(rows).to_csv(a.output / 'main_summary.tsv', sep='\t', index=False)
    arrays, selected = {}, []
    for name in ['Bdnf', 'Mpp7']:
        r = metadata.loc[name]
        if isinstance(r, pd.DataFrame):
            raise ValueError('Locus name is not unique: ' + name)
        lo, hi = int(r.start0) - 400000, int(r.end0) + 400000
        grid = np.arange(lo // 10000 * 10000, (hi + 9999) // 10000 * 10000, 10000)
        arrays[name + '__B1_starts'] = grid
        arrays[name + '__B1_log2'] = curves(ch, st, raw, lookup, r.chrom, grid, 9)
        selected.append(dict(gene_name=name, chromosome=r.chrom, start0=r.start0, end0=r.end0))
    np.savez_compressed(a.output / 'supplement.npz', **arrays)
    pd.DataFrame(selected).to_csv(a.output / 'supplement_genes.tsv', sep='\t', index=False)
    np.savez_compressed(a.output / 'raw_10kb.npz', chromosome=ch, start=st, condition_depth=raw)
    np.save(a.output / 'window_bins.npy', WINDOWS)
    shutil.copyfile(a.gene_metadata, a.output / 'gene_metadata.tsv')
    shutil.copyfile(a.exons, a.output / 'exons.bed')


if __name__ == '__main__':
    main()

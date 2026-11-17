"""Construct the TTX B1 grid and LAD methylation tables from generated measurements."""
from pathlib import Path
from io import StringIO
import argparse
import numpy as np
import pandas as pd


def write_grid(signals, decoded, output):
    values = np.load(signals)
    states = np.load(decoded)['states']
    chromosomes = values['chromosome']
    keep = np.isin(chromosomes, ['chr' + str(i) for i in range(1, 20)])
    if states.shape != values['condition_depth'].shape:
        raise ValueError('B1 values and HMM states must use the same bins and libraries')
    signal = values['condition_depth'][keep]
    groups = [[0, 1], [2, 3, 4], [5, 6, 7]]
    means = np.column_stack([signal[:, group].mean(axis=1) for group in groups])
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / 'B1_bins.npz', chromosome=chromosomes[keep],
                        start=values['start'][keep], means=means, clean_bp=values['clean_bp'][keep])
    np.savez_compressed(output / 'TTX_bins.npz', signal=signal[:, :2], state=states[keep, :2].max(axis=1))


def write_tables(grid, baseline, output):
    values = np.load(grid / 'B1_bins.npz')
    ttx = np.load(grid / 'TTX_bins.npz')
    methylation = np.load(baseline / 'baseline_counts.npz')
    chrom, start, bp = (values[key] for key in ['chromosome', 'start', 'clean_bp'])
    b1, state = ttx['signal'], ttx['state']
    meth, depth, sites = (methylation[key] for key in ['meth', 'depth', 'sites'])
    if meth.shape != (len(start), 3) or depth.shape != meth.shape or sites.shape != meth.shape:
        raise ValueError('Baseline CpG counts must use the B1 grid and three samples')
    with np.errstate(divide='ignore', invalid='ignore'):
        per = 100 * meth / depth
    valid = (sites >= 3).all(axis=1) & np.isfinite(per).all(axis=1)
    table = pd.DataFrame(dict(chromosome=chrom, start0=start, end0=start + 30000,
                              TTX_L9=b1[:, 0], TTX_L12=b1[:, 1], TTX_mean=values['means'][:, 0],
                              TTX_state=state, eligible=valid, mCG_mean=per.mean(axis=1)))
    for j, sample in enumerate(['A0', 'B0', 'C0']):
        table[sample + '_mCG'] = per[:, j]
        table[sample + '_sites'] = sites[:, j]
        table[sample + '_coverage_sum'] = depth[:, j]
    rows = []
    for chromosome in dict.fromkeys(chrom):
        indices = np.flatnonzero((chrom == chromosome) & (state > 0))
        for group in np.split(indices, np.flatnonzero(np.diff(start[indices]) != 30000) + 1):
            if not len(group):
                continue
            left, right = int(start[group[0]]), int(start[group[-1]]) + 30000
            b1_mean = np.average(b1[group], weights=bp[group], axis=0)
            m, n, c = meth[group].sum(axis=0), depth[group].sum(axis=0), sites[group].sum(axis=0)
            with np.errstate(divide='ignore', invalid='ignore'):
                percent = 100 * m / n
            row = dict(LAD_id=f'TTX_LAD_{len(rows) + 1:04d}', chromosome=chromosome,
                       start0=left, end0=right, length_bp=right - left, clean_bp=int(bp[group].sum()),
                       bins=len(group), B1_TTX_mean=float(b1_mean.mean()),
                       mCG_baseline_mean=float(percent.mean()), B1_L9=float(b1_mean[0]), B1_L12=float(b1_mean[1]))
            for j, sample in enumerate(['A0', 'B0', 'C0']):
                row.update({sample + '_mCG': float(percent[j]), sample + '_methylated': int(m[j]),
                            sample + '_coverage': int(n[j]), sample + '_CpG_records': int(c[j])})
            rows.append(row)
    output.mkdir(parents=True, exist_ok=True)
    # Preserve the TSV reading precision used for the plotted quantile groups.
    table = pd.read_csv(StringIO(table.to_csv(sep='\t', index=False)), sep='\t')
    table.to_csv(output / 'B1_mCG_bins.tsv.gz', sep='\t', index=False)
    pd.DataFrame(rows).to_csv(output / 'TTX_LADs.tsv', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    grid = commands.add_parser('grid')
    grid.add_argument('--signals', type=Path, required=True)
    grid.add_argument('--decoded', type=Path, required=True)
    grid.add_argument('--output', type=Path, required=True)
    table = commands.add_parser('tables')
    table.add_argument('--grid', type=Path, required=True)
    table.add_argument('--baseline', type=Path, required=True)
    table.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'grid':
        write_grid(args.signals, args.decoded, args.output)
    else:
        write_tables(args.grid, args.baseline, args.output)

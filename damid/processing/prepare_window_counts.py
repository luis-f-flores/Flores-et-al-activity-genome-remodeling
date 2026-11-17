#!/usr/bin/env python3
"""Prepare S1 compact inputs from the original R GATC and 400-kb count files.

The original per-feature table uses zero-based start and one-based end.
Reconstruct the source midpoint exactly, and require the site counts to
sum to every supplied 400-kb count before exporting count/site arrays.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'

import argparse
from pathlib import Path
import shutil

import numpy as np
import pandas as pd



def aggregate_sites(path, coords):
    index = {(row.chrom, int(row.start) // 400000): i for i, row in coords.iterrows()}
    total = np.zeros(len(coords), dtype=np.int64)
    positive = np.zeros(len(coords), dtype=np.int64)
    for chunk in pd.read_csv(path, sep='\t', header=None, names=['chrom', 'start', 'end', 'count'],
                             chunksize=500000, dtype={'chrom': 'category', 'start': 'int64', 'end': 'int64', 'count': 'int64'}):
        if (chunk['count'] < 0).any() or (chunk.end <= chunk.start).any():
            raise ValueError('Invalid GATC coordinates/counts')
        for chrom, group in chunk.groupby('chrom', observed=True):
            midpoint = group.start.to_numpy() + (group.end.to_numpy() - group.start.to_numpy() - 1) // 2
            local = midpoint // 400000
            dest = np.array([index.get((chrom, int(value)), -1) for value in local])
            valid = dest >= 0
            values = group['count'].to_numpy()[valid]
            np.add.at(total, dest[valid], values)
            np.add.at(positive, dest[valid], (values > 0).astype(np.int64))
    return total, positive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', type=Path, required=True, help='Directory with sample-400kb.counts.txt and sample-gatc.counts.txt[.gz]')
    parser.add_argument('--output', type=Path, required=True, help='A new compact-input directory')
    args = parser.parse_args()
    coords = pd.read_csv(args.counts / '6A-DBM-400kb.counts.txt', sep='\t', header=None, names=['chrom', 'start', 'end', 'count'])[['chrom', 'start', 'end']]
    arrays, libraries = {}, []
    for hours, time in ((6, '6 h'), (24, 'Recovery')):
        for condition, code in (('Mannitol', 'M'), ('KCl', 'K')):
            for rep in ('A', 'B'):
                for role, prefix in (('Fusion', 'DB'), ('Dam-only', 'D')):
                    sample = f'{hours}{rep}-{prefix}{code}'
                    bins = pd.read_csv(args.counts / f'{sample}-400kb.counts.txt', sep='\t', header=None,
                                       names=['chrom', 'start', 'end', 'count'])
                    if not bins[['chrom', 'start', 'end']].equals(coords):
                        raise ValueError(f'400-kb grid mismatch: {sample}')
                    path = args.counts / f'{sample}-gatc.counts.txt.gz'
                    if not path.exists():
                        path = path.with_suffix('')
                    count, positive = aggregate_sites(path, coords)
                    np.testing.assert_array_equal(count, bins['count'].to_numpy())
                    arrays[sample + '_count'], arrays[sample + '_positive'] = count, positive
                    libraries.append(dict(sample=sample, time=time, condition=condition, replicate=rep,
                                          role=role, total_assigned_counts=int(count.sum())))
                    print(sample, 'all site/bin counts agree', flush=True)
    args.output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(args.output / 'site_support_400kb.npz', **arrays)
    pd.DataFrame(libraries).to_csv(args.output / 'library_support.tsv', sep='\t', index=False)
    coords.to_csv(args.output / 'coordinates_400kb.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()

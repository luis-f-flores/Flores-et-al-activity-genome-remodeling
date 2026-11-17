"""Count deposited unique read coordinates in complete autosomal mm9 5-kb windows."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def count_reads(source, sizes_file, output, tag):
    sizes = dict(line.split() for line in sizes_file.read_text().splitlines())
    sizes = {chrom: int(length) for chrom, length in sizes.items()}
    chromosomes = [f'chr{i}' for i in range(1, 20)]
    records, offsets = [], {}
    for chrom in chromosomes:
        offsets[chrom] = len(records)
        records.extend((chrom, start, start + 5000) for start in range(0, sizes[chrom] // 5000 * 5000, 5000))
    grid = pd.DataFrame(records, columns=['chromosome', 'start0', 'end0'])
    counts = np.zeros(len(grid), dtype=np.int64)
    total = autosomal = 0
    for chunk in pd.read_csv(source, sep=r'\s+', header=None, names=['chr', 'strand', 'pos'],
                             dtype={'chr': 'category', 'strand': 'category', 'pos': 'int64'}, chunksize=1000000):
        total += len(chunk)
        assert chunk.strand.isin(['+', '-']).all()
        for chrom, group in chunk.groupby('chr', observed=True):
            if chrom not in offsets:
                continue
            assert group.pos.min() >= 0 and group.pos.max() < sizes[chrom]
            autosomal += len(group)
            indices = group.pos.to_numpy() // 5000
            indices = indices[indices < sizes[chrom] // 5000]
            counts += np.bincount(offsets[chrom] + indices, minlength=len(grid))
    output.mkdir(parents=True, exist_ok=True)
    grid.to_csv(output / 'mm9_5kb_windows.tsv', sep='\t', index=False)
    np.savez_compressed(output / f'{tag}_mm9_5kb.npz', counts=counts,
                        total_unique_reads=total, autosomal_reads=autosomal)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reads', type=Path, required=True)
    parser.add_argument('--sizes', type=Path, required=True)
    parser.add_argument('--tag', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    count_reads(args.reads, args.sizes, args.output, args.tag)

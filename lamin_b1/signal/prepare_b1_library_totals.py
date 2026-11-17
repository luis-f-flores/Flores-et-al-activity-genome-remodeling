"""Combine host and fly fragment totals used for B1 spike-in normalization."""
from pathlib import Path
import argparse
import bisect
import numpy as np
import pandas as pd

SAMPLES = ['L4', 'L6', 'L8', 'L2', 'L5', 'L7', 'L9', 'L12', 'L11', 'L14', 'L16', 'L10', 'L13', 'L15']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['reference', 'counts', 'fly', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    ref = np.load(args.reference)
    fly = pd.read_csv(args.fly, sep='\t').set_index('sample_id')
    usable = {}
    for chrom in ref['chromosomes']:
        chrom = str(chrom)
        starts = ref[chrom + '_mask_start'].tolist()
        ends = ref[chrom + '_mask_end'].tolist()
        actual = ref[chrom + '_actual_bp']
        mask_free = []
        for k, size in enumerate(actual):
            left, right = k * 5000, k * 5000 + int(size)
            j = bisect.bisect_left(starts, right) - 1
            mask_free.append(not (j >= 0 and ends[j] > left))
        usable[chrom] = (ref[chrom + '_ACGT_bp'] / actual >= .8) & np.array(mask_free)
    rows = []
    for sample in SAMPLES:
        counts = np.load(args.counts / f'{sample}.counts.npz')
        host = sum(int((counts[chrom] + counts[chrom + '_nonACGT_midpoints'])[mask].sum())
                   for chrom, mask in usable.items())
        fly_count = (int(fly.loc[sample, 'fly_nonmitochondrial_nonduplicate_proper_pairs_MAPQ10'])
                     if sample in fly.index else 0)
        rows.append(dict(sample_id=sample, host_usable_fragments=host, fly_eligible_fragments=fly_count))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()

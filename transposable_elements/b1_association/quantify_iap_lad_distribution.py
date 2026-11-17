"""Project IAP sequence onto the assessable bases of the B1 genomic grid."""
from pathlib import Path
from collections import defaultdict
from bisect import bisect_right
import argparse
import csv
import gzip
import numpy as np

CHROMS = ['chr' + str(i) for i in range(1, 20)] + ['chrX', 'chrY']


def merged(intervals):
    output = []
    for start, end in sorted(intervals):
        if output and start <= output[-1][1]: output[-1][1] = max(output[-1][1], int(end))
        else: output.append([int(start), int(end)])
    return output


def project(catalog, signals, references):
    refs, lengths = {}, {}
    for path in references:
        a = dict(np.load(path)); refs.update(a)
        lengths.update({str(c): int(n) for c, n in zip(a['chromosomes'], a['lengths'])})
    annotations = defaultdict(list)
    with gzip.open(catalog, 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['is_iap'].upper() == 'TRUE' and row['salmonte_target'] != 'RLTR1IAP_MM':
                annotations[row['chrom']].append((int(row['start_0based']), int(row['end_0based_exclusive'])))
    annotations = {c: merged(values) for c, values in annotations.items()}
    projected, union = {}, np.zeros(len(signals['start']), np.int64)
    for chrom in CHROMS:
        use = np.flatnonzero(signals['chromosome'] == chrom)
        masks = merged(list(zip(refs[chrom+'_mask_start'], refs[chrom+'_mask_end'])) +
                       list(zip(refs[chrom+'_nonACGT_start'], refs[chrom+'_nonACGT_end'])))
        mask_ends = [end for start, end in masks]
        clean = []
        for i in use:
            start = int(signals['start'][i]); end = min(start+30000, lengths[chrom])
            cursor = start; j = bisect_right(mask_ends, start)
            while j < len(masks) and masks[j][0] < end:
                lo, hi = masks[j]
                if cursor < lo: clean.append((cursor, min(lo, end)))
                cursor = max(cursor, hi); j += 1
            if cursor < end: clean.append((cursor, end))
        cs, ce = np.array(clean, np.int64).T
        cumulative = np.r_[0, np.cumsum(ce-cs)]
        widths = signals['clean_bp'][use]
        if cumulative[-1] != widths.sum(): raise ValueError('Assessable sequence lengths differ from B1 grid.')
        intervals = []
        for start, end in annotations.get(chrom, []):
            lo = np.searchsorted(ce, start, side='right'); hi = np.searchsorted(cs, end, side='left')
            for j in range(lo, hi):
                intervals.append((int(cumulative[j]+max(start, cs[j])-cs[j]),
                                  int(cumulative[j]+min(end, ce[j])-cs[j])))
        intervals = merged(intervals)
        starts, ends = np.array(intervals, np.int64).T if intervals else (np.array([], np.int64), np.array([], np.int64))
        boundaries = np.r_[0, np.cumsum(widths)]
        for start, end in intervals:
            lo = np.searchsorted(boundaries[1:], start, side='right')
            hi = np.searchsorted(boundaries[:-1], end, side='left')
            union[use[lo:hi]] += np.minimum(end, boundaries[lo+1:hi+1])-np.maximum(start, boundaries[lo:hi])
        projected[chrom+'_IAP_start'] = starts
        projected[chrom+'_IAP_end'] = ends
        projected[chrom+'_bin_width'] = widths
    if (union > signals['clean_bp']).any(): raise ValueError('IAP union exceeds assessable sequence.')
    return projected, union


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['catalog', 'signals', 'autosomal-reference', 'sex-reference', 'output']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=True)
    projected, union = project(a.catalog, np.load(a.signals), [a.autosomal_reference, a.sex_reference])
    np.savez_compressed(a.output/'IAP_projection.npz', **projected)
    np.savez_compressed(a.output/'IAP_union_bp_by_bin.npz', all_iap=union)


if __name__ == '__main__':
    main()

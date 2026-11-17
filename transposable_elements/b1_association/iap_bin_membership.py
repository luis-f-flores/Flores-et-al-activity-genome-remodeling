"""Map the original IAP catalogue onto supported, unmasked 30-kb B1 bins."""
from pathlib import Path
import argparse,bisect,csv,gzip,re
import numpy as np


def construct(catalog, signal_arrays, mask_paths, output):
    a = dict(np.load(signal_arrays))
    refs = {}
    for path in mask_paths:
        refs.update(dict(np.load(path)))
    lookup = {(str(c), int(s)): i for i, (c, s) in enumerate(zip(a['chromosome'], a['start']))}
    masks = {}
    for chrom in set(a['chromosome']):
        intervals = sorted((int(s), int(e)) for kind in ['mask', 'nonACGT']
                           for s, e in zip(refs[chrom + '_' + kind + '_start'], refs[chrom + '_' + kind + '_end']))
        merged = []
        for start, end in intervals:
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(end, merged[-1][1])
            else:
                merged.append([start, end])
        masks[chrom] = ([s for s, _ in merged], [e for _, e in merged])
    targets = ['IAPEY4_I', 'IAPA_MM', 'IAPEZI']
    normalize = lambda name: re.sub('[^A-Za-z0-9]', '', name).lower()
    names = {normalize(name): name for name in targets}
    members = {name: set() for name in targets}
    all_iap = np.zeros(len(a['start']), bool)
    with gzip.open(catalog, 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['is_iap'].upper() != 'TRUE':
                continue
            chrom, start, end = row['chrom'], int(row['start_0based']), int(row['end_0based_exclusive'])
            if chrom not in masks:
                continue
            starts, ends = masks[chrom]
            for bin_start in range(start // 30000 * 30000, ((end - 1) // 30000 + 1) * 30000, 30000):
                index = lookup.get((chrom, bin_start))
                if index is None:
                    continue
                lo, hi = max(start, bin_start), min(end, bin_start + 30000)
                j, masked = bisect.bisect_right(ends, lo), 0
                while j < len(starts) and starts[j] < hi:
                    masked += max(0, min(hi, ends[j]) - max(lo, starts[j])); j += 1
                if hi - lo - masked <= 0:
                    continue
                all_iap[index] = True
                target = names.get(normalize(row['salmonte_target']))
                if target:
                    members[target].add(index)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / 'memberships.npz', **{key: np.array(sorted(value), int) for key, value in members.items()})
    np.savez_compressed(output / 'bins.npz', chromosome=a['chromosome'], start=a['start'],
                        unmasked_bp=a['clean_bp'], B1=a['condition_depth'], IP_ids=a['IP_ids'],
                        TTX_state=a['states'][:, :2].max(axis=1), catalog_IAP=all_iap)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog', type=Path, required=True)
    parser.add_argument('--signal-arrays', type=Path, required=True)
    parser.add_argument('--mask', type=Path, required=True, action='append')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    construct(args.catalog, args.signal_arrays, args.mask, args.output)

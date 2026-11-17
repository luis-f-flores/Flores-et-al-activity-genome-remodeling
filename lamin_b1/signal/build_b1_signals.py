"""Aggregate supported 30-kb B1/Input measurements from 5-kb fragment counts."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np

AUTOSOMES = ['chr' + str(i) for i in range(1, 20)]
CHROMOSOMES = AUTOSOMES + ['chrX', 'chrY']
IP_IDS = ['L9', 'L12', 'L11', 'L14', 'L16', 'L10', 'L13', 'L15']
INPUT_IDS = ['L4', 'L6', 'L8', 'L2', 'L5', 'L7']
SAMPLE_IDS = INPUT_IDS + IP_IDS
INPUT_POOLS = [[0], [0], [1, 2], [1, 2], [1, 2], [3, 4, 5], [3, 4, 5], [3, 4, 5]]


def build(autosomal_reference, sex_reference, counts_dir, sex_counts_dir, library_summary, output):
    auto = dict(np.load(autosomal_reference))
    sex = dict(np.load(sex_reference))
    per = {s: dict(np.load(counts_dir / (s + '.counts.npz'))) for s in SAMPLE_IDS}
    xy = {s: dict(np.load(sex_counts_dir / (s + '.XY_counts.npz'))) for s in SAMPLE_IDS}
    with library_summary.open() as handle:
        qc = {r['sample_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    totals = {s: sum(int(per[s][c].sum()) for c in AUTOSOMES) for s in SAMPLE_IDS}
    h = np.array([totals[s] for s in IP_IDS])
    t = np.array([totals[s] for s in INPUT_IDS])
    fly = np.array([int(qc[s]['fly_eligible_fragments']) for s in IP_IDS])
    initial_host = np.array([int(qc[s]['host_usable_fragments']) for s in IP_IDS])
    if min(h.min(), t.min(), fly.min(), initial_host.min()) <= 0:
        raise ValueError('Normalization totals must be positive.')
    anchor = float(np.exp(np.log(initial_host / fly).mean()))
    parts, chrom, start, bp = [], [], [], []
    for c in CHROMOSOMES:
        ref = auto if c in AUTOSOMES else sex
        source = per if c in AUTOSOMES else xy
        length = int(ref['lengths'][list(ref['chromosomes']).index(c)])
        n = length // 30000
        counts = np.column_stack([source[s][c][:n * 6].reshape(n, 6).sum(1) for s in SAMPLE_IDS])
        clean = ref[c + '_clean_bp'][:n * 6].reshape(n, 6).sum(1)
        support = (clean > 0) & (counts[:, :6].min(1) >= 10) & (counts[:, 6:].min(1) > 0)
        parts.append(counts[support])
        chrom.extend([c] * int(support.sum()))
        start.extend((np.flatnonzero(support) * 30000).tolist())
        bp.extend(clean[support].tolist())
    counts = np.concatenate(parts)
    chrom, start, bp = np.array(chrom), np.array(start), np.array(bp)
    inp = np.column_stack([counts[:, pool].sum(1) / t[pool].sum() for pool in INPUT_POOLS])
    depth = (counts[:, 6:] / h) / inp
    calibrated = (counts[:, 6:] / fly) / (inp * anchor)
    if not (np.isfinite(depth).all() and (depth > 0).all()):
        raise ValueError('B1/Input measurements must be finite and positive.')
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / 'condition_30kb.npz', counts=counts, chromosome=chrom,
                        start=start, clean_bp=bp, fly_linear=calibrated, sample_ids=SAMPLE_IDS)
    np.savez_compressed(output / 'shared_hmm_input.npz', chromosome=chrom, start=start,
                        end=start + 30000, clean_bp=bp, x=depth, IP_ids=IP_IDS)
    np.savez_compressed(output / 'b1_signal.npz', chromosome=chrom, start=start, clean_bp=bp,
                        IP_ids=IP_IDS, counts=counts[:, 6:], library_sizes=h,
                        condition_depth=np.log2(depth), condition_fly=np.log2(calibrated))
    normalization = dict(calibration_anchor=anchor, autosomal_IP_totals=h.tolist(),
                         autosomal_Input_totals=t.tolist(), fly_fragments=fly.tolist())
    (output / 'normalization.json').write_text(json.dumps(normalization, indent=2) + '\n')
    print(f'Wrote {len(start):,} supported 30-kb bins for {len(IP_IDS)} B1 libraries.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--autosomal-reference', type=Path, required=True)
    parser.add_argument('--sex-reference', type=Path, required=True)
    parser.add_argument('--counts-dir', type=Path, required=True)
    parser.add_argument('--sex-counts-dir', type=Path, required=True)
    parser.add_argument('--library-summary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build(args.autosomal_reference, args.sex_reference, args.counts_dir, args.sex_counts_dir,
          args.library_summary, args.output)


if __name__ == '__main__':
    main()

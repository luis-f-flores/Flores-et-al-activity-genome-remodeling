"""Count paired fragments in autosomal 5-kb bins after sequence and genomic masks."""
from pathlib import Path
import argparse
import bisect
import csv
import subprocess
import numpy as np
import pysam

CHROMOSOMES = ['chr' + str(i) for i in range(1, 20)]
BASE = 5000


def overlap(starts, ends, start, end):
    j = bisect.bisect_left(starts, end) - 1
    return j >= 0 and ends[j] > start


def clean_sequence(seq, starts, ends):
    seq = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    acgt = (seq == 65) | (seq == 67) | (seq == 71) | (seq == 84)
    good = acgt.copy()
    for a, b in zip(starts, ends):
        good[a:b] = False
    return acgt, good


def write_table(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def prepare_reference(fasta, mask, output):
    masks = {chrom: [] for chrom in CHROMOSOMES}
    for line in mask.read_text().splitlines():
        if not line or line.startswith(('#', 'track', 'browser')):
            continue
        chrom, start, end = line.split()[:3]
        if chrom in masks:
            masks[chrom].append((int(start), int(end)))
    arrays, rows, lengths = {}, [], []
    with pysam.FastaFile(str(fasta)) as genome:
        for chrom in CHROMOSOMES:
            length = genome.get_reference_length(chrom)
            lengths.append(length)
            merged = []
            for start, end in sorted(masks[chrom]):
                if not 0 <= start < end <= length:
                    raise ValueError(f'Mask lies outside {chrom}: {start}-{end}')
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            starts = np.array([x for x, _ in merged], dtype=np.int64)
            ends = np.array([y for _, y in merged], dtype=np.int64)
            acgt, good = clean_sequence(genome.fetch(chrom), starts, ends)
            bins = np.arange(0, length, BASE)
            acgt_bp = np.add.reduceat(acgt, bins).astype(np.int64)
            clean_bp = np.add.reduceat(good, bins).astype(np.int64)
            changes = np.diff(np.concatenate(([False], ~acgt, [False])).astype(np.int8))
            arrays.update({chrom + '_nonACGT_start': np.flatnonzero(changes == 1),
                           chrom + '_nonACGT_end': np.flatnonzero(changes == -1),
                           chrom + '_clean_bp': clean_bp, chrom + '_actual_bp': np.minimum(BASE, length - bins),
                           chrom + '_ACGT_bp': acgt_bp, chrom + '_mask_start': starts, chrom + '_mask_end': ends})
            rows.append(dict(chromosome=chrom, total_bp=length, clean_bp=int(clean_bp.sum()),
                             nonACGT_bp=int((~acgt).sum()), mask_removed_ACGT_bp=int(acgt.sum() - good.sum()),
                             nonempty_clean_5kb_bins=int((clean_bp > 0).sum())))
    arrays.update(chromosomes=np.array(CHROMOSOMES), lengths=np.array(lengths))
    np.savez_compressed(output / 'clean_reference_bins.npz', **arrays)
    write_table(output / 'reference_summary.tsv', rows)


def count_sample(sample, bam, reference, output):
    subprocess.run(['samtools', 'quickcheck', '-v', str(bam)], check=True)
    ref = np.load(reference, allow_pickle=False)
    arrays = {}
    totals = dict(premask=0, fragment_mask_excluded=0, midpoint_nonACGT_excluded=0, retained=0)
    with pysam.AlignmentFile(str(bam), 'rb', threads=1) as reads:
        for chrom, length in zip(ref['chromosomes'], ref['lengths']):
            chrom, length = str(chrom), int(length)
            clean_bp = ref[chrom + '_clean_bp']
            counts = np.zeros(len(clean_bp), dtype=np.int64)
            non_acgt = np.zeros(len(clean_bp), dtype=np.int64)
            starts, ends = ref[chrom + '_mask_start'].tolist(), ref[chrom + '_mask_end'].tolist()
            ns, ne = ref[chrom + '_nonACGT_start'].tolist(), ref[chrom + '_nonACGT_end'].tolist()
            for read in reads.fetch(chrom):
                flag = read.flag
                if flag & 66 != 66 or flag & (4 | 8 | 256 | 512 | 1024 | 2048):
                    continue
                if read.next_reference_id != read.reference_id or read.template_length == 0:
                    continue
                if not read.has_tag('MQ'):
                    raise ValueError(f'{sample}: a properly paired read lacks mate mapping quality.')
                if min(read.mapping_quality, read.get_tag('MQ')) < 10:
                    continue
                totals['premask'] += 1
                size = abs(read.template_length)
                start = min(read.reference_start, read.next_reference_start)
                end = start + size
                if not 0 <= start < end <= length:
                    raise ValueError(f'{sample}: fragment lies outside {chrom}.')
                midpoint = start + size // 2
                k = midpoint // BASE
                if overlap(starts, ends, start, end):
                    totals['fragment_mask_excluded'] += 1
                    continue
                if overlap(ns, ne, midpoint, midpoint + 1):
                    non_acgt[k] += 1
                    totals['midpoint_nonACGT_excluded'] += 1
                    continue
                if clean_bp[k] <= 0:
                    raise ValueError('Counted midpoint lies in a bin with no unmasked sequence.')
                counts[k] += 1
                totals['retained'] += 1
            arrays[chrom] = counts
            arrays[chrom + '_nonACGT_midpoints'] = non_acgt
    assert totals['premask'] - totals['fragment_mask_excluded'] - totals['midpoint_nonACGT_excluded'] == totals['retained']
    np.savez_compressed(output / f'{sample}.counts.npz', **arrays)
    write_table(output / f'{sample}.counts.tsv', [dict(sample_id=sample, **totals)])
    print(f'{sample}: counted {totals["retained"]:,} fragments.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    reference = sub.add_parser('reference')
    for name in ['fasta', 'mask', 'output']:
        reference.add_argument('--' + name, type=Path, required=True)
    count = sub.add_parser('count')
    count.add_argument('--sample', required=True)
    for name in ['bam', 'reference', 'output']:
        count.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.command == 'reference':
        prepare_reference(args.fasta, args.mask, args.output)
    else:
        count_sample(args.sample, args.bam, args.reference, args.output)


if __name__ == '__main__':
    main()

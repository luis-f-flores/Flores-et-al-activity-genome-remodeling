"""Generate 50-bp RNA CPM tracks and three condition-mean tracks from BAMs."""
from pathlib import Path
import argparse
import csv
import subprocess
import numpy as np
import pysam
import pyBigWig

SAMPLES = ['Basal_A', 'Basal_B', 'Basal_C', 'Bic_A', 'Bic_B', 'Bic_C',
           'Bic_nimo_A', 'Bic_nimo_B', 'Bic_nimo_C']
GROUPS = {'TTX': SAMPLES[:3], 'Bic': SAMPLES[3:6], 'Bic_Nimo': SAMPLES[6:]}
FLAG = 2820
BIN = 50


def coverage_command(bam, output, retained, threads):
    return ['bamCoverage', '--bam', str(bam), '--outFileName', str(output),
            '--outFileFormat', 'bigwig', '--binSize', str(BIN),
            '--scaleFactor', repr(1e6 / retained), '--normalizeUsing', 'None',
            '--samFlagExclude', str(FLAG), '--numberOfProcessors', str(threads)]


def chromosome_aliases(chroms, reference):
    canonical = {str(i): 'chr' + str(i) for i in range(1, 20)}
    canonical.update(X='chrX', Y='chrY', MT='chrM')
    aliases = {}
    for name, mapped in canonical.items():
        source = name if name in chroms else mapped
        if source not in chroms or mapped not in reference or chroms[source] != reference[mapped]:
            raise ValueError('RNA/reference chromosome length mismatch: ' + mapped)
        aliases[source] = mapped
    if len({aliases.get(k, k) for k in chroms}) != len(chroms):
        raise ValueError('Chromosome aliases are not unique.')
    return aliases


def rename_track(source, output, chroms, aliases):
    temporary = output.with_name(output.name + '.partial')
    with pyBigWig.open(str(source)) as src, pyBigWig.open(str(temporary), 'w') as dst:
        if src.chroms() != chroms:
            raise ValueError('Coverage track and BAM reference differ.')
        dst.addHeader([(aliases.get(k, k), n) for k, n in chroms.items()])
        for chrom, size in chroms.items():
            for lo in range(0, size, 5000000):
                hi = min(size, lo + 5000000)
                intervals = src.intervals(chrom, lo, hi) or []
                starts, ends, values = [], [], []
                for start, end, value in intervals:
                    start, end = max(start, lo), min(end, hi)
                    if start < end:
                        starts.append(start); ends.append(end); values.append(value)
                if starts:
                    dst.addEntries([aliases.get(chrom, chrom)] * len(starts), starts,
                                   ends=ends, values=values)
    temporary.replace(output)


def condition_mean(inputs, output):
    handles = [pyBigWig.open(str(path)) for path in inputs]
    temporary = output.with_name(output.name + '.partial')
    try:
        chroms = handles[0].chroms()
        if any(handle.chroms() != chroms for handle in handles):
            raise ValueError('RNA tracks use different chromosome dictionaries.')
        with pyBigWig.open(str(temporary), 'w') as dst:
            dst.addHeader(list(chroms.items()))
            for chrom, size in chroms.items():
                for lo in range(0, size, 500000):
                    hi = min(size, lo + 500000)
                    starts = np.arange(lo, hi, BIN)
                    ends = np.minimum(starts + BIN, hi)
                    rows = []
                    for handle in handles:
                        row = np.full(len(starts), np.nan)
                        for start, end, value in handle.intervals(chrom, lo, hi) or []:
                            start, end = max(start, lo), min(end, hi)
                            row[(start-lo)//BIN:(end-lo+BIN-1)//BIN] = value
                        if not np.isfinite(row).all():
                            raise ValueError(f'RNA track has missing coverage: {chrom}:{lo}-{hi}')
                        rows.append(row)
                    values = np.mean(rows, axis=0)
                    cuts = np.r_[0, np.flatnonzero(np.diff(values) != 0) + 1, len(values)]
                    a, b = starts[cuts[:-1]].tolist(), ends[cuts[1:]-1].tolist()
                    dst.addEntries([chrom] * len(a), a, ends=b, values=values[cuts[:-1]].tolist())
        temporary.replace(output)
    finally:
        for handle in handles:
            handle.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam-dir', type=Path, required=True)
    parser.add_argument('--reference-fai', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    with args.reference_fai.open() as handle:
        reference = {row[0]: int(row[1]) for row in csv.reader(handle, delimiter='\t')}
    args.output.mkdir(parents=True, exist_ok=True)
    totals = []
    for sample in SAMPLES:
        bam = args.bam_dir / (sample + '.bam')
        subprocess.run(['samtools', 'quickcheck', '-v', str(bam)], check=True)
        with pysam.AlignmentFile(str(bam), 'rb') as handle:
            if handle.header.to_dict().get('HD', {}).get('SO') != 'coordinate' or not handle.has_index():
                raise ValueError('Coordinate-sorted, indexed BAM required: ' + str(bam))
            chroms = dict(zip(handle.references, handle.lengths))
        aliases = chromosome_aliases(chroms, reference)
        retained = int(subprocess.check_output(['samtools', 'view', '-c', '-F', str(FLAG), str(bam)], text=True))
        if retained <= 0:
            raise ValueError('No retained primary read ends: ' + sample)
        native = args.output / (sample + '.native.bw')
        output = args.output / (sample + '.CPM50.mm10.bw')
        subprocess.run(coverage_command(bam, native, retained, args.threads), check=True)
        rename_track(native, output, chroms, aliases)
        native.unlink()
        totals.append((sample, retained, 1e6 / retained))
    with (args.output / 'retained_read_counts.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample', 'retained_read_ends', 'CPM_scale'])
        writer.writerows(totals)
    for group, samples in GROUPS.items():
        condition_mean([args.output / (sample + '.CPM50.mm10.bw') for sample in samples],
                       args.output / (group + '.meanCPM50.mm10.bw'))


if __name__ == '__main__':
    main()

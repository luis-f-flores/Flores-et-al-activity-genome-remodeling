"""Write 90-kb smoothed B1/Input tracks on the supported 30-kb grid."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import pyBigWig

IPS = ['L9', 'L12', 'L11', 'L14', 'L16', 'L10', 'L13', 'L15']
INPUTS = ['L4', 'L6', 'L8', 'L2', 'L5', 'L7']
POOLS = [[0], [0], [1, 2], [1, 2], [1, 2], [3, 4, 5], [3, 4, 5], [3, 4, 5]]
CHROMS = ['chr' + str(i) for i in range(1, 20)] + ['chrX', 'chrY']
WIDTH = 30000

def smooth_three(ch, st, values):
    """Equal-weight linear mean at interior consecutive bins; edges unchanged."""
    interior=np.zeros(len(st), dtype=bool)
    interior[1:-1]=(ch[:-2]==ch[1:-1]) & (ch[1:-1]==ch[2:]) & (st[1:-1]-st[:-2]==WIDTH) & (st[2:]-st[1:-1]==WIDTH)
    out=values.copy()
    ix=np.flatnonzero(interior)
    out[ix]=(values[ix-1]+values[ix]+values[ix+1])/3
    return out, interior

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['counts', 'normalization', 'fai', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = np.load(args.counts)
    if data['sample_ids'].tolist() != INPUTS + IPS:
        raise ValueError('Unexpected Input/IP sample order')
    chrom, start, counts = data['chromosome'], data['start'], data['counts']
    norm = json.loads(args.normalization.read_text())
    ip_total = np.array(norm['autosomal_IP_totals'])
    input_total = np.array(norm['autosomal_Input_totals'])
    fly_total = np.array(norm['fly_fragments'])
    denominator = np.column_stack([counts[:, pool].sum(1) / input_total[pool].sum() for pool in POOLS])
    depth = (counts[:, 6:] / ip_total) / denominator
    fly_factor = ip_total / fly_total / norm['calibration_anchor']
    values = {'depth': depth, 'fly': depth * fly_factor}
    sizes = {row[0]: int(row[1]) for row in csv.reader(args.fai.open(), delimiter='\t')}
    rows = []
    for scale, linear in values.items():
        if not np.isfinite(linear).all() or not (linear > 0).all():
            raise ValueError('B1/Input ratios must be finite and positive')
        smoothed, _ = smooth_three(chrom, start, linear)
        signal = np.log2(smoothed)
        for i, sample in enumerate(IPS):
            path = args.output / f'{sample}_condition_{scale}_log2_smooth90kb.bw'
            with pyBigWig.open(str(path), 'w') as bw:
                bw.addHeader([(c, sizes[c]) for c in CHROMS])
                for c in CHROMS:
                    use = np.flatnonzero(chrom == c)
                    bw.addEntries([c] * len(use), start[use].tolist(),
                                  ends=np.minimum(start[use] + WIDTH, sizes[c]).tolist(),
                                  values=signal[use, i].tolist())
            rows.append(dict(sample=sample, scale=scale, path=path.name))
    with (args.output / 'tracks.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=['sample', 'scale', 'path'], delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


if __name__ == '__main__':
    main()

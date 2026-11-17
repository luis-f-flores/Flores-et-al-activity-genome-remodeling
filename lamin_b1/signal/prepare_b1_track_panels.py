"""Calculate autosomal B1 spread and the chromosome 19 track measurements."""
from pathlib import Path
import argparse
import csv
import numpy as np
import pandas as pd

GROUPS = [[0, 1], [2, 3, 4], [5, 6, 7]]


def smooth_three(chrom, start, linear):
    """Average three consecutive supported bins, retaining segment edges."""
    out = linear.copy()
    interior = ((chrom[:-2] == chrom[1:-1]) & (chrom[1:-1] == chrom[2:]) &
                (np.diff(start)[:-1] == 30000) & (np.diff(start)[1:] == 30000))
    ix = np.flatnonzero(interior) + 1
    out[ix] = (linear[ix - 1] + linear[ix] + linear[ix + 1]) / 3
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['signals', 'decoded', 'genes', 'fai', 'output']:
        p.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    v = np.load(a.signals)
    ch, st, raw = v['chromosome'], v['start'], v['condition_depth']
    auto = np.isin(ch, ['chr' + str(i) for i in range(1, 20)])
    means = np.column_stack([raw[:, group].mean(1) for group in GROUPS])
    np.savez_compressed(a.output / 'D_bin_means.npz', chromosome=ch[auto], start=st[auto],
                        means=means[auto], clean_bp=v['clean_bp'][auto])
    # bigWig stores track values as float32 before replicate averaging.
    smooth = np.log2(smooth_three(ch, st, 2 ** raw)).astype(np.float32).astype(float)
    calls = np.column_stack([np.load(a.decoded)['states'][:, group].max(1) for group in GROUPS])
    whole = ch == 'chr19'
    selected = whole & (st >= 10800000) & (st < 21600000)
    sizes = {r[0]: int(r[1]) for r in csv.reader(a.fai.open(), delimiter='\t')}
    np.savez_compressed(a.output / 'chr19_tracks.npz', starts=st[selected], sig=smooth[selected],
                        means=np.column_stack([smooth[selected][:, group].mean(1) for group in GROUPS]),
                        state=calls[selected], ws=st[whole],
                        wm=np.column_stack([smooth[whole][:, group].mean(1) for group in GROUPS]),
                        length=sizes['chr19'])
    genes = pd.read_csv(a.genes, sep='\t', header=None)
    genes = genes[(genes[0] == 'chr19') & (genes[1] < 21600000) & (genes[2] > 10800000)]
    if genes.shape[1] > 6:
        genes[3] = genes[6]
    genes.iloc[:, :6].to_csv(a.output / 'genes.bed', sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()

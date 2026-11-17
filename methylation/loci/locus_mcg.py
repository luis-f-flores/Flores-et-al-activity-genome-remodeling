"""Read CpG methylation at individual positions and average measured replicates."""
from pathlib import Path
import numpy as np
import pandas as pd
import pyBigWig


def extract(bigwigs, windows):
    samples, means = [], []
    for gene, chrom, lo, hi, time in windows:
        for group, suffix in [('0 h', '0'), ('Mannitol', f'{time}M'), ('KCl', f'{time}K')]:
            arrays = []
            for rep in 'ABC':
                sample = rep + suffix
                with pyBigWig.open(str(Path(bigwigs) / (sample + '.bw'))) as track:
                    values = np.asarray(track.values(chrom, lo, hi))
                    intervals = track.intervals(chrom, lo, hi) or []
                if not all(end - start == 1 for start, end, value in intervals):
                    raise ValueError('The mCG tracks must contain single-position intervals.')
                measured = values[np.isfinite(values)]
                if not np.all((measured >= 0) & (measured <= 100)):
                    raise ValueError('The mCG track values must be percentages.')
                arrays.append(values)
                for start, end, value in intervals:
                    samples.append(dict(gene=gene, chr=chrom, position0=start, position1=start+1,
                                        sample=sample, group=group, mCG=value))
            values = np.array(arrays)
            support = np.isfinite(values).sum(axis=0)
            averages = np.divide(np.nansum(values, axis=0), support,
                                 out=np.full(hi-lo, np.nan), where=support > 0)
            for index in np.flatnonzero(support):
                means.append(dict(gene=gene, chr=chrom, position0=lo+index, position1=lo+index+1,
                                  group=group, mCG=averages[index], n_replicates=int(support[index])))
    return pd.DataFrame(samples), pd.DataFrame(means)

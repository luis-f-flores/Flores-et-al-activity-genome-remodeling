"""Test IAP sequence-density enrichment across shared B1 states by chromosome rotations."""
from pathlib import Path
import argparse
import csv
import numpy as np

CONDITIONS = ['TTX', 'Bicuculline', 'Bicuculline + nimodipine']
GROUPS = [[0, 1], [2, 3, 4], [5, 6, 7]]
CHROMS = ['chr' + str(i) for i in range(1, 20)] + ['chrX', 'chrY']
N = 9999
SEED = 20261015


def make_map(w,st):
    cum=np.r_[0,np.cumsum(w)];k=np.r_[0,np.flatnonzero(st[1:]!=st[:-1])+1,len(w)]
    bounds=cum[k];labels=np.eye(3,dtype=np.int64)[st[k[:-1]]]
    pre=np.vstack([np.zeros(3,np.int64),np.cumsum(np.diff(bounds)[:,None]*labels,axis=0)])
    return bounds,labels,pre
def integral(x,bounds,labels,pre):
    L=bounds[-1];turn,pos=np.divmod(x,L);i=np.searchsorted(bounds[1:],pos,side='right')
    return turn[...,None]*pre[-1]+pre[i]+(pos-bounds[i])[...,None]*labels[i]
def shifted_overlap(s,e,shifts,m):
    return (integral(e[None,:]+shifts[:,None],*m)-integral(s[None,:]+shifts[:,None],*m)).sum(axis=1)

def write_table(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['signals', 'decoded', 'projection', 'output']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=True)
    signal = np.load(a.signals); states = np.load(a.decoded)['states']; projected = np.load(a.projection)
    null = np.zeros((N, len(CHROMS), 3, 3), np.int64)
    observed = np.zeros((len(CHROMS), 3, 3), np.int64); territory = np.zeros_like(observed)
    shifts = np.zeros((N, len(CHROMS)), np.int64); rng = np.random.default_rng(SEED)
    for ci, chrom in enumerate(CHROMS):
        use = np.flatnonzero(signal['chromosome'] == chrom)
        widths = projected[chrom+'_bin_width']; starts = projected[chrom+'_IAP_start']; ends = projected[chrom+'_IAP_end']
        np.testing.assert_array_equal(widths, signal['clean_bp'][use])
        shifts[:, ci] = rng.integers(0, int(widths.sum()), size=N)
        for j, group in enumerate(GROUPS):
            state = states[use][:, group].max(1); mapping = make_map(widths, state)
            observed[ci, j] = shifted_overlap(starts, ends, np.array([0]), mapping)[0]
            territory[ci, j] = mapping[2][-1]
            for offset in range(0, N, 128):
                null[offset:offset+128, ci, j] = shifted_overlap(starts, ends, shifts[offset:offset+128, ci], mapping)
            if not np.all(null[:, ci, j].sum(1) == (ends-starts).sum()):
                raise ValueError('Circular shifts did not conserve IAP sequence coverage.')
        print(chrom, 'complete', flush=True)
    np.savez_compressed(a.output/'IAP_spatial_null.npz', null_bp=null, observed_bp=observed,
                        territory_bp=territory, shifts_bp=shifts, chromosomes=np.array(CHROMS))
    distribution = [dict(condition=condition, state=k, IAP_bp=int(observed[:, j, k].sum()),
                         territory_bp=int(territory[:, j, k].sum()))
                    for j, condition in enumerate(CONDITIONS) for k in range(3)]
    write_table(a.output/'IAP_state_distribution.tsv', distribution)
    obs, terr, nu = observed.sum(0), territory.sum(0), null.sum(1)
    results = []
    for j, condition in enumerate(CONDITIONS):
        for label, numerator, denominator in [('LAD / non-LAD', [1, 2], [0]), ('High / intermediate', [2], [1])]:
            value = (obs[j, numerator].sum()/terr[j, numerator].sum())/(obs[j, denominator].sum()/terr[j, denominator].sum())
            values = (nu[:, j, numerator].sum(1)/terr[j, numerator].sum())/(nu[:, j, denominator].sum(1)/terr[j, denominator].sum())
            results.append(dict(scope='all', condition=condition, endpoint=label, density_ratio=value,
                                P_empirical=(1+np.count_nonzero(values >= value))/(N+1), null_mean=values.mean(),
                                null_q025=np.quantile(values, .025), null_q975=np.quantile(values, .975), permutations=N))
    pp = np.array([r['P_empirical'] for r in results]); order = np.argsort(pp); qq = np.empty_like(pp)
    qq[order] = np.minimum(1, np.minimum.accumulate((pp[order]*len(pp)/np.arange(1, len(pp)+1))[::-1])[::-1])
    for row, q in zip(results, qq): row['q_BH'] = q
    write_table(a.output/'IAP_enrichment_tests.tsv', results)


if __name__ == '__main__':
    main()

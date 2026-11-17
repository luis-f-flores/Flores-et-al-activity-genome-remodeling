#!/usr/bin/env python3
"""Recalculate the published continuous DamID measurements in Fig. S01G-N.

Starting point: 16 site-supported 400-kb count arrays. Full-depth
2-Mb sliding windows supply G and I-N; fixed nonoverlapping 2-Mb windows and
seeded assigned-count subsampling supply H. All calculations
are descriptive. The two biological pairs per treatment are not paired
across treatments, and overlapping windows are not independent replicates.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_INPUTS = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'damid_continuous'
TIMES = ((6, '6 h'), (24, 'Recovery'))
PAIRS = (('Mannitol', 'M', 'A'), ('Mannitol', 'M', 'B'),
         ('KCl', 'K', 'A'), ('KCl', 'K', 'B'))
SEED, DRAWS = 2026111301, 1000


def read_table(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


def window_indices(coords):
    """Every complete autosomal 2-Mb window, stepped by 400 kb."""
    regions, indices = [], []
    for chrom in (f'chr{i}' for i in range(1, 20)):
        local = coords[coords.chrom == chrom].sort_values('start')
        ix = local.index.to_numpy()
        for start in range(len(ix) - 4):
            use = ix[start:start + 5]
            starts = coords.loc[use, 'start'].to_numpy()
            ends = coords.loc[use, 'end'].to_numpy()
            if np.all(ends - starts == 400000) and np.all(np.diff(starts) == 400000):
                regions.append((chrom, int(starts[0]), int(ends[-1])))
                indices.append(use)
    if not regions:
        raise ValueError('No complete 2-Mb windows on the original 400-kb grid')
    return pd.DataFrame(regions, columns=['chrom', 'start0', 'end0']), np.array(indices)


def continuous_windows(coords, counts, libraries):
    base, indices = window_indices(coords)
    ranks, support = [], []
    for hours, time in TIMES:
        enrichments, eligibility, fields = [], [], []
        for condition, code, rep in PAIRS:
            fusion, dam = f'{hours}{rep}-DB{code}', f'{hours}{rep}-D{code}'
            fc = counts[fusion + '_count'][indices].sum(axis=1)
            dc = counts[dam + '_count'][indices].sum(axis=1)
            fp = counts[fusion + '_positive'][indices].sum(axis=1)
            dp = counts[dam + '_positive'][indices].sum(axis=1)
            signal = np.log2((fc / libraries.loc[fusion, 'total_assigned_counts'] * 1e6 + 1) /
                             (dc / libraries.loc[dam, 'total_assigned_counts'] * 1e6 + 1))
            eligible = (fc >= 10) & (dc >= 10) & (fp >= 10) & (dp >= 10) & (fc + dc >= 25)
            enrichments.append(signal)
            eligibility.append(eligible)
            row = base.copy()
            row['time'], row['condition'], row['replicate'] = time, condition, rep
            row['fusion_count'], row['Dam_count'] = fc, dc
            row['fusion_sites'], row['Dam_sites'] = fp, dp
            row['log2_enrichment'], row['pair_eligible'] = signal, eligible
            fields.append(row)
        signals = np.array(enrichments).T
        eligible = np.all(eligibility, axis=0)
        row = base.copy()
        row['time'], row['eligible'] = time, eligible
        for i, (condition, _, rep) in enumerate(PAIRS):
            row[condition + '_' + rep] = signals[:, i]
        row['mean_change'] = signals[:, 2:].mean(axis=1) - signals[:, :2].mean(axis=1)
        row['minimum_reads'] = np.min([f[['fusion_count', 'Dam_count']].min(axis=1) for f in fields], axis=0)
        row['minimum_sites'] = np.min([f[['fusion_sites', 'Dam_sites']].min(axis=1) for f in fields], axis=0)
        row['eligible_n'] = int(eligible.sum())
        row['percentile'], row['rank_from_gain'] = np.nan, np.nan
        row.loc[eligible, 'percentile'] = 100 * (row.loc[eligible, 'mean_change'].rank(method='average') - .5) / eligible.sum()
        row.loc[eligible, 'rank_from_gain'] = row.loc[eligible, 'mean_change'].rank(ascending=False, method='min')
        ranks.append(row)
        support.extend(fields)
    return pd.concat(ranks, ignore_index=True), pd.concat(support, ignore_index=True)


def depth_matched_sd(windows, support, libraries):
    """Calculate depth-matched signal spread on fixed 2-Mb tiles.

    The remainder category retains all assigned counts outside complete
    autosomal tiles, including sex chromosomes and chromosome-end remainders.
    Sampling order, unchanged lowest-depth libraries, and NumPy's marginals
    algorithm are preserved. The published SD used a float32 saved-signal
    round trip before float64 SD calculation; that precision is reproduced.
    """
    rng = np.random.default_rng(SEED)
    points, ratios, depth_support, draws_by_time = [], [], [], {}
    fixed = windows[windows.start0.mod(2000000).eq(0)].copy()
    for hours, time in TIMES:
        base = fixed[fixed.time == time].reset_index(drop=True)
        mask = base.eligible.to_numpy(dtype=bool)
        targets = {}
        for role, prefix in (('Fusion', 'DB'), ('Dam-only', 'D')):
            names = [f'{hours}{rep}-{prefix}{code}' for code in ('M', 'K') for rep in ('A', 'B')]
            targets[role] = int(libraries.loc[names, 'total_assigned_counts'].min())
        per_pair_sd = []
        for condition, code, rep in PAIRS:
            sub = support[(support.time == time) & (support.condition == condition) &
                          (support.replicate == rep) & support.start0.mod(2000000).eq(0)].reset_index(drop=True)
            if not sub[['chrom', 'start0', 'end0']].equals(base[['chrom', 'start0', 'end0']]):
                raise ValueError('Fixed-window count grids differ')
            pair = []
            for role, prefix, field in (('Fusion', 'DB', 'fusion_count'), ('Dam-only', 'D', 'Dam_count')):
                name = f'{hours}{rep}-{prefix}{code}'
                total = int(libraries.loc[name, 'total_assigned_counts'])
                target = targets[role]
                observed = sub[field].to_numpy(dtype=np.int64)
                partition = np.r_[observed, total - int(observed.sum())]
                if partition.min() < 0 or partition.sum() != total:
                    raise ValueError('Invalid whole-library subsampling partition')
                if target == total:
                    draws = np.broadcast_to(partition, (DRAWS, len(partition))).copy()
                else:
                    draws = rng.multivariate_hypergeometric(partition, target, size=DRAWS, method='marginals')
                if (draws < 0).any() or (draws > partition).any() or not np.all(draws.sum(axis=1) == target):
                    raise ValueError('Subsample count bounds or depth target failed')
                pair.append(draws[:, :-1] / target * 1e6)
                depth_support.append(dict(time=time, condition=condition, replicate=rep, library=name,
                                          role=role, original_assigned_counts=total,
                                          matched_assigned_counts=target, fixed_eligible_windows=int(mask.sum())))
            signal = np.log2((pair[0] + 1) / (pair[1] + 1)).astype(np.float32).astype(np.float64)
            sd = signal[:, mask].std(axis=1, ddof=1)
            per_pair_sd.append(sd)
            points.append(dict(time=time, condition=condition, replicate=rep, value=float(np.median(sd))))
        draws_by_time[f'{hours}_SD_draws'] = np.array(per_pair_sd)
        values = np.median(per_pair_sd, axis=1)
        ratios.append(dict(time=time, ratio_of_plotted_means=float(values[2:].mean() / values[:2].mean())))
    return fixed, pd.DataFrame(points), pd.DataFrame(ratios), pd.DataFrame(depth_support), draws_by_time


def calculate(inputs=DEFAULT_INPUTS, annotation=None):
    inputs = Path(inputs)
    coords = read_table(inputs / 'coordinates_400kb.tsv')
    libraries = read_table(inputs / 'library_support.tsv').set_index('sample')
    with np.load(inputs / 'site_support_400kb.npz', allow_pickle=False) as archive:
        counts = {name: archive[name] for name in archive.files}
    required = [f'{h}{rep}-{prefix}{code}' for h, _ in TIMES for _, code, rep in PAIRS for prefix in ('DB', 'D')]
    if not libraries.index.is_unique or set(libraries.index) != set(required):
        raise ValueError('Expected exactly the 16 publication libraries')
    if coords.duplicated(['chrom', 'start', 'end']).any():
        raise ValueError('Duplicate 400-kb coordinates')
    for name in required:
        count, positive = counts[name + '_count'], counts[name + '_positive']
        if count.shape != (len(coords),) or positive.shape != count.shape:
            raise ValueError(f'Count/site grid mismatch: {name}')
        if (count < 0).any() or (positive < 0).any() or (positive > count).any():
            raise ValueError(f'Invalid counts or positive-feature support: {name}')
        if int(count.sum()) != int(libraries.loc[name, 'total_assigned_counts']) or count.sum() <= 0:
            raise ValueError(f'Whole-library assigned-count total mismatch: {name}')
    windows, support = continuous_windows(coords, counts, libraries)
    fixed, points, ratios, depth, draws = depth_matched_sd(windows, support, libraries)
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from select_rank_examples import select
    if annotation is None:
        annotation = Path(os.environ.get('GB_DATA_ROOT', '../inputs')) / 'references/genes.bed'
    callouts, genes = select(windows, annotation)
    summary = dict(seed=SEED, subsamples=DRAWS, numpy_version=np.__version__,
                   eligible_sliding_windows={time: int(windows[(windows.time == time) & windows.eligible].shape[0]) for _, time in TIMES},
                   eligible_fixed_windows={time: int(fixed[(fixed.time == time) & fixed.eligible].shape[0]) for _, time in TIMES},
                   SD_points=points.to_dict('records'), SD_ratios=ratios.to_dict('records'))
    tables = {'windows_2Mb.tsv.gz': windows, 'window_read_support.tsv.gz': support,
              'fixed_2Mb_windows.tsv.gz': fixed, 'SD_points.tsv': points, 'SD_ratios.tsv': ratios,
              'depth_matching.tsv': depth, 'rank_callouts.tsv': callouts, 'gene_windows.tsv': genes}
    return summary, tables


def write_results(output, summary, tables):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(output / name, sep='\t', index=False, float_format='%.15g')
    (output / 'RESULTS.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=DEFAULT_INPUTS)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--annotation', type=Path)
    args = parser.parse_args()
    summary, tables = calculate(args.inputs, args.annotation)
    write_results(args.output, summary, tables)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

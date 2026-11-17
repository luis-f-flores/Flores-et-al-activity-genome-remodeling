"""Measure activity-response gene coverage and B1 signal across shared LAD calls."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

GROUPS = [[0, 1], [2, 3, 4], [5, 6, 7]]
NAMES = ['TTX', 'Bicuculline', 'Bicuculline + nimodipine']
SETS = ['rPRG', 'dPRG', 'SRG']
ALIASES = {'Fam46a': 'Tent5a', 'Gltscr2': 'Nop53'}


def merge(rows):
    out = []
    for a, b in sorted(rows):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(b, out[-1][1]))
        else:
            out.append((int(a), int(b)))
    return out


def interval_lookup(rows):
    result = {}
    for chrom, g in rows.groupby('chr', sort=False):
        intervals = merge(zip(g.s, g.e))
        result[chrom] = (np.array([a for a, b in intervals]), np.array([b for a, b in intervals]))
    return result


def overlap(chrom, left, right, lookup):
    if chrom not in lookup:
        return 0
    starts, ends = lookup[chrom]
    i = np.searchsorted(ends, left, side='right')
    j = np.searchsorted(starts, right, side='left')
    return int(np.maximum(0, np.minimum(ends[i:j], right) - np.maximum(starts[i:j], left)).sum())


def memberships(genes, lists):
    rows, indices = [], {}
    for name in SETS:
        symbols = [s.strip() for s in (lists / (name + '_gene_list.txt')).read_text().splitlines() if s.strip()]
        resolved = [ALIASES.get(s, s) for s in symbols]
        if len(resolved) != len(set(resolved)):
            raise ValueError('Gene alias collision: ' + name)
        indices[name] = np.flatnonzero(genes.gene.isin(resolved))
        found = set(genes.iloc[indices[name]].gene)
        rows.extend(dict(gene_set=name, original_symbol=s, coordinate_symbol=t,
                         coordinate_matched=t in found, alias_applied=s != t) for s, t in zip(symbols, resolved))
    return pd.DataFrame(rows), indices


def permutations(flags, indices, names):
    rows, nulls = [], {}
    for k, name in enumerate(SETS):
        ix = indices[name]
        if not len(ix):
            raise ValueError('No coordinate-matched genes: ' + name)
        rng = np.random.default_rng(20261012 + k)
        null = np.array([flags[rng.choice(len(flags), len(ix), replace=False)].sum(0) for _ in range(10000)])
        nulls[name] = null
        observed = flags[ix].sum(0)
        for j, condition in enumerate(names):
            values = 100 * null[:, j] / len(ix)
            x = 100 * observed[j] / len(ix)
            extreme = int(np.count_nonzero(abs(values - values.mean()) >= abs(x - values.mean()) - 1e-12))
            rows.append(dict(gene_set=name, condition=condition, genes=len(ix), LAD_genes=int(observed[j]),
                             observed_percent=x, null_mean=values.mean(), null_sd=values.std(),
                             null_q025=np.quantile(values, .025), null_q975=np.quantile(values, .975),
                             P_unadjusted=(extreme + 1) / 10001, draws=10000, seed=20261012 + k))
    return pd.DataFrame(rows), nulls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['signals', 'decoded', 'coordinates', 'gene-lists', 'autosomal-reference', 'sex-reference', 'rna', 'output']:
        parser.add_argument('--' + key, type=Path, required=True)
    a = parser.parse_args(); a.output.mkdir(parents=True, exist_ok=True)
    v = np.load(a.signals); states = np.load(a.decoded)['states']
    g = pd.read_csv(a.coordinates, sep='\t', names=['chr', 's', 'e', 'gene'])
    if not g.gene.is_unique or not (g.e > g.s).all():
        raise ValueError('Gene coordinates must be unique, positive-length intervals.')
    membership, indices = memberships(g, a.gene_lists)
    membership.to_csv(a.output / 'activity_gene_membership.tsv', sep='\t', index=False)
    state = np.column_stack([states[:, group].max(1) for group in GROUPS])
    bins = pd.DataFrame(dict(chr=v['chromosome'], s=v['start'], e=v['start'] + 30000))
    support_lookup = interval_lookup(bins)
    fractions, high = np.zeros((len(g), 3)), np.zeros((len(g), 3))
    support = np.array([overlap(r.chr, r.s, r.e, support_lookup) / (r.e - r.s) for r in g.itertuples()])
    for j, name in enumerate(['TTX', 'Bic', 'Bic_Nimo']):
        lad = interval_lookup(bins[state[:, j] >= 1]); high_lad = interval_lookup(bins[state[:, j] == 2])
        for i, r in enumerate(g.itertuples()):
            fractions[i, j] = overlap(r.chr, r.s, r.e, lad) / (r.e - r.s)
            high[i, j] = overlap(r.chr, r.s, r.e, high_lad) / (r.e - r.s)
        g[name + '_LAD_fraction'] = fractions[:, j]; g[name + '_high_fraction'] = high[:, j]
        g[name + '_LAD75'] = fractions[:, j] >= .75
    g['scored_span_fraction'] = support
    g.to_csv(a.output / 'activity_background_coverage.tsv', sep='\t', index=False)
    tested, nulls = permutations(fractions >= .75, indices, NAMES)
    tested.to_csv(a.output / 'activity_LAD_permutation.tsv', sep='\t', index=False)
    np.savez_compressed(a.output / 'activity_permutation_nulls.npz', **nulls)
    selected = pd.concat([g.iloc[indices[name]].assign(gene_set=name) for name in SETS])
    selected.to_csv(a.output / 'activity_LAD_coverage.tsv', sep='\t', index=False)
    masks = dict(np.load(a.autosomal_reference)); masks.update(dict(np.load(a.sex_reference)))
    bad, index = {}, {}
    for chrom in np.unique(v['chromosome']):
        bad[chrom] = interval_lookup(pd.DataFrame([(chrom, x, y) for key in ['mask', 'nonACGT']
            for x, y in zip(masks[chrom + '_' + key + '_start'], masks[chrom + '_' + key + '_end'])], columns=['chr', 's', 'e']))
        index[chrom] = np.flatnonzero(v['chromosome'] == chrom)
    rna = pd.read_csv(a.rna, sep='\t')
    rows = []
    for name in SETS:
        for i in indices[name]:
            r = g.iloc[i]; ii = index.get(r.chr, np.array([], dtype=int))
            use, weights = [], []
            for j in ii[(v['start'][ii] < r.e) & (v['start'][ii] + 30000 > r.s)]:
                left, right = max(r.s, int(v['start'][j])), min(r.e, int(v['start'][j]) + 30000)
                weight = right - left - overlap(r.chr, left, right, bad[r.chr])
                if weight > 0:
                    use.append(j); weights.append(weight)
            values = np.average(v['condition_depth'][use], axis=0, weights=weights) if use else np.full(8, np.nan)
            means = [values[group].mean() for group in GROUPS]
            row = dict(gene_set=name, gene_name=r.gene, chromosome=r.chr, start0=r.s, end0=r.e,
                       gene_length_bp=r.e-r.s, scored_span_fraction=support[i], assessed_clean_bp=sum(weights))
            for sid, value in zip(v['IP_ids'], values): row['B1_' + str(sid)] = value
            for j, label in enumerate(['TTX', 'Bic', 'Bic_Nimo']): row['B1_' + label + '_mean'] = means[j]
            for contrast, label, ref, target in [('Bic_vs_TTX', 'BicvsTTX', 0, 1),
                    ('Bic_Nimo_vs_TTX', 'Bic_NimovsTTX', 0, 2), ('Bic_Nimo_vs_Bic', 'Bic_NimovsBic', 1, 2)]:
                row[label + '_B1_delta'] = means[target] - means[ref]
                matches = rna[(rna.contrast == contrast) & (rna.gene_name == r.gene)]
                if len(matches) == 1:
                    for key in ['log2FoldChange', 'pvalue', 'padj']: row[label + '_' + key] = matches.iloc[0][key]
            rows.append(row)
    measured = pd.DataFrame(rows)
    measured.to_csv(a.output / 'activity_B1_RNA.tsv', sep='\t', index=False)
    summaries = []
    for name in SETS:
        d = measured[measured.gene_set == name]
        means = d[['B1_' + str(s) for s in v['IP_ids']]].mean().to_numpy()
        baseline = means[:2].mean()
        summaries.extend(dict(gene_set=name, IP=str(s), genes=len(d), mean_body_log2=means[j],
                              delta_from_TTX_mean=means[j]-baseline) for j, s in enumerate(v['IP_ids']))
    pd.DataFrame(summaries).to_csv(a.output / 'activity_B1_class_replicates.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()

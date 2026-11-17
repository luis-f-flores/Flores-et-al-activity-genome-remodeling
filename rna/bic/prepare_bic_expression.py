"""Calculate Bic RNA PCA, activity scores and DEG display tables."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

SAMPLES = ['Basal_A', 'Basal_B', 'Basal_C', 'Bic_A', 'Bic_B', 'Bic_C',
           'Bic_nimo_A', 'Bic_nimo_B', 'Bic_nimo_C']
CONTRASTS = {'Bic_vs_TTX': 'BicvsTTX', 'Bic_Nimo_vs_Bic': 'Bic_NimovsBic',
             'Bic_Nimo_vs_TTX': 'Bic_NimovsTTX'}
CLASSES = {'rapidPRGs': 'rPRG', 'delayedPRGs': 'dPRG', 'SRGs': 'SRG'}
ALIASES = {'Fam46a': 'Tent5a', 'Gltscr2': 'Nop53',
           '5430416O09Rik': 'Spaar', '9430020K01Rik': 'Jcad'}


def table(path):
    frame = pd.read_csv(path, sep='\t', low_memory=False)
    if not frame.gene_id.is_unique:
        raise ValueError('Duplicate gene IDs: ' + str(path))
    return frame.set_index('gene_id')


def write(frame, path):
    frame.to_csv(path, sep='\t', index=False)


def expression(inputs, gene_lists, output):
    counts = table(inputs / 'gene_count.xls')
    raw = counts[SAMPLES].to_numpy(float)
    if not (np.isfinite(raw).all() and (raw >= 0).all() and (raw == np.floor(raw)).all()):
        raise ValueError('The gene count matrix must contain nonnegative integer counts.')
    keep = (raw >= 10).sum(axis=1) >= 3
    reference = raw[keep & (raw > 0).all(axis=1)]
    geometric = np.exp(np.log(reference).mean(axis=1))
    size = np.median(reference / geometric[:, None], axis=0)
    size /= np.exp(np.log(size).mean())
    write(pd.DataFrame({'sample': SAMPLES, 'size_factor': size}), output / 'count_size_factors.tsv')
    logged = np.log2(raw / size + 1)
    filtered = logged[keep]
    top = np.argsort(np.var(filtered, axis=1, ddof=1))[-min(500, len(filtered)):]
    matrix = filtered[top].T
    matrix -= matrix.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    scores = u[:, :3] * singular[:3]
    variance = 100 * singular**2 / np.sum(singular**2)
    pca = pd.DataFrame(scores, columns=['pc1', 'pc2', 'pc3'])
    pca['sample'] = SAMPLES
    pca['group'] = ['TTX'] * 3 + ['Bic'] * 3 + ['Bic+Nimo'] * 3
    write(pca, output / 'RNA_PCA.tsv')
    (output / 'RNA_PCA_variance.json').write_text(json.dumps(
        {f'PC{i+1}': float(variance[i]) for i in range(3)}, indent=2) + '\n')

    classifications, tables, totals = [], {}, {}
    for contrast, stem in CONTRASTS.items():
        d = table(inputs / (stem + '_deg.xls'))
        tables[stem] = d
        fc, q = d.log2FoldChange, d.padj
        valid = np.isfinite(fc) & np.isfinite(q) & (q >= 0) & (q <= 1)
        up, down = valid & (fc >= .37) & (q <= .01), valid & (fc <= -.37) & (q <= .01)
        z = d[['gene_name', 'log2FoldChange', 'pvalue', 'padj']].reset_index()
        z.insert(0, 'contrast', contrast)
        z['classification'] = np.select([up, down, valid], ['up', 'down', 'not_DEG'], default='unavailable')
        classifications.append(z)
        totals[stem] = dict(DEGs=int((up | down).sum()), up=int(up.sum()), down=int(down.sum()))
    write(pd.concat(classifications, ignore_index=True), output / 'RNA_DEG_classification.tsv')
    (output / 'RNA_DEG_counts.json').write_text(json.dumps(totals, indent=2) + '\n')

    # Subsequent tables read the exported size factors with the original TSV precision.
    size = pd.read_csv(output / 'count_size_factors.tsv', sep='\t').size_factor.to_numpy()
    log = pd.DataFrame(np.log2(raw / size + 1), index=counts.index, columns=SAMPLES)
    common = set(counts.index[keep])
    for d in tables.values():
        common &= set(d.index[np.isfinite(d.log2FoldChange) & np.isfinite(d.padj)])
    ids = sorted(common)
    diagnostic = pd.DataFrame(index=pd.Index(ids, name='gene_id'))
    diagnostic['gene_name'] = counts.loc[ids, 'gene_name']
    for sample in SAMPLES:
        diagnostic[sample + '_log2_normalized_count_plus1'] = log.loc[ids, sample]
    for name, samples in [('TTX', SAMPLES[:3]), ('Bic', SAMPLES[3:6]), ('Bic+Nimo', SAMPLES[6:])]:
        diagnostic[name + '_SD'] = log.loc[ids, samples].std(axis=1, ddof=1)
    for stem, d in tables.items():
        d = d.loc[ids]
        for column in ['log2FoldChange', 'pvalue', 'padj']:
            diagnostic[stem + '_' + column] = d[column]
        diagnostic[stem + '_DEG'] = (d.log2FoldChange.abs() >= .37) & (d.padj <= .01)
    write(diagnostic.reset_index(), output / 'RNA_diagnostic_genes.tsv')

    lookup = counts.assign(symbol_key=counts.gene_name.astype(str).str.lower()).reset_index()
    activity, summary = [], []
    for name, short in CLASSES.items():
        symbols = [s.strip() for s in (gene_lists / (short + '_gene_list.txt')).read_text().splitlines() if s.strip()]
        values = []
        for symbol in symbols:
            resolved = ALIASES.get(symbol, symbol)
            match = lookup[lookup.symbol_key == resolved.lower()]
            if len(match) != 1:
                raise ValueError(f'Expected one count-table gene for {symbol}: {len(match)}')
            gid = match.iloc[0].gene_id
            v = log.loc[gid].to_numpy(float)
            sd = v.std(ddof=1)
            if not np.isfinite(v).all() or sd <= 0:
                raise ValueError('Activity gene has no finite expression spread: ' + symbol)
            z = (v - v.mean()) / sd
            values.append(z)
            activity.append(dict(original_symbol=symbol, gene_name=resolved, gene_id=gid,
                                 gene_set=name, **dict(zip(SAMPLES, z))))
        for j, sample in enumerate(SAMPLES):
            summary.append(dict(gene_set=name, sample=sample, n_genes=len(values),
                                mean_gene_z_score=float(np.asarray(values)[:, j].mean())))
    genes = pd.DataFrame(activity)
    if not genes.gene_id.is_unique:
        raise ValueError('Activity gene classes contain repeated gene IDs.')
    write(genes, output / 'activity_gene_z_scores.tsv')
    write(pd.DataFrame(summary), output / 'activity_scores_173.tsv')


def heatmap(inputs, expression_dir, annotation, output):
    classifications = pd.read_csv(expression_dir / 'RNA_DEG_classification.tsv', sep='\t')
    deg = classifications[(classifications.contrast == 'Bic_vs_TTX') & classifications.classification.isin(['up', 'down'])]
    deg = pd.concat([deg[deg.classification == direction].sort_values('gene_id')
                     for direction in ['up', 'down']], ignore_index=True)
    counts = table(inputs / 'gene_count.xls').loc[deg.gene_id, SAMPLES].to_numpy(float)
    size = pd.read_csv(expression_dir / 'count_size_factors.tsv', sep='\t').set_index('sample').loc[SAMPLES, 'size_factor'].to_numpy()
    normalized = counts / size
    logged = np.log2(normalized + 1)
    z = (logged - logged.mean(axis=1, keepdims=True)) / logged.std(axis=1, ddof=1)[:, None]
    order, trees = [], []
    for direction in ['up', 'down']:
        ix = np.flatnonzero(deg.classification.to_numpy() == direction)
        distances = np.clip(pdist(z[ix], metric='correlation'), 0, 2)
        tree = linkage(distances, method='average', optimal_ordering=True)
        leaves = leaves_list(tree)
        order.extend(ix[leaves])
        trees.append(dict(direction=direction, input_gene_ids=deg.gene_id.iloc[ix].tolist(),
                          linkage=tree.tolist(), leaf_order=leaves.tolist(),
                          ordered_gene_ids=deg.gene_id.iloc[ix[leaves]].tolist()))
    order = np.asarray(order)
    deg = deg.iloc[order].reset_index(drop=True)
    counts, normalized, z = counts[order], normalized[order], z[order]
    member = pd.read_csv(annotation, sep='\t')
    member = member[member.contrast == 'Bic_vs_TTX'].set_index('gene_id').reindex(deg.gene_id)
    in_lad = member.baseline_LAD.astype(str).str.lower() == 'true'
    deg['TTX_LAD_display'] = np.where(in_lad, 'LAD', np.where(member.annotation_status == 'unique', 'outside', 'unmapped'))
    effects = classifications[classifications.contrast == 'Bic_Nimo_vs_Bic'].set_index('gene_id').reindex(deg.gene_id)
    deg['Bic_Nimo_vs_Bic_log2FC'] = effects.log2FoldChange.to_numpy()
    deg['Bic_Nimo_vs_Bic_padj'] = effects.padj.to_numpy()
    deg['heatmap_row'] = np.arange(len(deg))
    for j, sample in enumerate(SAMPLES):
        deg[sample + '_raw_count'] = counts[:, j]
        deg[sample + '_norm_count'] = normalized[:, j]
        deg[sample + '_z'] = z[:, j]
    write(deg, output / 'RNA_DEG_heatmap.tsv')
    (output / 'RNA_heatmap_linkage.json').write_text(json.dumps(trees, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['expression', 'heatmap'])
    for name in ['inputs', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--gene-lists', type=Path)
    parser.add_argument('--expression', type=Path)
    parser.add_argument('--annotation', type=Path)
    a = parser.parse_args()
    if a.stage == 'expression' and a.gene_lists is None:
        parser.error('expression requires --gene-lists')
    if a.stage == 'heatmap' and (a.expression is None or a.annotation is None):
        parser.error('heatmap requires --expression and --annotation')
    a.output.mkdir(parents=True, exist_ok=True)
    if a.stage == 'expression':
        expression(a.inputs, a.gene_lists, a.output)
    else:
        heatmap(a.inputs, a.expression, a.annotation, a.output)


if __name__ == '__main__':
    main()

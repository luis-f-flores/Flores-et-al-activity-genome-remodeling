"""Prepare B1 panel tables from computed signals, LAD summaries and gene analyses."""
from pathlib import Path
import argparse
import json
import shutil
import numpy as np
import pandas as pd
from scipy.stats import rankdata, gaussian_kde

GROUPS = [[0, 1], [2, 3, 4], [5, 6, 7]]
NAMES = ['TTX', 'Bicuculline', 'Bicuculline + nimodipine']
POOLS = [[0], [0], [1, 2], [1, 2], [1, 2], [3, 4, 5], [3, 4, 5], [3, 4, 5]]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['signals', 'counts', 'decoded', 'summary', 'profiles', 'output']:
        p.add_argument('--' + key, type=Path, required=True)
    for key in ['pairwise', 'directional', 'b1-go', 'rna-go', 'context', 'activity', 'baseline-tpm']:
        p.add_argument('--' + key, type=Path)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=True)
    v = np.load(a.signals); counts = np.load(a.counts)['counts']; states = np.load(a.decoded)['states']
    ids, bp = v['IP_ids'], v['clean_bp']
    x = v['condition_depth']; matrix = x.T - x.mean(1)[None, :]
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    scores = u[:, :2] * s[:2]
    for k in range(2):
        if scores[np.argmax(abs(scores[:, k])), k] < 0:
            scores[:, k] *= -1
    pd.DataFrame(dict(scheme='condition_depth', scope='all', sample=ids, PC1=scores[:, 0], PC2=scores[:, 1])).to_csv(
        a.output / 'bin_PCA_scores.tsv', sep='\t', index=False)
    pd.DataFrame(dict(scheme='condition_depth', scope='all', component=np.arange(1, len(s)+1),
                      percent=100*s*s/sum(s*s))).to_csv(a.output / 'bin_PCA_variance.tsv', sep='\t', index=False)
    # Ranking the integer-count ratios avoids floating-point tie changes from scaling.
    denominator = np.column_stack([counts[:, pool].sum(1) for pool in POOLS])
    ranks = rankdata(counts[:, 6:] / denominator, axis=0)
    corr = np.corrcoef(ranks.T)
    pd.DataFrame([dict(scheme='condition_depth', scope='all', sample1=ids[i], sample2=ids[j], Spearman=corr[i, j])
                  for i in range(8) for j in range(8)]).to_csv(a.output / 'bin_correlations.tsv', sep='\t', index=False)
    np.savez_compressed(a.output / 'LAD_replication.npz', states=states, clean_bp=bp, IP_ids=ids)
    call = np.column_stack([states[:, group].max(1) for group in GROUPS])
    rows = [dict(condition=name, state=label, clean_bp=int(bp[call[:, j] == k].sum()),
                 Mb=bp[call[:, j] == k].sum()/1e6, percent_scored_genome=100*bp[call[:, j] == k].sum()/bp.sum())
            for j, name in enumerate(NAMES) for k, label in enumerate(['Non-B1-associated', 'Intermediate', 'High'])]
    pd.DataFrame(rows).to_csv(a.output / 'condition_state_extent.tsv', sep='\t', index=False)
    for name in ['domain_PCA_scores.tsv', 'domain_PCA_variance.tsv']:
        shutil.copyfile(a.summary / name, a.output / name)
    shutil.copyfile(a.profiles / 'profile_arrays.npz', a.output / 'LAD_profiles.npz')
    d = pd.read_csv(a.profiles / 'inside_outside.tsv', sep='\t')
    d = d[d.reference == 'Union'].copy()
    d['condition'] = d.condition.replace({'Bic': NAMES[1], 'Bic+Nimo': NAMES[2]})
    d = d.rename(columns={'inside_log2': 'LAD_mean_log2', 'outside_log2': 'nonLAD_mean_log2',
                          'contrast_ratio': 'inside_outside_ratio'})
    d[['scheme', 'sample', 'condition', 'LAD_mean_log2', 'nonLAD_mean_log2', 'inside_outside_ratio']].to_csv(
        a.output / 'calibration.tsv', sep='\t', index=False)
    for source, name in [(a.pairwise, 'pairwise_B1_tests.tsv'), (a.directional, 'B1_RNA_OR.tsv'),
                         (a.context, 'LAD_context.tsv')]:
        if source:
            shutil.copyfile(source, a.output / name)
    for source, name in [(a.b1_go, 'B1_GO_CC.tsv'), (a.rna_go, 'RNA_GO_CC.tsv')]:
        if source:
            d = pd.read_csv(source, sep='\t')
            if 'ontology' in d:
                d = d[d.ontology == 'CC']
            d.to_csv(a.output / name, sep='\t', index=False)
    if a.activity:
        for name in ['activity_B1_RNA.tsv', 'activity_B1_class_replicates.tsv',
                     'activity_LAD_coverage.tsv', 'activity_LAD_permutation.tsv']:
            shutil.copyfile(a.activity / name, a.output / name)
    if a.baseline_tpm:
        d = pd.read_csv(a.baseline_tpm, sep='\t', float_precision='round_trip')
        rows = []
        for state in range(3):
            values = d.loc[d.state == state, 'log2_TTX_TPM_plus1'].to_numpy()
            y = np.linspace(values.min(), values.max(), 100)
            density = gaussian_kde(values)(y); width = density/density.max()*.25
            q1, median, q3 = np.percentile(values, [25, 50, 75]); iqr = q3-q1
            rows.append(dict(state=state, n=len(values), y=y.tolist(), width=width.tolist(), q25=q1,
                             median=median, q75=q3, whislo=float(min(values[values >= q1-1.5*iqr])),
                             whishi=float(max(values[values <= q3+1.5*iqr]))))
        (a.output / 'TTX_TPM_state_violin.json').write_text(json.dumps(rows, indent=2)+'\n')


if __name__ == '__main__':
    main()

"""Recompute the published IAP-bin and baseline-matched non-IAP B1 summaries."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

TARGETS = ['IAPEY4_I', 'IAPA_MM', 'IAPEZI']


def calculate(inputs, output):
    inputs, output = Path(inputs), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    a = dict(np.load(inputs / 'bins.npz'))
    membership = np.load(inputs / 'memberships.npz')
    baseline = a['B1'][:, :2].mean(axis=1)
    caliper = 0.1 * baseline.std(ddof=1)
    weights = a['unmasked_bp'].astype(float)
    matches, match_rows = {}, []
    focal = sorted(set(np.concatenate([membership[k] for k in TARGETS]).tolist()))
    for chrom in sorted(set(a['chromosome'])):
        for state in [0, 1, 2]:
            eligible = np.flatnonzero((a['chromosome'] == chrom) &
                                      (a['TTX_state'] == state) & ~a['catalog_IAP'])
            candidates = eligible[np.lexsort((eligible, baseline[eligible]))]
            values = baseline[candidates]
            for target in (i for i in focal if a['chromosome'][i] == chrom and a['TTX_state'][i] == state):
                position = np.searchsorted(values, baseline[target])
                nearby = candidates[max(0, position - 5):min(len(values), position + 5)]
                ranked = sorted(nearby, key=lambda i: (abs(baseline[i] - baseline[target]), int(i)))
                selected = np.array([i for i in ranked if abs(baseline[i] - baseline[target]) <= caliper][:5], int)
                matches[target] = selected
                for control in selected:
                    match_rows.append({'target_bin': target, 'control_bin': int(control)})
    observed_matches = pd.DataFrame(match_rows).sort_values(['target_bin', 'control_bin']).reset_index(drop=True)
    rows = []
    for target in TARGETS:
        ids = membership[target]
        matched = np.array([i for i in ids if len(matches[i])], int)
        original = np.average(a['B1'][ids], weights=weights[ids], axis=0)
        target_mean = np.average(a['B1'][matched], weights=weights[matched], axis=0)
        control_mean = np.average(np.array([a['B1'][matches[i]].mean(axis=0) for i in matched]),
                                  weights=weights[matched], axis=0)
        for j, sample in enumerate(a['IP_ids']):
            rows.append({'target': target, 'sample': sample, 'bins': len(ids),
                         'matched_bins': len(matched), 'target_raw': original[j],
                         'target_matched_subset': target_mean[j], 'matched_control': control_mean[j],
                         'target_minus_matched': target_mean[j] - control_mean[j]})
    result = pd.DataFrame(rows)
    result = result.set_index(['target', 'sample'])
    result.reset_index().to_csv(output / 'IAP_B1_per_IP.tsv', sep='\t', index=False)
    observed_matches.to_csv(output / 'matched_bins.tsv.gz', sep='\t', index=False)
    return result.reset_index()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(calculate(args.inputs, args.output).to_string(index=False))

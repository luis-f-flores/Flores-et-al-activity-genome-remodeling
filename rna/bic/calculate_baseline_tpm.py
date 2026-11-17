"""Calculate Bic-experiment TPM and baseline RNA abundance by B1 state."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

GROUPS = {'TTX': ['Basal_A', 'Basal_B', 'Basal_C'],
          'Bic': ['Bic_A', 'Bic_B', 'Bic_C'],
          'Bic_Nimo': ['Bic_nimo_A', 'Bic_nimo_B', 'Bic_nimo_C']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['counts', 'states', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    counts = pd.read_csv(args.counts, sep='\t', float_precision='round_trip').set_index('gene_id')
    samples = sum(GROUPS.values(), [])
    if not counts.index.is_unique or not (counts.gene_length > 0).all():
        raise ValueError('Counts require unique gene IDs and positive lengths.')
    if not np.isfinite(counts[samples]).all().all() or not (counts[samples] >= 0).all().all():
        raise ValueError('Counts must be finite and nonnegative.')
    rpk = counts[samples].divide(counts.gene_length / 1000, axis=0)
    tpm = rpk.divide(rpk.sum(axis=0), axis=1) * 1e6
    if not np.allclose(tpm.sum(), 1e6, rtol=0, atol=1e-7):
        raise ValueError('TPM column totals are inconsistent.')
    means = pd.DataFrame({group: tpm[ids].mean(axis=1) for group, ids in GROUPS.items()})
    baseline = pd.read_csv(args.states, sep='\t', float_precision='round_trip')
    baseline['mean_TTX_TPM'] = means.loc[baseline.gene_id, 'TTX'].to_numpy()
    baseline['log2_TTX_TPM_plus1'] = np.log2(baseline.mean_TTX_TPM + 1)
    args.output.mkdir(parents=True, exist_ok=True)
    means.to_csv(args.output/'Bic_mean_TPM.tsv', sep='\t')
    baseline.to_csv(args.output/'baseline_TPM.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()

"""Assemble IAP-family RNA, B1 matching and LAD-distribution panel tables."""
from pathlib import Path
import argparse
import pandas as pd

TARGETS = ['IAPEY4_I', 'IAPA_MM', 'IAPEZI']
GROUPS = {'TTX': ['L9', 'L12'], 'Bicuculline': ['L11', 'L14', 'L16'],
          'Bicuculline + nimodipine': ['L10', 'L13', 'L15']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['family', 'b1', 'states', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    effects = pd.read_csv(args.family/'CROSS_PARADIGM_FAMILY_EFFECTS.tsv', sep='\t')
    columns = ['Name', 'KCl_Mann_6h', 'Bic_TTX', 'new_measured', 'KCl_6h_measured']
    effects = effects[columns].copy()
    for column in ['new_measured', 'KCl_6h_measured']:
        if effects[column].dtype != bool:
            raise ValueError('Expected boolean family eligibility: ' + column)
        effects[column] = effects[column].map({True: 'TRUE', False: 'FALSE'})
    effects.to_csv(args.output/'TE_family_effects.tsv', sep='\t', index=False)
    rna = pd.read_csv(args.family/'IAP_REPLICATE_VALUES.tsv', sep='\t')
    b1 = pd.read_csv(args.b1/'IAP_B1_per_IP.tsv', sep='\t')
    conditions = {sample: group for group, samples in GROUPS.items() for sample in samples}
    rows = []
    for target in TARGETS:
        for row in rna[rna.Name == target].itertuples():
            condition = ('TTX' if row.sample_id.startswith('Basal_') else
                         'Bicuculline + nimodipine' if row.sample_id.startswith('Bic_nimo_') else
                         'Bicuculline' if row.sample_id.startswith('Bic_') else None)
            if condition is None:
                raise ValueError('Unknown RNA condition: ' + row.sample_id)
            rows.append(dict(target=target, assay='RNA', region_type='SalmonTE family',
                             sample=row.sample_id, condition=condition, value=row.CPM,
                             units='input-pair CPM'))
        for kind, field in [('IAP-containing regions', 'target_matched_subset'),
                            ('Matched non-IAP regions', 'matched_control')]:
            for row in b1[b1.target == target].itertuples():
                rows.append(dict(target=target, assay='B1', region_type=kind,
                                 sample=row.sample, condition=conditions[row.sample],
                                 value=getattr(row, field), units='log2(B1/Input)'))
    frame = pd.DataFrame(rows)
    keys = ['target', 'assay', 'region_type', 'sample']
    if frame.duplicated(keys).any():
        raise ValueError('Duplicate per-sample IAP panel rows.')
    frame.to_csv(args.output/'IAP_display_values.tsv', sep='\t', index=False)
    for name in ['IAP_state_distribution.tsv', 'IAP_enrichment_tests.tsv']:
        table = pd.read_csv(args.states/name, sep='\t')
        table.to_csv(args.output/name, sep='\t', index=False)


if __name__ == '__main__':
    main()

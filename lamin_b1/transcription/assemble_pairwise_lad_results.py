"""Join pairwise LAD tests to the corresponding RNA and GO gene sets."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

CONTRASTS = ['Bic_vs_TTX', 'Bic_Nimo_vs_Bic', 'Bic_Nimo_vs_TTX']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairwise', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = args.pairwise / 'tables'
    tables = args.output / 'tables'
    go_dir = args.output / 'GO_B1'
    tables.mkdir(parents=True, exist_ok=True)
    go_dir.mkdir(parents=True, exist_ok=True)
    go = pd.read_csv(source / '64_gene_GO_links.tsv', sep='\t')
    rna = pd.read_csv(source / '64_gene_RNA_links.tsv', sep='\t')
    all_go, all_rna, all_tests = [], [], []
    for contrast in CONTRASTS:
        tests = pd.read_csv(source / f'64_{contrast}_condition_depth_tests.tsv', sep='\t').rename(
            columns={'row_index': 'pair_region_index', 'delta_z': 'B1_delta_z', 'P': 'B1_P', 'q': 'B1_q'})
        tests['contrast'] = contrast
        tests['scheme'] = 'condition_depth'
        tests['B1_group'] = np.where(tests.B1_q < .05,
                                    np.where(tests.B1_delta_z < 0, 'Lower B1', 'Higher B1'), 'Not detected')
        all_tests.append(tests)
        all_go.append(go[go.contrast == contrast].merge(tests, on=['contrast', 'pair_region_index'], validate='many_to_one'))
        all_rna.append(rna[rna.contrast == contrast].merge(tests, on=['contrast', 'pair_region_index'], validate='many_to_one'))
    pd.concat(all_go).to_csv(go_dir / '64_GO_membership.tsv', sep='\t', index=False)
    pd.concat(all_rna).to_csv(tables / 'pairwise_gene_RNA.tsv', sep='\t', index=False)
    pd.concat(all_tests).to_csv(tables / 'pairwise_tests.tsv', sep='\t', index=False)
    print('Wrote pairwise LAD tests with RNA and GO gene links.')


if __name__ == '__main__':
    main()

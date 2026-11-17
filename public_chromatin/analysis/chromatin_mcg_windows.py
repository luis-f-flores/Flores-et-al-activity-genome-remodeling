"""Combine original-grid CpG counts and paired DNMT3A/Input read counts."""
from pathlib import Path
import argparse
from io import StringIO
import numpy as np
import pandas as pd


def calculate(inputs, output, methylation=None):
    inputs, output = Path(inputs), Path(output)
    methylation = Path(methylation) if methylation else inputs
    grid = pd.read_csv(inputs / 'mm9_5kb_windows.tsv', sep='\t')
    valid = np.ones(len(grid), dtype=bool)
    for sample in ['A0', 'B0', 'C0', '1wk', '2wk']:
        data = np.load(methylation / f'{sample}_mm9_windows.npz')
        with np.errstate(divide='ignore', invalid='ignore'):
            grid[sample] = 100 * data['meth'] / data['depth']
        valid &= (data['sites'] >= 3) & np.isfinite(grid[sample])
    grid['Baseline'] = grid[['A0', 'B0', 'C0']].mean(axis=1)
    for rep in [1, 2, 3]:
        chip = np.load(inputs / f'DNMT3A_ChIP_{rep}_mm9_5kb.npz')
        control = np.load(inputs / f'Input_{rep}_mm9_5kb.npz')
        with np.errstate(divide='ignore', invalid='ignore'):
            grid[f'DNMT3A_rep{rep}'] = np.log2((chip['counts'] / int(chip['total_unique_reads'])) /
                                             (control['counts'] / int(control['total_unique_reads'])))
        # Round-trip each enrichment column through CSV before stable rank grouping.
        column = f'DNMT3A_rep{rep}'
        grid[column] = pd.read_csv(StringIO(grid[[column]].to_csv(index=False)))[column]
        valid &= np.isfinite(grid[f'DNMT3A_rep{rep}'])
    grid['eligible'] = valid
    # Round-trip the assembled table through CSV to preserve its ranking precision.
    grid = pd.read_csv(StringIO(grid.to_csv(index=False)))
    output.mkdir(parents=True, exist_ok=True)
    grid.to_csv(output / 'mCG_DNMT3A_windows.tsv.gz', sep='\t', index=False)
    print(f'{valid.sum():,} jointly supported windows')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--methylation', type=Path)
    args = parser.parse_args()
    calculate(args.inputs, args.output, args.methylation)
